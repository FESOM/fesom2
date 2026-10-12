#!/usr/bin/env python3
"""dump_diff.py — compare FESOM3 vs FESOM2 gid-keyed substep dumps.

Reads the little-endian stream record format written by src/infra/mod_dump.F90
and the FESOM2 oracle shim (fesom_dump_shim.F90):

    int32 step | int32 substep_id | int32 probe_gid | int32 nlevels |
    char[24] field_name | float64 values[nlevels]

Groups records by (step, substep_id, gid, field_name) and reports the per-record
max|delta|, the FIRST diverging substep (Fortran->Fortran gate is max|delta|=0),
and a magnitude histogram separating real divergence from FP noise.

Usage:
    dump_diff.py A B [--threshold T] [--glob]   # compare files (or prefixes with --glob)
    dump_diff.py A B --ignore-substep=2          # skip a substep id (repeatable; e.g. the
                                                 # M2-dead SW_AB=2 that FESOM3 does not emit)
    dump_diff.py A B --ignore-step               # compare the SAME physical state reached via a
                                                 # different per-segment step count (restart gate:
                                                 # straight step N vs split seg-2 step N-K). Each
                                                 # side must be windowed to ONE step; records re-key
                                                 # on (substep, gid, name) so the labels align.
    dump_diff.py --selftest                      # self-check (no external data)

Exit code 0 = match (all |delta| <= threshold), 1 = divergence (for ctest).
"""
import sys, os, struct, glob, math

HDR = struct.Struct("<4i24s")        # step, substep, gid, nlevels, name[24]
SUBSTEP_NAMES = {
    0: "INIT", 1: "PRESSURE_BV", 2: "SW_AB", 3: "PGF", 4: "MIXING",
    5: "VEL_RHS", 6: "VISC_FILTER", 7: "IMPL_VISC", 8: "SSH_RHS", 9: "SSH_SOLVE",
    10: "UPDATE_VEL", 11: "HBAR", 12: "ETA_N", 13: "ALE", 14: "GM_BOLUS",
    15: "TRACERS", 16: "THICKNESS",
}


def parse_file(path):
    """Return dict {(step, substep, gid, name): tuple(values)}."""
    out = {}
    with open(path, "rb") as f:
        data = f.read()
    off, n = 0, len(data)
    while off < n:
        if off + HDR.size > n:
            raise ValueError(f"{path}: truncated header at byte {off}")
        step, substep, gid, nlev, raw = HDR.unpack_from(data, off)
        off += HDR.size
        name = raw.split(b"\x00")[0].decode("ascii", "replace").strip()
        need = 8 * nlev
        if off + need > n:
            raise ValueError(f"{path}: truncated values at byte {off}")
        vals = struct.unpack_from("<%dd" % nlev, data, off)
        off += need
        out[(step, substep, gid, name)] = vals
    return out


def parse_tree(arg, use_glob):
    """Parse one file or all per-rank files matching a prefix (merge ranks)."""
    files = sorted(glob.glob(arg + ".*")) if use_glob else [arg]
    if not files:
        raise FileNotFoundError(f"no files for '{arg}'" + (".*" if use_glob else ""))
    merged = {}
    for p in files:
        merged.update(parse_file(p))
    return merged


def collapse_step(records, label):
    """Re-key a parsed dump so its step field is dropped (set to 0), for comparing the SAME
    physical state reached via a different per-segment step count — the restart gate compares the
    straight run's final step N against the split run's seg-2 final step N-K. Each side MUST be
    windowed to a single step (FESOM_DUMP_MINSTEP==FESOM_DUMP_MAXSTEP); the single-step guard makes
    masking impossible (two distinct steps cannot silently collapse onto one another)."""
    steps = sorted({k[0] for k in records})
    if len(steps) != 1:
        raise ValueError(
            f"--ignore-step: '{label}' carries {len(steps)} distinct step(s) {steps}; window the "
            "dump to ONE step (FESOM_DUMP_MINSTEP==FESOM_DUMP_MAXSTEP) so the alignment is unambiguous")
    return {(0, ss, gid, name): v for (st, ss, gid, name), v in records.items()}


def max_abs_diff(va, vb):
    if len(va) != len(vb):
        return math.inf
    return max((abs(a - b) for a, b in zip(va, vb)), default=0.0)


def max_abs(va):
    return max((abs(a) for a in va), default=0.0)


def compare(a, b, threshold, rel_floor=0.0):
    keys = sorted(set(a) | set(b))
    only_a = [k for k in keys if k not in b]
    only_b = [k for k in keys if k not in a]
    diffs = []  # (step, substep, gid, name, maxabsdiff)
    eff = {}    # per-key effective threshold = max(abs threshold, rel_floor*|ref|)
    for k in keys:
        if k in a and k in b:
            diffs.append((k[0], k[1], k[2], k[3], max_abs_diff(a[k], b[k])))
            eff[k] = max(threshold, rel_floor * max_abs(a[k])) if rel_floor > 0.0 else threshold
    # first diverging (step, substep) above the (possibly relative) effective threshold.
    # rel_floor>0 (np>1 restart with FESOM3_SHARED_SYNC=0 only) admits the redundant-element ~1-ULP
    # roundoff of unsynced shared-element copies (see zarr_diff.output_cmp doc) while still catching
    # any real regression (>= ~1e-4). With the sync (the default) np>1 restarts are exact: floor 0.
    first = None
    for step, substep, gid, name, d in sorted(diffs):
        if d > eff[(step, substep, gid, name)]:
            first = (step, substep)
            break
    # histogram of magnitudes
    buckets = {"=0": 0, "<=1e-12": 0, "<=1e-6": 0, "<=1e-3": 0, ">1e-3": 0}
    worst = 0.0
    for *_, d in diffs:
        worst = max(worst, d)
        if d == 0.0:
            buckets["=0"] += 1
        elif d <= 1e-12:
            buckets["<=1e-12"] += 1
        elif d <= 1e-6:
            buckets["<=1e-6"] += 1
        elif d <= 1e-3:
            buckets["<=1e-3"] += 1
        else:
            buckets[">1e-3"] += 1
    return diffs, only_a, only_b, first, buckets, worst


def report(diffs, only_a, only_b, first, buckets, worst, threshold):
    print(f"records compared: {len(diffs)}   worst |delta| = {worst:.3e}   "
          f"threshold = {threshold:.3e}")
    print("  magnitude histogram: " +
          "  ".join(f"{k}:{v}" for k, v in buckets.items()))
    if only_a:
        print(f"  WARNING: {len(only_a)} key(s) only in A")
    if only_b:
        print(f"  WARNING: {len(only_b)} key(s) only in B")
    if first is None and not only_a and not only_b:
        print("MATCH: all records within threshold.")
        return 0
    if first is not None:
        step, substep = first
        sname = SUBSTEP_NAMES.get(substep, f"substep#{substep}")
        print(f"DIVERGENCE: first at step {step}, substep {substep} ({sname}).")
        worst_here = [d for d in diffs if d[0] == step and d[1] == substep]
        for s, ss, gid, name, d in sorted(worst_here, key=lambda x: -x[4])[:5]:
            if d > threshold:
                print(f"    gid={gid:>7} field={name:<16} |delta|={d:.3e}")
    return 1


def selftest():
    import tempfile
    # Build a synthetic dump, write it, perturb a copy, verify detection.
    recs = []
    for step in (1, 2):
        for substep in (1, 3, 15):
            for gid in (1001, 2000):
                nlev = 3
                vals = [float(step * 100 + substep * 10 + gid + lv) for lv in range(nlev)]
                recs.append((step, substep, gid, "density", vals))

    def write(path, perturb_key=None, eps=1e-2):
        with open(path, "wb") as f:
            for step, substep, gid, name, vals in recs:
                v = list(vals)
                if perturb_key == (step, substep, gid):
                    v[1] += eps
                f.write(HDR.pack(step, substep, gid, len(v), name.encode()))
                f.write(struct.pack("<%dd" % len(v), *v))

    d = tempfile.mkdtemp(prefix="dumpdiff_")
    fa, fb = os.path.join(d, "a"), os.path.join(d, "b")
    ok = True

    # 1. identical files -> MATCH (exit 0)
    write(fa); write(fb)
    a, b = parse_file(fa), parse_file(fb)
    rc = report(*compare(a, b, 0.0), 0.0)
    print("[selftest] identical files ->", "PASS" if rc == 0 else "FAIL"); ok &= (rc == 0)

    # 2. perturbed copy -> DIVERGENCE flagged at the perturbed substep
    write(fb, perturb_key=(2, 3, 2000))
    res = compare(parse_file(fa), parse_file(fb), 0.0)
    rc = report(*res, 0.0)
    first = res[3]
    print("[selftest] perturbed copy ->",
          "PASS" if (rc == 1 and first == (2, 3)) else "FAIL")
    ok &= (rc == 1 and first == (2, 3))

    # 3. round-trip parse: values survive write->read exactly
    rt = all(parse_file(fa)[(s, ss, g, n)] == tuple(v)
             for s, ss, g, n, v in recs)
    print("[selftest] round-trip parse ->", "PASS" if rt else "FAIL"); ok &= rt

    # 4. --ignore-step: two single-step dumps with DIFFERENT step labels but identical state must
    #    MATCH after collapse (restart gate: straight step N vs split seg-2 step N-K).
    one_a = {(7, 1, 1001, "T"): (1.0, 2.0), (7, 15, 2000, "S"): (3.0,)}
    one_b = {(2, 1, 1001, "T"): (1.0, 2.0), (2, 15, 2000, "S"): (3.0,)}
    ca, cb = collapse_step(one_a, "A"), collapse_step(one_b, "B")
    rc = report(*compare(ca, cb, 0.0), 0.0)
    print("[selftest] ignore-step aligns N vs N-K ->", "PASS" if rc == 0 else "FAIL"); ok &= (rc == 0)
    #    a real divergence must still surface after collapse
    one_b_bad = dict(one_b); one_b_bad[(2, 15, 2000, "S")] = (3.5,)
    rc = report(*compare(collapse_step(one_a, "A"), collapse_step(one_b_bad, "B"), 0.0), 0.0)
    print("[selftest] ignore-step still detects divergence ->", "PASS" if rc == 1 else "FAIL"); ok &= (rc == 1)
    #    the single-step guard rejects a multi-step dump (masking is impossible)
    try:
        collapse_step({(1, 1, 1001, "T"): (1.0,), (2, 1, 1001, "T"): (1.0,)}, "multi")
        guard = False
    except ValueError:
        guard = True
    print("[selftest] ignore-step single-step guard ->", "PASS" if guard else "FAIL"); ok &= guard

    for p in (fa, fb):
        os.remove(p)
    os.rmdir(d)
    print("[selftest] ALL PASS" if ok else "[selftest] FAILURE")
    return 0 if ok else 1


def main(argv):
    if "--selftest" in argv:
        return selftest()
    args = [a for a in argv if not a.startswith("--")]
    threshold = 0.0
    rel_floor = 0.0
    use_glob = "--glob" in argv
    ignore_step = "--ignore-step" in argv
    ignore_substeps = set()
    for a in argv:
        if a.startswith("--threshold="):
            threshold = float(a.split("=", 1)[1])
        if a.startswith("--rel-floor="):
            rel_floor = float(a.split("=", 1)[1])
        if a.startswith("--ignore-substep="):
            ignore_substeps.add(int(a.split("=", 1)[1]))
    if len(args) != 2:
        print(__doc__)
        return 2
    a = parse_tree(args[0], use_glob)
    b = parse_tree(args[1], use_glob)
    if ignore_substeps:
        # key = (step, substep, gid, name); drop records of ignored substeps from both sides
        a = {k: v for k, v in a.items() if k[1] not in ignore_substeps}
        b = {k: v for k, v in b.items() if k[1] not in ignore_substeps}
        print("ignoring substep id(s): " + ", ".join(str(s) for s in sorted(ignore_substeps)))
    if ignore_step:
        a = collapse_step(a, args[0])
        b = collapse_step(b, args[1])
        print("ignoring step field (single-step alignment; compare same state, different step count)")
    if rel_floor > 0.0:
        print(f"np>1 rel-eps floor = {rel_floor:.0e} (redundant-element ~1-ULP roundoff, FESOM3_SHARED_SYNC=0)")
    return report(*compare(a, b, threshold, rel_floor), threshold)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
