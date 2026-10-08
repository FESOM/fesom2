#!/usr/bin/env python3
"""Compare two horizontal-advection dumps (FESOM2 oracle vs FESOM3) field-by-field
for the M1.1 operator byte-gate.

Both files use the self-describing format written by FESOM2
src/fesom_advhor_dump.F90 and FESOM3 src/infra/mod_advhor_dump.F90:

  char8  "FADVHDMP" | int32 nod2D,elem2D,edge2D,nl
  records: char24 name | int32 dtype(0=real64,1=int32) | int32 d1,d2 | d1*d2 values

The gated outputs are adv_flux_hor_{upw1,muscl} + del_ttf_advhoriz_{upw1,muscl}
(M1.1 horizontal), adv_flux_ver_{upw1,qr4c} + del_ttf_advvert_{upw1,qr4c} (M1.2
vertical; the node tendency, oce_ale_tracer.F90:240), and the M1.3 FCT (Zalesak)
path: fct_LO, the clipped adv_flux_{hor,ver}_fct, the limiter internals (fct_ttf_
max/min bounds, fct_plus/minus factors, pre-clip *_fct_ho fluxes) and the final
del_ttf_{advhoriz,advvert}_fct, plus the standalone adv_flux_hor_mfct and the
node thicknesses hnode/hnode_new. The M1.4 assembled-step fields gate the REAL
driver: valuesAB (init_tracers_AB AB(2) interpolation), del_ttf_{advhoriz,advvert,}
_step (do_oce_adv_tra FCT path + the del_ttf accumulation) and ..._step{non,_non}
(the non-FCT do_zero_flux path + per-tracer order knobs). The intermediates (ttf,
ttfAB, vel, helem, edge_up_dn_tri, tr_xy, wvel,
zbar_3d_n, Z_3d_n, area) are dumped too so any divergence is localised.

Usage:  advhor_diff.py <fesom2.bin> <fesom3.bin> [--signal 0.0]
Exit 0 iff every common field matches within the signal threshold (default 0 ==
byte-identical) and shapes/dims agree.
"""
import sys
import struct
import numpy as np

MAGIC = b"FADVHDMP"


def read_dump(path):
    with open(path, "rb") as f:
        buf = f.read()
    if buf[:8] != MAGIC:
        raise ValueError(f"{path}: bad magic {buf[:8]!r}")
    off = 8
    nod2D, elem2D, edge2D, nl = struct.unpack_from("<4i", buf, off)
    off += 16
    dims = dict(nod2D=nod2D, elem2D=elem2D, edge2D=edge2D, nl=nl)
    fields = {}
    n = len(buf)
    while off < n:
        name = buf[off:off + 24].decode("ascii", "replace").strip()
        off += 24
        dtype, d1, d2 = struct.unpack_from("<3i", buf, off)
        off += 12
        count = d1 * d2
        if dtype == 0:
            data = np.frombuffer(buf, dtype="<f8", count=count, offset=off)
            off += 8 * count
        elif dtype == 1:
            data = np.frombuffer(buf, dtype="<i4", count=count, offset=off)
            off += 4 * count
        else:
            raise ValueError(f"{path}: unknown dtype {dtype} for {name}")
        fields[name] = (dtype, d1, d2, data)
    return dims, fields


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    signal = 0.0
    for a in sys.argv[1:]:
        if a.startswith("--signal"):
            signal = float(a.split("=", 1)[1]) if "=" in a else 0.0
    if len(args) < 2:
        print(__doc__)
        sys.exit(2)
    p2, p3 = args[0], args[1]
    d2, f2 = read_dump(p2)
    d3, f3 = read_dump(p3)

    print(f"FESOM2: {p2}")
    print(f"FESOM3: {p3}")
    print(f"dims  FESOM2={d2}")
    print(f"dims  FESOM3={d3}")
    ok = True
    if d2 != d3:
        print("FAIL: mesh dims differ")
        ok = False

    names = list(f2.keys())
    extra3 = [k for k in f3 if k not in f2]
    print(f"\n{'field':24s} {'shape':>14s} {'dtype':>6s} {'max|delta|':>16s}  status")
    print("-" * 74)
    for name in names:
        if name not in f3:
            print(f"{name:24s} {'-':>14s} {'-':>6s} {'MISSING in F3':>16s}  FAIL")
            ok = False
            continue
        t2, a1, a2, v2 = f2[name]
        t3, b1, b2, v3 = f3[name]
        shape = f"{a1}x{a2}"
        if (t2, a1, a2) != (t3, b1, b2):
            print(f"{name:24s} {shape:>14s} {t2:>6d} shape/dtype F3={b1}x{b2}  FAIL")
            ok = False
            continue
        if v2.shape != v3.shape:
            print(f"{name:24s} {shape:>14s} {t2:>6d} {'len mismatch':>16s}  FAIL")
            ok = False
            continue
        if t2 == 0:
            d = np.abs(v2 - v3)
            mx = float(d.max()) if d.size else 0.0
            status = "PASS" if mx <= signal else "FAIL"
        else:
            nmis = int((v2 != v3).sum())
            mx = float(nmis)
            status = "PASS" if nmis == 0 else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"{name:24s} {shape:>14s} {t2:>6d} {mx:>16.6e}  {status}")

    if extra3:
        print(f"\n(only in FESOM3, ignored: {extra3})")

    print("\n" + ("TRACER-ADVECTION BYTE-GATE (M1.1 horiz + M1.2 vert + M1.3 FCT + M1.4 step): PASS (max|delta|=0)" if ok else
                  "TRACER-ADVECTION BYTE-GATE (M1.1 horiz + M1.2 vert + M1.3 FCT + M1.4 step): FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
