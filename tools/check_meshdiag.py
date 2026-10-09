#!/usr/bin/env python3
"""Check the hybrid-mesh layout of fesom.mesh.diag.zarr (hybrid-mesh plan
docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14). Standard library only (the ctest runs it with
the system python3, which has no numpy); the store is uncompressed Zarr v2 (mod_io_meshdiag).

usage: check_meshdiag.py <mesh_dir> <work_dir> <fesom_meshdiagdump> <np[,np...]> [launcher ...]

Runs `launcher -np <np> fesom_meshdiagdump` for every np (FESOM3_MESH_DIR=<mesh_dir>,
FESOM3_MESHDIAG_OUT=<work_dir>/md_np<np>.zarr), then checks every store:
  M1 face_nodes / face_edges / face_links / gradient_sca_x / gradient_sca_y have shape (4, elem2D),
     dims (n4, elem); fill_value -1 / -1 / -999 (the UGRID _FillValue), start_index 1, the
     topology names face_edges / face_links
  M2 face_nodes == the elements of elem2d.out (vertex sets; 3 or 4 columns, a triangle repeats
     its 1st id or has 3 columns): a triangle's 4th entry is -1, a quad has 4 distinct ids
  M3 face_edges: edge k joins face nodes k and k+1 (cyclic over the element's nv), checked
     against edge_nodes; a triangle's 4th entry is -1
  M4 face_links: across edge k the neighbour is the other element of edge_face_links (-999 on a
     boundary edge), and the link is mutual; a triangle's 4th entry is -999
  M5 gradient_sca_x/y: a triangle's 4th coefficient is 0; the coefficients sum to 0 (gradient of
     a constant) to roundoff
and compares every store with the first one (the store is partition-independent): integer arrays
identical (except nod_part / elem_part), real arrays to 1e-12 relative.
Exit 0 iff all checks pass; a missing mesh prints "ERROR: mesh not found" and exits 1 (the
meshdiag_* ctests are registered DISABLED at configure time when the mesh is absent).
"""
import json
import os
import shutil
import struct
import subprocess
import sys

DT = {"<i4": ("i", 4), "<f8": ("d", 8), "<f4": ("f", 4)}


def read_zarr(store, name):
    """Return (meta, attrs, data) of an uncompressed C-order Zarr v2 array; data is a flat list
    (C order, i.e. for shape (R, N) the entry (r, n) is data[r*N + n])."""
    d = os.path.join(store, name)
    meta = json.load(open(os.path.join(d, ".zarray")))
    attrs = json.load(open(os.path.join(d, ".zattrs"))) if os.path.exists(os.path.join(d, ".zattrs")) else {}
    if meta["compressor"] is not None or meta["order"] != "C":
        raise SystemExit(f"FAIL {name}: compressed or Fortran-order chunks are not supported")
    code, size = DT[meta["dtype"]]
    shape, chunks = meta["shape"], meta["chunks"]
    sep = meta.get("dimension_separator", ".")
    fill = meta["fill_value"]
    if len(shape) == 1:
        out = [fill] * shape[0]
        for c in range((shape[0] + chunks[0] - 1) // chunks[0]):
            p = os.path.join(d, str(c))
            if not os.path.exists(p):
                continue
            vals = struct.unpack("<%d%s" % (chunks[0], code), open(p, "rb").read())
            lo = c * chunks[0]
            n = min(chunks[0], shape[0] - lo)
            out[lo:lo + n] = vals[:n]
        return meta, attrs, out
    R, N = shape
    if chunks[0] != R:
        raise SystemExit(f"FAIL {name}: chunked along the first dim")
    C = chunks[1]
    out = [fill] * (R * N)
    for c in range((N + C - 1) // C):
        p = os.path.join(d, "0" + sep + str(c))
        if not os.path.exists(p):
            continue
        vals = struct.unpack("<%d%s" % (R * C, code), open(p, "rb").read())
        lo = c * C
        n = min(C, N - lo)
        for r in range(R):
            out[r * N + lo:r * N + lo + n] = vals[r * C:r * C + n]
    return meta, attrs, out


def read_elem2d(mesh_dir):
    with open(os.path.join(mesh_dir, "elem2d.out")) as f:
        ne = int(f.readline().split()[0])
        elems = []
        for _ in range(ne):
            ids = [int(x) for x in f.readline().split()]
            if len(ids) == 4 and ids[3] == ids[0]:
                ids = ids[:3]
            elems.append(ids)
    return elems


class Checker:
    def __init__(self):
        self.nfail = 0

    def ok(self, cond, msg):
        if not cond:
            self.nfail += 1
            print("  FAIL " + msg)
        return cond


def check_store(store, elems, chk):
    ne = len(elems)
    a = {}
    for nm, fill in (("face_nodes", -1), ("face_edges", -1), ("face_links", -999),
                     ("gradient_sca_x", None), ("gradient_sca_y", None)):
        meta, attrs, data = read_zarr(store, nm)
        a[nm] = data
        chk.ok(meta["shape"] == [4, ne], f"M1 {nm} shape {meta['shape']} != [4, {ne}]")
        chk.ok(attrs.get("_ARRAY_DIMENSIONS") == ["n4", "elem"], f"M1 {nm} dims {attrs.get('_ARRAY_DIMENSIONS')}")
        if fill is not None:
            chk.ok(meta["fill_value"] == fill, f"M1 {nm} fill_value {meta['fill_value']} != {fill}")
            chk.ok(attrs.get("start_index") == 1, f"M1 {nm} start_index {attrs.get('start_index')}")
    _, tattrs, _ = read_zarr(store, "fesom_mesh")
    chk.ok(tattrs.get("face_edge_connectivity") == "face_edges"
           and tattrs.get("face_face_connectivity") == "face_links", "M1 topology face_edges/face_links names")
    meta, _, en = read_zarr(store, "edge_nodes")
    ned = meta["shape"][1]
    _, _, efl = read_zarr(store, "edge_face_links")
    fn, fe, fl = a["face_nodes"], a["face_edges"], a["face_links"]
    gx, gy = a["gradient_sca_x"], a["gradient_sca_y"]
    col = lambda arr, k, e: arr[k * ne + e]
    bad = {k: 0 for k in ("M2", "M3", "M4", "M5")}
    ntri = nquad = nlink = nbnd = 0
    gmax = 0.0
    for e in range(ne):
        ids = elems[e]
        nv = len(ids)
        if nv == 3:
            ntri += 1
        else:
            nquad += 1
        f = [col(fn, k, e) for k in range(4)]
        if sorted(f[:nv]) != sorted(ids) or (nv == 3 and f[3] != -1) or len(set(f[:nv])) != nv:
            bad["M2"] += 1
        for k in range(4):
            ed, nb = col(fe, k, e), col(fl, k, e)
            if k >= nv:
                if ed != -1:
                    bad["M3"] += 1
                if nb != -999:
                    bad["M4"] += 1
                continue
            if not 1 <= ed <= ned or {en[ed - 1], en[ned + ed - 1]} != {f[k], f[(k + 1) % nv]}:
                bad["M3"] += 1
                continue
            t1, t2 = efl[ed - 1], efl[ned + ed - 1]
            if nb == -999:
                nbnd += 1
                if t1 != e + 1 or t2 != -999:
                    bad["M4"] += 1
            else:
                nlink += 1
                if {t1, t2} != {e + 1, nb} or e + 1 not in [col(fl, j, nb - 1) for j in range(4)]:
                    bad["M4"] += 1
        cx = [col(gx, k, e) for k in range(4)]
        cy = [col(gy, k, e) for k in range(4)]
        if nv == 3 and (cx[3] != 0.0 or cy[3] != 0.0):
            bad["M5"] += 1
        for c in (cx, cy):
            r = abs(sum(c)) / max(sum(abs(x) for x in c), 1e-300)
            gmax = max(gmax, r)
            if r > 1e-12:
                bad["M5"] += 1
    print(f"  {os.path.basename(store)}: {ne} faces ({ntri} triangles, {nquad} quads), {nlink} face links, "
          f"{nbnd} boundary sides; max |sum of gradient coefficients|/sum|c| = {gmax:.2e}")
    for k, v in bad.items():
        chk.ok(v == 0, f"{k}: {v} bad faces in {store}")
    if ntri:
        chk.ok(all(col(fn, 3, e) == -1 for e in range(ne) if len(elems[e]) == 3), "M2 triangle 4th == -1")


def compare_stores(s0, s1, chk):
    names = sorted(n for n in os.listdir(s0) if os.path.isdir(os.path.join(s0, n)))
    for nm in names:
        if nm in ("nod_part", "elem_part"):
            continue
        m0, _, d0 = read_zarr(s0, nm)
        m1, _, d1 = read_zarr(s1, nm)
        if not chk.ok(m0["shape"] == m1["shape"], f"compare {nm}: shapes differ"):
            continue
        if m0["dtype"] == "<i4":
            nd = sum(1 for x, y in zip(d0, d1) if x != y)
            chk.ok(nd == 0, f"compare {nm}: {nd} integer entries differ ({s1} vs {s0})")
        else:
            ref = max((abs(x) for x in d0), default=0.0)
            d = max((abs(x - y) for x, y in zip(d0, d1)), default=0.0) / max(ref, 1e-300)
            chk.ok(d <= 1e-12, f"compare {nm}: max|d|/max|ref| = {d:.2e} > 1e-12")
    print(f"  {os.path.basename(s1)} == {os.path.basename(s0)} ({len(names) - 2} arrays compared)")


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    mesh_dir, work, exe, nps = sys.argv[1:5]
    launcher = sys.argv[5:]
    if not os.path.exists(os.path.join(mesh_dir, "elem2d.out")):
        print(f"ERROR: mesh not found: no {mesh_dir}/elem2d.out")
        return 1
    nps = [int(x) for x in nps.split(",")]
    elems = read_elem2d(mesh_dir)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    chk = Checker()
    stores = []
    for n in nps:
        st = os.path.join(work, f"md_np{n}.zarr")
        env = dict(os.environ, FESOM3_MESH_DIR=mesh_dir, FESOM3_MESHDIAG_OUT=st)
        cmd = (launcher + ["-np", str(n)] if launcher else []) + [exe]
        r = subprocess.run(cmd, env=env, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True)
        if not chk.ok(r.returncode == 0 and "MESHDIAGDUMP OK" in r.stdout, f"np {n}: driver failed"):
            print(r.stdout[-3000:])
            continue
        print(f"np {n}:")
        check_store(st, elems, chk)
        stores.append(st)
    for st in stores[1:]:
        compare_stores(stores[0], st, chk)
    if chk.nfail == 0:
        shutil.rmtree(work, ignore_errors=True)
        print("MESHDIAG CHECK PASS")
        return 0
    print(f"MESHDIAG CHECK FAILED ({chk.nfail}); stores kept in {work}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
