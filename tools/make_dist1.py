#!/usr/bin/env python3
"""Hand-craft a SINGLE-RANK (dist_1) FESOM2 partition for a mesh.

METIS cannot emit a 1-PE partition. The FESOM3 model does not need one (at 1 rank it
builds its own identity partition, synthesize_1rank); dist_1 is for tools/mesh/check_dist.py
and the FESOM2 tools and 1-rank oracle runs (per-node/element accumulation order ==
FESOM3's global order; HANDOFF "Geometry byte-gate" + LESSONS L8). For np=1 the partition
is trivial: every node/elem/edge is owned, there is no halo, and the communication
structure is empty. The FESOM2 part mirrors FESOM2's `save_dist_mesh` np=1 output,
validated against the proven pi dist_1 template; format.out and the empty appended
sections (below) are always written.

Format (from the FESOM2 reader oce_mesh.F90 read_mesh):
  rpart.out          : line1 = npes(1); line2 = count = nod2D; then the contiguous
                       node->global identity map 1..nod2D (one per line).
  my_list00000.out   : n(0); myDim_nod2D, eDim_nod2D(0), <node list>;
                       myDim_elem2D, eDim_elem2D(0), eXDim_elem2D(0), <elem list>;
                       myDim_edge2D, eDim_edge2D(0), <edge list>.
                       Each SCALAR on its own line; each LIST list-directed (any
                       wrapping is safe — each read() starts a fresh record, and a
                       list read spans records until it has all its values).
  com_info00000.out  : rank (0), then empty np=1 comm structs (rPEnum=sPEnum=0, blank
                       zero-size arrays, rptr=sptr=1) — MESH-INDEPENDENT (the pi
                       template is checked to be rank + 3 of them).
Hybrid meshes (fesom_ini's dist_N/format.out, version 1): the marker format.out is
written and the appended sections are empty — my_list gets eXDim_nod2D(0) and its
(blank) list, com_info gets two more empty structs (com_nod2D_full, com_edge2D).

Usage: tools/make_dist1.py <mesh_dir> [pi_template_dir]
"""
import os, sys

PI_TEMPLATE = "/home/a/a270088/port2/fesom2/tests/data/MESHES/pi/dist_1"
WRAP = 30  # values per line in the identity lists (wrap-safe; avoids MB-long lines)
FORMAT_VERSION = 1  # dist_N/format.out (appended sections present)
EMPTY_COM = "0\n\n1\n\n0\n\n1\n\n"  # one empty comm struct (rPEnum .. slist)


def first_int(path):
    with open(path) as f:
        return int(f.read().split()[0])


def wrapped(n, width=WRAP):
    """Yield lines of '1 2 ... n', `width` values each (ends exactly at value n)."""
    out, buf = [], []
    for i in range(1, n + 1):
        buf.append(str(i))
        if len(buf) == width:
            out.append(" ".join(buf)); buf = []
    if buf:
        out.append(" ".join(buf))
    return out


def build_my_list(nod2D, elem2D, edge2D):
    lines = ["0",
             str(nod2D), "0", *wrapped(nod2D),
             str(elem2D), "0", "0", *wrapped(elem2D),
             str(edge2D), "0", *wrapped(edge2D),
             "0", ""]                       # eXDim_nod2D, its (empty) list
    return "\n".join(lines) + "\n"


def build_rpart(nod2D):
    return "1\n" + str(nod2D) + "\n" + "\n".join(str(i) for i in range(1, nod2D + 1)) + "\n"


def validate_against_pi(pi_dir):
    """Re-parse the proven pi my_list and assert our format model is correct."""
    toks = open(os.path.join(pi_dir, "my_list00000.out")).read().split()
    it = iter(toks)
    nxt = lambda: int(next(it))
    assert nxt() == 0, "my_list line1 (n) should be 0"
    myd_n, ed_n = nxt(), nxt(); assert ed_n == 0
    nod = [nxt() for _ in range(myd_n)]
    assert nod == list(range(1, myd_n + 1)), "node list not identity"
    myd_e, ed_e, exd_e = nxt(), nxt(), nxt(); assert ed_e == 0 and exd_e == 0
    elem = [nxt() for _ in range(myd_e)]
    assert elem == list(range(1, myd_e + 1)), "elem list not identity"
    myd_g, ed_g = nxt(), nxt(); assert ed_g == 0
    edge = [nxt() for _ in range(myd_g)]
    assert edge == list(range(1, myd_g + 1)), "edge list not identity"
    try:
        next(it); raise AssertionError("trailing tokens in pi my_list")
    except StopIteration:
        pass
    print(f"  [validate] pi my_list OK: nod2D={myd_n} elem2D={myd_e} edge2D={myd_g}")
    com = open(os.path.join(pi_dir, "com_info00000.out")).read()
    assert com == "0\n" + 3 * EMPTY_COM, "pi com_info is not rank + 3 empty structs"
    print("  [validate] pi com_info OK: rank + 3 empty comm structs")
    return myd_n, myd_e, myd_g


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mesh_dir = sys.argv[1].rstrip("/")
    pi_dir = sys.argv[2] if len(sys.argv) > 2 else PI_TEMPLATE

    print(f"make_dist1: target mesh = {mesh_dir}")
    validate_against_pi(pi_dir)

    nod2D = first_int(os.path.join(mesh_dir, "nod2d.out"))
    elem2D = first_int(os.path.join(mesh_dir, "elem2d.out"))
    edge2D = first_int(os.path.join(mesh_dir, "edgenum.out"))
    print(f"  mesh dims: nod2D={nod2D} elem2D={elem2D} edge2D={edge2D}")

    out = os.path.join(mesh_dir, "dist_1")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "rpart.out"), "w") as f:
        f.write(build_rpart(nod2D))
    with open(os.path.join(out, "my_list00000.out"), "w") as f:
        f.write(build_my_list(nod2D, elem2D, edge2D))
    # com_info is mesh-independent for np=1: rank 0, then the empty com_nod2D,
    # com_elem2D, com_elem2D_full (FESOM2) and com_nod2D_full, com_edge2D (appended)
    with open(os.path.join(out, "com_info00000.out"), "w") as f:
        f.write("0\n" + 5 * EMPTY_COM)
    with open(os.path.join(out, "format.out"), "w") as f:
        f.write(f"{FORMAT_VERSION}\n")
    print(f"  wrote {out}/{{rpart.out, my_list00000.out, com_info00000.out, format.out}}")


if __name__ == "__main__":
    main()
