#!/usr/bin/env python3
"""Split selected quads of a FESOM mesh into two triangles (makes a mixed mesh).

Selected are the quads whose centroid latitude (mean of the vertex latitudes) lies in
the band [LATMIN, LATMAX], plus every K-th of the remaining quads (element order).
A selected quad (1,2,3,4) is cut along its shorter diagonal (chord on the sphere):
1-3 -> (1,2,3), (1,3,4); 2-4 -> (1,2,4), (2,3,4). Both keep the clockwise vertex order.
At the coast the shorter diagonal can leave a triangle whose only neighbour is its twin
(a dead end, on which FESOM's setup stops); such a quad is cut along the other diagonal,
or left whole if both diagonals make a dead end.
The two triangles replace the quad in place and are written as 4 columns (4th = 1st).
Nodes and the depth file are unchanged; elements that are triangles already stay as they are.

Usage: split_quads.py <mesh_dir> <out_dir> [--band LATMIN LATMAX] [--every K]
  --band  : latitude band where every quad is split (default: none)
  --every : split every K-th quad outside the band (default 0: none)
"""
import argparse

import numpy as np

import mesh_io


def split(m, band=None, every=0):
    """Mesh with the selected quads split; returns (mesh, number of quads split,
    number cut along the longer diagonal, number left whole)."""
    nv = m.nnodes()
    quad = nv == 4
    sel = np.zeros(m.e2d, dtype=bool)
    if band is not None:
        clat = m.lat[m.elem].mean(axis=1)
        sel = quad & (clat >= band[0]) & (clat <= band[1])
    if every > 0:
        rest = np.nonzero(quad & ~sel)[0]
        sel[rest[every - 1::every]] = True

    p = mesh_io.xyz(m.lon, m.lat)[m.elem]
    d13 = np.linalg.norm(p[:, 0] - p[:, 2], axis=1)
    d24 = np.linalg.norm(p[:, 1] - p[:, 3], axis=1)
    # a triangle needs a neighbour across one of its two quad sides (sides 1-2, 2-3,
    # 3-4, 4-1); splitting a neighbour does not change which sides are shared
    sh = mesh_io.shared_sides(m)
    ok13 = (sh[:, 0] | sh[:, 1]) & (sh[:, 2] | sh[:, 3])
    ok24 = (sh[:, 0] | sh[:, 3]) & (sh[:, 1] | sh[:, 2])
    pref13 = d13 <= d24
    use13 = ok13 & (pref13 | ~ok24)
    longer = int((sel & (use13 != pref13) & (ok13 | ok24)).sum())
    whole = int((sel & ~ok13 & ~ok24).sum())
    sel &= ok13 | ok24

    e = m.elem
    t1 = np.where(use13[:, None], e[:, [0, 1, 2, 0]], e[:, [0, 1, 3, 0]])
    t2 = np.where(use13[:, None], e[:, [0, 2, 3, 0]], e[:, [1, 2, 3, 1]])

    # element k -> 1 or 2 output rows, in place
    nrow = np.where(sel, 2, 1)
    start = np.concatenate([[0], np.cumsum(nrow)[:-1]])
    out = np.empty((nrow.sum(), 4), dtype=np.int64)
    out[start[~sel]] = e[~sel]
    out[start[sel]] = t1[sel]
    out[start[sel] + 1] = t2[sel]
    return mesh_io.Mesh(m.ctok, m.flag, out, 4, m.depth), int(sel.sum()), longer, whole


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--band", nargs=2, type=float, metavar=("LATMIN", "LATMAX"),
                    help="split every quad in this latitude band (default: none)")
    ap.add_argument("--every", type=int, default=0, metavar="K",
                    help="split every K-th quad outside the band (default 0: none)")
    a = ap.parse_args()
    m = mesh_io.read_mesh(a.mesh_dir)
    q0 = int((m.nnodes() == 4).sum())
    s, nsplit, longer, whole = split(m, a.band, a.every)
    mesh_io.write_mesh(s, a.out_dir)
    print("split_quads: %d of %d quads split: %d elements (%d quads, %d triangles) -> %s"
          % (nsplit, q0, s.e2d, int((s.nnodes() == 4).sum()), int((s.nnodes() == 3).sum()),
             a.out_dir))
    print("  %d cut along the longer diagonal, %d left whole (avoiding dead-end triangles)"
          % (longer, whole))
    print("SPLIT quads_in=%d elems_in=%d split=%d elems_out=%d" % (q0, m.e2d, nsplit, s.e2d))


if __name__ == "__main__":
    main()
