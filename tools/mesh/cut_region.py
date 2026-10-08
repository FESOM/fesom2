#!/usr/bin/env python3
"""Cut a lon/lat box out of a FESOM mesh (small test meshes).

Kept are the elements with all vertices inside the box [LONMIN, LONMAX] x [LATMIN, LATMAX]
(longitudes compared modulo 360, so a box may cross the dateline). Of these only the
largest connected piece (elements connected over sides) is kept, and dead ends are
removed as in remove_deadends.py. Nodes are renumbered (original order kept), the
boundary flag is recomputed as "lies on a boundary edge" -- the open-ocean sides of the
box become coast like the land boundary -- and the depth file is subset to the kept nodes.

Usage: cut_region.py <mesh_dir> <out_dir> --lon LONMIN LONMAX --lat LATMIN LATMAX
"""
import argparse

import numpy as np

import mesh_io


def cut(m, lon, lat):
    """Cut-out of m (see module doc); returns (mesh, number of pieces dropped)."""
    inside = (np.mod(m.lon - lon[0], 360.0) <= lon[1] - lon[0]) \
        & (m.lat >= lat[0]) & (m.lat <= lat[1])
    c = mesh_io.subset(m, np.all(inside[m.elem], axis=1))
    if c.e2d == 0:
        raise SystemExit("cut_region: no element inside the box")
    ncomp, label = mesh_io.components(c)
    big = np.argmax(np.bincount(label))
    c = mesh_io.subset(c, label == big)
    c, _ = mesh_io.remove_deadends(c)
    c.flag = mesh_io.boundary_nodes(c).astype(np.int64)
    return c, ncomp - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mesh_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--lon", nargs=2, type=float, required=True, metavar=("LONMIN", "LONMAX"))
    ap.add_argument("--lat", nargs=2, type=float, required=True, metavar=("LATMIN", "LATMAX"))
    a = ap.parse_args()
    if not 0 < a.lon[1] - a.lon[0] <= 360 or a.lat[1] <= a.lat[0]:
        ap.error("need LONMIN < LONMAX <= LONMIN+360 and LATMIN < LATMAX")
    m = mesh_io.read_mesh(a.mesh_dir)
    c, dropped = cut(m, a.lon, a.lat)
    mesh_io.write_mesh(c, a.out_dir)
    nv = c.nnodes()
    print("cut_region: lon [%g, %g] lat [%g, %g]: %d nodes, %d elements (%d quads, "
          "%d triangles), %d smaller piece(s) dropped -> %s"
          % (a.lon[0], a.lon[1], a.lat[0], a.lat[1], c.n2d, c.e2d, int((nv == 4).sum()),
             int((nv == 3).sum()), dropped, a.out_dir))


if __name__ == "__main__":
    main()
