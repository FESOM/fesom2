#!/usr/bin/env python3
"""Remove dead-end elements (fewer than two neighbours) from a FESOM mesh.

FESOM2's setup (find_elem_neighbors_ini) stops on elements with fewer than two
neighbours. This tool removes them, iterating until none are left (a removal can make
a neighbour a dead end), drops the nodes no longer used, renumbers nodes and elements
(original order kept) and subsets the depth file (depth.out / aux3d.out) to the kept
nodes. Nodes of a removed element that stay in the mesh get boundary flag 1 (the
removed cell is land now). nemo_meshq.py stays a faithful port; this is a separate step.
depth.out layer counts are kept as they are; a node whose count is now above the max over
its remaining elements (the model's bottom invariant) is reported (check_mesh.py fails).

Usage: remove_deadends.py <mesh_dir> <out_dir>
"""
import sys

import numpy as np

import mesh_io


def main(src, dst):
    m = mesh_io.read_mesh(src)
    n0, e0 = m.n2d, m.e2d
    m, removed = mesh_io.remove_deadends(m)
    mesh_io.write_mesh(m, dst)
    print("remove_deadends: %d element(s) and %d node(s) removed: %d nodes, %d elements -> %s"
          % (removed, n0 - m.n2d, m.n2d, m.e2d, dst))
    lc = mesh_io.layer_counts(m)
    if lc is not None and np.any(lc[0] > lc[1]):
        bad = np.nonzero(lc[0] > lc[1])[0]
        print("  WARNING: %d node(s) with a layer count above the max over their elements "
              "(bottom invariant; counts left unchanged): %s%s"
              % (len(bad), ' '.join(str(i + 1) for i in bad[:10]), ' ...' if len(bad) > 10 else ''))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
