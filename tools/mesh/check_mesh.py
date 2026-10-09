#!/usr/bin/env python3
"""Validity checks of a FESOM mesh with 3- or 4-column elem2d.out (triangles, quads, mixed).

Errors (exit status 1):
  - element not clockwise (seen from outside the sphere) or not convex;
  - degenerate element: zero side, collinear corner, or a repeated vertex other than a
    triangle's 4th == 1st;
  - node not used by any element;
  - edge shared by more than two elements;
  - dead end: element with fewer than two neighbours (FESOM's setup stops on them);
  - depth.out only: node layer count above the max over its elements of the element
    count (min over its vertices) -- the model's bottom invariant (node levels come from
    depth.out unchanged, so fesom_ini would pass the count on and the model stop).
Information: element types, edges (boundary edges), connected components, nodes whose
boundary flag differs from "lies on a boundary edge", total spherical area.

The second-to-last line is a machine-readable summary
  SUMMARY nodes=.. elems=.. quads=.. tris=.. edges=.. bedges=.. ccw=.. nonconvex=..
          degenerate=.. unused=.. multiedge=.. deadends=.. layers=.. area=<m^2, %.15e>
and the last line PASS or FAIL.

Usage: check_mesh.py <mesh_dir>
"""
import sys

import numpy as np

import mesh_io

TOL = 1e-9                          # |sin(turn angle)| below this: collinear corner


def corner_turns(m):
    """sin of the turn angle at every vertex, positive = counter-clockwise turn, (e2d, 4)
    (slot 4 of a triangle = +inf, so it never counts), and a flag for zero sides."""
    nv = m.nnodes()
    p = mesh_io.xyz(m.lon, m.lat)[m.elem]                 # (e2d, 4, 3)
    s = np.full(m.elem.shape, np.inf)
    zero = np.zeros(m.e2d, dtype=bool)
    for n in (3, 4):
        sel = nv == n
        q = p[sel, :n]
        e_in = q - np.roll(q, 1, axis=1)                  # p_j - p_{j-1}
        e_out = np.roll(q, -1, axis=1) - q                # p_{j+1} - p_j
        t = np.einsum('ijk,ijk->ij', np.cross(e_in, e_out), q)
        l_in = np.linalg.norm(e_in, axis=2)
        l_out = np.linalg.norm(e_out, axis=2)
        zero[sel] = np.any((l_in == 0) | (l_out == 0), axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            s[sel, :n] = np.nan_to_num(t / (l_in * l_out))
    return s, zero


def spherical_area(m):
    """Area of every element on the unit sphere, positive for clockwise elements:
    triangles (1,2,3) + (1,3,4) (the second vanishes for a triangle)."""
    p = mesh_io.xyz(m.lon, m.lat)[m.elem]

    def tri(a, b, c):
        num = np.einsum('ij,ij->i', a, np.cross(b, c))
        den = 1 + np.einsum('ij,ij->i', a, b) + np.einsum('ij,ij->i', b, c) \
            + np.einsum('ij,ij->i', c, a)
        return -2 * np.arctan2(num, den)

    return tri(p[:, 0], p[:, 1], p[:, 2]) + tri(p[:, 0], p[:, 2], p[:, 3])


def check(m):
    """Dictionary of counts (see module doc) for the mesh m."""
    r = {}
    nv = m.nnodes()
    r['nodes'], r['elems'] = m.n2d, m.e2d
    r['quads'], r['tris'] = int((nv == 4).sum()), int((nv == 3).sum())

    # a repeated vertex other than the triangle's 4th == 1st
    e = m.elem
    rep = (e[:, 0] == e[:, 1]) | (e[:, 1] == e[:, 2]) | (e[:, 2] == e[:, 3]) \
        | (e[:, 0] == e[:, 2]) | (e[:, 1] == e[:, 3])

    s, zero = corner_turns(m)
    bad = rep | zero | np.any(np.abs(s) < TOL, axis=1)
    cw = np.all(np.where(np.isinf(s), -TOL, s) <= -TOL, axis=1)
    ccw = np.all(np.where(np.isinf(s), TOL, s) >= TOL, axis=1)
    r['degenerate'] = int(bad.sum())
    r['ccw'] = int((~bad & ccw).sum())
    r['nonconvex'] = int((~bad & ~cw & ~ccw).sum())

    used = np.zeros(m.n2d, dtype=bool)
    used[e.ravel()] = True
    r['unused'] = int((~used).sum())

    ed, side_edge, count = mesh_io.edges(m)
    r['edges'], r['bedges'] = len(ed), int((count == 1).sum())
    r['multiedge'] = int((count > 2).sum())
    nb = mesh_io.neighbour_count(m, side_edge, count)
    r['deadends'] = int((nb < 2).sum())
    r['components'] = mesh_io.components(m, side_edge, count)[0]
    onb = mesh_io.boundary_nodes(m, ed, count)
    r['flag_mismatch'] = int(((m.flag != 0) != onb).sum())
    r['area'] = float(spherical_area(m).sum()) * mesh_io.R_EARTH ** 2
    lc = mesh_io.layer_counts(m)
    above = np.nonzero(lc[0] > lc[1])[0] if lc is not None else np.zeros(0, dtype=np.int64)
    r['layers'] = len(above)

    r['_lists'] = {'degenerate': np.nonzero(bad)[0],
                   'ccw': np.nonzero(~bad & ccw)[0],
                   'nonconvex': np.nonzero(~bad & ~cw & ~ccw)[0],
                   'unused': np.nonzero(~used)[0],
                   'deadends': np.nonzero(nb < 2)[0],
                   'layers': above}
    return r


ERRORS = ('ccw', 'nonconvex', 'degenerate', 'unused', 'multiedge', 'deadends', 'layers')


def main(d):
    m = mesh_io.read_mesh(d)
    r = check(m)
    print('check_mesh %s' % d)
    print('  nodes %d, elements %d (quads %d, triangles %d), elem2d.out %d columns'
          % (r['nodes'], r['elems'], r['quads'], r['tris'], m.ncol))
    print('  edges %d, boundary edges %d, connected components %d'
          % (r['edges'], r['bedges'], r['components']))
    print('  boundary flag differs from "on a boundary edge" at %d node(s)' % r['flag_mismatch'])
    print('  total area %.6e m^2 (%.6f of the sphere)'
          % (r['area'], r['area'] / (4 * np.pi * mesh_io.R_EARTH ** 2)))
    fail = False
    for k, what in (('ccw', 'counter-clockwise elements'),
                    ('nonconvex', 'non-convex elements'),
                    ('degenerate', 'degenerate elements'),
                    ('unused', 'unused nodes'),
                    ('multiedge', 'edges in more than two elements'),
                    ('deadends', 'elements with fewer than two neighbours'),
                    ('layers', 'nodes above the bottom invariant')):
        n = r[k]
        ids = r['_lists'].get(k)
        more = ''
        if n and ids is not None:
            more = ': ' + ' '.join(str(i + 1) for i in ids[:10]) + (' ...' if n > 10 else '')
        print('  %-40s %8d %s%s' % (what, n, 'ERROR' if n else 'ok', more))
        fail = fail or n > 0
    print('SUMMARY ' + ' '.join('%s=%d' % (k, r[k]) for k in
                                ('nodes', 'elems', 'quads', 'tris', 'edges', 'bedges')
                                + ERRORS) + ' area=%.15e' % r['area'])
    print('FAIL' if fail else 'PASS')
    return 1 if fail else 0


if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0 if len(sys.argv) == 2 else 2)
    sys.exit(main(sys.argv[1]))
