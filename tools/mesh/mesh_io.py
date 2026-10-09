"""Shared helpers of the mesh tools: read/write a FESOM mesh directory, edges, neighbours.

A mesh directory holds
  nod2d.out  : n2d, then  id lon lat flag  (flag: 1 = land-boundary node)
  elem2d.out : e2d, then 3 or 4 node ids per element (3 columns -> 4th := 1st;
               a triangle in 4 columns repeats its 1st vertex as the 4th)
and optionally a depth file, which is subset line by line (kept lines byte-identical):
  depth.out  : nl, nl levels, n2d node depths, n2d layer counts  (nemo_meshq)
  aux3d.out  : nl, nl zbar, n2d node depths                      (FESOM2)
Coordinates are kept as the original text tokens, so kept nodes are written unchanged.
"""
import glob
import os
import shutil

import numpy as np
import scipy.sparse
import scipy.sparse.csgraph

R_EARTH = 6367500.0                 # FESOM r_earth (m)
DEPTH_FILES = ('depth.out', 'aux3d.out')
# fesom_ini outputs, derived from elem2d.out: stale once the elements change
INI_OUTPUTS = ('edges.out', 'edge_tri.out', 'edgenum.out', 'nlvls.out', 'elvls.out',
               'elvls_raw.out')


class Mesh:
    """lon/lat in degrees (float), ctok: (n2d, 2) coordinate text tokens, flag (n2d,) int,
    elem: (e2d, 4) 0-based node indices, ncol: 3 or 4 (columns of elem2d.out),
    depth: {file name: (header lines, node line blocks)} for the depth files present."""

    def __init__(self, ctok, flag, elem, ncol=4, depth=None):
        self.ctok = ctok
        self.lon = ctok[:, 0].astype(float)
        self.lat = ctok[:, 1].astype(float)
        self.flag = flag
        self.elem = elem
        self.ncol = ncol
        self.depth = depth if depth is not None else {}

    @property
    def n2d(self):
        return len(self.flag)

    @property
    def e2d(self):
        return len(self.elem)

    def nnodes(self):
        """Vertices per element: 3 (4th == 1st) or 4."""
        return np.where(self.elem[:, 3] == self.elem[:, 0], 3, 4)


def read_mesh(d):
    tok = open(os.path.join(d, 'nod2d.out')).read().split()
    n2d = int(tok[0])
    t = np.array(tok[1:1 + 4 * n2d]).reshape(n2d, 4)
    if not np.array_equal(t[:, 0].astype(np.int64), np.arange(1, n2d + 1)):
        raise SystemExit('%s/nod2d.out: node ids are not 1..n2d' % d)
    ctok = t[:, 1:3].copy()
    flag = t[:, 3].astype(np.int64)

    with open(os.path.join(d, 'elem2d.out')) as f:
        e2d = int(f.readline())
        first = f.readline()
        ncol = len(first.split())
        rest = f.read()
    if ncol not in (3, 4):
        raise SystemExit('%s/elem2d.out: %d columns (3 or 4 expected)' % (d, ncol))
    elem = np.array((first + rest).split(), dtype=np.int64)
    if elem.size != ncol * e2d:
        raise SystemExit('%s/elem2d.out: %d values for %d elements of %d columns'
                         % (d, elem.size, e2d, ncol))
    elem = elem.reshape(e2d, ncol) - 1
    if ncol == 3:
        elem = np.hstack([elem, elem[:, :1]])

    depth = {}
    for name in DEPTH_FILES:
        p = os.path.join(d, name)
        if not os.path.exists(p):
            continue
        lines = open(p).read().split('\n')
        if lines and lines[-1] == '':
            lines = lines[:-1]
        nl = int(lines[0])
        head, body = lines[:1 + nl], lines[1 + nl:]
        nblk = 2 if name == 'depth.out' else 1
        if len(body) != nblk * n2d:
            raise SystemExit('%s: %d node lines, %d expected (one value per line)'
                             % (p, len(body), nblk * n2d))
        depth[name] = (head, [np.array(body[k * n2d:(k + 1) * n2d], dtype=object)
                              for k in range(nblk)])
    return Mesh(ctok, flag, elem, ncol, depth)


def write_mesh(m, d):
    """Writes nod2d.out, elem2d.out (4 columns) and the depth files of m into d. fesom_ini
    outputs already in d (edge and level files, dist_*/) belong to the old elements and are
    removed: fesom_ini keeps existing edge files, so it would pair them with new dist_N."""
    os.makedirs(d, exist_ok=True)
    stale = [p for p in [os.path.join(d, f) for f in INI_OUTPUTS] if os.path.exists(p)] \
        + sorted(glob.glob(os.path.join(d, 'dist_*')))
    for p in stale:
        if os.path.isdir(p):
            shutil.rmtree(p)
        else:
            os.remove(p)
    if stale:
        print('  write_mesh: removed the stale fesom_ini output in %s: %s (rerun fesom_ini)'
              % (d, ' '.join(os.path.basename(p) for p in stale)))
    with open(os.path.join(d, 'nod2d.out'), 'w') as f:
        f.write('%9i \n' % m.n2d)
        f.write(''.join('%9i %8s %8s %8i\n' % (i + 1, x, y, g) for i, (x, y, g) in
                        enumerate(zip(m.ctok[:, 0].tolist(), m.ctok[:, 1].tolist(),
                                      m.flag.tolist()))))
    with open(os.path.join(d, 'elem2d.out'), 'w') as f:
        f.write('%9i \n' % m.e2d)
        f.write(''.join('%9i %9i %9i %9i\n' % tuple(q) for q in (m.elem + 1).tolist()))
    for name, (head, blocks) in m.depth.items():
        with open(os.path.join(d, name), 'w') as f:
            f.write('\n'.join(head) + '\n')
            for b in blocks:
                f.write('\n'.join(b.tolist()) + '\n')


def subset(m, keep_elem):
    """Mesh with the elements keep_elem (bool or index array, order kept) and only the
    nodes they use (original order kept, renumbered); flags and depth lines follow."""
    elem = m.elem[keep_elem]
    used = np.zeros(m.n2d, dtype=bool)
    used[elem.ravel()] = True
    new = np.full(m.n2d, -1, dtype=np.int64)
    new[used] = np.arange(used.sum())
    depth = {name: (head, [b[used] for b in blocks])
             for name, (head, blocks) in m.depth.items()}
    return Mesh(m.ctok[used], m.flag[used], new[elem], m.ncol, depth)


def edge_key(a, b, n2d):
    """Orientation-free key of the node pairs (a, b): min * n2d + max."""
    return np.minimum(a, b) * np.int64(n2d) + np.maximum(a, b)


def side_edges(elem, n2d):
    """edges() of the elements elem (e2d, 4) over n2d nodes."""
    a, b = elem, np.roll(elem, -1, axis=1)
    key = edge_key(a, b, n2d)
    real = a != b
    uk, inv, count = np.unique(key[real], return_inverse=True, return_counts=True)
    side_edge = np.full(a.shape, -1, dtype=np.int64)
    side_edge[real] = inv
    ed = np.stack([uk // n2d, uk % n2d], axis=1)
    return ed, side_edge, count


def edges(m):
    """Unique edges over the element sides 1-2, 2-3, 3-4, 4-1 (the degenerate side 4-1 of a
    triangle is dropped). Returns (ed (ned, 2) sorted node pairs, side_edge (e2d, 4) edge
    index of each side or -1, count (ned,) number of elements per edge)."""
    return side_edges(m.elem, m.n2d)


def shared_sides(m, side_edge=None, count=None):
    """(e2d, 4) bool: the side is shared with another element (its edge is in 2)."""
    if side_edge is None:
        _, side_edge, count = edges(m)
    shared = np.zeros(side_edge.shape, dtype=bool)
    real = side_edge >= 0
    shared[real] = count[side_edge[real]] >= 2
    return shared


def neighbour_count(m, side_edge=None, count=None):
    """Number of elements sharing a side with each element (sides of edges in 2 elements)."""
    return shared_sides(m, side_edge, count).sum(axis=1)


def boundary_mask(n2d, ed, count):
    """Nodes (of n2d) on an edge of ed that belongs to one element only (count == 1)."""
    on = np.zeros(n2d, dtype=bool)
    on[ed[count == 1].ravel()] = True
    return on


def boundary_nodes(m, ed=None, count=None):
    """Nodes on an edge that belongs to one element only."""
    if ed is None:
        ed, _, count = edges(m)
    return boundary_mask(m.n2d, ed, count)


def components(m, side_edge=None, count=None):
    """Connected components of the elements over shared sides: (ncomp, label (e2d,))."""
    if side_edge is None:
        _, side_edge, count = edges(m)
    real = side_edge >= 0
    el = np.broadcast_to(np.arange(m.e2d)[:, None], side_edge.shape)[real]
    eg = side_edge[real]
    g = scipy.sparse.coo_matrix((np.ones(len(el)), (el, eg)),
                                shape=(m.e2d, len(count))).tocsr()
    adj = g @ g.T
    return scipy.sparse.csgraph.connected_components(adj, directed=False)


def layer_counts(m):
    """depth.out only (else None): (c, mx), c = the node layer counts, mx = max over a
    node's elements of the element count (min over its vertices). The model's bottom
    invariant is c == mx; c >= mx always, so a node with c > mx violates it."""
    if 'depth.out' not in m.depth:
        return None
    c = np.array([int(s) for s in m.depth['depth.out'][1][1]], dtype=np.int64)
    emin = c[m.elem].min(axis=1)
    mx = np.zeros(m.n2d, dtype=np.int64)
    for k in range(m.elem.shape[1]):
        np.maximum.at(mx, m.elem[:, k], emin)
    return c, mx


def xyz(lon, lat):
    """Unit vectors of lon/lat (degrees)."""
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def remove_deadends(m, verbose=True):
    """FESOM2's rule: drop elements with fewer than two neighbours until none are left.
    Nodes of removed elements that stay in the mesh become land-boundary nodes (flag 1),
    since the removed cell is land now. Returns (mesh, number removed)."""
    removed = 0
    it = 0
    while True:
        nb = neighbour_count(m)
        dead = nb < 2
        nd = int(dead.sum())
        if nd == 0:
            return m, removed
        it += 1
        if verbose:
            print('  remove_deadends: pass %d removes %d element(s)' % (it, nd))
        touched = np.zeros(m.n2d, dtype=bool)
        touched[m.elem[dead].ravel()] = True
        m = subset(Mesh(m.ctok, np.where(touched, 1, m.flag), m.elem, m.ncol, m.depth), ~dead)
        removed += nd
