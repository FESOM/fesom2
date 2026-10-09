#!/usr/bin/env python3
"""Structural checks of fesom_ini output (edges, levels, dist_N) for triangle, quad and
mixed meshes. Everything is recomputed from nod2d.out / elem2d.out and the partition
vector (dist_N/rpart.out) and compared with the files:

Mesh files (once per mesh directory):
  - edges.out: the edges are exactly the element sides (1-2, 2-3, 3-4, 4-1; a triangle's
    degenerate side 4-1 dropped), each once; edgenum.out: edge2D and edge2D_in;
  - edge order: edges 1..edge2D_in lie in two elements, edge2D_in+1..edge2D in one;
  - edge_tri.out: the element(s) holding the edge (2nd <= 0 on the boundary), the 1st
    left of edges(1) -> edges(2) on the sphere and the 2nd right of it (the model's
    convention: the SSH stiffness matrix is symmetric only with it);
  - levels: with depth.out (NEMO) nlvls.out == layer counts + 1 and elvls.out == min over
    the element's vertices; otherwise (FESOM2 level finder) nlvls.out == max of elvls.out
    over the elements around the node.
dist_N (per number of parts):
  - every node owned by exactly one rank (owned list == nodes of part r, ascending);
  - every element owned once: in the owned list of the rank of its 1st vertex (its main
    owner); owned element list == elements with an owned vertex (ascending);
  - halo lists complete: halo nodes == vertices of owned elements not owned; eDim
    elements == side neighbours of owned elements without an owned vertex; eXDim ==
    remaining vertex neighbours; owned edges == edges with an owned node (ascending),
    eDim edges == sides of owned elements without an owned node;
  - com_info: node / element / full-element halos symmetric between ranks (what r receives
    from p == what p sends to r, same order), received from their owner (node: part,
    element: part of the 1st vertex), receive lists in local numbering.
  - appended sections (hybrid meshes; only if dist_N/format.out, version 1, is present):
    my_list eXDim_nod2D == vertices of the local elements 1..myDim+eDim+eXDim not among
    the nodes 1..myDim+eDim, grouped by owner, ascending within; com_nod2D_full receives
    every halo node (eDim + eXDim) exactly once, from its owner; com_edge2D receives every
    held edge whose unique owner (owner of edges(1,edge)) is another rank exactly once, from
    that owner -- incl. the edges that also lie in the receiving rank's 1..myDim_edge2D
    (owned edge ranges overlap: counted); both symmetric between ranks, ordered by global
    id. Without format.out the files must end after the FESOM2 part.

Usage: check_dist.py <mesh_dir> <npart> [<npart> ...]
  exit status 0 if every check passes, 1 otherwise; last line "check_dist: PASS|FAIL".
"""
import os
import sys

import numpy as np

import mesh_io


class Checker:
    def __init__(self):
        self.nfail = 0

    def check(self, cond, what):
        if cond:
            print('  ok: ' + what)
        else:
            print('  FAILED: ' + what)
            self.nfail += 1
        return cond


def ints(path):
    with open(path) as f:
        return np.fromstring(f.read(), dtype=np.int64, sep=' ')


class Tokens:
    """Sequential reader of the list-directed integer files of dist_N."""

    def __init__(self, path):
        self.v = ints(path)
        self.i = 0

    def take(self, n=None):
        if n is None:
            self.i += 1
            return int(self.v[self.i - 1])
        self.i += n
        return self.v[self.i - n:self.i]

    def done(self):
        return self.i == len(self.v)


def read_com(t):
    """One communication structure: {peer: list} for receive and send (local numbering)."""
    out = []
    for _ in range(2):                                  # receive, send
        npe = t.take()
        pe = t.take(npe)
        ptr = t.take(npe + 1)
        lst = t.take(int(ptr[-1]) - 1)
        out.append({int(p): lst[ptr[k] - 1:ptr[k + 1] - 1] for k, p in enumerate(pe)})
    return out


def check_mesh_files(c, d, m, nv, side_edge, count):
    print('== mesh files of %s (%d nodes, %d elements: %d quads, %d triangles)'
          % (d, m.n2d, m.e2d, int((nv == 4).sum()), int((nv == 3).sum())))
    e2d_n, e2d_in = ints(os.path.join(d, 'edgenum.out'))[:2]
    ed = ints(os.path.join(d, 'edges.out')).reshape(-1, 2) - 1
    et = ints(os.path.join(d, 'edge_tri.out')).reshape(-1, 2)
    nsides = len(count)
    c.check(len(ed) == e2d_n == nsides and len(et) == e2d_n,
            'edge count: edges.out %d, edgenum.out %d, element sides %d'
            % (len(ed), e2d_n, nsides))
    key = mesh_io.edge_key(ed[:, 0], ed[:, 1], m.n2d)
    skey = mesh_io.edge_key(m.elem, np.roll(m.elem, -1, axis=1), m.n2d)
    uk = np.unique(skey[side_edge >= 0])
    c.check(len(np.unique(key)) == len(key) and np.array_equal(np.sort(key), uk),
            'edges.out == element sides, each once')
    # global edge number (edges.out order) of every element side, -1 for a triangle's 4-1
    order = np.argsort(key)
    elem_edges = np.full(m.elem.shape, -1, dtype=np.int64)
    real = side_edge >= 0
    elem_edges[real] = order[np.searchsorted(key[order], skey[real])]
    nel = np.bincount(elem_edges[real], minlength=len(key))
    c.check(int((count == 2).sum()) == e2d_in and np.all(nel[:e2d_in] == 2)
            and np.all(nel[e2d_in:] == 1),
            'edge order: %d interior edges first (edge2D_in = %d), then %d boundary edges'
            % (int((nel == 2).sum()), e2d_in, int((nel == 1).sum())))
    # edge_tri: the elements holding the edge
    el = np.broadcast_to(np.arange(m.e2d)[:, None], m.elem.shape)[real]
    eg = elem_edges[real]
    o = np.lexsort((el, eg))
    first = np.full(len(key), -1, dtype=np.int64)
    second = np.full(len(key), -1, dtype=np.int64)
    eg_o, el_o = eg[o], el[o]
    start = np.r_[True, eg_o[1:] != eg_o[:-1]]
    first[eg_o[start]] = el_o[start]
    second[eg_o[~start]] = el_o[~start]
    t1, t2 = et[:, 0] - 1, et[:, 1] - 1
    ok_in = np.all(np.sort(np.stack([t1, t2], 1)[:e2d_in], 1)
                   == np.stack([first, second], 1)[:e2d_in])
    ok_bd = np.all(t1[e2d_in:] == first[e2d_in:]) and np.all(et[e2d_in:, 1] <= 0)
    c.check(ok_in and ok_bd, 'edge_tri.out == elements holding each edge (2nd <= 0 on the boundary)')
    # orientation: the element centre (vertex mean) against the great circle of the edge
    p = mesh_io.xyz(m.lon, m.lat)
    nv = m.nnodes()
    cen = sum(np.where((k < nv)[:, None], p[m.elem[:, k]], 0.0) for k in range(4)) / nv[:, None]
    nrm = np.cross(p[ed[:, 0]], p[ed[:, 1]])
    inner = et[:, 1] > 0
    left1 = np.einsum('ij,ij->i', nrm, cen[t1]) > 0
    right2 = np.einsum('ij,ij->i', nrm[inner], cen[t2[inner]]) < 0
    c.check(bool(left1.all() and right2.all()),
            'edge_tri.out: 1st element left of edges(1) -> edges(2), 2nd right (%d, %d wrong)'
            % (int((~left1).sum()), int((~right2).sum())))

    nlv = ints(os.path.join(d, 'nlvls.out'))
    elv = ints(os.path.join(d, 'elvls.out'))
    c.check(len(nlv) == m.n2d and len(elv) == m.e2d, 'nlvls.out / elvls.out lengths')
    if os.path.exists(os.path.join(d, 'depth.out')):
        head, blocks = m.depth['depth.out']
        cnt = np.array(blocks[1].tolist(), dtype=np.int64)
        c.check(np.array_equal(nlv, cnt + 1), 'nlvls.out == depth.out layer counts + 1')
        emin = nlv[m.elem].min(1)                       # (triangle: 4th == 1st)
        c.check(np.array_equal(elv, emin), 'elvls.out == min of nlvls.out over the vertices')
    else:
        nmax = np.zeros(m.n2d, dtype=np.int64)
        np.maximum.at(nmax, m.elem.ravel(), np.repeat(elv, 4))
        c.check(np.array_equal(nlv, nmax), 'nlvls.out == max of elvls.out over the elements around')
    return ed, elem_edges


def check_dist(c, d, npart, m, nv, side_edge, count, ed, elem_edges):
    dd = os.path.join(d, 'dist_%d' % npart)
    print('== %s' % dd)
    rp = ints(os.path.join(dd, 'rpart.out'))
    n2d, e2d = m.n2d, m.e2d
    c.check(rp[0] == npart and rp[1:1 + npart].sum() == n2d and len(rp) == 1 + npart + n2d,
            'rpart.out: %d parts, %d nodes' % (rp[0], rp[1:1 + npart].sum()))
    offs = np.r_[0, np.cumsum(rp[1:1 + npart])]
    part = np.searchsorted(offs, rp[1 + npart:], side='left') - 1     # contiguous id -> part
    en = m.elem
    owner = part[en[:, 0]]                                  # main owner of every element
    # side neighbours (-1 = none)
    real = side_edge >= 0
    el = np.broadcast_to(np.arange(e2d)[:, None], en.shape)[real]
    eg = side_edge[real]
    o = np.lexsort((el, eg))
    eg_o, el_o = eg[o], el[o]
    pair = np.nonzero(eg_o[1:] == eg_o[:-1])[0]
    nbr_a, nbr_b = el_o[pair], el_o[pair + 1]

    nod_owned = np.zeros(n2d, dtype=np.int64)
    el_owned_main = np.zeros(e2d, dtype=bool)
    ed_owned = np.zeros(len(ed), dtype=np.int64)
    recv = {}                                               # (kind, r, p) -> global ids
    send = {}
    bad = {'nodes': 0, 'halo nodes': 0, 'elements': 0, 'eDim elements': 0,
           'eXDim elements': 0, 'edges': 0, 'eDim edges': 0, 'files': 0,
           'receive lists': 0, 'owners': 0, 'eXDim nodes': 0, 'nod_full': 0, 'edge com': 0}
    fpath = os.path.join(dd, 'format.out')
    fmt = int(ints(fpath)[0]) if os.path.exists(fpath) else 0
    if fmt:
        c.check(fmt == 1, 'format.out: version %d' % fmt)
    else:
        print('  info: no format.out (FESOM2 layout only): appended sections not checked')
    uowner = part[ed[:, 0]]                                 # unique owner of every edge
    n_overlap = 0
    for r in range(npart):
        t = Tokens(os.path.join(dd, 'my_list%05d.out' % r))
        if t.take() != r:
            bad['files'] += 1
        mydn, edn = t.take(), t.take()
        nodes = t.take(mydn + edn) - 1
        myde, ede, exde = t.take(), t.take(), t.take()
        elems = t.take(myde + ede + exde) - 1
        mydg, edg = t.take(), t.take()
        edges_l = t.take(mydg + edg) - 1
        if fmt:
            exdn = t.take()
            xnodes = t.take(exdn) - 1
        if not t.done():
            bad['files'] += 1
        own = part == r
        if not np.array_equal(nodes[:mydn], np.nonzero(own)[0]):
            bad['nodes'] += 1
        nod_owned[nodes[:mydn]] += 1
        e_mine = own[en].any(1)
        if not np.array_equal(elems[:myde], np.nonzero(e_mine)[0]):
            bad['elements'] += 1
        el_owned_main[elems[:myde][owner[elems[:myde]] == r]] = True
        touched = np.zeros(n2d, dtype=bool)
        touched[en[e_mine].ravel()] = True
        if not np.array_equal(np.sort(nodes[mydn:]), np.nonzero(touched & ~own)[0]):
            bad['halo nodes'] += 1
        nb = np.zeros(e2d, dtype=bool)
        nb[nbr_b[e_mine[nbr_a]]] = True
        nb[nbr_a[e_mine[nbr_b]]] = True
        nb &= ~e_mine
        if not np.array_equal(np.sort(elems[myde:myde + ede]), np.nonzero(nb)[0]):
            bad['eDim elements'] += 1
        nx = touched[en].any(1) & ~e_mine & ~nb
        if not np.array_equal(np.sort(elems[myde + ede:]), np.nonzero(nx)[0]):
            bad['eXDim elements'] += 1
        g_own = own[ed].any(1)
        if not np.array_equal(edges_l[:mydg], np.nonzero(g_own)[0]):
            bad['edges'] += 1
        ed_owned[edges_l[:mydg]] += 1
        se = np.unique(elem_edges[e_mine])
        se = se[(se >= 0) & ~g_own[se.clip(0)]]
        if not np.array_equal(np.sort(edges_l[mydg:]), se):
            bad['eDim edges'] += 1

        t = Tokens(os.path.join(dd, 'com_info%05d.out' % r))
        if t.take() != r:
            bad['files'] += 1
        cn = read_com(t)
        ce = read_com(t)
        cf = read_com(t)
        if fmt:
            cnf = read_com(t)
            cg = read_com(t)
        if not t.done():
            bad['files'] += 1
        # receive lists: local numbering, owner = peer
        rl = np.concatenate([v for v in cn[0].values()]) if cn[0] else np.zeros(0, np.int64)
        if not np.array_equal(rl, np.arange(mydn + 1, mydn + edn + 1)):
            bad['receive lists'] += 1
        rl = np.concatenate([v for v in ce[0].values()]) if ce[0] else np.zeros(0, np.int64)
        if not np.array_equal(rl, np.arange(myde + 1, myde + ede + 1)):
            bad['receive lists'] += 1
        for p, v in cn[0].items():
            g = nodes[v - 1]
            bad['owners'] += int(np.any(part[g] != p))
            recv[('node', r, p)] = g
        for p, v in cn[1].items():
            send[('node', r, p)] = nodes[v - 1]
        for kind, cs in (('elem', ce), ('elem_full', cf)):
            for p, v in cs[0].items():
                g = elems[v - 1]
                bad['owners'] += int(np.any(owner[g] != p))
                recv[(kind, r, p)] = g
            for p, v in cs[1].items():
                send[(kind, r, p)] = elems[v - 1]
        full = np.sort(np.concatenate([elems[v - 1] for v in cf[0].values()])) \
            if cf[0] else np.zeros(0, np.int64)
        if not np.array_equal(full, np.sort(elems[myde:])):
            bad['receive lists'] += 1
        if fmt:
            # double-halo nodes: vertices of all local elements not among the local nodes
            vx = np.setdiff1d(np.unique(en[elems]), nodes)
            vx = vx[np.lexsort((vx, part[vx]))]                 # by owner, ascending within
            if not np.array_equal(xnodes, vx):
                bad['eXDim nodes'] += 1
            alln = np.concatenate([nodes, xnodes])
            nloc = mydn + edn + exdn
            # com_nod2D_full: every halo node once, from its owner, ordered by global id
            rl = np.sort(np.concatenate([v for v in cnf[0].values()])) if cnf[0] \
                else np.zeros(0, np.int64)
            if not np.array_equal(rl, np.arange(mydn + 1, nloc + 1)):
                bad['nod_full'] += 1
            for p, v in cnf[0].items():
                g = alln[v - 1]
                if np.any(part[g] != p) or np.any(np.diff(g) <= 0):
                    bad['nod_full'] += 1
                recv[('nod_full', r, p)] = g
            for p, v in cnf[1].items():
                if np.any(v < 1) or np.any(v > mydn):
                    bad['nod_full'] += 1
                send[('nod_full', r, p)] = alln[np.clip(v, 1, nloc) - 1]
            # com_edge2D: every held edge of another unique owner once, from that owner
            nedl = mydg + edg
            want = np.nonzero(uowner[edges_l] != r)[0] + 1
            rl = np.sort(np.concatenate([v for v in cg[0].values()])) if cg[0] \
                else np.zeros(0, np.int64)
            if not np.array_equal(rl, want):
                bad['edge com'] += 1
            for p, v in cg[0].items():
                g = edges_l[v - 1]
                if np.any(uowner[g] != p) or np.any(np.diff(g) <= 0):
                    bad['edge com'] += 1
                recv[('edge', r, p)] = g
            for p, v in cg[1].items():
                g = edges_l[np.clip(v, 1, nedl) - 1]
                if np.any(v < 1) or np.any(v > nedl) or np.any(uowner[g] != r):
                    bad['edge com'] += 1
                send[('edge', r, p)] = g
            # owned edge ranges overlap: edges in 1..myDim_edge2D of another unique owner
            ovl = np.nonzero(uowner[edges_l[:mydg]] != r)[0] + 1
            n_overlap += len(ovl)
            if not np.all(np.isin(ovl, rl)):
                bad['edge com'] += 1

    c.check(bad['files'] == 0, 'my_list/com_info files: rank ids and lengths')
    c.check(np.all(nod_owned == 1) and bad['nodes'] == 0,
            'every node owned by exactly one rank (owned lists == rpart.out)')
    c.check(bool(el_owned_main.all()) and bad['elements'] == 0,
            'every element owned once (in the list of the rank of its 1st vertex); '
            'owned elements == elements with an owned vertex')
    c.check(np.all(ed_owned >= 1) and bad['edges'] == 0,
            'every edge owned (owned edges == edges with an owned node)')
    c.check(bad['halo nodes'] == 0, 'halo nodes == vertices of owned elements, not owned')
    c.check(bad['eDim elements'] == 0 and bad['eXDim elements'] == 0,
            'halo elements: eDim == side neighbours, eXDim == other vertex neighbours')
    c.check(bad['eDim edges'] == 0, 'halo edges == sides of owned elements without an owned node')
    c.check(bad['receive lists'] == 0 and bad['owners'] == 0,
            'receive lists in local numbering, received from the owner')
    if fmt:
        c.check(bad['eXDim nodes'] == 0,
                'eXDim_nod2D == vertices of the local elements not among the local nodes, '
                'grouped by owner')
        c.check(bad['nod_full'] == 0,
                'com_nod2D_full: every halo node (eDim + eXDim) received once, from its owner')
        c.check(bad['edge com'] == 0,
                'com_edge2D: every held edge of another unique owner received once, from it '
                '(%d edges in overlapping owned ranges)' % n_overlap)
    asym = [k for k in set(recv) | set(send)
            if not np.array_equal(recv.get(k, []), send.get((k[0], k[2], k[1]), []))]
    c.check(not asym, 'node / element / full-element%s halos symmetric between ranks '
                      '(%d receive lists%s)' % (' / full-node / edge' if fmt else '', len(recv),
                                                '' if not asym else ', asymmetric: %s' % asym[:5]))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    d = sys.argv[1].rstrip('/')
    c = Checker()
    m = mesh_io.read_mesh(d)
    nv = m.nnodes()
    ed_u, side_edge, count = mesh_io.edges(m)
    ed, elem_edges = check_mesh_files(c, d, m, nv, side_edge, count)
    for np_ in sys.argv[2:]:
        check_dist(c, d, int(np_), m, nv, side_edge, count, ed, elem_edges)
    print('check_dist: %s' % ('PASS' if c.nfail == 0 else 'FAIL (%d)' % c.nfail))
    sys.exit(0 if c.nfail == 0 else 1)


if __name__ == '__main__':
    main()
