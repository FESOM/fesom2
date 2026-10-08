#!/usr/bin/env python3
"""Python port of nemo_meshq.m: NEMO ORCA025 (1/4 deg) grid -> FESOM quad mesh.

Faithful, statement-by-statement port of the MATLAB script (Steps 1-6) that makes
nod2d.out, elem2d.out (4 columns, quads SW-NW-NE-SE) and depth.out from the NEMO
files mesh_hgr.nc, mesh_zgr.nc and mask.nc. The goal is byte-identical output, so
MATLAB semantics are emulated:
  - float32 inputs stay float32 and every single/double mix is computed in float32
    (MATLAB "single wins"): the f-point corrections, the southern extension row and
    the `levels` recursion;
  - mbathy stays integer (int16);
  - order-dependent loops stay sequential (single-cell land removal modifies bathy
    while scanning; quad creation and node compaction define the numbering);
  - MATLAB arrays are 1-based and column-major: indices are converted to 0-based,
    the MATLAB index is kept in a comment at each manual fix, node ids keep the
    MATLAB column-major numbering nny + Ny*(nnx-1).
Plotting statements of the .m file are omitted (no effect on the output). Dead-end
elements are NOT removed here (see remove_deadends.py) so that the port stays exact.

Usage: nemo_meshq.py <nemo_dir> <out_dir>
  nemo_dir: directory with mesh_hgr.nc, mesh_zgr.nc, mask.nc
  out_dir : where nod2d.out, elem2d.out, depth.out are written
"""
import os
import sys

import numpy as np
import netCDF4

f32 = np.float32
PI = f32(np.pi)


def ncread(path, name):
    """Raw variable (no auto mask/scale), leading singleton time dim dropped.

    netCDF4 returns (t, y, x) / (y, x), i.e. the transpose of MATLAB's ncread,
    which is what the .m file takes with its trailing ' anyway."""
    with netCDF4.Dataset(path) as d:
        v = d.variables[name]
        v.set_auto_maskandscale(False)
        a = np.array(v[:])
    if a.ndim == 3 and a.shape[0] == 1:
        a = a[0]
    return a


def main(nemo_dir, out_dir):
    # ==========
    # Step 1. Read and clean
    # ==========
    Nx = 1442
    Ny = 1021
    Nz = 46
    hgr = os.path.join(nemo_dir, 'mesh_hgr.nc')
    lon = ncread(hgr, 'nav_lon')        # lon, lat on T points, (Ny, Nx) float32
    lat = ncread(hgr, 'nav_lat')
    lev = ncread(hgr, 'nav_lev')        # standard midlevels, float32
    lon2 = ncread(hgr, 'glamf')         # f points
    lat2 = ncread(hgr, 'gphif')

    # Extend by one line in south, to close cells.
    lon2 = np.vstack([lon2[0:1, :], lon2])
    dy = lat[1, 0] - lat[0, 0]                           # float32
    lat2 = np.vstack([lat2[0:1, :] - dy, lat2])
    Ny = Ny + 1

    bathy = ncread(os.path.join(nemo_dir, 'mesh_zgr.nc'), 'mbathy')    # int16
    mask = ncread(os.path.join(nemo_dir, 'mask.nc'), 'tmaskutil')      # int8, unused
    del mask

    # f-point corrections on rows 1020/1021 (rows 1015:1022 of lat2/lon2 in xx2/yy2;
    # end-2 = row 1020 -> local 5, end-1 = row 1021 -> local 6). All float32.
    rows = slice(1014, 1022)                              # MATLAB 1015:1022
    xx2 = np.cos(lat2[rows, :] * PI / f32(180)) * np.cos(lon2[rows, :] * PI / f32(180))
    yy2 = np.cos(lat2[rows, :] * PI / f32(180)) * np.sin(lon2[rows, :] * PI / f32(180))

    # ----------------
    # error correction
    a1 = np.arange(379, 384)                              # MATLAB [380,381,382,383,384]
    a2 = np.arange(1058, 1063)                            # MATLAB [1059:1063]
    k = np.arange(5, dtype=f32)                           # MATLAB [0:4]

    def interp(arr, r, a):
        return arr[r, a[0]] + k * (arr[r, a[4]] - arr[r, a[0]]) / f32(4)

    for r in (5, 6):                                      # end-2 (%1020), end-1 (%1021)
        xx2[r, a1] = interp(xx2, r, a1)
        yy2[r, a1] = interp(yy2, r, a1)
        xx2[r, a2] = interp(xx2, r, a2)
        yy2[r, a2] = interp(yy2, r, a2)

    aa = np.concatenate([a1, a2])
    for r, row in ((5, 1019), (6, 1020)):                 # lon2(1020,:), lon2(1021,:)
        lonc = np.arctan2(yy2[r, aa], xx2[r, aa]) * f32(180) / PI
        latc = np.arccos(np.sqrt(yy2[r, aa] ** 2 + xx2[r, aa] ** 2)) * f32(180) / PI
        lon2[row, aa] = lonc
        lat2[row, aa] = latc
    # end error correction
    # ----------------

    # =========
    # Step 2. Read/adjust topo
    # =========
    # remove land features that are one cell in size (sequential: bathy is
    # modified while scanning, x outer, y inner as in the .m loops)
    b = bathy.tolist()
    for j in range(1, Nx - 1):                            # nnx=2:Nx-1
        for i in range(1, Ny - 2):                        # nny=2:Ny-2
            if (b[i][j] <= 0 and b[i][j + 1] > 0 and b[i][j - 1] > 0
                    and b[i + 1][j] > 0 and b[i - 1][j] > 0):
                b[i][j] = b[i][j + 1]
    bathy = np.array(b, dtype=bathy.dtype)
    del b

    bathy[664, 1125:1128] = bathy[663, 1125:1128]         # (665,1126:1128)=(664,..) Gibraltar
    bathy[687:690, 1261] = bathy[687:690, 1262]           # (688:690,1262)=(..,1263) Bosforus
    bathy[684, 1252:1257] = bathy[683, 1252:1257]         # (685,1253:1257)=(684,..) Dardanelles
    bathy[550, 1320] = bathy[550, 1321]                   # (551,1321)=(551,1322)
    bathy[549, 1321] = bathy[549, 1322]                   # (550,1322)=(550,1323) Red Sea
    bathy[960, 1000:1005] = bathy[961, 1000:1005]         # (961,1001:1005)=(962,..)
    bathy[962, 1009] = bathy[962, 1008]                   # (963,1010)=(963,1009)
    bathy[924, 886] = bathy[924, 887]                     # (925,887)=(925,888)
    bathy[923, 887] = bathy[923, 888]                     # (924,888)=(924,889)
    bathy[975, 939] = bathy[975, 938]                     # (976,940)=(976,939)
    bathy[972:975, 938] = bathy[972:975, 937]             # (973:975,939)=(..,938)
    bathy[974, 939] = bathy[974, 938]                     # (975,940)=(975,939) CAA
    bathy[751, 284] = bathy[751, 283]                     # (752,285)=(752,284) Sakhalin
    bathy[745:749, 283] = bathy[745:749, 282]             # (746:749,284)=(..,283)

    # ========
    # Step 3. Node indices, depth and node numbers:
    # ========
    # node_ind: 1 at the 4 corners of every land cell (order-independent)
    node_ind = np.zeros(lon2.shape)
    land = bathy[0:Ny - 1, 1:Nx] <= 0                     # bathy(nny,nnx+1)<=0, nnx=1:Nx-1, nny=1:Ny-1
    for di in (0, 1):
        for dj in (0, 1):
            node_ind[di:Ny - 1 + di, dj:Nx - 1 + dj][land] = 1

    # node_depth(nny,nnx)=max(max(bathy(nny-1:nny,nnx:nnx+1))), nnx=1:Nx-1, nny=2:Ny-1
    node_depth = np.zeros(lon2.shape)
    node_depth[1:Ny - 1, 0:Nx - 1] = np.maximum.reduce([
        bathy[0:Ny - 2, 0:Nx - 1], bathy[1:Ny - 1, 0:Nx - 1],
        bathy[0:Ny - 2, 1:Nx], bathy[1:Ny - 1, 1:Nx]])
    node_depth[0, :] = node_depth[1, :]
    node_depth[-1, :] = node_depth[-2, :]
    node_depth[:, Nx - 1] = node_depth[:, 1]              # (:,Nx)=(:,2)

    # node_num(nny,nnx)=nny+Ny*(nnx-1): MATLAB (1-based, column-major) ids
    ii, jj = np.meshgrid(np.arange(1, Ny + 1), np.arange(1, Nx + 1), indexing='ij')
    node_num = ii + Ny * (jj - 1)
    # Fold: one row of cells maps half onto the other; remove the extra nodes.
    node_num[Ny - 2, 1:721] = node_num[Ny - 3, 1440:720:-1]   # (end-1,2:721)=(end-2,1441:-1:722)
    node_num[Ny - 2, 0] = node_num[Ny - 3, 1]                 # (end-1,1)=(end-2,2)
    node_num[:, Nx - 2] = node_num[:, 0]                      # (:,Nx-1)=(:,1) cyclicity

    # ========
    # Step 4. Quadrilateral elements (sequential: defines the element numbering)
    # ========
    bl = bathy.tolist()
    nn = node_num.tolist()
    tri = []
    for j in range(1, Nx - 1):                            # nnx=2:Nx-1
        for i in range(0, Ny - 3):                        # nny=1:Ny-3
            if bl[i][j] > 0:
                tri.append((nn[i][j - 1], nn[i + 1][j - 1], nn[i + 1][j], nn[i][j]))
    i = Ny - 3                                            # nny=Ny-2
    for j in range(1, 721):                               # nnx=2:721
        if bl[i][j] > 0:
            tri.append((nn[i][j - 1], nn[i + 1][j - 1], nn[i + 1][j], nn[i][j]))
    e2d = len(tri)
    tri = np.array(tri, dtype=np.int64).T                 # (4, e2d), MATLAB ids

    # ===========
    # Step 5. Remove unreferenced nodes, create nodal arrays
    # ===========
    node_used = np.zeros(Nx * Ny, dtype=np.int64)
    node_used[tri.ravel() - 1] = 1
    n2d = int(node_used.sum())

    # column-major traversal nnx=1:Nx, nny=1:Ny == linear id order nny+Ny*(nnx-1)
    lon2f = lon2.ravel(order='F')
    lat2f = lat2.ravel(order='F')
    indf = node_ind.ravel(order='F')
    depf = node_depth.ravel(order='F')
    xcoord = np.zeros(n2d, dtype=f32)                     # grown from single in MATLAB
    ycoord = np.zeros(n2d, dtype=f32)
    indices = np.zeros(n2d)
    depth = np.zeros(n2d)
    used = node_used.tolist()
    counter = 0
    for lin in range(Nx * Ny):
        if used[lin] == 1:
            xcoord[counter] = lon2f[lin]
            ycoord[counter] = lat2f[lin]
            indices[counter] = indf[lin]
            depth[counter] = depf[lin]
            counter = counter + 1
    counter = 0
    for lin in range(Nx * Ny):
        if used[lin] == 1:
            counter = counter + 1
            used[lin] = counter
    node_used = np.array(used, dtype=np.int64)
    tri = node_used[tri - 1]

    # lev contains layer mid-depth; levels the depths of layer bottoms (float32
    # arithmetic: lev is single, levels a double array holding single values)
    levels = np.zeros(Nz + 1)
    levels[1] = f32(2) * lev[0]
    for nz in range(2, Nz + 1):                           # nz=3:Nz+1
        lp = f32(levels[nz - 1])
        levels[nz] = f32(2) * (lev[nz - 1] - lp) + lp

    # ===========
    # Step 6. Ascii output
    # ===========
    os.makedirs(out_dir, exist_ok=True)
    ids = range(1, n2d + 1)
    with open(os.path.join(out_dir, 'nod2d.out'), 'w') as fid:
        fid.write('%9i \n' % n2d)
        fid.write(''.join('%9i %8.4f %8.4f %8i\n' % (n, x, y, m) for n, x, y, m in
                          zip(ids, xcoord.tolist(), ycoord.tolist(), indices.astype(np.int64).tolist())))

    with open(os.path.join(out_dir, 'elem2d.out'), 'w') as fid:
        fid.write('%9i \n' % e2d)
        fid.write(''.join('%9i %9i %9i %9i\n' % tuple(q) for q in tri.T.tolist()))

    dint = depth.astype(np.int64)
    with open(os.path.join(out_dir, 'depth.out'), 'w') as fid:
        fid.write('%9i \n' % (Nz + 1))
        fid.write(''.join('%8.4f\n' % v for v in levels.tolist()))
        fid.write(''.join('%8.4f\n' % v for v in levels[dint].tolist()))   # depth of layers (m)
        fid.write(''.join('%9i \n' % v for v in dint.tolist()))            # number of layers

    print('nemo_meshq: %d nodes, %d quads, %d levels -> %s' % (n2d, e2d, Nz + 1, out_dir))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
