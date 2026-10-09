#!/usr/bin/env python3
"""Diagnostics of the short runs on the cut-out meshes (hybrid-mesh plan Task 17,
tools/run_cutout_short.sh). Needs numpy, scipy and zarr (module load python3).

  cutout_diag.py selftest
  cutout_diag.py compare <ref_restart_dir> <new_restart_dir> <dt> <early_steps> <tol_early>
                         <tol_ocean> <tol_mix> <tol_ice> [csv]
      1 vs N: every checkpoint folder fesom.<YYYY>.<DDD>.<SSSSS> present in both restart
      directories, every field store in it: d = max|a - b| / max|a| (a = reference; a field that is
      0 everywhere in a is compared absolutely, d = max|b|). Two tiers: EARLY = every field of the
      checkpoints within the first <early_steps> steps against tol_early; the whole run in the
      classes MIX = tke (threshold-driven), ICE = area hice hsnow uice vice t_skin sigma11/12/22,
      OCEAN = everything else against tol_ocean / tol_mix / tol_ice. Prints a table field x
      checkpoint, the per-field maximum, and exits 1 if a tier exceeds its tolerance. csv: optional
      file for the full table.
  cutout_diag.py checker <mesh_dir> <output_dir> <var> [<var> ...] [--csv file]
      time series of the checkerboard and interface-noise measures (below) of the 2-D output
      fields <var> (zarr stores <output_dir>/<var>.fesom.<YYYY>.zarr, one record per step).

CHECKERBOARD MEASURES of a node field f (all weights are areas; element areas planar in the local
lon/lat tangent plane, A_n^Q = sum of A_e/4 over the quads at node n, so sum_n A_n^Q = sum_Q A_e):
  V_Q    the vertices of the quads (all nodes on the NEMO cut-out; the quad patches incl. their
         triangle-quad interface nodes on the mixed one); var = sum_{V_Q} A_n^Q (f - fbar)^2 with
         fbar the A^Q-weighted mean over V_Q.
  E_loc  = sum_Q A_e q_e^2 / var,  q_e = (f1 - f2 + f3 - f4)/4 -- each quad's projection onto its
         own (+1,-1,+1,-1) pattern. Sign-free, so it also sees a checkerboard whose amplitude or
         phase varies in space; 0 for a field linear on parallelograms, 1 for a pure checkerboard.
  E_glob = sum_C p_C^2 sum_C A^Q c'^2 / var: projection onto the GLOBAL 2-colouring c = +-1 of the
         quad-side graph (BFS per connected component C; c' = c minus its A^Q-mean on C, p_C =
         sum_C A^Q c' f / sum_C A^Q c'^2). This is the exact null mode of the SSH stiffness on an
         all-quad mesh (Task 9). 1 for a pure checkerboard. A component with an odd cycle (e.g.
         an odd ring of quads around an island) has no 2-colouring: it is left out and counted.
INTERFACE NOISE (Task 11: a linear tracer is not advected exactly at triangle-quad interface nodes):
  r_n = f_n - P_n(x_n), P_n the least-squares plane through the values at the edge neighbours of n
  (local tangent coordinates): exact (r = 0) for a linear field on ANY neighbourhood, so the
  irregular stencil at an interface node is not mistaken for noise; O(h^2 f'') for a smooth field,
  2 for a checkerboard. At interior nodes (not on a boundary edge), split into the classes tri (all
  elements triangles), quad (all quads) and iface (both); R_class = rms_class(r) / std(f), std
  over all nodes (area-weighted).
"""
import glob
import os
import sys

import numpy as np
import scipy.sparse
import scipy.sparse.csgraph

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "mesh"))
import mesh_io  # noqa: E402

MIX = {"tke"}
ICE = {"area", "hice", "hsnow", "uice", "vice", "t_skin", "sigma11", "sigma12", "sigma22"}
COORDS = {"lon", "lat", "nz", "nz1", "time"}


# ----------------------------------------------------------------------------- geometry
class Geo:
    """Element areas, node shares, quad set, quad-side colouring and node classes of a mesh
    (lon/lat in degrees, elem (e2d, 4) 0-based, triangle 4th == 1st)."""

    def __init__(self, lon, lat, elem):
        self.n2d = len(lon)
        self.elem = elem
        nv = np.where(elem[:, 3] == elem[:, 0], 3, 4)
        self.quad = nv == 4
        # planar element areas in the local tangent plane (cyclic lon wrap)
        lo, la = np.radians(lon)[elem], np.radians(lat)[elem]
        lo = lo - lo[:, :1]
        lo = (lo + np.pi) % (2 * np.pi) - np.pi
        x = mesh_io.R_EARTH * lo * np.cos(la.mean(axis=1))[:, None]
        y = mesh_io.R_EARTH * (la - la[:, :1])
        a = 0.5 * np.abs(np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1))
        self.area_e = a
        # node shares: all elements (A_e/nv over the unique vertices) and quads only (A_e/4)
        self.area_n = np.zeros(self.n2d)
        for k in range(4):
            uniq = (k < 3) | self.quad
            np.add.at(self.area_n, elem[uniq, k], a[uniq] / nv[uniq])
        q = elem[self.quad]
        self.area_q = np.zeros(self.n2d)
        for k in range(4):
            np.add.at(self.area_q, q[:, k], self.area_e[self.quad] / 4.0)
        self.vq = self.area_q > 0
        # quad-side graph -> 2-colouring per connected component
        a_, b_ = q.ravel(), np.roll(q, -1, axis=1).ravel()
        g = scipy.sparse.coo_matrix((np.ones(a_.size), (a_, b_)), shape=(self.n2d, self.n2d)).tocsr()
        g = ((g + g.T) > 0).astype(np.int8)
        ncomp, lab = scipy.sparse.csgraph.connected_components(g, directed=False)
        self.colour = np.zeros(self.n2d)
        self.comp = np.full(self.n2d, -1)
        self.odd_components = 0
        self.components = 0
        for c in range(ncomp):
            nodes = np.nonzero((lab == c) & self.vq)[0]
            if nodes.size == 0:
                continue
            self.components += 1
            order, pred = scipy.sparse.csgraph.breadth_first_order(g, nodes[0], directed=False)
            depth = np.zeros(self.n2d, dtype=np.int64)
            for v in order[1:]:          # BFS order: the predecessor's depth is already set
                depth[v] = depth[pred[v]] + 1
            col = np.where(depth[order] % 2 == 0, 1.0, -1.0)
            self.colour[order] = col
            # bipartite iff every quad side joins opposite colours
            gi, gj = g[order].nonzero()
            if np.any(self.colour[order[gi]] == self.colour[gj]):
                self.odd_components += 1
                self.colour[order] = 0.0
                continue
            self.comp[order] = c
        # node classes and edge-neighbour graph
        ed, _, count = mesh_io.side_edges(elem, self.n2d)
        bnd = mesh_io.boundary_mask(self.n2d, ed, count)
        has_tri = np.zeros(self.n2d, dtype=bool)
        has_quad = np.zeros(self.n2d, dtype=bool)
        for k in range(4):
            has_tri[elem[~self.quad, k]] = True
            has_quad[elem[self.quad, k]] = True
        self.cls = {"tri": has_tri & ~has_quad & ~bnd, "quad": has_quad & ~has_tri & ~bnd,
                    "iface": has_tri & has_quad & ~bnd}
        # linear-exact neighbour prediction: row n of W evaluates the least-squares plane through
        # the edge neighbours of n at x_n (interior nodes only; other rows stay 0)
        e1, e2 = ed[:, 0], ed[:, 1]
        adj = scipy.sparse.coo_matrix((np.ones(2 * len(ed)), (np.r_[e1, e2], np.r_[e2, e1])),
                                      shape=(self.n2d, self.n2d)).tocsr()
        rlon, rlat = np.radians(lon), np.radians(lat)
        wr, wc, wv = [], [], []
        for n in np.nonzero(~bnd)[0]:
            nbr = adj.indices[adj.indptr[n]:adj.indptr[n + 1]]
            dl = (rlon[nbr] - rlon[n] + np.pi) % (2 * np.pi) - np.pi
            X = np.column_stack([np.ones(nbr.size), dl * np.cos(rlat[n]), rlat[nbr] - rlat[n]])
            wr += [n] * nbr.size
            wc += nbr.tolist()
            wv += np.linalg.pinv(X)[0].tolist()
        self.W = scipy.sparse.coo_matrix((wv, (wr, wc)), shape=(self.n2d, self.n2d)).tocsr()

    def measures(self, f):
        """{"E_loc", "E_glob", "R_tri", "R_quad", "R_iface"} of the node field f (NaN if undefined)."""
        out = {}
        w = self.area_q[self.vq]
        fq = f[self.vq]
        var = np.sum(w * (fq - np.sum(w * fq) / np.sum(w)) ** 2)
        q = self.elem[self.quad]
        qe = (f[q[:, 0]] - f[q[:, 1]] + f[q[:, 2]] - f[q[:, 3]]) / 4.0
        out["E_loc"] = np.sum(self.area_e[self.quad] * qe ** 2) / var if var > 0 else np.nan
        eg = 0.0
        for c in np.unique(self.comp[self.comp >= 0]):
            m = self.comp == c
            wa = self.area_q[m]
            cc = self.colour[m] - np.sum(wa * self.colour[m]) / np.sum(wa)
            den = np.sum(wa * cc ** 2)
            if den > 0:
                eg += np.sum(wa * cc * f[m]) ** 2 / den
        out["E_glob"] = eg / var if var > 0 else np.nan
        sd = np.sqrt(np.sum(self.area_n * (f - np.sum(self.area_n * f) / np.sum(self.area_n)) ** 2)
                     / np.sum(self.area_n))
        r = f - self.W @ f
        for k, m in self.cls.items():
            out["R_" + k] = np.sqrt(np.mean(r[m] ** 2)) / sd if (m.any() and sd > 0) else np.nan
        return out


# ----------------------------------------------------------------------------- compare
def ckpt_folders(d):
    return {os.path.basename(p): p for p in glob.glob(os.path.join(d, "fesom.*.*.*"))
            if os.path.isdir(p) and not p.endswith(".tmp")}


def read_store(p):
    import zarr
    g = zarr.open(p, mode="r")
    k = [n for n in g.array_keys() if n not in COORDS]
    return np.asarray(g[k[0]][:], dtype=np.float64)


def compare(dref, dnew, dt, early_steps, tol_early, tol_ocean, tol_mix, tol_ice, csv=None):
    fa, fb = ckpt_folders(dref), ckpt_folders(dnew)
    tags = sorted(set(fa) & set(fb))
    if not tags:
        print(f"  no common checkpoints in {dref} and {dnew}")
        return 1
    t0 = None
    res = {}                                     # field -> [d per tag]
    hours = []
    for tag in tags:
        _, y, day, sec = tag.split(".")
        t = (int(day) - 1) * 86400 + int(sec)
        t0 = t if t0 is None else t0
        hours.append(t)
        names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(fa[tag], "*.zarr")))
        for nm in names:
            pb = os.path.join(fb[tag], nm + ".zarr")
            if not os.path.isdir(pb):
                print(f"  {tag}: {nm} missing in {dnew}")
                return 1
            a, b = read_store(os.path.join(fa[tag], nm + ".zarr")), read_store(pb)
            if a.shape != b.shape:
                print(f"  {tag}: {nm} shapes differ {a.shape} {b.shape}")
                return 1
            ma = np.max(np.abs(a))
            d = float(np.max(np.abs(b - a)) / ma) if ma > 0 else float(np.max(np.abs(b)))
            res.setdefault(nm, []).append(d)
    steps = [int(round(h / dt)) for h in hours]
    show = list(range(len(tags)))
    if len(show) > 10:                           # at most 10 columns on screen (all in the csv)
        show = sorted(set(np.linspace(0, len(tags) - 1, 10).round().astype(int)))
    every = steps[1] - steps[0] if len(steps) > 1 else steps[0]
    print(f"  d = max|N - 1| / max|1| per field and checkpoint (step of the run; ckpt every {every} steps)")
    print(f"  {'step':<12}" + "".join(f"{steps[i]:9d}" for i in show) + "      max  class")
    worst = {"OCEAN": 0.0, "MIX": 0.0, "ICE": 0.0}
    for nm in sorted(res, key=lambda k: -max(res[k])):
        cl = "MIX" if nm in MIX else "ICE" if nm in ICE else "OCEAN"
        worst[cl] = max(worst[cl], max(res[nm]))
        print(f"  {nm:<12}" + "".join(f"{res[nm][i]:9.1e}" for i in show)
              + f"  {max(res[nm]):9.2e}  {cl}")
    if csv:
        with open(csv, "w") as f:
            f.write("field," + ",".join(str(s) for s in steps) + "\n")
            for nm in sorted(res):
                f.write(nm + "," + ",".join(f"{v:.3e}" for v in res[nm]) + "\n")
    tol = {"OCEAN": tol_ocean, "MIX": tol_mix, "ICE": tol_ice}
    early = [i for i, s in enumerate(steps) if s <= early_steps]
    d_early = max([res[nm][i] for nm in res for i in early], default=0.0)
    fail = d_early > tol_early
    print(f"  EARLY (steps <= {early_steps}, all fields) max d {d_early:.3e}  (TOL {tol_early:g})"
          + ("   ABOVE" if fail else ""))
    for cl in ("OCEAN", "MIX", "ICE"):
        bad = worst[cl] > tol[cl]
        fail |= bad
        print(f"  {cl:<5} max d {worst[cl]:.3e}  (TOL {tol[cl]:g})" + ("   ABOVE" if bad else ""))
    return 1 if fail else 0


# ----------------------------------------------------------------------------- checker
def checker(mesh_dir, out_dir, names, csv=None):
    import zarr
    m = mesh_io.read_mesh(mesh_dir)
    geo = Geo(m.lon, m.lat, m.elem)
    name = os.path.basename(mesh_dir.rstrip("/"))
    ntri, nquad, nif = (int(geo.cls[k].sum()) for k in ("tri", "quad", "iface"))
    print(f"  mesh {name}: {m.n2d} nodes, {int(geo.quad.sum())} quads, {int((~geo.quad).sum())} triangles; "
          f"quad-side graph: {geo.components} component(s), {geo.odd_components} without a 2-colouring "
          f"(left out of E_glob); interior nodes tri/quad/iface = {ntri}/{nquad}/{nif}")
    keys = ("E_loc", "E_glob", "R_tri", "R_quad", "R_iface")
    rows = []
    summary = {}
    for nm in names:
        stores = sorted(glob.glob(os.path.join(out_dir, nm + ".fesom.*.zarr")))
        if not stores:
            print(f"  no output store for {nm} in {out_dir}")
            return 1
        recs, times = [], []
        for p in stores:
            g = zarr.open(p, mode="r")
            recs.append(np.asarray(g[nm][:], dtype=np.float64))
            times.append(np.asarray(g["time"][:], dtype=np.float64))
        x, t = np.concatenate(recs), np.concatenate(times)
        ser = {k: np.full(len(t), np.nan) for k in keys}
        for i in range(len(t)):
            for k, v in geo.measures(x[i]).items():
                ser[k][i] = v
        for i in range(len(t)):
            rows.append((nm, t[i], [ser[k][i] for k in keys]))
        summary[nm] = (t, ser)
    # on screen: every 6 h of model time, first and last record
    for nm, (t, ser) in summary.items():
        print(f"  {nm}: time series (h, every 6 h)")
        print(f"    {'hour':>8}" + "".join(f"{k:>11}" for k in keys))
        h = t / 3600.0
        pick = [i for i in range(len(t)) if i == 0 or i == len(t) - 1 or abs(h[i] / 6 - round(h[i] / 6)) < 1e-9]
        for i in pick:
            print(f"    {h[i]:8.1f}" + "".join(f"{ser[k][i]:11.3e}" for k in keys))
        n = len(t)
        day = max(1, int(round(86400.0 / (t[1] - t[0])))) if n > 1 else 1
        for k in ("E_loc", "E_glob", "R_iface", "R_quad"):
            s = ser[k]
            if np.all(np.isnan(s)):
                continue
            first, last = np.nanmean(s[:day]), np.nanmean(s[-day:])
            ratio = last / first if first > 0 else np.nan
            print(f"    {nm} {k:<7} first-day mean {first:.3e}  last-day mean {last:.3e}  ratio {ratio:.2f}"
                  f"  max {np.nanmax(s):.3e}")
    if csv:
        with open(csv, "w") as f:
            f.write("var,time_s," + ",".join(keys) + "\n")
            for nm, tt, v in rows:
                f.write(f"{nm},{tt:.1f}," + ",".join(f"{z:.4e}" for z in v) + "\n")
    return 0


# ----------------------------------------------------------------------------- selftest
def grid_mesh(ni, nj, split=None, dx=0.25):
    """Structured lon/lat quad grid (ni x nj nodes, clockwise quads SW-NW-NE-SE as nemo_meshq);
    split: boolean (ni-1, nj-1) mask of quads cut into two triangles."""
    lon, lat = np.meshgrid(np.arange(ni) * dx - 10.0, np.arange(nj) * dx + 2.0, indexing="ij")
    idx = np.arange(ni * nj).reshape(ni, nj)
    el = []
    for i in range(ni - 1):
        for j in range(nj - 1):
            sw, nw, ne, se = idx[i, j], idx[i, j + 1], idx[i + 1, j + 1], idx[i + 1, j]
            if split is not None and split[i, j]:
                el += [[sw, nw, ne, sw], [sw, ne, se, sw]]
            else:
                el.append([sw, nw, ne, se])
    return lon.ravel(), lat.ravel(), np.array(el, dtype=np.int64), idx


def selftest():
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(f"  {msg:<70} " + ("ok" if cond else "FAILED"))
        ok &= bool(cond)

    lon, lat, el, idx = grid_mesh(24, 18)
    g = Geo(lon, lat, el)
    ii, jj = np.meshgrid(np.arange(24), np.arange(18), indexing="ij")
    cb = np.where((ii + jj) % 2 == 0, 1.0, -1.0).ravel()
    mm = g.measures(cb)
    check(abs(mm["E_loc"] - 1) < 1e-2 and abs(mm["E_glob"] - 1) < 1e-12,
          f"pure checkerboard: E_loc {mm['E_loc']:.4f}, E_glob {mm['E_glob']:.15f}")
    lin = 0.3 + 2.0 * lon - 1.5 * lat
    mm = g.measures(lin)
    check(mm["E_loc"] < 1e-24 and mm["E_glob"] < 1e-3,
          f"linear field: E_loc {mm['E_loc']:.1e}, E_glob {mm['E_glob']:.1e}")
    smooth = np.sin(np.radians(lon) * 40) * np.cos(np.radians(lat) * 30)
    e = [g.measures(smooth + eps * cb)["E_glob"] for eps in (0.0, 1e-3, 1e-2)]
    check(e[0] < e[1] < e[2] and abs(e[2] / e[1] - 100) < 5,
          f"smooth + eps*checkerboard: E_glob grows as eps^2 ({e[0]:.1e} {e[1]:.1e} {e[2]:.1e})")
    mod = cb * (1 + 0.5 * np.sin(np.radians(lon) * 20))          # modulated checkerboard
    mm = g.measures(smooth + 1e-2 * mod)
    check(mm["E_loc"] > 0.5 * g.measures(smooth + 1e-2 * cb)["E_loc"],
          f"modulated checkerboard seen by E_loc ({mm['E_loc']:.2e})")
    rr = g.measures(cb)
    check(rr["R_quad"] > 1.0 and np.isnan(rr["R_tri"]) and np.isnan(rr["R_iface"]),
          f"all-quad grid: node classes quad only, checkerboard roughness {rr['R_quad']:.2f}")
    # mixed grid: the left half split into triangles -> quad patch on the right, interface column
    split = np.zeros((23, 17), dtype=bool)
    split[:11, :] = True
    lon, lat, el, idx = grid_mesh(24, 18, split)
    g = Geo(lon, lat, el)
    check(g.cls["iface"].sum() == 16 and g.cls["tri"].sum() > 0 and g.cls["quad"].sum() > 0,
          f"mixed grid: {g.cls['iface'].sum()} interface nodes (16 expected)")
    cbq = np.where(g.vq, cb, 0.0)
    mm = g.measures(cbq)
    check(abs(mm["E_glob"] - 1) < 1e-12 and g.components == 1 and g.odd_components == 0,
          f"mixed grid: checkerboard on the quad patch E_glob {mm['E_glob']:.15f}")
    mm = g.measures(0.3 + 2.0 * lon - 1.5 * lat)
    check(mm["E_loc"] < 1e-24 and max(mm["R_tri"], mm["R_quad"], mm["R_iface"]) < 1e-12,
          f"mixed grid: linear field E_loc {mm['E_loc']:.1e}, "
          f"R tri/quad/iface {mm['R_tri']:.1e} {mm['R_quad']:.1e} {mm['R_iface']:.1e}")
    q2 = g.measures(smooth)
    check(q2["R_iface"] < 0.1 and q2["R_quad"] < 0.1,
          f"mixed grid: smooth field R quad/iface {q2['R_quad']:.1e} {q2['R_iface']:.1e} (O(h^2))")
    # odd ring of 7 quads around a hole: no 2-colouring -> left out of E_glob
    ang = np.arange(7) * 2 * np.pi / 7
    lon = np.r_[np.cos(ang), 2 * np.cos(ang)]
    lat = np.r_[np.sin(ang), 2 * np.sin(ang)]
    el = np.array([[k, 7 + k, 7 + (k + 1) % 7, (k + 1) % 7] for k in range(7)], dtype=np.int64)
    g = Geo(lon, lat, el)
    mm = g.measures(np.arange(14, dtype=float))
    check(g.odd_components == 1 and mm["E_glob"] == 0.0,
          f"odd ring of 7 quads: {g.odd_components} component without 2-colouring, E_glob 0")
    print("cutout_diag selftest: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main(argv):
    if len(argv) >= 1 and argv[0] == "selftest":
        return selftest()
    if len(argv) >= 9 and argv[0] == "compare":
        return compare(argv[1], argv[2], float(argv[3]), int(argv[4]), float(argv[5]),
                       float(argv[6]), float(argv[7]), float(argv[8]), argv[9] if len(argv) > 9 else None)
    if len(argv) >= 4 and argv[0] == "checker":
        csv = None
        if "--csv" in argv:
            k = argv.index("--csv")
            csv = argv[k + 1]
            argv = argv[:k] + argv[k + 2:]
        return checker(argv[1], argv[2], argv[3:], csv)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
