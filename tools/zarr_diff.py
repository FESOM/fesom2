#!/usr/bin/env python3
"""zarr_diff.py — read-back / compare tool for the FESOM3 hand-rolled Zarr writer (M9).

Modes (added as the M9 stages land):
  --roundtrip <store>   Stage 0: open the zarrsmoke store (codec none) with zarr AND xarray, assert
                        the toy arrays match their generator formula exactly.
  --lz4 <store>         Stage 0 Task 0.3: same value checks on the lz4+consolidated store — read via
                        consolidated metadata, assert the compressor is LZ4 (numcodecs framing OK).

Python gate env: /work/ab0995/a270088/mambaforge/bin/python3 (xarray/zarr/numpy).
Exit 0 = pass, 1 = fail.
"""
import argparse
import sys

import numpy as np
import zarr
import xarray as xr


def _fail(msg):
    print("FAIL:", msg)
    sys.exit(1)


# --- rotated-grid vector transform: an INDEPENDENT numpy port of mod_mesh_rotate (Task 2.5 gate) ----
# Reproduces vector_r2g so the geographic-frame vector output can be checked against a reference that
# does NOT share the Fortran code. Euler angles MUST match fesom_outputsmoke's read_mesh(...50,15,-90).
ROT_ALPHA, ROT_BETA, ROT_GAMMA = 50.0, 15.0, -90.0


def _r2g_matrix(alpha_deg, beta_deg, gamma_deg):
    al, be, ga = np.radians([alpha_deg, beta_deg, gamma_deg])
    M = np.empty((3, 3))
    M[0, 0] = np.cos(ga) * np.cos(al) - np.sin(ga) * np.cos(be) * np.sin(al)
    M[0, 1] = np.cos(ga) * np.sin(al) + np.sin(ga) * np.cos(be) * np.cos(al)
    M[0, 2] = np.sin(ga) * np.sin(be)
    M[1, 0] = -np.sin(ga) * np.cos(al) - np.cos(ga) * np.cos(be) * np.sin(al)
    M[1, 1] = -np.sin(ga) * np.sin(al) + np.cos(ga) * np.cos(be) * np.cos(al)
    M[1, 2] = np.cos(ga) * np.sin(be)
    M[2, 0] = np.sin(be) * np.sin(al)
    M[2, 1] = -np.sin(be) * np.cos(al)
    M[2, 2] = np.cos(be)
    return M


def _g2r(glon, glat, M):
    xg = np.cos(glat) * np.cos(glon); yg = np.cos(glat) * np.sin(glon); zg = np.sin(glat)
    xr = M[0, 0] * xg + M[0, 1] * yg + M[0, 2] * zg
    yr = M[1, 0] * xg + M[1, 1] * yg + M[1, 2] * zg
    zr = M[2, 0] * xg + M[2, 1] * yg + M[2, 2] * zg
    rlat = np.arcsin(zr)
    rlon = np.where((xr == 0) & (yr == 0), 0.0, np.arctan2(yr, xr))
    return rlon, rlat


def _vector_r2g(u, v, glon, glat, M):
    """rotated -> geographic vector (tlon=u zonal, tlat=v meridional). flag_coord=1 (geographic
    position): the same physical rotation the writer applies with flag_coord=0 at the rotated coords."""
    rlon, rlat = _g2r(glon, glat, M)
    txg = -v * np.sin(rlat) * np.cos(rlon) - u * np.sin(rlon)
    tyg = -v * np.sin(rlat) * np.sin(rlon) + u * np.cos(rlon)
    tzg = v * np.cos(rlat)
    txr = M[0, 0] * txg + M[1, 0] * tyg + M[2, 0] * tzg
    tyr = M[0, 1] * txg + M[1, 1] * tyg + M[2, 1] * tzg
    tzr = M[0, 2] * txg + M[1, 2] * tyg + M[2, 2] * tzg
    vout = -np.sin(glat) * np.cos(glon) * txr - np.sin(glat) * np.sin(glon) * tyr + np.cos(glat) * tzr
    uout = -np.sin(glon) * txr + np.cos(glon) * tyr
    return uout, vout


# generator formulas (must match src/drivers/fesom_zarrsmoke.F90)
EXP1 = np.array([(i + 1) * 1.5 for i in range(7)], dtype=np.float64)
EXP2 = np.array([[(i + 1) * 100 + (j + 1) for j in range(3)] for i in range(5)], dtype=np.float32)
EXP3 = np.array([[(i + 1) * 100 + (j + 1) for j in range(3)] for i in range(5)], dtype=np.int32)


def _check_values(g):
    a1 = g["arr1d_f8"][:]
    if a1.dtype != np.dtype("<f8"):
        _fail(f"arr1d_f8 dtype {a1.dtype} != <f8")
    d = np.max(np.abs(a1 - EXP1))
    print(f"  arr1d_f8  shape={a1.shape} dtype={a1.dtype}  max|Δ|={d:.3e}")
    if d != 0.0:
        _fail(f"arr1d_f8 max|Δ|={d} != 0")

    a2 = g["arr2d_f4"][:]
    if a2.dtype != np.dtype("<f4") or a2.shape != (5, 3):
        _fail(f"arr2d_f4 dtype/shape {a2.dtype}/{a2.shape}")
    d = np.max(np.abs(a2.astype(np.float64) - EXP2.astype(np.float64)))
    print(f"  arr2d_f4  shape={a2.shape} dtype={a2.dtype}  max|Δ|={d:.3e}")
    if d != 0.0:
        _fail(f"arr2d_f4 max|Δ|={d} != 0 (C-order transpose / partial-chunk / lz4-framing bug?)")

    a3 = g["arr2d_i4"][:]
    if a3.dtype != np.dtype("<i4") or not np.array_equal(a3, EXP3):
        _fail(f"arr2d_i4 mismatch:\n{a3}\n!=\n{EXP3}")
    print(f"  arr2d_i4  shape={a3.shape} dtype={a3.dtype}  exact match")


def roundtrip(store_path):
    print(f"[roundtrip none] store = {store_path}")
    g = zarr.open_group(store_path, mode="r")
    if g["arr1d_f8"].compressor is not None:
        _fail("none store: arr1d_f8 compressor should be null")
    _check_values(g)

    ds = xr.open_zarr(store_path, consolidated=False)
    for name, dims in (("arr1d_f8", ("x",)), ("arr2d_f4", ("d0", "d1")), ("arr2d_i4", ("d0", "d1"))):
        if name not in ds or tuple(ds[name].dims) != dims:
            _fail(f"xarray: {name} missing or wrong dims")
    print(f"  xarray open_zarr OK  vars={list(ds.data_vars)}  group attrs={dict(g.attrs)}")
    print("ROUNDTRIP PASS (max|Δ|=0)")


def roundtrip_lz4(store_path):
    print(f"[roundtrip lz4 + consolidated] store = {store_path}")
    # open via CONSOLIDATED metadata (.zmetadata) — exercises zarr_consolidate
    g = zarr.open_consolidated(store_path, mode="r")
    comp = g["arr1d_f8"].compressor
    if comp is None or comp.codec_id != "lz4":
        _fail(f"lz4 store: compressor is {comp}, expected LZ4")
    print(f"  compressor = {comp}")
    _check_values(g)  # decoding lz4 correctly proves the numcodecs framing

    ds = xr.open_zarr(store_path, consolidated=True)
    d = float(np.max(np.abs(ds["arr2d_f4"].values.astype(np.float64) - EXP2.astype(np.float64))))
    if d != 0.0:
        _fail("xarray (consolidated+lz4) arr2d_f4 values differ from generator")
    print(f"  xarray open_zarr(consolidated=True) OK  vars={list(ds.data_vars)}")
    print("ROUNDTRIP-LZ4 PASS (max|Δ|=0)")


# mesh.diag vars FESOM3 does NOT emit / can't be value-compared:
#   face_edges, face_links, gradient_vec_* : FESOM3 never builds elem_edges/elem_neighbors/gradient_vec
#   nod_part, elem_part                    : partition descriptors (differ by rank count, by design)
#   fesom_mesh                             : topology marker (dummy value)
#   ulevels*, zbar_*_surface               : cavity/partial-cell, OFF in this config
MESHDIAG_SKIP = {"fesom_mesh", "nod_part", "elem_part",
                 "face_edges", "face_links", "gradient_vec_x", "gradient_vec_y",
                 "ulevels", "ulevels_nod2D", "zbar_n_surface", "zbar_e_surface"}


def meshdiag(zarr_path, nc_path, ftol=0.0):
    """Compare the FESOM3 mesh.diag Zarr store vs the FESOM2 fesom.mesh.diag.nc, over the emitted
    subset, in canonical order (both are global-id ordered). Raw values (mask_and_scale=False):
    integer/connectivity exact; floats max|Δ| <= ftol (chase 0)."""
    print(f"[meshdiag] zarr = {zarr_path}\n           ref  = {nc_path}")
    f3 = xr.open_zarr(zarr_path, consolidated=True, mask_and_scale=False)
    f2 = xr.open_dataset(nc_path, mask_and_scale=False)

    nfail = 0
    checked = 0
    # iterate ALL variables (data_vars + coords like nz/nz1 which share their dim name)
    for v in sorted(f3.variables):
        if v in MESHDIAG_SKIP:
            continue
        if v not in f2.variables:
            print(f"  - {v:18s} SKIP (not in FESOM2 ref)")
            continue
        a3 = np.asarray(f3[v].values)
        a2 = np.asarray(f2[v].values)
        if a3.shape != a2.shape:
            print(f"  ! {v:18s} SHAPE {a3.shape} != {a2.shape}")
            nfail += 1
            continue
        checked += 1
        if np.issubdtype(a2.dtype, np.integer):
            nbad = int(np.count_nonzero(a3.astype(np.int64) != a2.astype(np.int64)))
            print(f"  {'ok ' if nbad==0 else 'BAD'} {v:18s} int   mismatches={nbad}")
            if nbad:
                nfail += 1
        else:
            d = float(np.max(np.abs(a3.astype(np.float64) - a2.astype(np.float64))))
            ok = d <= ftol
            print(f"  {'ok ' if ok else 'BAD'} {v:18s} float max|Δ|={d:.3e}")
            if not ok:
                nfail += 1

    print(f"[meshdiag] compared {checked} vars; {nfail} failures")
    if nfail:
        _fail(f"meshdiag: {nfail} var(s) mismatch")
    print("MESHDIAG PASS")


import os


def _find_store(out_dir, var):
    """Locate <var>.fesom.<year>.zarr in out_dir (the per-variable-per-year store)."""
    cands = sorted(d for d in os.listdir(out_dir)
                   if d.startswith(var + ".fesom.") and d.endswith(".zarr"))
    if not cands:
        _fail(f"no store for {var} in {out_dir}")
    return os.path.join(out_dir, cands[0])


def output(out_dir, frame="geographic"):
    """Stage-2 Task 2.1 round-trip: verify the fesom_outputsmoke per-variable-per-year stores.
    Each value at (record k, canonical node g, 1-based) must equal the driver's generator formula
    (partition-INDEPENDENT), with a CF time coord and embedded lon/lat. No FESOM2 oracle needed —
    self-consistency proves the writer (registry, growing time dim, append, decomp, coords).
    The fld_u/fld_v vector pair (Task 2.5) is checked in `frame`: native => raw formula (max|Δ|=0);
    geographic => the numpy r2g reference (tol, independent of the Fortran rotation code)."""
    print(f"[output] dir = {out_dir}  vec_frame = {frame}")
    # generator formulas — MUST match src/drivers/fesom_outputsmoke.F90
    #   fld_a (snap) = g + 0.25k ; fld_b (snap) = 0.5g - k ; fld_m (mean of g-1,g,g+1) = g
    def gen(v, nrec, N):
        gg = np.arange(1, N + 1, dtype=np.float64)
        kk = np.arange(nrec, dtype=np.float64)[:, None]
        if v == "fld_a":
            return gg[None, :] + 0.25 * kk
        if v == "fld_b":
            return gg[None, :] * 0.5 - kk
        return np.broadcast_to(gg[None, :], (nrec, N))  # fld_m

    coords = {}
    nrec_full = None
    for v in ("fld_a", "fld_b", "fld_m"):
        sp = _find_store(out_dir, v)
        # decode_times=False: keep `time` as raw seconds (else xarray -> datetime64; the fact that it
        # CAN decode already proves the CF units/calendar are valid).
        ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=False, decode_times=False)
        if v not in ds:
            _fail(f"{v} not in {sp}")
        da = ds[v]
        if tuple(da.dims) != ("time", "nod2"):
            _fail(f"{v} dims {da.dims} != (time, nod2)")
        nrec, N = da.shape
        nrec_full = nrec
        d = float(np.max(np.abs(da.values.astype(np.float64) - gen(v, nrec, N))))
        print(f"  {v:6s} shape=({nrec},{N}) dtype={da.dtype}  max|Δ|={d:.3e}")
        if d != 0.0:
            _fail(f"{v} values max|Δ|={d} != 0 (time-append / C-order / decomp bug?)")
        # CF time coord
        t = np.asarray(ds["time"].values, dtype=np.float64)
        dt = float(np.max(np.abs(t - 3600.0 * np.arange(nrec)))) if nrec else 0.0
        if dt != 0.0:
            _fail(f"{v} time max|Δ|={dt} != 0")
        tu = ds["time"].attrs.get("units", "")
        if "since" not in tu or "calendar" not in ds["time"].attrs:
            _fail(f"{v} time missing CF units/calendar (units={tu!r})")
        # embedded lon/lat coords
        for c in ("lon", "lat"):
            if c not in ds:
                _fail(f"{v} store missing {c}")
            arr = np.asarray(ds[c].values)
            if arr.shape != (N,) or not np.all(np.isfinite(arr)):
                _fail(f"{v} {c} shape={arr.shape}/finite check failed")
            coords.setdefault(c, []).append(arr)
        print(f"    time[0..]={t[:min(3,nrec)]}  units={tu!r} cal={ds['time'].attrs.get('calendar')}")
    for c in ("lon", "lat"):
        if not np.array_equal(coords[c][0], coords[c][1]):
            _fail(f"{c} differs between field stores (coord embed inconsistent)")

    # 3-D field fld_3: value(g, L) = g + (L+1) at valid levels; below-bottom masked -> NaN.
    sp = _find_store(out_dir, "fld_3")
    ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=True, decode_times=False)  # mask fill->NaN
    if "fld_3" not in ds:
        _fail(f"fld_3 not in {sp}")
    da = ds["fld_3"]
    if tuple(da.dims) != ("time", "nz1", "nod2"):
        _fail(f"fld_3 dims {da.dims} != (time, nz1, nod2)")
    nrec, nz1, N = da.shape
    vals = da.values.astype(np.float64)
    Lz = np.arange(1, nz1 + 1, dtype=np.float64)[None, :, None]  # 1-based layer
    gg = np.arange(1, N + 1, dtype=np.float64)[None, None, :]
    expect = np.broadcast_to(gg + Lz, (nrec, nz1, N))
    fin = np.isfinite(vals)
    if not fin.any():
        _fail("fld_3 fully masked (no valid levels)")
    nan_below = int(np.count_nonzero(~fin))
    d = float(np.max(np.abs(vals[fin] - expect[fin]))) if fin.any() else 0.0
    print(f"  fld_3  shape=({nrec},{nz1},{N}) dtype={da.dtype}  valid max|Δ|={d:.3e}  "
          f"masked(NaN)={nan_below}  (below-bottom CF mask)")
    if d != 0.0:
        _fail(f"fld_3 valid-level max|Δ|={d} != 0")
    if nan_below == 0:
        _fail("fld_3 has no NaN -> below-bottom masking did not fire")
    # vertical coord present, monotonic positive-down
    if "nz1" not in ds:
        _fail("fld_3 store missing nz1 vertical coord")
    zc = np.asarray(ds["nz1"].values)
    if zc.shape != (nz1,) or not np.all(np.diff(zc) > 0):
        _fail(f"nz1 coord shape/monotonic check failed (shape={zc.shape})")

    _output_vector(out_dir, frame)                                    # node velocity vector
    _output_vector(out_dir, frame, "fld_eu", "fld_ev", hdim="elem")   # Task 2.7 element vector
    _check_freq2(out_dir, nrec_full)
    _check_elem(out_dir)                                              # Task 2.7 element scalars
    print("OUTPUT PASS (max|Δ|=0)")


def _check_freq2(out_dir, nrec_full):
    """fld_f2 (Task 2.6): a 2-D snapshot on a freq=2 STEP cadence. It must have floor(nrec_full/2)
    records (proves the per-field step_event fired) and every value == gid+1000 (k-independent)."""
    sp = _find_store(out_dir, "fld_f2")
    ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=False, decode_times=False)
    if "fld_f2" not in ds:
        _fail(f"fld_f2 not in {sp}")
    da = ds["fld_f2"]
    nrec, N = da.shape
    exp = nrec_full // 2 if nrec_full else 0
    if nrec != exp:
        _fail(f"fld_f2 nrec={nrec} != floor(nrec_full/2)={exp} (per-field freq=2 step_event failed)")
    gg = np.arange(1, N + 1, dtype=np.float64) + 1000.0
    d = float(np.max(np.abs(da.values.astype(np.float64) - gg[None, :]))) if nrec else 0.0
    print(f"  fld_f2 shape=({nrec},{N})  records=floor({nrec_full}/2)={exp}  max|Δ|={d:.3e}  (per-field freq=2)")
    if d != 0.0:
        _fail(f"fld_f2 values max|Δ|={d} != 0")


def _output_vector(out_dir, frame, name_x="fld_u", name_y="fld_v", hdim="nod2", vtol=1e-9):
    """A 3-D vector pair (nl-1 layers), node OR element. native => store == raw generator (max|Δ|=0);
    geographic => store == numpy vector_r2g(raw) at the store's EMBEDDED lon/lat (node coords, or elem
    centroid for elements) — |Δ| <= vtol, an independent reference — AND the rotation is non-vacuous."""
    su = _find_store(out_dir, name_x); sv = _find_store(out_dir, name_y)
    du = xr.open_zarr(su, consolidated=False, mask_and_scale=True, decode_times=False)
    dv = xr.open_zarr(sv, consolidated=False, mask_and_scale=True, decode_times=False)
    if name_x not in du or name_y not in dv:
        _fail(f"{name_x}/{name_y} missing")
    u, v = du[name_x], dv[name_y]
    if tuple(u.dims) != ("time", "nz1", hdim):
        _fail(f"{name_x} dims {u.dims} != (time, nz1, {hdim})")
    nrec, nz1, N = u.shape
    uu = u.values.astype(np.float64); vv = v.values.astype(np.float64)
    # raw generator — MUST match src/drivers/fesom_outputsmoke.F90: u=0.001g+0.5L, v=-0.002g+0.25L+1
    g = np.arange(1, N + 1, dtype=np.float64)
    Lz = np.arange(1, nz1 + 1, dtype=np.float64)
    u_raw = 0.001 * g[None, :] + 0.5 * Lz[:, None]            # (nz1, N), 1-based g & L
    v_raw = -0.002 * g[None, :] + 0.25 * Lz[:, None] + 1.0
    if frame == "native":
        eu2, ev2, bar = u_raw, v_raw, 0.0
    else:
        lon = np.radians(np.asarray(du["lon"].values, dtype=np.float64))[None, :]   # (1, N) geographic
        lat = np.radians(np.asarray(du["lat"].values, dtype=np.float64))[None, :]
        eu2, ev2 = _vector_r2g(u_raw, v_raw, lon, lat, _r2g_matrix(ROT_ALPHA, ROT_BETA, ROT_GAMMA))
        bar = vtol
        chg = float(max(np.max(np.abs(eu2 - u_raw)), np.max(np.abs(ev2 - v_raw))))
        if chg < 1e-3:
            _fail(f"{name_x}/{name_y} geographic rotation vacuous (max|rot-raw|={chg:.3e}) — not applied?")
    eu = np.broadcast_to(eu2, (nrec, nz1, N)); ev = np.broadcast_to(ev2, (nrec, nz1, N))
    finu = np.isfinite(uu)
    if not finu.any():
        _fail(f"{name_x} fully masked (no valid levels)")
    dU = float(np.max(np.abs(uu[finu] - eu[finu])))
    dV = float(np.max(np.abs(vv[finu] - ev[finu])))
    nan_below = int(np.count_nonzero(~finu))
    print(f"  {name_x}/{name_y} frame={frame:10s} shape=({nrec},{nz1},{N})  max|Δu|={dU:.3e} "
          f"max|Δv|={dV:.3e}  masked(NaN)={nan_below}")
    if dU > bar or dV > bar:
        _fail(f"{name_x}/{name_y} frame={frame} max|Δ|=({dU:.3e},{dV:.3e}) > {bar:.0e} (rotation/pairing bug?)")
    if nan_below == 0:
        _fail(f"{name_x} no NaN -> below-bottom masking did not fire")


def _check_elem(out_dir):
    """Task 2.7 ELEMENT scalars: fld_e2 (2-D, value g+500 — integer-valued, float32-exact) + fld_e3
    (3-D on FULL levels nz, value g+100*L at valid levels, below-bottom NaN). g = canonical elem id."""
    sp = _find_store(out_dir, "fld_e2")
    ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=False, decode_times=False)
    da = ds["fld_e2"]
    if tuple(da.dims) != ("time", "elem"):
        _fail(f"fld_e2 dims {da.dims} != (time, elem)")
    nrec, N = da.shape
    g = np.arange(1, N + 1, dtype=np.float64)
    d = float(np.max(np.abs(da.values.astype(np.float64) - (g[None, :] + 500.0))))
    print(f"  fld_e2 shape=({nrec},{N}) dtype={da.dtype}  max|Δ|={d:.3e}  (element 2-D)")
    if d != 0.0:
        _fail(f"fld_e2 values max|Δ|={d} != 0 (element decomp / coords bug?)")
    # embedded elem-centroid lon/lat present + finite
    for c in ("lon", "lat"):
        if c not in ds or np.asarray(ds[c].values).shape != (N,) or not np.all(np.isfinite(ds[c].values)):
            _fail(f"fld_e2 missing/!finite element coord {c}")

    sp = _find_store(out_dir, "fld_e3")
    ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=True, decode_times=False)
    da = ds["fld_e3"]
    if tuple(da.dims) != ("time", "nz", "elem"):
        _fail(f"fld_e3 dims {da.dims} != (time, nz, elem)")
    nrec, nz, N = da.shape
    vals = da.values.astype(np.float64)
    Lz = np.arange(1, nz + 1, dtype=np.float64)[None, :, None]
    gg = np.arange(1, N + 1, dtype=np.float64)[None, None, :]
    expect = np.broadcast_to(gg + 100.0 * Lz, (nrec, nz, N))
    fin = np.isfinite(vals)
    d = float(np.max(np.abs(vals[fin] - expect[fin]))) if fin.any() else 0.0
    nan_below = int(np.count_nonzero(~fin))
    print(f"  fld_e3 shape=({nrec},{nz},{N})  valid max|Δ|={d:.3e}  masked(NaN)={nan_below}  "
          f"(element 3-D, full-levels mask)")
    if d != 0.0:
        _fail(f"fld_e3 valid max|Δ|={d} != 0")
    if nan_below == 0:
        _fail("fld_e3 no NaN -> element below-bottom masking did not fire")


def restart(ckpt_root):
    """Restart Stage 3 (Task 3.2) round-trip: verify a checkpoint folder written by mod_io_restart.
    ckpt_root is the dir holding the fesom.<YYYY>.<DDD>.<SSSSS>/ folder. Asserts the per-field SNAPSHOT
    store (eta_n.zarr) opens in xarray with the right dims (NO time dim), carries _ARRAY_DIMENSIONS,
    has embedded finite lon/lat sized to the entity, value(g)==g (canonical, partition-independent),
    and that checkpoint.json parses with the expected keys. Self-consistency — no FESOM2 oracle."""
    import json
    print(f"[restart] root = {ckpt_root}")
    folders = sorted(d for d in os.listdir(ckpt_root)
                     if d.startswith("fesom.") and os.path.isdir(os.path.join(ckpt_root, d)))
    if not folders:
        _fail(f"no fesom.<YYYY>.<DDD>.<SSSSS> checkpoint folder in {ckpt_root}")
    ck = os.path.join(ckpt_root, folders[-1])
    print(f"  checkpoint folder = {folders[-1]}")

    # --- per-field snapshot store (eta_n): a single-variable single-entity store, NO time dim ---
    sp = os.path.join(ck, "eta_n.zarr")
    if not os.path.isdir(sp):
        _fail(f"eta_n.zarr missing in {ck}")
    # raw _ARRAY_DIMENSIONS check straight off disk (the ushow/xarray contract)
    with open(os.path.join(sp, "eta_n", ".zattrs")) as fh:
        za = json.load(fh)
    if za.get("_ARRAY_DIMENSIONS") != ["nod2"]:
        _fail(f"eta_n _ARRAY_DIMENSIONS {za.get('_ARRAY_DIMENSIONS')} != ['nod2']")
    ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=False)
    if "eta_n" not in ds:
        _fail(f"eta_n not in {sp}")
    da = ds["eta_n"]
    if tuple(da.dims) != ("nod2",):
        _fail(f"eta_n dims {da.dims} != (nod2,) — a snapshot store must have NO time dim")
    (N,) = da.shape
    g = np.arange(1, N + 1, dtype=np.float64)
    d = float(np.max(np.abs(da.values.astype(np.float64) - g)))
    print(f"  eta_n  shape=({N},) dtype={da.dtype}  max|Δ|={d:.3e}  (value(g)==g)")
    if d != 0.0:
        _fail(f"eta_n values max|Δ|={d} != 0 (decomp / C-order / coord-embed bug?)")
    for c in ("lon", "lat"):
        if c not in ds:
            _fail(f"eta_n store missing {c}")
        arr = np.asarray(ds[c].values)
        if arr.shape != (N,) or not np.all(np.isfinite(arr)):
            _fail(f"{c} shape={arr.shape}/finite check failed (expected ({N},), all finite)")
    print(f"  lon/lat embedded: shape=({N},) finite=True")

    # --- checkpoint.json provenance manifest ---
    with open(os.path.join(ck, "checkpoint.json")) as fh:
        meta = json.load(fh)
    need = {"format_version", "year", "day", "time_sec", "globalstep", "fesom_git", "npes_wrote"}
    miss = need - set(meta)
    if miss:
        _fail(f"checkpoint.json missing keys {miss}")
    if meta["year"] != 2000 or meta["day"] != 1 or float(meta["time_sec"]) != 3600.0:
        _fail(f"checkpoint.json clock fields wrong: {meta}")
    print(f"  checkpoint.json OK: {meta}")
    print("RESTART PASS (max|Δ|=0, snapshot store opens + checkpoint.json parses)")


def restart_state(ckpt_root, ab_order=2, tke=True):
    """Restart Stage 3 (Task 3.4) full-state verify: a checkpoint written by restart_register_state must
    contain EVERY expected oce+ice store (incl. EVP sigma + the real(MP) hbar/hnode) with the right
    entity x level-kind shape, embedded finite lon/lat, finite data, element stores ~2x the node count,
    and each store's value matches its class formula (g for 2-D, g+0.5*L for 3-D; L 1-based) — so the
    MP->WP staging of hbar/hnode is proven lossless. Self-consistency — no FESOM2 oracle."""
    import json
    print(f"[restart-state] root = {ckpt_root}  ab_order={ab_order}  tke={tke}")
    folders = sorted(d for d in os.listdir(ckpt_root)
                     if d.startswith("fesom.") and os.path.isdir(os.path.join(ckpt_root, d))
                     and not d.endswith(".tmp"))
    if not folders:
        _fail(f"no fesom.<YYYY>.<DDD>.<SSSSS> checkpoint folder in {ckpt_root}")
    ck = os.path.join(ckpt_root, folders[-1])
    print(f"  checkpoint folder = {folders[-1]}")

    # expected stores: name -> (entity 'nod2'|'elem', kind '2d'|'nz1'|'nz')
    exp = {}
    for n in ("eta_n", "d_eta", "hbar", "ssh_rhs_old", "area", "hice", "hsnow", "uice", "vice", "t_skin"):
        exp[n] = ("nod2", "2d")
    exp["hnode"] = ("nod2", "nz1")
    for tr in ("temp", "salt"):
        exp[tr] = ("nod2", "nz1"); exp[tr + "_AB"] = ("nod2", "nz1"); exp[tr + "_M1"] = ("nod2", "nz1")
        if ab_order == 3:
            exp[tr + "_M2"] = ("nod2", "nz1")
    for n in ("w", "w_expl", "w_impl"):
        exp[n] = ("nod2", "nz")
    if tke:
        exp["tke"] = ("nod2", "nz")
    for n in ("u", "v", "urhs_AB", "vrhs_AB"):
        exp[n] = ("elem", "nz1")
    if ab_order == 3:
        exp["urhs_AB3"] = ("elem", "nz1"); exp["vrhs_AB3"] = ("elem", "nz1")
    for n in ("sigma11", "sigma12", "sigma22"):
        exp[n] = ("elem", "2d")

    on_disk = sorted(s[:-5] for s in os.listdir(ck) if s.endswith(".zarr"))
    miss = sorted(set(exp) - set(on_disk))
    extra = sorted(set(on_disk) - set(exp))
    if miss:
        _fail(f"checkpoint missing expected stores: {miss}")
    if extra:
        _fail(f"checkpoint has unexpected stores: {extra}")
    print(f"  {len(on_disk)} stores present, all expected ({len(exp)} required)")

    n_node = {}
    n_elem = {}
    nfail = 0
    for name in sorted(exp):
        entity, kind = exp[name]
        sp = os.path.join(ck, name + ".zarr")
        # raw _ARRAY_DIMENSIONS straight off disk (the ushow/xarray contract: (vdim, hdim) order)
        with open(os.path.join(sp, name, ".zattrs")) as fh:
            ad = json.load(fh).get("_ARRAY_DIMENSIONS")
        ds = xr.open_zarr(sp, consolidated=False, mask_and_scale=False)
        da = ds[name]
        if kind == "2d":
            want_dims = (entity,)
        else:
            want_dims = (kind, entity)
        if tuple(da.dims) != want_dims or list(ad) != list(want_dims):
            print(f"  ! {name:12s} dims {tuple(da.dims)} / _ARRAY_DIMENSIONS {ad} != {want_dims}")
            nfail += 1
            continue
        # value formula: g (2-D) or g+0.5*L (3-D), g & L 1-based canonical
        if kind == "2d":
            (N,) = da.shape
            g = np.arange(1, N + 1, dtype=np.float64)
            d = float(np.max(np.abs(da.values.astype(np.float64) - g)))
            nlev = None
        else:
            nlev, N = da.shape
            g = np.arange(1, N + 1, dtype=np.float64)[None, :]
            Lz = np.arange(1, nlev + 1, dtype=np.float64)[:, None]
            d = float(np.max(np.abs(da.values.astype(np.float64) - (g + 0.5 * Lz))))
        finite = bool(np.all(np.isfinite(da.values)))
        # embedded lon/lat sized to the entity, all finite
        coord_ok = True
        for c in ("lon", "lat"):
            arr = np.asarray(ds[c].values) if c in ds else np.array([])
            if arr.shape != (N,) or not np.all(np.isfinite(arr)):
                coord_ok = False
        # vertical coord present + monotonic positive-down for 3-D
        vert_ok = True
        if kind != "2d":
            if kind not in ds:
                vert_ok = False
            else:
                zc = np.asarray(ds[kind].values)
                vert_ok = (zc.shape == (nlev,) and bool(np.all(np.diff(zc) > 0)))
        ok = (d == 0.0) and finite and coord_ok and vert_ok
        shp = f"({N},)" if kind == "2d" else f"({nlev},{N})"
        print(f"  {'ok ' if ok else 'BAD'} {name:12s} {entity}/{kind:3s} shape={shp:13s} "
              f"max|Δ|={d:.1e} finite={finite} coords={coord_ok} vcoord={vert_ok}")
        if not ok:
            nfail += 1
        (n_node if entity == "nod2" else n_elem)[name] = N

    # partition-/entity-count sanity: element stores ~2x the node count (triangular mesh, as in M9)
    Nn = max(n_node.values()); Ne = max(n_elem.values())
    ratio = Ne / Nn
    print(f"  node count={Nn}  elem count={Ne}  elem/node={ratio:.3f}  (expect ~2x)")
    if not (1.7 <= ratio <= 2.3):
        _fail(f"elem/node ratio {ratio:.3f} outside [1.7, 2.3] (entity decomp wrong?)")
    if len(set(n_node.values())) != 1:
        _fail(f"node stores disagree on horizontal size: {sorted(set(n_node.values()))}")
    if len(set(n_elem.values())) != 1:
        _fail(f"element stores disagree on horizontal size: {sorted(set(n_elem.values()))}")

    if nfail:
        _fail(f"restart-state: {nfail} store(s) failed shape/value/coord checks")
    print("RESTART-STATE PASS (all stores present, shapes+levels+values+coords OK, elem~2x node)")


def output_cmp(d1, d2, rel_floor=0.0):
    """Partition-independence / restart self-consistency: two output dirs hold value-identical
    stores. Default rel_floor=0.0 => strict max|Δ|=0 (np=1, and the exact-by-construction cases).

    rel_floor>0 admits a per-store RELATIVE machine-epsilon floor for the np>1 restart gate with
    FESOM3_SHARED_SYNC=0 ONLY: at np>1, FESOM2's any-node element ownership makes 562 boundary
    elements (CORE2 dist_2) redundantly owned; without the sync their uv copies diverge by ~1 ULP
    from order-dependent RHS sums, and the partition-INDEPENDENT checkpoint must dedup them to one
    value, so a same-np resume perturbs the other rank by ~1 ULP and propagates ~2e-16 into temp/w
    (FESOM2 behaves the same). With the sync (the default; docs/plans/completed/20261012-shared-element-sync.md)
    every owner takes the designated owner's copy, the one the writer keeps, and np>1 restarts are
    exact (floor 0). The floor is set far below any real regression (the t_skin/d_eta seeds were
    ~1e-4 >> 1e-12) so the gate stays meaningful. Stores are classified =0 (exact) / ~eps
    (<=floor) / BAD (>floor); a BAD at >floor still FAILS."""
    print(f"[output-cmp] {d1}  vs  {d2}" + (f"  (np>1 rel-eps floor={rel_floor:.0e})" if rel_floor else ""))
    stores = sorted(s for s in os.listdir(d1) if s.endswith(".zarr"))
    if not stores:
        _fail(f"no .zarr stores in {d1}")
    nfail = 0
    neps = 0
    for s in stores:
        ds1 = xr.open_zarr(os.path.join(d1, s), consolidated=False, mask_and_scale=False, decode_times=False)
        ds2 = xr.open_zarr(os.path.join(d2, s), consolidated=False, mask_and_scale=False, decode_times=False)
        for v in sorted(ds1.variables):
            a1 = np.asarray(ds1[v].values)
            a2 = np.asarray(ds2[v].values)
            if a1.shape != a2.shape:
                print(f"  ! {s}:{v} shape {a1.shape} != {a2.shape}")
                nfail += 1
                continue
            a1 = a1.astype(np.float64); a2 = a2.astype(np.float64)
            d = float(np.max(np.abs(a1 - a2))) if a1.size else 0.0
            scale = float(np.max(np.abs(a1))) if a1.size else 0.0
            # relative residual vs the store's own magnitude; rel=d if scale==0 (forces exact for all-zero)
            rel = d / scale if scale > 0.0 else d
            if d == 0.0:
                tag = "ok "
            elif rel_floor > 0.0 and rel <= rel_floor:
                tag = "eps"; neps += 1
            else:
                tag = "BAD"; nfail += 1
            print(f"  {tag} {s}:{v:8s} max|Δ|={d:.3e}  rel={rel:.2e}")
    if nfail:
        _fail(f"output-cmp: {nfail} mismatch(es)" + (f" (+{neps} at <=eps floor)" if neps else ""))
    if neps:
        print(f"OUTPUT-CMP PASS ({neps} store(s) at <=rel-eps floor {rel_floor:.0e} — np>1 "
              f"redundant-element roundoff without the shared-element sync; rest max|Δ|=0)")
    else:
        print("OUTPUT-CMP PASS (max|Δ|=0, partition-independent)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--roundtrip", metavar="STORE", help="Stage-0 zarrsmoke round-trip (codec none)")
    ap.add_argument("--lz4", metavar="STORE", help="Stage-0 lz4 + consolidated round-trip")
    ap.add_argument("--meshdiag", nargs=2, metavar=("ZARR", "NC"),
                    help="compare mesh.diag Zarr store vs FESOM2 fesom.mesh.diag.nc")
    ap.add_argument("--output", metavar="DIR", help="Stage-2 fesom_outputsmoke store verify (formula)")
    ap.add_argument("--restart", metavar="DIR", help="Stage-3 restart checkpoint verify (snapshot + json)")
    ap.add_argument("--restart-state", metavar="DIR", dest="restart_state",
                    help="Stage-3 Task 3.4 full-state checkpoint verify (all oce+ice stores)")
    ap.add_argument("--ab-order", type=int, default=2, dest="ab_order",
                    help="Adams-Bashforth order for --restart-state expected set (2 or 3; default 2)")
    ap.add_argument("--no-tke", action="store_true", dest="no_tke",
                    help="--restart-state: do NOT expect the optional tke store")
    ap.add_argument("--output-cmp", nargs=2, metavar=("DIR1", "DIR2"), dest="output_cmp",
                    help="compare two output dirs (partition-independence)")
    ap.add_argument("--frame", choices=("geographic", "native"), default="geographic",
                    help="vector frame for --output's fld_u/fld_v check (default geographic)")
    ap.add_argument("--ftol", type=float, default=0.0, help="float tolerance for --meshdiag (default 0)")
    ap.add_argument("--rel-floor", type=float, default=0.0, dest="rel_floor",
                    help="per-store RELATIVE eps floor for --output-cmp (default 0=strict; the np>1 "
                         "restart gate uses ~1e-12 only with FESOM3_SHARED_SYNC=0, for the "
                         "redundant-element roundoff of unsynced copies)")
    args = ap.parse_args()

    if args.roundtrip:
        roundtrip(args.roundtrip)
    elif args.lz4:
        roundtrip_lz4(args.lz4)
    elif args.meshdiag:
        meshdiag(args.meshdiag[0], args.meshdiag[1], ftol=args.ftol)
    elif args.output:
        output(args.output, frame=args.frame)
    elif args.restart:
        restart(args.restart)
    elif args.restart_state:
        restart_state(args.restart_state, ab_order=args.ab_order, tke=not args.no_tke)
    elif args.output_cmp:
        output_cmp(args.output_cmp[0], args.output_cmp[1], rel_floor=args.rel_floor)
    else:
        ap.error("no mode selected (use --roundtrip / --lz4 / --meshdiag / --output / --output-cmp / "
                 "--restart / --restart-state)")


if __name__ == "__main__":
    main()
