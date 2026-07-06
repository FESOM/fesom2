# Double-precision mesh-metric computation for single-precision FESOM (plan)

**Status:** design / not yet implemented. Companion to the SSH-CG mixed-precision fix on this
branch (`b1c39557`, `3fc098a3`) and to `docs/geometric-precision-plan.md`.

## Problem

In a single-precision build (`-DUSE_SINGLE_PRECISION`, `WP=4`) FESOM computes every static
mesh metric by **differencing single-precision coordinates**:

- `src/oce_mesh.F90:2211-2218` — `a = coord_nod2D(:,n2)-coord_nod2D(:,n1)`; `elem_area = ½|a×b|`
- `src/oce_mesh.F90:2541`, `:2551-2576` — `edge_dxdy`, `edge_cross_dxdy`
- `src/oce_mesh.F90:2627-2645` — `gradient_sca` (= Δcoord/area), `gradient_vec`, `metric_factor`, `elem_cos`

Node lon/lat (degrees, O(180)) are read into a **`WP` buffer** (`src/oce_mesh.F90` `rbuff`,
broadcast `MPI_WP`; `x=rbuff*rad`), so each node position is already SP-rounded (~1 m) **before**
any differencing. The relative error in an edge vector / area / gradient is therefore
≈ `0.4 m / edge_length`:

| edge length | metric relative error (SP-differenced) |
|---|---|
| 10 km | 8e-5 |
| 5 km (ng5) | 1.5e-4 (area ≈ 3e-4) |
| **1 km** | **8e-4** |
| 500 m | 1.5e-3 |
| 100 m | 8e-3 |

This is the source of the benign `#### MASS MATRIX PROBLEM` warnings (SP area inconsistency),
and it grows sharply below ~1 km — so it must be addressed before finer meshes. It is
**secondary** to the SSH-CG-solve NaN already fixed on this branch (CORE2, with larger cells,
passes), but it is a cheap, permanent, future-proofing correctness fix.

**Goal:** compute the static horizontal metrics in **double precision at setup**, store them in
the existing `WP` arrays. DP build stays bit-identical; SP build gets metrics accurate to `WP`
(~1e-7) instead of ~1e-4, independent of resolution.

## Approach — compute-in-DP / store-WP (zero ripple)

Because `coord_nod2D` is already SP-rounded, DP local temporaries alone do **not** help — the
DP values must come from the read.

1. **Read coordinates in DP.** In the mesh-read loop (`src/oce_mesh.F90:~340-426`): make the read
   buffer `rbuff`, the `x/y` scalars and the `rad` factor `real64`; broadcast with
   `MPI_DOUBLE_PRECISION`. Fill a **transient DP coordinate array** `coord_nod2D_dp(2, myDim+eDim)`
   (owned+halo, exactly like `coord_nod2D` at `:421-425` — no halo exchange, same as today) **and**
   down-convert into the existing `WP` `coord_nod2D` for runtime.
   - `coord_nod2D_dp`: `real64` allocatable field on `t_mesh` (`src/MOD_MESH.F90`), allocated at
     read, **deallocated at the end of mesh-metric setup** (not needed at runtime).

2. **DP-compute the metrics.** In `mesh_areas` (`src/oce_mesh.F90:2162-2355`) and the edge/gradient
   setup (`:2451-2650`): make the local temporaries (`a`, `b`, `deltaX*/deltaY*`, `dfactor`,
   `ax`, `ay`, and the `r_earth`/`r_earth**2` scalings) `real64`, source coordinates from
   `coord_nod2D_dp`, do all differencing/cross-products/scalings in `real64`, then store the
   result into the existing `WP` arrays: `elem_area`, `area`, `areasvol(_inv)`, `edge_dxdy`,
   `edge_cross_dxdy`, `gradient_sca`, `gradient_vec`, `metric_factor`, `elem_cos`.

3. **Deallocate `coord_nod2D_dp`** once metric setup completes.

### What deliberately does NOT change
- **Metric arrays stay `real(kind=WP)`** → every existing halo exchange
  (`exchange_elem(elem_area)` `:2220`; `exchange_nod(area/areasvol)` `:2322-2323`;
  `exchange_elem(metric_factor/elem_cos)` `:2526-2527`; …) keeps resolving to the `WP` overload.
  **No new DP exchange overloads** (the `exchange_*` generics are `WP`/integer only).
- **`coord_nod2D`/`geo_coord_nod2D` stay `WP`** at runtime → **zero ripple** through the ~30
  runtime consumers (grid rotation `gen_modules_rotate_grid.F90`, forcing, mixing, tides, ice,
  tracers, dynamics, diagnostics, output). Promoting the coordinates themselves to `real64` was
  rejected precisely because it breaks every `real(kind=WP)` `r2g`/`g2r`/`vector_g2r` consumer.
- **DP build (`WP=8`)**: `real64==WP`, `coord_nod2D_dp==coord_nod2D`, `MPI_DOUBLE_PRECISION==MPI_WP`
  → **bit-identical**.

## Out of scope (future work)
- **ALE per-timestep vertical arrays** `hnode`/`helem`/`zbar_3d_n`/`Z_3d_n` are recomputed every
  step (`vert_vel_ale`, `update_thickness_ale` `src/oce_ale.F90:1143-1249`) — evolving `WP` *state*,
  not static geometry. The cumulative `zbar_3d_n` sum over ~70 levels is a possible future
  DP-accumulation target.
- **Rotated meshes** (`force_rotation`, `g2r` `gen_modules_rotate_grid.F90:89`, `WP` args): the
  rotation stays `WP`. Add a DP `g2r` path or accept residual WP-rotation error; non-rotated
  global meshes (e.g. ng5) are unaffected.

## Files
- `src/MOD_MESH.F90` — add `real64` allocatable `coord_nod2D_dp` to `t_mesh`.
- `src/oce_mesh.F90` — DP read (`~340-426`), `mesh_areas` (`2162-2355`), edge/gradient setup
  (`2451-2650`); allocate/deallocate `coord_nod2D_dp`.

## Verification
1. **DP regression (bit-identical):** default DP build, regenerate a pi/core2 mesh diag, diff
   `elem_area`/`gradient_sca` against a pre-change DP diag → identical.
2. **SP accuracy gain:** SP build; `elem_area` should now match the DP reference to ~1e-6 (was
   ~1e-4). Frugal check only — read 1-D `elem_area`; do **not** load 3-D `nod_area` fully (OOMs on ng5).
3. **Symptom check:** SP CORE2 run — `#### MASS MATRIX PROBLEM` warnings should largely disappear.
4. **ng5 (Levante):** rebuild Intel SP; confirm no regression and observe any additional benefit
   (independent of the CG fix already under test).
