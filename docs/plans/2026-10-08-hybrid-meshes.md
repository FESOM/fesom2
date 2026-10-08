# Hybrid (mixed triangle + quad) meshes

## Overview
- FESOM3 runs on genuinely mixed meshes: triangles and quadrilaterals coexisting in one mesh,
  stored with 4 vertices per element (a triangle repeats its 1st vertex as the 4th), clockwise,
  as in the `qq` prototype (`/home/a/a270029/qq`, example only — not copied).
- Setup stays as in FESOM2: run `fesom_ini` (ported into this repo, mixed-mesh capable) to
  partition the mesh and write `dist_N/` + edge files, then run the model.
- First test meshes: NEMO ORCA025 (all quads, `/work/ab0995/a270029/meshes/nemo025q`), CORE2 and
  pi (all triangles), and a mixed mesh made by splitting part of the NEMO quads into triangles.
- Scope: everything at once — setup tools, mesh geometry, all ocean kernels, sea ice, and two new
  exchanges (edges; double halo for nodes) for future developments.
- Triangle runs are NOT required to stay bit-identical (user decision): kernels use unified
  4-slot formulas as in `qq`; triangle regression is gated at roundoff level.

## Context (from discovery)
- **Mesh input today**: `src/mesh/mod_mesh_read.F90` (`read_mesh` 1-rank, `read_mesh_local`
  multi-rank) reads `nod2d.out`, 3-column `elem2d.out`, `aux3d.out`, `nlvls.out` and the
  precomputed `edgenum.out` / `edges.out` / `edge_tri.out`; partition + halos come from FESOM2's
  `fesom_ini` (`dist_N/rpart.out`, `my_list*.out`, `com_info*.out`) via
  `src/infra/mod_partitioning.F90`. FESOM3 has no partitioner.
- **t_mesh** (`src/types/mod_mesh.F90`): `MAX_NV=4`, `elem2D_nnodes`, `gradient_sca(2*MAX_NV,·)`
  (only rows 1:6 filled, x 1:3 / y 4:6), `elem_edges` never built, `elem_neighbors` allocated `(3,·)`.
- **Geometry** (`src/mesh/mod_mesh_areas.F90`): centroid `/3`, triangle area, node area
  `elem_area/3`, P1 `gradient_sca`, `build_elem_adjacency` (3 sides), edge cross vectors (generic).
- **Triangle assumptions in kernels** (~25 modules): ~26 `elem2D_nodes(1:3,·)` slices, ~41
  `gradient_sca(1:3|4:6)` uses, ~24 `do …=1,3` vertex loops, ~30 `/3` averages — ocean dynamics
  (`oce_dyn_velrhs`, `oce_ale`, `oce_pgf`, `oce_dyn_ivertvisc`, `oce_dyn_vinv`), tracers
  (`oce_tracer_grad`, `oce_ale_tracer`, `oce_adv_tra_fct`, `oce_adv_tra_flux`, `oce_muscl_adv`,
  `oce_pressure_bv`, `oce_initial_state`), GM (`oce_fer_gm`), mixing (KPP, TKE, PP, MO conv),
  forcing (`mod_forcing_bulk`), ice (`mod_ice_dyn`, `mod_ice_fct`, `mod_ice_setup`,
  `mod_ice_oce_coupling` — latent shape bug at :105), IO (`mod_io_coords`, `mod_io_meshdiag`,
  `mod_io_means`), `infra/mod_geom_dump`, drivers (`fesom_pressuredump`, `fesom_conserve`).
- **SSH matrix** (`src/oce/oce_ssh_rhs.F90` `init_stiff_mat_ale`/`update_stiff_mat_ale`): sparsity
  from edges only, fixed `n_pos(12,·)` — misses quad diagonals. Ice mass matrix
  (`mod_ice_setup.F90:177`) shares it, with no not-found guard in its column search.
- **FESOM2 fesom_ini** (`/home/a/a270029/fesom2.0/src`): `fvom_init.F90` (1978 lines; reads
  4-column elements, still triangle-only at 314/334/403/443/858/1068/1509/1591/1600 and in
  `test_tri_ini`, edge building, `find_levels`), `gen_comm.F90` (659; has q=3/4 branches),
  `oce_local.F90` (`save_dist_mesh`), `fort_part.c`, bundled METIS 5.1.0 (`lib/metis-5.1.0`).
- **NEMO mesh**: `nemo_meshq.m` (manual fixes: f-point corrections rows 1020/1021, single-cell
  land removal, strait edits, tripolar fold, cyclic mapping), outputs `nod2d.out` (915 956 nodes,
  4th column = land-boundary flag), `elem2d.out` (902 454 quads, SW-NW-NE-SE = clockwise),
  `depth.out` (47 positive levels, positive node depths = levels(count+1), per-node layer counts,
  min 3).
- **Patterns**: tests via `add_fesom_test(<name> <nranks>)` in `test/CMakeLists.txt`; build with
  `./configure.sh`; ctest; core2 runs need `ulimit -s unlimited`; plans move to
  `docs/plans/completed/`; lessons in `docs/LESSONS.md`; handoff in `docs/HANDOFF.md`.

## Development Approach
- **testing approach**: TDD where an oracle exists (Python port vs `.m` outputs, ported fesom_ini
  vs FESOM2 fesom_ini, geometry identities, exchange round trips, matrix identities); regular
  (code first, then tests) for kernel conversions, gated by the existing self-tests,
  conservation tests and roundoff-level regression.
- complete each task fully before moving to the next; small focused changes
- **CRITICAL: every task MUST include new/updated tests** for code changes in that task
- **CRITICAL: all tests must pass before starting next task** — no exceptions
- **CRITICAL: update this plan file when scope changes during implementation**
- build in a scratch build directory (`cmake -S . -B <scratch>`), never in `build_intel_dp`
  while the user's jobs run from it
- triangle meshes must keep working after every task (ctest green); quad/mixed meshes are
  exercised as self-tests on cut-outs from Task 7 on and become runnable after Task 16 (first
  runs in Task 17)
- triangle regression against today's model (binary of commit 71e1d49, `tools/run_tri_regression.sh`,
  created in Task 7) is re-run after Tasks 7, 10, 14 and 16, so structural errors are caught early
- no pull requests until the user says; commits without Claude attribution; `git add` explicit paths

## Testing Strategy
- **unit/self-tests** (Fortran programs registered with `add_fesom_test`, np 1/2[/8]):
  geometry identities, exchange round trips, SSH and ice matrix identities, kernel tests
  already in the suite (adapted from bitwise-vs-old to tolerances where the old bitwise
  reference was the triangle-specific arithmetic).
- **tool tests**: Python (`tools/mesh/`) byte comparisons and mesh validity checks, run from a
  small shell/pytest-free script so they work on Levante's python3.
- **integration**: conservation tests (`fesom_conserve`, heat/salt 1e-12) extended to the NEMO
  cut-out and mixed meshes; short NEMO / mixed runs (1 vs N ranks); CORE2 regression and
  throughput.
- no e2e/UI tests in this project.

## Progress Tracking
- mark completed items with `[x]` immediately when done
- add newly discovered tasks with ➕ prefix
- document issues/blockers with ⚠️ prefix
- update plan if implementation deviates from original scope

## Solution Overview
1. **Mesh generation tools (Python)**: faithful port of `nemo_meshq.m`; a quad-splitting tool
   producing the mixed test mesh.
2. **Setup (`fesom_ini` in FESOM3)**: FESOM2's partitioner/halo/edge writer ported and fixed for
   mixed meshes; same file formats plus appended sections (edge communication, double-halo
   nodes) behind a format version.
3. **Model mesh layer**: 4-slot connectivity, mixed-mesh geometry (centroid, area, kite weights,
   Green–Gauss `gradient_sca(8)`), conventions without zero indices (`edge_tri`,
   `elem_neighbors`). `gradient_vec` stays deferred (no kernel uses it yet); when it comes, it
   must respect the "missing neighbour = element itself" convention.
4. **Kernels**: unified 4-slot formulas everywhere (element mean = Σ `elem_w`·f, gradients over
   4 slots); explicit vertex count (`elem2D_nnodes`) only where a repeated vertex would matter.
5. **SSH matrix** with node-to-all-element-vertices sparsity (quad diagonals).
6. **Sea ice**: strain rates/stress via the ocean gradient arrays; FCT mass matrix with the
   square-bilinear quad pattern (S/36)(4,2,1,2), lumped S/4.
7. **Exchanges**: `exchange_edge` (unique edge owner) and `exchange_nod_full` (double halo).

## Technical Details
**File formats**
- `elem2d.out`: count, then 3 or 4 ints per line (column count detected on read; 3 columns →
  4th := 1st). Triangle detection: `elem2D_nodes(4,el) == elem2D_nodes(1,el)`.
- Depth file: FESOM2's namelist semantics are kept (`use_depthfile='aux3d'` → `aux3d.out`:
  nl, zbar(nl), node depth); a new value `use_depthfile='nemo'` reads `depth.out` (positive
  levels and depths, followed by per-node layer counts). Then `nlevels_nod2D = count + 1` is
  taken as given — no further level adjustment, exactly as `nemo_meshq.m` defines it.
  fesom_ini writes `nlvls.out`; the model takes node levels ONLY from `nlvls.out` (single
  source) and reads zbar / node depth from `aux3d.out`, or `depth.out` if `aux3d.out` is absent.
- `dist_N`: all existing files keep their exact FESOM2 layout (FESOM2 and `tools/make_dist1.py`
  must still read them); a separate marker file `dist_N/format.out` (version) announces the
  appended sections: `my_list*.out` + `eXDim_nod2D` and its global ids; `com_info*.out` +
  `com_nod2D_full` and `com_edge2D`. Without `format.out` the model reads as today and the new
  exchanges stop with a clear message if called. `make_dist1.py` and the built-in 1-rank
  synthesis (`mod_partitioning.F90` `synthesize_1rank`) provide empty new lists.
- Edge ownership: FESOM2's `myDim_edge2D` holds every edge with at least one owned node
  (`gen_comm.F90:629-639`), so shared edges are in the owned range of two ranks. `com_edge2D`
  sends from the UNIQUE owner (owner of `edges(1,edge)`) to every other rank that holds the
  edge, including ranks where it lies in 1..myDim_edge2D.
- `edge_tri.out` keeps FESOM2's encoding on disk; the model converts a missing 2nd element to
  the 1st one at read time.

**Mesh arrays (t_mesh)**
- `elem2D_nodes(4,·)`, `elem2D_nnodes(·)` (3/4), `elem_edges(4,·)` (1-2, 2-3, 3-4, 4-1;
  triangle 4th = 1st), `elem_neighbors(4,·)` (missing = the element itself), `edge_tri(2,·)`
  (missing 2nd = 1st), `edge_cross_dxdy(4,·)` (3:4 = 0 for a missing 2nd element),
  `gradient_sca(8,·)` (x 1:4, y 5:8; triangle slot 4 = 0), `elem_w(4,·)`.
- Boundary tests: every `el(2) > 0` / `el(2) <= 0` becomes `el(2) /= el(1)` (or is dropped where
  the zero cross vector already makes the term vanish); every `nb <= 0` neighbour test becomes
  `nb == el`. These changes land in the SAME task as the convention (Task 7), otherwise boundary
  edges would count twice.
- Node arrays keep size `nNodL = myDim+eDim`; `owned_bounds` additionally returns
  `nNodX = myDim+eDim+eXDim` (fields needing the double halo are allocated to it).

**Geometry formulas** (radians with `elem_cos`, then `r_earth`; cyclic wrap on all vertices)
- centroid: area-weighted centroids of triangles (1,3,4) and (1,2,3); area: their sum
  (triangle: the two coincide with the plain triangle formulas).
- kite weight at vertex j: ½·[A(j, j−1, c) + A(j, j+1, c)] / A_elem (= mid-edge kite);
  triangle (1/3, 1/3, 1/3, 0).
- node area: Σ over elements and 4 slots of `elem_area · elem_w(k)`.
- `gradient_sca`: Green–Gauss, coefficient of vertex j ∝ (x_{j+1} − x_{j−1}) / (2A) over the
  element's `nv` vertices; slot 4 = 0 for triangles.
- element levels: `nlevels(el)` = min over vertices, `ulevels(el)` = max over vertices.

**Known property to monitor**: the Green–Gauss quad gradient is blind to the checkerboard
pattern (+1,−1,+1,−1) on parallelogram quads, so on all-quad meshes the SSH operator and the
pressure gradient have a checkerboard null mode (CG still converges thanks to the area/dt
diagonal). It is tested and documented (Task 9) and monitored in the short runs (Task 17).

**Ice FCT element mass** — triangle: S/6 diagonal, S/12 off-diagonal; quad: (S/36)·(4,2,1,2)
rotated per row; lumped mass = row sums (S/3 resp. S/4 per vertex).

## What Goes Where
- **Implementation Steps** (`[ ]`): code, tests, tools, docs in this repo.
- **Post-Completion**: long simulations (user), comparison of NEMO/mixed climatologies, the PR.

## Implementation Steps

### Task 1: Python port of nemo_meshq.m

**Files:**
- Create: `tools/mesh/nemo_meshq.py`
- Create: `tools/mesh/test_nemo_meshq.sh`

- [ ] copy the three `.m` outputs (`nod2d.out`, `elem2d.out`, `depth.out`) to a reference
      directory (`/work/ab0995/a270029/meshes/nemo025q/ref_matlab/`) before anything writes there
- [ ] write `test_nemo_meshq.sh`: runs the port into a scratch dir and `cmp`s the three files
      against the reference (TDD: fails first)
- [ ] port Steps 1–6 statement by statement (netCDF reads, f-point corrections at rows
      1020/1021, single-cell land removal, every strait edit, node indices/depth, fold +
      cyclic node mapping, quad creation incl. the half fold row, unused-node removal, levels,
      ASCII output with the same `%9i`/`%8.4f` formats); keep MATLAB 1-based indices visible in
      comments next to each manual fix
- [ ] run the test — byte-identical output must pass before Task 2

### Task 2: Mixed test mesh by splitting quads

**Files:**
- Create: `tools/mesh/split_quads.py`
- Create: `tools/mesh/check_mesh.py`
- Create: `tools/mesh/test_split_quads.sh`

- [ ] `check_mesh.py`: validity checks for any 3/4-column mesh — every element clockwise
      (both triangles of a quad), no degenerate element, every node used, triangle 4th = 1st,
      total spherical area reported
- [ ] `split_quads.py`: split selected quads (a latitude band + every k-th quad elsewhere,
      parameters on the command line) along the shorter diagonal into two clockwise triangles
      written as 4 columns; nodes and depth file copied unchanged
- [ ] tests: `check_mesh` passes on the NEMO mesh and on the split mesh; total area preserved
      to roundoff; element count = quads + split count; a deliberately flipped element is
      reported (error case)
- [ ] produce the mixed mesh in `/work/ab0995/a270029/meshes/nemo025mix/`
- [ ] run tests — must pass before Task 3

### Task 3: Port fesom_ini into FESOM3 (triangle meshes unchanged)

**Files:**
- Create: `src/ini/` (from FESOM2 `fvom_init.F90`, `gen_comm.F90`, `oce_local.F90` save part,
  `fort_part.c`, and the FESOM2 modules they need: `MOD_MESH`, `MOD_PARTIT`, `g_config`,
  `o_PARAM`, `g_rotate_grid`)
- Create: `lib/metis-5.1.0/` (bundled as in FESOM2, same IDXTYPEWIDTH/REALTYPEWIDTH and flags)
- Modify: `CMakeLists.txt` (separate target `fesom_ini` with its OWN module directory, so FESOM2's
  `mod_mesh`/`mod_partit` never clash with FESOM3's)
- Create: `tools/run_ini_compare.sh`

- [ ] bring the sources over with FESOM2's namelist semantics (`namelist.config`: meshpath,
      n_part, `use_depthfile`, cavity options) and make them compile in this build
- [ ] write `run_ini_compare.sh`: run FESOM2's and FESOM3's `fesom_ini` on pi (2, 4 parts) and
      CORE2 (8, 128 parts) and `cmp` every output file (TDD); if METIS output differs only
      because of build differences, fall back to structural equality via `check_dist.py`
      (Task 4) and record why
- [ ] fix until identical on triangle meshes
- [ ] run the comparison — must pass before Task 4

### Task 4: fesom_ini mixed-mesh capable

**Files:**
- Modify: `src/ini/*` (orientation test, edge building, element neighbours, find_levels,
  depth-file reading, halo/communication over 4 vertices)
- Create: `tools/mesh/check_dist.py`

- [ ] orientation test for triangles and quads (as qq `test_elem`: both sub-triangles clockwise,
      swap 2↔4 if anticlockwise, stop on a self-intersecting quad)
- [ ] edge building and edge_tri over `nv` sides; all `elem2D_nodes(1:3,·)` slices and
      `DO q=1,3` loops → 4 slots / `nv` (fvom_init 314/334/403/443/858/1068/1509/1591/1600)
- [ ] depth: `use_depthfile='nemo'` reads `depth.out` (positive values); trailing layer counts
      taken as given, no further level adjustment; `nlvls.out` written from them
- [ ] halo/communication lists and METIS graph over all element vertices (quad diagonals
      adjacent)
- [ ] `check_dist.py`: every element owned once, halos symmetric between ranks, node/element/edge
      lists complete, `nlvls.out` equal to the depth-file counts
- [ ] tests: re-run `run_ini_compare.sh` (triangles still identical); `check_dist.py` on NEMO and
      the mixed mesh for 2, 8, 128 parts; a mesh with a flipped element is fixed/reported (error case)
- [ ] run tests — must pass before Task 5

### Task 5: fesom_ini writes edge communication and double-halo nodes

**Files:**
- Modify: `src/ini/*` (save_dist_mesh + new list builders)
- Modify: `tools/mesh/check_dist.py`, `tools/make_dist1.py`

- [ ] marker file `dist_N/format.out` (version); existing files keep their FESOM2 layout
- [ ] `eXDim_nod2D`: vertices of elements 1..myDim+eDim+eXDim not among 1..myDim+eDim nodes;
      appended to `my_list*.out`; `com_nod2D_full` appended to `com_info*.out`
- [ ] `com_edge2D`: unique owner = owner of `edges(1,edge)`, sends to every other rank holding the
      edge (also where it lies in that rank's 1..myDim_edge2D); appended to `com_info*.out`
- [ ] `make_dist1.py` writes the marker and empty new lists
- [ ] tests: `check_dist.py` verifies the new lists (send/receive symmetric, every halo edge and
      double-halo node has exactly one owner; 4-rank case where owned edge ranges overlap) on pi,
      CORE2, NEMO, mixed; FESOM2's model still reads a new `dist_N` (existing files unchanged)
- [ ] run tests — must pass before Task 6

### Task 6: Model reads mixed meshes and the extended dist_N

**Files:**
- Modify: `src/mesh/mod_mesh_read.F90` (`read_mesh`, `read_mesh_local`)
- Modify: `src/infra/mod_partitioning.F90` (incl. `synthesize_1rank`), `src/types/mod_partit.F90`
- Modify: `src/types/mod_mesh.F90`
- Modify: `test/test_mesh.F90`, `test/test_partit.F90`

- [ ] 3/4-column `elem2d.out` (column detection), `elem2D_nnodes` from 4th = 1st, element levels
      min/max over 4 slots; node levels only from `nlvls.out`
- [ ] zbar / node depth from `aux3d.out`, or `depth.out` if `aux3d.out` is absent (positive accepted)
- [ ] `edge_tri` keeps the on-disk encoding in this task (convention change is Task 7)
- [ ] `format.out` present → read `eXDim_nod2D` list, `com_nod2D_full`, `com_edge2D`;
      `synthesize_1rank` provides empty lists
- [ ] tests: read NEMO and mixed meshes at np 1/2 (counts, nnodes histogram, levels vs
      `nlvls.out`); old pi `dist_N` still read; a malformed element line stops with a message
      (error case)
- [ ] run ctest — must pass before Task 7

### Task 7: Mixed-mesh geometry, conventions and analytic mixed patch

**Files:**
- Modify: `src/types/mod_mesh.F90` (`elem_w`, `elem_edges`, `elem_neighbors(4)`, `gradient_sca(8)`;
  binary mesh dump read/write layout)
- Modify: `src/mesh/mod_mesh_areas.F90`, `src/mesh/mod_mesh_rotate.F90` (cyclic wrap),
  `src/mesh/mod_mesh_analytic.F90` (4-slot edges, new conventions, mixed-patch generator)
- Modify: every `gradient_sca(1:3|4:6)` call site → `(1:4|5:8)` (mechanical, all modules)
- Modify: every `edge_tri` boundary test and `elem_neighbors` missing test (same task!):
  `oce_ale.F90:216,359`, `oce_adv_tra_hor.F90:73,171`, `oce_adv_tra_flux.F90:72`,
  `oce_ale_tracer.F90:302`, `oce_ssh_rhs.F90:321`, `oce_adv_tra_driver.F90:118`,
  `oce_adv_tra_fct.F90:182,255`, `oce_dyn_vinv.F90:100,225,293`, `oce_dyn_velrhs.F90:245`,
  `mod_mesh_areas.F90:339`, `fesom_advhordump.F90:253` (and any further hit of the grep)
- Modify: `src/infra/mod_geom_dump.F90`, gate tools under `tools/` (see below)
- Create: `test/test_geom_mixed.F90`, `tools/run_tri_regression.sh`

- [ ] geometry: centroid, area, kite weights, node area, `gradient_sca` (Green–Gauss),
      `elem_edges`, `elem_neighbors` (missing = self), `edge_tri` missing 2nd = 1st with
      `edge_cross_dxdy(3:4) = 0`, as in Technical Details
- [ ] boundary tests converted in the same change (`el(2) /= el(1)`, `nb == el`), grep shows none left
- [ ] `mod_mesh_analytic`: edges over `nv` sides, new conventions, generator for a small mixed patch
- [ ] gate migration decision recorded in HANDOFF: FESOM2-oracle byte gates (`tools/run_*_gate*.sh`,
      `geom_diff.py`, `mod_geom_dump`) are retired or remapped to the new layout and listed
- [ ] `test_geom_mixed` (TDD): on pi, a NEMO-mixed cut-out and the analytic mixed patch —
      Σ node area = Σ element area; Σ `elem_w` = 1; `gradient_sca` exact for linear fields on
      triangles and quads; cross-edge vectors close around each element; triangle geometry equals
      the previous formulas to roundoff
- [ ] adapt existing tests whose bitwise reference was the old triangle arithmetic to tolerances
- [ ] `run_tri_regression.sh` (pi, CORE2: a few steps vs the 71e1d49 binary, roundoff level) — run
- [ ] run ctest — must pass before Task 8

### Task 8: Edge exchange and double-halo node exchange

**Files:**
- Modify: `src/infra/mod_halo.F90`, `src/infra/mod_part_bounds.F90` (`nNodX`)
- Create: `test/test_halo_ext.F90`

- [ ] `exchange_edge` (2D, 3D) over `com_edge2D`; `exchange_nod_full` (2D, 3D) over
      `com_nod2D_full`, same buffer/Isend/Irecv pattern as `exchange_nod`
- [ ] `owned_bounds` returns `nNodX`; coordinates/levels of double-halo nodes available
- [ ] tests: round trip (fill owned with global id, exchange, compare) on pi, CORE2, NEMO at np 2,
      4, 8, incl. edges owned on two ranks; calling without `format.out` stops with the message
      (error case)
- [ ] run ctest — must pass before Task 9

### Task 9: SSH stiffness matrix for mixed meshes

**Files:**
- Modify: `src/oce/oce_ssh_rhs.F90` (`init_stiff_mat_ale`, `update_stiff_mat_ale`)
- Create: `test/test_ssh_stiff.F90`

- [ ] sparsity = node plus all vertices of all elements around it (dynamic, no fixed 12)
- [ ] assembly over edges with 4 coefficients (`gradient_sca(1:4|5:8)`, `npos(4)`)
- [ ] tests: matrix symmetric; row sums = `area/dt` term; quad diagonal entries present and
      nonzero; triangle-mesh matrix equals the previous one to roundoff; CG converges on the NEMO
      cut-out; checkerboard (+1,−1,+1,−1) on a parallelogram quad patch is in the null space of the
      gradient part (documented property, asserted so it can't change silently)
- [ ] run ctest — must pass before Task 10

### Task 10: Ocean dynamics kernels

**Files:**
- Modify: `src/oce/oce_dyn_velrhs.F90`, `src/oce/oce_ale.F90`, `src/oce/oce_pgf.F90`,
  `src/oce/oce_dyn_ivertvisc.F90`, `src/oce/oce_dyn_vinv.F90`, `src/oce/oce_dyn_visc.F90`
- Modify: `test/test_vinv.F90`, `test/test_ivertvisc.F90`
- Create: `test/test_pgf_mixed.F90`

- [ ] element means via `elem_w`, vertex loops over 4 slots / `nv`
- [ ] `oce_pgf`: per-vertex arrays and loops (`idx`, `dx10`, `do ni`) to 4 slots
- [ ] `oce_dyn_vinv`: KE and vorticity node↔element with `elem_w`; `rvo_upwind` over
      `elem_neighbors(1:4)` and `nv` side normals
- [ ] tests: existing exact tests (uniform flow, linear shear) on triangle and mixed patches;
      `test_pgf_mixed`: resting stratified ocean gives zero pressure-gradient force on the mixed patch
- [ ] `run_tri_regression.sh` — roundoff level
- [ ] run ctest — must pass before Task 11

### Task 11: Tracer kernels and MUSCL up/downwind elements for quads

**Files:**
- Modify: `src/oce/oce_tracer_grad.F90`, `src/oce/oce_ale_tracer.F90`, `src/oce/oce_adv_tra_fct.F90`,
  `src/oce/oce_adv_tra_flux.F90`, `src/oce/oce_adv_tra_hor.F90`, `src/oce/oce_adv_tra_driver.F90`,
  `src/oce/oce_muscl_adv.F90`, `src/oce/oce_pressure_bv.F90`, `src/oce/oce_initial_state.F90`
- Modify: `test/test_muscl_onthefly.F90`
- Create: `test/test_updown_mixed.F90`

- [ ] element gradients/means over 4 slots; FCT element bounds over 4 slots
- [ ] `find_up_downwind_triangles` → elements: at vertex k the angle spanned by the two
      adjacent sides (next/previous vertex), triangles and quads (qq `oce_muscl_adv.F90:56-255`)
- [ ] `oce_pressure_bv` node smoothing with `elem_w` and node area
- [ ] tests: `test_updown_mixed` — up/downwind element found for every interior edge of the mixed
      patch and the edge direction lies in its vertex angle; MUSCL equals the reference rule;
      linear tracer advected exactly on uniform flow (interior)
- [ ] run ctest — must pass before Task 12

### Task 12: Mixing and GM/Redi

**Files:**
- Modify: `src/oce/oce_fer_gm.F90`, `src/oce/oce_mixing_kpp.F90`, `src/oce/oce_mixing_tke.F90`,
  `src/oce/oce_ale_mixing_pp.F90`, `src/oce/oce_mo_conv.F90`

- [ ] element means via `elem_w`
- [ ] tests: KPP/TKE/PP/GM conservation variants (`fesom_conserve`) on the mixed cut-out
- [ ] run ctest — must pass before Task 13

### Task 13: Forcing and ice–ocean coupling

**Files:**
- Modify: `src/forcing/mod_forcing_bulk.F90`, `src/ice/mod_ice_oce_coupling.F90`

- [ ] element means via `elem_w`; fix the `elnodes(3)` shape bug in the coupling (:105)
- [ ] tests: forcing/coupling self-tests on CORE2 (roundoff) and the mixed cut-out
- [ ] run ctest — must pass before Task 14

### Task 14: IO, mesh diagnostics and dump drivers

**Files:**
- Modify: `src/io/mod_io_coords.F90`, `src/io/mod_io_meshdiag.F90`, `src/io/mod_io_means.F90`,
  `src/drivers/fesom_pressuredump.F90`, `src/drivers/fesom_conserve.F90`,
  `src/drivers/fesom_advhordump.F90`, `fesom_geomdump`, `fesom_meshdiagdump`, `fesom_evpdump`,
  `fesom_icefctdump` drivers

- [ ] element coordinates = centroids; mesh diagnostic output with 4-column connectivity
      (dimension `n4`), `face_edges`/`face_links` from the now-built `elem_edges`/`elem_neighbors`
- [ ] dump drivers follow the gate decision of Task 7 (retired or remapped)
- [ ] grep check: no `elem2D_nodes(1:3`, `gradient_sca(1:3`, `(4:6`, element `/3` average left in
      `src/` (documented exceptions only, e.g. the spline coefficient `oce_adv_tra_ver.F90:134`)
- [ ] tests: meshdiag file has 4 columns and triangle 4th = 1st; `run_tri_regression.sh` — roundoff
- [ ] run ctest — must pass before Task 15

### Task 15: Sea-ice dynamics (EVP / mEVP)

**Files:**
- Modify: `src/ice/mod_ice_dyn.F90`
- Create: `test/test_ice_dyn_mixed.F90`

- [ ] strain rates over 4 slots of `gradient_sca`; metric term with `elem_w` mean
- [ ] `stress2rhs`: node shares `elem_area·elem_w(k)` and `gradient_sca(k|k+4)`
- [ ] tests: `test_ice_dyn_mixed` — zero strain rate for uniform drift and rigid rotation on
      triangle and mixed patches; CORE2 ice step at roundoff vs before
- [ ] run ctest — must pass before Task 16

### Task 16: Sea-ice FCT and mass matrix

**Files:**
- Modify: `src/ice/mod_ice_setup.F90` (`ice_mass_matrix_fill`, loops at :187/:192),
  `src/ice/mod_ice_fct.F90`
- Create: `test/test_ice_mass.F90`

- [ ] `ice_mass_matrix_fill`: triangle S/6, S/12; quad (S/36)(4,2,1,2) rotated per row; on the
      new SSH pattern; not-found column stops with an error
- [ ] per-element coefficient table `icoef(4,4)` from the element mass pattern for
      `ice_TG_rhs` and `ice_fem_fct`; lumped mass = row sums (S/3, S/4); limiter loops over `nv`
- [ ] ⚠️ if any Taylor–Galerkin term is not a plain mass-matrix entry, stop and ask the user
      before choosing its quad form
- [ ] tests: row sums = lumped mass, total = total area; constant ice field stays constant under
      advection; mass conserved on the mixed patch; missing column triggers the error (error case)
- [ ] `run_tri_regression.sh` — roundoff
- [ ] run ctest — must pass before Task 17

### Task 17: Conservation and short runs on NEMO and mixed meshes

**Files:**
- Modify: `test/CMakeLists.txt` (NEMO cut-out / mixed variants of `fesom_conserve`)
- Create: `tools/run_nemo_short.sh`

- [ ] conservation tests (heat/salt 1e-12) on the NEMO cut-out and the mixed cut-out, np 1/2
- [ ] short NEMO and mixed runs (a few days) at 1 vs N ranks: stable, conserving, N vs 1 within
      roundoff; with ice on; checkerboard energy of SSH and of a tracer monitored
- [ ] run ctest + short runs — must pass before Task 18

### Task 18: CORE2 throughput and example job for the NEMO mesh

**Files:**
- Create: `work/job_levante_nemo025` (new file; `work/job_levante` untouched)
- Modify: `docs/HANDOFF.md`

- [ ] CORE2 throughput with 4-slot loops vs ~14 min per model year, reported
- [ ] document the workflow: `nemo_meshq.py` → (optional `split_quads.py`) → `fesom_ini` → model
- [ ] example job script for NEMO 0.25° next to `job_levante`; dry-run its setup part (no submission)

### Task 19: Verify acceptance criteria
- [ ] mixed meshes read, partitioned and run (NEMO, mixed, CORE2, pi)
- [ ] all decisions of the Overview/Technical Details implemented; grep checks clean
- [ ] full suite: `./configure.sh` in a scratch build + `ctest` (with `ulimit -s unlimited`)
- [ ] tool tests (`tools/mesh/test_*.sh`), `run_ini_compare.sh`, `run_tri_regression.sh` pass

### Task 20: [Final] Update documentation
- [ ] `docs/HANDOFF.md`: mixed-mesh section (formats, arrays, conventions, setup, exchanges,
      retired gates, checkerboard property)
- [ ] `docs/LESSONS.md`: lessons learned
- [ ] move this plan to `docs/plans/completed/`

## Post-Completion
*Items requiring manual intervention or external systems - no checkboxes, informational only*

**Manual verification**
- long NEMO 0.25° and mixed-mesh simulations (run by the user), climate-level comparison with
  CORE2 and with NEMO reference fields
- visual check of the mixed mesh and of the tripolar fold region

**External**
- pull request only when the user says so, after testing
