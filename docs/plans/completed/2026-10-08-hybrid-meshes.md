# Hybrid (mixed triangle + quad) meshes

## Overview
- FESOM3 runs on genuinely mixed meshes: triangles and quadrilaterals coexisting in one mesh,
  stored with 4 vertices per element (a triangle repeats its 1st vertex as the 4th), clockwise,
  as in the `qq` prototype (`/home/a/a270029/qq`, example only — not copied).
- Setup stays as in FESOM2: run `fesom_ini` (FESOM2's `mesh_part` subproject, ported into this
  repo and made mixed-mesh capable) to partition the mesh and write `dist_N/` + edge files, then
  run the model.
- First test meshes: NEMO ORCA025 (all quads, `/work/ab0995/a270029/meshes/nemo025q`), CORE2 and
  pi (all triangles), and a mixed mesh made by splitting part of the NEMO quads into triangles.
- Scope: everything at once — setup tools, mesh geometry, all ocean kernels, sea ice, and two new
  exchanges (edges; double halo for nodes) for future developments.
- Triangle runs are NOT required to stay bit-identical once the geometry/kernels change (user
  decision): kernels use unified 4-slot formulas as in `qq`; triangle regression is gated at
  roundoff level from Task 7 on (bitwise before that).

## Context (from discovery, verified 2026-10-08)
- **Mesh input today**: `src/mesh/mod_mesh_read.F90` (`read_mesh` 1-rank, `read_mesh_local`
  multi-rank) reads `nod2d.out`, 3-column `elem2d.out`, `aux3d.out` (only `zbar` matters — node
  `depth` is stored but never used by the model), `nlvls.out` (node levels) and the precomputed
  `edgenum.out` / `edges.out` / `edge_tri.out`; partition + halos come from FESOM2's `fesom_ini`
  (`dist_N/rpart.out`, `my_list*.out`, `com_info*.out`) via `src/infra/mod_partitioning.F90`.
  FESOM3 has no partitioner.
- **Mesh directories are read-only for us**: pi (`/home/a/a270088/port2/fesom2/tests/data/MESHES/pi`)
  and CORE2 (`/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2`) belong to another account, and
  `fesom_ini` writes its outputs INTO the mesh directory → all setup runs use writable copies
  under `/work/ab0995/a270029/meshes/`.
- **t_mesh** (`src/types/mod_mesh.F90`): `MAX_NV=4`, `elem2D_nnodes`, `gradient_sca(2*MAX_NV,·)`
  (only rows 1:6 filled, x 1:3 / y 4:6; allocated for OWNED elements only), `elem_edges` never
  built, `elem_neighbors` allocated `(3,·)`.
- **Geometry** (`src/mesh/mod_mesh_areas.F90`): centroid `/3`, triangle area, node area
  `elem_area/3`, P1 `gradient_sca`, `build_elem_adjacency` (3 sides), edge cross vectors
  (generic, but `:385,:394` decide the boundary case on `el2 > 0`).
- **Triangle assumptions in kernels** (~25 modules, swept with grep): ~26 `elem2D_nodes(1:3,·)`
  slices, ~41 `gradient_sca(1:3|4:6)` uses, ~24 `do …=1,3` vertex loops, ~30 `/3` element
  averages, AND a second kind that the element-mean rule does not cover — **node gathers with a
  hidden 1/3**, `Σ_e A_e·x_e / (3·area)`: `oce_ale.F90:96-109` (`UVnode`),
  `oce_dyn_velrhs.F90:209-231` (`w·du/dz`, `/(3·hnode)`), `oce_ale_tracer.F90:456-468`
  (`tr_xynodes`), `oce_pressure_bv.F90` `compute_sigma_xy` (:534-578) and `smooth_nod`
  (:314-393), `mod_mesh_areas.F90:445,:512-521`, `mod_ice_dyn.F90:100-110`.
- **Boundary tests on `edge_tri(2,·)` / `elem_neighbors`** (all must change together with the
  convention): `oce_ale.F90:216,359`, `oce_adv_tra_hor.F90:73,171`, `oce_adv_tra_flux.F90:72`,
  `oce_ale_tracer.F90:302`, `oce_ssh_rhs.F90:165,239` (`el(i) < 1`) and `:321`,
  `oce_adv_tra_driver.F90:118`,
  `oce_adv_tra_fct.F90:182,255`, `oce_dyn_vinv.F90:100,225,294`, `oce_dyn_velrhs.F90:253`,
  `mod_mesh_areas.F90:339,385,394`, `mod_mesh_analytic.F90:144`, `mod_io_meshdiag` face links,
  `fesom_advhordump.F90:253`; tests `test_vinv.F90:104,123,151,197,519`, `test_bottom.F90:295`,
  `test_muscl_onthefly.F90:49,251`. Viscosity (`oce_dyn_visc.F90:127-165`) and the ice boundary
  conditions (`mod_ice_setup.F90:67-71`, `mod_ice_dyn.F90:415,650`) instead rely on the global
  edge ORDER (interior edges first, `edge2D_in`).
- **SSH matrix** (`src/oce/oce_ssh_rhs.F90` `init_stiff_mat_ale`/`update_stiff_mat_ale`): sparsity
  from edges only, fixed `n_pos(12,·)` — misses quad diagonals. The ice FCT (mass matrix and its
  neighbour lists) uses the same CSR rows, so the new pattern carries over; the column search in
  `ice_mass_matrix_fill` (`mod_ice_setup.F90:177`) has no not-found guard. With Green–Gauss
  coefficients the matrix stays symmetric on quads (dual-face normals ∝ gradient coefficients).
- **Ice FCT** (`mod_ice_fct.F90`, `mod_ice_setup.F90`): divides by the node area `mesh%area` as
  its lumped mass in the low-order step (:238-262), the high-order iteration (:184-190) and the
  antidiffusive fluxes (:343-355); `fct_fluxes(elem,3)` (`mod_ice_setup.F90:137`); hard-coded
  P1 factors (`/12`, `icoef(3,3)` = 1 off-diagonal, −2 diagonal). TG rhs (:96-113) has three
  terms: advection `S·dx(n)·(um+u_q)/12` (= gradient × P1 mass moment), diffusion
  `S·diff·(dx_n dx_q + dy_n dy_q)`, and the TG term with the element-mean velocity `um/3`.
- **FESOM2 setup tool** (`/home/a/a270029/fesom2.0`): a standalone CMake subproject `mesh_part/`
  builds `fesom_ini` from 11 FESOM2 sources (`MOD_MESH`, `oce_modules`, `gen_modules_config`,
  `gen_modules_partitioning`, `gen_modules_rotate_grid`, `fvom_init`, `oce_local`, `gen_comm`,
  `MOD_READ/WRITE_BINARY_ARRAYS`, `MOD_PARTIT`) + `fort_part.c` + METIS 5.1.0 (`lib/`, 1.7 MB);
  it reads `namelist.config` (paths, geometry, run_config, machine). ALREADY quad-aware:
  `stiff_mat_ini` (METIS graph incl. quad diagonals, :1609-1660), the node-neighbour lists of
  `find_edges_ini` (:390-470), `gen_comm` (`merge(3,4,…)`). Still triangle-only: `test_tri_ini`
  (:300-336), `find_levels`/`find_levels_cavity` element loops (:858, :1068), `elem_center`
  (:1500-1519), `find_elem_neighbors_ini` (:1520-1600, incl. its "≥ 2 neighbours" stop and the
  `nod_in_elem2D` build).
- **NEMO mesh** (`nemo_meshq.m`; manual fixes: f-point corrections rows 1020/1021, single-cell
  land removal, strait edits, tripolar fold, cyclic mapping): `nod2d.out` (915 956 nodes, 4th
  column = land-boundary flag), `elem2d.out` (902 454 quads, SW-NW-NE-SE), `depth.out`
  (47 positive levels, positive node depths = levels(count+1), per-node layer counts, min 3).
  Checked: all quads convex and clockwise on the sphere, nearly rectangular (smallest
  corner-turn ratio 0.92); every edge has ≤ 2 elements (27 624 boundary edges); **40 quads have
  only one neighbour** (dead ends) — FESOM's setup stops on them. Inputs are float32
  (`nav_lon/lat/lev`, `glamf/gphif`), int16 (`mbathy`), int8 (`tmaskutil`).
- **Python**: the system `python3` has no numpy; tools run under `module load python3`
  (numpy 1.23.4, netCDF4 1.6.0).
- **Patterns**: tests via `add_fesom_test(<name> <nranks>)` in `test/CMakeLists.txt` (mesh via
  `FESOM3_MESH_DIR`, skipped when a mesh is absent); build with `./configure.sh`; ctest; core2
  runs need `ulimit -s unlimited`; plans move to `docs/plans/completed/`; lessons in
  `docs/LESSONS.md`; handoff in `docs/HANDOFF.md`.

## Development Approach
- **testing approach**: TDD where an oracle exists (Python port vs `.m` outputs, ported fesom_ini
  vs existing `dist_N`, geometry identities, exchange round trips, matrix identities); regular
  (code first, then tests) for kernel conversions, gated by the existing self-tests,
  conservation tests and roundoff-level regression.
- complete each task fully before moving to the next; small focused changes
- **CRITICAL: every task MUST include new/updated tests** for code changes in that task
- **CRITICAL: all tests must pass before starting next task** — no exceptions
- **CRITICAL: update this plan file when scope changes during implementation**
- build in a scratch build directory (`cmake -S . -B <scratch>`), never in `build_intel_dp`
  while the user's jobs run from it
- triangle meshes keep working after every task (ctest green). Through Task 6 triangle runs stay
  BIT-IDENTICAL (reading/setup only); from Task 7 on, `tools/run_tri_regression.sh` compares pi
  and CORE2 against a frozen binary of the commit before Task 7 (roundoff level) after Tasks 7,
  10, 14 and 16
- quad/mixed meshes are exercised as self-tests on the small cut-out meshes from Task 6 on and
  are fully runnable after Task 16 (first runs in Task 17)
- production note: rebuilding `build_intel_dp` after Task 7 changes the trajectory of a running
  triangle experiment (e.g. results17) at roundoff level — keep a copy of the current binary if
  such a run must continue unchanged
- no pull requests until the user says; commits without Claude attribution; `git add` explicit paths

## Testing Strategy
- **unit/self-tests** (Fortran programs registered with `add_fesom_test`, np 1/2[/8]):
  geometry identities, exchange round trips, SSH and ice matrix identities, kernel tests
  already in the suite (adapted from bitwise-vs-old to tolerances where the old bitwise
  reference was the triangle-specific arithmetic).
- **test meshes** (created in Task 2/5, all under `/work/ab0995/a270029/meshes/`): writable copies
  of pi and CORE2 with regenerated `dist_N`; small cut-outs of the cleaned NEMO mesh and of the
  mixed mesh (a few thousand nodes) with `dist_1/2/8` — fast enough for ctest on a login node.
- **tool tests**: Python (`tools/mesh/`) byte comparisons and mesh validity checks, run from a
  small shell script under `module load python3`.
- **integration**: conservation tests (`fesom_conserve`, heat/salt 1e-12) on the cut-outs;
  1-vs-N rank runs on the cut-outs; CORE2 regression and throughput. Full-size NEMO/mixed runs
  are the user's (job script prepared).
- no e2e/UI tests in this project.

## Progress Tracking
- mark completed items with `[x]` immediately when done
- add newly discovered tasks with ➕ prefix
- document issues/blockers with ⚠️ prefix
- update plan if implementation deviates from original scope

## Solution Overview
1. **Mesh tools (Python)**: faithful port of `nemo_meshq.m`; dead-end cleanup (separate tool, so
   the port stays faithful); quad splitting for the mixed mesh; region cut-outs for tests.
2. **Setup (`fesom_ini` in FESOM3)**: FESOM2's `mesh_part` subproject ported as is, then made
   mixed-mesh capable; same file formats plus appended sections (edge communication,
   double-halo nodes) announced by a separate marker file.
3. **Model mesh layer**: 4-slot connectivity, mixed-mesh geometry (centroid, area, kite weights,
   Green–Gauss `gradient_sca(8)`), per-(node, element) kite weights for node gathers,
   conventions without zero indices (`edge_tri`, `elem_neighbors`). `gradient_vec` stays
   deferred (no kernel uses it yet); when it comes, it must respect "missing neighbour =
   element itself".
4. **Kernels**: unified 4-slot formulas everywhere (element mean = Σ `elem_w`·f, node gather =
   Σ `nod_in_elem2D_w`·A·x / area, gradients over 4 slots); explicit vertex count
   (`elem2D_nnodes`) where a repeated vertex would matter (plain counts, unique-vertex loops).
5. **SSH matrix** with node-to-all-element-vertices sparsity (quad diagonals).
6. **Sea ice**: strain rates/stress via the ocean gradient arrays and kite weights; FCT with the
   square-element mass contributions (S/36)(4,2,1,2) and its own lumped mass (S/4 per quad
   vertex) in all FCT steps (user's specification).
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
  source) and reads zbar from `aux3d.out`, or from `depth.out` if `aux3d.out` is absent.
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
- Edge order: interior edges first, boundary edges last (`edge2D_in` in `edgenum.out`) — kept
  for mixed meshes (viscosity and ice boundary conditions depend on it).
- `edge_tri.out` keeps FESOM2's encoding on disk; the model converts a missing 2nd element to
  the 1st one at read time.
- Output (mesh diagnostics / UGRID `face_nodes`): 4 columns, a triangle's 4th entry is the
  UGRID `_FillValue` −1 (user decision); internally the model repeats the 1st vertex.

**Mesh arrays (t_mesh)**
- `elem2D_nodes(4,·)`, `elem2D_nnodes(·)` (3/4), `elem_edges(4,·)` (1-2, 2-3, 3-4, 4-1;
  triangle 4th = 1st), `elem_neighbors(4,·)` (missing = the element itself), `edge_tri(2,·)`
  (missing 2nd = 1st), `edge_cross_dxdy(4,·)` (3:4 = 0 for a missing 2nd element),
  `gradient_sca(8,·)` (x 1:4, y 5:8; triangle slot 4 = 0), `elem_w(4,·)` (kite weights),
  `nod_in_elem2D_w(MAX_ADJACENT,·)` (kite weight of node n in its k-th element — for node
  gathers; owned nodes only need owned elements).
- Boundary tests: every `el(2) > 0` / `el(2) <= 0` / `el(i) < 1` becomes `el(2) /= el(1)` (or is
  dropped where the zero cross vector already makes the term vanish); every `nb <= 0` neighbour
  test becomes `nb == el`. These changes land in the SAME task as the convention (Task 7),
  otherwise boundary edges would count twice.
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
- `rvo_upwind` inward side normals: orientation reference = element centroid (not "the third
  vertex", which is ambiguous for quads).

**Known property — checkerboard mode**: the Green–Gauss quad gradient is blind to the pattern
(+1,−1,+1,−1) on parallelogram quads; ORCA025 quads are nearly rectangular, so on the all-quad
mesh the SSH operator and the pressure gradient have an (almost exact) checkerboard null mode.
The velocity divergence cannot drive it (it is the adjoint of the gradient), but freshwater
fluxes and nonlinear terms can. Asserted in a test (Task 9), measured in the short runs
(Task 17); ⚠️ if it grows, stop and ask the user (remedy is a design decision).

**Ice FCT** (user's specification)
- consistent mass matrix, element contributions — triangle: S/6 diagonal, S/12 off-diagonal;
  quad (regular-square simplification): row n1 (S/36)·(4, 2, 1, 2) for (n1, n2, n3, n4), row n2
  (2, 4, 2, 1), and so on.
- lumped mass = row sums: S/3 per triangle vertex, S/4 per quad vertex; it REPLACES `mesh%area`
  in all FCT steps (low-order step, high-order iteration, antidiffusive fluxes). The deviation
  from the scalar node area used in the ice dynamics is accepted (user).
- antidiffusive coefficients = consistent − lumped per element: triangle (−2, 1, 1)·S/12, quad
  (−5, 2, 1, 2)·S/36 (rotated per row); `fct_fluxes(elem, 4)`.
- TG rhs: advection term `S·dx(n)·(um+u_q)/12` → `dx(n)·(M_e u)_q` with the element mass
  contributions above; diffusion term unchanged in form (4-slot gradients); TG second-order
  term with the kite-weighted element-mean velocity instead of `um/3`. For triangles all
  three equal today's formulas up to roundoff.
- ice mass is then conserved in the lumped-mass sense (Σ lumped·m), which is what the tests check.

## What Goes Where
- **Implementation Steps** (`[ ]`): code, tests, tools, docs in this repo.
- **Post-Completion**: full-size and long simulations (user), climatology comparisons, the PR.

## Implementation Steps

### Task 1: Python port of nemo_meshq.m

**Files:**
- Create: `tools/mesh/nemo_meshq.py`
- Create: `tools/mesh/test_nemo_meshq.sh`

- [x] copy the three `.m` outputs (`nod2d.out`, `elem2d.out`, `depth.out`) to a reference
      directory (`/work/ab0995/a270029/meshes/nemo025q/ref_matlab/`) before anything writes there
- [x] write `test_nemo_meshq.sh` (runs under `module load python3`): runs the port into a scratch
      dir and `cmp`s the three files against the reference (TDD: fails first)
- [x] port Steps 1–6 statement by statement, emulating MATLAB semantics: float32 inputs stay
      float32 and every single/double mix is computed in float32 (MATLAB "single wins") — the
      `levels` recursion and the f-point corrections in particular; `mbathy` stays integer;
      1-based → 0-based indices with the MATLAB index kept in a comment at each manual fix
- [x] keep the order-dependent loops sequential (single-cell land removal modifies `bathy` while
      scanning; quad creation and node compaction define the numbering)
- [x] ASCII output with the same formats (`'%9i \n'` headers incl. the trailing space,
      `'%9i %8.4f %8.4f %8i\n'`, `'%9i %9i %9i %9i\n'`, `'%8.4f\n'`, `'%9i \n'`)
- [x] acceptance: `elem2d.out` and `depth.out` byte-identical; `nod2d.out` byte-identical except
      possibly the ~20 fold points corrected with single-precision cos/atan2/acos (MATLAB vs
      numpy may differ by one unit in the last place there) — any such difference listed and
      ≤ 1e-4 — result: elem2d.out and depth.out identical; nod2d.out differs in ONE line
      (node 711389, fold point of MATLAB row 1020/1021: port -83.6347, MATLAB -83.6348, i.e.
      1 float32 ulp of the atan2 result straddling the print rounding; correctly-rounded
      cos/atan2/acos variants were tried and match worse)
- [x] run the test — must pass before Task 2

### Task 2: Mesh tools — dead-end cleanup, quad splitting, cut-outs

**Files:**
- Create: `tools/mesh/check_mesh.py`, `tools/mesh/remove_deadends.py`,
  `tools/mesh/split_quads.py`, `tools/mesh/cut_region.py`
- Create: `tools/mesh/test_mesh_tools.sh`
- ➕ Create: `tools/mesh/mesh_io.py` (shared read/write, edges, neighbours, dead-end removal)

- [x] `check_mesh.py`: validity checks for any 3/4-column mesh — every element clockwise and
      convex on the sphere, no degenerate element, every node used, triangle 4th = 1st, every
      edge in ≤ 2 elements, count of elements with < 2 neighbours, total spherical area
- [x] `remove_deadends.py`: removes elements with fewer than two neighbours, iterating until none
      are left, drops orphan nodes, renumbers, keeps the depth file consistent (FESOM2's rule;
      `nemo_meshq.py` stays faithful) → cleaned NEMO mesh `/work/ab0995/a270029/meshes/nemo025q_fesom/`
      — result: 6 passes remove 46 elements (40 + 6 new dead ends), 90 nodes: 915 866 nodes,
      902 408 quads; nodes of removed elements get boundary flag 1 (the cell is land now)
- [x] `split_quads.py`: split selected quads (a latitude band + every k-th quad elsewhere,
      parameters on the command line) along the shorter diagonal into two clockwise triangles
      written as 4 columns; nodes and depth file unchanged → `/work/ab0995/a270029/meshes/nemo025mix/`
      — result (`--band -5 5 --every 5`): 216 956 quads split → 685 452 quads + 433 912 triangles;
      ➕ at the coast the shorter diagonal can leave a triangle whose only neighbour is its twin
      (838 cases = dead ends) → such quads are cut along the other diagonal (a quad is left whole
      only if both diagonals make a dead end: 0 cases)
- [x] `cut_region.py`: cut a lon/lat box out of a mesh (elements fully inside), renumber, recompute
      the boundary flag, subset the depth file → small test meshes `nemo025q_cut/`, `nemo025mix_cut/`
      (a few thousand nodes, with land and open-ocean boundaries closed as coast) — result: box
      lon −20…−5, lat −2…12 (Gulf of Guinea / West African coast): 2 509 nodes; `nemo025q_cut`
      2 391 quads, `nemo025mix_cut` 609 quads + 3 564 triangles; the cut keeps the largest
      connected piece and removes dead ends; flag = node on a boundary edge
- [x] tests: `check_mesh` passes on the cleaned, mixed and cut meshes (0 dead ends); the original
      NEMO mesh reports exactly the 40 dead ends; split preserves total area to roundoff and
      element count = quads + split count; a deliberately flipped element is reported (error case)
      — `tools/mesh/test_mesh_tools.sh [nemo_dir] [out_root]` (shared helpers in
      `tools/mesh/mesh_io.py`); also a bow-tie quad and a repeated vertex are reported
- [x] run tests — must pass before Task 3 — PASS (~35 s; area rel. change of the split 6.9e-16);
      the persistent test meshes were made with out_root `/work/ab0995/a270029/meshes`

### Task 3: Port FESOM2's mesh_part (fesom_ini) into FESOM3, triangle meshes unchanged

**Files:**
- Create: `mesh_part/` (FESOM2's subproject: its `CMakeLists.txt` and the 11 Fortran sources +
  `fort_part.c` copied under `mesh_part/src/`, with their own module directory so FESOM2's
  `mod_mesh`/`mod_partit` never clash with FESOM3's)
- Create: `lib/metis-5.1.0/` (as in FESOM2, same compile definitions incl. `METISRANDOMSEED`)
- Modify: `CMakeLists.txt` (`add_subdirectory(mesh_part)`, target `fesom_ini`), `configure.sh`
- Create: `tools/run_ini_compare.sh`

- [x] copy the subproject; keep FESOM2's `namelist.config` semantics (paths, geometry,
      run_config, machine: n_part); build it in the scratch build — sources from
      `/home/a/a270029/fesom2.0` (incl. its 5 `associate_*.h` headers) unchanged; private module
      dir `<build>/mesh_part/module`; METIS added `EXCLUDE_FROM_ALL` (only libmetis is built);
      FESOM2's mesh_part Fortran flags and C definitions; `configure.sh` passes `CC`/`CXX`;
      binary `<build>/bin/fesom_ini`
- [x] `run_ini_compare.sh` (TDD): make writable copies of pi and CORE2 (`nod2d.out`,
      `elem2d.out`, `aux3d.out` + a `namelist.config`) under `/work/ab0995/a270029/meshes/`, run
      the ported `fesom_ini` for pi (2, 8 parts) and CORE2 (8, 128 parts) and `cmp` every output
      (`edges.out`, `edge_tri.out`, `edgenum.out`, `nlvls.out`, `elvls.out`, `dist_N/*`) against the
      existing files the model uses today; second reference: FESOM2's own `mesh_part` built in
      scratch; if METIS output differs only because of the build, fall back to structural
      equality via `check_dist.py` (Task 4) and record why — writable copies
      `/work/ab0995/a270029/meshes/pi_w/` and `core2_w/` (+ CORE2 2 parts); FESOM2's build in
      `/work/ab0995/a270029/fesom3_hybrid/f2_meshpart_build/`, its runs on separate copies under
      `/work/ab0995/a270029/fesom3_hybrid/ini_compare/`
- [x] fix until identical on triangle meshes — nothing to fix: EVERY output (edges, levels,
      `dist_2/8/128`) is byte-identical to FESOM2's own build; against the existing mesh files:
      `edges.out`, `edge_tri.out`, `edgenum.out` (pi, CORE2), CORE2 `elvls/nlvls.out.20260528_regenerated`
      + `elvls_raw.out`, CORE2 `dist_2` and `dist_128` identical.
      ⚠️ Not reproducible (made by OLDER fesom_ini versions, not by the build): pi `dist_2/8` and
      `elvls/nlvls.out` (3 element levels differ), CORE2 `dist_8` (2023) and its current
      `elvls/nlvls.out` (2021 originals, hash = the mesh git LFS). These partitions are almost
      perfectly node-balanced (CORE2 dist_8: 15857/15858 per part), which suggests METIS was called
      without the node-level weights; a scratch probe with the reference levels, and one with unit weights, as
      weights did not reproduce them either (unit weights: pi dist_2 part sizes match, labels
      swapped). Reported as `info` by the script, structure to be checked by `check_dist.py` (Task 4)
- [x] run the comparison — must pass before Task 4 — `tools/run_ini_compare.sh` PASS (~30 s);
      ctest 48/48 in the scratch build

### Task 4: fesom_ini mixed-mesh capable

**Files:**
- Modify: `mesh_part/src/fvom_init.F90` (and `gen_comm.F90` / `oce_local.F90` only if an audit
  finds a gap)
- Create: `tools/mesh/check_dist.py`
- ➕ Create: `tools/run_ini_mixed.sh` (fesom_ini + check_dist on the NEMO meshes, error cases)

- [x] `test_tri_ini` → orientation test for triangles and quads (as qq `test_elem`: both
      sub-triangles clockwise, swap 2↔4 if anticlockwise, stop on a self-intersecting quad)
      — sub-triangles (1,2,3), (1,3,4) in (rotated) lon/lat; flips reported (first 10 + count),
      fixed in memory only; no element of the NEMO meshes is flipped or self-intersecting
- [x] `find_levels` / `find_levels_cavity` element loops, `elem_center`,
      `find_elem_neighbors_ini` (4 sides, `elem_edges` 1-2/2-3/3-4/4-1; the "≥ 2 neighbours"
      stop kept for all element types) and the `nod_in_elem2D` build → `nv` vertices
      — element-mean divisors stay literals (`/3.0` triangles, `/4.0` quads: LESSONS L7,
      `-no-prec-div`); `elem_center` = vertex mean (only orients `edge_tri`)
- [x] keep (audit only) the parts that are already quad-aware: `stiff_mat_ini` METIS graph with
      quad diagonals, `find_edges_ini` node neighbours, `gen_comm`
      — audit: no gap (`gen_comm` node/element halos and `mymesh`, `oce_local` `save_dist_mesh`
      edge lists all use `merge(3,4,…)` / 4 slots; `find_edges_ini` edges and `elem_edges` generic)
- [x] depth: `use_depthfile='nemo'` reads `depth.out` (positive values); trailing layer counts
      taken as given, no further level adjustment; `nlvls.out` written from them
      — `elvls.out` = min over the vertices (no isolated-cell iteration, no thresholding)
- [x] edge order interior-first / boundary-last and `edge2D_in` correct for mixed meshes
- [x] `check_dist.py`: every element owned once, halos symmetric between ranks, node/element/edge
      lists complete, `nlvls.out` equal to the depth-file counts, edge order + `edge2D_in`
      — also edges == element sides, `edge_tri`, receive lists from the owner; triangle meshes:
      `nlvls` == max of `elvls` around the node
- [x] tests: `run_ini_compare.sh` still identical for triangles; `check_dist.py` on the cleaned
      NEMO, mixed and cut meshes for 2, 8, 128 parts; the original NEMO mesh stops with the
      dead-end message (error case); a flipped element is fixed and reported
      — `tools/run_ini_mixed.sh`: cut meshes 2/8 parts, full NEMO + mixed 2/8/128; also
      `pi_w`, `core2_w` and (read-only) the old pi `dist_2/8` and CORE2 `dist_8` of Task 3:
      structurally correct; error cases: dead end, flipped quad + triangle (outputs identical
      to the unflipped mesh), bow-tie quad stops, check_dist detects a corrupted `dist_2`
- [x] run tests — must pass before Task 5 — `run_ini_compare.sh` PASS (triangles still
      byte-identical), `run_ini_mixed.sh` PASS (~3 min), scratch ctest 48/48

### Task 5: fesom_ini writes edge communication and double-halo nodes

**Files:**
- Modify: `mesh_part/src/oce_local.F90` (save_dist_mesh) + new list builders
- Modify: `tools/mesh/check_dist.py`, `tools/make_dist1.py`
- ➕ Modify: `tools/run_ini_compare.sh`, `tools/run_ini_mixed.sh` (appended sections, dist_1/4,
  reader and corruption tests)

- [x] marker file `dist_N/format.out` (version); existing files keep their FESOM2 layout
      — `format.out` = version 1 (written with `rpart.out`); the FESOM2 part of every
      `my_list*.out` / `com_info*.out` is byte-identical, the new lists are appended
      (module `dist_ext` in `oce_local.F90`)
- [x] `eXDim_nod2D`: vertices of elements 1..myDim+eDim+eXDim not among 1..myDim+eDim nodes;
      appended to `my_list*.out`; `com_nod2D_full` appended to `com_info*.out`
      — eXDim list grouped by owner, ascending within (local numbers myDim+eDim+1..);
      `com_nod2D_full` receives every halo node (eDim + eXDim) from its owner; per peer both
      sides ordered by global id
- [x] `com_edge2D`: unique owner = owner of `edges(1,edge)`, sends to every other rank holding the
      edge (also where it lies in that rank's 1..myDim_edge2D); appended to `com_info*.out`
      — holders = ranks of the vertices of the elements on the edge; local edge number =
      position in `myList_edge2D`; per peer ordered by global id
- [x] `make_dist1.py` writes the marker and empty new lists
- [x] regenerate `dist_1/2/8` (+ 128 for CORE2, NEMO) in all writable test meshes
      — pi_w 1/2/8, core2_w 1/2/8/128, nemo025q_fesom + nemo025mix 1/2/8/128, the cut-outs
      1/2/4/8 (dist_1 by `make_dist1.py`; made by `run_ini_compare.sh` / `run_ini_mixed.sh`)
- [x] tests: `check_dist.py` verifies the new lists (send/receive symmetric, every halo edge and
      double-halo node has exactly one owner; 4-rank case where owned edge ranges overlap) on
      pi, CORE2, NEMO, mixed and the cut-outs; FESOM2's reader still reads a new `dist_N`
      (existing files unchanged)
      — `run_ini_mixed.sh`: dist_4 of the cut-outs has 148 (quad) / 239 (mixed) edges in
      overlapping owned ranges, all received from the unique owner; corrupted eXDim list,
      `com_nod2D_full`, `com_edge2D` detected. `run_ini_compare.sh`: `rpart.out` identical and
      `my_list`/`com_info` = FESOM2 build's file + appended part; the model's existing
      FESOM2-sequence reader (`test_partit`, `test_halo`, np 2/8) reads the new pi_w dist_N
- [x] run tests — must pass before Task 6 — `run_ini_compare.sh` PASS, `run_ini_mixed.sh` PASS,
      `test_nemo_meshq.sh`, `test_mesh_tools.sh` PASS, scratch ctest 48/48

### Task 6: Model reads mixed meshes and the extended dist_N (triangles bit-identical)

**Files:**
- Modify: `src/mesh/mod_mesh_read.F90` (`read_mesh`, `read_mesh_local`)
- Modify: `src/infra/mod_partitioning.F90` (incl. `synthesize_1rank`), `src/types/mod_partit.F90`
- Modify: `src/types/mod_mesh.F90`
- Modify: `test/test_mesh.F90`, `test/test_partit.F90`, `test/CMakeLists.txt`
- ➕ Create: `test/test_mesh_read.F90`; Modify: `test/test_halo.F90`, `test/test_types.F90`
- ➕ Modify: `tools/mesh/cut_region.py`, `tools/mesh/test_mesh_tools.sh` (see ⚠️ below)

- [x] 3/4-column `elem2d.out` (column detection), `elem2D_nnodes` from 4th = 1st, element levels
      min/max over 4 slots; node levels only from `nlvls.out`
      — `read_elem_line` (both paths): columns from the first line, every line must match,
      3 columns → 4th := 1st, ids checked in 1..nod2D; triangles now store slot 4 = slot 1
      (was 0); `enforce_cw_orientation` turns a quad by swapping 2↔4 (triangles unchanged)
- [x] zbar from `aux3d.out`, or `depth.out` if `aux3d.out` is absent (positive accepted)
      — `open_depth_file`; depth.out node depths stored negated (aux3d sign; never used)
- [x] `edge_tri` keeps the on-disk encoding in this task (convention change is Task 7)
- [x] `format.out` present → read `eXDim_nod2D` list (+ coordinates/levels of those nodes),
      `com_nod2D_full`, `com_edge2D`; `synthesize_1rank` provides empty lists
      — `partit%dist_format` (0 = FESOM2 layout), `eXDim_nod2D`, ids appended to
      `myList_nod2D`; `coord_nod2D`, `geo_coord_nod2D`, `bc_index_nod2D`, `depth`,
      `nlevels_nod2D` sized myDim+eDim+eXDim (everything element-derived stays nNodL);
      1-rank: dist_format = 1 with empty structs
- [x] tests: read the NEMO and mixed cut-outs at np 1/2 (counts, nnodes histogram, levels vs
      `nlvls.out`); old pi `dist_N` still read; a malformed element line stops with a message
      (error case); pi and CORE2 runs BIT-IDENTICAL to the binary before this task
      — `test_mesh_read` (+ `test_partit`) on nemo025q_cut, nemo025mix_cut, pi_w at np 1/2,
      original pi dist_2, malformed case; bit-identity vs 9b0d3a6: production driver
      (production config, 4 steps, all-node dumps + every-step output) on pi np 1/2, pi_w np
      2/8, CORE2 np 2/8, core2_w np 8, and `fesom_conserve` (zstar+Redi+TKE, 20 steps) on pi
      np 1/2, pi_w np 2: every dump/output file and the log identical.
      ⚠️ The cut-outs violated the model's bottom invariant (23 / 17 nodes on the cut kept
      a layer count above every element left around them → `read_mesh` stops):
      `cut_region.py` now lowers such counts (and the node depth) to the deepest remaining
      element; cut-outs regenerated, `run_ini_mixed.sh` re-run (PASS)
- [x] freeze the binary of this task's commit as the triangle-regression reference
      — `/work/ab0995/a270029/fesom3_hybrid/frozen_task6/` (all 54 binaries + README)
- [x] run ctest — must pass before Task 7 — scratch ctest 64/64; `run_ini_compare.sh`,
      `run_ini_mixed.sh`, `test_mesh_tools.sh`, `test_nemo_meshq.sh` PASS

### Task 7: Mixed-mesh geometry, conventions, node-gather weights, analytic mixed patch

**Files:**
- Modify: `src/types/mod_mesh.F90` (`elem_w`, `nod_in_elem2D_w`, `elem_edges`,
  `elem_neighbors(4)`, `gradient_sca(8)`; binary mesh dump read/write layout)
- Modify: `src/mesh/mod_mesh_areas.F90`, `src/mesh/mod_mesh_rotate.F90` (cyclic wrap),
  `src/mesh/mod_mesh_analytic.F90` (4-slot edges, new conventions, mixed-patch generator)
- Modify: every `gradient_sca(4:6)` call site → `(5:7)` (pure index shift; `(1:3)` and 3-vertex
  lists stay — triangle semantics unchanged; the 4-slot conversion is done per kernel later)
- Modify: every boundary test listed in Context (src AND test files), in the same change
- Modify: `src/infra/mod_geom_dump.F90`, gate tools under `tools/` (see below)
- Create: `test/test_geom_mixed.F90`, `tools/run_tri_regression.sh`

- [x] geometry: centroid, area, kite weights `elem_w`, node area, `gradient_sca` (Green–Gauss),
      `elem_edges`, `elem_neighbors` (missing = self), `edge_tri` missing 2nd = 1st with
      `edge_cross_dxdy(3:4) = 0`, `compute_mesh_resolution` with kite weights
      — `mod_mesh_areas`: centroid/lon-lat area from the sub-triangles (1,3,4)+(1,2,3) (qq),
      `elem_area` = that area × `elem_cos`; kite weights normalised by the sum of the 4 side
      triangles (triangles exactly 1/3); Green–Gauss with the longitudes relative to vertex 1
      (one planar embedding, as FESOM2); `edge_tri` converted at read time (both paths) and in
      the analytic builder; `elem_w`/`nod_in_elem2D_w` added to the binary mesh layout.
      `mod_mesh_rotate` needed no change (the existing wrap helpers suffice)
- [x] `nod_in_elem2D_w(k,n)`: kite weight of node n in element `nod_in_elem2D(k,n)`
      — owned nodes computed, halo nodes by `exchange_nod` (the halo slot order is the owner's)
- [x] boundary tests converted in the same change (`el(2) /= el(1)`, `nb == el`); grep shows none left
      — all sites of Context + `test_bottom`, `test_vinv`, `test_muscl_onthefly`, meshdiag
      `edge_face_links`; the grep finds only the convention builders themselves (`edge_tri.out`
      encoding at read, the analytic builder) and the `-1` invariant check
- [x] `mod_mesh_analytic`: edges over `nv` sides, new conventions, generator for a small mixed patch
      — `generate_mixed_patch` (quads + triangle pairs, jittered interior nodes, convex);
      `build_edges` over nv sides, boundary edge_tri(2)=edge_tri(1), interior edges first
      (`edge2D_in`); triangles store slot 4 = slot 1
- [x] gate migration decision recorded in HANDOFF: FESOM2-oracle byte gates (`tools/run_*_gate*.sh`,
      `geom_diff.py`, `mod_geom_dump`) are retired or remapped to the new layout and listed
      — HANDOFF "Gate migration for hybrid meshes": all FESOM2-oracle gates + oracle-dump helpers
      retired (listed); `mod_geom_dump` remapped (FESOM2 6-row `gradient_sca` from slots 1:3/5:7);
      FESOM3-only gates (restart, output, conserve, ctest) kept
- [x] `test_geom_mixed` (TDD): on pi, the NEMO-mixed cut-out and the analytic mixed patch —
      Σ node area = Σ element area; Σ `elem_w` = 1 per element; Σ_k `nod_in_elem2D_w`·A_e = node
      area; `gradient_sca` exact for linear fields on triangles and quads; cross-edge vectors
      close around each element; triangle geometry equals the previous formulas to roundoff
      — np 1/2 on `nemo025mix_cut` (+ the patch at np 1) and `pi_w` (rotated 50/15/-90; the
      unrotated pi has elements around the geographic pole); closure tested as
      A·g_j = rot(d_j − d_{j−1}) with the stored `edge_cross_dxdy`; conventions checked too;
      errors ≤ 1.4e-14 (G1–G3), ≤ 8e-15 (G4), ≤ 1.3e-13 (G5), triangles vs P1 ≤ 8e-15
- [x] adapt existing tests whose bitwise reference was the old triangle arithmetic to tolerances
      — none needed: every existing test passes with the new geometry (only the convention
      changes above were needed)
- [x] `run_tri_regression.sh` (pi, CORE2: a few steps vs the frozen Task-6 binary, roundoff) — run
      — pi np 1/2 + CORE2 np 8 production config 4 steps, `fesom_conserve` pi np 1/2: PASS. Metric
      max|new−ref|/max|ref| per field: step 1 before the first SSH solve ≤ 1e-9 (measured ≤
      1.3e-11), all substeps/outputs ≤ 1e-2 (the CG stops at relative residual 1e-5; measured pi
      ≤ 6e-9, CORE2 ≤ 5.4e-4 in Kv/Av), conservation 1e-12
- [x] run ctest — must pass before Task 8 — scratch ctest 68/68; `run_ini_compare.sh`,
      `run_ini_mixed.sh`, `test_mesh_tools.sh`, `test_nemo_meshq.sh` PASS
- [x] ➕ follow-up: tighten `run_tri_regression.sh` to roundoff and explain the pi np1≠np2
      difference — env knob `FESOM3_SSH_SOLTOL`/`FESOM3_SSH_MAXITER` in `oce_ssh_solve.F90`
      (default unchanged); reference rebuilt as `frozen_task6_soltol` (5b2a7f3 + only the knob,
      byte-identical to `frozen_task6` with the knob unset); all runs solve to 1e-14; tiered
      metric `max|Δ|/max|ref|` (pi / core2): PRE 1e-9, EARLY (steps 1–3 pi, 1 core2) 1e-10 / 1e-8,
      OCEAN 1e-6 / 2e-4, ICE 1e-5 / 2e-3, MIX (Kv, Av) 1e-6 / 1e-2; measured pi EARLY ≤ 9.2e-13,
      OCEAN ≤ 3.3e-8, CORE2 EARLY 2.1e-9 (ice), OCEAN 4.4e-6 dumps / 2.2e-5 float32 u/v. np1≠np2 on
      pi = IC extrapolation (`extrap_nod3D` Gauss–Seidel with per-rank halo snapshot, FESOM2 design,
      not fixed); new `decomp:pi_w:2/8` cases compare np 1 vs np N from one checkpoint (≤ 2.5e-11).
      HANDOFF + LESSONS L60

### Task 8: Edge exchange and double-halo node exchange

**Files:**
- Modify: `src/infra/mod_halo.F90`, `src/infra/mod_part_bounds.F90` (`nNodX`)
- Create: `test/test_halo_ext.F90`
- ➕ Modify: `test/CMakeLists.txt`, `tools/run_ini_compare.sh` (dist_4 of pi_w/core2_w),
  `tools/run_ini_mixed.sh` (check_dist on them)

- [x] `exchange_edge` (2D, 3D) over `com_edge2D`; `exchange_nod_full` (2D, 3D) over
      `com_nod2D_full`, same buffer/Isend/Irecv pattern as `exchange_nod`
      — real 2D/3D wrappers on the existing broadcast-only cores (`core_2D_r`, `core_3D_r`);
      `dist_format == 0` stops (MPI_Abort) with "<list> is not available: dist_N has no
      format.out"; 1-rank synthesis = empty lists = no-op
- [x] `owned_bounds` returns `nNodX`
      — new OPTIONAL last argument `nNodX` (myDim+eDim+eXDim_nod2D; npes 1: `nod2D`), the
      139 existing calls unchanged
- [x] tests: round trip (fill owned with global id, exchange, compare) on the regenerated pi and
      CORE2 copies and the NEMO cut-out at np 2, 4, 8, incl. edges owned on two ranks; calling
      without `format.out` stops with the message (error case)
      — `test_halo_ext`: nodes owned = gid → `exchange_nod_full` 2D/3D → all nNodX slots ==
      `myList_nod2D`; edges set only on the unique owner (owner of `edges(1,·)` from
      edges.out, node owners gathered) → `exchange_edge` 2D/3D → all myDim+eDim edges == gid;
      receive lists only halo slots / foreign-owned edges; owned-range edges of another owner
      counted (> 0 required; core2_w np4 1190, mixed cut np4 227 = check_dist). ctest: np 1
      (synthesis), nemo025q_cut 2/4/8, nemo025mix_cut 4, pi_w 2/4/8, core2_w 2/4/8, original
      pi dist_2 (no format.out) stops (PASS regex). dist_4 of pi_w/core2_w generated by
      `run_ini_compare.sh` (== FESOM2 build; CORE2 dist_4 == the existing one, now gated),
      checked in `run_ini_mixed.sh`
- [x] run ctest — must pass before Task 9 — scratch ctest 80/80 (one earlier full run had a
      transient `test_mesh_np1` failure, passed on rerun and in 3 further full runs);
      `run_ini_compare.sh`, `run_ini_mixed.sh` PASS (pi_w/core2_w dist_1/2/8/128 unchanged)

### Task 9: SSH stiffness matrix for mixed meshes

**Files:**
- Modify: `src/oce/oce_ssh_rhs.F90` (`init_stiff_mat_ale`, `update_stiff_mat_ale`)
- Create: `test/test_ssh_stiff.F90`
- ➕ Modify: `test/CMakeLists.txt`, `src/mesh/mod_mesh_analytic.F90` (`generate_quad_patch`,
  FESOM edge orientation in `build_edges`)

- [x] sparsity = node plus all vertices of all elements around it (dynamic, no fixed 12); the
      diagonal stays the first CSR entry (Jacobi preconditioner relies on it)
      — row = self, edge neighbours in edge order (as before), then the remaining vertices of
      the elements around the node (`nod_in_elem2D`, marker array); row length counted first
      (≤ 1 + #edges + #elements); triangle meshes keep the identical CSR
- [x] assembly over edges with 4 coefficients (`gradient_sca(1:4|5:8)`, `npos(4)`)
      — init and `update_stiff_mat_ale`; scalar scatter loop over the 4 slots (a triangle's
      npos(4)=npos(1) with coefficient 0 — a vector-subscript assignment would be many-one)
- [x] tests: matrix symmetric; row sums = `area/dt` term; quad diagonal entries present and
      nonzero; triangle-mesh matrix equals the previous one to roundoff; CG converges on the NEMO
      cut-out; checkerboard (+1,−1,+1,−1) on a parallelogram quad patch is in the null space of
      the stiffness part (documented property, asserted so it can't change silently)
      — `test_ssh_stiff` S1–S8 on nemo025q_cut np1/2, nemo025mix_cut np1/2/8, pi_w np1/2, core2_w
      np8 (+ analytic mixed and parallelogram patches at np1): sym ≤ 1.4e-13 (CORE2), row sums
      ≤ 3.4e-16, quad diagonal ≥ 1.4e-2 of the row's largest off-diagonal, triangles (pi_w,
      core2_w) CSR identical and values BITWISE equal to the previous assembly, update with
      dhe = −H/2 adds K/2 (≤ 3e-15), CG 224 its (np1) to 2.5e-13. Checkerboard: |K f|/|K| ≈ 1e-16
      on the parallelogram patch AND on the (2-colourable) NEMO quad cut-out — the Green–Gauss
      coefficient of vertex j ∝ x_{j+1}−x_{j−1}, so c₁ = −c₃, c₂ = −c₄ on ANY quad: the mode is an
      exact null mode of K on every all-quad mesh, not only on parallelograms (input for Task 17).
      ➕ `build_edges` of the analytic generator now orients edges as FESOM (edge_tri(1) left of
      edges(1)→edges(2)); the old lo→hi order made the patch matrices non-symmetric
- [x] run ctest — must pass before Task 10 — scratch ctest 88/88; `run_tri_regression.sh` 7/7

### Task 10: Ocean dynamics kernels

**Files:**
- Modify: `src/oce/oce_dyn_velrhs.F90`, `src/oce/oce_ale.F90`, `src/oce/oce_pgf.F90`,
  `src/oce/oce_dyn_ivertvisc.F90`, `src/oce/oce_dyn_vinv.F90`, `src/oce/oce_dyn_visc.F90`
- Modify: `test/test_vinv.F90`, `test/test_ivertvisc.F90`
- Create: `test/test_pgf_mixed.F90`
- ➕ Modify: `test/CMakeLists.txt`, `CMakeLists.txt` + `tools/run_tri_regression.sh` (MPI launch
  flag, see the ⚠️ finding below), `docs/LESSONS.md` (L61)

- [x] element means via `elem_w`; gradients over 4 slots; vertex loops over 4 slots / `nv`
      — `update_vel`, `compute_vel_rhs` (SSH gradient, `gradient_sca(1:4|5:8)`), `dhe` and the
      zstar `helem` (`sum(elem_w*…)`), `impl_vert_visc_ale` face means of `w_i`, the scalar
      momentum vertice→element mean; `oce_dyn_visc` needed no change (edge-based, boundary by
      `edge2D_in`, `len = sqrt(sum(elem_area(el)))`)
- [x] node gathers via `nod_in_elem2D_w`: `UVnode` (`oce_ale.F90:96-109`), `w·du/dz`
      (`oce_dyn_velrhs.F90:209-231`) — `Σ A_e·w·x / area` (the hidden `3·` removed)
- [x] `oce_pgf`: per-vertex arrays and loops (`idx`, `dx10`, `do ni`, `sum(drho_dz)/3`) to 4 slots
      — full-cell and Shchepetkin; vertex `drho_dz` mean with `elem_w`; FESOM2's
      `(/1.0,1.0,1.0/)*z_n` → the scalar `z_n`
- [x] `oce_dyn_vinv`: KE and vorticity node↔element with `elem_w`; `rvo_upwind` over
      `elem_neighbors(1:4)` and `nv` side normals oriented by the centroid
      — KE scatter as a scalar 4-slot loop (a triangle's repeated slot); side normals over 4
      slots (a triangle's side 4 is degenerate → normal 0, never weighted), inward = toward the
      area-weighted centroid of (1,3,4)+(1,2,3)
- [x] tests: existing exact tests (uniform flow, linear shear) on triangle and mixed patches;
      `test_pgf_mixed`: resting stratified ocean gives zero pressure-gradient force on the mixed patch
      — `test_vinv`/`test_ivertvisc` generalised to 4 slots + `FESOM3_TEST_PATCH=mixed|quad` /
      `FESOM3_TEST_ROTATE=0`; new ctests: vinv on the mixed and parallelogram patches, NEMO-mixed
      cut np1/2, NEMO-quad cut np2; ivertvisc on the mixed patch and NEMO-mixed cut np2. Finding:
      the median-dual circulation with element-constant velocity is exact for linear fields on
      triangles and PARALLELOGRAMS only (V2 3.9e-15 on the skewed parallelogram patch, 3.7e-14 on
      the NEMO cut-outs = lon/lat rectangles; 6e-2 on the jittered general quads, first order) —
      V2/V9a assert roundoff at "exact" nodes, first-order consistency (< 0.1) elsewhere; V9a
      normalised by the vorticity term (a pure shear's tendency is itself 0 on planar meshes).
      `test_pgf_mixed` (NEMO-mixed cut np1/2 + patch, NEMO-quad cut np2, pi_w np1): P1 resting
      linfs PGF = 0, P2 linear pressure exact, P3 Shchepetkin on tilted levels with linear
      stratification = 0, P4 horizontal density gradient exact, N1 uniform-velocity node gather
      exact, N2 `update_vel` SSH gradient exact — all ≤ 5e-13; mutation check: the old `/3`
      forms fail P3 (0.12) and N1 (0.33)
- [x] `run_tri_regression.sh` — roundoff level — 7/7 PASS, levels as after Task 7 (pi PRE 6.3e-15,
      EARLY 9.2e-13, OCEAN 1.0e-7; CORE2 EARLY 2.1e-9, OCEAN 2.5e-5, ICE 1.2e-4, MIX 2.4e-4;
      decomp np1 vs np2/8 OCEAN 6.4e-10; conservation ≤ 3e-14)
- [x] run ctest — must pass before Task 11 — scratch ctest 99/99
- ⚠️ (fixed in the test launchers) the login-node Open MPI 4.1.2 `vader` single-copy mechanism
  (`knem`, default) delivered STALE halo payloads (right count, data of an earlier exchange) when
  the per-call send buffer reused a heap address — seen as `test_vinv` failures on the NEMO quad
  cut np2; `--mca btl_vader_single_copy_mechanism none` added to `FESOM3_MPI_LAUNCH_FLAGS` and
  `run_tri_regression.sh` (LESSONS L61). Batch (UCX) runs are not affected; other login-node
  scripts with `--mca btl self,vader` need the flag too
- [x] ➕ follow-up: root cause of the stale halo payload — code bug (non-blocking MPI on a
      temporary / freed buffer) or MPI library? — **Open MPI 4.1.2 bug, code clean.** vader's knem
      registration cache (keyed by virtual address) misses glibc heap trims (`brk`), so a buffer
      re-allocated at the same address is read from the old pinned pages. Standalone reproducer
      `tools/mpi_knem_stale_repro.F90` + `tools/run_mpi_knem_stale_repro.sh` (correct contiguous
      Isend/Irecv/Waitall): 788/800 stale under knem with heap trimming; 0 with single-copy
      `none`/`emulated`, `mpi_leave_pinned 0`, no trim, or `pml ucx`; `test_vinv` quad cut np2
      likewise passes under knem once trimming or the cache is off. `mod_halo` audit + an
      `-check arg_temp_created` build: no argument temporaries anywhere (test_vinv, prod pi np2,
      CORE2 np8). Workaround kept (also added to `run_conserve_pi.sh`, `run_bottom_audit.sh`),
      comments in `CMakeLists.txt`/`env.sh`/`run_tri_regression.sh` corrected, L61 rewritten (L35
      was the same bug), HANDOFF: production (UCX + `env.sh`) not affected, not a cause of the
      dotty `w`. ctest 99/99, `run_tri_regression.sh` 7/7

### Task 11: Tracer kernels and MUSCL up/downwind elements for quads

**Files:**
- Modify: `src/oce/oce_tracer_grad.F90`, `src/oce/oce_ale_tracer.F90`, `src/oce/oce_adv_tra_fct.F90`,
  `src/oce/oce_adv_tra_flux.F90`, `src/oce/oce_adv_tra_hor.F90`, `src/oce/oce_adv_tra_driver.F90`,
  `src/oce/oce_muscl_adv.F90`, `src/oce/oce_pressure_bv.F90`, `src/oce/oce_initial_state.F90`
- Modify: `test/test_muscl_onthefly.F90`
- Create: `test/test_updown_mixed.F90`
- ➕ Modify: `test/CMakeLists.txt`; `src/mesh/mod_mesh_read.F90` (comment: subroutine rename)

- [x] element gradients/means over 4 slots; FCT element bounds over 4 slots
      — `tracer_gradient_elements` `gradient_sca(1:4|5:8)`; FCT a2 max/min over 4 slots (a
      triangle's repeated slot leaves max/min unchanged); `oce_adv_tra_hor/flux/driver` needed
      no code change (edge-based, boundary tests from Task 7) — comments triangle → element
- [x] node gathers via `nod_in_elem2D_w`: `tr_xynodes`, `compute_sigma_xy`, `smooth_nod`
      — `tr_xynodes` `Σ A·w·g / areasvol` (hidden `/3` removed); `compute_sigma_xy` weights
      and `vol` = kite shares `elem_area·nod_in_elem2D_w`, gradients over 4 slots; `smooth_nod`
      = kite-share-weighted mean of the `elem_w` element means (triangles: FESOM2's
      `ΣA(a1+a2+a3)/(3ΣA)` to roundoff)
- [x] `find_up_downwind_triangles` → elements: at vertex k the angle spanned by the two adjacent
      sides (next/previous vertex), triangles and quads (qq `oce_muscl_adv.F90:56-255`); the
      multi-rank `coord_elem(2,3,·)` / `e_nodes(3,·)` exchange becomes 4 vertices
      — renamed `find_up_downwind_elements`; one helper `elem_sides` (nv from slot 4 == slot 1,
      b = next, c = previous) for the 1-rank and multi-rank paths; the exact "x along a side"
      compares (`ab == ax`, `ax == 0`) became `|·| < 1e-10` rad (on structured quad grids the
      edge direction is often collinear with a neighbour's side; qq uses 0.01); the array name
      `edge_up_dn_tri` kept (t_tracer_work binary layout, drivers)
- [x] `oce_initial_state` extrapolation averages over the element's UNIQUE vertices (`nv`)
- [x] tests: `test_updown_mixed` — up/downwind element found for every interior edge of the mixed
      patch and the edge direction lies in its vertex angle; MUSCL equals the reference rule;
      linear tracer advected exactly on uniform flow (interior)
      — `test_updown_mixed` (NEMO-mixed cut np1/2/8 + analytic mixed and parallelogram patches
      at np1, NEMO-quad cut np2, pi_w rotated np2): U1 element found at every interior edge node
      (interior = vertex angles sum to 2π; 0 missing), U2 found element around the anchor and the
      direction inside its vertex angle by cross products (0 outside), L1 linear tracer: element
      gradient exact (≤ 1.2e-13) and the MUSCL face values = the exact edge-midpoint value on
      every edge with both elements, both flow signs (≤ 2.2e-12), L2 node tendency = −u·∇T at
      all-triangle (≤ 8e-14) and all-parallelogram (≤ 4.3e-14) interior nodes; MUSCL == reference
      rule: `test_muscl_onthefly` part R (bitwise) now also on the unrotated NEMO-mixed cut
      np1/2 and NEMO-quad cut np2 (`FESOM3_TEST_ROTATE=0`; synth_cavity over 4 slots)
- ⚠️ finding (scheme property, not a kernel bug; relevant for Task 17): the median-dual edge
  scheme is linearly exact only where the dual cell is all triangles (Galerkin-P1 equivalence)
  or all rectilinear parallelograms. At triangle–quad interface nodes and on general quads the
  node tendency of a linear tracer in uniform flow has an O(1) LOCAL relative error although
  the MUSCL face values are exact: 0.22 at the NEMO-mixed interface nodes, 0.24 on the jittered
  mixed patch (the face value is taken at the edge midpoint, but the dual face runs
  midpoint→centroid; the errors cancel edge by edge only for those two shapes). Asserted
  bounded (< 0.5) in L2; a remedy (e.g. face values at the dual-segment midpoints) is a design
  decision for the user
- [x] run ctest — must pass before Task 12 — scratch ctest 107/107; `run_tri_regression.sh` 7/7
      (pi PRE 6.3e-15, EARLY 9.2e-13, OCEAN 1.0e-7 as after Task 10; CORE2 EARLY 2.1e-9, OCEAN
      3.9e-5, ICE 1.2e-4, MIX 2.4e-4; decomp OCEAN 6.4e-10; conservation <= 4.3e-16 totals)

### Task 12: Mixing and GM/Redi

**Files:**
- Modify: `src/oce/oce_fer_gm.F90`, `src/oce/oce_mixing_kpp.F90`, `src/oce/oce_mixing_tke.F90`,
  `src/oce/oce_ale_mixing_pp.F90`, `src/oce/oce_mo_conv.F90`

- ➕ Modify: `src/drivers/fesom_conserve.F90` (helem invariant = `elem_w` mean), `test/CMakeLists.txt`

- [x] element means via `elem_w` (incl. `fer_uv`, `oce_fer_gm.F90:296-318`)
      — `fer_gamma2vel` `sum(elem_w*(Γ_nz−Γ_nz+1))/helem` (the `onethird` constant gone);
      KPP `viscAE`, TKE `Av`, PP `Av = mix_coeff_PP·sum(elem_w·Kv²) + A_ver`; `mo_convect`
      maxval of the convective ramp over all 4 slots (a triangle's repeated slot leaves the max
      unchanged); `elnodes(MAX_NV)` everywhere, comments "3 corner nodes" → vertices. Not
      changed: the GM/Redi resolution scaling `(area(n)·2/refres²)^(1/2)` (FESOM2 tuning; the
      factor 2 = element/node area ratio of triangles — on all-quad meshes it overstates the
      resolution by √2; a user decision, recorded in the progress file); the KPP blmc smoothing
      already goes through `smooth_nod` (Task 11)
- [x] tests: KPP/TKE/PP/GM conservation variants (`fesom_conserve`) on the mixed cut-out
      — 10 new ctests `fesom_conserve_mixcut_{PP,KPP,TKE,GM,REDI}_np{1,2}` (zstar, 20 steps,
      PHC IC on the closed NEMO-mixed cut-out, tol 1e-12): measured |drift| ≤ 1.3e-14 heat/salt,
      ≤ 9.3e-15 volume; the driver's helem invariant now checks the `elem_w` mean (with the old
      `/3` check it is wrong on quads)
- [x] run ctest — must pass before Task 13 — scratch ctest 117/117; `run_tri_regression.sh` 7/7
      (pi PRE 4.3e-11 (Kv; was 6.3e-15 — the TKE/PP Av element mean now goes through elem_w,
      thresholds amplify), EARLY 3.2e-12, OCEAN 7.4e-8; CORE2 EARLY 7.5e-9, OCEAN 2.8e-5, ICE
      1.0e-4, MIX 6.6e-4; decomp OCEAN 6.4e-10 dumps/7.0e-9 float32 w; conservation totals
      ≤ 4.3e-16)

### Task 13: Forcing and ice–ocean coupling

**Files:**
- Modify: `src/forcing/mod_forcing_bulk.F90`, `src/ice/mod_ice_oce_coupling.F90`
- ➕ Create: `test/test_coupling_mixed.F90`; Modify: `test/CMakeLists.txt`

- [x] element means via `elem_w`; fix the `elnodes(3)` shape bug in the coupling (:105)
      — `oce_fluxes_mom` and `forcing_stress_surf`: `stress_surf = sum(elem_w*stress_node_surf)`
      over the 4 slots, `elnodes(MAX_NV)` (the coupling assigned the 4-slot `elem2D_nodes(:,e)`
      to `elnodes(3)`); comments updated. The rest of the forcing (bulk coefficients, wind/ice
      stress, JRA55 interpolation) and `oce_fluxes` are per node — no change
- [x] tests: forcing/coupling self-tests on CORE2 (roundoff) and the mixed cut-out
      — `test_coupling_mixed` (node fields = hashes of the global node id): C1 uniform node stress
      → same element stress (the former /3 mean gives 1/3 error on quads, asserted as mutation
      check); C2 momentum input Σ_e A_e·stress_surf == Σ_n area·stress_node_surf with partial ice
      and ice-ocean drag (each element counted once, on the owner of its 1st vertex); C3 triangles
      == the former sum/3 (≤ 3e-16 rel.); C4 `forcing_stress_surf` == `oce_fluxes_mom` at a_ice = 0
      (bitwise, np 1); C5 `oce_fluxes` balancing (∫virtual_salt = ∫relax_salt = 0, freshwater
      balance, ocean_area = Σ areasvol) — NEMO-mixed cut np1/2 (+ analytic mixed patch at np1),
      NEMO-quad cut np2, pi_w np1, core2_w np8 (rotated): all ≤ 9e-15. CORE2 roundoff:
      `run_tri_regression.sh` 7/7 (its prod cases run JRA55 + sea ice + this coupling; CORE2
      EARLY 2.1e-9, OCEAN 2.5e-5, ICE 1.2e-4, MIX 2.5e-4; pi PRE 6.3e-15, EARLY 9.2e-13, OCEAN
      1.1e-7; decomp OCEAN 2.5e-9; conservation totals ≤ 3.5e-16). A forcing+ice RUN on the
      cut-out waits for Tasks 15/16 (ice dynamics/FCT still triangle-only) — Task 17
- [x] run ctest — must pass before Task 14 — scratch ctest 122/122 (wiped + reconfigured)

### Task 14: IO, mesh diagnostics and dump drivers

**Files:**
- Modify: `src/io/mod_io_coords.F90`, `src/io/mod_io_meshdiag.F90`, `src/io/mod_io_means.F90`,
  `src/drivers/fesom_pressuredump.F90`, `src/drivers/fesom_conserve.F90`,
  `src/drivers/fesom_advhordump.F90`, `fesom_geomdump`, `fesom_meshdiagdump`, `fesom_evpdump`,
  `fesom_icefctdump` drivers
- ➕ Modify: `src/mesh/mod_mesh_areas.F90` (`elem_center` public), `src/infra/mod_geom_dump.F90`,
  `src/drivers/fesom_advhordump_mr.F90` (status note), `test/CMakeLists.txt`,
  `tools/run_tri_regression.sh` (element-vector comparison), `docs/HANDOFF.md`;
  Create: `tools/check_meshdiag.py`

- [x] element coordinates = centroids; mesh diagnostic output with 4-column `face_nodes`
      (triangle 4th = UGRID `_FillValue` −1), `face_edges`/`face_links` from the now-built
      `elem_edges`/`elem_neighbors`
      — `io_coords_compute` (means + restart stores, and the r2g point of element vectors) uses
      `mod_mesh_areas elem_center` (area-weighted, cyclic-wrapped; FESOM2's plain unwrapped 3-vertex
      mean was on the wrong side of the globe for the elements across the rotated dateline: pi 77,
      CORE2 566; all others agree to < 1e-10 deg). Meshdiag: `face_nodes`, `face_edges` (fill −1),
      `face_links` (missing neighbour and triangle 4th −999 = FESOM2's fill), `gradient_sca_x/y` (4
      coefficients, triangle 4th 0), all dim `n4`, Zarr `fill_value` declared; topology gains
      `face_edge_connectivity`/`face_face_connectivity`; global ids (partition-independent)
- [x] dump drivers follow the gate decision of Task 7 (retired or remapped)
      — none deleted, each a diagnostic dump with its byte gate retired (status note in each header,
      table in HANDOFF): `mod_geom_dump` FESOM2 layout incl. `edge_tri(2)=0`, triangle-only (stops
      on quads); `fesom_meshdiagdump` = hybrid layout; `fesom_pressuredump` 4-slot `elem_w` means
      (KPP `viscAE`, unstable-element count); advhordump(_mr)/evpdump/icefctdump have no
      triangle code of their own; `fesom_conserve` already converted (Task 12)
- [x] grep check: no `elem2D_nodes(1:3`, `gradient_sca(1:3`, `(5:7`, element `/3` average or
      hidden-1/3 node gather left in `src/` (documented exceptions only, e.g. the spline
      coefficient `oce_adv_tra_ver.F90:134`, KPP's cube roots)
      — left outside the ice: `mod_geom_dump` (FESOM2 oracle layout, triangle-only guard),
      `mod_mesh_analytic:122-138` (the generator builds triangles, slot 4 set after),
      `mod_mesh_areas:208-209` + `oce_dyn_vinv:270-271` (sub-triangle centroid `/(3(s1+s2))`),
      `mod_mesh_areas:283` (triangle kite weights 1/3), spline `oce_adv_tra_ver:134-135`, KPP cube
      roots. Pending Task 15/16 (not converted here; the `mod_ice_dyn` ones cleared in Task 15, the `mod_ice_fct`/`mod_ice_setup` ones in Task 16): `mod_ice_dyn.F90` (:160-167 strain rates,
      :219/:514 val3, :269/:468 `elnodes(3)`, :342-361 msum/asum/aa/elevation gradient, :530-535,
      :570, :596-598), `mod_ice_fct.F90` (:73, :96-102, :292 `icoef(3,3)`, :335, :414, :454, :474,
      :488, :502), `mod_ice_setup.F90` (:179, :186)
- [x] tests: meshdiag file has 4 columns with −1 for triangles; `run_tri_regression.sh` — roundoff
      — ctests `meshdiag_{mixed,quad,pi_w}` (`tools/check_meshdiag.py`, stdlib python: runs
      `fesom_meshdiagdump` at np 1,2,8 / 1,2 / 1,2; M1 shapes/dims/fill/topology, M2 face_nodes ==
      elem2d.out vertex sets with −1 4th for triangles, M3 face_edges join nodes k,k+1 per
      edge_nodes, M4 face_links == other element of edge_face_links and mutual, M5 triangle 4th
      gradient 0 and Σ coefficients ≈ 0 (≤ 1.9e-16); np N store == np 1 (ints exact, reals
      1e-12); mutation check: a repeated 1st vertex or swapped links fail). `run_tri_regression.sh`
      7/7 — element vectors compared by component away from the moved (rotated-dateline)
      centroids and by speed on them (documented in its header; the moved set must be > 1 deg
      jumps and ≤ 5 %): pi PRE 6.3e-15, EARLY 9.2e-13, OCEAN 1.1e-7 (u/v elsewhere; speed at the 77
      moved 2.7e-8); CORE2 EARLY 2.1e-9, OCEAN 2.5e-5 (speed at the 566 moved 6.6e-6), ICE 1.2e-4,
      MIX 2.5e-4; decomp OCEAN ≤ 2.5e-9; conservation ≤ 3.5e-16
- [x] run ctest — must pass before Task 15 — scratch ctest 125/125

### Task 15: Sea-ice dynamics (EVP / mEVP)

**Files:**
- Modify: `src/ice/mod_ice_dyn.F90`
- Create: `test/test_ice_dyn_mixed.F90`
- ➕ Modify: `test/CMakeLists.txt`, `tools/run_tri_regression.sh` (CORE2 EARLY tolerance, header)

- [x] strain rates over 4 slots of `gradient_sca`; metric term with the `elem_w` mean
      — `stress_tensor` and the fused mEVP loop: `gradient_sca(1:4|5:8)` over `elem2D_nodes(:,el)`,
      metric terms `metric_factor·sum(elem_w·v)` (FESOM2 `sum(v)/3`)
- [x] `stress2rhs` (EVP and the fused mEVP loop): node shares `elem_area·elem_w(k)` and
      `gradient_sca(k|k+4)`; SSH-tilt terms (`aa = g·elem_area/3`, mEVP `bb`) → `elem_area·elem_w(k)`;
      `msum`/`asum` via `elem_w`
      — scatter loops over `MAX_NV` (a triangle's slot 4 has gradient and kite weight 0, adds
      exactly 0), metric share `sigma·elem_w(k)·metric_factor`; `any(m_ice/a_ice <= 0)` over 4
      slots (repeated slot unchanged); `elnodes(MAX_NV)`, `val3`/`meancos` gone; the mod_ice_dyn
      grep hits listed under Task 14 are cleared. `stress_tensor`/`stress2rhs` made public (test)
- [x] ocean-velocity node gather (`mod_ice_dyn.F90:100-110`) via `nod_in_elem2D_w`
      — `Σ A·w·UV / Σ A·w` over the wet elements (FESOM2 `Σ A·UV / Σ A`, no hidden 1/3: the
      normalisation by the wet shares keeps the cavity partial-patch behaviour)
- [x] tests: `test_ice_dyn_mixed` — zero strain rate for uniform drift and rigid rotation on
      triangle and mixed patches; CORE2 ice step at roundoff vs before
      — `test_ice_dyn_mixed` (NEMO-mixed cut np1/2 + analytic triangle, mixed and parallelogram
      patches at np1, NEMO-quad cut np2, pi_w np1 and core2_w np8 rotated): D1 uniform drift /
      rigid rotation zero strain (≤ 1.7e-15), linear deformation exact (≤ 3.7e-15); on the sphere
      rigid rotation about the pole: the metric terms cancel du/dy (max|eps|/Ω 3.6e-4 NEMO-mixed,
      1.5e-7 NEMO-quad, 3.3e-2 pi, 7.4e-3 CORE2 vs 0.10–0.50 without metric — discretisation
      error, identical to FESOM2 on triangles); triangles == former P1 strain formulas ≤ 1.1e-16;
      D2 `stress2rhs` = exact negative adjoint of the strain rates (Σ u·R = −Σ A σ:ε, ≤ 2.3e-17)
      and == the former formula at all-triangle nodes (0); D3 SSH tilt EVP/mEVP: Σ area·rhs_a =
      −gΣAη_x (global) and rhs_a = −g·a for a linear η on the patches (≤ 2.6e-15); D4 uniform
      strain is force-free at interior nodes after one subcycle, EVP and mEVP (≤ 1.4e-13; stress
      861/5.7 Pa); D5 ocean2ice uniform velocity exact, triangles == former gather (≤ 4.2e-16).
      Mutation: the former 3-slot scatter fails D2/D3/D4 on the mixed and quad patches and the
      NEMO-mixed cut (13 failures), triangle patch unaffected. CORE2 ice step: `fesom_evpdump`
      (CORE2 1 rank, 120 subcycles, prescribed forcing) new vs a build of 5b7dae1 (before this
      task): EVP uice 1.7e-6, sigma ≤ 2.0e-4; mEVP uice 3.8e-11, sigma ≤ 5.6e-9 — the same floor
      as the frozen Task-6 binary vs 5b7dae1 (identical ice code, Task-7 geometry roundoff only):
      EVP uice 8.8e-7, sigma ≤ 1.0e-4, mEVP 1.5e-11 / 2.1e-9 — roundoff amplified by the stiff
      EVP; srfoce_u/v ≤ 7.4e-16. `run_tri_regression.sh` 7/7 after raising the CORE2 EARLY
      tolerance 1e-8 → 1e-7: CORE2 EARLY T 2.2e-8 (was 2.1e-9; roundoff now enters inside every EVP
      subcycle), OCEAN 3.1e-5, ICE 8.9e-5, MIX 7.9e-4; the same case with mEVP in both runs gives
      EARLY T 1.1e-15, ICE 6.9e-6 (amplification by the standard EVP, not a kernel error); pi PRE
      6.3e-15, EARLY 9.2e-13, OCEAN 1.1e-7, ICE 6.6e-7; decomp OCEAN ≤ 4.4e-10; conservation totals
      ≤ 3.5e-16
- ⚠️ finding (pre-existing FESOM2 property, not a kernel bug): uniform ice AT REST is not a
  numerically stable state of the standard EVP (delta_min-regularised viscosity): on the plain
  triangle patch roundoff grows to O(1e-2) m/s within 120 subcycles; mEVP has zero stress at rest.
  It is the amplifier behind the CORE2 ICE regression tier (and the EVP dump floor above)
- [x] run ctest — must pass before Task 16 — scratch ctest 130/130 (wiped + reconfigured)

### Task 16: Sea-ice FCT and mass matrix (user's specification)

**Files:**
- Modify: `src/ice/mod_ice_setup.F90` (`ice_mass_matrix_fill`, loops at :187/:192,
  `fct_fluxes` allocation :137), `src/ice/mod_ice_fct.F90`, `src/types` ice work type (lumped mass)
- Create: `test/test_ice_mass.F90`
- ➕ Modify: `test/CMakeLists.txt`

- [x] `ice_mass_matrix_fill`: triangle S/6, S/12; quad (S/36)(4,2,1,2) rotated per row; on the new
      SSH pattern; not-found column stops with an error
      — element contributions `elem_area·ice_mcoef(n,q,nv)` over the `nv` unique vertices; the
      coefficient tables `ice_mcoef`/`ice_acoef(MAX_NV,MAX_NV,3:4)` are parameters in `mod_ice_fct`
      (triangle slot 4 = 0, one code path for both element types); a column missing from the CSR
      row → message (global ids) + `MPI_Abort`
- [x] FCT lumped mass = row sums (S/3, S/4) stored once; it replaces `mesh%area` in the
      low-order step, the high-order iteration and the antidiffusive fluxes
      — `ice%work%fct_lumpedmass` (nodes) = row sums of `fct_massmatrix`, exchanged to the halo
      (the fluxes of owned elements divide by it at halo vertices); a row without wet elements
      (all cavity) keeps the node area (FESOM2 behaviour, keeps the divides finite)
- [x] antidiffusive coefficients consistent − lumped per element (triangle (−2,1,1)·S/12, quad
      (−5,2,1,2)·S/36); `fct_fluxes(elem, 4)`; limiter loops over `nv`
      — `ice_acoef`; flux computation, ±sums, limiting and the update loop over
      `elem2D_nnodes` (a 4-slot limiter on triangles would take the min at vertex 1 twice:
      mutation check fails F5 by 0.15–0.25)
- [x] TG rhs as in Technical Details (advection via the element mass contributions, diffusion
      with 4-slot gradients, TG term with the kite-weighted mean velocity)
      — `dt·(dx(n)·(M_e u)_q + dy(n)·(M_e v)_q) − dt·S·diff·∇n·∇q − ½dt²·S·(ū·∇n)(ū·∇q)` with
      `ū = Σ elem_w·u`, over 4 slots (a triangle's slot 4 adds exactly 0)
- [x] tests: row sums = lumped mass, total = total area; a constant ice field stays constant
      (low and high order); Σ lumped·m conserved on the mixed patch; triangles equal today's FCT to
      roundoff; missing column triggers the error (error case)
      — `test_ice_mass` (NEMO-mixed cut np1/2 + analytic triangle, mixed and parallelogram
      patches at np1, NEMO-quad cut np2, pi_w np1 and core2_w np8 rotated): F0 tables ==
      cyclic-distance coefficients, row sums 1/3, 1/4, acoef = mcoef − diag(rowsum) (≤ 2.8e-17);
      F1 `fct_massmatrix` == an independent assembly (0); F2 row sums == lumped == Σ S/nv on owned
      and halo nodes (≤ 5.3e-16), Σ lumped == total wet area (≤ 1.4e-14), triangles: lumped ==
      node area (≤ 6.7e-16); F3 constant field at rest with ice_diff = 10: TG rhs ≤ 2e-18, low
      order, high order and final == c (≤ 3.4e-16), uniform drift on the patches: nodes ≥ 3 rings
      from the boundary stay c (≤ 3.2e-16); F4 Σ lumped·f conserved for a/m/m_snow with a
      non-uniform velocity and limiter-active fields (≤ 2.3e-15, fields change by up to 0.25);
      F5 triangles: TG rhs, high order, low order, final == the former FCT (FESOM2 formulas, in
      the test) ≤ 7.4e-16 (pi_w, core2_w np8, triangle patch); ctest
      `test_ice_mass_missing_col_np1` (a CSR column removed → the message, regex pair).
      `fesom_icefctdump` CORE2 1 rank new vs the pre-task binary: rhs ≤ 6.2e-16, a/m/m_snow ≤
      1.8e-15 (EVP), ≤ 2.4e-15 (mEVP)
- [x] `run_tri_regression.sh` — roundoff
      — 7/7: pi PRE 6.3e-15, EARLY 9.2e-13, OCEAN 1.1e-7, ICE 5.9e-7, MIX 4.4e-9; CORE2 PRE
      2.5e-11, EARLY 2.2e-8, OCEAN 3.1e-5, ICE 9.7e-5, MIX 7.8e-4 (the Task-15 levels); decomp
      OCEAN ≤ 8.8e-10; conservation ≤ 2.8e-14
- [x] run ctest — must pass before Task 17 — scratch ctest 136/136 (wiped + reconfigured)

### Task 17: Conservation and 1-vs-N runs on the cut-out meshes

**Files:**
- Modify: `test/CMakeLists.txt` (NEMO and mixed cut-out variants of `fesom_conserve`)
- Create: `tools/run_cutout_short.sh`
- ➕ Modify: `src/drivers/fesom_conserve.F90` (sea-ice option `FESOM3_CONSERVE_ICE`)
- ➕ Create: `tools/cutout_diag.py` (1-vs-N checkpoint comparison, checkerboard measures, self-test)

- [x] conservation tests (heat/salt 1e-12) on the NEMO and mixed cut-outs, np 1/2, with ice on
      — `fesom_conserve` had no sea ice: new presence knob `FESOM3_CONSERVE_ICE` runs, every step
      as in the production runloop, ocean2ice → EVP/mEVP (`FESOM3_WHICHEVP`) → ice TG rhs → ice FCT
      → `oce_fluxes_mom` (ice-ocean drag into `stress_surf`/`stress_node_surf`); the cut-outs are
      tropical (a cold start has no ice), so the ice is SEEDED analytically (a_ice
      0.45·(1 + sin(k·lon)·cos(k·lat)), 5° wavelength, ice-free troughs so the limiter acts; m_ice,
      m_snow ∝ a_ice) and driven by a uniform wind-on-ice stress (0.1, −0.05) N/m²; no cut_off /
      thermodynamics (they exchange heat/salt/water with the ocean — a budget, not a conserved
      quantity). Gated: heat, salt, volume and the ice masses Σ fct_lumpedmass·f (a, m, m_snow)
      within 1e-12. 6 ctests `fesom_conserve_{q,mix}cut_ice_np{1,2}` + `_ice_mevp_np2`, zstar with
      the production ocean physics the driver supports (TKE, GM + Redi, splines, vector invariant,
      rvo_upwind 1, w split), 20 steps: |drift| ≤ 1.2e-14 (ocean) and ≤ 2.8e-15 (ice), max|uice|
      0.35–0.39 m/s. Mutation (ice low-order step with `mesh%area` instead of the lumped mass): ice
      drift 1.3e-6–7e-6 → CONSERVATION VIOLATED
- [x] short runs (a few days) on the cut-outs at 1 vs 2/8 ranks: stable, conserving, N vs 1
      within roundoff; checkerboard energy of SSH and of a tracer measured and reported
      — `tools/run_cutout_short.sh` (6/6 PASS, ~25 min): 5 days (240 steps) on both cut-outs, np 1
      vs 2 and 8, each restarted from ONE np-1 0-step checkpoint (L60), SSH soltol 1e-14.
      **prod** (= the run_tri_regression prod physics: linfs, JRA55 1958, EVP ice, GM + Redi, TKE,
      splines, vinv, rvo_upwind 1, SW penetration, w split): stable on all np (max|eta| 0.27 m at
      the cold-start adjustment, max|uv| 0.67 m/s quad / 0.65 mixed; 10-day run 0.75 m/s), ice-free
      (tropics, a_ice = 0 — the ice is exercised in the cons case). 1 vs N per field every 6 h
      from the restart checkpoints (double, `cutout_diag.py compare`): first day ≤ 1.9e-10 (all
      fields, both meshes, np 2 and 8); 5 days: quad cut-out v 1.0e-4, tke 7.1e-5, T 1.5e-6, S
      6.4e-7, eta 1.8e-6; mixed cut-out v 3.8e-9, tke 2.2e-9, T 6.5e-11, S 2.7e-11, eta 4.0e-11. The growth
      on the quad cut-out starts at day ~3.9 near (−7°E, −1°N), where the equatorial flow runs
      against the artificial closed southern wall of the cut (np 2 vs np 8 grows the same way:
      1e-10 → 9e-6), amplified by TKE's thresholds: **prodpp** (same with PP instead of TKE, the
      clean diagnostic like mEVP for the EVP) stays ≤ 1.1e-7 over 5 days on both meshes. A
      second TKE effect: its thresholds turn roundoff into a local O(1e-9) jump within hours,
      and whether it fires depends on the arithmetic (np, build, dump probes on/off): one
      exploratory np-1 run gave 4.4e-9 at step 12 (6e-12 without TKE; same in a -O0 -check
      bounds,uninit build, no out-of-bounds or uninitialised read) and 7.9e-3 after 5 days. Gates:
      EARLY (first day, all fields) 1e-6, whole run 5e-2 (OCEAN/MIX), margins ≥ 200×
      resp. 6× over the worst exploratory run. **cons** (fesom_conserve with ice, 5 days, np
      1/2/8): heat/salt/volume |drift| ≤ 1.3e-14, ice masses ≤ 8.2e-15, max|uice| 1.12 m/s.
      Checkerboard (`cutout_diag.py checker`, ssh/sst/sss snapshots every step of the np-1 prod
      run; E_glob = projection on the global 2-colouring of the quad-side graph = the exact SSH
      null mode, E_loc = per-quad (+1,−1,+1,−1) energy, both normalised by the variance over the
      quad vertices; R_tri/R_quad/R_iface = rms of the linear-exact neighbour residual / std at
      triangle-only, quad-only and interface nodes), daily means of the 10-day run: **SSH** E_glob
      no growth (quad cut-out 1e-16…4e-15, mixed 1e-7…3e-6 fluctuating), E_loc rises during the
      spin-up (day 1→3: 8e-5 → 1.8e-4 quad, 1.1e-4 → 2.4e-4 mixed) and then saturates/declines
      (day 10: 1.3e-4 / 2.8e-4). **Tracers** grow over the whole 10 days: SST E_loc 7.3e-6 →
      4.3e-5 (quad) and 9.2e-6 → 6.1e-5 (mixed), SSS E_loc 9.8e-6 → 1.2e-4 / 3.0e-5 → 3.3e-4
      (E_glob stays ≤ 7e-7); on the mixed cut-out the grid-scale residual grows at quad and
      interface nodes more than at triangle nodes (SSS day 1 → 10: R_quad 9.9e-3 → 4.4e-2,
      R_iface 9.5e-3 → 5.7e-2, R_tri 1.2e-2 → 2.0e-2; SST R_iface 1.0e-2 → 2.7e-2, R_tri 9.0e-3 →
      1.5e-2) — consistent with the Task-11 finding (linear tracer not advected exactly at
      interface nodes and irregular quads). Time series: `cutout_short{,_10d}/prod_*/checker.csv`
- ⚠️ the tracer checkerboard / grid-scale energy GROWS (SST ×6, SSS ×11–12 in E_loc over 10
      days, faster at quad and triangle–quad interface nodes than at triangle nodes); the SSH
      null mode does not grow. Not an exponential instability within 10 days, but no saturation
      either — a remedy (e.g. hourglass control or the dual-segment face values of Task 11) is
      the user's design decision; nothing implemented
- [x] ⚠️ if the checkerboard energy grows, stop and ask the user before continuing
      (user ordered uninterrupted execution — reported, not asked)
- [x] run ctest + short runs — must pass before Task 18 — scratch ctest 142/142 (wiped +
      reconfigured), `run_cutout_short.sh` 6/6, `run_tri_regression.sh conserve:pi:1
      conserve:pi:2` PASS (the only cases using the changed driver; ice off = unchanged path)

### Task 18: CORE2 throughput, setup workflow and example job for the NEMO mesh

**Files:**
- Create: `work/job_levante_nemo025` (new file; `work/job_levante` untouched)
- Modify: `docs/HANDOFF.md`

- [x] CORE2 throughput with 4-slot loops vs ~14 min per model year, reported (user's run or a
      short timing run on the cut-out scaled — no production submission by us)
      — no user run of the hybrid code yet → login-node timing, old `frozen_task6` (5b2a7f3, pre
      4-slot) vs new scratch build (ec099ea), original CORE2, production driver with the job_levante
      physics, DEFAULT SSH tolerance, no output/restart, loop timers over steps 21..N (setup and
      first steps excluded), interleaved. np 16 (4 + 4 runs, 200 steps): LOOP TOTAL old min 456.5 /
      median 458.0 ms/step, new 492.8 / 499.1 → **new/old 1.079 (min) / 1.090 (median)**; np 8 (3 + 3,
      120 steps): 850.9 / 853.9 vs 926.9 / 931.9 → **1.089 / 1.091**; run-to-run noise 0.8–3.3 %.
      Only component above 20 %: sea ice ×1.29 (88 → 114 ms/step at np 16); `perf` (all ranks): EVP
      subcycles 1.33× on the most ice-loaded rank, the others wait in MPI; ocean ×1.01–1.09.
      **Scaled: 14 min × 1.08–1.09 ≈ 15.1–15.3 min per model year on 512 ranks** (indicative; np 8–16
      on a shared login node). NEMO 0.25° ESTIMATE (mesh counts 4.3–8.4× the CORE2 work per step,
      900 s step): 2.2–4.3 h per model year on 512 ranks, 1.1–2.1 h on 1024. Not optimised. Details
      and scripts: `docs/HANDOFF.md` "CORE2 throughput of the 4-slot code",
      `/work/ab0995/a270029/fesom3_hybrid/t18_timing/`
- [x] document the workflow: `nemo_meshq.py` → `remove_deadends.py` → (optional `split_quads.py`)
      → `fesom_ini` → model
      — `docs/HANDOFF.md` "Mixed-mesh setup workflow": exact commands, inputs/outputs, fesom_ini
      (build location, namelist.config, files written, keeps old edge/level files, exit status),
      checkers (check_mesh, check_dist, check_meshdiag, make_dist1), model env, runtime/memory on
      ORCA025. Every command run end to end in a scratch dir (`t18_workflow/run_workflow.sh`):
      nemo_meshq 5 s/0.7 GB, remove_deadends 12 s, split 5 s, cut 3 s, check_mesh 4 s/0.9 GB,
      fesom_ini 1024 parts 83 s/0.3 GB (max 16 neighbour ranks, limit 32), check_dist 1024 PASS 73 s;
      results identical to the persistent test meshes except the known fold point (nod2d.out, and one
      split diagonal in nemo025mix); model smoke run on the regenerated mixed cut-out (np 4, 8 steps)
- [x] example job script for NEMO 0.25° next to `job_levante`; dry-run its setup part (no submission)
      — `work/job_levante_nemo025` (untracked like `work/`): job_levante physics copied; mesh
      nemo025q_fesom (nemo025mix in a comment), 8 nodes × 128 = 1024 ranks, 1 h; setup checks mesh,
      depth file, binary, NPART == SLURM_NTASKS, edge/level files + complete dist_NPART, runs
      `srun -n 1 -N 1 fesom_ini` if missing; `DRYRUN=1` prints, `DRYRUN_INI=1` also runs fesom_ini.
      USER: choose — STEP_PER_DAY 96, run length 1 month, nodes/wall time, RUNID/output tree.
      Dry run as is: `setup: ERROR: missing …/build_intel_dp/bin/fesom_ini (rebuild build_intel_dp
      from bottom-at-vertices)` + would-run fesom_ini line for dist_1024 (exit 1, correct: the user's
      binary predates the branch); with the scratch build + a mesh with dist_1024: `ok dist_1024:
      rpart.out, format.out (version 1), 1024 my_list + com_info` … `would run: srun -l -n 1024
      …/fesom_lifecycle_native_mr` `DRY RUN done`; `NPART=8 DRYRUN_INI=1` on a scratch nemo025mix:
      fesom_ini 21 s, check_dist PASS (logs `t18_workflow/t18_dryrun_*.log`)

### Task 19: Verify acceptance criteria
All runs below use ONE build: a clean clone of the committed HEAD 98b5164
(`git clone` → `/work/ab0995/a270029/fesom3_hybrid/t19_clone`, `./configure.sh --compiler intel
--precision dp --clean --build`), logs `/work/ab0995/a270029/fesom3_hybrid/t19_*.log`.
- [x] mixed meshes read, partitioned and run on the cut-outs (NEMO, mixed) and full CORE2/pi
      — **partitioned**: `run_ini_mixed.sh` (clone `fesom_ini` + `check_dist.py` PASS on
      nemo025q_cut / nemo025mix_cut dist_1/2/4/8, nemo025q_fesom / nemo025mix dist_1/2/8/128, pi_w
      1/2/4/8, core2_w 1/2/4/8/128, old pi dist_2/8 + CORE2 dist_8; all error cases detected) and
      `run_ini_compare.sh` (pi/CORE2 == FESOM2's build) both PASS; the 1 036 files of the six test
      meshes (dist_N, format.out, edges, levels) regenerated BYTE-IDENTICAL (md5 before/after,
      `t19_meshes_md5_{before,after}.txt`). **Read**: ctests `test_mesh_read_*`/`test_partit_*` on
      both cut-outs and pi_w np 1/2, original pi dist_2, malformed case. **Run** (np 1 and np > 1):
      cut-outs — `run_cutout_short.sh` 6/6 (prod/prodpp/cons, 5 days, np 1 vs 2/8: stable, max|eta|
      0.27 m, max|uv| ≤ 0.70 m/s; 1 vs N first day ≤ 1.9e-10, whole run quad prod 1.0e-4 (TKE, as in
      Task 17), mixed prod 3.8e-9, prodpp ≤ 1.1e-7; cons drift ocean ≤ 1.3e-14, ice ≤ 8.2e-15) + 22
      `fesom_conserve_{mixcut,qcut}_*` ctests np 1/2; pi — `run_tri_regression.sh` prod:pi:1/2,
      decomp:pi_w:2/8 (np 1 vs np N from one checkpoint ≤ 8.8e-10), conserve:pi:1/2; CORE2 —
      prod:core2:8 (original mesh) and the extra case decomp:core2_w:8 (np 1 vs np 8 on the
      fesom_ini CORE2 from one checkpoint, 4 steps: PRE 0, EARLY 2.6e-12, OCEAN 2.7e-5, ICE 7.3e-5,
      MIX 2.2e-5 — PASS at the CORE2 tolerances; `t19_tri_regression_core2w.log`). **Full NEMO
      0.25°** (optional item): production physics (job_levante set, linfs, JRA55, EVP, TKE, GM/Redi,
      900 s step), np 8 on the existing dist_8 (no dist_16 written into the shared meshes), cold start,
      192 steps = 2 days each (`t19_full/run_full.sh`, logs `t19_full/*/log`): nemo025q_fesom 13.4 min,
      nemo025mix 14.5 min wall (shared with the other runs), finite, exit 0; max|eta| 2.55 / 2.52 m
      (step 106), max|uv| 2.42 / 2.60 m/s (steps 38 / 34, cold-start adjustment), then bounded
      (day 2: 1.8 m, 2.1–2.2 m/s) — global cold-start values, not judged further (full-size runs are
      the user's)
- [x] all decisions of the Overview/Technical Details implemented; grep checks clean
      — decision → code (all in `src/` unless noted): 4-slot storage, triangle 4th = 1st,
      3/4-column detection `mod_mesh_read.F90:122` (`read_elem_line`), quad turn 2↔4 `:228`;
      depth `aux3d.out`/`depth.out` `:169`, node levels only from `nlvls.out` `:100,:478`, element
      levels min/max over vertices `:314-315`; fesom_ini ported `mesh_part/`, `CMakeLists.txt:136`,
      quad orientation `mesh_part/src/fvom_init.F90:303-360`, `use_depthfile='nemo'` `:726-748`,
      ≥ 2-neighbour stop kept (run_ini_mixed dead-end case); `format.out` + appended sections
      `mesh_part/src/oce_local.F90:14-283` (`dist_ext`, unique edge owner), read
      `infra/mod_partitioning.F90:157-190`, empty lists in `synthesize_1rank` `:86` and
      `tools/make_dist1.py:22`; exchanges `infra/mod_halo.F90:160-203` (`exchange_edge`,
      `exchange_nod_full`, `require_dist_format` stop), `owned_bounds(..., nNodX)`
      `infra/mod_part_bounds.F90:53-68`; edge order interior-first checked by
      `tools/mesh/check_dist.py:117`; `edge_tri` on-disk kept, missing 2nd → 1st at read
      `mod_mesh_read.F90:93,:535`; geometry `mesh/mod_mesh_areas.F90`: centroid/area `:181-263`,
      kite `elem_w` `:267-299`, `nod_in_elem2D_w` `:301-324`, node area `:512`, Green–Gauss
      `gradient_sca(8)` `:326-356`, `elem_edges`/`elem_neighbors` (missing = self) `:358-418`,
      `edge_cross_dxdy(3:4)=0` `:467`, resolution `:556`; `gradient_vec` deferred
      (`types/mod_mesh.F90:79`, not computed); kernels — element means `elem_w` and node gathers
      `nod_in_elem2D_w` in `oce_ale.F90:103`, `oce_dyn_velrhs.F90:158-328`, `oce_pgf.F90:105`,
      `oce_dyn_vinv.F90:195` (rvo_upwind centroid normals `:256-271`), `oce_dyn_ivertvisc.F90:75`,
      `oce_ale_tracer.F90:465`, `oce_pressure_bv.F90:347,376,578`, `oce_muscl_adv.F90:57`
      (`find_up_downwind_elements`), `oce_fer_gm.F90:293`, KPP/TKE/PP, `oce_mo_conv`,
      `mod_forcing_bulk.F90:10`, `mod_ice_oce_coupling.F90:27`; unique-vertex loops over
      `elem2D_nnodes` only where a repeat matters (`oce_initial_state.F90:280`,
      `oce_ssh_rhs.F90:162`, `mod_ice_fct.F90:142-552`, `mod_ice_setup.F90:202`); SSH sparsity
      `oce_ssh_rhs.F90:87-315`; ice dynamics `mod_ice_dyn.F90:70-203` (ocean2ice `:115`); ice FCT
      `mod_ice_fct.F90:51-155` (`ice_mcoef`/`ice_acoef`, TG `ubar` `:154`), lumped mass, `fct_fluxes(:,MAX_NV)`,
      missing-column stop `mod_ice_setup.F90:145-246`; output `io/mod_io_coords.F90:33` (centroids),
      `io/mod_io_meshdiag.F90:19-115` (`face_nodes` fill −1); checkerboard asserted
      `test/test_ssh_stiff.F90:17,66,82`, measured `tools/cutout_diag.py`; Python tools
      `tools/mesh/{nemo_meshq,remove_deadends,split_quads,cut_region,check_mesh}.py`. No decision
      without code. **Greps** (src/): `elem2D_nodes(1:3` → only `mod_geom_dump` (FESOM2 oracle,
      triangle-only guard) and `mod_mesh_analytic:122-138` (generator, slot 4 set at :143);
      `gradient_sca(1:3|4:6|5:7` → only `mod_geom_dump:65`; `do …=1,3` → none over vertices (ice
      tracers, smoke/restart drivers, 3 smoothing passes `mod_mesh_areas:576`); `elnodes(3)` → none;
      `enodes(3)` → `oce_adv_tra_flux.F90:41` (edge nodes, only 1:2 used — harmless FESOM2 leftover);
      `/3` element means / hidden 1/3 → only documented exceptions (sub-triangle centroid
      `/(3(s1+s2))` `mod_mesh_areas:209-210`, `oce_dyn_vinv:270-271`; spline `oce_adv_tra_ver:134-135`;
      KPP cube roots), the rest are comments; zero-index boundary tests (`el(2)>0`, `el(i)<1`,
      `nb<=0`, `edge_tri(2,·)>0`) → only the convention builders (`mod_mesh_read:93`,
      `mod_mesh_analytic:240-243`) and the `-1` invariant (`mod_mesh_areas:409`), none in test/;
      triangle/quad branches (`nnodes==3`, `nv==3`) → none in kernels: only the geometry setup
      (exact triangle kite 1/3 `mod_mesh_areas:282`, adjacency slot 4 := 1 `:413`), the CW turn
      (`mod_mesh_read:253`, swap 2↔3 vs 2↔4), the analytic generator and fesom_ini (literal
      `/3.0`,`/4.0`, L7); no TODO/FIXME added by this plan (the 5 in `mesh_part/` are FESOM2's,
      the 2 in `src/` pre-date it); no debug prints/stops added (all added `write`s are error
      messages, the SSH-knob echo or the conserve driver's report). Fixed: stale comment in
      `mod_mesh_read.F90` `assert_bottom_invariant` (cited the old `1/(3*vol)` and old line numbers;
      comment only, no test)
- [x] full suite: `./configure.sh` in a scratch build + `ctest` (with `ulimit -s unlimited`)
      — clean clone of HEAD 98b5164: `./configure.sh --compiler intel --precision dp --clean
      --build` built everything from the committed tree (65 targets, 62 binaries incl. `fesom_ini`,
      no missing file), `ctest -j 8` **142/142** in 49 s, none skipped (`t19_clone_build.log`,
      `t19_ctest.log`); scratch build with the comment fix: 142/142 (`t19_scratch_ctest.log`)
- [x] tool tests (`tools/mesh/test_*.sh`), `run_ini_compare.sh`, `run_tri_regression.sh` pass
      — `test_nemo_meshq.sh` PASS (1 fold line, as Task 1), `test_mesh_tools.sh` PASS,
      `run_ini_compare.sh` PASS, `run_ini_mixed.sh` PASS, `run_tri_regression.sh` (all 7 default
      cases) 7/7 — pi PRE 6.3e-15 EARLY 9.2e-13 OCEAN 1.1e-7 ICE 5.9e-7 MIX 4.4e-9; CORE2 PRE
      2.5e-11 EARLY 2.2e-8 OCEAN 3.1e-5 ICE 9.7e-5 MIX 7.8e-4; decomp OCEAN ≤ 8.8e-10; conserve drift
      ≤ 2.8e-14 (= the Task-16 levels) — plus decomp:core2_w:8 1/1, `run_cutout_short.sh` 6/6
      (logs `t19_{nemo_meshq,mesh_tools,ini_compare,ini_mixed,tri_regression,
      tri_regression_core2w,cutout_short}.log`)

### Task 20: [Final] Update documentation
- [x] `docs/HANDOFF.md`: mixed-mesh section (formats, arrays, conventions, setup workflow,
      exchanges, retired gates, checkerboard property, ice FCT lumped mass)
      — one section "Hybrid (mixed triangle + quad) meshes": new subsections File formats, Mesh
      arrays and conventions, Geometry, Kernels, Exchanges and halo rings, SSH stiffness matrix and
      the checkerboard null mode, Sea ice, Tests and tools; the Task-18 subsections (setup workflow,
      CORE2 throughput, NEMO job) kept in it; "Gate migration for hybrid meshes" MOVED into it from
      "Gate recipes (archived)" (pointer left there; CORE2 EARLY 1e-7 and the end-of-plan levels
      added); new "Known findings / open questions for the user" (15 items). "Where we are" status
      line; the SSH-solver line (Jacobi preconditioner `28af61c`, `FESOM3_SSH_SOLTOL`/`MAXITER`) under
      the M2.12c-2 milestone
- [x] `docs/LESSONS.md`: lessons learned — L62–L69 (MATLAB precision port, 4-slot slot-use
      classification, frozen reference + amplifier control runs, checkerboard null mode, exact
      ingredients vs assembled tendency, FCT lumped mass, wrapped element embedding, verify setup
      outputs against the reader)
- [x] move this plan to `docs/plans/completed/` (moved by the orchestrator after the review phases)

## Post-Completion
*Items requiring manual intervention or external systems - no checkboxes, informational only*

**Manual verification (user)**
- full-size NEMO 0.25° and mixed-mesh short runs (job script from Task 18), then long runs;
  climate-level comparison with CORE2 and with NEMO reference fields; watch the checkerboard
  diagnostics on the full mesh
- visual check of the mixed mesh and of the tripolar fold region

**External**
- pull request only when the user says so, after testing
