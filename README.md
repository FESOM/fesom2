# FESOM_F

FESOM_F is a clean-architecture **Fortran** reimplementation of the FESOM2 ocean + sea-ice model.
It began as a line-for-line transcription of FESOM2's numerical kernels and reproduced a reference
FESOM2 configuration **bit-for-bit** up to milestone M10; the surrounding "plumbing" (data types, the
parallel layer, configuration, I/O) is redesigned for clarity and future flexibility. The guiding rule
was *"first reproduce exactly, then optimize."* Since v1 the model has moved beyond FESOM2 on purpose
(see below).

**Main achievements beyond FESOM2:**

- **Bottom at vertices.** The vertex column defines the bottom and each element takes the depth of
  its shallowest vertex, so no velocity touches topography; tracer cells are straight prisms with the
  same area at every level.
- **Vector-invariant momentum advection** (`FESOM3_MOMADV_OPT=1`) with an upwind-blended face
  vorticity, together with a Courant-dependent explicit/implicit split of the vertical velocity
  (`FESOM3_WSPLIT`).
- **Quad and hybrid (triangle + quad) meshes.** One set of 4-slot kernels runs triangle, quad and
  mixed meshes, from mesh generation and partitioning (`fesom_ini`) to output. FESOM_F runs on the NEMO
  ORCA025 quad mesh, and triangle meshes give the same results as before at roundoff level.

**Status:** v1 is feature-complete. Milestones M0–M10 are done and tagged (`m0`, `m1`, `m2-mvp`,
`m3` … `m10`): dynamical core, sea-ice EVP/mEVP, GM/Redi, KPP, zstar ALE, TKE, multi-year production
runs, Zarr output, and restart/checkpoint. Through M10 the model was byte-exact against FESOM2.

### Changes since v1 (September–October 2026)

The results now differ from FESOM2 on purpose. Bit-identity with the FESOM2 oracle is retired;
conservation, mesh-delta, regression and self-consistency gates take its place
([TESTING.md](TESTING.md) §6a).

- **Bottom at vertices.**
  - `nlvls.out` is authoritative: the vertex column defines the bottom, and each element takes the
    depth of its shallowest vertex (3 or 4), so velocities touching topography vanish by construction.
  - The scalar-cell area is 1-D and the same at every level, and the z* stretch covers the full column.
  - New conservation and no-leakage gates (`fesom_conserve`), including the Redi path.
  - See [IMPLEMENTATION.md](IMPLEMENTATION.md) §Bottom at vertices.
- **Time step:** `FESOM3_STEP_PER_DAY` sets the number of steps per day.
- **Forcing:**
  - The time bracket is reset at a year rollover, which fixes the leap-year crash.
  - Three fixes synced from FESOM2 2.7.4–2.8.0.
- **Vertical mixing:**
  - Cubic-spline N² and shear at the interfaces (`FESOM3_N2_SPLINES`, `FESOM3_SHEAR_SPLINES`), and a
    ramped convective adjustment.
  - `N2smth_v` is removed.
- **Momentum advection:** a vector-invariant form (`FESOM3_MOMADV_OPT=1`). The face vorticity is
  upwind-blended (`FESOM3_RVO_UPWIND`), and implicit w_i advection is in advective form.
- **w split (`FESOM3_WSPLIT`):** a smooth, Courant-dependent explicit/implicit split of the vertical
  velocity, with implicit FCT tracer advection.
- **MUSCL:** up/downwind gradients are computed on the fly (reference RK3 rule), and the stored
  per-edge gradients are gone. The CORE2 MFCT step is about 9 % faster.
- **SSH solver:**
  - Jacobi preconditioner and fewer halo exchanges.
  - Tolerance and iteration limit can be set with `FESOM3_SSH_SOLTOL` and `FESOM3_SSH_MAXITER`.
- **Hybrid meshes (triangles and quads):**
  - Python mesh tools (`tools/mesh/`): NEMO quad mesh, checks, splitting, cut-outs.
  - `fesom_ini`, a port of FESOM2's mesh_part with METIS 5.1.0 that handles mixed meshes.
  - Unified 4-slot kernels with kite weights for ocean dynamics, tracers, mixing and GM, forcing and
    ice coupling, EVP/mEVP and the sea-ice FCT.
  - Triangle-only meshes give the same results at roundoff level. See §5.
- **Shared-element velocity sync:**
  - A border element belongs to every rank that owns one of its vertices. Each step, `update_vel` now
    sends one designated copy to the other owners (`FESOM3_SHARED_SYNC=0` turns this off).
  - This fixes the volume leak of the earlier multi-rank runs.
  - A run resumed from a restart now matches a straight run bit for bit at np > 1.
- **Zarr I/O:** the chunk writer and reader no longer copy through `transfer()`. This fixes the
  end-of-run restart crash on large meshes.

Production runs on Levante with these changes conserve volume to roundoff: CORE2 on 512 ranks (about
95 model years per day) and the NEMO ORCA025 quad mesh on 1280 ranks (about 29 model years per day).

See **[IMPLEMENTATION.md](IMPLEMENTATION.md)** for the architecture and design decisions, and
**[TESTING.md](TESTING.md)** for how it is validated against FESOM2 and how it scales.

> **A note on the name.** The project is called **FESOM_F** in this documentation, but it lives in a
> repository still named `fesom3`, and you will see that working name in real identifiers: the
> environment-variable prefix `FESOM3_`, build directories like `build_intel_dp`, driver binaries
> `fesom_*`, and git tags `m0…m10`. Those are the actual names — use them verbatim. Only the prose
> uses "FESOM_F".

---

## 1. Requirements

FESOM_F currently builds and runs on the **Levante** supercomputer (DKRZ). Everything is provided
through environment modules, loaded automatically by `env.sh`:

| | Anchor (bit-identity) | Portability |
|---|---|---|
| Compiler | `intel-oneapi-compilers/2022.0.1` | `gcc/11.2.0` |
| MPI | `openmpi/4.1.2-intel-2021.5.0` | `openmpi/4.1.2-gcc-11.2.0` |
| netCDF | `netcdf-c/4.8.1` + `netcdf-fortran/4.5.3` | same |

- **netCDF-Fortran is required** (forcing and initial-condition reads); the build aborts if
  `nf-config` is not on the `PATH`.
- **liblz4 is optional** — it enables the LZ4 output compressor; without it the writer is
  uncompressed-only.
- **A C compiler is required** for the mesh partitioner `fesom_ini` (`mesh_part/`, with METIS 5.1.0
  vendored in `lib/metis-5.1.0/`); `env.sh` sets `CC=mpicc CXX=mpicxx` and `configure.sh` passes them
  to CMake. Builds that bypass `configure.sh` must configure a C compiler themselves.
- **The Python mesh tools** (`tools/mesh/*.py`, `tools/make_dist1.py`) need numpy/scipy:
  `module load python3` (the system `python3` has neither).
- The **Intel + double-precision** build is the *anchor*: its compiler flags are copied verbatim
  from FESOM2 v2.7.3, which is what makes bit-identical results possible. The GNU build is for
  portability testing.

---

## 2. Build

`configure.sh` sources `env.sh` for you (loading the modules and the Levante MPI workaround), runs
CMake, and builds:

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor build
```

Binaries land in `build_intel_dp/bin/` — the model drivers and the mesh partitioner `fesom_ini`
(§5) — and Fortran modules in `build_intel_dp/module/`. A binary built before the hybrid-mesh work
cannot read quad meshes and has no `fesom_ini`: rebuild before using a NEMO/quad mesh. The portability
build is the same with `--compiler gnu` (→ `build_gnu_dp/`):

```bash
./configure.sh --compiler gnu --precision dp --clean --build
```

`configure.sh` options: `--compiler intel|gnu`, `--precision dp|sp|hp`, `--debug` (builds into a
separate `build_<compiler>_<precision>_debug/` so it never overwrites the anchor), `--clean`,
`--build`, `--jobs N`. The build directory name encodes the compiler and precision, so several builds
coexist.

> Build and run on a Levante node — `env.sh` only knows how to set up Levante and will stop on any
> other host.

---

## 3. Quick check

Run the self-test suite on a login node:

```bash
cd build_intel_dp && ctest --output-on-failure
```

For a stronger check, run the conservation gate — an unforced 20-step pi run that verifies total
heat and salt content hold to round-off under `zstar`, that no velocity survives in a partly-land
prism, and that the layer thicknesses stay consistent:

```bash
bash tools/run_conserve_pi.sh      # expect: GATE OK
bash tools/run_bottom_audit.sh     # mesh-delta audit on pi + core2
```

The FESOM2 oracle byte-gates (`run_step_gate.sh` and friends) are retired — see
[TESTING.md](TESTING.md) §6.

What these checks mean — and the full validation story — is in [TESTING.md](TESTING.md).

---

## 4. Running a simulation

The production driver is **`fesom_lifecycle_native_mr`** (multi-rank, computes its own air–sea
fluxes from atmospheric forcing). A run is configured through `FESOM3_*` environment variables; the
example below mirrors `tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`:

```bash
source env.sh intel                         # modules + Levante MPI workaround

POOL=/pool/data/AWICM/FESOM2
export FESOM3_MESH_DIR="$POOL/MESHES_FESOM2.1/core2"
export FESOM3_IC_FILE="$POOL/INITIAL/phc3.0/phc3.0_winter.nc"
export FESOM3_FORCING_DIR="$POOL/FORCING/JRA55-do-v1.4.0"
export FESOM3_FORCING=JRA55
export FESOM3_RUNOFF_FILE="$POOL/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc"
export FESOM3_SSS_FILE="$POOL/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc"
export FESOM3_START_CLOCK="0 1 1958"        # "timeofday day year"
export FESOM3_NSTEPS=35040                   # 2 model years at 48 steps/day

# physics (the default FESOM2 production configuration)
export FESOM3_WHICH_ALE=zstar
export FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SW_PENE=1 FESOM3_CHL_SWEENEY=1

srun -l -n 864 ./build_intel_dp/bin/fesom_lifecycle_native_mr > run.log 2>&1
```

The physics options added since v1 are off by default (the shared-element sync is on). The current
CORE2 and NEMO production runs use `FESOM3_WHICH_ALE=linfs` and add:

```bash
export FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0     # vector-invariant momentum advection
export FESOM3_WSPLIT=1                               # Courant-dependent explicit/implicit w split
export FESOM3_SHEAR_SPLINES=1 FESOM3_N2_SPLINES=1    # spline shear and N² (set both)
export FESOM3_STEP_PER_DAY=48                        # 1800 s on CORE2; 72 (1200 s) on NEMO 0.25°
```

> **On/off switches are presence-based.** `FESOM3_FER_GM`, `FESOM3_REDI`, `FESOM3_MIX_TKE`,
> `FESOM3_MIX_KPP`, `FESOM3_SW_PENE`, `FESOM3_CHL_SWEENEY`, `FESOM3_WSPLIT`, `FESOM3_SHEAR_SPLINES`,
> `FESOM3_N2_SPLINES` and `FESOM3_STEP_DIAG` are on whenever the variable is set and not empty, also
> for `=0`; leave it unset to turn the option off. `FESOM3_MOMADV_OPT`, `FESOM3_RVO_UPWIND`,
> `FESOM3_WSPLIT_MAXCFL`/`_MINCFL` and `FESOM3_SHARED_SYNC` take values.

**Where the knobs are defined.** The authoritative, always-current list of run variables is the
environment block at the top of **`src/drivers/fesom_lifecycle_native_mr.F90`** (each is read with
`get_environment_variable`); the SSH solver knobs and `FESOM3_SHARED_SYNC` in the table below are read
where they act (the solver, `src/infra/mod_shared_elem.F90`). The `tools/*.sbatch` and `tools/run_lifecycle_*.sh` scripts are working
examples. The most common variables:

| Variable | Meaning |
|---|---|
| `FESOM3_MESH_DIR` | mesh directory (must contain `dist_<NP>/`) |
| `FESOM3_IC_FILE` | initial T/S netCDF |
| `FESOM3_FORCING_DIR`, `FESOM3_FORCING` | atmospheric forcing directory and set (`JRA55`/`CORE2`) |
| `FESOM3_RUNOFF_FILE`, `FESOM3_SSS_FILE` | runoff and sea-surface-salinity restoring climatologies |
| `FESOM3_START_CLOCK` | cold-start clock, `"timeofday day year"` |
| `FESOM3_NSTEPS` | number of time steps (or use `FESOM3_RUN_LENGTH` + `FESOM3_RUN_UNIT`) |
| `FESOM3_STEP_PER_DAY` | time steps per day (default `48` = 1800 s; CORE2 uses 48, finer meshes need more) |
| `FESOM3_WHICHEVP` | sea-ice solver: `0` = EVP, `1` = mEVP |
| `FESOM3_WHICH_ALE` | vertical coordinate: `linfs` (default) or `zstar`. `zlevel` is not ported, and the driver does not reject it: do not use it |
| `FESOM3_FER_GM`, `FESOM3_REDI` | Gent–McWilliams bolus advection, Redi isopycnal diffusion |
| `FESOM3_MIX_KPP`, `FESOM3_MIX_TKE`, `FESOM3_SW_PENE` | vertical mixing and shortwave penetration |
| `FESOM3_SHEAR_SPLINES`, `FESOM3_N2_SPLINES` | cubic-spline shear and N² at the interfaces (set both) |
| `FESOM3_MOMADV_OPT` | momentum advection: `2` (default) = flux form, `1` = vector invariant |
| `FESOM3_RVO_UPWIND` | upwind blend of the face vorticity in [0, 1] (default `0`; vector invariant only) |
| `FESOM3_WSPLIT`, `FESOM3_WSPLIT_MAXCFL`, `FESOM3_WSPLIT_MINCFL` | Courant-dependent explicit/implicit w split; explicit CFL cap (default `1.0`) and onset (default half the cap) |
| `FESOM3_OUTPUT`, `FESOM3_OUTPUT_EVERY` | output directory (enables output) and step cadence |
| `FESOM3_RESTART`, `FESOM3_RESTART_LENGTH`, `FESOM3_RESTART_UNIT` | checkpoint directory and cadence |
| `FESOM3_SHARED_SYNC` | shared border elements: `1` (default) = every step the other owners take the designated owner's `uv` (no volume leak, exact np>1 restarts); `0` = off, the old unsynced copies (A/B comparisons only). Read in `src/infra/mod_shared_elem.F90` |
| `FESOM3_STEP_DIAG` | per-step global diagnostic line (max\|eta\|, max\|uv\|, ice, Tmax, Smax, and the shared-element numbers `divsum=`, `dUV=`, `dAB=`, `dsig=`); two extra collectives per step, off by default |
| `FESOM3_SSH_SOLTOL`, `FESOM3_SSH_MAXITER` | SSH solver tolerance and iteration cap — test/diagnostic only, read in `src/oce/oce_ssh_solve.F90` (defaults `1e-5`, `2000`) |

> **Configuration today vs. intent.** The production drivers are configured by these `FESOM3_*`
> variables plus in-code defaults. A write-once configuration *module* (`mod_config`, reading a
> single `namelist.config`) already exists, but the lifecycle drivers do not read it yet;
> consolidating run configuration into one namelist file is the intended future interface.

**Restart.** Set `FESOM3_RESTART=<dir>` to write checkpoints (cadence via `FESOM3_RESTART_LENGTH` /
`FESOM3_RESTART_UNIT`). On the next launch, if a complete checkpoint is found under that directory
(a `restart.latest` pointer to a folder containing `checkpoint.json`), the run resumes from it
byte-exactly; otherwise it cold-starts. See [IMPLEMENTATION.md](IMPLEMENTATION.md) §I/O.

---

## 5. Using different meshes

A mesh is a directory holding the global mesh files plus one sub-directory of partition files per MPI
rank count:

```
core2/
  nod2d.out  elem2d.out  aux3d.out  nlvls.out  edges.out  edge_tri.out  edgenum.out
  elvls.out                                       # present but NOT read (see below)
  dist_1/  dist_2/  dist_8/  dist_864/   ...      # one per rank count (partition files)
```

**Triangle, quad and mixed meshes.** `elem2d.out` may have 3 or 4 node ids per line (the column count
is taken from the first line; a triangle in a 4-column file repeats its 1st vertex). Meshes made from
NEMO grids carry a `depth.out` instead of `aux3d.out`; the model reads the levels from `depth.out`
when `aux3d.out` is absent.

> **`elvls.out` is no longer read.** FESOM_F places the bottom at **vertices**: `nlvls.out`
> (the per-vertex level count) is authoritative and each element's bottom is *derived* as the
> shallowest of its 3 or 4 vertices, so an element is wet only where all of its corners
> are. The file can stay in the mesh directory; nothing opens it. See
> [IMPLEMENTATION.md](IMPLEMENTATION.md) §Bottom at vertices.

Point `FESOM3_MESH_DIR` at this directory. The driver **auto-selects `dist_<NP>/`** from the number
of MPI ranks (`srun -n <NP>`); a single-rank run reads the global `nod2d.out`/`elem2d.out` directly.
So switching meshes — or switching rank counts on the same mesh — is just changing `FESOM3_MESH_DIR`
and `-n`, provided the matching `dist_<NP>/` exists.

These meshes are used routinely:

- **pi** — a small test mesh (~3k nodes), shipped with the FESOM2 oracle at
  `/home/a/a270088/port2/fesom2/tests/data/MESHES/pi` (with `dist_1/`, `dist_2/`, `dist_8/`). This is
  the default when `FESOM3_MESH_DIR` is unset, and what the byte-gates and `ctest` use — it runs in
  seconds on a login node.
- **CORE2** — the global ~1° triangle production mesh (126,858 nodes, 244,659 triangles, 47 levels)
  at `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2`, with partitions up to `dist_864`.
- **NEMO 0.25°** — the ORCA025 grid as a quad mesh, `nemo025q_fesom` (915,866 nodes, 902,408 quads,
  46 levels), and the mixed test mesh `nemo025mix` (685,452 quads and 433,912 triangles on the same
  nodes), at `/work/ab0995/a270029/meshes/`. Production runs use 1280 ranks and a 1200 s step.

To run on a new mesh, create its `dist_<NP>/` partition for your target rank count with the in-tree
partitioner **`build_*/bin/fesom_ini`** (FESOM2's `mesh_part`, ported): run it in a work directory
whose `namelist.config` sets `MeshPath`, `n_part=<NP>` and, for a NEMO mesh with `depth.out`,
`use_depthfile="nemo"`. It writes the edge and level files and `dist_<NP>/` into the mesh directory,
including `dist_<NP>/format.out`, which announces the appended halo lists needed by
`exchange_edge`/`exchange_nod_full`. Then point `FESOM3_MESH_DIR` at the mesh directory.

A mesh made from a NEMO grid goes through
`tools/mesh/nemo_meshq.py` → `remove_deadends.py` → [`split_quads.py` for a mixed mesh,
`cut_region.py` for a regional cut-out] → `check_mesh.py` → `fesom_ini` → `check_dist.py`. The exact
commands, the `namelist.config` and the pitfalls are in `docs/HANDOFF.md`, section "Hybrid (mixed
triangle + quad) meshes" → "Mixed-mesh setup workflow".

---

## 6. Output and reading results

Set `FESOM3_OUTPUT=<dir>` to enable output. Variables and writer options are read from `namelist.io`
in the run directory (see `config/namelist.io`); without it, a default variable set is written at the
`FESOM3_OUTPUT_EVERY` step cadence. Output is written as **Zarr** stores — one per variable per
period:

```
<dir>/ssh.fesom.1958.zarr
<dir>/temp.fesom.1958.zarr
<dir>/fesom.mesh.diag.zarr        # mesh geometry / connectivity (UGRID-style)
```

The stores are partition-independent (the same regardless of rank count) and open directly with
xarray:

```python
import xarray as xr
ds = xr.open_zarr("ssh.fesom.1958.zarr")   # CF time axis, embedded lon/lat
```

> **Format change (hybrid meshes).** `fesom.mesh.diag.zarr` now holds 4-column per-face arrays
> (dim `n4`: `face_nodes`, `face_edges`, `face_links`, `gradient_sca_x/y`). A triangle's 4th slot is
> the `_FillValue` (−1; −999 in `face_links`; coefficient 0 in `gradient_sca_x/y`), so xarray returns
> the integer arrays as float with NaN unless opened with `mask_and_scale=False`. Element `lon`/`lat` are now the cyclic-wrapped centroid (elements across
> the rotated dateline moved). Details: `docs/HANDOFF.md` → "Gate migration for hybrid meshes" →
> "Output format changes".

---

## 7. Repository layout

```
src/            Fortran source
  params/         precision, constants, configuration            (6 modules)
  types/          mesh / partition / dynamics / tracers / ice    (5)
  infra/          MPI halo layer, partitioning, clock, timers    (10)
  mesh/           mesh read, geometry, rotation                  (4)
  oce/            ocean physics kernels                          (27)
  ice/            sea-ice dynamics, thermodynamics, coupling     (6)
  forcing/        atmospheric forcing read + bulk fluxes         (3)
  io/             Zarr output, decomposition, restart            (8)
  step/           timestep orchestration                         (2)
  drivers/        program entry points + byte-gate drivers       (33)
mesh_part/      fesom_ini sources (FESOM2's mesh partitioner, mixed-mesh capable)
lib/metis-5.1.0/  vendored METIS used by fesom_ini
cmake/          compiler-flag and helper modules
config/         namelist.io (output configuration)
env/            per-host module/setup files (levante.dkrz.de)
test/           ctest self-tests
tools/          byte-gate scripts, oracle runner, diff tools, sbatch jobs
  mesh/           Python mesh tools (nemo_meshq, remove_deadends, split_quads,
                  cut_region, check_mesh, check_dist) + test_*.sh
docs/           HANDOFF, LESSONS, plans (internal working logs)
configure.sh    build helper          env.sh   environment loader
```

---

## 8. Documentation map

- **[IMPLEMENTATION.md](IMPLEMENTATION.md)** — architecture, design decisions, differences from FESOM2.
- **[TESTING.md](TESTING.md)** — validation approach, byte-gate methodology, results, scalability.
- **`docs/HANDOFF.md`** — the detailed internal development log (milestone-by-milestone state).
- **`docs/LESSONS.md`** — accumulated gotchas and subtle reproducibility traps.
- **`docs/plans/`** — per-milestone implementation plans (and `completed/`).
