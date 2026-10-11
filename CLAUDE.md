# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commit attribution policy

Only humans may be represented as authors or contributors.

Do not add AI tools or models as authors, co-authors, contributors, or
signatories in commit messages or pull-request descriptions. This includes
Claude, Codex, ChatGPT, Gemini, Copilot, Cursor, Aider, and Devin.

Never add `Co-authored-by`, `Generated-by`, `Assisted-by`, or similar
attribution referring to an AI system.

Before creating commits, ensure the repository hooks are installed:

    pre-commit install

Never bypass repository hooks with `--no-verify`.

If the hook rejects a commit, remove the prohibited attribution instead of
bypassing the hook.

## Build

FESOM2 is Fortran (plus a small C++ thread library) built with CMake and MPI; NetCDF (C and Fortran) is required.

    ./configure.sh ubuntu                       # platform name selects env/<platform>/shell
    ./configure.sh ubuntu -DENABLE_OPENMP=ON -DBUILD_TESTING=ON -DBUILD_MESHDIAG=ON -DBUILD_MESHPARTITIONER=ON
    cd build && make install -j                 # incremental rebuild

- `configure.sh` sources `env.sh` (guesses the platform from the hostname if none is given), runs cmake in `${BUILD_DIR:-build}` with the source dir as install prefix, so the binary lands in `bin/fesom.x`. Only arguments starting with `-` are forwarded to cmake.
- Options worth knowing (`CMakeLists.txt`, `src/CMakeLists.txt`): `ENABLE_OPENMP`, `OPENMP_REPRODUCIBLE`, `USE_SINGLE_PRECISION`, `ASYNCHRONOUS_IO_THREADS` (OFF everywhere; without it the async reader/writer threads run inline), `FESOM_WITH_XIOS` (forces async threads off), `FESOM_PROFILING` (ON; writes `fesom.stats` into ResultPath), `ENABLE_MULTIO`, `RECOM_COUPLED`, `USE_ICEPACK`. Presets in `CMakePresets.json`.

## Tests

    cd build
    ctest -N                                        # list
    ctest -R integration_pi_mpi2 --output-on-failure  # single test
    ctest -L conservation                           # needs -DENABLE_CONSERVATION_TESTS=ON and python3

- Tests live in `tests/` (integration, meshpartitioner, openmp, unit, perf); the harness is `cmake/FesomTesting.cmake` (`add_fesom_test`, `add_fesom_test_with_options`), which patches the `config/namelist.*` files into `build/tests/integration/<test>/` and runs `fesom.x` there. Pass/fail is decided by `cmake/CheckFesomRun.cmake` (exit code, success markers, required output files).
- Bundled data in `tests/data` (pi, pi_cavity, soufflet meshes; CORE2 forcing for 1948). Only `dist_2` and `dist_8` partitions of pi are bundled; other counts need `bin/fesom_meshpart`. The JRA55 test registers only if `tests/data/FORCING/JRA55` exists.
- The bundled CORE2 time axis (0, 6, 12, …) is read as days, because FESOM ignores the units attribute and divides by `nm_nc_freq` (1 in the test namelist).
- Conservation tests are judged by `tools/check_conservation.py` from the CONS log lines.
- Legacy pFUnit tests under `test/` are used only when `tests/` is absent.
- CI (`.github/workflows/`): `fesom2_ctest_runner.yml` runs ctest in `ghcr.io/fesom/fesom2_docker:fesom2_ci-master`; `fesom2_main.yml` runs `mkrun pi test_pi -m docker`, then `fcheck .` against `setups/test_pi/setup.yml`. Update reference values with `setups/update_truth_values.sh`.
- The only pre-commit hook is `tools/check_ai_attribution.py` (commit-msg). There is no Fortran linter.

## Architecture

- `src/fesom_main.F90` → `fesom_init` / `fesom_runloop` / `fesom_finalize` in `src/fesom_module.F90`. Init order: MPI (`MPI_THREAD_MULTIPLE` only with async I/O threads), `par_init`, namelists (`setup_model`), mesh, dynamics/tracer/arrays, `forcing_setup`, `ice_setup`, restart read. Each step: `update_atm_forcing` → ice → `oce_fluxes` → `oce_timestep_ale` → diagnostics, `output`, restart.
- State lives in derived types held by the module object `f`: `t_mesh` (`MOD_MESH.F90`), `t_partit` (`MOD_PARTIT.F90`: `mype`, `npes`, `MPI_COMM_FESOM`, `myDim_*`/`eDim_*` owned and halo sizes), `t_dyn`, `t_tracer`, `t_ice`. Routines unpack them with `#include "associate_part_def.h"` / `associate_mesh_def.h` plus the matching `_ass.h`.
- Domain decomposition is precomputed (`dist_N/` next to the mesh, from `fesom_meshpart`); `n_levels`/`n_part` in `namelist.config` describe it. Node loops run over `1, myDim_nod2D+eDim_nod2D` (owned + halo), with `!$OMP PARALLEL DO` around them.
- Timers: the `runtime …` table at the end of the log comes from `MPI_Wtime` stamps in `fesom_module.F90`; finer sections use `fesom_profiler_start/end` (`src/fesom_profiler.F90`, see `docs/developer_documentation/profiler.rst`). The profiler reduces sections across ranks by slot order, so every rank must enter every section, in the same order.
- I/O:
  - output: `io_meandata.F90`, parallel NetCDF4/MPI-IO when `&io_parallel parallel_write` is set
  - restarts: `io_restart.F90`, `io_fesom_file.F90`, with `io_redistribute.F90` for reader/writer schedules (`redist_build_scatter`/`redist_scatter`: a fixed schedule, readers own contiguous global blocks, one Alltoallv)
  - plus XIOS / MultIO interfaces
- Atmospheric forcing (`namelist.forcing`, `&nam_sbc`): `update_atm_forcing` (`gen_forcing_couple.F90`) → `sbc_do` → `getcoeffld` (`gen_surface_forcing.F90`).
  - `getcoeffld` runs once per new record per field. One I/O rank per field (`next_io_rank`) reads the record through `forcing_provider_async_module` (look-ahead block reader, `PREFETCH_SIZE`) and sends it out. With `forcing_distribution='bcast'` that is the full grid; with `'scatter'` each rank gets only its bilinear-stencil box.
  - The record is then interpolated to nodes into the `coef_a`/`coef_b` slope form. `data_timeinterp` evaluates it every step.
  - Monthly SSS, chlorophyll and runoff go through `gen_modules_read_NetCDF.F90` (rank 0 reads, then broadcast).
- Docs: `docs/` (Sphinx; `docs/forcing_configuration.rst`, `docs/getting_started/getting_started.rst`, `docs/developer_documentation/`).
