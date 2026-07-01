# FESOM profiler scalar-diagnostics facility (as built)

Branch `sp-cg-diagnostics` (off `support_sp`).

## What it is

A small, reusable way to record arbitrary **per-call scalar diagnostics** (solver iterations,
residuals, CFL numbers, event counts) and have them reported, MPI-reduced, into `fesom.stats`
next to the timing table. It mirrors the `fesom_profiler` patterns (name→index registry,
per-entry accumulation, reduce-then-print) but accumulates *values* instead of wall time.
Storage is double precision (`PROF_WP`) regardless of the model `WP`, so SP and DP diagnostics
are directly comparable.

## Files

- **`src/fesom_diagnostics.F90`** (new; auto-globbed by `src/CMakeLists.txt`): type
  `diag_stat{name,sum,sumsq,min,max,count}`, registry `diags(MAX_DIAGS=64)`,
  `find_or_create_diag`, `diag_init/enabled/reset`, generic **`diag_count(name,value)`**
  (specifics for `integer`, `real(real32)`, `real(real64)` — no cast at call sites in SP or DP),
  and **`diag_report(unit,comm,rank)`** (collective; reuses the profiler's reduce pattern).
- **`src/fesom_profiler.F90`**: `use fesom_diagnostics`; `diag_init` in init, `diag_reset` in
  reset, and `diag_report` **before** `print_detailed_report` so the block prints *above* the
  timing table (parser-safe — `compare_stats.py` keys on the timing rows).
- **`src/solver.F90`**: records `ssh_cg_iters` (and `ssh_cg_nonconv`) at the end of
  `ssh_solve_cg`, guarded by `#if defined(FESOM_PROFILING)`; plus an always-on `mype==0`
  non-convergence WARNING at `maxiter` (previously a silent stall).

## Usage

```fortran
use fesom_diagnostics
call diag_count("ssh_cg_iters", iter)          ! integer
call diag_count("bt_max_cfl", real(cfl, 8))    ! real
```
Guard hot call sites with `#if defined(FESOM_PROFILING)` for zero cost when profiling is off.
Solver-agnostic: with split-explicit subcycling `ssh_cg_iters` simply won't fire; that path can
record its own scalars (e.g. `bt_max_cfl`) with the same one-liner.

## Output (top of `fesom.stats`)

```
================ PER-CALL DIAGNOSTICS (scalar counters) ================
                        Name     Mean/call     Min      Max      Std     Count
ssh_cg_iters                       177.3750  162.0000  195.0000  8.7384   512
```
`Mean/call = Σsum/Σcount`; Min/Max across ranks & calls; Std from sumsq.

## Validation

- Builds clean SP & DP with `-DFESOM_PROFILING=ON`; module dependency order handled by CMake.
  Without profiling the `diag_count` calls compile out (DP numerics bit-unchanged; only an
  always-on non-convergence warning is added, which does not alter results).
- pi mpi2 and core2 mpi16 (1-day) both emit the block. **First result** (core2, 127k nodes):
  DP mean 181.0 iters/solve, SP 177.4, zero non-convergence — used to overturn the
  "SP does more CG iterations" hypothesis (see `docs/sp-scaling-diagnosis.md`).

## Optional follow-ups

- Extend `scripts/pi_sp_dp/compare_stats.py` + `scripts/levante_sp_dp/collect_scaling.py` to
  parse the PER-CALL DIAGNOSTICS block so diagnostics flow into the scaling analysis.
- Add more diagnostics as needed (`ssh_cg_resid`, `bt_max_cfl`, blow-up counts).
