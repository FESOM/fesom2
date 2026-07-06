# Known bug: wrong `sprod` in the reproducible-OpenMP SSH CG solver

**Status:** confirmed, not yet fixed (documented as a future-fix plan).
**Scope:** general correctness bug — affects **double *and* single precision**. Independent of the
single-precision workstream, though it compounds SP fragility.
**Activation:** only when built with `-DOPENMP_REPRODUCIBLE=ON` (which defines
`__openmp_reproducible`, `src/CMakeLists.txt:314-315`). Default builds are **not** affected — they
take the correct branch.

## Location
`src/solver.F90`, subroutine `ssh_solve_cg`, the "Scalar products for beta" block, the
`#else` (`__openmp_reproducible`) branch — currently lines 255-256:

```fortran
#else
    sprod(1) = sum(rr(1:myDim_nod2D) * zz(1:myDim_nod2D))   ! r . z
    sprod(1) = sum(rr(1:myDim_nod2D) * rr(1:myDim_nod2D))   ! r . r   <-- BUG: writes sprod(1) again
#endif
```

## What the code should compute
Each CG iteration needs two independent scalars:
- `sprod(1) = r . z` (preconditioned inner product) — used by `be = sprod(1)/s_old` and
  `s_old = sprod(1)` (lines 268-269).
- `sprod(2) = r . r` (residual norm squared) — used by the exit test
  `if (sqrt(sprod(2)/nod2D) < rtol) exit` (line 265).

The **default** (`#if !defined(__openmp_reproducible)`) branch does this correctly:
```fortran
sprod(1:2)=0.0_WP
DO row=1, myDim_nod2D
   sprod(1)=sprod(1)+rr(row)*zz(row)   ! r . z
   sprod(2)=sprod(2)+rr(row)*rr(row)   ! r . r
END DO
```

## The defect
In the reproducible branch, `sprod(1)` is assigned twice (the second assignment overwrites `r.z`
with `r.r`) and **`sprod(2)` is never assigned** on this path (nor zeroed — it is an uninitialised
local; the `sprod(1:2)=0` init lives only in the other branch). The subsequent
`MPI_Allreduce(sprod, 2, ...)` then also reduces the garbage `sprod(2)` across ranks.

## Consequences (general — DP and SP)
1. **Wrong `beta`.** `be = sprod(1)/s_old = (r.r)/(r_old.z_old)` instead of the required
   `(r.z)/(r_old.z_old)`. Preconditioned CG must use `r.z`; using `r.r` breaks conjugacy of the
   search directions → **slower convergence, stalling, or divergence** of the SSH solve.
2. **Error propagation.** `s_old = sprod(1) = r.r` feeds the next iteration's `al = s_old/s_aux`
   numerator, which should be `r.z` — compounding the error.
3. **Broken convergence test.** `sqrt(sprod(2)/nod2D) < rtol` reads uninitialised/stale `sprod(2)`:
   - **false early exit** → an *unconverged* SSH field silently feeds the barotropic/baroclinic
     dynamics (wrong physics, hard to notice), or
   - **never exits** → runs to the iteration cap every solve (wasted compute); and `sqrt` of a
     garbage-negative value yields `NaN`, making `NaN < rtol` always false — i.e. a guaranteed
     max-iterations solve.

Because reproducible-OpenMP builds are used specifically to certify bit-for-bit results across
thread counts, this is exactly the configuration where a silently-degraded solver is most harmful.

## SP relevance (cross-reference)
Separately, `ssh_solve_cg` computes all its dot products and step scalars
(`s_old, s_aux, al, be, sprod`) in `WP` and reduces them with `MPI_WP`. In single precision at
high resolution the SP global reductions lose enough precision to drive `s_aux = p.Ap` toward zero,
so `al = s_old/s_aux` blows up to Inf/NaN across the whole solution vector — the systemic
`vert_vel_ale` "Nan in Wvel" at step 1 (seen on ng5 for GNU and on CORE2 for Intel). The planned
mixed-precision fix (DP inner products + `MPI_DOUBLE_PRECISION`) touches this same block; the
`sprod(2)` typo should be corrected at the same time (but is a distinct, precision-independent bug).

## Fix plan
One-line correction:
```fortran
    sprod(2) = sum(rr(1:myDim_nod2D) * rr(1:myDim_nod2D))
```
Validation:
- Build with `-DOPENMP_REPRODUCIBLE=ON` and confirm the SSH solver converges in a sane iteration
  count and that results match a default (non-reproducible) build to solver tolerance.
- Ideally add a CI/local check that exercises the `OPENMP_REPRODUCIBLE=ON` path so this branch is
  not left untested again.
