# SSH solver optimization — analysis & options (for a future dedicated branch)

Scope note: this is a **design/analysis note for later work**, not implemented here. The
`sp-cg-diagnostics` branch is deliberately kept to the diagnostics facility + the `soltol`
namelist option. The preconditioner / communication-avoiding work should live on its own branch
(suggested name `ssh-solver-precond`). Context: `docs/sp-scaling-diagnosis.md` showed the SSH CG
solve is the strong-scaling drag, at an SP≈DP iteration count (~180/solve on core2), so the win
is *fewer collectives* — dominated by *fewer iterations*.

## Allreduces per CG iteration: it's 2, and that's already optimal

Reading `src/solver.F90` `ssh_solve_cg` (loop `:201–284`):

- **Per iteration = 2 `MPI_Allreduce`:**
  - `:226` — `p·Ap` (1 scalar) → needed for `alpha`, *after* the mat-vec `App` (`:206-209`).
  - `:264` — one **2-element** reduce of `(r·z, r·r)` → needed for `beta` and the convergence
    test; `r·r` is folded into the same collective for free.
- **Plus 2 one-time setup reduces** before the loop: `:157` (‖rhs‖²) and `:195` (initial `r·z`).

So a solve costs `2 + 2·niter` reductions. **2/iteration is the classic-PCG minimum** — you need
`p·Ap` to form `alpha` before updating `r`, and `r·z` afterward to form `beta`; they are separated
by the mat-vec and cannot be trivially merged. (An earlier note that said "4/iteration" miscounted
by including the two setup reduces.)

The only sub-2 option is a **communication-avoiding / pipelined CG** (Chronopoulos–Gear): restructure
the recurrences to a single reduce/iteration and overlap it with the mat-vec. Real but nontrivial,
with a mild numerical-stability cost, and it only buys ~2×. Lower priority than cutting iterations.

## Preconditioner: the current one is weak (root cause of ~180 iters)

`ssh_solve_preconditioner` (`src/solver.F90:31–95`) is a **single symmetrized-Jacobi sweep**
(MITgcm, JGR 1997): `M⁻¹` stored as diagonal `1/a_r` and off-diagonals
`-0.5·(a_i/a_r)/(a_r + a_diag_i)`, applied as one sparse mat-vec (`zz = pr_values·rr`, `:176,245`).
This is only marginally stronger than diagonal scaling, so ~180 iterations for a global elliptic
SSH solve is expected. Reducing this iteration count is the highest-value scaling lever (it cuts
the `2·niter` collectives that dominate at high task counts, for both precisions).

### Options, ranked by fit to FESOM's structure (CSR mat-vec + `exchange_nod` halo, OpenMP/`!$ACC`)

1. **Chebyshev polynomial preconditioner — recommended.** Apply `M⁻¹ ≈ p_k(A)` as `k` SSH-operator
   mat-vecs with Chebyshev coefficients from eigenvalue bounds. Fit: the apply is **pure mat-vec +
   halo** (the ops that already scale), adds **zero new global reductions**, and is OpenMP/GPU
   friendly (no sequential triangular solve). Bounds are cheap: `λ_max` via a few power iterations
   once (store in `solverinfo`), `λ_min` via a Gershgorin/diagonal estimate. `k=2–4` typically cuts
   iterations 2–3×. Net high-np effect: fewer CG iterations → fewer collectives, extra work being
   scalable mat-vecs.
2. **k-sweep symmetrized SSOR/Jacobi.** Trivial extension of the current code (loop the existing
   `pr_values` apply `k` times, `exchange_nod` between sweeps). Small effort, modest iteration
   reduction, but each apply adds `k` halo exchanges (latency) → smaller net high-np win than (1).
3. **Block-Jacobi ILU(0)** (per-rank incomplete factorization). Stronger per iteration, purely
   local (no extra halo), but the triangular solves are **sequential** → fights the OpenMP/`!$ACC`/
   GPU parallelism the code is built around. Good on paper, awkward here.
4. **Multigrid / multilevel** on the SSH operator. Strongest (iteration count ~O(1),
   resolution-independent — attractive for ng5/high-res), but a large unstructured-mesh effort.
   Long-term.

**Recommendation:** prototype the **Chebyshev** preconditioner behind a namelist switch (default =
current symmetrized-Jacobi), validate SP-vs-DP fields on pi + core2-1day, and measure `ssh_cg_iters`
(now available) and SYPD across task counts. It composes independently with R3 (split-explicit
subcycling), which removes the elliptic solve's collectives entirely.

## Related levers (from `docs/sp-scaling-diagnosis.md`)
- R3 — split-explicit barotropic subcycling: removes all ~360 elliptic collectives/timestep.
- R1b — communication-avoiding CG (2→1 reduce/iter): smaller, independent.
- R1c — DP-accumulated dot-products: accuracy only; does **not** change iteration count (measured).
