# SP strong-scaling on core2: diagnosis and improvement plan

## The observation

core2 strong scaling, SP (`WP=4`) vs DP (`WP=8`), 128-core Levante nodes. Efficiency is
relative to the smallest count (64 tasks); ideal = 1.

| ntasks | SP SYPD | DP SYPD | SP/DP | SP eff | DP eff |
|---|---|---|---|---|---|
| 64  | 11.60 | 6.40 | 1.81 | 1.00 | 1.00 |
| 72  | 13.48 | 7.42 | 1.82 | 1.03 | 1.03 |
| 128 | 25.34 | 14.54 | 1.74 | 1.09 | 1.14 |
| 144 | 27.53 | 16.35 | 1.68 | 1.06 | 1.14 |
| 256 | 48.87 | 30.41 | 1.61 | 1.05 | 1.19 |
| 288 | 53.80 | 33.94 | 1.59 | 1.03 | 1.18 |
| 432 | 77.05 | 51.57 | 1.49 | 0.98 | 1.19 |
| 512 | 87.89 | 60.10 | 1.46 | 0.95 | 1.17 |
| 864 | 115.39 | 89.03 | 1.30 | 0.74 | 1.03 |

Two things to explain: **(1)** SP/DP falls 1.81 → 1.30 with task count; **(2)** SP efficiency
rolls off to 0.74 while DP stays superlinear (~1.1–1.2) almost to 864.

## What it means (performance model)

Two effects are expected and largely fundamental; a third is SP-specific and fixable.

**(a) DP's superlinear efficiency is a cache effect.** DP arrays are twice the bytes, so DP is
more memory-bandwidth-bound. As tasks increase, each rank's subdomain shrinks and begins to fit
in cache, so DP keeps speeding up *faster than linear* (eff 1.1–1.2). SP arrays are already half
the size, so SP is *less* memory-bound and has a *smaller* cache dividend — it has already
"cashed in" the benefit DP is still collecting as domains shrink. That is why DP efficiency sits
above SP efficiency across the whole sweep.

**(b) The falling ratio and SP rolloff are Amdahl's law on latency-bound communication.** SP's
~1.8× advantage at low task counts is the compute (bandwidth-halving) win. Communication does
**not** get faster in SP at high task counts: the halo exchanges and the solver's global
reductions become *latency-bound* (dominated by message count and collective latency, not by
bytes), and SP's half-size messages don't reduce latency. Because SP compute is ~2× faster, the
communication fraction becomes dominant at a *lower* task count for SP — so SP hits the scaling
wall first, its efficiency falls, and the SP/DP ratio converges toward 1 while DP (slower
compute, larger cache dividend) keeps scaling.

**(c) An SP-specific pathology in the SSH solver steepens SP's rolloff** — see below. This is
the part worth fixing.

## Root cause in the code

The scaling runs used `use_ssh_se_subcycl = .false.`
(`config/namelist.dyn:67`, and the run namelists), so FESOM used the **preconditioned CG SSH
solver** `ssh_solve_cg` in `src/solver.F90` — *not* the halo-only split-explicit subcycling.

1. **CG does 4 `MPI_Allreduce` per iteration** (`src/solver.F90:153,191,221,259`), each a scalar
   in `MPI_WP`. Global reductions are the classic latency-bound, non-scaling cost; their count
   is `4 × (CG iterations)` per timestep, on top of ~4 other per-step collectives (item 3). At
   high task counts these dominate the barotropic part of the step for both precisions.

2. **The CG inner products and residual accumulate in `WP`** — real(4) under SP:
   `s_old = s_old + rhs*rhs` with `s_old` real(WP) (`solver.F90:143-147`), reduced in `MPI_WP`
   (`:153`), and the convergence test `sqrt(sprod(2)/nod2D) < rtol` with
   `rtol = soltol*sqrt(s_old/nod2D)`, `soltol = 1e-5` (`src/MOD_DYN.F90:18`) — all `WP`.
   Summing ~126,858 nodal products in real(4) gives a residual-norm noise floor of order
   `ε·√N ≈ 1.2e-7 · 356 ≈ 4e-5`, which is **larger than the tolerance `1e-5`**. So in SP the
   convergence test sits *below its own round-off noise floor*: near the tolerance the CG
   residual estimate is noisy, so the solver stalls / behaves erratically and tends to take
   **extra iterations** (bounded by `maxiter = 2000`). Every extra iteration is another 4
   `MPI_Allreduce` — latency-bound collectives that cost proportionally more at high task
   counts. Net effect: **SP's `oce_ssh_solve` share grows faster with task count than DP's,
   steepening SP's efficiency rolloff.** This is the SP-specific mechanism.

3. **Other per-timestep collectives** (latency-bound, precision-independent): the blow-up check,
   1× `MPI_Allreduce` (`src/write_step_info.F90:671`), and the reference-profile sums, 3× in
   `MPI_WP` (`src/oce_ale_pressure_bv.F90:3354-3356`).

4. **Halo exchange** (`src/gen_halo_exchange.F90`): ~123 nonblocking point-to-point exchanges
   per timestep, **one Isend + one Irecv per neighbouring rank** (neighbours capped at 32,
   `src/MOD_PARTIT.F90:21`). SP buffers are correctly half-size via `MPI_WP`, but message
   *count / latency* is precision-independent → latency-bound at high task counts. No
   per-timestep communication is hard-forced to DP (CVMix is off in SP; only profiler / CMOR /
   coupling use `MPI_DOUBLE_PRECISION`, none per step).

5. **Partitioning** (`src/fort_part.c:107`) minimises **edge-cut** (`METIS_OBJTYPE_CUT`,
   recursive), which bounds communication *volume* but not neighbour *count*, so the message
   count (latency) stays high at large task counts.

**On "should SP use a different partition?"** No. The mesh topology is identical, so there is no
SP-specific partition, and SP messages already halve via `MPI_WP`. The partition *objective* is a
lever (R4) that helps both precisions — and disproportionately SP, because SP is
communication-bound sooner — but there is no per-precision partition to make.

## Corroboration from the downloaded runs (`core2_sp_dp/core2_scaling`)

Parsing all 18 `results/fesom.stats` (via `pi_sp_dp/compare_stats.py::parse_stats`) confirms the
hypothesis and cleanly separates the two mechanisms.

**SSH-solve share of total runtime rises with task count and is larger for SP:**

| tasks | ssh% SP | ssh% DP |
|---|---|---|
| 64 | 4.0 | 3.4 |
| 128 | 5.4 | 3.7 |
| 256 | 8.9 | 6.2 |
| 512 | 15.1 | 11.2 |
| 864 | 25.5 | 19.9 |

**Per-component strong-scaling efficiency (vs 64 tasks) — SSH-solve collapses; the rest scales
(DP superlinearly):**

| tasks | ssh eff SP | ssh eff DP | rest eff SP | rest eff DP |
|---|---|---|---|---|
| 64 | 1.00 | 1.00 | 1.00 | 1.00 |
| 128 | 0.80 | 1.05 | 1.11 | 1.14 |
| 256 | 0.47 | 0.66 | 1.11 | 1.22 |
| 512 | 0.25 | 0.36 | 1.07 | 1.28 |
| 864 | 0.12 | 0.18 | 0.95 | 1.24 |

Key facts from the data:
- **The non-solver bulk scales well** — DP *superlinearly* (1.14→1.28, the cache effect),
  SP near-linear (1.07–1.11, smaller cache dividend because SP is less memory-bound). Confirms
  performance-model points (a) and (b).
- **The CG SSH solve is anti-scaling** — efficiency falls to 0.12–0.18, and its *absolute* wall
  time **increases** from 512→864 tasks (SP 23.9→30.9 s), the signature of a
  global-reduction/latency-bound component that does not scale.
- **It hits SP harder** — worse per-component efficiency at every count (0.12 vs 0.18 at 864)
  and a larger runtime share (25.5% vs 19.9%). The tell that SP takes ≥ DP iterations: absolute
  SSH time *converges* (SP 30.9 ≈ DP 31.1 s at 864) despite SP's ~2× faster local matvec — SP's
  speed advantage in that kernel is entirely consumed.
- **Attribution:** SP overall efficiency at 864 is 0.74 while its non-SSH part runs at 0.95, so
  the SSH solver accounts for essentially all of SP's rolloff. If SSH scaled like the rest,
  SP@864 would be ~0.93 efficiency (≈145 SYPD vs 115, ~26% more throughput) — the target of the
  collective-reduction fixes (R1a/R1b/R3).

Caveats: the run logs do **not** print CG iteration counts (only timestep-loop barrier lines),
so "SP does more iterations" was, at this point, only implied by the converging absolute times.
The `run.log` BENCHMARK block also shows the SSH solve is load-imbalanced at low task counts
(dp_n64: mean 41 s, max 80 s) but becomes uniform-and-slow at 864 (all ranks blocked on
reductions), consistent with collective-latency dominance rather than compute imbalance at scale.

### MEASURED: CG iteration counts (local core2 1-day, NP=16) — hypothesis overturned

With the new `ssh_cg_iters` diagnostic (see `docs/profiling-diagnostics-plan.md`), a local core2
1-day run (127k nodes) gives:

| build | mean CG iters/solve | range | non-convergence |
|---|---|---|---|
| DP | 181.0 | 164–203 | 0 |
| SP | 177.4 | 162–195 | 0 |

**SP does *not* do more iterations — marginally fewer, and never fails to converge.** The
pessimistic `ε_sp·√N ≈ 4×10⁻⁵ > soltol=1e-5` estimate did **not** materialize: the
Jacobi-preconditioned *relative* residual (a ratio) reaches 1e-5 in SP. On the small pi mesh SP
and DP also matched (12.3 vs 12.3). So:

- **The SP scaling penalty is not extra CG iterations.** The equal SSH-solve wall time at 864
  tasks is explained by *equal iterations × latency-identical scalar `MPI_Allreduce`* (a
  1-scalar reduce costs the same for `MPI_REAL` and `MPI_DOUBLE_PRECISION` — it's pure latency).
  So the SSH solve is the scaling drag because it is **collective-latency-bound at a fixed
  iteration count**, and since SP compute is ~2× faster, that fixed collective cost is a larger
  *fraction* of SP's runtime → SP efficiency rolls off first. This is pure Amdahl, not an SP
  numerical pathology.
- **~180 CG iterations per SSH solve is itself very high** for a Jacobi-preconditioned system —
  that is 180 × (4 `MPI_Allreduce`) ≈ 720 latency-bound collectives per timestep, for both
  precisions. Cutting that count is the real lever.

This reprioritizes the fixes below.

## Recommendations (ranked — updated after measuring iteration counts)

The scaling drag is the **number of latency-bound global reductions per timestep**
(≈180 CG iters × 4 `MPI_Allreduce` ≈ 720), at an SP-vs-DP-identical iteration count. The levers
therefore target *reducing collectives*, not SP round-off. They help both precisions and, because
SP is comm-bound sooner, help SP's rolloff most.

### R1a — Fewer CG iterations via a better preconditioner  *(biggest lever)*
~180 Jacobi-PCG iterations/solve is high. A stronger preconditioner (block-Jacobi / additive
Schwarz / a few-level multigrid on the SSH operator) that cuts iterations ~2–4× directly removes
that fraction of the per-step collectives for **both** precisions. Larger effort, largest payoff.

### R1b — Fewer reductions *per* CG iteration
The current CG issues 4 `MPI_Allreduce`/iteration (`solver.F90:153,191,221,259`); a standard CG
needs 2, and communication-avoiding variants (Chronopoulos–Gear / pipelined CG) need **1** and
overlap it with the mat-vec. Reducing 4→1–2 cuts the dominant high-np cost with the same numerics.

### R1c — (was R1) DP-accumulated inner products — *accuracy only, not scaling*
Making the CG dot-products/residual `real(8)` **does not** change the iteration count (measured
SP≈DP), so it does **not** help scaling. Keep it only as an optional SP *accuracy* improvement
(tighter residual), low priority.

### R2 — Cut the per-timestep collective count
Fuse the three reference-profile `MPI_Allreduce` (`oce_ale_pressure_bv.F90:3354-3356`) into a
single length-`3*(nl-1)` reduction; consider making the blow-up check
(`write_step_info.F90:671`) non-blocking or every-N-steps. Helps both precisions at high task
counts.

### R3 — enable split-explicit barotropic subcycling  *(largest structural win)*
`use_ssh_se_subcycl = .true.` (`namelist.dyn`, `se_BTsteps = 50`) replaces the
reduction-heavy CG with halo-only subcycling — **no `MPI_Allreduce` in the barotropic loop** —
removing *all* ~720 elliptic-solve collectives/timestep at a stroke. Biggest scaling win for both
precisions (and it moots R1a/R1b). It changes the barotropic numerics, so it needs its own
SP-vs-DP validation and `se_BTsteps` re-tuning — but given the measurement (the CG collectives,
not SP round-off, are the drag), this is the highest-leverage direction, not merely strategic.

### R4 — Partition objective for latency
Try `METIS_OBJTYPE_VOL` with `METIS_PartGraphKway` (`src/fort_part.c:107`) to reduce
interface-node / neighbour count (message count), trading a little load balance. Helps latency at
high task counts for both precisions, SP more.

### R5 — Practical scheduling *(immediate, no code)*
SP's efficiency sweet spot is at **lower task count than DP**. On core2, SP at 128–256 runs at
~1.05–1.09 efficiency (48.9 SYPD at 256); pushing SP to 864 wastes cores (0.74). Run SP jobs at
fewer cores than DP for best core-hour efficiency and reserve high task counts for DP.

## Verification

1. **Confirm the mechanism from data already in hand (no new runs).** Extend
   `scripts/levante_sp_dp/collect_scaling.py` — it already parses `fesom.stats` sections via
   `scripts/pi_sp_dp/compare_stats.py::parse_stats` (populates `rep.sections`) — to plot, across
   task count, the absolute time and %Total of `oce_ssh_solve`, `oce_tracer_solve`,
   `ice_timestep`, and `output`/`restart`, for SP and DP. **Expected signature:** the
   `oce_ssh_solve` share rises with task count and rises *faster for SP* than DP. (The example
   `fesom.stats` already exposes a `>oce_ssh_solve` row.)
2. **Confirm extra SP iterations (optional, tiny change).** Print the CG exit `iter`
   (`solver.F90`, after the `Do iter = 1, maxiter` loop) and do one short SP and one DP run at
   512 tasks; compare mean iterations per timestep. Expect SP > DP today.
3. **After R1.** Re-run task counts {256, 512, 864} in SP; expect the SP `oce_ssh_solve` share
   and CG iteration count to fall to ≈ DP, SP efficiency at 864 to recover toward ~0.9, and the
   SP/DP ratio to hold nearer ~1.8 further out. The DP build must remain bit-for-bit unchanged
   (regression guard). Re-plot with `collect_scaling.py`.

## Scope / sequencing

Start with verification step 1 (the `collect_scaling.py` section breakdown) to confirm the
`oce_ssh_solve` signature from the runs already done. Then implement **R1** (contained,
reversible, DP-neutral) as the concrete deliverable. R2/R4 are small follow-ons; R3 is a separate
numerics study; R5 is immediate operational guidance. Relatedly, the geometry-precision work in
`docs/geometric-precision-plan.md` is orthogonal (accuracy, not scaling) but shares the theme of
keeping selected DP where round-off matters — the CG reductions here are exactly such a place.
