# FESOM_F — Testing and Validation

This document describes how FESOM_F is tested, how it is compared against FESOM2, and what the
results are — correctness and scalability. For the architecture being tested, see
[IMPLEMENTATION.md](IMPLEMENTATION.md); for building and running, see [README.md](README.md).

> Naming, as elsewhere: the project is **FESOM_F** in prose; real identifiers (`FESOM3_*`, `fesom_*`,
> tags `m0…m10`) keep the repository's working name `fesom3`.

---

## Part I — Overview

FESOM_F is validated against one demanding standard: **`max|Δ| = 0`** — bit-for-bit identical results
to FESOM2. Because the numerical kernels are transcribed verbatim (see
[IMPLEMENTATION.md](IMPLEMENTATION.md)), the two codes should agree to the last bit, and the tests
assert exactly that. A non-zero difference is treated as a defect to be found and removed, not a
tolerance to be accepted.

Validation rests on two pillars:

1. **Automated self-tests** (`ctest`) — 153 cases over the infrastructure (types, partitioning, halo
   exchange, I/O decomposition, restart redistribution, …), the kernels and conservation, on triangle,
   quad and mixed meshes, run on 1 to 8 ranks, on both the Intel and GNU toolchains.
2. **Byte-gates against the FESOM2 oracle** — an instrumented FESOM2 v2.7.3 and FESOM_F are run on
   identical input and their internal fields are compared substep by substep, demanding `max|Δ| = 0`.
   Gates exist for every subsystem, on single and multiple ranks, on a small test mesh and on the
   global CORE2 mesh.

A third property is checked throughout: **partition independence** — results, output, and restart
files do not depend on the number of MPI ranks.

**Headline results.** Milestones M0–M10 are byte-exact against FESOM2 on 1 rank and on up to 864
ranks, for both sea-ice solver variants. A full-physics run reproduces FESOM2 bit-for-bit over a
model year (1.1 million compared records), and a free-running two-year run is physically stable.
Restart resumes byte-exactly. On performance, FESOM_F is at or ahead of FESOM2 per timestep at scale.

---

## Part II — Reference

### 1. Philosophy: byte-exact validation and its floors

"Byte-exact" (`max|Δ| = 0`) means two independently compiled programs produce identical bits for
every compared value. This is achievable only because the kernels are transcribed verbatim and the
build and runtime are matched to FESOM2 (compiler flags, reduction order, floating-point mode — see
[IMPLEMENTATION.md](IMPLEMENTATION.md) §11). It is a stronger and more useful bar than "agrees to
N digits": it catches the class of subtle plumbing bugs that a tolerance would hide.

There are, however, a few honest **reproducibility floors** — places where exact equality is either
physically meaningless or inherited from FESOM2:

- **Denormals and flush-to-zero.** Bit-identity depends on the floating-point *runtime* mode, not
  just the source. A 2-year run once diverged because FESOM_F kept a denormal value FESOM2 flushed to
  zero; matching the mode (FTZ on) restored equality (LESSONS L51).
- **Reduction order.** Global sums are bit-identical only because FESOM_F mirrors FESOM2's exact MPI
  reduction; this is engineered, not assumed (LESSONS L6).
- **Shared (multiply-owned) elements at multiple ranks — resolved 2026-10-12.** FESOM2 assigns an
  element to a rank if *any* of its nodes is owned, so on CORE2/dist_2 **562 elements** are "owned" by
  more than one rank (24,798 at dist_512). Their `uv` copies used to differ by ~1 ULP from
  order-dependent sums: the restart gate admitted that as a 1e-12 floor, and in long CORE2 runs the
  `rvo_upwind` blend amplified it into a volume and salt leak (LESSONS L70). Since the shared-element
  sync (`FESOM3_SHARED_SYNC=1`, the default) every owner takes the designated owner's copy each step —
  the copy the partition-independent stores keep — so np>1 restarts are exact. The floor applies only
  with `FESOM3_SHARED_SYNC=0` — see §7–8.

### 2. The FESOM2 oracle

The reference ("oracle") is **FESOM2 v2.7.3** (git SHA `9271ae92`), at
`/home/a/a270088/port2/fesom2`, instrumented with dump shims that write its internal fields at each
substep. It is run with `tools/run_oracle_pi.sh`, which sets up a run on the small **pi** mesh and
emits per-rank, gid-keyed dumps. The oracle is deterministic: a fresh run reproduces a committed
reference fixture bit-for-bit, which is what lets it serve as a fixed target.

### 3. The dump / gid-keyed comparison

Both codes write the same binary record (`src/infra/mod_dump.F90` and the FESOM2 shim), one record
per probed field at each substep:

```
int32 step | int32 substep_id | int32 probe_gid | int32 nlevels | char[24] field_name | float64 values[nlevels]
```

The key is **gid** — the *global* node or element id — so a value can be matched across runs
regardless of how the mesh was partitioned. `substep_id` names the point in the timestep
(`PRESSURE_BV`, `MIXING`, `SSH_SOLVE`, `TRACERS`, …). `tools/dump_diff.py` groups records by
`(step, substep, gid, field)`, reports the per-record `max|Δ|`, the first diverging substep, and a
magnitude histogram that separates real divergence from floating-point noise. Its options drive the
gates:

- `--glob` merges the per-rank dump files (multi-rank comparison keyed on gid),
- `--threshold T` sets the tolerance (default `0` — byte-exact),
- `--ignore-substep` skips a substep a given configuration does not compute,
- `--ignore-step` compares the *same physical state* reached via different per-segment step counts —
  the restart gate (straight run vs. checkpoint-and-resume).

### 4. Gate methodology

Every byte-gate follows the same three steps:

1. **Run the oracle** — instrumented FESOM2, dumping its fields.
2. **Run FESOM_F** on the *identical* input state, dumping the same fields.
3. **Compare** with `dump_diff.py`, requiring `max|Δ| = 0`.

For example, `tools/run_step_gate.sh` runs the real FESOM2 `oce_timestep_ale` and the FESOM_F
`step_oce` on identical single-rank pi state and compares every per-substep node dump (density,
pressure, `N²`, `Kv`, `ssh_rhs`, `d_eta`, `hbar`, `eta_n`, `w`, `T`, `S`, layer thickness). Because
each kernel consumes the previous kernel's live output, a clean `max|Δ| = 0` proves the whole
assembled step — data flow and dispatch — is byte-faithful, not just the kernels in isolation.

### 5. The `ctest` self-test suite

`test/CMakeLists.txt` registers **167 ctest cases**, run on the Intel and GNU builds. The core
infrastructure cases:

| Test | Ranks | Checks |
|---|---|---|
| `test_params`, `test_types` | 1 | constants/precision parameters; derived-type definitions |
| `test_partit` | 1, 2, 8 | partition invariants (owned nodes partition exactly; elements/edges boundary-redundant) |
| `test_halo` | 1, 2, 8 | halo exchange identity; stale-halo corruption probe |
| `test_dump` | 1 | the gid-keyed dump binary round-trip |
| `test_mesh` | 1 | mesh infrastructure; depth-independent scalar-cell area |
| `test_bottom` | 1, 2, 8 | the bottom-at-vertices contract (§6a); element bounds derived from vertex columns, edge interval, `hnode` column sums, `edge_len` units |
| `fesom_conserve_zstar` | 1, 2 | conservation + no-leakage gate (§6a) |
| `test_io_decomp` | 1, 2, 8 | output redistribution to writer subset, chunking, partial chunks |
| `test_io_decomp_gather` | 1, 2 | the inverse (restart-read) redistribution |
| `test_io_posix` | 1 | the `bind(C)` POSIX shims (rename/fsync/rmtree) |
| `test_vector_rotate` | 1 | vector rotation invertibility (output velocity to geographic) |
| `test_io_means` | 1 | output writer lifecycle + per-field cadence |
| `dump_diff_selftest` | — | the comparison tool itself (synthetic perturbation detection) |
| `fesom_analytic` | 1, 2 | end-to-end analytic driver runs clean |
| `test_shared_elem` (+ `_corrupt`, `_badknob`, `_off`) | pi 1/2/8, pi_w 2/4, cut-outs 2, core2_w 8 | shared elements: the list `com_elem2D_shared` vs an allgather oracle (S1), its exchange and the copy-spread probe (S2), the mechanism through `update_vel` + `compute_hbar_ale` with the sync on and off (S3: divergence-sum ratio ~1e-17 vs 1e-5..4e-4); the error stops and `FESOM3_SHARED_SYNC=0` |
| `fesom_conserve_linfs_vinv_shsync` | 2 | linfs, vector invariant with `rvo_upwind` 1, sea ice: divergence-sum ratio < 1e-12 and zero shared-copy spreads in every step (`FESOM3_DIVSUM_TOL`) |

Kernel tests from the bottom-at-vertices and momentum work (`test_vert_spline`, `test_vinv`,
`test_muscl_onthefly`, `test_ivertvisc`, `test_wsplit`, `test_wimpl_tra`, `test_minlayers`, the
`fesom_conserve_zstar_*` variants) complete the pre-hybrid set. **Hybrid (triangle + quad) meshes**
add (details: `docs/HANDOFF.md` → "Hybrid (mixed triangle + quad) meshes" → "Tests and tools of the
plan"):

| Test | Meshes | Checks |
|---|---|---|
| `test_mesh_read` (+ `_malformed[2-5]`) | quad/mixed cut-outs, pi_w | 3/4-column `elem2d.out`, orientation (flipped quad/triangle turned clockwise), levels, `dist_N/format.out`; each malformed-file case stops |
| `test_partit` (hybrid cases) | quad/mixed cut-outs, pi_w | appended `dist_N` sections, `format.out` required (`FESOM3_EXPECT_FORMAT=1`) |
| `test_geom_mixed` | analytic patches, pi_w | centroids, kite weights `elem_w` (vs independent shoelace kites), Green–Gauss `gradient_sca`, node areas, `mesh_resolution` |
| `test_halo_ext` (+ `_noformat`) | cut-outs, pi_w, core2_w at np 2/4/8 | `exchange_edge` / `exchange_nod_full` round trips; FESOM2-layout `dist_N` error |
| `test_ssh_stiff` | patches, mixed cut-out, pi_w, core2_w | SSH matrix pattern and values (quadratic field on the parallelogram patch), symmetry, triangle bit-identity, quad checkerboard null mode |
| `test_pgf_mixed`, `test_updown_mixed`, `test_coupling_mixed` | patches, cut-outs, pi_w, core2_w | pressure gradient, density gradient / smoother / GM bolus / IC extrapolation, MUSCL up/downwind elements, surface-stress gather |
| `test_ice_dyn_mixed`, `test_ice_mass` (+ `_missing_col`) | patches, quad cut-out, pi_w, core2_w | EVP strain/stress adjoint, ice FCT mass matrix, lumped mass, conservation and bounds |
| `fesom_ini_{aux3d,nodepth,onelem,count1,nonconvex}` | copies of the quad cut-out | `fesom_ini` level finder on quads, its error stops (`tools/check_fesom_ini.sh`) |
| mixed/quad cases of `test_vinv`, `test_ivertvisc`, `test_muscl_onthefly` | patches, cut-outs | the converted ocean kernels |
| `meshdiag_{mixed,quad,pi_w}` | cut-outs, pi_w | 4-column mesh-diagnostics store (`tools/check_meshdiag.py`) |
| `fesom_conserve_mixcut_*`, `fesom_conserve_{q,mix}cut_ice_*` | cut-outs, np 1/2 | conservation (PP/KPP/TKE/GM/Redi; with sea ice, EVP and mEVP) |

The mesh-dependent hybrid cases read meshes from the hard-coded root `/work/ab0995/a270029/meshes`
(`test/CMakeLists.txt`). Whether a mesh exists is decided at configure time
(`if(EXISTS <meshes>/<mesh>/nod2d.out)`): a case whose mesh is absent is registered as DISABLED
(ctest prints "Not Run (Disabled)"); a missing `dist_N`/`com_info` file in an existing mesh
directory FAILS. "167/167 passed" therefore holds on a filesystem that has these meshes. The
shared-element guards live in existing tests: `test_vinv` V13 (`uv_rhsAB` copies bitwise equal with
the upwind blend on) and `test_ice_dyn_mixed` D6 (EVP/mEVP sigma copies bitwise equal), at np>1.

These run in under a minute on a login node (`ctest -j 8`):

```bash
cd build_intel_dp && ctest --output-on-failure
```

### 6. The byte-gate catalog

> **Bottom at vertices retires the oracle-dependent gates.** FESOM_F now derives each element's
> bottom from its vertex columns instead of reading `elvls.out`, which makes the model
> deliberately different from FESOM2 — so every gate that compares against the **oracle** is no
> longer meaningful and is retired: `run_step_gate.sh`, `run_geom_gate*.sh`,
> `run_lifecycle_*_gate_*`, `run_advhor_gate*`, `run_pressure_gate*`, `run_ic*_gate*`,
> `run_ice*_gate*`.
>
> The **self-consistency** gates are not oracle-dependent and remain the live net:
> `run_conserve_pi.sh` (§6a), `run_restartroundtrip.sh`, `run_output_gate.sh` (partition
> independence), the restart straight-vs-resume gates, and `run_bottom_audit.sh` (§6a).
> The table below documents the gates as they were built; treat the oracle rows as history.

Beyond the self-tests, byte-gates cover each subsystem as it was built. They live in `tools/` as
`run_*_gate*.sh`; representative ones:

| Group | Example gate | Mesh / ranks | What it asserts |
|---|---|---|---|
| Kernels | `run_pressure_gate.sh`, `run_step_gate.sh`, `run_forcing_gate.sh` | pi, 1 | EOS/pressure, the whole assembled step, forcing read |
| CORE2 setup | `run_geom_gate_core2.sh`, `run_ic_gate_core2.sh` | CORE2, 1 | mesh geometry; 3-D initial conditions |
| Lifecycle | `run_lifecycle_gate_core2.sh`, `run_lifecycle_fullynative_gate_core2.sh`, `…_kpp_…`, `…_zstar_…` | CORE2, 1 | multi-step coupled runs per physics configuration |
| Multi-rank | `run_step_gate_multirank.sh`, `run_lifecycle_fullynative_gate_multirank.sh`, `run_lifecycle_jra55_gate_multirank.sh` | pi/CORE2, 2–864 | the same fields match per gid across rank counts |
| Output | `run_output_gate.sh`, `run_meshdiag_gate.sh` | pi, 1/2/8 | Zarr stores round-trip; mesh diagnostics match FESOM2 |
| Restart | `run_restart_gate_core2.sh`, `run_restart_gate_multirank.sh`, `run_restartroundtrip.sh` | CORE2, 1/2/8 | resume reproduces the straight-through trajectory |

Each gate is `max|Δ| = 0` at its stated anchor. As built, the np>1 restart gate admitted the
multi-rank element floor of §1/§7; since the shared-element sync it is strict at every np (the floor
only with `FESOM3_SHARED_SYNC=0`). The HANDOFF log lists the exact field/record counts per gate.

### 6a. The bottom-at-vertices gates

With the oracle gates retired, two FESOM_F-internal checks carry the numerical validation of the
vertex-bottom scheme.

**`tools/run_conserve_pi.sh` — conservation and no-leakage.** Runs `fesom_conserve` (an unforced
20-step pi run) at `zstar` np 1/2, `zstar`+Redi np 1, and `linfs` np 1, checking four per-step
invariants: total heat/salt content, `UV(:,nz,e) == 0` for
`nz >= nlevels(e)` (no velocity in a partly-land prism), `helem == mean(hnode over elnodes)` across
the element's full range, and finiteness of `UV`/`w`/`hnode`/tracers. Registered as ctest at np 1
and 2.

The conservation part gates on **zstar only**, at a relative tolerance of `1e-12`:

| mode | np | heat drift | salt drift |
|---|---|---|---|
| zstar | 1 | `-5.97e-15` | `-1.04e-15` |
| zstar | 2 | `+1.34e-14` | `-1.43e-14` |
| zstar + Redi | 1 | `-1.07e-15` | `-6.91e-16` |
| linfs | 1 | `-2.17e-04` | `-1.09e-06` |

**The Redi case is not optional coverage.** `diff_ver_part_redi_expl` is the only caller of the
`tr_xynodes` node-average whose denominator the depth-independent area changed, and it runs only
under `if (Redi)` (`oce_ale_tracer.F90:225`), which defaults to `.false.`. Without `FESOM3_REDI` that
code path never executes.

**And conservation alone does not validate that denominator.** The Redi tendency is a telescoping
flux divergence — `del_ttf += (vd_flux(nz)-vd_flux(nz+1))*dt/areasvol` — so the column total is
conserved whatever `tr_xynodes` contains. The denominator is pinned separately by `test_bottom`'s
*Redi: node-averaged gradient is exact on the WET area*: fed a constant element gradient (what a
linear tracer field produces), averaging over the elements that actually contribute returns it to
`4.4e-16`, while dividing by `areasvol` is off by up to `8.6e-01` — 86 % at the worst pi node.

**`linfs` is not tracer-conserving and cannot be gated on drift.** The linear free surface freezes
`hnode`, so the surface vertical advective flux `-w*T*area` at `nzmin` (`oce_adv_tra_ver.F90:66`) is
a real source/sink with no thickness change to balance it — the textbook reason `zlevel`/`zstar`
exist. Its drift is also non-monotone, i.e. a free-surface adjustment transient. So `linfs` runs for
the invariants only, with drift reported and never gated. `zstar` is the production mode, so the
gate covers the path that matters.

**`tools/run_bottom_audit.sh` — mesh delta.** Dumps the built geometry on pi and core2 and checks
the element-bottom shift against numbers measured from the mesh files *before* any code was
written, so the code answers back to the prediction rather than defining it:

| mesh | elements changed | volume | stagnant bottom cells |
|---|---|---|---|
| pi | 651 / 5839 | `+1.55 %` | 0 |
| core2 | 10465 / 244659 | `+0.38 %` | 0 |

Zero stagnant bottom cells is a hard requirement, not a quality metric: three unguarded divides
(`oce_ale.F90:88`, `oce_ale.F90:377`, `oce_pressure_bv.F90:310`) produce NaN without it, which is
why the mesh setup also asserts it at runtime.

### 6b. The hybrid-mesh gates

With every kernel on 4-slot formulas, triangle results changed at roundoff, so bit-identity with the
pre-hybrid code is replaced by tiered roundoff gates (rationale: `docs/HANDOFF.md` → "Gate migration
for hybrid meshes"):

- **`tools/run_tri_regression.sh`** — runs a FROZEN reference build (the last bit-identical state)
  and the current build side by side on pi and CORE2 for a few steps (SSH solved to `1e-14`) and
  compares per field in tiers PRE / EARLY / OCEAN / ICE / MIX / ICE_EARLY with per-mesh tolerances (the
  table is in the script header); `decomp:` cases compare np 1 vs np N restarted from one checkpoint,
  `conserve:` cases run `fesom_conserve`.
- **`tools/run_cutout_short.sh`** — 5-day production (`prod`, `prodpp`) and conservation-with-ice
  (`cons`) runs on the quad and mixed NEMO cut-outs, np 1 vs 2/8 from one checkpoint, plus the
  checkerboard time series (`tools/cutout_diag.py`).
- **`tools/run_ini_compare.sh`** — the ported `fesom_ini` against FESOM2's own `mesh_part` build on
  pi and CORE2 (byte-identical outputs).
- **`tools/run_ini_mixed.sh`** — `fesom_ini` on the NEMO quad/mixed meshes, every output checked by
  `tools/mesh/check_dist.py`, plus error cases (dead ends, flipped and bow-tie elements, corrupted
  partitions).
- The Python mesh tools have their own suites: `tools/mesh/test_nemo_meshq.sh`,
  `tools/mesh/test_mesh_tools.sh`.

### 7. Multi-rank and partition independence

Correctness across ranks is tested in two complementary ways:

- **Per-gid byte-gates** at `dist_2`, `dist_8`, and on up to `dist_864` (CORE2). Because dumps are
  keyed on global id, the multi-rank result is compared against the single-rank (and the oracle)
  value for the same entity.
- **Partition independence**: the output gate's `zarr_diff.py --output-cmp` asserts that the Zarr
  stores written at np=1, np=2, and np=8 are byte-identical (`max|Δ| = 0`) for every variable. The
  restart path adds a cross-rank check: a checkpoint written at one rank count restores exactly at
  another (`C2 ≡ C8`).

Shared (multiply-owned) border elements (§1) are kept equal by the shared-element sync: every step,
`update_vel` gives each of them the designated owner's `uv`, the copy the stores keep, and guards pin
the other per-element fields (`test_vinv` V13, `test_ice_dyn_mixed` D6). The multi-rank restart gate
therefore asserts strict `max|Δ| = 0` at every np. Before the sync (and still with
`FESOM3_SHARED_SYNC=0`) the copies differed by ≤1 ULP and the gate admitted a relative ε floor of
`1e-12` at np>1. The global divergence sum is tested at np>1 too: `|Σ ssh_rhs_old| / Σ|ssh_rhs_old|`
stays at roundoff (`test_shared_elem` S3, `fesom_conserve` with `FESOM3_DIVSUM_TOL`,
`tools/run_conserve_pi.sh`, and the `divsum=` field of `FESOM3_STEP_DIAG`).

### 8. Output and restart validation

**Output (Zarr).** `tools/zarr_diff.py` validates stores several ways: `--roundtrip` (values match a
known generator formula), `--lz4` (compressed read-back), `--meshdiag` (the mesh-diagnostics store
matches the FESOM2 oracle), `--output`/`--output-cmp` (per-variable correctness and
partition-independence), and `--frame` (vector rotation to geographic vs. native). The output gate
also sweeps the writer knobs (chunk sizes, compression, writer count) and asserts every variant is
value-identical to the default.

**Restart.** The restart gate runs the model straight through N steps, and separately runs it K
steps, checkpoints, starts a *fresh process*, reads the checkpoint, and runs the remaining N−K steps
— then compares the full state (including sea ice and the EVP stress tensor `sigma`) with
`dump_diff.py --ignore-step`. Results:

- **Same rank count:** byte-exact resume — `max|Δ| = 0` for the entire state (FESOM_F serializes
  `sigma`, which FESOM2's own restart omits, so even the ice stress matches), at np=1 and, since the
  shared-element sync, at np>1 (CORE2 np 2, mid-run and month-boundary splits, default and
  vector-invariant physics). Only with `FESOM3_SHARED_SYNC=0` does the ≤1-ULP element floor of §1/§7
  apply at np>1.
- **Different rank count:** the state is *restored* exactly (partition-independent storage,
  `C2 ≡ C8`), after which the run continues on a physically valid — but, because element ownership
  differs, not bit-identical — trajectory.

### 9. Results

**Milestones (all byte-exact, `max|Δ| = 0`, against FESOM2 v2.7.3 at the stated anchor):**

| Milestone | Tag | Subsystem |
|---|---|---|
| M0 | `m0` | foundation: types, partitioning, halo, dump, mesh |
| M1 | `m1` | tracer advection (FCT) |
| M2 | `m2-mvp` | dynamical core + multi-rank whole-step |
| M3 | `m3` | sea ice EVP/mEVP + air–sea coupling |
| M4 | `m4` | GM/Redi |
| M5 | `m5` | KPP vertical mixing + shortwave penetration |
| M6 | `m6` | zstar vertical coordinate |
| M7 | `m7` | TKE vertical mixing |
| M8 | `m8` | multi-year production runs (CORE2, JRA55-do) |
| M9 | `m9` | Zarr output |
| M10 | `m10` | restart / checkpoint |

Each milestone is byte-exact on 1 rank and multi-rank, for both EVP variants.

**Full-physics headline run (M8).** The full default configuration (zstar + TKE + GM/Redi +
shortwave penetration + native sea ice, JRA55-do forcing from 1958) on CORE2 at **dist_864** (864
ranks) reproduces FESOM2 bit-for-bit over a **full model year — 1,138,800 compared records**,
spanning ~17,280 steps (3.4× past the day-107 mark where the flush-to-zero issue, §1/LESSONS L51,
once appeared). A free-running **two-year run (35,040 steps)** with no oracle is physically stable —
no NaN/Inf, bounded global peaks, a correct seasonal ice cycle — confirming the model is not merely
byte-faithful but well-behaved.

**Compiler / precision matrix.**

| Build | Status |
|---|---|
| Intel + double | **anchor** — byte-exact vs FESOM2; all gate results |
| GNU + double | portability — all 21 `ctest` cases pass; pure-data fields bit-identical to Intel (only rotation-derived coordinates differ at ~1e-13, expected libm differences) |
| single precision | builds (scaffolding, §IMPLEMENTATION 6); not a byte-gate target |

### 10. Scalability and performance

Performance is measured, not estimated, with a **two-point method**: run from cold start for 100 and
for 600 steps and take per-step cost as `(T₆₀₀ − T₁₀₀)/500`, which removes fixed setup cost. FESOM_F
reports per-component timings directly (`src/infra/mod_timer.F90`, in ms/step); FESOM2 reports
seconds-per-run, so a fair comparison must divide by the step count first.

The one component that was initially slower than FESOM2 was **surface forcing**: reading the
DEFLATE-compressed JRA55-do netCDF re-opened the file and read both time brackets on every rank at
every crossing — 15.74 ms/step vs FESOM2's 7.03 (2.24×, +8.7 ms/step). Switching to a persistent file
handle with a double buffer (`FESOM3_FORCING_PERSIST`, default on; `src/forcing/mod_forcing_read.F90`)
fixed it. Measured at **dist_512** (steady state, two-point):

| Component | FESOM2 | FESOM_F (before) | FESOM_F (after) |
|---|---:|---:|---:|
| forcing | 7.03 | 17.71 | **5.78** (3.07× faster than before; now below FESOM2) |
| ocean step | 30.24 | — | **28.66** (0.95×, FESOM_F faster) |
| **loop total** | 47.81 | 58.26 | **44.49** |

So after the fix FESOM_F's whole timestep at 512 ranks (44.49 ms/step) is below FESOM2's (47.81) —
its all-ranks forcing read beats FESOM2's rank-0-read-plus-broadcast on this machine. The full
investigation, including the ruled-out hypotheses, is in
`docs/plans/2026-06-28-forcing-perf-investigation.md`. (Systematic scaling/throughput work remains
future — see [IMPLEMENTATION.md](IMPLEMENTATION.md) §12.)

### 11. Reproducing the tests

On a Levante login node:

```bash
./configure.sh --compiler intel --precision dp --clean --build   # build the anchor
cd build_intel_dp && ctest -j 8 --output-on-failure              # 153 self-tests
cd .. && bash tools/run_conserve_pi.sh                           # conservation gate → GATE OK
bash tools/run_tri_regression.sh                                 # triangle roundoff gate (§6b)
```

(`ulimit -s unlimited` first: the CORE2 ctest cases need it.) The FESOM2 byte-gates such as
`tools/run_step_gate.sh` are retired (§6).

Multi-rank and CORE2 gates (and the dist_512/dist_864 SLURM jobs) are the `tools/run_*_gate*.sh` and
`tools/*.sbatch` scripts. The development log `docs/HANDOFF.md` records the exact command and
expected record count for each.

---

*See also:* [README.md](README.md) (build & run) · [IMPLEMENTATION.md](IMPLEMENTATION.md)
(architecture) · `docs/HANDOFF.md` (development log) · `docs/LESSONS.md` (reproducibility traps) ·
`docs/plans/`.
