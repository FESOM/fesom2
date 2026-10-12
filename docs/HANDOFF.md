# FESOM3 — Handoff (durable state across sessions)

Single source of truth for "where are we / what's next". Update at the end of every task.
Full pre-M2.12 milestone detail + the per-gate recipes live in [`HANDOFF-archive.md`](HANDOFF-archive.md).

## Where we are

- **Shared-element velocity sync ✅ DONE (2026-10-12, plan `plans/completed/20261012-shared-element-sync.md`,
  on `bottom-at-vertices`):** a border element is owned by every rank that owns one of its vertices;
  `update_vel` now gives every such element one owner for `uv` (the highest vertex-owner rank, whose copy
  the restart writer keeps), every step (`FESOM3_SHARED_SYNC`, default 1). Fixes the linfs volume and
  salt leak of the CORE2 runs results7/11/17/21 and makes np>1 restarts exact. np 1 is bit-identical,
  np>1 changes at roundoff. Next for the user: rebuild the production binary and resume results17 or
  results21 with `FESOM3_STEP_DIAG=1` → section "Shared elements: one owner for `uv`" → "Production check".
- **Hybrid (mixed triangle + quad) meshes ✅ DONE (2026-10-09, plan `plans/completed/2026-10-08-hybrid-meshes.md`,
  commits `d4c6f8d`…`d07f44f` + the docs commit, on `bottom-at-vertices`, unpushed):** model, `fesom_ini`
  (FESOM2's `mesh_part`, ported) and the Python mesh tools run triangle, quad and mixed meshes (NEMO
  ORCA025 → `nemo025q_fesom`, `nemo025mix`); unified 4-slot kernels, so triangle runs are no longer
  bit-identical (roundoff gate `tools/run_tri_regression.sh`); ctest 153/153 (after the review fixes); full NEMO 0.25° 2-day runs
  at np 8 finite. Next for the user: rebuild `build_intel_dp` from this branch, first full-size NEMO run
  (`work/job_levante_nemo025`), the open design questions → section "Hybrid (mixed triangle + quad)
  meshes" → "Known findings / open questions for the user".
- **Milestone:** M2 (minimal ocean dynamical core) — **multi-rank MVP byte-match ACHIEVED (tag `m2-mvp`).**
  **→ M3 (sea ice EVP) ✅ COMPLETE 2026-06-23 (all M3a–M3f DONE, `max|Δ|=0` 1-rank AND multi-rank) — ✅ COMMITTED
  (`b591153`, tag `m3`).**
  **→ M4 (GM/Redi) ✅ COMPLETE 2026-06-24 (all M4a–M4f DONE, `max|Δ|=0` 1-rank AND multi-rank) — ✅ COMMITTED (`7a088de`, tag `m4`).**
  M4a (producers) + M4b (GM diffusivity/streamfunction/bolus-velocity) + M4c (GM bolus into advection) + M4d (Redi
  isopycnal diffusion) + M4e (GM+Redi in the FORCED/fully-native lifecycle) + M4f (multi-rank) all byte-gated
  `max|Δ|=0` vs FESOM2 (CORE2 1-rank AND dist_2/dist_8, both whichEVP). Plan: `docs/plans/2026-06-23-m4-gm-redi.md`.
  **→ M5 (KPP vertical mixing + sw_pene) ✅ COMPLETE 2026-06-24 (all M5a–M5d DONE, `max|Δ|=0`
  1-rank AND multi-rank) — tag `m5`.** Scoped + decomposed M5a–M5d (plan: `docs/plans/2026-06-24-m5-kpp.md`).
  **M5a (isolated KPP module gate, 95 pressure-gate fields) ✅ DONE; M5b (wire KPP into `step_oce` +
  UNFORCED CORE2 lifecycle) ✅ DONE: `max|Δ|=0` KPP-alone 195+325 AND KPP+GM+Redi 195.**
  M5a-1 (init + wscale) + M5a-2 (ri_iwmix) + M5a-3 (prestep dVsq/ustar/Bo + dbsfc + bldepth) +
  M5a-4 (blmix_kpp + enhance + combine + viscAE average) all byte-exact on BOTH pi AND CORE2, all
  FOUR first try; M5b assembled the real `oce_mixing_kpp_driver` + the `mix_scheme_nmb==1` dispatch
  (one new physics step beyond M5a: `smooth_blmc` — LESSONS **L45**).
  **M5c (forced/fully-native lifecycle + sw_pene) ✅ DONE (first try, LESSONS L46): the PRODUCTION
  CORE2 column physics (KPP + GM + Redi + shortwave penetration) byte-exact in the fully-native
  forced lifecycle — `max|Δ|=0` at 195 (3-step) AND 325 (5-step), BOTH whichEVP**
  (`tools/run_lifecycle_kpp_native_gate_core2.sh`). KEY FINDING: of the two tracer-TDMA terms the
  plan named, only `sw_3d` is LIVE in work_core — the `ghats` nonlocal flux is gated off by
  `use_kpp_nonlclflx=.false.` (absent from work_core, like MLD1_ind/K_hor in M4); ghats was
  transcribed anyway + byte-gated by a `KPP_NONLCL=1` variant (`max|Δ|=0`, non-vacuous 7.3e-5).
  `use_sw_pene` consolidated to its single `mod_config` home (a latent duplicate with `mod_param_phys`).
  **M5d (multi-rank) ✅ DONE (first try, the M4f lesson): the production KPP+sw_pene+GM+Redi native
  lifecycle byte-matches FESOM2 `max|Δ|=0` at MULTI-RANK — CORE2 dist_2 + dist_8, BOTH whichEVP,
  195 (3-step) + 325 (5-step)** (`tools/run_lifecycle_kpp_native_gate_multirank.sh`). PURE WIRING:
  the KPP module was optional-`partit` from the start (`oce_mixing_kpp_driver` does owned-loops +
  `exchange_nod(blmc/Kv_double/ghats/viscA)` + `smooth_blmc`); the only NEW work was the
  `fesom_lifecycle_native_mr` KPP/sw_pene block (LOCAL-sized) + `cal_shortwave_rad`'s optional
  `partit` (loop owned+halo `nNodL`, since the local mesh's `mesh%nod2D` = the GLOBAL count). NO
  regression (ctest 13/13 + 1-rank production + pi M2.10c forcing-gate `sw_3d` all `max|Δ|=0`).
  **→ M6 IN PROGRESS (2026-06-24): M6a = zstar ALE.** Plan: `docs/plans/2026-06-24-m6-zstar.md`
  (decomposed M6a-1..4). **DIRECTION (user 2026-06-24):** target is `which_ALE='zstar'` (the canonical
  `config/namelist.config` default — the `work_core` copy I'd read was customized to zlevel; my brief
  zlevel/"paper-parity" detour was WRONG). **DROP "paper-parity" framing** — the goal is a faithful 1:1
  port of FESOM2 with BETTER architecture; the paper is only a conversion-methodology example. **NO
  aEVP** (mEVP/whichEVP=1 already done in M3). **TKE → its own plan** (M7), ported from the FORTRAN
  oracle but with C-port experience (NO CVMix, NO diagnostics-array bloat). **Production runs (months/
  years, compare to Fortran ~byte-identical) → its own plan** (M8).
  **M6a-1 (Shchepetkin full-free-surface PGF `pressure_force_4_zxxxx_shchepetkin`) ✅ DONE** — byte-exact
  `max|Δ|=0` (CORE2 pressure gate `pgf_x_shchep`/`pgf_y_shchep`, ~3e-5; shared by zlevel/zstar). En
  route, fixed a **pre-existing m5 BASELINE CRASH**: the pressure shim left `use_sw_pene=.true.` into the
  tracer solve → M5c `sw_3d` term read UNALLOCATED `dyn%work%sw_3d` (Release tolerated it until a heap
  re-layout → SIGSEGV; debug `-check all` pinned it). FIX mirrors the oracle (`use_sw_pene=.false.`
  before the tracer-solve loop). LESSON **L47**.
  **M6a-2 (zstar thickness machinery, UNFORCED lifecycle) ✅ DONE (first try, 2026-06-24):** `max|Δ|=0`
  vs FESOM2, **195 (3-step) + 325 (5-step)** records — SSH/hnode now EVOLVE (max|eta_n| 0.35→0.68), so
  the whole zstar free surface is exercised: the Shchepetkin PGF + the NEW `update_stiff_mat_ale`
  (per-step SSH-stiffness 2nd-term update via the lagged `dhe`, accumulating to track the surface) +
  the `vert_vel_ale` zstar stretch (`dd=dhbar/H` proportional over all full-prism layers) + the
  `update_thickness_ale` commit (`hnode_new→hnode`, recompute `zbar_3d_n`/`Z_3d_n` bottom-up).
  `tools/run_lifecycle_zstar_gate_core2.sh` (oracle `WHICH_ALE=zstar` env + FESOM3_WHICH_ALE=zstar; NO
  oracle SOURCE change — zstar is native FESOM2). zstar is SIMPLER than zlevel (no local-zstar fallback/
  min_hnode/CFL_z; all scalar divides → no L29 trap). No regression: linfs lifecycle 195 + step-65 +
  ctest 13/13 all `max|Δ|=0`.
  **M6a-3 (forced/native zstar + real freshwater flux) ✅ DONE (2026-06-24): the FULL PRODUCTION zstar
  config byte-matches FESOM2 `max|Δ|=0` — zstar + KPP + GM + Redi + sw_pene + native ice/forcing,
  BOTH whichEVP, 195 (3-step) + 325 (5-step); reduced (PP) variants too.** `use_virt_salt=.false.`
  (`is_nonlinfs=1`) for zstar: TWO water_flux terms went live — `compute_ssh_rhs_ale`
  (`ssh_rhs-=alpha*water_flux*areasvol`, FESOM2 `:2122`) AND `compute_hbar_ale`
  (`ssh_rhs_old-=water_flux*areasvol`, `:2262`) — plus `vert_vel_ale` `Wvel(1)-=water_flux` and the
  real-salt path. KEY WIRING BUG fixed: `obudget` writes `atm%real_salt_flux` (ice-brine `rsf=fwice*Sice`,
  the `.not. use_virt_salt` branch — M3 pre-wrote it) but the native driver passed a LOCAL zeroed
  `real_salt_flux` to `step_oce` → S diverged ~2e-3; fix = pass `atm%real_salt_flux`. bc_surface /
  oce_fluxes / obudget all ALREADY had the `is_nonlinfs`/`use_virt_salt` branches (M3 anticipated zstar),
  so M6a-3 was the 2 water_flux terms + the real_salt_flux wiring + the driver flags.
  `tools/run_lifecycle_zstar_native_gate_core2.sh` (+ WHICH_ALE env in the forced oracle runner +
  the fullynative gate). No oracle SOURCE change. No regression: linfs native (reduced + production) +
  linfs/zstar unforced + step-65 + ctest 13/13 all `max|Δ|=0`.
  **M6a-4 (multi-rank zstar) ✅ DONE (first try, 2026-06-24) → M6a (zstar ALE) COMPLETE.** The full
  production zstar config byte-matches FESOM2 `max|Δ|=0` at MULTI-RANK: CORE2 dist_2 + dist_8,
  reduced (PP) AND production (KPP+GM+Redi+sw_pene), BOTH whichEVP, 3-step. PURE WIRING (the M4f/M5d
  lesson): the zstar kernels were optional-`partit` from the start (update_stiff_mat_ale owned-edge +
  vert_vel_ale/update_thickness_ale owned loops + the exchanges in place; `dhe` read at owned el(i)
  only — M2.12c-2 invariant, no halo extension); the only new work was the `fesom_lifecycle_native_mr`
  zstar block (FESOM3_WHICH_ALE + use_virt_salt/is_nonlinfs + `atm%real_salt_flux`, mirror of M6a-3) +
  the MR gate (`run_lifecycle_zstar_native_gate_multirank.sh`). No regression: linfs MR native dist_2 +
  ctest 13/13. **RESUME at M7 = TKE** (its own plan — port from the FORTRAN oracle with C-port
  experience: NO CVMix, NO diagnostics-array bloat; `mix_scheme='TKE'`, prognostic `tke`, `work_tke_dump`/
  `work_linfs_tke` exist; decompose like KPP). Then M8 = production months/years runs (netCDF output +
  restart write + run-length-in-years + Fortran comparison, ~byte-identical). **Tag `m6` at the M6a
  commit.** The WHOLE default FESOM config (zstar + KPP + GM + Redi + sw_pene, native ice/forcing) is
  now byte-exact 1-rank AND multi-rank, both whichEVP.
  **→ M7 (TKE vertical mixing — Gaspar'90 + Blanke-Delecluse mixing length, library-free) ✅ COMPLETE
  2026-06-25 (all M7a–M7e DONE, `max|Δ|=0` 1-rank AND multi-rank, both whichEVP) — tag `m7`.** Plan:
  `docs/plans/2026-06-24-m7-tke.md`. `mix_scheme='cvmix_TKE'` (mix_scheme_nmb==5) → a prognostic-TKE
  Av/Kv producer transcribed DIRECTLY from the FESOM2 CVMix files (`cvmix_tke.F90` integrate_tke +
  `cvmix_utils_addon.F90` solve_tridiag + `gen_modules_cvmix_tke.F90` calc_cvmix_tke) into a clean
  `src/oce/oce_mixing_tke.F90` — NO CVMix library, NO tke_type put/get, NO diagnostics-output bloat
  (the C-port interpretation). `tke_only=.true.` (IDEMIX/Langmuir/Dirichlet branches all dead);
  `tke_mxl_choice==2` (the only ported mixing length). The `&param_tke` DOUBLES are baked in
  `mod_param_phys` (⚠️ `tke_cd=3.75`, the namelist-over-codedefault trap; used at the Neumann surface BC).
  **M7a-1 (isolated column CONTROLLED REPLAY):** `fesom_tkereplay` ingests the oracle `tke_dump` per-step
  INPUTS → integrate_tke → emits OUTPUTS, byte-identical to the oracle (13 tags × 3 steps `max|Δ|=0`).
  The literals transcribed VERBATIM (`6.6`, `**(3./2.)`, `sqrt(2D0)`, reciprocal-multiply Thomas) are
  byte-exact under the `-r8` anchor — the C-port's 6.6/pow/fmax traps are ABSENT in Fortran↔Fortran.
  **M7b (unforced linfs+TKE lifecycle):** the prognostic `tke` recurrence (`dyn%work%tke` step→step) +
  the Av/Kv producer + mo_convect — `max|Δ|=0` 195+325. **M7c (forced/fully-native, both whichEVP):** the
  surface WIND term (`forc_tke_surf=|stress_node_surf|/density_0` → Neumann `(cd*forc**(3./2.))/dzt(1)`)
  + sw_pene + GM + Redi — `max|Δ|=0` 195+325 BOTH whichEVP. **M7d (multi-rank dist_2/dist_8):** PURE WIRING
  (the M4f/M5d lesson) — `calc_cvmix_tke` is optional-`partit` from the start (owned-loop +
  `exchange_nod(tke_Kv)`/`exchange_nod(tke_Av)` BEFORE the node→elem 3-vertex average; `tke` is NEVER
  exchanged — the recurrence is partition-local) — `max|Δ|=0`. **M7e (zstar+TKE production):** the FULL
  default config with TKE — `max|Δ|=0` 1-rank AND multi-rank, both whichEVP (zstar `hnode`/`Z_3d_n` are
  time-varying, read live by `calc_cvmix_tke` — a pure config flip; TKE does not touch the ALE thickness
  machinery). NET consumer-side change vs M6 was ZERO (Av element / Kv node UNCHANGED; no ghats; `sw_3d`
  already M5c) — markedly simpler than M5. No regression: KPP / zstar+KPP / ctest 13/13 all `max|Δ|=0`.
  **`tke` restart serialization = M8** (the FIRST stateful mixing field — `dyn%work%tke`; not in
  `write_t_dyn_work`). **Tag `m7`. RESUME at M8** = production months/years runs (netCDF output + restart
  write + run-length-in-years + Fortran comparison ~byte-identical) — its own plan.
  **→ M8 (production long simulations — months/years, CORE2/JRA55 multinode) ✅ COMPLETE 2026-06-27
  (tag `m8`).** Plan: `docs/plans/2026-06-25-m8-long-simulations.md`. M8a (clock) + M8b (forcing
  record/day rollover, `forcing_sbc_do`) + M8c Steps 1+2 (JRA55-do-v1.4.0 production forcing-read gate +
  month-climatology SSS/Sweeney-chl + year-file rollover) all `max|Δ|=0` (the short-gate byte bar holds
  by induction). **M8c Step 3 (2-year headline, JRA55 1958-1959, dist_864) BYTE-EXACT `max|Δ|=0`:** an
  earlier day-107 divergence was root-caused to a **runtime FTZ (flush-to-zero) mismatch** — the FESOM3
  process ran with the MXCSR FTZ bit OFF (kept a denormal `m_snow` ~1e-309 m) while FESOM2 had it ON
  (flushed to 0.0), flipping the `if(hsn>0)` ice-albedo branch. FIXED in `fesom_lifecycle_native_mr.F90`
  via `ieee_set_underflow_mode(gradual=.false.)` (startup + per-step; FESOM3-only, no FESOM2 change; **L51**).
  Confirmed byte-exact over a full model year (1,138,800 records, steps 1→~17280, 3.4× past day-107) +
  the 5100-step FTZ gate + the separately-gated year-rollover. (Full 2-yr clean re-run pending only because
  the heavyweight oracle leg crashed on its own output at the year boundary — not an F3 divergence.) M8f:
  full no-regression green (ctest 13/13 Intel+GNU + production MR gate both whichEVP `max|Δ|=0`) after
  backing out the fw_ bisect scaffolding from both codes. **netCDF output + restart (incl. `tke`
  serialization) = M9.**
  **→ PERF side-investigation (2026-06-28, ✅ RESOLVED): F3 was at parity/faster everywhere EXCEPT
  surface forcing (`2.24×` F2, `+8.7 ms/step`). Root cause (measured via `mod_timer` sub-timers +
  microbench): the JRA55-do files are DEFLATE-compressed netCDF-4/HDF5 (one slice/chunk), and the legacy
  `forcing_getcoeffld` did `nc_open` (HDF5 metadata re-parse) + read BOTH brackets (~800 KB inflate each) +
  `nc_close` every crossing on every rank. FIX = FESOM2's persistent handle + double-buffer (reopen only on
  year change; reuse the cached `t_indx` slice, read only the new `t_indx_p1`) — `src/forcing/mod_forcing_read.F90`,
  toggle `FESOM3_FORCING_PERSIST` (default ON). dist_512: forcing `17.7→5.78 ms/step` (3.07×), open→0, read
  halved; **5.78 < F2's 7.03 ⇒ F3 now beats F2 on forcing AND overall** (LOOP `58.3→44.5`). Byte-exact
  (`…jra55_gate_multirank.sh 8 {12,48}` ⇒ `|Δ|=0`; `ctest` 13/13). Why all-ranks-read beats F2's
  rank-0+`Bcast`: the Levante broadcast is slow (two-copy vader, KNEM off) — confirms the "do NOT broadcast"
  call. Full writeup: `docs/plans/2026-06-28-forcing-perf-investigation.md` (RESOLVED section at top).**
  **M8a (clock + run-length driver) ✅ COMPLETE 2026-06-25 (first try):** ported `src/infra/mod_clock.F90`
  (`g_clock` VERBATIM — `clock`/`clock_init`/`check_fleapyr`/`is_fleapyr` + `clock_nsteps`=`get_run_steps`;
  `r_restart` rehomed to `mod_clock`; `clock_finish`/`clock_newyear`/`use_transit` deferred to M9) and wired
  it into `fesom_lifecycle_native_mr`: `step_per_day=48`/`mod_config%dt=1800` so `clock` advances per step;
  `FESOM3_START_CLOCK` ("t d y", default "0 1 1948") → rank-0 writes a 2-identical-line `.clock`, barrier,
  all ranks `clock_init`; `nsteps=clock_nsteps(partit)` (`FESOM3_NSTEPS` short-gate override kept); `call
  clock` at the top of each step. Forcing path UNCHANGED (still hard-codes the day-1 rdate) ⇒ byte-NEUTRAL
  within day 1. Gates `max|Δ|=0`: unit `fesom_clocktest` (`tools/run_clocktest.sh` — hand-computed table for
  day/month/year rollovers + `clock_nsteps` s/d/m/y, Intel AND GNU, 1+2 rank); zstar+TKE MR 325 records
  dist_2 (BOTH whichEVP) + dist_8, 5 steps; no-regression ctest 13/13 (Intel+GNU dp) + step-65 1-rank.
  **M8b (forcing record/day rollover — `forcing_sbc_do`) ✅ COMPLETE 2026-06-26:** ported
  `gen_surface_forcing.F90:1524-1567` as `forcing_sbc_do` (per-field running rdate from `mod_clock` with
  the `-dt/2` half-step; crossing test `rdate>nc_time(t_indx_p1) .and. nc_time(t_indx)<nc_time(ntime)` ⇒
  re-fire `forcing_getcoeffld` + wind g2r re-rotate) MINUS the year branch (M8c) and the leap special-case
  (dead for `'none'`). `forcing_getcoeffld` now PERSISTS the per-field bracket `t_indx`/`t_indx_p1` in
  `t_ffile` (`nc_Ntime`==existing `ntime`); driver cold-start coef build seeded from the clock
  (`rdate_cold` NO half-step, SSS at `i=month`) — required REORDERING the M8a clock block to BEFORE the
  forcing setup. **GATES `max|Δ|=0` (both whichEVP, CORE2 full-year pool):** login `dist_2`+`dist_8` 48
  steps = 3120 records with **18 `getcoeffld` re-fires** (6-hourly winds/q/Tair 1/2→2/3→3/4→4/5 at steps
  ~7/19/31/43 + daily radiation 1/2 at ~25); **PRODUCTION `dist_864` (7 nodes/864 ranks) — the FIRST
  fesom3 SLURM-batch inter-node byte-gate — 3120 records, both whichEVP** (also dist_128 batch). No
  regression: ctest 13/13 + zstar+TKE MR stub day-1 (325 records, dist_2 both whichEVP + dist_8). NEW
  `tools/run_lifecycle_long_gate_multirank.sh` + `…_dist864.sbatch`; `run_lifecycle_forced_core2.sh` + the
  MR gate parametrized with `FORCING_OVERRIDE`/`START_CLOCK` + `srun`/204800-stack SLURM branches. TWO
  latent harness bugs fixed en route — see **L48** (`grep|head` SIGPIPE + `set -e`; `ulimit -s unlimited`
  raise on a hard-capped compute node).
  **M8c (JRA55-do production forcing — the REAL target, not the CORE2 substitute) ✅ Steps 1+2 COMPLETE 2026-06-26:**
  The production atmosphere is **JRA55-do-v1.4.0** (3-hourly, `gregorian` + `include_fleapyear=.true.`, `tmid=0`,
  fields `uas/vas/huss/rsds/rlds/tas/prra/prsn`, start **1958** — NOT CORE2/1948; M2.10→M8b used the CORE2 NCAR
  substitute). **Step 1 (forcing-read gate)** byte-matches the 8 atmosphere fields + bulk from 1958 — exposed and
  fixed the **calendar NUL-byte trap** (`nc_get_att_text` returned `'gregorian'//char(0)` ⇒ `trim()=='gregorian'`
  silently FALSE ⇒ noleap branch ⇒ wrong slice ⇒ `ssh_rhs`~1e4; CORE2's `'NOLEAP'` had masked it for 8 milestones —
  **L49**). `max|Δ|=0`: `dist_2` (both whichEVP) + `dist_8`; CORE2 no-regression 48-step `dist_2` (3120 records) still
  `max|Δ|=0`. **Step 2 (rollover)** ports the **year-file branch** (`forcing_sbc_do` `yearnew/=yearold` ⇒
  `forcing_read_grid` per field + `force_newcoeff`; grid re-entrancy via `if(allocated)deallocate`) + **SSS monthly
  read-ahead** + **Sweeney chl monthly** (`roll_monthly_clim` in the driver: `update_monthly_flag` last-instant trigger,
  `i=month/(mstep>1→+1)/wrap`; chl needed `nc_get_att_dp` optional-`stat` + `miss=-99` for the attribute-less Sweeney
  file — **L50**). `max|Δ|=0` at `dist_2`: Jan→Feb month gate (SSS read-ahead fires step 2), Dec31'58→Jan1'59 YEAR
  gate (atmosphere 1958→1959 rollover + Dec→Jan SSS), chl-Sweeney month gate — all non-vacuous (stdout `YEAR ROLLOVER`
  / `roll_monthly_clim` markers). Files: `mod_io_netcdf.F90` (NUL sanitize + att `stat`), `mod_forcing_other.F90`
  (no-fill sentinel), `mod_forcing_read.F90` (year branch + re-entrancy), `fesom_lifecycle_native_mr.F90` (forc_set
  JRA55/CORE2 wiring, `roll_monthly_clim`, chl-Sweeney), runners (`run_lifecycle_jra55_gate_multirank.sh` +
  `FORCING_SET`/`CHL_SWEENEY` knobs). **NOTE: the M8 plan's "Goal & target config" + M8c task still say
  CORE2/1948/noleap/chl-`'None'` — STALE (pre-JRA55-correction); the `work_core` namelist (SSS `'CORE2'` monthly, chl
  `'Sweeney'` monthly) is ground truth.**
  **M8c Step 3 (2-year headline) ✅ COMPLETE 2026-06-27 — BYTE-EXACT `max|Δ|=0` (day-107 divergence root-caused + FIXED).**
  An earlier divergence at step 5095 (day ~107) was traced to a **runtime flush-to-zero (FTZ) mismatch between the two
  PROCESSES** — NOT a port bug. FESOM3 ran with the MXCSR FTZ bit OFF, keeping a **denormal `m_snow`** (~1.3e-309 m =
  physically-zero snow, node 119505) where FESOM2 (FTZ ON) flushed it to exactly 0.0; that flipped the ice-albedo branch
  `if(hsn>0)` (snow 0.85 vs ice 0.65) → `t_skin` diverges 2.6% → freshwater flux → global `net` (`integrate_nod_2D`) →
  SSH spill at all 126858 nodes. Confirmed emergent, not transcription: every snow/ice routine AND all *compile* flags
  are byte-identical (`m_ice`/`a_ice` stay `|Δ|=0` through 5094); the difference was the per-process FP runtime mode,
  probed via `ieee_get_underflow_mode` (F2 `gradual=.false.`, F3 `gradual=.true.`). **FIX:** FESOM3 forces FTZ on with
  `ieee_set_underflow_mode(gradual=.false.)` (startup after `par_init` + re-asserted per-step) in
  `fesom_lifecycle_native_mr.F90` — no FESOM2 change, physically correct (**L51**; bisected via
  `tools/onset_allnode.py` denormal-filtered + `extract_gid.py`; memory `m8c-day107-divergence`). **Byte-exact confirmed:**
  the 5100-step FTZ gate `max|Δ|=0` (331500 records, past day-107) AND the 2-yr gate overlap `max|Δ|=0` over **1,138,800
  records spanning steps 1→~17280 — a full model year, 3.4× past the day-107 flip** — plus the separately-gated year
  rollover (Step 2). (Full 2-yr clean re-run pending only because the heavyweight oracle leg crashed on its own NetCDF
  output at the year boundary — an oracle-side I/O crash, not an F3 divergence; lean re-run with oracle output disabled
  in flight.) Corroborating physics: the free-running F3 stability run
  (`tools/run_lifecycle_2yr_freerun_f3_dist864.sbatch`, job 25927003, GLOBAL `MPI_MAX` per-step diagnostics) completed
  **35040 steps, 0 NaN/Inf/FPE**, global peaks bounded with NO secular drift (`max|eta|` 1.8–2.0 m, `max|uv|` 1.5–2.9 m/s,
  `Tmax` 30–32.5 °C, `Smax` ~41.04 flat), correct ice **seasonal cycle** (`a_ice`→1.0, `m_ice` 2.0→6.1→4.6 m), and the
  `1958→1959` rollover + all 24 monthly SSS/chl read-aheads fired. Known follow-up: the year-rollover crossing-test
  dereferences a stale `t_indx_p1` — safe for 2920→2920 (1958→1959) but must reset before a leap-OUT crossing
  (2928→2920, e.g. 1960→1961) for runs past 1960.
  (M3 scoped 2026-06-22 into M3a–M3f.) M3a (ice foundation + cold-start IC + FCT mass
  matrix) ✅ DONE 2026-06-22 (`max|Δ|=0`, 4 fields, CORE2 1-rank, `tools/run_ice_gate_core2.sh`). M3b (ocean2ice + EVP
  dynamics) ✅ DONE 2026-06-22 (`max|Δ|=0`, 7 fields × BOTH whichEVP=0 standard-EVP AND whichEVP=1 mEVP, CORE2 1-rank,
  `tools/run_evp_gate_core2.sh`). M3c (ice FCT advection) ✅ DONE 2026-06-22 (`max|Δ|=0`, 6 fields = `rhs_a/m/ms` +
  post-advection `a_ice/m_ice/m_snow`, × BOTH whichEVP=0 AND whichEVP=1, CORE2 1-rank, `tools/run_icefct_gate_core2.sh`;
  first try). M3d (thermodynamics: `cut_off` + `thermodynamics` + `therm_ice`/`budget`/`obudget`/`flooding`) ✅ DONE
  2026-06-23 (`max|Δ|=0`, 6 fields = post-thermo `a_ice/m_ice/m_snow` + `flx_h`/`flx_fw` + `t_skin`, × BOTH whichEVP=0
  AND whichEVP=1, CORE2 1-rank, `tools/run_icethermo_gate_core2.sh`; first try, L38). M3e (`oce_fluxes` coupling-out —
  THE PAYOFF) ✅ DONE 2026-06-23 (`max|Δ|=0`, 5 fields = `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` +
  `stress_surf`, × BOTH whichEVP=0 AND whichEVP=1, CORE2 1-rank, `tools/run_iceflux_gate_core2.sh`; first try, L39 —
  the prescribed M2.11c-2 flux dump is now NATIVE). **M3f-1 (`ice_timestep` assembly) + M3f-2 (native-flux coupled
  multi-step forced lifecycle) ✅ DONE 2026-06-23 (first try, L40): `max|Δ|=0` on the 195-record (3-step) AND
  325-record (5-step) NODE substep set, BOTH whichEVP=0 AND whichEVP=1 — the native sea-ice + air-sea budget
  (`ocean2ice → ice_timestep → oce_fluxes_mom → oce_fluxes`) drives the multi-step coupled CORE2 ocean step
  byte-exactly; `tools/run_lifecycle_native_gate_core2.sh`. M3f-3a (native CORE2 NCAR forcing READ) + M3f-3c (native
  NCAR bulk Ch/Ce + stress_atmoce + wind-on-ice stress_atmice) ✅ DONE 2026-06-23 (L41) — 14/16 atmospheric arrays
  native, `max|Δ|=0` self-checked in the coupled lifecycle, both EVP variants. **M3f-3b (runoff + SSS monthly
  climatology — the last 2/16) ✅ DONE 2026-06-23 (first try, L42) → M3f-3 COMPLETE: the FULLY-NATIVE CORE2 lifecycle
  (NO prescribed atmosphere — the whole air-sea forcing computed in-driver) is `max|Δ|=0` vs FESOM2, 195 records (3-step)
  AND 325 records (5-step), BOTH whichEVP=0 and whichEVP=1 (`tools/run_lifecycle_fullynative_gate_core2.sh`). **M3f-4
  (multi-rank) ✅ DONE 2026-06-23 (first try, L43) → M3 COMPLETE: the FULLY-NATIVE lifecycle byte-matches FESOM2
  `max|Δ|=0` at MULTI-RANK (CORE2 dist_2 + dist_8, BOTH whichEVP, 195 + 325 records) —
  `tools/run_lifecycle_fullynative_gate_multirank.sh`. The M3a–e kernels were already MR-ready (M2.12 optional-partit);
  the only new work was `mesh%ocean_area` allreduce + the `mod_forcing_bulk` nNodL loop bound + the
  `fesom_lifecycle_native_mr` driver. Tag `m3`.**
  **CURRENT STATUS (2026-06-22):** the whole M2 dynamical core is BYTE-COMPLETE — `max|Δ|=0` vs FESOM2 on pi 1-rank
  (per-kernel + whole-step) AND across multi-step CORE2 (unforced + forced lifecycle, incl. the free-surface CG) AND
  **the WHOLE ASSEMBLED ocean step now byte-matches FESOM2 at MULTI-RANK** (pi dist_2 + dist_8). **M2.12 (multi-rank,
  MVP exit) is ✅ DONE: M2.12a (geometry) + M2.12b (tracer advection) + M2.12c-1 (pre-SSH DYNAMICS → ssh_rhs) + c-2
  (SSH stiffness + free-surface CG → d_eta) + c-3 (post-SSH ALE update + tracer SOLVE + the WHOLE step via step_oce,
  all 65 substep records `max|Δ|=0`).** The whole-model multi-rank byte-match — the architectural MVP — holds at dist_2
  AND dist_8 with no 1-rank regression. **MVP CLOSED 2026-06-22: baseline re-confirmed after a clean rebuild (all 8
  byte-gates `max|Δ|=0`/green) → M2.12b–c committed → tag `m2-mvp` created.** **PRODUCTION VALIDATION DONE 2026-06-22:
  the multi-rank FREE-RUNNING CORE2 lifecycle (dt=1800, reduced-M2) is now BYTE-EXACT vs FESOM2 (195 records
  `max|Δ|=0`); it exposed two latent multi-rank bugs — a UV halo over-exchange (FIXED, eDim) and a `tr_xy` halo "drift"
  that turned out to be an OpenMPI `vader` KNEM single-copy corruption of large (>131 KB) messages, fixed by
  `env.sh` exporting `OMPI_MCA_btl_vader_single_copy_mechanism=none` (LESSONS L35; not a FESOM3 bug).** Next milestone:
  M3 (sea ice EVP + `oce_fluxes` air-sea budget). See "Next task".

A compact list; full per-milestone detail + the gate recipes are archived in
[`HANDOFF-archive.md`](HANDOFF-archive.md).

### Completed milestones (all `max|Δ|=0` vs FESOM2 at the stated anchor)
- **M0** Foundation — tag `m0`. CMake + params/types/partitioning/halo/dump/mesh/step; 13/13 ctest (Intel+GNU dp).
- **M1** Tracer advection — tag `m1`, pi 1-rank. Geometry + horiz/vert advection (upw1/MUSCL/QR4C) + FCT (Zalesak) limiter + assembled step (`tools/run_advhor_gate.sh`).
- **M2.1–M2.9** Minimal dynamical core — pi 1-rank. EOS/`hpressure`/N² + PGF + Coriolis/AB2 `vel_rhs` + momentum advection + biharmonic viscosity + implicit vertical viscosity (TDMA) + SSH (CSR stiffness + preconditioned CG) + ALE update + PP mixing + `mo_convect` + tracer solve, AND the whole `step_oce` assembly (`tools/run_pressure_gate.sh` 57 fields + `tools/run_step_gate.sh` 65 records).
- **M2.10** Forcing — pi 1-rank, the FIRST netCDF I/O. Read + NCAR bulk coeffs/wind stress + SW penetration (`tools/run_forcing_gate.sh` 20 fields).
- **M2.11a** CORE2 geometry — CORE2 1-rank, 19 fields (`tools/run_geom_gate_core2.sh`); closed the CW-swap caveat (≈244654/244659 swaps).
- **M2.11b** Initial conditions (`do_ic3d`) — CORE2 1-rank, the FIRST 3D netCDF read (`tools/run_ic_gate_core2.sh`).
- **M2.11c** Multi-step CORE2 lifecycle — the whole dynamical core byte-exact over multiple steps, UNFORCED + FORCED (`tools/run_lifecycle_gate_core2.sh` + `…_forced_…`, 195 records each). The CORE2 free-surface-CG "reproducibility floor" was SOLVED — a vectorised preconditioner divide, one `!DIR$ NOVECTOR` (LESSONS **L29**).
- **M2.12a** Multi-rank local-mesh remap + geometry gate — pi dist_2 + dist_8, 19 fields every rank (`tools/run_geom_gate_multirank.sh`); the FIRST multi-rank byte-gate (LESSONS **L30**). Detail in "Next task".
- **M2.12b** Multi-rank tracer advection gate — pi dist_2 + dist_8, 13 fields every rank (`tools/run_advhor_gate_multirank.sh`): the `find_neighbors` halo dance + `exchange_elem_full(tr_xy/elem_area)` + `exchange_nod(fct_LO/fct_plus/minus)` + owned/halo loop bounds, all byte-exact; `del_ttf` (FCT + MUSCL) `max|Δ|=0` (LESSONS **L31**). Detail in "Next task".
- **M2.12c-1** Multi-rank pre-SSH DYNAMICS chain — pi dist_2 + dist_8, `max|Δ|=0` (25 records: density/pressure/bvfreq/Kv/ssh_rhs) (`tools/run_stepdyn_gate_multirank.sh`). Lifted the whole chain `compute_vel_nodes → pressure_bv(+smooth_nod) → pressure_force_4_linfs → oce_mixing_pp → mo_convect → compute_vel_rhs(+momentum_adv_scalar) → viscosity_filter(visc_filt_bidiff) → impl_vert_visc_ale → compute_ssh_rhs_ale` to multi-rank via the M2.12b OPTIONAL-`partit` pattern: `mod_part_bounds.local_dims` (owned/halo/full bounds) + the FESOM2 exchanges (`exchange_nod` UVnode/bvfreq-per-smooth-sweep/UVnode_rhs/ssh_rhs, `exchange_elem` U_c/V_c) + the visc interior-edge test on the GLOBAL id (`myList_edge2D(ed)>edge2D_in`). New `exchange_nod` rank-3 node-block variant in `mod_halo`. **NO 1-rank regression** (step 65 + advhor 40 + lifecycle 195 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L32**. Detail in "Next task".
- **M2.12c-2** Multi-rank SSH stiffness + free-surface CG — pi dist_2 + dist_8, `max|Δ|=0` on `d_eta` (30 records = the 25 c-1 + 5 `d_eta`; `tools/run_stepdyn_gate_multirank.sh`, substep 9 now compared). The FIRST multi-rank ITERATIVE solver (38 CG iters on 8 ranks, byte-exact). Lifted `init_stiff_mat_ale` + `solve_ssh_ale` (`ssh_solve_preconditioner`/`ssh_solve_cg`) via the same OPTIONAL-`partit` pattern: owned loops (`nNodO`), arrays `nNodL`, `exchange_nod(diag_values)` in the precond + `exchange_nod(pp/rr/x)` in the CG + new `allreduce_sum` (in `mod_halo`) for the dot-products; `rtol`/convergence denominators stay GLOBAL (`mesh%nod2D`). **KEY FINDING: NO mesh-infra extension was needed** — the HANDOFF/L32 hypothesis (extend `elem2D_nodes`/`gradient_sca` to the halo + `enforce_cw` on owned+eDim) was WRONG; the stiffness owned rows assemble FULLY LOCALLY (verified partition invariant: for owned edges both triangles are owned — FESOM2's own `elem2D_nodes`/`gradient_sca` are owned-only too). **NO 1-rank regression** (step 65 + pressure 57 + advhor MR + 13/13 ctest all `max|Δ|=0`). See LESSONS **L33**. Detail in "Next task".
  - **SSH solver now (`28af61c`, 2026-10-08): diagonal (Jacobi) preconditioner** — `pr_values` = 1/diag of the owned rows (length nNodO) instead of FESOM2's MITgcm symmetrised-Jacobi matrix, so M⁻¹r is local and the per-iteration `exchange_nod(rr)` (and the setup exchange of `diag_values`) is gone; one `exchange_nod(pp)` + two allreduces per iteration remain; CG loops fused. NOT bit-identical to FESOM2 any more (header of `src/oce/oce_ssh_solve.F90`). Env knobs (hybrid plan, `e4dd8d2`): `FESOM3_SSH_SOLTOL` (soltol, default 1e-5) and `FESOM3_SSH_MAXITER` (default 2000), read once at the first solve (`ssh_solve_env_overrides`); `tools/run_tri_regression.sh` solves to 1e-14 with them (L60).
- **M2.12c-3** Multi-rank post-SSH ALE update + tracer SOLVE + the WHOLE assembled step — pi dist_2 + dist_8, `max|Δ|=0` on ALL **65 substep records** (density/pressure/bvfreq/Kv/ssh_rhs/d_eta/hbar/eta_n/hnode_new/w/T/S/hnode; `tools/run_step_gate_multirank.sh`). **The whole multi-rank ocean step byte-matches FESOM2 — the architectural MVP.** Lifted `update_vel`/`compute_hbar_ale`/`update_eta_n`/`vert_vel_ale`(+`compute_CFLz`/`compute_Wvel_split`)/`update_thickness_ale` (`oce_ale.F90`) + `solve_tracers_ale`/`diff_tracers_ale`/`diff_part_hor_redi`/`diff_ver_part_impl_ale` (`oce_ale_tracer.F90`; the FCT advection was M2.12b) via the OPTIONAL-`partit` pattern (`owned_bounds`/`local_dims` + the FESOM2 exchanges), and threaded the optional `partit` through `mod_step_oce::step_oce` (absent ⇒ 1-rank verbatim; the 1-rank callers `fesom_stepdump`/`fesom_lifecycle` unchanged). New whole-step MR driver `fesom_stepfull_mr` (= `step_oce` through `partit`, gid-keyed per-rank dumps); the FESOM2 oracle is UNCHANGED (`fesom_step_dump` npes>1 already drives the whole `oce_timestep_ale` + dumps all substeps). pi ships only dist_1/2/8 (no dist_32; dist_2+dist_8 = the c-1/c-2 coverage). **NO 1-rank regression** (step 65 + pressure 57 + advhor MR dist_2/8 + stepdyn MR dist_2/8 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L34**. Detail in "Next task".

## Build

```bash
./configure.sh --compiler intel --precision dp --clean --build   # anchor
./configure.sh --compiler gnu   --precision dp --clean --build   # portability
cd build_intel_dp && ctest --output-on-failure                   # self-tests
```
Anchor = Intel + DP + FESOM2-v2.7.3-exact flags (see docs/LESSONS.md L1). Build dirs:
`build_<compiler>_<precision>/` (Release); `--debug` builds into `build_<compiler>_<precision>_debug/`
(kept SEPARATE so a Debug binary never clobbers the Release anchor and silently breaks the byte-gates —
the L14 footgun, now fixed in configure.sh). Login-node runs of 1–8 ranks are fine for self-tests.
Debug is for `-check all` OOB/FPE only, NOT byte-comparable to the Release oracle (L10).

## Oracle — PROVEN RUNNABLE (2026-06-19) ✅

The whole byte-gate pipeline is validated end-to-end:
- `bin/fesom.x` (FESOM2 v2.7.3, git SHA 9271ae92) is prebuilt. The node dump shim is
  already wired into `oce_ale.F90`.
- **`tools/run_oracle_pi.sh [run_dir] [nsteps] [nranks]`** runs pi in ~0.25 s and writes
  `<prefix>.<rank>` dumps (density/pressure/bvfreq/sw_alpha-beta/Kv/ssh_rhs/d_eta/hbar/
  eta_n/w/T/S/hnode, substeps 1–16). Quirks it handles: fresh `fesom.clock` = `0 1 1948`
  twice; `run_length_unit='s'` = steps; absolute ClimateDataPath.
- The dump byte-format is **identical** to FESOM3 `mod_dump`; `tools/dump_diff.py` parses
  real FESOM2 dumps (150 records) and self-compares `max|Δ|=0`. A fresh run byte-matches
  the committed fixture `test/refdata/pi_oracle_default/dump.*` → FESOM2 is deterministic.
- M1 (DONE) used the end-of-`ocean_setup` shim `fesom_advhor_dump.F90` (prescribe inputs, call
  the real kernels, dump, stop before forcing). **M2.1 reuses that shim pattern** for EOS/pressure
  (see "Next task"); later M2 kernels that run past forcing need the reduced-M2 namelist (below).

## Canonical references

- **Algorithm oracle + byte-gate target:** FESOM2 v2.7.3 `/home/a/a270088/port2/fesom2/src/`
  (tag `fesom2.7.3-cport-instr`). Transcribe FROM here, gate `max|Δ|=0` AGAINST here.
  Run via `tools/run_oracle_pi.sh` (proven).
- **Structural template (structure only, NOT math/flags):**
  `/home/a/a270088/fesom3/design_refs/tracer_dwarf/lib/`.
- **Reduced M2 oracle namelist:** `mix_scheme='PP'`, `Fer_GM=.false.`, `Redi=.false.`,
  `which_ale='linfs'`, `opt_visc=7` (NOT the shipped KPP/GM CORE2 namelist).

## Environment (Levante)

- Toolchains via modules (see `env/levante.dkrz.de/shell.{intel,gnu}`): intel-oneapi
  2022.0.1 + openmpi 4.1.2-intel; gcc 11.2.0 + openmpi 4.1.2-gcc. netCDF loaded (only
  needed M2.10+). Login `gfortran` is 8.5.0 with no MPI — always build via `configure.sh`.


## Gate recipes (archived)

The detailed per-gate recipes — the 1-rank FESOM2 oracle setup (hand-crafted `dist_1`, the
end-of-`ocean_setup` prescribe-and-stop shim, the FADVHDMP dump format), the M1 advection
kernel-gate recipe, the M2.1–M2.9a dynamics-kernel recipe (+ every per-kernel extension), the
M2.10 forcing gate, and the M2.11b IC gate — are in [`HANDOFF-archive.md`](HANDOFF-archive.md).
The live "Oracle" section above + the "M2.12 entry notes" below cover what the current work needs.
The hybrid-mesh gate migration (which FESOM2-oracle gates are retired, the roundoff gate
`tools/run_tri_regression.sh` that replaces them) is in "Hybrid (mixed triangle + quad) meshes" →
"Gate migration for hybrid meshes" below.

## Shared elements: one owner for `uv` (plan `plans/completed/20261012-shared-element-sync.md`, 2026-10-12)

**Ownership rule (FESOM2's, unchanged).** A rank owns every element that has at least one vertex it
owns, so an element on a partition border is owned by two or three ranks (CORE2: 562 at dist_2, 24,798
of 244,659 at dist_512). Each owner updates its own copy of the element fields.

**The CORE2 leak (9–12 Oct 2026).**
- Leaking: results7, results11, results17 from July 1965, results21. Global salt +1.4e-5 to +2.3e-5
  psu/yr, global-mean SSH drifting. Clean: results_orig, results, results2, results3, results5,
  results17 up to June 1965 (salt −1.6e-6 to −2.1e-6 psu/yr, mean SSH flat).
- Volume rules differ by code: f2 runs (results_orig, results_2.0) `nod_area(k,n) × hnode`, f3 runs
  `area(n) × hnode` (`src/mesh/mod_mesh_areas.F90:515`). Report: artifact `Gdj6DMQ4n2BQVZ4Si6cAUD` v2,
  figures and scripts in `/work/ab0995/a270029/fesom3_hybrid/runs_report/core2/`.
- Mechanism: the viscosity filter adds edge contributions in each rank's own edge order, so the `uv`
  copies of a shared element differed by roundoff; the `rvo_upwind` blend (results7 on, results5 before
  it) can turn that into O(1); the two ranks then computed different fluxes across their common edges and
  `ssh_rhs_old` no longer summed to zero. Evidence: in the leaking checkpoints the stored u, v disagree
  with the stored `ssh_rhs_old` at 86–107 nodes, all on dist_512 borders next to shared elements (clean
  checkpoints agree to 4e-11 everywhere); a restart re-synced the copies and paused the leak for about
  six months. LESSONS L70.

**The fix.**
- `partit%com_elem2D_shared` (`src/types/mod_partit.F90`), built by `build_com_elem2D_shared`
  (`src/infra/mod_shared_elem.F90`) at the end of `read_mesh_local` from the owned elements and the
  vertex owners in `com_nod2D` (no file format change, no MPI for the lists). The **designated owner**
  is the highest rank among the vertex owners: `decomp_redistribute_2d_r`
  (`src/io/mod_io_decomp.F90:220-237`) writes in rank order, so this is the copy restart and output keep.
  Blocks are sorted by global id with partners ascending. One-time checks: an `MPI_Alltoall` of the
  pairwise counts, a global-id self-check through the exchange, and the rule itself (an owned element
  without an owned vertex stops). Every stop goes through `die` (`MPI_Abort`).
- Startup log on rank 0: `shared-element sync: ON, <n> shared elements, <c> extra copies, max per rank:
  recv <r>, send <s>`. Production dist_512 (checked offline): `24798 shared elements, 25512 extra
  copies, max per rank: recv 139, send 136`.
- `update_vel` (`src/oce/oce_ale.F90`) calls `exchange_elem_shared(UV)` right after the owned update and
  before `exchange_elem(UV)`, every step, under the subset timer "(of which) shared sync" of
  "dynamics ssh". `FESOM3_SHARED_SYNC`: unset or `1` = on, `0` = off (the old unsynced copies, for A/B
  comparisons; the list is still built), any other value stops.
- Only `uv` is synced. `uv_rhsAB` and the EVP/mEVP stresses are computed from equal inputs (nodal fields
  their owners compute and exchange, per-element geometry from the same global coordinates) and stay
  bitwise equal; guards `test_vinv` V13 and `test_ice_dyn_mixed` D6. If a guard ever fails, the plan's
  "Fields" fallback says where to sync that field.

**Diagnostics.**
- `FESOM3_STEP_DIAG=1` (production driver): the step line gets, after `Smax`, `divsum=` (the ratio
  `|Σ ssh_rhs_old| / Σ|ssh_rhs_old|` over owned nodes; a linfs leak shows above ~1e-12; zstar adds the
  freshwater term, so there it is only a report), `dUV=`, `dAB=`, `dsig=` (max|copy − designated copy|
  of `uv`, `uv_rhsAB`, sigma; all 0 with the sync). `tools/run_cutout_short.sh` keeps its field numbers.
- `fesom_conserve`: one end-of-run line `shared elements (sync ON|OFF): max divsum= … dUV= … dAB= …
  dsig= …`; `FESOM3_DIVSUM_TOL` (unset = report only) stops the run if the max ratio exceeds it or, with
  the sync on, any spread is nonzero.

**Tests and gates.**
- `test_shared_elem`: S1 the list against an allgather oracle, S2 the exchange and the probes, S3 the
  mechanism through `update_vel` + `compute_hbar_ale` (ratio 9e-20 to 2e-17 with the sync, 7e-6 to
  4e-4 without, on pi np 1/2/8, pi_w np 2/4, the NEMO cut-outs np 2, core2_w np 8), plus the error cases
  (corrupt owner map, `FESOM3_SHARED_SYNC=2`) and the OFF log. `fesom_conserve_linfs_vinv_shsync_np2`
  (linfs, vinv, `rvo_upwind` 1, sea ice, `FESOM3_DIVSUM_TOL=1e-12`). ctest 167/167.
- `tools/run_conserve_pi.sh` SHARED block (linfs + ice, np 2 and 8): max ratio 1.6e-16 to 3.1e-16,
  spreads 0. Sync off: `uv` spread 6.9e-18 after 20 steps, already amplified to 5.8e-11 in `uv_rhsAB`.
- Bitwise against the frozen pre-change build (`run_tri_regression.sh` with all `TOL_*=0`): np 1 pi and
  CORE2, and np 2 pi with `FESOM3_SHARED_SYNC=0`, every printed d = 0. The default
  `run_tri_regression.sh` cases pass 7/7; the np>1 differences stay at their earlier levels.
- Restart: `tools/run_restart_gate_core2.sh` is strict at every np now (floor 1e-12 only with
  `FESOM3_SHARED_SYNC=0`): CORE2 np 1 and np 2, mid-run and month-boundary splits, checkpoint (full
  state incl. ice, sigma, velocity, AB) and live dumps all `max|Δ| = 0`, with the default physics and
  with `FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0`. Before the sync the np 2 boundary case needed the
  floor (temp 2.2e-16). `run_restart_gate_multirank.sh` (np 2 → np 8): `max|Δ| = 0`.
- Cost (CORE2 np 8, 48 steps, production physics of `tools/prod_env.sh`, no output, login node): "(of
  which) shared sync" 0.51 ms/step (max 0.66) of a 974 ms loop, 0.05 %; with `FESOM3_SHARED_SYNC=0` the
  loop took 960 ms (login-node noise is larger than the sync). The same runs: sync on, `divsum` <= 1.4e-15
  and `dUV = dAB = dsig = 0` in all 48 steps; sync off, `dUV` grows to 3.9e-16 and `dAB` to 5.8e-11 (the
  ratio is still at roundoff after one day).

**Production check (the user's runs).** Rebuild the production binary when you decide. Resume
results17 from `fesom.1967.001.00000` (or results21 from `fesom.1959.001.00000`) with the same settings,
at least into autumn 1967 (the leak paused for about six months after each restart). With
`FESOM3_STEP_DIAG=1` expect `divsum` around 1e-15 to 1e-13 and `dUV = dAB = dsig = 0`; the global-mean
SSH should stay at its restart value (0.0888 m for results17) and the salt trend return to about −2e-6
psu/yr. Optional A/B from the same restart: `FESOM3_SHARED_SYNC=0` (`dUV` > 0, the leak returns) and
`FESOM3_RVO_UPWIND=0` (the trigger).

## Hybrid (mixed triangle + quad) meshes (plan `plans/completed/2026-10-08-hybrid-meshes.md`, 2026-10-08/09)

The model, its setup tool `fesom_ini` and the Python mesh tools handle triangle, quad and mixed meshes
in ONE code path: every element has 4 vertex slots (a triangle repeats its 1st vertex), clockwise, and
every kernel uses unified 4-slot formulas — there is no triangle/quad branch in any kernel. Commits
`d4c6f8d` (Task 1) … `d07f44f` (Task 19) on `bottom-at-vertices`, `hybrid:` commits, one or more per
task (the messages list what changed), unpushed. Triangle results are NOT bit-identical to FESOM2 or to the
pre-hybrid code any more (roundoff, gated by `tools/run_tri_regression.sh`, see "Gate migration for
hybrid meshes" below). Scratch work, logs and the frozen reference binaries of the plan:
`/work/ab0995/a270029/fesom3_hybrid/` (task log `progress-hybrid-meshes.txt`); test meshes in
`/work/ab0995/a270029/meshes/` (list under the setup workflow). Order of the subsections: file formats →
arrays and conventions → geometry → kernels → exchanges → SSH matrix and checkerboard → sea ice → tests
and tools → setup workflow → CORE2 throughput → NEMO job → gate migration → open findings. Lessons
from the plan: `LESSONS.md` L60–L69.

### File formats

- **`elem2d.out`**: count, then 3 OR 4 node ids per line; the column count is taken from the first
  line and every line must match (`read_elem_line`, `src/mesh/mod_mesh_read.F90:122`; a malformed line →
  message + error stop, ctest `test_mesh_read_malformed_np1`); 3 columns → 4th := 1st. Triangle ⇔
  `elem2D_nodes(4,e) == elem2D_nodes(1,e)`. `enforce_cw_orientation` (`:228`) turns an anticlockwise
  quad by swapping 2↔4 (a triangle 2↔3, as before); `fesom_ini` checks the same in `test_tri_ini`
  (sub-triangles (1,2,3), (1,3,4); flips reported and fixed in memory, a bow-tie/degenerate quad stops).
- **Depth**: `fesom_ini`'s `namelist.config` keeps FESOM2's `use_depthfile='aux3d'` (`aux3d.out`: nl,
  zbar, node depths) and gains `use_depthfile='nemo'`: `depth.out` = nl, the positive levels, the
  positive node depths (= levels(count+1)), then the per-node layer counts → `nlvls.out` = count + 1 as
  given, `elvls.out` = min over the vertices, no thresholding or isolated-cell iteration. The model
  takes node levels ONLY from `nlvls.out` and zbar from `aux3d.out`, or from `depth.out` when
  `aux3d.out` is absent (`open_depth_file`, `:169`; depth.out node depths stored negated, unused).
  Element levels: `nlevels(e)` = min, `ulevels(e)` = max over the vertices.
- **`nod2d.out`** of the NEMO meshes: 4th column = land-boundary flag (as `nemo_meshq.m`).
- **Edge files** `edges.out`, `edge_tri.out`, `edgenum.out` (edge2D, edge2D_in): FESOM2 layout,
  interior edges first (viscosity and the ice boundary conditions rely on that order); `edge_tri.out`
  keeps FESOM2's on-disk encoding (missing 2nd element ≤ 0) — the model converts it at read time
  (`mod_mesh_read.F90:93` 1-rank, `:535` multi-rank).
- **`dist_N/`**: `rpart.out`, `my_list%05d.out`, `com_info%05d.out` keep their exact FESOM2 layout
  (FESOM2 and `tools/make_dist1.py` still read them). `dist_N/format.out` (one integer, version 1 =
  `DIST_FORMAT_VERSION` of module `dist_ext`, `mesh_part/src/oce_local.F90:14-283`) announces the
  APPENDED sections: `my_list` + `eXDim_nod2D` and the global ids of the double-halo nodes (vertices of
  the local elements 1..myDim+eDim+eXDim_elem2D that are not among the nodes 1..myDim+eDim_nod2D;
  grouped by owner, ascending within; local numbers myDim+eDim+1, …); `com_info` + `com_nod2D_full`
  then `com_edge2D`, each in FESOM2's 8-record com-struct layout, per peer both sides ordered by global
  id. The model reads them when `format.out` exists (`partit%dist_format`, 0 = FESOM2 layout; an
  unknown version stops) in `read_dist_partition` (`src/infra/mod_partitioning.F90:125`); the 1-rank
  synthesis `synthesize_1rank` (`:86`) and `make_dist1.py` provide empty lists (dist_format = 1). An old
  `dist_N` without `format.out` still runs; only the two new exchanges then stop with a message.
- **Output** (Task 14): `face_nodes`/`face_edges` with 4 columns and the UGRID `_FillValue` −1 for a
  triangle's 4th, element coordinates = wrapped centroids — see "Gate migration for hybrid meshes" →
  "Output format changes".

### Mesh arrays and conventions (`src/types/mod_mesh.F90`, `src/types/mod_partit.F90`)

- `MAX_NV = 4`, `MAX_ADJACENT = 32` (`mod_mesh.F90:30-31`). `elem2D_nodes(MAX_NV,·)` clockwise,
  triangle 4th = 1st; `elem2D_nnodes(·)` = 3 or 4 — used only where a repeated vertex would matter
  (unique-vertex loops, plain counts: IC extrapolation, SSH assembly index set, ice mass matrix and FCT
  limiter).
- `elem_edges(MAX_NV,·)` (sides 1-2, 2-3, 3-4, 4-1) and `elem_neighbors(MAX_NV,·)` (neighbour across
  side k), owned elements, built in `build_elem_adjacency` (`src/mesh/mod_mesh_areas.F90:358`): a
  triangle's slot 4 repeats slot 1; a MISSING neighbour is the element ITSELF (no zero index anywhere).
- `edge_tri(2,·)`: a missing 2nd element = the 1st, and `edge_cross_dxdy(3:4,e) = 0` there, so a term
  weighted by the cross vector vanishes without a test. Boundary tests are `edge_tri(2,e) /=
  edge_tri(1,e)` and `nb == el`; every former `el(2) > 0` / `el(i) < 1` / `nb <= 0` test (src and test)
  was converted in the same change as the convention (Task 7), otherwise boundary edges count twice.
- `gradient_sca(2*MAX_NV,·)`: x-coefficients in rows 1:4, y in 5:8 (FESOM2: 6 rows, x 1:3 / y 4:6); a
  triangle's slot 4 = 0; OWNED elements only. `gradient_vec` is still not computed; when it comes it
  must respect "missing neighbour = self".
- `elem_w(MAX_NV,·)` (owned elements): kite (median-dual) weight of each vertex, Σ = 1, a triangle
  exactly (1/3, 1/3, 1/3, 0). `nod_in_elem2D_w(MAX_ADJACENT,·)`: kite weight of node n in its k-th
  element `nod_in_elem2D(k,n)`, computed for owned nodes and halo-exchanged (`compute_nod_weights`,
  `:301`). Both are in the binary mesh layout.
- **The two rules every kernel follows:** element mean = `sum(elem_w(:,e)*f(elem2D_nodes(:,e)))`;
  node gather = `Σ_k elem_area(e_k)·nod_in_elem2D_w(k,n)·x(e_k) / area(n)` (FESOM2's `Σ A·x /
  (3·area)` hid the 1/3 in the 3; every such gather was converted); gathers over WET elements
  only normalise by the wet kite shares `Σ A·w` instead (`compute_sigma_xy`, `smooth_nod`, ice
  `ocean2ice` — FESOM2's partial-patch normalisation).
- Node arrays keep size nNodL = myDim+eDim. The mesh FILE arrays (`coord_nod2D`, `geo_coord_nod2D`,
  `bc_index_nod2D`, `depth`, `nlevels_nod2D`) and `myList_nod2D` hold myDim+eDim+eXDim_nod2D entries;
  `owned_bounds(…, nNodX=)` (optional, `src/infra/mod_part_bounds.F90:53`) returns that size (npes 1:
  nod2D).

### Geometry (`src/mesh/mod_mesh_areas.F90`, Task 7)

- Centroid and area from the sub-triangles (1,3,4) + (1,2,3) (qq), area-weighted and cyclic-wrapped
  (`elem_center`, `:181`, public — also the output element coordinate); `elem_area` = lon/lat area ×
  `elem_cos` (cosine of the centroid latitude). For a triangle both sub-triangles coincide with the
  plain formulas.
- Kite weights (`compute_elem_weights`, `:267`): ½·[A(j, j−1, c) + A(j, j+1, c)] normalised by the
  sum of the 4 side triangles (= the element area for a convex quad, so Σ `elem_w` = 1 to roundoff);
  node area = Σ `elem_area·elem_w` (`compute_node_areas`, `:472`).
- `gradient_sca` Green–Gauss (`compute_gradient_sca`, `:326`): coefficient of vertex j ∝
  (x_{j+1} − x_{j−1})/(2A) over the nv vertices, longitudes relative to vertex 1 (one planar embedding,
  as FESOM2). Exact for linear fields on every convex element; for triangles = P1 to roundoff
  (`test_geom_mixed` G4, G7).
- `compute_mesh_resolution` (`:556`) with kite-weighted means. `mod_mesh_rotate.F90` needed no change.
- Analytic test patches (`src/mesh/mod_mesh_analytic.F90`): `generate_mixed_patch` (quads + triangle
  pairs, jittered interior nodes), `generate_quad_patch` (parallelograms, optional skew); `build_edges`
  orients edges as FESOM (edge_tri(1) left of edges(1)→edges(2)) and puts interior edges first.

### Kernels: unified 4-slot formulas (Tasks 10–16)

Every kernel loops over the 4 slots with `gradient_sca(1:4|5:8)`, `elem_w` means and `nod_in_elem2D_w`
gathers; a triangle's slot 4 contributes exactly 0 (gradient 0, weight 0) or repeats a value where only
max/min is taken (FCT bounds, `mo_convect`). Converted: ocean dynamics (`oce_ale.F90`
`compute_vel_nodes`, `update_vel`, `dhe`, zstar `helem`; `oce_dyn_velrhs.F90`; `oce_pgf.F90` incl. the
Shchepetkin per-vertex loops over MAX_NV; `oce_dyn_ivertvisc.F90`; `oce_dyn_vinv.F90`: KE and vorticity
with `elem_w`, `rvo_upwind` side normals oriented toward the centroid, a triangle's degenerate side 4
has normal 0); tracers (`oce_tracer_grad.F90`, `oce_ale_tracer.F90` `tr_xynodes`, `oce_adv_tra_fct.F90`
bounds, `oce_pressure_bv.F90` `compute_sigma_xy`/`smooth_nod`, `oce_initial_state.F90` extrapolation
over the unique vertices); the MUSCL up/downwind search `find_up_downwind_elements`
(`src/oce/oce_muscl_adv.F90:57`; helper `elem_sides`: b = next, c = previous vertex for both types; the
exact collinearity compares became |·| < 1e-10 rad; the array name `edge_up_dn_tri` was kept); mixing
and GM (`oce_fer_gm.F90` `fer_uv`, KPP `viscAE`, TKE/PP `Av`, `oce_mo_conv.F90`); forcing and coupling
element stress (`mod_forcing_bulk.F90`, `mod_ice_oce_coupling.F90` — its `elnodes(3)` shape bug went
with it). The edge-based kernels (`oce_dyn_visc.F90`, `oce_adv_tra_hor/flux/driver.F90`) needed only
the boundary-test conversion. Exceptions left in `src/` (Task 19 greps, all documented in the plan):
`mod_geom_dump` (FESOM2 layout, stops on quads), the analytic generator, the sub-triangle centroid
`/(3(s1+s2))`, the exact triangle kite 1/3, the spline and the KPP cube roots; `fesom_ini` keeps literal
`/3.0`, `/4.0` divisors (L7).

### Exchanges and halo rings (Task 8, `src/infra/mod_halo.F90:159-203`)

- **Rings.** Nodes: 1..myDim (owned) | eDim_nod2D (FESOM2 halo) | eXDim_nod2D (NEW double halo: the
  remaining vertices of all local elements, i.e. of the eDim edge-sharing and the eXDim node-only
  element rings, L56). Edges: FESOM2's 1..myDim_edge2D holds EVERY edge with an owned node, so a shared
  edge lies in the owned range of two ranks; then eDim_edge2D.
- **`exchange_edge`** (2D/3D real, over `com_edge2D`): the UNIQUE owner = owner of `edges(1,edge)` sends
  to every other rank holding the edge, also where it lies in that rank's owned range.
  **`exchange_nod_full`** (2D/3D real, over `com_nod2D_full`): every halo node eDim + eXDim from its
  owner; arrays sized nNodX. Both run on the existing broadcast-only cores (`core_2D_r`/`core_3D_r`) like
  `exchange_nod`; `require_dist_format` aborts with "<list> is not available: dist_N has no format.out"
  on a FESOM2-layout `dist_N`; at 1 rank the lists are empty (no-op). No kernel uses them yet (added for
  future developments). Test `test_halo_ext` (round trips at np 2/4/8 on pi_w, core2_w and the cut-outs,
  edges in two owned ranges > 0 required, no-format error case).

### SSH stiffness matrix and the checkerboard null mode (Task 9)

- `init_stiff_mat_ale`/`update_stiff_mat_ale` (`src/oce/oce_ssh_rhs.F90:87`, `:252`): CSR row = self
  FIRST (the Jacobi preconditioner reads the diagonal there), the edge neighbours in edge order (as
  before), then the remaining vertices of the elements around the node (quad diagonals); the row length
  is counted first (no fixed `n_pos(12,·)`); scatter over the 4 slots. Triangle meshes keep the
  IDENTICAL CSR and bitwise-equal values (`test_ssh_stiff` S5). The ice mass matrix reuses this pattern.
  The solver itself (diagonal preconditioner, `FESOM3_SSH_SOLTOL`/`FESOM3_SSH_MAXITER`) is described
  under the M2.12c-2 milestone above.
- **Checkerboard property.** Green–Gauss gives c₁ = −c₃ and c₂ = −c₄ on ANY quad, so the pattern
  (+1,−1,+1,−1) has zero gradient in every quad: on a 2-colourable all-quad mesh the checkerboard is an
  EXACT null mode of the stiffness part of the SSH operator (|K f|/|K| ≈ 1e-16 on the parallelogram
  patch AND on the NEMO quad cut-out — not only on parallelograms, as the plan first assumed) and of the
  pressure gradient; only the mass (area/dt) part controls it. Asserted in `test_ssh_stiff` S8 so that
  it cannot change silently. Measured over 10 days on the cut-outs (`tools/cutout_diag.py checker`): the
  SSH null mode does not grow, the tracer grid-scale energy does (open findings below).

### Sea ice (Tasks 15–16)

- **EVP/mEVP** (`src/ice/mod_ice_dyn.F90`): strain rates over `gradient_sca(1:4|5:8)` with
  kite-weighted metric terms (`metric_factor·sum(elem_w·v)`); `stress2rhs` node shares
  `elem_area·elem_w(k)` (the exact negative adjoint of the strain rates, `test_ice_dyn_mixed` D2); SSH
  tilt `g·elem_area·elem_w(k)`; `ocean2ice` gather `Σ A·w·UV / Σ A·w` over the wet elements.
- **FCT with the user's mass matrix** (`src/ice/mod_ice_fct.F90:73-98` tables,
  `src/ice/mod_ice_setup.F90:189` `ice_mass_matrix_fill`): consistent mass matrix from the element
  contributions `elem_area·ice_mcoef(·,·,nv)` — triangle S/6 diagonal, S/12 off-diagonal; quad
  (regular-square simplification) row n1 (S/36)·(4, 2, 1, 2) for (n1, n2, n3, n4), rotated per row. **Lumped mass**
  `ice%work%fct_lumpedmass` = row sums (S/3 per triangle vertex, S/4 per quad vertex), exchanged to the
  halo; it REPLACES `mesh%area` in ALL FCT steps (low order, the high-order iteration, the antidiffusive
  fluxes); a row without wet elements keeps the node area. Antidiffusive operator `ice_acoef` =
  consistent − lumped: triangle (−2, 1, 1)·S/12, quad (−5, 2, 1, 2)·S/36; `fct_fluxes(elem, MAX_NV)`;
  flux, sum, limiter and update loops over `elem2D_nnodes` (a 4-slot limiter would visit vertex 1 twice).
  TG rhs: `dx(n)·(M_e u)_q`, diffusion over the 4 gradient slots, TG term with ū = Σ `elem_w`·u. A CSR
  column missing from the SSH pattern → message with global ids + MPI_Abort (ctest
  `test_ice_mass_missing_col_np1`). Ice mass is conserved as Σ `fct_lumpedmass`·m; on quads the lumped
  mass S/4 differs from the kite node area the dynamics use (accepted by the user). For triangles
  lumped = node area and every term equals FESOM2 to roundoff (`test_ice_mass` F5).

### Tests and tools of the plan

- **ctest** (153 in all, ~50 s with `-j 8`): `test_mesh_read` (+ flipped elements, 5 malformed
  `elem2d.out` cases), `test_geom_mixed` (+ quad kite weights / centroid vs independent shoelace
  formulas, `mesh_resolution`), `test_halo_ext`, `test_ssh_stiff` (+ quadratic-field values on the
  parallelogram patch), `test_pgf_mixed` (+ `compute_sigma_xy`, `smooth_nod`, `fer_gamma2vel`,
  `extrap_nod3D`), `test_updown_mixed`, `test_coupling_mixed`, `test_ice_dyn_mixed`, `test_ice_mass`
  (+ FCT bounds), `test_partit_*` (`FESOM3_EXPECT_FORMAT=1`), `fesom_ini_*` (`tools/check_fesom_ini.sh`:
  FESOM2 level finder on quads, the 'nemo' error stops, non-convex quad),
  `meshdiag_{mixed,quad,pi_w}`, `fesom_conserve_mixcut_*` (PP/KPP/TKE/GM/REDI) and
  `fesom_conserve_{q,mix}cut_ice_*` (driver knobs `FESOM3_CONSERVE_ICE`, `FESOM3_WHICHEVP`), plus mixed
  and quad cases of `test_vinv`, `test_ivertvisc`, `test_muscl_onthefly`. The hybrid tests share
  `test/test_util.F90` (mesh setup, check/finish, reductions). Kernel tests pick the mesh with
  `FESOM3_MESH_DIR`, the analytic patches with `FESOM3_TEST_PATCH=mixed|quad`, and `FESOM3_TEST_ROTATE`
  (0: NEMO cut-outs unrotated; 1: pi_w rotated — the unrotated pi has elements around the geographic
  pole). The launcher flags `FESOM3_MPI_LAUNCH_FLAGS` (`CMakeLists.txt:157`) include the knem workaround
  (L61). Every mesh-dependent hybrid ctest reads its mesh from the HARD-CODED root
  `/work/ab0995/a270029/meshes` (`test/CMakeLists.txt`, not a cache variable); a case whose mesh is
  absent at configure time (`if(EXISTS <meshes>/<mesh>/nod2d.out)`) is registered DISABLED ("Not Run
  (Disabled)"), and a missing `dist_N`/`com_info` file in an existing mesh directory FAILS — so
  "153/153" holds only on a filesystem with these meshes (this account's).
- **Scripts:** `tools/mesh/test_nemo_meshq.sh`, `tools/mesh/test_mesh_tools.sh`,
  `tools/run_ini_compare.sh` (fesom_ini == FESOM2's own mesh_part build on pi/CORE2),
  `tools/run_ini_mixed.sh` (NEMO meshes + `tools/mesh/check_dist.py` + error cases),
  `tools/run_tri_regression.sh` (triangle roundoff gate; tiers incl. ICE_EARLY = the sea-ice output
  of the first steps), `tools/run_cutout_short.sh` (5-day prod/prodpp/cons runs on both
  cut-outs, np 1 vs 2/8 from one checkpoint, prodpp gated at the tight TOL_PP = 1e-5, checkerboard
  time series via `tools/cutout_diag.py`; its ICE class is vacuous -- the tropical cut-outs stay
  ice-free and ice seeded into the checkpoint melts within 2 steps, so a 1-vs-N sea-ice compare on
  quads needs a run that keeps ice, which no driver switch provides yet),
  `tools/check_meshdiag.py`, `tools/run_mpi_knem_stale_repro.sh` (L61). Shared: `tools/prod_env.sh`
  (the production physics env of the two run scripts), `tools/lib_ini.sh` (the run_ini harness).

### Mixed-mesh setup workflow (triangle, quad and mixed meshes)

Every command below was run as written on 2026-10-09 (scratch copy, script
`/work/ab0995/a270029/fesom3_hybrid/t18_workflow/run_workflow.sh`, log `t18_workflow.log` there);
times and peak memory are login-node values for the full ORCA025 mesh. The Python tools need numpy:
`source /usr/share/Modules/init/bash; module load python3; hash -r` (the system `python3` has none).

1. **NEMO → FESOM quad mesh** — `python3 tools/mesh/nemo_meshq.py <nemo_dir> <out_dir>`, e.g.
   `python3 tools/mesh/nemo_meshq.py /work/ab0995/a270029/meshes/nemo025q <scratch>/nemo025q_py`.
   Reads `mesh_hgr.nc`, `mesh_zgr.nc`, `mask.nc`; writes `nod2d.out` (4th column = land-boundary
   flag), `elem2d.out` (4 columns, SW-NW-NE-SE), `depth.out` (47 positive levels, node depths, layer
   counts). Faithful port of `nemo_meshq.m` (identical to the MATLAB output except one fold point,
   1 float32 ulp; `tools/mesh/test_nemo_meshq.sh`). 5 s, 0.7 GB. The output still holds 40 dead-end
   quads, which fesom_ini rejects.
2. **Dead ends** — `python3 tools/mesh/remove_deadends.py <in_dir> <out_dir>`: iterative, 46
   elements / 90 nodes removed → 915 866 nodes, 902 408 quads (`nemo025q_fesom`). 12 s, 0.9 GB.
   `depth.out` layer counts are kept; a node left above the bottom invariant is reported (none here).
   The mesh tools write into `<out_dir>` and remove fesom_ini outputs already there (edge/level files,
   `dist_*`), which belong to the old elements.
3. **Mixed mesh (optional)** — `python3 tools/mesh/split_quads.py <in_dir> <out_dir> --band -5 5
   --every 5`: every quad with centroid latitude in the band plus every 5th other quad split into two
   triangles (shorter diagonal; the longer one where the shorter leaves a dead-end triangle) →
   `nemo025mix`, 685 452 quads + 433 912 triangles, same nodes and `depth.out`. 5 s, 0.8 GB.
4. **Regional cut-out for tests (optional)** — `python3 tools/mesh/cut_region.py <in_dir> <out_dir>
   --lon -20 -5 --lat -2 12` → `nemo025q_cut` / `nemo025mix_cut` (2 509 nodes; largest connected
   piece, dead ends removed, boundary flag recomputed, layer counts at the cut lowered to keep the
   model's bottom invariant). 3 s.
5. **Check every mesh** — `python3 tools/mesh/check_mesh.py <mesh_dir>`: prints a `SUMMARY` line
   (nodes, elems, quads, tris, edges, bedges, ccw, nonconvex, degenerate, unused, multiedge,
   deadends, layers, area; `layers` = `depth.out` nodes above the bottom invariant) and PASS/FAIL;
   exit 1 on any error incl. dead ends (the raw `nemo_meshq.py` output fails with `deadends=40`).
   4 s, 0.9 GB.
6. **Partition: `fesom_ini`** (FESOM2's `mesh_part`, ported: sources `mesh_part/src/`, METIS 5.1.0 in
   `lib/metis-5.1.0/`). Built together with the model by `./configure.sh` / `cmake --build` as
   `<build>/bin/fesom_ini` (`build_intel_dp/bin/fesom_ini` once the user rebuilds from this branch;
   scratch: `/work/ab0995/a270029/fesom3_hybrid/build/bin/fesom_ini`). Serial; run it in a work
   directory holding `namelist.config`:
   ```
   &paths
   MeshPath="<mesh_dir>/"
   /
   &geometry
   use_depthfile="nemo"
   /
   &run_config
   /
   &machine
   n_levels=1
   n_part=<N>
   /
   ```
   (`use_depthfile="nemo"` = NEMO `depth.out`; FESOM triangle meshes leave it out = `'aux3d'`; the
   rotation defaults 50/15/−90 equal the model's hard-coded ones), then
   `cd <work_dir> && <build>/bin/fesom_ini > fesom_ini.out 2>&1` (inside a batch job
   `srun -n 1 -N 1 …`). It STOPs with exit status 0 also on errors — judge by `Total time` in the log
   and `dist_<N>/rpart.out`. It stops on a self-intersecting or non-convex quad and on a `depth.out`
   layer count outside 2..nl−1 (the model needs ≥ 2 layers per column). It writes INTO `MeshPath`: `edges.out`, `edge_tri.out`, `edgenum.out`
   (edge2D, edge2D_in), `nlvls.out`, `elvls.out`, `elvls_raw.out` and `dist_<N>/` = `rpart.out`,
   `my_list%05d.out`, `com_info%05d.out` in the FESOM2 layout plus the appended hybrid sections
   (double-halo nodes, `com_nod2D_full`, `com_edge2D`) announced by `dist_<N>/format.out` (version 1).
   It KEEPS existing edge/level files: delete `edges.out edge_tri.out edgenum.out elvls.out
   elvls_raw.out nlvls.out dist_*` first when the mesh changed (`clean()` in `tools/run_ini_mixed.sh`).
   One run per N. ORCA025: 128 parts 31 s; 1024 parts 83 s (METIS 42 s, writing 34 s), peak 0.3 GB;
   at most 16 neighbour ranks per communication list for 128 and 1024 parts (the model's
   `com_struct` holds `MAX_NEIGHBOR_PARTITIONS = 32`; `fesom_ini` stops above it for every list before
   filling it, `read_com` before reading it — count with
   `/work/ab0995/a270029/fesom3_hybrid/t18_workflow/max_neighbours.py <mesh_dir> <N>` for very high
   N). `dist_1`: `python3 tools/make_dist1.py <mesh_dir>` (only for `check_dist.py`/FESOM2 tools; the
   model builds its 1-rank partition itself).
7. **Check the partition** — `python3 tools/mesh/check_dist.py <mesh_dir> <N> [<N> …]` (edges, edge
   order, levels, ownership, halos, all communication lists incl. the appended sections) → `PASS`;
   1024 parts of ORCA025: 73 s, 0.9 GB.
8. **Model** — `FESOM3_MESH_DIR=<mesh_dir>/` and exactly N ranks (the driver reads
   `dist_<number of ranks>/`). `elem2d.out` with 3 or 4 columns is detected; zbar comes from
   `aux3d.out`, or from `depth.out` when `aux3d.out` is absent; node levels only from `nlvls.out`.
   No other mesh setting exists (rotation hard-coded); IC, forcing, runoff, SSS and chlorophyll are
   regular-grid files interpolated at run time. Only `FESOM3_STEP_PER_DAY` is resolution-dependent.
   Smoke test: the mixed cut-out made by the commands above, production physics, 4 ranks, 8 steps.
9. **Mesh diagnostics (optional)** — `FESOM3_MESHDIAG_OUT=<store> fesom_meshdiagdump`, checked by
   `tools/check_meshdiag.py <mesh_dir> <work_dir> <fesom_meshdiagdump> <np[,np…]> [launcher …]`
   (ctests `meshdiag_{mixed,quad,pi_w}`).

Regression scripts of the workflow: `tools/mesh/test_nemo_meshq.sh`, `tools/mesh/test_mesh_tools.sh
[nemo_dir] [out_root]` (steps 2–5; `out_root=/work/ab0995/a270029/meshes` made the persistent test
meshes), `tools/run_ini_mixed.sh` (step 6–7 on all NEMO meshes + error cases), `tools/run_ini_compare.sh`
(triangle meshes byte-identical to FESOM2's own mesh_part build). Test meshes in
`/work/ab0995/a270029/meshes/`: `nemo025q_fesom`, `nemo025mix` (dist_1/2/8/128), `nemo025q_cut`,
`nemo025mix_cut` (dist_1/2/4/8), `pi_w` (dist_1/2/4/8), `core2_w` (dist_1/2/4/8/128) — the NEMO ones
were made from the MATLAB output in `nemo025q/`; from the Python port, `elem2d.out`/`depth.out` and all
cut-out files are identical, `nod2d.out` differs in the one fold point and `nemo025mix` in the two
triangles of the quad at that point (its diagonal choice flips). `work/README.md` ("FESOM3 has no
partitioner") predates `fesom_ini` (the user's file, not edited).

### CORE2 throughput of the 4-slot code (login-node timing, 2026-10-09)

- **Setup:** old = `/work/ab0995/a270029/fesom3_hybrid/frozen_task6` (commit `5b2a7f3`, before the
  4-slot kernels, triangle results bit-identical to the pre-hybrid code), new = scratch build at
  `ec099ea`; both Release `-O3`. Original CORE2 mesh (`dist_16`, `dist_8`), production driver with the
  `job_levante` physics (linfs, JRA55, EVP, GM+Redi, TKE, splines, vinv, rvo_upwind 1, SW pene,
  w split), DEFAULT SSH tolerance, no output, no restart, `FESOM3_TIMING_EVERY=20`; cost per step =
  the model's loop timers (`LOOP TOTAL`) over steps 21..N from the cumulative reports (setup, the
  first 20 steps and I/O excluded). Old/new interleaved; np 16: 4 runs each of 200 steps, np 8: 3
  each of 120 steps. Scripts and logs: `/work/ab0995/a270029/fesom3_hybrid/t18_timing/`
  (`run_one.sh`, `series.sh`, `parse.py`, `t18_timing_np{16,8}.txt`, `log_*.txt`).
- **Result:** `LOOP TOTAL` ms/step — np 16: old min 456.5 / median 458.0, new min 492.8 / median
  499.1 → **new/old 1.079 (min) / 1.090 (median)**; np 8: old 850.9 / 853.9, new 926.9 / 931.9 →
  **1.089 / 1.091**. Noise: run-to-run spread of the loop time 2.5 % (old) / 3.3 % (new) at np 16,
  0.8 % / 1.7 % at np 8; 20-step blocks give the same ratio (1.08–1.10). Whole-run wall clock agrees
  (np 16 old 93.7–95.7 s, new 101.0–104.3 s).
- **Per component** (np 16 median, ms/step old → new): ice 88.4 → 113.8 (**×1.29**), GM/Redi
  36.1 → 39.0 (×1.08), dynamics u,v,w 75.5 → 82.2 (×1.09), mix/pres/EOS 52.0 → 55.1 (×1.06), tracer
  148.5 → 149.8 (×1.01), SSH 32.8 → 33.5, forcing/fluxes unchanged; np 8 the same pattern (ice ×1.29).
- **Hot spot (only component > 20 %): the sea-ice EVP.** `perf` over all 16 ranks (80 steps,
  `t18_perf_cmp.txt`): `evpdynamics` (stress_tensor/stress2rhs inlined, 120 subcycles over 4 slots)
  1.24× summed over ranks and 1.33× on the most ice-loaded rank, which sets the pace — the other ranks
  wait in MPI (`mca_btl_vader_component_progress` +19 %); ice FCT `ice_fem_fct` 1.28×, `ice_TG_rhs`
  1.62× (small in absolute terms). Ocean: `smooth_nod` now a separate routine (+0.7 % of the step),
  `fer_gamma2vel` 1.72×, `momentum_adv_vinv` 1.15×, `impl_vert_visc_ale` 1.19×. All compute outside
  MPI: ×1.048. Not optimised (Task 18 reports only).
- **Scaled:** 14 min per model year on 512 ranks (results17, pre-hybrid binary) × 1.08–1.09 ≈
  **15.1–15.3 min per model year**. Indicative only: login-node ratios at np 8–16 (shared node, MPI
  over shared memory); at 512 ranks communication weighs more and the ice-rank imbalance may differ.

### NEMO 0.25° example job: `work/job_levante_nemo025`

- New file next to `work/job_levante` (untracked like everything in `work/`; never `git add`ed).
  Same header style, modules (`env.sh intel`), limits and FESOM3_* physics as `job_levante` (copied);
  changed: mesh `/work/ab0995/a270029/meshes/nemo025q_fesom` (comment: `nemo025mix` as alternative),
  8 nodes × 128 = 1024 ranks (`NPART`, must match `dist_1024`), 1 h, output tree
  `results_nemo025/<mesh>_1958_1m`.
- **Setup part:** checks `nod2d.out`, `elem2d.out`, the depth file (picks `use_depthfile` 'aux3d' or
  'nemo'), the model binary, `SLURM_NTASKS == NPART`, and that the edge/level files and a complete
  `dist_<NPART>/` (rpart, format.out, NPART my_list + com_info) exist; if not, writes the
  `namelist.config` above into `work/ini_<mesh>_<NPART>/` and runs `srun -l -n 1 -N 1 fesom_ini`
  there, then re-checks. `DRYRUN=1 bash job_levante_nemo025` runs the checks and prints the would-be
  `namelist.config`, the fesom_ini and the srun lines and the FESOM3 environment, launches nothing;
  `DRYRUN=1 DRYRUN_INI=1` also runs fesom_ini on the login node (prepares `dist_N` before
  submission). Overrides `MESH`, `NPART`, `BIN`, `INIDIR` for tests.
- **Dry runs (2026-10-09):** as is → stops with `ERROR: missing build_intel_dp/bin/fesom_ini` (the
  user's binary of 2026-10-08 predates the hybrid branch; **rebuild `build_intel_dp` from this branch
  before use** — the old model binary cannot read quad meshes) and `dist_1024` missing → prints the
  fesom_ini call. With `BIN` = scratch build and a scratch mesh that has `dist_1024`: all checks ok,
  prints `srun -l -n 1024 …/fesom_lifecycle_native_mr`. With `NPART=8 DRYRUN_INI=1` on a scratch
  `nemo025mix`: fesom_ini ran (21 s), `check_dist.py` PASS afterwards.
- **USER: choose** (conservative placeholders): `FESOM3_STEP_PER_DAY=96` (900 s; ORCA025 cells
  shrink to ~10 km at high latitudes; the cut-outs ran 10 days stable at 48), run length 1 month
  (`FESOM3_RUN_LENGTH=1`, `FESOM3_RUN_UNIT=m`), nodes / wall time (8 nodes, 1 h), `RUNID` / output
  tree. Restart settings as `job_levante` (yearly; a shorter run writes its end-of-run checkpoint).
- **Wall-time ESTIMATE (not measured):** work per step relative to CORE2 from the mesh counts —
  nodes 7.2× (915 866 vs 126 858), wet node-levels 8.4× (31.1 M vs 3.7 M), element-levels 4.3×
  (quad mesh; mixed 5.3×), edges 4.9× (mixed 5.5×) → factor 4.3–8.4; × 2 steps (900 s) × 1.09 (4-slot
  cost) × 14 min per CORE2 year on 512 ranks → **2.2–4.3 h per model year on 512 ranks, 1.1–2.1 h on
  1024 ranks (8 nodes)**, assuming linear scaling (1024 ranks still hold 3.6× more nodes per rank than
  CORE2 on 512) → one month ≈ 6–11 min + startup (mesh read and IC for 0.9 M nodes, minutes) +
  fesom_ini 1.5 min on the first run; 1 h leaves ≥ 4× margin; a model year on 8 nodes: ask for ~3 h.
  Unknowns: more CG iterations of the SSH solve at 0.25° (CG ≈ 4 % of the CORE2 step), a larger ice share
  (dense ORCA grid in the Arctic; the EVP is the 4-slot hot spot), startup time at 1024 ranks.

### Gate migration for hybrid meshes (decision, 2026-10-08, plan `plans/completed/2026-10-08-hybrid-meshes.md` Task 7)

From Task 7 on the mesh geometry is the unified 4-slot one (area-weighted centroid, kite weights
`elem_w`, Green–Gauss `gradient_sca(8)` with x in 1:4 / y in 5:8, `edge_tri(2) == edge_tri(1)` and
`elem_neighbors == self` on the boundary). For triangles it equals FESOM2's P1 geometry only up to
roundoff, so no FESOM2-oracle byte gate (`max|Δ|=0`) can pass any more.
- **Retired** (kept in `tools/` as the record of how M0–M10 were proven; not to be run as gates):
  every FESOM2-oracle byte gate — `run_geom_gate{,_core2,_multirank}.sh` (+ `geom_diff.py` as a gate),
  `run_advhor_gate{,_multirank}.sh`, `run_pressure_gate{,_core2}.sh`, `run_step_gate{,_multirank}.sh`,
  `run_stepdyn_gate_multirank.sh`, `run_forcing_gate.sh`, `run_ic_gate_core2.sh`,
  `run_ice_gate_core2.sh`, `run_evp_gate_core2.sh`, `run_icefct_gate_core2.sh`,
  `run_iceflux_gate_core2.sh`, `run_icethermo_gate_core2.sh`, `run_tke_replay_gate_core2.sh`,
  `run_meshdiag_gate.sh`, and all `run_lifecycle_*gate*` scripts (`.sh` and `.sbatch`), plus the
  oracle-dump helpers they call (`run_*dump*.sh`, `run_oracle_pi.sh`, `run_kppinit_oracle.sh`,
  `run_tkedump_oracle.sh`).
- **Remapped:** `src/infra/mod_geom_dump.F90` still writes FESOM2's 6-row `gradient_sca` (x = slots
  1:3, y = slots 5:7), so `fesom_geomdump` + `tools/geom_diff.py --signal=<tol>` remain usable as a
  roundoff-level diagnostic. The other dump drivers build and run unchanged; Task 14 retires or
  remaps each of them.
- **Dump drivers (Task 14):** none deleted — each stays a buildable diagnostic dump whose FESOM2
  byte gate is retired (a status note heads every driver). `fesom_geomdump`/`mod_geom_dump`:
  REMAPPED to the FESOM2 layout (3 vertices, 6-row `gradient_sca`, `edge_tri(2) = 0` on boundary
  edges as FESOM2 stores it), triangle meshes only — a mesh with quads stops with a message.
  `fesom_meshdiagdump`: REMAPPED to the hybrid store layout (next bullet), checked by the ctests
  `meshdiag_{mixed,quad,pi_w}` (`tools/check_meshdiag.py`); `zarr_diff.py --meshdiag` vs FESOM2 is
  retired. `fesom_pressuredump`: REMAPPED — its own node→element means (KPP `viscAE`, the unstable-
  element count) run over the 4 slots with `elem_w` like the kernels. `fesom_advhordump{,_mr}`,
  `fesom_evpdump`, `fesom_icefctdump`: no triangle-specific code of their own (they run the
  converted kernels; the ice kernels convert in Tasks 15/16), unchanged. `fesom_conserve` is a kept
  ctest gate (its `helem` invariant was converted in Task 12).
- **Output format changes (Task 14) — check notebooks/scripts that read FESOM3 stores:**
  `fesom.mesh.diag.zarr` (`FESOM3_MESHDIAG`) now has 4-column per-face arrays, dim `n4`:
  `face_nodes` (a triangle's 4th = `_FillValue` −1, i.e. Zarr `fill_value` −1), the new
  `face_edges` (edge k joins nodes k, k+1; triangle 4th −1) and `face_links` (neighbour across edge
  k; missing neighbour AND a triangle's 4th = −999, FESOM2's value), `gradient_sca_x/y` (4
  coefficients, triangle 4th = 0); the topology names `face_edge_connectivity`/
  `face_face_connectivity`. xarray decodes a declared `_FillValue`, so with the default
  `mask_and_scale=True` these integer arrays come back as float with NaN (`mask_and_scale=False`
  keeps the integers). The element coordinates `lon`/`lat` of the field-output and restart stores
  (and the point where element vectors u/v/bolus are r2g-rotated) are now the area-weighted,
  cyclic-wrapped centroid (`elem_center`) instead of FESOM2's plain unwrapped 3-vertex mean: equal
  to < 1e-10 deg everywhere except the elements across the ROTATED dateline (pi 77, CORE2 566),
  where the old coordinate was on the wrong side of the globe and the old geographic u/v were
  rotated with the wrong angle. `run_tri_regression.sh` compares element vectors there by speed
  (see its header). The user's notebooks in `work/` read FESOM2's `fesom.mesh.diag.nc` — unaffected.
- **Kept** (FESOM3-only, unaffected): `run_restart_gate_core2.sh`, `run_restart_gate_multirank.sh`
  (split vs straight), `run_output_gate.sh` (writer vs generator), `run_conserve_pi.sh` and the
  ctest suite.
- **Replacement:** `tools/run_tri_regression.sh` (Task 7, tightened in the Task-7 follow-up) — run
  after Tasks 7, 10, 14 and 16. Reference = `/work/ab0995/a270029/fesom3_hybrid/frozen_task6_soltol`:
  commit `5b2a7f3` (the last bit-identical state) + ONLY the SSH-solver tolerance knob (README there;
  with the knob unset it is byte-identical to the original `frozen_task6`, which is kept). Cases:
  `prod:pi:1`, `prod:pi:2`, `prod:core2:8` (production driver, job_levante physics, 4 steps, all
  substep node dumps in double + float32 output every step, ref vs new), `decomp:pi_w:2`,
  `decomp:pi_w:8` (new build np 1 vs np N, both restarted from one 0-step np-1 checkpoint — see the
  next bullet), `conserve:pi:1/2` (`fesom_conserve` zstar+Redi+TKE 20 steps, drift and final totals
  ≤ 1e-12).
  - **Solver:** every prod/decomp run sets `FESOM3_SSH_SOLTOL=1e-14` (new env knob in
    `oce_ssh_solve.F90`, also `FESOM3_SSH_MAXITER`; unset = production soltol 1e-5 / maxiter 2000;
    ~3× CG iterations). With 1e-5 the iterate differences swamp the roundoff: the same binary at 1e-5
    vs 1e-14 differs by 3.5e-4 in `d_eta` at step 1 on pi — that is why the first gate could only be 1e-2.
  - **Metric** per field `d = max|a−b|/max|a|` over nodes/elements, levels, steps/records; tiers and
    classes (pi / core2 tolerance): PRE = step 1 before the first solve (1e-9 / 1e-9; measured 6e-15 /
    1.3e-11); EARLY = every double dump except Kv over steps 1..3 on pi, step 1 on core2 (1e-10 / 1e-8,
    core2 raised to 1e-7 in Task 15 for the EVP-amplified roundoff; measured 9.2e-13 / 2.1e-9); OCEAN =
    all steps of the double dumps (T, S, eta, d_eta, w, ssh_rhs,
    density, …) and all float32 output records of temp/salt/ssh/sst/sss/u/v/w/unod/vnod/bolus
    (1e-6 / 2e-4, float32 +2^-21; measured 3.3e-8 dumps, 1.0e-7 u/v / 4.4e-6 dumps, 2.2e-5 u/v); ICE =
    a_ice/m_ice/m_snow (1e-5 / 2e-3; measured 3.9e-7 / 1.2e-4); MIX = Kv, Av (1e-6 / 1e-2; measured
    2.2e-9 / 2.4e-4). decomp measured: PRE 0, EARLY ≤ 4.3e-13, OCEAN ≤ 2.5e-11 (dumps).
  - **Why not tighter:** the remaining floor is amplification of roundoff by stiff/discontinuous
    physics, not roundoff: on CORE2 the sea ice (EVP at rest, Δ_min-regularised viscosity) lifts T at
    ice-covered nodes to 2e-9 already in step 1 and to ~1e-6..1e-5 after 4 steps; on pi the jump
    comes at step 4 (first JRA55 record refresh moves the ice; the same jump appears in np1 vs npN
    from a common state and with the summer IC); TKE thresholds drive Kv/Av. Steps 1–3 on pi are
    pure roundoff (≤ 1e-12), hence the EARLY tier.
  - **Levels at the end of the plan** (Task 19, clean clone of `98b5164`, all 4-slot kernels in): pi
    PRE 6.3e-15, EARLY 9.2e-13, OCEAN 1.1e-7, ICE 5.9e-7, MIX 4.4e-9; CORE2 PRE 2.5e-11, EARLY 2.2e-8
    (T at ice nodes; with mEVP in both runs 1.1e-15), OCEAN 3.1e-5, ICE 9.7e-5, MIX 7.8e-4; decomp
    OCEAN ≤ 8.8e-10; conservation drift ≤ 2.8e-14. Extra case `decomp:core2_w:8` (np 1 vs np 8 on full
    CORE2 from one checkpoint): EARLY 2.6e-12, OCEAN 2.7e-5, ICE 7.3e-5, MIX 2.2e-5. Element vectors at
    the moved dateline centroids are compared by speed (Task 14, header of the script).
- **Decomposition dependence on the original pi (resolved, by design, not a bug):** a COLD START
  differs between np 1 and np N (frozen and pre-hybrid 9b0d3a6 binaries alike; density 4.7e-3 at step
  1, pi_w too: 5.5e-3) because the INITIAL CONDITION does. `oce_initial_state::extrap_nod3D` (FESOM2
  `gen_support.F90` extrap_nod3D, transcribed for byte-identity) fills land-masked PHC points by a
  Gauss–Seidel neighbour average over OWNED nodes, reading halo values that are refreshed only after
  each full sweep — the result depends on the partition and node order. Proven on pi: all 513
  (node, level) density differences at step 1 (34 nodes) are points that `getcoeffld` leaves as dummy
  (2175 of 3140 nodes need extrapolation at some level); nothing else differs. From a common state
  (one checkpoint written by a 0-step np-1 cold start, restarted at np 1/2/8, tight solver) the model
  is decomposition-independent to roundoff (original pi and pi_w, frozen and current build: step 1
  ≤ 1.2e-14, steps 2–3 ≤ 2e-13). Not changed: a partition-independent extrapolation (Jacobi sweep,
  or extrapolate on the global mesh) would change every production IC and break the FESOM2 byte
  heritage — a design decision for the user. Consequence: compare np 1 vs np N only from a common
  restart (as `decomp:` does), never from two cold starts.

### Known findings / open questions for the user

Collected from the plan's ⚠️ notes and the `[deviation]` lines of the task log. None is fixed; the ones
marked **decide** need a design decision.

1. **Decomposition-dependent cold start** (Task 7 follow-up, L60) — **decide**: `extrap_nod3D` (FESOM2's
   Gauss–Seidel IC extrapolation with per-rank halo snapshots) makes np 1 and np N start from different
   T/S at the extrapolated points; from one checkpoint the model is decomposition-independent to
   roundoff. A partition-independent sweep would change every production IC (details in "Gate
   migration for hybrid meshes").
2. **Linear tracers are not advected exactly at triangle–quad interface nodes and on general quads**
   (Task 11) — **decide**: the median-dual edge scheme is linearly exact only where the dual cell is all
   triangles or all parallelograms; elsewhere the node tendency of a linear tracer in uniform flow has an
   O(1) LOCAL relative error (0.22 NEMO-mixed interface nodes, 0.24 jittered patch) although the MUSCL
   face values are exact. Possible remedy: face values at the dual-segment midpoints.
3. **Vorticity of a linear shear** is exact only on triangles and parallelograms (Task 10): 6e-2 (first
   order) on the jittered general quads; the NEMO quads are lon/lat rectangles (exact).
4. **Tracer checkerboard / grid-scale energy grows** (Task 17 ⚠️) — **decide**: 10-day np-1 runs on the
   cut-outs, E_loc SST ×6, SSS ×11–12, faster at quad and interface nodes than at triangle nodes; not
   exponential within 10 days but not saturating either. The SSH null mode does not grow. Remedies
   (hourglass control, the dual-segment face values of item 2) are design decisions. Watch it on the
   full mesh (`tools/cutout_diag.py checker` works on any mesh with output).
5. **GM/Redi resolution scaling** `(area(n)·2/refscalresol²)^(1/order)` in `init_Redi_GM`
   (`src/oce/oce_fer_gm.F90:103`, `:126`) — **decide**: the factor 2 is the triangle element/node area
   ratio; on all-quad meshes it overstates the resolution by √2 (larger `fer_scal`). FESOM2 tuning, left
   unchanged.
   - **Same question for the viscosity length on quads** (code review 2026-10-09) — **decide**:
     `visc_filt_bidiff` (`src/oce/oce_dyn_visc.F90:132`, `:168`) uses `len =
     sqrt(sum(elem_area(el)))` of the edge's elements; a quad has about twice a triangle's area at
     the same node spacing, so `len` is ≈ √2 larger, and the `/elem_area` divisor and the 4-side
     stencil change the effective coefficient too. On all-quad and mixed meshes the effective
     viscosity therefore differs by an O(1) factor from what the triangle-tuned `visc_gamma0/1/2`
     give. Formula left unchanged (physics tuning); an element-type-independent resolution length
     would be the alternative.
6. **Output element coordinates changed at the rotated-dateline elements** (Task 14): pi 77, CORE2 566
   elements now carry the wrapped centroid; their lon/lat and the geographic u/v/bolus_u/bolus_v change
   by O(1) (the old values were wrong). Notebooks comparing old and new output see this.
7. **Standard EVP is unstable for uniform ice at rest** (Task 15, pre-existing FESOM2 property): roundoff
   grows to O(1e-2) m/s within 120 subcycles on the triangle patch; it is the amplifier behind the CORE2
   ICE regression tier. mEVP (zero stress at rest) is the clean diagnostic.
8. **TKE thresholds amplify roundoff** (Task 17): 1-vs-N on the quad cut-out reaches 1e-4 (v) after 5
   days, and a threshold can turn roundoff into a local O(1e-9) jump within hours depending on build,
   np and dump probes; with PP instead of TKE ≤ 1.1e-7. Compare decompositions with PP, or expect this.
9. **Open MPI 4.1.2 knem stale payloads on the login node** (Task 10, L61): library bug, workaround
   `--mca btl_vader_single_copy_mechanism none` in the test launchers; batch runs (UCX, `env.sh`) are not
   affected (see "Open notes / risks").
10. **Cost of the 4-slot code**: CORE2 ×1.08–1.09 per step (≈ 15.1–15.3 min per model year on 512
    ranks), the sea-ice EVP ×1.29 is the hot spot; not optimised.
11. **At most `MAX_NEIGHBOR_PARTITIONS = 32` neighbour ranks per communication list**
    (`src/types/mod_partit.F90:12`, `mesh_part/src/MOD_PARTIT.F90`; `rPE`/`sPE` hold 32, `rptr`/`sptr` 33
    in both): `fesom_ini` stops ("Increase MAX_NEIGHBOR_PARTITIONS") for all five lists when more than 32
    receive or send PEs are counted, before the fixed arrays are filled (`gen_comm.F90` for com_nod2D,
    com_elem2D, com_elem2D_full; `set_pe_ptr` in `oce_local.F90` for com_nod2D_full, com_edge2D; no ctest,
    it needs a > 32-neighbour partition); the model's `read_com` dies with the list name and count before
    reading (ctest `test_partit_many_neighbours_*`). ORCA025 at 128 and 1024 parts: at most 16. Count before going
    to very high N (`max_neighbours.py`, setup workflow step 6).
12. **Full NEMO 0.25° cold start** (Task 19, production physics, 900 s, np 8, 2 days): finite, but
    max|eta| 2.55 / 2.52 m and max|uv| 2.42 / 2.60 m/s (quad / mixed mesh) within the first day, bounded
    afterwards (day 2: 1.8 m, 2.1–2.2 m/s) — global cold-start adjustment on a mesh never run before; to
    be judged with the user's longer runs.
13. **`configure.sh` now passes `CC`/`CXX`** (`-DCMAKE_C_COMPILER`/`-DCMAKE_CXX_COMPILER`, for METIS and
    `fort_part.c` of `fesom_ini`); builds that bypass `configure.sh` need a C compiler configured.
    **Rebuild `build_intel_dp` from this branch** before any NEMO run (the old binary cannot read quad
    meshes and has no `fesom_ini`); a rebuild moves a running triangle experiment at roundoff level —
    keep a copy of its binary if it must continue unchanged.
14. **Old reference partitions not reproducible** (Task 3): pi `dist_2/8` + its `elvls/nlvls.out` and
    CORE2 `dist_8` + its current levels come from older fesom_ini versions; the ported fesom_ini equals
    FESOM2's own current build byte for byte, and the old files pass `check_dist.py`.
15. **Still open from before the plan — dotty `w` near topography** in the user's CORE2 run: not
    explained by the MUSCL change nor by the knem bug (production uses UCX); candidate one-change tests
    from the same restart (user's runs): `FESOM3_MOMADV_OPT=2`, `FESOM3_WSPLIT` unset,
    `FESOM3_RVO_UPWIND=0`.

## Next task

### ✅ RESTART / CHECKPOINT — **COMPLETE & tagged `m10` (2026-06-29).** Partition-independent Zarr checkpoints (full prognostic state incl. sea ice + EVP sigma), written periodically + at end-of-run, read at startup on `r_restart`, resuming byte-exactly. Built on the M9 Zarr stack (`mod_io_zarr` + `mod_io_decomp` + new `decomp_gather` inverse); new `src/io/mod_io_restart.F90` (writer+reader, atomic tmp→rename finalize + `restart.latest` + keep-N prune), `mod_io_posix.F90` (bind(C) rename/unlink/fsync/rmtree — no `execute_command_line`), `mod_io_coords.F90` (shared coord helper). Lifecycle wiring in `fesom_lifecycle_native_mr` (detect/read/AB-guard/write/cadence + `clock_finish`).
- **Gate (6.1) split-vs-straight self-consistency:** np=1 mid+boundary **exactly `max|Δ|=0`** over the WHOLE state (ice+sigma+velocity+AB+tracers); np=2 mid exact; np=2 boundary **≤1 ULP** (temp 2.2e-16) — a **FESOM2-inherent** limit (any-node element ownership → 562 redundantly-owned boundary elements not halo-synced; FESOM2's own nc restart dedups identically; a live sync would break M0–M9 byte-identity → user-accepted caveat; `tools/run_restart_gate_core2.sh`, `--rel-floor` 0 at np=1 / 1e-12 at np>1). **Resolved 2026-10-12:** the same unsynced copies leaked volume in CORE2 runs; the shared-element sync (section "Shared elements: one owner for `uv`") keeps them equal, np>1 restarts are now exact and the gate's floor is 0 (1e-12 only with `FESOM3_SHARED_SYNC=0`). Found+fixed 2 missing carried fields the gate caught: `ice%thermo%t_skin` + `dyn%d_eta`.
- **Gate (6.2) cross-partition storage+restore:** `C2(np2) ≡ C8(np8)` **strict `max|Δ|=0`** all 27 stores (`tools/run_restart_gate_multirank.sh`, np8 reads C2, zero steps, re-writes C8).
- **No-regression (6.3):** ctest **Intel 21/21 + GNU 21/21** (sweep caught & fixed a GNU-only cpp `/*`-in-comment build break); M9 zarrsmoke/meshdiag/output gates GREEN; production MR lifecycle byte-gate np2 **195 records `max|Δ|=0`**.
- **Plan:** `docs/plans/completed/2026-06-28-restart-checkpoint.md` (all boxes ticked; Claim-scope refinement documents the np>1 ULP). Memory: `restart-tskin-carried-ice-state`, `restart-np-element-1ulp-inherent`.
- **NEXT MILESTONE: OPEN — FESOM3 v1 is feature-complete & byte-exact vs FESOM2 (`m0…m10`: dyn core, sea ice EVP/mEVP, GM/Redi, KPP, zstar, TKE, long runs, output, restart).** No more core physics to byte-port. The genuine remaining directions (each a phase, user picks): **(a) scientific/production validation** — multi-decade real runs (JRA55-do from 1958), SST/SSS RMS + deep-drift + ice area/volume vs obs & FESOM2 (the "is it a good climate run" check the byte-gates do NOT cover; Post-Completion manual HPC); **(b) mixed precision** — flip WP→real32, protect cancellation islands (foundation laid); **(c) performance/SYPD profiling + optimization** (forcing-perf was a first taste); **(d) beyond-v1 physics** — IDEMIX, backscatter (uke/iwe), BGC/REcoM, ice-shelf cavities, icebergs, coupling (OIFS/OASIS/icepack) — fresh ports, not byte-ports. Decide with the user.

### ✅ M9 (Zarr OUTPUT) — **COMPLETE & tagged `m9` (commit `b631399`, 2026-06-28).** All stages (scalars + 3-D + NODE + ELEMENT vectors/scalars) + Task 2.6 namelist/knobs + F1 full no-regression sweep, all byte-gated. The lone open plan box = the MANUAL `ushow` display smoke. (Restart was the milestone that followed — now DONE, see above.)

**Plan:** `docs/plans/completed/2026-06-28-m9-zarr-output.md` (checkboxes ticked through F1; moved to completed/). Design SETTLED.

**DONE (F1 done, all byte-gated, ctest 18/18 Intel+GNU, no regression):**
- **Stage 0 — `mod_io_zarr`** (Zarr v2): JSON + C-order transpose + partial-chunk fill-pad + codecs `none`/**lz4** +
  `.zmetadata`; `zarr_write_chunk` 1d/2d/**3d**/int + `zarr_rewrite_zarray` (grow time dim). Gate `run_zarrsmoke.sh`.
- **Stage 1 — `mod_io_decomp`** (canonical `MPI_Alltoallv` → writer subset, no rank-0 gather; ctest `test_io_decomp`
  np 1/2/8) + **`mod_io_meshdiag`** + `fesom_meshdiagdump` → **`fesom.mesh.diag.zarr` 18 vars `max|Δ|=0` vs FESOM2**,
  partition-indep (`run_meshdiag_gate.sh` np 1/2/8); wired into lifecycle (`FESOM3_MESHDIAG`).
- **Stage 2 SCALARS — `mod_io_means`** (Tasks 2.1/2.2/2.3/2.4): per-variable-per-year stores
  `<dir>/<name>.fesom.<YYYY>.zarr` with growing time dim (chunk_time=1 append), CF `time`+calendar, embedded
  `lon`/`lat`(+`nz`/`nz1` for 3-D). **2-D node scalars** ssh/sst/sss/a_ice/m_ice/m_snow + **3-D** temp/salt(nz1)/w(nz),
  **snapshot AND mean** (accumulate+divide in OUTPUT precision = FESOM2 io_meandata semantics), 3-D below-bottom
  `NC_FILL` mask (nlevels-based). Wired into lifecycle (`FESOM3_OUTPUT`/`FESOM3_OUTPUT_EVERY`). Gate
  `run_output_gate.sh` (`fesom_outputsmoke` + `zarr_diff.py --output/--output-cmp`): synthetic fields == generator
  formula `max|Δ|=0`, **partition-indep dist_2≡dist_8≡1-rank**, mean divide + 3-D mask verified. No FESOM2 oracle
  (state byte-exact thru M8 + writer self-consistency).
- **Stage 2 VECTORS — Task 2.5 (`unod`/`vnod` + r2g rotation):** ported `vector_r2g` VERBATIM (FESOM2
  `gen_modules_rotate_grid.F90:164-202`) into `mod_mesh_rotate.F90` (rotated→geo, inverse of the byte-gated
  `vector_g2r`); `means_define_vector3d` pairs two node3d fields, rotated together at write. **NODE-based**
  `unod`/`vnod` (`dyn%uvnode`) — FESOM2's actual default velocity output (every `namelist.io` uses it, NEVER elem
  `u`/`v`); reuses `means_define_node3d`, no new decomp. **FESOM2 ORDER matched:** `io_r2g` rotates the accumulated
  SUM (`io_meandata.F90:2265`) BEFORE `compute_means` divides (:2335) → `write_vector_3d` rotates-sum-then-divides
  (`<f4` promote→rotate→demote→divide = io_r2g r4 branch). `vec_frame` geographic(default)|native = FESOM2
  `vec_autorotate`, via `FESOM3_VEC_FRAME`. Gates: round-trip ctest `test_vector_rotate` (`g2r∘r2g==id` 8.9e-15);
  `run_output_gate.sh` native==raw `max|Δ|=0`, geographic==INDEPENDENT numpy `vector_r2g` ref `≈2e-13`≪1e-9
  (non-vacuous), partition-indep both frames. ctest 16→17.
- **Task 2.6 (namelist.io knobs + chunk_time/chunk_vert/lz4/n_writers/filesplit + per-field cadence):** writer is
  now CONFIGURABLE. `mod_io_zarr` gained `zarr_read_chunk` (1d/2d/3d + lz4 `LZ4_decompress_safe`) for the
  `chunk_time>1` **read-modify-write** append (record t -> time-chunk t/ct slot mod(t,ct)). `mod_io_means`
  REFACTORED to PER-FIELD record tracking: each field keeps its own `t` + period store; `means_begin`/`means_write`
  replaced by **`means_output(io, istep, t_means_clock)`** evaluating the FESOM2-ported events (`event_due`:
  annual/monthly/daily/hourly/step from `gen_events.F90`) per field. Global knobs (chunk_time, chunk_vert
  [ceil(nz/cv) vert chunks], compressor [data-array lz4], filesplit [y|m -> `<name>.fesom.<YYYY>[_<MM>]`]) via
  `means_init` args + `FESOM3_*` env. **`means_read_namelist`** parses `&nml_general` + `&nml_list` (FESOM2
  derived-type-array trick: 1 `t_io_entry` per 5 flat values, +5th `op` mean|snap). Lifecycle now LIST-DRIVEN
  (`register_output_var` dispatch by id; env-fallback default set kept; accumulate EVERY step via `means_has`).
  GATES: `run_output_gate.sh` knob sweep — chunk_time={1,2,3} + chunk_vert + lz4 + n_writers=2 all value-identical
  to base `max|Δ|=0`, combo partition-indep np1≡np2, per-field freq=2 (`fld_f2`) => floor(nrec/2) records;
  `test_io_means` ctest pins the namelist parse Intel+GNU (17→18). No regression (zarrsmoke + meshdiag GREEN; GNU
  writer round-trips, pure-data fields bit-identical to Intel, only rotation-derived lon/lat/unod/vnod differ
  ~1e-13 = expected ifort/gfortran libm).
- **Task 2.7 (ELEMENT output — user-requested):** added the element decomp `De` + owned elem-centroid coords to
  `mod_io_means`: ROTATED centroid `sum(coord_nod2D(1:2,elem2D_nodes(1:3,e)))/3` (FESOM2 io_r2g:3004 VERBATIM,
  flag_coord=0) for the vector r2g + the GEOGRAPHIC centroid `r2g(rotated)` (deg) for the embed + `mesh%nlevels`
  for the element mask. New `means_define_elem2d`/`means_define_elem3d`/`means_define_vector3d_elem`; the field
  registration refactored to a shared `add_field` core + node/elem wrappers; the write path resolves the entity
  context (`D`/`nO`/`nlevown`/`rlon`/`rlat`) once in `write_one_field` (pointer into target `io`) + threads it
  through the writers/`emit_chunks`/`put_static`/`open_field_store`; `def_field_store` uses per-field `hdim`
  ('nod2'|'elem'). Lifecycle registers **u/v** (`dyn%uv`), **Av** (`dyn%work%Av` full levels), **bolus_u/v**
  (`dyn%fer_uv`, Fer_GM-guarded) via the list dispatch + default set + `means_has`-guarded accumulate;
  `config/namelist.io` gains u/v/Av rows. GATE (`run_output_gate.sh`, np 1/2/8): elem `native==raw` `max|Δ|=0`,
  `geographic ==` numpy elem-centroid r2g ref (`≈4.8e-13`≪1e-9), elem scalars `fld_e2`/`fld_e3` `max|Δ|=0`,
  partition-indep `dist_2≡dist_8≡1` `max|Δ|=0`, element nlevels mask fires. ctest 18/18 Intel+GNU; no regression.

**KEY LESSONS (reuse for F1):**
- Writer path (decomp_redistribute → zarr_write_chunk; store-create ordering rank0-define+barrier+writers-write) is
  PROVEN — reuse verbatim. POSIX `mkdir` (NOT `execute_command_line`; fork post-MPI_Init SEGFAULTS).
- xarray: `decode_times=False` to compare raw `time` seconds (else →datetime64); `mask_and_scale=False` for exact;
  3-D fill = `NC_FILL_DOUBLE` 9.9692099683868690e36 (=FESOM2), xarray masks →NaN. 2-D scalars use `fill_value:null`.
- mean = `means_accumulate` every step (sum) / `means_write` divide-by-count + reset, IN output precision (FESOM2
  io_meandata.F90:2107/2335). Snapshot = the count=1 special case.
- **Deferred:** `face_edges`/`face_links`/`gradient_vec` (FESOM3 never builds the sources). FESOM2 oracle runs SLOW.

**VERIFY current state (run first to confirm green before building on it):**
```
cd /home/a/a270088/fesom3 && source env.sh intel
cmake build_intel_dp && cmake --build build_intel_dp -j 8
(cd build_intel_dp && ctest)            # expect 17/17 (incl. test_vector_rotate)
tools/run_zarrsmoke.sh                  # Stage 0
tools/run_meshdiag_gate.sh 1            # Stage 1 (also: 2, 8)
tools/run_output_gate.sh                # Stage 2: np 1 2 8, scalars+3D+VECTORS, geographic+native, partition-indep
```
**Working tree is UNCOMMITTED** (the user commits): new `src/io/mod_io_{zarr,decomp,meshdiag,means}.F90`,
`src/drivers/fesom_{zarrsmoke,meshdiagdump,outputsmoke}.F90`, `test/test_io_decomp.F90`, `test/test_vector_rotate.F90`,
`test/test_io_means.F90` (namelist parser), `tools/{run_zarrsmoke,run_meshdiag_gate,run_output_gate}.sh`,
`tools/zarr_diff.py`, `config/namelist.io`; modified `CMakeLists.txt`, `test/CMakeLists.txt`,
`src/mesh/mod_mesh_rotate.F90` (vector_r2g), `src/io/mod_io_{zarr,means}.F90` (Task 2.6: read-chunk + per-field
events/knobs/namelist), `src/drivers/fesom_lifecycle_native_mr.F90` (list-driven output), plan, HANDOFF.

**Task 2.6 (namelist.io knobs + chunk_time/chunk_vert/lz4/n_writers/filesplit + per-field cadence) ✅ DONE** — see
the Task 2.6 DONE bullet above + the plan. New/modified: `mod_io_zarr` (`zarr_read_chunk`), `mod_io_means` (per-field
`means_output`/events/RMW/`means_read_namelist`), `fesom_lifecycle_native_mr` (list-driven), `fesom_outputsmoke` +
`tools/{zarr_diff.py,run_output_gate.sh}`, `test/test_io_means.F90` (+`test/CMakeLists.txt`), `config/namelist.io`.

**Task 2.7 (ELEMENT output) DONE 2026-06-28** — see the Task 2.7 DONE bullet above + the plan (checkboxes ticked).
New/modified: `mod_io_means` (element decomp `De` + elem-centroid coords + `means_define_elem2d/3d` +
`means_define_vector3d_elem` + the shared `add_field` core + entity-context-resolving write path), `fesom_lifecycle_native_mr`
(u/v/Av/bolus_u/v registration + accumulate), `fesom_outputsmoke` + `tools/{zarr_diff.py,run_output_gate.sh}` (elem
gate: native==raw + geographic==numpy elem-centroid r2g + partition-indep, all `max|Δ|=0`/≈4.8e-13), `config/namelist.io`.

**F1 (full no-regression + gate sweep) DONE 2026-06-28:** ctest 18/18 Intel+GNU; zarrsmoke + meshdiag + output gate
(np 1/2/8, native+geographic, Task 2.6 knobs, ELEMENT) all GREEN; production MR lifecycle byte-gate
(`run_lifecycle_fullynative_gate_multirank.sh`, CORE2 forced np=2) `max|Δ|=0` **both whichEVP** (130 records each,
worst |Δ|=0) — element-output wiring proven inert/additive; **real-lifecycle output integration DONE** — the same np=2
run with `FESOM3_OUTPUT`+`FESOM3_OUTPUT_EVERY=1` wrote all 14 default-set stores from the REAL loop (node
ssh/sst/sss/a_ice/m_ice/m_snow/temp/salt/w/unod/vnod **+ element u/v/Av**), xarray opens each (dims, CF time, finite
node+elem-centroid coords, 2 records, 244659 elem ≈ 2× node), bolus correctly absent (Fer_GM off).

**M9 DONE & tagged `m9`** (commit `b631399`, 2026-06-28): Task 2.6 + 2.7 + F1 all byte-gated; the working tree
(`mod_io_zarr`/`mod_io_means`/both drivers/`config/namelist.io`/`tools{zarr_diff.py,run_output_gate.sh}`/
`test{test_io_means.F90,CMakeLists.txt}` + docs) committed, plan moved to `docs/plans/completed/`, tag created. The
only open M9 item is the MANUAL `ushow` display smoke (the xarray/zarr round-trip is the automated proxy).

**RESTART / CHECKPOINT = DONE & tagged `m10`** (2026-06-29) — see the `✅ RESTART / CHECKPOINT` section at the top of
"Next task" for the full result. The brainstorm-prep + plan are at `docs/plans/2026-06-28-restart-handoff.md` and
`docs/plans/completed/2026-06-28-restart-checkpoint.md` (forks decided: F-A/F-B = canonical partition-independent Zarr
reusing the M9 stack; F-C = serialize EVP `sigma` for true `max|Δ|=0` on ice; cadence periodic+end; clock_finish
ported; field set = the FESOM2 superset + the gate-found carried `t_skin`/`d_eta`). Gates GREEN (6.1 split-vs-straight
np=1 exact / np>1 ≤1-ULP FESOM2-inherent; 6.2 cross-np restore `C2≡C8`; 6.3 no-regression Intel+GNU ctest 21/21 +
M9 gates + production lifecycle). **NEXT MILESTONE: OPEN — v1 is feature-complete & byte-exact (`m0…m10`); no core
physics left to port.** Remaining directions (each a phase; see the top of "Next task" for the full list): scientific
multi-decade validation, mixed-precision (WP→real32), performance/SYPD, or beyond-v1 physics (IDEMIX/backscatter/BGC/
cavities/coupling). Decide with the user.

M0–M10 COMPLETE & tagged (m0…m10); forcing-perf ✅ RESOLVED (`708fdf4`/`ce2288a`/`1ef9516`; L133).

**M9 = model OUTPUT as hand-rolled Zarr v2** — xarray-readable, openable by the user's `ushow` (`/home/a/a270088/ushow`),
UGRID-1.0, + a `fesom.mesh.diag.zarr` analog. **Settled design:**
- **Write engine:** canonical global-id ordering (partition-INDEPENDENT files) + uniform chunks (configurable
  `time/vert/horiz`) + **distributed chunk-writers** (configurable `n_writers` subset; one `MPI_Alltoallv` redistributes
  compute→canonical-chunk layout, level-by-level for 3D; **NO single-rank gather**) + pluggable codec (`none` default,
  **lz4** first). Hand-rolled (NCZarr exists on Levante but we chose control/parallelism). The user explicitly rejected
  BOTH gather-to-root AND partition-dependent files → the distributed-writer redistribution is the resolution.
- **4 modules** (`src/io/mod_io_*`, dwarf-aligned, auto-GLOB-built): `mod_io_zarr` (v2 writer: JSON + C-order chunks +
  codec), `mod_io_decomp` (the Alltoallv redistribution + writer subset), `mod_io_meshdiag` (UGRID mesh.diag),
  `mod_io_means` (def_stream registry + accumulate/average + interval write). Driver gets 2 hooks (meshdiag at setup;
  `output(istep)` **after `step_oce` ~line 730, ABOVE the `if(.not.step_diag) cycle` at line 740**); physics untouched.
- **On disk:** per-variable-per-year stores (`temp.fesom.1964.zarr`, … — pyfesom/ushow-compatible) + separate
  `fesom.mesh.diag.zarr`; coords embedded, connectivity in mesh.diag; yearly split (configurable).
- **Fields (v1):** node `T,S,ssh,sst,sss,a_ice,m_ice,m_snow` + element `u,v` + node `w`; mean+snapshot per-var; per-var
  freq; float32 default; vectors default **geographic** (must **port `vector_r2g`** — only `vector_g2r` exists), knob to native.
- **Gates:** (1) mesh.diag vs FESOM2 `fesom.mesh.diag.nc`; (2) field values vs FESOM2 output **chasing `max|Δ|=0`** via
  byte-aligned averaging; (3) partition-independence `dist_2 ≡ dist_8`; (4) python round-trip; (5) ushow smoke. Python
  gate env: `/work/ab0995/a270088/mambaforge/bin/python3` (xarray/zarr/numpy).

**Build order:** Task 0.1 (`mod_io_zarr` scaffold) → 0.2 (chunk write + round-trip) → 0.3 (lz4) → Stage 1 (mesh.diag
1.2a/1.2b → MR 1.3) → Stage 2 (fields 2.1→2.6). **1-rank gate FIRST then multi-rank, every stage** (optional-`partit`).
Plan-review caught: `gradient_vec` is unbuilt in F3 (deferred from mesh.diag); no invertible `vector_r2g` (port it from
FESOM2 `gen_modules_rotate_grid.F90`); the line-740 `cycle` would skip output (hook goes above it); CF `time:calendar`
needed (L49-sensitive); pin `chunk_time=1` for v1.

**Memory:** [[project-fesom3-implementation-state]], [[stay-close-to-fortran]], [[project-levante-mpi-knem-gotcha]],
[[project-forcing-is-jra55-not-core2]] (calendar sensitivity for the CF `time` coord).

---

**M2.10 forcing DONE + M2.11a geometry DONE + M2.11b initial conditions DONE + M2.11c lifecycle DONE (unforced AND
forced) — the multi-step CORE2 lifecycle is `max|Δ|=0` (the "CG floor" was SOLVED 2026-06-21; the forced re-verify
CONFIRMED 2026-06-21: 195 records, worst |Δ|=0).** **All of M2 is now byte-exact end-to-end on CORE2.**

**→ M2.12 (multi-rank, MVP exit) IN PROGRESS.** Decomposed into M2.12a/b/c in the plan
(`docs/plans/2026-06-18-fesom3-architecture.md` Task M2.12) — see "M2.12 entry notes" below. **Key finding: the
multi-rank foundation ALREADY EXISTS** (`mod_halo.F90` real MPI exchange_nod/elem + `mod_partitioning` myList/com,
tested 1/2/8-rank), so M2.12 is faithful transcription of FESOM2's `read_mesh`/`find_neighbors`/`mesh_areas` (with
their halo exchanges) gated PER-RANK vs same-partition FESOM2 — not a from-scratch parallel build.
- **M2.12a — local-mesh remap + per-rank GEOMETRY byte-gate — ✅ DONE (2026-06-21).** `max|Δ|=0` on ALL 19 geometry
  fields, every rank, on pi **dist_2 + dist_8** (`tools/run_geom_gate_multirank.sh`). `read_mesh_local` (global→local
  scatter + owned `nod_in_elem2D` + `enforce_cw`(owned)) + partition-aware `compute_geometry` (owned centers
  `exchange_elem`'d for owned-edge `edge_cross_dxdy`; owned-node areas local). See LESSONS L30.
- **M2.12b — multi-rank TRACER ADVECTION byte-gate — ✅ DONE (2026-06-21).** `max|Δ|=0` on ALL 13 advection fields,
  every rank, pi **dist_2 + dist_8** vs same-partition FESOM2 (`tools/run_advhor_gate_multirank.sh`). Lifted the
  advection subtree (init_tracers_AB → tracer_gradient_elements → muscl_adv_init/find_up_downwind_triangles/
  fill_up_dn_grad → adv_tra_hor/ver → oce_tra_adv_fct → do_oce_adv_tra → advect_tracer) to multi-rank via an
  **OPTIONAL `partit`** (absent ⇒ the proven 1-rank path VERBATIM; present+npes>1 ⇒ owned/halo bounds + the FESOM2
  exchanges). The explicit-shape dummies + per-element arithmetic are UNCHANGED so codegen — and the byte-match incl.
  any vectorised divide (L29) — is preserved; only loop bounds change + guarded exchanges are added. New machinery:
  `src/infra/mod_part_bounds.F90` (`owned_bounds`/`is_multirank`); the `find_neighbors` halo dance
  (`complete_nod_in_elem_halo` in `read_mesh_local` — fills halo-node `nod_in_elem2D`, re-localised through eXDim);
  `exchange_elem_full` (2D-r/2D-i/3D-r on `com_elem2D_full`) + `core_blk_r` in `mod_halo`; `compute_geometry` halo
  exchanges (`elem_area` full + `area`/`areasvol`/`*_inv`); `exchange_elem(tr_xy)` (full) + `exchange_nod(fct_LO)` +
  `exchange_nod(fct_plus/minus)`; `find_up_downwind_triangles` `coord_elem`/`e_nodes` (full-halo, global-id match);
  `nboundary_lay` owned+halo, NO exchange (partition-local — matches FESOM2). Gate: `tools/run_advhordump_multirank.sh`
  (oracle `fesom_advhor_dump.F90` npes>1 branch) + `fesom_advhordump_mr` (F3) → `advhor_diff.py` per rank.
  ⚠️ `edge_up_dn_grad` is captured PRE-`do_oce_adv_tra` (FESOM2 REUSES it as the FCT `AUX` scratch → bignumber=1e3).
  **No regression:** 1-rank geom(19)/advhor(40)/pressure(57)/step(65) + CORE2 lifecycle(195) + 13/13 ctest all
  `max|Δ|=0`/green. See LESSONS L31.
- **M2.12c-1 — multi-rank pre-SSH DYNAMICS chain — ✅ DONE (2026-06-22).** `max|Δ|=0` on pi dist_2 + dist_8 (25
  records density/pressure/bvfreq/Kv/ssh_rhs, `tools/run_stepdyn_gate_multirank.sh`). Lifted EVERY dynamics kernel
  up to `compute_ssh_rhs_ale` with the M2.12b optional-`partit` pattern (`mod_part_bounds.local_dims` for the
  owned/halo/full bounds + the exact FESOM2 exchanges, arithmetic UNCHANGED). New machinery: `mod_halo` rank-3
  `exchange_nod` node-block variant (UVnode/UVnode_rhs `(2,nl-1,nod)`); `local_dims` (adds `nEdgeL`/`nElemL` to
  `owned_bounds`). The exchanges added at the FESOM2 sites: `exchange_nod(UVnode)` (compute_vel_nodes),
  `smooth_nod` per-sweep `exchange_nod(bvfreq)`, `exchange_nod(UVnode_rhs)` (momentum_adv_scalar),
  `exchange_elem(U_c/V_c)` (visc_filt_bidiff, between its two Laplacian sweeps), `exchange_nod(ssh_rhs)`
  (compute_ssh_rhs_ale). visc interior-edge test lifted to the GLOBAL id (`myList_edge2D(ed)>edge2D_in`). Gate:
  FESOM3 `fesom_stepdump_mr` (dynamics chain only, no step_oce — the un-lifted CG is never reached) + the FESOM2
  `fesom_step_dump` shim extended to npes>1 (owned-prescribe + `exchange_elem(UV)`, runs the REAL `oce_timestep_ale`
  whose built-in per-rank gid-keyed `dump_shim_record_node` emit the probes; the post-ssh_rhs substeps 9/11/12/13/15/16
  are `--ignore-substep`'d since FESOM3 stops at ssh_rhs). Key finding: c-1 needs NO mesh-infra change — every kernel
  reads `elem2D_nodes`/`gradient_sca` only at OWNED elements (the rest read UV/UV_rhs/helem/elem_area at the M2.12b-
  exchanged full halo). See LESSONS **L32**.
- **M2.12c-2 — multi-rank SSH stiffness + free-surface CG — ✅ DONE (2026-06-22).** `max|Δ|=0` on `d_eta`, pi dist_2
  + dist_8 (30 records, `tools/run_stepdyn_gate_multirank.sh` with substep 9 un-ignored; 38 CG iters on 8 ranks). The
  FIRST multi-rank ITERATIVE solver. **KEY FINDING — the scoped "mesh-infra extension" (this bullet's prior text +
  L32) was WRONG and NOT needed.** A `dist_N`-file check (for each owned edge, is `edge_tri ≤ myDim_elem2D`? → 0 halo
  on every rank) proved the stiffness owned rows assemble FULLY LOCALLY. Invariants (pi dist_2/8; universal in FESOM2
  since its own `elem2D_nodes`/`gradient_sca` are owned-only, `oce_mesh.F90:497`/`:2466`): (i) every edge incident to
  an owned node is owned; (ii) both triangles of an owned edge are owned; (iii) owned-element nodes ≤ `nNodL`. So
  `init_stiff_mat_ale` lifts with ONLY bounds (`nod2D→nNodO`, `edge2D→nEdgeO`; `n_num(nNodL)`, `n_pos(12,nNodO)`;
  `ssh_stiff%dim=mesh%nod2D` global) — NO `elem2D_nodes`/`gradient_sca` halo, NO `enforce_cw` on halo, NO exchange in
  the assembly. The CG (`oce_ssh_solve`): arrays `nNodL`, owned loops `nNodO`, `exchange_nod(diag_values)` in the
  precond + `exchange_nod(pp/rr/x)` per-iter + NEW `mod_halo::allreduce_sum` (scalar+vec, `MPI_IN_PLACE`/`MPI_SUM`/
  `MPI_COMM_FESOM` = the oracle `solver.F90`) for the owned-partial-sum dot-products; `rtol`/convergence denominators
  GLOBAL (`mesh%nod2D`). Optional `partit` absent ⇒ 1-rank VERBATIM. The cross-rank `MPI_Allreduce` byte-matches
  (same OpenMPI+comm-size+op+8-byte-type → deterministic tree, L6); the L29 NOVECTOR precond divide carries over. Gate:
  `fesom_stepdump_mr` extended (`init_stiff_mat_ale`+`solve_ssh_ale`, dump `d_eta` at `DUMP_SUBSTEP_SSH_SOLVE`=9), the
  oracle already dumps `d_eta` at substep 9 (NO oracle change). **NO 1-rank regression** (step 65 + pressure 57 + advhor
  MR dist_2 + 13/13 ctest all `max|Δ|=0`). See LESSONS **L33**.
- **M2.12c-3 — multi-rank post-SSH ALE update + tracer SOLVE + the WHOLE step — ✅ DONE (2026-06-22).** `max|Δ|=0` on
  ALL **65 substep records** (density/pressure/bvfreq/Kv/ssh_rhs/d_eta/hbar/eta_n/hnode_new/w/T/S/hnode), pi dist_2 +
  dist_8 (`tools/run_step_gate_multirank.sh`, only SW_AB id=2 ignored = the 1-rank `run_step_gate.sh` set). **The whole
  multi-rank ocean step byte-matches FESOM2 — the architectural MVP.** Lifted EVERY post-SSH kernel via the M2.12b/c
  optional-`partit` pattern (bounds from `owned_bounds`/`local_dims`, arithmetic UNCHANGED): `update_vel` (owned
  elements + `exchange_elem_full(UV)` — rank-3, the eDim+eXDim superset of FESOM2's `exchange_elem`); `compute_hbar_ale`
  (ssh_rhs_old zero owned+halo, edge-div over owned edges, hbar update owned, then `exchange_nod(hbar)` ALWAYS — outside
  the linfs guard — so the owned-element `dhe` reads valid halo hbar); `update_eta_n` (owned+halo, no exchange);
  `vert_vel_ale` (Wvel zero owned+halo, edge-div owned, cumsum/divide owned, `exchange_nod(w)` + `exchange_nod(hnode_new)`)
  + `compute_CFLz`/`compute_Wvel_split` (owned+halo); `update_thickness_ale` (linfs no-op, helem unchanged + already
  full-halo-valid ⇒ FESOM2's `exchange_elem(helem)` is value-neutral and skipped, L33 "no machinery you don't need").
  Tracer SOLVE (`oce_ale_tracer.F90`): `solve_tracers_ale` (per-tracer `exchange_nod(values)` after the diffusion solve +
  the salinity clamp over owned+halo; `tr_xy` sized to the LOCAL `nElemF`); `diff_tracers_ale` ALE-reconstruct owned;
  `diff_part_hor_redi` owned edges (owned-node del_ttf complete, invariant i; `Ki`/`areasvol` read at the edge's halo
  node are prescribed/M2.12b-exchanged); `diff_ver_part_impl_ale` owned-node TDMA. Threaded the optional `partit` through
  `mod_step_oce::step_oce` (absent ⇒ 1-rank VERBATIM; the 1-rank callers `fesom_stepdump`/`fesom_lifecycle` unchanged).
  New whole-step MR driver `src/drivers/fesom_stepfull_mr.F90` (= `step_oce` through `partit`, local-sized state + tracer
  machinery + forcing inputs, gid-keyed per-rank dumps); the FESOM2 oracle is UNCHANGED (`fesom_step_dump` npes>1 already
  drives the whole `oce_timestep_ale` + its built-in dump_shim emits all substeps). The byte-gate on the rich analytic
  state (non-trivial T/S/UV/SSH at 2 AND 8 ranks) subsumes the suggested rest-at-rest / gravity-wave / `stale_halo_max_*`
  sanity probes — matching FESOM2's exact owned loop bounds + exchanges makes any stale halo unobservable in the owned
  dumps. **NO 1-rank regression** (step 65 + pressure 57 + advhor MR dist_2/8 + stepdyn MR dist_2/8 + 13/13 ctest all
  `max|Δ|=0`). See LESSONS **L34**.
- **MVP CLOSED 2026-06-22** — the whole-model multi-rank byte-match (M2.12) is committed and tagged `m2-mvp` (baseline
  re-confirmed after a clean rebuild: all 8 byte-gates `max|Δ|=0`/green, see the gate list in "RESUME HERE").
- **PRODUCTION VALIDATION (multi-rank FREE-RUNNING lifecycle) — ✅ DONE, BYTE-EXACT.** Built the multi-rank analog
  of the 1-rank CORE2 lifecycle (`fesom_lifecycle_mr` + `tools/run_lifecycle_gate_multirank.sh`, CORE2 dist_2, dt=1800,
  reduced-M2) — the FIRST test of *multi-rank × free-running × many-steps* (the m2-mvp gates are all single-step
  prescribe-and-stop, which MASK live-multi-step halo bugs). It found two latent bugs, **both now fixed**:
  1. ✅ **UV halo blow-up FIXED** — `update_vel` exchanged UV over `com_elem2D_full` (eDim+eXDim); FESOM2 uses
     `com_elem2D` (eDim). Wrong UV halo → viscosity → blow-up. Fix (`oce_ale.F90`): eDim `exchange_elem(UV)`. Model now
     STABLE + matches FESOM2's physical values to all printed digits.
  2. ✅ **`tr_xy` halo "drift" SOLVED — it was an OpenMPI `vader` KNEM single-copy bug, NOT a FESOM3 bug.** The
     exchange logic is correct; KNEM (the `vader` shared-memory single-copy path; CMA is blocked by
     `kernel.yama.ptrace_scope=3`) CORRUPTS messages >~131 KB on levante. A cross-rank buffer trace localized it to the
     MPI transport (`sbuf` packed correct → `rbuf` corrupt from byte 134656 → unpack faithful). The CORE2 `tr_xy`
     block exchange (478 KB) is the FIRST live F3 message above the threshold (every pi gate stays under it; FESOM2's
     `MPI_TYPE_INDEXED` path dodges KNEM) — which is why it only surfaced in the big-mesh free-running lifecycle.
     **Fix: `env.sh` exports `OMPI_MCA_btl_vader_single_copy_mechanism=none`** (forces byte-faithful copy-in/copy-out).
     Verified: CORE2 dist_2 lifecycle gate **195 records `max|Δ|=0`**; no 1-rank/pi regression. Full writeup +
     meta-lessons: [`HANDOFF-multirank-lifecycle-trxy.md`](HANDOFF-multirank-lifecycle-trxy.md) + LESSONS **L35**.
- **ACTIVE NEXT: M3 — sea ice (EVP).** The prescribed `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`stress_surf`
  become the real air-sea budget (sea ice EVP + advection + thermo → `oce_fluxes`). **Scoped 2026-06-22** (oracle map +
  call-sequence pass — see "M3 entry notes" below). M3 mirrors M2: byte-gate each ice kernel against FESOM2 v2.7.3,
  `max|Δ|=0`. Standard EVP path only (`whichEVP=0`; NO icepack/meltponds/cavity/oasis-yac). Decomposed:
  - **M3a — ice foundation + cold-start IC — ✅ DONE (2026-06-22).** `max|Δ|=0` on ALL 4 fields (`ice_a_ice`,
    `ice_m_ice`, `ice_m_snow`, `ice_massmatrix` [870146 CSR entries]), CORE2 1-rank (`tools/run_ice_gate_core2.sh`).
    Rewrote `src/types/mod_ice.F90` (full `t_ice`: `data(1:3)`=a_ice/m_ice/m_snow + work EVP stress/strain + thermo +
    rheology params, ALL initialisers VERBATIM from `MOD_ICE.F90` incl. un-suffixed singles); new
    `src/ice/mod_ice_setup.F90` (`ice_setup`/`ice_allocate`/`ice_mass_matrix_fill`/`ice_initial_state`, optional-`partit`
    pattern). **`ice_mass_matrix_fill` rides the M2.6 `ssh_stiff` CSR** — instead of FESOM2's `nn_num`/`nn_pos` it scans
    `ssh_stiff%colind_loc(rowptr_loc(row):..)` (the SAME ordered neighbour list `init_stiff_mat_ale` fills), and it
    byte-matched → confirms the FESOM3 CSR column order == FESOM2 `nn_pos`. Oracle: NEW `port2/fesom2/src/fesom_ice_dump.F90`
    (env `FESOM_ICE_DUMP`, dumps after `ice_setup`, stops) wired in `fesom_module.F90` after `ice_setup`; runner
    `tools/run_icedump_core2.sh` (the forced-lifecycle init path, use_ice=.true.). NEW driver `src/drivers/fesom_icedump.F90`.
    Restart binary I/O for ice DEFERRED (M3 is cold-start; gates use mod_dump). **No regression** (13/13 ctest +
    step-65 `max|Δ|=0`). ⚠️ CORE2 `namelist.ice` (the oracle's effective config) DIFFERS from the `MOD_ICE.F90` type
    defaults: `ice_gamma_fct=0.5` (not 0.25), `ice_diff=0.0` (NO artificial ice diffusion, not 10.0), `albw=0.1` (not
    0.066); and the namelist parses literals as DOUBLES (so `cd_oce_ice=0.0055` etc. = `_WP` doubles, NOT single→WP) —
    bake the namelist doubles into the M3b driver config (these don't affect the SST-sign IC / geometry mass matrix).
  - **M3b — ocean2ice + EVP dynamics — ✅ DONE (2026-06-22).** `max|Δ|=0` on ALL 7 fields (`ice_srfoce_u`/`srfoce_v`
    [ocean2ice], `ice_uice`/`vice`, `ice_sigma11`/`12`/`22`), CORE2 1-rank, for **BOTH `whichEVP=0` (standard EVP) AND
    `whichEVP=1` (modified EVP / mEVP)** (`tools/run_evp_gate_core2.sh`; mEVP added per user request). New
    `src/ice/mod_ice_dyn.F90`: `ocean2ice` (srfoce_u/v/temp/salt/ssh from the ocean surface; node-avg of `UV(:,1,:)` +
    `exchange_nod`) + std-EVP `EVPdynamics` (`stress_tensor` strain→stress + `stress2rhs` ∇·σ + 120-subcycle velocity
    update [ocean drag + Coriolis + ∇·σ + wind-on-ice stress] + edge BC + per-subcycle `exchange_nod(u_ice/v_ice)`) +
    `EVPdynamics_m` (mEVP, Bouillon-2013 α/β pseudo-time iteration on `uice_aux`/`vice_aux`, inlined
    `stress_tensor_m`+`stress2rhs_m`+`ssh2rhs`, the `bc_index_nod2D` node BC) + `EVPdynamics_solve` dispatcher
    (whichEVP SELECT CASE; aEVP=2 deferred). Optional-`partit` pattern throughout (1-rank verbatim; exchanges guarded by
    `is_multirank`). Gate driver `src/drivers/fesom_evpdump.F90` (`FESOM3_WHICHEVP` env) prescribes the analytic ocean
    forcing (surface UV / hbar / stress_atmice — identical formulas to the oracle, from byte-identical rotated coords),
    runs ocean2ice + the EVP solve, dumps. Oracle: `port2/fesom2/src/fesom_ice_dump.F90::evp_dump_write` (env
    `FESOM_EVP_DUMP`, dispatches on namelist `whichEVP`, runs the REAL `ocean2ice`+`EVPdynamics`/`EVPdynamics_m`, dumps,
    stops; wired in `fesom_module.F90` after `ice_dump_write`) + `tools/run_evpdump_core2.sh [run] [out] [whichEVP]`.
    **M3a-deferred items built:** `bc_index_nod2D` (the mEVP node boundary mask) + `uice_aux`/`vice_aux` allocation in
    `mod_ice_setup` (`ice_allocate`/`build_bc_index_nod2D`). ⚠️ **`bc_index_nod2D` lives in `t_ice`, NOT `t_mesh`** —
    adding any component to `t_mesh` (which has type-bound `write(unformatted)`) trips an ifort `#7976` generic-
    serialization cascade across every `write_bin_array` call; `t_ice` has no type-bound I/O so it's safe. Value is
    byte-identical (0/1 mask). ⚠️ `cd_oce_ice`/`delta_min` set to the namelist DOUBLES (`0.0055_WP`/`1.0e-11_WP`) in the
    driver, NOT the `t_ice` single→WP defaults (the M3a precision note). mEVP uses `rdt=ice_dt` (NOT `ice_dt/steps`),
    `alpha_evp`/`beta_evp`=250 (exact). **No regression** (step 65 + M3a ice 4 + 13/13 ctest all `max|Δ|=0`).
  - **M3c — ice FCT advection — ✅ DONE (2026-06-22; first try).** `max|Δ|=0` on ALL 6 fields (`ice_rhs_a`/`rhs_m`/
    `rhs_ms` from `ice_TG_rhs` + post-advection `ice_a_ice`/`m_ice`/`m_snow` from `ice_fct_solve`), CORE2 1-rank, for
    **BOTH `whichEVP=0` AND `whichEVP=1`** feeding the advection (`tools/run_icefct_gate_core2.sh`). New
    `src/ice/mod_ice_fct.F90` (`ice_TG_rhs` Taylor–Galerkin rhs + `ice_fct_solve` = `ice_solve_high_order` 3-iter Jacobi
    → `ice_solve_low_order` → `ice_fem_fct` Zalesak limiter ×3), optional-`partit` pattern. **KEY: the neighbour-list
    basis is the ssh_stiff CSR (`rowptr_loc`/`colind_loc`), NOT FESOM3's own `nn_pos`** — M3a proved `fct_massmatrix`
    byte-matches FESOM2's `mass_matrix` position-by-position, so `colind_loc(clo:clo2)` IS FESOM2's `nn_pos(1:nn_num,row)`
    as an ordered list (incl. self), and `sum(fct_massmatrix(clo:clo2)*field(colind_loc(clo:clo2)))` reproduces
    `sum(mass_matrix*field(nn_pos))` byte-for-byte (same operands + reduction order). Driver
    `src/drivers/fesom_icefctdump.F90` (extends `fesom_evpdump`: ocean2ice + EVP → `uice/vice`, then `ice_TG_rhs` +
    `ice_fct_solve`); oracle `fesom_ice_dump.F90::ice_fct_dump_write` (env `FESOM_FCT_DUMP`); runners
    `run_icefctdump_core2.sh`/`run_icefct_gate_core2.sh`. ⚠️ `ice_diff=0.0`/`ice_gamma_fct=0.5` set by the driver to the
    CORE2 namelist values (NOT the `t_ice` 10.0/0.25 defaults; both exactly representable, no ULP subtlety). Oracle
    `ENABLE_OPENMP=OFF` ⇒ the `ice_fem_fct` `!$OMP ORDERED`/`ATOMIC` flux scatters run serially in element order = the
    FESOM3 serial path (accumulation byte-matches). `ice_fct_solve` is called as an external (not in
    `ice_fct_interfaces`); `ice_TG_rhs` IS. **No regression** (step 65 + M3a ice 4 + M3b evp 7×2 + M3c fct 6×2 + 13/13
    ctest all `max|Δ|=0`; the byte-exact `rhs` independently re-confirms M3b — the `uice/vice` feeding `ice_TG_rhs` were
    byte-identical). See LESSONS **L37**.
  - **M3d — ice thermodynamics — ✅ DONE (2026-06-23; first try, L38).** `max|Δ|=0` on ALL 6 fields (post-thermo
    `ice_a_ice`/`m_ice`/`m_snow` + `ice_flx_h`/`flx_fw` + `ice_t_skin`), CORE2 1-rank, for **BOTH `whichEVP=0` AND
    `whichEVP=1`** feeding the advection upstream (`tools/run_icethermo_gate_core2.sh`). New `src/ice/mod_ice_thermo.F90`
    (`cut_off` hmin/Armin clamp + `thermodynamics` per-node driver + `therm_ice` 0-layer Semtner/Hibler-1984 ice-class
    growth + `budget` 5-iter Newton-Raphson surface temp + `obudget` open-ocean growth+evaporation + `flooding` snow→ice
    + `TFrez` Millero-1978 + solar-zenith/albw helpers [ported for linkage, dead at `open_water_albedo=0`]),
    optional-`partit` (`exchange_nod(ustar)` guarded). **Atmospheric forcing + thermo flux diagnostics bundled into
    `t_atmflux`** (FESOM3 explicit-dataflow vs FESOM2's `g_forcing_arrays`/`g_forcing_param`/`o_arrays` globals); the gated
    outputs live in `t_ice`. Added `h_cutoff`/`hpdf` to `t_ice_thermo` (read unconditionally by `therm_ice`, dead at
    `new_iclasses=.false.`). Driver `src/drivers/fesom_icethermodump.F90` (extends `fesom_icefctdump`: + prescribed
    analytic atmosphere + `cut_off` + `thermodynamics`); oracle `fesom_ice_dump.F90::ice_thermo_dump_write` (env
    `FESOM_THERMO_DUMP`, wired in `fesom_module.F90` after `ice_fct_dump_write`). **KEY (L38): the `&ice_therm` params are
    read by the oracle as namelist DOUBLES (con=2.1656/consn/hmin/Armin=0.01/emiss=0.97/albedos/albw=0.1 [not 0.066]/
    Sice/h0/c_melt) — NOT the `t_ice_thermo` single→WP defaults; the driver re-sets them as `_WP` doubles. Likewise
    `Ch_atm_ice=Ce_atm_ice=0.00175_WP` (namelist.forcing doubles), `ref_sss_local=.true.` (⇒ rsss=S_oc),
    `use_virt_salt=.true.` (linfs), `l_snow=.true.`. cc=rhowat·4190 / cl=rhoice·3.34e5 are exactly representable
    (ice_init recompute = the single→WP default). All thermo arithmetic is SCALAR per-node ⇒ no L29 vectorised-divide
    trap.** No regression (step 65 + M3c fct 6×2 + 13/13 ctest all `max|Δ|=0`; M3c re-run since the shared oracle
    `fesom_ice_dump.F90` was edited).
  - **M3e — `oce_fluxes` coupling-out (THE PAYOFF) — ✅ DONE (2026-06-23; first try, L39).** `max|Δ|=0` on ALL 5 fields
    (`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` [nod2D] + `stress_surf` [2,elem2D]), CORE2 1-rank, for **BOTH
    `whichEVP=0` AND `whichEVP=1`** (`tools/run_iceflux_gate_core2.sh`). **These 5 fields ARE the proven M2.11c-2
    `fesom_flux_dump` set — the `fesom_lifecycle` PRESCRIBED M3-gap is now produced NATIVELY.** New
    `src/ice/mod_ice_oce_coupling.F90`: `oce_fluxes_mom` (ice-ocean drag blended with atm-ocean stress by ice
    concentration on nodes → averaged to `stress_surf` on elements) + `oce_fluxes` (the reduced air-sea budget:
    `heat_flux=-flx_h`, `water_flux=-flx_fw` then globally balanced, `virtual_salt`/`relax_salt` each globally balanced
    to zero net via `integrate_nod`/`ocean_area`) + a private `integrate_nod_2D` (the EXACT FESOM2 `gen_support.F90:318`
    sequential reduction, NOT `sum()`). **KEY (L39): `mesh%ocean_area` — the balancing divisor — was an un-gated landmine:**
    computed since M2.11a as `sum(area(1,1:nNodO))` but never consumed, so never gated; `oce_fluxes` is its first consumer,
    and `sum()` ≠ FESOM2's sequential `areasvol` loop is a 1-ULP trap that would drift every balanced field. Fixed
    `mod_mesh_areas.F90` to the faithful FESOM2 `oce_mesh.F90:2385` loop (geometry-neutral — `ocean_area`/`areawithcav`
    are the only changed scalars, no gate consumes them but M3e). `dens_flux` (MOC diagnostic, needs `sw_alpha`/`sw_beta`
    EOS) DEFERRED (out of M3 scope, not a surface BC). Driver `src/drivers/fesom_icefluxdump.F90` (extends
    `fesom_icethermodump`: + prescribed analytic `stress_atmoce_x/y` + `Ssurf`); oracle
    `fesom_ice_dump.F90::ice_flux_dump_write` (env `FESOM_OCEFLUX_DUMP`, after `ice_thermo_dump_write`; + `wr_r2` for the
    2-D `stress_surf`); runners `run_icefluxdump_core2.sh`/`run_iceflux_gate_core2.sh`. `t_atmflux` extended with the
    oce_fluxes I/O (`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`heat_flux_in`/`stress_node_surf`/
    `stress_atmoce_x/y`/`Ssurf`/`surf_relax_S`=1.929e-06). **No regression** (M3e 5×2 + M3d 6×2 + step 65 1-rank + step 65
    multirank dist_2 + pressure 57 + 13/13 ctest all `max|Δ|=0`).
  - **M3f — whole ice step + forced lifecycle (NATIVE fluxes) + multi-rank.** IN PROGRESS — M3f-1 + M3f-2 ✅ DONE
    (2026-06-23, first try, L40); M3f-3 (native CORE2 forcing read) + M3f-4 (multi-rank) remain.
    - **M3f-1 — `ice_timestep` assembly — ✅ DONE.** New `src/ice/mod_ice_step.F90::ice_timestep` = the faithful
      FESOM2 `ice_setup_step.F90:96` chain `EVPdynamics_solve → ice_TG_rhs → ice_fct_solve → cut_off → thermodynamics`
      (CMIP6 `dyngr*` / `h_ice` diagnostics + cavity cleans omitted — unconsumed in the reduced config). Refactored
      `fesom_icefluxdump` to call it; M3e gate re-run `max|Δ|=0` (byte-neutral, both whichEVP).
    - **M3f-2 — native-flux coupled multi-step forced lifecycle — ✅ DONE (first try).** New driver
      `src/drivers/fesom_lifecycle_native.F90` (= the proven M2.11c-2 ocean init + `step_oce`, but the prescribed
      `FESOM3_FLUX_FILE` is REPLACED by the native chain `ocean2ice → ice_timestep → oce_fluxes_mom → oce_fluxes`,
      runloop order ocean2ice→[atm]→ice_timestep→oce_fluxes_mom→oce_fluxes→step_oce). The per-step POST-bulk
      atmospheric forcing (16 nod2D arrays: shortwave/longwave/Tair/shum/prec_rain/prec_snow/runoff/u_wind/v_wind/
      Ch-Ce_atm_oce/stress_atmoce_x-y/stress_atmice_x-y/Ssurf) is PRESCRIBED from a NEW oracle dump
      (`port2/fesom2/src/fesom_atmflux_dump.F90`, env `FESOM_ATMFLUX_DUMP`, wired in `fesom_module.F90` next to
      `flux_dump_record`) so the native sea-ice + air-sea coupling is isolated from the CORE2 forcing read (M3f-3).
      **Result: `max|Δ|=0` on the 13 NODE substeps × 5 probes — 195 records (3 steps) AND 325 records (5 steps) — for
      BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP).** The driver's built-in per-step native-vs-oracle flux
      self-check (`FESOM3_FLUX_FILE`) reads `max|Δ(hf,wf,vs,rs,ss)|=0` at EVERY step. This is the FIRST test of the
      MULTI-STEP ice evolution (sigma elastic memory + t_skin carry-over + `values_old` — the M3a-e gates were all
      single-step from cold start) AND the FIRST native ice↔ocean coupled run. Gate
      `tools/run_lifecycle_native_gate_core2.sh [run] [nsteps] [whichEVP]` (oracle `run_lifecycle_forced_core2.sh`
      extended with `FESOM_ATMFLUX_DUMP` + a `whichEVP` namelist patch). **No regression** (step 65 1-rank + ctest 13/13
      + M2.11c-2 prescribed-flux forced lifecycle 195 + M3e flux 5×2 all `max|Δ|=0`). See LESSONS **L40**.
    - **M3f-3 — native CORE2 forcing read + bulk. M3f-3a + M3f-3c ✅ DONE 2026-06-23 (first try, L41); M3f-3b remains.**
      **14/16 atmospheric arrays are now NATIVE** (`max|Δ|=0`); only runoff + Ssurf (the monthly climatology) stay
      prescribed.
      - **M3f-3a — native CORE2 NCAR forcing READ — ✅ DONE.** The 8 NCAR fields (shortwave/longwave/Tair/shum/
        prec_rain/prec_snow/u_wind/v_wind from u_10/v_10/q_10/t_10/ncar_rad/ncar_precip) read + bilinear-interp +
        g2r-rotated via the M2.10a `mod_forcing_read` (already byte-proven on these EXACT files, on pi), now on the
        CORE2 mesh + the per-step `rdate` advance (clock `0 1 1948`, dt=1800, `sbc_do:1527` half-step;
        `rdate(n)=julday(1948,1,1,noleap)+(2n-1)·900/86400`; cold-start getcoeffld + per-step timeinterp, NO
        re-trigger within forcing day 1). Driver `src/drivers/fesom_forcing_core2.F90` self-checks the 8 fields vs the
        M3f-2 `atmflux_f2` dump: `max|Δ|=0`, 3 steps.
      - **M3f-3c — native NCAR bulk + wind stress — ✅ DONE.** `Ch_atm_oce_arr`/`Ce_atm_oce_arr` (`forcing_bulk_ncar`,
        M2.10b, fed the LIVE `ice%srfoce_temp/u/v` from ocean2ice — same inputs the oracle's
        `ncar_ocean_fluxes_mode` reads), `stress_atmoce_x/y` (`forcing_wind_stress`, M2.10b), and the NEW
        wind-on-ice `stress_atmice_x/y` (`forcing_ice_stress` in `mod_forcing_bulk`:
        `Cd_atm_ice·rhoair·|uw-uice|·(uw-uice)`, Cd_atm_ice=0.0012 const, fed the previous-step `ice%uice/vice`).
        Validated IN the lifecycle via a per-step self-check in `fesom_lifecycle_native` (env `FESOM3_FORCING_DIR`):
        the native NCAR read + bulk + 2 stresses are recomputed each step and compared to the prescribed atm —
        `max|Δ|=0` on all 14 fields, BOTH whichEVP, while the step itself still USES the prescribed values (so the
        195-record gate cannot regress). No regression: M2.10 forcing gate (pi, `mod_forcing_bulk` changed) PASS,
        ctest 13/13, the 195-record lifecycle gate MATCH.
      - **M3f-3b — runoff + SSS monthly climatology read (the last 2/16) — ✅ DONE (2026-06-23; first try, L42).**
        New `src/forcing/mod_forcing_other.F90` (`read_other_NetCDF` + `interp_2d_field`, faithful from FESOM2
        `gen_modules_read_NetCDF.F90:6` / `gen_interpolation.F90:145`) — a DISTINCT reader from `mod_forcing_read`
        (one 2D slice, raw-grid dummy fill, direct per-node bilinear, NO time-interp). runoff: `read_other_NetCDF
        ('Foxx_o_roff', rec 1, check_dummy=.false.)` → missing→0, then `/1000` (kg/s/m²→m/s); Ssurf: `('SALT', rec
        month=1, check_dummy=.true.)` → 30-neighbour fill. Both vertices (`do_onvert=.true.`; centroid path
        `error stop`-guarded). +2 real64 netCDF wrappers (`nc_get_slice_dp`/`nc_get_att_dp` in `mod_io_netcdf`). Wired
        into `fesom_lifecycle_native` (read ONCE at setup — both constant for a Jan run); self-check `max|d clim
        (runoff,Ssurf)|=0`. **THEN dropped the prescribed atm entirely → M3f-3 COMPLETE: the FULLY-NATIVE lifecycle**
        (`FESOM3_ATMFLUX_FILE` absent ⇒ compute all 16 in-driver via `apply_native_forcing`) is `max|Δ|=0` vs FESOM2,
        195 records (3-step) + 325 records (5-step), BOTH whichEVP, `tools/run_lifecycle_fullynative_gate_core2.sh`.
        ⚠️ TRAP (L42): the standalone 1-rank drivers do NOT populate `partit%myDim_nod2D` (it stays 0) — passing
        `partit` to an optional-`partit` routine whose loop bound is `partit%myDim_*` gives `num=0` → SILENT all-zero
        output (no crash). FIX: omit `partit` at 1-rank (count from `mesh%nod2D`), exactly like `ocean2ice`/
        `ice_timestep`. M3f-4 will pass partit (dims set by `read_dist_partition`).
    - **M3f-4 — multi-rank — ✅ DONE (2026-06-23; first try, L43) → M3 COMPLETE, tag `m3`.** `max|Δ|=0` on the per-rank
      gid-keyed NODE substeps — **195 records (3-step) AND 325 records (5-step), BOTH whichEVP=0 (std EVP) AND
      whichEVP=1 (mEVP), on CORE2 dist_2 AND dist_8** (`tools/run_lifecycle_fullynative_gate_multirank.sh [np] [nsteps]
      [whichEVP]`). The MULTI-RANK fully-native lifecycle (`src/drivers/fesom_lifecycle_native_mr.F90` = the
      `fesom_lifecycle_mr` MR mesh/state/dump scaffold ⊕ the `fesom_lifecycle_native` ice/atm/forcing setup + runloop,
      at LOCAL sizes, `, partit` on every kernel call) drives the whole coupled sea-ice + air-sea budget + atmosphere
      (16/16 native arrays) byte-exactly at multi-rank with NO prescribed input. **KEY (L43): the M3a–e kernels were
      ALREADY MR-ready** (each transcribed with the M2.12 optional-`partit` pattern; `ice_setup` sizes via `local_dims`,
      `ocean2ice`/EVP/FCT exchange their fields, `cut_off`+`thermodynamics` loop owned+halo so outputs stay halo-valid,
      `oce_fluxes` does the cross-rank `integrate_nod_2D`). The EVP needs **NO sigma element-halo** — `stress2rhs`
      scatters owned-elements→owned-nodes and the owned-node-completeness invariant (M2.12a) makes the owned rhs
      complete (FESOM2 `ice_EVP.F90` has only `exchange_nod(U_ice,V_ice)`). Only TWO things were genuinely missing:
      (1) `mesh%ocean_area`/`ocean_areawithcav` were a LOCAL-only owned sum → added `allreduce_sum` at npes>1
      (`mod_mesh_areas.F90`, FESOM2 `oce_mesh.F90:2389` `MPI_AllREDUCE`) — the M3e flux-balance divisor; (2) the 3
      `mod_forcing_bulk` routines looped `mesh%nod2D` (GLOBAL) → optional `partit` + `owned_bounds` → loop `nNodL`
      (FESOM2 computes the bulk over `myDim+eDim`, no exchange). The NCAR read needed only `frc%nnod=nNodL`;
      `read_other_NetCDF` gets `partit` at MR (1-rank omits it, the L42 trap). Oracle: `run_lifecycle_forced_core2.sh`
      now takes `[np]` (FESOM2 auto-reads dist_<np>; the 1-rank-only flux/atmflux dumps skip at npes/=1, the gid-keyed
      `dump_shim` is the gate). **No regression:** ctest 13/13 + 1-rank fully-native 195 + iceflux 5×2 + forcing pi +
      step-65 1-rank + step-65 MR dist_2 all `max|Δ|=0`. ⚠️ multi-rank levante needs the `env.sh` KNEM flag (L35).

**↳ RESUME HERE (next session): M9 = the production I/O harness** (netCDF mean/snapshot output +
restart write/read + mesh.nc) — M8 (production long runs) is ✅ COMPLETE (tag `m8`). M8 delivered the
clock + forcing record/day/month/year rollover + run-length-in-years driver, all `max|Δ|=0` at `dist_864`;
the 2-year JRA55 1958-1959 headline is **byte-exact `max|Δ|=0`** — the day-107 divergence was root-caused to a
runtime FTZ (flush-to-zero) mismatch (FESOM3 process kept a denormal `m_snow` that FESOM2 flushed to 0.0) and
FIXED via `ieee_set_underflow_mode(gradual=.false.)` in the driver (**L51**); confirmed byte-exact over a full
model year (1.14M records, 3.4× past day-107) + the FTZ gate + the year-rollover gate. **Next = M9** (see the
"M9 scope" section of `docs/plans/2026-06-25-m8-long-simulations.md`):
port `io_meandata.F90` (the FIRST netCDF-WRITE gate — `mod_io_netcdf.F90` only READS) + `io_restart.F90`
incl. **`tke` serialization** (the first stateful mixing field) + `clock_finish` (the `.clock` write,
M8a-deferred) + `io_mesh_info`; gate = restart-reproducibility (`max|Δ|=0` straight-through vs split-restart)
+ field-by-field vs the F2 output. **Known M8 follow-up (only matters past 1960):** the year-rollover
crossing-test dereferences a stale `t_indx_p1` — safe for the 2920→2920 (1958→1959) headline, but must reset
before a leap-OUT crossing (2928→2920, e.g. 1960→1961).

**DIRECTION (user 2026-06-24):** target = `which_ALE='zstar'` (canonical `config/namelist.config`; the
`work_core` copy was customized to zlevel — my brief zlevel/paper-parity detour was wrong, corrected).
**Drop "paper-parity" framing** (goal = faithful 1:1 port + better architecture; paper = example only).
**No aEVP.** **TKE = M7, own plan, from the FORTRAN oracle with C-port experience (no CVMix, no
diagnostics-array bloat).** **Production runs (months/years, ~byte-identical vs Fortran) = M8, own plan.**

### M6 progress (zstar ALE) — Plan: `docs/plans/2026-06-24-m6-zstar.md`.
- **Decomposition:** M6a-1 (Shchepetkin PGF) → M6a-2 (zstar thickness machinery, unforced) → M6a-3
  (forced/native + real freshwater flux, both whichEVP) → M6a-4 (multi-rank). Then M7 TKE / M8 production.
- **Config (canonical `config/namelist.config`):** `which_ALE='zstar'`, `which_pgf='shchepetkin'`
  (DEFAULT, `oce_modules.F90:190`), `use_partial_cell/use_cavity/use_floatice=.false.` ⇒ cavity/floatice/
  partial-cell branches DEAD. `which_ALE/='linfs' ⇒ use_virt_salt=.false.` (`oce_setup_step.F90:120`) ⇒
  real freshwater flux replaces virtual salt (M6a-3). `update_stiff_mat_ale` runs per-step (`:3921`).
  ⚠️ zstar is SIMPLER than zlevel (no local-zstar fallback / min_hnode / lzstar_lev / CFL_z) — it
  stretches proportionally over all full-prism layers (`hnode*(1+η/H)`), all scalar divides (no L29).
- **M6a-1 ✅ DONE (byte-exact, first try):** NEW `pressure_force_4_zxxxx_shchepetkin` in `src/oce/oce_pgf.F90`
  (Shchepetkin density-Jacobian PGF, FESOM2 `oce_ale_pressure_bv.F90:2104-2339`; self-contained, NO
  hpressure; shared by zlevel/zstar). CORE2 pressure gate `pgf_x_shchep`/`pgf_y_shchep` `max|Δ|=0`,
  non-vacuous ~3.1e-5/2.8e-5, 97 fields (env `FESOM3_PGF_SHCHEP`/`FESOM_PGF_SHCHEP`). At rest the
  Shchepetkin ≈ linfs to ~1e-17 (flat layers ⇒ `dz_dx≈0`), so the η≠0 stress test was M6a-2 (passed).
  Committed `801a68f`.
- **⚠️ PRE-EXISTING m5 BASELINE BUG fixed in M6a-1 (L47):** the pressure shim left `use_sw_pene=.true.`
  into the M2.9a tracer solve → the M5c `sw_3d` term read UNALLOCATED `dyn%work%sw_3d` (Release tolerated
  it until a heap re-layout → SIGSEGV; debug `-check all`: `oce_ale_tracer.F90:665`). FIX: mirror the
  oracle — `use_sw_pene=.false.` before the tracer-solve loop. The M5a-3 META-LESSON realized.
- **M6a-2 ✅ DONE (byte-exact, FIRST try, 2026-06-24):** the zstar free surface byte-matches FESOM2
  `max|Δ|=0` — **195 (3-step) + 325 (5-step)** in the UNFORCED CORE2 lifecycle, SSH/hnode EVOLVING
  (max|eta_n| 0.35→0.68). Three new pieces (all optional-`partit`, owned loops, M6a-4-ready):
  (1) **`update_stiff_mat_ale`** (NEW, `oce_ssh_rhs.F90`; FESOM2 `:1892-2001`) — per-step SSH-stiffness
  2nd-term update, ADD `-dhe(elem)*(gradient_sca·edge_cross_dxdy)*factor` over owned edges; byte-exact
  because the CSR (rowptr_loc/colind_loc/values) is M2.6-proven + same edge order; `dhe` is LAGGED (prev
  step's `compute_hbar_ale`; step-1 dhe=0); additions ACCUMULATE to track the moving surface.
  (2) **`vert_vel_ale` zstar branch** (`oce_ale.F90`; FESOM2 `:2755-2821`) — `dd=(hbar-hbar_old)/H`
  (H over `nlevels_nod2D_min` full prisms), `Wvel(nz)-=(zbar_3d_n(nz)-zbar_3d_n(nzmax))*dd/dt`,
  `hnode_new(nz)=hnode(nz)+(zbar_3d_n(nz)-zbar_3d_n(nz+1))*dd` (water_flux term deferred to M6a-3).
  (3) **`update_thickness_ale` zstar commit** (`oce_ale.F90`; FESOM2 `:1378-1436`) — `hnode_new→hnode`,
  recompute `zbar_3d_n`/`Z_3d_n` BOTTOM-UP from the fixed anchor, `helem`=node-mean, `exchange_elem`.
  Wired into `mod_step_oce`: PGF dispatch (`which_ALE=='linfs'?linfs:Shchepetkin`, pass `density_m_rho0`)
  + `if(which_ALE/='linfs') call update_stiff_mat_ale` before `compute_ssh_rhs_ale`. Driver:
  `FESOM3_WHICH_ALE=zstar` in `fesom_lifecycle.F90`; oracle `run_lifecycle_core2.sh` gained `WHICH_ALE`
  env (NO oracle source change — zstar is native FESOM2). Gate `tools/run_lifecycle_zstar_gate_core2.sh`.
  No regression: linfs lifecycle 195 + step-65 + ctest 13/13 all `max|Δ|=0`.
- **M6a-3 ✅ DONE (2026-06-24): the FULL PRODUCTION zstar config byte-exact.** Forced/fully-native zstar
  (`use_virt_salt=.false.`, `is_nonlinfs=1`) `max|Δ|=0` vs FESOM2 — **195 (3-step) + 325 (5-step), BOTH
  whichEVP, reduced (PP) AND production (KPP+GM+Redi+sw_pene)** + every per-step native-flux self-check
  `=0`. `tools/run_lifecycle_zstar_native_gate_core2.sh` (+ `WHICH_ALE` env in `run_lifecycle_forced_core2.sh`
  / `run_lifecycle_fullynative_gate_core2.sh`). Three real-freshwater pieces (the M3-anticipated branches
  in bc_surface/oce_fluxes/obudget were already there — these were the GAPS the gate found, fixed in order):
  - **TWO water_flux terms, not one:** `compute_ssh_rhs_ale` `ssh_rhs-=alpha*water_flux*areasvol(nzmin)`
    (FESOM2 `oce_ale.F90:2122-2134`, open-ocean branch; the substep-8 divergence) AND `compute_hbar_ale`
    `ssh_rhs_old-=water_flux*areasvol` (`:2262`) AND `vert_vel_ale` `Wvel(nzmin)-=water_flux` (`:2809`).
    All threaded as an OPTIONAL `water_flux` arg (absent ⇒ linfs/pressure-shim skip; guarded
    `which_ALE/='linfs'`), so linfs is byte-neutral.
  - **the `atm%real_salt_flux` wiring bug:** `obudget` computes the ice-brine salt flux `rsf=fwice*Sice`
    into `atm%real_salt_flux` (the `if(.not. use_virt_salt)` branch — M3 already wrote it, matches the
    oracle verbatim), but `fesom_lifecycle_native` passed a LOCAL zeroed `real_salt_flux` to `step_oce`
    → S diverged ~2e-3 at substep 15. FIX: pass `atm%real_salt_flux`. (linfs ignores it via is_nonlinfs=0.)
  - **driver flags:** `FESOM3_WHICH_ALE=zstar` ⇒ `which_ALE='zstar'` + `atm%use_virt_salt=.false.`
    (drives both `oce_fluxes`'s levitating-ice branch AND `obudget`'s `rsf`) + `is_nonlinfs=1.0`
    (bc_surface uses `real_salt_flux` + the `sval*water_flux` advective-heat term). NO oracle SOURCE change.
  No regression: linfs native (reduced + production) + linfs/zstar UNFORCED + step-65 + ctest 13/13 `max|Δ|=0`.
- **M6a-4 ✅ DONE (first try, 2026-06-24) → M6a COMPLETE.** The production zstar config byte-matches
  FESOM2 `max|Δ|=0` at MULTI-RANK — CORE2 dist_2 + dist_8, reduced (PP) AND production (KPP+GM+Redi+
  sw_pene), BOTH whichEVP. PURE WIRING (M4f/M5d): the zstar kernels are optional-`partit` (update_stiff
  owned-edge loop; vert_vel_ale/update_thickness_ale owned loops; the exchanges — `exchange_elem(helem)`,
  `exchange_nod(Wvel/hnode_new/ssh_rhs/ssh_rhs_old)` — already in place; `dhe` read at owned el(i) only,
  no halo extension, the M2.12c-2 owned-edge invariant). Only NEW work: the `fesom_lifecycle_native_mr`
  zstar block (FESOM3_WHICH_ALE + `atm%use_virt_salt`/`is_nonlinfs` + `atm%real_salt_flux`, LOCAL-sized,
  mirror of M6a-3) + `WHICH_ALE` thread in `run_lifecycle_fullynative_gate_multirank.sh` +
  `run_lifecycle_zstar_native_gate_multirank.sh`. No oracle SOURCE change. No regression: linfs MR
  native dist_2 + ctest 13/13. **Tag `m6`.**

### M5 progress (KPP vertical mixing + sw_pene) — scoped + started 2026-06-24. Plan: `docs/plans/2026-06-24-m5-kpp.md`.
- **Decomposition:** M5a (producers, isolated gate; sub-stepped a1–a4 mirroring the C-port K1–K8) →
  M5b (wire into step + UNFORCED lifecycle) → M5c (forced/native + sw_pene + ghats) → M5d (multi-rank). Tag `m5`.
- **KEY STRUCTURAL FINDINGS (from reading the actual `oce_ale_mixing_kpp.F90` + the dispatch + the validated C port):**
  - `oce_ale.F90:3713` dispatch: `oce_mixing_KPP(Av, Kv_double, ...)` — **`Av` (viscAE) stays element-based** (KPP
    averages node→elem internally, `minmix=3e-3` floor) ⇒ `impl_vert_visc` UNCHANGED; then `Kv(:,n)=Kv_double(:,n,1)`
    (T-channel) + `mo_convect` ⇒ the tracer vertical-diffusion TDMA is **UNCHANGED** (same single `Kv` as PP). The
    S-channel `Kv_double(:,:,2)` is computed-but-dead in CORE2 (only a gated output). So KPP integration = the KPP
    module itself + `dbsfc` (in pressure_bv) + sw_pene wiring + TWO new tracer terms (ghats nonlocal + sw_3d).
  - `dbsfc` is filled in `oce_ale_pressure_bv.F90:332-339` (a NEW output of the already-ported pressure_bv) — in
    FESOM3 a `dyn%work%dbsfc` field; `ghats`/`blmc` are READ by the tracer solve ⇒ live in `dyn%work` (not KPP-locals).
  - `cal_shortwave_rad` (sw_3d producer) is ALREADY ported + byte-proven (M2.10c); M5 wires it into the lifecycle +
    adds the TDMA consumer (`oce_ale_tracer.F90:991-996` sw_3d term, `:900-945` ghats term). `ddmix` DEFERRED
    (`double_diffusion=.false.`). Config DOUBLES: `Ricr=0.3`, `concv=1.6`, `visc_sh_limit`/`diff_sh_limit=5e-3`,
    `A_ver=1e-4`, `K_ver=1e-5`, `Kv0_const=.true.`.
- **HARNESS (no oracle rebuild needed):** the oracle `build/lib/libfesom.so` (Jun-23 M4 build) ALREADY carries the
  `FESOM_KPP_DUMP_DIR` instrumentation (`kpp_dump_*` symbols; env-gated off by default → harmless). `bin/fesom.x`
  (Apr-23, stale) is a thin stub that loads that `.so` at runtime — the gates use `build/bin/fesom.x`. **NEW
  `tools/run_kppinit_oracle.sh`** runs CORE2 1-rank with `mix_scheme=KPP` (proven forced-runner setup, no PP
  downgrade) → emits `kpp_init_rank0.txt` + `kpp_wscale_rank0.txt` **AND all the per-node K3/K5/K6/K7/K8 dumps**
  (`kpp_dump_s1_{ri_*,prestep,dVsq,dbsfc,blmc_*,diffK*,ghats,viscA,viscAE,dkm1,bldepth}_rank0.txt`) — the M5a-2/3/4
  reference data is ALREADY GENERATED at `/scratch/a/a270088/kppinit_oracle/kpp_oracle/`. KPP+GM+Redi+linfs runs
  cleanly in the oracle at 1-rank.
- **M5a-1 ✅ DONE (byte-exact, first try):** NEW `src/oce/oce_mixing_kpp.F90` (`oce_mixing_kpp_init` builds Vtc/cg/
  deltaz/deltau + the wmt/wst lookup tables; `wscale`) + NEW driver `src/drivers/fesom_kppdump.F90` (emits the two
  txt files in the oracle's exact es24.16 format, Ricr=0.3/concv=1.6). FESOM3 vs oracle `kpp_init_rank0.txt`
  (429949 lines) + `kpp_wscale_rank0.txt` (20302 lines) are **IDENTICAL** (`diff` clean = `max|Δ|=0`; es24.16 =
  exact double round-trip). Validates the constant arithmetic, the `**1/3`/`1/4`/`1/2` table build (same Intel libm),
  + the wscale INT()/clamp/interp. NEW-FILE clean rebuild OK; no existing source touched (no regression possible).
- **M5a-2 ✅ DONE (byte-exact, first try):** `ri_iwmix` (interior Ri-mixing → `viscA`/`diffK` T&S) appended to the
  CORE2 pressure gate (option (b), the M2-M4 house method — RESOLVED the harness question). It reuses the pressure
  shim's already-prescribed strong-shear `uvnode` (FESOM3:556 == oracle:476, "Ri factor spans [0,~0.7]") + smoothed
  `bvfreq`. **`ri_viscA`/`ri_diffKt`/`ri_diffKs` all `max|Δ|=0`; all 70 pre-existing pressure fields still PASS;
  ctest 13/13** (`tools/run_pressure_gate_core2.sh` now 73 fields). Non-trivial: `ri_viscA` max=5.10e-3 (=visc_sh_limit
  +A_ver ⇒ frit=1 reached), `ri_diffKt/s` max=5.01e-3 — the full shear-instability shape range exercised. KEY: the L29
  divide is kept scalar via local POINTERS (`vA/dK` like `oce_mixing_pp`), and `AMAX1`/`AMIN1` kept VERBATIM (REAL
  intrinsics — Fortran↔Fortran is exact; the C port's 1e-9 drift was from mapping them to fmax). `ri_iwmix` loops
  OWNED nodes (`nNodO`, the oracle's `myDim_nod2D`); optional-`partit` from the start (M5d-ready). Changed: FESOM3
  `oce_mixing_kpp.F90` (+ri_iwmix), `mod_param_phys.F90` (+`diff_sh_limit`), `fesom_pressuredump.F90` (call+dump);
  oracle `oce_ale_mixing_kpp.F90` (`public ri_iwmix`), `fesom_pressure_dump.F90` (call+dump). Both rebuilt; the
  pressure gate's 70-field no-regression confirms the oracle rebuild didn't shift the dynamical core.
- **M5a-3 ✅ DONE (byte-exact, first try, 2026-06-24):** the KPP driver prestep (dVsq/ustar/Bo) + `dbsfc`(pressure_bv)
  + `bldepth` (the C-port HIGHEST-RISK routine: bulk-Ri accumulation + sw interp + ekman/monob limit) appended to BOTH
  pressure gates. **All 11 new fields `max|Δ|=0` on pi (3140) AND CORE2 (126858): `kpp_stress`/`kpp_sw3d` (inputs) +
  `kpp_dVsq`/`kpp_dbsfc`/`kpp_ustar`/`kpp_Bo` (prestep+dbsfc) + `kpp_hbl`/`kpp_kbl`/`kpp_bfsfc`/`kpp_stable`/`kpp_caseA`
  (bldepth outputs).** NON-VACUOUS: hbl spans 0.13→1180 m, BOTH forcing branches fire (57485 stable / 69373 unstable
  nodes on CORE2). Implementation:
  - **`bldepth`** ported into `src/oce/oce_mixing_kpp.F90` (explicit-dataflow args, optional-`partit`, owned-node loops;
    `smooth_hbl=.false.` ⇒ no exchange). All divides are scalar (inner nz-loop EXITs at the Rib crossing + calls
    `wscale`) ⇒ NO L29 SIMD trap; `SIGN`/`AMIN1` kept VERBATIM. Calls the M5a-1 `wscale`/tables.
  - **`dbsfc`** added as an OPTIONAL output of `pressure_bv` (`oce_pressure_bv.F90`): a separate `!DIR$ NOVECTOR` loop
    (scalar divide, matching the oracle's db_max-reduction density loop), guarded by `present(dbsfc)` so the
    lifecycle/step callers (dbsfc absent) are byte-neutral. The gate driver passes it on the smoothed `pressure_bv`
    call. Faithful to FESOM2 `oce_ale_pressure_bv.F90:326-339` (`-g*(rho_surf-rho_full)/rho_full`).
  - **Harness (the M5a-2 method, option (b)):** BOTH shims prescribe surface fluxes (`heat_flux`=200·.../`water_flux`/
    `stress_node_surf`) + a physical decaying `sw_3d` (exp(zbar/20), surface ~2.4e-5 K m/s) ANALYTICALLY (byte-identical
    formula, M3b style; NOT `cal_shortwave_rad` — that wiring is M5c) so Bo spans both signs. Compute dVsq/ustar/Bo
    inline (the oce_mixing_KPP prestep, verbatim both sides), then call `bldepth`. Oracle: `bldepth`+`dVsq`/`ustar`/
    `bfsfc`/`stable`/`caseA`/`kbl` made PUBLIC; `dbsfc` filled by enabling `mix_scheme_nmb=1` + `oce_mixing_kpp_init`
    (allocates the KPP arrays) around the smoothed `pressure_bv`; `sw_3d` allocated in-shim (forcing-init is later).
  - ⚠️ **META-LESSON (uninit-memory gate fragility):** the oracle shim allocated `fer_tapfac` (M4a) WITHOUT init →
    its below-bottom region was uninitialised garbage that COINCIDENTALLY matched FESOM3's 0.0. My new
    `oce_mixing_kpp_init` (allocated ~11 arrays before the `fer_tapfac` alloc) shifted the heap → `fer_tapfac`
    below-bottom became 1.0 on pi (CORE2 still 0.0). FIX: `fer_tapfac = 0.0_WP` after the alloc (matches FESOM3,
    deterministic). The consumed fields (`slope_tapered`/`neutral_slope`/`sigma_xy`) all stayed `max|Δ|=0` — it was
    a dump-only below-bottom artifact, not a physics bug. **Lesson: any shim array that's dumped+gated must be
    explicitly initialised below-bottom (don't rely on fresh-alloc memory); a later allocation can shift the heap.**
  - **No regression:** pi pressure gate (now incl. fer_tapfac fix) + CORE2 pressure gate (84 fields) + 1-rank step-65
    + CORE2 lifecycle 195 + ctest 13/13 all `max|Δ|=0`/green. FESOM3 changed: `oce_mixing_kpp.F90` (+bldepth),
    `oce_pressure_bv.F90` (+dbsfc), `fesom_pressuredump.F90` (M5a-3 block). Oracle (uncommitted): `oce_ale_mixing_kpp.F90`
    (public bldepth+arrays), `fesom_pressure_dump.F90` (M5a-3 block + fer_tapfac init) + `libfesom.so` rebuilt.
- **M5a-4 ✅ DONE (byte-exact, first try, 2026-06-24) → M5a COMPLETE (isolated KPP module gate closed).** `blmix_kpp`
  (BL mixing coeffs `blmc(3)` + `dkm1(3)` + nonlocal `ghats`, eqn 10/11/20: T uses diffK(1)→blmc(2), S diffK(2)→blmc(3),
  mom diffK→blmc(1)) + `enhance` (kbl-1 interface blend) + the driver tail (combine: within-BL `max(.,blmc)`, outside
  `ghats=0`; node→elem viscAE average + `minmix=3e-3` surface floor) appended to BOTH pressure gates. **All 11 new
  fields `max|Δ|=0` on pi (3140) AND CORE2 (126858): `kpp_blmc1/2/3` + `kpp_dkm1m/t/s` + `kpp_ghats` + final
  `kpp_viscA`/`kpp_diffKt`/`kpp_diffKs` + `kpp_viscAE`** (the pressure gates now total **95 fields** each). Implementation:
  - **`blmix_kpp` + `enhance`** ported into `src/oce/oce_mixing_kpp.F90` (explicit-dataflow args — `blmc`/`ghats`/`dkm1`
    passed in, vs the oracle's module arrays — optional-`partit`, owned-node loops; `blmc` pre-zeroed over owned+halo
    like the oracle). They read the M5a-2 `viscA_ri`/`diffK_ri` (interior) + the M5a-3 `hbl`/`kbl`/`bfsfc`/`stable`/
    `caseA`/`ustar` — all byte-proven, so blmix consumes them directly. **NO L29 NOVECTOR needed** (first try): every
    blmix divide (`/dthick`, `/(hbl+epsln)`, `/(wm+epsln)`, `/(ws+epsln)`) is SCALAR per-node — the inner nz-loop EXITs
    at kbl + calls `wscale`, and enhance's `delta` is one scalar divide; `AMIN1`/`MIN`/`INT`/`ABS`/`SIGN` kept VERBATIM.
  - **Combine + node→elem viscAE average** done INLINE in BOTH shims (the M5a-3 prestep pattern), operating on SEPARATE
    copies `viscA_fin`/`diffK_fin` so the M5a-2 `ri_viscA`/`ri_diffKt`/`ri_diffKs` dumps still echo the interior (the
    combine overwrites viscA/diffK in place in the real driver; the copies isolate the two gate points). `viscAE` is the
    `pp_Av`-style element field (`SUM(viscA_fin(nz,elnodes))/3`, gid-order byte-exact at 1-rank); `kpp_ghats` is the
    post-combine nonlocal flux (== post-enhance here: blmix only writes nz<kbl, combine only zeros nz≥kbl which were
    already 0). Oracle: `blmix_kpp`/`enhance`/`dkm1` made PUBLIC (like `bldepth`); the shim calls them on `viscA_ri`/
    `diffK_ri` (`blmc`/`ghats`/`dkm1` already allocated+zeroed by the M5a-3 `oce_mixing_kpp_init`).
  - **No regression:** pi pressure gate (95) + CORE2 pressure gate (95) + ctest 13/13 + step-65 1-rank (65) all
    `max|Δ|=0`/green; the heap-shift-sensitive `tsol_*` + `fer_tapfac` (L:M5a-3) stayed byte-exact on BOTH meshes.
    FESOM3 changed: `oce_mixing_kpp.F90` (+blmix_kpp/+enhance), `fesom_pressuredump.F90` (M5a-4 block + 11 dumps).
    Oracle (uncommitted): `oce_ale_mixing_kpp.F90` (public blmix_kpp/enhance/dkm1), `fesom_pressure_dump.F90` (M5a-4
    block + 11 dumps) + `libfesom.so` rebuilt.
- **M5b ✅ DONE (2026-06-24) — KPP wired into `step_oce`; UNFORCED CORE2 lifecycle byte-exact.** `max|Δ|=0` KPP-alone
  195 (3-step) + 325 (5-step) AND KPP+GM+Redi 195 (`tools/run_lifecycle_kpp_gate_core2.sh [run] [nsteps]`, env
  `FER_GM=1`/`REDI=1`). NEW `oce_mixing_kpp_driver` (renamed from `oce_mixing_KPP` to dodge the case-insensitive clash
  with the module name) + KPP `dyn%work` arrays (allocated only when KPP, not serialized) + `mix_scheme_nmb==1` dispatch
  in `mod_step_oce` (pressure_bv `dbsfc=`; `sw_alpha_beta` fires for `Fer_GM.or.Redi.or.is_kpp`; `Kv=Kv_double(:,:,1)`
  over `nNodL` then `mo_convect`; `Av`/viscAE element-based → `impl_vert_visc` UNCHANGED). `step_oce` gained an OPTIONAL
  `stress_node_surf(2,nod)` (KPP ustar; unforced⇒0, passed unallocated for PP⇒seen absent). The reduced-M2 drivers now
  pin `mix_scheme_nmb=2` (the new dispatch defaults to KPP=1). `FESOM3_MIX_KPP` env in `fesom_lifecycle` (mix_scheme_nmb=1
  + `use_sw_pene=.false.` + work_core DOUBLES + alloc + `oce_mixing_kpp_init`); oracle `run_lifecycle_core2.sh` `MIX_KPP=1`
  KEEPS `mix_scheme='KPP'` (NO oracle source change). **ROOT-CAUSE (L45): the ONE step beyond M5a-4 was
  `smooth_blmc=.true.`** (oracle `:439-449`: 3 area-weighted `smooth_nod` sweeps × the 3 `blmc` channels, BEFORE the
  combine) — the M5a-4 isolated gate called `blmix`/`enhance`/combine DIRECTLY so never exercised it; omitting it diverged
  ~1e-3 in Kv at the kbl-1 BL level. Fixed by reusing the M2.1 byte-proven `smooth_nod` (made public in `oce_pressure_bv`;
  == FESOM2 `gen_support.F90` `smooth_nod3D`). **No regression:** PP step-65 1-rank + dist_2 + PP lifecycle 195 +
  CORE2/pi pressure (95) + ctest 13/13 all `max|Δ|=0`. *(Historical M5b plan follows.)*
- **M5c ✅ DONE (2026-06-24, first try, L46) — forced/fully-native lifecycle + sw_pene = THE production CORE2 column
  physics byte-exact.** `max|Δ|=0` on the 13 NODE substeps × 5 probes: **195 (3-step) AND 325 (5-step), BOTH whichEVP**,
  with KPP + GM + Redi + shortwave-penetration ALL on in the fully-native forced lifecycle
  (`tools/run_lifecycle_kpp_native_gate_core2.sh [run] [nsteps] [whichEVP]`). The per-step native-vs-oracle flux
  self-check reads `max|Δ(hf,wf,vs,rs,ss)|=0` (so the `cal_shortwave_rad` heat_flux modification matches too).
  - **KEY FINDING (L46): only `sw_3d` is live in work_core; the `ghats` nonlocal flux is DEAD.** The plan named two new
    tracer-TDMA terms, but the oracle (`oce_ale_tracer.F90:892`) guards ghats with `if (use_kpp_nonlclflx)` FIRST —
    and `use_kpp_nonlclflx` defaults `.false.`, is never set `.true.` in the oracle, and is ABSENT from work_core (like
    MLD1_ind/K_hor in M4). So the production M5c term is `sw_3d` alone.
  - **`sw_3d` term** added to `diff_ver_part_impl_ale` (`oce_ale_tracer.F90`, FESOM2 `:991-996`, guarded `use_sw_pene
    .and. id==1`); reads `dyn%work%sw_3d`. **`cal_shortwave_rad`** (M2.10c, byte-proven) wired into the native runloop
    right after `oce_fluxes` (it IS the last line of FESOM2 `oce_fluxes`, `ice_oce_coupling.F90:873`; the oracle flux
    dump `fesom_module.F90:746` records `heat_flux` AFTER it) — fills `dyn%work%sw_3d` + adds the visible band back to
    `atm%heat_flux`. chl = constant 0.1 (work_core `chl_data_source='None'`, `chl_const=0.1`); albw = `ice%thermo%albw`
    (0.1). The native driver gained `FESOM3_MIX_KPP`/`FESOM3_SW_PENE` (+ the M5b KPP `dyn%work` alloc block); passes
    `atm%stress_node_surf` (oce_fluxes_mom output) to `step_oce` for KPP ustar.
  - **⚠️ `use_sw_pene` was a latent DUPLICATE.** It already lived in `mod_config` (the g_config home, read by the KPP
    `bldepth` + `fesom_lifecycle` + `fesom_pressuredump`); I'd briefly added a second to `mod_param_phys`. Consolidated
    to the single `mod_config` source, default flipped `.true.`→`.false.` ("unset=off"): the gate drivers configure
    in-code (no namelist read) and `dyn%work%sw_3d` is allocated only with sw_pene, so a `.true.` default would fire the
    term against unallocated memory in every PP driver.
  - **ghats transcribed + gated anyway** (faithful, complete port; guarded by `use_kpp_nonlclflx`): the `KPP_NONLCL=1`
    variant (`tools/run_lifecycle_kppnonlcl_native_gate_core2.sh`, injects `use_kpp_nonlclflx=.true.` into the oracle
    `&tracer_phys` + the FESOM3 `FESOM3_KPP_NONLCL` env) is `max|Δ|=0` BOTH whichEVP and NON-VACUOUS (ghats moves T/S by
    up to 7.3e-5 over 3 steps vs ghats-off) — so the dead term is byte-verified, not an ungated landmine (L39). `ref_sss`
    /`ref_sss_local` added to `mod_param_phys` (the ghats salinity branch `rsss`).
  - **No regression:** ctest 13/13 + M5b unforced KPP (195) + M4e GM+Redi native (195) + M3f baseline fully-native (195)
    all `max|Δ|=0`. Changed: `oce_ale_tracer.F90` (sw_3d + ghats terms), `mod_param_phys.F90` (+use_kpp_nonlclflx/ref_sss),
    `mod_config.F90` (use_sw_pene default), `fesom_lifecycle_native.F90` (KPP/sw_pene wiring), the forced-oracle +
    fully-native gate runners (+MIX_KPP/SW_PENE/KPP_NONLCL env), + 2 new gate scripts. No oracle SOURCE change.
- **M5d ✅ DONE (2026-06-24, first try — the M4f lesson) → M5 COMPLETE.** The production KPP+sw_pene+GM+Redi
  fully-native FORCED lifecycle byte-matches FESOM2 `max|Δ|=0` at MULTI-RANK: CORE2 **dist_2 + dist_8, BOTH whichEVP,
  195 (3-step) + 325 (5-step)** — all 8 (dist × EVP × nsteps) cells green (`tools/run_lifecycle_kpp_native_gate_multirank.sh
  [np] [nsteps] [whichEVP]`). PURE WIRING: the KPP module (`oce_mixing_kpp.F90`) was written optional-`partit` from the
  START — `oce_mixing_kpp_driver` already does owned-loop bounds + `exchange_nod(blmc 1/2/3, Kv_double 1/2, ghats, viscA)`
  + the 3-sweep `smooth_blmc` (per-sweep `smooth_nod` exchange), and the sw_3d/ghats tracer terms loop `nNodO` (owned).
  The only NEW work: (i) the `fesom_lifecycle_native_mr` KPP/sw_pene/nonlcl config block (LOCAL-sized `nNodL`/`nElemF`,
  mirror of the M5c 1-rank block + the M4f GM/Redi block) + `cal_shortwave_rad` after `oce_fluxes` (partit) +
  `stress_node_surf=atm%stress_node_surf` to `step_oce`; (ii) **`cal_shortwave_rad` gained an optional `partit`** —
  ⚠️ in the LOCAL mesh `mesh%nod2D` = the GLOBAL count (`read_mesh_local`), so it can't be the loop bound; now loops
  `nNodL` (owned+halo = oracle `myDim+eDim`) via `owned_bounds` (absent partit ⇒ `nNodL=mesh%nod2D`, 1-rank/pi verbatim).
  The gate runners `run_lifecycle_fullynative_gate_multirank.sh` + the forced oracle thread MIX_KPP/SW_PENE/KPP_NONLCL.
  No regression (ctest 13/13 + 1-rank production 195 + pi M2.10c forcing-gate `sw_3d`/`chl`/`heat_flux_sw` all `max|Δ|=0`).
  Changed: `fesom_lifecycle_native_mr.F90`, `oce_shortwave_pene.F90` (optional partit), the MR gate runner, + the new
  `run_lifecycle_kpp_native_gate_multirank.sh`. No oracle SOURCE change.
- **M5b PLAN (done above):** the isolated KPP module is fully byte-proven — now ASSEMBLE the driver + wire it into the step.
  Build the real `oce_mixing_KPP` driver in `oce_mixing_kpp.F90` (prestep dVsq/ustar/Bo + `ri_iwmix` → `bldepth` →
  `blmix_kpp` → `enhance` → combine → node→elem viscAE average — all four sub-kernels already byte-proven, so this is
  assembly, not new physics) taking `viscAE`(elem)+`Kv_double`(node,ntr); add KPP arrays (`Kv_double`/`ghats`/`blmc`/
  `hbl`/`kbl`/`dbsfc`/`dVsq`/`ustar`/`Bo`/`bfsfc`/`stable`/`caseA`/`dkm1`) to `dyn%work`, allocated only when KPP.
  Dispatch in `mod_step_oce` (`oce_ale.F90:3711`): after the producers call `oce_mixing_KPP`, then `Kv(:,n)=
  Kv_double(:,n,1)` (single T-channel Kv → the tracer TDMA is UNCHANGED), then `mo_convect` (already ported). `Av`=
  viscAE stays element-based → `impl_vert_visc` UNCHANGED. `dbsfc` is already the optional pressure_bv output (M5a-3);
  wire it through `dyn%work`. Gate the UNFORCED CORE2 lifecycle `max|Δ|=0` (unforced ⇒ zero surface flux + sw ⇒
  ghats=0, sw_3d=0 ⇒ tests the KPP Kv/Av integration + mo_convect-after-KPP). M5c adds the ghats + sw_3d tracer terms
  (forced/native); M5d multi-rank (the kernels are optional-`partit` from the start → pure wiring, the M4f lesson).

**M4 (GM/Redi) ✅ COMPLETE + COMMITTED (tag `m4`).** Byte-exact `max|Δ|=0` 1-rank AND multi-rank (CORE2 dist_2/dist_8),
unforced AND forced/fully-native, GM-only AND GM+Redi, BOTH whichEVP. Gates: `run_lifecycle_gm_gate_core2.sh` (M4c),
`run_lifecycle_redi_gate_core2.sh` (M4d), `run_lifecycle_gmredi_native_gate_core2.sh` (M4e),
`run_lifecycle_gmredi_native_gate_multirank.sh` (M4f). M4 enabled by env: FESOM3 drivers read `FESOM3_FER_GM`/
`FESOM3_REDI`; oracle runners read `FER_GM=1`/`REDI=1` to KEEP work_core `Fer_GM`/`Redi=.true.`.

**M4f KEY FINDING (the M4 capstone):** the M4 routines were ALREADY MULTI-RANK READY — every kernel was transcribed
(M4a–M4d) with the M2.12/M3 optional-`partit` pattern (owned-loop bounds + the FESOM2 halo exchanges), and
`mod_step_oce::step_oce` already threads the optional `partit` to all of them. So M4f was a pure WIRING task (the
GM/Redi config + LOCAL-sized array allocation into `fesom_lifecycle_native_mr` + the gate env) — the partit-present
path that M4a–M4d wrote but never gated just worked, FIRST try, NO halo extension, NO `dist_N` invariant check needed
(the GM/Redi producers reuse the SAME mesh ops — `nod_in_elem2D`+`gradient_sca` over owned elements, owned-edge loops —
the M2.12 dynamics already proved invariant). `sw_alpha_beta`/`tracer_gradient_z` compute over owned+halo (no exchange);
`compute_sigma_xy`/`compute_neutral_slope`/`init_Redi_GM`/`fer_solve_Gamma`/`fer_gamma2vel`/`vert_vel_ale` exchange
their outputs; the Redi diff terms loop owned edges/nodes reading exchanged/halo-valid inputs.

**M5 = KPP vertical mixing (the real production CORE2 config: `mix_scheme='KPP'` instead of the reduced-M2 PP).**
The production `work_core` runs KPP + Fer_GM + Redi + sw_pene + `which_ALE='zlevel'`. M4 added GM+Redi; M5 swaps PP→KPP
and turns on `use_sw_pene` (the shortwave penetration term in the tracer TDMA, currently `.false.` to dodge the
unallocated `sw_3d`). Mirror the M4 method: port `oce_ale_mixing_kpp.F90` (the KPP boundary-layer + interior scheme),
gate `max|Δ|=0` vs FESOM2 with `mix_scheme='KPP'`, then the forced/fully-native lifecycle, then multi-rank. The oracle
has work-dirs staged for the future targets (`work_zstar_kpp`/`work_tke_dump`). See "Milestone ordering" below + the
plan doc `2026-06-18-fesom3-architecture.md`.

⚠️ multi-rank levante needs the `env.sh` KNEM flag (L35). ⚠️ Gate every consumed intermediate (L29). ⚠️ **Read the
ACTUAL FESOM2 `.F90` + `work_core` namelists, NOT the plan summaries** (they were wrong on `MLD1_ind` + `K_hor` in M4).
M3 ✅ COMMITTED (`b591153`, tag `m3`); M4 ✅ COMMITTED (`7a088de`, tag `m4`); the M2 baseline stays at `m2-mvp`.
**M5 ✅ COMPLETE + COMMITTED — tag `m5` (`5c5a03d`). M5a+M5b (`30760c3`), M5c forced/native + sw_pene (`2fd05b8`),
M5d multi-rank + M5 capstone (`5c5a03d`, tag `m5`). The M2 baseline stays at `m2-mvp`.**

**Milestone ordering (plan doc `2026-06-18-fesom3-architecture.md`):** M4 = GM/Redi → M5 = KPP + production
multi-year (paper-parity) → M6 = beyond-paper (zstar/zlevel ALE, TKE, aEVP). mEVP (whichEVP=1) was already
delivered in M3. **The real production CORE2 config (`work_core`) runs KPP + Fer_GM + Redi + sw_pene + which_ALE='zlevel'**
— the reduced-M2 gate sed-downgrades ALL of these (linfs/PP/no-GM/no-Redi/no-sw_pene, `run_lifecycle_forced_core2.sh:44-53`);
the oracle has work-dirs staged for the future targets (`work_zstar_dump`/`work_zstar_kpp`/`work_zstar_tke`/`work_tke_dump`).

**M4 = GM/Redi (turn ON `Fer_GM`+`Redi` over the existing PP+linfs core; gate `max|Δ|=0` vs FESOM2 with the matching
config).** Decomposed M4a–M4f (full detail + oracle map + SIMD-divide watch-list in the plan file):
- **M4a ✅ DONE** — producers (`mesh_resolution` compute, `sw_alpha_beta`, `compute_sigma_xy`,
  `compute_neutral_slope`→`slope_tapered`/`fer_tapfac`); feed-forward, the safest first gate. (`MLD1_ind` turned out
  DEAD under work_core — see M4b.)
- **M4b ✅ DONE** — GM diffusivity + streamfunction + bolus velocity (new `oce_fer_gm.F90`: `init_Redi_GM`,
  `fer_solve_Gamma` TDMA, `fer_gamma2vel`; `fer_w` in `vert_vel_ale`). `Fer_GM=T`, `Redi=F`. 70-field pressure gate
  `max|Δ|=0`, first try.
- **M4c ✅ DONE** — GM bolus into tracer advection (add/subtract `fer_uv`/`fer_w` in `solve_tracers_ale` + GM chain
  wired into `step_oce`); unforced lifecycle gate `max|Δ|=0` (195+325 records). `Fer_GM=T`, `Redi=F`.
- **M4d ✅ DONE** — Redi isopycnal diffusion (`Ki` real field; K13/K23 in `diff_part_hor_redi`, NEW
  `diff_ver_part_redi_expl` K31/K32, K33 in `diff_ver_part_impl_ale`). Gated `Fer_GM=T, Redi=T` (the plan's `Fer_GM=F`
  is vacuous — work_core `K_hor=0`); unforced lifecycle `max|Δ|=0` (195+325). Merges with M4e-unforced.
- **M4e ✅ DONE** — both on in the FORCED / fully-native lifecycle (production tracer physics, + ice + native forcing).
  `Fer_GM=T`, `Redi=T`. Byte-exact 195+325 records, BOTH whichEVP (`tools/run_lifecycle_gmredi_native_gate_core2.sh`).
- **M4f ✅ DONE** — multi-rank (the M4 routines were ALREADY MR-ready; pure wiring + gating). CORE2 dist_2/dist_8,
  195+325 records, BOTH whichEVP, `max|Δ|=0` (`tools/run_lifecycle_gmredi_native_gate_multirank.sh`). Tag `m4`.

**M4 ✅ COMPLETE + COMMITTED (tag `m4`). M4a producers + M4b GM diffusivity/streamfunction/bolus-velocity + M4c
GM-bolus-into-advection + M4d Redi isopycnal diffusion + M4e GM+Redi in the FORCED/fully-native lifecycle + M4f
multi-rank ✅ DONE + byte-gated 2026-06-23/24 (`max|Δ|=0`, CORE2 1-rank AND dist_2/dist_8). ACTIVE NEXT = M5 (KPP).**
- **M4a** — `mesh_resolution` (`mod_mesh_areas.F90::compute_mesh_resolution`, rides the geom gate → 20 fields) +
  `sw_alpha_beta`/`compute_sigma_xy`/`compute_neutral_slope` (`oce_pressure_bv.F90`): the 6 fields
  `sw_alpha`/`sw_beta`/`sigma_xy`/`neutral_slope`/`slope_tapered`/`fer_tapfac` appended to the CORE2 pressure gate.
  Both harnesses force the `work_core` taper config around the producer calls (`Fer_GM=Redi=Redi_Ktaper=.true.`,
  `ODM95_Scr=0.2e-2` — the work_core override, not the 1e-2 default), then restore; the oracle snapshots
  `slope_tapered` before the shim zeros it at `:483`.
- **M4b** — new `src/oce/oce_fer_gm.F90` (`init_Redi_GM` + `fer_solve_Gamma` TDMA + `fer_gamma2vel`, explicit-dataflow
  + optional-`partit`) + `fer_w` added to `vert_vel_ale` (`oce_ale.F90`, `if(Fer_GM)`-guarded → byte-neutral at
  `Fer_GM=.false.`). The 6 fields `fer_K`/`fer_c`/`fer_scal`/`fer_gamma`/`fer_uv`/`fer_w` appended to the CORE2
  pressure gate → **70 PASS / 0 FAIL** (FIRST try, NO NOVECTOR needed). FESOM3 `fesom_pressuredump` + oracle
  `fesom_pressure_dump` each force `Fer_GM=T, Redi=F` + work_core GM config (`K_GM_max=1000`/`min=2`/`cm=3`/`cmin=0.1`/
  `resscalorder=2`/`rampmax=rampmin=-1`/`Ktaper=F`/`scaling_resolution=T`/`scaling_GMzexp=T`/`zref=500`/`smin=0.6`)
  around the GM calls + the existing `vert_vel_ale` (so it fills `fer_w`), then restore; fer arrays init'd to the
  FESOM2 `oce_setup_step.F90:962-967` values (`fer_K=500`/`fer_c=1`/`fer_scal=0`/`fer_gamma=0`, `fer_uv/fer_w=0`) so
  below-bottom byte-matches. **KEY: `MLD1_ind` is DEAD under work_core** (read ONLY in `init_Redi_GM`'s
  `scaling_Ferreira` branch; work_core uses `scaling_GMzexp` → GMzexp `exp(-|zbar_3d_n|/zref)`, needs no MLD1_ind/
  bvref) → NOT implemented (L33, the plan's "fold into M4b" assumed Ferreira); `fer_tapfac`/`neutral_slope` likewise
  unconsumed by GM (`K_GM_Ktaper=F`, `scaling_FESOM14=F`). Only the work_core GM path is ported — `init_Redi_GM`
  `error stop`s on `scaling_Ferreira/Rossby/GINsea/FESOM14/K_GM_Ktaper/Redi` (the **Redi Ki path is M4d**, which
  removes Redi from the guard + passes a `Ki` arg). `zbar_n_bot(n)` (absent in FESOM3's mesh) ==
  `mesh%zbar_3d_n(nlevels_nod2D(n),n)` (FESOM2 `oce_ale.F90:550`); `fer_solve_Gamma` mirrors the oracle's
  `tr => fer_gamma(:,:,n)` POINTER (dodges the L29 SIMD-divide trap). No regression (pressure 70 + step-65 1-rank +
  step-65 MR dist_2 + ctest 13/13 `max|Δ|=0`).
- **M4c** — GM bolus wired INTO the step + lifecycle (`mod_step_oce`: producers `sw_alpha_beta`→`compute_sigma_xy`
  after PGF, then `init_Redi_GM`→`fer_solve_Gamma`→`fer_gamma2vel` after `update_eta_n` before `vert_vel_ale`; bolus
  add/subtract in `solve_tracers_ale` — once around the tracer loop over owned+halo, `uv += fer_uv`/`w_e += fer_w`/
  `w += fer_w` before advection, subtracted before the salinity clamp; all `if(Fer_GM)`). GM aux fields added to
  `t_dyn_work` (non-serialized; allocated only when `Fer_GM`). The UNFORCED CORE2 lifecycle with GM ON is **`max|Δ|=0`
  — 195 records (3-step) AND 325 records (5-step)** (`tools/run_lifecycle_gm_gate_core2.sh`, FIRST try). Gate: the
  oracle base namelist IS work_core (GM params already correct) so `run_lifecycle_core2.sh` gained `FER_GM=1` to KEEP
  `Fer_GM=.true.` (Redi sed off); FESOM3 `fesom_lifecycle` gained `FESOM3_FER_GM` (sets `Fer_GM` + work_core GM config
  + allocates/inits the GM arrays). NO oracle source change. `compute_neutral_slope` SKIPPED (Redi/Ktaper-only,
  unconsumed at M4c → M4d). `fer_w` into BOTH `w` and `w_e` (faithful: the non-FCT HO vertical advection reads `w_e`; it
  equals `w` when the split is off). No
  regression: GM-OFF lifecycle 195 `|Δ|=0` (the `Fer_GM`-guarded edits are byte-neutral), step-65 1-rank `|Δ|=0`,
  ctest 13/13.
- **M4d** — Redi isopycnal diffusion, gated `Fer_GM=T + Redi=T` (production GM+Redi) UNFORCED lifecycle: **`max|Δ|=0`,
  195 (3-step) + 325 (5-step) records** (`tools/run_lifecycle_redi_gate_core2.sh`, FIRST try, NO NOVECTOR). **KEY (from
  the actual Fortran): the plan's "Fer_GM=F, Redi=T" isolating config is VACUOUS** — work_core `namelist.tra` sets
  `K_hor=0` so `init_Redi_GM` F1 `Ki(nzmin)=K_hor**(..)=0`; the only non-zero `Ki` is the `if(Redi.and.Fer_GM)` coupling
  `Ki=max(fer_scal*Redi_Kmax,K_GM_min)` (needs Fer_GM=T). So M4d gates BOTH on (merges plan M4d+M4e-unforced); M4c
  already pinned GM-only so a pass isolates the Redi terms. Implemented: `init_Redi_GM` Redi path (Redi removed from the
  guard; `Redi_Kmax<=0` sync; F1/coupling/F2 Ki + `Redi_Ktaper` sqrt-split; `Ki`/`fer_tapfac` optional args);
  `compute_neutral_slope` into `step_oce` (`if(Fer_GM.or.Redi)`); `diff_part_hor_redi` K13/K23 (`Tz`/`SxTz`/`SyTz`,
  `replace_all` across the 5 level-ranges); NEW `diff_ver_part_redi_expl` K31/K32; `diff_ver_part_impl_ale` K33
  (`slope_tapered(3)^2*Ki` augments `Kv` in the tridiag); `tracer_gradient_z`→`tr_z` (new `t_tracer_work` field).
  Redi aux (`Ki`/`neutral_slope`/`slope_tapered`/`fer_tapfac`) in `t_dyn_work`, allocated only when Redi; ALL Redi reads
  `if(Redi)`-guarded so the off configs (arrays unallocated) are byte-neutral (`-(Kv+0)*..==-Kv*..` by IEEE +0). Gate:
  `run_lifecycle_core2.sh` `REDI=1` keeps work_core `Redi=.true.`; `fesom_lifecycle` `FESOM3_REDI` enables it. NO oracle
  source change; Redi non-vacuous (max|uv| differs from GM-only). No regression: M4c GM-only 195 `|Δ|=0`, GM-OFF 195
  `|Δ|=0`, step-65 `|Δ|=0`, ctest 13/13.
- **M4e** — GM+Redi in the FORCED / fully-native lifecycle (production tracer physics + ice + native forcing): **`max|Δ|=0`,
  195 (3-step) + 325 (5-step) records, BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP)**
  (`tools/run_lifecycle_gmredi_native_gate_core2.sh`, FIRST try, 2026-06-24). **Pure WIRING — NO new physics code.** The
  M4a–M4d kernels (`step_oce`'s `Fer_GM`/`Redi`-guarded chain) are unchanged; M4e just enables them in the native run:
  (1) FESOM3 `fesom_lifecycle_native.F90` got the `FESOM3_FER_GM`/`FESOM3_REDI` env reads + work_core GM+Redi config +
  `dyn%work` GM/Redi array allocation block, COPIED VERBATIM from `fesom_lifecycle.F90`; `step_oce` (1-rank, no `partit`)
  reads the module flags + consumes the `dyn%work`/`fer_uv`/`tr_z` arrays → the GM+Redi chain runs on the FORCED state
  with no further change (first GM/Redi exercise on a non-zero-flux ocean). (2) Oracle `run_lifecycle_forced_core2.sh`
  gained `FER_GM`/`REDI` env (mirror `run_lifecycle_core2.sh`) — `=1` KEEPS work_core `Fer_GM`/`Redi=.true.`; NO oracle
  source change. (3) `run_lifecycle_fullynative_gate_core2.sh` gained `FER_GM`/`REDI` env (passes to oracle + exports
  `FESOM3_FER_GM`/`FESOM3_REDI`); thin wrapper `run_lifecycle_gmredi_native_gate_core2.sh` (FER_GM=1 REDI=1). The native
  forcing + flux self-checks print `max|Δ|=0` at every step under GM+Redi. No regression: M3f GM/Redi-OFF fully-native
  195 `|Δ|=0` (the shared-oracle-runner + gate-script edits are byte-neutral OFF), M4d unforced GM+Redi 195 `|Δ|=0`,
  ctest 13/13. Build incremental (only `fesom_lifecycle_native` relinked — no core-physics/new-file change).
- **M4f** — multi-rank GM+Redi (production tracer physics at MULTI-RANK): **`max|Δ|=0`, CORE2 dist_2 AND dist_8,
  195 (3-step) + 325 (5-step) records, BOTH whichEVP** (`tools/run_lifecycle_gmredi_native_gate_multirank.sh`,
  FIRST try, 2026-06-24). **Pure WIRING — NO new physics, NO halo extension, NO `dist_N` invariant check needed.** The
  M4 routines were ALREADY MR-ready (transcribed M4a–M4d with the M2.12/M3 optional-`partit` pattern; `step_oce` threads
  `partit` to all of them — the partit-present path was written but never gated). The exchanges already in place:
  `compute_sigma_xy`→`exchange_nod(sigma_xy)` (the rank-3 superset of the oracle's 3×`MPI_BARRIER` per-component dance =
  same bytes), `compute_neutral_slope`→`exchange_nod(neutral_slope/slope_tapered)`, `init_Redi_GM`→`exchange_nod(fer_c/
  fer_K/Ki)`, `fer_solve_Gamma`→`exchange_nod(fer_gamma)`, `fer_gamma2vel`→`exchange_elem(fer_uv)`, `vert_vel_ale`→
  `exchange_nod(fer_w)`. `sw_alpha_beta`+`tracer_gradient_z` compute owned+halo (no exchange — the compute-at-halo idiom).
  The Redi diff terms loop OWNED edges/nodes reading exchanged `slope_tapered`/`Ki` + halo-valid `tr_z` (reusing the
  M2.12c-3 owned-edge invariant). Only edit: the GM/Redi config + LOCAL-sized (nNodL/nElemF) array allocation block into
  `fesom_lifecycle_native_mr.F90` (`fer_uv` sized nElemF like `dyn%uv`) + the gate env. No regression: GM/Redi-OFF MR
  fully-native dist_2 195 `|Δ|=0`, ctest 13/13.
**M4 ✅ COMMITTED (tag `m4`, M4a–M4f together).** FESOM3 changed: `oce_fer_gm.F90` (NEW)/`oce_ale.F90`/
`oce_ale_tracer.F90`/`oce_pressure_bv.F90`/`oce_tracer_grad.F90`/`mod_dyn.F90`/`mod_tracer.F90`/`mod_step_oce.F90`/
`mod_mesh_areas.F90`/`mod_geom_dump.F90`/`fesom_pressuredump.F90`/`fesom_lifecycle.F90`/`fesom_lifecycle_native.F90`/
`fesom_lifecycle_native_mr.F90` + `tools/run_lifecycle_core2.sh`/`run_lifecycle_gm_gate_core2.sh`/
`run_lifecycle_redi_gate_core2.sh`/`run_lifecycle_forced_core2.sh`/`run_lifecycle_fullynative_gate_core2.sh`/
`run_lifecycle_fullynative_gate_multirank.sh`/`run_lifecycle_gmredi_native_gate_core2.sh` (NEW)/
`run_lifecycle_gmredi_native_gate_multirank.sh` (NEW). Oracle (UNCOMMITTED working-tree, by design):
`fesom_pressure_dump.F90`/`fesom_geom_dump.F90` + `libfesom.so` rebuilt.
**NEXT: M5 — KPP vertical mixing + production multi-year (paper-parity).** Plan:
`docs/plans/2026-06-23-m4-gm-redi.md` (M4 log) + `2026-06-18-fesom3-architecture.md` (milestone ordering). The PRIOR
resume target follows.

**↳ (prior) M3f — whole ice step + forced lifecycle (NATIVE fluxes) + multi-rank, tag `m3`.**
M3a + M3b + M3c + M3d + M3e are ALL DONE — the ice foundation, the EVP dynamics (ocean2ice + standard EVP + mEVP), the
ice FCT advection (`ice_TG_rhs` + `ice_fct_solve`), the ice thermodynamics (`cut_off` + `thermodynamics` →
`flx_h`/`flx_fw`/`t_skin`), AND the air-sea coupling-out (`oce_fluxes_mom` + `oce_fluxes` →
`heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`/`stress_surf`) all byte-match FESOM2 `max|Δ|=0`. The M3 oracle
machinery is PROVEN + EXTENDED (`fesom_ice_dump` shim now has `ice_dump_write` [M3a], `evp_dump_write` [M3b],
`ice_fct_dump_write` [M3c], `ice_thermo_dump_write` [M3d], AND `ice_flux_dump_write` [M3e]; runners
`run_ice*`/`run_evp*`/`run_icefct*`/`run_icethermo*`/`run_iceflux*`). M2 stays closed (tag `m2-mvp`); no-regression
baseline holds (13/13 ctest + step-65 1-rank/multirank + M3d-thermo-6×2 + M3e-flux-5×2 all `max|Δ|=0`). **M3f assembles
the native ice step into the lifecycle and drops the prescribed fluxes:** build a real `ice_timestep` (ocean2ice →
EVPdynamics_solve → ice_TG_rhs → ice_fct_solve → cut_off → thermodynamics) + `oce_fluxes_mom` + `oce_fluxes`, wire them
into a new forced lifecycle driver (the analog of `fesom_lifecycle` but with NATIVE fluxes replacing the
`FESOM3_FLUX_FILE` prescription), and byte-gate the multi-step forced CORE2 lifecycle vs `run_lifecycle_forced_*`
(195+ records). Then multi-rank: thread the optional `partit` + add the ice halo exchanges (the momentum-stress /
surface-flux halo reads in `oce_fluxes_mom`/`oce_fluxes` — `a_ice`/`uice`/`vice`/`srfoce`/`stress_atmoce` at the halo;
`oce_fluxes_mom` writes owned elements only so its node loop needs valid halo `stress_node_surf`; the `integrate_nod`/
`ocean_area` `allreduce_sum` is already wired guarded by `is_multirank` — and `mesh%ocean_area` needs the cross-rank
`allreduce_sum` added in `mod_mesh_areas.F90`, currently local-only, see its TODO comment). ⚠️ The real atmospheric
forcing read (NCAR bulk → `stress_atmoce_x/y` + `Ch/Ce_atm_oce_arr` + prec/runoff) must feed the native thermo/fluxes
— the M2.10 forcing read is proven, but the bulk-formula `stress_atmoce` assembly (`gen_bulk`/`fesom_forcing_dump.F90`
path) is NOT yet ported (M3f scope). ⚠️ **Any multi-rank run on levante needs the env.sh MPI flag** (sourcing `env.sh`).

Regression sanity (re-run any time — e.g. after a rebuild; all `max|Δ|=0`/green as of 2026-06-24):
```bash
./configure.sh --compiler intel --precision dp --build           # incremental (edits-only M4b–M4e);
                                                                 # --clean only when NEW files were added (L19)
bash tools/run_step_gate.sh                      # 65 records, worst |Δ|=0 (1-rank)
bash tools/run_step_gate_multirank.sh 2          # 65 records, worst |Δ|=0
bash tools/run_pressure_gate.sh                  # pi 1-rank, 95 fields, max|Δ|=0 (incl. M5a-2/3/4 kpp_* + fer_tapfac fix)
bash tools/run_pressure_gate_core2.sh            # CORE2 97 fields incl. M4a/b fer_* + M5a-2 ri_* + M5a-3/4 kpp_* + M6a-1 pgf_*_shchep (Shchepetkin zlevel PGF), max|Δ|=0
bash tools/run_lifecycle_gate_core2.sh           # GM/Redi-OFF lifecycle 195 records (no-regression baseline)
bash tools/run_lifecycle_zstar_gate_core2.sh 3   # M6a-2: zstar ALE unforced lifecycle 195, max|Δ|=0 (CORE2 1-rank)
bash tools/run_lifecycle_zstar_native_gate_core2.sh 3 0  # M6a-3: forced/native zstar 195, max|Δ|=0 (reduced; add FER_GM=1 REDI=1 MIX_KPP=1 SW_PENE=1 for production)
FER_GM=1 REDI=1 MIX_KPP=1 SW_PENE=1 bash tools/run_lifecycle_zstar_native_gate_multirank.sh 2 3 0  # M6a-4: PRODUCTION zstar MR dist_2 195, max|Δ|=0 (needs env.sh KNEM, L35)
bash tools/run_lifecycle_kpp_gate_core2.sh 3     # M5b: KPP-alone unforced lifecycle 195, max|Δ|=0 (CORE2 1-rank)
FER_GM=1 REDI=1 bash tools/run_lifecycle_kpp_gate_core2.sh 3  # M5b: KPP+GM+Redi 195, max|Δ|=0 (production combo)
bash tools/run_lifecycle_gm_gate_core2.sh        # M4c: GM bolus lifecycle 195 records, max|Δ|=0 (CORE2 1-rank)
bash tools/run_lifecycle_redi_gate_core2.sh      # M4d: GM+Redi lifecycle 195 records, max|Δ|=0 (CORE2 1-rank)
bash tools/run_lifecycle_fullynative_gate_core2.sh        # M3f: GM/Redi-OFF fully-native lifecycle 195 (no-regr)
bash tools/run_lifecycle_gmredi_native_gate_core2.sh 3 0  # M4e: GM+Redi fully-native lifecycle 195, max|Δ|=0 (whichEVP 0)
bash tools/run_lifecycle_gmredi_native_gate_core2.sh 3 1  # M4e: ... whichEVP 1 (mEVP)
bash tools/run_lifecycle_gmredi_native_gate_multirank.sh 2 3 0  # M4f: GM+Redi MR fully-native dist_2 195, max|Δ|=0
bash tools/run_lifecycle_gmredi_native_gate_multirank.sh 8 3 1  # M4f: ... dist_8 whichEVP 1 (needs env.sh KNEM flag, L35)
bash tools/run_icefct_gate_core2.sh              # M3c: 6 fields × whichEVP 0/1, max|Δ|=0 (CORE2 1-rank)
bash tools/run_iceflux_gate_core2.sh             # M3e: 5 fields × whichEVP 0/1, max|Δ|=0 (CORE2 1-rank)
cd build_intel_dp && ctest                       # 13/13
```
**Before assuming any kernel needs a halo extension, run the cheap `dist_N`-file invariant check first (L33).** The
optional-`partit`-absent path is byte-for-byte the proven 1-rank code, so 1-rank callers never regress.

✅ **RESOLVED — the free-surface CG "reproducibility floor" was a bug, now fixed (2026-06-21; LESSONS L29).** The L28
"iterative-solver floor" diagnosis was WRONG (a CG with byte-identical `A`/`b`/`x0`/kernels cannot manufacture
divergence from iteration count). Root cause: the preconditioner off-diagonal divide in
`oce_ssh_solve.F90::ssh_solve_preconditioner` was **auto-vectorised** (`divpd`) by the compiler while the FESOM2 oracle
compiles it SCALAR (`divsd`) — packed vs scalar division differ ~1 ULP under `-no-prec-div -fimf-use-svml`, so the
**un-gated `pr_values` array** drifted in 1299/870146 entries, seeding `z=M⁻¹r` and surfacing in the CG residual at
~iter 6 (CORE2 136 iters; pi's 37 iters stayed below the bit). The asymmetry: the oracle writes via a LOCAL POINTER
(compiler can't disprove aliasing → scalar), F3 wrote the DERIVED-TYPE COMPONENT (provably distinct → vectorised).
**Fix: one `!DIR$ NOVECTOR`** on the precond loop → scalar `divsd`, `pr_values` `max|Δ|=0`. **Verified (production,
un-instrumented binaries):** `tools/run_pressure_gate_core2.sh` PASS `max|Δ|=0` on all 28 fields incl.
`d_eta`/`eta_n`/`uv_upd`/`hbar`; and the previously-blocked `tools/run_lifecycle_gate_core2.sh` now MATCHes — **195
records (13 substeps × 5 probes × 3 steps), worst |Δ| = 0**. The whole dynamical core is byte-exact across multiple
steps on CORE2. Method + generalisable lessons (gate every consumed intermediate; a divide byte-matches only when
operands AND SIMD width match) in LESSONS L29. *Historical scoping notes + the DONE sub-gates follow.*

M2.11 was decomposed into byte-gated sub-tasks (mirroring M2.10a/b/c); scoping done **2026-06-20** (3 oracle-investigation
passes — all input paths verified: CORE2 mesh `/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2/`, IC `phc3.0_winter.nc`):

- **M2.11a — CORE2 geometry byte-gate — ✅ DONE (2026-06-20; `max|Δ|=0`, 19 fields, `tools/run_geom_gate_core2.sh`).**
  Scaled mesh-read ~40× (CORE2 `nod2D=126858`/`elem2D=244659`/`edge2D=371644`/`edge2D_in=362333`/`nl=48`) and **CLOSED
  the L8 CW-swap deferral**: FESOM3 made **244654/244659** `enforce_cw_orientation` swaps (≈100%, mesh stored CCW; pi
  had 0) and the post-swap `elem2D_nodes`/centroid/`elem_area`/`gradient_sca`/`area` are all `max|Δ|=0` vs the oracle's
  runtime `test_tri`. **`tools/make_dist1.py`** hand-crafts dist_1 (validated vs the pi template); `fesom_geomdump`
  needed NO code change (just `FESOM3_MESH_DIR`); the geom shim stops at mesh_setup (~2 s, mesh-only). See LESSONS L26.
  ↳ **ACTIVE NEXT: M2.11b.** *(historical scoping notes follow.)* **No `dist_1` existed** (only dist_2/4/…/512) →
  hand-crafted from the pi template
  (`tests/data/MESHES/pi/dist_1/`: `rpart.out`=npes/nod2D/identity-list, `my_list00000.out`=identity,
  `com_info00000.out`=empty halo w/ blank lines). Run the existing FESOM2 geom-dump shim (`fesom_geom_dump.F90`, fires
  at `mesh_setup`, STOPS before forcing → NO lifecycle/output-crash exposure) + `fesom_geomdump` on CORE2 1-rank;
  extend `tools/run_geom_gate.sh` for CORE2 → `max|Δ|=0`. **Safe deferrals confirmed:** min `nlevels=5` (nodes+elems)
  → NO single-layer columns (L18 TDMA OOB guard NOT needed); **cavity + partial-cell OFF**.
- **M2.11b — initial conditions (`do_ic3d`) — ✅ DONE (2026-06-21; `max|Δ|=0`, 3 fields [Z_3d_n/ic_temp/ic_salt], CORE2
  1-rank, `tools/run_ic_gate_core2.sh`).** Ported `oce_initial_state`→`do_ic3d` (`gen_ic3d.F90:493`) into
  `src/oce/oce_initial_state.F90` = netCDF 3D-double read (`nc_get_var3d_dp` in `mod_io_netcdf` + reused
  `forcing_binarysearch`/bilinear from `mod_forcing_read`) + vertical LINEAR interp onto `Z_3d_n` + `extrap_nod3D`
  (`gen_support.F90:400`; fill `dummy=1e10` by Gauss-Seidel neighbour-average — the partition-order step, deterministic at
  1-rank) + Kelvin guard + `insitu2pot` (`ptheta` Bryden-1973 RK4 + `atg`, in `oce_pressure_bv.F90`). **IC file CONFIRMED:
  `phc3.0_winter.nc`** (360×180×33, the work_core production file — NOT woa18; the test harness path also has it). phc3.0 land
  = NaN (no `_FillValue`) → missing-value mask = `ieee_is_nan`. `idlist=2,1` ⇒ salt→data(2) FIRST, temp→data(1) SECOND,
  `t_insitu=.true.`. Oracle shim `fesom_ic_dump.F90` dumps the live `Tclim`/`Sclim` (no prescribe). See LESSONS L27.
- **M2.11c-1 — unforced multi-step lifecycle + CORE2 dynamics-kernel gate — ✅ DONE (2026-06-21; whole dynamical core
  incl. CG `d_eta` `max|Δ|=0` over multiple steps post-L29).** Built: FESOM3 `src/drivers/fesom_lifecycle.F90` (cold-start CORE2 + `do_ic3d`
  IC + N-step `step_oce` loop, multi-step AB2 in place), oracle `tools/run_lifecycle_core2.sh` (REAL multi-step
  `oce_timestep_ale` + built-in dump_shim, `use_ice=.false.` unforced — `forcing_setup` is then a no-op, all surface
  fluxes 0), `tools/run_lifecycle_gate_core2.sh` (compares 13 NODE substeps × 5 probes × N steps), AND the NEW
  `tools/run_pressuredump_core2.sh`+`tools/run_pressure_gate_core2.sh` (the M2.1-M2.9 per-kernel gate on CORE2 —
  `fesom_pressuredump` now env-configurable `FESOM3_STEP_PER_DAY=48`). Oracle fixes (uncommitted): `io_meandata.F90::output`
  + `io_restart.F90::write_initial_conditions` both `if(partit%npes==1) return` (the 1-rank `io_gather init_nod2D_lists`
  bug hits BOTH output and restart-write); `use_sw_pene=.false.` in the run dir (else `sw_3d` unallocated → segfault, L28).
  **After the L29 fix the CORE2 pressure gate is `max|Δ|=0` on all 28 dynamical-core fields INCLUDING `d_eta`, and the
  multi-step `run_lifecycle_gate_core2.sh` MATCHes (195 records, worst |Δ|=0).** Debug `-check all` clean. See LESSONS L29.
- **M2.11c-2 — FORCED lifecycle — ✅ DONE (2026-06-21; forced dynamical core byte-exact on CORE2, CG floor resolved by L29, re-verify CONFIRMED).**
  The REAL forced FESOM2 lifecycle (`use_ice=.true.`, real CORE2 NCAR forcing at the 1948 stubs + pool runoff/SSS, the ice
  EVP + `oce_fluxes` producing the air-sea fluxes) runs cleanly at 1-rank on CORE2 — **the ice-at-1-rank unknown is
  DE-RISKED** (`tools/run_lifecycle_forced_core2.sh`). Built: oracle `port2/fesom2/src/fesom_flux_dump.F90` (per-step
  full-field dump of `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` [nod2D] + `stress_surf` [2,elem2D] BEFORE
  `oce_timestep_ale`, env `FESOM_FLUX_DUMP`, wired in `fesom_module.F90` runloop), FESOM3 `fesom_lifecycle` reads them per
  step (`FESOM3_FLUX_FILE`) and prescribes into `step_oce` (the M2.5/M2.8 prescribe-the-unsourced-input pattern), gate
  `tools/run_lifecycle_forced_gate_core2.sh`. **Result identical to unforced — now byte-exact on ALL substeps:** the
  re-verify ran 2026-06-21 POST-L29 and `tools/run_lifecycle_forced_gate_core2.sh` MATCHes — **195 records (13 substeps ×
  5 probes × 3 steps), worst |Δ| = 0** (the pre-L29 `d_eta` drift 5.5e-17 → 3.5e-6 is GONE; the NOVECTOR fix makes the
  forced CG byte-exact exactly like the unforced). FESOM3 step diagnostics match the oracle to all printed digits
  (step 1 uv=0.23283/eta=0.34914 … step 3 uv=0.45194/eta=0.67866). So **the M3-gap flux prescription is byte-exact.**
  Calendar: CORE (noleap) forcing REQUIRES `include_fleapyear=.false.`; `use_sw_pene=.false.` (matches the ported step_oce
  — no `sw_3d` term). Debug `-check all` clean (incl. the flux-read path).
- **SCOPE (LESSONS L25):** `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` air-sea budget = **M3** (the `obudget`
  in `ice_thermo_oce.F90` + `oce_fluxes`) → prescribed from the oracle dump at M2.11c. The SSS/runoff/chl climatology
  reads + the JRA55 (gregorian) forcing variant are deferred (CORE2 NCAR stubs suffice). M1's multi-rank advection gate
  + the local-mesh remap ride **M2.12**.

- **The reduced-M2 oracle namelist was NEVER needed** (through M2.9b). Every M2.x gate runs the prescribe-and-stop
  shim at end-of-`ocean_setup`, which drives the REAL kernels (up to the FULL `oce_timestep_ale` at M2.9b) on
  PRESCRIBED state — `oce_timestep_ale` reads forcing ARRAYS (set by the shim), not files, so the 1-rank forcing-file
  hang (L8) is bypassed. The shim FORCES the reduced-M2 dispatch (`mix_scheme_nmb=2`/`Fer_GM=.false.`/`Redi=.false.`/
  `opt_visc=7`/`use_wsplit=.false.`) over the shipped pi namelist (which keeps all GM/Redi/KPP arrays allocated).
  M2.10 forcing is where reading real files + running PAST the shim's stop first becomes necessary.
- **`use_wsplit=.true.` / `adv_tra_vert_impl` was an unported kernel** (do_oce_adv_tra's FCT implicit-vertical-advection
  correction; guarded with `error stop`); M2.9b forced `use_wsplit=.false.` (M1.4 precedent). **Ported 2026-10-02**
  together with the smooth Shchepetkin split and run end to end in the `wsplit` conserve configs — see
  "Implicit vertical advection: the smooth `use_wsplit` split" under "Momentum advection options".

**M1 multi-rank gate (folded into M2.12b):** the local-mesh remap (global→local
numbering/connectivity/geometry, com-structs) is ✅ DONE (M2.12a — geometry byte-gated `max|Δ|=0` on dist_2+dist_8;
the multi-rank `enforce_cw_orientation` swap path is CLOSED). **M2.12b** then lifts M1.1–M1.4 advection to multi-rank
+ adds the `find_neighbors` halo `nod_in_elem2D` dance, and **M2.12c** gates the WHOLE model on 1/8/32-rank. For
M2.12b re-confirm the M1.1–M1.4 kernels' dropped halo exchanges (tr_xy/edge_up_dn_grad/fct_LO/fct_plus_minus/
del_ttf) + loop bounds (myDim vs myDim+eDim) are correctly lifted; the gate target is the post-exchange OWNED values
on the SAME partition (L8 accumulation-order caveat). Still deferred to a richer mesh: the FCT
`AUX`/`edge_up_dn_grad`-scratch cavity caveat (L11).

## M2.12 entry notes (scoped 2026-06-21; M2.12a + M2.12b DONE — the rest carries to M2.12c)

The full decomposition (M2.12a/b/c) lives in `docs/plans/2026-06-18-fesom3-architecture.md` Task M2.12.
**M2.12a (local-mesh remap + geometry gate) AND M2.12b (multi-rank tracer advection gate) are ✅ DONE** (both
`max|Δ|=0` pi dist_2+dist_8, see "Next task"). **Active next = M2.12c** (whole-model: lift the dynamics + tracer
SOLVE, threading the optional `partit` exactly like M2.12b). The notes below are now mostly reference; the
M2.12b-flagged items (the `find_neighbors` halo dance, `exchange_elem_full`, the area/elem_area exchanges) are DONE.

- **Foundation EXISTS — do NOT rebuild it.** `src/infra/mod_halo.F90` is a real multi-rank MPI exchange
  (manual pack → Isend/Irecv → `MPI_Waitall` → broadcast-only owner→halo unpack; `exchange_nod` 2D/3D-real +
  2D-int, `exchange_elem` 2D/3D-real), and `src/infra/mod_partitioning.F90::read_dist_partition` already loads
  `myList_nod2D/elem2D/edge2D` (local→global) + `com_nod2D`/`com_elem2D`/`com_elem2D_full`. Both tested
  1/2/8-rank (`test_partit`, `test_halo`). FESOM2 uses precompiled MPI_TYPE_INDEXED datatypes vs FESOM3's
  manual pack — same bytes moved, no arithmetic → byte-identical. **The geometry pipeline
  (`compute_geometry`/`setup_vertical`/`enforce_cw_orientation`) is already partition-agnostic** (loops over
  `mesh%nod2D`/`elem2D`/`edge2D`, reads only local arrays) → it runs UNCHANGED once the local arrays are built.
- **Gate rule (L8) — the load-bearing constraint.** Compare FESOM3 `dist_N` vs FESOM2 `dist_N` PER-RANK on
  OWNED entries (1..myDim). NOT vs 1-rank global: both `dist_N` runs share the same `myList` order → local idx
  i ↔ same global id AND the per-node area sum order matches → byte-identical; the 1-rank global uses a
  different element permutation (non-associative FP) → would mismatch at ~ULP.
- **Algorithm (FESOM2 `port2/fesom2/src/oce_mesh.F90`):** ⟶DONE in M2.12a (`read_mesh_local`): `read_mesh` :212
  (scatter — implemented with a full-size inverse map, not FESOM2's chunked `mapping`; identical result); edges read
  :1787 (scatter + negative-localize); owned `nod_in_elem2D` build (`find_neighbors` :2021-2049, OWNED part);
  `setup_vertical` `exchange_nod(nlevels_nod2D_min)` :1669. `elem2D_nodes` is `(MAX_NV,myDim_elem2D)` — OWNED only.
  ⟶DONE in M2.12b: the `find_neighbors` HALO dance — `exchange_nod(num)` :2051, pack global elem-ids :2057
  → `exchange_nod` → re-localize through **eXDim** :2066-2073 (`complete_nod_in_elem_halo` in `read_mesh_local`,
  completes halo `nod_in_elem2D` for MUSCL/advection); and `mesh_areas` :2162 — `compute_geometry` now does
  `exchange_nod(area/areasvol/*_inv)` :2322 + `exchange_elem_full(elem_area)` :2220 at npes>1 (the owned-node area
  accumulation stays owned-only, the exchange fills the halo). `elem_neighbors`/`elem_edges` still SKIPPED (not
  needed for advection; may be needed by the M2.12c dynamics — re-check when lifting `momentum_adv`/visc).
- **Infra note (M2.12a + M2.12b):** the owned-entry GEOMETRY gate (M2.12a) did NOT need `exchange_elem_full` — the
  fix there was precompute owned element centers + `exchange_elem` them. **M2.12b ADDED `exchange_elem_full`**
  (`com_elem2D_full`, eDim+eXDim; 2D-r/2D-i/3D-r variants + `core_blk_r`) — needed for `elem_area` AND `tr_xy` AND
  the `coord_elem`/`e_nodes` of `find_up_downwind_triangles` (all read at eXDim through a halo node's element list).
- **Oracle (⟶DONE):** `port2/fesom2/src/fesom_geom_dump.F90` extended to npes>1 (per-rank LOCAL OWNED dump,
  `<path>.<mype5>`); `libfesom.so` rebuilt. **Gate (⟶DONE):** `tools/run_geom_gate_multirank.sh [np]`
  (FESOM2 dist_N dump vs FESOM3 dist_N dump, per-rank owned, `geom_diff.py`); `max|Δ|=0` on dist_2 + dist_8.
  CLOSED the deferred multi-rank `enforce_cw_orientation` swap path. Lessons in **L30**.
- **Build discipline:** clean-rebuild after adding NEW files (L19); `ulimit -s unlimited` for the dump writer (L20).

## M3 entry notes (scoped 2026-06-22)

The full M3 decomposition (M3a–M3f) is in "Next task". M3 is sea ice (standard EVP) → the air-sea budget. The
FESOM2 source is the same v2.7.3 oracle: `ice_EVP.F90` / `ice_fct.F90` / `ice_thermo_oce.F90` / `ice_oce_coupling.F90`
/ `ice_setup_step.F90` / `MOD_ICE.F90`. **M3a–M3e ✅ DONE (all `max|Δ|=0`); active next = M3f** (assemble the native
ice step + `oce_fluxes` into a forced lifecycle, drop the prescribed fluxes, then multi-rank; tag `m3`).

- **Call sequence (FESOM2 `fesom_module.F90` runloop):** `ocean2ice` (:677, set ice's view of the ocean surface) →
  `ice_timestep` (:707 = `ice_setup_step.F90:96`: `EVPdynamics`(:215) → `ice_TG_rhs`(:258) → `ice_fct_solve`(:261) →
  `cut_off`(:295) → `thermodynamics`(:314)) → `oce_fluxes_mom`(:715) → `oce_fluxes`(:716). Ice runs at `ice_dt =
  ice_ave_steps·dt` (=`dt` here; `Tevp_inv=3/ice_dt`); `evp_rheol_steps=120` subcycles; the σ stress tensor is
  PROGNOSTIC (elastic memory persists across steps — already a `t_ice` field in `mod_ice.F90`).
- **Oracle = the proven forced lifecycle.** `tools/run_lifecycle_forced_core2.sh` ALREADY runs the REAL FESOM2 ice
  EVP+advection+thermo+oce_fluxes at 1-rank on CORE2 (`use_ice=.true.`, real NCAR forcing) and dumps the air-sea fluxes
  via `fesom_flux_dump.F90`. M3 per-kernel gates extend this: add a prescribe-and-stop `fesom_ice_dump` shim (mirror
  `fesom_step_dump`) that runs `ice_setup` + one `ice_timestep` + `oce_fluxes` on the do_ic3d-IC'd ocean and dumps each
  intermediate (`mod_dump` byte format), per-kernel gated like M2.1–M2.9. The flux dump is the M3e/M3f target.
- **Reuses (don't rebuild):** `ice_mass_matrix_fill` rides the **M2.6 `ssh_stiff` CSR** (rowptr/nn_num/nn_pos); the
  ice FCT advection reuses the M1/M2.12b tracer-FCT pattern; the cold-start IC reads `SST` from the M2.11b `do_ic3d`
  state + `geo_coord_nod2D`; dumps reuse `mod_dump`/`mod_binary_arrays`. `t_ice` skeleton already in
  `src/types/mod_ice.F90` (sigma11/12/22 prognostic, a_ice/m_ice/m_snow, uice/vice, work) — flesh it out from `MOD_ICE.F90`.
- **Config (reduced, like reduced-M2):** standard EVP only — `whichEVP=0`; NO `__icepack`/`use_meltponds`/`use_cavity`/
  `__oasis`/`__yac`/`__oifs`; `num_itracers=3` (a_ice/m_ice/m_snow). Namelists `&ice_dyn` + `&ice_therm` (defaults in
  `MOD_ICE.F90`: `pstar=30000`, `ellipse=2`, `c_pressure=20`, `delta_min=1e-11`, `evp_rheol_steps=120`,
  `ice_gamma_fct=0.25`, `ice_diff=10`, `cd_oce_ice=5.5e-3`, `Sice=4`, `h0=0.5`, densities rhoice=910/rhosno=290/...).
- **Gate rule (same L8):** per-rank owned dumps vs same-partition FESOM2; ice fields are nod2D (a/m/uv/flx) or elem2D
  (sigma). Watch (plan §M3): `bc_index_nod2D` multi-rank safety, NaN-vs-0 ice masking in diagnostics, ice_dt sync.

## Open notes / risks

- **Stale MPI payloads on the login node = an Open MPI 4.1.2 bug, NOT a halo-code bug (2026-10-08,
  hybrid plan Task 10 follow-up, LESSONS L61; supersedes the "corrupts >131 KB" reading of L35).**
  `pml ob1` + `btl vader` + knem single-copy caches knem registrations by virtual address and misses
  glibc heap trims (`brk`); a buffer re-allocated at the same address is then read from the OLD pinned
  pages (right count, correct first 32 KB, payload of an earlier message). Standalone proof:
  `tools/run_mpi_knem_stale_repro.sh` (stale only under knem + heap trim; clean with single-copy
  `none`/`emulated`, `mpi_leave_pinned 0`, or `pml ucx`). `mod_halo` was audited: contiguous buffers
  alive until `MPI_Waitall`, exact request counts, no overlapping exchanges, and an
  `-check arg_temp_created` build creates no argument temporary in `test_vinv` or the production
  driver (pi np 2, CORE2 np 8). **Production is not affected**: `work/job_levante` sources `env.sh`
  → `shell.intel` sets `OMPI_MCA_pml=ucx` (vader unused; UCX tested clean with its knem/cma
  transports), and `env.sh` has also exported `OMPI_MCA_btl_vader_single_copy_mechanism=none` since
  M2.12. Hence this bug **cannot explain the dotty `w` near topography** in the CORE2 batch runs.
  Exposed are only login-node `mpiexec --mca btl self,vader` runs WITHOUT `env.sh` and without the
  flag (ctest gets it from `FESOM3_MPI_LAUNCH_FLAGS`; `run_tri_regression.sh`, `run_conserve_pi.sh`,
  `run_bottom_audit.sh` pass it; the other `tools/run_*.sh` source `env.sh`).
- The FESOM2 oracle is built + PROVEN (M0.7 geometry + M1.1–M1.4 advection + M2.1–M2.9a
  pressure/EOS/N²/PGF/vel_rhs/momadv/viscosity/ivertvisc/SSH/ALE-update/PP-mixing/convective-adjustment/tracer-solve gates
  all `max|Δ|=0`, **57 fields**, `tools/run_pressure_gate.sh`) + M2.9b the WHOLE assembled `oce_timestep_ale`
  (`tools/run_step_gate.sh`, **65 records** vs the REAL step's built-in dumps, `--ignore-substep=2`).
  Self-tests run standalone (13/13 ctest). M2+ kernel gates extend the proven 1-rank end-of-`ocean_setup` shim
  pattern (`fesom_pressure_dump.F90` drives the isolated kernels; `fesom_step_dump.F90` drives the REAL
  `compute_vel_nodes` + `oce_timestep_ale` on prescribed state). The prescribe-and-stop shim STILL sufficed for
  M2.9b — the assembled step reads forcing ARRAYS (set by the shim), not files, so the reduced-M2 namelist + the
  forcing-file read were never needed. **M2.10 forcing is next** (the FIRST netCDF I/O), see "Next task".
- **⚠️ CLEAN-rebuild rule (L19):** when a milestone adds NEW `src/**/*.F90` files (not just edits), run
  `./configure.sh --compiler intel --precision dp --clean --build` before trusting the byte-gate — an
  incremental build re-runs CMake configure (GLOB `CONFIGURE_DEPENDS`) and can link objects against a MIX
  of stale/fresh `.mod` interfaces, ULP-drifting EVERY field (the M2.6 first-gate red herring).
- **M2.12 (multi-rank, MVP exit) is ✅ DONE** (a/b/c-1/c-2/c-3). The local-mesh remap + the whole-model multi-rank
  byte-match (all 65 substeps, dist_2+dist_8) landed cleanly. The deferred cavity caveat (L11) + dist_32 (pi has no
  such partition) ride a richer mesh; M3 is next.
- **M2.12c-1/c-2/c-3 uncommitted artifacts:** FESOM3 — drivers `src/drivers/fesom_stepdump_mr.F90` (c-1/c-2: dynamics
  chain → d_eta, direct kernel calls) + **`src/drivers/fesom_stepfull_mr.F90` (c-3 NEW: the WHOLE step via `step_oce`
  through `partit`)**; gates `tools/run_stepdyn_gate_multirank.sh` (c-1/c-2) + **`tools/run_step_gate_multirank.sh`
  (c-3 NEW: all 65 substeps, ignore only SW_AB id=2)**; `mod_halo.F90` (c-1 rank-3 `exchange_nod`; c-2 `allreduce_sum`),
  `mod_part_bounds.F90` (`local_dims`); the c-1 lifted dynamics kernels (`oce_pressure_bv`/`oce_pgf`/`oce_ale_mixing_pp`/
  `oce_mo_conv`/`oce_dyn_velrhs`/`oce_dyn_visc`/`oce_dyn_ivertvisc`), c-2 `oce_ssh_rhs.F90`+`oce_ssh_solve.F90`, and
  **c-3 `oce_ale.F90` (update_vel/compute_hbar_ale/update_eta_n/vert_vel_ale/compute_CFLz/compute_Wvel_split/
  update_thickness_ale optional partit) + `oce_ale_tracer.F90` (solve_tracers_ale/diff_tracers_ale/diff_part_hor_redi/
  diff_ver_part_impl_ale optional partit) + `mod_step_oce.F90` (optional partit threaded to every kernel)**. FESOM2
  oracle — **UNCHANGED for c-3** (the c-1 `fesom_step_dump.F90` npes>1 extension — dropped the `npes/=1` return; nodes
  owned+halo, element velocity OWNED + `exchange_elem(UV)`; ⚠️ do NOT prescribe UV at eXDim elements,
  `elem2D_nodes(1)` can be beyond the eDim node halo → `coord_nod2D` OOB — already drives the whole `oce_timestep_ale`
  and dumps ALL substeps). Rebuild the oracle with `make -C port2/fesom2/build fesom.x` only if its tree changed.
- The FESOM2 oracle shim edits (`fesom_advhor_dump.F90`, `fesom_pressure_dump.F90`, `fesom_step_dump.F90` [now npes>1],
  `fesom_forcing_dump.F90`, `fesom_ic_dump.F90` [M2.11b], `oce_setup_step.F90`, `fesom_module.F90`, geom/ale shims) **+
  the M2.10a `io_netcdf_workaround_module.F90` np=1 fix** live as UNCOMMITTED working-tree instrumentation in
  `port2/fesom2` (not committed there, by design); `libfesom.so` must be rebuilt (`make -C port2/fesom2/build
  fesom.x`, then re-run `cmake .` first if a NEW shim file was added so the GLOB picks it up) after editing a
  shim. The `next_io_rank` np=1 fix is behavior-preserving (sequential I/O, value-identical) and arguably a real
  bug fix (the code's own TODO admits the recursion never ends at 1 rank).

## MUSCL horizontal tracer advection: on-the-fly up/downwind reconstruction (reference RK3 rule)

FESOM2 stores the per-edge up/downwind tracer gradient `edge_up_dn_grad(4, nl-1, nEdge)`,
rebuilt every step by `fill_up_dn_grad` (triangle gradient on the levels shared by both edge
nodes, a Miura node-averaged gradient on the others and at edges without both triangles).
FESOM3 follows the reference RK3 code (`qq/oce_stepRK3.F90`, `t_hor_adv_muscl_RK3`) instead:

- `init_tracers_AB` builds `twork%tr_xy` (elemental gradient of `values`, full element halo),
  persistent scratch; `solve_tracers_ale` reuses it for the horizontal diffusion / Redi.
- `adv_tra_hor_muscl(…, tr_xy, edge_up_dn_tri, …)`: side 1 takes `tr_xy` of
  the upwind triangle, side 2 that of the downwind triangle, on the levels where that triangle
  is wet; a missing or dry triangle gives a **zero** increment. No stored array, no node
  average. Kept from FESOM3 (not from the reference): `edge_dxdy` with the folded metric, edge
  flux output for FCT, cavities (`ulevels`), the gradient of `values` (not of the advected `valuesAB`, as FESOM2).
- **Not bit-identical to FESOM2 any more**: identical on every (edge, level) where both edge
  triangles exist and the level is shared by both edge nodes; different at coasts and
  bathymetry / cavity steps (core2: 0.8 M of 10.9 M edge-levels, max |Δflux| 0.4 % of max
  |flux|; pi 1.3-1.5 %).
- `fill_up_dn_grad` is removed (2026-10-08) together with the FESOM2-comparison parts S/T of
  `test_muscl_onthefly` and the `edge_up_dn_grad` dump record; the numbers below are from
  before the removal.

Net: `test_muscl_onthefly` (np 1/2; pi as read and with synthesised cavities): part R = the
production kernels equal the stored-array kernels fed with the reference rule, bitwise, both
velocity signs; part S = equal to FESOM2 bitwise on the shared levels, different elsewhere
(teeth); part T = timing. core2 (login node, one MFCT call incl. the gradient preparation):
np 4 42 vs 74 ms, np 1 180 vs 316 ms (~43 % faster than FESOM2's fill + kernel); memory
-535 MiB in total (~1 MiB/rank at 512). Run the test on core2 with `ulimit -s unlimited`.

**Single MUSCL kernel (2026-10-08)**: as in the reference RK3 code there is ONE MUSCL kernel,
`adv_tra_hor_muscl` = FESOM2's `adv_tra_hor_mfct` (no boundary clamp; where the up/downwind
triangle is missing or dry the face value is ⅔·T_up + ⅓·T_dn). FESOM2's clamped MUSCL and
everything of `nboundary_lay` (init, `t_tracer_work` field, dump record) are removed;
`tra_adv_hor` is `UPW1` or `MUSCL` (`'MFCT'` is no longer accepted). The production FCT
path is unchanged bitwise (only the name). Without FCT the scheme runs but is not monotone.

## Momentum advection options (`momadv_opt`)

| value | scheme | where |
|---|---|---|
| 1 | vector invariant, `(curl u + f) x u + grad(u^2/2) + w du/dz` | `src/oce/oce_dyn_vinv.F90` |
| 2 | flux/scalar form (default, all drivers pin it) | `momentum_adv_scalar`, `src/oce/oce_dyn_velrhs.F90` |

Both only ADD into `UV_rhsAB`, dispatched as siblings in `compute_vel_rhs`; `momadv_opt==2`
is unreachable from the `==1` branch, so enabling the new scheme cannot perturb the old one.

`FESOM3_MOMADV_OPT` (**value**-based, not presence-based — `momadv_opt` is 1 or 2; default
2) is read by `fesom_conserve` and `fesom_lifecycle_native_mr`; `tools/run_conserve_pi.sh`
gates it at np=1 and np=8. The other production drivers pin `momadv_opt = 2`.

`rvo_upwind` (`FESOM3_RVO_UPWIND`, value in [0,1], default 0) blends the face
relative-vorticity reconstruction toward the 1-2 upwind edge-neighbours, weighted by the
inflow `d+|d|` through each face (inward normals rebuilt from the element's OWN vertices --
per-edge geometry is owned-edges-only). 0 = plain vertex average, structurally bit-identical;
1 = fully upwind. Convex (max principle), energy-neutral (`u.[(f+zeta) x u] == 0` for any
zeta), linear in the coefficient. Net: `test_vinv` V8-V12 + the `momadv-vinv-upw` gate
configs at np=1/8. `mesh%elem_neighbors` is built once in `build_elem_adjacency`
(`compute_geometry`); `elem_edges` remains unbuilt.

Revalidated after the fact (critical re-review of `6e38cf6`): the one open question was
the halo ring — `exchange_elem` fills only eDim while `elem_neighbors` may hold halo ids.
Settled, not patched: edge-neighbours of an owned element are in eDim by construction
(FESOM2 `gen_comm.F90`: `com_elem2D` = edge-sharing, `com_elem2D_full` = node-sharing),
verified at np=2/8 and guarded by `test_vinv` A5 (production fill/exchange/read pattern
with a sentinel). `exchange_elem_full` is NOT needed here. See L56.

**There is no FESOM2 oracle for `momadv_opt==1`** (the v2.7.3 branch aborts: "not adapted
mom_adv advection typ for ALE"). Its net is `test/test_vinv.F90` — analytic cases that are
exact to round-off — plus the conservation gate. See L54 and
`docs/plans/completed/2026-10-01-momadv-vector-invariant.md`.

### Implicit vertical advection: the smooth `use_wsplit` split

`w = w_e + w_i` face by face in `compute_Wvel_split` (end of `vert_vel_ale`, owned+halo,
`src/oce/oce_ale.F90`): `w_i = f·w`, `w_e = w − w_i`, with the implicit share `f` from
`oce_wsplit::wsplit_implicit_fraction(CFL_z, wsplit_mincfl, wsplit_maxcfl)`, the C¹
limiting function of Shchepetkin (2015, Ocean Modelling 91, Sec. 3.1 / Fig. 9; NEMO
`wAimp`). With `D = Cu_max − Cu_min`, `F = 4·Cu_max·D`, `Cu_cut = 2·Cu_max − Cu_min`:
`f = 0` for `Cu ≤ Cu_min`; `f = x²/(F + x²)`, `x = Cu − Cu_min`, on the bend;
`f = (Cu − Cu_max)/Cu` for `Cu ≥ Cu_cut`, so the explicit Courant number `Cu·(1−f)`
saturates EXACTLY at `wsplit_maxcfl` (FESOM2's cap, reached without the kink — the
parameter keeps its FESOM2 meaning). `f` and `f'` are continuous at both joints, `Cu_e` is
non-decreasing, `f → 1`. FESOM2's hard switch (`w_e = w·C/CFL_z` above `C`) is the
degenerate `Cu_min = Cu_max` limit (`F = 0`): admissible to `wsplit_check_params`, used
only as the positive control of the smoothness test, NOT a mode (decision 2026-10-02).
The argument is FESOM2's `CFL_z` (`compute_CFLz`, unchanged):
`|w(nz)|·dt·(1/h(nz−1) + 1/h(nz))`, the interface flux counted against BOTH adjacent
cells (one term at the surface), so `wsplit_maxcfl = 1` means `|w|·dt/h ≈ 0.5` per cell
on uniform layers.

Parameters: `t_dyn%wsplit_mincfl` (type default 0.5 — an ABSOLUTE number, `0.5·maxcfl`
only at the default cap) and `%wsplit_maxcfl` (1.0), both serialized (`test_types` round
trip); admissible iff `maxcfl > 0` and `0 ≤ mincfl ≤ maxcfl` (`wsplit_check_params`).
Validation is SETUP-time and the kernels are stateless (house convention): every driver
that enables `use_wsplit` calls `wsplit_assert_params(mincfl, maxcfl)` once — the
obligation is documented in `t_dyn` next to the parameters; `fesom_conserve`,
`fesom_lifecycle_native_mr` and `fesom_pressuredump` do — so a driver that lowers the cap
without the onset cannot run the discontinuous pair in silence; the ≥ 2-layer column
precondition of the implicit solves is a mesh invariant asserted once, collectively, in
`read_mesh` (`assert_min_layers`; ctest `test_minlayers_np1` trips it). Env hooks in
`fesom_conserve` and `fesom_lifecycle_native_mr`:
`FESOM3_WSPLIT` (presence → `use_wsplit`), `FESOM3_WSPLIT_MAXCFL` (value, default 1.0),
`FESOM3_WSPLIT_MINCFL` (value, default `0.5·maxcfl`, so a config that sets only the cap
keeps the type's onset/cap ratio; an unparsable value is an `error stop`, not the
default), `FESOM3_WSPLIT_EXPECT_SPLIT=N` (value; `fesom_conserve` only: `error stop`
unless ≥ N owned faces split at the last step and ≥ 1 has `f ≥ 0.5`). Off by default
everywhere: `fesom_lifecycle`, `_mr`, `_native`, `fesom_stepdump*`, `fesom_stepfull_mr`,
`fesom_advhordump_mr` pin `use_wsplit=.false.`; `fesom_conserve` and
`fesom_lifecycle_native_mr` default to `.false.` unless `FESOM3_WSPLIT` is set. With it
off: `w_e == w`, `w_i == +0` bitwise, the previous off path verbatim, every pre-existing
gate drift row unchanged. Switching it on in `job_levante` and choosing `wsplit_mincfl`
is the user's call.

Where `w_i` acts: (1) momentum — the TDMA `impl_vert_visc_ale` in upwind flux form, with
the advective-form diagonal correction for `momadv_opt==1` (L57; `test_ivertvisc` I1-I5);
momentum advection (`momentum_adv_scalar`, `momentum_adv_vinv`) reads `w_e`. (2) non-FCT
tracers — `do_wimpl` in the vertical-diffusion TDMA `diff_ver_part_impl_ale`. (3) FCT
tracers — FESOM2's sequence in `do_oce_adv_tra` (`oce_adv_tra_driver.F90`): low-order
upwind step with `w_e`, then `adv_tra_vert_impl(dt, w_i, fct_LO)` (per-column upwind
backward-Euler TDMA on `fct_LO`, column sums `= hnode_new`, f3 1-D area rule; ≥ 2 layers
is a stated precondition, asserted at setup, not in the kernel; `oce_adv_tra_ver.F90`),
then the LO upwind flux RECOMPUTED with the full
`w` so the antidiffusive flux is `HO(w) − LO(w)`. (`HO(w) − LO(w_e)` would double-count
the `w_i` transport and still telescope — a conservation gate cannot see it.)

Nets: `test/test_wsplit.F90` W (function: cap error ≤ 4.4e-17, joints exact, C¹ to
1.6e-6, `Cu_e` monotone to 1.23 eps·Cu), S (split identities on pi: `w_e + w_i = w` to
1 ulp, `w_i = f·w` bitwise, every `CFL_z` class populated), S5 (the HALO contract: inputs
from the global node index, then `cfl_z/w_e/w_i` at every halo node == the owner's
bitwise against an `exchange_nod`'d copy — 1429 faces at np 2, 0 mismatches; positive
control with local-index inputs 1429/1429 mismatch. No consumer exchanges `w_e/w_i` and
the conserve gate is blind to a wrong halo value, so this is the only net on it), X
(transition smoothness of the non-FCT chain QR4C + `do_wimpl`, `Kv = 0`, quadratic `T`,
601-point `Cu` scan over `[0, 3·maxcfl]` on one representative column per
`(ulevels, nlevels)` class: second-difference ratio `max|D2|/(ΔCu·max|D1|)` = 0.96 / 1.17
for `w ≷ 0` against the bound 5 (a C⁰ linear ramp would pass 10 for `w < 0`); positive
control = the degenerate `(1,1)` function: 100.3 / 64.1, its kink `D2 = 1.0025·J·ΔCu` at
the closed-form jump `J = cq·h²·(1+Cu_c)/2`). Both tests synthesise cavity columns
(`ulevels = 3` on every 10th global node) so the `nzmin`-based layouts are not tested at
`nzmin = 1` only. `test/test_wimpl_tra.F90` C1-C6 (identity bitwise; constancy of uniform
`T` through explicit(`w_e`) + implicit(`w_i`) with `hnode_new` advanced by the full `w`
≤ 5.3e-16 for three parameter sets, both signs, open surface, incl. a 2-layer column;
conservation ≤ 1e-15; boundedness at `CFL_z` 10.8 / 14.0 against an uncapped positive
control overshooting by 0.95 / 5.6 of the range; the fully implicit closed form; C6 the
driver `do_oce_adv_tra` BITWISE equal to the test's own assembly at np 1/2 for BOTH signs
of `w` — the assembly carries `oce_tra_adv_fct`, inactive for `w > 0` (every factor 1
on 95469 nonzero faces, closed form of the whole step to 1.4e-14 K) and active for
`w < 0` — the wrong `HO(w) − LO(w_e)` differing by 29.9 / 58.8 K·m = the predicted
double-counted `w_i` transport to 5e-14 for `w > 0`, by 3.9 / 1.4 K·m with the limiter
for `w < 0`); ctest `test_minlayers_np1` trips the setup-time < 2-layer refusal.
End to end: `fesom_conserve` prints for every config one last-step owned-face statistics
line (max `CFL_z`; faces with `w_i ≠ 0` = `CFL_z > mincfl`; capped faces
`CFL_z ≥ Cu_cut = 2·maxcfl − mincfl`; `f ≥ 0.5` = `CFL_z ≥ 2·maxcfl`; max `f`). The
20-step pi cold start has max `CFL_z` = 2.76e-2, so the default cap 1.0 splits NOTHING;
the gate configs `wsplit` / `momadv-vinv-wsplit` / `GM+wsplit` use
`FESOM3_WSPLIT_MAXCFL=0.005` (onset 0.0025, chosen over 0.01 which left 26 faces at
`f ≥ 0.5`): 8632 faces split, 1130 capped, 500 at `f ≥ 0.5`, max `f` 0.818 (np 8: 8628 /
496; vinv: 8647 / 497), drift heat/salt/vol ~1e-14 (`tools/run_conserve_pi.sh` np 1/8 +
vinv np 1 + GM np 1; ctest `fesom_conserve_zstar_{wsplit,vinv_wsplit}_np{1,2}` with
`EXPECT_SPLIT=1000`). The guards are ctests too: `fesom_conserve_zstar_wsplit_vacuous_np1`
(cap 0.5 → "vacuous config", no drift report; also pins the driver default
`mincfl = 0.5·maxcfl`), `fesom_conserve_wsplit_badparams_np1` ((0.9, 0.5) rejected),
`fesom_conserve_wsplit_badvalue_np1` (`MAXCFL=0.5x` rejected). The conserve tracers are
FCT, so the gate runs path (3) and the momentum TDMA; path (2) is pinned at unit level by
X3 (the chain in the fully implicit limit == the closed-form upwind row recursion).

Retired byte-gate fields: `fesom_pressuredump` still forces `use_wsplit=.true.` (cap 1.0)
so `compute_Wvel_split` runs there (its diagnostic counts faces with `w_i ≠ 0`), but it
no longer dumps `w_split_e/w_split_i` (nothing read them); `tools/pressure_diff.py` keeps
the FESOM2-side names in `SKIP_FIELDS` (printed as SKIP, never compared against FESOM2's
hard split; the FESOM2 byte-gates are retired since bottom-at-vertices anyway). That
driver calls `diff_tracers_ale` only, never `do_oce_adv_tra`, so it does not reach
`adv_tra_vert_impl`. Plan: `docs/plans/completed/2026-10-02-wsplit-smooth.md`; see L57,
L58.


