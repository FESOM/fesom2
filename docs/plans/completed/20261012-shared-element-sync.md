# Shared-element velocity sync: one owner for uv

## Overview
- **Problem.** The model uses FESOM2's partition rule: a rank owns every element that has at least one vertex it owns. Elements on a partition border therefore belong to two or three ranks. dist_512 on CORE2 has 24,798 such shared elements out of 244,659.
  - Each owner updates its own copy of the element velocity `dyn%uv` in `update_vel`. Nothing makes the copies equal again: `exchange_elem(UV)` refreshes only the eDim halo (`src/oce/oce_ale.F90:158`).
  - When the copies differ, the two ranks compute different fluxes across their common edges. The barotropic divergence `ssh_rhs_old` (`src/oce/oce_ale.F90:199-229`) then no longer sums to zero over the globe.
  - With the linear free surface this creates or destroys water. The global-mean SSH drifts, which it must not, and the salt trend becomes +1.4e-5 to +2.3e-5 psu/yr instead of about −2e-6.
- **Why the copies differ** (code reading, 12 Oct 2026; confirmed by a plan review):
  - The viscosity filter adds each edge's contribution to both neighbouring elements in the rank's own local edge order (`src/oce/oce_dyn_visc.F90`, loops over `nEdgeL`). The owners sum in different orders, so their `UV_rhs`, and after `update_vel` their `uv`, differ by roundoff every step.
  - The upwind blend of the vector-invariant scheme (`src/oce/oce_dyn_vinv.F90:302-314`) switches when the inflow `d = u·n` through an edge changes sign. A roundoff difference in `uv` can therefore become an O(1) difference in `uv_rhsAB`, and then in `uv`. results5 (1 Oct, before `rvo_upwind`) is clean; results7 (2 Oct, with it) leaks.
  - Everything that writes `uv_rhsAB` works per element on nodal fields that each node's owner computes and exchanges:
    - the Coriolis term (`src/oce/oce_dyn_velrhs.F90:119-120`);
    - vorticity, KE and `omega_e` (`oce_dyn_vinv.F90:149, 201, 236`);
    - `w_e`;
    - the nodal tendency of `momentum_adv_scalar` (`oce_dyn_velrhs.F90:319`).
  - The EVP and mEVP stresses `sigma11/12/22` (`src/ice/mod_ice_dyn.F90:165-201, 607-629`) are also per element. They read `u_ice`, exchanged every subcycle (`:439-442, 677-680`), and `a_ice`/`m_ice`, exchanged after the FCT.
  - So once `uv` is equal on all owners, the copies of `uv_rhsAB` and sigma should stay bitwise equal without their own sync.
  - The 9 Oct restart gate supports this: at np 2 (boundary case, `gate_bisect/gate_HEADzarr_full.log`) only u, v, urhs_AB, vrhs_AB, temp, w and w_expl differ, while sigma is exact.
- **Evidence of the leak** (CORE2 checkpoints, 9–12 Oct 2026):
  - results7, results11, results17 from July 1965, and results21 leak up to 2.4–10 Sv.
  - In those checkpoints the stored u, v disagree with the stored `ssh_rhs_old` at 86 to 107 nodes. All of them sit on dist_512 borders and touch shared elements. Clean checkpoints agree everywhere to 4e-11.
  - A restart re-syncs the copies, because every owner reads the same value from file. results17's mean SSH stays flat for about six months after each restart.
- **Fix.** Give every shared element one designated owner: the highest rank among the owners of its vertices. Right after the owned update in `update_vel`, every other owner takes the designated owner's `uv`. This runs every step and does what a restart does today.
  - `compute_hbar_ale`, `vert_vel_ale`, the tracer advection and the next step's `compute_vel_rhs` all read `uv` after that point. The border fluxes then agree on both ranks, so the leak is gone whatever the trigger.
  - Only `uv` is synced. Guard tests check that the `uv_rhsAB` and sigma copies stay bitwise equal. If a guard fails, that field gets the same sync (fallback in Task 4).
  - The designated owner holds the copy the restart and output writers already keep. np>1 restarts therefore become exact, and the restart gate's np>1 floor of 1e-12 drops to 0.
  - np 1 has no shared elements, so its results stay bit-identical.
- **qq prototype** (`/home/a/a270029/qq/oce_comm.F90:340-560, 1040-1110`): it uses the same ownership rule and never syncs the shared copies either. Its step driver has commented-out `MPI_AllREDUCE` checks of `sum(Wvel(1,:))` (`oce_stepRK3.F90:35-89`), the same mass check this plan adds as a diagnostic.

## Context (from discovery)
- **Partition data:** `t_partit` and `com_struct` are in `src/types/mod_partit.F90:12-60`, with a fixed `MAX_NEIGHBOR_PARTITIONS = 32`.
  - `read_dist_partition` (`src/infra/mod_partitioning.F90:136-232`) reads `com_nod2D`. Its receive list is sized eDim (`:216, 249`): it covers every eDim halo node with local indices, grouped by owner rank. Every rank therefore knows the owner of all its owned and eDim nodes without any communication.
  - `partit%part` holds cumulative node counts per rank, not an owner map.
  - `set_empty_comms` (`:111-128`) clears the lists at 1 rank.
  - `die` (`:272-277`) prints a message and calls `MPI_Abort`. It is private.
  - `write_partit` and `read_partit` are used only by the round trip in `test/test_types.F90:90-129`.
- **Local mesh:** `read_mesh_local` (`src/mesh/mod_mesh_read.F90:413-571`).
  - `read_mesh` (`:28`) calls it only when npes > 1 (`:40-45`).
  - All 42 `read_mesh` calls pass program-level variables, and `test_util::setup_test_mesh` takes `partit` as `inout` (`test/test_util.F90:35`).
  - It allocates `elem2D_nodes` and `elem2D_nnodes` for the owned elements only, in local node numbers (`:500-516`).
  - `partit` is `intent(in)` in both routines today.
- **Halo layer** (`src/infra/mod_halo.F90`):
  - The broadcast-only cores `core_2D_r`, `core_2D_i`, `core_3D_r` and `core_blk_r` (`:239-381`) return at once on an empty list.
  - The exchange-and-compare `probe` (`:396-410`, behind `stale_halo_max_elem`) is the pattern for the copy-spread diagnostic.
  - `allreduce_max` is real only (`:42-44`).
- **Where `uv` is written:** only in `update_vel` (`src/oce/oce_ale.F90:114-159`: owned loop `:146-157`, `exchange_elem(UV)` `:158`).
  - The GM bolus add and subtract (`src/oce/oce_ale_tracer.F90:142, 172`) uses the exchanged `fer_uv` (`oce_fer_gm.F90:318`), so it keeps the copies equal.
  - Step order (`src/step/mod_step_oce.F90:212-286`): `compute_vel_rhs` → `viscosity_filter` → `impl_vert_visc_ale` → `compute_ssh_rhs_ale` → `solve_ssh_ale` → `update_vel` → `compute_hbar_ale` → `vert_vel_ale` → `solve_tracers_ale`.
  - `update_vel` runs inside the timer region `TMR_OCE_SSH` (`:228-244`).
- **`uv_rhsAB`:** shape `(AB_order-1, 2, nl-1, elem)`, i.e. `(1, 2, nl-1, elem)` for AB2. It stores one level only (`src/drivers/fesom_lifecycle_native_mr.F90:360`).
  - It is written only in `compute_vel_rhs` and in the momentum advection it calls: momentum-advection call at `oce_dyn_velrhs.F90:132-136`, AB blend at `:140-148`.
- **Which copy the writer keeps (verified):** `decomp_redistribute_2d_r` (`src/io/mod_io_decomp.F90:220-237`) writes the received values in rank order without dedup, so the last write wins.
  - The **highest rank's copy** therefore ends up in restart and output files.
  - A restart read gives every owner that copy.
- **Knobs, diagnostics, timers:**
  - The production knobs are read in the env block of `src/drivers/fesom_lifecycle_native_mr.F90` (around `:200-300`). README calls that block the authoritative list.
  - The SSH solver knobs are read by the solver itself and listed in README's table.
  - The `FESOM3_STEP_DIAG` block (`:1023-1039`) does one `MPI_MAX` of 6 values and prints one line on rank 0. `tools/run_cutout_short.sh:131` reads fields 5, 7 and 9 of that line.
  - `src/infra/mod_timer.F90` keeps a fixed registry of named timers, with subset timers such as "(of which) ssh solve".
- **Tests:**
  - The ctests are registered in `test/CMakeLists.txt` (`add_fesom_test`, `add_mesh_test`). Error cases pair `PASS_REGULAR_EXPRESSION` with `FAIL_REGULAR_EXPRESSION "did NOT stop"` (`test_ice_mass_missing_col_np1`). `fesom_conserve` already has ctests (`:278-291`).
  - Meshes:
    - FESOM2 pi: the FESOM2 dist layout, like production dist_512; dist_2 and dist_8.
    - pi_w, core2_w, nemo025mix_cut and nemo025q_cut: the fesom_ini layout; dist_2, dist_4 and dist_8.
  - `test_vinv` runs at np 2 on pi and on the cut-outs. Its fields are analytic, equal on every copy. V8–V12 cover `rvo_upwind`, and V11 uses u = CU·lat² (`:569-575`).
  - `test_ice_dyn_mixed` runs at np 2 on the cut-outs and np 8 on core2_w, with EVP and mEVP. Its fields are hashes of the global id.
  - Gates:
    - `tools/run_restart_gate_core2.sh`: `BUILD=`; the floor is set at `:85` (comment `:82-84`).
    - `tools/run_restart_gate_multirank.sh`: `BUILD=`.
      - **Its run dir defaults to `/scratch/a/a270088/restart_gate_multirank` (`:21`), which it deletes with `rm -rf` (`:32`).**
      - Both restart gates default `F3` to a270088's worktree (`:14`, core2 `:30`).
      - Always pass `F3` and the run dir.
    - `tools/run_tri_regression.sh`: `REF_BIN`, `NEW_BIN`, `ROOT` and `TOL_*` can be overridden. Float32 outputs get a 2⁻²¹ allowance even with a tolerance of 0 (`:264`), and `:283` greps `step 20 ` in the conserve logs.
    - `tools/run_conserve_pi.sh` (`BIN=`, which defaults to `build_intel_dp`, so always set it):
      - its heredoc config loop runs at np 1 (`:88-97`);
      - multi-rank cases have their own blocks (`:71`, `:116-141`).
- **Build:** sources are found with `CONFIGURE_DEPENDS` globs (`CMakeLists.txt:100-111`), so a new module under `src/infra/` needs no CMake change. The same globs mean a build always compiles the tree it points at.

## Development Approach
- **Testing approach:** Regular (code first, then tests in the same task).
- Complete each task fully before moving to the next; small, focused changes.
- **Every task includes new or updated tests; all tests pass before the next task starts.**
- Update this plan when scope changes.
- np 1 has no shared elements, so its results must stay bit-identical. Tasks 3 and 7 check this explicitly.
- np>1 results change at roundoff level, on purpose. `FESOM3_SHARED_SYNC=0` must reproduce the old np>1 results bitwise; Task 3 checks this.
- Build only in scratch dirs, never in `build_intel_dp`:
  - `S=/work/ab0995/a270029/fesom3_hybrid/shared_sync`.
  - `$S/src_base`: a **frozen copy** of the tree before Task 1. A build of the live tree would pick up the edits through the globs.
  - `$S/build_base`: the baseline build, from `$S/src_base`.
  - `$S/build`: the work, from the repository.
  - Rebuild all targets after any flags change.
- git is not on PATH: use `/sw/spack-levante/git-2.31.1-25ve7r/bin/git`.
- ctest: `cd $S/build && bash -c 'ulimit -s unlimited; ctest -j 8'`.

## Testing Strategy
- **New ctest `test_shared_elem`** runs at np 1; np 2 and 8 on pi (the production dist layout); np 2 and 4 on pi_w; np 2 on nemo025mix_cut and nemo025q_cut (quads); np 8 on core2_w. It has three parts:
  - S1: the list, checked against an independent oracle.
  - S2: the exchange and the copy-spread probe.
  - S3: the mechanism, through the real kernels `update_vel` and `compute_hbar_ale`, with the sync on and off.
- **Error cases**, at np 2 on pi (the builder runs only when npes > 1): an invalid knob value, and a vertex with no known owner.
- **Guard ctests:**
  - `test_vinv` V13: `uv_rhsAB` copies bitwise equal with `rvo_upwind = 1`, on V11's field with varying vorticity.
  - `test_ice_dyn_mixed` D6: sigma copies bitwise equal, with EVP and with mEVP.
- **Full model:**
  - A new `fesom_conserve` ctest at np 2: linfs, vector-invariant scheme, `rvo_upwind = 1`, sea ice, `FESOM3_DIVSUM_TOL=1e-12`.
  - `run_conserve_pi.sh` adds np 2 and 8 configs, also with the scalar momentum advection (momadv 2, the default).
  - Each requires a divergence-sum ratio below 1e-12 and copy spreads of exactly 0 in every step.
- **Gates:**
  - The CORE2 restart gate becomes strict at np>1.
  - The multirank restart gate and the triangle regression keep their tolerances.
  - Bitwise checks against `$S/build_base` with `run_tri_regression.sh`: np 1, and np 2 with the sync off. Every printed d must be 0, because the script admits a 2⁻²¹ allowance on float32 output even at tolerance 0.
- **Cost:** a subset timer reports the sync's share of the step. The share is reported, not gated.
- This project has no UI or e2e tests.

## Progress Tracking
- Mark completed items `[x]` immediately; new tasks get ➕, blockers ⚠️.
- Update this plan if the implementation deviates from it.

## Solution Overview
- **Approach A (chosen): a list built at run time, with its own exchange.**
  - At the end of `read_mesh_local`, each rank builds `partit%com_elem2D_shared` from data it already has: its owned elements, their vertices and the owner of every vertex (from `com_nod2D`).
  - A shared element's owner set is the set of its vertex owners. The designated owner is the **highest rank** in that set, which holds the copy the restart and output writers keep.
  - The designated owner sends the element to every other owner. Both sides build their lists from the same owner sets and sort them by global element id, so they agree on the order without extra communication. A one-time self-check proves it.
  - The list is built even when the sync is off, so the copy-spread diagnostic can still show the drift.
- **Rejected:**
  - **B, merging into the UV halo exchange.** It would change which owner feeds each halo element, even with the sync off. A second list would also have to stay consistent with `com_elem2D`.
  - **C, identical copies through a fixed summation order in the viscosity filter.** It changes np 1 results (FESOM2 bit-identity). The next order-dependent kernel would bring the leak back, and copies that already differ would never be repaired.
  - **fesom_ini writes the list into `dist_N`.** No help for the pool dist_512 used in production, and it changes the file format.
  - **Compute each element on one owner only.** That means refactoring every element loop.
- **Sync point:** `update_vel`, after the owned update and before `exchange_elem(UV)`.
  - The eDim halo copies then carry the synced value too, whichever owner `com_elem2D` names as their sender.
  - `compute_ssh_rhs_ale` reads the unsynced `UV_rhs` (`oce_ssh_rhs.F90:345-358`). Its copies differ only by roundoff, which affects the SSH solve but not the `hbar` budget.
- **Frequency:** every step.
  - `FESOM3_SHARED_SYNC=1` is the default. `0` turns the sync off, giving the old behaviour for A/B comparisons.
  - There is no every-N option. `ssh_rhs_old` is computed from `uv` right after `update_vel`, so every step with unequal copies leaks; N>1 would only slow the leak.
- **Fields:** `uv` only.
  - `uv_rhsAB` and sigma are pinned by the guard tests and, in the full model, by the copy-spread diagnostic.
  - Fallback if a guard fails:
    - `uv_rhsAB(1,:,:,:)`: sync right after the momentum-advection call and before the AB blend (`oce_dyn_velrhs.F90:137-140`). Both the blend and the stored history then use the synced value.
    - sigma: sync at the end of `EVPdynamics` and `EVPdynamics_m`, once per ice step. This needs a rank-1 `exchange_elem_shared_2D_r`, which is added only then.
- **Diagnostic:** with `FESOM3_STEP_DIAG`, the step line also prints three things:
  - the divergence-sum ratio `|Σ ssh_rhs_old| / Σ|ssh_rhs_old|` over owned nodes; a leak shows as a ratio above about 1e-12;
  - the copy spread of `uv`, which shows the drift when the sync is off;
  - the copy spreads of `uv_rhsAB` and sigma; a broken guard shows as a nonzero spread.
- **Cost:** one extra element exchange per step.
  - At 512 ranks a rank receives and sends about 50 elements, each 2 × (nl−1) values. A step already does many exchanges.
  - A subset timer "(of which) shared sync" under "dynamics ssh" reports the sync's share of the step. The expectation is well below 1 %; the number is reported, not gated.

## Technical Details
- **Data in `t_partit`:**
  - `type(com_struct) :: com_elem2D_shared`.
    - Receive blocks: this rank's non-designated copies (local element numbers in `1..myDim_elem2D`), one block per designated owner.
    - Send blocks: the elements this rank is the designated owner of, one block per other owner.
  - `logical :: sync_shared_elem = .true.`: set from the knob by the builder. Tests set it directly.
  - Not serialized: `write_partit`/`read_partit` stay as they are. Check that `test_types`' round trip compares only the serialized fields. `set_empty_comms` clears the new list.
- **Builder** `build_com_elem2D_shared(mesh, partit)`: new module `src/infra/mod_shared_elem.F90`, called at the end of `read_mesh_local`.
  - `partit` becomes `intent(inout)` in `read_mesh` and `read_mesh_local`.
  - Re-entrant: on entry it deallocates and resets the old lists. The corrupt-owner test calls it a second time.
  - Every error stop goes through `die` (made public in `mod_partitioning`), so a failing rank aborts the job instead of leaving the others waiting in a collective.
  1. **Knob (Task 3).** `FESOM3_SHARED_SYNC` unset or `1` → on; `0` → off; any other value → `die('FESOM3_SHARED_SYNC must be 0 or 1')`. Value-based, like `FESOM3_MOMADV_OPT`.
  2. **Owner map** `own(1:myDim_nod2D+eDim_nod2D)`. Owned nodes get `mype`. For each receive block p of `com_nod2D`, set `own(rlist(rptr(p):rptr(p+1)-1)) = rPE(p)`.
  3. **Owner sets.** For each owned element e, S is the set of distinct `own(v)` over `v = elem2D_nodes(1:elem2D_nnodes(e), e)`.
     - If a vertex lies outside `1..myDim+eDim_nod2D` or has no owner, stop with `shared-element list: element <global id> has a vertex with no known owner`.
     - If |S| > 1, the designated owner is d(e) = max(S). If d is mype, add e to the send block of every r in S \ {mype}. Otherwise add e to the receive block of d.
  4. **Lists.**
     - Sort each block by global id (`myList_elem2D(e)`), with partners in ascending order.
     - Fill `rPE/rptr/rlist` and `sPE/sptr/slist` with `rptr(1) = sptr(1) = 1`, as `read_com` does.
     - More partners than `MAX_NEIGHBOR_PARTITIONS` stops with the existing "Increase MAX_NEIGHBOR_PARTITIONS" message. The partners are a subset of the `com_nod2D` partners, so a dist that passed `read_com` cannot trigger this; the check stays as a guard.
  5. **Self-check, once.**
     - Put `myList_elem2D` of the owned elements into an integer array and send it through `exchange_elem_shared`.
     - Every received slot must equal its own global id; otherwise stop, naming the rank pair.
     - This proves both sides list the same elements in the same order.
  6. **Log on rank 0:** `shared-element sync: ON|OFF, <n> shared elements, <c> extra copies, max per rank: recv <r>, send <s>`.
     - n is counted once, at the designated owner; n and c are summed over ranks.
     - r and s are maxima, through a new `allreduce_max_i0` in the generic `allreduce_max`.
- **Exchange** (`mod_halo`): a generic `exchange_elem_shared`. Both variants are no-ops on an empty list, and so at np 1.
  - `_blk_r(arr(:,:,:))` on `core_blk_r`, for `uv (2,nl-1,elem)` and the section `uv_rhsAB(1,:,:,:)`.
  - `_2D_i(arr(:))` on `core_2D_i`, for the self-check.
  - A rank-1 real variant is added only if the sigma fallback needs it.
- **Copy-spread probe** (`mod_halo`): a generic `shared_copy_max_elem(arr, partit)` for rank-1 and rank-3 real arrays.
  - It returns this rank's max|arr − designated copy| over its receive slots. The value is local, so the step diagnostic keeps a single `MPI_MAX`.
  - Rank 1 reuses `probe`; rank 3 gets a new block probe on `core_blk_r`.
- **`update_vel`** (`oce_ale.F90:158`):
  ```fortran
  if (is_multirank(partit)) then
      if (partit%sync_shared_elem) then                 ! one owner per shared element
          call timer_start(TMR_OCE_SHSYNC)
          call exchange_elem_shared(UV, partit)
          call timer_stop(TMR_OCE_SHSYNC)
      end if
      call exchange_elem(UV, partit)                     ! FESOM2 :172 (eDim com_elem2D)
  end if
  ```
  - The ifs are nested because `partit` is optional and Fortran does not short-circuit.
  - `TMR_OCE_SHSYNC` is a new subset timer, "(of which) shared sync", whose parent is `TMR_OCE_SSH`, like `TMR_OCE_SOLVE`.
- **Step diagnostic** (`fesom_lifecycle_native_mr.F90:1023-1039`):
  - `dstat(7) = shared_copy_max_elem(dyn%uv)`.
  - `dstat(8) = shared_copy_max_elem(dyn%uv_rhsAB(1,:,:,:))`.
  - `dstat(9)` = the largest of the three sigma spreads, or 0 when sigma is not allocated.
  - The `MPI_MAX` grows from 6 to 9 values.
  - One `MPI_SUM` of `[Σ ssh_rhs_old(1:nNodO), Σ |ssh_rhs_old(1:nNodO)|]` gives the ratio.
  - The fields `divsum=`, `dUV=`, `dAB=` and `dsig=` are appended **after `Smax`**, so `run_cutout_short.sh` keeps its field numbers.
  - With linfs the ratio is expected below 1e-12. With zstar, `ssh_rhs_old` also carries `−water_flux·areasvol`, so the ratio is only reported.
- **`fesom_conserve`:**
  - It computes the same four numbers every step, but **prints no new per-step line**: `run_tri_regression.sh:283` greps `step 20 `.
  - Their run maxima go into one end-of-run line, inside the tail that `run_conserve_pi.sh` shows.
  - New knob `FESOM3_DIVSUM_TOL` (unset = report only). It stops with an error if the max ratio exceeds the tolerance, or if the sync is on and any spread is nonzero.
- **Restart gate:** the np>1 floor (`run_restart_gate_core2.sh:85`) becomes 0 while the sync is on (the default). It stays at 1e-12 with `FESOM3_SHARED_SYNC=0`.

## What Goes Where
- **Implementation Steps:** code, ctests, gate scripts and docs in this repository.
- **Post-Completion:** production reruns, which the user runs and configures.

## Implementation Steps

### Task 1: Build the shared-element list

**Files:**
- Create: `src/infra/mod_shared_elem.F90`
- Modify: `src/types/mod_partit.F90`
- Modify: `src/infra/mod_halo.F90` (integer exchange wrapper, `allreduce_max_i0`)
- Modify: `src/infra/mod_partitioning.F90` (`set_empty_comms`, `die` public)
- Modify: `src/mesh/mod_mesh_read.F90`
- Create: `test/test_shared_elem.F90`
- Modify: `test/CMakeLists.txt`

- [x] Before changing any code, freeze the baseline and build it:
  - Export `HEAD` with `git archive` into `$S/src_base`.
  - Copy the uncommitted `src/io/mod_io_zarr.F90` over it.
  - Build `$S/build_base` from `$S/src_base`. Tasks 3 and 7 compare against it.
- [x] Add `com_elem2D_shared` and `sync_shared_elem` (default `.true.`) to `t_partit`. Clear the list in `set_empty_comms`, and make `die` public.
- [x] Add the integer variant `exchange_elem_shared_2D_i` (generic `exchange_elem_shared`) and `allreduce_max_i0` to `mod_halo`.
- [x] Implement `build_com_elem2D_shared` (Technical Details, builder steps 2–6), re-entrant and stopping through `die`. Call it at the end of `read_mesh_local`, and make `partit` `intent(inout)` in `read_mesh` and `read_mesh_local`.
- [x] Create `test_shared_elem` with S1. The independent oracle: allgather the owned global node ids into a global owner map, then recompute every owned element's owner set and designated owner. Check that:
  - an element is in the receive block of its designated owner exactly when it is shared and this rank is not that owner;
  - it is in the send block for every other owner exactly when this rank is the designated owner;
  - blocks are ascending by global id and partners are ascending;
  - every shared element has exactly one designated owner (a per-element count, summed over ranks);
  - at np 1 the list is empty.
- [x] Add the error-case ctest at np 2 on pi. With `FESOM3_TEST_SHARED_CORRUPT=1` the test clears the receive blocks of `com_nod2D` and calls the builder again. The test passes on `PASS_REGULAR_EXPRESSION "has a vertex with no known owner"`. It fails on `FAIL_REGULAR_EXPRESSION "did NOT stop"`, which the test prints if the call returns.
- [x] Register `test_shared_elem` at np 1; np 2 and 8 on pi; np 2 and 4 on pi_w; np 2 on nemo025mix_cut and nemo025q_cut; np 8 on core2_w.
- [x] Run ctest: all pass, including `test_types`.
- ➕ Additions to the builder (12 Oct): before the self-check, one `MPI_Alltoall` of the per-partner send counts checks that both sides of every pair agree; a mismatch would otherwise hang the self-check exchange instead of stopping with a message. It also stops if an owned element has no owned vertex (the partition rule the design relies on). S1 checks the rule from the other side too: every vertex owner owns the element.
- Measured: pi 35 shared elements at np 2, 226 at np 8; core2_w 2,099 at np 8 (2,102 extra copies).

### Task 2: Block exchange and the copy-spread probe

**Files:**
- Modify: `src/infra/mod_halo.F90`
- Modify: `test/test_shared_elem.F90`

- [x] Add `exchange_elem_shared_blk_r` to the generic `exchange_elem_shared`.
- [x] Add the generic `shared_copy_max_elem`: rank 1 through `probe`, rank 3 through a new block probe on `core_blk_r`. It returns the local maximum.
- [x] S2 block test: fill a `(2,nl-1,elem)` array and the `(1,:,:,:)` section of a `(1,2,nl-1,elem)` array with `1000*mype + global id`. Check that:
  - before the exchange, the block probe (allreduce max) is positive at np>1;
  - after it, every shared copy holds `1000*d(e) + global id` and the probe is 0;
  - non-shared owned slots and halo slots are bitwise unchanged.
- [x] S2 rank-1 probe test: plant a known difference on one receive slot of a rank-1 array. The probe must return exactly that difference on that rank and 0 on the others.
- [x] S2 at np 1: the arrays stay bitwise unchanged and the probes return 0.
- [x] Run ctest: all pass.

### Task 3: Sync `uv` in `update_vel`, with the `FESOM3_SHARED_SYNC` knob

**Files:**
- Modify: `src/infra/mod_shared_elem.F90`
- Modify: `src/infra/mod_timer.F90`
- Modify: `src/oce/oce_ale.F90`
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90` (env-block comment only)
- Modify: `test/test_shared_elem.F90`
- Modify: `test/CMakeLists.txt`

- [x] In the builder, read the knob (builder step 1) into `partit%sync_shared_elem`, and print ON or OFF in the log line.
- [x] Add the subset timer `TMR_OCE_SHSYNC`, "(of which) shared sync", under `TMR_OCE_SSH`.
- [x] In `update_vel`, call `exchange_elem_shared(UV)` before `exchange_elem(UV)`, timed and as shown in Technical Details. Add a comment saying why.
- [x] Add `FESOM3_SHARED_SYNC` to the env-block comment of the production driver, so that block stays the full list of run variables.
- [x] S3 mechanism test, through the real kernels:
  - Set `uv` and `uv_rhs` from hashes of the global id, so they are equal on every copy. Add δ = 1e-3 m/s to `uv_rhs` on this rank's receive slots (the non-designated copies). Set `d_eta = 0`.
  - Allocate `hbar`, `hbar_old`, `dhe`, `ssh_rhs_old` and `helem` the way `fesom_conserve` and `test_vinv` do.
  - Call `update_vel`, then `compute_hbar_ale`. Compute r = |Σ ssh_rhs_old| / Σ|ssh_rhs_old| over owned nodes (allreduce).
  - Sync on: r < 1e-12 and `shared_copy_max_elem(uv)` is 0.
  - Sync off (`partit%sync_shared_elem = .false.`): r > 1e-8 and the probe ≈ δ.
  - At np 1: r < 1e-12 either way.
- [x] Error and OFF cases as ctests at np 2 on pi:
  - `FESOM3_SHARED_SYNC=2` must match `"FESOM3_SHARED_SYNC must be 0 or 1"`, paired with `FAIL_REGULAR_EXPRESSION "did NOT stop"`.
  - `FESOM3_SHARED_SYNC=0` must still build the list and log `shared-element sync: OFF`.
- [x] np 1 bit-identity:
  - Command: `REF_BIN=$S/build_base/bin NEW_BIN=$S/build/bin ROOT=$S/tri TOL_PRE=0 TOL_EARLY=0 TOL_OCEAN=0 TOL_ICE=0 TOL_MIX=0 TOL_ICE_EARLY=0 tools/run_tri_regression.sh prod:pi:1 prod:core2:1`.
  - Every printed d must be 0. A PASS alone is not enough, because float32 output gets a 2⁻²¹ allowance.
- [x] OFF bit-identity: run the same command with `FESOM3_SHARED_SYNC=0` exported and the case `prod:pi:2`. Every printed d must be 0.
- [x] Run ctest: all pass.
- Measured (12 Oct): S3 ratio with the sync on 9e-20 to 2e-17 on every mesh, off 7e-6 (pi np 2) to 4e-4 (pi_w np 4), uv spread off = δ exactly. np 1 bit-identity: pi and CORE2, 98 printed d, all 0. OFF bit-identity: pi np 2, 49 printed d, all 0; the startup log says `shared-element sync: OFF, 35 shared elements`.
- ➕ The timer name field (`TNAME`) is now 24 characters wide, so "(of which) shared sync" fits; the report prints the existing names unchanged.

### Task 4: Guard tests for the `uv_rhsAB` and sigma copies

**Files:**
- Modify: `test/test_vinv.F90`
- Modify: `test/test_ice_dyn_mixed.F90`

- [x] Add `test_vinv` V13 at np>1 (V12 is taken: "blend linear in rvo_upwind", `:642-658`):
  - Use V11's velocity field, u = CU·lat², v = 0, which has varying vorticity. `set_shear` has constant vorticity, so the blend would be inert there.
  - Call `relative_vorticity` and `momentum_adv_vinv` with `rvo_upwind = 1.0` and `uv_rhsAB = 0`.
  - The allreduce max of `shared_copy_max_elem(dyn%uv_rhsAB(1,:,:,:))` must be exactly 0.
  - So the guard is not vacuous, also require that on some shared element the blend changes `uv_rhsAB` by more than 1e-6 relative to `rvo_upwind = 0`.
- [x] Add `test_ice_dyn_mixed` D6 at np>1:
  - Use a hashed ice velocity: mEVP sigma stays 0 at rest. Run at least 2 subcycles.
  - Run `EVPdynamics_solve` with EVP, then with mEVP. After each, the allreduce max of `shared_copy_max_elem` of `sigma11`, `sigma12` and `sigma22` must be exactly 0.
  - So the check is not vacuous, require nonzero sigma on some shared element.
- [x] (not needed) If a guard fails: ➕ add the sync for that field at the place given under Fields in Solution Overview, gated by `partit%sync_shared_elem`. For sigma, also add `exchange_elem_shared_2D_r` with an S2-style test. Record the cause in this plan and rerun the guard.
- [x] Run ctest: all pass.
- Measured (12 Oct): both guards hold, no fallback sync. V13: spread 0, blend change on the shared elements 0.98 (pi np 2), 0.71 (mixed), 0.79 (quad). D6: spread 0 with EVP and mEVP; max|sigma| on the shared elements 1.0e3 / 17 Pa (mixed np 2), 3.6e2 / 6.4 Pa (quad np 2), 3.8e3 / 66 Pa (core2_w np 8).

### Task 5: Leak and copy-spread diagnostics in the drivers

**Files:**
- Modify: `src/drivers/fesom_lifecycle_native_mr.F90`
- Modify: `src/drivers/fesom_conserve.F90`
- Modify: `tools/run_conserve_pi.sh`
- Modify: `test/CMakeLists.txt`

- [x] In `fesom_lifecycle_native_mr`, extend the `FESOM3_STEP_DIAG` block as described in Technical Details, with the new fields after `Smax`.
- [x] In `fesom_conserve`, add the per-step numbers (no new per-step line), the end-of-run line with their maxima, and `FESOM3_DIVSUM_TOL`.
- [x] Add multi-rank configs to `run_conserve_pi.sh`, in their own blocks like the `for np` loop at `:71`, since the heredoc loop runs at np 1 only. Common settings: pi, linfs, `FESOM3_CONSERVE_ICE=1`, 20 steps.
  - Sync on, `FESOM3_DIVSUM_TOL=1e-12` (a hard gate: ratio < 1e-12, spreads 0):
    - `FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0` at np 2 and 8;
    - `FESOM3_MOMADV_OPT=2` (scalar advection, the default) at np 2.
  - Sync off (`FESOM3_SHARED_SYNC=0`), vector-invariant scheme, np 2: report only. Record whether `dUV`/`dAB` > 0, i.e. whether the copies drift without the sync.
- [x] Register the ctest `fesom_conserve_linfs_vinv_shsync_np2` (pi, the first sync-on config) following `fesom_conserve_zstar_np*` (`test/CMakeLists.txt:278-291`).
- [x] Run `BIN=$S/build/bin/fesom_conserve tools/run_conserve_pi.sh` (all configs): it must be green. Write the measured ratios and spreads into the script comment.
- [x] Run ctest: all pass (167/167).
- Measured (12 Oct, pi, 20 steps): max ratio 1.6e-16 (vinv np 2), 3.1e-16 (vinv np 8), 1.8e-16 (scalar np 2), spreads 0. Sync off: `dUV` 6.9e-18, `dAB` 5.8e-11, i.e. the copies drift and the blend amplifies the drift by ~1e7 within 20 steps; the ratio stays 1.6e-16 this early. Every other config of the script reports ratios <= 3.8e-16 and zero spreads.
- ➕ `run_conserve_pi.sh`: every driver run prints one more line, so every `tail` grew by one. The script needs `source env.sh intel` first (mpirun).
- ➕ README also lists `FESOM3_STEP_DIAG` with its new fields.

### Task 6: Gates and cost

**Files:**
- Modify: `tools/run_restart_gate_core2.sh`
- Modify: `tools/run_restart_gate_multirank.sh` (comment only)
- Modify: `tools/zarr_diff.py`, `tools/dump_diff.py` (help text only)

- [x] In `run_restart_gate_core2.sh`, make the np>1 floor (`:85`) 0 unless `FESOM3_SHARED_SYNC=0`. Update the comment at `:82-84` and the GREEN message at `:151`.
- [x] Update the "562 redundantly-owned" comments in `run_restart_gate_multirank.sh:8-10`, `zarr_diff.py:573, 612, 637-639` and `dump_diff.py:105, 254`.
- [x] Run `F3=/home/a/a270029/fesom3.0 BUILD=$S/build tools/run_restart_gate_core2.sh $S/rg_core2 1 2`: np 1 and np 2 must be exact (max|Δ| = 0) in both the mid-run and the boundary case.
- [x] Run the same again with `FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0` exported and the run dir `$S/rg_core2_vinv`: it must be exact too.
- [x] (not needed: np 2 exact) If np 2 is not exact: locate the differing store and elements first, and do not set the floor to 0 before the cause is understood. ⚠️ Record it here.
- [x] Run `F3=/home/a/a270029/fesom3.0 BUILD=$S/build tools/run_restart_gate_multirank.sh $S/rg_mr`. **Never run it without a run dir**: the default is another user's scratch dir, which the script deletes.
- [x] Run `NEW_BIN=$S/build/bin ROOT=$S/tri tools/run_tri_regression.sh` (default cases and tolerances). Both gates must be green. Record the np>1 differences.
- [x] Cost, reported and not gated:
  - Run `fesom_lifecycle_native_mr` on CORE2 at np 8 for 48 steps, with the production physics from `tools/prod_env.sh` and no output or restart.
  - Record the share of "(of which) shared sync" in the FINAL timer report (expected well below 1 %).
  - As a cross-check, record the step time with `FESOM3_SHARED_SYNC=0`.
- Measured (12 Oct): CORE2 restart gate np 1 and np 2, mid-run and boundary, checkpoint and live dumps `max|Δ| = 0` with the floor at 0, default physics and vinv + `rvo_upwind` 1. Multirank gate np 2 → 8 `max|Δ| = 0`. Cost: "(of which) shared sync" 0.51 ms/step (max 0.66) of a 974 ms loop on CORE2 np 8, 0.05 %; the sync-off run took 960 ms per step (login-node noise exceeds the sync). Same runs: sync on `divsum` <= 1.4e-15, spreads 0 in all 48 steps; sync off `dUV` up to 3.9e-16, `dAB` up to 5.8e-11.
- Triangle regression (default cases and tolerances, vs `frozen_task6_soltol`): 7/7 pass. np>1 differences at the levels measured before the sync: pi np 2 OCEAN 9.7e-8, ICE 3.9e-7; CORE2 np 8 EARLY 2.2e-8, OCEAN 3.1e-5, ICE 9.7e-5, MIX 7.8e-4; decomp pi_w np 2 / 8 OCEAN 2.2e-11 / 8.8e-10; conserve np 2 final salt 3.5e-16.

### Task 7: Verify acceptance criteria
- [x] At np>1 with the sync on, the `uv` copies are equal after every `update_vel` (S3, and `dUV = 0` in the conserve configs).
- [x] The `uv_rhsAB` and sigma copies stay bitwise equal (Task 4 guards, `dAB = dsig = 0` in the conserve configs). If not, the fallback sync is in place and recorded.
- [x] At np>1 with the sync on, the divergence-sum ratio is below 1e-12, including after a deliberate perturbation. With the sync off it is clearly nonzero (S3).
- [x] Rerun the np 1 bit-identity and OFF bit-identity commands from Task 3, now that Task 5 has edited the drivers: every printed d is 0.
- [x] The restart gate is exact at np 2, with the default and with the vector-invariant settings (Task 6).
- [x] The sync's share of the step is recorded (Task 6).
- [x] Full ctest passes: `cd $S/build && bash -c 'ulimit -s unlimited; ctest -j 8'`.
- Measured (12 Oct): np 1 rerun after Task 5, pi and CORE2: 98 printed d, all 0; OFF rerun pi np 2: 49 printed d, all 0. CORE2 np 8 production physics, 48 steps: `divsum` <= 1.4e-15, `dUV = dAB = dsig = 0`. Offline check of production dist_512: the ownership rule holds (no owned element without an owned vertex, every vertex owner owns the element), 24,798 shared elements, 25,512 extra copies, max per rank recv 139 / send 136.

### Task 8: [Final] Update documentation
- [x] `docs/HANDOFF.md`:
  - Add a shared-element section: the ownership rule, the designated owner (highest rank = the writer's copy), the sync point, the knob, the diagnostics and the guards.
  - Add the CORE2 leak history (results7/11/17/21).
  - Correct the M10 restart note (`:869`): the "user-accepted caveat" is resolved, and np>1 restarts are exact.
- [x] `docs/LESSONS.md`, new lesson L70:
  - A shared element is a duplicated prognostic variable.
  - An order-dependent sum makes the copies differ by roundoff; a branchy operator makes the difference O(1); duplicated border fluxes then turn it into a mass leak.
  - Test the global divergence sum at np>1, not only conservation at np 1.
- [x] `TESTING.md` §1, §7 and §8:
  - The redundant-element floor is gone with the sync; it applies only with `FESOM3_SHARED_SYNC=0`.
  - Add the new ctests and the strict np>1 restart gate.
- [x] `IMPLEMENTATION.md`: add `com_elem2D_shared`, `exchange_elem_shared` and `shared_copy_max_elem` to the parallel-layer description.
- [x] `README.md`: add `FESOM3_SHARED_SYNC` to the run-variables table. Note that it is read in `mod_shared_elem`, the way the SSH solver reads its own knobs.
- [x] `docs/HANDOFF-next-session.md`: replace the "implement the plan" item with the result.
- [x] Move this plan to `docs/plans/completed/`.
- Notes: `docs/HANDOFF-next-session.md` stays untracked, as before. Code comments, tests and tools cite the plan at its `completed/` path. `TESTING.md` §5 and §6 were updated too (ctest count 167, the strict restart gate).

## Post-Completion
*For the user: runs and run settings are the user's.*

**Production verification:**
- Rebuild the production binary when you decide. The work builds only in `$S/build`.
- Resume results17 from `fesom.1967.001.00000` (or results21 from `fesom.1959.001.00000`) with the new binary and the same settings, at least into autumn 1967. The leak paused for about six months after each restart, so earlier months prove nothing.
- With `FESOM3_STEP_DIAG=1`, for a short segment or the whole run, the step line should show:
  - `divsum` around 1e-15 to 1e-13;
  - `dUV = dAB = dsig = 0`.
- Expect the global-mean SSH to stay at its restart value (0.0888 m for results17), and the global salt trend to return to about −2e-6 psu/yr.
- The startup log at dist_512 should report 24,798 shared elements, the same as the offline count.
- Optional: from the same restart, run one segment with `FESOM3_SHARED_SYNC=0`; `dUV` > 0 should show the drift, and the leak should return. Run one with `FESOM3_RVO_UPWIND=0` to confirm the trigger.
- The FINAL timer report shows the sync's share of the step at 512 ranks.

**Report:**
- Once the production check is done, the CORE2 runs page (artifact `Gdj6DMQ4n2BQVZ4Si6cAUD`) can get a short note.
