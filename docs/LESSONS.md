# FESOM3 — Lessons (mistakes + fixes, never repeated)

Append-only log of non-obvious gotchas. Each entry: what bit us, why, the fix.

## L1 — Anchor compile flags come from FESOM2 v2.7.3, NOT the tracer_dwarf

The byte-gate oracle is **FESOM2 v2.7.3** (`/home/a/a270088/port2/fesom2/src/CMakeLists.txt`).
The tracer_dwarf is a *structural* reference whose flags differ and were never byte-gated.
Differences that would silently break bit-identity if we copied the dwarf:

- dwarf adds `-no-prec-sqrt -ip` to **ifx** too; FESOM2 restricts them to **ifort classic**.
- dwarf adds `-march=core-avx2 -mtune=core-avx2` on Levante Intel; **FESOM2 v2.7.3 COMMENTS
  THAT LINE OUT** (`src/CMakeLists.txt:356`). So the FESOM2 Levante Intel build targets the
  SSE2 baseline (no AVX2, **no FMA contraction**). Adding `-march=core-avx2` enables FMA →
  different bits. Our anchor therefore uses **no `-march`** on Levante.

Faithful anchor (Intel classic, Levante), from FESOM2 v2.7.3 `src/CMakeLists.txt:335-340`:
```
-O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero
-no-wrap-margin -fpe0 -fpp  -no-prec-sqrt -ip
```
GNU (FESOM2 `:396`,`:434`): `-O3 -ffloat-store -finit-local-zero -finline-functions
-fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp
-march=znver3 -mtune=znver3 -ftree-vectorize -flto`.

Codified in `cmake/fesom_flags.cmake`. `-fpe0` was assumed to imply flush-to-zero (FTZ) of
denormals — a bit-affecting behavior we must match, so keep it. **CAVEAT (see L51):** `-fpe0`
traps FP *exceptions* but does NOT pin the runtime MXCSR FTZ/DAZ bits — those are per-process and
can differ or be reset by a library, which is exactly what bit the 2-yr headline at day 107. FESOM3
asserts FTZ explicitly via `ieee_set_underflow_mode(gradual=.false.)`. `-init=zero` zeroes locals
(prevents spurious `-fpe0` traps on uninitialised reads). When a real byte-gate runs
(M0.7+), re-verify the actual FESOM2 build's flags from its build dir, not just the CMake.
**VERIFIED at the M1 geometry gate:** both build with ifort 2021.5.0 via the same mpif90,
and the FESOM2 build flags (`build/src/CMakeFiles/fesom.dir/flags.make`) are byte-for-byte
the anchor list above. Same compiler + same flags ⇒ scalar SVML transcendentals match.

## L7 — `-no-prec-div`: divide by a LITERAL, not a runtime `real(nv)` (arity trap)

The single hardest bug in the M1 geometry gate. The anchor flag `-no-prec-div` lets ifort
implement `x/y` as `x * recip(y)`. For a **compile-time constant** divisor (`x/3.0_WP`)
ifort folds an exact `1/3` and the result is one value; for a **runtime** divisor
(`x/real(nv,WP)`, nv from `elem2D_nnodes`) it emits a hardware reciprocal approximation
that differs by **1 ULP**. That 1 ULP in a centroid latitude `sum(lat)/3` flowed into
`cos`/`tan` (steep near the rotated pole → up to ~25% of elements differ), then
`elem_area`, `gradient_sca`, `area`, `edge_cross_dxdy` — the whole geometry failed the
gate. **Fix:** divide by the literal `3.0_WP`/`3.0_MP` (what FESOM2 writes), NOT
`real(nv,WP)`. This is the plan's arity caveat made concrete: the `/3 → /n_vert`
generalization is NOT universally bit-safe. v1 is triangles only (nv==3); quad support
must branch on nv and divide by the matching literal. **Non-causes ruled out** (cost a
day if not): the intrinsic `sum()` vs an explicit accumulation loop are bit-identical here;
auto-vectorization of `cos` was NOT the cause (`!DIR$ NOVECTOR` had zero effect — removed);
coords + `elem2D_nodes` were proven `max|Δ|=0` first, which is what localised it to the
divide. Lesson: when a geometry byte-gate shows ~1-ULP transcendental diffs, suspect a
runtime divisor under `-no-prec-div` before anything else.

## L8 — 1-rank FESOM2 oracle (hand-crafted dist_1) + geometry byte-faithfulness

To byte-gate anything that accumulates per node/element (`area`, later FCT scatter), the
oracle must run **1-rank**, because FESOM2 builds `nod_in_elem2D` in LOCAL element order
(oce_mesh.F90 find_neighbors), so a multi-rank run reorders the per-node `Σ elem_area/3`
→ ULP diffs vs FESOM3's global order. METIS can't produce `dist_1`; hand-craft it
(see HANDOFF "Geometry byte-gate"). FESOM2's reader accepts blank lines for zero-size
halo arrays (zero-trip `read(*,*)` skips a record). 1-rank forcing init hangs on a login
node — irrelevant: dump geometry at `mesh_setup` and STOP before forcing.
Geometry transcription gotchas the gate caught (all now `max|Δ|=0`):
- `elem_center`/`edge_center` wrap arithmetic must be VERBATIM (FESOM2 wraps each lon vs
  `amin`; edge_center shifts a(1)/b(1) asymmetrically). An "equivalent" rewrite diverged.
- **Vertex (CW) order — gated only on the 0-swap case so far.** pi has 0/5839 elements
  needing the clockwise swap (it's pre-oriented), so the geometry gate did NOT exercise
  `enforce_cw_orientation`'s reorder path. soufflet has 228/5700, CORE2 ~certainly >0.
  Within-element node order is FP-order-sensitive (centroid `sum(lat)/3`; `gradient_sca`
  column→vertex map; elem_area is swap-invariant). FESOM3 `enforce_cw_orientation` is
  byte-identical to FESOM2 runtime `test_tri` (oce_mesh.F90:1706-1726: same b/c, same
  `trim_cyclic` on comp 1, same `r=b1*c2-b2*c1`, same `r>0`→swap 2,3) operating on the
  same elem2d.out + byte-proven coords ⇒ identical swaps by construction. **But empirically
  confirm on a swap mesh (soufflet, or the M2.11 single-rank CORE2 gate) before trusting it.**
- `elem_area`/`area`/`areasvol` accumulate UNSCALED then `*r_earth²` ONCE at the end
  (a single deferred scaling), matching mesh_areas — not per-element scaling.
- `edge_tri.out` stores **-999** for boundary (no 2nd triangle); FESOM2 `load_edges`
  does `where(edge_tri<0) edge_tri=0`. Replicate, or the boundary marker leaks.

## L6 — Oracle is runnable; how to drive it (proven 2026-06-19)

Running the prebuilt instrumented FESOM2 on pi to produce reference dumps works and is
fast (~0.25 s, 2 steps, 2 ranks). Gotchas, all handled by `tools/run_oracle_pi.sh`:
- `fesom.clock` must be non-empty even for a fresh start: two lines `0 1 1948` (the
  forcing year; columns = seconds, day-of-year, year). Empty → `severe (24): end-of-file`.
- `run_length_unit='s'` means STEPS (not seconds) — set `run_length=2` for a 2-step run.
- `ClimateDataPath` in the shipped namelist is relative (`..//test/input/global/`); from a
  scratch run dir it must be absolute. Forcing paths were already absolute.
- Login-node MPI needs `--mca pml ob1 --mca btl self,vader` (UCX/IB settings fail).
- A fresh run byte-matches a prior run's dump (`max|Δ|=0`) → FESOM2 is deterministic, so
  `max|Δ|=0` byte-gating is well-posed. The shipped pi config runs KPP/GM (substep MIXING
  dumps `Kv`); the **reduced M2 config** (PP, no GM/Redi, linfs, opt_visc=7) is a separate
  namelist to apply when byte-gating M2.

## L5 — Mesh geometry (M0.7): ocean-only meshes, units, deferred byte-gate

- **pi is OCEAN-only** (no land elements). Σ elem_area = ocean area ≈ 3.40e14 m² ≈
  0.67·(4πr²), NOT the full sphere. Don't assert full-sphere area; assert the ocean
  fraction (~0.5–0.85). pi extent is global (lon 0–360, lat −78..89), nl=48.
- **FESOM2's geometry convention is radians × r_earth**: coords stored in radians,
  `elem_area *= r_earth²`. A Cartesian analytic mesh must therefore store coords as
  *radians-like* (physical metres / r_earth) so the shared pipeline yields physical
  m². Test a linear field with the physical position = coord·r_earth.
- **gradient_sca assumes CW node order** — run CW enforcement (test_tri) BEFORE
  computing geometry, on file AND analytic meshes. (pi files are pre-oriented: 0 swaps.)
- **FESOM2 geometry byte-gate is DEFERRED to M1**: it needs an instrumented-FESOM2
  run dumping reference elem_area/gradient_sca/areas, which wasn't producible this
  session. M0.7 is gated by self-consistency (CW, ocean-area, gradient annihilates
  constants, gradient reproduces a linear field, adjacency). Deferred geometry not yet
  transcribed: gradient_vec (M2 momentum), mesh_resolution smoothing (M4 GM), and the
  multi-rank mesh remap (M2.12). Edge left/right orientation on the analytic mesh is
  not FESOM2's (analytic only; pi reads edge_tri from file).

## L4 — Type-design choices (M0.3) deviating from the plan letter

Documented here per the plan's "update scope when implementation deviates":
- **Field names kept FESOM2-verbatim** (`elem2D_nodes`, `nod_in_elem2D`, `nlevels`,
  `gradient_sca`, ...) so kernels transcribe line-for-line. The plan's `elem_nodes`/
  `elem_nnodes` naming → realized as `elem2D_nodes(MAX_NV,:)` + new `elem2D_nnodes(:)`.
- **node→element adjacency is DENSE** (`nod_in_elem2D(MAX_ADJACENT,:)` + `_num`), not
  CSR as D2 suggested. FESOM2 and the dwarf both use dense; matching it preserves
  kernel iteration order for byte-identity. (CSR was a memory aspiration, not in the refs.)
- **Types are lean/growable**: dropped beyond-v1 fields (cavity, iceberg, OASIS, split-
  explicit `se_*`, backscatter/UKE, energy diagnostics `ke_*`, DVD). Aux 3D WORK fields
  (density/N²/Kv/Av/hpressure/sw_α-β) are added to `t_dyn_work` at M2.1 with verified
  FESOM2 names rather than guessed now. GM `fer_uv/fer_w` declared, serialized at M4.
- **Serialization = derived-type `read/write(unformatted)`** (DTIO) over an `access='stream'`
  unit; size-prefixed via `mod_binary_arrays`. Only WP-real variants exist — at DP/SP
  MP==WP so MP mesh arrays serialize through them; HP compiles IO out. Full restart-file
  orchestration (`mod_io_restart`) deferred to M2.11.

FESOM2's GNU Levante flags include `-flto`, but FESOM2 links one big executable;
FESOM3 links a static `libfesom3.a` + separate driver/test executables. `-flto` puts
LTO objects in the archive, and gfortran 8.5's `ar`/`ranlib` need the LTO plugin
(`gcc-ar`/`gcc-ranlib`) plus `-flto` at link — the default tools fail with
"plugin needed to handle lto object" and unresolved module symbols. GNU is the
portability build, not the bit-identity anchor (that's Intel), and `-flto` can reduce
FP reproducibility, so we omit it. Revisit only if GNU↔FESOM2-GNU byte-identity is wanted.

## L2 — Precision: `-r8` promotion AND explicit `_WP` kinds move together

FESOM2 relies on `-r8`/`-fdefault-real-8`: unsuffixed literals like `1.0` are **double**.
Transcribe literals verbatim (keep unsuffixed where FESOM2 leaves them unsuffixed) AND build
with `-r8`, or the bits differ. The cpp define `USE_SINGLE_PRECISION` (sets `WP=4`) and the
`-r4` flag are flipped together by `configure.sh --precision sp`. `MP = max(WP,4)`.

## L9 — Prior gates validate hidden orderings; the kernel-gate recipe (M1.1, PASSED first try)

M1.1 (horizontal tracer advection: upw1 + MUSCL + scatter) hit `max|Δ|=0` vs FESOM2 on
the FIRST gate run — no bug-hunt. Two reasons, both reusable:

- **A passing gate transitively proves orderings the next gate depends on.** M1.1's
  `find_up_downwind_triangles` picks a triangle by `atan2` comparisons, and `fill_up_dn_grad`
  area-weight-AVERAGES `tr_xy` over `nod_in_elem2D(:,node)` — both sensitive to the ORDER of
  `nod_in_elem2D`. That order was never dumped directly, but the **geometry gate's `area`
  match (`area(nz,n)=Σ_k elem_area(nod_in_elem2D(k,n))/3`, an FP sum) already pinned it**:
  if the k-order differed, `area` would have differed. So `edge_up_dn_tri` (integer indices)
  and `edge_up_dn_grad` matched by construction. Lesson: when designing a gate, list the
  order-sensitive inputs and check whether an earlier `max|Δ|=0` already constrains them —
  it often does, and that's why faithful transcription "just works."

- **The kernel-gate oracle recipe** (now proven, reuse for M1.2+): wire an env-gated,
  npes==1 dump shim into FESOM2 at the LATEST setup point that has what you need but is
  still before the 1-rank forcing hang (M1.1 = end of `ocean_setup`, after
  `init_thickness_ale`). The shim PRESCRIBES analytic inputs from the byte-proven coords
  (so FESOM2 and FESOM3 generate identical `ttf`/`vel` independently — no captured-input
  file) and calls the **real FESOM2 kernels**, then STOPS. Gate every intermediate, not
  just the final target, so a divergence localises itself. Pin scalar knobs (dt, num_ord)
  as shared constants in BOTH codes, not from a namelist.

- **`helem` at init is `which_ale`-independent.** Both linfs and zstar `init_thickness_ale`
  reduce to `helem(nz,e)=zbar(nz)-zbar(nz+1)` (full cells) at the initial state because
  `hbar=eta=0` there. So a kernel gate that dumps at init can compute `helem` from `zbar`
  directly without importing ALE thickness evolution (deferred to M2.7). Verify it as a
  dumped field anyway (M1.1 did — `max|Δ|=0`), since partial cells/cavity would break it.

## L10 — Vertical advection (M1.2): QR4C wants ALE 3D depths, not static Z/zbar (PASSED first try)

M1.2 (vertical advection: upw1 + QR4C + scatter) hit `max|Δ|=0` vs FESOM2 on the FIRST gate
run, like M1.1. Reusable specifics:

- **QR4C divides by the *per-node ALE* mid-depths `Z_3d_n`/`zbar_3d_n`, NOT the static 1-D
  `Z`/`zbar`.** FESOM2 builds them in `init_ale` (oce_ale.F90:531-566). At the initial state
  (no cavity, full cells) they reduce to: `zbar_3d_n(nz,n)=zbar(nz)` for all nz; interior
  `Z_3d_n(nz,n)=Z(nz)`; and the surface/bottom interfaces use
  `Z_3d_n(nzmin)=zbar(nzmin)+(zbar(nzmin+1)-zbar_n_srf)/2` and
  `Z_3d_n(nzmax-1)=zbar(nzmax-1)+(zbar_n_bot-zbar(nzmax-1))/2` with `zbar_n_srf=zbar(nzmin)`,
  `zbar_n_bot=zbar(nzmax)`. The FESOM2 shim used the LIVE `mesh%Z_3d_n` (init_ale ran earlier
  in ocean_setup); the FESOM3 driver REPLICATES that formula. Both `zbar_3d_n` and `Z_3d_n`
  were dumped as gated fields (`max|Δ|=0`) so the replication is proven, not assumed. NOTE the
  geom gate did NOT cover `Z`/`zbar` (it dumps level *counts* only); the advhor gate now does
  (via the 3-D arrays). FESOM3 `read_mesh` builds `Z=0.5*(zbar(k)+zbar(k+1))` in ONE step vs
  FESOM2's two steps (`Z=zbar(k)+zbar(k+1); Z=0.5*Z`) — byte-identical because ×0.5 is exact.

- **`/2` is byte-exact under `-no-prec-div`; `/3` is not (cf. L7).** The boundary `Z_3d_n`
  `/2` and the QR4C `/3.0_WP` literal divisors fold to exact/identical constants. And QR4C's
  `qc/qu/qd = (ttf diff)/(Z_3d_n(k)-Z_3d_n(k+1))` are RUNTIME divisors — fine here because
  BOTH codes divide by the *same* runtime operand (byte-identical `Z_3d_n`), so the reciprocal
  approximation matches. L7's trap was a runtime-vs-LITERAL *mismatch* between the two codes,
  not runtime division per se. Transcribe the divisor expression verbatim and operands match.

- **The "D=2 shallow-column double-write" is moot on pi** (its min column is `nlevels_nod2D=5`
  = 4 layers; `sort -n nlvls.out | head -1`). The QR4C 2-layer double-write (2nd-layer and
  bottom-1 both hitting interface `nzmin+1`, which nets the flux there to ~0) only triggers at
  `nlevels_nod2D=3`; a 1-layer column (`=2`) would read `ttf(0)` OOB. pi has neither, so
  faithful statement-order transcription is automatically correct. Confirmed with a Debug
  (`-check all`) run of `fesom_advhordump` (EXIT 0, no OOB/FPE). Re-examine for meshes with
  shallower minima (soufflet, CORE2) — there the double-write/OOB order matters and must match.

- **Debug build is NOT byte-comparable to the Release oracle.** A `-check all` (`-O0`, no
  `-no-prec-div`) FESOM3 dump differs from the `-O3` FESOM2 oracle on EVERY FP-computed field
  (horizontal too, and even `area` at ~3.7e-4) — a uniform flag effect, not a bug. Use Debug
  ONLY to catch OOB/uninit/FPE; byte-gate Release-vs-Release (anchor flags both sides, L1).

## L11 — FCT (M1.3): config knobs, the ttf-gradient quirk, the bignumber bottom-fill (PASSED first try)

M1.3 (FCT/Zalesak limiter `oce_tra_adv_fct` + the `use_lo` scatter) hit `max|Δ|=0` vs FESOM2 on
the FIRST gate run, like M1.1/M1.2. The 15 new fields (fct_LO, clipped adf_h/adf_v, fct_ttf_max/
min, fct_plus/minus, del_ttf_*_fct, hnode/hnode_new, standalone MFCT) all matched. Specifics:

- **The pi FCT config is MFCT/QR4C/FCT with `opth=0.0`, `optv=1.0`** (namelist.tra tracer 1/2:
  `1,'MFCT','QR4C','FCT ',0.,1.` → hor.Ord=0 ⇒ 3rd-order MFCT, vert.Ord=1 ⇒ 4th-order QR4C), NOT
  the `num_ord=0.75` M1.1/M1.2 used for the standalone MUSCL/QR4C dumps. Pin opth/optv as shared
  constants in BOTH codes. The blended num_ord term is already proven (M1.1 MUSCL@0.75); a
  standalone `adv_flux_hor_mfct`@0.75 dump finally gates the MFCT kernel (no-clamp variant) fully.

- **`edge_up_dn_grad` is the gradient of `values` (ttf), NOT `valuesAB` (ttfAB).** FESOM2
  `init_tracers_AB` (oce_tracer_mod.F90:126-127) has the `valuesAB` gradient call COMMENTED OUT
  and uses `values`; then `do_oce_adv_tra` reconstructs `valuesAB` (ttfAB) in MUSCL/MFCT/QR4C
  using that ttf-gradient (note line 131 "WHY NOT AB HERE? DSIDOREN!"). So the gate reuses
  `eudg=grad(ttf)` (already built for M1.1) for the `MFCT(ttfAB)` call — do NOT recompute
  grad(ttfAB). LO fluxes use ttf, HO fluxes use ttfAB, both with eudg=grad(ttf).

- **The `a2` bignumber fill clobbers EVERY element's bottom layer.** `do nz=nu1,nl1-1` (real
  bounds) then `if(nl1<=nl-1) do nz=nl1,nl-1` with `nl1=nlevels(elem)-1`, so layer
  `nlevels(elem)-1` (the element's own bottom layer) is set to ∓bignumber, leaving the deepest
  node layer's admissible increment unconstrained. An off-by-one (`nl1` vs `nl1-1`) here would
  diverge — pi has shallow elements so the fill is exercised and gated.

- **AUX scratch: FESOM2 reuses `edge_up_dn_grad`; FESOM3 uses a LOCAL.** Byte-identical on pi
  because `a2` writes `AUX(1:2, ulevels(elem):nl-1, elem)` for every element and (with
  ulevels==1, no cavity) `a3` only reads written entries — the uninitialised local never reaches
  a result. On a CAVITY mesh `a3` reads `AUX(:,nz,elem)` at `nz<ulevels(elem)` (FESOM2 = stale
  edge_up_dn_grad there); unreproducible with a fresh array → **re-gate FCT on a cavity mesh
  (M2+).** In the SHIM, pass a SEPARATE `aux_scratch` (not edge_up_dn_grad) to the real
  `oce_tra_adv_fct` so edge_up_dn_grad survives intact for its own dump record. `dmax1/dmin1`→
  generic `max/min` (identical IEEE at WP=8). exchange_nod(fct_plus,fct_minus) is a 1-rank no-op.

- **`hnode/hnode_new` at init = `zbar(nz)-zbar(nz+1)`** (node analog of helem, full cells;
  `hnode_new=hnode`, oce_ale.F90:1217), built like helem and gated. Used by the LO vertical
  update, the b2 limiter divisor, and the `use_lo` scatter reconstruction. Runtime divisors
  `areasvol`/`hnode_new` match FESOM2 (same byte-identical operands; cf. L7/L10).

- **A strong FCT gate needs ttfAB sharp+large vs the smooth ttf, or the limiter never clips.**
  The b2 increment is `flux·dt/(areasvol·hnode_new)` with `dt/(areasvol·hnode_new) ~ 1.8e-10` on
  pi, so for smooth realistic-magnitude fields the antidiffusive increment sits far below the
  admissible bound ⇒ `ae≡1` ⇒ b1/b2/b3 are computed but multiplied by 1, masking transcription
  bugs in the limiting *selection*. A high-wavenumber, amplitude-20/10 ttfAB (vs amplitude-1 ttf)
  drove clipping at ~41% of nodes (fct_plus/minus min=0), changing ~14% of the antidiffusive
  fluxes — so the sign-based clip selection was genuinely exercised. Verify clipping fraction
  post-run (count fct_plus<1); byte-identity holds regardless, but gate STRENGTH needs it.

## L12 — Assembled step (M1.4): gate the REAL driver, not an inline copy (PASSED first try)

M1.4 (the assembled advection step: `do_oce_adv_tra` + `init_tracers_AB` + the `adv_tracers_ale`
`del_ttf += advhoriz+advvert` accumulation) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run.
Unlike M1.1–M1.3 (the shim INLINED the orchestration), M1.4 drives FESOM2's OWN `init_tracers_AB`
+ `do_oce_adv_tra`. The driver-gate recipe (reuse for every future assembled-routine gate):

- **Run the real FESOM2 routine in the oracle shim AFTER the inline-kernel records are dumped.**
  The driver overwrites the shared `tracers%work` (adv_flux_*, fct_LO, fct_*, edge_up_dn_grad as
  FCT scratch) and `tracers%data(1)`. Putting the driver section LAST (just before `close(u)`)
  needs ZERO edits to the proven M1.1–M1.3 dump — purely additive, regression-safe (re-run the
  whole gate: the 34 prior fields must still PASS, which they did). The alternative (driver-first)
  works only because the inline section rebuilds work from scratch; last-is-simpler.

- **Prescribe the driver's INPUTS, not its intermediates, or the new numeric never runs.** M1.3
  prescribed `valuesAB`(ttfAB) directly; that never exercises FESOM2's AB interpolation. M1.4
  prescribes `values`(smooth ttf) + `valuesold(1)`(sharp ttfAB) so `init_tracers_AB` COMPUTES
  `valuesAB = -(0.5+ε)·valuesold + (1.5+ε)·values` (ε=0.1) and the gate compares that.

- **The AB-offset `+` did NOT trip the L7 literal-vs-runtime trap — but gate it anyway.** Feared:
  FESOM2's `epsilon` is a runtime module var (`(1.5_WP+epsilon)` computed live), so if FESOM3
  folds `1.5+0.1→1.6` at compile time the bits could differ under fast-math (as `/` does under
  `-no-prec-div`). They DON'T: `1.5d0+0.1d0` rounds to exactly the literal `1.6d0` (and
  `0.5+0.1→0.6d0`), so `valuesAB` matched byte-for-bit. `+`/`*` reassociation is far less fragile
  than `/` here. Defensive choice kept regardless: `ab_epsilon` is a non-parameter module var in
  mod_config (mirrors FESOM2), so the compiler can't fold even in principle. Lesson: a runtime
  scalar constant in a SUM is usually fold-safe, but the only proof is the oracle gate.

- **pi runs `use_wsplit=.true.`; the gate must FORCE `.false.` on both sides.** With w-split on,
  the FCT path calls `adv_tra_vert_impl` (implicit vertical, not ported until 2026-10-02 — L58) +
  recomputes the LO vertical on full `w`. Forcing `dynamics%use_wsplit=.false.` in the shim (= the
  FESOM3 `t_dyn` default) gates the matched EXPLICIT path (`w==w_e`, both passed the single
  prescribed wvel). Do not gate against the production `.true.` setting — at the time it would have
  hit the unported kernel's FESOM3 `error stop` guard. The implicit w-split path is a separate gate.

- **Non-FCT (`do_zero_flux`) + per-tracer order knobs closed in the same gate.** A second config
  on tracer 1 (MUSCL/QR4C/NON, ph=pv=0.75) exercises the `do_zero_flux=.true.` dispatch (HO scheme
  applied directly to valuesAB, scatter WITHOUT use_lo) and the tra_adv_ph/pv knobs the M1.3 gate
  fixed at 0/1. FCT vs non-FCT `del_ttf` differ by max 26 yet BOTH byte-match FESOM2 — proving the
  branch selection, not just one path. (The HANDOFF had flagged this gap for M1.4.)

- **WP pointers onto MP work arrays bind by kind VALUE.** `do_oce_adv_tra` associates `real(WP)`
  pointers (`ttf`, `fct_LO`, `adv_flux_*`, ...) with the `real(MP)` `tracers%work`/`tracers%data`
  components and passes them to the WP kernels; legal because MP==WP==8 at the DP/SP anchor (the
  pointer/argument match is on the kind integer 8, not the parameter NAME). The `tracers`/`w`/`we`
  dummies need `target`. Revisit only for FP16 (WP=2, MP=4) — a deferred precision decision.

## L13 — Pressure/EOS/N² (M2.1): split EOS + the horizontal smoother both byte-match (PASSED first try)

M2.1 (`pressure_bv`: split Jackett-McDougall EOS → `density_m_rho0`, top-down `hpressure`, N²
`bvfreq` + the horizontal `smooth_nod` sweep) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run, like
all of M1. The 10 fields (3 inputs T/S/density_ref, 3 depths, density_m_rho0, hpressure, bvfreq raw
+ smoothed) all matched. Reusable specifics:

- **`density_ref == density_0` on pi.** `use_density_ref` defaults `.false.`, so `arrays_init` sets
  `density_ref = density_0 = 1030` and `init_ref_density` is NOT called (`oce_setup_step.F90:218`).
  The plan's "subtract the `density_ref(nz,node)` ARRAY, not the scalar" is therefore bit-TRIVIAL on
  pi (the array is a constant 1030), but transcribe the array form anyway — a cavity/`use_density_ref`
  mesh (M2.11) makes it non-constant. N² divides by the SCALAR `density_0`; the PGF anomaly subtracts
  the ARRAY. Two different density_0 uses in one routine — don't conflate.

- **The horizontal `smooth_nod` (mass-matrix sweep) byte-matched by faithful transcription** — the one
  horizontally-coupled step in M2. Same reason as M1 (L9): it weights by `elem_area` (geom-proven) and
  iterates `nod_in_elem2D(:,n)` / `elem2D_nodes(:,el)` (order area-gate-proven). The `1/(3*Σarea)`
  runtime reciprocal matches because both codes divide the same byte-identical operand (L7/L10). 1-rank:
  the per-cycle `exchange_nod(bvfreq)` is a no-op, dropped (lifted at M2.12). With `N2smth_hidx=1` only
  the first sweep runs (the `do q=1,N_smooth-1` loop is empty). The smoother changed 100% of valid
  entries on the prescribed T/S, so it is a NON-vacuous gate (verify: `max|bvfreq_raw - bvfreq|>0`).

- **Dump bvfreq BOTH pre- and post-smoothing (`bvfreq_raw` + `bvfreq`) to localize.** A two-call toggle
  (`pressure_bv` with `N2smth_h=.false.` then `.true.`) separates an EOS/N²-difference bug (shows in
  bvfreq_raw) from a smoother bug (shows in bvfreq-but-not-raw). Same "gate every intermediate" rule as
  M1. density_m_rho0/hpressure are identical across the two calls (smoothing only touches bvfreq).

- **The caller pre-zeros the three outputs; the kernel does NOT.** FESOM2 `pressure_bv` writes only
  valid levels `nz=nzmin..nzmax-1` (hpressure/density) / `..nzmax` (bvfreq) and leaves the below-bottom
  entries untouched (not zeroed at allocation — heap, not `-init=zero` stack). So the dump's below-bottom
  region is indeterminate unless zeroed. BOTH the FESOM2 shim and the FESOM3 driver `=0` the outputs
  before calling → deterministic 0 there → the rectangular dump compares clean. Faithful: the FESOM3
  kernel is `intent(inout)` and also doesn't zero (matches FESOM2); the zeroing is a gate-harness concern.

- **Force the gate knobs in the shim; KPP/GM vs reduced-M2 does NOT change the EOS fields.** EOS density
  is a pure function of T/S/Z (Jackett-McDougall), so the shipped KPP/GM pi namelist gives the same
  density/hpressure/bvfreq as reduced-M2 — no need to assemble the reduced namelist for M2.1. The shim
  pins only what `pressure_bv` reads: `state_equation=1`, `which_ale='linfs'` (routes hpressure through
  this routine; zstar/zlevel compute it in `pressure_force_4_zxxxx`, M2.x), `N2smth_v=.false.`/
  `N2smth_hidx=1`, `ldiag_dMOC=.false.` (skip density_dmoc), `mix_scheme_nmb=-1` (`mixing_kpp=.false.`
  → dbsfc untouched). MLD1/2/3 + dbsfc + dMOC are KPP-only diagnostics that do NOT feed the gated fields
  (density_m_rho0 is computed BEFORE rho_surf/dbsfc1; the MLD logic only READS bvfreq/rhopot) — omitted.

- **Shared-lib rebuild: `make fesom.x` rebuilds `libfesom.so` but does NOT relink the exe — and that's
  fine.** The shim lives in `libfesom.so`; `fesom.x` resolves it at runtime (`ldd` → build/lib64). After
  editing a shim, `make -C build fesom.x` updates the `.so` (timestamp moves) while the exe stays old;
  running it loads the new `.so`. Confirm with `nm -D libfesom.so | grep <shim_symbol>`. (Same as M1.4.)

## L14 — Hydrostatic PGF (M2.2) byte-matched; and the configure.sh `--debug` Release-clobber footgun

M2.2 (`oce_pgf.F90` `pressure_force_4_linfs_fullcell`: the `gradient_sca` contraction of the M2.1
`hpressure` → element `pgf_x`/`pgf_y`) hit `max|Δ|=0` vs FESOM2 on the FIRST gate run, like all of M1
+ M2.1. The PGF gate is now the 2-field tail (pgf_x/pgf_y) of the same `tools/run_pressure_gate.sh`
(12 fields total). Specifics:

- **PGF is the SAME contraction shape as M1.1's `tracer_gradient_elements`, on a pre-gated input.**
  `pgf_{x,y}(nz,elem) = Σ_k gradient_sca({1:3,4:6},elem)·hpressure(nz,elnodes_k)/density_0`. Both
  operands were already `max|Δ|=0`: `gradient_sca` from the geometry gate, `hpressure` from the M2.1
  pass (computed in the smoothed pressure_bv call — hpressure is identical across the raw/smoothed
  calls, only bvfreq changes, L13). So faithful transcription "just works" (the L9 transitive-gate
  pattern again). `density_0` is a RUNTIME divisor but byte-identical on both sides (=1030), so the
  `-no-prec-div` reciprocal matches (L7/L10). Transcribe VERBATIM: the `/density_0` is INSIDE the sum
  (each of the 3 products divided, then summed), not `sum(...)/density_0`. The shim calls the REAL
  FESOM2 `pressure_force_4_linfs_fullcell` directly (the `pressure_force_4_linfs` dispatcher is pure
  namelist branching — no numerics — so it is not gated, cf. M1.4 which DID gate `do_oce_adv_tra`
  because that dispatch had order-knob numerics). Caller pre-zeros pgf_x/pgf_y (the kernel writes only
  `ule..nle`; FESOM2 leaves below-bottom at its allocation-zero), same as the M2.1 output treatment.
  Gate STRENGTH check: pgf is 66% non-zero at ~1e-5 m/s² (realistic PGF accel) on the prescribed T/S.

- **`configure.sh --debug` overwrote the Release anchor build dir → a spurious `area=3.7e-4` gate FAIL.**
  The footgun: `configure.sh` named the build dir `build_<compiler>_<precision>` with NO build-type
  suffix, so `--debug` reconfigured **build_intel_dp itself to CMAKE_BUILD_TYPE=Debug** and rebuilt its
  binaries at `-O0`/no-`-no-prec-div`. Running a parallel `--debug` build alongside the M1 advhor gate
  meant the gate executed the freshly-clobbered DEBUG `fesom_advhordump` → every FP field diverged,
  with `area` at the tell-tale `3.66e-4` (the exact L10 Debug-vs-Release signature). It LOOKED like an
  M2.2 regression in M1; it was a build-dir collision. **Fix:** `configure.sh` now appends `_debug` to
  the build dir for Debug builds (`build_intel_dp_debug`), so Debug never touches the Release anchor.
  **Diagnosis tell:** if a previously-`max|Δ|=0` gate suddenly shows `area≈3.7e-4` + diffs on EVERY FP
  field, you are running a Debug binary — check `CMAKE_BUILD_TYPE` in the build dir's CMakeCache.txt and
  the binary mtime/size (Debug `fesom_advhordump` ~3.1 MB vs Release ~1.75 MB) BEFORE suspecting a code
  regression. (The real M2.2 gates, run against the restored Release anchor, are all `max|Δ|=0`.)

## L15 — Partial vel_rhs (M2.3): the Debug `-check all` catches what Release + a passing gate hide

M2.3 (`compute_vel_rhs` Coriolis + AB2 + PGF + SSH-gradient, the non-advection part; momentum
advection deferred to M2.4) byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (19 fields).
But the Release gate PASSED **with a latent non-conformance** that only the Debug `-check all` run
exposed — the single most reusable lesson here:

- **`elem2D_nodes` is `(MAX_NV=4, elem2D)` in FESOM3 (quad-capable, L4), `(3, elem2D)` in FESOM2.**
  The faithful-looking transcription `elnodes = mesh%elem2D_nodes(:,elem)` (FESOM2 writes exactly
  that) assigns a 4-extent RHS into the 3-extent `elnodes(3)`. In **Release** ifort silently copies
  the LHS extent (3) — reading `elem2D_nodes(1:3)`, the correct triangle nodes — so the gate is
  `max|Δ|=0`. In **Debug** (`-check all`) it traps: `severe (408) Shape mismatch: extent of dim 1 of
  ELNODES is 3 and ... MESH is 4`. **Fix:** slice explicitly `elem2D_nodes(1:3,elem)` (what `oce_pgf`
  already did). Lesson: a green Release byte-gate does NOT prove conformance; ALWAYS run the Debug
  `-check all` pass (L10/L13/L14 said so for OOB/FPE — add *array-shape* to that list). Any FESOM2
  `(:,elem)` on `elem2D_nodes`/`gradient_sca`/other `MAX_NV`-dim arrays must become `(1:3,elem)` in
  FESOM3.

- **`coriolis` byte-matches via `r2g` transitively — a NEW geometry field gated for free (L9 again).**
  `coriolis(elem)=2·omega·sin(lat_geo)` where `lat_geo` is `r2g` applied to the rotated element
  centroid (`elem_center`). The geometry gate proved `coord_nod2D` (built by `g2r` with the rotation
  matrix), and `r2g` reuses the SAME matrix with the inverse operations — so `coriolis` is `max|Δ|=0`
  by construction, no new transcendental risk. `elem_center` uses the same literal `/3.0_WP` (L7) it
  uses for `elem_cos`. Guard the cartesian/analytic path (no `r2g` on cartesian coords → asin>1 NaN
  under `-fpe0`): for `cartesian=.true.` fill from the stored lat directly (ungated, benign).

- **Two-call gate covers BOTH AB ff branches.** `compute_vel_rhs` has a `save :: lfirst`; the first
  call uses `ff=1.0` (Euler start, `lfirst .and. .not. r_restart`), every later call `ff=ab2=1.6`.
  Driving the REAL routine TWICE in the shim (dump `uv_rhs_eul` then `uv_rhs_ab2`) exercises both —
  they differ over 66% of entries, so the ff selection is non-vacuous. The 2nd call naturally consumes
  the 1st call's `UV_rhsAB` (= Coriolis) as its "previous", a self-consistent step (no reset needed).
  FESOM3 takes `lfirst` as an explicit argument (`.true.` then `.false.`) since it has no `save` state;
  `r_restart` is folded into it (v1 has no restart — M2.11).

- **Forcing the gated config in the shim: three knobs that bite.** (1) `dynamics%ldiag_ke` **defaults
  to `.true.`** in FESOM2 `MOD_DYN` — leave it and `compute_vel_rhs` writes the `ke_*_AB` diag arrays
  (extra work, and FESOM3 has no `ke_*`); force `.false.`. (2) `momadv_opt` is `2` on pi (calls
  `momentum_adv_scalar`); force `0` so the partial gate isolates Coriolis+AB2+PGF (the FESOM3 kernel
  simply omits the momadv call). (3) `dt` is pinned to a shared constant (`1800.0_WP`) in BOTH codes,
  **not** pi's namelist `dt=86400/36=2400 s` — the gate only needs both sides equal in the final
  `dt*(…)/elem_area` scaling, and pinning removes the namelist dependency (L9).

- **The `use_pice=0` path lets a minimal fake `ice` satisfy the real routine.** `compute_vel_rhs`
  unconditionally associates `m_ice => ice%data(2)%values` / `m_snow => ice%data(3)%values` at entry,
  but with `which_ale='linfs'` → `use_pice=0` they are never read. `ice_setup` runs AFTER `ocean_setup`
  (where the shim fires), so no live `ice` exists yet — the shim builds a throwaway `type(t_ice)` with
  just `data(2:3)%values` allocated. The pointer targets exist; the values are irrelevant. (FESOM3 v1
  has no ice at all, so its kernel drops the `p_ice` term entirely — byte-identical since p_ice=0.)

- **Same transitive-gate reason M2.3 "just worked":** every operand was already byte-pinned —
  `coriolis` (r2g/geometry), `elem_area`/`gradient_sca` (geometry gate), `pgf_x`/`pgf_y` (M2.2),
  `ab_epsilon`-derived `ab1=-0.6`/`ab2=1.6` (L12 fold-safe), and the runtime `/elem_area` divisor is
  byte-identical on both sides (L7/L10/L14). The only genuinely new arithmetic is the AB blend and the
  SSH-gradient `sum(gradient_sca·(-g·eta))`, both simple `+`/`*`/`sum` — fold-safe.

## L16 — Momentum advection (M2.4): the full UV_rhs gate; OpenMP-off oracle; gate the intermediate, not the perturbation (PASSED first try)

M2.4 (`momentum_adv_scalar` → the FULL `UV_rhs`: Coriolis + AB2 + PGF + SSH-grad + momentum advection)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (21 fields), like all of M1 + M2.1–M2.3.
`momentum_adv_scalar` lives in `src/oce/oce_dyn_velrhs.F90` (private, called by `compute_vel_rhs` when
`momadv_opt==2` — mirroring FESOM2's own `oce_ale_vel_rhs.F90` layout, NOT the plan's separate
`oce_dyn_momadv.F90`). Reusable specifics:

- **The FESOM2 oracle is built `ENABLE_OPENMP=OFF`, so momadv's `omp_set_lock(partit%plock)` +
  `!$OMP ORDERED` paths compile OUT.** momadv is the FIRST gated kernel with OpenMP locks/ordered
  reductions (the M1/M2.3 kernels had `!$OMP DO` but no locks). With OpenMP off they are dead code and
  the per-node edge-scatter runs SERIALLY in natural edge order (1..edge2D) — which is exactly FESOM3's
  serial loop. So the order-sensitive accumulation matches by construction (the L8/L9 argument: edges
  order is geometry-gate-proven, `nod_in_elem2D` order is area-gate-proven). Check `ENABLE_OPENMP` in
  the oracle `build/CMakeCache.txt` before trusting a parallel-reduction kernel's serial transcription;
  if it were ON with `__openmp_reproducible`, the `!$OMP ORDERED` would still serialize, but plain
  OpenMP would reorder → ULP diffs (re-gate would need the reproducible build or 1-thread).

- **Gate the operator's own intermediate, not just its (small) contribution to the dominant field.**
  momadv adds only ~`elem_area·avg(uvnode_rhs)` ≈ 8 to a Coriolis `UV_rhsAB` of ≈1e4 (≈0.07%), and ~1%
  of the final `uv_rhs`. A small perturbation to a dominant field is a WEAK gate *if you only watch the
  dominant field*. Fix: dump `uvnode_rhs` (the post-`areasvol_inv` nodal advection) as its OWN record —
  it is 68.9% non-zero with BOTH signs (the edge-scatter `+nod(1)` / `−nod(2)` branches), gated at its
  native ~1e-6 magnitude. So a bug in EITHER the vertical (`w·du/dz`) OR horizontal (`u·du/dx`) pass
  shows in `uvnode_rhs` directly, independent of how little it moves `uv_rhs`. (The EXACT `max|Δ|=0`
  gate would catch any ULP diff regardless — but the separate intermediate LOCALISES it, the L9/L13
  "gate every intermediate" rule.) And because the gate runs the REAL FESOM2 vertical+horizontal passes,
  FESOM3 matching `uvnode_rhs` proves BOTH passes byte-identical *by construction* (a no-op vertical
  pass would leave `uvnode_rhs` horizontal-only ≠ FESOM2 → FAIL).

- **Prescribe the velocity fully-defined on `1:nl-1` to keep the `0·below-bottom` products clean.** The
  horizontal scatter reads `un2(nz)*UV(:,nz,el2)` for `nz` up to `max(nl1,nl2)` — i.e. BELOW el2's
  bottom, where `un2(nz)=0` (zero-filled) but `UV(:,nz,el2)` is whatever sits there. `0·finite=0`
  (clean) but `0·NaN=NaN`. The driver/shim prescribe `UV` for ALL `nz=1..nl-1` (every element, every
  level), so the below-bottom reads are finite zeros-of-the-analytic-formula, never uninitialised heap →
  no NaN under `-fpe0`. (Same defensive-prescription reason M1.2 set `wvel` on the full column.)

- **`w_e` prescribed like M1.2's `wvel`** — `1e-4·sin(2·lon)·cos(lat)·cos(0.3·nz)`, sign varying in
  space AND depth so the `w·du/dz` finite difference `wu(nz)−wu(nz+1)` is non-trivial. Dumped as an
  input record (`w_e`) so the prescription itself is gated, like `wvel`/`uv_in`.

- **The momadv divisors follow the established rules.** `/(3._WP·hnode(nz,n))` is a RUNTIME divisor but
  byte-identical operands on both sides → reciprocal matches (L7/L10). The vertice→element `/3.0_WP` is
  a LITERAL divisor (exact 1/3 fold, byte-safe; triangles only — the L7 arity caveat). `elem2D_nodes`
  accessed PER-COMPONENT (`(1,el)`/`(2,el)`/`(3,el)`, scalar) not as a `(:,el)` slice → no MAX_NV=4-vs-3
  shape mismatch (the L15 trap; Debug `-check all` confirmed clean).

## L17 — Biharmonic viscosity (M2.4, opt_visc=7): no new geometry; large |du| to exercise the flow-aware max (PASSED first try)

M2.4 biharmonic viscosity (`viscosity_filter(7)` → `visc_filt_bidiff`, FESOM2 `oce_dyn.F90:591-744`)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (24 fields now), like all of M1 + M2.1-M2.4.
It lives in `src/oce/oce_dyn_visc.F90` (mirroring FESOM2's `oce_dyn.F90` — a SEPARATE file/operator run
AFTER `compute_vel_rhs`, NOT inside it like momadv). Reusable specifics:

- **opt_visc=7 needs NO new geometry — the HANDOFF's `gradient_vec` worry was unfounded for THIS scheme.**
  `visc_filt_bidiff` is a pure edge-based ∇² applied TWICE (a biharmonic): pass 1 scatters the across-edge
  velocity jump `u1=UV(el1)-UV(el2)` × a flow-aware coeff into the element field `U_c`/`V_c`; pass 2
  scatters the across-edge jump of `U_c` into `UV_rhs/elem_area`. It reads only `edge_tri`, `elem_area`,
  `ulevels`/`nlevels`, `edge2D_in` (+ `dynamics%uv`/`uv_rhs`/`work%u_c`/`v_c`) — every one already gated.
  `gradient_vec` is still deferred but is NOT required here (it would be for a vector-Laplacian scheme; the
  edge-difference biharmonic sidesteps it). So the L9 transitive-gate pattern again → `max|Δ|=0` first run.

- **First gated kernel to use `edge2D_in` (interior-edge-only loop; free slip on boundary edges).** FESOM2
  `if(myList_edge2D(ed)>edge2D_in) cycle`. At 1-rank `myList_edge2D` is identity AND the mesh orders interior
  edges first (`1..edge2D_in`, boundary `edge2D_in+1..edge2D` — `fvom_init` writes them so, line 554-555), so
  the lift is `do ed=1,edge2D; if(ed>edge2D_in) cycle`. `edge2D_in` is byte-pinned transitively: same
  `edgenum.out` (FESOM3 reads line 2), same edge ordering proven by the geometry gate (`edge_tri`/
  `edge_cross_dxdy` were `max|Δ|=0`). For interior edges BOTH `edge_tri(1/2,ed)>0`, so no `el2<0` guard is
  needed (unlike momadv, which looped ALL edges). `edge_tri(:,ed)` is the size-2 left/right slot — NOT a
  MAX_NV dim, so no L15 shape trap, and `elem2D_nodes` is never touched. Debug `-check all` EXIT 0.

- **The biharmonic is TWO SEPARATE edge loops with `U_c`/`V_c` between — the serial structure replaces the
  `!$OMP BARRIER`.** Pass 1 must FULLY fill `U_c` before pass 2 reads `U_c(el1)-U_c(el2)`; FESOM2 enforces
  this with an `!$OMP BARRIER` between the loops, FESOM3's two sequential `do ed` loops get it for free.
  `exchange_elem(U_c/V_c)` between them is a 1-rank no-op (lifted at M2.12). `U_c`/`V_c` are zeroed over all
  elements first. The OpenMP-off oracle (L16, re-confirmed `ENABLE_OPENMP=OFF` in the oracle CMakeCache)
  compiles out the `omp_set_lock`/`!$OMP ORDERED` scatter → serial edge order → matches FESOM3's serial loop.

- **Gate STRENGTH: the flow-aware coeff `max(γ0, max(γ1·|du|, γ2·|du|²))` needs LARGE |du|, or only γ0 fires.**
  The pi gammas are `γ0=0.003` (background; OVERRIDES the `t_dyn` default 0.03 — pin the namelist value!),
  `γ1=0.1`, `γ2=0.285`, `γ0_h=γ1_h=0` (pure biharmonic). The winners cross at `|du|=γ0/γ1=0.03` (γ0→γ1) and
  `|du|=γ1/γ2=0.351` (γ1→γ2). At pi-realistic velocities (~0.1 m/s) `|du|≲0.02` so ONLY γ0 (background) is
  ever selected → the γ1/γ2 branches are computed-but-not-selected, masking a transcription bug there (the
  L11 weak-gate trap). Fix: bump the SHARED prescribed `UV` amplitude to 2.0/1.5 m/s (a strong-current
  stress test; the operator gate has no stability constraint) so `|du|` spans all three — measured
  selection **γ0/γ1/γ2 = 20.4/78.6/1.0%** of edge-levels (γ2 on ~2000 levels, a clear margin). A
  driver-side, NON-gated diagnostic counts the actual `max` winner to verify this. NOTE: γ2 is a near-dead
  path in PRODUCTION pi (the namelist even labels γ2 "only used for opt_visc=5/8") — it fires only on
  synthetic |du|>0.351; the gate exercises it deliberately. Bumping UV re-gated M2.3/M2.4 with new values
  (automatic — identical formula on both sides; all prior fields stayed `max|Δ|=0`).

- **Gate the operator's own intermediate (the L16 rule).** Dump `visc_u_c`/`visc_v_c` (the first-stage
  Laplacian, the element field pass 1 builds) as their OWN records, so a pass-1 bug localises there while a
  pass-2 bug shows in `uv_rhs_visc`-but-not-`visc_u_c`. The post-viscosity `uv_rhs_visc` is the gate target;
  viscosity only READS `UV` and ADDS its increment into the incoming `uv_rhs` (= `uv_rhs_ab2`, already gated),
  so running it after the 2nd `compute_vel_rhs` is well-defined.

- **The viLapl (harmonic-addition) term is transcribed but contributes 0 on pi (γ_h=0) — a deferred sub-gate.**
  `viLapl=dt·max(γ0_h, γ1_h·|du|)·len`; with γ0_h=γ1_h=0 it is exactly 0 (no FPE: `0·sqrt(...)`=0). The
  γ_h>0 combined harmonic+biharmonic path is ungated on pi (like the cavity/zstar/use_density_ref deferrals);
  re-gate it when a config with γ_h>0 is exercised. All divisors follow the established rules: `len=
  sqrt(sum(elem_area))`, `-dt·sqrt(...)`, final `/elem_area(el)` are RUNTIME divisors byte-identical on both
  sides (L7/L10/L14). Local intermediates are `real(WP)`; FESOM2 hard-codes them `real(kind=8)` — byte-equal
  at the DP anchor (WP=8), a deferred SP-precision nuance.

## L18 — Implicit vertical viscosity (M2.5): a sequential TDMA byte-matches by pure transitivity (PASSED first try)

M2.5 implicit vertical viscosity (`impl_vert_visc_ale`, FESOM2 `oce_ale.F90:3315-3526`) byte-matched FESOM2
`max|Δ|=0` on the FIRST Release gate run (28 fields now), like all of M1 + M2.1-M2.4. It lives in a new
`src/oce/oce_dyn_ivertvisc.F90`, a SEPARATE operator run AFTER `viscosity_filter` (FESOM2 `oce_ale.F90:3874`).
Reusable specifics:

- **A tridiagonal Thomas solve is a strictly SEQUENTIAL recurrence — there is NO summation/scatter order
  ambiguity to defend (unlike the edge/node scatters of M2.4/M2.4-visc).** Per element the column is
  independent; forward elimination then back substitution are recurrences with one well-defined order. So the
  gate reduces ENTIRELY to "are the operands byte-identical?" — and they are: `UV`, the incoming `uv_rhs` (=
  the M2.4-visc-gated `uv_rhs_visc`), the prescribed `w_i`/`Av`/`stress_surf`, the geometry-proven
  `helem`/`zbar_e_bot`/`ulevels`/`nlevels`, and `C_d`/`density_0`/`dt`. The L9 transitive-gate pattern again →
  `max|Δ|=0` first run. The many runtime divisors (`zinv`, the Z/zbar differences, the Thomas `1/b`, `1/m`)
  are byte-identical operands on both sides, so the `-no-prec-div` reciprocal matches (L7/L10/L14).

- **Inputs not yet produced by the core (`Av` from PP mixing M2.8, `stress_surf` from forcing M2.10) → pass
  as EXPLICIT dummy args, not as type fields.** Honest about provenance, minimal diff, defers the
  type-placement decision; when M2.8/M2.10 land the caller just sources them (`dyn%work%Av`, `forcing%...`) and
  the kernel is unchanged. FESOM2 pulls both from `o_ARRAYS`; the shim sets those module arrays directly.
  Prescribe `Av` STRICTLY POSITIVE (a real viscosity: `[1e-3,1.4e-2]` m²/s here) over all `nz=1..nl` (the TDMA
  reads `Av(nzmax)`), and `stress_surf` sign-varying (~0.1 N/m², both signs → the surface BC is exercised
  fully).

- **New mesh field `zbar_e_bot` + the (missing) `helem` — both full-cell trivial; serialization is
  unallocated-safe.** The pressure driver built `hnode` but NOT `helem`; M2.5 needs both. Full cells (pi,
  `use_partial_cell=.false.`): `helem(nz,e)=zbar(nz)-zbar(nz+1)` (the element analog of `hnode`),
  `zbar_e_bot(e)=zbar(nlevels(e))` (FESOM2 `init_bottom_elem_thickness`). The kernel rebuilds the per-column
  `zbar_n`/`Z_n` bottom-up from these. Added `zbar_e_bot` to `t_mesh` + write/read serialization;
  `write_bin_array` writes `0` for an unallocated array and `read_bin_array` leaves it unallocated, so the
  ctest mesh round-trip (where `zbar_e_bot` is never built) stays green — 13/13.

- **The single-layer-column OOB is a REAL benign-OOB in FESOM2 that pi avoids — CHECK THE MESH (`elvls.out`)
  before assuming `-check all` clean.** For `nlevels(elem)==2` the "last row" reads `Z_n(nzmax-2)=Z_n(0)` and
  `UV(:,nzmax-2)=UV(:,0)` — values that are overwritten by the "first row" block / multiplied by `a(top)=0`,
  so the Release result is byte-correct (FESOM2's Release oracle tolerates it), but `-check all` would TRAP.
  pi has NO such columns (min element `nlevels=5` per `elvls.out` — i.e. ≥4 layers), so the EXACT transcription
  is both byte-identical AND `-check all` clean (EXIT 0, confirmed). A `nlevels==2` shelf column (CORE2) needs
  a guard matched to FESOM2's behaviour — deferred to M2.11. Determining this up front (read `elvls.out`)
  avoided a guess.

- **`w_i` (implicit vertical velocity) prescribed sign-varying to exercise both upwind branches; the TDMA must
  be shown NON-vacuous.** A DISTINCT analytic formula from `w_e` (`sin(lon)·cos(2·lat)·cos(0.25·nz)`) so its
  sign varies in space AND depth → `min(0,wu)`/`max(0,wu)` (+ the `wd` pair) all fire. Driver diagnostic
  (non-gated): `w_i>0/<0 = 75360/75360` (both signs) and `max|d(uv_rhs)|=0.985` (the solve substantially
  transforms the rhs, not a near-identity) — the L11 weak-gate guard.

- **No auto-generated `*_interface` module for `impl_vert_visc_ale` (only the `_vtransp` variant) → declare an
  explicit interface block in the shim.** FESOM2's interface-gen tool skips it because it is declared in an
  explicit interface block inside `oce_ale.F90:269`; the other gated routines (`pressure_bv`,
  `compute_vel_rhs`, `visc_filt_bidiff`) DO get `<name>_interface.mod`. So the shim copies the interface block
  verbatim. Pin `C_d=0.0025` in the shim (the routine reads it from `o_PARAM`; matches the FESOM3
  `mod_param_phys` default).

- **Gate target is `UV_rhs`, NOT `UV`.** The routine OVERWRITES `UV_rhs` with the Thomas solution; `UV` is
  read-only (the velocity update `UV += UV_rhs` happens later in the timestep). The HANDOFF "post-solve UV"
  note was imprecise — the dumped gate target is `uv_rhs_ivv`. Out-of-pi branches dropped: `ldiag_ke`
  `ke_wind`/`ke_drag` (no `ke_*` in FESOM3 `t_dyn`, as `compute_vel_rhs`), and the `toy_ocean`
  dbgyre/neverworld2/soufflet constant-`C_d` friction (pi `toy_ocean=.false.` → quadratic bottom drag
  `-C_d·|UV_bottom|`).

## L19 — SSH stiffness + CG solve (M2.6): the first iterative solver byte-matches; the incremental-build ULP-drift footgun (PASSED after a clean rebuild)

M2.6 SSH (CSR stiffness `init_stiff_mat_ale` + `compute_ssh_rhs_ale` + preconditioned-CG `solve_ssh_ale`,
FESOM2 `oce_ale.F90:1584/2012/3272` + `solver.F90`) byte-matched FESOM2 `max|Δ|=0` on all 4 new fields
(`ssh_stiff_diag`, `ssh_Aeta`, `ssh_rhs`, `d_eta` → **32 fields** total). Built `src/oce/oce_ssh_rhs.F90`
(stiffness + rhs) + `src/oce/oce_ssh_solve.F90` (precond + CG). The FIRST iterative solver in the port.

- **THE FOOTGUN (cost me one red gate): after ADDING NEW module files, an incremental FESOM3 build is
  byte-UNRELIABLE — clean-rebuild before gating.** The first gate run FAILED with a *uniform ~few-ULP
  relative* delta on EVERY field — including prescribed trig-derived inputs (`temp` 7e-15, `coriolis` 5e-20)
  while integer/geometry fields (`zbar`, `hnode`) stayed exactly 0. That split is the L10 "uniform flag/
  codegen effect" signature, NOT a kernel bug. Root cause: CMake's `GLOB CONFIGURE_DEPENDS` detected the 2
  new `oce_ssh_*.F90`, re-ran `cmake` configure mid-incremental-build, and the partial recompile linked
  objects built against a MIX of stale/fresh `.mod` interfaces → ULP-drifted FP codegen vs the FESOM2
  oracle. `./configure.sh --clean --build` (full recompile, consistent `.mod` set) restored `max|Δ|=0` on
  all 32 fields on the next run. **Rule: when a milestone adds NEW source files (not just edits), do a CLEAN
  Release rebuild before trusting the byte-gate.** (Edits-only incremental builds stayed byte-exact through
  M2.1–M2.5; it is specifically NEW files + the reconfigure that drift.)

- **A CG solve byte-matches by the same L9 transitivity as the M2.5 TDMA — given byte-identical operands the
  iteration is a deterministic recurrence.** Same matrix + rhs + x0 ⇒ identical per-iteration scalars ⇒ the
  convergence test (`sqrt(Σrr²/nod2D) < soltol·sqrt(Σb²/nod2D)`, soltol=1e-5) fires at the SAME iteration
  (here **37**) ⇒ identical `d_eta`. So a byte-matching `d_eta` proves the entire 37-iteration recurrence
  (matvecs + dot-products + preconditioner) is bit-identical. The operands are all pre-gated: `ssh_stiff`
  (geometry-proven `edge_tri`/`edges`/`gradient_sca`/`edge_cross_dxdy`/`zbar_e_bot`/`areasvol` + the edge-
  scatter accumulation order), `ssh_rhs` (same edge order + the M2.5-gated `uv_rhs`/`helem`), `g`/`dt`/`α`/`θ`.

- **REDUCTION FORM is part of the bits.** The FESOM2 oracle is `ENABLE_OPENMP=OFF` and `__openmp_reproducible`
  undefined, so its dot-products compile to the explicit serial `DO row; s=s+…; END DO` (the `!$OMP REDUCTION`
  is an inert comment). Transcribe THAT serial form for `s_old`/`s_aux`/`sprod` — NOT a `sum()` intrinsic
  (the compiler may reduce it in a different order). The matrix-vector products DO use `sum()` over a CSR
  slice, exactly as FESOM2 — both codes `sum()` the same slice ⇒ byte-identical. (Same OpenMP-off reasoning
  as L16/L17's serial edge scatter.)

- **dt for the stiffness is the pi NAMELIST dt (=86400/36=2400 s), NOT the shim's `dt=1800`, and must be
  computed the SAME way (not a `2400.0` literal).** `init_stiff_mat_ale` runs at `ocean_setup` (oce_setup_step
  :140) where dt is still the namelist value, BEFORE the shim overrides it to 1800 for the vel_rhs gates. dt
  enters M2.6 ONLY via the stiffness (`factor=g·dt·α·θ`, mass `areasvol/dt`); compute_ssh_rhs has no dt. So
  the FESOM3 driver builds the stiffness with `dt_ssh = 86400._WP/real(36,WP)` (FESOM2's `gen_model_setup`
  formula) — a literal `2400.0_WP` could differ by an ULP under `-no-prec-div` (the L7 trap) and break the
  whole matrix.

- **`α=θ=1.0` on pi (FESOM2 default, namelist doesn't set them) ⇒ two simplifications:** the linfs water-flux
  term `(1-α)·ssh_rhs_old` VANISHES (so `ssh_rhs_old` is multiplied by 0 — set it to 0, don't bother
  prescribing it), and `factor = g·dt`. Force `α=θ=1.0` in BOTH shim and driver defensively (the stiffness
  already used 1.0; compute_ssh_rhs re-reads them).

- **linfs builds the stiffness ONCE and never updates it** (`oce_ale.F90:3921` calls `update_stiff_mat_ale`
  only for NON-linfs). So the single `init_stiff_mat_ale` assembly with the unperturbed depth IS the
  production matrix on pi/reduced-M2 — no `update_stiff_mat_ale` port needed for M2.

- **1-rank collapses the global-numbering machinery to the identity.** The CG uses `colind_loc`/`rowptr_loc`
  (LOCAL CSR); the FESOM2 global remap (`rpart.out` mapping, the per-PE nza offset, `myList_nod2D`) is
  identity at 1-rank and is dropped — build the local CSR and set `colind=colind_loc`, `rowptr=rowptr_loc`.
  `exchange_nod`(diag/rr/pp/x) and the `MPI_Allreduce`(s_old/s_aux/sprod) are no-ops (the local serial sum
  IS the global sum), dropped. The new multi-rank bit-identity risk — the cross-rank dot-product reduction
  ORDER — is an M2.12 concern.

- **Gate the matrix with a MATVEC, not just the diagonal.** `ssh_stiff_diag` (the row-diagonal = mass +
  self-stiffness) localises the assembly, but `ssh_Aeta = A·eta_n` (full CSR matvec against the prescribed
  `eta_n`) exercises EVERY non-zero — the off-diagonals are the bulk of the operator. Both `max|Δ|=0`.

- **The `Σ ssh_rhs` telescoping check is ~1e-13 RELATIVE, not absolute.** `ssh_rhs` is an edge-divergence
  (`+/-（c1+c2)` scattered into the two edge nodes) so `Σ_nodes` telescopes to 0 in exact arithmetic; the FP
  roundoff scales with the `ssh_rhs` magnitude. With UV bumped to 2.0/1.5 m/s (the M2.4-visc strength test)
  `ssh_rhs ~ 1e8`, so `Σ = -3.7e-5` (≈1e-13 relative) — consistent, non-vacuous. The plan's "~1e-13" was the
  small-UV estimate.

- **The prescribe-and-stop shim STILL works for M2.6** (the worry it might need to "run past forcing" was
  unfounded): `init_stiff_mat_ale` runs at `ocean_setup` line 140, BEFORE the end-of-`ocean_setup` shim, so
  the matrix is already assembled; the shim just prescribes `UV`/post-TDMA `UV_rhs`/`d_eta`=0/`ssh_rhs_old`=0
  and drives the REAL `compute_ssh_rhs_ale` + `solve_ssh_ale`. The reduced-M2 namelist is still NOT needed.
  No auto-gen `*_interface` module for `compute_ssh_rhs_ale`/`solve_ssh_ale` (internal interface blocks) →
  declare explicit interfaces in the shim (as for `impl_vert_visc_ale`, L18).

- **`zbar_e_srf` = `zbar(ulevels(elem))` = `zbar(1)` = 0 on pi (no cavity).** init_stiff's H factor is
  `(zbar_e_bot − zbar_e_srf)`; computed inline from `zbar`+`ulevels` (no new mesh field — faithful to the
  FESOM2 non-cavity default `oce_ale.F90:525`). The CG `nod2D` divisor (`rtol`, exit test) is the global
  `nod2D=3140`. Debug `-check all` clean (RC=0): `n_pos(12,nod2D)` is wide enough (pi max node degree < 11),
  no CSR-slice OOB.

## L20 — ALE velocity/SSH/thickness-W update (M2.7): the post-CG tail byte-matches; in-place overwrite of prescribed inputs; the Debug I/O stack-temp overflow (PASSED first try)

M2.7 (the linfs post-CG tail: `update_vel` + `compute_hbar_ale` + the `eta_n` blend + `vert_vel_ale`
→ `compute_CFLz` + `compute_Wvel_split`) byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run
(all 11 new fields → **43 fields**), like all of M1 + M2.1-M2.6. Built `src/oce/oce_ale.F90` (FESOM3
consolidates `update_vel` — FESOM2 keeps it in `oce_dyn.F90` — with the `oce_ale.F90` ALE routines into
one M2.7 module). Reusable specifics:

- **The whole post-CG chain is `max|Δ|=0` by pure L9 transitivity — every operand was already pinned.**
  `update_vel` is a `gradient_sca` contraction of `-g·θ·dt·d_eta` (d_eta from the M2.6 CG, gradient_sca
  geometry-gated) added to the post-TDMA `UV_rhs` (M2.5; the CG solve never touches it) into `UV`. The
  `compute_hbar_ale` + `vert_vel_ale` edge-divergences reuse the SAME edge order + `helem`/`edge_cross_dxdy`/
  `areasvol`/`area` the M2.6 `compute_ssh_rhs_ale` already gated. So nothing genuinely new arithmetic-wise →
  faithful transcription "just works" (the recurring L9 pattern). `hbar` is prescribed identically on both
  sides. `dt=1800`/`θ=1`/`α=1` pinned (as M2.3-M2.6); the OpenMP-off oracle (L16/L17/L19) serial edge order
  matches the serial loops.

- **Three prescribed inputs are OVERWRITTEN in place by the M2.7 kernels — save copies BEFORE the chain or
  the INPUT dump records become wrong.** `update_vel` overwrites `UV` (but `uv_in` was already saved at
  M2.3); the `eta_n` blend overwrites `dynamics%eta_n`; `compute_Wvel_split` (inside `vert_vel_ale`)
  overwrites `dynamics%w_e`/`w_i` — the SAME arrays prescribed as the M2.4/M2.5 inputs. Since the dump fires
  at the END (after all compute), the `eta_n`/`w_e`/`w_i` INPUT records would otherwise dump the M2.7 output.
  Fix: save `eta_n_in`/`w_e_in`/`w_i_in` right before the chain (they are fully consumed by M2.4/M2.5 first),
  dump those for the input records, and dump the post-M2.7 values as NEW records (`eta_n_upd`/`w_split_e`/
  `w_split_i`). This is correct physics: in the real timestep `vert_vel_ale` recomputes `w_e`/`w_i` for the
  NEXT step's momadv/ivertvisc. Same save-the-input pattern as `uv_in` (M2.3) and the M2.1 raw/smoothed split.

- **`compute_CFLz` keeps its TWO-statement form for byte-reproducibility (the L19 "reduction form is part of
  the bits" rule, made explicit by FESOM2).** `CFL_z(nz)=CFL_z(nz)+c1` then `CFL_z(nz+1)=c2` — FESOM2's own
  comment says this exact split (vs. folding both into one accumulate) is "for the sake of reproducibility …
  (rounding error)". Transcribe it verbatim; do NOT combine. `c1`/`c2` are scalars here (a DIFFERENT `c1`
  from the `vert_vel_ale` per-level array).

- **linfs collapses `vert_vel_ale` to W-only.** `which_ale='linfs'` ⇒ the `zlevel`/`zstar` thickness-
  redistribution branches are NOT taken, so `hnode_new` stays = `hnode` (its init value — gate it anyway to
  confirm the linfs path leaves it untouched) and there is no surface Wvel/hnode correction; `compute_hbar_ale`'s
  water-flux term (`.not. linfs`) vanishes. So `vert_vel_ale` is just: zero W → edge-scatter `div(UV·h)` →
  cumsum bottom-up → `/area`. Fer_GM (`fer_UV`/`fer_Wvel`) and ldiag_ke (`ke_*`) branches dropped (no GM / no
  `ke_*` in v1, as `compute_vel_rhs`/`impl_vert_visc_ale`). The full-free-surface thickness evolution gets its
  own later gate.

- **Non-vacuity needs `use_wsplit=.true.` (pi production) + a large UV — else `compute_Wvel_split` is the
  trivial `Wvel_e=Wvel` branch.** The split only deviates where `CFL_z > wsplit_maxcfl` (=1.0 on pi). With the
  M2.4-visc-bumped UV (2.0/1.5 m/s) the W divergence is large enough that `CFL_z>1` on **13253** (nz,node) →
  the `dd`/`Wvel_e`/`Wvel_i` formula is genuinely exercised (a driver diagnostic counts it; the L11/L17 weak-
  gate guard). Force `use_wsplit=.true.`/`wsplit_maxcfl=1.0` in BOTH shim and driver (the namelist sets them,
  but assert — defensive like the M2.6 α=θ=1). `max|w|=0.042` m/s (physical), `max|uv_upd|=2.5` (the SSH-grad
  correction + UV_rhs visibly moved UV) — all non-vacuous.

- **The "`K_v⁻` deformation bound" plan bullet was a MISLABEL — there is no Kv bound in the post-CG ALE path.**
  `Kv` (vertical diffusivity) is produced by PP mixing (M2.8), not the velocity/SSH/thickness update. Grep of
  the `oce_ale.F90` step tail (`update_vel`→`compute_hbar_ale`→`eta_n`→`vert_vel_ale`) confirmed no Kv/
  deformation bound there. Folded the bullet into M2.8.

- **THE DEBUG FOOTGUN (not a kernel bug): the I/O dump writer SEGFAULTs under Debug `-check all` on the default
  8 MB stack — `ulimit -s unlimited` fixes it.** The compute ran CLEAN under `-check all` (every M2.7 diagnostic
  printed → all six kernels completed, no OOB/shape/FPE; the L15 `elem2D_nodes(1:3,·)` slices are correct). The
  crash was at `mod_advhor_dump.F90:68` `write(u) real(a, real64)` on the FIRST big `wr_r3` (`uv_in`, an
  EXISTING record) — ifort `-O0` puts the ~4.4 MB `real(...,MP)` array temporary on the stack, overflowing the
  8192 KB default (the extra M2.7 records pushed cumulative pressure over the edge). `ulimit -s unlimited` →
  EXIT 0, dump written. Release (`-O3`) is unaffected (different temp handling). **Diagnosis tell:** a segfault
  in `wr_r{2,3}`/`real(a,real64)` AFTER all compute diagnostics printed is a stack-temp overflow, NOT a kernel
  fault — raise the stack, don't hunt the kernel. (Byte-gate is Release-vs-Release; Debug is OOB-only, L10.)

## L21 — PP vertical mixing (M2.8): three sequential passes byte-match by transitivity; prescribe `uvnode` (not `UV`) to isolate the add + control the Ri-factor range (PASSED first gate run)

M2.8 (Pacanowski-Philander Richardson-number mixing: `oce_mixing_pp` → `Kv` on nodes, `Av` on elements)
byte-matched FESOM2 `max|Δ|=0` on the FIRST Release gate run (3 new fields `uvnode`/`pp_Kv`/`pp_Av` →
**46 fields**), like all of M1 + M2.1-M2.7. Built `src/oce/oce_ale_mixing_pp.F90` (`oce_mixing_pp` +
`Kv0_background_qiang` + `Kv0_background`, mirroring the oracle filename). Reusable specifics:

- **Three SEQUENTIAL passes, byte-match by pure L9 transitivity.** Pass 1 (nodes): `Kv := factor =
  shear/(shear + 5·max(N²,0) + 1e-14)`, the inverse-Richardson factor `1/(1+5·Ri)` with `shear =
  |d(uvnode)/dz|²`. Pass 2 (elements): `Av = mix_coeff_PP·mean₃(factor²) + A_ver`. Pass 3 (nodes): `Kv =
  mix_coeff_PP·factor³ + K_ver`. The ordering is LOAD-BEARING — `Kv` is the scratch that holds the factor
  (pass 1) AND the final diffusivity (pass 3), so the `Av` elem loop MUST run between them (Av reads factor²,
  Kv overwrites with factor³). No scatter/accumulation (`sum(Kv(nz,elnodes)**2)` is a fixed 3-element array
  order, elnodes geom-gated) → no order ambiguity → every operand pinned (`bvfreq` M2.1, `uvnode`/`Z_3d_n`
  geometry) ⇒ `max|Δ|=0` first run. The recurring L9 pattern (like M2.5's TDMA).

- **Prescribe `dyn%uvnode` DIRECTLY (a new input), do NOT compute it from the prescribed `UV`.** PP reads the
  nodal velocity `dynamics%uvnode` (the area-weighted elem→node average `compute_vel_nodes` fills each step,
  `oce_dyn.F90:177`). Two wins from prescribing it rather than running `compute_vel_nodes`: (1) **isolation** —
  the existing `UV` (and hence all M2.3-M2.7 records) is UNTOUCHED, so M2.8 is a pure additive gate (no
  re-gate churn); (2) **non-vacuity control** — the PP factor depends on the VERTICAL shear `Δuvnode/Δz`, and
  the existing `UV`'s depth term (`-0.005·nz`, tuned so it cancels in the across-EDGE viscosity `du`) gives a
  shear ~3e-8 → factor ~1e-3 (technically non-zero in dp, but a weak gate). A bespoke `uvnode` with a strong
  vertical shear (`1.2·cos(lat)sin(lon)·sin(0.5·nz)` etc., per-layer jump ~0.6 m/s) makes the factor span
  **[5.8e-7, 0.93]** (62.4% of node-levels > 0.1) → `f²`/`f³` exercised across their full range, `max|Kv|`=8e-3
  (800× the `K_ver`=1e-5 background), `max|Av|`=8.7e-3 (87× `A_ver`). This is the M2.5 precedent (prescribe the
  not-yet-sourced input: `Av`/`stress_surf` there): `compute_vel_nodes` is gated with the step at M2.9, where
  `impl_vert_visc_ale` will also source its `Av` from `dyn%work%Av` instead of the prescribed analytic one.

- **`A_ver` is a pi NAMELIST override (1e-4), NOT the `mod_param_phys` default (1e-3) — force it.** The pi
  `namelist.oce` sets `A_ver=1.e-4`; FESOM3's `mod_param_phys` default is `0.001`. PP's `Av=…+A_ver` reads it,
  so the driver+shim must FORCE `A_ver=1e-4` (+ `mix_coeff_PP=0.01`, `K_ver=1e-5`, `Kv0_const=.true.`, all pi
  defaults). Proof `A_ver` was unused by M2.1-M2.7: M2.5 gated `max|Δ|=0` with the oracle at `A_ver=1e-4` and
  FESOM3 at its 1e-3 default — if `impl_vert_visc_ale` read `A_ver` they would have diverged. (`Kv0_const=.false.`
  Qiang lat/depth background + the cavity `nzmin>1` path are transcribed but UNGATED on pi — M2.11.)

- **PP OVERWRITES the shared `Av`; save the M2.5 prescribed `Av` first (the L20 save-the-input pattern).** In
  the FESOM2 shim `Av`/`Kv` are o_ARRAYS globals: PP clobbers the SAME `Av` the M2.5 section prescribed, and the
  dump fires at the END → the `Av` record would dump PP's output. Save `Av_in = Av` before PP, dump `Av_in` for
  the M2.5 record, dump post-PP `Av`/`Kv` as NEW records `pp_Av`/`pp_Kv`. (FESOM3 has no conflict: PP writes
  `dyn%work%Av`, distinct from the driver's local M2.5 `Av` — but the DUMPED values match.) Both sides pre-zero
  `Kv`/`Av` before PP (it writes only `nz∈[nzmin+1,nzmax-1]`; surface/bottom/below stay a deterministic 0).

- **`target` on the `dyn` dummy is required** for the `UVnode=>dyn%uvnode` / `Kv=>dyn%work%Kv` pointer aliases
  (ifort `#6796`) — FESOM2 declares `dynamics` `target` too; mirror it. The only compile error; otherwise the
  transcription + Debug `-check all` were clean first try (the L15 `elem2D_nodes(1:3,elem)` slice avoids the
  MAX_NV=4 trap; FESOM2's `elem2D_nodes(:,elem)` is `(3,·)` on its runtime path so its `:` is already 1:3).

## L22 — Convective adjustment (M2.8b): the FIRST gate to deliberately perturb the shared T/S (an unstable band) for non-vacuity — the cascade auto-re-verifies (PASSED first gate run)

M2.8b (`mo_convect`: where N²<0, floor `Kv` nodes / `Av` elements to `instabmix_kv`=0.1 — the convective
overturning proxy, run after PP, FESOM2 `oce_ale.F90:3729`) byte-matched FESOM2 `max|Δ|=0` on the FIRST gate
run (2 new fields `moc_Kv`/`moc_Av` → **48 fields**). Built `src/oce/oce_mo_conv.F90` (`mo_convect`). Specifics:

- **The convective branch is UNTESTABLE on a stable T/S — perturb the shared input to make `bvfreq<0`, and the
  cascade re-verifies for free.** The standing prescribed T/S (T↓, S↑ with depth) is statically STABLE
  everywhere → `bvfreq≥0` → the `if (bvfreq<0) Kv=max(Kv,0.1)` floor NEVER fires → a vacuous gate (passes
  trivially, a transcription bug in the condition/`max` uncaught). Fix: add a LOCALIZED unstable band to T — a
  warm subsurface lens `+8·max(0,cos lat·cos lon)·min(nz,8)/8` whose downward warming over the top 8 levels
  overcomes the −0.20·nz(T)/+0.03·nz(S) stabilisation in the warm hemisphere → `bvfreq<0` on **5165 node-levels
  / 9936 elem-levels**, the floor genuinely fires (`max|ΔKv|`=0.090, `max|ΔAv|`=0.096, i.e. ~0.01→0.1). This is
  the FIRST gate to deliberately change the shared T/S — it CASCADES (T→density→hpressure→pgf→vel_rhs→…→every
  downstream field, 40+ records), but **all M2.1-M2.8 records re-verify `max|Δ|=0` automatically** (both sides
  use the identical T) — the cascade is harmless, only the new `moc_Kv`/`moc_Av` are added. Keep it BOUNDED
  (≤+8 °C → T≤~28 °C) and S>0 so the EOS `sqrt(s)` stays well-posed; verify no NaN. (Side effect, noted not
  re-committed: the unstable region clamps PP's `5·max(N²,0)=0` → the M2.8 Ri factor now reaches exactly **1.0**,
  vs the [5.8e-7, 0.93] of the M2.8 commit's stable T/S — `pp_Kv`/`pp_Av` re-verify regardless.)

- **`mo_convect` overwrites `Kv`/`Av` in place → save the PP output before calling it (the L20/L21 pattern).**
  `pp_Kv`/`pp_Av` (the M2.8 records) must dump the PRE-adjustment PP output, but the dump fires at the end after
  `mo_convect` has floored `Kv`/`Av`. Save `pp_Kv_save`/`pp_Av_save` right after PP, dump those for `pp_Kv`/
  `pp_Av` and the live post-adjustment `Kv`/`Av` as NEW records `moc_Kv`/`moc_Av`. The difference localises the
  gate: `pp_*` isolates PP, `moc_*` isolates the convective floor.

- **Port only the reachable branch; OMIT (don't stub) the forcing-coupled `use_momix` block.** FESOM2
  `mo_convect` has three enhancements: `use_instabmix` (convective, M2-live), `use_windmix` (near-surface, uses
  only params+`Kv`/`Av` → transcribe guarded-off) and `use_momix` (TB04 Monin-Obukhov). The momix block reads
  forcing/ice fields not ported yet (`water_flux`/`heat_flux`/`stress_node_surf`/ice/`mo`/`mixlength`/`mo_length`/
  `pmlktmo`) — they don't EXIST, so it cannot compile → OMIT it entirely (a comment marks it deferred to M2.10),
  not a dead `if`-guard. Force `use_momix=.false.` in the gate. **Oracle-side gotcha:** the FESOM2 `mo_convect`
  signature is `(ice, partit, mesh)` and it pointer-assigns `u_ice=>ice%uice` / `v_ice=>ice%vice` /
  `a_ice=>ice%data(1)%values` UNCONDITIONALLY (before the `use_momix` guard), so the shim's `ice_dummy` must have
  those three allocated (size `nod2D`) even though `use_momix=.false.` never dereferences them — extend the M2.3
  `ice_dummy` (which only had `data(2)`/`data(3)` for `compute_vel_rhs`'s `m_ice`/`m_snow`). FESOM3's ported
  `mo_convect(dyn, mesh)` drops `ice` entirely (no momix → no ice).

## L23 — Tracer-solve assembly (M2.9a): the FIRST `Kv` consumer (sourced LIVE), `NaN·0` poisoning from Redi-off slope reads, module-procedure vs free-subroutine symbol mangling (PASSED first gate run)

M2.9a (the tracer diffusion solve: `diff_tracers_ale` = `diff_part_hor_redi` horizontal diffusion + ALE
reconstruct + `diff_ver_part_impl_ale` implicit vertical-diffusion TDMA + `bc_surface`) byte-matched FESOM2
`max|Δ|=0` on the FIRST gate run (9 new `tsol_*` records → **57 fields**). Built into `src/oce/oce_ale_tracer.F90`
(mirrors FESOM2's file, NOT the plan's `oce_solve_tracers.F90` — the L16 layout precedent). Specifics:

- **The TDMA is the FIRST consumer of a PP output — source `dyn%work%Kv` LIVE (post-`mo_convect`), don't
  prescribe it.** Every prior M2 gate PRESCRIBED its diffusivity/viscosity (M2.5 `Av`, M2.8 `uvnode`→`Kv`); here
  `diff_ver_part_impl_ale` reads the `Kv` that PP+`mo_convect` already wrote into `o_ARRAYS Kv` / FESOM3
  `dyn%work%Kv` earlier in the SAME gate. Append M2.9a after the M2.8b section so `Kv` is live; both sides consume
  `Kv∈[0, 0.1]` (the `mo_convect` floor) — non-vacuous AND it realises the integration the step (M2.9b) needs.
  Gate `max|dT|`=1.26 °C / `max|dS|`=0.11 confirms the solve actually moves T/S.

- **On pi the tracer-solve diffusion reduces to the implicit vertical TDMA ALONE; transcribe THAT path, defer the
  rest.** `K_hor=0` (horizontal diffusion `diff_part_hor_redi` is a no-op in production), `i_vert_diff=.true.`
  (skip the explicit `diff_ver_part_expl_ale`), T/S use `'FCT'` → `do_wimpl=.false.` (the implicit vertical
  ADVECTION terms are off — FCT carries advection explicitly via the full `w`), `Redi=.false.` → `isredi=0` kills
  every isoneutral `slope_tapered`/`Ki` term, `mix_scheme='PP'` → `use_kpp_nonlclflx=.false.` (no nonlocal flux),
  `smooth_bh_tra=.false.` (no biharmonic). So the FESOM3 TDMA = `Kv` vertical diffusion + the `bc_surface` row;
  the do_wimpl advection block is transcribed-but-unexercised, the rest deferred. **linfs makes the ALE
  reconstruct's `del_ttf += T·(hnode-hnode_new)` term vanish** (`hnode_new==hnode`) → `T* = T + del_ttf/hnode`
  (the zstar non-trivial case waits for M2.11). To exercise `diff_part_hor_redi` non-vacuously anyway, PRESCRIBE
  `Ki>0` (the L21 prescribe-the-unsourced-input pattern: `Ki` needs `mesh_resolution`, deferred to M4) and source
  `tr_xy` from the gated M1.1 `tracer_gradient_elements` (re-run per tracer; `do_oce_adv_tra` normally sets it).

- **Redi-off reads `slope_tapered`/`tr_z` BEFORE multiplying by `isredi=0` → `NaN·0=NaN` poisons the oracle.**
  The shipped namelist has `Fer_GM=.true.`/`Redi=.true.` so `slope_tapered`/`tr_z` are ALLOCATED but not yet
  FILLED at end-of-`ocean_setup`; FESOM2 computes `Fx=Kh·(Tx + SxTz·isredi)` where `SxTz=Σ(Tz·slope_tapered)/2` —
  if `slope_tapered`/`tr_z` are uninitialised garbage/NaN, `SxTz·0 = NaN` and `Fx=NaN`, while FESOM3 (which OMITS
  the slope terms for Redi-off) computes a finite `Fx=Kh·Tx` → the gate would MISMATCH (or NaN-propagate). Fix:
  the shim explicitly `slope_tapered=0`/`tr_z=0` (guarded by `allocated`) before the solve, so `SxTz=0` cleanly
  and both sides agree on `Fx=Kh·Tx` (`Tx+0.0·0.0 == Tx` to the bit). General rule: when transcribing a
  feature-OFF branch that DROPS a `×flag` term FESOM2 still evaluates, zero the dropped term's operands in the
  oracle so its `×0` is a clean finite 0, not `NaN·0`.

- **Oracle-side: distinguish module procedures from free subroutines — the symbol mangling differs.** The shim
  drives REAL FESOM2 routines via either a `use` of an auto-generated `*_interface` module (e.g.
  `diff_tracers_ale_interface`) OR an explicit `interface` block (for free subroutines with no `*_interface`
  module, like `impl_vert_visc_ale`/`mo_convect`). `tracer_gradient_elements` is a MODULE PROCEDURE of `o_tracers`
  → its symbol is mangled (`o_tracers_mp_tracer_gradient_elements_`), so an explicit external `interface` block
  emits a call to the UNMANGLED `tracer_gradient_elements_` → **`undefined symbol` at load time** (links fine,
  fails at runtime `symbol lookup error`). Fix: `use o_tracers, only: tracer_gradient_elements` (let the module
  provide the mangled name). Check `module … contains` membership before choosing interface-block vs `use`.

- **Two transcription gotchas the build caught.** (1) Fortran specification-expression order: a dummy whose
  array bound references ANOTHER dummy (`tr_xy(2, mesh%nl-1, …)`) must have that dummy (`mesh`) DECLARED EARLIER
  in the spec part — Intel `#6158`/`#6415`; reorder the type-decls (the arg LIST order is independent and stays
  matching the call). (2) `real_salt_flux` lives in `g_forcing_arrays`, NOT `o_ARRAYS` with its sibling surface
  fluxes `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` — `use` the right module per symbol. Keep FESOM2's
  `zinv1=zinv2` backup bookkeeping in the TDMA verbatim (`1/dz` of the layer above) rather than recomputing
  `1/(Z_n(nz-1)-Z_n(nz))` inline — byte-identical, but faithful transcription removes reviewer doubt (L7).

## L24 — `step_oce` ASSEMBLY (M2.9b): the whole ocean step byte-matches the REAL `oce_timestep_ale` via its own built-in dumps; physical UV to avoid the `check_blowup` stop (PASSED first gate run)

M2.9b assembles every M2.1-M2.9a leaf kernel into one driver (`mod_step_oce::step_oce`, mirroring FESOM2
`oce_timestep_ale` + the `compute_vel_nodes` the main loop runs just before it) with **LIVE data flow** — each
kernel reads the PREVIOUS kernel's output, not a prescribed input — and gates the WHOLE step substep-by-substep.
- **The oracle dumps come for free from the REAL driver.** FESOM2's `oce_timestep_ale` already has built-in
  `dump_shim_record_node` calls at every substep, and FESOM3's `mod_dump` is byte-format-identical to that shim
  (M0.6). So the gate just prescribes the clean state + forcing arrays + reduced-M2 config and calls the REAL
  `compute_vel_nodes` + `oce_timestep_ale` ONCE — no new oracle dump code. FESOM3's `step_oce` emits the same
  records inline (same substep IDs/names), `dump_diff.py` compares (13 NODE fields × 5 probes = 65 records).
  This is the **mod_dump / `dump_diff.py`** gate vehicle (the run_oracle path), DISTINCT from the per-kernel
  FADVHDMP / `pressure_diff.py` gate — both coexist; the step gate is its own `run_step_gate.sh`.
- **Use the pi namelist (KPP/GM/Redi) + FORCE the reduced-M2 dispatch in the shim** rather than authoring a
  reduced-M2 namelist. The shipped pi namelist guarantees ALL the GM/Redi/KPP arrays (`sw_alpha`/`sigma_xy`/
  `neutral_slope`/…) are ALLOCATED at `ocean_setup`; the shim then just redirects the DISPATCH (`mix_scheme_nmb=2`,
  `Fer_GM=.false.`, `Redi=.false.`, `opt_visc=7`, …). The dead-in-M2 producers (`sw_alpha_beta`/`compute_sigma_xy`/
  `compute_neutral_slope`) still RUN inside `oce_timestep_ale` (FESOM3's `step_oce` OMITS them), so they dump the
  SW_AB substep (id=2) that FESOM3 never emits → pass `--ignore-substep=2` to `dump_diff.py` (a clean, documented
  deferral) rather than faking the records. The "reduced-M2 namelist running past forcing" the HANDOFF anticipated
  was NOT needed — the prescribe-and-stop shim drives `oce_timestep_ale` standalone (it reads forcing ARRAYS, not
  files; the 1-rank forcing-file hang L8 is bypassed exactly as for M2.1-M2.9a).
- **`max|Δ|=0` first run by L9 transitivity** — every operand was already byte-pinned per-kernel; M2.9b only adds
  the WIRING (uvnode← `compute_vel_nodes`, Av← PP into `impl_vert_visc_ale`, Kv← PP into the tracer TDMA, the full
  `solve_tracers_ale` loop, the linfs `update_thickness_ale` no-op). The data flow being correct is exactly what a
  whole-step substep match proves. `compute_vel_nodes` is validated transitively: PP's Kv (a dumped substep)
  consumes it, so Kv matching ⇒ uvnode matched.
- **`update_thickness_ale` is a no-op for linfs** (FESOM2 `oce_ale.F90:1226`): neither the zlevel nor the zstar
  thickness-redistribution branch is taken (hnode/helem/zbar_3d_n/Z_3d_n fixed), only `exchange_elem(helem)` —
  a 1-rank no-op. The `hnode` THICKNESS dump equals the init value on both sides.
- **`use_wsplit=.false.` is FORCED (M1.4 precedent).** `do_oce_adv_tra`'s `use_wsplit=.true.` path needs the FCT
  implicit vertical-advection correction `adv_tra_vert_impl` (a distinct kernel, unported until 2026-10-02 — L58 —
  and until then guarded with an explicit `error stop`). With `.false.` the explicit/implicit split is trivial (`w_e=w`, `w_i=0`) AFTER `vert_vel_ale`, but
  `impl_vert_visc_ale` runs BEFORE `compute_Wvel_split` so it still consumes the PRESCRIBED `w_i` ≠ 0 (vertical
  momentum advection IS exercised). The w-split itself is gated at M2.7. Porting `adv_tra_vert_impl` +
  `use_wsplit=.true.` is a scoped follow-up.
- **Physical prescribed UV (0.50/0.40 m/s), NOT the viscosity gate's 2.0/1.5.** The strong-current stress test
  drives the single-step elevation `eta_n ~13 m`, past FESOM2's `±10 m` `check_blowup` guard → `oce_timestep_ale`
  STOPs ("eta_n become NaN or <-10,>10") AFTER all substep dumps fire (the gate still PASSES — dumps complete
  before the stop, verified by `hnode` matching). But an abnormal-termination + "NaN" log makes a permanent gate
  look fragile. Reducing UV to a physical 0.50/0.40 keeps `eta_n ∈ [-4.6, 1.9]` so the REAL step COMPLETES cleanly
  (the shim's own `stop` fires). The gate's correctness rests on the dump COMPARISON (all 13 substeps present +
  matched), not the exit code (`run_stepdump_pi.sh` tolerates the exit but checks dump existence + completeness),
  so a blowup would NOT silently pass — but a clean run is the trustworthy default. Branch coverage for the strong
  flow (the γ0/γ1/γ2 viscosity selection) stays at M2.4, not here.

## L25 — Forcing READ (M2.10a): the FIRST netCDF I/O — the 1-rank async-netcdf infinite-recursion fix, per-node partition-independence, the noleap `julday=365·yyyy` (710820 not 2.4e6) two-stage time-interp cancellation (PASSED first gate run)

The first kernel that reads real files (CORE2 NCAR stubs `test/input/global/{u_10,v_10,q_10,ncar_rad,t_10,ncar_precip}.1948.nc`,
192×94, 5 records, ~360 KB each). All 8 atmospheric fields (`u_wind`/`v_wind`/`Tair`/`shum`/`shortwave`/`longwave`/
`prec_rain`/`prec_snow`) `max|Δ|=0` vs the REAL FESOM2 `sbc_do` on pi 1-rank. Built `src/io/mod_io_netcdf.F90`
(thin `use netcdf` wrapper) + `src/forcing/mod_forcing_read.F90` (the read+interp) + `vector_g2r` into
`mod_mesh_rotate` + driver `src/drivers/fesom_forcingdump.F90`; gate `tools/run_forcing_gate.sh`.

- **The L8 "1-rank forcing hang" is a REAL infinite recursion, now FIXED in the oracle.** Empirically: a 1-rank
  full run hangs in `g_sbf::getcoeffld → io_netcdf_workaround::next_io_rank → mpi_topology`. Root cause:
  `next_io_rank_helper` (its own TODO admits it) infinite-recurses when the only rank IS `SEQUENTIAL_IO_RANK`(0) —
  the async-NetCDF IO-rank selector. Patched `port2/fesom2/src/io_netcdf_workaround_module.F90`: at `partit%npes==1`
  short-circuit to sequential I/O on rank 0 (`async_netcdf_allowed=.false.`) — **VALUE-IDENTICAL** (only changes
  which rank reads+bcasts, not the bytes). After the patch the full model runs PAST forcing and through a clean
  first step (`FDBG step=1 uv=0.13 ssh=0.10 Smax=37.4`, no NaN); it then crashes at OUTPUT I/O
  (`io_meandata::init_nod2d_lists`, a SEPARATE 1-rank gather bug) — irrelevant, the prescribe-and-stop shim stops
  before output. This is an UNCOMMITTED FESOM2 working-tree change (like the dump shims).
- **Forcing is PARTITION-INDEPENDENT per node** (each node's value is a pure function of its own geo-coords + the
  global file data + the model time — NO cross-node accumulation), unlike `area`/FCT scatter. So a multi-rank
  oracle would byte-match too; the np=1 fix simply keeps the 1-rank anchor (consistent with every prior gate).
- **The oracle shim** (`port2/fesom2/src/fesom_forcing_dump.F90`) is wired right AFTER `forcing_setup` in
  `fesom_module.F90:315` (NOT end-of-`ocean_setup` — `sbc_ini` only runs inside `forcing_setup`, which is AFTER
  `ocean_setup` and gated on `use_ice=.true.`, true on pi). It pins the model time, drives the REAL `sbc_do`, maps
  `atmdata→arrays` exactly as `gen_forcing_couple::update_atm_forcing:681-694`, dumps (FADVHDMP), stops.
- **`julday` for `noleap`/`none`/`365_days` = `365·yyyy`** (the else branch, `gen_surface_forcing.F90:1887`), so the
  pi CORE2 time axis is `365·1948=710820`-scale — NOT the ~2.43e6 NR Julian-Day-Number of the JRA gregorian case
  (that's where the HANDOFF's "2.4e6 cancellation" figure comes from). Same cancellation PRINCIPLE, smaller magnitude.
- **The time interp is a TWO-STAGE form that must NOT be collapsed**: `coef_a=(data2-data1)/Δt`;
  `coef_b=data1-coef_a·nc_time(t_indx)`; `atmdata=rdate·coef_a+coef_b`. Both `rdate·coef_a` (~710820·a) and `coef_b`
  (~−710820·a) are large and nearly cancel to the O(1) physical value. Refactoring to the algebraically-equal
  `data1+coef_a·(rdate−nc_time)` rounds the large `coef_b` differently → drift. (FESOM2 getcoeffld:1015-1016 +
  data_timeinterp:1041.)
- **Two rdates**: the coefficients are built at the COLD-START rdate (`nc_sbc_ini:643` — NO half-step, clock-init
  day 1 / sec 0) and the data is evaluated at the per-step rdate (`sbc_do:1527` — WITH the `−dt/86400/2` half-step,
  pinned day 1 / sec 43200 / dt 2400). The gate driver replicates BOTH. (For pi both land in the same window
  → `t_indx=1`, so a single-rdate read would also match, but the two-phase form is faithful and general.)
- **Raw slices stay `real(4)`** (the on-disk type) and promote to WP inside the bilinear weight expression — exactly
  FESOM2 (`real(4) sbcdata · real(WP) wgt`); `real4→real8` is exact so it's byte-identical either way, but keeping
  `real(4)` reads the on-disk bytes verbatim. The periodic-lon halo (`nlon+2`, mirror cols 1↔nlon-1 / nlon↔2, then
  `ic_cyclic` ±360) makes lon monotonic for `binarysearch`; lat-flip only if descending (CORE2 ascends → no flip).
- **netCDF link via `nf-config`** in CMakeLists.txt (the SAME spack `netcdf-fortran-4.5.3` the oracle links → byte
  reads) + **`-Wl,-rpath`** to both netcdf-fortran/-c lib dirs so binaries run without `LD_LIBRARY_PATH` (the spack
  modules set compile paths but not runtime). `use netcdf` (the F90 `.mod`) is compatible across intel 2021.5.0
  (netcdf build) ↔ 2022.0.1 (FESOM3) — no need for the F77 `netcdf.inc`.
- **The dump drivers do NOT call `set_partition`**, so `partit%myDim_nod2D` stays 0 — loop forcing over `mesh%nod2D`
  (== `myDim+eDim` at 1-rank), the same convention every FESOM3 kernel uses (`oce_pressure_bv.F90:72`
  `do node=1,mesh%nod2D`). The driver sets `frc%nnod=mesh%nod2D` and the read routines loop `1,frc%nnod`.
- **SCOPE**: M2.10a is the READ only. `heat_flux`/`water_flux` come from the air-sea **obudget** in
  `ice_thermo_oce.F90` (+ `oce_fluxes`), which is **M3 (ice/thermo)** per the plan ("M3 → thermo → oce_fluxes") —
  even for open water. So M2.10b (bulk transfer coeffs `Cd/Ch/Ce` + wind `stress_surf`) + M2.10c (SW penetration
  `sw_3d`) finish M2.10's producible fields; `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` stay prescribed
  until M3. `gen_bulk_formulae.F90` computes ONLY the transfer coefficients (NOT stress/heat/evap).
- **M2.10b (bulk + wind stress, 17 fields, PASSED first gate run).** `ncar_ocean_fluxes_mode` (the LIVE NCAR routine:
  Large&Yeager 2004 + **Large-2009 drag** [the `u10**6` term + the 33 m/s → `2.34e-3` cap, NOT the commented-out
  L-Y2004 6a form], `n_itts=5`, a 3-height Monin-Obukhov stability iteration) → `Cd`/`Ch`/`Ce` — gated against the REAL
  routine. Then `stress_atmoce=Cd·(ρ_air·|Δu|)·Δu` (`Δu=u_wind−(1−Swind)·u_w`, `Swind=0`) and the node→elem
  `stress_surf(elem)=sum(stress_node_surf(elnodes))/3` (`a_ice=0` ⇒ `stress_node_surf=stress_atmoce`). **Byte traps —
  transcribe the un-suffixed default-real literals VERBATIM** (the L16 family): `inc_ratio=1.0e-4`, `inv_rhoair=1./1.3`,
  `tmelt=273.15`, `rhoair=1.3` (MOD_ICE type-defaults / the gen_bulk local) — a `_WP` suffix would round differently;
  also `(ustar*ustar)` not `ustar**2`, `atan(1.0_WP)` kept a runtime call, `test=abs(cd−cd_prev)/(cd+1.0e-8_WP)`.
  **SST + surface ocean velocity are PRESCRIBED** (a fresh `type(t_ice)` dummy supplies the thermo TYPE-DEFAULTS —
  `ice%thermo` is a non-allocatable component so `inv_rhoair`/`tmelt`/`rhoair` auto-initialize; only `srfoce_temp/u/v`
  need allocation). The wind stress + node→elem are inlined in BOTH sides (trivial formulas, L9-transitive — `Cd` is the
  gated substance). Non-vacuous: `Cd∈[5.5e-5,1.7e-2]`, the wide Tair(−44..31)−SST(−1..19) range fires both
  stability branches. `elem2D_nodes(1:3,elem)` slice (FESOM3 MAX_NV=4; FESOM2's is (3,·)) — the L15 trap. The bulk
  consumes the LIVE M2.10a `u_wind`/`v_wind`/`Tair`/`shum` (a real read→bulk assembly).
- **M2.10c (shortwave penetration, 20 fields total, PASSED first gate run).** `cal_shortwave_rad`
  (`oce_shortwave_pene.F90`, Morel&Antoine 1994 / Sweeney 2005): `swsurf=(1−albw)·shortwave·0.54`; `heat_flux+=swsurf`
  (the visible band is REMOVED from the +upward non-solar `heat_flux` and redeposited as `sw_3d`); chl floor `0.02`;
  the v1/v2/sc1/sc2 polynomials in `c=log10(chl)`; the two-exponential `sw_3d(k)=swsurf·(v1·exp(zbar/sc1)+v2·exp(zbar/
  sc2))` over `zbar_3d_n` (the per-node ALE depth, pressure-gate-proven) with the `aux<1e-5`/`k==nzmax` cutoff;
  `swsurf/=vcpw` (W/m²→K·m/s, `vcpw=4.2e6` — exactly representable so `_WP` is moot). Consumes the LIVE M2.10a
  `shortwave`; `chl`/`heat_flux`/`a_ice=0` prescribed. **`albw=0.066` forced to `0.066_WP` on BOTH sides** (the shim
  sets `ice%thermo%albw=0.066_WP`) so the un-suffixed default-real ambiguity can't bite. The shim's dummy ice needs
  `ice%data(1)%values` (= `a_ice`) allocated (the routine pointer-assigns it); the REAL `cal_shortwave_rad` is driven
  via an explicit interface. Non-vacuous: the chl floor fires on 476 polar nodes; `sw_3d` decays to 0 over 48 levels.
- **Build/CMake recap (M2.10):** netCDF added to FESOM3 via `nf-config` (include+`--flibs`) + `-Wl,-rpath` to the
  netcdf-fortran/-c lib dirs (self-contained binaries; no `LD_LIBRARY_PATH`). `use netcdf` (the F90 `.mod`) works
  across the intel 2021.5.0 (netcdf) ↔ 2022.0.1 (FESOM3) minor-version gap. `src/io/` + `src/forcing/` were already in
  `FESOM3_LIB_DIRS` (auto-globbed once created). One forcing gate (`run_forcing_gate.sh`) now covers all 20 fields
  (read→bulk→SW pene, a real producer chain on the LIVE read); `pressure_diff.py` is the generic comparator (its
  "PRESSURE/EOS…" trailer is cosmetic).

## L26 — CORE2 geometry byte-gate (M2.11a): the L8 CW-swap deferral CLOSED; the proven geom pipeline scales 40× with NO code change (PASSED first run)

M2.11a (the FIRST gate on the CORE2 mesh — nod2D=126858 / elem2D=244659 / edge2D=371644 / edge2D_in=362333 / nl=48,
~40× pi) byte-matched FESOM2 `max|Δ|=0` on the FIRST run, on **all 19 geometry fields** (coord/geo_coord/elem2D_nodes
+ elem_area/elem_cos/metric_factor/gradient_sca + edge_dxdy/edge_cross_dxdy + area/areasvol/(inv) + edges/edge_tri +
the level arrays). Run it: `tools/run_geom_gate_core2.sh`. Reusable specifics:

- **The L8 CW-orientation-swap deferral is now EMPIRICALLY CLOSED.** pi has 0/5839 swaps, so `enforce_cw_orientation`'s
  reorder path (`mod_mesh_read.F90:150-174`: `r=b1·c2−b2·c1`, swap nodes 2,3 when `r>0`) was byte-identical-to-FESOM2-
  test_tri *by construction* but NEVER exercised. CORE2 is stored CCW → FESOM3 reports **CW swaps = 244654/244659**
  (≈100%), and the post-swap `elem2D_nodes` (integer node order) AND every geometry quantity that depends on it
  (centroid `sum/3`, `elem_area`, `gradient_sca`, the per-node `area=Σelem_area/3` accumulation) are `max|Δ|=0` vs the
  oracle's runtime `test_tri`. The reorder is the SAME b/c/trim_cyclic/r>0 logic operating on the SAME byte-proven
  coords, so it matches — but it is now *gated*, not just *argued*. This was the single biggest M2.11 risk; it's gone.
- **The whole geom pipeline scaled 40× with ZERO FESOM3 code change.** `fesom_geomdump` already parameterized the mesh
  dir (`FESOM3_MESH_DIR`); `read_mesh`/`compute_geometry`/`enforce_cw_orientation` are all allocatable-sized and
  mesh-agnostic. So the FESOM3 side is just `FESOM3_MESH_DIR=<core2> mpirun -n1 fesom_geomdump` — no new code, no
  rebuild. The diff was exactly 0 (NOT the L10/L19 uniform-ULP signature), confirming the M2.10 build is sound for
  geometry and an incremental rebuild was unnecessary here.
- **The single-rank `dist_1` is hand-crafted + MESH-INDEPENDENT in its hard part.** FESOM2 always reads a partition;
  METIS can't make dist_1, so `tools/make_dist1.py` mirrors `save_dist_mesh` np=1 (validated by re-parsing the proven pi
  dist_1): `rpart.out`=npes(1)/count(nod2D)/identity-map 1..nod2D; `my_list00000.out`=`0` then per-domain
  `myDim`/`eDim`(0)[/`eXDim`(0)]/identity-list (each SCALAR on its OWN line — list-directed READ starts a new record per
  statement so a scalar read consumes a whole line; LISTS may wrap freely since a list read spans records until full);
  `com_info00000.out`=empty np=1 comm structs (rPEnum=sPEnum=0, blank zero-size arrays, rptr=sptr=1) — this file is
  **byte-identical across meshes** (no nod2D/elem2D content), so COPY pi's verbatim. The oracle read it first try
  ("rpart is read / myLists are read / communication arrays are read"). Lives in `<core2>/dist_1/` (the mesh dir, like
  pi's), NOT in the FESOM3 repo (external artifact; the generator is the tracked deliverable).
- **The geom-dump shim STOPS at mesh_setup (MPI_FINALIZE+stop, `fesom_geom_dump.F90:88`), BEFORE ocean_setup/forcing**,
  so the CORE2 oracle run is ~2 s and reads ONLY the mesh (no IC/forcing — bypasses every downstream 1-rank issue).
  CORE2 `test_tri`=3.5 ms, `load_edges`=0.38 s, 246 MB dump. Fast even at 40× pi.
- **Safe deferrals CONFIRMED on CORE2** (read the mesh, don't guess): min `nlevels=5` for nodes AND elems → NO
  single-layer columns (the L18 `impl_vert_visc_ale` benign-OOB guard is NOT needed, same as pi); cavity + partial-cell
  OFF. So the M2.11b/c kernels inherit pi's safe assumptions. (The deferred cavity/partial-cell + the FCT AUX-scratch
  cavity caveat (L11) still await a cavity mesh — `pi_cavity` exists in tests/data, a later gate.)

## L27 — Initial conditions `do_ic3d` (M2.11b): the 3D-climatology IC byte-matches on CORE2 1-rank; the order-dependent `extrap_nod3D` sweep is deterministic at 1-rank; insitu→potential RK4 (PASSED first gate run)

M2.11b ported the 3D-climatology initial conditions (`oce_initial_state`→`do_ic3d`) and byte-matched FESOM2 `max|Δ|=0`
on the FIRST gate run — all **3 fields** (`Z_3d_n` input, `ic_temp` potential-T, `ic_salt`) on CORE2 1-rank, in ~19 s
total. New: `src/oce/oce_initial_state.F90` (`do_ic3d` + `nc_readGrid` + `nc_ic3d_ini` + `getcoeffld` + `extrap_nod3D`),
`insitu2pot`/`ptheta`/`atg` in `oce_pressure_bv.F90`, `nc_get_var3d_dp` in `mod_io_netcdf.F90`, driver
`src/drivers/fesom_icdump.F90`, oracle shim `port2/fesom2/src/fesom_ic_dump.F90`. Run it: `tools/run_ic_gate_core2.sh`.
Reusable specifics:

- **The IC reads phc3.0_winter.nc (360×180×33), NOT woa18.** Confirmed from `work_core/namelist.tra` `&tracer_init3d`:
  `n_ic3d=2`, `idlist=2,1`, `filelist=2× phc3.0_winter.nc`, `varlist='salt','temp'`, `t_insitu=.true.`. So **salt (ID 2)
  is read FIRST into `tracers%data(2)%values`, temp (ID 1) SECOND into `data(1)%values`** (the `idlist` ordering, not the
  data-slot ordering), then `data(1)` is converted in-situ→potential. `ClimateDataPath=/pool/.../INITIAL/phc3.0/` — FESOM3
  reads the SAME pool file the oracle does (`FESOM3_IC_FILE`).
- **phc3.0 land = NaN (no `_FillValue`), so the missing-value mask reduces to `ieee_is_nan`.** The file has no
  `_FillValue`/`missing_value`; valid data is in [-2.1, 41.4]. FESOM2's `nf_inq_var_fill` returns the *default*
  `NF_FILL_DOUBLE`≈9.97e36, but no point equals it, and the `v>dummy`(=1e10) range check subsumes that branch anyway. So
  FESOM3 detects missing with `ieee_is_nan(v) .or. v<-0.99·dummy .or. v>dummy` (drop the `==FILL_VALUE` term — provably
  vacuous here). `nf90_get_var` (F90) does **no auto fill/scale masking** (base-library, like F77 `nf_get_vara_double`),
  so NaN land bits pass through verbatim for the `ieee_is_nan` test — bytes identical to the oracle.
- **The IC data is read as real(8) (double), NOT real(4) like the M2.10 forcing.** FESOM2 `getcoeffld` uses
  `nf_get_vara_double` into a `real(WP)` cube; FESOM3 `nc_get_var3d_dp` reads into `real64` then promotes to WP (==real64
  at the anchor). The bilinear weights `(x2-x)·(y2-y)/denom` + the vertical-interp `cf_a·Z+cf_b` are per-node-independent
  WP arithmetic on byte-identical operands (`geo_coord_nod2D/rad` geometry-proven by M2.11a, `Z_3d_n` linfs full-cell =
  reference mid-depth at init, `nc_lon/lat/depth` read identically) → `max|Δ|=0` by L9 transitivity. `forcing_binarysearch`
  (mod_forcing_read) is the SAME `d=1e-9` bisection the IC needs — reused verbatim (no second copy).
- **`extrap_nod3D` is THE order-dependent step, and it's deterministic at 1-rank.** Its Gauss-Seidel surface sweep
  (`do while(success)`) reads/writes a `work_array` in NODE order, accumulating valid neighbours through
  `nod_in_elem2D`/`elem2D_nodes` in their stored order; the outer `do while(glob_max>0.99·dummy)` repeats until no surface
  dummy remains; then a downward vertical fill. At 1-rank `exchange_nod` is a no-op and `eDim=0`, so the node order +
  `nod_in_elem2D` order match FESOM2's global order (geometry-gate proven, L9) → byte-identical. This is the one step that
  would diverge at multi-rank (deferred to M2.12). `elem2D_nodes(1:3,el)` slice (NOT `(:,el)`) avoids the L15 MAX_NV trap.
- **`insitu2pot` uses the 1-D `Z(nz)` (not `Z_3d_n`) for the pressure proxy** `pp=abs(Z(nz))` (FESOM2 `:3113`, a
  partial-cell stability choice at init). `ptheta`/`atg` are verbatim Bryden-1973 RK4 with `_WP`-suffixed literals; mutating
  the by-reference dummies `t,p` is harmless (the caller passes scalar temporaries re-read each iteration).
- **The oracle shim is trivial: just dump the LIVE tracers.** `oce_initial_state` runs at `oce_setup_step.F90:253` (early in
  ocean_setup); by the end-of-ocean_setup dump point `data(1)/(2)%values` already hold the final potential-T/S (lines 256-259
  only copy values→valuesold). So `ic_dump_write` PRESCRIBES NOTHING — it dumps `data(1)/(2)%values` + `mesh%Z_3d_n` and
  stops (no kernel call, unlike the pressure shim). Wired FIRST (before `advhor_dump_write`) so the IC is pristine.
- **Oracle namelist override: `which_ALE 'zlevel'→'linfs'`** (+ `mix_scheme→PP`, `Fer_GM/Redi→.false.`, the reduced-M2 set).
  At init (eta=0, full cells) linfs and zlevel give byte-identical `Z_3d_n`, but forcing linfs matches FESOM3's linfs build
  exactly AND lightens ocean_setup. **CORE2 ocean_setup completes at 1-rank in ~19 s** (init_stiff on 126858 nodes + the
  extrap sweep) — no hang; the dump stops before forcing_setup (the L8 hang is downstream). `Z_3d_n` is gated as an input,
  so the linfs-build assumption is verified, not assumed (`max|Δ|=0`).
- **Range-print red herring:** the oracle's `do_ic3d` "global min salt" (5.628) loops only WET levels
  (`ulevels_nod2D:nlevels_nod2D-1`), while the FESOM3 driver's `minval` spans the full array incl. the bottom-zeroed 0.0 →
  the two *printed* S-mins differ (5.628 vs 0.0). This is NOT a data mismatch — the byte-gate compares the full arrays
  (incl. the bottom zeros) and is `max|Δ|=0`. Don't chase printed summary stats; trust the field-by-field diff.
- **Declaration-order trap (Intel #6415):** a dummy array whose bounds reference another dummy
  (`arr(mesh%nl-1, mesh%nod2D)`) must be declared AFTER that dummy (`type(t_mesh) :: mesh`) — Intel flags "name conflicts
  with prior uses" otherwise. `pressure_bv`/`insitu2pot` already declare `mesh` first; mirror that order.

## L28 — Full lifecycle (M2.11c): the FIRST real multi-step run; the whole dynamical core byte-matches on CORE2; the free-surface CG hits the iterative-solver reproducibility floor

> ⚠️ **SUPERSEDED by L29 (2026-06-21).** The "CG reproducibility floor" diagnosed below was NOT a floor — it was an
> auto-vectorised preconditioner divide (`pr_values`, an un-gated CG operand). Fixed with one `!DIR$ NOVECTOR`;
> `d_eta` and the full multi-step CORE2 lifecycle are now `max|Δ|=0`. The localization facts below are accurate;
> only the "irreducible floor" CONCLUSION (and the "OPEN DECISION" it framed) are wrong. Read L29.

M2.11c built the FIRST real time-stepping run (not a prescribe-and-stop shim): FESOM3's
`src/drivers/fesom_lifecycle.F90` (cold-start CORE2 mesh + `do_ic3d` phc3.0 IC + N-step
runloop calling `mod_step_oce::step_oce`) vs the REAL FESOM2 multi-step lifecycle
(`tools/run_lifecycle_core2.sh`, the built-in per-substep `dump_shim` over N steps). Result:
the **entire ported dynamical core is byte-identical (`max|Δ|=0`) on the 40× CORE2 mesh**, with
ONE exception — the free-surface CG solve `d_eta`, which sits at the iterative-solver
reproducibility floor. Reusable specifics:

- **Multi-step AB2 needs NO extra bookkeeping — persist the state and pass `lfirst=(n==1)`.** The
  velocity Adams-Bashforth array `uv_rhsAB` rotates IN PLACE inside `compute_vel_rhs` (reads the
  previous step's slot, overwrites with this step's Coriolis+momadv); the tracer AB history
  `valuesold` rotates IN PLACE inside `init_tracers_AB` (`valuesold(1)=values` each step, AFTER
  computing `valuesAB`). So the driver just loops `step_oce(n, dt, n==1, …)` on the SAME persistent
  `dyn`/`tracers`/`mesh` — do NOT reset `uv_rhsAB`/`valuesold` between steps. Cold start: `valuesold=values`
  (=IC), `uv_rhsAB=0`, `eta_n/UV/w=0` — both sides identical, so step 1 is Euler (`ff=1.0`), steps 2+ AB2.
- **Unforced first gate (`use_ice=.false.`) cleanly isolates the lifecycle from forcing/ice.** In FESOM2
  `forcing_setup` is internally guarded by `if (use_ice)` (`gen_forcing_init.F90:43`) and the whole
  per-step forcing/flux block (`fesom_module.F90:673-715`) is `if(use_ice)` — so `use_ice=.false.` skips
  ALL forcing reads + the air-sea budget; the surface flux arrays (`heat_flux`/`water_flux`/`virtual_salt`/
  `relax_salt`/`stress_surf`) stay at their `arrays_init` zeros (`oce_setup_step.F90:980-997`). FESOM3
  prescribes them all = 0. This gates the lifecycle + multi-step evolution + `step_oce` on CORE2 with the
  real `do_ic3d` IC, with zero forcing complexity. (M2.11c-2 forced: prescribe the oracle's per-step
  fluxes — but it inherits the same CG floor below.)
- **Two CORE2 1-rank crashes the unforced run hits — both the SAME `io_gather init_nod2D_lists` bug.**
  (1) `output()` (`io_meandata.F90`) and (2) the restart write `write_initial_conditions` → `ini_ocean_io`
  (`io_restart.F90`, called UNCONDITIONALLY on the first step's init, regardless of `restart_length_unit='off'`)
  both reach `io_gather::init_nod2D_lists`, which derefs an unallocated `remPtr_nod2D` / a size-0
  `rank0List_nod2D` on the sole rank (the gather machinery is multi-rank-only). Fix = a 1-line
  `if (partit%npes==1) return` at the top of EACH (uncommitted oracle instrumentation, like the dump shims).
  `finalize_output` is then a safe no-op. The gate needs neither output nor restart.
- **`use_sw_pene` MUST be OFF when `use_ice=.false.` — else a segfault, not a numeric mismatch.** With
  `use_ice=.false.`, `cal_shortwave_rad` (called only inside the skipped `oce_fluxes`) never runs, so `sw_3d`
  (allocated only in `gen_forcing_init`) stays UNallocated; but the tracer TDMA `diff_ver_part_impl_ale`
  derefs `sw_3d` under `if (use_sw_pene .and. ID==1)` → null-pointer segfault at step 1. `use_sw_pene=.false.`
  also matches the ported `step_oce` (M2.9 never consumes `sw_3d`). The CORE2 production namelist has
  `use_sw_pene=.true.` → override it in the reduced-M2 run dir.
- **THE CG REPRODUCIBILITY FLOOR (the one non-`max|Δ|=0` field).** `tools/run_pressure_gate_core2.sh` (the
  M2.1-M2.9 per-kernel pressure gate, run on CORE2) DEFINITIVELY localizes it: **all 27 dynamical-core fields
  are `max|Δ|=0` on CORE2** — `density`/`hpressure`/`bvfreq`/`pgf`/`coriolis`/`uv_rhs*`/`visc`/`uv_rhs_ivv` +
  the FULL SSH stiffness matrix `ssh_stiff_diag` (diagonal) AND `ssh_Aeta` (the matvec `A·eta_n`, which
  exercises every CSR nonzero) AND `ssh_rhs` + the M2.8 `pp_Kv`/`pp_Av` + M2.9 tracer-solve `tsol_*`. ONLY
  `d_eta` (the preconditioned-CG solution) diverges, at `~4e-14` (the strong-UV analytic gate; `~1.4e-16`
  ≈ 1 ULP at the physical cold-start lifecycle step 1), and everything DOWNSTREAM of it inherits the seed
  (`uv_upd`/`hbar`/`eta_n`/`w`). So with **byte-identical A (matvec+diag), b (`ssh_rhs`), and x0=0, the
  136-iteration CG nevertheless accumulates a sub-ULP/iteration rounding difference** between the two
  separately-linked binaries (oracle `libfesom.so` whole-model vs FESOM3 `libfesom3.a`), both Intel 2021.5.0,
  same anchor flags (`-fp-model precise -no-prec-div -ip`), `__openmp_reproducible` NOT defined (both use the
  serial DO-loop dot-products), `ENABLE_OPENMP=OFF`. Ruled out: the matvec (replacing FESOM3's `sum()` CSR
  matvec with an explicit sequential DO-loop changed the result by **0** — FESOM3's `sum()` is already
  sequential and equals the oracle's, consistent with `ssh_Aeta` passing); the dot-products (serial DO, same
  range `myDim==nod2D` at 1-rank); the preconditioner (built from byte-identical `values`+`diag_values`);
  iteration count (`d_eta` Δ `4e-14` ≪ the `soltol=1e-5` tolerance → same count). The residual ULP almost
  certainly comes from `-ip` instruction-scheduling of the long CG recurrence differing across the two link
  units. On **pi** (3140 nodes, 37 CG iters; M2.6/M2.9b) this floor was below the last bit → `max|Δ|=0`; on
  **CORE2** (126858 nodes, 136 iters) it surfaces. This is the documented [[project-bit-identity-reality]]
  floor for a global iterative solve. **Consequence:** a multi-step CORE2 `max|Δ|=0` gate is BLOCKED by this
  seed (it amplifies chaotically: `~1e-16` → `~3e-6` over 3 steps). The achievable+achieved CORE2 gate is the
  per-kernel `run_pressure_gate_core2.sh` (`max|Δ|=0` on the whole dynamical core; `d_eta` at the floor).
  **OPEN DECISION for the user:** accept the CG floor (declare M2.11c closed "to the iterative-solver floor")
  vs a deeper fix (e.g., a deterministic-reduction CG, or matching the oracle's exact `-ip` codegen). See
  HANDOFF "Next task".
- **M2.11c-2 (FORCED lifecycle): the ice-at-1-rank unknown is DE-RISKED + the M3-gap flux prescription is byte-exact.**
  The REAL forced FESOM2 lifecycle (`use_ice=.true.`, real CORE2 NCAR forcing at the 1948 stubs + pool runoff/SSS, the ice
  EVP + `oce_fluxes` air-sea budget) runs CLEANLY at 1-rank on CORE2 (`tools/run_lifecycle_forced_core2.sh`). Two setup
  facts: (1) **CORE (noleap) forcing REQUIRES `include_fleapyear=.false.`** — FESOM2 stops with a calendar-consistency error
  otherwise (the noleap `julday=365·yyyy`, L25); the work_core default is `.true.` (for the JRA gregorian production forcing).
  (2) **`use_sw_pene=.false.`** for the gate so `cal_shortwave_rad` is skipped → `heat_flux` is the raw obudget value and the
  tracer TDMA has no `sw_3d` term — matching the ported `step_oce` (M2.9). (Forcing paths: `make_full_path` prepends
  ClimateDataPath only for RELATIVE paths, so absolute stub/pool paths are used verbatim.) The forced gate
  (`tools/run_lifecycle_forced_gate_core2.sh`): oracle dumps the per-step fluxes via `fesom_flux_dump.F90` (full-field
  `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt` + `stress_surf`, BEFORE `oce_timestep_ale`); FESOM3 `fesom_lifecycle`
  reads them (`FESOM3_FLUX_FILE`) and prescribes them into `step_oce` (the M2.5/M2.8 prescribe-the-unsourced-input pattern).
  **Result identical to unforced:** every pre-CG substep `max|Δ|=0` — INCLUDING `ssh_rhs`, which now consumes the prescribed
  wind stress (so the stress prescription is byte-exact). *(Pre-L29 only `d_eta`+downstream diverged at the CG floor, first
  divergence `5.5e-17`; **POST-L29 the forced gate is `max|Δ|=0` on ALL substeps** — 195 records, re-verified 2026-06-21. See L29.)*
  So the forced dynamical core is byte-exact on CORE2 with REAL air-sea fluxes; the M3 ice/thermo budget is the only piece
  prescribed (it's genuinely M3). Debug `-check all` clean (incl. the flux-read). This is the M2-MVP capstone with real forcing.

## L29 — The CORE2 "CG reproducibility floor" (L28) was NOT a floor: an auto-vectorised preconditioner divide. RESOLVED, `max|Δ|=0`

**L28 is SUPERSEDED.** L28 declared the CORE2 free-surface CG `d_eta` divergence (`~4e-14`) an irreducible
iterative-solver reproducibility floor ("all CG inputs byte-identical; the 136-iter recurrence accumulates
sub-ULP/iter across two link units"). **That conclusion was wrong.** A CG with byte-identical `A`, `b`, `x0` and
byte-identical per-iteration kernels is a deterministic recurrence — iteration count does not manufacture
divergence. The divergence meant ONE per-iteration operation differed, and it did: **the preconditioner array
`pr_values` was NOT byte-identical** (it is built inside `ssh_solve_preconditioner` and was never gated).

**Root cause.** FESOM3's precond off-diagonal `K_ri = -0.5*(a_ri/a_rr)/(a_rr+a_ii)` was **auto-vectorised**
(packed `divpd`/`mulpd`) by the compiler; the FESOM2 oracle compiles the SAME source SCALAR (`divsd`). Packed
and scalar division differ by **~1 ULP** under the anchor flags (`-no-prec-div -fimf-use-svml`, SSE2/no-FMA on
Levante) — they are NOT the same micro-op. So `1299 / 870146` `pr_values` entries drifted by `~6.6e-24`. That
seed enters `z = M⁻¹ r` every iteration; it is sub-ULP in the early CG dot-products (so `sum(rhs²)`,
`sum(r0·z0)`, and iters 1–5 all byte-matched), then **surfaces in the residual `sum(r·z)` at ~iter 6** and
cascades through `β → p → x` to `~4e-14` by iter 136. On **pi** (3140 nodes, 37 iters) the seed stayed below the
last bit → `max|Δ|=0` (why pi never saw it); **CORE2** (126858 nodes, 136 iters) surfaced it.

**Why the asymmetry (identical source + flags, different codegen).** The oracle writes the result via a LOCAL
POINTER `pr_values(...)` (`solver.F90:81`); the compiler can't prove it doesn't alias the `ssh_stiff%values`
reads → it keeps the loop SCALAR. FESOM3 wrote the DERIVED-TYPE COMPONENT `ssh_stiff%pr_values(...)`, provably
distinct from the `%values` component → the compiler auto-VECTORISES the divide. Same arithmetic, different SIMD
width, 1-ULP-different result.

**Fix (one line):** `!DIR$ NOVECTOR` on the precond off-diagonal loop in `oce_ssh_solve.F90::ssh_solve_preconditioner`,
forcing the oracle's scalar `divsd`. After it, the F3 precond disassembles to `divsd:3 divpd:0` — identical to the
oracle — and `pr_values` is `max|Δ|=0`. **Verified:** `tools/run_pressure_gate_core2.sh` PASS `max|Δ|=0` on all 28
fields INCLUDING `d_eta`/`eta_n`/`uv_upd`/`hbar`; and the previously-blocked **multi-step gate**
`tools/run_lifecycle_gate_core2.sh` now MATCHes — **195 records (13 substeps × 5 probes × 3 steps), worst |Δ| = 0**.
The whole dynamical core is byte-exact across multiple steps on CORE2. **The FORCED gate too (re-verified 2026-06-21):**
`tools/run_lifecycle_forced_gate_core2.sh` (real CORE2 forcing + ice + prescribed air-sea fluxes) now MATCHes —
**195 records, worst |Δ| = 0** — where pre-L29 its `d_eta` drifted `5.5e-17 → 3.5e-6` over 3 steps. So the NOVECTOR
fix makes BOTH the unforced and forced free-surface CG byte-exact; M2 is byte-identical end-to-end on CORE2.

**How it was localised (reusable method).** (1) A per-iteration scalar trace of `s_old/s_aux/al/sprod(1)/sprod(2)`
written to a file from BOTH binaries (env-guarded, dumped AFTER the loop so the loop codegen is unperturbed). The
diff showed iters 1–5 byte-identical, first harmful divergence `sprod(1)=sum(r·z)` at iter 6 (sprod(2)=sum(r·r)
diverged harmlessly at iter 3 — it only feeds the convergence test, not the recurrence). (2) `sprod(1)` differing
with `sprod(2)` matching ⇒ `rr` identical, `zz` differs ⇒ the preconditioner. (3) An array dump of `pr_values`
(the un-gated build output) confirmed `max|Δ|=6.6e-24` while gated `values` was `0`. (4) `objdump` of both
`ssh_solve_preconditioner` showed `divpd:9` (F3) vs `divpd:0` (oracle). Disassembly is ground truth — the scalar
trace localised the iteration, the array dump localised the array, the disasm localised the instruction.

**Generalisable lessons.**
- **Gate every intermediate that a downstream kernel consumes, not just the named outputs.** `pr_values` was the
  one CG operand never dumped; it hid the bug for an entire milestone. The companion gate dumped `ssh_stiff_diag`,
  `ssh_Aeta`, `ssh_rhs` — but not the preconditioner it builds. (cf. [[feedback-tick-plan-checkboxes]] discipline.)
- **A floating-point divide (or any op) byte-matches only when the OPERANDS *and the SIMD width* match.** Packed
  vs scalar `divpd`/`divsd` differ ~1 ULP under `-no-prec-div`. `-fp-model precise` stops reduction REASSOCIATION
  but NOT auto-vectorisation of element-wise divides. Watch for auto-vectorised divides/reciprocals/sqrt in any
  ported kernel; force the oracle's width with `!DIR$ NOVECTOR` (or match its pointer-vs-component access form).
- **An iterative solver is NOT a reproducibility excuse.** "It's a long recurrence across two binaries" was a
  plausible-but-false story (L28). Same inputs + same kernels ⇒ same output, period. When a recurrence diverges,
  bisect it to the first differing scalar, then to the array, then to the instruction — don't accept a "floor".
- **A sub-ULP seed can hide for many iterations.** `pr_values` drift was invisible in `sum(rhs²)`, `sum(r0·z0)`,
  and iters 1–5; value-dependent rounding only surfaced it at iter 6. Matching a few early probes ≠ byte-identity.

## L30 — Multi-rank local-mesh remap (M2.12a): the geometry pipeline is partition-agnostic; gate per-rank vs same-partition

The first multi-rank byte-gate. FESOM3's local mesh at npes>1 byte-matches FESOM2 per-rank on all 19 geometry
fields (pi dist_2 + dist_8), via `tools/run_geom_gate_multirank.sh`. Reusable lessons:

- **The geometry math is partition-agnostic — only the BOUNDS are local.** `compute_geometry` reads only local
  arrays (`coord_nod2D`, `elem2D_nodes`, `edges`, `edge_tri`) and never global ids. So multi-rank = build the
  LOCAL mesh arrays (scatter the global files through the `myList_*` inverse map) + run the SAME geometry with
  local loop bounds. A `local_bounds(mesh, partit, ...)` helper returns the mesh global counts at npes==1 (the
  proven 1-rank path is byte-for-byte unchanged) and `partit%myDim/eDim[/eXDim]` at npes>1.
- **Gate PER-RANK vs SAME-PARTITION FESOM2, NOT vs 1-rank global (L8).** Both `dist_N` runs read the same
  `my_list<rank>.out`, so local index i ↔ the same global id AND the per-node area accumulation order is the same
  → byte-identical. The 1-rank global uses a different element permutation; non-associative FP would diverge ~ULP.
  Dump per-rank OWNED slices (`1..myDim_*`); at npes==1 owned==global so the existing 1-rank gate is untouched.
- **Owned-entry geometry needs ONLY the element-CENTER halo exchange — not the area or nod_in_elem2D machinery.**
  The one place an OWNED quantity reads a HALO value is `edge_cross_dxdy`: an owned edge can border a halo (eDim)
  element whose center it needs, and `elem2D_nodes` is stored OWNED-only (so `elem_center` can't run on a halo
  element). FESOM2 solves this by precomputing owned element centers and `exchange_elem`-ing them
  (oce_mesh.F90:2528-2529); `edge_cross_dxdy` then reads the exchanged center array, never `elem_center` on a
  halo. Everything else owned is local: owned-node `area` sums over owned adjacent elements (the partition
  guarantees an owned node's full element-neighbourhood is owned), so FESOM2's `exchange_nod(area)` and the
  `find_neighbors` nod_in_elem2D halo dance only fix HALO entries (which the owned gate doesn't compare). Defer
  them until the dynamics actually consume halos (M2.12b) — don't build machinery the current gate can't exercise.
- **File-read arrays fill owned+halo for free; computed arrays need the exchange.** Because the remap reads the
  GLOBAL files, `coord_nod2D`/`nlevels`/`nlevels_nod2D`/`depth` can be scattered to BOTH owned and halo local
  slots directly (no exchange). Only COMPUTED geometry (centers/elem_cos) needs a halo exchange. `elem2D_nodes`
  is the exception — stored owned-only, since a halo element's nodes may not all be local.
- **MP vs WP in the exchange.** Mesh arrays are `MP=max(WP,4)`; `exchange_elem` takes `WP`. At dp/sp `MP==WP` so
  it compiles, but route `elem_cos` (MP) through a `WP` scratch to stay correct if a future build has `MP/=WP`
  (NVHPC half). Element centers are kept `WP` and exchanged directly.
- **`enforce_cw_orientation` is per-element deterministic → consistent across ranks** (a function of the element's
  node coords only), so the swap decision is identical on every rank with no communication. Closes the deferred
  multi-rank CW-swap caveat (pi had 0 swaps at 1-rank; the swap path was first exercised on CORE2 at M2.11a).

## L31 — Multi-rank tracer advection (M2.12b): lift bounds + exchanges, NOT the arithmetic; gate every consumed intermediate

The whole tracer-advection subtree byte-matches FESOM2 per-rank on pi dist_2 + dist_8 (13 fields,
`tools/run_advhor_gate_multirank.sh`). The reusable lessons:

- **Lift a byte-proven kernel to multi-rank by changing ONLY the loop bounds + adding the halo exchanges — never
  the array dummies or the per-element arithmetic.** Keeping the explicit-shape dummies (`ttf(nl-1, mesh%nod2D)`
  etc.) and the exact statements means the COMPILER EMITS THE SAME CODE, so the byte-match vs FESOM2 — including
  any auto-vectorised divide (the L29 `divpd`/`divsd` trap) — is preserved by construction. The actual arrays are
  allocated to LOCAL sizes; the dummy's "global" trailing extent is just metadata never exceeded (loops stay
  `≤ owned`, the leading dims that drive the stride are identical), so it's runtime-safe in Release.
- **OPTIONAL `partit` contains the blast radius.** Make `partit` an optional arg on every advection routine:
  absent ⇒ the 1-rank path runs VERBATIM (global bounds, no exchanges) so the proven 1-rank gates can't regress
  and the existing 1-rank callers (`solve_tracers_ale` etc.) need NO change; present+npes>1 ⇒ owned/halo bounds +
  exchanges. A tiny `mod_part_bounds.owned_bounds(mesh,…,partit)` / `is_multirank(partit)` (npes==1 ⇒ global)
  centralises it. Beats threading mandatory `partit` through the whole step (which ripples into the dynamics).
  NB: Fortran `.and.` does NOT short-circuit — guard with `is_multirank()` (a `present()` check inside), never
  `present(partit) .and. partit%npes>1`.
- **The `find_neighbors` halo dance is the prerequisite the geometry gate couldn't exercise.** MUSCL
  `fill_up_dn_grad` / `find_up_downwind_triangles` loop OWNED edges but read a halo node's FULL element list
  (`nod_in_elem2D`), which reaches the **eXDim** second element-halo layer. So after building owned `nod_in_elem2D`
  you must: `exchange_nod(num)` → per slot pack local→GLOBAL ids, `exchange_nod`, store back → re-localise every
  entry GLOBAL→local through the full-halo inverse map (`imap_elem`, 1..myDim+eDim+eXDim). The eXDim halo
  GUARANTEES a halo node's element list is fully local, so the re-localise never hits a 0.
- **Halo coverage is chosen by ARRAY SIZE in FESOM2; replicate with an explicit `exchange_elem_full`.** FESOM2's
  `exchange_elem` picks `com_elem2D` (eDim) vs `com_elem2D_full` (eDim+eXDim) from `ubound(arr)≤myDim+eDim`.
  Things read at eXDim — `elem_area`, `tr_xy`, and the `coord_elem`/`e_nodes` of the upwind/downwind search — need
  the FULL halo (`com_elem2D_full`). Exchanging the SCALED `elem_area` post-accumulation == FESOM2's
  scale-then-broadcast (a per-element multiply commutes with the owner→halo copy).
- **A reused scratch array is NOT a gateable field at the wrong time.** FESOM2's `oce_tra_adv_fct` reuses
  `edge_up_dn_grad` as its `AUX` scratch (filled with `bignumber=1e3` below the bottom); FESOM3 uses a separate
  allocatable. So `edge_up_dn_grad` only matches BEFORE `do_oce_adv_tra` — capture it right after
  `init_tracers_AB`/`fill_up_dn_grad`, not after the FCT step. (The mismatch is a dump-timing artifact, not a bug:
  every real field — `del_ttf`, `fct_LO`, `fct_plus/minus` — was already `max|Δ|=0`.) General rule: when gating an
  intermediate, dump it at the point its value is live, before any in-place scratch reuse downstream.
- **`nboundary_lay` is computed partition-locally with NO exchange — and that's correct for the gate.** FESOM2's
  `muscl_adv_init` builds it from OWNED edges only (a `min` over each node's owned-edge set) and never exchanges.
  A halo node's value is therefore the partition-local min — possibly not the global min — but FESOM3 computes it
  the SAME way on the SAME partition, so they match byte-for-byte. Don't "fix" it with an exchange; match FESOM2.
- **Prescribe gate inputs halo-consistently.** Node fields are set at owned+halo from the (geometry-gated) rotated
  coords directly; the element velocity has no local `elem2D_nodes` at halo elements, so prescribe it at OWNED
  elements then `exchange_elem_full` — identical on both codes (same coords, same partition) without needing a
  halo exchange of the inputs themselves.

## L32 — Multi-rank pre-SSH dynamics (M2.12c-1): the whole dynamics RHS lifts mechanically; the byte-match is bounds + exchanges, the work is knowing which halos a kernel actually reads

The whole pre-SSH dynamics chain (`compute_vel_nodes → pressure_bv(+smooth_nod) → pressure_force_4_linfs →
oce_mixing_pp → mo_convect → compute_vel_rhs(+momentum_adv_scalar) → viscosity_filter(visc_filt_bidiff) →
impl_vert_visc_ale → compute_ssh_rhs_ale`) byte-matches FESOM2 per-rank on pi dist_2 + dist_8 (25 records
density/pressure/bvfreq/Kv/ssh_rhs, `tools/run_stepdyn_gate_multirank.sh`). The reusable lessons:

- **The M2.12b recipe scales to the whole dynamics with ZERO new ideas:** optional `partit` (absent ⇒ 1-rank path
  VERBATIM, so the proven gates can't regress and `step_oce`'s existing 1-rank callers need no change) + a bounds
  helper (`mod_part_bounds.local_dims` returns `nNodO/nNodL/nEdgeO/nEdgeL/nElemO/nElemL/nElemF`; global counts when
  partit is absent) + the FESOM2 exchanges at the FESOM2 sites, arithmetic/dummies UNCHANGED. Nine kernels, ~1 hour,
  byte-exact first run on dist_2 AND dist_8. Don't reason about physical "completeness" of halo values — replicate
  FESOM2's exact loop bounds + exchanges and the bits follow (the L31 principle, confirmed at scale).
- **The load-bearing skill is knowing which halo each array read needs — read FESOM2's bound, don't guess.** Node
  kernels split: the EOS/mixing/convection node loops run owned+halo (`myDim+eDim`) because a later same-kernel
  ELEMENT loop reads the just-computed node field at the element's 3 corners (halo), and computing the halo in-place
  beats an exchange (uvnode/bvfreq are already halo-valid from upstream). The accumulate-into-node kernels
  (momentum_adv_scalar, compute_ssh_rhs_ale, smooth_nod) loop OWNED edges/nodes, scatter into owned+halo, then
  `exchange_nod` the result (every edge incident to an owned node is owned ⇒ the owned-node value is complete; the
  exchange only fixes the halo). visc_filt_bidiff is the one owned+HALO-edge kernel (`nEdgeL`), with `exchange_elem`
  between its two Laplacian sweeps and the interior-edge test on the GLOBAL id `myList_edge2D(ed)>edge2D_in` (the
  1-rank `ed>edge2D_in` is WRONG at MR — local index vs global threshold).
- **c-1 needs NO mesh-infra change because no pre-SSH kernel reads `elem2D_nodes`/`gradient_sca` at a halo element.**
  They read those only at OWNED elements (an owned element's nodes are local; `compute_vel_rhs`/`pgf`/mixing/convection
  all loop owned elements). The halo reads are UV/UV_rhs/helem/elem_area/edge_cross_dxdy — UV is `exchange_elem_full`'d,
  helem/elem_area are full-halo (M2.12b), edge_cross_dxdy is owned-edge-only. So the local mesh's owned-only
  `elem2D_nodes` (M2.12a) suffices. The SSH STIFFNESS (c-2) is the first kernel that DOES read `elem2D_nodes`/
  `gradient_sca` at eDim-halo triangles of owned edges → it needs the mesh-infra extension. Identify this boundary
  before lifting: it tells you exactly when the cheap bounds-only lift ends and the infra work begins.
- **`exchange_nod` needed a rank-3 node-block variant.** `UVnode`/`UVnode_rhs` are `(2,nl-1,nod)` — the existing
  `exchange_nod` handled only rank-1/rank-2 node fields. Added `exchange_nod_blk_r` (the existing `core_blk_r` on
  `com_nod2D`): the leading two dims are a contiguous per-node block, exactly like `exchange_elem_full`'s 3D variant.
- **Reuse the gid-keyed dump for the MR gate — it is already per-rank.** `mod_dump` writes `<prefix>.<mype5>` and the
  rank owning a probe gid writes it (`resolve_probes` over `myDim`). So a per-rank gid-keyed gate needs NO new dump
  code: each global probe is owned by the same rank on both codes (same `dist_N`/myList), and `dump_diff.py` matches
  by gid. The oracle ran the REAL `oce_timestep_ale` at npes>1 (extend the `fesom_step_dump` shim: drop the `npes/=1`
  return; prescribe nodes owned+halo, element velocity OWNED + `exchange_elem` — prescribing UV at an eXDim element
  SEGFAULTS, its `elem2D_nodes(1)` can be a node beyond the eDim halo so `coord_nod2D` OOB); FESOM3 dumps only the
  pre-ssh_rhs substeps, so `--ignore-substep` the oracle's post-ssh substeps (9/11/12/13/15/16 + SW_AB 2).

## L33 — Multi-rank SSH stiffness + free-surface CG (M2.12c-2): the scoped "mesh-infra extension" was NOT needed; a cross-rank iterative solver byte-matches

The SSH stiffness assembly (`init_stiff_mat_ale`) + the preconditioned CG (`solve_ssh_ale`) byte-match FESOM2 per-rank
on pi dist_2 + dist_8 — `max|Δ|=0` on `d_eta` (30 records incl. the c-1 chain; `tools/run_stepdyn_gate_multirank.sh`,
substep 9 un-ignored). The FIRST multi-rank iterative solver, with cross-rank `MPI_Allreduce` dot-products. Lessons:

- **VERIFY a scoped invariant empirically before building the machinery it implies.** Both the HANDOFF and L32 scoped
  c-2 as needing "the first mesh-infra extension": extend `elem2D_nodes` to the halo + `exchange_elem(gradient_sca)` +
  `enforce_cw` on owned+eDim, because `init_stiff_mat_ale` "reads `elem2D_nodes`/`gradient_sca`/`zbar_e_bot` at the
  eDim-halo triangles of owned edges". **That hypothesis was WRONG.** A 20-line Python check over the `dist_N` files
  (for each OWNED edge, is `edge_tri(:,ed) ≤ myDim_elem2D`?) found **0** halo triangles on every rank — the owned rows
  assemble **FULLY LOCALLY**, no exchange, no infra extension. The tell I should have trusted first: **FESOM2's OWN
  `elem2D_nodes` AND `gradient_sca` are allocated owned-only** (`oce_mesh.F90:497` `(3,myDim_elem2D)`, `:2466`
  `(6,myDim_elem2D)`), yet its `init_stiff` reads them at `edge_tri(:,ed)` for owned `ed` — which is only safe (no OOB)
  if owned edges have owned triangles. The oracle's allocation IS the proof of the invariant. The three partition
  invariants (verified pi dist_2/8, universal in FESOM2 since the owned-only alloc is unconditional): (i) every edge
  incident to an owned node is owned ⇒ looping owned edges visits every contribution to an owned row; (ii) both
  triangles of an owned edge are owned ⇒ the element-array reads stay in the owned-only arrays; (iii) an owned
  element's nodes are within owned+halo (`≤ nNodL`, `maxlocnode==nNodL` exactly) ⇒ `n_num(elnodes)` in bounds and the
  CSR may have HALO columns. Don't infer "needs a halo" from "reads an element array at `edge_tri`"; check whether
  `edge_tri` of an OWNED edge ever leaves the owned set. (c-1's `compute_ssh_rhs_ale` looks identical and ALSO only
  reads owned elements for owned edges — the c-1 driver's `nElemF`-sized UV + `exchange_elem_full` was defensive
  over-provisioning driven by the owned+HALO-edge visc kernel, NOT by ssh_rhs.)
- **A cross-rank iterative solver byte-matches — the L29 corollary, confirmed at MR.** Same `A`/`b`/`x0`/`M⁻¹` +
  byte-identical per-iteration reductions ⇒ byte-identical recurrence. The only new MR risk is the reduction ORDER of
  the dot-products, and `MPI_Allreduce(MPI_SUM)` is byte-identical between the two codes: same OpenMPI (4.1.2-intel),
  same comm size, same op, same 8-byte type ⇒ the library picks the same reduction tree ⇒ same operation order (the
  same determinism FESOM2's own reproducibility rests on, L6). **38 CG iters on 8 ranks, `max|Δ|=0` on `d_eta`** — the
  L29 NOVECTOR precond divide carries over unchanged. So "it's a long cross-rank recurrence" is NOT a reproducibility
  excuse any more than "it's a long recurrence" was (L28→L29).
- **The lift is pure bounds+exchange (the L31/L32 recipe again).** `nod2D→nNodO` for owned loops; arrays
  `rr/zz/pp/App` + `diag_values` sized `nNodL`; `exchange_nod(diag_values)` in the precond (a halo node's diagonal
  lives on its owner — the off-diag `K_ri` reads it), `exchange_nod(pp)` before `A·p` and `exchange_nod(rr)` before
  `M⁻¹r` (halo COLUMNS of the mat-vec), `allreduce_sum` after each owned-partial-sum dot-product. The convergence/`rtol`
  denominators stay **GLOBAL** (`mesh%nod2D`, which the local-mesh remap sets to the global count). Keep the explicit
  `DO row; s=s+…` dot-product form (NOT `sum()`): the oracle is `ENABLE_OPENMP=OFF`/`__openmp_reproducible` undefined,
  so its `#if !defined(__openmp_reproducible)` reduction compiles to that serial loop (L16). Optional `partit` absent
  ⇒ the 1-rank path runs VERBATIM (no exchange, no allreduce: the local sum already IS the global sum).
- **All ranks reporting the SAME CG iter count is the non-vacuity + consistency check.** The exit test
  `sqrt(sprod(2)/nod2D) < rtol` uses the allreduce'd `sprod(2)` and the global `nod2D`, so every rank computes the same
  test and exits at the same iteration (38). A solver where ranks disagreed on the iter count would mean a non-global
  convergence test or a drifting reduction — a red flag even before checking `max|Δ|`.

## L34 — Multi-rank WHOLE STEP (M2.12c-3): the post-SSH tail lifts with the same recipe; threading an optional partit through the assembly closes the multi-rank MVP

The post-SSH ALE update (`update_vel`/`compute_hbar_ale`/`update_eta_n`/`vert_vel_ale`(+`compute_CFLz`/
`compute_Wvel_split`)/`update_thickness_ale`) + the tracer SOLVE (`solve_tracers_ale`/`diff_tracers_ale`/
`diff_part_hor_redi`/`diff_ver_part_impl_ale`) lifted to multi-rank with the M2.12b/c optional-`partit` recipe, and
threading the optional `partit` through `mod_step_oce::step_oce` made the WHOLE assembled ocean step byte-match FESOM2
per-rank on pi dist_2 + dist_8 — `max|Δ|=0` on all 65 substep records (`tools/run_step_gate_multirank.sh`). **The whole
multi-rank model is byte-identical to FESOM2 — the architectural MVP.** The reusable lessons:

- **Thread an optional through the assembly, not an `if(present)` ladder.** `step_oce` calls ~16 kernels; the lift is
  one `type(t_partit), intent(in), optional :: partit` on `step_oce` + `, partit` appended to every kernel call.
  Fortran passes a non-present optional actual as "absent" to the callee, so `call kernel(..., partit)` does the right
  thing whether `partit` is present or not — no branching, and the 1-rank callers (`fesom_stepdump`/`fesom_lifecycle`)
  that omit `partit` are untouched (their `step_oce` runs every kernel's 1-rank path VERBATIM). One keyword caveat:
  when the optional sits AFTER another optional in the callee (`solve_ssh_ale(dyn,mesh,n_iter,partit)`), call it by
  keyword — `call solve_ssh_ale(dyn, mesh, partit=partit)`.
- **The accumulate-then-exchange pattern dominates the post-SSH tail, and the one easy-to-miss exchange is the one
  OUTSIDE a linfs guard.** `compute_hbar_ale` skips the water_flux term AND its `exchange_nod(ssh_rhs_old)` for linfs —
  but the `exchange_nod(hbar)` right after sits OUTSIDE that guard and fires ALWAYS, because the next loop (`dhe` over
  owned elements) reads `hbar` at the element's 3 nodes, which span the halo. Read the FESOM2 control flow to the
  closing `endif`: an exchange that looks like it belongs to the gated branch may be unconditional. `update_eta_n`
  loops owned+halo (no exchange — both operands are already halo-valid); `vert_vel_ale` is owned-edge-scatter +
  owned-cumsum then `exchange_nod(w)`/`exchange_nod(hnode_new)`; `compute_CFLz`/`compute_Wvel_split` are owned+halo
  (so the NEXT step's momentum advection reads valid halo `w_e`/`w_i`); `diff_ver_part_impl_ale` is a per-column TDMA
  with no halo coupling (owned-node loop, no exchange); `solve_tracers_ale` does `exchange_nod(values)` after each
  tracer's solve + clamps salinity over owned+halo.
- **`exchange_elem` has no rank-3 variant — rank-3 `UV` uses `exchange_elem_full`, a SUPERSET of FESOM2's eDim
  `exchange_elem`, and that is safe because the owned dumps never depend on the extra halo.** FESOM2's `update_vel`
  ends with `exchange_elem(UV)` (eDim); FESOM3's `exchange_elem` interface is rank-1/rank-2 only, so the rank-3
  `(2,nl-1,elem)` UV must go through `exchange_elem_full` (eDim+eXDim). Before accepting the superset, confirm the
  gated fields depend only on OWNED values: every owned-edge kernel reads `UV` at the edge's triangles, which are
  owned (invariant ii), so the eDim/eXDim halo of `UV` is never read into an owned dump — the extra eXDim refresh is
  invisible to the gate (and strictly *more* correct for the next step's halo reads). Same logic retires
  `update_thickness_ale`'s `exchange_elem(helem)` for linfs: helem is unchanged AND already full-halo-valid, so the
  exchange is value-neutral — skip it (L33 "no machinery you don't need").
- **In the partitioned mesh, `mesh%nod2D`/`elem2D`/`edge2D` hold the GLOBAL counts — any work array a caller sizes
  from them must be re-sized to the LOCAL count.** `solve_tracers_ale` allocated `tr_xy(2,nl-1,mesh%elem2D)`; at
  multi-rank that is the global element count (wrong + huge). Size it `nElemF` (local, via `local_dims`). The kernels'
  explicit-shape dummies keep the global declared bound (preserving codegen / the L29 divide), but the ACTUAL
  allocation in the driver/caller must be local — the dummy maps onto the local storage and only owned/local indices
  are ever accessed.
- **The whole-step gate cost almost nothing in new infra because the oracle was already a whole-step driver.** The
  c-1 `fesom_step_dump` npes>1 extension drives the REAL `oce_timestep_ale` and its built-in dump_shim emits ALL
  substeps — so c-3 needed ZERO oracle change. Only the FESOM3 side needed work: the post-SSH lift + a whole-step MR
  driver (`fesom_stepfull_mr` = `step_oce` through `partit`, the multi-rank analog of the 1-rank `fesom_stepdump`).
  The gate is the 1-rank `run_step_gate.sh` ignore set (only SW_AB id=2) applied per-rank.
- **A byte-gate on a rich analytic state subsumes the "physical sanity" probes (rest-at-rest, gravity wave,
  `stale_halo_max_*`).** `max|Δ|=0` vs FESOM2 across all 65 substeps on non-trivial T/S/UV/SSH at 2 AND 8 ranks is a
  strictly stronger statement than any single physical scenario, and because we replicate FESOM2's exact owned loop
  bounds + exchanges, a stale halo cannot affect an owned dump — it would have to change an owned value to be
  observable, and it doesn't. The gate IS the stale-halo probe.
- **The whole multi-rank port was bounds + exchanges end-to-end — never arithmetic (L31→L34, confirmed at MODEL
  scale).** From geometry (c-a) through advection (c-b), the dynamics RHS (c-1), the iterative SSH solver (c-2), to the
  ALE update + tracer solve + the assembled step (c-3): not one arithmetic line changed, the codegen (and the L29
  vectorised divide) is preserved verbatim, and every gate was `max|Δ|=0` on the first or second run. The discipline
  that produced this: transcribe FESOM2's loop bounds + halo exchanges exactly, gate per-rank vs the same partition
  (L8), keep the optional-`partit`-absent path byte-for-byte the proven 1-rank code so nothing can regress, and verify
  a scoped partition invariant empirically before building machinery for it (L33).

## L35 — Production validation (multi-rank FREE-RUNNING lifecycle): the single-step prescribe-and-stop gate MASKS live-multi-step halo bugs; UV over-exchange FIXED, `tr_xy` "exchange" was a KNEM single-copy MPI bug — RESOLVED, `max|Δ|=0`

The m2-mvp byte-match was proven by SINGLE-STEP prescribe-and-stop gates (`run_step_gate_multirank.sh` etc.): the
whole state is prescribed at owned+halo each run, so a wrong/stale HALO is invisible — the prescribed halo is always
correct, and the step runs once. **That cannot validate halo MAINTENANCE across steps.** The production-validation
**multi-rank FREE-RUNNING lifecycle** (`fesom_lifecycle_mr` + `run_lifecycle_gate_multirank.sh`, CORE2 dist_2,
dt=1800 — the multi-rank analog of the byte-exact 1-rank L29 lifecycle) was the FIRST test of *multi-rank ×
free-running × many-steps*, and it immediately exposed two latent bugs the single-step gates masked:

- **UV halo over-exchange — FIXED.** `update_vel` exchanged UV over `com_elem2D_full` (eDim+eXDim); FESOM2 uses
  `com_elem2D` (eDim, `exchange_elem(UV)`). The extra eXDim fill diverged from FESOM2's stale-0 eXDim and (compounded
  by the same exchange pathology as the open bug below) corrupted the UV halo → viscosity blew it up (eta→2708 by
  step 3) while FESOM2 stayed at eta≈0.54. Fix: match FESOM2's eDim exchange. Model then STABLE + matches FESOM2 to
  all printed digits. **Lesson: match FESOM2's exact halo WIDTH per field — UV is eDim, `tr_xy` is full; the extra
  halo is not "harmless superset" when a downstream kernel reads it.**

- **`tr_xy` halo "drift" was an OpenMPI `vader` KNEM single-copy bug — RESOLVED.** The element tracer-gradient
  `tr_xy`, exchanged over `com_elem2D_full` before the MUSCL `fill_up_dn_grad`, got a WRONG halo at step 2+ (OWNED
  byte-exact; HALO 456/636 wrong, F2≈1e-6 vs F3≈1e-11), drifting T/S → gross step-3 `ssh_rhs`. The exchange LOGIC was
  provably correct, and it was — the corruption was BELOW our code, in the MPI transport. **Bisection (cross-rank
  buffer trace in `core_blk_r`):** rank1's packed `sbuf` was correct, rank0's unpack `rbuf→arr` was faithful, but the
  bytes *in transit* `sbuf→rbuf` were correct for the first **134656 bytes** (16832 reals) then garbage — a contiguous
  truncation, not a shift or a logic error. **Root cause:** OpenMPI 4.1.2's `vader` (shared-memory) BTL uses a
  single-copy mechanism for large messages; on levante it falls back to **KNEM** (CMA needs ptrace, blocked by
  `kernel.yama.ptrace_scope=3`), and KNEM here CORRUPTS messages above the eager limit. Proof: same binary, same
  buffers — under `--mca btl self,tcp` `total-diff=0`; under `--mca btl_vader_single_copy_mechanism none` (forces the
  copy-in/copy-out two-copy path) `total-diff=0`; only default `vader`+KNEM corrupts. (Mechanism found in L61: the knem registration cache misses glibc heap trims -> stale pages.) **Fix:** `env.sh` exports
  `OMPI_MCA_btl_vader_single_copy_mechanism=none` (levante branch) → the CORE2 dist_2 free-running lifecycle is
  `max|Δ|=0` (195 records), with no 1-rank/pi regression (ctest 13/13, step gates 65/65 `=0`). RULED OUT (correctly, in
  the end): eXDim halo ordering (`myList_elem2D` is read VERBATIM from the same `my_list*.out` as FESOM2 — identical
  owned AND halo), `core_blk_r` impl, com corruption, the gid-test (it was telling the truth — our logic is right).
  **Why it hid until now:** the bug needs a message >~131 KB; every pi-mesh gate's halos stay well under it, and FESOM2's
  `MPI_TYPE_INDEXED` exchange uses a different `vader` path that dodges KNEM. CORE2's `tr_xy` block exchange (94×636 =
  478 KB) is the FIRST live F3 message above the threshold — only the big-mesh free-running lifecycle could expose it.

**Meta-lessons:** (1) a single-step prescribe-and-stop byte-gate is necessary but NOT sufficient — only a free-running
multi-step gate validates halo maintenance; and only a BIG-MESH one exercises the message sizes that trip a buggy MPI
single-copy path. (2) When OWNED is byte-exact and the exchange logic is provably correct yet the HALO is still wrong,
suspect the LAYER BELOW your code (the MPI transport): bisect pack→transport→unpack with a cross-rank buffer trace
keyed by gid; here `sbuf`(correct)→`rbuf`(corrupt)→`arr`(faithful) pointed straight at MPI. (3) "Provably impossible by
every verified fact" usually means a fact lives one layer down — verify the transport itself, not just your indices.
(4) The clean gid-exchange self-test (`g3(owned)=gid`) verifies the slist/rlist MAPPING and was NOT misleading — it
correctly said the logic was right; trust it for what it tests and look elsewhere (not "it must be a subtle logic bug
the test can't see"). (5) `myList_elem2D` (owned AND halo) is byte-identical to FESOM2 by construction — both read the
same `dist_N/my_list*.out` verbatim — so a position-keyed dump diff IS gid-aligned; don't chase an "ordering" ghost.

## L36 — Sea-ice EVP dynamics (M3b): ocean2ice + standard EVP + mEVP byte-match at once; a new `t_mesh` component trips an ifort type-bound-I/O cascade (put ice fields in `t_ice`)

M3b ported `ocean2ice` + the EVP momentum solve and byte-matched FESOM2 `max|Δ|=0` on 7 fields
(`ice_srfoce_u/v` + `uice/vice` + `sigma11/12/22`), CORE2 1-rank, for **BOTH `whichEVP=0` (standard EVP) AND
`whichEVP=1` (mEVP)** — `tools/run_evp_gate_core2.sh`, both modes in one run. New `src/ice/mod_ice_dyn.F90`
(`ocean2ice`, `EVPdynamics` [std], `EVPdynamics_m` [mEVP], `EVPdynamics_solve` dispatcher), same optional-`partit`
pattern as the M2.12 kernels (1-rank verbatim; `is_multirank`-guarded `exchange_nod`). Gate driver
`src/drivers/fesom_evpdump.F90` prescribes analytic surface UV / `hbar` / `stress_atmice` (= the oracle, from the
byte-identical rotated `coord_nod2D`), runs the coupling + EVP, dumps; oracle is `fesom_ice_dump.F90::evp_dump_write`
(env `FESOM_EVP_DUMP`, dispatches on namelist `whichEVP`, runs the REAL FESOM2 routines).

- **The big gotcha (reusable): do NOT add a component to `t_mesh`.** FESOM3's `t_mesh` has type-bound
  `write(unformatted)`/`read(unformatted)` (`write_t_mesh` calls `write_bin_array` on each component). Adding ANY new
  allocatable component (tried `bc_index_nod2D`, both `WP` and `MP`) makes ifort emit error **#7976** ("allocatable
  dummy may only be associated with an allocatable actual") for EVERY `write_bin_array` call in `write_t_mesh` — a
  generic-resolution/type-bound-I/O cascade, not a real error in the added line. **Fix:** put ice-owned mesh-derived
  fields in `t_ice` (which has NO type-bound I/O), not `t_mesh`. FESOM2 keeps `bc_index_nod2D` in `mesh%`; we keep it
  in `ice%` — the VALUE is byte-identical (0/1 boundary-node mask), so the byte-gate is unaffected. Storage location is
  a free FESOM3 design choice; faithful = same value + same arithmetic, not same struct.
- **mEVP specifics** (`ice_maEVP.F90` `EVPdynamics_m`, the INLINED form — NOT the standalone `stress_tensor_m`/
  `stress2rhs_m`, which are unused in that path): `rdt = ice_dt` (the FULL step, NOT `ice_dt/evp_rheol_steps` like
  std-EVP); `det2=1/(1+alpha_evp)`, `det1=alpha_evp*det2`; iterates an AUXILIARY velocity `uice_aux/vice_aux` (init =
  `uice/vice`, copied back at the end); velocity BC is the node mask `bc_index_nod2D(i)` (in `det`) PLUS the edge
  zeroing; `pressure = pressure_fac/(delta+delta_min)` with `pressure_fac=det2*pstar*msum*exp(...)`. std-EVP and mEVP
  group the metric/elevation terms differently (`metric_factor*sum/3` vs `sum*val3*metric_factor`; `9.81_WP*A/3` vs
  `g*val3*A`) — transcribe EACH verbatim against ITS own FESOM2 source; do NOT unify them.
- **Namelist-double precision (M3a note paid off):** `cd_oce_ice`/`delta_min` are read by the oracle's namelist as
  DOUBLES (`0.0055`/`1.0e-11`), but the `t_ice` defaults are un-suffixed single→WP (`5.5e-3`/`1.0e-11`) — ~1 ULP off.
  The driver overrides `ice%cd_oce_ice=0.0055_WP`/`ice%delta_min=1.0e-11_WP`. `alpha_evp`/`beta_evp`=250, `pstar`=30000,
  `ellipse`=2, `c_pressure`=20, `theta_io`=0 are integer-valued → exactly representable → no override needed.
- `bc_index_nod2D` + `uice_aux`/`vice_aux` (both M3a-deferred) built in `mod_ice_setup` (`ice_allocate` +
  `build_bc_index_nod2D`, unconditional like FESOM2's "also for whichEVP==0"). `ice_setup` `mesh` stays `intent(in)`
  (the mask is in `ice`, not `mesh`). **No regression** (step 65 + M3a ice 4 + 13/13 ctest all `max|Δ|=0`).
- **Meta:** when the user asks for an extra variant (here mEVP), porting both variants of a dispatched kernel at once
  is cheap — they share `ocean2ice`, the prescribe driver, the oracle shim, and the gate harness; only the inner solve
  differs. Gate both by parameterising `whichEVP` end-to-end (driver env + oracle namelist patch + one runner loop).

## L37 — Sea-ice FCT advection (M3c): the whole Zalesak limiter transcribes verbatim once the neighbour-list basis is pinned; an M3a invariant pays its dividend downstream

M3c ported the ice FCT advection (`ice_TG_rhs` Taylor-Galerkin rhs + `ice_fct_solve` = `ice_solve_high_order` 3-iter
Jacobi → `ice_solve_low_order` → `ice_fem_fct` Zalesak limiter ×3 tracers) and byte-matched FESOM2 `max|Δ|=0` on **6
fields** (`rhs_a/m/ms` + post-advection `a_ice/m_ice/m_snow`), CORE2 1-rank, for **BOTH `whichEVP=0` AND `whichEVP=1`**
feeding the advection — `tools/run_icefct_gate_core2.sh`, **first try, no iteration.** New `src/ice/mod_ice_fct.F90` +
driver `src/drivers/fesom_icefctdump.F90` (extends `fesom_evpdump`: ocean2ice + EVP → `uice/vice`, then `ice_TG_rhs` +
`ice_fct_solve`) + oracle `fesom_ice_dump.F90::ice_fct_dump_write` (env `FESOM_FCT_DUMP`).

- **The key enabler — an upstream invariant, re-used not re-derived.** FESOM2's FCT indexes the (row, neighbour)
  mass-matrix slot via `mesh%nn_num`/`nn_pos` (`mass_matrix(clo:clo2)` aligned with `nn_pos(1:cn,row)`). FESOM3 has its
  OWN `mesh%nn_pos` (built by `muscl_adv_init`) but it may order neighbours differently — so DON'T use it. Instead use
  the **ssh_stiff CSR** (`rowptr_loc`/`colind_loc`), the same ordered list `ice_mass_matrix_fill` scanned to build
  `fct_massmatrix`. M3a PROVED `fct_massmatrix` byte-matches FESOM2's `mass_matrix` position-by-position ⇒
  `colind_loc(rowptr_loc(row):rowptr_loc(row+1)-1)` IS FESOM2's `nn_pos(1:nn_num,row)` as an ORDERED list (incl. self).
  So `sum(fct_massmatrix(clo:clo2)*field(colind_loc(clo:clo2)))` reproduces FESOM2's `sum(mass_matrix(clo:clo2)*
  field(nn_pos))` with the SAME operands AND the SAME reduction order → byte-identical. The cluster min/max (`maxval`/
  `minval` over the neighbour SET) only needs the same set, so the CSR basis serves there too. **Meta: when a kernel
  needs a neighbour list, reach for the basis a PRIOR gate already pinned (the CSR), not a parallel structure with its
  own ordering — the FCT solve then depends only on M3a, and inherits its byte-proof for free.**
- **`ice_diff=0.0` / `ice_gamma_fct=0.5`** are the CORE2 `namelist.ice` values, NOT the `t_ice` defaults (10.0 / 0.25);
  the driver overrides `ice%ice_diff`/`ice%ice_gamma_fct` (the FCT reads them). Both are exactly representable (unlike
  `cd_oce_ice=0.0055`) so no ULP subtlety — but they MUST be the namelist values or the diffusion/limiter blend differs
  by a large margin. `ice_diff=0` zeroes the TG diffusion term, but `scale_area=2.0e8` is still kept correct so
  `sqrt(elem_area/scale_area)` is finite (else `0*NaN=NaN`).
- **`ENABLE_OPENMP=OFF` in the oracle (verified in `CMakeCache.txt`)** ⇒ `ice_fem_fct`'s `!$OMP ORDERED`/`ATOMIC` flux
  scatters run serially in element order = FESOM3's serial path, so the node-accumulation FP order matches (same as
  M3b's `stress2rhs`). Whenever a ported kernel has a scatter-accumulate, confirm the oracle's OpenMP state — a
  threaded non-`__openmp_reproducible` build would NOT byte-match.
- **Dump rhs AND the final tracers** (the pressure-gate "dump-inputs-too" habit): `rhs_a/m/ms` localise a `ice_TG_rhs`
  bug, the post-FCT `a_ice/m_ice/m_snow` localise an `ice_fct_solve` bug. Both `=0` ⇒ both kernels independently
  confirmed. The advection is non-trivial (`a_ice` 0.900 → 0.913 — the antidiffusive flux pushed past the IC max, so
  the limiter clamp path is genuinely exercised). `ice_fct_solve` is called as an external (no interface) exactly like
  `ice_timestep` (it is NOT in `ice_fct_interfaces`); `ice_TG_rhs` IS, so `use ice_fct_interfaces, only: ice_TG_rhs`.
- **No regression:** step 65 + M3a ice 4 + M3b evp 7×2 + M3c (6×2 whichEVP) + 13/13 ctest all `max|Δ|=0`/green (the
  M3b evp gate was re-run explicitly since the shared oracle `fesom_ice_dump.F90` was edited; the byte-exact `rhs` also
  independently re-confirms M3b — the `uice/vice` feeding `ice_TG_rhs` were byte-identical).

## L38 — Sea-ice thermodynamics (M3d): a Newton-iteration growth model byte-matches verbatim; the win is the namelist-DOUBLE precision discipline + all-scalar arithmetic (no L29 trap)

M3d ported the ice thermodynamics (`cut_off` hmin/Armin clamp + `thermodynamics` per-node driver + `therm_ice` 0-layer
Semtner/Hibler-1984 ice-class growth + `budget` 5-iter Newton-Raphson surface temperature + `obudget` open-ocean
growth+evaporation + `flooding` snow→ice + `TFrez`) and byte-matched FESOM2 `max|Δ|=0` on **6 fields** (post-thermo
`a_ice/m_ice/m_snow` + `flx_h/flx_fw` + `t_skin`), CORE2 1-rank, for **BOTH `whichEVP=0` AND `whichEVP=1`** feeding the
advection upstream — `tools/run_icethermo_gate_core2.sh`, **first try, no iteration.** The dumped state is rich
(`flx_h` to ~6 kW/m², `t_skin` to ~15 °C, surface T −1.9→30 °C ⇒ freezing+melting+snow-melt+flooding all exercised). New
`src/ice/mod_ice_thermo.F90`; driver `fesom_icethermodump.F90` (extends `fesom_icefctdump`); oracle
`fesom_ice_dump.F90::ice_thermo_dump_write` (env `FESOM_THERMO_DUMP`). Reusable specifics:

- **A 5-iteration Newton-Raphson is not a reproducibility risk when its operands are byte-pinned** (cf. the CG L29
  scare). `budget` iterates `t = t + (A1+A2+C)/A3` 5× per ice class, and the per-class `t` even *carries across* the
  iclasses loop (a faithful quirk: each `budget` call mutates the shared skin-temp `t`). Byte-matched on the FIRST gate
  because every operand is pinned (prescribed atmosphere + M3b/M3c-proven ice state + do_ic3d-proven `srfoce_temp/salt`
  + geometry-proven `geo_coord`) and the whole subtree is **SCALAR per-node arithmetic** — there is no array loop to
  vectorise, so the L29 packed-vs-scalar divide trap *cannot* arise. Transcribe verbatim; it just works (the L9 pattern
  at full thermodynamic complexity).

- **The `&ice_therm` namelist is the precision battleground (extends the M3a note).** The oracle's `ice_init` reads
  `&ice_therm` → parses each literal as a **DOUBLE**, overwriting the `t_ice_thermo` single→WP type defaults. So
  `con=2.1656`, `consn=0.31`, `hmin=Armin=0.01`, `emiss_ice=emiss_wat=0.97`, the four albedos, and `albw=0.1` (CORE2,
  *not* the 0.066 LY2004 default) all differ ~1 ULP between the namelist-double and the single→WP default — the driver
  MUST re-set them as `_WP` doubles (the namelist-parse == `_WP`-literal identity, proven through M3a–M3c, makes this
  byte-exact). Same for `Ch_atm_ice=Ce_atm_ice=0.00175_WP` (namelist.forcing `&forcing_exchange_coeff` doubles).
  Params NOT in the namelist match for free: `h_ml=2.5_WP` (suffixed both sides); `rho*/inv_rho*/clh*/tmelt/boltzmann/
  cpair` (MOD_ICE single→WP defaults, untouched); and `cc=rhowat·4190 / cl=rhoice·3.34e5` recompute in `ice_init` to
  values that are *exactly representable* (4294750 / 303940000 are integers below the round-off) so the double recompute
  == the single→WP literal product. Recipe: override every namelist-listed param in the driver; trust the rest.

- **Config switches that change the math must be matched, not just the literals.** `which_ALE='linfs'` ⇒
  `use_virt_salt=.true.` (the `(rsss-Sice)/rsss` freshwater branch + the flooding `+iflice·…/rsss` term), and
  `ref_sss_local=.true.` (namelist.tra) ⇒ `rsss = S_oc` (the per-node surface salinity, *not* the global `ref_sss`).
  `l_snow=.true.` ⇒ `rain=prec_rain, snow=prec_snow, evap_in=0`. Get any of these wrong and `flx_fw` diverges while the
  heat-only fields still match — gate every output (we dump `flx_h` AND `flx_fw` AND `t_skin`) so the failure localises.

- **FESOM3 has no global forcing state — bundle it.** FESOM2's `thermodynamics` reads/writes ~22 `g_forcing_arrays`/
  `o_arrays` module arrays. The FESOM3 port passes them as ONE `t_atmflux` argument (per-node inputs + the thermo flux
  diagnostics + the scalar config), keeping the explicit-dataflow architecture; the GATED outputs
  (`a/m_ice/m_snow/flx_h/flx_fw/t_skin`) stay in `t_ice`. The diagnostic outputs (`real_salt_flux/fw_ice/fw_snw/hf_Q*`)
  are write-only here but become M3e's `oce_fluxes` inputs — `t_atmflux` is their natural home, not throwaway scratch.

- **Dead branches still need their symbols.** `obudget`'s `open_water_albedo>0` solar-zenith block is dead in the gated
  config (=0) but must link — port `compute_solar_zenith_angle`/`albw_taylor`/`albw_briegleb` anyway, and supply local
  `daynew/timenew` placeholders for the (never-taken) `g_clock` reads. (Also: the oracle pointer-assigns into an
  `intent(in)` `ithermp%albw` in that block — ifort allows it; the FESOM3 transcription mirrors it, dead, harmless.)

- **No regression:** step 65 + M3c (6×2 whichEVP) + 13/13 ctest all `max|Δ|=0`/green (the M3c gate re-run since the
  shared oracle `fesom_ice_dump.F90` was edited — adding `ice_thermo_dump_write` left the existing shims byte-identical).

## L39 — Sea-ice → ocean coupling-out (M3e, THE PAYOFF): the air-sea budget byte-matches first try; a precomputed reduction scalar that was never consumed is an un-gated landmine (the L29 pattern, pre-empted)

M3e ported the air-sea coupling-out (`oce_fluxes_mom` ice-ocean+atm-ocean momentum stress → `stress_surf`, and
`oce_fluxes` heat/freshwater/salt budget → `heat_flux`/`water_flux`/`virtual_salt`/`relax_salt`) from
`ice_oce_coupling.F90` and byte-matched FESOM2 `max|Δ|=0` on **5 fields**, CORE2 1-rank, for **BOTH `whichEVP=0` AND
`whichEVP=1`** — `tools/run_iceflux_gate_core2.sh`, **first try.** These 5 fields ARE the proven M2.11c-2
`fesom_flux_dump` set: until M3e the ocean step READ them from a FESOM2 dump (the "M3-gap" prescription); now FESOM3
produces them natively. New `src/ice/mod_ice_oce_coupling.F90`; driver `fesom_icefluxdump.F90` (extends
`fesom_icethermodump`); oracle `fesom_ice_dump.F90::ice_flux_dump_write` (env `FESOM_OCEFLUX_DUMP`). Reusable specifics:

- **A precomputed scalar that has NEVER been consumed is an un-gated landmine — the L29 pattern, caught BEFORE it bit.**
  `mesh%ocean_area` is the divisor in every flux-balancing step (`net = integrate_nod(field)/ocean_area`). FESOM3 had
  computed it as `sum(mesh%area(1,1:nNodO))` since M2.11a — never gated, because nothing read it until M3e. FESOM2
  computes it as an explicit *sequential* `do`-loop over `areasvol(ulevels,n)` (`oce_mesh.F90:2385`). A `sum()`
  intrinsic and a sequential loop are NOT guaranteed bit-identical even under `-fp-model precise` (precise stops
  *reassociation*, but the intrinsic's base evaluation order is processor-dependent), and a 1-ULP `ocean_area` would
  drift EVERY balanced field. Fix: rewrite `mod_mesh_areas.F90`'s `ocean_area`/`ocean_areawithcav` as the faithful
  FESOM2 loop. **Generalisable: when a long-dormant precomputed value becomes load-bearing, re-derive it against the
  oracle's exact arithmetic before trusting it — "computed but unconsumed" == "ungated" (cf. L29's `pr_values`).**

- **`integrate_nod` is an explicit sequential loop, never `sum()`/`dot_product()`.** Transcribe FESOM2
  `gen_support.F90:318` verbatim: `lval=0; do row=1,nNodO: lval=lval+data(row)*areasvol(ulevels_nod2D(row),row)`, then
  (multi-rank) `allreduce_sum`. The summation order is the byte-match — a reduction intrinsic would invite the compiler
  to reorder. At 1-rank the allreduce is identity, so the local partial sum IS the result.

- **`oce_fluxes` is bookkeeping, not physics — the byte-match risk is concentrated, not diffuse.** Strip the disabled
  features (no `__icepack`/`use_cavity`/`use_icebergs`/`lwiso`/`use_landice_water`/`use_age_tracer`/`__oasis`) and the
  ~600-line routine collapses to ~80: sign-flips (`heat_flux=-flx_h`), two globally-balanced salt fluxes
  (`virtual_salt`, `relax_salt`), and a globally-balanced freshwater flux. The ONLY ULP-sensitive operations are the
  three `integrate_nod` calls and the `/ocean_area` divides — every per-node summand is a byte-pinned M3d/do_ic3d
  product. Pin the global integrals and the rest is exact. Gate the `use_virt_salt=.true.` (linfs) path; `ref_sss_local`
  ⇒ `rsss=salt(1,n)`; the additive order of the freshwater `flux(n)` is transcribed left-to-right verbatim.

- **dens_flux (the MOC diagnostic, `oce_fluxes:691`) is the one piece you DON'T port at M3e.** It needs `sw_alpha`/
  `sw_beta` (EOS) + `vcpw` and feeds only the MOC diagnostic — neither in M3 scope, and it is NOT a surface BC. Skipping
  it does not touch the 5 gated fields (the oracle computes it into a zeroed `sw_alpha`/`sw_beta` ⇒ `dens_flux=0`, never
  dumped). Don't thread EOS state through M3e just to reproduce a deferred diagnostic.

- **Same precision/Intel discipline as M3d.** `surf_relax_S=1.929e-06` is a `namelist.tra` DOUBLE (NOT the o_PARAM
  default `10/(60*3600*24)`); driver re-sets `1.929e-06_WP`. `density_0=1030.0_WP` (mod_constants == o_PARAM, exact). The
  new analytic prescriptions (`stress_atmoce_x/y`, `Ssurf`) must be byte-identical functions on both sides. And (L27
  again) declare the `mesh` dummy BEFORE `stress_surf(2,mesh%elem2D)` or ifort #6415 aborts — the only iteration this
  gate needed was that one declaration reorder.

- **No regression:** M3e 5×2 + M3d 6×2 (shared oracle re-run) + step 65 (1-rank) + step 65 (multirank dist_2,
  confirming the `mod_mesh_areas` change is geometry-neutral) + pressure 57 + 13/13 ctest all `max|Δ|=0`/green.

## L40 — Native-flux coupled forced lifecycle (M3f-1 + M3f-2): the multi-step ice↔ocean coupling byte-matches first try; isolate the forcing READ from the ice/flux COUPLING by prescribing the post-bulk atmosphere

M3f-1 assembled `ice_timestep` (`src/ice/mod_ice_step.F90` = the FESOM2 `ice_setup_step.F90:96` chain
`EVPdynamics_solve → ice_TG_rhs → ice_fct_solve → cut_off → thermodynamics`) and M3f-2 wired it into a NATIVE-flux
forced lifecycle (`src/drivers/fesom_lifecycle_native.F90`): the proven M2.11c-2 ocean init + `step_oce`, but the
prescribed `FESOM3_FLUX_FILE` (heat_flux/water_flux/virtual_salt/relax_salt/stress_surf) REPLACED by the live chain
`ocean2ice → ice_timestep → oce_fluxes_mom → oce_fluxes`. **`max|Δ|=0` on the 195-record (3-step) AND 325-record
(5-step) NODE substep set, BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP) — first try, no iteration.** The chain is
the M3a–M3e leaf kernels, all already byte-proven single-step; M3f-2 is the FIRST multi-step + coupled test.

- **Prescribe the post-bulk ATMOSPHERE, not the fluxes — isolate the ice/flux COUPLING from the forcing READ.** The
  established prescribe-the-input discipline applied at the next layer up: instead of prescribing the 5 surface fluxes
  (M2.11c-2), prescribe the 16 per-step POST-`update_atm_forcing` atmospheric arrays (shortwave/longwave/Tair/shum/
  prec_rain/prec_snow/runoff/u_wind/v_wind/Ch-Ce_atm_oce_arr/stress_atmoce_x-y/stress_atmice_x-y/Ssurf) and let FESOM3
  compute the fluxes NATIVELY. New oracle shim `port2/fesom2/src/fesom_atmflux_dump.F90` (env `FESOM_ATMFLUX_DUMP`,
  wired in `fesom_module.F90` next to `flux_dump_record`, i.e. after the ice step) dumps them; FESOM3 reads them and
  runs the native chain. This makes the byte-gate localize cleanly: a forcing-read bug (M3f-3) is now SEPARATE from a
  coupling bug (M3f-2). Module homes: g_forcing_arrays (the 11 atm fields), o_ARRAYS (stress_atmoce_x/y, Ssurf),
  `ice%` (stress_atmice_x/y).
- **A coupled fixed-point byte-matches when every link does.** The ice reads the LIVE ocean (ocean2ice reads
  `dyn%uv(:,1,:)` node-avg + `tracers(1:2,1)` + `mesh%hbar`), produces fluxes, the ocean step reads the fluxes and
  updates the ocean, the ice reads the updated ocean next step… Step 1: both sides cold-start identical (ocean IC =
  do_ic3d, ice IC = SST-sign cold start, UV=0). Each link is byte-proven (ocean state via M2.11c, ice/flux via
  M3a–e), so the whole loop stays byte-locked across steps. The 5-step run (worst |Δ|=0) confirms it does not drift.
- **A built-in per-step native-vs-oracle flux self-check is the cheapest localizer.** The driver optionally reads the
  M2.11c-2 `flux_f2` dump (`FESOM3_FLUX_FILE`) and prints `max|Δ(hf,wf,vs,rs,ss)|` each step BEFORE `step_oce`
  consumes them. `=0` at every step pinpoints the coupling as byte-exact independently of the 195-record substep gate
  (and it immediately diagnosed the mEVP harness bug below). When a coupled gate can compare an intermediate against a
  known-good oracle dump, do it inline — don't wait for the end-to-end gate to fail.
- **The prognostic ice state IS the new test.** `values_old` (a_ice/m_ice/m_snow saved INSIDE `thermodynamics` before
  the update — `mod_ice_thermo.F90:250`), the EVP `sigma11/12/22` elastic memory (carries across model steps, was 0 at
  the single-step gate), the `t_skin` Newton initial guess, and `uice/vice` all persist in the live `ice` object. The
  multi-step gate is the first to exercise their carry-over; first-try `=0` means the single-step transcriptions were
  faithful to the carry-over too.
- **A fixed-namelist oracle must be parameterized to match a FESOM3 switch.** The forced-lifecycle oracle
  (`run_lifecycle_forced_core2.sh`) runs whatever `whichEVP` its `namelist.ice` holds (=0). FESOM3 `whichEVP=1` then
  "mismatched" — NOT a port bug, the oracle was still running std-EVP. The self-check flagged it instantly (heat_flux
  off by 1824 W/m² at step 1). Fix: the runner now `sed`-patches `whichEVP` in `namelist.ice` to a passed arg, the
  gate threads it to BOTH sides → mEVP also `max|Δ|=0`. Lesson: when gating a dispatched kernel against a namelist-
  driven oracle, parameterize the oracle's namelist from the same switch — never assume the oracle default matches.
- **`ice_timestep` is sequencing only — omit the unconsumed diagnostics.** The CMIP6 dynamical-growth-rate diagnostics
  (`ice%thermo%dyngr*`, `ice_setup_step:203`/`:307`) and the post-thermo `h_ice/h_snow` effective thicknesses
  (`:320`) are write-only output; nothing in the reduced config reads them, so omitting them from the assembled
  routine is byte-neutral for the gated ocean substeps (the gate proved it). Faithful = same consumed dataflow, not
  every diagnostic line.
- **No regression:** step 65 (1-rank) + ctest 13/13 + M2.11c-2 prescribed-flux forced lifecycle 195 + M3e flux 5×2 all
  `max|Δ|=0`/green (the shared oracle `fesom_module.F90` + `fesom_atmflux_dump.F90` edits are byte-neutral: the
  `atmflux_dump_record` call is env-gated + npes==1, so the prescribed-flux path is unchanged).

## L41 — Native CORE2 forcing read + bulk (M3f-3a/c): a proven READ machinery re-validates on a 40× mesh by getting the per-step clock right; validate a live-input kernel by self-check inside the coupled run

M3f-3a/c made 14/16 atmospheric arrays NATIVE in the coupled lifecycle (`max|Δ|=0`, both whichEVP): the 8 NCAR fields
(M2.10a `mod_forcing_read`) + the NCAR bulk Ch/Ce (M2.10b `forcing_bulk_ncar`) + `stress_atmoce` (M2.10b
`forcing_wind_stress`) + the NEW wind-on-ice `stress_atmice` (`forcing_ice_stress`). Only runoff + Ssurf (the monthly
climatology) stay prescribed (M3f-3b). Driver `fesom_forcing_core2.F90` (standalone NCAR-read self-check vs the M3f-2
`atmflux_f2` dump) + `fesom_lifecycle_native` extended with a `FESOM3_FORCING_DIR` per-step native-forcing self-check.

- **A byte-proven READ machinery re-validates on a new mesh for free — the only new thing is the per-step clock.**
  M2.10a's `mod_forcing_read` was gated on pi against the SAME CORE2 NCAR stub files; the CORE2 read differs ONLY in
  `geo_coord_nod2D` (where the bilinear lands) — the netCDF read, time-axis transform, bilinear, and g2r rotation are
  identical. So M3f-3a is just: drive the proven routines on CORE2 + advance `rdate` per step. The load-bearing detail
  is the clock: FESOM2 `clock()` (`gen_modules_clock.F90:38`) adds `dt` at the TOP of each step (`fesom_module:673`),
  so at step n `timenew=n·dt`, `daynew=1` (n<48); cold-start `getcoeffld` uses `rdate_cold=julday(1948,1,1,noleap)`
  (NO half-step, `nc_sbc_ini:643`), per-step `timeinterp` uses `rdate(n)=julday+(daynew-1)+timenew/86400 - dt/2/86400`
  (the `-dt/2` half-step, `sbc_do:1528`) → `rdate(n)=710820+(2n-1)·900/86400`. For a short run (within forcing day 1)
  the `getcoeffld` re-trigger (`sbc_do:1561`, `rdate>nc_time(t_indx_p1)`) never fires, so the cold-start coefficients
  apply every step — `getcoeffld` once + `timeinterp` per step is byte-faithful. (Multi-day runs crossing a forcing
  interval need the re-trigger — an M3f-3 longer-run deferral.) First try `max|Δ|=0`, 3 steps.
- **Validate a kernel with LIVE inputs by a self-check INSIDE the coupled run — don't try to gate it standalone.** The
  NCAR bulk needs the ocean surface state (`ice%srfoce_temp/u/v`, set by `ocean2ice`) and `stress_atmice` needs the
  previous-step `ice%uice/vice` — neither is in any dump, so a standalone gate can't feed them. Instead: in the
  lifecycle, after `ocean2ice` (srfoce live) and before `ice_timestep` (uice still the prev step's — exactly when
  FESOM2 `update_atm_forcing` runs), recompute the native forcing and print `max|Δ|` vs the prescribed atm. The step
  itself keeps USING the prescribed values (the proven M3f-2 path), so the 195-record gate CANNOT regress from a native
  bug — the self-check independently proves byte-equality. The oracle's `ncar_ocean_fluxes_mode`
  (`gen_bulk_formulae.F90:170`) reads `ice%srfoce_u/v/temp`; feeding the native bulk the same live arrays → `max|Δ|=0`
  on Ch/Ce/stress, both EVP variants (the mEVP `uice` flows into `stress_atmice` byte-exactly too).
- **`stress_atmice` is the one genuinely-new forcing arithmetic** (`gen_forcing_couple.F90:759`): a single combined
  FESOM2 loop sets BOTH `stress_atmoce` (drag `Cd_atm_oce_arr` from the bulk, wind relative to ocean via `Swind=0`) and
  `stress_atmice` (CONSTANT `Cd_atm_ice=0.0012` namelist drag — `AOMIP_drag_coeff=.false.` so no `cal_wind_drag_coeff`
  — wind relative to ice). Splitting it into `forcing_wind_stress` (M2.10b, oce) + new `forcing_ice_stress` is
  byte-identical (per-node independent arithmetic); each carries the `ulevels>1` cavity guard the combined loop had.
- **No regression:** M2.10 forcing gate (pi, `mod_forcing_bulk` recompiled) PASS + ctest 13/13 + the 195-record native
  lifecycle gate MATCH (both whichEVP) all `max|Δ|=0`/green.

## L42 — Runoff + SSS climatology read (M3f-3b) → the FULLY-NATIVE lifecycle (M3f-3 complete): a distinct read machinery byte-matches first try; the standalone-driver `partit%myDim=0` trap makes an optional-partit routine SILENTLY output all-zeros (no crash)

M3f-3b ported the last 2/16 atmospheric arrays — runoff + sea-surface-salinity restoring (`Ssurf`) — then dropped the
prescribed atmosphere entirely: the FULLY-NATIVE CORE2 lifecycle (whole air-sea forcing computed in-driver: 8 NCAR +
NCAR bulk + 2 stresses + runoff + Ssurf → ocean2ice → ice_timestep → oce_fluxes → step_oce) is now `max|Δ|=0` vs FESOM2
— **195 records (3-step) AND 325 records (5-step), BOTH whichEVP=0 and whichEVP=1**. New `src/forcing/mod_forcing_other.F90`
(`read_other_NetCDF` + `interp_2d_field`, from FESOM2 `gen_modules_read_NetCDF.F90:6` / `gen_interpolation.F90:145`) +
two real64 netCDF wrappers (`nc_get_slice_dp`/`nc_get_att_dp` in `mod_io_netcdf`); new gate
`tools/run_lifecycle_fullynative_gate_core2.sh` (NO `FESOM3_ATMFLUX_FILE`). FESOM2 oracle UNCHANGED (reuses the M3f-2
`atmflux_f2` dump for the self-check + the `run_lifecycle_forced_core2.sh` 195-record reference).

- **The runoff/SSS read is a DISTINCT machinery — do NOT shoehorn it into `mod_forcing_read`.** Where `mod_forcing_read`
  (NCAR) is bilinear-index + two-stage time-interp coefficients over bracketing slices, `read_other_NetCDF` is: read ONE
  2D slice, fill missing/land values ON THE RAW REGULAR GRID, then a direct per-node bilinear with NO time interpolation
  (CORE runoff is time-constant; the SSS climatology is read once per month — `update_monthly_flag=mstep==1` for a short
  January run, so a single record `i=month=1`). runoff: `read_other_NetCDF('Foxx_o_roff', rec 1, check_dummy=.false.)`
  → missing→0, then `/1000` (kg/s/m²→m/s); Ssurf: `read_other_NetCDF('SALT', rec 1, check_dummy=.true.)` → a 30-neighbour
  expanding-box average fills land. Both `do_onvert=.true.` (vertices); the centroid path is unexercised → guarded with
  `error stop` (don't ship ungated code).
- **Why it's byte-exact with NO new gate subtlety:** the raw-grid dummy fill (a deterministic serial loop over the global
  `ncdata(lon,lat)`) AND the per-node bilinear are PARTITION-INDEPENDENT (each owned node interpolates from the full global
  raw grid; no cross-rank reduction). So owned values are byte-identical at any partition — like `mod_forcing_read`, every
  rank reads the file directly (the oracle's read-on-0 + BCast moves the same bytes; no BCast needed). The on-disk fields
  are `float`; netCDF converts float→real64 on read EXACTLY as the oracle's `real64 ncdata`. **The `missing_value`
  attribute (`1.e30f` runoff / `-99.f` SALT) is ALSO float → reading it into real64 (`nc_get_att_dp`) makes the `==miss`
  equality match the (also float→real64) data values bit-for-bit** — read the attribute the same way the data is read, or
  the mask drifts. `geo_coord_nod2D/rad` is the same deg conversion already proven in M3f-3a.
- **THE TRAP (cost the only real debugging): a standalone 1-rank driver does NOT populate `partit%myDim_nod2D` — it stays
  0.** These hand-built drivers (`fesom_lifecycle`, `fesom_lifecycle_native`, …) call `par_init` + `read_mesh` but the
  ocean/ice kernels are invoked WITHOUT `partit` (the M2.12 optional-`partit`-ABSENT path loops `mesh%nod2D`); nobody ever
  reads `partit%myDim_nod2D`, so it's never set (the `synthesize_1rank` that would set it isn't on this path). I wrote
  `read_other_NetCDF` to derive the node count as `present(partit) ? partit%myDim_nod2D+eDim_nod2D : mesh%nod2D` (the
  M2.12 idiom) and PASSED `partit` → `num=0` → the interp loop `do n=1,num` no-ops → **`model_2Darray` came out all-zeros
  with NO crash and NO error** (the netCDF read, the fill, everything upstream was correct; only the final write was
  empty). Presented as "runoff AND Ssurf both EXACTLY 0" (diff = the full prescribed value). Bisected with 3 debug prints
  (read `[min,max]` OK → fill `[min,max]` OK → `num=0`, `haspartit=T`). **FIX: at 1-rank OMIT the optional `partit`** (so
  `num=mesh%nod2D`) — exactly how the same driver calls `ocean2ice`/`ice_timestep`. **Generalizable: an optional-`partit`
  routine whose loop bound comes from `partit%myDim_*` will SILENTLY produce zero-length output if a caller passes a
  partit whose dims are unset — a degenerate loop is not an error. Either omit partit at 1-rank, or assert `num>0`. M3f-4
  (multi-rank) WILL pass partit, where `read_dist_partition` sets the dims.**
- **Completion method — self-check, then flip.** M3f-3b first gated the 2 new fields by a per-step self-check INSIDE the
  prescribe-mode lifecycle (`max|d clim(runoff,Ssurf)|=0` while the step still USES the prescribed values, so the
  195-gate can't regress — the L41 pattern). Then the driver was made dual-mode: `FESOM3_ATMFLUX_FILE` present ⇒ prescribe
  + self-check; absent ⇒ FULLY NATIVE (a `compute_native_forcing` helper feeds either the self-check print or an
  `apply_native_forcing` that writes `atm%*`/`ice%stress_atmice`). Flipping to fully-native and re-running the 195/325-gate
  is the end-to-end proof — the per-step flux self-check (`max|d(hf,wf,vs,rs,ss)|=0`) confirms the native atmosphere drives
  identical air-sea fluxes. **No regression:** ctest 13/13 + forcing pi + step-65 1-rank + step-65 MR dist_2 + pressure 57
  + iceflux 5×2 all `max|Δ|=0`/green (the `mod_io_netcdf` additions don't touch existing readers; the kernels are untouched).

## L43 — Multi-rank fully-native lifecycle (M3f-4): the whole coupled sea-ice + air-sea + atmosphere step byte-matches at multi-rank FIRST TRY; the M2.12/M3 kernels were already MR-ready — the only NEW work was two loop-bound/reduction fixes

M3f-4 built the MULTI-RANK analog of the fully-native lifecycle (`fesom_lifecycle_native_mr` + `tools/run_lifecycle_
fullynative_gate_multirank.sh`) and byte-matched FESOM2 `max|Δ|=0` on the per-rank gid-keyed NODE substeps — **195
records (3-step) AND 325 records (5-step), BOTH whichEVP=0 (std EVP) AND whichEVP=1 (mEVP), on CORE2 dist_2 AND
dist_8**, ALL on the FIRST gate run. This closes M3 (tag `m3`): the native sea-ice EVP + advection + thermo + air-sea
budget AND the whole atmosphere (NCAR read + bulk + 2 stresses + runoff + Ssurf, 16/16 arrays) drive the multi-step
coupled CORE2 ocean step byte-exactly at multi-rank with NO prescribed input. The reusable lessons:

- **The hard part was an AUDIT, not new code.** The single most valuable step was verifying — kernel by kernel against the
  FESOM2 oracle — which halos each ice/coupling routine reads and which it already exchanges, BEFORE writing anything.
  The finding: every M3a–e kernel was ALREADY multi-rank-correct because each was transcribed with the M2.12 optional-
  `partit` pattern (owned/halo bounds + the FESOM2 exchanges) and gated for arithmetic at 1-rank. `ice_setup` sizes via
  `local_dims` (nNodL/nElemL); `ocean2ice` exchanges `u_w/v_w` (srfoce_temp/salt/ssh computed over owned+halo, no
  exchange); the EVP exchanges only `uice/vice` per subcycle; `ice_fct_solve` exchanges its solve intermediates +
  final `a/m/m_snow`; `cut_off` + `thermodynamics` LOOP over owned+halo (so their outputs — incl. `values_old`,
  `flx_h`, `flx_fw`, `t_skin` — stay halo-valid with no tail exchange); `oce_fluxes`'s `integrate_nod_2D` does the
  cross-rank `allreduce_sum`. Only TWO things were genuinely missing.
- **The EVP needs NO sigma (element) halo — the owned-node-completeness invariant carries it.** `stress_tensor` writes
  sigma11/12/22 at OWNED elements; `stress2rhs` reads sigma at OWNED elements (its element loop is `1..myDim_elem2D`)
  and scatters into nodes; the velocity update reads `u_rhs_ice` at OWNED nodes only. Because an owned node's complete
  element-neighbourhood is owned (the M2.12a area invariant), the owned-node rhs is complete from the owned-element
  scatter — so sigma is never read at a halo element and FESOM2 (`ice_EVP.F90`) has NO `exchange_elem(sigma)`, only
  `exchange_nod(U_ice,V_ice)`. Verified by grepping the oracle's exchange calls (std EVP: 1 velocity exchange; mEVP:
  the same on `u_ice_aux/v_ice_aux`). The prognostic sigma elastic memory persists owned-only across steps — correct.
  **Meta: before adding an element-halo exchange, check whether the consuming kernel only reads OWNED elements; the
  scatter-into-owned-nodes pattern needs no element halo (L33's "verify the invariant before building machinery").**
- **Missing #1 — `mesh%ocean_area` was a local-only sum (the M3e landmine, now closed at MR).** `compute_node_areas`
  built `ocean_area`/`ocean_areawithcav` as `vol`/`vol2` summed over OWNED nodes (correct at 1-rank where local==global,
  but at npes>1 that is only the rank's partial sum). `oce_fluxes` divides the flux-balance net by `ocean_area`, so a
  partial divisor would scale every balanced field wrong on every rank. Fix (`mod_mesh_areas.F90`): at npes>1
  `allreduce_sum(vol/vol2)` — FESOM2 `oce_mesh.F90:2389` `MPI_AllREDUCE(MPI_SUM, MPI_DOUBLE_PRECISION)` over the SAME
  owned partial sums, byte-identical (L6 deterministic tree). This was the `mod_mesh_areas.F90` TODO; L39 had already
  flagged `ocean_area` as un-gated until `oce_fluxes` consumed it — M3f-4 is where its MR form first matters.
- **Missing #2 — the bulk forcing looped over `mesh%nod2D` (GLOBAL), the classic MR loop-bound trap.** The three
  `mod_forcing_bulk` routines (`forcing_bulk_ncar`/`forcing_wind_stress`/`forcing_ice_stress`) looped `do i=1,mesh%nod2D`
  — at npes>1 that is the global count (wrong + OOB on the local arrays). FESOM2 computes the bulk + stresses over
  `myDim+eDim` (`gen_forcing_couple.F90:703/738`) with NO exchange (the per-node bulk is partition-independent given
  halo-valid inputs). Fix: optional `partit` + `owned_bounds` → loop `nNodL`; absent/npes==1 ⇒ `mesh%nod2D` (1-rank
  verbatim). The NCAR read (`mod_forcing_read`) needed only `frc%nnod = nNodL` (it already loops `frc%nnod`); and
  `read_other_NetCDF` already derived `myDim+eDim` from `present(partit)` (the L42 trap) — so at MR, PASS `partit`
  (`read_dist_partition` sets `myDim`), at 1-rank OMIT it. **Sweep `grep "do .*=.*mesh%(nod2D|elem2D)"` over every
  module on the new MR path — a loop bounded by the global mesh count is the first thing to break, and it fails
  silently (OOB read of garbage, or short loop) rather than loudly.**
- **The oracle's 1-rank-only dumps skip cleanly at npes/=1 — the gid-keyed `dump_shim` is the MR validation.** The
  `fesom_flux_dump` / `fesom_atmflux_dump` shims `print "skipped (npes/=1)"` and return, so there is no per-step flux
  self-check at multi-rank (that was the 1-rank isolation aid). The REAL gate is the built-in per-substep `dump_shim`
  (`<prefix>.<mype5>`, keyed by `myList_nod2D` gid) which IS multi-rank — both codes share `dist_<NP>/myList`, so each
  global probe id is owned by the same rank on both, and `dump_diff.py --glob` matches per-rank by gid (the L8/L34
  same-partition rule). The forced FESOM2 lifecycle runs cleanly at npes>1 (the 1-rank `io_meandata`/`next_io_rank`
  workarounds are inactive — they patched 1-rank-only bugs; the normal `io_gather` runs). The MR driver is fully-native
  ONLY (no prescribe path — the per-rank atm dump the prescribe path would need isn't produced at MR).
- **Building the MR driver = `fesom_lifecycle_mr` (MR mesh/state/dump scaffold) ⊕ `fesom_lifecycle_native` (ice/atm/
  forcing setup + runloop), at LOCAL sizes.** Take the `set_partition`+`read_mesh`(local remap)+`nNodO/nNodL/nElemF`
  bounds + gid-keyed `dump_init` + `do_ic3d(...,partit)`/`muscl_adv_init(...,partit)`/`init_stiff_mat_ale(...,partit)`
  from the ocean MR driver; take the `ice_setup`/`alloc_atm`/native-forcing-read + the `ocean2ice → apply_native_forcing
  → ice_timestep → oce_fluxes_mom → oce_fluxes → step_oce` runloop from the native driver; size EVERY array `nNodL`/
  `nElemF` and append `, partit` to every kernel call. No arithmetic, no new exchange logic in the driver — the kernels
  carry it. First try, `max|Δ|=0`. **No regression:** ctest 13/13 + 1-rank fully-native 195 + iceflux 5×2 + forcing pi
  + step-65 1-rank + step-65 MR dist_2 all `max|Δ|=0`/green.

## L44 — GM/Redi in the forced/native lifecycle + multi-rank (M4e + M4f): once each kernel is transcribed MR-ready, "turn the feature on in the next driver" is pure WIRING — the integration + multi-rank milestones cost a config block, not new physics

M4e (GM+Redi in the FORCED / fully-native 1-rank lifecycle) and M4f (the same at MULTI-RANK, CORE2 dist_2/dist_8) each
byte-matched FESOM2 `max|Δ|=0` on the FIRST gate run — 195 (3-step) + 325 (5-step) records, BOTH whichEVP, with NO new
physics code, NO halo extension, and (M4f) NO `dist_N` invariant check needed. M4 (GM bolus + Redi isopycnal diffusion)
is now complete end-to-end (1-rank + multi-rank, unforced + forced/fully-native, GM-only + GM+Redi, both EVP; tag `m4`).
The reusable lessons:

- **The integration milestone (M4e) is a config + allocation block, copied verbatim.** All the GM/Redi *physics* landed
  in M4a–M4d inside `step_oce`'s `Fer_GM`/`Redi`-guarded chain (producers → `init_Redi_GM`/`fer_solve_Gamma`/
  `fer_gamma2vel` → the Redi diff terms), reading the module flags + the `dyn%work` GM/Redi arrays. So "enable GM+Redi in
  the native forced lifecycle" = paste the `FESOM3_FER_GM`/`FESOM3_REDI` env reads + the work_core config + the
  `dyn%work` array `allocate`/init block from `fesom_lifecycle.F90` into `fesom_lifecycle_native.F90`; `step_oce` (called
  1-rank, no `partit`) does the rest. The oracle side is symmetric: the forced runner already runs the real GM+Redi
  `oce_timestep_ale` — it just needed a `FER_GM`/`REDI` env to KEEP work_core `Fer_GM`/`Redi=.true.` instead of the
  reduced-M2 sed-off (mirror the env already added to the unforced runner). NO oracle source change. The new test value:
  it is the FIRST GM/Redi exercise on a FORCED (non-zero surface-flux) ocean — the native heat/freshwater/salt budget
  perturbs the surface T/S → density → bolus/slope, and over multiple steps GM/Redi feeds back into the ice via the
  evolving surface state ocean2ice reads. That whole coupled loop byte-matched.
- **The multi-rank milestone (M4f) is the SAME wiring on the MR driver — because each kernel was transcribed MR-ready
  from the start.** Every M4 routine was written (M4a–M4d) with the M2.12/M3 optional-`partit` pattern (owned-loop bounds
  + the FESOM2 halo exchanges), and `mod_step_oce::step_oce` already threads its optional `partit` to all of them. The
  partit-present path therefore EXISTED but had never been exercised until M4f — and it just worked. The ONLY edit was the
  GM/Redi config + array-allocation block into `fesom_lifecycle_native_mr` at LOCAL sizes (`nNodL`/`nElemF`; `fer_uv`
  sized `nElemF` to match `dyn%uv`, since the bolus reads it over owned+halo `nElemL`). This is L43's lesson generalised:
  **if you pay the optional-`partit` tax at transcription time, the multi-rank gate is a driver-wiring task, not a
  porting task.** (M4f reused L43's `fesom_lifecycle_native_mr` scaffold unchanged — only the GM/Redi block was added.)
- **The `dist_N` invariant check (L33) was UNNECESSARY here, and recognising that saved a step.** The HANDOFF/plan
  insisted on running the cheap invariant check before assuming the GM/Redi halos suffice. But the GM/Redi producers use
  the EXACT mesh operations the M2.12 dynamics already proved invariant: `compute_sigma_xy` accumulates element gradients
  over `nod_in_elem2D(n)` for owned `n` (owned node ⇒ owned element-neighbourhood, M2.12a) exactly like
  `momentum_adv_scalar`; the Redi diff terms (`diff_part_hor_redi`/`diff_ver_part_redi_expl`/`diff_ver_part_impl_ale`)
  loop OWNED edges/nodes exactly like the M2.12c-3 off-path `diff_part_hor_redi`. Reusing a proven invariant ≠ needing to
  re-verify it. The gate IS the check: a pass on dist_2 AND dist_8 confirms it. (Run the explicit check when a kernel uses
  a NEW access pattern; skip it when it reuses one already gated MR.)
- **Two MR idioms the M4 producers lean on, both byte-exact.** (1) *Compute-at-halo instead of exchange:* `sw_alpha_beta`
  and `tracer_gradient_z` loop owned+halo (`nNodL`) and recompute the per-node function at halo nodes from halo-valid
  inputs (the tracer `values` are exchanged at the end of each tracer's solve; geometry is valid everywhere) — so the
  halo outputs are byte-identical to the owner's, with NO exchange. Valid precisely because the function is pure per-node
  scalar arithmetic. (2) *A rank-3 `exchange_nod` = the oracle's per-component `MPI_BARRIER` aux dance:* `compute_sigma_xy`
  exchanges `sigma_xy(2,nl-1,nNodL)` with one rank-3 `exchange_nod`; FESOM2 does a 3×`MPI_BARRIER` + per-component
  `exchange_nod(aux)`. Same owner→halo bytes moved, byte-identical result — the barriers are synchronization, not data.
- **No-regression for an integration/MR milestone = prove the OFF path is byte-neutral + the SHARED harness edits are
  byte-neutral.** Since M4e/M4f touched no core-physics source (only the two native drivers + the gate scripts), the
  kernel/step gates are provably unchanged (the binaries didn't relink). The two real risks: (a) the driver's new
  GM/Redi block must be skipped when `FESOM3_FER_GM`/`REDI` are unset (it is — `use_fer_gm`/`use_redi` false ⇒ block not
  entered ⇒ arrays unallocated ⇒ `step_oce`'s `Fer_GM`/`Redi`-guards short-circuit); (b) the SHARED oracle runner
  (`run_lifecycle_forced_core2.sh`, used by M3f's gates too) must still sed `Fer_GM`/`Redi` OFF by default. Confirmed by
  re-running the GM/Redi-OFF fully-native gate (1-rank AND dist_2) `max|Δ|=0` + the M4d unforced GM+Redi gate + ctest 13/13.

## L45 — Wiring KPP into the step (M5b): an isolated-kernel gate that calls sub-kernels DIRECTLY skips the driver's own glue steps — the live driver exposed the missing `smooth_blmc` as a ~1e-3 Kv drift

**Context.** M5a byte-proved the whole KPP module (init/wscale/ri_iwmix/bldepth/blmix_kpp/enhance) `max|Δ|=0` on pi AND
CORE2 — 95 pressure-gate fields, all four sub-steps first try. M5b assembled the real `oce_mixing_kpp_driver` (prestep →
ri_iwmix → bldepth → blmix → enhance → combine → node→elem viscAE average) + the `mix_scheme_nmb==1` dispatch in
`step_oce`, and gated the UNFORCED CORE2 lifecycle. The KPP-alone gate diverged at step 1, substep 4 (MIXING) — `Kv`
`|Δ|=9.6e-4`, "close but not byte-exact," ~12/65 records, localized to specific nodes.

**The hunt (a textbook bisection).** The PP path was byte-exact and every M5a kernel was byte-proven *individually*, so
deterministic code + identical inputs ⇒ identical output. A divergence therefore means an INPUT differs. Bisecting:
- Compared `Kv` (post-`mo_convect`, the dumped field) vs `Kv_double` (pre-`mo_convect`): `Kv_double` MATCHED at the probe
  node, `Kv` didn't → suspected `mo_convect`. But `mo_convect`'s namelist + inputs matched, so that was a red herring.
- Printed the FULL column: the divergence was a single level — node 1001, **nz=3 = kbl-1**: oracle `diffKt(3)=9.726e-4`,
  FESOM3 `1e-5`. The `kbl-1` level is exactly what `enhance` writes, so I re-read `enhance` (byte-identical) and worked
  the algebra BACKWARD: with the dumped `dkm1=0`, `caseA=1`, `diffK=1e-5`, the oracle's 9.7e-4 needs `δ≈4.59`, but the
  oracle's OWN `hbl=10.498`/`zbar` give `δ≈0.05` → 9.5e-6. **The oracle's final `blmc` was inconsistent with its own
  `enhance` inputs** ⇒ a step exists BETWEEN enhance and the final dump that I wasn't replicating.
- That step is `smooth_blmc` (oracle `oce_ale_mixing_kpp.F90:439-449`, default `.true.`): 3 area-weighted `smooth_nod`
  sweeps × the 3 `blmc` channels, run AFTER enhance, BEFORE the combine. A node's `blmc` becomes a neighbour-smoothed
  value — which is why node 1001's local 9.5e-6 became 9.7e-4 after smoothing.

**Root cause + fix.** The M5a-4 isolated gate called `blmix_kpp`/`enhance`/combine DIRECTLY (the pressure shim, both
sides) — it never ran the driver's `smooth_blmc`, so the smoothing was invisible until the live driver. Fix: add the
`smooth_blmc` block to the driver, reusing the **M2.1 byte-proven** `smooth_nod` (the bvfreq smoother, made `public` in
`oce_pressure_bv`; it IS FESOM2 `gen_support.F90::smooth_nod3D`). Result: KPP-alone 195+325 + KPP+GM+Redi 195 all
`max|Δ|=0`, no regression.

**Generalisable.**
- **An isolated-kernel gate proves the KERNELS, not the DRIVER.** When the gate harness calls sub-kernels directly
  (`call blmix; call enhance; <inline combine>`), it skips any glue the production driver runs between/around them
  (`smooth_blmc`, exchanges, in-place-vs-copy, floors). RE-READ the real driver top-to-bottom when assembling it and
  list EVERY step between the proven calls — don't assume "assembly of proven pieces" is complete.
- **A new feature's gate may exercise a code path its component gates never hit.** M5a-3/a-4 gated bldepth/blmix with
  `use_sw_pene=.true.` + nonzero fluxes; the UNFORCED M5b run has `ustar=0`/`bfsfc=0`/`use_sw_pene=.false.` — a different
  branch set. (Here that wasn't the bug, but it's why "all kernels green" didn't imply "driver green.")
- **Bisect a "close but not exact" drift to a SINGLE (node, level), then work the local formula BACKWARD.** Solving
  `enhance` for the δ the oracle's output implies, and finding it contradicted the oracle's own geometry, is what proved
  the difference was UPSTREAM of enhance (a missing step) rather than a transcription error IN enhance — far faster than
  re-diffing every kernel.
- **The oracle's built-in instrumentation is the fastest bisection tool.** `FESOM_KPP_DUMP_DIR`/`FESOM_KPP_DUMP_STEP=1`
  dumped every KPP intermediate (dVsq/dbsfc/ri_*/bldepth/blmc_*/dkm1/final) gid-keyed es24.16 — comparing per-node-per-
  level localised the divergence to one cell in minutes vs instrumenting the F3 side blind.
- **House-keeping traps cleared in passing:** the driver subroutine `oce_mixing_KPP` clashed (case-insensitive) with the
  module name `oce_mixing_kpp` (ifort #6645) → renamed `oce_mixing_kpp_driver`. And adding a `mix_scheme_nmb==1` dispatch
  means every reduced-M2 driver must now PIN `mix_scheme_nmb=2` (the o_PARAM default is 1=KPP) or it silently reroutes the
  PP gates to KPP. Passing an UNALLOCATED allocatable to an optional dummy = "absent" (F2008) is the clean way to make one
  shared `step_oce` call serve both PP (no `stress_node_surf`) and KPP.

## L46 — Shortwave penetration + the DEAD ghats nonlocal flux (M5c): of the two tracer-TDMA terms the plan named, only `sw_3d` is live in work_core — `ghats` is gated off by `use_kpp_nonlclflx=.false.`; and `use_sw_pene` already existed in the right module

**Context.** M5c's RESUME note (and the M5 plan) listed TWO new terms for `diff_ver_part_impl_ale`: the `sw_3d`
shortwave heating AND the `ghats` KPP nonlocal counter-gradient flux. Reading the ACTUAL oracle
(`oce_ale_tracer.F90:892-996`) showed the ghats block is guarded by `if (use_kpp_nonlclflx)` **FIRST**, then
`mix_scheme_nmb==1` second. `use_kpp_nonlclflx` defaults `.false.` (`oce_modules.F90:168`, `&tracer_phys`), is never
assigned `.true.` anywhere in the oracle source, and is **ABSENT from every work_core namelist** → the ghats term is
DEAD in the production CORE2 config (exactly like `MLD1_ind`/`K_hor` in M4). The plan's "mix_scheme==1" was the
secondary test; the real gate is the namelist flag. So the production M5c work is `sw_3d` alone.

**Result.** Production gate `run_lifecycle_kpp_native_gate_core2.sh` (KPP+sw_pene+GM+Redi, fully-native forced
lifecycle) `max|Δ|=0` at 195 (3-step) AND 325 (5-step), BOTH whichEVP — FIRST try, no NOVECTOR. `cal_shortwave_rad`
(M2.10c, byte-proven) wired into the runloop right after `oce_fluxes` (it IS the last line of FESOM2 `oce_fluxes`,
`ice_oce_coupling.F90:873`); the oracle flux dump (`fesom_module.F90:746`) records `heat_flux` AFTER it, so the
per-step flux self-check reads `0` for the visible-band-modified `heat_flux` too. The ghats term was still transcribed
faithfully and gated by a dedicated `KPP_NONLCL=1` variant (`use_kpp_nonlclflx=.true.` injected into the oracle
`&tracer_phys` + the FESOM3 `FESOM3_KPP_NONLCL` env): `max|Δ|=0` BOTH whichEVP, and NON-VACUOUS (ghats changes T/S
by up to 7.3e-5 over 3 steps vs ghats-off) — so the otherwise-dead term is byte-verified, not an L39 landmine.

**Generalisable.**
- **"Read the actual code, not the plan summary" — the FLAG edition.** Before implementing a term the plan names,
  grep its gating flag's DEFAULT *and* its presence in the live `work_core` namelist. A term guarded by a flag that
  defaults off and is unset is DEAD — implementing it as "active production work" wastes effort and (worse) ships
  ungated code. Confirm liveness with `grep -rn "<flag> *= *\.true\." src/` (does anything turn it on?) + `grep -rn
  "<flag>" work_core/` (does the config set it?).
- **Before adding a config flag, grep for it — it may already exist in the correct module.** `use_sw_pene` already
  lived in `mod_config` (the g_config analog — the faithful home, since FESOM2 `use_sw_pene` is `g_config`) and was
  ALREADY read by the KPP `bldepth` + `fesom_lifecycle` + `fesom_pressuredump`. Adding a SECOND `use_sw_pene` to
  `mod_param_phys` made two flags; the production gate passed only because both happened to be `.true.`. A duplicate
  that "works by coincidence" is a latent split-brain bug — one source of truth, in the module that matches the FESOM2
  origin group.
- **An in-code default of `.true.` is a footgun when the consumer dereferences a conditionally-allocated array.**
  `dyn%work%sw_3d` is allocated ONLY when sw_pene is enabled, but the byte-gate drivers configure physics in-code (they
  do not read `namelist.config`), so a `.true.` default would fire the `sw_3d` term against unallocated memory in every
  PP/reduced driver. Flip such gate-driver defaults to `.false.` ("unset = off") and let the few feature drivers opt in
  explicitly — fail-safe beats faithful-to-the-namelist-default here.
- **Gate the dead code too.** When a faithful port includes a term that's off in production, add a cheap variant that
  turns it on (both oracle + F3) and prove it byte-exact + non-vacuous, rather than shipping untested transcription.

**M5d (multi-rank) — pure wiring, first try (the L44/M4f lesson again).** The whole production KPP+sw_pene+GM+Redi
native lifecycle byte-matched at MULTI-RANK (CORE2 dist_2 + dist_8, BOTH whichEVP, 195+325) on the FIRST try — because
the KPP module was written optional-`partit` from the start (`oce_mixing_kpp_driver` already does owned-loops +
`exchange_nod(blmc/Kv_double/ghats/viscA)` + the 3-sweep `smooth_blmc`), so M5d was a config block in
`fesom_lifecycle_native_mr` (LOCAL-sized `nNodL`) + one signature fix. The ONE new MR gotcha:
- **In the LOCAL mesh, `mesh%nod2D` = the GLOBAL node count** (`read_mesh_local: mesh%nod2D = nNodG`), NOT the local
  owned/owned+halo count. A 1-rank-correct routine that loops `do n=1,mesh%nod2D` over LOCAL-sized arrays is therefore
  OUT OF BOUNDS at multi-rank. `cal_shortwave_rad` had exactly this (it predates the optional-`partit` discipline, M2.10c)
  — fixed by adding an optional `partit` and looping `nNodL` from `owned_bounds` (absent ⇒ `nNodL=mesh%nod2D`, the 1-rank
  /pi path verbatim). When lifting an OLDER routine to MR, grep it for bare `mesh%nod2D`/`mesh%elem2D` loop bounds — those
  are the global-count landmines; the routines written since M2.12 already use `owned_bounds`/`local_dims`.

## L47 — The Shchepetkin zlevel PGF (M6a-1) byte-matches first try; a latent m5 shim crash (unallocated `sw_3d`) surfaced only when the heap re-laid-out — the M5a-3 META-LESSON realized

**M6a-1 (the full-free-surface Shchepetkin & McWilliams density-Jacobian PGF, `pressure_force_4_zxxxx_shchepetkin`)
transcribed byte-exact first try** (CORE2 pressure gate, `pgf_x_shchep`/`pgf_y_shchep` `max|Δ|=0`, non-vacuous ~3e-5).
Two reusable points:
- **A new algorithm that REDUCES to the proven one in the gate's regime still needs a real byte-gate, but know what it
  does/doesn't test.** At rest (η=0) with flat full cells, `Z_3d_n(nz,:)` is horizontally constant, so the Shchepetkin
  surface-slope correction `dz_dx = Σ gradient_sca·Z_3d_n ≈ 0` (sum of a constant times shape-fn derivatives = 0 to ULP)
  and the scheme collapses to the hydrostatic `drho_dx·helem·g/ρ0` integral = the linfs PGF (the two agreed to ~1e-17).
  The isolated gate is byte-exact vs the oracle (validates EVERY line incl. the Newton-interp `drho_dz` + the length-3
  array divides — the L29 watch), but because `drho_dz`'s contribution is `·dz_dx ≈ 0` at rest, a tiny `drho_dz`
  divergence would be MASKED here. The η≠0 regime (lifecycle, where the surface tilts ⇒ `dz_dx ≠ 0`) is the strong test
  — defer the divide-stress verification there, and SAY SO in the handoff rather than claiming the rest gate covers it.
- **A Release build can run for a whole milestone on an unallocated-allocatable read, then segfault later when the heap
  re-lays-out — exactly the M5a-3 META-LESSON, now with teeth.** The CORE2/pi pressure shim set `use_sw_pene=.true.` (M5a-3,
  for KPP `bldepth`) but never restored it before the M2.9a tracer solve, so the M5c `sw_3d` term in
  `diff_ver_part_impl_ale` (`if(use_sw_pene .and. id==1)`) dereferenced UNALLOCATED `dyn%work%sw_3d`. At m5-commit the
  Release descriptor happened to point at readable zeros (term added 0 ≡ oracle's term-off) so the gate was green; a later
  run/rebuild shifted the heap so the descriptor pointed at unmapped memory → SIGSEGV (Release backtrace named only
  `diff_ver_part_impl_ale`; the DEBUG `-check all` pinned it: "Attempt to fetch from allocatable variable SW_3D when it is
  not allocated", `oce_ale_tracer.F90:665`). FIX: mirror the oracle (`fesom_pressure_dump.F90:653`) — restore
  `use_sw_pene=.false.` right before the tracer-solve loop (the `bldepth` above already ran and uses a LOCAL `sw3d_kpp`
  ARG, not the module flag). Generalisable: (i) when a shim flips a global flag for one kernel, restore it before the next
  kernel that reads it — mirror the oracle's flag lifecycle exactly, don't leave it set; (ii) if a "green" gate suddenly
  segfaults with NO source change, suspect an uninitialised/unallocated read whose luck ran out on a heap re-layout —
  reach for the DEBUG `-check all` build first (it converts the silent UB into a named, located error).

## L48 — The FIRST SLURM-batch gate (M8b `dist_864`) exposed TWO latent `set -e` footguns in the shared harness — both invisible to the login-node mpirun gates, both about a non-fatal command made fatal by `set -euo pipefail`

**M8b's port (`forcing_sbc_do`) was byte-exact at `dist_2`/`dist_8` on the FIRST clean run — but two HARNESS bugs,
dormant since fesom3 had only ever gated on the login node with `mpirun`, surfaced the moment the run got LONGER
(48 steps) and moved to a SLURM compute node (`srun`).** Both were `set -e` + a command that returns non-zero in a
context the login gates never hit. The byte-gate diff (`max|Δ|=0`) was never in question; the harness just aborted
before the comparison ran (oracle "completes" but F3 never launches → no `run_f3.log`).

- **`grep … | head` SIGPIPE (long runs).** The oracle runner's LAST line was a diagnostic
  `grep -E '…' run.log | head` (`run_lifecycle_forced_core2.sh`). With `set -o pipefail`, once the run is long
  enough that `grep` produces enough matching lines (~50 at 48 steps; ~19 at 16 steps did NOT trip it) that `head`
  closes the pipe after 10 BEFORE `grep` finishes writing, `grep` dies with SIGPIPE → exit 141 → pipefail makes the
  pipeline 141 → `set -e` aborts the oracle runner → the gate's `bash oracle | tail -4` propagates 141 → the gate
  aborts after `[1/3]`, before `[2/3]`. The legacy gates ran 3-5 steps (<10 matches) so it NEVER fired. **Symptom:
  `GATE_EXIT=141`.** FIX: `… | head || true` (a diagnostic's exit must never gate the run).
- **`ulimit -s unlimited` RAISE on a hard-capped compute node.** `env/levante.dkrz.de/shell.intel` (sourced by
  `env.sh`) did a bare `ulimit -s unlimited`. On the LOGIN node the hard stack limit is unlimited so this succeeds;
  on a SLURM COMPUTE node the default is `soft=8192 hard=unlimited`, BUT my batch script first ran `ulimit -s 204800`
  which sets BOTH soft AND hard to 204800 — after which `ulimit -s unlimited` tries to RAISE the hard limit and gets
  "Operation not permitted" (rc=1). Under the oracle's `set -e`, `source env.sh` (which sources shell.intel) returns
  non-zero → the oracle aborts at `source env.sh`, BEFORE `mkdir` → **missing run dir, ~1 s, no output at all.** The
  debugging trap: env.sh sources shell.intel with `>/dev/null 2>&1`, so the failing line is invisible; and a DIRECT
  oracle run (no pre-`ulimit -s 204800`) WORKED, masking it. FIXES (both applied): (a) `shell.intel`/`shell.gnu`
  `ulimit -s unlimited 2>/dev/null || true` (best-effort — the env must survive a caller's `set -e` when the hard
  limit is already capped); (b) do NOT pre-cap the hard limit in the batch script — let the gate runners set the
  204800 dump-writer stack per leg AFTER sourcing env.sh (lowering soft+hard is always allowed; raising is not).
- **META (reusable for M8c/d/e and any future batch gate):** a Fortran↔Fortran byte-port is "done" when the DIFF is
  zero, but a gate is a SHELL PIPELINE wrapped in `set -euo pipefail`, and every `cmd | head`, `[ ] && cmd`, bare
  `ulimit`/`module`/`hostname -f` in it (or in `env.sh`) is a latent non-zero that a longer run or a different node
  can make fatal. When a batch gate dies in seconds with no diff output, the bug is almost always the HARNESS, not
  the port — re-run the failing leg with `bash -x` (un-redirected) on the actual compute node to find the aborting
  line, and prefer `|| true` on every non-essential diagnostic. The batch harness itself
  (`run_lifecycle_long_gate_dist864.sbatch` + the `srun`/`SLURM_JOB_ID`-gated branches in the runners) is NEW
  infrastructure for M8 — reuse it, it now works at 864 ranks.

## L49 — The FIRST JRA55-do byte-gate (M8c Step 1): a netCDF calendar attribute's trailing NUL made `trim(cal)=='gregorian'` silently FALSE — and the CORE2 substitute had hidden the trap for the whole M2.10→M8b port

**Every fesom3 forcing gate from M2.10 to M8b used the CORE2 NCAR substitute (`noleap`); the FIRST gate on the REAL
production atmosphere — JRA55-do-v1.4.0 (`gregorian` + `include_fleapyear=.true.`, 3-hourly, `tmid=0`) — diverged hard
(`ssh_rhs` ~1.0e4 at step 1).** The port logic was correct; the bug was a single trailing NUL byte on a string read
from netCDF.

- **Bisection (reusable recipe).** A per-field `atm_selfcheck` (reads the F2 `atmflux` STREAM dump, prints `max|Δ|`
  of each native atmosphere/bulk field over owned nodes) localized the divergence to the atmosphere READ, not the
  ice/flux coupling. Then: range diagnostics → argmax node → raw on-disk slice → time-axis scale → calendar
  attribute. The slice index was wrong: `rdate` on the gregorian Julian-Day scale (~2436205) vs `nc_time` built on
  the NOLEAP scale (~714684) ⇒ `binarysearch` capped to the wrong record ⇒ garbage wind/heat.
- **Root cause — the NUL byte.** `nc_get_att_text` returned `calendar = 'gregorian'//char(0)` (len 10): netCDF
  C-string attributes carry a trailing NUL. Fortran `trim()` strips trailing SPACES, NOT NUL — so `forcing_julday`'s
  `trim(calendar)=='gregorian'` was `.false.` ⇒ the routine fell to the `else` (noleap) branch ⇒ `julday=365*yyyy`
  instead of the gregorian Numerical-Recipes formula ⇒ wrong `nc_time` ⇒ wrong slice ⇒ `ssh_rhs` blew up.
- **Why it hid for 8 milestones.** CORE2's calendar string is `'NOLEAP'`, which `lc()`-lowercases to `'noleap'` and
  takes the SAME `else` branch *whether or not* the NUL is stripped — so the substitute dataset NEVER exercised the
  `=='gregorian'` comparison. The bug was invisible until the real JRA55 data finally drove the other branch. (The
  concrete cost of the "we never used CORE2" reality the user flagged: a substitute that drives a DIFFERENT branch of
  a shared routine masks any bug on the branch it skips.)
- **FIX (one place, byte-NEUTRAL for CORE2).** In `nc_get_att_text`, after the get, blank every control byte:
  `do k=1,len(txt); if (iachar(txt(k:k)) < 32) txt(k:k)=' '; end do`. The CORE2 no-regression gate (48 steps,
  `dist_2`, 3120 records) stayed `max|Δ|=0`; JRA55 then byte-matched `dist_2` (both whichEVP) + `dist_8`, all
  `max|Δ|=0`.
- **META.** Any string read from a C library (netCDF attribute, `units`, a variable `name`) can carry a trailing
  `char(0)` that `trim()` will not remove, so every `trim(x)=='literal'` is a latent silent-false. When a string
  compare "should" match but doesn't, dump `len_trim(x)` and the `iachar` of each byte before trusting it. And widen
  the lens: when a port is gated only with a SUBSTITUTE input, the gate proves only the branches that input drives —
  the first run on the real production data is itself a milestone, not a formality.

## L50 — Multi-year forcing rollover (M8c Step 2): year-file + month-climatology read-ahead byte-match across the SAME read routines as the cold start; cross every boundary CHEAPLY by starting at 23:00 of the last day, and an attribute-less file needs a no-fill sentinel

**With M8b's record/day crossing already byte-exact, M8c Step 2 closed the remaining JRA55 production rollovers — the
atmosphere YEAR-file (1958→1959), the SSS monthly climatology, and the Sweeney chlorophyll monthly climatology — all
`max|Δ|=0` at `dist_2` (Jan→Feb month gate, Dec31→Jan1 year gate, and a chl-Sweeney month gate).** The port is small
because every rollover re-drives a read routine the cold start already proved; the work is the TRIGGER and re-entrancy.

- **Year-file branch (`forcing_sbc_do`).** Port `gen_surface_forcing.F90:1510-1519` verbatim: `if (yearnew/=yearold)`
  ⇒ for each field `call forcing_read_grid(frc, fld, yearnew)` (re-opens `<var>.<yearnew>.nc`, re-anchors `nc_time`
  to the new year) + `force_newcoeff=.true.`, then OR `force_newcoeff` into the per-field crossing test so
  `getcoeffld` re-fires that step. `yearnew`/`yearold` come straight from `mod_clock` (the clock sets `yearold=yearnew`
  at the top of `clock`, then bumps `yearnew` on the day-365→1 wrap) — no separate state to track. The bilinear source
  indices (`forcing_build_bilin`) are NOT rebuilt: the grid is identical year-to-year (FESOM2 re-reads grid via
  `nc_readTimeGrid` but never re-runs `nc_sbc_ini`). `getcoeffld`'s fresh `binarysearch(rdate)` resets the stale
  prev-year `t_indx` to the new year's early-January bracket.
- **`forcing_read_grid` re-entrancy (prereq).** The cold-start grid read had a bare `allocate(nc_lon/nc_lat/nc_time)`;
  a second call aborts (ifort: already allocated). Guard each with `if (allocated) deallocate` first. NOTE `ntime`
  CAN change across a leap boundary (2920 non-leap vs 2928 leap) so the arrays must be freed, not reused — 1958→1959
  is 2920→2920 (no leap), but the headline must not assume constant `ntime`.
- **Monthly read-ahead (`roll_monthly_clim`, driver).** Port `gen_surface_forcing.F90:1590-1618`:
  `update_monthly_flag = (day_in_month==num_day_in_month(fleapyear,month) .and. timenew==86400._WP) .or. mstep==1`;
  on a hit read the NEXT month `i=month; if(mstep>1) i=i+1; if(i>12) i=1` via the SAME `read_other_NetCDF('SALT'/'chl',
  i)` the cold start used (so the slice bytes match). The trigger is the LAST instant of the month (read-ahead, NOT the
  month-change step) — `timenew==86400._WP` is EXACT because `dt=1800` divides 86400 and any clock-file start time is a
  multiple of `dt`. Keyed off the driver step counter `n` as `mstep`. The cold-start SSS/chl read at `i=month` stays
  (feeds the diagnostic + is byte-identical to the redundant `mstep==1` re-read).
- **Attribute-less file needs a no-fill sentinel (chl).** `Sweeney_2005.nc`'s `chl` variable has NO `missing_value`
  attribute; FESOM2's `read_other_NetCDF` captures `status` and IGNORES it, leaving `miss` undefined — harmless there
  only because the chl field is fully extrapolated to `[.001,.6]` with no `-99`/sentinels, so the fill loop never
  fires. FESOM3's `nc_get_att_dp` error-STOPS on an absent attribute, so add an optional `stat` and set `miss=-99`
  when absent: the hardcoded `==-99` check stays the only trigger ⇒ no fill ⇒ byte-identical to FESOM2's effective
  no-fill (the chl gate confirmed `max|Δ|=0`). (PHC2 SSS HAS `missing_value=-99` ⇒ `stat==0`, `miss=-99` anyway —
  the change is a no-op for the path that was already gated.)
- **Cheap boundary-crossing gates.** To cross a month/year boundary in 2 steps instead of a full day, start at 23:00
  of the last day: `START_CLOCK="82800 31 1958"` crosses Jan→Feb at step 3 with the SSS read-ahead firing at step 2;
  `START_CLOCK="82800 365 1958"` crosses 1958→1959 (the year-file rollover fires at step 3, the Dec→Jan SSS read-ahead
  at step 2). Each gate is non-vacuous (the `YEAR ROLLOVER`/`roll_monthly_clim` stdout markers prove the path fired)
  AND `max|Δ|=0` — a passing gate that never crossed the boundary is the failure mode to guard against.
- **META.** A multi-year run adds NO new physics — only triggers and re-entrancy around already-proven reads. Gate
  each new trigger at the EARLIEST step it fires by SEEDING the clock next to the boundary (don't run a real year),
  and print a stdout marker per fired path so a green gate can't secretly be a no-op. The production config is whatever
  the REAL work-dir namelist says (`work_core`: SSS `'CORE2'` monthly, chl `'Sweeney'` monthly) — NOT what an earlier
  plan assumed (the M8 scope doc still said chl `'None'`/constant from the CORE2-substitute era; the work-dir
  namelist, not the doc, is ground truth).

## L51 — The 2-year headline (M8c Step 3) diverged at day 107 because the FESOM3 PROCESS ran with flush-to-zero (FTZ) OFF while FESOM2 had it ON — root cause FOUND and FIXED, byte-exact

The full-production 2-yr gate (JRA55-do 1958, zstar+TKE+GM+Redi+sw_pene, dist_864) diverged first at
**step 5095 (day ~107), SSH_RHS, at all 126858 nodes at once** — the fingerprint of a GLOBAL scalar
(the freshwater `net` from `integrate_nod_2D`) spilling into `water_flux` everywhere. Bisected to the
origin with a windowed all-owned-node `FESOM_DUMP_ALL` (`tools/onset_allnode.py`) + per-gid extract
(`tools/extract_gid.py`).

- **The symptom (branch flip).** At node gid 119505, `m_snow` = **0.0 exactly in FESOM2** but a
  **denormal 1.34e-309 in FESOM3** (a snow thickness of ~1e-309 m = physically zero). The ice `budget`
  picks albedo via `if (hsn.gt.0.0_WP)` (snow ~0.85 vs ice ~0.65); the denormal `hsn>0` flips F3 to the
  snow branch → absorbed-SW shift → `t_skin` diverges **2.6%** (−0.4350 vs −0.4614) →
  evap/sublimation/`thdgr` → freshwater `flux` → global `net` → `water_flux` at every node → SSH spill.
- **ROOT CAUSE — a runtime FP-environment mismatch, NOT a transcription bug.** Every snow/ice routine is
  byte-identical F3↔F2 and the *compile* flags match (`-fp-model precise -no-prec-div -fimf-use-svml
  -init=zero -no-prec-sqrt -ip …`, both with `-fpe0`). But the **MXCSR flush-to-zero (FTZ) bit differed
  between the two PROCESSES at runtime**: FESOM2 ran FTZ **ON** (underflow → exactly 0.0), FESOM3 FTZ
  **OFF** (underflow → denormal kept). So where FESOM2 flushed the underflowing snow residual to 0.0,
  FESOM3 retained 1.34e-309 — and the `>0` test did the rest. Probed directly with
  `ieee_get_underflow_mode`: F2 `gradual=.false.` (flush-to-zero), F3 `gradual=.true.` (denormals kept).
- **Why `-fpe0` did NOT guarantee it.** `-fpe0` traps on FP *exceptions*; it does not pin the FTZ/DAZ
  MXCSR bits for the life of the process — those can start cleared and/or be reset by a linked library's
  init. (This SUPERSEDES the older note near the top of this file that "`-fpe0` implies FTZ" — it does
  not by itself; the runtime mode must be asserted explicitly.)
- **THE FIX (FESOM3 only — no FESOM2 change).** Force FTZ on with the standard intrinsic in
  `src/drivers/fesom_lifecycle_native_mr.F90`: `use, intrinsic :: ieee_arithmetic` →
  `if (ieee_support_underflow_control(1.0_WP)) call ieee_set_underflow_mode(gradual=.false.)` once at
  startup (after `par_init`) AND re-asserted each step (to survive any library MXCSR reset). Physically
  correct — a ~1e-309 m snow thickness IS zero. No tolerance, no reference change.
- **Confirmed byte-exact.** With the fix: FTZ gate (5100 steps, past day-107) `max|Δ|=0` over 331500
  records; and the 2-yr gate overlap `max|Δ|=0` over **1,138,800 records spanning steps 1→~17280 — a full
  model year (day 1→~360)**, 3.4× past the day-107 flip point, through the whole seasonal cycle. The
  year-rollover is separately byte-gated (M8c Step 2). (Full 2-yr clean re-run pending only because the
  heavyweight oracle leg crashed on ITS OWN output at the year boundary — not an F3 divergence.)
- **Diagnostic leftovers worth keeping.** (1) The first *byte* divergence in the window was a denormal
  red herring (`m_snow`=2.75e-308 vs 3.66e-308 at node 42847 — underflow noise absorbed by any normal-mag
  sum); `onset_allnode` needs a denormal floor (`max(|a|,|b|) < 1e-100 ⇒ ignore`) or it reports the wrong
  origin. The damage was never the denormal *arithmetic* — it was the discrete `>0` branch downstream.
  (2) `m_ice`/`a_ice` stayed `|Δ|=0` throughout — only the near-zero `m_snow` ever reached the denormal regime.
- **META.** Byte-identical source + identical *compile* flags do NOT guarantee bit-identity: the
  **runtime** FP environment (MXCSR FTZ/DAZ) is per-process and can differ, and it only bites once a field
  underflows into the denormal range and meets a discrete `if (x>0)` test. When a long gate diverges where
  short gates passed, suspect (a) a near-zero field crossing a branch — dump the *integrand/inputs upstream
  of the global reduction* (all 5 probes lighting at once = a global scalar) — AND (b) a runtime FP-mode
  mismatch: probe `ieee_get_underflow_mode`/MXCSR on BOTH processes, not just the source.

## L52 — Chasing the "F3 is ~1.5× slower" perf gap: I guessed twice and was wrong twice; the only real gap is FORCING (2.2×), ocean+ice are at parity/faster, and the blunders were a units misread + validating at the wrong scale

The story is fully captured in `docs/plans/2026-06-28-forcing-perf-investigation.md` (fresh-session handoff).
The transferable lessons:
- **MEASURE, DON'T GUESS — and convert units before comparing two codes.** I "diagnosed" the forcing read
  pattern and shipped a byte-exact rank-0-read+`MPI_Bcast` "fix" before measuring it. The A/B then showed it
  **+60% WORSE at dist_128** and **neutral at dist_512** → reverted (`83615d7`→`d8ad992`). Then I compared
  FESOM2's `runtime ocean: 18.14` against FESOM3's `28.66` and "found" a broad 2× ocean slowdown — but F2's
  block is **seconds for the whole run**, F3's timer is **ms/step**: `18.14 s / 600 = 30.24 ms/step`, so F3
  ocean is actually *faster*. ALWAYS `/nsteps*1000` the F2 numbers first. With correct units the picture is
  simple and matches the user's memory: ocean/ice parity-or-faster, **forcing alone is 2.24× (+8.7 ms/step)**.
- **Validate perf at the SCALE the gap lives.** dist_128 on one node masks I/O via the OS page cache (all
  ranks read from cache after the first); the forcing read cost is a multi-node phenomenon. A single-node A/B
  sent me the wrong way. Trust dist_512; treat dist_128 as a smoke test only.
- **The page cache beats a "redundant-read elimination" broadcast on-node**, and the ~800 KB slice broadcast
  is slow two-copy vader anyway (KNEM single-copy disabled in env.sh for the corruption bug). Moving *where*
  the read happens was neutral at scale → the forcing cost is NOT the per-rank read transfer; the real
  suspects are F2's persistent file handle + double-buffer slice cache (read 1 slice/crossing, F3 reads 2 +
  reopens every crossing) — but that too must be measured (sub-time open/read/interp) before fixing.
- **Byte-gating needs no compute partition** (login-node `run_lifecycle_jra55_gate_multirank.sh`); only the
  *perf measurement* needs scale. So correctness work can proceed during a Levante maintenance window; the
  speed confirmation waits for `compute`.
- **The permanent `mod_timer` breakdown (L-infra) is what made this legible** — it localizes the gap to one
  component and is directly comparable to F2's `rtime_*` (after the unit conversion).

## L53 — Bottom at vertices: the min∘max asymmetry that decided the bathymetry, the "pure alias" that was a wetness guard, and the regression net that was a step counter

**Context.** Moving the bottom from elements to vertices (`nlevels(e) = min` over the element's
vertex columns, `elvls.out` no longer read). Four things bit, and the first three are the kind that
look settled until you measure them.

**1. `min ∘ max` and `min ∘ min` are not two flavours of the same choice.** The element reduction is
fixed at `min` by the physics — velocity must vanish at a land corner — so it is not negotiable.
That makes `blayer = max over adjacent elvls` compose as `min ∘ max`, a morphological **closing**
where the dilation and erosion cancel; and `blayer = min over adjacent elvls` compose as `min ∘ min`,
a **double erosion** with nothing to undo it. Measured on core2: `+0.38 %` volume, 0 stagnant cells
versus `-9.27 %` volume, 3966 nodes whose deepest cell has no wet adjacent element and can only
exchange vertically. The "conservative" option was the destructive one. **Always compose the two
reductions on paper before picking the vertex rule, and measure both on the real mesh.**

**2. "It becomes a pure alias" is a claim about values, so measure it.** `nlevels_nod2D_min` looked
like it would collapse onto `nlevels_nod2D` under the vertex bottom — both mean "the range where the
cell is a full prism", surely. On pi they differ by **mean -4.0 levels, min -30, at 88% of nodes**.
Worse, at `oce_muscl_adv.F90:303` the array is not a range at all but the **only wetness guard** on
```fortran
edge_up_dn_grad(1:2,nz,edge) = tr_xy(1, nz, edge_up_dn_tri(:,edge))
```
which has no wetness test of its own and reads `tr_xy` — an `intent(out)` array written only over
`ulevels(elem)..nlevels(elem)-1`, so uninitialised heap below an element's bottom. Substituting the
node column would have read past the up/downwind triangles' bottoms on ~85% of pi's interior edges.
**Before retiring an array as redundant, grep its consumers for the one that is using it as a guard
rather than as a bound.**

**3. A latent inequality flipped and took `helem` with it.** `update_thickness_ale` rebuilt `helem`
one layer short of the element's bottom, which was safe because `nlevels_nod2D_min(n) <= nlevels(e)`
meant no node ever stretched that layer. The vertex bottom reverses it to
`nlevels(e) <= nlevels_nod2D(n)`, so the deepest layer *is* stretched at the element's deeper nodes
while `helem` stays frozen — silently putting the SSH/continuity budget and the tracer volume out of
agreement under zstar. **When a scheme inverts a producer, enumerate the inequalities that held
between the old and new quantities; the loop bounds that quietly relied on them are the bugs.**

**4. The regression net was a step counter.** `fesom_analytic` was written into the plan as "the
single most valuable free check" — flat-bottom, so bit-identical across the change. It runs no
physics: `model_step` is `model%nsteps_done = model%nsteps_done + 1` (`mod_step_oce`'s real step is
`model_ocean_step`, which that driver never calls). Five tasks were going to gate on it. **Open the
driver and read its step function before calling anything a regression net.**

**What replaced it.** A conservation gate (`tools/run_conserve_pi.sh`) built and baselined *before*
any behaviour changed, checking total heat/salt content plus three invariants per step. Building it
first is what made every later task measurable — and it earned its keep immediately: reverting just
the `helem` loop bound reproduces
`helem VIOLATED elem 1 nz 18 helem/mean= 6.0000000000000000E+01 6.0000195317207513E+01`.

**5. A gate that never runs the code path proves nothing — check the feature flags.** The
`tr_xynodes` denominator was the one change flagged as a physics judgement, yet every gate run
against it had `Redi = .false.`: `diff_ver_part_redi_expl` is its only caller, it is guarded by
`if (Redi)` (`oce_ale_tracer.F90:225`), and the driver inherited the `.false.` default because it
was derived from a *reduced-config* lifecycle driver that had no GM/Redi knobs. Zero coverage, and
nothing failed to say so. **After building a gate, list the `if (feature)` guards around the code
you changed and confirm the gate actually turns them on.**

And when it did run: conservation still could not tell right from wrong, because the Redi tendency
is a telescoping flux divergence (`del_ttf += (vd_flux(nz)-vd_flux(nz+1))*dt/areasvol`) that
conserves for *any* `tr_xynodes`. The discriminating test is a value test, not a budget test — feed
a constant element gradient and require the node average to return it: wet-area `4.4e-16`, `areasvol`
off by `8.6e-01`. **A conservative operator launders a wrong coefficient; budget checks cannot see
inside one.**

**And one thing that is not a trap.** `linfs` drifts `-2.2e-04` in heat over 20 unforced pi steps
and that is correct behaviour, not a leak: the linear free surface freezes `hnode`, so the surface
vertical advective flux `-w*T*area` at `nzmin` (`oce_adv_tra_ver.F90:66`) is a real source/sink with
no thickness change to balance it. The drift is non-monotone — a free-surface adjustment transient.
Gate conservation on **zstar** (round-off, `~1e-15`) and run `linfs` for the invariants only.

## L54 — A ported scheme with no oracle needs an EXACT analytic test, not a plausible one

`momadv_opt==1` (vector-invariant momentum advection, `oce_dyn_vinv`) has **no FESOM2
oracle**: the v2.7.3 branch exists but aborts with *"not adapted mom_adv advection typ for
ALE"*, and the only implementation (`/home/a/a270029/qq/oce_vinv_mom_adv.F90`) is pre-ALE
research code that was never tested upstream. So the usual byte-gate is unavailable and
correctness has to be designed in.

What worked, and what did not:

- **A tautology is worse than no test.** The first design gated on
  `u . [(f+zeta) x u] == 0`. That is identically — in fact *bitwise* — zero for ANY value
  of `zeta`, because `(f+zeta)*zhat x u` is perpendicular to `u` by construction. It would
  have reported `0.0` forever and manufactured confidence. Before trusting an invariant,
  check it can actually FAIL.
- **An exact analytic case is worth more than a global budget.** `u = alpha*y` on a
  `cartesian=.true.` mesh gives `zeta == -alpha` to relative 5e-13 over 81729 node-levels.
  That one case pins the magnitude, the sign, AND the `area(nz,n) -> area(n)` deviation,
  because `zeta = circulation/area` only reproduces `-alpha` if `area(n)` equals the sum of
  `elem_area/3` over the node's WET adjacent elements.
- **Build analytic fields in `coord_nod2D`, not `geo_coord_nod2D`.** Meshes are read with
  `force_rotation=.true.`, so the geometry (`elem_area`, `edge_cross_dxdy`, element
  centres) lives in the ROTATED frame. A field linear in geographic latitude is not linear
  there; the test failed worst at the rotated pole (mean `zeta/(-alpha)` = 0.836) and the
  cause looked like a metric bug in the operator for some time.
- **`cartesian=.true.` is a legitimate unit-test device.** It sets `elem_cos = 1` AND
  `ay = 1` in `compute_elem_area`, so the metric stays self-consistent and the median-dual
  contour closes; the discrete curl then becomes exact for a linear field. On the spherical
  metric `edge_cross_dxdy(1:2)` carries `elem_cos(el1)` and `(3:4)` carries
  `elem_cos(el2)` at the same edge midpoint, so no exactness assertion is possible.
- **Design the test so the known bug fails it.** `qq` sets `uvert(nl1+1) = 0`, which under
  bottom-at-vertices leaves a spurious bottom tendency (`nlevels(e)` is the `min` over the
  element's nodes, while `vert_vel_ale` zeroes `Wvel` only at `nlevels_nod2D(n)`). The
  originally planned test — "`w == 0` gives zero contribution" — is blind to it, since
  every vertical term is proportional to `w`. The test that sees it is **uniform `u` with
  NONZERO `w`**: 5.4e-2 relative with `qq`'s bottom face, 2.6e-16 with the fix.

## L55 — Element-local geometry beats per-edge arrays at MR; per-rank invariants beat global tallies

Three lessons from the rvo_upwind face-vorticity blend (`oce_dyn_vinv`):

- **Per-edge geometry is owned-edges-only** (`edge_dxdy`/`edge_cross_dxdy`/`edge_len` over
  `nEdgeO`), and an owned element can carry ONE halo edge — so an element-wise consumer of
  per-edge geometry breaks at MR. The safe pattern: rebuild edge vectors from the element's
  OWN vertices (`coord_nod2D` + `trim_cyclic` + its `elem_cos`), and take neighbour VALUES
  through `exchange_elem`-filled element fields rather than neighbour vertex lists
  (`elem2D_nodes` is owned-only). Bonus: the inward normal oriented "toward the third
  vertex" needs no edge_tri left/right convention at all, and `test_vinv` V10 PROVES it
  geometrically (neighbour centroid on the outward side, same metre metric — a positive
  diagonal rescaling does not preserve dot-product signs, so never mix metrics in such a
  test).
- **FESOM element ownership OVERLAPS at partition seams** (`myDim` elements = elements
  touching an owned node), so summing per-owned-element counts across ranks double-counts
  boundary-strip items (pi dist_2: 6 boundary edges counted twice). Write per-rank
  invariants; a cross-rank tally mismatch there is NOT a bug. A ground-truth probe against
  the global mesh files settles such alarms in minutes — run it before redesigning.
- **A stencil-widening feature moves "interior" one ring inward in constant-field tests.**
  The blend reads the NEIGHBOURS' vertex vorticity, so constant-zeta inertness holds only
  where the element AND its neighbours avoid boundary nodes (`elem_int2`). The V9 failures
  at 0.89 relative were the mask, not the physics.

## L56 — The two element halo rings: eDim is the EDGE-sharing ring, eXDim the node-only ring; a receiver-side reading of the halo code raised a false seam alarm that the sender-side definition plus a 2/8-rank probe settled

**Context.** Revalidating `6e38cf6` (rvo_upwind), reading only our side of the halo
machinery: `exchange_elem` fills `com_elem2D` (the eDim block) and never the eXDim block;
`build_elem_adjacency` stores whatever `edge_tri` holds, which can be a halo id; the one
test that walks neighbours (V11) skips every halo neighbour as "no local data". From that
alone the inference was: an owned element's far neighbour across its halo edge may be
eXDim, `omega_e(nz,nb)` stays at its initialised 0 there, the blend silently treats that
neighbour as zero vorticity, and results depend on the partition at seams. No bounds error
would ever fire (the array spans `nElemF`), so nothing in the suite could see it.

**What settled it.** The DEFINITION of the rings lives in the partitioner, not in the
receiver: FESOM2 `gen_comm.F90` builds `com_elem2D` as the non-owned elements "only those
sharing an edge" with an owned element, and `com_elem2D_full` as "all neighbors" (sharing
a node). So an edge-neighbour of an owned element is in eDim *by construction*; an eXDim
element is never across an edge. A probe at np=2 and np=8 (now permanent as `test_vinv`
A5): 17–35 halo-neighbour slots per rank, every one inside the eDim ring, every one
holding its owner's value after one standard `exchange_elem`. No defect. Switching to
`exchange_elem_full` "to be safe" would have bought nothing except larger receive lists
on every momentum step.

**Lessons.**
- **Choose the exchange by the stencil, not by fear.** Element across an EDGE
  (`elem_neighbors`, `edge_tri`) → `exchange_elem`. All elements around a NODE
  (`nod_in_elem2D` patches, MUSCL's `tmp` field) → `exchange_elem_full`. The ring names
  are not "near/far"; they are "edge/node".
- **A halo alarm from the receiver's code alone is a hypothesis, not a finding.** The
  sender-side definition (partitioner) is the ground truth; read it before changing the
  exchange. Same move as L55's global-mesh probe, and again it took minutes.
- **A blind test does not justify changing production — it justifies a test.** V11 was
  blind by design; A5 reuses the exact production pattern (owned fill → one
  `exchange_elem` → read at every `elem_neighbors` id) with a sentinel, and PRINTS the
  number of slots that pass only because the exchange filled them. That printed count is
  the test's teeth; at npes>1 it must be >0 on every rank, which A5 asserts.
- **A guard built on partition metadata must handle the serial identity partition
  explicitly.** In a test that never calls `set_partition`, `partit%myDim_elem2D` is 0 and
  `myList_elem2D` is unallocated at npes==1; the first A5 draft crashed serially for
  exactly that reason while passing at np=2/8.

## L57 — An implicit operator inherits the FORM of the explicit scheme it completes: flux form for a flux-form scheme, advective form for a vector-invariant one — the split `w_i` was silently the wrong one

**Context.** `use_wsplit` moves the part of `w` above the vertical CFL limit into the
momentum TDMA (`impl_vert_visc_ale`), transcribed from FESOM2 in upwind FLUX form:
`a += min(0,wu)·zinv, b += max(0,wu)·zinv − min(0,wd)·zinv, c −= max(0,wd)·zinv`, row
sums `(wu − wd)·zinv`. Applied to a uniform `u` that is `−u·∂w_i/∂z`. The scalar scheme
is flux form everywhere, so the piece is balanced by the horizontal `u·∇·u` through
continuity. The vector-invariant scheme is advective form everywhere (`ζ×u + ∇KE`;
Block C builds `w_e ∂u/∂z` as `δ(w_e u) − u δw_e`), and nothing balances the implicit
`u ∂w_i/∂z`: wherever the split is active, a uniform column with divergent `w_i` is
accelerated — O(1) relative to the vertical advection there, invisible everywhere else,
and invisible in every production driver because they all pin `use_wsplit=.false.`.

**The fix is a row-sum correction, not a new operator.** Subtract `(wu − wd)·zinv` from
the diagonal (the bottom row drops `wd`, so there `wu·zinv`): the operator becomes the
upwind ADVECTIVE form with the inflow-face velocity, row sums vanish, constants are
preserved implicitly, diagonal dominance is kept, and the Thomas sweep is unchanged.
Gated on `momadv_opt == 1`, so the FESOM2-pinned path is structurally untouched.

**Lessons.**
- **Check the row sums.** Flux form ⇔ row sums `= δw·zinv`; advective form ⇔ row sums
  `= 0`. That one line of algebra tells which form an implicit operator is, faster than
  reading its stencil. A uniform-`u` test with divergent `w_i` is the executable version.
- **Operator splitting must split one operator.** Explicit `w_e` part and implicit `w_i`
  part only add up to the scheme if both have the same form. When a new explicit scheme
  is dropped in (momadv_opt=1), every implicit completion it relies on (here the TDMA)
  must be re-derived, not inherited.
- **Pin both branches with closed forms.** A single-face `w_i` gives exact per-column
  solutions for BOTH forms (I3/I4); the flux-form closed form is the regression pin that
  the untouched path really is untouched, which the gate on `momadv_opt==2` alone
  cannot say (it only sees conservation).
- **Run the gate config before claiming the gate.** The planned `wsplit` conserve
  configs (with a split-face count so they could not be vacuous) died before the first
  step: `use_wsplit=.true.` error-stops in `do_oce_adv_tra` because the tracer-side
  `adv_tra_vert_impl` is unported. The momentum fix is therefore pinned at unit level
  only (I1-I4) and is dormant in every driver; the end-to-end gate belongs to the
  `adv_tra_vert_impl` port, and the driver hook was reverted rather than left as a
  knob that always aborts. (Superseded 2026-10-02: `adv_tra_vert_impl` is ported, the
  hook is back as `FESOM3_WSPLIT`, and the `wsplit` / `momadv-vinv-wsplit` gate configs
  run — L58.)

## L58 — "Smooth vs hard" CFL switching is ONE parameter, not two schemes; the FCT low-order explicit/implicit upwind pair is EXACTLY split-invariant, so a smoothness test needs operators that differ at first order; the constancy test is the consistency proof of an operator split; and a cold-start gate needs a MEASURED cap or it is vacuous

**Context.** `use_wsplit` (FESOM2 `compute_Wvel_split`) sends the part of `w` above a
vertical Courant cap to the implicit upwind solves (momentum TDMA, tracer-diffusion TDMA,
FCT `adv_tra_vert_impl`). FESOM2's split is a hard switch (`w_e = w·C/CFL_z` above
`C = wsplit_maxcfl`): `d(w_e)/d(CFL_z)` jumps from 1 to 0 at `C`, so a face oscillating
around the threshold flips between high-order explicit and partly first-order implicit
treatment. The port (`docs/plans/completed/2026-10-02-wsplit-smooth.md`, commits
431e3e7..bf73965) replaced it by Shchepetkin's (2015, Ocean Modelling 91, Sec. 3.1) C¹
limiting function `oce_wsplit::wsplit_implicit_fraction(Cu, Cu_min, Cu_max)` — `f = 0`
below `Cu_min`, the bend `x²/(F + x²)` with `x = Cu − Cu_min`, `F = 4·Cu_max·(Cu_max −
Cu_min)`, and `f = (Cu − Cu_max)/Cu` above `Cu_cut = 2·Cu_max − Cu_min` — and ported the
one missing consumer, the FCT `adv_tra_vert_impl`, so `use_wsplit=.true.` runs end to end
for the first time (L57's dormant momentum correction is live). There is no FESOM2 oracle
for the smooth split; the tests are the specification (L54).

**What settled it.**
- *One parameter.* At `Cu_min = Cu_max`, `F = 0` and `Cu_cut = Cu_max`: the function
  collapses to FESOM2's switch in exact arithmetic. The hard switch is the corner `(1,1)`
  of the parameter space, not a second scheme — it needs no mode flag and serves as the
  positive control of the smoothness test (`test_wsplit` X2). Plumbing the new parameter
  through the old `dd` branch first (TDD stage 2) failed S1/S2 by 2 ulp / non-bitwise and
  X1 by 10x — the same numbers the final positive control reproduces.
- *Split invariance of the LO pair.* Explicit upwind with `w_e` followed by implicit
  upwind with `w_i` across ONE face leaves the donor cell at
  `T₂·(h₂ − dt·w_e)/(h₂ − dt·w + dt·w_i) = T₂` and the receiver at `T₁h₁ + dt·w·T₂`, for
  ANY `f`. The FCT low-order pair therefore has no first-order kink to show: the smoothness
  test the plan first put there would have passed for the hard switch too (a tautology,
  L54). The kink is first order only where the explicit and implicit operators differ at
  first order — the non-FCT path, QR4C (centred) on `w_e` plus the upwind TDMA on `w_i`,
  where the tendency's derivative with respect to `w` jumps by `J = [UPW − QR4C](T)` for
  a `T` with curvature. `test_wsplit` X runs there (`Kv = 0`, quadratic `T`, 601-point
  `Cu` scan): second-difference ratio `max|D2|/(ΔCu·max|D1|)` 0.96 / 1.17 (`w ≷ 0`) for
  the smooth function, 100.3 / 64.1 for the `(1,1)` control (bound 10), and the control's
  kink equals the closed form `J = cq·h²·(1+Cu_c)/2` per unit `Cu` to 0.25 %
  (`D2/(J·ΔCu) = 1.0025`; the factor `(1+Cu_c)` is the explicit pre-advection
  `T* = T + w_c·P_q T` the TDMA then acts on — the shorthand `cq·h²/2` would have read
  2.005, and the 0.25 % is the predicted `ΔCu/(1+Cu_c)` curvature of the implicit branch).
- *The constancy test is the proof.* `hnode_new` is advanced by the FULL `w`; for uniform
  `T` the explicit LO step with `w_e` gives
  `T* = T·(1 + dt·δw_i·area/(areasvol·hnode_new))` and the implicit step solves
  `(hnode_new + dt·δw_i·area/areasvol)·T^{n+1} = hnode_new·T*`, so `T^{n+1} = T` exactly,
  for any split with `w_e + w_i = w`. One test (`test_wimpl_tra` C2/C2b: ≤ 5.3e-16 for
  three parameter sets, both signs, non-uniform layers, open surface) pins the split
  identity, the flux-form coefficients, the unsigned surface row and the `hnode_new`
  usage at once.
- *The wiring bug a conservation gate cannot see.* The antidiffusive flux must be
  `HO(w) − LO(w)`, the LO flux RECOMPUTED on the full `w` after the implicit step;
  `HO(w) − LO(w_e)` double-counts the `w_i` transport yet still telescopes, so heat and
  salt conserve to round-off either way. `test_wimpl_tra` C6 calls `do_oce_adv_tra` itself
  and compares with the test's own assembly from the public parts: bitwise equal at
  np 1/2, while the wrong assembly differs by 29.9 / 58.8 K·m — exactly the predicted
  double-counted transport `(−w_i(nz)·T(nz) + w_i(nz+1)·T(nz+1))·area·dt/areasvol` to
  5e-14 — and the whole step matches a closed form to 1.4e-14 K.
- *A cap must be measured.* The 20-step pi cold start from rest has max `CFL_z` =
  2.76e-2 (owned faces above 0.005 / 0.01 / 0.02 / 0.04: 2568 / 499 / 26 / 0). At the
  default cap 1.0 a `wsplit` conserve config splits NOTHING and would have passed as a
  no-op; with `FESOM3_WSPLIT_EXPECT_SPLIT` set it exits 1 (the positive control, now the
  ctest `fesom_conserve_zstar_wsplit_vacuous_np1`; the guard is value-based: ≥ N split
  faces and ≥ 1 at `f ≥ 0.5`). The gate cap 0.005 (onset 0.0025) gives 8632 split faces,
  1130 on the capped branch `CFL_z ≥ Cu_cut = 2·maxcfl − mincfl = 0.0075` (where `Cu_e`
  saturates — NOT at `maxcfl` itself, where a face is still on the bend), 500 at
  `f ≥ 0.5`, max `f` 0.82; 0.01 would have left 26 faces at `f ≥ 0.5`, a margin a small
  dynamics shift erases. `fesom_conserve` prints one statistics line for every config
  with exactly those three counts (`w_i ≠ 0` is by construction the count above
  `mincfl`, `f ≥ 0.5` the count above `2·maxcfl`; a first version printed both members
  of each pair and labelled the `> maxcfl` count "the cap" — the review caught it).

**Lessons.**
- **Before adding a "mode", check whether it is a corner of the parameter space.** The
  hard switch is `Cu_min = Cu_max`; "smooth vs hard" is one knob. A degenerate corner is a
  free positive control and needs no flag.
- **A smoothness (or any sensitivity) test needs operators that differ at FIRST order.**
  Work out what the kink is proportional to before choosing the path: the obvious pair
  (the LO explicit/implicit upwind pair) is exactly invariant and would have passed for
  the broken scheme. Derive the closed-form jump and require the positive control to
  REPRODUCE it (1.0025·J), not merely to fail the bound.
- **A free-stream/constancy identity is the consistency proof of an operator split.** It
  holds for every split that sums to `w` and only with the right `hnode_new`, so it fails
  on any coefficient, surface-row or thickness mistake. Run it first; a conservation test
  proves less (column sums alone).
- **A conservation gate is blind to telescoping wiring bugs.** Assert the ASSEMBLY against
  an independent assembly from the public parts — bitwise where the parts are the same
  calls — and make the wrong assembly's difference equal its predicted value.
- **A gate config that may fire nothing must prove it fired.** Measure the quantity the
  switch keys on (max `CFL_z` 2.8e-2 here), set the threshold from the measurement with a
  margin, print the counts, and make an empty count a FAILURE (`EXPECT_SPLIT`).
- **Tolerance lessons from the review.** (i) Test C⁰ at a joint by evaluating BOTH branch
  formulas at the joint and requiring agreement (`test_wsplit` W4): a value check of the
  function alone cannot tell continuity from a lucky branch selection. (ii) At a joint
  whose closed-form slope is 0 a relative tolerance is meaningless; bound the one-sided
  difference quotients ABSOLUTELY (`|q| ≤ 2h/F`, W5). (iii) Above the cap `1 − f =
  Cu_max/Cu` is formed by cancellation, so the error of `Cu·(1−f)` scales as `eps·Cu`,
  not with ulps of `Cu_max`: the first "4 ulp of `Cu_max`" bound failed at `Cu ≈ 17.5`
  (measured 1.23 eps·Cu), the error-model bound `4·eps·Cu` passes with a 3x margin (W6).

## L59 — The cost was the RULE, not the storage: keeping FESOM2's node-averaged fallback made the on-the-fly MUSCL barely faster; adopting the reference RK3 rule (triangle gradient or zero) made it ~43 % faster

**Context.** Goal: drop FESOM2's stored per-edge MUSCL gradient `edge_up_dn_grad` as the
reference RK3 code does (`qq/oce_stepRK3.F90`). First attempt: keep results bit-identical,
so the Miura node average FESOM2 uses at coasts / bathymetry steps became its own per-node
array `gnod`. Measured on core2: the first version was ~20 % SLOWER than the stored array;
after shaping the lookup (one branch per level) and computing `gnod` only where read it was
~9 % faster. The user then asked for the reference implementation proper: triangle gradient
where the up/downwind triangle is wet, ZERO increment otherwise -- no node average at all.
That kernel is ~43 % faster than FESOM2's fill + kernel (np 4: 42 vs 74 ms), and the
results differ from FESOM2 only at coast/step levels (core2: max |Δflux| 0.4 % of max).

**Lessons.**
- **Profile what the stored array encodes, not just that it is stored.** The expensive part
  was the node-averaged fallback (gathers over the node's elements, per edge in FESOM2 /
  per node in the bit-identical port), not the array itself.
- **Bit-identity is a choice with a price; state it as one.** Equivalence pinned the old
  rule's cost into the new code. When the user wants the reference scheme, test it EXACTLY
  instead: the production kernels must equal the stored-array kernels fed with the new
  rule's array (bitwise), and must equal FESOM2 bitwise where the two rules coincide, and
  differ elsewhere (teeth) -- test_muscl_onthefly parts R and S.
- **pi is too small to decide performance:** the 13 MiB array sits in cache there; core2
  changed the ranking of the intermediate versions. Time on the production mesh (the test
  binary takes `FESOM3_MESH_DIR`; run it with `ulimit -s unlimited`).
- **A bit-identity test on synthesised cavities caught a wrong claim within minutes:**
  FESOM2's `nl12 = 0` at boundary edges makes segment (D) start at level 1, above a cavity
  roof (the velocity is 0 there).

## L60 — A "roundoff" regression gate must make the solver converge to roundoff, and a cross-np comparison must start from one common state: the CG tolerance and the IC extrapolation each hid the roundoff behind 1e-4..1e-2 differences

**Context.** Hybrid-mesh plan, Task 7 follow-up. The first triangle regression gate (new
geometry vs frozen pre-Task-7 binary) could only be 1e-2 after the first SSH solve, and the
frozen binary itself differed between pi np 1 and np 2 at step 1 (density 4.7e-3).

**Lessons.**
- **An iterative solver stopped at a relative residual turns roundoff into tolerance-sized
  noise.** The SSH CG stops at soltol = 1e-5: the SAME binary at 1e-5 vs 1e-14 differs by
  3.5e-4 in d_eta on pi. Solve both runs to 1e-14 (env `FESOM3_SSH_SOLTOL`, ~3x iterations)
  and pi steps 1-3 agree to <= 1e-12. Before trusting a "roundoff" gate, run the reference
  against itself with the tolerance changed: that shows the noise floor the gate can see.
- **What is left is physics amplifying roundoff, and it is local in time and space.** CORE2
  T differs by 2e-9 after step 1 only at ice-covered nodes (EVP at rest is stiff); pi jumps at
  step 4, the first forcing record refresh. Find the event, then tier the gate (early steps
  tight, later steps / ice / TKE-driven fields at their measured floors) instead of one loose
  tolerance for everything.
- **A cold start is not partition-independent when the IC is built by an order-dependent
  sweep.** FESOM2's `extrap_nod3D` Gauss-Seidel fills land-masked IC points with per-rank halo
  snapshots, so np 1 and np N start from different T/S at exactly the extrapolated points
  (pi: 513 (node, level), all of them dummy after the bilinear step). Check the step-0 fields
  before suspecting halo bugs; compare decompositions from ONE checkpoint (restart at np 1 and
  np N) -- then np 1/2/8 agree to 1e-14 at step 1.


## L61 — Login-node MPI (Open MPI 4.1.2, ob1 + vader + knem) delivered STALE payloads: its knem registration cache misses glibc heap trims (an MPI-library bug, not a halo-code bug)

**Context.** Hybrid-mesh plan, Task 10. `test_vinv` on the NEMO quad cut-out at np 2 failed
V9b/V11/V12 with garbage-sized errors while np 1, the mixed cut-out and pi passed; a
bounds-checked build was clean and the failure was deterministic. Per-call logging of
`core_3D_r` showed every count matching and the sender's buffer correct, yet the receiver got
the payload of the exchange FOUR calls earlier (the previous KE exchange).

**Root cause (Task 10 follow-up, proven).** `pml ob1` sends a message above 32 KB
(`btl_vader_rndv_eager_limit`) as a copied 32 KB prefix plus a single-copy RGET of the rest; with
knem the receiver's kernel copies from a knem region the SENDER registered. vader caches these
registrations in an rcache keyed by virtual address (`mpi_leave_pinned` = auto -> on, because the
memory patcher claims hook support) and relies on memory hooks to invalidate them. glibc's
`free()` shrinks the heap top with `brk` (trim threshold), the next `allocate` re-grows it at the
SAME virtual address with fresh pages -- and the cache still hands out the old region, whose
pages knem keeps pinned with the data of the last message sent from there. Right count, correct
first 32 KB, stale tail. Evidence:
- `tools/run_mpi_knem_stale_repro.sh` (standalone, `tools/mpi_knem_stale_repro.F90`: contiguous
  allocatable buffers passed by first element, `MPI_Waitall` before `deallocate`, one message per
  pair and call): with `MALLOC_TRIM_THRESHOLD_=0 MALLOC_TOP_PAD_=0` 788 of 800 receives are stale
  under knem; 0 without the trim, 0 with `single_copy_mechanism none` or `emulated`, 0 with knem
  and `--mca mpi_leave_pinned 0` (no cache), 0 under `pml ucx` (`UCX_TLS` mm,knem / mm,cma,
  `UCX_RNDV_THRESH=1024`). An explicit `madvise(MADV_DONTNEED)` from user code IS caught (0 stale);
  strace shows the FESOM test does 45 `brk` and no `madvise` per rank.
- `test_vinv` nemo025q_cut np 2 under knem: FAILS by default, PASSES (1e-16) with
  `MALLOC_TRIM_THRESHOLD_=100000000000`, with `--mca mpi_leave_pinned 0`, with `--mca memory
  ^patcher`, with single-copy `none`/`emulated`, and under `pml ucx` (incl. knem, rndv 1 KB).
  (`cma` is not offered on Levante: `kernel.yama.ptrace_scope`; the value is ignored -> knem.)
- The code is clean: every non-blocking call in `src/` is in `mod_halo`'s four cores (contiguous
  local buffers passed by element, alive until after `MPI_Waitall(nreq, ...)` over exactly the
  posted requests, no exchange overlaps another, so the shared tag is matched in order). A build
  with `-check arg_temp_created -warn interfaces` creates NO argument temporary anywhere while
  running `test_vinv` np 2 and the production driver on pi np 2 and CORE2 np 8.

**Lessons.**
- **"Right count, payload of an earlier message of the same size" has two classic causes**: a
  copy-in temporary at a non-blocking call (code) or a stale registration cache (library). Tell
  them apart with `-check arg_temp_created` and with the transport knobs that disable only the
  cache (`mpi_leave_pinned 0`, `MALLOC_TRIM_THRESHOLD_` huge) -- a code bug survives those.
- **Every login-node `mpiexec --mca pml ob1 --mca btl self,vader` run needs
  `--mca btl_vader_single_copy_mechanism none`** (or `OMPI_MCA_btl_vader_single_copy_mechanism=none`,
  which `env.sh` has exported since M2.12 for the SAME bug: L35's `tr_xy` "correct for 134656
  bytes, then garbage" was this stale tail). The ctest launcher (`FESOM3_MPI_LAUNCH_FLAGS`),
  `tools/run_tri_regression.sh`, `run_conserve_pi.sh` and `run_bottom_audit.sh` pass it
  explicitly because they run without `env.sh`. `mpi_leave_pinned 0` is an
  equivalent alternative. It hits ANY point-to-point or collective payload above 32 KB, not only
  halos; reusing persistent buffers in `mod_halo` would hide it there only.
- **Batch runs are not affected**: `env.sh`/`env/levante.dkrz.de/shell.intel` select `pml ucx`
  (vader unused, `btl self`), and UCX's own memory hooks catch the trim (tested intra-node with
  its knem and cma transports; inter-node IB not testable on a login node).
- **It hid for months** because whether a buffer lands on trimmed-and-regrown heap pages depends
  on the allocation history of the whole program (2-D fields, pi np 2 rarely hit it); earlier
  "transient MPI failures" of login-node ctest runs without `env.sh` may have been this.


## L62 — A faithful port of a MATLAB mesh script is only byte-checkable if it keeps MATLAB's precision rules: float32 inputs stay float32 and every single/double mix is computed in single ("single wins")

**Context.** Hybrid-mesh plan, Task 1 (`tools/mesh/nemo_meshq.py`, commit `d4c6f8d`). The MATLAB
outputs of `nemo_meshq.m` (`nod2d.out`, `elem2d.out`, `depth.out`, saved to
`/work/ab0995/a270029/meshes/nemo025q/ref_matlab/` before anything could overwrite them) were the
byte oracle. The NEMO inputs are float32 (`nav_lon/lat/lev`, `glamf/gphif`), int16 (`mbathy`) and int8.

**Lessons.**
- **Port the arithmetic types, not only the statements.** In MATLAB `single op double` is single, so
  the `levels` recursion and the f-point corrections run in float32; numpy promotes to float64 and
  prints different 4th decimals. Emulating "single wins" made `elem2d.out` and `depth.out`
  byte-identical and left ONE line of `nod2d.out` (node 711389, a fold point: -83.6347 vs -83.6348).
- **The last ulp of a transcendental function is not portable -- list it, do not chase it.** That
  line is 1 float32 ulp of `atan2` straddling the print rounding; correctly-rounded cos/atan2/acos
  variants and double-precision semantics matched WORSE. The test (`tools/mesh/test_nemo_meshq.sh`)
  allows <= 20 listed coordinate lines differing by <= 1e-4 and requires the rest byte-identical.
- **Save the oracle first and keep order-dependent loops sequential.** Land removal edits `bathy`
  while scanning; quad creation and node compaction define the numbering -- vectorising them changes
  the mesh, not just the speed. Order-independent loops may be vectorised.


## L63 — One code path for triangles and quads (4 slots, triangle 4th = 1st) removes every branch only if each USE of the repeated slot is classified: weighted sums and max/min are safe, counts and per-vertex updates are not -- and the boundary convention changes everywhere in one commit

**Context.** Hybrid-mesh plan, Tasks 6-16. Storage: `elem2D_nodes(4,·)` with a triangle repeating
its 1st vertex; degenerate-slot conventions so that no kernel needs `if (nv == 3)`: gradient slot 4 =
0, kite weight `elem_w(4)` = 0, a missing neighbour = the element itself, a missing `edge_tri(2)` =
`edge_tri(1)` with `edge_cross_dxdy(3:4) = 0`, a triangle's side 4 has zero length (normal 0).

**What it took.**
- **Weighted sums are safe** (slot 4 adds exactly 0: gradients, `elem_w` means, the SSH and ice-TG
  scatters, `stress2rhs`), **max/min are safe** (the repeated value changes nothing: FCT bounds,
  `mo_convect`), but **anything that counts or updates per vertex is not**: the ice FCT limiter over 4
  slots visits vertex 1 twice (mutation check: `test_ice_mass` F5 fails by 0.15-0.25), the IC
  extrapolation would average vertex 1 twice. Those loops run over `elem2D_nnodes` (Task 16, Task 11).
- **The boundary tests must change in the SAME commit as the convention** (Task 7): every
  `el(2) > 0`, `el(i) < 1`, `nb <= 0` in src AND tests became `el(2) /= el(1)` / `nb == el`;
  half-converted, a boundary edge counts twice. The grep afterwards must find only the convention
  builders.
- **A grep for element means (`/3`, `1:3`) does not find node gathers.** FESOM2 writes the node gather
  as `sum(A_e*x_e)/(3*area)` -- the 1/3 is the triangle kite weight hidden in a `3*` (`UVnode`,
  `w*du/dz`, `tr_xynodes`, `compute_sigma_xy`, `smooth_nod`, the node area itself). The discovery
  grep had to look for `3*` and `/(3` too; they became `nod_in_elem2D_w` gathers.
- **Invariant checks in drivers are formulas too**: `fesom_conserve` asserted `helem == sum(hnode)/3`,
  wrong on quads once zstar `helem` is the `elem_w` mean (Task 12); the coupling assigned the 4-slot
  `elem2D_nodes(:,e)` to `elnodes(3)` (shape bug, Task 13).

**How to apply.** For a unified-slot data structure write the per-use rule down (sum / extremum /
count) and review each loop against it; add a mutation test that runs the "wrong" variant on a quad
or mixed patch -- the triangle patch alone does not see these bugs.


## L64 — When bit-identity is retired, freeze the reference at the last bit-identical commit, and raise a tolerance only after a control run without the amplifier (mEVP for EVP, PP for TKE)

**Context.** Hybrid-mesh plan, Tasks 6-17. From Task 7 on triangle results change at roundoff, so
the FESOM2 byte gates were retired and `tools/run_tri_regression.sh` compares against a frozen
binary. Twice a level had to move: CORE2 EARLY 2.1e-9 -> 2.2e-8 when the EVP became 4-slot (Task 15),
and the cut-out 1-vs-N reached 1e-4 after 5 days with TKE (Task 17).

**Lessons.**
- **Freeze before you change the arithmetic, and change the reference only by the test knob.**
  `frozen_task6` = all binaries of `5b2a7f3` (last bit-identical state, proven by every dump and log
  of pi/CORE2 runs being byte-identical to the previous commit); `frozen_task6_soltol` = the same tree
  + ONLY the `FESOM3_SSH_SOLTOL` patch, checked byte-identical to `frozen_task6` with the knob unset.
  Build everything in scratch trees; the user's `build_intel_dp` runs production jobs directly from
  `bin/` and was never touched.
- **A tolerance increase needs a control that removes the suspected amplifier.** Same case with mEVP
  in both runs: CORE2 EARLY 1.1e-15 (vs EVP 2.2e-8) -- the converted ice code is right, the standard
  EVP (unstable for ice at rest, Task 15) amplifies roundoff that now enters inside each of the 120
  subcycles. Same case with PP instead of TKE: 1-vs-N <= 1.1e-7 over 5 days (vs 1e-4) -- TKE
  thresholds switch on roundoff. Only then raise the tier, keeping the same margin (x5 in Task 15).
- **Keep the per-mesh tolerance table and the measured levels in the gate's header** and update them
  with each change (`run_tri_regression.sh` header, HANDOFF "Gate migration for hybrid meshes"); a
  gate whose measured level is not recorded cannot tell a regression from noise.


## L65 — An exact null mode passes every consistency test: test the operator's kernel directly -- the Green-Gauss checkerboard is a null mode on EVERY all-quad mesh, not only on parallelograms

**Context.** Hybrid-mesh plan, Task 9 (`test_ssh_stiff` S8) and Task 17. The plan expected an
"almost exact" checkerboard null mode of the SSH operator on the nearly rectangular ORCA quads.

**What the kernel test showed.** The Green-Gauss coefficient of vertex j is proportional to
x_{j+1} - x_{j-1}, so c1 = -c3 and c2 = -c4 on ANY quad: the pattern (+1,-1,+1,-1) has zero gradient
in every quad, and on any 2-colourable all-quad mesh it is an EXACT null mode of the stiffness part
(|K f|/|K| = 1.8e-16 on the NEMO quad cut-out, not just on the parallelogram patch). Linear-exactness,
symmetry, row-sum and convergence tests (S1-S7) all pass for such an operator. Only the mass part
controls the mode. In the 10-day cut-out runs the SSH null mode stayed flat, but the TRACER grid-scale
energy grew (SST x6, SSS x11-12, faster at quad and interface nodes) -- a design question for the user.

**How to apply.** For every new discrete operator on a new element type, apply it to the candidate
spurious patterns (checkerboard, stripes, per-element alternation) on a real mesh of that class and
assert the result, so the property is known and cannot change silently; then measure the same
patterns in a run (`tools/cutout_diag.py checker`: global 2-colouring projection, per-quad energy,
linear-exact neighbour residual by node class).


## L66 — Exact ingredients do not make an exact scheme on mixed meshes: test the ASSEMBLED node tendency per node class (all-triangle, all-parallelogram, interface, general quad)

**Context.** Hybrid-mesh plan, Tasks 10-11. The plan asked for "linear tracer advected exactly" and
"exact linear shear" on mixed patches. Every ingredient was exact: element gradients (<= 1.2e-13), the
MUSCL face values at the edge midpoints on every edge (<= 2.2e-12), the element velocities.

**What happened.** The median-dual edge scheme assembles the node tendency from face values at the
edge MIDPOINT, while the dual face runs midpoint -> centroid; the errors cancel edge by edge only when
the dual cell is all triangles (Galerkin-P1 equivalence) or all parallelograms. At triangle-quad
interface nodes and on general quads the node tendency of a linear tracer in uniform flow has an O(1)
local relative error (0.22 NEMO-mixed interface, 0.24 jittered patch). Same for the vorticity of a
linear shear: exact on triangles and parallelograms, 6e-2 (first order) on jittered quads.

**Lessons.**
- **Split exactness tests by node class** and assert roundoff only where the scheme is exact
  (`test_updown_mixed` L2, `test_vinv` V2/V9a), with a first-order bound elsewhere -- one global
  tolerance either fails for a scheme property or hides a kernel bug.
- **Test patches must contain the irregular case on purpose**: the NEMO quads are lon/lat rectangles,
  so the real meshes alone never show it; the jittered `generate_mixed_patch` does.
- A scheme property found this way is a finding for the user (remedy = design decision, e.g. face
  values at the dual-segment midpoints), not something to "fix" inside a conversion task.


## L67 — A consistent mass matrix brings its own lumped mass: an FCT must divide by the row sums in EVERY step and conservation must be measured with the same weights -- on quads they are not the geometric node area

**Context.** Hybrid-mesh plan, Tasks 16-17 (the user's specification). Ice FCT element mass
contributions: triangle S/6 diagonal, S/12 off-diagonal; quad (S/36)(4,2,1,2) rotated per row. Row
sums: S/3 per triangle vertex (= the kite node area), S/4 per quad vertex (NOT the kite share of a
non-square quad). FESOM2 divided by the node area `mesh%area` in the low-order step, the high-order
iteration and the antidiffusive fluxes.

**Lessons.**
- **Store the row sums once (`ice%work%fct_lumpedmass`) and use them in all FCT steps**; the
  antidiffusive operator is consistent - lumped (`ice_acoef`: triangle (-2,1,1)/12, quad (-5,2,1,2)/36),
  as compile-time tables per vertex count, so both element types share one code path and no runtime
  `/12`, `/36` appears (L7). Exchange it to the halo: owned elements divide by it at halo vertices.
- **Conserved is sum(lumped*m), not sum(area*m).** The conservation gate (`fesom_conserve` with
  `FESOM3_CONSERVE_ICE`, seeded ice) measures the lumped sums; the mutation "low-order step with
  `mesh%area`" drifts by 1e-6 to 7e-6 and fails, the correct code drifts <= 2.8e-15. For triangles
  lumped == node area, so triangle meshes could never have shown the difference.
- **Check the CSR column search**: the mass matrix lives on the SSH pattern, and a missing column used
  to be silent; it now stops with global ids (`test_ice_mass_missing_col_np1`).


## L68 — Element geometry in lon/lat needs ONE cyclic-wrapped embedding per element (relative to vertex 1): per-pair trimming and FESOM2's plain vertex mean are wrong across the rotated dateline and at the pole

**Context.** Hybrid-mesh plan, Tasks 7 and 14. The 4-slot geometry (`src/mesh/mod_mesh_areas.F90`)
computes centroids, areas, kite weights and Green-Gauss gradients in rotated lon/lat.

**What happened.**
- Trimming each vertex-pair difference separately gave inconsistent embeddings for elements spanning
  > 180 deg in longitude; taking all longitudes relative to vertex 1 (one planar embedding, as FESOM2's
  gradient did) keeps the gradient exact for linear fields and consistent with `elem_area`.
- FESOM2's output element coordinate was the plain unwrapped mean of the 3 rotated vertices: on the
  wrong side of the globe for the 77 (pi) / 566 (CORE2) elements across the rotated dateline, and the
  element velocities there were rotated to geographic with the wrong angle. Output now uses
  `elem_center` (area-weighted, wrapped); the regression compares element vectors there by speed.
- The UNROTATED pi has elements around the geographic pole, where no lon/lat embedding is meaningful:
  geometry tests (`test_geom_mixed`) read pi rotated (`FESOM3_TEST_ROTATE=1`, as the model does);
  the NEMO cut-outs are read unrotated.

**How to apply.** One helper for the element centroid used by geometry AND output; test geometry on
a mesh that crosses the (rotated) dateline; when an output coordinate changes, check that the change
is a > 1 deg jump at a few % of elements and < 1e-10 deg elsewhere before calling it roundoff.


## L69 — A setup tool's exit status and the reader's array sizes are not checks: verify every derived mesh and partition against the MODEL's invariants and limits

**Context.** Hybrid-mesh plan, Tasks 2-6 and 18.

**What happened.**
- `fesom_ini` (FESOM2's `mesh_part`) STOPs with exit status 0 also on errors; only `Total time` in its
  log and `dist_N/rpart.out` tell success. It also KEEPS existing edge/level files, so a changed mesh
  needs them deleted first. `tools/mesh/check_dist.py` re-derives ownership, halos, edge order and
  every communication list (incl. the appended hybrid ones) independently.
- The model reads com lists into fixed `MAX_NEIGHBOR_PARTITIONS = 32` arrays and `read_com` has no
  guard; fesom_ini's new lists are allocatable and unchecked too. ORCA025 needs <= 16 at 128 and 1024
  parts -- fine, but only because it was counted (Task 18).
- Derived meshes broke model invariants the generating tool did not know: the cut-outs left 23/17
  nodes with a layer count above every remaining element (the model's bottom invariant stops in
  `read_mesh`; Task 6 fixed `cut_region.py`), and splitting coastal quads along the shorter diagonal
  made 838 dead-end triangles that FESOM's setup rejects (`split_quads.py` takes the other diagonal).

**How to apply.** After producing a mesh or partition, run `check_mesh.py`, `check_dist.py` AND the
model's own reader on it (ctest `test_mesh_read`); count neighbours against the reader's fixed limits
before using very high part counts; do not trust exit codes of ported Fortran setup tools.


## L70 — A shared element is a duplicated prognostic variable: an order-dependent sum makes the copies differ by roundoff, a branchy operator makes it O(1), and the duplicated border fluxes turn it into a mass leak -- test the global divergence sum at np>1

**Context.** CORE2 runs results7/11/17/21 (2-12 Oct 2026); plan
`docs/plans/completed/20261012-shared-element-sync.md`.

**What happened.**
- FESOM2's partition rule gives an element to every rank that owns one of its vertices; dist_512 has
  24,798 such shared elements out of 244,659. Each owner updated its own copy of `uv` in `update_vel`,
  and nothing made the copies equal again: `exchange_elem(UV)` refreshes only the eDim halo.
- The viscosity filter adds each edge's contribution in the rank's own local edge order, so the copies
  differed by roundoff every step. The `rvo_upwind` blend of the vector-invariant scheme switches when
  the inflow `u.n` through an edge changes sign, which can turn a roundoff difference into O(1).
- The two ranks then computed different fluxes across their common edges, `ssh_rhs_old` no longer
  summed to zero over the globe, and the linear free surface created water: the global-mean SSH
  drifted and the salt trend was +1.4e-5 to +2.3e-5 psu/yr instead of about -2e-6. results5 (before
  `rvo_upwind`) was clean, results7 (with it) leaked; each restart re-synced the copies (every owner
  reads the one written value) and paused the leak for about six months.
- No gate saw it. The conservation gates ran at np 1 (no shared elements) or at np 2 with zstar and
  without the blend for 20 steps, where the roundoff never got amplified. The restart gate's np>1
  "FESOM2-inherent 1-ULP floor" was the same defect seen from the other side and had been accepted as
  a caveat.
- Fix: one designated owner per shared element (the highest vertex-owner rank, whose copy the restart
  writer keeps); `update_vel` syncs `uv` every step. The `uv_rhsAB` and sigma copies then stay bitwise
  equal by themselves (guards `test_vinv` V13, `test_ice_dyn_mixed` D6). `test_shared_elem` S3: the
  divergence-sum ratio is ~1e-17 with the sync and 7e-6 to 4e-4 without it; np>1 restarts are exact.

**How to apply.** Treat every field stored on multiply-owned entities as duplicated state: compute it
once and send it, or test that the copies stay bitwise equal. At np>1, test the global divergence sum
`|sum ssh_rhs_old| / sum|ssh_rhs_old|` and the copy spreads (`FESOM3_STEP_DIAG`, `fesom_conserve`
`FESOM3_DIVSUM_TOL`), not only conservation at np 1. When a gate needs a "harmless ULP floor" at
np>1, find the duplicated state behind it before accepting it.
