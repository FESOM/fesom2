#!/usr/bin/env bash
# M2.11c-1 oracle: run the instrumented FESOM2 (build/bin/fesom.x -> build/lib64/libfesom.so)
# on the CORE2 mesh at SINGLE RANK through its REAL multi-step lifecycle (fesom_init ->
# fesom_runloop(N) -> fesom_finalize), UNFORCED (use_ice=.false. => forcing_setup is a
# no-op, gen_forcing_init.F90:43; all surface fluxes stay 0, arrays_init-zeroed). The
# built-in dump_shim (FESOM_DUMP_FILE) emits the 13 NODE substeps per step inside the REAL
# oce_timestep_ale; dump_shim_finalize closes the file in fesom_finalize. The 1-rank output()
# crash is bypassed by the npes==1 early-return added to io_meandata.F90 (uncommitted shim).
#
# This is the FIRST true time-stepping oracle run (not a prescribe-and-stop shim): it gates
# FESOM3's lifecycle driver (src/drivers/fesom_lifecycle) for the multi-step state evolution
# (AB2 velocity + tracer rotation, state carry-over) on CORE2 with the do_ic3d phc3.0 IC.
#
# Reduced-M2 overrides on work_core (the proven M2 oracle config):
#   namelist.config: which_ALE zlevel->linfs, use_ice .true.->.false., yearnew 1958->1948,
#                    run_length 50->N, ResultPath ./   (step_per_day=48 -> dt=1800 kept)
#   namelist.oce:    mix_scheme KPP->PP, Fer_GM/Redi .true.->.false.
#   namelist.dyn/tra already match (opt_visc=7, momadv_opt=2, use_wsplit=.false., AB_order=2,
#                    MUSCL/QR4C/FCT, i_vert_diff=.true., K_hor=0, use_instabmix, ...)
# Needs the hand-crafted CORE2 dist_1/ (tools/make_dist1.py, M2.11a).
#
#   tools/run_lifecycle_core2.sh [run_dir] [dump_prefix] [nsteps] [nranks]
# [nranks] (default 1) selects CORE2 dist_<nranks>; npes>1 is the production-validation
# multi-rank oracle (the built-in dump_shim is gid-keyed -> per-rank <DUMP>.<mype5> dumps;
# the npes==1 io_meandata early-return is inactive so the normal multi-rank io_gather runs).
set -euo pipefail

F2=/home/a/a270088/port2/fesom2
RUN="${1:-/scratch/a/a270088/lifecycle_core2}"
DUMP="${2:-$RUN/life_f2}"
NSTEPS="${3:-3}"
NP="${4:-1}"
# M4c/M4d: FER_GM=1 / REDI=1 keep the work_core Fer_GM/Redi=.true. (GM bolus / Redi isopycnal
# diffusion ON) instead of the reduced-M2 sed-off. Default 0 = the proven GM/Redi-off gate.
FER_GM="${FER_GM:-0}"
REDI="${REDI:-0}"
# M5b: MIX_KPP=1 KEEPS work_core mix_scheme='KPP' (the production boundary-layer scheme)
# instead of the reduced-M2 sed-to-PP. use_sw_pene stays .false. (unforced => sw_3d=0).
# Default 0 = the proven PP gate.
MIX_KPP="${MIX_KPP:-0}"
# M7b: MIX_TKE=1 sets mix_scheme='cvmix_TKE' (the prognostic TKE producer) and supplies the
# namelist.cvmix &param_tke group (from work_tke_dump). use_sw_pene stays .false. (unforced =>
# sw_3d=0). Default 0 = the PP/KPP gate. (MIX_TKE takes precedence over MIX_KPP.)
MIX_TKE="${MIX_TKE:-0}"
# M6a-2: WHICH_ALE selects the vertical coordinate (default linfs = the reduced-M2 gate).
# WHICH_ALE=zstar keeps the production zstar free surface (Shchepetkin PGF + thickness stretch).
WHICH_ALE="${WHICH_ALE:-linfs}"

source /home/a/a270088/fesom3/env.sh intel >/dev/null 2>&1

rm -rf "$RUN"; mkdir -p "$RUN"
cp "$F2"/work_core/namelist.* "$RUN"/
# M7b: cvmix_TKE needs the &param_tke namelist group (tke_cd=3.75, the namelist-over-codedefault).
[ "$MIX_TKE" = 1 ] && cp "$F2"/work_tke_dump/namelist.cvmix "$RUN"/
ln -sf "$F2"/build/bin/fesom.x "$RUN"/fesom.x      # loads build/lib64/libfesom.so (output 1-rank fix)
printf '0 1 1948\n0 1 1948\n' > "$RUN"/fesom.clock # cold start, year 1948 (unforced: year irrelevant)

python3 - "$RUN" "$NSTEPS" "$FER_GM" "$REDI" "$MIX_KPP" "$WHICH_ALE" "$MIX_TKE" <<'PY'
import re, sys
run, nsteps, fer_gm, redi, mix_kpp, which_ale, mix_tke = sys.argv[1:8]
# namelist.config
p = run + '/namelist.config'; s = open(p).read()
s = re.sub(r"ResultPath\s*=\s*'[^']*'", "ResultPath       = './'", s, count=1)
s = re.sub(r"which_ALE\s*=\s*'zlevel'", f"which_ALE          = '{which_ale}'", s, count=1)
s = re.sub(r"use_ice\s*=\s*\.true\.",   "use_ice                  = .false.", s, count=1)
# use_sw_pene OFF: with use_ice=.false. cal_shortwave_rad (inside oce_fluxes) never runs,
# so sw_3d (forcing_init-allocated) stays UNallocated; the tracer TDMA derefs it under
# use_sw_pene=.true. -> segfault. .false. matches the ported step_oce (no sw_3d consumer).
s = re.sub(r"use_sw_pene\s*=\s*\.true\.","use_sw_pene              = .false.", s, count=1)
s = re.sub(r"yearnew\s*=\s*1958",       "yearnew = 1948", s, count=1)
s = re.sub(r"run_length\s*=\s*\d+",     f"run_length        = {nsteps}", s, count=1)
open(p, 'w').write(s)
# namelist.oce  (reduced-M2 dynamics)
p = run + '/namelist.oce'; s = open(p).read()
# M7b: MIX_TKE=1 -> mix_scheme='cvmix_TKE' (prognostic TKE). Else M5b: keep 'KPP' when
# MIX_KPP=1, otherwise reduce to 'PP'.
if mix_tke == '1':
    s = re.sub(r"mix_scheme\s*=\s*'KPP'", "mix_scheme         = 'cvmix_TKE'", s, count=1)
elif mix_kpp != '1':
    s = re.sub(r"mix_scheme\s*=\s*'KPP'", "mix_scheme         = 'PP'", s, count=1)
# M4c/M4d: keep work_core Fer_GM/Redi=.true. when FER_GM=1 / REDI=1; else reduce them off.
if fer_gm != '1':
    s = re.sub(r"Fer_GM\s*=\s*\.true\.",  "Fer_GM             = .false.", s, count=1)
if redi != '1':
    s = re.sub(r"Redi\s*=\s*\.true\.",    "Redi               = .false.", s, count=1)
open(p, 'w').write(s)
PY

cd "$RUN"
export FESOM_DUMP_FILE="$DUMP"          # built-in per-substep node dumps (multi-step)
export FESOM_DUMP_MAXSTEPS="$NSTEPS"
ulimit -s unlimited
echo "run_lifecycle_core2: $NP rank(s), $NSTEPS steps (unforced) -> ${DUMP}.<rank>"
mpirun --mca pml ob1 --mca btl self,vader --oversubscribe -n "$NP" ./fesom.x > "$RUN/run.log" 2>&1 || true
if ls "${DUMP}".* >/dev/null 2>&1; then
    echo "run_lifecycle_core2: done -> ${DUMP}.00000 ($(stat -c%s "${DUMP}".00000) bytes)"
    grep -E 'FESOM Run|fesom_dump_shim|complete|blowup|NaN|Error|error' "$RUN/run.log" | head -20 || true
else
    echo "run_lifecycle_core2: NO DUMP — see $RUN/run.log"; tail -40 "$RUN/run.log"; exit 1
fi
