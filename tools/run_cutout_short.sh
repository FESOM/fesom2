#!/usr/bin/env bash
# SHORT RUNS ON THE CUT-OUT MESHES (hybrid-mesh plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 17):
# the first multi-day runs of the mixed-mesh model, on the regional NEMO cut-outs (lon -20..-5,
# lat -2..12, closed boundaries; nemo025q_cut = all quads, nemo025mix_cut = quads + triangles).
#
# Per mesh three cases, each at np 1 and at every np in NPS:
#   prod   the production driver fesom_lifecycle_native_mr with the work/job_levante physics (the
#          prod case of tools/run_tri_regression.sh: linfs, JRA55 from 1958-01-01, EVP sea ice, GM +
#          Redi, TKE, splines, vector-invariant momentum with rvo_upwind 1, SW penetration, w split),
#          NDAYS days. All np restart from ONE checkpoint written by a 0-step np-1 cold start: a cold
#          start is not decomposition-independent (the IC extrapolation into land-masked PHC points,
#          LESSONS L60), a common state is. The SSH is solved to FESOM3_SSH_SOLTOL=$SSH_SOLTOL so that
#          the 1-vs-N difference measures roundoff and not the CG tolerance. Checks:
#            STABLE   every run finishes; max|eta| <= ETA_MAX and max|uv| <= UV_MAX at every step
#                     (FESOM3_STEP_DIAG, finite)
#            1 vs N   every CKPT_STEPS steps the full prognostic state (restart checkpoints, double):
#                     d = max|N - 1| / max|1| per field (tools/cutout_diag.py compare). EARLY tier: all
#                     fields over the first day <= TOL_EARLY; whole run: classes OCEAN / MIX (tke) /
#                     ICE <= TOL_OCEAN / TOL_MIX / TOL_ICE. The whole-run tolerances are loose on
#                     purpose: TKE's thresholds turn roundoff into a local O(1e-9) jump within hours
#                     (the switch fires or not depending on the arithmetic: np, build, dumps on/off),
#                     and from day ~3 the flow against the artificial closed wall at the equator
#                     (near -7E, -1N) amplifies it by ~e every few hours (see prodpp)
#            CHECKERBOARD  ssh, sst, sss of the np-1 run (snapshots every step, double): the
#                     checkerboard energies E_loc / E_glob and the interface-noise ratios R_tri /
#                     R_quad / R_iface (definitions in tools/cutout_diag.py) as time series. REPORTED,
#                     not gated: growth is a design question for the user (plan Technical Details).
#   prodpp prod with PP vertical mixing instead of TKE (FESOM3_MIX_TKE unset), 1 vs N only: the
#          diagnostic that removes the threshold amplifier, as mEVP does for the EVP (Task 15). Its
#          whole run is gated at the tight TOL_PP (OCEAN, MIX and ICE classes)
#   cons   fesom_conserve, unforced, zstar + the production ocean physics the driver supports (TKE,
#          GM + Redi, splines, vector invariant, rvo_upwind 1, w split) with SEA ICE ON (seeded ice,
#          wind on ice, EVP + FCT advection + ice-ocean drag; FESOM3_CONSERVE_ICE), NDAYS days: the
#          driver's gate |relative drift| <= TOL_CONS of heat, salt, volume and the three ice masses.
#          The cut-outs are tropical: the prod runs stay ice-free (ice module on, a_ice = 0), so the
#          converted ice kernels run coupled to the ocean only here. linfs (production) does not
#          conserve tracers (surface w*T term), hence a separate zstar case for the conservation.
#
# Measured 2026-10-08 (5 days, np 1 vs 2 and 8): EARLY <= 1.9e-10 on both meshes; whole run prod
# quad cut-out v 1.0e-4, tke 7.1e-5, T 1.5e-6 (growth from day ~3.9 at the equatorial wall),
# mixed cut-out v 3.8e-9; prodpp <= 1.1e-7 on both; one exploratory np-1 run in which the TKE
# switch fired early reached 7.9e-3 (quad). cons |drift| <= 1.3e-14 (ocean), 8.2e-15 (ice).
# The ICE class of 1 vs N is vacuous here (a_ice = 0 in prod and prodpp). Ice seeded into the
# common checkpoint does not help: the tropical forcing melts it within 2 steps (a_ice 0.48 -> 0 by
# step 3, review fixes 2026-10-09) and the melt event lifts OCEAN 1 vs N to 1e-2. A non-vacuous
# ice compare needs a run that keeps ice (no switch for the ice thermodynamics exists).
#
# usage: tools/run_cutout_short.sh [mesh ...]       (default: nemo025q_cut nemo025mix_cut)
#   env: BIN (default /work/ab0995/a270029/fesom3_hybrid/build/bin), ROOT (run directory, default
#        /work/ab0995/a270029/fesom3_hybrid/cutout_short), NDAYS (default 5), NPS (default "2 8"),
#        CKPT_STEPS (default 12 = 6 h), SSH_SOLTOL (default 1e-14), 1 vs N tolerances (see the
#        measured levels above) TOL_EARLY (default 1e-6), TOL_OCEAN / TOL_MIX (prod whole run,
#        default 5e-2), TOL_ICE (prod whole run, default 1e-6), TOL_PP (prodpp whole run,
#        default 1e-5; prodpp: measured <= 1.1e-7), TOL_CONS (default 1e-12),
#        ETA_MAX (default 3 m), UV_MAX (default 3 m/s), CASES (default "prod prodpp cons"),
#        KEEP=1 keeps the checkpoints (default: removed after the comparison; logs, csv tables and
#        the np-1 output stay), COMPARE_ONLY=1
# Exit 0 iff every check passes. Login-node sized (np <= 8; ~25 min for the defaults).
set -uo pipefail
HYB=/work/ab0995/a270029/fesom3_hybrid
BIN=${BIN:-$HYB/build/bin}
ROOT=${ROOT:-$HYB/cutout_short}
NDAYS=${NDAYS:-5}
NPS=${NPS:-"2 8"}
CKPT_STEPS=${CKPT_STEPS:-12}
SSH_SOLTOL=${SSH_SOLTOL:-1e-14}
TOL_EARLY=${TOL_EARLY:-1e-6}
TOL_OCEAN=${TOL_OCEAN:-5e-2}
TOL_MIX=${TOL_MIX:-5e-2}
TOL_ICE=${TOL_ICE:-1e-6}
TOL_PP=${TOL_PP:-1e-5}
TOL_CONS=${TOL_CONS:-1e-12}
ETA_MAX=${ETA_MAX:-3}
UV_MAX=${UV_MAX:-3}
CASES=${CASES:-"prod prodpp cons"}
KEEP=${KEEP:-0}
COMPARE_ONLY=${COMPARE_ONLY:-0}
SRC=$(cd "$(dirname "$0")/.." && pwd)
# btl_vader_single_copy_mechanism none: Open MPI 4.1.2's knem registration cache misses glibc heap
# trims and then delivers the stale payload of an earlier message (LESSONS L61)
MPI="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpiexec --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe"
source "$SRC/tools/prod_env.sh"          # prod_env / prod_ocean_env (work/job_levante physics)
IC=$PROD_IC
DT=1800
NSTEPS=$((NDAYS * 86400 / DT))
meshes=("$@")
[ ${#meshes[@]} -eq 0 ] && meshes=(nemo025q_cut nemo025mix_cut)

source "$SRC/env.sh" intel >/dev/null 2>&1
ulimit -s unlimited
PY="source /usr/share/Modules/init/bash; module load python3; hash -r; python3"

case_env() {   # case_env <prod|prodpp>: the prod physics, SSH solved tightly; PP mixing for prodpp
    prod_env $SSH_SOLTOL
    [ "$1" != prod ] && unset FESOM3_MIX_TKE
    return 0
}

run_prod() {   # run_prod <mesh> <case_dir> [prodpp]: checkpoint, then np 1 and NPS from it
    local mesh=$1 cd_=$2 kind=${3:-prod} n
    rm -rf "$cd_"; mkdir -p "$cd_/ckpt" "$cd_/w0"
    # output of the np-1 run: ssh/sst/sss snapshots every step, double
    printf "&nml_list\n  io_list = 'ssh   ', 1, 's', 8, 'snap',\n            'sst   ', 1, 's', 8, 'snap',\n            'sss   ', 1, 's', 8, 'snap',\n/\n" > "$cd_/namelist.io"
    (
        cd "$cd_/w0" || exit 1
        case_env $kind
        export FESOM3_MESH_DIR=/work/ab0995/a270029/meshes/$mesh FESOM3_NSTEPS=0 FESOM3_RESTART=$cd_/ckpt \
            FESOM3_RESTART_IN=$cd_/ckpt FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1 \
            FESOM3_NAMELIST_IO=$cd_/none
        $MPI -np 1 "$BIN/fesom_lifecycle_native_mr" > log 2>&1
    ) || { echo "  checkpoint run failed (see $cd_/w0/log)"; return 1; }
    [ -s "$cd_/ckpt/restart.latest" ] || { echo "  no checkpoint written (see $cd_/w0/log)"; return 1; }
    for n in 1 $NPS; do
        mkdir -p "$cd_/np$n" "$cd_/ck$n"
        (
            cd "$cd_/np$n" || exit 1
            case_env $kind
            export FESOM3_MESH_DIR=/work/ab0995/a270029/meshes/$mesh FESOM3_NSTEPS=$NSTEPS \
                FESOM3_RESTART=$cd_/ck$n FESOM3_RESTART_IN=$cd_/ckpt FESOM3_RESTART_UNIT=s \
                FESOM3_RESTART_LENGTH=$CKPT_STEPS FESOM3_NAMELIST_IO=$cd_/none
            [ $n -eq 1 ] && [ $kind = prod ] && export FESOM3_OUTPUT=$cd_/out1 FESOM3_NAMELIST_IO=$cd_/namelist.io
            $MPI -np $n "$BIN/fesom_lifecycle_native_mr" > log 2>&1
        ) || { echo "  np $n run failed (see $cd_/np$n/log)"; return 1; }
        grep -q "RESTART MODE" "$cd_/np$n/log" || { echo "  np $n did not restart (see $cd_/np$n/log)"; return 1; }
    done
}

check_stable() {   # check_stable <log>: max|eta|, max|uv| over all steps from the STEP_DIAG lines
    awk -v em=$ETA_MAX -v um=$UV_MAX -v ns=$NSTEPS '
        /fesom_lifecycle_native_mr: step [0-9]+ / {
            n++; e=$5+0; u=$7+0; a=$9+0
            if ($5 !~ /^[0-9.Ee+-]+$/ || $7 !~ /^[0-9.Ee+-]+$/) bad=1
            if (e>me) {me=e; se=$3}; if (u>mu) {mu=u; su=$3}; if (a>ma) ma=a; le=e; lu=u }
        END { printf "    %d steps; max|eta| %.3f m (step %d), max|uv| %.3f m/s (step %d), last %.3f / %.3f; max a_ice %g\n",
                     n, me, se, mu, su, le, lu, ma
              exit !(n==ns && !bad && me<=em && mu<=um) }' "$1"
}

check_prod() {   # check_prod <mesh> <case_dir> [prodpp]
    local mesh=$1 cd_=$2 kind=${3:-prod} n rc=0 tols="$TOL_OCEAN $TOL_MIX $TOL_ICE"
    [ $kind = prodpp ] && tols="$TOL_PP $TOL_PP $TOL_PP"
    for n in 1 $NPS; do
        echo "  STABLE np $n:"
        check_stable "$cd_/np$n/log" || { echo "    ABOVE the bounds (ETA_MAX $ETA_MAX, UV_MAX $UV_MAX) or not finished"; rc=1; }
    done
    for n in $NPS; do
        echo "  1 vs $n (restart checkpoints every $CKPT_STEPS steps):"
        bash -c "$PY '$SRC/tools/cutout_diag.py' compare '$cd_/ck1' '$cd_/ck$n' $DT $((86400 / DT)) $TOL_EARLY $tols '$cd_/d_1vs$n.csv'" || rc=1
    done
    if [ $kind = prod ]; then
        echo "  CHECKERBOARD (np 1, reported, not gated):"
        bash -c "$PY '$SRC/tools/cutout_diag.py' checker /work/ab0995/a270029/meshes/$mesh '$cd_/out1' ssh sst sss --csv '$cd_/checker.csv'" || rc=1
    fi
    if [ "$KEEP" != 1 ]; then rm -rf "$cd_/ckpt"; for n in 1 $NPS; do rm -rf "$cd_/ck$n"; done; fi
    return $rc
}

run_cons() {   # run_cons <mesh> <case_dir>: fesom_conserve with ice, np 1 and NPS
    local mesh=$1 cd_=$2 n
    rm -rf "$cd_"
    for n in 1 $NPS; do
        mkdir -p "$cd_/np$n"
        (
            cd "$cd_/np$n" || exit 1
            prod_ocean_env
            export FESOM3_MESH_DIR=/work/ab0995/a270029/meshes/$mesh FESOM3_IC_FILE=$IC FESOM3_NSTEPS=$NSTEPS \
                FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_ICE=1 FESOM3_WHICHEVP=0 FESOM3_CONSERVE_TOL=$TOL_CONS
            $MPI -np $n "$BIN/fesom_conserve" > log 2>&1
        )
    done
}

check_cons() {   # check_cons <case_dir>
    local cd_=$1 n rc=0
    for n in 1 $NPS; do
        if grep -q "sea ice ON" "$cd_/np$n/log" && grep -q "CONSERVATION OK" "$cd_/np$n/log"; then
            echo "  np $n: conservation OK (TOL_CONS $TOL_CONS), relative drift over $NSTEPS steps:"
        else
            echo "  np $n: conservation gate FAILED or ice off (see $cd_/np$n/log)"; rc=1
        fi
        grep -A4 "relative drift" "$cd_/np$n/log" | tail -4 | sed "s/^/   /"
        grep "ice step $NSTEPS " "$cd_/np$n/log" | sed "s/.*max|uice|/    max|uice| at the last step/"
    done
    return $rc
}

mkdir -p "$ROOT"
nfail=0
bash -c "$PY '$SRC/tools/cutout_diag.py' selftest" > "$ROOT/selftest.log" 2>&1 || {
    cat "$ROOT/selftest.log"; echo "run_cutout_short: cutout_diag selftest FAILED"; exit 1; }
for mesh in "${meshes[@]}"; do
    for kind in $CASES; do
        cd_=$ROOT/${kind}_${mesh}
        echo "== $kind:$mesh ($NDAYS days = $NSTEPS steps, np 1 vs $NPS)"
        if [ "$COMPARE_ONLY" != 1 ]; then
            case $kind in prod*) run_prod $mesh "$cd_" $kind ;; *) run_cons $mesh "$cd_" ;; esac || { echo "  FAIL $kind:$mesh (run)"; nfail=$((nfail+1)); continue; }
        fi
        case $kind in prod*) check_prod $mesh "$cd_" $kind ;; *) check_cons "$cd_" ;; esac \
            && echo "  PASS $kind:$mesh" || { echo "  FAIL $kind:$mesh"; nfail=$((nfail+1)); }
    done
done
ncase=$((${#meshes[@]} * $(echo $CASES | wc -w)))
echo "run_cutout_short: $((ncase-nfail))/$ncase cases pass (NDAYS=$NDAYS, SSH_SOLTOL=$SSH_SOLTOL," \
     "TOL_EARLY=$TOL_EARLY TOL_OCEAN=$TOL_OCEAN TOL_MIX=$TOL_MIX TOL_ICE=$TOL_ICE TOL_PP=$TOL_PP" \
     "TOL_CONS=$TOL_CONS)"
[ $nfail -eq 0 ]
