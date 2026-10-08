#!/usr/bin/env bash
# TRIANGLE REGRESSION at ROUNDOFF level (hybrid-mesh plan docs/plans/2026-10-08-hybrid-meshes.md,
# from Task 7 on). The mixed-mesh geometry and the unified 4-slot kernels change triangle results at
# roundoff level, so the bitwise identity of Task 6 is replaced by this script: it runs the FROZEN
# binaries of commit 5b2a7f3 (the last bit-identical state, /work/ab0995/a270029/fesom3_hybrid/
# frozen_task6) and the current build side by side on pi and CORE2 for a few steps and compares.
#
# Cases (all from a cold start, PHC winter IC):
#   prod:<mesh>:<np>      production driver fesom_lifecycle_native_mr with the work/job_levante physics
#                         (linfs, JRA55, mEVP off, GM+Redi, TKE, splines, vector-invariant momentum with
#                         rvo_upwind 1, SW penetration, w split), NSTEPS steps (default 4), node dumps
#                         of every substep (FESOM_DUMP_ALL) and the field output every step
#   conserve:<mesh>:<np>  fesom_conserve, zstar + Redi + TKE, 20 steps: both must pass the driver's
#                         conservation gate (|relative heat/salt/volume drift| <= TOL_CONS = 1e-12) and
#                         their final heat/salt totals must agree to TOL_CONS (relative)
#
# Metric (prod), per field: d = max|new - ref| / max|ref|, the max over all nodes (elements), levels
# and records compared. Two levels, because the SSH solver is where "roundoff" stops being roundoff:
#   PRE   step 1 BEFORE the first SSH solve -- density, pressure, bvfreq (PRESSURE_BV), Kv (MIXING),
#         ssh_rhs (SSH_RHS; it carries the whole momentum right-hand side: pressure gradient,
#         momentum advection, viscosity, GM, ...): pure roundoff, d <= TOL_PRE (default 1e-9).
#   ALL   every dumped substep of every step (+ d_eta, hbar, eta_n, hnode(_new), w, T, S) and the
#         last record of every output store both runs wrote (temp, salt, ssh, u, v, w, unod, vnod,
#         Av, bolus_u/v, sst, sss, a_ice, m_ice, m_snow; float32): d <= TOL (default 1e-2).
#         The CG stops at a relative residual soltol = 1e-5, so a roundoff change of the matrix or
#         the rhs moves eta by up to ~1e-5 of its range at the first solve (measured CORE2 step 1:
#         1.9e-5, at a few nodes of one small basin); threshold physics (TKE, ice) then turns that
#         into larger local differences (CORE2, 4 steps: up to ~2e-3). A kernel bug (a wrong weight,
#         a boundary edge counted twice) shows at step 1 at O(1e-3..1) and is caught by PRE.
# The d of every field is printed so the margin is visible.
#
# usage: tools/run_tri_regression.sh [case ...]
#   default cases: prod:pi:1 prod:pi:2 prod:core2:8 conserve:pi:1 conserve:pi:2
#   env: NEW_BIN (default /work/ab0995/a270029/fesom3_hybrid/build/bin), REF_BIN (frozen_task6),
#        ROOT (run directory, default /work/ab0995/a270029/fesom3_hybrid/tri_regression),
#        NSTEPS (prod steps, default 4), TOL_PRE, TOL, TOL_CONS, COMPARE_ONLY=1 (reuse the runs)
# Exit 0 iff every case passes. Login-node sized (np <= 8; ~3 min for the default cases).
set -uo pipefail
HYB=/work/ab0995/a270029/fesom3_hybrid
NEW_BIN=${NEW_BIN:-$HYB/build/bin}
REF_BIN=${REF_BIN:-$HYB/frozen_task6}
ROOT=${ROOT:-$HYB/tri_regression}
NSTEPS=${NSTEPS:-4}
TOL_PRE=${TOL_PRE:-1e-9}
TOL=${TOL:-1e-2}
TOL_CONS=${TOL_CONS:-1e-12}
COMPARE_ONLY=${COMPARE_ONLY:-0}
SRC=$(cd "$(dirname "$0")/.." && pwd)
MPI="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpiexec --mca pml ob1 --mca btl self,vader --oversubscribe"
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
IC=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
cases=("$@")
[ ${#cases[@]} -eq 0 ] && cases=(prod:pi:1 prod:pi:2 prod:core2:8 conserve:pi:1 conserve:pi:2)

source "$SRC/env.sh" intel >/dev/null 2>&1
ulimit -s unlimited

mesh_dir() {
    case $1 in
        pi)    echo /home/a/a270088/port2/fesom2/tests/data/MESHES/pi ;;
        core2) echo /pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2 ;;
        *)     echo /work/ab0995/a270029/meshes/$1 ;;
    esac
}

run_one() {   # run_one <kind> <mesh> <np> <bin_dir> <run_dir>
    local kind=$1 mesh=$2 np=$3 bin=$4 rd=$5
    rm -rf "$rd"; mkdir -p "$rd"
    (
        cd "$rd" || exit 1
        export FESOM3_MESH_DIR=$(mesh_dir $mesh) FESOM3_IC_FILE=$IC
        if [ "$kind" = prod ]; then
            export FESOM3_FORCING_DIR=$POOL FESOM3_FORCING=JRA55 FESOM3_RUNOFF_FILE=$POOL/CORE2_runoff.nc \
                FESOM3_SSS_FILE=$POOL/PHC2_salx.nc FESOM3_START_CLOCK="0 1 1958" FESOM3_WHICH_ALE=linfs \
                FESOM3_WHICHEVP=0 FESOM3_FER_GM=1 FESOM3_REDI=1 FESOM3_MIX_TKE=1 FESOM3_SHEAR_SPLINES=1 \
                FESOM3_N2_SPLINES=1 FESOM3_MOMADV_OPT=1 FESOM3_RVO_UPWIND=1.0 FESOM3_SW_PENE=1 \
                FESOM3_CHL_SWEENEY=1 FESOM3_STEP_DIAG=1 FESOM3_WSPLIT=1 FESOM3_NSTEPS=$NSTEPS \
                FESOM3_OUTPUT=$rd/out FESOM3_OUTPUT_EVERY=1 FESOM3_NAMELIST_IO=$rd/none \
                FESOM_DUMP_FILE=$rd/dn FESOM_DUMP_ALL=1 FESOM_DUMP_MAXSTEPS=$NSTEPS
            case $mesh in core2*) export FESOM3_STEP_PER_DAY=48 ;; esac
            $MPI -np $np "$bin/fesom_lifecycle_native_mr" > log 2>&1
        else
            export FESOM3_NSTEPS=20 FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL=$TOL_CONS \
                FESOM3_REDI=1 FESOM3_MIX_TKE=1
            $MPI -np $np "$bin/fesom_conserve" > log 2>&1
        fi
    )
}

compare_prod() {   # compare_prod <ref_run_dir> <new_run_dir>
    bash -c "source /usr/share/Modules/init/bash; module load python3; hash -r; \
             python3 - '$1' '$2' '$TOL_PRE' '$TOL' '$SRC/tools'" <<'EOF'
import sys, os, glob
import numpy as np
import zarr
dref, dnew, tol_pre, tol, tools = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
sys.path.insert(0, tools)
from dump_diff import parse_tree
PRE_SUBSTEPS = (1, 4, 8)          # PRESSURE_BV, MIXING, SSH_RHS: step 1, before the first SSH solve

def dist(ra, rb):
    return float(np.max(np.abs(rb - ra)) / max(np.max(np.abs(ra)), 1e-300))

# (1) substep node dumps, keyed (step, substep, gid, name); all ranks merged
a = parse_tree(os.path.join(dref, "dn"), True)
b = parse_tree(os.path.join(dnew, "dn"), True)
if set(a) != set(b):
    print(f"  dump record sets differ: ref {len(a)} new {len(b)}"); sys.exit(1)
pre, alld = {}, {}
for nm in sorted({k[3] for k in a}):
    keys = [k for k in a if k[3] == nm]
    alld["dump:" + nm] = dist(np.concatenate([np.asarray(a[k]) for k in keys]),
                              np.concatenate([np.asarray(b[k]) for k in keys]))
    kp = [k for k in keys if k[0] == 1 and k[1] in PRE_SUBSTEPS]
    if kp:
        pre[nm] = dist(np.concatenate([np.asarray(a[k]) for k in kp]),
                       np.concatenate([np.asarray(b[k]) for k in kp]))

# (2) output stores: the last record of every float array both runs wrote; fill values masked
def stores(d):
    out = {}
    for p in glob.glob(os.path.join(d, "**", "*.zarr"), recursive=True):
        try:
            g = zarr.open(p, mode="r")
        except Exception:
            continue
        for k in g.array_keys():
            arr = g[k]
            if arr.dtype.kind != "f" or k in ("lon", "lat", "nz", "nz1", "time"):
                continue
            x = np.asarray(arr[-1] if arr.ndim > 1 else arr[:], dtype=np.float64)
            if arr.fill_value is not None:
                x[x == np.float64(np.asarray(arr.fill_value, dtype=arr.dtype))] = np.nan
            out[k] = x
    return out
sa, sb = stores(os.path.join(dref, "out")), stores(os.path.join(dnew, "out"))
common = sorted(set(sa) & set(sb))
for k in common:
    ra, rb = sa[k], sb[k]
    if ra.shape != rb.shape or not np.array_equal(np.isfinite(ra), np.isfinite(rb)):
        print(f"  output {k}: shapes or masks differ"); sys.exit(1)
    m = np.isfinite(ra)
    alld["out:" + k] = dist(ra[m], rb[m]) if m.any() else 0.0

fail = False
print(f"  PRE (step 1 before the first SSH solve), TOL_PRE {tol_pre:g}:")
for k, d in sorted(pre.items(), key=lambda x: -x[1]):
    print(f"    {k:<16s} {d:.3e}" + ("   ABOVE" if d > tol_pre else ""))
    fail |= d > tol_pre
print(f"  ALL ({len(alld) - len(common)} dump fields x all steps, {len(common)} output arrays), TOL {tol:g}:")
for k, d in sorted(alld.items(), key=lambda x: -x[1]):
    print(f"    {k:<16s} {d:.3e}" + ("   ABOVE" if d > tol else ""))
    fail |= d > tol
if not pre or not common:
    print("  nothing compared (missing dumps or output)"); fail = True
sys.exit(1 if fail else 0)
EOF
}

compare_conserve() {   # compare_conserve <case_dir>
    local ok=0
    for v in ref new; do
        if ! grep -q "CONSERVATION OK" "$1/$v/log"; then echo "  $v: conservation gate FAILED"; ok=1; fi
    done
    [ $ok -ne 0 ] && return 1
    grep -A3 "relative drift" "$1/new/log" | sed "s/^/  new /"
    # "fesom_conserve: step 20  heat= <h>  salt= <s>": fields 5 and 7 (ref), 12 and 14 (new)
    paste -d' ' <(grep "step 20 " "$1/ref/log" | tail -1) <(grep "step 20 " "$1/new/log" | tail -1) | \
        awk -v tol="$TOL_CONS" '{h=($12-$5)/$5; s=($14-$7)/$7; if(h<0)h=-h; if(s<0)s=-s;
            printf "  final totals new vs ref: heat %.3e  salt %.3e (TOL_CONS %g)\n", h, s, tol;
            exit !(h<=tol && s<=tol)}'
}

nfail=0
for c in "${cases[@]}"; do
    kind=${c%%:*}; rest=${c#*:}; mesh=${rest%%:*}; np=${rest##*:}
    cd_=$ROOT/${kind}_${mesh}_np${np}
    echo "== $c"
    if [ "$COMPARE_ONLY" != 1 ]; then
        run_one $kind $mesh $np "$REF_BIN" "$cd_/ref"; rr=$?
        run_one $kind $mesh $np "$NEW_BIN" "$cd_/new"; rn=$?
        if [ $rr -ne 0 ] || [ $rn -ne 0 ]; then
            echo "  run failed (ref rc=$rr, new rc=$rn): see $cd_/{ref,new}/log"; nfail=$((nfail+1)); continue
        fi
    fi
    if [ $kind = prod ]; then compare_prod "$cd_/ref" "$cd_/new"; else compare_conserve "$cd_"; fi \
        && echo "  PASS $c" || { echo "  FAIL $c"; nfail=$((nfail+1)); }
done
echo "run_tri_regression: $((${#cases[@]}-nfail))/${#cases[@]} cases pass" \
     "(TOL_PRE=$TOL_PRE, TOL=$TOL, TOL_CONS=$TOL_CONS)"
[ $nfail -eq 0 ]
