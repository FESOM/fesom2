#!/usr/bin/env bash
# TRIANGLE REGRESSION at ROUNDOFF level (hybrid-mesh plan docs/plans/completed/2026-10-08-hybrid-meshes.md,
# from Task 7 on). The mixed-mesh geometry and the unified 4-slot kernels change triangle results at
# roundoff level, so the bitwise identity of Task 6 is replaced by this script: it runs the FROZEN
# reference binaries and the current build side by side on pi and CORE2 for a few steps and compares.
#
# Reference: /work/ab0995/a270029/fesom3_hybrid/frozen_task6_soltol = commit 5b2a7f3 (the last
# bit-identical state) + ONLY the FESOM3_SSH_SOLTOL/FESOM3_SSH_MAXITER env knob of the SSH CG
# (src/oce/oce_ssh_solve.F90; unset => soltol 1e-5 as before; with the knob unset its results are
# byte-identical to frozen_task6). See the README there.
#
# Cases (all from a cold start, PHC winter IC):
#   prod:<mesh>:<np>      production driver fesom_lifecycle_native_mr with the work/job_levante physics
#                         (linfs, JRA55, EVP sea ice, GM+Redi, TKE, splines, vector-invariant momentum
#                         with rvo_upwind 1, SW penetration, w split), NSTEPS steps (default 4), node
#                         dumps of every substep (FESOM_DUMP_ALL, double) and the field output every
#                         step (float32): REF binary vs NEW build at the same np.
#   decomp:<mesh>:<np>    DECOMPOSITION check of the NEW build: np 1 vs np <np>, both restarted from
#                         ONE checkpoint written by a 0-step np-1 cold start. A cold start cannot be
#                         compared across np: the IC extrapolation into land-masked PHC points
#                         (oce_initial_state extrap_nod3D, FESOM2's Gauss-Seidel sweep with a per-rank
#                         halo snapshot) is decomposition-dependent BY DESIGN (pi: 34 nodes / 513
#                         levels differ, density up to 4.7e-3 relative). From a common state np1 == npN
#                         to roundoff. Use fesom_ini meshes (pi_w, core2_w: dist_N with format.out).
#   conserve:<mesh>:<np>  fesom_conserve, zstar + Redi + TKE, 20 steps: both must pass the driver's
#                         conservation gate (|relative heat/salt/volume drift| <= TOL_CONS = 1e-12) and
#                         their final heat/salt totals must agree to TOL_CONS (relative)
#
# SSH SOLVER. Every prod/decomp run solves the SSH system to FESOM3_SSH_SOLTOL=$SSH_SOLTOL (default
# 1e-14; production uses 1e-5). With 1e-5 the CG iterate differences (~1e-5 of eta at the first
# solve) swamp the roundoff and a gate can only be ~1e-2; with 1e-14 the solve itself stays at
# roundoff (pi step 1: d_eta 4e-15) at ~3x the CG iterations.
#
# Metric, per field: d = max|a - b| / max|a| over all nodes/elements, levels and records compared.
# Fields are checked in tiers and classes, each with its own tolerance (pi / core2 meshes; meshes
# whose name starts with "pi" use the pi column, every other mesh the core2 column; env overrides
# N_EARLY TOL_PRE TOL_EARLY TOL_OCEAN TOL_ICE TOL_MIX TOL_ICE_EARLY apply to all meshes):
#   tier/class  fields                                                    pi      core2
#   PRE         step 1 before the first SSH solve: density, pressure,     1e-9    1e-9
#               bvfreq, Kv, ssh_rhs (the whole momentum rhs)
#   EARLY       every double dump field except Kv, every substep of       1e-10   1e-7
#               steps 1..N_EARLY (pi 3, core2 1)
#   OCEAN       double dumps of all steps except Kv (T, S, eta_n, hbar,   1e-6    2e-4
#               d_eta, ssh_rhs, w, density, pressure, bvfreq, hnode) and
#               all output records of temp salt ssh sst sss u v w unod
#               vnod bolus_u bolus_v (float32: + 2^-21, 4 float32 ulps)
#   ICE         output records of a_ice m_ice m_snow (float32 + 2^-21)    1e-5    2e-3
#   MIX         Kv (all steps), Av (output) -- TKE thresholds             1e-6    1e-2
#   ICE_EARLY   the ICE output records of steps 1..N_EARLY (float32       1e-6    1e-5
#               + 2^-21): an ice-only kernel error shows at step 1, below
#               the ICE / MIX levels that the later amplification needs
# Measured 2026-10-08 (Task-7 build vs reference, 4 steps; decomp = np1 vs np2/np8 on pi_w):
#   PRE   pi <= 6.3e-15, CORE2 1.3e-11 (ssh_rhs); decomp exactly 0
#   EARLY pi <= 9.2e-13 (w, steps 1-3), decomp <= 4.3e-13; CORE2 step 1 <= 2.1e-9 (T at
#         ice-covered nodes: the sea ice, stiff at rest -- Delta_min-regularised EVP viscosity --
#         amplifies the geometry roundoff into ice velocity/thickness and the ice-ocean fluxes)
#   OCEAN pi <= 3.3e-8 dumps (d_eta), 1.0e-7 float32 u/v -- the jump comes at step 4, the first
#         JRA55 record refresh, which moves the ice (same with the summer IC and in np1 vs npN on
#         the original pi); decomp <= 2.5e-11 dumps, 6.4e-10 output; CORE2 <= 4.4e-6 dumps (T),
#         2.2e-5 float32 u/v
#   ICE   pi <= 3.9e-7, CORE2 1.2e-4;  MIX pi <= 2.2e-9, decomp 1.9e-11, CORE2 2.4e-4 (TKE)
#   ICE_EARLY (2026-10-09, review fixes): pi steps 1-3 and decomp exactly 0 in float32, CORE2
#         step 1 <= 1.4e-6 (m_ice)
# ELEMENT COORDINATES (Task 14). The element coordinate of the output (and the point at which the
# element vectors u/v and bolus_u/bolus_v are r2g-rotated) is now the area-weighted, cyclic-wrapped
# centroid (mod_mesh_areas elem_center); the reference took the plain mean of the 3 rotated vertices
# without the wrap, which lies on the wrong side of the globe for the elements across the rotated
# dateline (pi 77, CORE2 566 elements; every other element agrees to < 1e-10 deg). For a vector pair
# the comparison therefore (a) finds those elements as the ones whose embedded lon/lat moved by more
# than 1e-6 deg -- every moved element must have moved by more than 1 deg (the wrap fix, not a
# drift) and at most 5 % may move, else FAIL; (b) compares the components on all other elements
# and the rotation-invariant speed hypot(u, v) on the moved ones (normalised by max|speed|), both
# in the field's usual tier. Scalars are not rotated and are compared everywhere. Measured after
# Task 14: pi u/v 1.1e-7 elsewhere, speed 2.7e-8 at the moved elements; CORE2 2.5e-5 / 6.6e-6.
# SEA-ICE DYNAMICS (Task 15). The EVP/mEVP kernels now use the 4-slot formulas (kite-weighted metric
# means and node shares), so roundoff enters INSIDE each of the 120 EVP subcycles, not only through
# the geometry. The standard EVP amplifies it: CORE2 EARLY T 2.2e-8 at step 1 (was 2.1e-9), so the
# CORE2 EARLY tolerance went 1e-8 -> 1e-7 (same ~5x margin); the other tiers stayed at their levels
# (OCEAN 3.1e-5, ICE 8.9e-5, MIX 7.9e-4). Evidence that this is amplification, not a kernel error:
# the same case with mEVP (FESOM3_WHICHEVP=1 in both runs; mEVP runs the same converted ocean2ice,
# SSH-tilt, ice-strength and metric code) gives EARLY T 1.1e-15 and ICE 6.9e-6; triangles equal the
# former EVP formulas to <= 1.1e-16 element by element (test_ice_dyn_mixed on core2_w np8).
#
# Why these levels: below them the differences are no longer roundoff but its amplification by
# stiff / discontinuous physics (sea ice above all, TKE), so they are the tightest defensible ones
# (margins ~10-100x). For scale: the SAME binary with the production SSH tolerance 1e-5 vs 1e-14
# differs by 3.5e-4 in d_eta at step 1 on pi (EARLY fails by 6 orders); a kernel error (a wrong
# weight, a boundary edge counted twice) shows at O(1e-3..1). The d of every field is printed.
#
# usage: tools/run_tri_regression.sh [case ...]
#   default cases: prod:pi:1 prod:pi:2 prod:core2:8 decomp:pi_w:2 decomp:pi_w:8 conserve:pi:1
#                  conserve:pi:2
#   env: NEW_BIN (default /work/ab0995/a270029/fesom3_hybrid/build/bin), REF_BIN (frozen_task6_soltol),
#        ROOT (run directory, default /work/ab0995/a270029/fesom3_hybrid/tri_regression),
#        NSTEPS (prod/decomp steps, default 4), SSH_SOLTOL (default 1e-14), TOL_* (see above),
#        TOL_CONS, COMPARE_ONLY=1 (reuse the runs)
# Exit 0 iff every case passes. Login-node sized (np <= 8; ~5 min for the default cases).
set -uo pipefail
HYB=/work/ab0995/a270029/fesom3_hybrid
NEW_BIN=${NEW_BIN:-$HYB/build/bin}
REF_BIN=${REF_BIN:-$HYB/frozen_task6_soltol}
ROOT=${ROOT:-$HYB/tri_regression}
NSTEPS=${NSTEPS:-4}
SSH_SOLTOL=${SSH_SOLTOL:-1e-14}
TOL_CONS=${TOL_CONS:-1e-12}
COMPARE_ONLY=${COMPARE_ONLY:-0}
SRC=$(cd "$(dirname "$0")/.." && pwd)
# btl_vader_single_copy_mechanism none: Open MPI 4.1.2's knem registration cache misses glibc heap
# trims and then delivers the stale payload of an earlier message (LESSONS L61,
# tools/run_mpi_knem_stale_repro.sh)
MPI="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpiexec --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe"
source "$SRC/tools/prod_env.sh"          # prod_env: the production physics (work/job_levante)
IC=$PROD_IC
cases=("$@")
[ ${#cases[@]} -eq 0 ] && cases=(prod:pi:1 prod:pi:2 prod:core2:8 decomp:pi_w:2 decomp:pi_w:8 \
                                 conserve:pi:1 conserve:pi:2)

source "$SRC/env.sh" intel >/dev/null 2>&1
ulimit -s unlimited

mesh_dir() {
    case $1 in
        pi)    echo /home/a/a270088/port2/fesom2/tests/data/MESHES/pi ;;
        core2) echo /pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2 ;;
        *)     echo /work/ab0995/a270029/meshes/$1 ;;
    esac
}

tolerances() {   # tolerances <mesh> -> "PRE N_EARLY EARLY OCEAN ICE MIX ICE_EARLY"
    local t
    case $1 in pi*) t=(1e-9 3 1e-10 1e-6 1e-5 1e-6 1e-6) ;; *) t=(1e-9 1 1e-7 2e-4 2e-3 1e-2 1e-5) ;; esac
    echo "${TOL_PRE:-${t[0]}} ${N_EARLY:-${t[1]}} ${TOL_EARLY:-${t[2]}} ${TOL_OCEAN:-${t[3]}}" \
         "${TOL_ICE:-${t[4]}} ${TOL_MIX:-${t[5]}} ${TOL_ICE_EARLY:-${t[6]}}"
}

run_one() {   # run_one <kind> <mesh> <np> <bin_dir> <run_dir> [nsteps] [extra env ...]
    local kind=$1 mesh=$2 np=$3 bin=$4 rd=$5 ns=${6:-$NSTEPS}
    shift 6 2>/dev/null || shift $#
    rm -rf "$rd"; mkdir -p "$rd"
    (
        cd "$rd" || exit 1
        export FESOM3_MESH_DIR=$(mesh_dir $mesh) FESOM3_IC_FILE=$IC
        for v in "$@"; do export "$v"; done
        if [ "$kind" = prod ]; then
            prod_env $SSH_SOLTOL
            export FESOM3_NSTEPS=$ns FESOM3_OUTPUT=$rd/out FESOM3_OUTPUT_EVERY=1 \
                FESOM3_NAMELIST_IO=$rd/none FESOM_DUMP_FILE=$rd/dn FESOM_DUMP_ALL=1 \
                FESOM_DUMP_MAXSTEPS=$ns
            case $mesh in core2*) export FESOM3_STEP_PER_DAY=48 ;; esac
            $MPI -np $np "$bin/fesom_lifecycle_native_mr" > log 2>&1
        else
            export FESOM3_NSTEPS=20 FESOM3_WHICH_ALE=zstar FESOM3_CONSERVE_TOL=$TOL_CONS \
                FESOM3_REDI=1 FESOM3_MIX_TKE=1
            $MPI -np $np "$bin/fesom_conserve" > log 2>&1
        fi
    )
}

compare_prod() {   # compare_prod <run_dir_a> <run_dir_b> <mesh>
    bash -c "source /usr/share/Modules/init/bash; module load python3; hash -r; \
             python3 - '$1' '$2' $(tolerances $3) '$SRC/tools'" <<'EOF'
import sys, os, glob
import numpy as np
import zarr
da, db = sys.argv[1], sys.argv[2]
n_early = int(sys.argv[4])
tol = dict(zip(("PRE", "EARLY", "OCEAN", "ICE", "MIX", "ICE_EARLY"), map(float, [sys.argv[3]] + sys.argv[5:10])))
sys.path.insert(0, sys.argv[10])
from dump_diff import parse_tree
PRE_SUBSTEPS = (1, 4, 8)          # PRESSURE_BV, MIXING, SSH_RHS: step 1, before the first SSH solve
F32 = 2.0 ** -21                  # 4 float32 ulps of max|a|: rounding of the output to float32
MIX_DUMP, MIX_OUT, ICE_OUT = {"Kv"}, {"Av"}, {"a_ice", "m_ice", "m_snow"}
VEC_PAIR = {"u": "v", "v": "u", "bolus_u": "bolus_v", "bolus_v": "bolus_u"}   # element vectors

def dist(ra, rb):
    return float(np.max(np.abs(rb - ra)) / max(np.max(np.abs(ra)), 1e-300))

# (1) substep node dumps (double), keyed (step, substep, gid, name); all ranks merged
a = parse_tree(os.path.join(da, "dn"), True)
b = parse_tree(os.path.join(db, "dn"), True)
if set(a) != set(b):
    print(f"  dump record sets differ: {len(a)} vs {len(b)}"); sys.exit(1)
res = {t: {} for t in tol}            # tier/class -> {field: d}
lim = {t: {} for t in tol}            # per-field limit (float32 allowance)
def put(tier, key, d, extra=0.0):
    res[tier][key] = d; lim[tier][key] = tol[tier] + extra
cat = lambda r, ks: np.concatenate([np.asarray(r[k]) for k in ks])
for nm in sorted({k[3] for k in a}):
    keys = [k for k in a if k[3] == nm]
    put("MIX" if nm in MIX_DUMP else "OCEAN", "dump:" + nm, dist(cat(a, keys), cat(b, keys)))
    ke = [k for k in keys if k[0] <= n_early]
    if ke and nm not in MIX_DUMP:
        put("EARLY", nm, dist(cat(a, ke), cat(b, ke)))
    kp = [k for k in keys if k[0] == 1 and k[1] in PRE_SUBSTEPS]
    if kp:
        put("PRE", nm, dist(cat(a, kp), cat(b, kp)))

# (2) output stores (float32): every record of every float array both runs wrote; fill values masked
def stores(d):
    out, elon = {}, {}
    for p in glob.glob(os.path.join(d, "**", "*.zarr"), recursive=True):
        try:
            g = zarr.open(p, mode="r")
        except Exception as e:
            print(f"  output store {p}: cannot open ({e})"); sys.exit(1)
        for k in g.array_keys():
            arr = g[k]
            if arr.dtype.kind != "f" or k in ("lon", "lat", "nz", "nz1", "time"):
                continue
            x = np.asarray(arr[:], dtype=np.float64)
            if arr.fill_value is not None:
                x[x == np.float64(np.asarray(arr.fill_value, dtype=arr.dtype))] = np.nan
            out[k] = x
            if "elem" in arr.attrs.get("_ARRAY_DIMENSIONS", []) and "lon" in g and "lat" in g:
                elon[k] = (np.asarray(g["lon"][:], dtype=np.float64), np.asarray(g["lat"][:], dtype=np.float64))
    return out, elon
(sa, ca), (sb, cb) = stores(os.path.join(da, "out")), stores(os.path.join(db, "out"))
if set(sa) != set(sb):
    print(f"  output variables differ: only in {da}: {sorted(set(sa) - set(sb))},"
          f" only in {db}: {sorted(set(sb) - set(sa))}"); sys.exit(1)
common = sorted(sa)
for k in common:
    ra, rb = sa[k], sb[k]
    if ra.shape != rb.shape or not np.array_equal(np.isfinite(ra), np.isfinite(rb)):
        print(f"  output {k}: shapes or masks differ"); sys.exit(1)
    m = np.isfinite(ra)
    mv = None
    if k in VEC_PAIR and VEC_PAIR[k] in common and k in ca and k in cb:
        dlon = np.abs((cb[k][0] - ca[k][0] + 180.0) % 360.0 - 180.0)
        dlat = np.abs(cb[k][1] - ca[k][1])
        mv = (dlon > 1e-6) | (dlat > 1e-6)              # element coordinate moved (Task 14)
        if not mv.any():
            mv = None
        elif (mv & (dlon <= 1.0) & (dlat <= 1.0)).any() or mv.sum() > 0.05 * mv.size:
            print(f"  output {k}: {int(mv.sum())} element coordinates moved, not all by > 1 deg"
                  " or more than 5 % (expected: only the elements across the rotated dateline)")
            sys.exit(1)
    if mv is None:
        d = dist(ra[m], rb[m]) if m.any() else 0.0
    else:
        keep = m & ~mv
        pa, pb = sa[VEC_PAIR[k]], sb[VEC_PAIR[k]]
        spa, spb = np.hypot(ra, pa), np.hypot(rb, pb)
        mm = m & mv
        d_comp = dist(ra[keep], rb[keep]) if keep.any() else 0.0
        d_speed = float(np.max(np.abs(spb[mm] - spa[mm])) / max(np.max(spa[m]), 1e-300)) if mm.any() else 0.0
        d = max(d_comp, d_speed)
        print(f"  output {k}: {int(mv.sum())} element coordinates moved (rotated-dateline centroids):"
              f" components elsewhere {d_comp:.3e}, speed there {d_speed:.3e}")
    if k in ICE_OUT:
        put("ICE", "out:" + k, d, F32)
        e = m.copy(); e[n_early:] = False          # the records of steps 1..N_EARLY
        put("ICE_EARLY", "out:" + k, dist(ra[e], rb[e]) if e.any() else 0.0, F32)
    elif k in MIX_OUT: put("MIX", "out:" + k, d, F32)
    else:              put("OCEAN", "out:" + k, d, F32)

fail = False
title = {"PRE": "step 1 before the first SSH solve",
         "EARLY": f"all substeps of steps 1..{n_early}, Kv excluded",
         "OCEAN": "all steps / records", "ICE": "sea ice, all records",
         "MIX": "threshold-driven mixing, all steps / records",
         "ICE_EARLY": f"sea ice, the records of steps 1..{n_early}"}
for t in ("PRE", "EARLY", "OCEAN", "ICE", "MIX", "ICE_EARLY"):
    print(f"  {t} ({title[t]}), TOL {tol[t]:g}" + (" (+2^-21 for float32 output)" if t in ("OCEAN", "ICE", "MIX", "ICE_EARLY") else "") + ":")
    for k, d in sorted(res[t].items(), key=lambda x: -x[1]):
        bad = d > lim[t][k]
        print(f"    {k:<16s} {d:.3e}" + ("   ABOVE" if bad else ""))
        fail |= bad
if not res["PRE"] or not common:
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

run_decomp() {   # run_decomp <mesh> <np> <case_dir>: np1 and np<np> from one np-1 checkpoint
    local mesh=$1 np=$2 cd_=$3 r1 rn rc
    rm -rf "$cd_/ckpt" "$cd_/out1" "$cd_/out$np"   # a stale checkpoint folder blocks the write
    mkdir -p "$cd_/ckpt" "$cd_/out1" "$cd_/out$np"
    run_one prod $mesh 1 "$NEW_BIN" "$cd_/w0" 0 FESOM3_RESTART=$cd_/ckpt FESOM3_RESTART_IN=$cd_/ckpt \
        FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1 || return 1
    [ -s "$cd_/ckpt/restart.latest" ] || { echo "  no checkpoint written (see $cd_/w0/log)"; return 1; }
    for n in 1 $np; do
        run_one prod $mesh $n "$NEW_BIN" "$cd_/np$n" $NSTEPS FESOM3_RESTART=$cd_/out$n \
            FESOM3_RESTART_IN=$cd_/ckpt FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1 || return 1
        grep -q "RESTART MODE" "$cd_/np$n/log" || { echo "  np $n did not restart (see $cd_/np$n/log)"; return 1; }
    done
}

nfail=0
for c in "${cases[@]}"; do
    kind=${c%%:*}; rest=${c#*:}; mesh=${rest%%:*}; np=${rest##*:}
    cd_=$ROOT/${kind}_${mesh}_np${np}
    echo "== $c"
    if [ "$COMPARE_ONLY" != 1 ]; then
        if [ $kind = decomp ]; then
            run_decomp $mesh $np "$cd_" || { echo "  FAIL $c (run)"; nfail=$((nfail+1)); continue; }
        else
            run_one $kind $mesh $np "$REF_BIN" "$cd_/ref"; rr=$?
            run_one $kind $mesh $np "$NEW_BIN" "$cd_/new"; rn=$?
            if [ $rr -ne 0 ] || [ $rn -ne 0 ]; then
                echo "  run failed (ref rc=$rr, new rc=$rn): see $cd_/{ref,new}/log"; nfail=$((nfail+1)); continue
            fi
        fi
    fi
    case $kind in
        prod)     compare_prod "$cd_/ref" "$cd_/new" $mesh ;;
        decomp)   compare_prod "$cd_/np1" "$cd_/np$np" $mesh ;;
        *)        compare_conserve "$cd_" ;;
    esac && echo "  PASS $c" || { echo "  FAIL $c"; nfail=$((nfail+1)); }
done
echo "run_tri_regression: $((${#cases[@]}-nfail))/${#cases[@]} cases pass" \
     "(SSH_SOLTOL=$SSH_SOLTOL, TOL_CONS=$TOL_CONS; per-mesh TOL table in the header)"
[ $nfail -eq 0 ]
