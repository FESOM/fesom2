#!/usr/bin/env bash
# Restart Stage 6 (Task 6.1) PRIMARY self-consistency byte-gate — the HEADLINE gate that proves the
# whole milestone: a SPLIT run (K steps -> checkpoint -> FRESH process -> restart_read -> N-K steps)
# reaches a final prognostic state BYTE-IDENTICAL (max|Δ|=0) to a STRAIGHT-THROUGH N-step run, sea
# ice + EVP sigma included, at np=1 AND np=2. It is also where Task 5.2's first-resumed-step AB guard
# is truly tested (a wrong guard => uv / uv_rhsAB diverge in comparison #1).
#
# TWO complementary comparisons, so the FULL state incl. ice is covered:
#   #1 CHECKPOINT compare (comprehensive — covers ice+sigma+velocity+AB): straight writes an
#      end-of-run checkpoint C_straight; split writes C_split at the same physical step. zarr_diff.py
#      --output-cmp compares ALL stores (eta_n, u/v, uv_rhsAB, w/w_e/w_i, temp/salt values/AB/M1,
#      ice area/hice/hsnow/uice/vice, EVP sigma11/12/22, + coords) => max|Δ|=0. This is the primary
#      full-state proof AND exercises Task 5.2 (the restored uv_rhsAB must be AB2-blended, not Euler).
#   #2 FESOM_DUMP_ALL compare (independent LIVE-state, catches any writer blind spot): both runs dump
#      the live node dyn/tracer state at their FINAL step; dump_diff.py --glob --ignore-step => max|Δ|=0.
#      This reads the live arrays DIRECTLY (not via the restart writer). NOTE: FESOM_DUMP_ALL (the
#      mod_step_oce dump set) covers NODE dyn/tracer fields only (density/pressure/bvfreq/Kv/ssh_rhs/
#      d_eta/hbar/eta_n/hnode_new/w/T/S/hnode) — it does NOT cover ice, sigma, or element velocity;
#      those are covered by comparison #1. The dump step labels differ (straight step N vs split
#      seg-2 step N-K) so --ignore-step aligns the single-step windows on (substep,gid,name).
#
# Splits run: a MID-RUN split (default 1948-01-01 cold start; K mid-run) AND a BOUNDARY split whose
# checkpoint lands ON the Jan->Feb month boundary (FESOM3_START_CLOCK near 86400 on day 31), so the
# forcing-from-clock resume + roll_monthly_clim read-ahead is actually exercised on resume (a real
# month crossing, asserted via the "slice 2" log). Same-np throughout (bit-identical resume is a
# same-np guarantee — cross-np is the Task 6.2 restore round-trip).
#
#   tools/run_restart_gate_core2.sh [run_dir] [NP...]     # default NP list: "1 2"
set -uo pipefail
F3="${F3:-/home/a/a270088/fesom3/.claude/worktrees/restart-checkpoint}"
BUILD="${BUILD:-$F3/build_intel_dp}"
PY="${PY:-/work/ab0995/a270088/mambaforge/bin/python3}"
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NCAR=/home/a/a270088/port2/fesom2/test/input/global
RUN="${1:-/scratch/a/a270088/restart_gate_core2}"
shift || true
NPS="${*:-1 2}"
N="${N:-4}"            # straight-through total steps
K="${K:-2}"            # split point (seg-1 steps; seg-2 runs N-K) — mid-run when 0<K<N
MID_START="0 1 1948"        # cold mid-run start (Jan 1; K lands mid-day -> no boundary)
BND_START="82800 31 1948"   # 23:00 Jan 31; seg-1 (K=2) ends ON the 86400 Jan->Feb boundary

source "$F3/env.sh" intel >/dev/null 2>&1
BIN="$BUILD/bin/fesom_lifecycle_native_mr"
[ -x "$BIN" ] || { echo "missing $BIN — build first (cmake --build $BUILD --target fesom_lifecycle_native_mr)"; exit 2; }
rm -rf "$RUN"; mkdir -p "$RUN"

export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE"
export FESOM3_FORCING_DIR="$NCAR"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc"
export FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"
export FESOM3_WHICHEVP=0          # EVP (sigma carries across steps; the hard case for F-C)
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe)
GREP='nod2D|INITIALISATION|RESTART MODE|RESTART RUN|state restored|registered|checkpoint written|final checkpoint|roll_monthly|done'
rc=0

# run one lifecycle segment. UNIT=off => exactly one (end-of-run) checkpoint at the segment's last step.
#   $1 NP  $2 restart_dir  $3 nsteps  $4 start_clock  $5 restart_in(""=use restart_dir, fresh=cold)
#   $6 dump_prefix(""=no dump)  $7 dump_step  $8 logfile
run_seg() {
    local np="$1" rdir="$2" nst="$3" sclk="$4" rin="${5:-}" dpfx="${6:-}" dstep="${7:-}" log="$8"
    local -a e=( -u FESOM_DUMP_ALL -u FESOM_DUMP_FILE -u FESOM_DUMP_MINSTEP -u FESOM_DUMP_MAXSTEPS
                 FESOM3_RESTART="$rdir" FESOM3_RESTART_IN="${rin:-$rdir}"
                 FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1
                 FESOM3_NSTEPS="$nst" FESOM3_START_CLOCK="$sclk" )
    if [ -n "$dpfx" ]; then
        mkdir -p "$(dirname "$dpfx")"
        e+=( FESOM_DUMP_ALL=1 FESOM_DUMP_FILE="$dpfx" FESOM_DUMP_MINSTEP="$dstep" FESOM_DUMP_MAXSTEPS="$dstep" )
    fi
    env "${e[@]}" "${MPIRUN[@]}" -n "$np" "$BIN" >"$log" 2>&1
    local s=$?
    grep -E "$GREP" "$log" | sed 's/^/      /'
    return $s
}

# straight + split + the two comparisons for one (label, start_clock) case at one np.
do_case() {
    local np="$1" label="$2" sclk="$3" want_boundary="$4"
    # Strict max|Δ|=0 at every np. A border element is owned by every rank that owns one of its
    # vertices (562 such elements at CORE2 dist_2); update_vel syncs uv to the designated owner's
    # copy (FESOM3_SHARED_SYNC=1, the default; docs/plans/completed/20261012-shared-element-sync.md), which is
    # the copy the restart writer keeps, so a resumed np>1 run continues from exactly the state of
    # the straight run. Only with FESOM3_SHARED_SYNC=0 (the unsynced copies, A/B only) is the old
    # ~1-ULP redundant-element floor admitted: the straight run keeps every owner's own copy while
    # the resumed one starts all owners from the written copy. 1e-12 is far below any real
    # regression (the t_skin/d_eta seeds were ~1e-4).
    local floor=0; [ "$np" -gt 1 ] && [ "${FESOM3_SHARED_SYNC:-1}" = 0 ] && floor=1e-12
    local base="$RUN/np$np/$label"
    local strt="$base/straight" splt="$base/split"
    rm -rf "$base"; mkdir -p "$strt" "$splt"
    local rem=$(( N - K ))
    echo "==================================================================================="
    echo "### np=$np  $label split   (N=$N, K=$K, seg-2=$rem steps; start='$sclk') ###"
    echo "==================================================================================="

    echo "--- STRAIGHT: $N steps (cold), end-of-run checkpoint + FESOM_DUMP_ALL@step$N ---"
    run_seg "$np" "$strt" "$N" "$sclk" "$strt" "$strt/dump/node" "$N" "$strt.log"
    local s_strt=$?

    echo "--- SPLIT seg-1: $K steps (cold) -> checkpoint@step$K (NO dump) ---"
    run_seg "$np" "$splt" "$K" "$sclk" "$splt" "" "" "$splt.seg1.log"
    local s_s1=$?
    echo "--- SPLIT seg-2: $rem steps (FRESH process, resumes) + FESOM_DUMP_ALL@step$rem ---"
    run_seg "$np" "$splt" "$rem" "$sclk" "$splt" "$splt/dump/node" "$rem" "$splt.seg2.log"
    local s_s2=$?

    # run sanity: all three exited 0; seg-2 must have actually resumed
    for tag in "straight:$s_strt" "seg1:$s_s1" "seg2:$s_s2"; do
        [ "${tag#*:}" = 0 ] || { echo "FAIL  ${tag%%:*} exit ${tag#*:}"; rc=1; }
    done
    grep -q 'THIS IS A RESTART RUN' "$splt.seg2.log" \
        && echo "PASS  seg-2 resumed (r_restart=.true.)" || { echo "FAIL  seg-2 did not resume"; rc=1; }

    # boundary case: the straight run must have actually crossed the Jan->Feb boundary (read-ahead
    # slice 2). A mid-run K never hits it; this proves the resume exercised roll_monthly_clim.
    if [ "$want_boundary" = yes ]; then
        if grep -q 'slice 2 ' "$strt.log" && grep -q 'slice 2 ' "$splt.seg2.log"; then
            echo "PASS  boundary crossed: roll_monthly_clim read-ahead fired (Feb slice 2) in straight AND split seg-2"
        else
            echo "FAIL  boundary NOT crossed (no 'slice 2' read-ahead) — start clock mistimed"; rc=1
        fi
    fi

    # ---- comparison #1: end-of-run CHECKPOINT, ALL stores (incl. ice+sigma+velocity+AB) ----
    local sf pf
    sf=$(cat "$strt/restart.latest" 2>/dev/null); pf=$(cat "$splt/restart.latest" 2>/dev/null)
    echo "--- CMP #1 checkpoint  straight=$sf  split=$pf ---"
    if [ -z "$sf" ] || [ "$sf" != "$pf" ]; then
        echo "FAIL  checkpoint folder mismatch ('$sf' vs '$pf')"; rc=1
    else
        echo "    CMD: $PY tools/zarr_diff.py --output-cmp $strt/$sf $splt/$pf --rel-floor $floor"
        "$PY" "$F3/tools/zarr_diff.py" --output-cmp "$strt/$sf" "$splt/$pf" --rel-floor "$floor"
        if [ $? = 0 ]; then echo "PASS  #1 checkpoint $([ "$floor" = 0 ] && echo 'max|Δ|=0' || echo '<=1-ULP (np>1 element-ownership floor)') (full state incl. ice+sigma+velocity+AB)"
        else echo "FAIL  #1 checkpoint diverged (see per-store max|Δ| above)"; rc=1; fi
    fi

    # ---- comparison #2: live FESOM_DUMP_ALL at the final step (independent of the restart writer) ----
    echo "--- CMP #2 live FESOM_DUMP_ALL (node dyn/tracer), straight step$N vs split seg-2 step$rem ---"
    echo "    CMD: $PY tools/dump_diff.py --glob --ignore-step --rel-floor=$floor $strt/dump/node $splt/dump/node"
    "$PY" "$F3/tools/dump_diff.py" --glob --ignore-step --rel-floor="$floor" "$strt/dump/node" "$splt/dump/node"
    if [ $? = 0 ]; then echo "PASS  #2 live dump $([ "$floor" = 0 ] && echo 'max|Δ|=0' || echo '<=1-ULP (np>1 floor)') (node dyn/tracer state)"
    else echo "FAIL  #2 live dump diverged (see above)"; rc=1; fi
    echo
}

for NP in $NPS; do
    do_case "$NP" mid "$MID_START" no
    do_case "$NP" bnd "$BND_START" yes
done

echo "==================================================================================="
[ "$rc" = 0 ] \
    && echo "run_restart_gate_core2: PRIMARY GATE GREEN — split == straight-through (np: $NPS; mid-run + boundary; ice+sigma). max|Δ|=0 exact at every np$([ "${FESOM3_SHARED_SYNC:-1}" = 0 ] && echo '; np>1 with FESOM3_SHARED_SYNC=0: <=1-ULP redundant-element floor (1e-12)')" \
    || echo "run_restart_gate_core2: PRIMARY GATE FAILED (np: $NPS) — see FAIL lines above"
exit "$rc"
