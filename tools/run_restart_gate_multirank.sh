#!/usr/bin/env bash
# Restart Stage 6 (Task 6.2) SECONDARY cross-partition gate — STORAGE + RESTORE round-trip.
#
# Proves F-B = partition-INDEPENDENT storage: a checkpoint written at np=2 is restored by a FRESH
# np=8 process and IMMEDIATELY re-written (ZERO model steps) — the re-written C8 must equal C2 exactly
# (max|Δ|=0), canonically. This is STORAGE+RESTORE, NOT cross-np EVOLUTION (which is physically non-
# bit-identical in FESOM — SSH-CG allreduce reduces in comm-size-dependent order; every project
# byte-gate is same-np). With zero steps there is NO order-dependent RHS, so even the 562 shared
# (multiply-owned) boundary elements round-trip EXACTLY: np=8 reads the canonical value into every owner, then
# re-dedups to the same canonical value => strict max|Δ|=0 (rel-floor 0), no element-ownership floor.
# (Since the shared-element sync, docs/plans/completed/20261012-shared-element-sync.md, the stepping runs keep
# the owners' copies equal too; this gate needs no step and is unchanged by it.)
#
#   tools/run_restart_gate_multirank.sh [run_dir] [NP_WRITE] [NP_READ]    # default 2 -> 8
set -uo pipefail
F3="${F3:-/home/a/a270088/fesom3/.claude/worktrees/restart-checkpoint}"
BUILD="${BUILD:-$F3/build_intel_dp}"
PY="${PY:-/work/ab0995/a270088/mambaforge/bin/python3}"
COREMESH=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2
ICFILE=/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc
POOL=/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0
NCAR=/home/a/a270088/port2/fesom2/test/input/global
RUN="${1:-/scratch/a/a270088/restart_gate_multirank}"
NPW="${2:-2}"          # np that WRITES C2
NPR="${3:-8}"          # np that RESTORES + re-writes C8
K="${K:-2}"            # cold steps before C2 is written
MID_START="0 1 1948"

source "$F3/env.sh" intel >/dev/null 2>&1
BIN="$BUILD/bin/fesom_lifecycle_native_mr"
[ -x "$BIN" ] || { echo "missing $BIN — build first"; exit 2; }
[ -d "$COREMESH/dist_$NPW" ] || { echo "missing $COREMESH/dist_$NPW"; exit 2; }
[ -d "$COREMESH/dist_$NPR" ] || { echo "missing $COREMESH/dist_$NPR"; exit 2; }
rm -rf "$RUN"; mkdir -p "$RUN"
C2="$RUN/c_np$NPW"; C8="$RUN/c_np$NPR"; mkdir -p "$C2" "$C8"

export FESOM3_MESH_DIR="$COREMESH" FESOM3_IC_FILE="$ICFILE" FESOM3_FORCING_DIR="$NCAR"
export FESOM3_RUNOFF_FILE="$POOL/CORE2_runoff.nc" FESOM3_SSS_FILE="$POOL/PHC2_salx.nc"
export FESOM3_WHICHEVP=0
ulimit -s unlimited
MPIRUN=(mpirun --mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe)
GREP='nod2D|INITIALISATION|RESTART RUN|state restored|registered|checkpoint written|final checkpoint'
rc=0

echo "=== WRITE: np=$NPW cold $K steps -> C2 (end-of-run checkpoint) ==="
env FESOM3_RESTART="$C2" FESOM3_RESTART_IN="$C2" FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1 \
    FESOM3_NSTEPS="$K" FESOM3_START_CLOCK="$MID_START" \
    "${MPIRUN[@]}" -n "$NPW" "$BIN" >"$RUN/write.log" 2>&1
s1=$?; grep -E "$GREP" "$RUN/write.log" | sed 's/^/   /'
[ "$s1" = 0 ] || { echo "FAIL  write seg exit $s1 (see $RUN/write.log)"; rc=1; }

echo "=== RESTORE: np=$NPR reads C2, ZERO steps -> re-writes C8 ==="
env FESOM3_RESTART="$C8" FESOM3_RESTART_IN="$C2" FESOM3_RESTART_UNIT=off FESOM3_RESTART_LENGTH=1 \
    FESOM3_NSTEPS=0 FESOM3_START_CLOCK="$MID_START" \
    "${MPIRUN[@]}" -n "$NPR" "$BIN" >"$RUN/restore.log" 2>&1
s2=$?; grep -E "$GREP" "$RUN/restore.log" | sed 's/^/   /'
[ "$s2" = 0 ] || { echo "FAIL  restore seg exit $s2 (see $RUN/restore.log)"; rc=1; }
grep -q 'THIS IS A RESTART RUN' "$RUN/restore.log" \
    && echo "PASS  np=$NPR resumed (r_restart=.true.)" || { echo "FAIL  np=$NPR did not resume"; rc=1; }

c2f=$(cat "$C2/restart.latest" 2>/dev/null); c8f=$(cat "$C8/restart.latest" 2>/dev/null)
echo "=== COMPARE  C2=$c2f  (np$NPW)  vs  C8=$c8f  (np$NPR), STRICT max|Δ|=0 ==="
if [ -z "$c2f" ] || [ -z "$c8f" ]; then
    echo "FAIL  missing checkpoint ('$c2f' / '$c8f')"; rc=1
else
    echo "    CMD: $PY tools/zarr_diff.py --output-cmp $C2/$c2f $C8/$c8f --rel-floor 0"
    "$PY" "$F3/tools/zarr_diff.py" --output-cmp "$C2/$c2f" "$C8/$c8f" --rel-floor 0
    if [ $? = 0 ]; then echo "PASS  C2 == C8 max|Δ|=0 (partition-independent STORAGE+RESTORE, np$NPW->np$NPR)"
    else echo "FAIL  C2 != C8 (see per-store max|Δ| above)"; rc=1; fi
fi

echo "==================================================================================="
[ "$rc" = 0 ] \
    && echo "run_restart_gate_multirank: SECONDARY GATE GREEN — C2(np$NPW) == C8(np$NPR) max|Δ|=0 (partition-independent restore)" \
    || echo "run_restart_gate_multirank: SECONDARY GATE FAILED — see FAIL lines above"
exit "$rc"
