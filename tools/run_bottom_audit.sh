#!/usr/bin/env bash
# MESH-DELTA AUDIT for the bottom-at-vertices change (docs/plans/
# 20260910-fesom3-bottom-at-vertices.md, Task 8).
#
# Dumps the fully-built FESOM3 mesh geometry on pi and core2 and checks how the element
# bottom moved when it stopped being read from elvls.out and started being derived as the
# min over the element's vertex columns. The expected numbers were measured from the mesh
# files BEFORE any code was written, so this is the code answering back to the prediction.
# A disagreement means one of the two is wrong -- stop and reconcile rather than adjust
# the expectation.
#
#   pi     651 elements changed of 5839,  +1.55 % volume, 0 stagnant cells
#   core2  10465 changed of 244659,       +0.38 % volume, 0 stagnant cells
#
# tools/bottom_delta.py needs numpy-free plain python3 but a modern one; the login-node
# /usr/bin/python3 is 3.6.8 and works. PY can be overridden.
#
# Usage:  bash tools/run_bottom_audit.sh [run_dir]
set -euo pipefail

F3="${F3:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${BIN:-$F3/build_intel_dp/bin/fesom_geomdump}"
RUN="${1:-${RUN:-/tmp/bottom_audit_$USER}}"
PY="${PY:-python3}"
PIMESH="${PIMESH:-/home/a/a270088/port2/fesom2/tests/data/MESHES/pi}"
COREMESH="${COREMESH:-/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2}"
MPIFLAGS="${MPIFLAGS:---mca pml ob1 --mca btl self,vader --mca btl_vader_single_copy_mechanism none --oversubscribe}"  # LESSONS L61

[ -x "$BIN" ] || { echo "run_bottom_audit: missing $BIN (build first)"; exit 1; }
mkdir -p "$RUN"

fail=0

audit() {
    local name="$1" mesh="$2" expc="$3" expv="$4"
    if [ ! -d "$mesh" ]; then
        echo "=== $name: mesh not available at $mesh -- SKIPPED ==="
        return 0
    fi
    echo "=== $name ==="
    FESOM3_MESH_DIR="$mesh" FESOM3_GEOM_OUT="$RUN/geom_$name.bin" \
        mpirun $MPIFLAGS -n 1 "$BIN" > "$RUN/geomdump_$name.log" 2>&1 || {
            echo "run_bottom_audit: geomdump FAILED for $name (see $RUN/geomdump_$name.log)"
            tail -5 "$RUN/geomdump_$name.log"; fail=1; return 0; }
    "$PY" "$F3/tools/bottom_delta.py" "$RUN/geom_$name.bin" "$mesh" \
        --expect-changed "$expc" --expect-volume "$expv" || fail=1
    echo
}

audit pi    "$PIMESH"   651   1.55
audit core2 "$COREMESH" 10465 0.38

if [ "$fail" -ne 0 ]; then
    echo "run_bottom_audit: AUDIT FAILED"
    exit 1
fi
echo "run_bottom_audit: AUDIT OK"
