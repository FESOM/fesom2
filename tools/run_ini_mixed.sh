#!/usr/bin/env bash
# fesom_ini on quad and mixed (triangle + quad) meshes (hybrid-mesh plan, Task 4).
#
# 1. Runs fesom_ini (use_depthfile='nemo': depth.out) on the writable NEMO test meshes under
#    MESH_ROOT and checks every output with tools/mesh/check_dist.py (edges, edge order and
#    edge2D_in, levels == depth.out counts, dist_N ownership, halos, communication):
#      nemo025q_cut, nemo025mix_cut (2, 8 parts), nemo025q_fesom, nemo025mix (2, 8, 128).
#    Old outputs (edges, levels, dist_*) are removed first: fesom_ini keeps existing
#    edge/level files instead of rewriting them.
# 2. check_dist.py on the triangle meshes: the writable copies pi_w, core2_w (made by
#    tools/run_ini_compare.sh) and, read-only, the existing dist_N of the original pi
#    (dist_2, dist_8) and CORE2 (dist_8) made by older fesom_ini versions.
# 3. Error cases (on copies under WORK_ROOT):
#    - the original NEMO mesh (40 dead-end quads) stops with the "Insufficient number of
#      neighbors" message and writes no dist_N;
#    - a mixed cut-out with one anticlockwise quad and one anticlockwise triangle: both are
#      fixed and reported, and every output equals that of the unflipped mesh;
#    - a self-intersecting (bow-tie) quad stops;
#    - check_dist.py fails on a dist_2 with an owned and a halo node exchanged.
# tools/run_ini_compare.sh (triangles byte-identical to FESOM2's build) is run separately.
#
# Usage: tools/run_ini_mixed.sh [build_dir]
#   build_dir defaults to /work/ab0995/a270029/fesom3_hybrid/build (holds bin/fesom_ini)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${1:-/work/ab0995/a270029/fesom3_hybrid/build}"
INI="$BUILD/bin/fesom_ini"
MESH_ROOT="${MESH_ROOT:-/work/ab0995/a270029/meshes}"
WORK_ROOT="${WORK_ROOT:-/work/ab0995/a270029/fesom3_hybrid/ini_mixed}"
NEMO_ORIG="$MESH_ROOT/nemo025q"
PI_SRC=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
CORE2_SRC=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2

[ -x "$INI" ] || { echo "run_ini_mixed: missing executable $INI" >&2; exit 2; }
ulimit -s unlimited
source /usr/share/Modules/init/bash
module load python3 >/dev/null 2>&1
hash -r

fail=0
ok()  { echo "  ok: $*"; }
bad() { echo "  FAILED: $*"; fail=1; }
OUTS="edges.out edge_tri.out edgenum.out elvls.out elvls_raw.out nlvls.out"
mkdir -p "$WORK_ROOT/log"

# clean <mesh_dir>: removes the fesom_ini outputs
clean() {
    local f
    for f in $OUTS; do rm -f "${1:?}/$f"; done
    rm -rf "${1:?}"/dist_*
}

# run_ini <mesh_dir> <npart> <log>: one fesom_ini run in its own working directory
# (fesom_ini STOPs with exit status 0: callers check the log and the outputs)
run_ini() {
    local mesh="$1" np="$2" log="$3" rd
    rd="$WORK_ROOT/run/$(basename "$mesh")_$np"
    mkdir -p "$rd"
    printf '&paths\nMeshPath="%s/"\n/\n&geometry\nuse_depthfile="nemo"\n/\n&run_config\n/\n&machine\nn_levels=1\nn_part=%s\n/\n' \
        "$mesh" "$np" > "$rd/namelist.config"
    (cd "$rd" && "$INI" > "$log" 2>&1) || true
}

# check <mesh_dir> <nparts...>
check() {
    local d="$1" name log; shift
    name=$(basename "$d"); log="$WORK_ROOT/log/check_${name}.log"
    if python3 "$HERE/mesh/check_dist.py" "$d" "$@" > "$log" 2>&1; then
        ok "check_dist $d ($*): $(grep -c '  ok:' "$log") checks"
    else
        bad "check_dist $d ($*):"; grep -v '  ok:' "$log"
    fi
}

echo "[1/3] fesom_ini + check_dist on the NEMO quad and mixed meshes"
for spec in nemo025q_cut:2,8 nemo025mix_cut:2,8 nemo025q_fesom:2,8,128 nemo025mix:2,8,128; do
    name=${spec%%:*}; parts=${spec#*:}; d="$MESH_ROOT/$name"
    clean "$d"
    for np in ${parts//,/ }; do
        log="$WORK_ROOT/log/${name}_$np.log"
        run_ini "$d" "$np" "$log"
        [ -f "$d/dist_$np/rpart.out" ] && grep -q 'Total time' "$log" \
            || bad "fesom_ini $name np=$np (log $log)"
    done
    check "$d" ${parts//,/ }
done

echo "[2/3] check_dist on the triangle meshes"
for spec in "$MESH_ROOT/pi_w:2,8" "$MESH_ROOT/core2_w:2,8,128" "$PI_SRC:2,8" "$CORE2_SRC:8"; do
    d=${spec%%:*}; parts=${spec#*:}
    if [ -f "$d/dist_${parts%%,*}/rpart.out" ]; then
        check "$d" ${parts//,/ }
    else
        bad "$d/dist_${parts%%,*} missing (run tools/run_ini_compare.sh first)"
    fi
done

echo "[3/3] error cases"
# original NEMO mesh: dead ends (inputs linked read-only, outputs land in the copy)
d="$WORK_ROOT/nemo025q_orig"; mkdir -p "$d"; clean "$d"
for f in nod2d.out elem2d.out depth.out; do ln -sf "$NEMO_ORIG/$f" "$d/$f"; done
run_ini "$d" 2 "$WORK_ROOT/log/deadend.log"
if grep -q 'Insufficient number of neighbors' "$WORK_ROOT/log/deadend.log" && [ ! -d "$d/dist_2" ]; then
    ok "original NEMO mesh stops on a dead end ($(grep -m1 'Insufficient' "$WORK_ROOT/log/deadend.log" | xargs))"
else
    bad "original NEMO mesh: no dead-end stop (log $WORK_ROOT/log/deadend.log)"
fi

# flipped quad + flipped triangle in a copy of the mixed cut-out
ref="$MESH_ROOT/nemo025mix_cut"; d="$WORK_ROOT/flip"; mkdir -p "$d"; clean "$d"
cp "$ref/nod2d.out" "$ref/depth.out" "$d/"
python3 - "$ref/elem2d.out" "$d/elem2d.out" <<'EOF'
import sys
lines = open(sys.argv[1]).read().split('\n')
el = [l.split() for l in lines[1:] if l.strip()]
q = next(i for i, e in enumerate(el) if e[0] != e[3])          # first quad: 2 <-> 4
t = next(i for i, e in enumerate(el) if e[0] == e[3])          # first triangle: 2 <-> 3
el[q][1], el[q][3] = el[q][3], el[q][1]
el[t][1], el[t][2] = el[t][2], el[t][1]
el[t][3] = el[t][0]
with open(sys.argv[2], 'w') as f:
    f.write(lines[0] + '\n' + ''.join('%9s %9s %9s %9s\n' % tuple(e) for e in el))
print('    flipped quad %d and triangle %d' % (q + 1, t + 1))
EOF
log="$WORK_ROOT/log/flip.log"
run_ini "$d" 2 "$log"
if grep -q 'vertices exchanged for quad' "$log" && grep -q 'vertices exchanged for triangle' "$log" \
   && grep -q 'made clockwise for *2 *elements' "$log"; then
    ok "flipped quad and triangle reported"
else
    bad "flipped elements not reported (log $log)"
fi
n=0
for f in $OUTS $(cd "$ref/dist_2" && ls | sed 's|^|dist_2/|'); do
    cmp -s "$d/$f" "$ref/$f" || n=$((n+1))
done
[ "$n" = 0 ] && ok "flipped mesh: all outputs identical to the unflipped mesh" \
             || bad "flipped mesh: $n output files differ from the unflipped mesh"

# bow-tie quad (nodes 2 <-> 3 of the first quad)
d="$WORK_ROOT/bowtie"; mkdir -p "$d"; clean "$d"
cp "$ref/nod2d.out" "$ref/depth.out" "$d/"
awk 'NR>1 && !done && $1!=$4 {t=$2; $2=$3; $3=t; done=1} {print}' "$ref/elem2d.out" > "$d/elem2d.out"
run_ini "$d" 2 "$WORK_ROOT/log/bowtie.log"
if grep -q 'self-intersecting or degenerate quad' "$WORK_ROOT/log/bowtie.log" && [ ! -d "$d/dist_2" ]; then
    ok "bow-tie quad stops"
else
    bad "bow-tie quad does not stop (log $WORK_ROOT/log/bowtie.log)"
fi

# check_dist must fail on a corrupted partition: owned node <-> halo node in rank 0
d="$WORK_ROOT/flip"
python3 - "$d/dist_2/my_list00000.out" <<'EOF'
import sys
tok = open(sys.argv[1]).read().split()
myd = int(tok[1])                      # tokens: rank, myDim, eDim, nodes...
tok[3], tok[3 + myd] = tok[3 + myd], tok[3]
open(sys.argv[1], 'w').write('\n'.join(tok) + '\n')
EOF
if python3 "$HERE/mesh/check_dist.py" "$d" 2 > "$WORK_ROOT/log/check_corrupt.log" 2>&1; then
    bad "check_dist passes a corrupted dist_2"
else
    ok "check_dist detects a corrupted dist_2 ($(grep -c FAILED "$WORK_ROOT/log/check_corrupt.log") failed checks)"
fi

if [ "$fail" = 0 ]; then echo "run_ini_mixed: PASS"; else echo "run_ini_mixed: FAIL"; exit 1; fi
