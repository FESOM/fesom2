#!/usr/bin/env bash
# fesom_ini on quad and mixed (triangle + quad) meshes (hybrid-mesh plan, Task 4).
#
# 1. Runs fesom_ini (use_depthfile='nemo': depth.out) on the writable NEMO test meshes under
#    MESH_ROOT and checks every output with tools/mesh/check_dist.py (edges, edge order and
#    edge2D_in, levels == depth.out counts, dist_N ownership, halos, communication):
#      nemo025q_cut, nemo025mix_cut (2, 4, 8 parts), nemo025q_fesom, nemo025mix (2, 8, 128),
#    plus dist_1 of each (tools/make_dist1.py). Old outputs (edges, levels, dist_*) are
#    removed first: fesom_ini keeps existing edge/level files instead of rewriting them.
#    Since Task 5 check_dist.py also checks the appended sections (dist_N/format.out:
#    double-halo nodes, com_nod2D_full, com_edge2D); the 4-part cut-outs must have edges
#    in overlapping owned ranges (an edge in 1..myDim_edge2D of two ranks).
# 2. check_dist.py on the triangle meshes: the writable copies pi_w (dist_1/2/4/8), core2_w
#    (dist_1/2/4/8/128) (made by tools/run_ini_compare.sh) and, read-only, the existing
#    dist_N of the original pi (dist_2, dist_8) and CORE2 (dist_8) made by older fesom_ini
#    versions (FESOM2 layout only).
# 3. Error cases (on copies under WORK_ROOT):
#    - the original NEMO mesh (40 dead-end quads) stops with the "Insufficient number of
#      neighbors" message and writes no dist_N;
#    - a mixed cut-out with one anticlockwise quad and one anticlockwise triangle: both are
#      fixed and reported, and every output equals that of the unflipped mesh;
#    - a self-intersecting (bow-tie) quad stops;
#    - check_dist.py fails on a dist_2 with an owned and a halo node exchanged, and on a
#      dist_4 with two double-halo nodes exchanged, a com_nod2D_full or a com_edge2D
#      receive entry duplicated.
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
source "$HERE/lib_ini.sh"          # python3 module, ok/bad/fail, OUTS, ini_namelist
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
    ini_namelist "$rd/namelist.config" "$mesh" "$np" nemo
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
for spec in nemo025q_cut:2,4,8 nemo025mix_cut:2,4,8 nemo025q_fesom:2,8,128 nemo025mix:2,8,128; do
    name=${spec%%:*}; parts=${spec#*:}; d="$MESH_ROOT/$name"
    clean "$d"
    for np in ${parts//,/ }; do
        log="$WORK_ROOT/log/${name}_$np.log"
        run_ini "$d" "$np" "$log"
        [ -f "$d/dist_$np/rpart.out" ] && grep -q 'Total time' "$log" \
            || bad "fesom_ini $name np=$np (log $log)"
    done
    python3 "$HERE/make_dist1.py" "$d" > "$WORK_ROOT/log/${name}_1.log" 2>&1 \
        || bad "make_dist1.py $name (log $WORK_ROOT/log/${name}_1.log)"
    check "$d" 1 ${parts//,/ }
    if [[ ",$parts," == *",4,"* ]]; then
        n=$(sed -n '/dist_4$/,/^==\|check_dist:/p' "$WORK_ROOT/log/check_${name}.log" \
            | sed -n 's/.*(\([0-9]*\) edges in overlapping owned ranges).*/\1/p')
        [ "${n:-0}" -gt 0 ] && ok "$name dist_4: $n edges in overlapping owned ranges, received from their unique owner" \
                            || bad "$name dist_4: no edges in overlapping owned ranges"
    fi
done

echo "[2/3] check_dist on the triangle meshes"
for spec in "$MESH_ROOT/pi_w:1,2,4,8" "$MESH_ROOT/core2_w:1,2,4,8,128" "$PI_SRC:2,8" "$CORE2_SRC:8"; do
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
elif grep -q 'FAILED: every node owned by exactly one rank' "$WORK_ROOT/log/check_corrupt.log"; then
    ok "check_dist detects a corrupted dist_2 ($(grep -c FAILED "$WORK_ROOT/log/check_corrupt.log") failed checks)"
else
    bad "check_dist fails a corrupted dist_2 for another reason (log $WORK_ROOT/log/check_corrupt.log)"
fi

# check_dist must fail on corrupted appended sections (copies of the mixed cut-out dist_4)
for what in xnodes nod_full edge; do
    d="$WORK_ROOT/corrupt_$what"; mkdir -p "$d"; clean "$d"
    for f in nod2d.out elem2d.out depth.out $OUTS; do cp "$ref/$f" "$d/"; done
    cp -r "$ref/dist_4" "$d/"
    python3 - "$d/dist_4" "$what" <<'EOF'
import sys
d, what = sys.argv[1], sys.argv[2]
def load(f):
    return open(f).read().split()
def save(f, tok):
    open(f, 'w').write('\n'.join(tok) + '\n')
if what == 'xnodes':                   # exchange the first two double-halo nodes of rank 0
    f = d + '/my_list00000.out'
    tok = load(f)
    i = 3 + int(tok[1]) + int(tok[2])                     # myDim_elem2D
    i += 3 + sum(int(x) for x in tok[i:i + 3])            # myDim_edge2D
    i += 2 + int(tok[i]) + int(tok[i + 1])                # eXDim_nod2D
    assert int(tok[i]) >= 2
    tok[i + 1], tok[i + 2] = tok[i + 2], tok[i + 1]
else:                                  # duplicate a receive entry of rank 1
    f = d + '/com_info00001.out'
    tok = load(f)
    i = 1
    for s in range(4 if what == 'nod_full' else 5):       # skip to the wanted structure
        start = i
        for _ in range(2):                                 # receive, send
            npe = int(tok[i]); ptr_end = int(tok[i + 2 * npe + 1]) if npe else 1
            i += 1 + npe + npe + 1 + ptr_end - 1
    npe = int(tok[start])
    r0 = start + 1 + npe + npe + 1                         # first rlist entry
    assert npe > 0 and int(tok[start + 2 * npe + 1]) > 2
    tok[r0 + 1] = tok[r0]
save(f, tok)
EOF
    case $what in
        xnodes)   want='FAILED: eXDim_nod2D == vertices of the local elements' ;;
        nod_full) want='FAILED: com_nod2D_full: every halo node' ;;
        edge)     want='FAILED: com_edge2D: every held edge' ;;
    esac
    if python3 "$HERE/mesh/check_dist.py" "$d" 4 > "$WORK_ROOT/log/check_corrupt_$what.log" 2>&1; then
        bad "check_dist passes a dist_4 with corrupted $what"
    elif grep -q "$want" "$WORK_ROOT/log/check_corrupt_$what.log"; then
        ok "check_dist detects corrupted $what in dist_4 ($(grep FAILED "$WORK_ROOT/log/check_corrupt_$what.log" | head -1 | cut -c1-70)...)"
    else
        bad "check_dist fails corrupted $what for another reason (log $WORK_ROOT/log/check_corrupt_$what.log)"
    fi
done

if [ "$fail" = 0 ]; then echo "run_ini_mixed: PASS"; else echo "run_ini_mixed: FAIL"; exit 1; fi
