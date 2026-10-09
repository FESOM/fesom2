#!/usr/bin/env bash
# Test of the mesh tools (check_mesh, remove_deadends, split_quads, cut_region) on the
# NEMO ORCA025 quad mesh made by nemo_meshq. Builds, under <out_root>,
#   nemo025q_fesom/  : dead ends removed           (remove_deadends.py)
#   nemo025mix/      : band 5S-5N + every 5th quad split into triangles (split_quads.py)
#   nemo025q_cut/, nemo025mix_cut/ : lon -20..-5, lat -2..12 (Gulf of Guinea, West
#                      African coast; cut_region.py)
# and checks:
#   - the original mesh is valid except for exactly 40 dead ends;
#   - the cleaned, mixed and cut meshes pass check_mesh (0 dead ends);
#   - split: element count = elements + split count, triangles = 2 x split count,
#     total area preserved to roundoff, nod2d.out and depth.out unchanged;
#   - cut-outs: a few thousand nodes, the mixed cut-out holds quads and triangles, every
#     node layer count = max over its elements (the model's bottom invariant);
#   - error cases: flipped quad, flipped triangle, bow-tie quad, repeated vertex and a node
#     layer count above the bottom invariant are reported.
#
# Usage: tools/mesh/test_mesh_tools.sh [nemo_dir] [out_root]
#   out_root defaults to a scratch dir (removed at exit); /work/ab0995/a270029/meshes makes the test meshes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEMO="${1:-/work/ab0995/a270029/meshes/nemo025q}"
SCR="$(mktemp -d "${TMPDIR:-/tmp}/mesh_tools.XXXXXX")"
trap 'rm -rf "$SCR"' EXIT
ROOT="${2:-$SCR}"                   # the scratch dir (removed at exit) unless given

source /usr/share/Modules/init/bash
module load python3 >/dev/null 2>&1
hash -r

fail=0
ok()  { echo "  ok: $*"; }
bad() { echo "  FAILED: $*"; fail=1; }
# value of key $2 in the SUMMARY/SPLIT line of output file $1
val() { tr ' ' '\n' < "$1" | sed -n "s/^$2=//p" | tail -1; }
# check_mesh <dir> <log>: runs check_mesh.py, returns its exit status
check() { python3 "$HERE/check_mesh.py" "$1" > "$2" 2>&1; }
# expect_pass <dir> <name>
expect_pass() {
    if check "$1" "$SCR/$2.log" && [ "$(val "$SCR/$2.log" deadends)" = 0 ]; then
        ok "$2 passes check_mesh ($(grep '^  nodes' "$SCR/$2.log" | sed 's/^  //'))"
    else
        bad "$2 check_mesh:"; cat "$SCR/$2.log"
    fi
}

echo "[1/5] check_mesh on the original mesh $NEMO"
if check "$NEMO" "$SCR/orig.log"; then
    bad "original mesh passes (40 dead ends expected)"
else
    n=$(val "$SCR/orig.log" deadends)
    others=0
    for k in ccw nonconvex degenerate unused multiedge; do
        [ "$(val "$SCR/orig.log" $k)" = 0 ] || others=1
    done
    if [ "$n" = 40 ] && [ $others = 0 ]; then
        ok "original mesh: exactly 40 dead ends, no other error"
    else
        bad "original mesh: $n dead ends (40 expected) or other errors"; cat "$SCR/orig.log"
    fi
fi

echo "[2/5] remove_deadends -> $ROOT/nemo025q_fesom"
python3 "$HERE/remove_deadends.py" "$NEMO" "$ROOT/nemo025q_fesom"
expect_pass "$ROOT/nemo025q_fesom" cleaned

echo "[3/5] split_quads -> $ROOT/nemo025mix"
python3 "$HERE/split_quads.py" "$ROOT/nemo025q_fesom" "$ROOT/nemo025mix" \
    --band -5 5 --every 5 | tee "$SCR/split.log"
expect_pass "$ROOT/nemo025mix" mixed
ein=$(val "$SCR/split.log" elems_in); ns=$(val "$SCR/split.log" split)
eout=$(val "$SCR/split.log" elems_out)
if [ "$ns" -gt 0 ] && [ "$eout" = $((ein + ns)) ] && [ "$(val "$SCR/mixed.log" elems)" = "$eout" ] \
   && [ "$(val "$SCR/mixed.log" tris)" = $((2 * ns)) ]; then
    ok "element count $eout = $ein + $ns split, triangles = 2 x $ns"
else
    bad "element counts: in $ein, split $ns, out $eout"
fi
if python3 -c 'import sys; a, b = map(float, sys.argv[1:]); print("    area %.15e -> %.15e, rel %.1e" % (a, b, abs(b-a)/a)); sys.exit(abs(b-a) > 1e-12*a)' \
     "$(val "$SCR/cleaned.log" area)" "$(val "$SCR/mixed.log" area)"; then
    ok "total area preserved to roundoff"
else
    bad "total area changed by the split"
fi
if cmp -s "$ROOT/nemo025q_fesom/nod2d.out" "$ROOT/nemo025mix/nod2d.out" \
   && cmp -s "$ROOT/nemo025q_fesom/depth.out" "$ROOT/nemo025mix/depth.out"; then
    ok "nod2d.out and depth.out unchanged by the split"
else
    bad "split changed nod2d.out or depth.out"
fi

echo "[4/5] cut_region -> $ROOT/nemo025q_cut, $ROOT/nemo025mix_cut"
for m in nemo025q_fesom:nemo025q_cut nemo025mix:nemo025mix_cut; do
    src=${m%%:*}; dst=${m##*:}
    python3 "$HERE/cut_region.py" "$ROOT/$src" "$ROOT/$dst" --lon -20 -5 --lat -2 12 \
        | tee "$SCR/${dst}_cut.log"
    expect_pass "$ROOT/$dst" "$dst"
    n=$(val "$SCR/$dst.log" nodes)
    if [ "$n" -ge 1000 ] && [ "$n" -le 10000 ]; then
        ok "$dst: $n nodes"
    else
        bad "$dst: $n nodes (a few thousand expected)"
    fi
    # the model's bottom invariant: node count = max over its elements of the element
    # count (min over the vertices), also for the nodes on the cut (check_mesh layers=)
    if [ "$(val "$SCR/$dst.log" layers)" = 0 ]; then
        ok "$dst: node layer counts = max over the elements ($(grep -o '[0-9]* node layer count' "$SCR/${dst}_cut.log" | cut -d' ' -f1) lowered by the cut)"
    else
        bad "$dst: bottom invariant violated"
    fi
done
if [ "$(val "$SCR/nemo025mix_cut.log" quads)" -gt 0 ] && [ "$(val "$SCR/nemo025mix_cut.log" tris)" -gt 0 ]; then
    ok "nemo025mix_cut holds quads and triangles"
else
    bad "nemo025mix_cut is not mixed"
fi

echo "[5/5] error cases on copies of nemo025mix_cut"
# mutate <name> <python expression on e (elements, 0-based, 4 columns) and q/t (first quad/triangle)>
mutate() {
    mkdir -p "$SCR/$1"
    cp "$ROOT/nemo025mix_cut/"{nod2d.out,depth.out} "$SCR/$1/"
    python3 - "$ROOT/nemo025mix_cut/elem2d.out" "$SCR/$1/elem2d.out" "$2" <<'EOF'
import sys
import numpy as np
tok = open(sys.argv[1]).read().split()
e = np.array(tok[1:], dtype=np.int64).reshape(-1, 4)
q = np.nonzero(e[:, 3] != e[:, 0])[0][0]
t = np.nonzero(e[:, 3] == e[:, 0])[0][0]
exec(sys.argv[3])
with open(sys.argv[2], 'w') as f:
    f.write('%9i \n' % len(e))
    f.write(''.join('%9i %9i %9i %9i\n' % tuple(r) for r in e.tolist()))
EOF
}
# expect_error <name> <key> <count>
expect_error() {
    if check "$SCR/$1" "$SCR/$1.log"; then
        bad "$1 passes check_mesh"
    elif [ "$(val "$SCR/$1.log" "$2")" = "$3" ]; then
        ok "$1 reported: $(grep ERROR "$SCR/$1.log" | sed 's/^  //; s/  */ /g')"
    else
        bad "$1: $2 = $(val "$SCR/$1.log" "$2"), $3 expected"; cat "$SCR/$1.log"
    fi
}
mutate flipped 'e[q, [1, 3]] = e[q, [3, 1]]; e[t, [1, 2]] = e[t, [2, 1]]'
expect_error flipped ccw 2
mutate bowtie 'e[q, [1, 2]] = e[q, [2, 1]]'
expect_error bowtie nonconvex 1
mutate repeated 'e[q, 2] = e[q, 1]'
expect_error repeated degenerate 1
# one node layer count raised above the max over its elements (depth.out block 2)
mutate deeper 'pass'
python3 - "$SCR/deeper/depth.out" <<'EOF'
import sys
lines = open(sys.argv[1]).read().split('\n')
nl = int(lines[0]); n2d = (len([l for l in lines if l != '']) - 1 - nl) // 2
i = 1 + nl + n2d                    # the layer count of node 1
lines[i] = '%9i ' % (int(lines[i]) + 1)
open(sys.argv[1], 'w').write('\n'.join(lines))
EOF
expect_error deeper layers 1

if [ "$fail" -ne 0 ]; then
    echo "FAIL"; exit 1
fi
echo "PASS"
