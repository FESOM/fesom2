#!/usr/bin/env bash
# Regression test of tools/mesh/nemo_meshq.py (Python port of nemo_meshq.m):
# runs the port into a scratch dir and compares the three outputs against the
# MATLAB outputs saved in ref_matlab/.
#   elem2d.out, depth.out : byte-identical
#   nod2d.out             : byte-identical, or differing only in the coordinates of
#                           at most 20 lines (the fold f-points corrected with
#                           single-precision cos/atan2/acos, row 1020/1021), each by
#                           <= 1e-4 degrees; every differing line is listed.
#
# Usage: tools/mesh/test_nemo_meshq.sh [nemo_dir] [ref_dir] [scratch_dir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEMO="${1:-/work/ab0995/a270029/meshes/nemo025q}"
REF="${2:-$NEMO/ref_matlab}"
OUT="${3:-$(mktemp -d "${TMPDIR:-/tmp}/nemo_meshq.XXXXXX")}"

source /usr/share/Modules/init/bash
module load python3 >/dev/null 2>&1
hash -r

mkdir -p "$OUT"
echo "[1/2] nemo_meshq.py $NEMO -> $OUT"
python3 "$HERE/nemo_meshq.py" "$NEMO" "$OUT"

echo "[2/2] compare against $REF"
fail=0
for f in elem2d.out depth.out; do
    if cmp -s "$OUT/$f" "$REF/$f"; then
        echo "  $f: identical"
    else
        echo "  $f: DIFFERS"; cmp "$OUT/$f" "$REF/$f" || true; fail=1
    fi
done

if cmp -s "$OUT/nod2d.out" "$REF/nod2d.out"; then
    echo "  nod2d.out: identical"
else
    python3 - "$OUT/nod2d.out" "$REF/nod2d.out" <<'EOF' || fail=1
import sys
a = open(sys.argv[1]).read().split('\n')
b = open(sys.argv[2]).read().split('\n')
if len(a) != len(b):
    sys.exit("  nod2d.out: line count differs (%d vs %d)" % (len(a), len(b)))
bad, worst = [], 0.0
for i, (la, lb) in enumerate(zip(a, b)):
    if la == lb:
        continue
    fa, fb = la.split(), lb.split()
    if i == 0 or len(fa) != 4 or len(fb) != 4 or fa[0] != fb[0] or fa[3] != fb[3]:
        sys.exit("  nod2d.out: line %d differs beyond coordinates:\n    %r\n    %r" % (i + 1, la, lb))
    # in units of the last printed decimal (%8.4f), so 1e-4 is not lost to roundoff
    d = 1e-4 * max(abs(round(1e4 * float(fa[c])) - round(1e4 * float(fb[c]))) for c in (1, 2))
    worst = max(worst, d)
    bad.append((i + 1, la, lb, d))
for ln, la, lb, d in bad:
    print("    line %d: port %r  matlab %r  |d|=%.4g" % (ln, la, lb, d))
ok = len(bad) <= 20 and worst <= 1e-4
print("  nod2d.out: %d coordinate line(s) differ, max |d| = %.4g -> %s"
      % (len(bad), worst, "OK" if ok else "FAIL"))
sys.exit(0 if ok else 1)
EOF
fi

if [ "$fail" -ne 0 ]; then
    echo "FAIL"; exit 1
fi
echo "PASS"
