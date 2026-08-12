#!/bin/bash
# build_variant.sh <frozen-name> [PROBE_FLAG ...] — atomic, provenance-stamped
# build of one probe variant from THIS worktree into
# /work/ab0995/a270088/sp_paper/frozen/<name>/{bin,lib64} + PROVENANCE.txt.
# Refuses to build from a dirty tree, uses a variant-specific BUILD_DIR, and
# sha256-stamps the artifacts, so a frozen binary can always be traced to
# (branch, commit, flags) and never mixed up with another variant.
set -euo pipefail
cd "$(dirname "$0")"
name=$1; shift
flags=""
for f in "$@"; do flags="$flags -D$f"; done
flags="${flags# }"
[ -z "$(git status --porcelain)" ] || { echo "ERROR: dirty tree — commit first"; exit 1; }
commit=$(git rev-parse --short=12 HEAD); branch=$(git rev-parse --abbrev-ref HEAD)
export BUILD_DIR=build_$name
# PRECIS=dp builds a double-precision variant (for DP-base probes); default sp
precis=${PRECIS:-sp}
case $precis in
  sp) pflag="-DUSE_SINGLE_PRECISION=ON" ;;
  dp) pflag="-DUSE_SINGLE_PRECISION=OFF" ;;
  *)  echo "ERROR: PRECIS must be sp or dp"; exit 1 ;;
esac
if [ -n "$flags" ]; then
  ./configure.sh $pflag "-DCMAKE_Fortran_FLAGS=$flags"
else
  ./configure.sh $pflag
fi
F=/work/ab0995/a270088/sp_paper/frozen/$name
mkdir -p $F/bin $F/lib64
cp bin/fesom.x $F/bin/
cp lib64/libfesom.so $F/lib64/ 2>/dev/null || cp lib/libfesom.so $F/lib64/
{ echo "variant:  $name"
  echo "branch:   $branch"
  echo "commit:   $commit"
  echo "cpp:      ${flags:-<none>} (+ USE_SINGLE_PRECISION=$([ $precis = sp ] && echo ON || echo OFF))"
  echo "built:    $(date -Is) on $(hostname) by $USER"
  echo "builddir: $PWD/$BUILD_DIR"
  sha256sum $F/bin/fesom.x $F/lib64/libfesom.so
} > $F/PROVENANCE.txt
echo "== $F frozen:"; cat $F/PROVENANCE.txt
