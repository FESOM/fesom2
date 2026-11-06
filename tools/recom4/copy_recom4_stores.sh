#!/bin/sh
# copy_recom4_stores.sh FROM TO
#
# RECOM4 M10 (#57): biology from another run. Copies RECOM4's part of the FESOM_F
# checkpoint FROM (a run with RECOM4 on) into the checkpoint TO (for example a
# physics spin-up), so that a run resumed from TO with REcoM_restart = .true.
# starts with FROM's biology and TO's physics and clock:
#   - the tracer stores tra0003 ... with their history (_AB, _M1, _M2),
#   - the benthic pools BenN, BenC, BenSi, BenCalc,
#   - recom4_restart_tag.
#
# FROM and TO are checkpoint folders (fesom.<tag>/, with checkpoint.json), not
# restart directories. TO is changed in place: copy it first to keep it.
# Both must be on the same mesh (checked on the node coordinates of temp). If
# TO was written with RECOM4 on, its restart tag must equal FROM's (the same
# configuration); the run checks the tag again when it reads TO.
#
# The tracers' AB history comes from FROM's run, with FROM's velocities; the
# first step after the restart uses it.

set -eu

fail() {
    echo "copy_recom4_stores: $*" >&2
    exit 1
}

[ $# -eq 2 ] || { echo "usage: $0 FROM_CHECKPOINT TO_CHECKPOINT" >&2; exit 2; }
from=${1%/}
to=${2%/}

for d in "$from" "$to"; do
    [ -f "$d/checkpoint.json" ] || fail "$d is not a checkpoint folder (no checkpoint.json)"
done
[ "$(cd "$from" && pwd -P)" != "$(cd "$to" && pwd -P)" ] || fail "FROM and TO are the same folder"
[ -d "$from/recom4_restart_tag.zarr" ] || fail "$from has no recom4_restart_tag.zarr: not written with RECOM4 on"

# same mesh: the node coordinates stored with temp must be identical
for c in lon lat; do
    [ -d "$from/temp.zarr/$c" ] && [ -d "$to/temp.zarr/$c" ] || \
        fail "temp.zarr/$c missing in FROM or TO: cannot check the mesh"
    diff -r "$from/temp.zarr/$c" "$to/temp.zarr/$c" > /dev/null || \
        fail "node coordinates ($c) differ: not the same mesh"
done

# same configuration, if TO was written with RECOM4 on
if [ -d "$to/recom4_restart_tag.zarr" ]; then
    diff -r "$from/recom4_restart_tag.zarr" "$to/recom4_restart_tag.zarr" > /dev/null || \
        fail "restart tags differ: FROM and TO are different RECOM4 configurations"
fi

# the stores to copy: tracers 3 onwards, the pools, the tag
list=""
for s in "$from"/tra[0-9][0-9][0-9][0-9]*.zarr; do
    [ -d "$s" ] || continue
    b=$(basename "$s")
    j=$(echo "$b" | sed 's/^tra\([0-9][0-9][0-9][0-9]\).*/\1/')
    [ "$j" -ge 3 ] || continue
    list="$list $b"
done
[ -n "$list" ] || fail "$from has no tracer stores tra0003 onwards"
for b in BenN BenC BenSi BenCalc recom4_restart_tag; do
    [ -d "$from/$b.zarr" ] || fail "$from has no $b.zarr"
    list="$list $b.zarr"
done

n=0
for b in $list; do
    rm -rf "${to:?}/$b"
    cp -R "$from/$b" "$to/$b"
    n=$((n + 1))
done
echo "copy_recom4_stores: $n stores copied from $from into $to"
