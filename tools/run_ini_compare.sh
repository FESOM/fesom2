#!/usr/bin/env bash
# Byte comparison of the ported fesom_ini (mesh_part/, FESOM2's mesh partitioner) on the
# triangle meshes pi and CORE2 (hybrid-mesh plan, Task 3).
#
# Writable copies of the read-only meshes (nod2d.out, elem2d.out, aux3d.out) are kept under
# MESH_ROOT (pi_w/, core2_w/); fesom_ini writes its outputs into them. Each run starts from
# the three input files only (old outputs are removed first: FESOM2's fesom_ini keeps
# existing edge/level files instead of rewriting them).
#
# Two references:
#   1. FESOM2's own mesh_part build (F2_INI), run on separate copies under WORK_ROOT:
#      EVERY output must be byte-identical (edges.out, edge_tri.out, edgenum.out,
#      elvls.out, elvls_raw.out, nlvls.out, dist_N/*). This is the port check.
#   2. The files in the original mesh directories. Byte-identical are the ones produced by
#      this FESOM2 version: edges.out, edge_tri.out, edgenum.out (pi, CORE2), CORE2
#      elvls/nlvls.out.20260528_regenerated, CORE2 dist_2 and dist_128. Made by older
#      fesom_ini versions and therefore NOT reproducible (reported, not gated):
#        - pi elvls.out/nlvls.out (3 element levels differ) and pi dist_2/dist_8, CORE2
#          dist_8 (2023), CORE2 elvls.out/nlvls.out (the 2021 originals): their partitions
#          are nearly perfectly node-balanced, which suggests METIS was called without the
#          node-level weights (PART_WEIGHTED); a scratch run with unit weights matches the pi
#          dist_2 part sizes but not the node lists (part labels swapped), one with the
#          existing levels as weights matches nothing either.
#      Their structure is checked by check_dist.py (Task 4).
#
# FESOM2's mesh_part is built once (scratch) with
#   cmake -S /home/a/a270029/fesom2.0/mesh_part -B <dir> -DCMAKE_BUILD_TYPE=Release \
#         -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc   (after source env.sh intel)
#
# Usage: tools/run_ini_compare.sh [build_dir]
#   build_dir defaults to /work/ab0995/a270029/fesom3_hybrid/build (holds bin/fesom_ini)
set -euo pipefail

BUILD="${1:-/work/ab0995/a270029/fesom3_hybrid/build}"
INI="$BUILD/bin/fesom_ini"
F2_INI="${F2_INI:-/work/ab0995/a270029/fesom3_hybrid/f2_meshpart_build/bin/fesom_meshpart}"
MESH_ROOT="${MESH_ROOT:-/work/ab0995/a270029/meshes}"
WORK_ROOT="${WORK_ROOT:-/work/ab0995/a270029/fesom3_hybrid/ini_compare}"
PI_SRC=/home/a/a270088/port2/fesom2/tests/data/MESHES/pi
CORE2_SRC=/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2

for x in "$INI" "$F2_INI"; do
    [ -x "$x" ] || { echo "run_ini_compare: missing executable $x" >&2; exit 2; }
done
ulimit -s unlimited

fail=0
ok()  { echo "  ok: $*"; }
bad() { echo "  FAILED: $*"; fail=1; }
OUTS="edges.out edge_tri.out edgenum.out elvls.out elvls_raw.out nlvls.out"

# prepare <src> <dst>: dst holds the three input files of src and no fesom_ini output
prepare() {
    local src="$1" dst="$2" f
    mkdir -p "$dst"
    for f in nod2d.out elem2d.out aux3d.out; do
        cmp -s "$src/$f" "$dst/$f" || cp "$src/$f" "$dst/$f"
    done
    for f in $OUTS; do rm -f "${dst:?}/$f"; done
    rm -rf "${dst:?}"/dist_*
}

# run_ini <exe> <mesh_dir> <npart> <log>: one fesom_ini run in its own working directory
run_ini() {
    local exe="$1" mesh="$2" np="$3" log="$4" rd
    rd="$WORK_ROOT/run/$(basename "$mesh")_$(basename "$exe")_$np"
    mkdir -p "$rd"
    printf '&paths\nMeshPath="%s/"\n/\n&geometry\n/\n&run_config\n/\n&machine\nn_levels=1\nn_part=%s\n/\n' \
        "$mesh" "$np" > "$rd/namelist.config"
    if ! (cd "$rd" && "$exe" > "$log" 2>&1); then
        bad "$exe on $(basename "$mesh") np=$np (log $log)"
    fi
}

# ndiff <dirA> <dirB> <files...>: echoes the number of differing (or missing) files
ndiff() {
    local a="$1" b="$2" n=0 f; shift 2
    for f in "$@"; do cmp -s "$a/$f" "$b/$f" || n=$((n+1)); done
    echo "$n"
}
dist_files() { (cd "$1/dist_$2" && ls); }

# compare_mesh <name> <src> <exact_vs_src_files> <exact_vs_src_dists> <nparts...>
compare_mesh() {
    local name="$1" src="$2" exact_files="$3" exact_dists="$4"; shift 4
    local dst="$MESH_ROOT/${name}_w" ref="$WORK_ROOT/f2/${name}" np f n files
    echo "== $name ($src) -> $dst"
    prepare "$src" "$dst"
    prepare "$src" "$ref"
    mkdir -p "$WORK_ROOT/log"
    for np in "$@"; do
        run_ini "$INI"    "$dst" "$np" "$WORK_ROOT/log/${name}_f3_$np.log"
        run_ini "$F2_INI" "$ref" "$np" "$WORK_ROOT/log/${name}_f2_$np.log"
    done
    # 1. port vs FESOM2's own build: everything identical
    n=$(ndiff "$dst" "$ref" $OUTS)
    [ "$n" = 0 ] && ok "$name: $OUTS identical to FESOM2's build" \
                 || bad "$name: $n of the edge/level files differ from FESOM2's build"
    for np in "$@"; do
        files=$(dist_files "$ref" "$np")
        n=$(ndiff "$dst/dist_$np" "$ref/dist_$np" $files)
        [ "$n" = 0 ] && ok "$name dist_$np: $(echo $files | wc -w) files identical to FESOM2's build" \
                     || bad "$name dist_$np: $n files differ from FESOM2's build"
    done
    # 2. port vs the files in the original mesh directory
    for f in $exact_files; do
        cmp -s "$dst/${f%%=*}" "$src/${f##*=}" && ok "$name: ${f%%=*} identical to existing ${f##*=}" \
                                               || bad "$name: ${f%%=*} differs from existing ${f##*=}"
    done
    for np in "$@"; do
        files=$(cd "$src/dist_$np" && ls)
        n=$(ndiff "$dst/dist_$np" "$src/dist_$np" $files)
        if [[ " $exact_dists " == *" $np "* ]]; then
            [ "$n" = 0 ] && ok "$name dist_$np: identical to existing" \
                         || bad "$name dist_$np: $n files differ from existing"
        else
            echo "  info: $name dist_$np: $n of $(echo $files | wc -w) files differ from existing (older fesom_ini, not gated)"
        fi
    done
    for f in elvls.out nlvls.out; do
        cmp -s "$dst/$f" "$src/$f" || echo "  info: $name $f differs from existing (older level finder, not gated)"
    done
}

compare_mesh pi "$PI_SRC" \
    "edges.out=edges.out edge_tri.out=edge_tri.out edgenum.out=edgenum.out" \
    "" 2 8
compare_mesh core2 "$CORE2_SRC" \
    "edges.out=edges.out edge_tri.out=edge_tri.out edgenum.out=edgenum.out elvls.out=elvls.out.20260528_regenerated nlvls.out=nlvls.out.20260528_regenerated elvls_raw.out=elvls_raw.out" \
    "2 128" 2 8 128

if [ "$fail" = 0 ]; then echo "run_ini_compare: PASS"; else echo "run_ini_compare: FAIL"; exit 1; fi
