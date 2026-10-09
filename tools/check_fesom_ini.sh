#!/usr/bin/env bash
# fesom_ini paths without another test (ctest fesom_ini_*, hybrid-mesh review fixes): runs
# fesom_ini on a COPY of a quad mesh (nod2d.out, elem2d.out, depth.out of a NEMO cut-out)
# in <work_dir>, never in the mesh directory itself.
#   aux3d     : the FESOM2 level finder on quads (aux3d.out = the first nl+1+nod2D lines of
#               depth.out, use_depthfile="aux3d"), 2 parts; checks that dist_2 is written
#               and that nlvls.out == max of elvls.out over the elements around each node
#               (4 vertices of a quad) -> prints "check_fesom_ini: aux3d OK"
#   nodepth   : use_depthfile="nemo" without depth.out          -> fesom_ini error message
#   onelem    : use_depthfile="nemo" with use_depthonelem=.true. -> fesom_ini error message
#   count1    : a depth.out layer count of 1 (the model needs 2)  -> fesom_ini error message
#   nonconvex : the 1st vertex of element 1 (a quad) moved inside the triangle of its
#               other three vertices (reflex corner at vertex 1)  -> fesom_ini error message
# The error cases print "check_fesom_ini: dist_2 written" if fesom_ini went on (the ctest
# FAIL regex). Standard tools only (awk), no numpy.
#
# Usage: tools/check_fesom_ini.sh <fesom_ini> <quad_mesh_dir> <work_dir> <case>
set -uo pipefail
INI=$1; SRC=$2; WD=$3; CASE=$4
M=$WD/mesh
rm -rf "$WD"; mkdir -p "$M" "$WD/run"
cp "$SRC/nod2d.out" "$SRC/elem2d.out" "$M/"
nl=$(head -1 "$SRC/depth.out" | tr -d ' ')
n2d=$(head -1 "$SRC/nod2d.out" | tr -d ' ')
depthfile=nemo; onelem=.false.
case $CASE in
    aux3d)     head -n $((1 + nl + n2d)) "$SRC/depth.out" > "$M/aux3d.out"; depthfile=aux3d ;;
    nodepth)   ;;
    onelem)    cp "$SRC/depth.out" "$M/"; onelem=.true. ;;
    count1)    awk -v i=$((2 + nl + n2d)) 'NR == i {print "        1 "; next} {print}' \
                   "$SRC/depth.out" > "$M/depth.out" ;;
    nonconvex) cp "$SRC/depth.out" "$M/"
               read -r v1 v2 v3 v4 < <(sed -n 2p "$SRC/elem2d.out")
               awk -v a=$v1 -v b=$v2 -v c=$v3 -v d=$v4 '
                   NR == FNR { if (FNR > 1) { x[$1] = $2; y[$1] = $3 }; next }
                   FNR > 1 && $1 == a { $2 = sprintf("%.10f", (x[b] + x[c] + x[d]) / 3)
                                        $3 = sprintf("%.10f", (y[b] + y[c] + y[d]) / 3) }
                   { print }' "$SRC/nod2d.out" "$SRC/nod2d.out" > "$M/nod2d.out" ;;
    *) echo "check_fesom_ini: unknown case $CASE"; exit 2 ;;
esac
printf '&paths\nMeshPath="%s/"\n/\n&geometry\nuse_depthfile="%s"\nuse_depthonelem=%s\n/\n&run_config\n/\n&machine\nn_levels=1\nn_part=2\n/\n' \
    "$M" "$depthfile" "$onelem" > "$WD/run/namelist.config"
(cd "$WD/run" && "$INI") 2>&1
[ -s "$M/dist_2/rpart.out" ] && written=1 || written=0
if [ "$CASE" != aux3d ]; then
    [ $written = 1 ] && echo "check_fesom_ini: dist_2 written"
    exit 0
fi
[ $written = 1 ] || { echo "check_fesom_ini: aux3d FAILED: no dist_2/rpart.out"; exit 1; }
# nlvls.out == max of elvls.out over the elements around each node (3 or 4 vertices)
if awk 'FILENAME == ARGV[1] { if (FNR > 1) { n = (NF == 4 && $4 != $1) ? 4 : 3
                                              for (k = 1; k <= n; k++) v[FNR-1, k] = $k
                                              nv[FNR-1] = n }; next }
        FILENAME == ARGV[2] { e = FNR; for (k = 1; k <= nv[e]; k++)
                                  if ($1 > mx[v[e, k]]) mx[v[e, k]] = $1; ne = e; next }
        { nn++; if ($1 != mx[FNR]) bad++ }
        END { printf "check_fesom_ini: %d elements, %d nodes, %d node levels != max over the elements\n",
                     ne, nn, bad; exit (bad > 0 || nn == 0) }' \
       "$M/elem2d.out" "$M/elvls.out" "$M/nlvls.out"; then
    echo "check_fesom_ini: aux3d OK"
else
    echo "check_fesom_ini: aux3d FAILED: nlvls.out"
    exit 1
fi
