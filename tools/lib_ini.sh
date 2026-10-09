# Shared harness of tools/run_ini_compare.sh and tools/run_ini_mixed.sh (sourced): python3
# with numpy (module), the ok/bad reporters (bad sets fail=1), the fesom_ini output files and
# the namelist of one fesom_ini run.
source /usr/share/Modules/init/bash
module load python3 >/dev/null 2>&1
hash -r

fail=0
ok()  { echo "  ok: $*"; }
bad() { echo "  FAILED: $*"; fail=1; }
OUTS="edges.out edge_tri.out edgenum.out elvls.out elvls_raw.out nlvls.out"

# ini_namelist <file> <mesh_dir> <npart> [use_depthfile]: namelist.config of a fesom_ini run
ini_namelist() {
    local geo=""
    [ $# -ge 4 ] && geo="use_depthfile=\"$4\"\n"
    printf "&paths\nMeshPath=\"%s/\"\n/\n&geometry\n${geo}/\n&run_config\n/\n&machine\nn_levels=1\nn_part=%s\n/\n" \
        "$2" "$3" > "$1"
}
