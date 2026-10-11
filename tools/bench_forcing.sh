#!/usr/bin/env bash
# Time the atmospheric forcing path on the pi mesh across rank and thread counts.
#
# Reuses the namelists of the integration_pi_mpi2 CTest (configure with
# -DBUILD_TESTING=ON and run that test once), partitions pi for every rank
# count asked for, and runs fesom.x once per (mode, ranks, threads). The
# profiler sections of the forcing update (fesom.stats, Max over ranks) and
# "runtime ice read forc" are collected into one CSV.
#
#   tools/bench_forcing.sh -b build -o /tmp/bench [-n "1 2 4 8 16 32"] [-t "1 2 4"]
#                          [-m "bcast scatter"] [-s 96] [-d 1] [-f FORCING_DIR] [-F NC_FREQ] [-r REPEAT]
#
# -f points at a forcing directory to use instead of the bundled CORE2 one,
# e.g. a longer series from tools/make_synthetic_forcing.py. FESOM ignores the
# units attribute of the time axis, so pass -F 24 for an axis in hours (the
# bundled CORE2 axis, 0 6 12 .., is read as days with the default nm_nc_freq=1).
set -euo pipefail

SRC_DIR=$(cd "$(dirname "$0")/.." && pwd)
BUILD_DIR=build
OUT_DIR=
NPS="1 2 4 8 16 32"
THREADS="1"
MODES="bcast scatter"
STEPS=96
DAYS=1
FORCING_DIR=
REPEAT=1
NC_FREQ=
MAX_CORES=$(nproc)

while getopts "b:o:n:t:m:s:d:f:F:r:c:h" opt; do
  case $opt in
    b) BUILD_DIR=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    n) NPS=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    m) MODES=$OPTARG ;;
    s) STEPS=$OPTARG ;;
    d) DAYS=$OPTARG ;;
    f) FORCING_DIR=$(cd "$OPTARG" && pwd) ;;
    F) NC_FREQ=$OPTARG ;;
    r) REPEAT=$OPTARG ;;
    c) MAX_CORES=$OPTARG ;;
    *) sed -n '2,18p' "$0"; exit 1 ;;
  esac
done
[[ -n $OUT_DIR ]] || { echo "need -o OUT_DIR" >&2; exit 1; }

BUILD_DIR=$(cd "$BUILD_DIR" && pwd)
FESOM=$BUILD_DIR/bin/fesom.x
MESHPART=$BUILD_DIR/bin/fesom_meshpart
TEMPLATE=$BUILD_DIR/tests/integration/integration_pi_mpi2
[[ -x $FESOM ]] || { echo "no $FESOM" >&2; exit 1; }
[[ -f $TEMPLATE/namelist.config ]] || { echo "run 'ctest -R integration_pi_mpi2' in $BUILD_DIR first" >&2; exit 1; }
mkdir -p "$OUT_DIR"
OUT_DIR=$(cd "$OUT_DIR" && pwd)

# setnml FILE KEY VALUE: replace "KEY = ..." keeping the rest of the file
setnml() {
  sed -i -E "s|^([[:space:]]*)$2[[:space:]]*=[^!]*|\1$2=$3 |" "$1"
}

# --- mesh with all the partitions we need -----------------------------------
MESH=$OUT_DIR/mesh_pi
mkdir -p "$MESH"
for f in "$SRC_DIR"/tests/data/MESHES/pi/*; do
  name=$(basename "$f")
  [[ $name == dist_* ]] && continue
  ln -sfn "$f" "$MESH/$name"
done
for np in $NPS; do
  [[ -f $MESH/dist_$np/rpart.out ]] && continue
  if [[ -f $SRC_DIR/tests/data/MESHES/pi/dist_$np/rpart.out ]]; then
    ln -sfn "$SRC_DIR/tests/data/MESHES/pi/dist_$np" "$MESH/dist_$np"
    continue
  fi
  pdir=$OUT_DIR/partition_$np
  mkdir -p "$pdir" "$MESH/dist_$np"
  cp "$SRC_DIR/config/namelist.config" "$pdir/"
  setnml "$pdir/namelist.config" MeshPath "'$MESH/'"
  setnml "$pdir/namelist.config" ResultPath "'$pdir/'"
  setnml "$pdir/namelist.config" n_levels 1
  setnml "$pdir/namelist.config" n_part "$np"
  (cd "$pdir" && "$MESHPART" > partition.log 2>&1) || { echo "partitioning for $np failed, see $pdir/partition.log" >&2; exit 1; }
done

# --- runs ---------------------------------------------------------------------
CSV=$OUT_DIR/bench_forcing.csv
SECTIONS="update_atm_forcing forc_sbc_do sbc_getcoeff sbc_read_distribute sbc_ioread sbc_bcast sbc_scatterv sbc_fill_land sbc_interp sbc_time_interp forc_wind_drag forc_ncar_fluxes"
if [[ ! -f $CSV ]]; then
  echo "mode,np,threads,rep,steps,wall,read_forc_max,$(echo $SECTIONS | tr ' ' ',')" > "$CSV"
fi

export OMPI_MCA_rmaps_base_oversubscribe=1
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe

for mode in $MODES; do
for np in $NPS; do
for nt in $THREADS; do
  (( np * nt <= MAX_CORES )) || continue
  for rep in $(seq 1 "$REPEAT"); do
    run=$OUT_DIR/run_${mode}_np${np}_t${nt}_r${rep}
    rm -rf "$run"; mkdir -p "$run/results"
    cp "$TEMPLATE"/namelist.* "$run/"
    setnml "$run/namelist.config" MeshPath "'$MESH/'"
    setnml "$run/namelist.config" ResultPath "'$run/results/'"
    setnml "$run/namelist.config" step_per_day "$STEPS"
    setnml "$run/namelist.config" run_length "$DAYS"
    setnml "$run/namelist.config" run_length_unit "'d'"
    setnml "$run/namelist.config" restart_length "$DAYS"
    setnml "$run/namelist.config" logfile_outfreq "$STEPS"
    [[ -n $FORCING_DIR ]] && setnml "$run/namelist.config" ForcingDataPath "'$FORCING_DIR/'"
    [[ -n $NC_FREQ ]] && setnml "$run/namelist.forcing" nm_nc_freq "$NC_FREQ"
    if grep -q "forcing_distribution" "$run/namelist.forcing"; then
      setnml "$run/namelist.forcing" forcing_distribution "'$mode'"
    else
      sed -i -E "s|^([[:space:]]*)&nam_sbc|&\n   forcing_distribution='$mode'|" "$run/namelist.forcing"
    fi
    # a fresh clock: the template's may have advanced past the start
    printf '0.0 1 1948\n0.0 1 1948\n' > "$run/results/fesom.clock"

    t0=$(date +%s.%N)
    if ! (cd "$run" && OMP_NUM_THREADS=$nt OMP_PROC_BIND=false \
          bash -c "ulimit -s unlimited 2>/dev/null; exec mpiexec -n $np --bind-to none '$FESOM'" > fesom.out 2> fesom.err); then
      echo "run failed: $run" >&2; continue
    fi
    t1=$(date +%s.%N)
    wall=$(echo "$t1 - $t0" | bc)
    readforc=$(awk -F: '/runtime ice read forc/ {split($2,a," "); print a[3]}' "$run/fesom.out" | tail -1)
    row="$mode,$np,$nt,$rep,$((STEPS*DAYS)),$wall,${readforc:-}"
    for s in $SECTIONS; do
      # columns: name mean min max calls ...; Max over ranks is the critical path
      v=$(awk -v s=">$s" '$1==s {print $4; exit}' "$run/results/fesom.stats")
      row="$row,${v:-0}"
    done
    echo "$row" | tee -a "$CSV"
  done
done
done
done
echo "results in $CSV"
