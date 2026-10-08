#!/usr/bin/env bash
# LESSONS L61 evidence: Open MPI 4.1.2 vader/knem registration cache vs glibc heap trimming.
# Builds tools/mpi_knem_stale_repro.F90 (correct, contiguous, fully waited Isend/Irecv) and runs it
# at np 2 on the login node under several transport settings. Expected on Levante (2026-10-08):
#   ob1+vader, knem (the default)          + heap trim   STALE (788 of 800 receives)
#   ob1+vader, knem                        no heap trim  ok
#   ob1+vader, single_copy none / emulated + heap trim   ok
#   ob1+vader, knem, mpi_leave_pinned 0    + heap trim   ok   (no knem registration cache)
#   ucx (mm,knem / mm,cma, rndv 1 KB)      + heap trim   ok
# "heap trim" = MALLOC_TRIM_THRESHOLD_=0 MALLOC_TOP_PAD_=0: glibc free() returns the heap top to the
# kernel with brk on every call, the next allocate re-grows it at the SAME virtual address with
# fresh pages. usage: tools/run_mpi_knem_stale_repro.sh [build_dir]   (default: a mktemp dir)
set -uo pipefail
OMPI=/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLD=${1:-$(mktemp -d)}
mkdir -p "$BLD"
"$OMPI/mpif90" -O3 "$SRC/mpi_knem_stale_repro.F90" -o "$BLD/mpi_knem_stale_repro" || exit 2
OB1="--mca pml ob1 --mca btl self,vader --oversubscribe"
UCX="--mca pml ucx --mca pml_ucx_tls any --mca pml_ucx_devices any --oversubscribe"
TRIM="-x MALLOC_TRIM_THRESHOLD_=0 -x MALLOC_TOP_PAD_=0"
run() {   # run <label> <mpiexec args ...>
    local label=$1; shift
    printf '%-46s ' "$label"
    "$OMPI/mpiexec" "$@" -np 2 "$BLD/mpi_knem_stale_repro" 2>&1 | grep 'mpi_knem_stale_repro:' | sed 's/.*: //'
}
run "ob1 vader knem (default), heap trim"      $OB1 $TRIM
run "ob1 vader knem (default), no heap trim"   $OB1
run "ob1 vader single_copy none, heap trim"    $OB1 $TRIM --mca btl_vader_single_copy_mechanism none
run "ob1 vader single_copy emulated, heap trim" $OB1 $TRIM --mca btl_vader_single_copy_mechanism emulated
run "ob1 vader knem, mpi_leave_pinned 0, trim" $OB1 $TRIM --mca mpi_leave_pinned 0
run "ucx mm,knem rndv 1KB, heap trim"          $UCX $TRIM -x UCX_TLS=self,mm,knem -x UCX_RNDV_THRESH=1024
run "ucx mm,cma rndv 1KB, heap trim"           $UCX $TRIM -x UCX_TLS=self,mm,cma -x UCX_RNDV_THRESH=1024
