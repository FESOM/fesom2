#!/usr/bin/env bash
# Source this to load the FESOM3 build environment for the current host.
#   source env.sh <compiler>      # compiler in {intel, gnu}
# Sets FC/CC/CXX (mpi wrappers), FESOM_PLATFORM_STRATEGY, and runtime env.
#
# Anchor toolchain = intel (matches FESOM2 v2.7.3). gnu is the portability build.

_compiler="${1:-intel}"

# Ensure `module` is defined (it is a shell function, not inherited by child
# shells — so a `./configure.sh` subprocess would otherwise lack it).
if ! type module >/dev/null 2>&1; then
    for _init in /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash \
                 /etc/profile.d/modules.sh; do
        # shellcheck disable=SC1090
        [ -f "$_init" ] && source "$_init" && break
    done
fi

if [ -n "$BASH_VERSION" ]; then
    _SOURCE="${BASH_SOURCE[0]}"
elif [ -n "$ZSH_VERSION" ]; then
    _SOURCE=${(%):-%N}
fi
_DIR="$( cd "$( dirname "${_SOURCE}" )" && pwd )"

_host="$(hostname -f)"
if [[ "$_host" =~ ^levante ]] || [[ "$_host" =~ \.lvt\.dkrz\.de$ ]]; then
    export FESOM_PLATFORM_STRATEGY="levante.dkrz.de"
    # Runtime MPI workaround (levante): the OpenMPI 4.1.2 vader BTL's single-copy path
    # (KNEM here — CMA is unavailable under kernel.yama.ptrace_scope=3) CORRUPTS large
    # shared-memory messages on this host: a contiguous nonblocking Isend/Irecv delivers
    # the eager prefix (32 KB) correctly, then fills the tail with STALE data. Root cause
    # (LESSONS L61, tools/run_mpi_knem_stale_repro.sh): the knem registration cache keys
    # on the virtual address and misses glibc heap trims (brk), so a buffer re-allocated
    # at the same address is read from the old pinned pages. This silently
    # broke FESOM3's manual-pack halo exchange (mod_halo core_blk_r) for the CORE2-scale
    # tr_xy element-block exchange (478 KB) — the FIRST live F3 message above the threshold;
    # every pi-mesh gate stayed under it, and FESOM2's MPI_TYPE_INDEXED exchange uses a
    # different vader path that dodges it. Forcing the copy-in/copy-out (two-copy) path is
    # byte-identical and only marginally slower. Verified: with this set the CORE2 multi-rank
    # lifecycle is byte-exact vs FESOM2; without it the tr_xy halo drifts at step 2.
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
else
    echo "env.sh: unknown host '$_host'; only levante.dkrz.de is configured." >&2
    return 1 2>/dev/null || exit 1
fi

_shellfile="${_DIR}/env/${FESOM_PLATFORM_STRATEGY}/shell.${_compiler}"
if [[ ! -e "${_shellfile}" ]]; then
    echo "env.sh: no shell file for compiler '${_compiler}': ${_shellfile}" >&2
    return 1 2>/dev/null || exit 1
fi
echo "env.sh: sourcing ${_shellfile} (platform=${FESOM_PLATFORM_STRATEGY})"
source "${_shellfile}"
