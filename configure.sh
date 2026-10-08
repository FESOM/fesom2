#!/usr/bin/env bash
# FESOM3 configure + build helper.
#
#   ./configure.sh [--compiler intel|gnu] [--precision dp|sp]
#                  [--debug] [--clean] [--build] [--jobs N] [-- <extra cmake args>]
#
# Produces build_<compiler>_<precision>/ (Release) or build_<compiler>_<precision>_debug/
# (--debug, kept separate so it never clobbers the Release anchor). Default = anchor
# (intel, dp, Release).
set -euo pipefail

compiler="intel"
precision="dp"
build_type="Release"
do_clean=0
do_build=0
jobs="$(nproc 2>/dev/null || echo 4)"
extra_cmake=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --compiler)  compiler="$2"; shift 2 ;;
        --precision) precision="$2"; shift 2 ;;
        --debug)     build_type="Debug"; shift ;;
        --clean)     do_clean=1; shift ;;
        --build)     do_build=1; shift ;;
        --jobs)      jobs="$2"; shift 2 ;;
        --)          shift; extra_cmake=("$@"); break ;;
        *) echo "configure.sh: unknown arg '$1'" >&2; exit 1 ;;
    esac
done

case "$precision" in
    dp) prec_flag="" ;;
    sp) prec_flag="-DUSE_SINGLE_PRECISION=ON" ;;
    hp) prec_flag="-DUSE_HALF_PRECISION=ON" ;;
    *) echo "configure.sh: precision must be dp|sp|hp" >&2; exit 1 ;;
esac

src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
build_dir="${src_dir}/build_${compiler}_${precision}"
# Debug builds get their OWN dir so they never clobber the Release anchor build.
# A Debug binary sitting in the Release dir silently breaks the byte-gates: at -O0
# without -no-prec-div it diverges from the FESOM2 -O3 oracle on every FP field
# (even `area` at ~3.7e-4 — see LESSONS L10/L14).
if [[ "${build_type}" == "Debug" ]]; then
    build_dir="${build_dir}_debug"
fi

# Load the toolchain (sets FC=mpif90, FESOM_PLATFORM_STRATEGY, runtime env).
# shellcheck disable=SC1090
source "${src_dir}/env.sh" "${compiler}"

if [[ ${do_clean} -eq 1 ]]; then
    echo "configure.sh: removing ${build_dir}"
    rm -rf "${build_dir}"
fi
mkdir -p "${build_dir}"

echo "configure.sh: cmake -> ${build_dir} (compiler=${compiler} precision=${precision} type=${build_type})"
cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DCMAKE_Fortran_COMPILER="${FC}" \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    ${prec_flag} \
    "${extra_cmake[@]}"

if [[ ${do_build} -eq 1 ]]; then
    echo "configure.sh: building with ${jobs} jobs"
    cmake --build "${build_dir}" -j "${jobs}"
fi
