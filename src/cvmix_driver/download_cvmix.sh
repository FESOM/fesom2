#!/usr/bin/env bash
set -euo pipefail

# Clone and build CVMix-src locally for FESOM2
# Usage:
#   ./get_cvmix.sh [--branch <branch_or_tag>] [--prefix <install_prefix>]
# Defaults:
#   branch: master
#   prefix: $(pwd)/CVMix-install
# Result:
#   Installs CVMix to <prefix>, producing:
#     <prefix>/lib/libcvmix.a
#     <prefix>/include/*.mod
#
# Then configure FESOM with CMake options:
#   -DUSE_EXTERNAL_CVMIX=ON \
#   -DCVMIX_ROOT=<prefix>
# or
#   -DCVMIX_LIB_DIR=<prefix>/lib -DCVMIX_MOD_DIR=<prefix>/include
REPO_URL="https://github.com/CVMix/CVMix-src.git"
BRANCH="master"
PREFIX="CVMix-install"
SRC_DIR="CVMix-src"
BUILD_DIR="CVMix-build"

# to extract the proper src and install dir for cmake, so the are always 
# synchronized with the download script
[[ "${1:-}" == "--print-prefix" ]] && {
  echo "$PREFIX"
  exit 0
}
[[ "${1:-}" == "--print-src_dir" ]] && {
  echo "$PREFIX"
  exit 0
}
[[ "${1:-}" == "--print-build_dir" ]] && {
  echo "$PREFIX"
  exit 0
}
PREFIX=$(pwd)/${PREFIX}
SRC_DIR=$(pwd)/${SRC_DIR}
BUILD_DIR=$(pwd)/${BUILD_DIR}


while [[ $# -gt 0 ]]; do
  case "$1" in
    --branch)
      BRANCH="$2"; shift 2;;
    --prefix)
      PREFIX="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

mkdir -p "$(dirname "$PREFIX")"

if [[ ! -d "$SRC_DIR/.git" ]]; then
  echo "Cloning CVMix-src ($BRANCH) ..."
  git clone --branch "$BRANCH" --depth 1 "$REPO_URL" "$SRC_DIR"
else
  echo "Updating existing CVMix-src ..."
  git -C "$SRC_DIR" fetch --depth 1 origin "$BRANCH"
  git -C "$SRC_DIR" checkout -q "$BRANCH"
  git -C "$SRC_DIR" reset --hard "origin/$BRANCH"
fi

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

cd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_MODULE_DIRECTORY="$PREFIX/include" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      "$SRC_DIR"

cmake --build . --config Release -- -j
cmake --install .

echo
echo "CVMix installed to: $PREFIX"
echo "  lib: $PREFIX/lib"
echo "  mods: $PREFIX/include"
echo
echo "Configure FESOM CMake with:"
echo "  -DUSE_EXTERNAL_CVMIX=ON -DCVMIX_ROOT=$PREFIX"
