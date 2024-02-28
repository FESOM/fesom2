#!/usr/bin/env bash

set -e
HERE=$PWD

SOURCE_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
BUILD_DIR=${BUILD_DIR:-build}

source env.sh # source this from your run script too
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}
cmake .. -DCMAKE_INSTALL_PREFIX=$HERE -DCMAKE_BUILD_TYPE=Debug ${CMAKE_ARGS} $@
  # not required when re-compiling
  # additional cmake arguments can be passed to configure.sh
  # this also includes fesom specific options in CMakeLists, can be used as -DFESOM_COUPLED=ON
make install -j`nproc --all`
