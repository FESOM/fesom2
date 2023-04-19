#!/usr/bin/env bash

set -e

source env.sh # source this from your run script too
mkdir build || true # make sure not to commit this to svn or git
cd build
cmake .. $@ -DCMAKE_BUILD_TYPE=Debug # not required when re-compiling
                                     # additional cmake arguments can be passed to configure.sh
				     # this also includes fesom specific options in CMakeLists, can be used as -DFESOM_COUPLED=ON
make install -j`nproc --all`
