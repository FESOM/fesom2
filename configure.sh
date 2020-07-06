#!/usr/bin/env bash

set -e

source env.sh # source this from your run script too
mkdir build || true # make sure not to commit this to svn or git
cd build
cmake .. # not required when re-compiling
sed -i -e 's/-lFALSE//g' src/CMakeFiles/fesom.dir/link.txt # workaround for cray
make install -j`nproc --all`
