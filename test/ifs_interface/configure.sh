#!/usr/bin/env bash

# custom build script in use at ECMWF

set -e

#cd ../../
source env.sh # source this from your run script too

HERE=$PWD
echo "+ mkdir -p build"
mkdir -p build
echo "+ cd build"
cd build
echo + cmake .. -DCMAKE_INSTALL_PREFIX=$HERE -DBUILD_SHARED_LIBS=ON -DENABLE_IFS_INTERFACE=ON ${CMAKE_ARGS} $@ # not required when re-compiling
cmake .. -DCMAKE_INSTALL_PREFIX=$HERE -DBUILD_SHARED_LIBS=ON -DENABLE_IFS_INTERFACE=ON ${CMAKE_ARGS} $@ # not required when re-compiling
echo "+ make install -j`nproc --all`"
make install -j`nproc --all`

