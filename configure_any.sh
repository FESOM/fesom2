#!/usr/bin/env bash

# custom build script in use at ECMWF

set -e

LIB=no
while getopts "l" OPT
do
  case "$OPT" in
    l) LIB=yes;;
  esac
done
shift $((OPTIND-1))

source env.sh # source this from your run script too

if [[ ${LIB} = yes ]]; then
  mkdir build.lib || true # build dir for library
  cd build.lib
  cmake -DBUILD_FESOM_AS_LIBRARY=ON -DFESOM_INSTALL_PREFIX=$PWD/../ -DBUILD_SHARED_LIBS=ON -DASYNCHRONOUS_IO_THREADS=OFF -DENABLE_OPENMP=ON .. # not required when re-compiling
  sed -i -e 's/-lFALSE//g' src/CMakeFiles/fesom.dir/link.txt # workaround for the moment on cray
else
  mkdir build || true # build dir for binary
  cd build
  cmake .. # not required when re-compiling
fi
make install -j`nproc --all`
