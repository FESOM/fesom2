#!/usr/bin/env bash

set -e

MYPATH=$PWD
cd ../
source env.sh # source this from your run script too
cd $MYPATH

mkdir build || true # make sure not to commit this to svn or git
cd build
cmake .. # not required when re-compiling
make
ln -s ../mesh_part/build/fesom_ini ../../bin/fesom_ini.x


