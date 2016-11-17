#!/bin/sh

cd lib/metis-5.1.0
make
cd ../../

cd lib/parms
make
cd ../../

cd src
make fesom_ini
make clean
make fesom
