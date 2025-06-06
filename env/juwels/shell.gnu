#!/usr/bin/bash
# Use this file to source the environment in your
# preprocessing or postprocessing scripts
module --force purge

#_______________________________________________________________________________
# module load Stages/2023
# module load Stages/2024
module load Stages/2025
 
#_______________________________________________________________________________
# use GCC + openMPI
module load GCC OpenMPI
export FC=gfortran
export F77=gfortran
export MPIFC=gfortran
export FCFLAGS=-free
export CC=cc
export CXX=c++

##_______________________________________________________________________________
## use INTEL + MPI
#module load Intel ParaStationMPI
## module load Intel ParaStationMPI/5.10.0-1-mt # -mt ... Multi-threaded Version
#export FC=mpifort
#export F77=mpifort
#export MPIFC=mpifort
#export FCFLAGS=-free
#export CC=mpicc
#export CXX=mpic++

#_______________________________________________________________________________
# Intel oneAPI Math Kernel Library
module load imkl

#_______________________________________________________________________________
# NetCDF (network Common Data Form) software libraries
module load netCDF netCDF-Fortran

#_______________________________________________________________________________
# Others libraries
module load CMake
module load git

#_______________________________________________________________________________
#module load Python
#module load Perl
#module load CDO
#module unload ecCodes
#module load libaec FFTW cURL

#_______________________________________________________________________________
# show list of now loaded module environment
module list
