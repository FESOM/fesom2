#!/bin/bash
#SBATCH --job-name=fesom2.0
#SBATCH --ntasks=36
#SBATCH --time=00:30:00
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out
#SBATCH --partition=batch

module load CMake Intel IntelMPI imkl netCDF netCDF-Fortran

set -x

ulimit -s unlimited

# determine JOBID
JOBID=`echo $SLURM_JOB_ID |cut -d"." -f1`

ln -s ../bin/fesom.x .           # cp -n ../bin/fesom.x
cp -n ../config/namelist.config  .
cp -n ../config/namelist.forcing .
cp -n ../config/namelist.oce     .
cp -n ../config/namelist.ice     .

date
srun --mpi=pmi2 ./fesom.x > "fesom2.0.out"
date

#qstat -f $PBS_JOBID
#export EXITSTATUS=$?
#if [ ${EXITSTATUS} -eq 0 ] || [ ${EXITSTATUS} -eq 127 ] ; then
#sbatch job_ollie
#fi
