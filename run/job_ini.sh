#!/bin/bash
#PBS -j oe
#PBS -N fesom_ini
#PBS -l walltime=01:00:00
#PBS -l nodes=150:ppn=24
#PBS -q bmq
#PBS -V

echo $NCPUS
export OMP_WAIT_POLICY=PASSIVE
export CRAY_OMP_CHECK_AFFINITY=TRUE
export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR
pwd

date
aprun -B ../bin/fesom_ini > fesom_ini.out
date

qstat -f $PBS_JOBID
