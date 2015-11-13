#!/bin/bash
#PBS -j oe
#PBS -N fvom
#PBS -l walltime=12:00:00
#PBS -l nodes=12:ppn=24
#PBS -q mppq
#PBS -A hbk00032
#PBS -V

echo $NCPUS
export OMP_WAIT_POLICY=PASSIVE
export CRAY_OMP_CHECK_AFFINITY=TRUE
export OMP_NUM_THREADS=1

cd /home/h/hbkdsido/fvom/bin/

date
aprun -n 288 fvom.x  > fvom.out
date

qstat -f $PBS_JOBID

#msub job.ll
