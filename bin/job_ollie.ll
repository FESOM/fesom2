#!/bin/bash
#SBATCH --job-name=fesom2.0
#SBATCH -p mpp
#SBATCH --ntasks=72
#SBATCH --time=0:30:00
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out

set -x

ulimit -s 1000000

module purge
module load intel.compiler
module load intel.mpi
module load netcdf/4.4.0_intel

# determine JOBID
JOBID=`echo $SLURM_JOB_ID |cut -d"." -f1`

#export I_MPI_FABRICS=shm:tmi

cd /home/ollie/pscholz/trunk/bin
date
#### srun --mpi=pmi2 --ntasks=72 ./fvom.x > "fvom-${SLURM_JOB_ID}.out"
srun --mpi=pmi2 --ntasks=72 ./fvom.x > "fvom.out"
date
# qstat -f $PBS_JOBID

#export EXITSTATUS=$?
#if [ ${EXITSTATUS} -eq 0 ] || [ ${EXITSTATUS} -eq 127 ] ; then
#sbatch job.ll_ollie
#fi

