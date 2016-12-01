#!/bin/bash
#SBATCH --job-name=fesom2.0
#SBATCH -p mpp
#SBATCH --ntasks=72
#SBATCH --time=1:00:00

set -x

ulimit -s 1000000

#module purge
#module load intel.compiler
#module load intel.mpi
#module load netcdf/4.4.0_intel

module purge
module load PrgEnv-cray
module swap mvapich2_cce cray-impi
module load netcdf
module load cray-libsci
module load craype-broadwell
module load slurm
module list
module load perftools-base
module load perftools
export LD_LIBRARY_PATH=/opt/cray/papi/default/lib:$LD_LIBRARY_PATH
#pat_build -O apa -f fvom.x
#pat_build -O 000305.apa
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/opt/slurm-1/15.08.8-1/lib/

# determine JOBID
JOBID=`echo $PBS_JOBID |cut -d"." -f1`
export I_MPI_FABRICS=shm:tmi
#export I_MPI_FABRICS=tcp

cd /home/ollie/dsidoren/fesom2.0/bin
date
srun --mpi=pmi2 -l -n 72 ./fvom.x+pat > fvom.out
date
qstat -f $PBS_JOBID

#export EXITSTATUS=$?
#if [ ${EXITSTATUS} -eq 0 ] || [ ${EXITSTATUS} -eq 127 ] ; then
#sbatch job.ll_ollie
#fi

