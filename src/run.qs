#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=112gb
####PBS -l select=1:ncpus=64:mem=192gb 
####PBS -l select=1:ncpus=32:mem=192gb 
####PBS -l select=1:ncpus=8:mem=56gb 
####PBS -l select=1:ncpus=8:mem=112gb 
####PBS -l select=1:ncpus=32
#PBS -N FVOM
#PBS -o /uv/user/kkorchuk/fvom/results
#PBS -e /uv/user/kkorchuk/fvom/results
eval `$MODULESHOME/bin/modulesinit`
#module load intel.compiler
module load mpt
###module load intel.compiler
date

###export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uv/soft/intel/composer_xe_2011_sp1.9.293/composer_xe_2011_sp1.9.293/mkl/lib/intel64/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uv/soft/intel/composer_xe_2013.1.117/composer_xe_2013/compiler/lib/intel64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uv/soft/intel/composer_xe_2013.1.117/composer_xe_2013/mkl/lib/intel64/
cd $PBS_O_WORKDIR
date
echo "***********CONFIG*************"
cat ../config/namelist.config
echo "***********OCE*************"
cat ../config/namelist.oce
echo "***********ICE*************"
cat ../config/namelist.ice
echo "***********ITERATIONS*************"
python reset_clock.py
mpirun -np $NCPUS ./fvom.x
date
qstat -f $PBS_JOBID
