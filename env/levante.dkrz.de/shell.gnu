# make the contents as shell agnostic as possible so we can include them with bash, zsh and others
export LC_ALL=en_US.UTF-8
export CPU_MODEL=AMD_EPYC_ZEN3

module --force purge

module load git

module load gcc/11.2.0-gcc-11.2.0

# both mpi below work
#module load intel-oneapi-mpi/2021.5.0-gcc-11.2.0
module load openmpi/4.1.2-gcc-11.2.0

# both below work not sure whats the diff?
#module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-gcc-11.2.0
#module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-gcc-11.2.0

module load netcdf-c/4.8.1-gcc-11.2.0
module load netcdf-fortran/4.5.3-gcc-11.2.0

export FC=mpif90 CC=mpicc CXX=mpicxx 

# following is only needed for libblas which is needed by params lib and often provided by lapack
#module load intel-oneapi-mkl/2022.0.1-gcc-11.2.0
# so use the LD_LIBRARY_PATH or other paths like prefix paths etc for cmake
#export LD_LIBRARY_PATH=/sw/spack-levante/intel-oneapi-mkl-2022.0.1-ttdktf/mkl/2022.0.1/lib/intel64:$LD_LIBRARY_PATH spack load intel-oneapi-mkl@2022.0.1%gcc@11.2.0

#other alternative blas
#spack load netlib-lapack@3.9.1%gcc@11.2.0

ulimit -s unlimited # without setting the stack size we get a segfault from the levante netcdf library at runtime
ulimit -c 0 # do not create a coredump after a crash

#OPENMPI specific runtime settings some deviation from ref:https://docs.dkrz.de/doc/levante/running-jobs/runtime-settings.html
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl=self
export OMPI_MCA_osc="pt2pt"
export UCX_IB_ADDR_TYPE=ib_global
# for most runs one may or may not want to disable HCOLL
export OMPI_MCA_coll="^ml,hcoll"
export OMPI_MCA_coll_hcoll_enable="0"
export HCOLL_ENABLE_MCAST_ALL="0"
export HCOLL_MAIN_IB=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
export UCX_UNIFIED_MODE=y
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_io="romio321"
export UCX_HANDLE_ERRORS=bt
