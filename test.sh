#!/bin/bash
# Run simples FESOM2 test in a container.
#
# With singularity on ollie
#
# module load singularity/3.9.1
# cd fesom2
# singularity exec /home/ollie/nkolduno/SINGULARITY/fesom_refactoring.sif ./test.sh
#
# With docker on Linux/Mac
# docker run -it -v "$(pwd)"/fesom2:/fesom/fesom2 koldunovn/fesom2_test:refactoring /bin/bash
# cd fesom2
# ./test.sh
#

set -e

machine="docker"
tests="test_pi"

function red_echo() {
    echo -e "\e[31m$1\e[0m"
}

function green_echo() {
    echo -e "\e[32m$1\e[0m"
}

for test in $tests; do

    ./configure.sh ubuntu
    echo $test
    mkrun pi $test -m $machine
    cd work_pi
    chmod +x job_docker_new
    ./job_docker_new
    fcheck . && green_echo "Test $test passed" || red_echo "Test $test failed"
    cd ../
done

echo "Performing UGRID Compliance Test"
pip install ugrid-checks xarray
# Needed to combine the mesh diag and output
python3 -c "import xarray as xr; ds1 = xr.open_dataset('./test/output_pi/sst.fesom.1948.nc'); ds2 = xr.open_dataset('./test/output_pi/fesom.mesh.diag.nc'); xr.merge([ds1, ds2]).to_netcdf('./test/output_pi/merged_ugrid_check.nc')"
ugrid-checker ./test/output_pi/merged_ugrid_check.nc && green_echo "UGRID Compliance Test Passed" || red_echo "UGRID Compliance Test Failed"
