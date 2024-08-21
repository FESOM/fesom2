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

for test in $tests; do

    ./configure.sh ubuntu
    echo $test
    mkrun pi $test -m $machine
    cd work_pi_one_timestep
    chmod +x job_docker_new
    ./job_docker_new_one_timestep
    fcheck .
    cd ../

done
