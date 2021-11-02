#!/bin/bash
set -e
cd ../

machine="docker"
tests="test_pi test_pi_linfs test_pi_zstar test_pi_partial test_pi_floatice test_pi_visc7 test_pi_zstar"

for test in $tests; do

./configure.sh ubuntu
echo $test
    mkrun pi $test -m $machine
    pwd
    cd work_pi
    chmod +x job_docker_new
    ./job_docker_new
    fcheck .
    cd ../

done

