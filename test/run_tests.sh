#!/bin/bash
set -e
cd ../

machine="docker"
tests="test_pi test_souf test_pi_linfs test_pi_zstar test_pi_partial test_pi_floatice test_pi_visc7"

./configure.sh ubuntu

for test in $tests; do

echo $test
    mkrun pi $test -m $machine
    pwd
    cd work_pi
    chmod +x job_docker_new
    ./job_docker_new
    echo "This was ${test}"
    fcheck .
    cd ../

done

othertest="test_lib_compiles"

for test in $othertest; do
	
    echo $othertest
    ./test/ifs_interface/configure_lib.sh -l

    FILE=./lib/libfesom.a
    if [ -f "$FILE" ]; then
	echo "$FILE compiled and linked."
    else
	echo "$FILE does not exist."
    fi
done
