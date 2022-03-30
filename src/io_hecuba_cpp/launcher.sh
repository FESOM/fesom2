# Launch Cassandra
cassandra -R &
export CONTACT_NAMES=127.0.0.1
export CQLSH_HOST=127.0.0.1
export EXECUTION_NAME=fesom2
cd ..
cd ..
# prepare test run
mkrun pi test_pi -m docker
cd work_pi/
cp /fesom_hecuba/src/io_hecuba_cpp/fesom_datamodel.yaml .
./job_docker_new