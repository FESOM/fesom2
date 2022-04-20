# Original source code used to create the POC:

https://gist.github.com/Luthaf/4df78ca52b3caf7fbe0e

Then you can build the example with make command,

# Setting up Hecuba needs 
python3 setup.py install --c_binding=/usr/local/hecuba

# compilation

in the cluster, in order to have the Hecuba includes you first need to 
load the Hecuba modules:

```
module use /apps/HECUBA/modulefiles
module load Hecuba/1.1
```
Note that in the release version, the C++ API is still not present so you 
have to use an special module specially created for such purpose:
```
module use /apps/HECUBA/modulefiles
module load Hecuba/1.0_api
```
```
# set global var where the includes and libs of Hecuba where compiled and installed
export HECUBA_ROOT=/usr/local/hecuba
g++ -o apitest   apitest.cpp  -std=c++11 -I ${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib  -lhfetch   -Wl,-rpath,${HECUBA_ROOT}/lib
```

currently implemented in V4 of the container

```
docker pull jberlindocker/eflows4hpc-wp5:wp5_dev_v4
```

# Running the POC - docker
to run the POC, the following commands should be executed:

```
# inside FESOM2 folder
rm -rf build/
./configure.sh ubuntu
cassandra -R
export CONTACT_NAMES=127.0.0.1
export EXECUTION_NAME=fesom2
cd fesom_hecuba/
mkrun pi test_pi -m docker
cd work_pi/
cp /fesom_hecuba/src/io_hecuba_cpp/*.yaml .
```
The following files should be presented in the working dir:

```
 fcheck_values.csv
 *.yaml (datamodel files)
 job_docker_new
 namelist.config
 namelist.cvmix
 namelist.forcing
 namelist.ice
 namelist.icepack
 namelist.io
 namelist.oce
```

to execute the model run:
```
./job_docker_new
```

To check the tables has been created correctly:

```
# set first the node (in the container is not needed since it points to 127.0.0.1)
export CQLSH_HOST=127.0.0.1
nodetool -h cassandra_node status
cqlsh -e "describe keyspaces"
cqlsh -e "SELECT table_name FROM system_schema.tables WHERE keyspace_name = 'fesom2';"
```
note that in the docker container, cassandro node is 127.0.0.1, to check the content just pick one of the list and execute the following command:

```
! cqlsh -e "select * from fesom2.<table_name>"

 storage_id                           | cluster_id | block_id | payload
--------------------------------------+------------+----------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 1c77c532-4473-8466-2f5a-f7854b82a5a8 |          0 |        0 | 0x000000000000f03f000000000000004000000000000008400000000000001040000000000000144000000000000018400000000000001c4000000000000020400000000000002240000000000000244000000000000026400000000000002840
```

cqlsh -e "SELECT * from fesom2.sst"
```
 lat | ts         | metrics
-----+------------+--------------------------------------
   0 | 1071644672 | e303b960-81dd-0f31-4acf-98871b76c417
```

# Miscellaneous

## install hecuba from a TAG version
``` 
   
    rm -rf hecuba/
    git clone https://github.com/bsc-dd/hecuba.git
    git fetch
    cd hecuba/
    git fetch
    git checkout tags/v1.1 -b v1.1 
    python3 setup.py install --c_binding=/usr/local/
     
   17  git log
   18  history
``` 

## update Hecuba installation (container)
```  
     cd hecuba/
     git pull
     rm -rf build/
     #if it doesnt exist....
     git clone https://github.com/jbeder/yaml-cpp.git
     git branch -a
     cd yaml-cpp/
     rm -rf build/
     mkdir build
     cmake -DYAML_BUILD_SHARED_LIBS=ON
     make install
     cd ..
     ll
     python3 setup.py install --c_binding=/usr/local/
```

------------------------------------- update docker container
update docker container (eflows4hpc-wp5:wp5_dev_v4 is the new version)
```
# get the id of the container currently running asuming the one you made the changes
docker ps
docker commit d7e7e076d01c jberlindocker/eflows4hpc-wp5:wp5_dev_v4
docker push jberlindocker/eflows4hpc-wp5:wp5_dev_v4
```

