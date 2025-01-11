.. _chap_getting_started:

Getting Started with FESOM2
***************************

This chapter describes several ways of getting started with FESOM2. First we show a minimum set of comands that will lead to a working setup on systems where FESOM2 is used activelly. We also have instructions for Docker/Singularity and Ubuntu.

TL;DR version for supported HPC systems
=======================================

Supported systems are: generic ``ubuntu``, ``albedo`` at AWI, ``levante`` at DKRZ, ``JURECA`` at JSC, ``HLRN``, ``Hazel Hen``, ``MareNostrum 4`` at BSC. During configuration the system will be recognised and apropriate environment variables and compiler options should be used.

::

    git clone https://github.com/FESOM/fesom2.git
    cd fesom2
    bash -l ./configure.sh

Create file ``fesom.clock`` in the output directory with the following content (if you plan to run with COREII forcing):

::

    0 1 1958
    0 1 1958

after that, one has to adjust the run script for the target system and run it:

::

    cd work
    sbatch job_albedo


Detailed steps of compiling and runing the code
===============================================

The following section assumes you are located on one of the supported HPC systems. To install FESOM2 on your local machine we recoment to use `Docker based installation`_ and read about `Necessary Ubuntu packages`_ if you decide not to use Docker. 

First thing is to checkout FESOM2 code from the repository. The code is developed in open repository on GitHub_. 

.. _GitHub: https://github.com/FESOM/fesom2/

Build model executable with Cmake
---------------------------------

Clone the GitHub repository with a git command:

::

    git clone https://github.com/FESOM/fesom2.git


The repository contains model code and two additional libraries: `Metis` (domain partitioner) and `Parms` (solver), necessary to run FESOM2. To build FESOM2 executable one have to compile Parms library and the code of the model (`src` folder). In order to build executable that is used for model domain partitioning (distribution of the model mesh between CPUs) one have to compile `Metis` library and also some code located in the src directory (see :ref:`partitioning`). Building of the model executable and the partitioner is usually done automatically with the use of CMake. If you going to build the code not on one of the supported platforms (ollie, DKRZ, HLRN, HAZELHEN, and BSC, general Ubuntu), you might need to do some (usually small) modifications described in `Adding new platform for compilation`_ section.

Change to the `fesom2` folder and execute:

::

    cd fesom2
    git checkout refactoring

As a good practice, if one wants to make modifications to the source code or any of the files, it is advisable to create a branch from refactoring:

::

    git checkout -b <my branch> refactoring

After confirming that the right FESOM2 branch is being used, compile the model with:

::
    
    bash -l ./configure.sh

In the best case scenario, your platform will be recognized and the Parms library and model executable will be built and copied to the bin directory. If something went wrong have a look at Troubleshooting_ section.

If you would like to select platform manually (which is necessary in the case of Ubuntu, for example), type:

::

    bash -l ./configure.sh ubuntu


Data and mesh files
-------------------

The FESOM2 repository contains only very small example meshes and data (in the ``test`` directory, see the note below). However, if you want to run realistic simulations, you ether have to have them on your system, or download an archive with sample data. There is a chance that your system already has some of the necesseary files, you can check it in the ``setups/paths.yml`` file. If not, the easiest way to start is to download example set from `DKRZ cloud`_  (12 Gb) by executing:

::

    curl https://swift.dkrz.de/v1/dkrz_035d8f6ff058403bb42f8302e6badfbc/FESOM2.0_tutorial/FESOM2_one_year_input.tar > FESOM2_one_year_input.tar

and untar:

::

    tar -xvf FESOM2_one_year_input.tar

You will have a folder named ``FESOM2_one_year_input`` that contains all the data you need to do initial run of the model. The `mesh` directory contains two meshes: ``pi`` and ``core2``. The ``pi`` mesh is very small global FESOM2 mesh, that can run relativelly fast even on a laptop. The ``CORE`` mesh is our 1 degree equivalent mesh and is used in many tuning and testing studies. Mesh folders already include several prepared partitionings (``dist_`` folders), so you don't have to worry about partitioning during your first steps with FESOM.

The ``input`` folder contains files with initial conditions (``phc3.0``) and atmospheric forcing (``JRA55``) for one year (1958).

.. note:: You can find more standard FESOM2 meshes in https://gitlab.awi.de/fesom . Download instructions are available in each mesh repository.


.. _DKRZ cloud: https://swiftbrowser.dkrz.de/download/FESOM2.0_tutorial/FESOM2_one_year_input.tar

.. note::  The FESOM2 distribution contains minimal set of data to run the model in the ``test`` directory, namelly ``pi`` and ``soufflet`` (channel) meshes, WOA13 initial conditions and CORE2 forcing data for one day. Those are mainly used for testing, and require a bit more involved modification of namelists. For more details see instructions on `Docker based installation`_. 


Preparing the run
------------------

You have to do several basic things in order to prepare the run. 

First, be aware of the files you need to modify according to your run configurations. Normally, those are:

- ``namelist.config``: inside of the ``config`` folder. In this file you can set several configurations, such as the path to your mesh, climatology and results, as well as run length, units and start year of your run. 

- ``namelist.forcing``: inside of the ``config`` folder. In this file you can set the path to your forcing files.

- ``job_<name-of-hpc>``: inside of the ``work`` folder. In this file you can set other important configurations, such as the time, tasks and tasks per node you allocate to your run.

The exact changes necessary to those file are indicated later in this documentation. Before doing so, create a directory to store your output. Usually, it is created in the model root directory:

::

    mkdir results

You might make a link to some other directory located on the part of the system where you have a lot of storage. 

In your results directory, create a file named ``fesom.clock`` (NOTE: if you change ``runid`` in ``namelist.config`` to something like ``runid=mygreatrun``, the file will be named ``mygreatrun.clock``).           

Inside the file you have to put two identical lines:

::

    0 1 1958
    0 1 1958

This is initial date of the model run, or the time of the `cold start` of your model. In case you want to start your run with a specific forcing from a specific year, substitute 1958 to the desired year. More detailed explanation of the clock file will be given in the `The clock file`_ section. 

The next step is to make some changes in the model configuration. All runtime options can be set in the namelists that are located in the config directory:

::

    cd ../config/

As mentioned before, in this directory, you will normally have to change two files: ``namelist.config`` and ``namelist.forcing``. Both of these files ask for paths to initial conditions. Normally, these paths can be found under ``./setups/paths.yml``. 

Changing namelist.config
========================

In ``namelist.config``, the options that you might want to change for your first FESOM2 run are:

- ``run_length``: length of the model run in run_length_unit (see below). 

- ``run_length_unit``: units of the run_length. Can be ``y`` (year), ``m`` (month), ``d`` (days), ``s`` (model steps).

.. note:: you might need to adjust the run time to the length of your run. In some setups and/or for some machines, if you set ``run_length`` to 10 and ``run_length_unit`` to ``y``, for example, the run time needs to be enough for a 10-year run at once.  

- ``yearnew``: define the same as the year in your ``fesom.clock``;

- ``MeshPath``: path to the mesh you would like to use (e.g. ``/youdir/FESOM2_one_year_input/mesh/pi/``, slash at the end is important!);

- ``ClimateDataPath``: path to the folder with the file with model temperature and salinity initial conditions (e.g. ``/youdir/FESOM2_one_year_input/input/phc3.0/``). The name of the file with initial conditions is defined in `namelist.oce`, but during first runs you probably don't want to change it;

- ``ResultPath``: path to your results folder. The output of the model will be stored there.

More detailed explanation of options in the ``namelist.config`` is in the section :ref:`chap_general_configuration`.

Changing namelist.forcing
=========================

In ``namelist.forcing``, the options you need to change for your first FESOM2 run depends on the forcing you decide to use to initialize your experiment. Please note that the year you initialize your experiment with needs to be included in the forcing data files.

In section ``&nam_sbc``, change the path of all the files to the path to the forcing you have chosen. For example, if you want to initialize your experiment with JRA55 forcing on ``levante``, the path to each fiel will be:

::

'/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/<name_of_variable>'

More detailed explanation of options in the ``namelist.forcing`` is in the section :ref:`chap_forcing_configuration`.

Running the model
-----------------

Change to the ``work`` directory. You should find several batch scripts that are used to submit model jobs to different HPC machines. The scripts also link ``fesom.x`` executable to the ``work`` directory and copy namelists with configurations from config folder.

.. note::
   Model executable, namelists and job script will be located in the same directory (usually ``work``).

If you are working on AWI's ``albedo`` supercomputer, you have to use ``job_albedo``, in other case use the job script for your specific platform, or try to modify one of the existing ones.

In the job file, the changes are done based on the HPC you are using. For ``levante``, you should adapt for example:

- ``#SBATCH --job-name``: name of your experiment; e.g. myexperiment_001;

- ``#SBATCH --ntasks-per-node``: number of cores per node. This number has to be divisible by the number of tasks. If you choose the ``ntasks``/4, for example, you will run your experiment with 4 nodes;

- ``#SBATCH --ntasks``: number of cores. This number has to be the same of your desired mesh partitioning. It is the ``xx`` number in your ``dist_xx`` mesh folder;

- ``#SBATCH --time``: be generous with your run time, in case you are running a longer simulation and the job is not being resubmmited after each time step;

- ``#SBATCH -A <account>``: define your project account.


On ``levante`` the submission of your job is done by executing the following command:

::

    sbatch job_levante

The job is then submitted. In order to check the status of your job on ollie you can execute:

::

    squeue -u yourusername

The output of the model run should appear in the ``results`` directory that you have specified in the ``namelist.config``. After the run is finished the ``fesom.clock`` file (or if you change your runid, ``runid.clock``)  will be updated with information about the time of your run's end, that allows running the next time portion of the model experiment by just resubmitting the job with ``sbatch job_ollie``.

Some files will also be stored on the work folder. Those are

- A file containing information about errors during job preparation and submission, usually containing ``err.out`` in its name;

- A file containing information about the job itself, such as duration, folders, etc, usually contining ``out.out`` in its name;

- A file containing information about the simulation, usually called ``fesom2-0.out``;
  
- A binary file ``fesom.x`` specific to that simulation;

- A copy of the namelists used to define the configurations of your run. 

In case your simulation crashes, usually the job error file or ``fesom2-0.out`` contain valuable information to either fix the issue causing the crash or to give the developers an idea of what can be done to help you.


Other things you need to know earlier on
========================================

The clock file
--------------

The clock file is located in your output directory (specified in ``ResultPath`` option of ``namelist.config``) and controls the time. At the start of a new experiment that we want to initialize from climatology (a so-called cold start), the ``fesom.clock`` file would usually look like this:

::

    0 1 1958
    0 1 1958

In this example, ``1958`` is the first available year of the atmospheric ``JRA55`` forcing. The two identical lines tell the model that this is the start of the experiment and that there is no restart file to be read. Also make sure that the ``yearnew`` option of the ``namelist.config`` is set to the year you would like the cold start to begin (1958 in this case).

Let's assume that we run the model with a timestep of 30 minutes (= 1800 seconds) for a full year (1948). After the run is successfully finished, the clock file will then automatically be updated and look like this:

::

    84600.0 365 1958
    0.0     1   1959

where the first row is the second of the day of the last time step of the model, and the second row gives the time when the simulation is to be continued. The first row indicates that the model ran for 365 days (in 1958) and 84600 seconds, which is ``1 day - 1`` FESOM timestep in seconds. In the next run, FESOM2 will look for restart files for the year 1958 and continue the simulation at the 1st of January in 1959.


Tricking FESOM2 into accepting existing restart files
-----------------------------------------------------
The simple time management of FESOM2 allows to easily trick FESOM2 to accept existing restart files. Let's assume that you have performed a full ``JRA55`` cycle until the year 2019 and you want to perform a second cycle, restarting from the last year of the first cycle. This can be done by (copying and) renaming the last year into:

::

    mv fesom.2019.ice.nc fesom.1957.ice.nc
    mv fesom.2019.oce.nc fesom.1957.oce.nc

by changing the clock file into:

::

    84600.0 365 1957
    0.0     1   1958
    
In case the second cycle starts again at the very first year (e.g. 1958 in ``JRA55``) of the forcing, namelist.config needs to be modified, otherwise the model will always perform a cold start in 1958 instead of restarting from the 1957 restart files:

::

    &clockinit
    timenew=0.0
    daynew=1
    yearnew=1957



.. _partitioning:

Build partitioner executable
----------------------------

First meshes you will use probably will come with several predefined partitionings (``dist_XXXX`` folders). However at some point you might need to create partitioning yourself. To do so you have to first compile the partitioner. First you change to the ``mesh_part`` directory:

::

    cd mesh_part

if you work on the one of the supported systems, you shoule be able to execute:

::

    bash -l ./configure.sh

or, in case of the Ubuntu, or other customly defined system:

::

    bash -l ./configure.sh ubuntu

The ``cmake`` should build the partitioner for you. If your system is not supported yet, have a look on how to add custom system in `Adding new platform for compilation`_. The executable ``fesom_ini.x`` should now be available in ``bin`` directory. Now you can proceed with `Running mesh partitioner`_.


Running mesh partitioner
------------------------

You have to do this step only if your mesh does not have partitioning for the desired number of cores yet. You can understand if the partitioning exists by the presence of the ``dist_XXXX`` folder(s) in your mesh folder, where XXX is the number of CPUs. If the folder contains files with partitioning, you can just skip this step.

Partitioning is going to split your mesh into pieces that correspond to the number of cores you going to request. Now FESOM2 scales until 300 vertices per core, further increase in the amount of cores will probably have relatively small effect.

In order to tell the partitioner how many cores you need the partitioning for, one has to edit ``&machine`` section in the ``namelist.config`` file (see also :ref:`chap_general_configuration`). There are two options: ``n_levels`` and ``n_part``. FESOM mesh can be partitioned with use of several hierarchy levels and ``n_levels`` define the number of levels while ``n_part`` the number of partitions on each hierarchy level. The simplest case is to use one level and ``n_part`` just equal to the number of cores and we recoment to use it at the beggining:

::

    n_levels=1
    n_part= 288

This will prepear your mesh to run on 288 computational cores.

In order to run the partitioner change to the ``work`` directory. You should find several batch scripts that are used to submit partitioner jobs to HPC machines (have ``_ini_`` in their names). The scripts also links ``fesom_ini.x`` executable to the ``work`` directory and copy namelists with configurations from ``config`` folder (for partitioner we actually need only ``namelist.config``, but scripts copy everything).

.. note::
   For the partitioner to run, the ``fesom_ini.x`` executable, configuration namelists (in particular ``namelist.config``) and job script have to be located in the same directory (usually ``work``).

If you are working on AWI's ``ollie`` supercomputer, you have to use ``job_ini_ollie``, in other case use the job script for your specific HPC platform, or try to modify one of the existing ones. For relativelly small meshes (up to 1M nodes) and small partitions it is usually fine just to run the partitioner on a login node (it is serial anyway), like this:

::

    ./fesom_ini.x

.. note::
   Make sure that you have the same enviroment that was used during compilation of ``fesom_ini.x``. Usually the easiest way to do this is to first (example for ``ollie`` platform)::

       source ../env/ollie/shell


   This file (``shell``) is used to setup the environment during the compilation of both ``fesom_ini.x`` and ``fesom.x``.

If you trying to partition large mesh, then on ``ollie`` for example the submission of your partitioning job is done by executing the following command:

::

    sbatch job_ini_ollie


Model spinup / Cold start at higher resolutions
-----------------------------------------------

Cold start of the model at high mesh resolutions with standard values for timestep and viscosity will lead to instabilities that cause the model to crash. If no restart files are available and a spinup has to be performed, the following changes should be made for the first month long simulation and then adjusted gradually over the next 6-8 months:

- First thing to try, that usually helps, is to set in the ``namelist.oce``::

    w_split=.true.

- Try to reduce the timestep in ``namelist.config``, for example to:

  ::

      step_per_day=720

  or even lower (e.g. value 1440 will lead to 1 minute timestep).

.. note::
   Make sure that for the high resolution runs (with mesh resolution over considerable portions of the domain finer than 25-10 km) you don't use the combination of default "Easy Backscatter" vescosity (``visc_option=5``) and ``easy_bs_return= 1.5``. This is true not only for the spinup, but for the whole duration of the run. The "Easy Backscatter" option works very good on low resolution meshes, but for high resolution meshes (eddy resolving) it makes more harm than good. If you would like to use ``visc_option=5`` for high resolution runs, put ``easy_bs_return= 1.0``.


- In ``namelist.oce`` make sure that ``visc_option`` is set to 7 or 5 (see also the note above about option 5) and increase ``gamma1`` to something like:

  ::

      gamma1=0.8


or even higher. After running for about a month try to reduce it. If you change the values of run lengh and restart output frequency (which you probably want to do during the spinup, to run for short periods), don't forget to change them back in the ``namelist.config``:

::

    run_length= 1
    run_length_unit='m'
    ...
    restart_length=1
    restart_length_unit='m'

Increase the timestep gradually. Very highly resolved meshes may require an inital timestep of one-two minutes or even less.

Adding new platform for compilation
-----------------------------------

In order to add a new platform for compilation, you simply have to specify the computational environment. In a simplest case this requires:

- To edit the ``env.sh`` file.
- To add a folder with the name of the platform to the ``env`` folder and put the ``shell`` file with enrionment setup.

In the ``env.sh`` file you have to add one more ``elif`` statement in to the ``if`` control stucture, where the platform (let's call it ``mynewhost``) is selected::

    elif [[  $LOGINHOST = mynewhost ]]; then
        STRATEGY="mynewhost"

As you can see in the ``env.sh`` file some host systems are authomatically identified by using regular expressions, but the simpliest way is just to explicitly provide the name of the host system.

The next step is to create additional folder in the ``env`` folder::

    mkdir ./env/mynewhost

and add a file name with the name ``shell`` to it. This file will be sourced before the compilation, so you can setup the environment (bash syntax) in it. Please have a look at the ``shell`` file in other folders for examples. Now you should be able to do::

    bash -l ./configure.sh mynewhost

to do the compilation.

If you are lucky this will be everything you need. However in more complicated cases one  had to adjust CMake files (``CMakeLists.txt`` located in folders), so the knowlege of CMake is required.

Change compiler options
-----------------------

Compiler options for FESOM2 code can be changed in the ``./src/CMakeLists.txt`` file. Currently the defenition of compiler options for Intel compiler looks like::

    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
        target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero)

At present only Intel and GNU compilers are supported, but the user can realtivelly easy add options by following the same pattern.


Troubleshooting
===============

Error ``can not determine environment for host:``
-------------------------------------------------

If you on Ubuntu system, add ``ubuntu`` as input parameter for ``configure.sh``:

::

    ./configure.sh ubuntu

Otherwise you have to add another system - have a look at `Adding new platform for compilation`_ section.

Model blows up
--------------

There could by many reasons for this, but the first thing to try is to reduce time step or/and increase model viscosity for short period of time. Have a look at `Model spinup / Cold start at higher resolutions`_ for instructions.


Docker based installation
=========================

The best way to run the model locally is to use Docker container. You obviously have to have Docker installed for your system. The Docker image we are going to use have all necessary libraries installed plus have the ``mkrun`` python script (`Docker file`_), that helps to create FESOM2 configurations. As a result of the steps below, you will run ``pi`` mesh for one day using data files that comes with the model.

.. _Docker file: https://github.com/FESOM/FESOM2_Docker/tree/master/fesom2_test

- Get the image::
    
    docker pull ghcr.io/fesom/fesom2_docker:fesom2_test_refactoring-master
    # if you use Mac Silicon (M1 M2 etc) use:
    docker pull --platform linux/amd64 ghcr.io/fesom/fesom2_docker:fesom2_test_refactoring-master

- Go to the folder with your version of fesom2 folder (NOT inside fesom2 folder, one up, the one you run ``git clone https://github.com/FESOM/fesom2.git`` in).
- Run::

    docker run -it -v "$(pwd)"/fesom2:/fesom/fesom2 koldunovn/fesom2_test:refactoring2 /bin/bash
    # if you use Mac Silicon (M1 M2 etc) use:
    docker run --platform linux/amd64 -it -v "$(pwd)"/fesom2:/fesom/fesom2 ghcr.io/fesom/fesom2_docker:fesom2_test_refactoring-master /bin/bash

- This should get you inside the container. You now can edit the files in your fesom2 folder (on host system), but run compule and run the model inside the container.
- When inside the container, to compile do:
  
  ::

    cd fesom2
    git checkout refactoring
    bash -l configure.sh ubuntu

- To prepare the run (this will do the test with pi mesh)::

    mkrun pi test_pi -m docker

- To run the model:

  ::
    
    cd work_pi/
    ./job_docker_new

As a next step you can modify the setup in ``work_pi`` to try different parameters. You can also follow the steps described in `Detailed steps of compiling and runing the code`_. To make your life a bit easier place ``FESOM2_one_year_input`` in the ``fesom2`` folder, so that the data are available inside the container. You also can generate setup that would use ``JRA55`` forcing, and adjust it - this will save you some time on editing ``namelist.forcing``, since original setup in ``work_pi`` folder use old ``CORE2`` forcing. 

  ::

    mkrun pi_jra55 test_pi -m docker -f JRA55

Necessary Ubuntu packages
=========================

Here is the list of packages you need to install on ``Ubuntu`` to compile and run FESOM2. Should work (with adjustments for package managers and names) for other linux distributions.

  ::

    apt-get -y install make gfortran gcc g++ libblas-dev libopenmpi-dev
    apt-get -y install cmake vim git libnetcdf-dev libnetcdff-dev libpmi2-pmix


