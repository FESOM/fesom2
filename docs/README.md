FESOM2 tutorial
==============

Very fast version
----------------
```bash
git clone https://gitlab.dkrz.de/FESOM/fesom2.git fesom2 
cd fesom2
bash -l configure.sh
```

create file fesom.clock in the output directory with the following content:

```bash
  0 1  1948
  0 1  1948
```

```bash
cd ../bin
sbatch job_ollie
```

Longer version
-------------

First thing is to checkout FESOM 2 code from the repository. Now it is located at DKRZ gitlab (https://gitlab.dkrz.de/). You have to have an active DKRZ account (more information on how to register [here](https://www.dkrz.de/up/my-dkrz/getting-started/account/DKRZ-user-account)). Before your account can be added to the members of the project you have to at least once log in to https://gitlab.dkrz.de/ with your DKRZ account so that it appears in the GitLab data base. To be added as a member of the FESOM project on GitLab you have to send an email to Dmitry Sidorenko (dmitry.sidorenko@awi.de) with the request. Please make sure that you are able to login to https://gitlab.dkrz.de/ before sending the request.

The project with FESOM2 is located here: 
https://gitlab.dkrz.de/FESOM/fesom2/

Clone the repository with git command:

```bash
git clone https://gitlab.dkrz.de/FESOM/fesom2.git
```
The repository contains model code and two additional libraries (`Metis` and `Parms`) necessary for running FESOM2. To build FESOM2 executable one have to compile `Parms` library and the code of the model (`src` folder). In order to build executable that is used for model domain partitioning (distribution of the model mesh between CPUs) one have to compile `Metis` library and also some code located in the `src` directory. The first step of building model executable is usually done automatically with the use of CMake, but it will work only on supported platforms (ollie, DKRZ, HLRN, and HAZELHEN). To build executable for partitioning and model executable on other, unsupported, platforms you can still use good old Makefiles for now, but it might change in the future in favor of complete CMake based build process. Both ways will be explained below. 

#### Build model executable with Cmake (default and recommended)

Go to the `fesom2` folder and execute:

```bash
bash -l configure.sh 
```
In the best case scenario, your platform will be recognized and the `Parms` library and model executable will be built and copied to the `bin` directory. If it does not work - try to use `make` (see below).

#### Build partitioner executable with make (the only way for now)

First, build `Metis` library. In the ideal case, the only things you have to do is:

```bash
cd lib/metis-5.1.0
make clean
make install
```

However, you might need modification of the build options in your `Makefile`. The `Makefile` consist of the two parts - generic one contained in the `Makefile` itself and the `Makefile.in` that have all platform/compiler specific options and is "imported" during the `make` command execution. There are several `Makefile.in_*` files present, and you have to copy the one that is suitable for your platform into the `Makefile.in`. 

After successful compilation of `Metis` you can build partitioner executable:
```bash
cd ../src/
make run_ini
```
The `Makefiles` have the same structure, and you might need to replace default `Makefile.ini` with the one that is suitable for your platform.

#### Build model executable with make (not recommended, but sometimes the only way)

If Cmake fails, your platform is not supported, or you would like to do something crazy, you still can use make.
First, we have to build `Parms` library:

```bash
cd ../parms/
make cleanall
make
```
note the `cleanall` instead of usual `clean`.

To build model executable:
```bash
cd ../../src/
make clean
make
```

Please see note related to the Makefiles in the previous section.

Preparing the run
----------------

You have to do several basic things in order to prepare the run. First, create a directory where results will be stored. Usually, it is created in the model root directory:

```bash
mkdir results
```

you might make a link to some other directory located on the part of the system where you have a lot of storage. In the `results` directory, you have to create `runid.clock` file, where `runid` is usually `fesom`, so the file is most often named `fesom.clock`. Inside the file you have to put two identical lines:
```
0 1 1948
0 1 1948
```
the meaning of this lines will be explained later. 

The next step is to make some changes in the model configuration. All runtime options can be set in the namelists that are located in the `config` directory:
```bash
cd ../config/
```

There are four configuration files, but we are only interested in the `namelist.config` for now. The options that you might want to change for your first FESOM2 run are:

* `run_length` length of the model run in `run_length_unit` (see below)
* `run_length_unit` units of the `run_length`. Can be `y`(year), `m`(month), `d`(days), `s`(model steps)
* `MeshPath` - path to the mesh you would like to use
* `ClimateDataPath` - path to the file with model temperature and salinity initial conditions
* `ForcingDataPath` - path to the forcing data

Running mesh partitioner
-----------------------
You have to do this step only if your mesh does not have partitioning for the desired amount of CPU yet. You can understand if the partitioning exists by the presence of the `dist_XXXX` folder in the folder with your mesh, where `XXX` is the number of CPUs. If the folder contains files with partitioning, you can just skip to the next step, otherwise see instructions in this section. 

Partitioning is going to split your mesh into pieces that correspond to the number of CPUs you going to request. Now FESOM2 scales until 300 nodes per CPU, further increase in the amount of CPU will probably have relatively small effect.

In order to tell the partitioner how many CPUs you need the partitioning for, one has to edit `&machine` section in the `namelist.config` file. There are two options: `n_levels` and `n_part`. FESOM mesh can be partitioned with use of several hierarchy levels and `n_levels` define the number of levels while `n_part` the number of partitions on each hierarchy level. The simplest case is to use one level and `n_part` just equal to the number of CPUs:
```bash
n_levels=1
n_part= 288           
```
If you would like to use more level, note that the product of levels should be equal to the number of CPUs, so for three levels and 288 CPUs:
```bash
n_levels=3
n_part= 2, 4, 36
```

Then change to the `work` directory. You should find several batch scripts that are used to submit model jobs to HPC machines. The scripts also link `fvom_ini.x` executable to the `work` directory and copy namelists with configurations from `config` folder. **NOTE, model executable, configuration namelists and job script have to be located in the same directory (usually `work`)**. 
 
If you are working on AWI's `ollie` supercomputer, you have to use `job_ini_ollie`, in other case use the job script for your specific HPC platform, or try to modify one of the existing ones.

On `ollie` the submission of your partitioning job is done by executing the following command:
```bash
sbatch job_ini_ollie
```
The job is then submitted to the supercomputer. In order to check the status of your job on `ollie` you can execute:
```bash
squeue -u yourusername
```


Running the model
----------------
Change to the `work` directory. You should find several batch scripts that are used to submit model jobs to HPC machines. The scripts also link `fesom.x` executable to the `work` directory and copy namelists with configurations from `config` folder. **NOTE, model executable, configuration namelists and job script have to be located in the same directory (usually `work`)**.  If you are working on AWI's `ollie` supercomputer, you have to use `job_ollie`, in other case use the job script for your specific HPC platform, or try to modify one of the existing ones.

On `ollie` the submission of your job is done by executing the following command:
```bash
sbatch job_ollie
```
The job is then submitted to the supercomputer. In order to check the status of your job on `ollie` you can execute:
```bash
squeue -u yourusername
```

Results of the model run should appear in the `results` directory that you have specified in the `namelist.config`. After the run is finished the `fesom.clock`file (or if you change your `runid`, `runid.clock`)  will be updated with information about the time of your run's end, that allows running the next time portion of the model experiment by just resubmitting the job with `sbatch job_ollie`. 


The clock file
----------------

The clock file is usually located in your output directory (e.g. `results`) and controls the time. 
At the start of a new experiment that we want to be initialized from climatology (a so-called cold start),
the `fesom.clock` file would usually look like this
```
0 1 1948
0 1 1948
```
In this example, 1948 is the first available year of the atmospheric `CORE2` forcing. The two identical
lines tell the model that this is the start of the experiment and that there is no restart file to be read.

Let's assume that we run the model with a timestep of 30 minutes (= 1800 seconds) for a full year (1948).
After the run is successfully finished, the clock file will then automatically be updated and look like this:
```
84600.0 365 1948
0.0     1   1949
```
where the first row is the second last time step of the model, and the second row gives the time where the simulation
is to be continued. The first row indicates that the model ran for 365 days (in 1948) and 84600 seconds, which is `1 day - 1 FESOM timestep`
in seconds. In the next job, FESOM2 will look for restart files for the year 1948 and continue the simulation at the 1st of January in 1949.

Since 1948 is a leap year (366 days), this is an exceptional case and the `fesom.clock` file after two full years
(1948--1949) would look like this:
```
84600.0 364 1949
0.0     1   1950
```
Note that dependent on the forcing data set (using a different calendar), a year could only have 360 days.

Tricking FESOM2 into accepting existing restart files
----------------

The simple time management of FESOM allows to easily trick FESOM to accept existing restart files. Let's assume that you have performed a full `CORE2` cycle until the year 2009 and you want
to perform a second cycle, restarting from the last year of the first cycle. This can be done by (copying and)
renaming the last year into:
```bash
mv fesom.2009.ice.nc fesom.1947.ice.nc
mv fesom.2009.oce.nc fesom.1947.oce.nc
```
By changing the clock file into
```
84600.0 364 1947
0.0     1   1948
```
the model will then restart from the last year of the first cycle, but using CORE2 forcing from 1948 onwards.


Todo
----

#### recomendations on model spin up


