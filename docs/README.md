FESOM2 tutorial
==============

Very fast version
----------------
```bash
svn checkout --username yourusername https://swrepo1.awi.de/svn/awi-fvom/branches/fvom_ALE/ fesom2 

cd fesom2/lib/metis-5.1.0
make clean
make

cd ../parms/
make cleanall
make

cd ../../src/
make clean
make

cd ../
mkdir results
cd results
```

creater file fesom.clock with the following content:

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

First thing is to checkout FESOM 2 code from repository. Now it is located at AWI FORGE https://swrepo1.awi.de/, that is hosted by AWI. You have to register there and ask Dmitry Sidorenko to add you to the members of the project, that have access to the code.

The project with FESOM2 located here: 
https://swrepo1.awi.de/scm/?group_id=61

At this page you will see an svn command, that will allow you to checkout the code using svn, something like this:

```bash
svn checkout --username koldunovn https://swrepo1.awi.de/svn/awi-fvom/trunk
```

At present we are working on the fvom_ALE brunch, so the command that you actually have to use should be:

```bash
svn checkout --username koldunovn https://swrepo1.awi.de/svn/awi-fvom/branches/fvom_ALE/ fesom2 
```
Before you can compile FESOM2 itself, you have to compile two libraries that are nessesary for work, and that come together with FESOM2. The libraries are:

 * metis (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) - that maps fesom computational domain to CPUs 

 * parms (http://www-users.cs.umn.edu/~saad/software/pARMS/) - library of parallel solvers.

It is better to use versions of libraries that come with fesom.

Metis compilation
----------------
```bash
cd lib/metis-5.1.0
make clean
make install
```
If you are lucky, everything is going to work. If not, you have to get in to the Makefile and see how to fix it. In all Makefiles that comes with fesom there is an include statement:
```Makefile
 include Makefile.in
```
The `Makefile.in` file contains platform specific settings, like name of your compiler, path to the libraries and so on. If you are trying to compile `metis`, or any other part of the fesom distribution not at the default platform (`ollie`), search the the directory for files like `Makefile.in_ubuntu` or `Makefile.in_dkrz`. If you find a file for your platform, just copy it to the `Makefile.in`:
```bash
cp Makefile.in_ubuntu Makefile
```
cross your fingures and hope it will work.

Parms compilation
----------------
THe procedure of parms compilation is very similar to the one for metis, except you have to use `make cleanall` instead of `make clean`:

```bash
cd ../parms/
make cleanall
make
```
You also have to use `Makefile.in` that is created for your platform.

FESOM2 compilation
-----------------
The compilation of FESOM2 is also quite standard if all defaults are working for you:

```bash
cd ../../src/
make clean
make
```
Once again pay attention to `Makefile.in_platform`.

Preparing the run
----------------

You have to do several basic things in order to prepare the run. First create directory where results will be stored. Usually it is created in the model root directory:

```bash
mkdir results
```

you might make a link to some other directory located on the part of the system where you have a lot of storage. In the `results` directory you have to create `runid.clock` file, where `runid` is usually `fesom`, so the file is most often named `fesom.clock`. Inside the file you have to put two identical lines:
```
0 1 1948
0 1 1948
```
the meaning of this lines will be explained later. 

The next step is make some changes in the model configuration. All runtime options can be se in the namelists that are located in the `config` directory:
```bash
cd ../config/
```
There are four configuration files, but we are only interested in the `namelist.config`. The options that you might whant to change for your first FESOM2 run are:

* `run_length` length of the model run in `run_length_unit` (see below)
* `run_length_unit` units of the `run_length`. Can be `y`(year), `m`(month), `d`(days), `s`(seconds)
* `MeshPath` - path to the mesh you would like to use
* `ClimateDataPath` - path to the file with model temperature and salinity initial conditions
* `ForcingDataPath` - path to the forcing data

Running the model
----------------
Change to the `bin` directory. You should find the your compiled executable `fvom.x` and several batch scripts, that one use to submit model execution job to HPC machines. If you are working on AWI's `ollie`, you have to use `job_ollie`, in other case use the job script for your specific HPC platform, or try to modify one of existing ones.

On `ollie` the submission of your job is done by executing the following command:
```bash
sbatch job_ollie
```
The job is then submitted to the supercomputer. In order to check the status of your job on `ollie` you can execute:
```bash
squeue -u yourusername
```

































