============================
GitHub Actions CI Test Suite
============================

We run a series of automatic tests to ensure the model compiles and runs with various features enabled. 

If you want to add a new test, you can add a new test configuration under ``.github/workflows``, using the main configuration as an example. Each test will:

* Compile FESOM in an Ubuntu container
* Prepare a minimal test run folder with appropriate namelists
* Run one day, with one restart
* Run a checksum comparison of produced output vs. precomputed "known truths"


You may need to add a specific mesh or specific input files. Some important locations are:

* ``./test/meshes`` for grids. Please try to use a PI mesh (or something similarly small)
* ``./test/input/`` for forcing files. If you need these copied into the work folder, you will need to add those steps to your GitHub Actions workflow defintion! An example is provided in the iceberg test, where we include some non-standard input files.


For questions regarding the test suite, please get in touch with paul.gierz@awi.de (GitHub @pgierz).