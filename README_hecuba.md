# Description of modifications to master source for hecuba IO template

Note additional directory in `src`:`src/hecuba_io/` containing a module for hecuba io `io_mod_hecuba.F90`. This module acts as interface between hecuba/cassandra and fesom.
As an example module contains a contains a subroutine `output_hecuba(istep)` that is modelled similarly to that of `src/io_meandata.F90`. To demo its working it currently writes 2 variables `sst` and `temp` into binary files for every timestep. Idea is if we can replace regular writes in the module with that of hecuba interface API  and write to cassandra. Once this can be done we can extend to more variables and have a namelist to control data writing and what needs to be written.

The module is used and subroutine is called in `src/fvom_main.F90` see modifications there which are wrapped around preprocessor option `USE_HECUBAIO` set in master CMakeLists.txt as `set(USE_HECUBAIO ON CACHE BOOL "compile fesom with hecuba interface to cassandra DB for IO.")`. this allows to enable and disable the interface safely without disturbing the model setup.

similarly same flag `USE_HECUBAIO` is used in `src/CMakeLists.txt` to include sources specific to interface for compiling together with the fesom2 executable. You may optionally want to add hecuba's C API library to this cmakelist to also compile it along fesom2 src either as gitsubmodule or just subdirectory. 

fesom2's original way to setup remains same, e.g., `./configure.sh ubuntu` 
