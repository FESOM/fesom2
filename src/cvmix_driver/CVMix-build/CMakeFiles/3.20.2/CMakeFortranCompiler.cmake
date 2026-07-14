set(CMAKE_Fortran_COMPILER "/albedo/soft/sw/spack-sw/intel-oneapi-mpi/2021.6.0-ypwmdz5/mpi/2021.6.0/bin/mpiifort")
set(CMAKE_Fortran_COMPILER_ARG1 " -qmkl")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "2021.6.0.20220226")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/albedo/soft/sw/spack-sw/intel-oneapi-mpi/2021.6.0-ypwmdz5/mpi/2021.6.0/include;/albedo/soft/sw/spack-sw/intel-oneapi-mkl/2022.1.0-7235vfh/mkl/2022.1.0/include/intel64/lp64;/albedo/soft/sw/spack-sw/intel-oneapi-mkl/2022.1.0-7235vfh/mkl/2022.1.0/include;/albedo/soft/sw/spack-sw/intel-oneapi-compilers/2022.1.0-kwhaebe/compiler/2022.1.0/linux/compiler/include/intel64;/albedo/soft/sw/spack-sw/intel-oneapi-compilers/2022.1.0-kwhaebe/compiler/2022.1.0/linux/compiler/include/icc;/albedo/soft/sw/spack-sw/intel-oneapi-compilers/2022.1.0-kwhaebe/compiler/2022.1.0/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/8/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core;iomp5;mpifort;mpi;dl;rt;pthread;mkl_intel_lp64;mkl_intel_thread;mkl_core;iomp5;ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/albedo/soft/sw/spack-sw/intel-oneapi-mpi/2021.6.0-ypwmdz5/mpi/2021.6.0/lib/release;/albedo/soft/sw/spack-sw/intel-oneapi-mpi/2021.6.0-ypwmdz5/mpi/2021.6.0/lib;/albedo/soft/sw/spack-sw/intel-oneapi-mkl/2022.1.0-7235vfh/mkl/2022.1.0/lib/intel64;/albedo/soft/sw/spack-sw/intel-oneapi-compilers/2022.1.0-kwhaebe/compiler/2022.1.0/linux/compiler/lib/intel64_lin;/albedo/soft/sw/spack-sw/intel-oneapi-compilers/2022.1.0-kwhaebe/compiler/2022.1.0/linux/lib;/albedo/soft/sw/spack-sw/intel-oneapi-mpi/2021.6.0-ypwmdz5/mpi/2021.6.0/libfabric/lib;/usr/lib/gcc/x86_64-redhat-linux/8;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
