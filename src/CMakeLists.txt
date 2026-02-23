cmake_minimum_required(VERSION 3.16)

project(fesom LANGUAGES C CXX Fortran VERSION 2.0.0)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake ${CMAKE_MODULE_PATH})

include(fesom_setup)

find_package(MPI REQUIRED)

#TODO: these machine specific changes must be done using cmake toolchain then here otherwise will make cmake another raps. 
if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
   set(FESOM_PLATFORM_STRATEGY $ENV{FESOM_PLATFORM_STRATEGY} CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
else()
   set(FESOM_PLATFORM_STRATEGY "notset" CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
endif()

# Cray compiler specific optimizations
if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
   # Use alternative cray-mpich library, about 5% faster with cray-mpich/7.7.3
   # Not available for modules cray-mpich > 7.7.3
   add_compile_options(-craympich-mt)
   # Make sure to also set these variables in the runtime environment:
   #    MPICH_MAX_THREAD_SAFETY=multiple # allows highest MPI thread level (i.e. MPI_THREAD_MULTIPLE)
   #    MPICH_CRAY_OPT_THREAD_SYNC=0 # the Cray MPICH library falls back to using the pthread_mutex-based thread-synchronization implementation
   #    MPICH_OPT_THREAD_SYNC=0 # seems to be a duplicate variable which also appears in some documentation instead of MPICH_CRAY_OPT_THREAD_SYNC (but this one brings a huge speed gain on aleph)
endif()

# Backward compatibility: convert old DISABLE_MULTITHREADING to new ASYNCHRONOUS_IO_THREADS
if(DEFINED DISABLE_MULTITHREADING)
   if(DISABLE_MULTITHREADING)
      set(ASYNCHRONOUS_IO_THREADS OFF CACHE BOOL "enable asynchronous I/O threads" FORCE)
   else()
      set(ASYNCHRONOUS_IO_THREADS ON CACHE BOOL "enable asynchronous I/O threads" FORCE)
   endif()
   message(STATUS "Converted legacy DISABLE_MULTITHREADING=${DISABLE_MULTITHREADING} to ASYNCHRONOUS_IO_THREADS=${ASYNCHRONOUS_IO_THREADS}")
endif()

#TODO: these machine specific block can easyly go in cmake toolchain.
if    (${FESOM_PLATFORM_STRATEGY} STREQUAL levante.dkrz.de )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL leo-dcgp )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL mn5-gpp )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL albedo )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL atosecmwf )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL juwels )
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)   
else()
   option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)
endif()
if (ASYNCHRONOUS_IO_THREADS)
   message(STATUS "ASYNCHRONOUS I/O THREADS @ ${FESOM_PLATFORM_STRATEGY}: ON")
else()
   message(STATUS "ASYNCHRONOUS I/O THREADS @ ${FESOM_PLATFORM_STRATEGY}: OFF")
endif()


option(ENABLE_OPENACC "compile with OpenACC support" OFF)
message(STATUS "ENABLE_OPENACC: ${ENABLE_OPENACC}")
option(DISABLE_OPENACC_ATOMICS "disable kernels using atomic statement for reproducible results" ON)
set(GPU_COMPUTE_CAPABILITY "cc80" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")
# set(GPU_FLAGS "${GPU_COMPUTE_CAPABILITY}" CACHE STRING "GPU arch for nvfortran compiler")
set(GPU_FLAGS "cuda12.2,${GPU_COMPUTE_CAPABILITY}" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")

option(ENABLE_OPENMP "build FESOM with OpenMP" OFF)
message(STATUS "ENABLE_OPENMP: ${ENABLE_OPENMP}")
if(${ENABLE_OPENMP})
   find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

option(USE_ICEPACK "Use ICEPACK" OFF)
message(STATUS "USE_ICEPACK: ${USE_ICEPACK}")


# option to trigger building a library version of FESOM
# we do not always build the library along with the executable to avoid having two targets here in the CMakeLists.txt
# two targets would allow e.g. setting different compiler options or preprocessor definition, which would be error prone.
set( ENABLE_IFS_INTERFACE_DEFAULT OFF )
if(DEFINED BUILD_FESOM_AS_LIBRARY AND NOT DEFINED ENABLE_IFS_INTERFACE)
  message(DEPRECATION "Use ENABLE_IFS_INTERFACE instead of BUILD_FESOM_AS_LIBRARY")
  set(ENABLE_IFS_INTERFACE_DEFAULT ${BUILD_FESOM_AS_LIBRARY})
endif()
option(ENABLE_IFS_INTERFACE   "Enable IFS interface" ${ENABLE_IFS_INTERFACE_DEFAULT})
if(DEFINED FESOM_ENABLE_IFS_INTERFACE)
  set(ENABLE_IFS_INTERFACE ${FESOM_ENABLE_IFS_INTERFACE}) # To distinguish option in a nested cmake project (bundle)
endif()
message(STATUS "ENABLE_IFS_INTERFACE: ${ENABLE_IFS_INTERFACE}")


set(ENABLE_MULTIO_DEFAULT OFF)
if(ENABLE_IFS_INTERFACE)
   if(USE_ICEPACK)
      message(FATAL_ERROR "Could not enable IFS interface as it is incompatible with USE_ICEPACK")
   endif()
   if( NOT DEFINED ENABLE_MULTIO )
     find_package(multio QUIET)
     if( multio_FOUND )
      set(ENABLE_MULTIO_DEFAULT ON)
     endif()
   endif()
endif()

option(ENABLE_MULTIO "Enable MultIO" ${ENABLE_MULTIO_DEFAULT})
if(DEFINED FESOM_ENABLE_MULTIO)
  set(ENABLE_MULTIO ${FESOM_ENABLE_MULTIO}) # To distinguish option in a nested cmake project (bundle)
endif()
message(STATUS "ENABLE_MULTIO: ${ENABLE_MULTIO}")
if(ENABLE_MULTIO)
  find_package(multio)
  if(multio_FOUND)
    message(STATUS "Found multio: ${multio_DIR}")
  else()
    message(FATAL_ERROR "Could not enable multio as multio could not be found")
  endif()
  if(NOT TARGET multio-fapi)
    message(FATAL_ERROR "Could not enable multio as multio was not compiled with 'multio-fapi' target")
  endif()
endif()


set(FESOM_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/../.." CACHE FILEPATH "directory where FESOM will be installed to via 'make install'")

# get our source files
set(src_home ${CMAKE_CURRENT_SOURCE_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt

# Link FESOM2 fortran code
file(GLOB sources_Fortran ${src_home}/*.F90)

# Link custom CVMix files and fesom2cvmix wrapper-interface files
if(CVMIX)
   file(GLOB sources_cvmix_addon ${src_home}/cvmix_driver/*.F90)
   list(APPEND sources_Fortran ${sources_cvmix_addon})
endif()

# Link ICEPACK code
if(USE_ICEPACK)
   file(GLOB sources_ice_pack ${src_home}/icepack_drivers/*.F90
                              ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
   list(APPEND sources_Fortran ${sources_ice_pack})
   
# Link IFS-Interface code   
elseif(ENABLE_IFS_INTERFACE) # ICEPACK + IFS_INTERFACE NOT SUPPORTED (YET)
   if(ENABLE_MULTIO)
      list(APPEND sources_Fortran ${src_home}/ifs_interface/iom.F90)
   endif()

   list(APPEND sources_Fortran ${src_home}/ifs_interface/ifs_interface.F90
                               ${src_home}/ifs_interface/ifs_modules.F90
                               ${src_home}/ifs_interface/ifs_notused.F90
                               ${src_home}/ifs_interface/mpp_io.F90)
endif()

# Link RECOM code
if(RECOM_COUPLED)  # ADD RECOM FILES
   file(GLOB sources_recom ${src_home}/int_recom/*.F90
                           ${src_home}/int_recom/recom/*.F90)
   list(APPEND sources_Fortran ${sources_recom})                           
endif()                             

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
set(FESOM_ORIGINAL_VERSION_FILE ${src_home}/fesom_version_info.F90)
set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
list(REMOVE_ITEM sources_Fortran ${FESOM_ORIGINAL_VERSION_FILE}) # we want to compile the generated file instead
list(APPEND sources_Fortran ${FESOM_GENERATED_VERSION_FILE})
set_source_files_properties(${FESOM_GENERATED_VERSION_FILE} PROPERTIES GENERATED TRUE)

add_custom_target(fesom_version_info-generated.F90 ALL
                  COMMENT "Determining ${PROJECT_NAME} git SHA ..."
                  COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

#if(FESOM_STANDALONE)
#   list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_driver.F90)
#endif()
list(REMOVE_ITEM sources_Fortran ${src_home}/fvom_init.F90 ${src_home}/oce_local.F90 ${src_home}/gen_comm.F90)
list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_main.F90)

find_package( NETCDF REQUIRED COMPONENTS C Fortran )
find_package( MPI REQUIRED COMPONENTS C Fortran )

# depends on the metis library
#add_subdirectory(../lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
#include_directories(../lib/metis-5.1.0/include)

add_subdirectory(async_threads_cpp)


# fesom library
add_library(${PROJECT_NAME} ${sources_Fortran} $<TARGET_OBJECTS:async_threads_cpp>)
add_dependencies(${PROJECT_NAME} fesom_version_info-generated.F90)

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:module/fesom>)
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)

# build CVMix library if CVMIX is ON
if (CVMIX)
   message(STATUS "CVMIX: ${CVMIX}")
   
   # read CVMix install and src directories from the download script, so that there
   # is only one place where they are defined
   execute_process(COMMAND bash ${src_home}/cvmix_driver/download_cvmix.sh --print-src_dir
                   OUTPUT_VARIABLE CVMIX_DIR_src OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND bash ${src_home}/cvmix_driver/download_cvmix.sh --print-prefix
                   OUTPUT_VARIABLE CVMIX_DIR_install OUTPUT_STRIP_TRAILING_WHITESPACE)  
   
   # set default CVMix installation path
   set(CVMIX_ROOT "${src_home}/cvmix_driver/${CVMIX_DIR_install}_${CMAKE_Fortran_COMPILER_ID}") 
   message(STATUS "set default CVMIX_ROOT=${CVMIX_ROOT}")
   
   # Where library is located --> fix in download_cvmix.sh that it as always installed 
   # in /lib and not in /lib64
   set(CVMIX_LIB "${CVMIX_ROOT}/lib/libcvmix.a")
   
   # Create parent directory early so CMake does not fail on missing include paths
   file(MAKE_DIRECTORY ${CVMIX_ROOT})
   file(MAKE_DIRECTORY ${CVMIX_ROOT}/include/cvmix)
   
   # Collect CVMix source files so that CMake knows what changes triggers a 
   # possible rebuild of the CVMix library
   file(GLOB_RECURSE CVMIX_SOURCES
         "${src_home}/cvmix_driver/${CVMIX_DIR_src}/src/*.F90"
         "${src_home}/cvmix_driver/${CVMIX_DIR_src}/src/*/*.F90")
        
   # Add script dependency (rebuild if script changes)
   set(CVMIX_DOWNLOAD "${src_home}/cvmix_driver/download_cvmix.sh")
   
   # Custom command: build CVMix if not yet built
   add_custom_command(
         OUTPUT ${CVMIX_LIB}
         COMMAND bash ${src_home}/cvmix_driver/download_cvmix.sh --prefix ${CVMIX_ROOT}
         WORKING_DIRECTORY ${src_home}/cvmix_driver
         DEPENDS ${CVMIX_SOURCES} ${CVMIX_DOWNLOAD} # trigger rebuild if source or download script change
         COMMENT "Downloading and building CVMix into ${CVMIX_ROOT}"
         VERBATIM)
   
   # check if ${CVMIX_ROOT}/lib/libcvmix.a exist or not, only if it does not 
   # exist build library
   add_custom_target(cvmix_autobuild DEPENDS ${CVMIX_LIB})
 
   # build cvmix external libray
   add_library(cvmix_external STATIC IMPORTED GLOBAL)
   set_target_properties(cvmix_external PROPERTIES
         IMPORTED_LOCATION ${CVMIX_LIB}
         INTERFACE_INCLUDE_DIRECTORIES "${CVMIX_ROOT}/include/cvmix")

   # cvmix_autobuild is a custom target that runs your download_cvmix.sh script 
   # and produces libcvmix.a. cvmix_external is an imported library, meaning CMake 
   # doesn’t build it — it just points to an existing file (libcvmix.a).
   add_dependencies(cvmix_external  cvmix_autobuild)
   
   # Link CVMix library to fesom, build chain becomes: fesom → cvmix_external → cvmix_autobuild
   target_link_libraries(${PROJECT_NAME} PRIVATE cvmix_external)
   
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
if(${ENABLE_OPENMP})
   target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(ENABLE_MULTIO)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __MULTIO)
   target_link_libraries(${PROJECT_NAME} PRIVATE multio-fapi)
endif()

if(ASYNCHRONOUS_IO_THREADS)
   target_compile_definitions(${PROJECT_NAME} PRIVATE ASYNCHRONOUS_IO_THREADS)
endif()

if(${FESOM_COUPLED})
   message(STATUS "FESOM is coupled")
   if(${OASIS_WITH_YAC})
      message(STATUS "Looking for OASIS with YAC lib")
      include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASISYAC.cmake)
   else()
      message(STATUS "Looking for OASIS")
      include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
   endif()
   target_compile_definitions(${PROJECT_NAME} PRIVATE __oasis)
endif()
if(${ASYNC_ICEBERGS})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __async_icebergs)
endif()
if(${OIFS_COUPLED})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __oifs)
endif()

if(USE_ICEPACK)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __icepack)
endif()

if(CVMIX)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __cvmix)
endif()

if(ENABLE_IFS_INTERFACE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __ifsinterface)
endif()

if(VERBOSE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()

if(OPENMP_REPRODUCIBLE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __openmp_reproducible)
endif()

if(${RECOM_COUPLED})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __recom USE_PRECISION=2 __3Zoo2Det __coccos)# __usetp)
endif()

if(${CISO_COUPLED})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __ciso)
endif()

if(${FESOM_PROFILING})
   target_compile_definitions(${PROJECT_NAME} PRIVATE FESOM_PROFILING)
endif()

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL  IntelLLVM )
   
   # Base compiler flags
   target_compile_options(${PROJECT_NAME} PRIVATE -O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero -no-wrap-margin -fpe0 -fpp)
   
   # compiler flags not supported by IntelLLVM
   if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
      target_compile_options(${PROJECT_NAME} PRIVATE -no-prec-sqrt -ip )
   endif()   
   
   if(ENABLE_OPENMP)
      target_compile_options(${PROJECT_NAME} PRIVATE -qopenmp)
   endif()   
   
   # platform specific compiler flags
   # add -fpe0 for RAPS environment
   if    (${FESOM_PLATFORM_STRATEGY} STREQUAL levante.dkrz.de )
      # PROCESSORS:  2x AMD 7763 CPU; 128 cores in total, 256 GB main memory
      # --> compiler flag performance test: core2, 512CPUs, 3mon. simulat.
      # | Intel/MPI | nothing                                                             | 144s |
      # | Intel/MPI | -O2                                                                 | 144s |
      # | Intel/MPI | -O3                                                                 | 144s |
      # |*Intel/MPI | -03 -march=core-avx2 -mtune=core-avx2 -qopt-zmm-usage=low           | 140.95s |
      # |           | -align array64byte -unroll-aggressive -qopt-malloc-options=2        |         |
      #target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2 -qopt-zmm-usage=low -align array64byte -unroll-aggressive -qopt-malloc-options=2)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL leo-dcgp )
      target_compile_options(${PROJECT_NAME} PRIVATE -xCORE-AVX512 -qopt-zmm-usage=high -align array64byte -ipo)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL mn5-gpp )
      target_compile_options(${PROJECT_NAME} PRIVATE -xCORE-AVX512 -qopt-zmm-usage=high -align array64byte -ipo)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL albedo)
      # PROCESSORS: 2x AMD Rome Epyc 7702 (64 cores each), 128 cores, 256 GB RAM
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -fPIC -qopt-malloc-options=2 -qopt-prefetch=5 -unroll-aggressive) # -g -traceback -check) #NEC mpi option
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL atosecmwf )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL juwels )
      # PROCESSOR: 2× Intel Xeon Platinum 8168 CPU, 2× 24 cores, 2.7 GHz, 96 (12× 8) GB DDR4, 2666 MHz
      # --> compiler flag performance test: core2, 192CPUs, 1mon. simulat.
      # | Intel/Para...MPI    | nothing                                                   | 309s |
      # | Intel/Para...MPI    | -O2                                                       | 114s |
      # | Intel/Para...MPI    | -O3                                                       | 115s |
      # | Intel/Para...MPI    | -O3 -xCORE-AVX512                                         | 117s |
      # | Intel/Para...MPI    | -O3 -xCORE-AVX512 -qopt-zmm-usage=high -align array64byte | 119s |
      # | Intel/Para...MPI    | -O2 -xCORE-AVX2                                           | 114s |
      # |*Intel/Para...MPI    | -O3 -xCORE-AVX2                                           | 112s |
      # | Intel/Para...MPI    | -O3 -xCORE-AVX2 -qopt-streaming-stores=always             | 126s |
      # | Intel/Para...MPI    | -O3 -xCORE-AVX2 -qopt-prefetch=5                          | 116s |
      # | Intel/Para...MPI    | -O3 -xCORE-AVX2 -funroll-loops                            | 113s |
      # |*Intel/Para...MPI-mt | -O3 -xCORE-AVX2                                           | 113s |
      target_compile_options(${PROJECT_NAME} PRIVATE -xCORE-AVX2)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -xHost)
   endif()
   
   # # debugging related compiler flags
   # target_compile_options(${PROJECT_NAME} PRIVATE -qopenmp -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)
   # target_compile_options(${PROJECT_NAME} PRIVATE -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)


elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
#    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
   if(${FESOM_PLATFORM_STRATEGY} STREQUAL ubuntu )
      message(STATUS "Allowing type mismatches on Ubuntu for CI Testing" )  # NOTE(PG): Would be nicer to grab the CI=True from the env variable
      target_compile_options(${PROJECT_NAME} PRIVATE -O2 -g -fbacktrace -ffloat-store -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fallow-argument-mismatch -cpp)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -O3 -ffloat-store -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp)
   endif()
   if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
      target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
   endif()
   
   # platform specific compiler flags
   if    (${FESOM_PLATFORM_STRATEGY} STREQUAL juwels )
      # PROCESSOR: 2× Intel Xeon Platinum 8168 CPU, 2× 24 cores, 2.7 GHz, 96 (12× 8) GB DDR4, 2666 MHz
      # --> compiler flag performance test: core2, 192CPUs, 1mon. simulat.
      # | GCC/openMPI | nothing                                                           | 390s |
      # | GCC/openMPI | -O2                                                               | 192s |
      # |*GCC/openMPI | -O3 -march=skylake-avx512 -mtune=skylake-avx512                   | 173s | chatgpt recomendation
      # |             | -mprefer-vector-width=512 -falign-loops=64 -falign-functions=64   |      |
      # |             | -falign-jumps=64  (chatGPT recomendation)                         |      |
      target_compile_options(${PROJECT_NAME} PRIVATE -march=skylake-avx512 -mtune=skylake-avx512 -mprefer-vector-width=512 -falign-loops=64 -falign-functions=64 -falign-jumps=64)      
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL levante.dkrz.de )
      # PROCESSORS:  2x AMD 7763 CPU; 128 cores in total, 256 GB main memory
      # --> compiler flag performance test: core2, 512CPUs, 3mon. simulat.
      # | GCC/openMPI | nothing                                                           | 659s |
      # | GCC/openMPI | -O2                                                               | 364s |
      # | GCC/openMPI | -O3                                                               | 316s |
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3                                   | 278s | chatgpt recomendation
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize                  | 278s | chatgpt recomendation
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -funroll-loops                    | 284s | chatgpt recomendation
      # |*GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize -flto            | 276s | chatgpt recomendation
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize -flto            | 281s | chatgpt recomendation
      # |             | -falign-functions=32 -falign-loops=32 -falign-jumps=32            |      |
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize -flto            | 281s | chatgpt recomendation
      # |             | -falign-functions=64 -falign-loops=64 -falign-jumps=64            |      |
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize -flto            | 280s | chatgpt recomendation
      # |             | -mno-fma4 -mavx2 -mfma                                            |      |
      # | GCC/openMPI | -O3 -march=znver3 -mtune=znver3 -ftree-vectorize -flto            | 280s | chatgpt recomendation
      # |             | -mcpu=znver3                                                      |      |
      target_compile_options(${PROJECT_NAME} PRIVATE -march=znver3 -mtune=znver3 -ftree-vectorize -flto)      
   else()
      #[[if(NOT (${FESOM_PLATFORM_STRATEGY} STREQUAL ubuntu))
         target_compile_options(${PROJECT_NAME} PRIVATE -native)
      endif()  ]]  
   endif()
   
   # # debugging related compiler flags
   #target_compile_options(${PROJECT_NAME} PRIVATE 
                           #-fallow-argument-mismatch 
                           #-fno-fast-math -fno-signed-zeros 
                           #-freciprocal-math 
                           #-mveclibabi=svml 
                           #-flto
                           #-pg 
                           #-fbacktrace 
                           #-fcheck=all,bounds 
                           #-finit-real=snan 
                           #-finit-integer=-9999 
                           #-fsanitize=undefined,address
                        #)
   # if use -fsanitize=undefined,address you also need ... PRIVATE -lubsan)           
   #target_link_libraries(${PROJECT_NAME} PRIVATE -lubsan)
                     

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
	#target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O3)
   target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O2 -hnoacc -M878)  #-hnoacc is a workaround for cray automatically activate -hacc, -M878 is to suppress ftn-878 warning
   if(${ENABLE_OPENMP})
      target_compile_options(${PROJECT_NAME} PRIVATE -homp)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -hnoomp)
   endif()
   if(${ENABLE_OPENACC}) # remove debug later
	   target_compile_options(${PROJECT_NAME} PRIVATE -hacc –hlist=a -g -DENABLE_OPENACC )
      if(${DISABLE_OPENACC_ATOMICS})
          target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_OPENACC_ATOMICS)
      endif()
   endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC )
   target_compile_definitions(${PROJECT_NAME} PRIVATE ENABLE_NVHPC_WORKAROUNDS)
   if(${ENABLE_OPENACC})
      target_compile_options(${PROJECT_NAME} PRIVATE
      $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=math_uniform,cuda12.2,cc80>
      $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs -Minfo=all -acc=verystrict -gpu=cc80>)
      set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc=verystrict -Mnofma -gpu=math_uniform,cuda12.2,cc80")
      set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc=verystrict -gpu=cc80")
      if(${DISABLE_OPENACC_ATOMICS})
         target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_OPENACC_ATOMICS)
      endif()
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE 
      $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma>
      $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs>)
   endif()
endif()

target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${YAC_Fortran_INCLUDE_DIRECTORIES} ${YAXT_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PRIVATE ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})


# OASIS-YAK
target_link_libraries(${PROJECT_NAME} PRIVATE ${YACUTILS_Fortran_LIBRARIES} ${YACCORE_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PRIVATE ${YAXT_Fortran_LIBRARIES} ${YAXTC_Fortran_LIBRARIES})


#set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

#target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${YAC_Fortran_INCLUDE_DIRECTORIES} ${YAXT_Fortran_INCLUDE_DIRECTORIES})
##target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES} ${OASIS_Fortran_LIBRARIES})

## OASIS-YAK
#target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_C ${YACUTILS_Fortran_LIBRARIES} ${YACCORE_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_C ${YAXT_Fortran_LIBRARIES} ${YAXTC_Fortran_LIBRARIES})

   
# fesom.x executable
add_executable(${PROJECT_NAME}.x ${src_home}/fesom_main.F90)
target_link_libraries(${PROJECT_NAME}.x PUBLIC ${PROJECT_NAME})

# fesom_meshdiag executable (optional)
if(BUILD_MESHDIAG)
    add_executable(fesom_meshdiag ${src_home}/fesom_meshdiag.F90)
    target_link_libraries(fesom_meshdiag PUBLIC ${PROJECT_NAME})
    target_link_libraries(fesom_meshdiag PRIVATE MPI::MPI_Fortran)
    
    # Install the mesh diagnostics utility
    install(TARGETS fesom_meshdiag
            RUNTIME DESTINATION bin
            COMPONENT tools)
    
    message(STATUS "Mesh diagnostics utility will be built as fesom_meshdiag")
endif()


get_target_property(FLAGS ${PROJECT_NAME} COMPILE_OPTIONS)
message(STATUS " --> Final Compile options for ${PROJECT_NAME}: ${FLAGS}")


### Export and installation

# Add meshdiag to additional targets if built
if(BUILD_MESHDIAG)
    set(additional_targets ${additional_targets} fesom_meshdiag)
endif()

fesom_export(TARGETS ${PROJECT_NAME} fesom.x ${additional_targets})
