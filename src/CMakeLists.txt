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

if(DEFINED ENV{ENABLE_ALEPH_CRAYMPICH_WORKAROUNDS}) # be able to set the initial cache value from our env settings for aleph, not only via cmake command
   option(ALEPH_CRAYMPICH_WORKAROUNDS "workaround for performance issues on aleph" ON)
else()
   option(ALEPH_CRAYMPICH_WORKAROUNDS "workaround for performance issues on aleph" OFF)
endif()
if(DEFINED ENV{ENABLE_ALBEDO_INTELMPI_WORKAROUNDS}) # be able to set the initial cache value from our env settings for aleph, not only via cmake command
   option(ALBEDO_INTELMPI_WORKAROUNDS "workaround for performance issues on albedo" ON)
else()
   option(ALBEDO_INTELMPI_WORKAROUNDS "workaround for performance issues on albedo" OFF)
endif()

#TODO: these machine specific changes must be done using cmake toolchain then here otherwise will make cmake another raps. 

if(ALEPH_CRAYMPICH_WORKAROUNDS)
   # todo: enable these options only for our targets
   add_compile_options(-craympich-mt) # alternative cray-mpich library, about 5 % faster with cray-mpich/7.7.3 on aleph, not available for modules cray-mpich > 7.7.3; todo: test compiling and performance with cray-mpich > 7.7.3
   # make sure to also set these variables in the runtime environment:
   #    MPICH_MAX_THREAD_SAFETY=multiple # allows highest MPI thread level (i.e. MPI_THREAD_MULTIPLE)
   #    MPICH_CRAY_OPT_THREAD_SYNC=0 # the Cray MPICH library falls back to using the pthread_mutex-based thread-synchronization implementation
   #    MPICH_OPT_THREAD_SYNC=0 # seems to be a duplicate variable which also appears in some documentation instead of MPICH_CRAY_OPT_THREAD_SYNC (but this one brings a huge speed gain on aleph)
   #add_compile_options(-DDISABLE_PARALLEL_RESTART_READ) # reading restarts is slow when doing it on parallel on aleph, switch it off for now
   add_compile_options(-DENABLE_ALEPH_CRAYMPICH_WORKAROUNDS)
endif()
if(ALBEDO_INTELMPI_WORKAROUNDS)
   add_compile_options(-DENABLE_ALBEDO_INTELMPI_WORKAROUNDS)
   add_compile_options(-DDISABLE_PARALLEL_RESTART_READ)
endif()

#TODO: this machine specific block can easyly go in  cmake toolchain.
if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
   message(STATUS "multithreading disabled for Levante") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  leo-dcgp )
   message(STATUS "multithreading disabled for Leonardo") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  mn5-gpp )
   message(STATUS "multithreading disabled for MN5 GPP") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  albedo )
   message(STATUS "multithreading disabled for Albedo") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL atosecmwf )
   message(STATUS "multithreading disabled for Atos ECMWF")
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
else()
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
endif()

option(ENABLE_OPENACC "compile with OpenACC support" OFF)
message(STATUS "ENABLE_OPENACC: ${ENABLE_OPENACC}")

set(NV_GPU_ARCH "cc80" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")

option(ENABLE_OPENMP "build FESOM with OpenMP" OFF)
message(STATUS "ENABLE_OPENMP: ${ENABLE_OPENMP}")
if(ENABLE_OPENMP)
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
file(GLOB sources_Fortran ${src_home}/*.F90)
if(USE_ICEPACK)
   file(GLOB sources_ice_pack ${src_home}/icepack_drivers/*.F90
                              ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
   list(APPEND sources_Fortran ${sources_ice_pack})
elseif(ENABLE_IFS_INTERFACE) # ICEPACK + IFS_INTERFACE NOT SUPPORTED (YET)
   if(ENABLE_MULTIO)
      list(APPEND sources_Fortran ${src_home}/ifs_interface/iom.F90)
   endif()

   list(APPEND sources_Fortran ${src_home}/ifs_interface/ifs_interface.F90
                               ${src_home}/ifs_interface/ifs_modules.F90
                               ${src_home}/ifs_interface/ifs_notused.F90
                               ${src_home}/ifs_interface/mpp_io.F90)
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
target_compile_definitions(${PROJECT_NAME} PRIVATE -DMETIS_VERSION=5 -DPART_WEIGHTED -DMETISRANDOMSEED=35243)

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:module/fesom>)
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_Fortran)

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
if(ENABLE_OPENMP)
   target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(ENABLE_MULTIO)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __MULTIO)
   target_link_libraries(${PROJECT_NAME} PRIVATE multio-fapi)
endif()

if(DISABLE_MULTITHREADING)
   target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_MULTITHREADING)
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

if(ENABLE_IFS_INTERFACE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __ifsinterface)
endif()

if(VERBOSE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()

if(OPENMP_REPRODUCIBLE)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __openmp_reproducible)
endif()

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
   target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -init=zero -no-wrap-margin -fpe0) # add -fpe0 for RAPS environment
   if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL leo-dcgp )
      target_compile_options(${PROJECT_NAME} PRIVATE -O3 -xCORE-AVX512 -qopt-zmm-usage=high -align array64byte -ipo)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL mn5-gpp )
      target_compile_options(${PROJECT_NAME} PRIVATE -O3 -xCORE-AVX512 -qopt-zmm-usage=high -align array64byte -ipo)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  albedo)
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -O3 -ip -fPIC -qopt-malloc-options=2 -qopt-prefetch=5 -unroll-aggressive) #NEC mpi option
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL atosecmwf )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -xHost)
   endif()
#    target_compile_options(${PROJECT_NAME} PRIVATE -g -traceback ) #-check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)
#    target_compile_options(${PROJECT_NAME} PRIVATE -qopenmp -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)
#    target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
#    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
   target_compile_options(${PROJECT_NAME} PRIVATE -O2 -g -ffloat-store -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
   if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
      target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
   endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
	#target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O3)
   target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O2 -hnoacc -M878)  #-hnoacc is a workaround for cray automatically activate -hacc, -M878 is to suppress ftn-878 warning
   if(ENABLE_OPENMP)
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
   target_compile_options(${PROJECT_NAME} PRIVATE -fast -fastsse -O3 -Mallocatable=95 -Mr8 -pgf90libs)
   if(ENABLE_OPENACC)
      # additional compiler settings
      target_compile_options(${PROJECT_NAME} PRIVATE -acc -ta=tesla:${NV_GPU_ARCH} -Minfo=accel)
      set(CMAKE_EXE_LINKER_FLAGS "-acc -ta=tesla:${NV_GPU_ARCH}")
   endif()
   if(ENABLE_OPENMP)
      target_compile_options(${PROJECT_NAME} PRIVATE -Mipa=fast)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -Mipa=fast,inline)
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

   
# fesom.x executable
add_executable(${PROJECT_NAME}.x ${src_home}/fesom_main.F90)
target_link_libraries(${PROJECT_NAME}.x PUBLIC ${PROJECT_NAME})


### Export and installation

fesom_export(TARGETS ${PROJECT_NAME} fesom.x)
