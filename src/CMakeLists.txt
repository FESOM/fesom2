cmake_minimum_required(VERSION 3.9)

project(fesom C Fortran)

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

if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
   message(STATUS "multithreading disabled for Levante") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  albedo )
   message(STATUS "multithreading disabled for Albedo") # multithreading suddenly produces an error, disable it until a fix is found. issue #413
   option(DISABLE_MULTITHREADING "disable asynchronous operations" ON)
else()
   option(DISABLE_MULTITHREADING "disable asynchronous operations" OFF)
endif()

option(ENABLE_OPENACC "compile with OpenACC support" OFF)
option(DISABLE_OPENACC_ATOMICS "disable kernels using atomic statement for reproducible results" ON)
set(GPU_COMPUTE_CAPABILITY "cc80" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")
set(GPU_FLAGS "cuda11.7,${GPU_COMPUTE_CAPABILITY}" CACHE STRING "GPU arch for nvfortran compiler (cc35,cc50,cc60,cc70,cc80,...)")

option(ENABLE_OPENMP "build FESOM with OpenMP" OFF)
if(${ENABLE_OPENMP})
   find_package(OpenMP REQUIRED)
endif()

# option to trigger building a library version of FESOM
# we do not always build the library along with the executable to avoid having two targets here in the CMakeLists.txt
# two targets would allow e.g. setting different compiler options or preprocessor definition, which would be error prone.
option(BUILD_FESOM_AS_LIBRARY "Build a library instead of an executable" OFF)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
if(${RECOM_COUPLED})
		file(GLOB sources_Fortran ${src_home}/*.F90 
                                          ${src_home}/int_recom/*.F90 
                                          ${src_home}/int_recom/recom/*.F90)

elseif(${USE_ICEPACK})
   file(GLOB sources_Fortran ${src_home}/*.F90
                             ${src_home}/icepack_drivers/*.F90
                             ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
else()
   file(GLOB sources_Fortran ${src_home}/*.F90)
   list(REMOVE_ITEM sources_Fortran ${src_home}/recom_nitrogenss.F90)    # REcoM
   list(REMOVE_ITEM sources_Fortran ${src_home}/recom_sinking.F90)       # REcoM
   list(REMOVE_ITEM sources_Fortran ${src_home}/recom_modules_diag.F90)  # REcoM
endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
file(GLOB sources_C ${src_home}/*.c)
list(REMOVE_ITEM sources_C ${src_home}/psolve_feom.c) # does the file still exist?

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
set(FESOM_ORIGINAL_VERSION_FILE ${src_home}/fesom_version_info.F90)
set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
list(REMOVE_ITEM sources_Fortran ${FESOM_ORIGINAL_VERSION_FILE}) # we want to compile the generated file instead
list(APPEND sources_Fortran ${FESOM_GENERATED_VERSION_FILE})
add_custom_command(OUTPUT 5303B6F4_E4F4_45B2_A6E5_8E2B9FB5CDC4 ${FESOM_GENERATED_VERSION_FILE} # the first arg to OUTPUT is a name for a file we never create to make sure this command will run on every re-build (let our file be the second arg, as the first file is inadvertently removed by make)
                   COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                   COMMENT "determining ${PROJECT_NAME} git SHA ...")

#if(${FESOM_STANDALONE})
#   list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_driver.F90)
#endif()
list(REMOVE_ITEM sources_Fortran ${src_home}/fvom_init.F90)
list(REMOVE_ITEM sources_C       ${src_home}/fort_part.c)

# depends on the metis library
#add_subdirectory(../lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
#include_directories(../lib/metis-5.1.0/include)
# depends on the parms library
add_subdirectory(../lib/parms ${PROJECT_BINARY_DIR}/parms)

add_subdirectory(async_threads_cpp)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)

if(${BUILD_FESOM_AS_LIBRARY})
  add_library(${PROJECT_NAME}_C ${sources_C})
else()
add_library(${PROJECT_NAME}_C ${sources_C})
endif()
target_compile_definitions(${PROJECT_NAME}_C PRIVATE PARMS USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(${PROJECT_NAME}_C parms) #metis

# create our binary (set its name to name of this project)
add_executable(${PROJECT_NAME} ${sources_Fortran})
target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS -DMETIS_VERSION=5 -DPART_WEIGHTED -DMETISRANDOMSEED=35243)
if(${DISABLE_MULTITHREADING})
   target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_MULTITHREADING)
endif()
if(${FESOM_COUPLED})
   include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __oasis)
endif()
if(${OIFS_COUPLED})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __oifs)
endif()
if(${RECOM_COUPLED})
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __recom USE_PRECISION=2)
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __recom USE_PRECISION=2 __coccos)
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __recom USE_PRECISION=2 __3Zoo2Det)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __recom USE_PRECISION=2 __coccos __3Zoo2Det)
endif()
if(${CISO_COUPLED})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __ciso)
endif()
if(${USE_ICEPACK})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __icepack)
endif()
if(${VERBOSE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()
if(${OPENMP_REPRODUCIBLE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __openmp_reproducible)
endif()
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
#if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
#	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero -no-wrap-margin)
##	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-#ftrapuv ) #-init=zero)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
   if(${BUILD_FESOM_AS_LIBRARY})
        target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -init=zero -no-wrap-margin -fpe0) # add -fpe0 for RAPS environment
   else()
        target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -init=zero -no-wrap-margin)
   endif()
   if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   elseif(${FESOM_PLATFORM_STRATEGY} STREQUAL  albedo)
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -O3 -ip -fPIC -qopt-malloc-options=2 -qopt-prefetch=5 -unroll-aggressive) #NEC mpi option
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
   if(${ENABLE_OPENMP})
      target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -homp)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -hnoomp)
   endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC )
   target_compile_definitions(${PROJECT_NAME} PRIVATE ENABLE_NVHPC_WORKAROUNDS)
   #target_compile_options(${PROJECT_NAME} PRIVATE -fast -fastsse -O3 -Mallocatable=95 -Mr8 -pgf90libs)
   target_compile_options(${PROJECT_NAME} PRIVATE -Mnofma -Mallocatable=95 -Mr8 -pgf90libs)
   if(${ENABLE_OPENACC})
      # additional compiler settings
      target_compile_options(${PROJECT_NAME} PRIVATE -acc -gpu=${GPU_FLAGS} -Minfo=accel)
      set(CMAKE_EXE_LINKER_FLAGS "-acc -gpu=${GPU_FLAGS}")
      if(${DISABLE_OPENACC_ATOMICS})
          target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_OPENACC_ATOMICS)
      endif()
   endif()
   if(${ENABLE_OPENMP})
      target_compile_options(${PROJECT_NAME} PRIVATE -Mipa=fast)
      #else()
      ## Mipa enables -O2, so it is disabled during gpu development for testing
       #target_compile_options(${PROJECT_NAME} PRIVATE -Mipa=fast,inline)
   endif()
endif()
if(${BUILD_FESOM_AS_LIBRARY})
  target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
  target_compile_options(${PROJECT_NAME}_C PRIVATE -fPIC)
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} async_threads_cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
if(${ENABLE_OPENMP} AND NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
   target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_Fortran_FLAGS}) # currently we only have OpenMP in the Fortran part
   target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_Fortran)
endif()

#set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom.x" CACHE FILEPATH "file path where the FESOM binary should be put")
#get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
#get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
#install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
set(FESOM_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/.." CACHE FILEPATH "directory where FESOM will be installed to via 'make install'")
if(${BUILD_FESOM_AS_LIBRARY})
   install(TARGETS ${PROJECT_NAME} DESTINATION "${FESOM_INSTALL_PREFIX}/lib")
   install(TARGETS ${PROJECT_NAME}_C DESTINATION "${FESOM_INSTALL_PREFIX}/lib")
else()
   set(FESOM_INSTALL_FILEPATH "${FESOM_INSTALL_PREFIX}/bin/fesom.x")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
endif()
