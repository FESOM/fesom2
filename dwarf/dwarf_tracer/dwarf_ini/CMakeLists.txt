cmake_minimum_required(VERSION 3.4)

project(fesom Fortran)

if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
   set(FESOM_PLATFORM_STRATEGY $ENV{FESOM_PLATFORM_STRATEGY} CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
else()
   set(FESOM_PLATFORM_STRATEGY "notset" CACHE STRING "switch to platform specific compile settings, this is usually determined via the env.sh script")
endif()

option(ENABLE_OPENACC "compile with OpenACC support" ON)
option(DISABLE_OPENACC_ATOMICS "disable kernels using atomic statement for reproducible results" OFF)
option(ENABLE_NVTX "enable use of nvtx tools" OFF)

option(ASYNCHRONOUS_IO_THREADS "enable asynchronous I/O threads" OFF)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
#if(${USE_ICEPACK})
#   file(GLOB sources_Fortran ${src_home}/*.F90
#                             ${src_home}/icepack_drivers/*.F90
#                             ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
#else()
file(GLOB sources_Fortran ${src_home}/*.F90)
#endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
#file(GLOB sources_C ${src_home}/*.c)

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
#set(FESOM_ORIGINAL_VERSION_FILE ${src_home}/fesom_version_info.F90)
#set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
#list(REMOVE_ITEM sources_Fortran ${FESOM_ORIGINAL_VERSION_FILE}) # we want to compile the generated file instead
#list(APPEND sources_Fortran ${FESOM_GENERATED_VERSION_FILE})
#add_custom_command(OUTPUT 5303B6F4_E4F4_45B2_A6E5_8E2B9FB5CDC4 ${FESOM_GENERATED_VERSION_FILE} # the first arg to OUTPUT is a name for a file we never create to make sure this command will run on every re-build (let our file be the second arg, as the first file is inadvertently removed by make)
#                   COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
#                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
#                   COMMENT "determining ${PROJECT_NAME} git SHA ...")

#if(${FESOM_STANDALONE})
#   list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_driver.F90)
#endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fvom_init.F90)
#list(REMOVE_ITEM sources_C       ${src_home}/fort_part.c)

# depends on the metis library
#add_subdirectory(../lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
#include_directories(../lib/metis-5.1.0/include)

#add_subdirectory(async_threads_cpp)

#include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)

#add_library(${PROJECT_NAME}_C ${sources_C})
#target_compile_definitions(${PROJECT_NAME}_C PRIVATE USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)

# create our binary (set its name to name of this project)
add_executable(${PROJECT_NAME} ${sources_Fortran})
#target_compile_definitions(${PROJECT_NAME} PRIVATE -DMETIS_VERSION=5 -DPART_WEIGHTED -DMETISRANDOMSEED=35243)
#if(ASYNCHRONOUS_IO_THREADS)
#   target_compile_definitions(${PROJECT_NAME} PRIVATE ASYNCHRONOUS_IO_THREADS)
#endif()
#if(${FESOM_COUPLED})
#   include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __oasis)
#endif()
#if(${OIFS_COUPLED})
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __oifs)
#endif()
#if(${USE_ICEPACK})
#   target_compile_definitions(${PROJECT_NAME} PRIVATE __icepack)
#endif()
if(${VERBOSE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
#       target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -init=zero -no-wrap-margin)
        target_compile_options(${PROJECT_NAME} PRIVATE 
                               $<$<CONFIG:DEBUG>:-r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit>
                               $<$<CONFIG:RELEASE>:-r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -ip -g -traceback>
                              )
   if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   else()
      target_compile_options(${PROJECT_NAME} PRIVATE -xHost)
   endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
        target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
                target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
        endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL NVHPC)
   if(${ENABLE_OPENACC})
       target_compile_options(${PROJECT_NAME} PRIVATE
                              $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma -Minfo=all -acc=verystrict -gpu=math_uniform,cuda11.7,cc80>
                              $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs -Minfo=all -acc=verystrict -gpu=cc80>
                             )
       set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-acc=verystrict -Mnofma -gpu=math_uniform,cuda11.7,cc80")
       set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-acc=verystrict -gpu=cc80")
       if(${DISABLE_OPENACC_ATOMICS})
          target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_OPENACC_ATOMICS)
       endif()
   else()
       target_compile_options(${PROJECT_NAME} PRIVATE 
                              $<$<CONFIG:DEBUG>:-Mallocatable=95 -Mr8 -pgf90libs -Mnofma>
                              $<$<CONFIG:RELEASE>:-Mallocatable=95 -Mr8 -pgf90libs>
                             )
   endif()
   if(${ENABLE_NVTX})
       target_link_libraries(${PROJECT_NAME} nvhpcwrapnvtx)
   endif()
endif()


#elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
#       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
#endif()
#target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} async_threads_cpp)
#set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom.x" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
