cmake_minimum_required(VERSION 3.4)

project(fesom C Fortran)

option(DISABLE_MULTITHREADING "disable asynchronous operations" OFF)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
if(${USE_ICEPACK})
   file(GLOB sources_Fortran ${src_home}/*.F90
                             ${src_home}/icepack_drivers/*.F90
                             ${src_home}/icepack_drivers/Icepack/columnphysics/*.F90)
else()
   file(GLOB sources_Fortran ${src_home}/*.F90)
endif()
#list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
file(GLOB sources_C ${src_home}/*.c)

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

add_library(${PROJECT_NAME}_C ${sources_C})
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
if(${USE_ICEPACK})
   target_compile_definitions(${PROJECT_NAME} PRIVATE __icepack)
endif()
if(${VERBOSE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE VERBOSE)
endif()
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -init=zero -no-wrap-margin)
#	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fimf-use-svml -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit ) #-ftrapuv ) #-init=zero)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
		target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
	endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -hnoomp)
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} async_threads_cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom.x" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
