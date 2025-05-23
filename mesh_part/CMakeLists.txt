cmake_minimum_required(VERSION 3.4)

project(fesom_ini C Fortran)
find_package(MPI REQUIRED)



# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}/../src)
set(sources_Fortran ${src_home}/MOD_MESH.F90 ${src_home}/oce_modules.F90 ${src_home}/gen_modules_config.F90 ${src_home}/gen_modules_partitioning.F90 ${src_home}/gen_modules_rotate_grid.F90 ${src_home}/fvom_init.F90 ${src_home}/oce_local.F90 ${src_home}/gen_comm.F90 ${src_home}/MOD_READ_BINARY_ARRAYS.F90 ${src_home}/MOD_WRITE_BINARY_ARRAYS.F90 ${src_home}/MOD_PARTIT.F90)
set(sources_C ${src_home}/fort_part.c)


add_subdirectory(${src_home}/../lib/metis-5.1.0             ${PROJECT_BINARY_DIR}/metis)
include_directories(${src_home}/../lib/metis-5.1.0/include)

add_library(${PROJECT_NAME}_C ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(${PROJECT_NAME}_C metis)

add_executable(${PROJECT_NAME} ${sources_Fortran})
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip)
#	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -g -traceback -check all,noarg_temp_created,bounds,uninit) 
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	target_compile_options(${PROJECT_NAME} PRIVATE -fdefault-real-8 -ffree-line-length-none)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
      target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
   endif()
endif()
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C MPI::MPI_Fortran)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
