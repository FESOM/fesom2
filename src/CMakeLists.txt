cmake_minimum_required(VERSION 3.13) # 20-11-2018
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")

project(fesom C Fortran)

option(DISABLE_MULTITHREADING "disable asynchronous operations" OFF)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_Fortran ${src_home}/*.F90)
#list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
file(GLOB sources_C ${src_home}/*.c)

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
set(FESOM_ORIGINAL_VERSION_FILE ${src_home}/fesom_version_info.F90)
set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
list(REMOVE_ITEM sources_Fortran ${FESOM_ORIGINAL_VERSION_FILE}) # we want to compile the generated file instead
list(REMOVE_ITEM sources_Fortran ${FESOM_GENERATED_VERSION_FILE}) # we don't want to include twice the generated file (included by the glob *.F90)
list(APPEND sources_Fortran ${FESOM_GENERATED_VERSION_FILE})
add_custom_command(OUTPUT 5303B6F4_E4F4_45B2_A6E5_8E2B9FB5CDC4 ${FESOM_GENERATED_VERSION_FILE} # the first arg to OUTPUT is a name for a file we never create to make sure this command will run on every re-build (let our file be the second arg, as the first file is inadvertently removed by make)
                   COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                   COMMENT "determining fesom git SHA ...")

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

add_library(fesom_C ${sources_C})
target_compile_definitions(fesom_C PRIVATE PARMS USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(fesom_C parms) #metis

# create our binary (set its name to name of this project)
add_executable(fesom ${sources_Fortran})
target_compile_definitions(fesom PRIVATE
  PARMS
  -DMETIS_VERSION=5 -DPART_WEIGHTED -DMETISRANDOMSEED=35243
  $<$<BOOL:${DISABLE_MULTITHREADING}>:DISABLE_MULTITHREADING>
)
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_%3CLANG%3E_COMPILER_ID
target_compile_options(fesom PRIVATE
  $<$<BOOL:${FESOM_COUPLED}>:__oasis>
  $<$<BOOL:${OIFS_COUPLED}>:__oifs>
  $<$<BOOL:${VERBOSE}>:VERBOSE>

  # or use $<COMPILE_LANG_AND_ID:language,compiler_ids>
  $<$<Fortran_COMPILER_ID:Intel>:
    -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero -no-wrap-margin
   >

   $<$<Fortran_COMPILER_ID:GNU>:
    -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none
    $<$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},10>:
      -fallow-argument-mismatch # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
     >
   >

  $<$<Fortran_COMPILER_ID:Cray>:
    -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64
   >

  $<$<Fortran_COMPILER_ID:PGI>:
    -fast -fastsse -O3 -Mallocatable=95 -Minline -Mr8 -Munroll -pgf90libs
    $<$<BOOL:${FESOM_OPENACC}>:
      -acc -ta=tesla:cc35
      -Minfo=accel
     >
   >
  )

target_link_options(fesom PRIVATE
  $<$<Fortran_COMPILER_ID:PGI>:
    -acc -ta=tesla:cc35
  >
)
if(${FESOM_COUPLED})
   if(${CRAY})
      include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS_cray.cmake)
   else()
      if(${OIFS_COUPLED})
        include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS_oifs.cmake)
      else()
        include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
      endif()
   endif()
endif()
target_include_directories(fesom PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(fesom PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(fesom PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(fesom fesom_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
target_link_libraries(fesom fesom_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
target_link_libraries(fesom async_threads_cpp)
set_target_properties(fesom PROPERTIES LINKER_LANGUAGE Fortran)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom.x" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/fesom DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
