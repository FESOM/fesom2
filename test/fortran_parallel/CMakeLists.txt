cmake_minimum_required(VERSION 3.4)

project(pfunittests_mpi.x Fortran)

set(PFUNIT_MPI_DIR CACHE PATH "path to the parallel pfunit installation")
if(NOT PFUNIT_MPI_DIR)
   message(FATAL_ERROR "PFUNIT_MPI_DIR not set")
endif()

# create the corresponding F90 file for all our pfunit files (*.pf) via pFUnitParser.py
file(GLOB pfunitparser_files ${CMAKE_CURRENT_LIST_DIR}/*.pf)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_CURRENT_LIST_DIR})
foreach(f ${pfunitparser_files})
   set(outfile ${CMAKE_CURRENT_BINARY_DIR}${f}.F90)
   add_custom_command(
      OUTPUT ${outfile}
      COMMAND python ${PFUNIT_MPI_DIR}/bin/pFUnitParser.py ${f} ${outfile}
      DEPENDS ${f}
   )
   list(APPEND sources_Fortran ${outfile})
endforeach()

file(GLOB unittest_files ${CMAKE_CURRENT_LIST_DIR}/*.F90)
list(APPEND sources_Fortran ${unittest_files})

list(APPEND sources_Fortran ${PFUNIT_MPI_DIR}/include/driver.F90)
find_library(PFUNIT_Fortran_LIBRARIES pfunit HINTS ${PFUNIT_MPI_DIR}/lib)

set(sources_under_test ${CMAKE_CURRENT_LIST_DIR}/../../fesom_cpl/toposervice.F90)
list(APPEND sources_Fortran ${sources_under_test})

add_executable(${PROJECT_NAME} ${sources_Fortran})
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
   target_compile_options(${PROJECT_NAME} PRIVATE -fpp -fp-model precise)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
   target_compile_options(${PROJECT_NAME} PRIVATE -cpp -ffree-line-length-none)
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR} ${PFUNIT_MPI_DIR}/mod)
target_link_libraries(${PROJECT_NAME} ${PFUNIT_Fortran_LIBRARIES})
target_compile_definitions(${PROJECT_NAME} PRIVATE USE_MPI)
