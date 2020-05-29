cmake_minimum_required(VERSION 3.16.5)

project(fesom2_fortran_tests_parallel Fortran)
set(LIB_TARGET ${PROJECT_NAME}-lib)

find_package(PFUNIT REQUIRED)
enable_testing()

add_library(${LIB_TARGET} ${CMAKE_CURRENT_LIST_DIR}/../../src/mpi_topology_module.F90)
target_include_directories(${LIB_TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR} ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${LIB_TARGET} ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES})
set_target_properties(${LIB_TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_pfunit_ctest (${PROJECT_NAME}
  TEST_SOURCES mpi_topology_module_tests.pf # must be a path relative to CMAKE_CURRENT_SOURCE_DIR, then the generated files will be placed in CMAKE_CURRENT_BINARY_DIR (see add_pfunit_ctest.cmake)
  LINK_LIBRARIES ${LIB_TARGET}
  MAX_PES 6
  )

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
     target_compile_options(${LIB_TARGET} PRIVATE -fpp -fp-model precise)
   target_compile_options(${PROJECT_NAME} PRIVATE -fpp -fp-model precise)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
     target_compile_options(${LIB_TARGET} PRIVATE -cpp -ffree-line-length-none)
   target_compile_options(${PROJECT_NAME} PRIVATE -cpp -ffree-line-length-none)
endif()
