cmake_minimum_required(VERSION 3.16.5)

project(fesom2_fortran_tests Fortran CXX)
set(LIB_TARGET ${PROJECT_NAME}-lib)

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindNETCDF.cmake)

find_package(PFUNIT REQUIRED)
enable_testing()

add_library(${LIB_TARGET} ${CMAKE_CURRENT_LIST_DIR}/../../src/forcing_provider_async_module.F90
                          ${CMAKE_CURRENT_LIST_DIR}/../../src/async_threads_module.F90
                          ${CMAKE_CURRENT_LIST_DIR}/../../src/forcing_provider_netcdf_module.F90
                          ${CMAKE_CURRENT_LIST_DIR}/../../src/forcing_lookahead_reader_module.F90
           )

add_subdirectory(../../src/async_threads_cpp ${PROJECT_BINARY_DIR}/async_threads_cpp)

target_include_directories(${LIB_TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR} ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${LIB_TARGET} ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES})
target_link_libraries(${LIB_TARGET} async_threads_cpp)
set_target_properties(${LIB_TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_pfunit_ctest (${PROJECT_NAME}
  TEST_SOURCES forcing_provider_module_tests.pf # must be a path relative to CMAKE_CURRENT_SOURCE_DIR, then the generated files will be placed in CMAKE_CURRENT_BINARY_DIR (see add_pfunit_ctest.cmake)
               forcing_provider_netcdf_module_tests.pf
               forcing_lookahead_reader_module_tests.pf
  LINK_LIBRARIES ${LIB_TARGET}
  )

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
     target_compile_options(${LIB_TARGET} PRIVATE -fpp -fp-model precise)
   target_compile_options(${PROJECT_NAME} PRIVATE -fpp -fp-model precise)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
     target_compile_options(${LIB_TARGET} PRIVATE -cpp -ffree-line-length-none)
   target_compile_options(${PROJECT_NAME} PRIVATE -cpp -ffree-line-length-none)
endif()

add_custom_command(
   TARGET ${PROJECT_NAME} POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/
          ${CMAKE_CURRENT_BINARY_DIR}/fixtures/)
