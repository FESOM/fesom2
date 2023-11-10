cmake_minimum_required(VERSION 3.4)

project(async_threads_cpp CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(ThreadsManagerFCMacros.h MACRO_NAMESPACE "ThreadsManagerFCMacros_" SYMBOLS init_ccall begin_ccall end_ccall)

add_library(${PROJECT_NAME} ${sources_CXX})

target_include_directories(${PROJECT_NAME}
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PUBLIC  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL Cray )
   target_compile_options(${PROJECT_NAME} PRIVATE -hstd=c++11)
else()
   target_compile_options(${PROJECT_NAME} PRIVATE -std=c++11)
endif()
if(${BUILD_FESOM_AS_LIBRARY})
  target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
  install(TARGETS ${PROJECT_NAME} DESTINATION "${FESOM_INSTALL_PREFIX}/lib")
endif()
