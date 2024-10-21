cmake_minimum_required(VERSION 3.16)

project(async_threads_cpp CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

set (CMAKE_CXX_STANDARD 11)

# get our source files
file(GLOB sources ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(ThreadsManagerFCMacros.h MACRO_NAMESPACE "ThreadsManagerFCMacros_" SYMBOLS init_ccall begin_ccall end_ccall)

add_library(${PROJECT_NAME} OBJECT ${sources})
target_include_directories(${PROJECT_NAME}
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
