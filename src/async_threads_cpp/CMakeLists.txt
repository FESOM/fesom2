cmake_minimum_required(VERSION 3.4)

project(async_threads_cpp CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(ThreadsManagerFCMacros.h MACRO_NAMESPACE "ThreadsManagerFCMacros_" SYMBOLS init_ccall begin_ccall end_ccall)

add_library(${PROJECT_NAME} ${sources_CXX})
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_LIST_DIR}
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL Cray )
   target_compile_options(${PROJECT_NAME} PRIVATE -hstd=c++11)
else()
   target_compile_options(${PROJECT_NAME} PRIVATE -std=c++11)
endif()
