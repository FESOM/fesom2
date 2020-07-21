cmake_minimum_required(VERSION 3.4)

project(forcing_provider_async CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_CXX ${src_home}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(ForcingReaderFCMacros.h MACRO_NAMESPACE "ForcingReaderFCMacros_" SYMBOLS init_ccall begin_ccall end_ccall)

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
