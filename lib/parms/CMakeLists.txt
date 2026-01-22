cmake_minimum_required(VERSION 3.4)

project(parms C)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB all_sources ${src_home}/src/*.c ${src_home}/src/DDPQ/*.c)

include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindBLAS.cmake")

# create our library (set its name to name of this project)
add_library(${PROJECT_NAME} ${all_sources})
target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS USE_MPI REAL=double DBL FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 HAS_BLAS)
target_include_directories(${PROJECT_NAME}
	PRIVATE ${src_home}/src/../include ${src_home}/src/include
	INTERFACE ${src_home}/src/../include
)
target_link_libraries(${PROJECT_NAME} INTERFACE ${BLAS_C_LIBRARIES} $ENV{UBUNTU_BLAS_LIBRARY})
if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
	target_compile_options(${PROJECT_NAME} PRIVATE -no-prec-div -no-prec-sqrt -fast-transcendentals -fp-model precise)
endif()
