cmake_minimum_required(VERSION 3.4)

project(parms C)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB all_sources ${src_home}/src/*.c ${src_home}/src/DDPQ/*.c)

find_package(BLAS) # standard way to find blas
if( NOT BLAS_FOUND) # try older way to find blas
   include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindBLAS.cmake")
endif()

# create our library (set its name to name of this project)
if(${BUILD_FESOM_AS_LIBRARY})
  add_library(${PROJECT_NAME} ${all_sources})
else()
  add_library(${PROJECT_NAME} ${all_sources})
endif()

if(${BLAS_FOUND})
  message("BLAS FOUND ---------------------------")
  target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS USE_MPI REAL=double DBL FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 HAS_BLAS)
  target_link_libraries(${PROJECT_NAME} INTERFACE ${BLAS_C_LIBRARIES} $ENV{UBUNTU_BLAS_LIBRARY})
else()
  message("BLAS NOT FOUND ***********************")
  target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS USE_MPI REAL=double DBL FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8)
  target_link_libraries(${PROJECT_NAME} INTERFACE $ENV{UBUNTU_BLAS_LIBRARY})
endif()

target_include_directories(${PROJECT_NAME}
	PRIVATE ${src_home}/src/../include ${src_home}/src/include
	INTERFACE ${src_home}/src/../include
)
if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
   target_compile_options(${PROJECT_NAME} PRIVATE -no-prec-div -no-prec-sqrt -fast-transcendentals -fp-model precise)
   
   if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
      target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
   endif()
endif()
if(${BUILD_FESOM_AS_LIBRARY})
  target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)
  install(TARGETS ${PROJECT_NAME} DESTINATION "${FESOM_INSTALL_PREFIX}/lib")
endif()
