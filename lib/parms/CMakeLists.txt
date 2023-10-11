cmake_minimum_required(VERSION 3.4)

project(parms C)

# Set location to look for find_package modules
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../../cmake )

# Set -fPIC flag etc.
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB all_sources ${src_home}/src/*.c ${src_home}/src/DDPQ/*.c)

# TODO: support legacy BLAS without LAPACK?
#include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindBLAS.cmake")
if( NOT LAPACK_LIBRARIES )
  # First check if we have the Cray libraries
  set( _cray_libsci_loaded $ENV{CRAY_LIBSCI_DIR} )
  if( _cray_libsci_loaded )
    # set( _CRAY_PRGENV $ENV{PE_ENV} )
    # string( TOLOWER "${_CRAY_PRGENV}" _cray_prgenv  )
    # set( LAPACK_LIBRARIES sci_${_cray_prgenv} )
  else()
    find_package(LAPACK REQUIRED) # cannot compile without! or BLAS
  endif()
endif()
if( NOT TARGET MPI::MPI_C )
  find_package(MPI REQUIRED COMPONENTS C)
endif()


# create our library (set its name to name of this project)
add_library(${PROJECT_NAME} ${all_sources})

target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS USE_MPI REAL=double DBL FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8)

target_include_directories(${PROJECT_NAME}
  PRIVATE ${src_home}/src/include
  PUBLIC  $<BUILD_INTERFACE:${src_home}/src/../include>
)

#target_link_libraries(${PROJECT_NAME} PRIVATE ${BLAS_C_LIBRARIES} $ENV{UBUNTU_BLAS_LIBRARY})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_C)



if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
    target_compile_options(${PROJECT_NAME} PRIVATE -no-prec-div -no-prec-sqrt -fast-transcendentals -fp-model precise)
    if(${FESOM_PLATFORM_STRATEGY} STREQUAL  levante.dkrz.de )
        target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -mtune=core-avx2)
    endif()
endif()
