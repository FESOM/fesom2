cmake_minimum_required(VERSION 3.16)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "setting default build type: Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

project(FESOM2.0)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Default to using shared libs")
set(TOPLEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(FESOM_COUPLED OFF CACHE BOOL "compile fesom standalone or with oasis support (i.e. coupled)")
set(OIFS_COUPLED OFF CACHE BOOL "compile fesom coupled to OpenIFS. (Also needs FESOM_COUPLED to work)")
set(CRAY OFF CACHE BOOL "compile with cray ftn")
set(USE_ICEPACK ON CACHE BOOL "compile fesom with the Iceapck modules for sea ice column physics.")
set(OPENMP_REPRODUCIBLE OFF CACHE BOOL "serialize OpenMP loops that are critical for reproducible results")
set(RECOM_COUPLED OFF CACHE BOOL "compile fesom including biogeochemistry, REcoM3")
set(CISO_COUPLED OFF CACHE BOOL "compile ciso coupled to REcoM3. RECOM_COUPLED has to be active")
set(USE_MULTIO OFF CACHE BOOL "Use MULTIO for IO, either grib or binary for now. This also means path to MULTIO installation has to provided using env MULTIO_INSTALL_PATH='..' and multio configuration yamls must be present to run the model with MULTIO")
set(OASIS_WITH_YAC OFF CACHE BOOL "Useing a version of OASIS compiled with YAC instead of SCRIP for interpolation?")
set(ASYNC_ICEBERGS ON CACHE BOOL "compile fesom with or without support for asynchronous iceberg computations")
set(VERBOSE OFF CACHE BOOL "toggle debug output")
#add_subdirectory(oasis3-mct/lib/psmile)
add_subdirectory(src)

foreach( _file fesom-config.cmake fesom-config-version.cmake fesom-targets.cmake )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_BINARY_DIR}/src/${_file} ${PROJECT_BINARY_DIR}/${_file}  )
endforeach()

# Define ${PROJECT_NAME}_DIR in PARENT_SCOPE so that a `find_package( <this-project> )` in a bundle
# will easily find the project without requiring a `HINT <this-project>_BINARY_DIR` argument
if( NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR )
    # Guard needed because PARENT_SCOPE cannot be used in top-level CMake project

    set( fesom_DIR ${fesom_DIR} PARENT_SCOPE )
endif()
