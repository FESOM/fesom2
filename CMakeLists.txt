cmake_minimum_required(VERSION 3.16)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "setting default build type: Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

project(FESOM2.0 LANGUAGES C CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

set(BUILD_SHARED_LIBS ON CACHE BOOL "Default to using shared libs")
set(TOPLEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(FESOM_COUPLED       OFF CACHE BOOL "compile fesom standalone or with oasis support (i.e. coupled)")
set(OIFS_COUPLED        OFF CACHE BOOL "compile fesom coupled to OpenIFS. (Also needs FESOM_COUPLED to work)")
set(USE_YAC             OFF CACHE BOOL "compile fesom with yac")
set(CRAY                OFF CACHE BOOL "compile with cray ftn")
set(USE_ICEPACK         OFF CACHE BOOL "compile fesom with the Iceapck modules for sea ice column physics.")
set(OPENMP_REPRODUCIBLE OFF CACHE BOOL "serialize OpenMP loops that are critical for reproducible results")
set(RECOM_COUPLED       OFF CACHE BOOL "compile fesom including biogeochemistry, REcoM3")
set(CISO_COUPLED        OFF CACHE BOOL "compile ciso coupled to REcoM3. RECOM_COUPLED has to be active")
set(USE_MULTIO          OFF CACHE BOOL "Use MULTIO for IO, either grib or binary for now. This also means path to MULTIO installation has to provided using env MULTIO_INSTALL_PATH='..' and multio configuration yamls must be present to run the model with MULTIO")
set(OASIS_WITH_YAC      OFF CACHE BOOL "Useing a version of OASIS compiled with YAC instead of SCRIP for interpolation?")
set(ASYNC_ICEBERGS      ON  CACHE BOOL "compile fesom with or without support for asynchronous iceberg computations")
set(VERBOSE             OFF CACHE BOOL "toggle debug output")
set(FESOM_PROFILING     OFF CACHE BOOL "Enable enhanced profiling system")
set(CVMIX               ON  CACHE BOOL "Download/Compile/Install/Link CVMix libray")

# CVMix is a fixed double-precision library (its kind is cvmix_r8 = real64). It is
# linked in BOTH double- and single-precision FESOM builds: the FESOM wrappers in
# src/cvmix_driver/ convert WP<->cvmix_r8 at the CVMix call boundary, so CVMix always
# runs in double precision. This keeps the available mixing schemes identical in SP and
# DP (so SP and DP configurations can be compared like-for-like) and leaves the default
# DP build numerically unchanged. CVMIX can still be turned off explicitly (-DCVMIX=OFF).

# Testing options
option(BUILD_TESTING "Build tests" OFF)

# Tools
option(BUILD_MESHPARTITIONER "Build mesh partitioning tools (fesom_meshpart executable)" OFF)
option(BUILD_MESHDIAG "Build mesh diagnostics utility (fesom_meshdiag executable)" OFF)

# NetCDF build option
option(BUILD_NETCDF "Download and build NetCDF libraries when system libraries are incompatible" OFF)

# YAC build option
option(BUILD_YAC "Download and build YAC, YAXT, and libfyaml libraries from source" OFF)

# Include ExternalProject module for building NetCDF if requested
if(BUILD_NETCDF)
    include(ExternalProject)
    include(cmake/BuildNetCDF.cmake)
endif()

# Include ExternalProject module for building YAC if requested
if(BUILD_YAC)
    include(ExternalProject)
    include(cmake/BuildYAC.cmake)
endif()

# Compile-time diagnostics. "strict" stops suppressing argument mismatches, which
# gfortran >=10 then reports as errors.
set(FESOM_BUILD_CHECKS "off" CACHE STRING "Compile-time diagnostics: off | warn | strict")
set_property(CACHE FESOM_BUILD_CHECKS PROPERTY STRINGS off warn strict)

set(FESOM_BUILD_CHECK_FLAGS "")
set(FESOM_ALLOW_ARG_MISMATCH ON)
if(NOT FESOM_BUILD_CHECKS STREQUAL "off")
   if(NOT FESOM_BUILD_CHECKS MATCHES "^(warn|strict)$")
      message(FATAL_ERROR "FESOM_BUILD_CHECKS must be one of: off, warn, strict (got '${FESOM_BUILD_CHECKS}')")
   endif()
   if(FESOM_BUILD_CHECKS STREQUAL "strict")
      set(FESOM_ALLOW_ARG_MISMATCH OFF)
   endif()

   if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      set(FESOM_BUILD_CHECK_FLAGS -Wall -Wextra -Wno-unused-dummy-argument -Wno-compare-reals)
   elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel OR CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
      set(FESOM_BUILD_CHECK_FLAGS -warn all,noexternals)
   else()
      message(WARNING "FESOM_BUILD_CHECKS=${FESOM_BUILD_CHECKS}: no flag mapping for compiler ${CMAKE_Fortran_COMPILER_ID}.")
   endif()
   message(STATUS "FESOM_BUILD_CHECKS=${FESOM_BUILD_CHECKS} (allow argument mismatch: ${FESOM_ALLOW_ARG_MISMATCH})")
endif()

# Runtime diagnostics for development builds. Appended after the production flags,
# so they win on -O and -finit-*; must be set here to reach mesh_part too.
set(FESOM_RUNTIME_CHECKS "off" CACHE STRING "Runtime diagnostics for development builds: off | checks | paranoid")
set_property(CACHE FESOM_RUNTIME_CHECKS PROPERTY STRINGS off checks paranoid)

set(FESOM_RUNTIME_CHECK_FLAGS "")
if(NOT FESOM_RUNTIME_CHECKS STREQUAL "off")
   if(NOT FESOM_RUNTIME_CHECKS MATCHES "^(checks|paranoid)$")
      message(FATAL_ERROR "FESOM_RUNTIME_CHECKS must be one of: off, checks, paranoid (got '${FESOM_RUNTIME_CHECKS}')")
   endif()

   if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      if(FESOM_RUNTIME_CHECKS STREQUAL "paranoid")
         set(FESOM_RUNTIME_CHECK_FLAGS -fcheck=all -fbacktrace
                                       -ffpe-trap=invalid,zero,overflow
                                       -finit-real=snan -finit-integer=-9999)
      else()
         set(FESOM_RUNTIME_CHECK_FLAGS -fcheck=bounds -fbacktrace)
      endif()
   elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel OR CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
      if(FESOM_RUNTIME_CHECKS STREQUAL "paranoid")
         set(FESOM_RUNTIME_CHECK_FLAGS -check all -traceback -fpe0 -init=snan,arrays)
      else()
         set(FESOM_RUNTIME_CHECK_FLAGS -check bounds -traceback)
      endif()
   else()
      message(WARNING "FESOM_RUNTIME_CHECKS=${FESOM_RUNTIME_CHECKS}: no flag mapping for compiler ${CMAKE_Fortran_COMPILER_ID}; building without runtime checks.")
   endif()

   if(FESOM_RUNTIME_CHECK_FLAGS)
      string(REPLACE ";" " " _fesom_checks_pretty "${FESOM_RUNTIME_CHECK_FLAGS}")
      message(STATUS "FESOM_RUNTIME_CHECKS=${FESOM_RUNTIME_CHECKS} -> ${_fesom_checks_pretty}")
   endif()
endif()

#add_subdirectory(oasis3-mct/lib/psmile)
add_subdirectory(src)

# Enable testing if requested
if(BUILD_TESTING)
    enable_testing()
    # Only add tests if the tests directory exists
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
        add_subdirectory(tests)
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
        # Fallback to existing test directory if tests/ doesn't exist yet
        message(STATUS "Using existing test/ directory for testing")
        add_subdirectory(test)
    else()
        message(WARNING "BUILD_TESTING is ON but no tests directory found")
    endif()
endif()

# Add mesh_part directory if requested
if(BUILD_MESHPARTITIONER)
    message(STATUS "Configuring mesh partitioner as part of the main build")
    add_subdirectory(mesh_part)
    
    # Install the mesh partitioner
    install(TARGETS fesom_meshpart
            RUNTIME DESTINATION bin
            COMPONENT tools)
    
    # Create a convenience target for building just the mesh partitioner
    add_custom_target(meshpart
        COMMENT "Building just the mesh partitioner"
        DEPENDS fesom_meshpart)
    
    message(STATUS "Mesh partitioner will be built as part of the main build")
endif()

foreach( _file fesom-config.cmake fesom-config-version.cmake fesom-targets.cmake )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_BINARY_DIR}/src/${_file} ${PROJECT_BINARY_DIR}/${_file}  )
endforeach()

# Add uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

# Define ${PROJECT_NAME}_DIR in PARENT_SCOPE so that a `find_package( <this-project> )` in a bundle
# will easily find the project without requiring a `HINT <this-project>_BINARY_DIR` argument
if( NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR )
    # Guard needed because PARENT_SCOPE cannot be used in top-level CMake project

    set( fesom_DIR ${fesom_DIR} PARENT_SCOPE )
endif()
