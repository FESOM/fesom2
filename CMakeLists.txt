cmake_minimum_required(VERSION 3.16)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "setting default build type: Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

project(FESOM2.0)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

set(BUILD_SHARED_LIBS ON CACHE BOOL "Default to using shared libs")
set(TOPLEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(FESOM_COUPLED OFF CACHE BOOL "compile fesom standalone or with oasis support (i.e. coupled)")
set(OIFS_COUPLED OFF CACHE BOOL "compile fesom coupled to OpenIFS. (Also needs FESOM_COUPLED to work)")
set(CRAY OFF CACHE BOOL "compile with cray ftn")
set(USE_ICEPACK OFF CACHE BOOL "compile fesom with the Iceapck modules for sea ice column physics.")
set(OPENMP_REPRODUCIBLE OFF CACHE BOOL "serialize OpenMP loops that are critical for reproducible results")
set(RECOM_COUPLED OFF CACHE BOOL "compile fesom including biogeochemistry, REcoM3")
set(CISO_COUPLED OFF CACHE BOOL "compile ciso coupled to REcoM3. RECOM_COUPLED has to be active")
set(USE_MULTIO OFF CACHE BOOL "Use MULTIO for IO, either grib or binary for now. This also means path to MULTIO installation has to provided using env MULTIO_INSTALL_PATH='..' and multio configuration yamls must be present to run the model with MULTIO")
set(OASIS_WITH_YAC OFF CACHE BOOL "Useing a version of OASIS compiled with YAC instead of SCRIP for interpolation?")
set(ASYNC_ICEBERGS ON CACHE BOOL "compile fesom with or without support for asynchronous iceberg computations")
set(VERBOSE OFF CACHE BOOL "toggle debug output")
set(FESOM_PROFILING OFF CACHE BOOL "Enable enhanced profiling system")

# Testing options
option(BUILD_TESTING "Build tests" OFF)

# Tools
option(BUILD_MESHPARTITIONER "Build mesh partitioning tools (fesom_meshpart executable)" OFF)
option(BUILD_MESHDIAG "Build mesh diagnostics utility (fesom_meshdiag executable)" OFF)

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
