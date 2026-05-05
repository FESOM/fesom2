# BuildYAC.cmake
# External project build for YAC (Yet Another Coupler)
# YAC is used for coupling FESOM2 with atmospheric models

include(ExternalProject)
include(GNUInstallDirs)

# Include dependencies
include(${CMAKE_CURRENT_LIST_DIR}/BuildLibfyaml.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/BuildYAXT.cmake)

# Set installation prefix for built YAC libraries
set(YAC_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/yac)

# Create directories early so CMake doesn't complain about non-existent paths
file(MAKE_DIRECTORY ${YAC_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${YAC_INSTALL_PREFIX}/lib)
file(MAKE_DIRECTORY ${YAC_INSTALL_PREFIX}/include)

# Find required dependencies
find_package(MPI REQUIRED)

# Use FESOM's FindNETCDF.cmake to detect NetCDF the same way FESOM does
# This ensures consistency and reuses the same NetCDF library
if(NOT NETCDF_C_FOUND OR NOT NETCDF_Fortran_FOUND)
    message(STATUS "YAC needs NetCDF - using FESOM's NetCDF detection...")
    include(${CMAKE_CURRENT_LIST_DIR}/FindNETCDF.cmake)
endif()

# Verify NetCDF is now available
if(NOT NETCDF_C_FOUND OR NOT NETCDF_Fortran_FOUND)
    message(FATAL_ERROR "YAC requires NetCDF, but NetCDF was not found. "
                        "Please ensure NetCDF is available or use -DBUILD_NETCDF=ON")
endif()

# Get NetCDF paths for YAC configuration
if(BUILD_NETCDF)
    # NetCDF was built by BuildNetCDF.cmake
    set(YAC_NETCDF_C_ROOT ${NETCDF_INSTALL_PREFIX})
    set(YAC_NETCDF_Fortran_ROOT ${NETCDF_INSTALL_PREFIX})
else()
    # NetCDF found by FindNETCDF.cmake
    # Extract root directories from include paths
    get_filename_component(YAC_NETCDF_C_ROOT ${NETCDF_C_INCLUDE_DIRECTORIES} DIRECTORY)
    get_filename_component(YAC_NETCDF_Fortran_ROOT ${NETCDF_Fortran_INCLUDE_DIRECTORIES} DIRECTORY)
endif()

# YAC configuration using autotools/configure
# YAC uses autotools build system and provides autogen.sh to run autoreconf
# Note: YAC only builds static libraries (.a), does not support --enable-shared
set(YAC_CONFIGURE_COMMAND
    bash -c "cd <SOURCE_DIR> && ./autogen.sh && CC=${MPI_C_COMPILER} FC=${MPI_Fortran_COMPILER} CFLAGS='-fPIC -I${YAXT_INSTALL_PREFIX}/include -I${LIBFYAML_INSTALL_PREFIX}/include -I${YAC_NETCDF_C_ROOT}/include' FCFLAGS='-fPIC -I${YAXT_INSTALL_PREFIX}/include -I${YAC_NETCDF_Fortran_ROOT}/include' LDFLAGS='-L${YAXT_INSTALL_PREFIX}/lib -L${LIBFYAML_INSTALL_PREFIX}/lib -L${YAC_NETCDF_C_ROOT}/lib -L${YAC_NETCDF_Fortran_ROOT}/lib' LIBS='-lyaxt -lyaxt_c -lfyaml -lnetcdf -lnetcdff' ./configure --prefix=${YAC_INSTALL_PREFIX} --with-yaxt-root=${YAXT_INSTALL_PREFIX} --with-fyaml-root=${LIBFYAML_INSTALL_PREFIX}"
)

# Determine number of parallel jobs
if(NOT DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
    include(ProcessorCount)
    ProcessorCount(N)
    if(N EQUAL 0)
        set(N 1)
    endif()
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
endif()

# Build YAC, for now master
ExternalProject_Add(yac-external
    DEPENDS yaxt-external libfyaml-external
    GIT_REPOSITORY https://gitlab.dkrz.de/dkrz-sw/yac.git
    GIT_TAG master
    GIT_SHALLOW ON
    PREFIX ${CMAKE_BINARY_DIR}/external/yac-build
    CONFIGURE_COMMAND ${YAC_CONFIGURE_COMMAND}
    BUILD_COMMAND make -j${CMAKE_BUILD_PARALLEL_LEVEL}
    BUILD_IN_SOURCE ON
    BUILD_BYPRODUCTS
        ${YAC_INSTALL_PREFIX}/lib/libyac.a
    INSTALL_COMMAND make install
    LOG_DOWNLOAD OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    USES_TERMINAL_DOWNLOAD ON
    USES_TERMINAL_CONFIGURE ON
    USES_TERMINAL_BUILD ON
    USES_TERMINAL_INSTALL ON
    STEP_TARGETS configure build install
)

# Create imported targets for the built YAC libraries
# YAC builds multiple static libraries that need to be linked together
add_library(YAC::yac_core STATIC IMPORTED GLOBAL)
set_target_properties(YAC::yac_core PROPERTIES
    IMPORTED_LOCATION ${YAC_INSTALL_PREFIX}/lib/libyac_core.a
)
add_dependencies(YAC::yac_core yac-external)

add_library(YAC::yac_utils STATIC IMPORTED GLOBAL)
set_target_properties(YAC::yac_utils PROPERTIES
    IMPORTED_LOCATION ${YAC_INSTALL_PREFIX}/lib/libyac_utils.a
)
add_dependencies(YAC::yac_utils yac-external)

add_library(YAC::yac_mtime STATIC IMPORTED GLOBAL)
set_target_properties(YAC::yac_mtime PROPERTIES
    IMPORTED_LOCATION ${YAC_INSTALL_PREFIX}/lib/libyac_mtime.a
)
add_dependencies(YAC::yac_mtime yac-external)

add_library(YAC::yac_mci STATIC IMPORTED GLOBAL)
set_target_properties(YAC::yac_mci PROPERTIES
    IMPORTED_LOCATION ${YAC_INSTALL_PREFIX}/lib/libyac_mci.a
)
add_dependencies(YAC::yac_mci yac-external)

# Find LAPACK (needed by YAC for linear algebra operations)
find_package(LAPACK QUIET)
if(NOT LAPACK_FOUND)
    # Fallback to system LAPACK/BLAS
    set(LAPACK_LIBRARIES "-llapack;-lblas" CACHE STRING "LAPACK libraries" FORCE)
endif()

# Main YAC library that depends on all components
add_library(YAC::yac STATIC IMPORTED GLOBAL)
set_target_properties(YAC::yac PROPERTIES
    IMPORTED_LOCATION ${YAC_INSTALL_PREFIX}/lib/libyac.a
    INTERFACE_INCLUDE_DIRECTORIES "${YAC_INSTALL_PREFIX}/include"
    INTERFACE_LINK_LIBRARIES "YAC::yac_core;YAC::yac_utils;YAC::yac_mtime;YAC::yac_mci;YAXT::yaxt;libfyaml::fyaml;${NETCDF_C_LIBRARIES};${NETCDF_Fortran_LIBRARIES};${LAPACK_LIBRARIES}"
)
add_dependencies(YAC::yac yac-external)

# Set cache variables that FindYAC.cmake and src/CMakeLists.txt expect
set(YAC_FOUND TRUE CACHE BOOL "YAC found" FORCE)
set(yac_DIR ${YAC_INSTALL_PREFIX} CACHE PATH "YAC installation directory" FORCE)
set(YAC_Fortran_INCLUDE_DIRECTORIES ${YAC_INSTALL_PREFIX}/include CACHE PATH "YAC Fortran include directories" FORCE)
set(YAC_LIBRARY ${YAC_INSTALL_PREFIX}/lib/libyac.a CACHE FILEPATH "YAC library" FORCE)
set(YAC_Fortran_LIBRARIES YAC::yac CACHE STRING "YAC Fortran library targets" FORCE)

# Create progress monitoring targets
add_custom_target(yac-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting YAC build (3/3) ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Downloading and building YAC from gitlab.dkrz.de..."
    COMMAND ${CMAKE_COMMAND} -E echo "    YAC depends on YAXT and libfyaml (already built)"
    VERBATIM
)

add_dependencies(yac-external-configure yac-progress)

# Create a convenience target
add_custom_target(build_yac DEPENDS yac-external)

# Add completion message
add_custom_target(yac-complete
    COMMAND ${CMAKE_COMMAND} -E echo "=== YAC build complete! ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    YAC and all dependencies have been built successfully"
    COMMAND ${CMAKE_COMMAND} -E echo "    Installation locations:"
    COMMAND ${CMAKE_COMMAND} -E echo "      libfyaml: ${LIBFYAML_INSTALL_PREFIX}"
    COMMAND ${CMAKE_COMMAND} -E echo "      YAXT:     ${YAXT_INSTALL_PREFIX}"
    COMMAND ${CMAKE_COMMAND} -E echo "      YAC:      ${YAC_INSTALL_PREFIX}"
    VERBATIM
)
add_dependencies(yac-complete yac-external)
add_dependencies(build_yac yac-complete)

message(STATUS "BUILD_YAC is enabled - YAC will be built from source")
message(STATUS "YAC installation prefix: ${YAC_INSTALL_PREFIX}")
message(STATUS "YAC will use:")
message(STATUS "  YAXT from:     ${YAXT_INSTALL_PREFIX}")
message(STATUS "  libfyaml from: ${LIBFYAML_INSTALL_PREFIX}")
message(STATUS "  NetCDF from:   ${YAC_NETCDF_C_ROOT}")
message(STATUS "NOTE: YAC build will start after YAXT and libfyaml complete...")

# Add YAXT and libfyaml library directories to RPATH
# This is necessary because they are shared libraries built from source
# Without this, executables will fail at runtime with "cannot open shared object file"
list(APPEND CMAKE_INSTALL_RPATH "${YAXT_INSTALL_PREFIX}/lib")
list(APPEND CMAKE_INSTALL_RPATH "${LIBFYAML_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}" PARENT_SCOPE)
