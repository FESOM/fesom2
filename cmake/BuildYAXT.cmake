# BuildYAXT.cmake
# External project build for YAXT (Yet Another eXchange Tool)
# YAXT is a dependency of YAC (Yet Another Coupler)

include(ExternalProject)
include(GNUInstallDirs)

# Set installation prefix for built YAXT libraries
set(YAXT_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/yaxt)

# Create directories early so CMake doesn't complain about non-existent paths
file(MAKE_DIRECTORY ${YAXT_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${YAXT_INSTALL_PREFIX}/lib)
file(MAKE_DIRECTORY ${YAXT_INSTALL_PREFIX}/include)

# Find required dependencies
find_package(MPI REQUIRED)

# YAXT uses autotools build system
# Need to run autoreconf before configure since git clone doesn't include configure script
set(YAXT_CONFIGURE_COMMAND
    bash -c "cd <SOURCE_DIR> && autoreconf -iv && CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} CFLAGS=-fPIC FCFLAGS=-fPIC ./configure --prefix=${YAXT_INSTALL_PREFIX} --enable-shared --disable-static --without-regard-for-quality"
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

# Build YAXT
ExternalProject_Add(yaxt-external
    GIT_REPOSITORY https://gitlab.dkrz.de/dkrz-sw/yaxt.git
    GIT_TAG master
    GIT_SHALLOW ON
    PREFIX ${CMAKE_BINARY_DIR}/external/yaxt-build
    CONFIGURE_COMMAND ${YAXT_CONFIGURE_COMMAND}
    BUILD_COMMAND make -j${CMAKE_BUILD_PARALLEL_LEVEL}
    BUILD_IN_SOURCE ON
    BUILD_BYPRODUCTS
        ${YAXT_INSTALL_PREFIX}/lib/libyaxt.so
        ${YAXT_INSTALL_PREFIX}/lib/libyaxt_c.so
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

# Create imported targets for the built YAXT libraries
# This ensures proper dependency propagation and RPATH handling

# YAXT C library imported target
add_library(YAXT::yaxt_c SHARED IMPORTED GLOBAL)
set_target_properties(YAXT::yaxt_c PROPERTIES
    IMPORTED_LOCATION ${YAXT_INSTALL_PREFIX}/lib/libyaxt_c.so
    INTERFACE_INCLUDE_DIRECTORIES ${YAXT_INSTALL_PREFIX}/include
)
add_dependencies(YAXT::yaxt_c yaxt-external)

# YAXT Fortran library imported target
add_library(YAXT::yaxt SHARED IMPORTED GLOBAL)
set_target_properties(YAXT::yaxt PROPERTIES
    IMPORTED_LOCATION ${YAXT_INSTALL_PREFIX}/lib/libyaxt.so
    INTERFACE_INCLUDE_DIRECTORIES ${YAXT_INSTALL_PREFIX}/include
    INTERFACE_LINK_LIBRARIES "YAXT::yaxt_c"
    IMPORTED_LINK_INTERFACE_LIBRARIES "YAXT::yaxt_c"
)
add_dependencies(YAXT::yaxt yaxt-external)

# Set cache variables for use by BuildYAC.cmake and src/CMakeLists.txt
set(YAXT_FOUND TRUE CACHE BOOL "YAXT found" FORCE)
set(YAXT_ROOT ${YAXT_INSTALL_PREFIX} CACHE PATH "YAXT installation directory" FORCE)
set(YAXT_INCLUDE_DIRS ${YAXT_INSTALL_PREFIX}/include CACHE PATH "YAXT include directories" FORCE)
set(YAXT_LIBRARIES "YAXT::yaxt;YAXT::yaxt_c" CACHE STRING "YAXT library targets" FORCE)
# Set variables that src/CMakeLists.txt expects
set(YAXT_Fortran_INCLUDE_DIRECTORIES ${YAXT_INSTALL_PREFIX}/include CACHE PATH "YAXT Fortran include directories" FORCE)
set(YAXT_Fortran_LIBRARIES YAXT::yaxt CACHE STRING "YAXT Fortran library target" FORCE)
set(YAXTC_Fortran_LIBRARIES YAXT::yaxt_c CACHE STRING "YAXT C library target" FORCE)

# Create progress monitoring targets
add_custom_target(yaxt-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting YAXT build ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Downloading and building YAXT from gitlab.dkrz.de..."
    VERBATIM
)

add_dependencies(yaxt-external-configure yaxt-progress)

# Create a convenience target
add_custom_target(build_yaxt DEPENDS yaxt-external)

# Add completion message
add_custom_target(yaxt-complete
    COMMAND ${CMAKE_COMMAND} -E echo "=== YAXT build complete! ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    YAXT has been built successfully at ${YAXT_INSTALL_PREFIX}"
    VERBATIM
)
add_dependencies(yaxt-complete yaxt-external)
add_dependencies(build_yaxt yaxt-complete)

message(STATUS "BUILD_YAC is enabled - YAXT will be built from source as a dependency of YAC")
message(STATUS "YAXT installation prefix: ${YAXT_INSTALL_PREFIX}")
message(STATUS "NOTE: YAXT build progress will be shown in real-time below...")
