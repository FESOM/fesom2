# BuildLibfyaml.cmake
# External project build for libfyaml (YAML parser library)
# libfyaml is a dependency of YAC (Yet Another Coupler)

include(ExternalProject)
include(GNUInstallDirs)

# Set installation prefix for built libfyaml
set(LIBFYAML_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/libfyaml)

# Create directories early so CMake doesn't complain about non-existent paths
file(MAKE_DIRECTORY ${LIBFYAML_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${LIBFYAML_INSTALL_PREFIX}/lib)
file(MAKE_DIRECTORY ${LIBFYAML_INSTALL_PREFIX}/include)

# libfyaml uses autotools/configure build system
# We need to configure it with appropriate flags
set(LIBFYAML_CONFIGURE_COMMAND
    <SOURCE_DIR>/bootstrap.sh &&
    <SOURCE_DIR>/configure
    --prefix=${LIBFYAML_INSTALL_PREFIX}
    --enable-shared
    --disable-static
    CC=${CMAKE_C_COMPILER}
    CFLAGS=-fPIC
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

# Build libfyaml
# Note: We only build the library, not the tools/utilities to avoid linking issues
ExternalProject_Add(libfyaml-external
    GIT_REPOSITORY https://github.com/pantoniou/libfyaml.git
    GIT_TAG master
    GIT_SHALLOW ON
    PREFIX ${CMAKE_BINARY_DIR}/external/libfyaml-build
    CONFIGURE_COMMAND ${LIBFYAML_CONFIGURE_COMMAND}
    BUILD_COMMAND make -j${CMAKE_BUILD_PARALLEL_LEVEL} -C src libfyaml.la
    BUILD_IN_SOURCE ON
    BUILD_BYPRODUCTS ${LIBFYAML_INSTALL_PREFIX}/lib/libfyaml.so
    INSTALL_COMMAND make -C src install-libLTLIBRARIES install-includeHEADERS && make install-pkgconfigDATA
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

# Create imported target for the built libfyaml library
add_library(libfyaml::fyaml SHARED IMPORTED GLOBAL)
set_target_properties(libfyaml::fyaml PROPERTIES
    IMPORTED_LOCATION ${LIBFYAML_INSTALL_PREFIX}/lib/libfyaml.so
    INTERFACE_INCLUDE_DIRECTORIES ${LIBFYAML_INSTALL_PREFIX}/include
)
add_dependencies(libfyaml::fyaml libfyaml-external)

# Set cache variables for use by BuildYAC.cmake
set(LIBFYAML_FOUND TRUE CACHE BOOL "libfyaml found" FORCE)
set(LIBFYAML_ROOT ${LIBFYAML_INSTALL_PREFIX} CACHE PATH "libfyaml installation directory" FORCE)
set(LIBFYAML_INCLUDE_DIRS ${LIBFYAML_INSTALL_PREFIX}/include CACHE PATH "libfyaml include directories" FORCE)
set(LIBFYAML_LIBRARIES libfyaml::fyaml CACHE STRING "libfyaml library target" FORCE)

# Create progress monitoring targets
add_custom_target(libfyaml-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting libfyaml build ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Downloading and building libfyaml from github.com/pantoniou/libfyaml..."
    VERBATIM
)

add_dependencies(libfyaml-external-configure libfyaml-progress)

# Create a convenience target
add_custom_target(build_libfyaml DEPENDS libfyaml-external)

# Add completion message
add_custom_target(libfyaml-complete
    COMMAND ${CMAKE_COMMAND} -E echo "=== libfyaml build complete! ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    libfyaml has been built successfully at ${LIBFYAML_INSTALL_PREFIX}"
    VERBATIM
)
add_dependencies(libfyaml-complete libfyaml-external)
add_dependencies(build_libfyaml libfyaml-complete)

message(STATUS "BUILD_YAC is enabled - libfyaml will be built from source as a dependency of YAC")
message(STATUS "libfyaml installation prefix: ${LIBFYAML_INSTALL_PREFIX}")
message(STATUS "NOTE: libfyaml build progress will be shown in real-time below...")
