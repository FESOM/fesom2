# BuildCVMix.cmake
# External project build for CVMix library
# Follows the pattern of BuildNetCDF.cmake for consistency

include(ExternalProject)
include(GNUInstallDirs)

# Set installation prefix for built CVMix library
set(CVMIX_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/cvmix)

# Create directories early so CMake doesn't complain about non-existent paths
file(MAKE_DIRECTORY ${CVMIX_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${CVMIX_INSTALL_PREFIX}/lib)
file(MAKE_DIRECTORY ${CVMIX_INSTALL_PREFIX}/include)
file(MAKE_DIRECTORY ${CVMIX_INSTALL_PREFIX}/include/cvmix)

# CVMix git repository and version
set(CVMIX_GIT_REPOSITORY "https://github.com/CVMix/CVMix-src.git")
set(CVMIX_GIT_TAG "v1.0.0")  # Can be made configurable via cache variable

# CRITICAL: Set common cmake args for CVMix build
# CVMix must be compiled with compatible Fortran flags to link properly with FESOM2
set(CVMIX_CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CVMIX_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_MODULE_DIRECTORY=${CVMIX_INSTALL_PREFIX}/include
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
)

# Pass CMAKE_Fortran_FLAGS to CVMix if set
# This ensures CVMix is compiled with the same flags as FESOM2 (e.g., -fdefault-real-8, -r8)
if(CMAKE_Fortran_FLAGS)
    list(APPEND CVMIX_CMAKE_ARGS
        "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
    )
    message(STATUS "CVMix will be built with CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif()

# Pass build-type-specific flags if they exist
if(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_BUILD_TYPE MATCHES "Release")
    list(APPEND CVMIX_CMAKE_ARGS
        "-DCMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}"
    )
endif()

if(CMAKE_Fortran_FLAGS_DEBUG AND CMAKE_BUILD_TYPE MATCHES "Debug")
    list(APPEND CVMIX_CMAKE_ARGS
        "-DCMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}"
    )
endif()

# Intel compiler: Ensure critical real/integer kind flags are passed
# FESOM2 uses -r8 (8-byte reals) and -i4 (4-byte integers) for Intel
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # Check if -r8 or similar flags are in CMAKE_Fortran_FLAGS
    # If not explicitly set, add them to ensure compatibility
    string(FIND "${CMAKE_Fortran_FLAGS}" "-r8" has_r8)
    string(FIND "${CMAKE_Fortran_FLAGS}" "-real-size" has_real_size)

    if(has_r8 EQUAL -1 AND has_real_size EQUAL -1)
        # Add Intel-specific flags for real/integer kinds if not already present
        set(CVMIX_EXTRA_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i4")
        list(APPEND CVMIX_CMAKE_ARGS
            "-DCMAKE_Fortran_FLAGS=${CVMIX_EXTRA_FLAGS}"
        )
        message(STATUS "CVMix (Intel): Adding -r8 -i4 flags for real/integer kind compatibility")
    endif()
endif()

# GNU compiler: Ensure critical real/double kind flags are passed
# FESOM2 uses -fdefault-real-8 -fdefault-double-8 for GNU
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    string(FIND "${CMAKE_Fortran_FLAGS}" "-fdefault-real-8" has_real8)

    if(has_real8 EQUAL -1)
        # Add GNU-specific flags for real/double kinds if not already present
        set(CVMIX_EXTRA_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
        list(APPEND CVMIX_CMAKE_ARGS
            "-DCMAKE_Fortran_FLAGS=${CVMIX_EXTRA_FLAGS}"
        )
        message(STATUS "CVMix (GNU): Adding -fdefault-real-8 -fdefault-double-8 flags for kind compatibility")
    endif()
endif()

# Cray compiler: Ensure critical real kind flags are passed
# FESOM2 uses -s real64 for Cray
if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    string(FIND "${CMAKE_Fortran_FLAGS}" "-s real64" has_real64)

    if(has_real64 EQUAL -1)
        set(CVMIX_EXTRA_FLAGS "${CMAKE_Fortran_FLAGS} -s real64")
        list(APPEND CVMIX_CMAKE_ARGS
            "-DCMAKE_Fortran_FLAGS=${CVMIX_EXTRA_FLAGS}"
        )
        message(STATUS "CVMix (Cray): Adding -s real64 flag for real kind compatibility")
    endif()
endif()

# NVIDIA/PGI compiler: Ensure critical real kind flags are passed
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    string(FIND "${CMAKE_Fortran_FLAGS}" "-r8" has_r8)

    if(has_r8 EQUAL -1)
        set(CVMIX_EXTRA_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
        list(APPEND CVMIX_CMAKE_ARGS
            "-DCMAKE_Fortran_FLAGS=${CVMIX_EXTRA_FLAGS}"
        )
        message(STATUS "CVMix (NVHPC): Adding -r8 flag for real kind compatibility")
    endif()
endif()

# Build CVMix library using ExternalProject
ExternalProject_Add(cvmix-build
    GIT_REPOSITORY ${CVMIX_GIT_REPOSITORY}
    GIT_TAG ${CVMIX_GIT_TAG}
    GIT_SHALLOW TRUE
    PREFIX ${CMAKE_BINARY_DIR}/external/cvmix-src
    CMAKE_ARGS ${CVMIX_CMAKE_ARGS}
    CMAKE_CACHE_ARGS
        -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    BUILD_ALWAYS OFF
    UPDATE_DISCONNECTED TRUE
    INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
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

# Create proper imported target for the built CVMix library
# This ensures proper dependency propagation
add_library(CVMix::cvmix STATIC IMPORTED GLOBAL)
set_target_properties(CVMix::cvmix PROPERTIES
    IMPORTED_LOCATION ${CVMIX_INSTALL_PREFIX}/lib/libcvmix.a
    INTERFACE_INCLUDE_DIRECTORIES ${CVMIX_INSTALL_PREFIX}/include/cvmix
)
add_dependencies(CVMix::cvmix cvmix-build)

# Set the CVMix variables that will be used by src/CMakeLists.txt
set(CVMIX_FOUND TRUE CACHE BOOL "CVMix library found" FORCE)
set(CVMIX_ROOT ${CVMIX_INSTALL_PREFIX} CACHE PATH "CVMix installation root" FORCE)
set(CVMIX_LIBRARY CVMix::cvmix CACHE STRING "CVMix library target" FORCE)
set(CVMIX_INCLUDE_DIR ${CVMIX_INSTALL_PREFIX}/include/cvmix CACHE PATH "CVMix include directory" FORCE)

# Create progress monitoring targets
add_custom_target(cvmix-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting CVMix build ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Repository: ${CVMIX_GIT_REPOSITORY}"
    COMMAND ${CMAKE_COMMAND} -E echo "    Version: ${CVMIX_GIT_TAG}"
    COMMAND ${CMAKE_COMMAND} -E echo "    Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}"
    COMMAND ${CMAKE_COMMAND} -E echo "    Build type: ${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} -E echo "    Installing to: ${CVMIX_INSTALL_PREFIX}"
    VERBATIM
)

add_dependencies(cvmix-build-configure cvmix-progress)

# Create a convenience target
add_custom_target(build_cvmix DEPENDS cvmix-build)

# Add completion message
add_custom_target(cvmix-complete
    COMMAND ${CMAKE_COMMAND} -E echo "=== CVMix build complete! ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Library: ${CVMIX_INSTALL_PREFIX}/lib/libcvmix.a"
    COMMAND ${CMAKE_COMMAND} -E echo "    Modules: ${CVMIX_INSTALL_PREFIX}/include/"
    VERBATIM
)
add_dependencies(cvmix-complete cvmix-build)
add_dependencies(build_cvmix cvmix-complete)

message(STATUS "CVMix will be built from source using ExternalProject")
message(STATUS "  Version: ${CVMIX_GIT_TAG}")
message(STATUS "  Installation prefix: ${CVMIX_INSTALL_PREFIX}")
message(STATUS "  Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
if(CMAKE_Fortran_FLAGS)
    message(STATUS "  Fortran flags: ${CMAKE_Fortran_FLAGS}")
endif()
