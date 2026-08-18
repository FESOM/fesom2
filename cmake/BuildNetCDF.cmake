# BuildNetCDF.cmake
# External project build for NetCDF-C and NetCDF-Fortran libraries
# when system libraries are incompatible with the selected compiler

include(ExternalProject)
include(GNUInstallDirs)

# Set installation prefix for built NetCDF libraries
set(NETCDF_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/external/netcdf)

# Create directories early so CMake doesn't complain about non-existent paths
file(MAKE_DIRECTORY ${NETCDF_INSTALL_PREFIX})
file(MAKE_DIRECTORY ${NETCDF_INSTALL_PREFIX}/lib)
file(MAKE_DIRECTORY ${NETCDF_INSTALL_PREFIX}/include)

# Find required dependencies
find_package(MPI REQUIRED)

# Find ZLIB with more robust search
find_package(ZLIB QUIET)
if(NOT ZLIB_FOUND)
    # Try pkg-config approach
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(ZLIB QUIET zlib)
    endif()
    
    # Manual search if pkg-config fails
    if(NOT ZLIB_FOUND)
        find_path(ZLIB_INCLUDE_DIRS zlib.h PATHS /usr/include /usr/local/include)
        find_library(ZLIB_LIBRARIES z PATHS /usr/lib /usr/lib64 /usr/local/lib /lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu)
        if(ZLIB_INCLUDE_DIRS AND ZLIB_LIBRARIES)
            set(ZLIB_FOUND TRUE)
        endif()
    endif()
endif()

if(NOT ZLIB_FOUND)
    message(FATAL_ERROR "ZLIB is required for NetCDF but could not be found. Please install zlib development packages (e.g., zlib1g-dev)")
endif()

# Find HDF5 with robust detection - prefer serial version to avoid MPI conflicts with Intel MPI
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    # Use serial HDF5 to avoid MPI library conflicts between Intel MPI and OpenMPI HDF5
    pkg_check_modules(HDF5 QUIET hdf5-serial)
    if(HDF5_FOUND)
        set(HDF5_C_LIBRARIES ${HDF5_LINK_LIBRARIES})
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
        set(HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS})
        message(STATUS "Found serial HDF5 via pkg-config (avoids MPI conflicts):")
        message(STATUS "  Libraries: ${HDF5_C_LIBRARIES}")
        message(STATUS "  Include dirs: ${HDF5_INCLUDE_DIRS}")
        message(STATUS "  Library dirs: ${HDF5_LIBRARY_DIRS}")
    else()
        message(WARNING "Serial HDF5 not found, falling back to system HDF5 (may cause MPI conflicts)")
        pkg_check_modules(HDF5 QUIET hdf5)
        if(HDF5_FOUND)
            set(HDF5_C_LIBRARIES ${HDF5_LINK_LIBRARIES})
            set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
            set(HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS})
            message(STATUS "Found HDF5 via pkg-config:")
            message(STATUS "  Libraries: ${HDF5_C_LIBRARIES}")
        endif()
    endif()
endif()

# Fallback to CMake's find_package if pkg-config failed
if(NOT HDF5_FOUND)
    find_package(HDF5 QUIET COMPONENTS C)
    if(HDF5_FOUND)
        message(STATUS "Found HDF5 via CMake: ${HDF5_C_LIBRARIES}")
    endif()
endif()

# Manual fallback if both approaches failed
if(NOT HDF5_FOUND)
    find_path(HDF5_INCLUDE_DIRS hdf5.h PATHS /usr/include/hdf5/openmpi /usr/include/hdf5/serial /usr/include)
    find_library(HDF5_C_LIBRARIES hdf5 PATHS /usr/lib/x86_64-linux-gnu/hdf5/openmpi /usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib)
    if(HDF5_INCLUDE_DIRS AND HDF5_C_LIBRARIES)
        set(HDF5_FOUND TRUE)
        message(STATUS "Found HDF5 manually: ${HDF5_C_LIBRARIES}")
    endif()
endif()

# Set common cmake args for both NetCDF-C and NetCDF-Fortran
set(NETCDF_COMMON_CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${NETCDF_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DBUILD_SHARED_LIBS=ON
    -DENABLE_TESTS=OFF
    -DENABLE_EXAMPLES=OFF
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
)

# Add compiler-specific flags for Intel compilers  
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # Find Intel OneAPI MPI library directory using latest symlink if available
    # Use environment variable if available, otherwise try common paths
    if(DEFINED ENV{INTEL_ONEAPI_ROOT})
        set(INTEL_ONEAPI_SEARCH_PATHS "$ENV{INTEL_ONEAPI_ROOT}/mpi/latest" "$ENV{INTEL_ONEAPI_ROOT}/mpi/2021.16")
    else()
        set(INTEL_ONEAPI_SEARCH_PATHS 
            /opt/intel/oneapi/mpi/latest 
            /opt/intel/oneapi/mpi/2021.16
            $ENV{HOME}/intel/oneapi/mpi/latest
            $ENV{HOME}/intel/oneapi/mpi/2021.16
            /usr/local/intel/oneapi/mpi/latest
            /usr/local/intel/oneapi/mpi/2021.16
        )
    endif()
    
    find_file(INTEL_MPI_LIB_DIR 
        NAMES lib
        PATHS ${INTEL_ONEAPI_SEARCH_PATHS}
        NO_DEFAULT_PATH
    )
    if(INTEL_MPI_LIB_DIR)
        get_filename_component(INTEL_MPI_LIB_PATH ${INTEL_MPI_LIB_DIR} DIRECTORY)
        set(INTEL_MPI_LIB_PATH "${INTEL_MPI_LIB_PATH}/lib")
        list(APPEND NETCDF_COMMON_CMAKE_ARGS
            -DCMAKE_C_FLAGS=-w
            -DCMAKE_CXX_FLAGS=-w
            -DCMAKE_EXE_LINKER_FLAGS=-L${INTEL_MPI_LIB_PATH}
            -DCMAKE_SHARED_LINKER_FLAGS=-L${INTEL_MPI_LIB_PATH}
            -DENABLE_PARALLEL=OFF
            -DENABLE_PARALLEL_TESTS=OFF
            -DENABLE_PNETCDF=OFF
        )
        message(STATUS "Added Intel MPI library path for NetCDF build: ${INTEL_MPI_LIB_PATH}")
    else()
        message(WARNING "Intel MPI library directory not found")
    endif()
endif()

# Add compiler-specific flags for NVIDIA compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    list(APPEND NETCDF_COMMON_CMAKE_ARGS
        -DCMAKE_C_FLAGS=-w
        -DCMAKE_CXX_FLAGS=-w
        -DENABLE_PARALLEL=OFF
        -DENABLE_PARALLEL_TESTS=OFF
        -DENABLE_PNETCDF=OFF
    )
    message(STATUS "Added NVIDIA compiler flags for NetCDF build")
endif()

# Add ZLIB paths (required)
list(APPEND NETCDF_COMMON_CMAKE_ARGS 
    -DZLIB_LIBRARY=${ZLIB_LIBRARIES}
    -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIRS}
)

# Add math library for NetCDF-C using environment-provided paths
if(DEFINED ENV{MATH_LIBRARY_PATHS})
    string(REPLACE ":" ";" MATH_LIB_PATHS "$ENV{MATH_LIBRARY_PATHS}")
    find_library(MATH_LIBRARY m PATHS ${MATH_LIB_PATHS} NO_DEFAULT_PATH)
    if(NOT MATH_LIBRARY)
        find_library(MATH_LIBRARY m PATHS ${MATH_LIB_PATHS})
    endif()
else()
    find_library(MATH_LIBRARY m)
endif()

if(MATH_LIBRARY)
    message(STATUS "Found math library for NetCDF: ${MATH_LIBRARY}")
    # Explicitly provide math library path to NetCDF build
    list(APPEND NETCDF_COMMON_CMAKE_ARGS 
        -DMATH_LIBRARY=${MATH_LIBRARY}
        -DCMAKE_REQUIRED_LIBRARIES=${MATH_LIBRARY}
    )
else()
    message(WARNING "Math library not found - trying fallback options")
    # Try standard system paths as fallback
    foreach(MATH_PATH /usr/lib/x86_64-linux-gnu /usr/lib /lib/x86_64-linux-gnu /lib /usr/lib64 /lib64)
        if(EXISTS "${MATH_PATH}/libm.so" OR EXISTS "${MATH_PATH}/libm.a")
            set(MATH_LIBRARY "${MATH_PATH}/libm.so")
            if(NOT EXISTS "${MATH_LIBRARY}")
                set(MATH_LIBRARY "${MATH_PATH}/libm.a")
            endif()
            message(STATUS "Found math library fallback: ${MATH_LIBRARY}")
            list(APPEND NETCDF_COMMON_CMAKE_ARGS 
                -DMATH_LIBRARY=${MATH_LIBRARY}
                -DCMAKE_REQUIRED_LIBRARIES=${MATH_LIBRARY}
            )
            break()
        endif()
    endforeach()
endif()

# Add HDF5 support for NetCDF-4 and compression capabilities

# Add HDF5 configuration if available
set(NETCDF_C_CMAKE_ARGS ${NETCDF_COMMON_CMAKE_ARGS})
if(HDF5_FOUND)
    # Find the actual HDF5 library file
    find_library(HDF5_ACTUAL_LIBRARY hdf5 PATHS ${HDF5_LIBRARY_DIRS} NO_DEFAULT_PATH)
    if(NOT HDF5_ACTUAL_LIBRARY)
        find_library(HDF5_ACTUAL_LIBRARY hdf5)
    endif()
    
    list(APPEND NETCDF_C_CMAKE_ARGS
        -DHDF5_C_LIBRARY=${HDF5_ACTUAL_LIBRARY}
        -DHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIRS}
        -DCMAKE_LIBRARY_PATH=${HDF5_LIBRARY_DIRS}
        -DENABLE_HDF5=ON
        -DUSE_HDF5=ON
        -DENABLE_NETCDF_4=ON
    )
    message(STATUS "NetCDF will use HDF5 library: ${HDF5_ACTUAL_LIBRARY}")
else()
    list(APPEND NETCDF_C_CMAKE_ARGS
        -DENABLE_HDF5=OFF
        -DUSE_HDF5=OFF
        -DENABLE_NETCDF_4=OFF
    )
endif()

# Build NetCDF-C first (required dependency for NetCDF-Fortran)
ExternalProject_Add(netcdf-c
    URL "https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz"
    URL_HASH SHA256=bc104d101278c68b303359b3dc4192f81592ae8640f1aee486921138f7f88cb7
    DOWNLOAD_NAME "netcdf-c-4.9.2.tar.gz"
    PREFIX ${CMAKE_BINARY_DIR}/external/netcdf-c
    CMAKE_ARGS
        ${NETCDF_C_CMAKE_ARGS}
        -DENABLE_DAP=OFF
        -DENABLE_BYTERANGE=OFF
        -DENABLE_PARALLEL_TESTS=OFF
        -DENABLE_PNETCDF=OFF
        -DENABLE_PLUGINS=OFF
        -DENABLE_FILTER_TESTING=OFF
        -DENABLE_NCZARR=OFF
        -DENABLE_NCZARR_FILTERS=OFF
        -DENABLE_LIBXML2=OFF
    BUILD_ALWAYS OFF
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

# Build NetCDF-Fortran (depends on NetCDF-C) - using older compatible version
ExternalProject_Add(netcdf-fortran
    DEPENDS netcdf-c
    URL "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.5.4.tar.gz"
    DOWNLOAD_NAME "netcdf-fortran-4.5.4.tar.gz"
    PREFIX ${CMAKE_BINARY_DIR}/external/netcdf-fortran
    CMAKE_ARGS
        ${NETCDF_COMMON_CMAKE_ARGS}
        -DnetCDF_DIR=${NETCDF_INSTALL_PREFIX}/lib/cmake/netCDF
        -DCMAKE_PREFIX_PATH=${NETCDF_INSTALL_PREFIX}
        -DENABLE_TESTS=OFF
        -DENABLE_FORTRAN_TYPE_CHECKS=OFF
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5
    BUILD_ALWAYS OFF
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

# Create proper imported targets for the built NetCDF libraries
# This ensures proper dependency propagation and RPATH handling

# NetCDF-C imported target
add_library(NetCDF::netcdf SHARED IMPORTED GLOBAL)
set(NETCDF_LINK_LIBS "m;z")
if(HDF5_FOUND AND HDF5_ACTUAL_LIBRARY)
    list(APPEND NETCDF_LINK_LIBS ${HDF5_ACTUAL_LIBRARY})
endif()
set_target_properties(NetCDF::netcdf PROPERTIES
    IMPORTED_LOCATION ${NETCDF_INSTALL_PREFIX}/lib/libnetcdf.so
    INTERFACE_INCLUDE_DIRECTORIES ${NETCDF_INSTALL_PREFIX}/include
    IMPORTED_LINK_INTERFACE_LIBRARIES "${NETCDF_LINK_LIBS}"
)
add_dependencies(NetCDF::netcdf netcdf-c)

# NetCDF-Fortran imported target
add_library(NetCDF::netcdff SHARED IMPORTED GLOBAL)  
set_target_properties(NetCDF::netcdff PROPERTIES
    IMPORTED_LOCATION ${NETCDF_INSTALL_PREFIX}/lib/libnetcdff.so
    INTERFACE_INCLUDE_DIRECTORIES ${NETCDF_INSTALL_PREFIX}/include
    INTERFACE_LINK_LIBRARIES "NetCDF::netcdf"
    IMPORTED_LINK_INTERFACE_LIBRARIES "NetCDF::netcdf"
)
add_dependencies(NetCDF::netcdff netcdf-fortran)

# Set the NetCDF variables that will be used by src/CMakeLists.txt
set(NETCDF_C_FOUND 1 CACHE BOOL "NetCDF-C found" FORCE)
set(NETCDF_Fortran_FOUND 1 CACHE BOOL "NetCDF-Fortran found" FORCE)
set(NETCDF_C_INCLUDE_DIRECTORIES ${NETCDF_INSTALL_PREFIX}/include CACHE PATH "NetCDF-C include directories" FORCE)
set(NETCDF_Fortran_INCLUDE_DIRECTORIES ${NETCDF_INSTALL_PREFIX}/include CACHE PATH "NetCDF-Fortran include directories" FORCE)

# Use the imported targets - this ensures proper dependency propagation
set(NETCDF_C_LIBRARIES NetCDF::netcdf CACHE STRING "NetCDF-C library target" FORCE)
set(NETCDF_Fortran_LIBRARIES NetCDF::netcdff CACHE STRING "NetCDF-Fortran library target" FORCE)

# Create progress monitoring targets
add_custom_target(netcdf-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting NetCDF-C build (1/2) ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Downloading, configuring, building, and installing NetCDF-C..."
    VERBATIM
)

add_custom_target(netcdf-fortran-progress
    COMMAND ${CMAKE_COMMAND} -E echo "=== Starting NetCDF-Fortran build (2/2) ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Building NetCDF-Fortran (depends on NetCDF-C)..."
    VERBATIM
)

# Add dependencies to show progress messages
add_dependencies(netcdf-c-configure netcdf-progress)
add_dependencies(netcdf-fortran-configure netcdf-fortran-progress)

# Create a convenience target that ensures NetCDF is built before FESOM
add_custom_target(build_netcdf DEPENDS netcdf-c netcdf-fortran)

# Add completion message
add_custom_target(netcdf-complete
    COMMAND ${CMAKE_COMMAND} -E echo "=== NetCDF build complete! ==="
    COMMAND ${CMAKE_COMMAND} -E echo "    Both NetCDF-C and NetCDF-Fortran have been built successfully."
    VERBATIM
)
add_dependencies(netcdf-complete netcdf-c netcdf-fortran)
add_dependencies(build_netcdf netcdf-complete)

message(STATUS "BUILD_NETCDF is enabled - NetCDF-C v4.9.2 and NetCDF-Fortran v4.5.4 will be built from source")
message(STATUS "NetCDF installation prefix: ${NETCDF_INSTALL_PREFIX}")
message(STATUS "NOTE: NetCDF build progress will be shown in real-time below...")
message(STATUS "      This may take several minutes depending on your system.")