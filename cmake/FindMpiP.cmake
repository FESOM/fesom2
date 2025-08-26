# FindMpiP.cmake - Find or build mpiP MPI profiling library
#
# Variables set by this module:
#   MpiP_FOUND        - True if mpiP was found/built successfully
#   MpiP_LIBRARIES    - Libraries to link against for mpiP
#   MpiP_LIBRARY_DIR  - Directory containing mpiP library

cmake_minimum_required(VERSION 3.16)

include(ExternalProject)
include(FindPackageHandleStandardArgs)

# Set default paths
set(MpiP_INSTALL_DIR "${CMAKE_BINARY_DIR}/external/mpip-install" CACHE PATH "mpiP install directory")

# First, try to find an existing mpiP installation
find_library(MpiP_LIBRARY
    NAMES mpiP
    HINTS ${MpiP_INSTALL_DIR}/lib
    PATHS /usr/local/lib /opt/local/lib
)

# If not found, download and build mpiP
if(NOT MpiP_LIBRARY)
    message(STATUS "mpiP not found, will download and build from GitHub")
    
    # Find MPI compiler (required by mpiP)
    find_package(MPI REQUIRED)
    if(NOT MPI_C_COMPILER)
        message(FATAL_ERROR "MPI C compiler not found. mpiP requires MPI.")
    endif()

    # Set up ExternalProject to download and build mpiP with debug symbols
    ExternalProject_Add(mpip_external
        PREFIX ${CMAKE_BINARY_DIR}/external/mpip
        
        # Download from GitHub
        GIT_REPOSITORY "https://github.com/LLNL/mpiP.git"
        GIT_TAG "3.5"
        GIT_SHALLOW TRUE
        
        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/mpip-src
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/mpip-src
        INSTALL_DIR ${MpiP_INSTALL_DIR}
        
        # Configure step with debug symbols and BFD/unwind support
        CONFIGURE_COMMAND 
            ${CMAKE_COMMAND} -E env 
                CC=${MPI_C_COMPILER} 
                CXX=${MPI_CXX_COMPILER} 
                FC=${MPI_Fortran_COMPILER}
                CFLAGS=-g\ -O2 
                CXXFLAGS=-g\ -O2
            <SOURCE_DIR>/configure 
                --prefix=<INSTALL_DIR>
                --enable-bfd
                --enable-unwind
        
        # Build step
        BUILD_COMMAND make -j${CMAKE_BUILD_PARALLEL_LEVEL}
        
        # Install step
        INSTALL_COMMAND make install
        
        # Update timestamps to avoid unnecessary rebuilds
        UPDATE_COMMAND ""
        
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )

    # Set the library path after build
    set(MpiP_LIBRARY "${MpiP_INSTALL_DIR}/lib/libmpiP.so" CACHE FILEPATH "mpiP library" FORCE)
    
    # Create an imported target that depends on the ExternalProject
    add_library(MpiP::MpiP SHARED IMPORTED GLOBAL)
    add_dependencies(MpiP::MpiP mpip_external)
    
    set_target_properties(MpiP::MpiP PROPERTIES
        IMPORTED_LOCATION "${MpiP_LIBRARY}"
    )
    
    # Indicate that we built mpiP from source and targets need to depend on mpip_external
    set(MpiP_BUILT_FROM_SOURCE TRUE CACHE BOOL "mpiP was built from source")
    set(MpiP_FOUND TRUE)
    
else()
    message(STATUS "Found existing mpiP installation")
    
    # Create imported target for existing installation
    add_library(MpiP::MpiP SHARED IMPORTED GLOBAL)
    set_target_properties(MpiP::MpiP PROPERTIES
        IMPORTED_LOCATION "${MpiP_LIBRARY}"
    )
    
    # Using existing installation, no external project dependency needed
    set(MpiP_BUILT_FROM_SOURCE FALSE CACHE BOOL "mpiP was built from source")
    set(MpiP_FOUND TRUE)
endif()

# Get library directory from library path
if(MpiP_LIBRARY)
    get_filename_component(MpiP_LIBRARY_DIR "${MpiP_LIBRARY}" DIRECTORY)
endif()

# Set standard variables
set(MpiP_LIBRARIES "${MpiP_LIBRARY}")

# Standard CMake package handling
find_package_handle_standard_args(MpiP
    REQUIRED_VARS MpiP_LIBRARY
)