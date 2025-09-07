#!/usr/bin/env bash

set -e

HERE=$PWD
SOURCE_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
BUILD_DIR=${BUILD_DIR:-build}

# Parse arguments
CLEAN_BUILD=false
ENV_ARGS=()

# Process arguments
for arg in "$@"; do
    if [ "$arg" = "--clean" ]; then
        CLEAN_BUILD=true
    else
        ENV_ARGS+=("$arg")
    fi
done

# Clean if requested - do this before sourcing env.sh as cleaning doesn't need environment setup
if [ "$CLEAN_BUILD" = true ]; then
    echo "Cleaning build directory and installed files..."
    if [ -d "${BUILD_DIR}" ]; then
        if [ -f "${BUILD_DIR}/CMakeCache.txt" ]; then
            # Only attempt make clean and uninstall if we have a valid build directory
            cd ${BUILD_DIR}
            echo "Running make clean..."
            make clean 2>/dev/null || true
            
            echo "Running make uninstall..."
            make uninstall 2>/dev/null || echo "No uninstall target available"
            
            cd ${HERE}
        fi
        
        # Remove the build directory for a completely fresh build
        echo "Removing build directory..."
        rm -rf ${BUILD_DIR}
    fi
    echo "Clean completed successfully"
    # Recreate build directory
    mkdir -p ${BUILD_DIR}
fi

# Source environment with the original arguments (minus --clean)
source env.sh "${ENV_ARGS[@]}"

mkdir -p ${BUILD_DIR} # make sure not to commit this to svn or git
cd ${BUILD_DIR}

# Don't pass environment arguments (like 'ubuntu') to cmake
# Only pass arguments that start with - (like -DENABLE_OPENMP=ON)
CMAKE_ONLY_ARGS=""
for arg in "${ENV_ARGS[@]}"; do
    if [[ "$arg" =~ ^- ]]; then
        CMAKE_ONLY_ARGS="$CMAKE_ONLY_ARGS $arg"
    fi
done

cmake ${SOURCE_DIR} -DCMAKE_INSTALL_PREFIX=$HERE -DCMAKE_BUILD_TYPE=Debug ${CMAKE_ARGS} ${CMAKE_ONLY_ARGS}
    # not required when re-compiling
    # additional cmake arguments can be passed to configure.sh
    # this also includes fesom specific options in CMakeLists, can be used as -DFESOM_COUPLED=ON
make install -j`nproc --all`
