#!/usr/bin/env bash

set -e

# Get the absolute path of the script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FESOM_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${SCRIPT_DIR}/build"
BIN_DIR="${FESOM_ROOT}/bin"
METIS_DIR="${FESOM_ROOT}/lib/metis-5.1.0"

# Create build directory
mkdir -p "${BUILD_DIR}"

# Create bin directory if it doesn't exist
mkdir -p "${BIN_DIR}"

# Source environment if env.sh exists
if [ -f "${FESOM_ROOT}/env.sh" ]; then
    echo "Sourcing environment from ${FESOM_ROOT}/env.sh"
    source "${FESOM_ROOT}/env.sh"
fi

# Check if METIS exists
if [ ! -d "${METIS_DIR}" ]; then
    echo "ERROR: METIS not found at ${METIS_DIR}"
    echo "Please ensure the METIS library is available in the FESOM2 lib directory."
    exit 1
fi

# Configure with CMake
echo "Configuring mesh partitioner build in ${BUILD_DIR}"
cd "${BUILD_DIR}"
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${FESOM_ROOT}" \
    -DFESOM_SOURCE_DIR="${FESOM_ROOT}" \
    -DMETIS_ROOT="${METIS_DIR}" \
    "$@"

# Check if configuration was successful
if [ $? -ne 0 ]; then
    echo "CMake configuration failed. Please check the error messages above."
    exit 1
fi

# Build the mesh partitioner
echo "Building mesh partitioner..."
make -j$(nproc)

if [ $? -ne 0 ]; then
    echo "Build failed. Please check the error messages above."
    exit 1
fi

# Create bin directory if it doesn't exist
mkdir -p "${BIN_DIR}"

# Create a symlink in the bin directory
ln -sf "${BUILD_DIR}/fesom_meshpart" "${BIN_DIR}/fesom_meshpart"

# Install the executable
make install

echo ""
echo "Build complete!"
echo "The fesom_meshpart executable is available at:"
echo "  ${BUILD_DIR}/fesom_meshpart"
echo "A symlink has been created at:"
echo "  ${BIN_DIR}/fesom_meshpart"


