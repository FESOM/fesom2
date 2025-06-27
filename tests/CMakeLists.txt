#===============================================================================
# tests/CMakeLists.txt - Main testing configuration for FESOM2
#===============================================================================

cmake_minimum_required(VERSION 3.16)

# Include our testing utilities
include(${CMAKE_SOURCE_DIR}/cmake/FesomTesting.cmake)

# Set up MPI for testing
setup_mpi_testing()

# Ensure the fesom target exists
if(NOT TARGET fesom)
    message(FATAL_ERROR "fesom target not found. Tests require the main executable to be built.")
endif()

# Testing configuration options
option(ENABLE_INTEGRATION_TESTS "Enable integration tests that run FESOM" ON)
option(ENABLE_MPI_TESTS "Enable MPI tests (requires mpiexec/mpirun)" ON)
set(TEST_TIMEOUT 600 CACHE STRING "Default timeout for tests in seconds")

# Validate test data directory
set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
if(NOT EXISTS "${TEST_DATA_DIR}")
    message(WARNING "Test data directory ${TEST_DATA_DIR} not found. Creating minimal structure.")
    file(MAKE_DIRECTORY "${TEST_DATA_DIR}/meshes")
    file(MAKE_DIRECTORY "${TEST_DATA_DIR}/forcing")
    
    # Create a minimal README for test data
    file(WRITE "${TEST_DATA_DIR}/README.md" 
"# FESOM2 Test Data

This directory should contain minimal test data for FESOM2 integration tests:

## Required Structure:
- `meshes/`: Small test mesh files
- `forcing/`: Minimal forcing/climate data files

## Setup:
Copy minimal test datasets from your FESOM2 data repository or create
symbolic links to existing test data.

Example:
```bash
# Link to existing test mesh
ln -s /path/to/fesom2/test/meshes/test_mesh ./meshes/test_mesh

# Link to minimal forcing data  
ln -s /path/to/minimal/forcing/data ./forcing/
```

The tests will automatically configure namelists to point to these directories.
")
endif()

# Print test configuration
message(STATUS "FESOM2 Testing Configuration:")
message(STATUS "  Integration tests: ${ENABLE_INTEGRATION_TESTS}")
message(STATUS "  MPI tests: ${ENABLE_MPI_TESTS}")
message(STATUS "  Test timeout: ${TEST_TIMEOUT}s")
message(STATUS "  Test data directory: ${TEST_DATA_DIR}")
message(STATUS "  MPI executable: ${MPIEXEC_EXECUTABLE}")

# Add integration tests if enabled
if(ENABLE_INTEGRATION_TESTS)
    add_subdirectory(integration)
endif()

# Add unit tests (for future expansion)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/unit/CMakeLists.txt")
    add_subdirectory(unit)
endif()

# Create a convenience target to run all tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running all FESOM2 tests"
    VERBATIM
)

# Create target to run only integration tests
add_custom_target(run_integration_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "integration_"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running FESOM2 integration tests"
    VERBATIM
)

# Create target to run only MPI tests
if(ENABLE_MPI_TESTS)
    add_custom_target(run_mpi_tests
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "_mpi"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running FESOM2 MPI tests"
        VERBATIM
    )
endif()
