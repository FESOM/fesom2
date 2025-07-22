#===============================================================================
# test/ecbundle/test_ecbundle_meshpart.cmake - ecbundle build test with meshpartitioner
#===============================================================================

# This script tests the ecbundle build process for FESOM2 with meshpartitioner enabled
# It uses the --with-meshpart=ON option from bundle.yml

message(STATUS "Starting ecbundle build test with meshpartitioner")
message(STATUS "Test directory: ${TEST_DIR}")
message(STATUS "Bundle file: ${BUNDLE_FILE}")

# Create test directory
file(MAKE_DIRECTORY "${TEST_DIR}")

# Step 1: Use ecbundle-create to checkout sources
message(STATUS "Step 1: Running ecbundle-create...")
execute_process(
    COMMAND ${ECBUNDLE_CREATE} --bundle ${BUNDLE_FILE}
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE CREATE_RESULT
    OUTPUT_VARIABLE CREATE_OUTPUT
    ERROR_VARIABLE CREATE_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

if(NOT CREATE_RESULT EQUAL 0)
    message(FATAL_ERROR "ecbundle-create failed with result ${CREATE_RESULT}")
    message(FATAL_ERROR "Output: ${CREATE_OUTPUT}")
    message(FATAL_ERROR "Error: ${CREATE_ERROR}")
endif()

message(STATUS "ecbundle-create completed successfully")
message(STATUS "Output: ${CREATE_OUTPUT}")

# Verify that source directory was created
if(NOT EXISTS "${TEST_DIR}/source")
    message(FATAL_ERROR "source directory was not created by ecbundle-create")
endif()

# Verify that FESOM sources are present
if(NOT EXISTS "${TEST_DIR}/source/fesom")
    message(FATAL_ERROR "FESOM sources not found in ${TEST_DIR}/source/fesom")
endif()

# Step 2: Use ecbundle-build to build with meshpartitioner enabled
message(STATUS "Step 2: Running ecbundle-build with --with-meshpart=ON...")
execute_process(
    COMMAND ${ECBUNDLE_BUILD} --with-meshpart=ON
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE BUILD_RESULT
    OUTPUT_VARIABLE BUILD_OUTPUT
    ERROR_VARIABLE BUILD_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

if(NOT BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR "ecbundle-build with meshpartitioner failed with result ${BUILD_RESULT}")
    message(FATAL_ERROR "Output: ${BUILD_OUTPUT}")
    message(FATAL_ERROR "Error: ${BUILD_ERROR}")
endif()

message(STATUS "ecbundle-build with meshpartitioner completed successfully")
message(STATUS "Output: ${BUILD_OUTPUT}")

# Step 3: Verify that the build produced the expected artifacts
message(STATUS "Step 3: Verifying build artifacts...")

# Check for build directory
if(NOT EXISTS "${TEST_DIR}/build")
    message(FATAL_ERROR "build directory was not created by ecbundle-build")
endif()

# Check for FESOM executable (should be in build/bin)
if(NOT EXISTS "${TEST_DIR}/build/bin/fesom.x")
    message(FATAL_ERROR "FESOM executable not found in ${TEST_DIR}/build/bin/fesom.x")
endif()

# Check for meshpartitioner executable (should be in build/bin)
if(NOT EXISTS "${TEST_DIR}/build/bin/fesom_meshpart")
    message(FATAL_ERROR "meshpartitioner executable not found in ${TEST_DIR}/build/bin/fesom_meshpart")
endif()

# Check that executables are actually executable
execute_process(
    COMMAND test -x "${TEST_DIR}/build/bin/fesom.x"
    RESULT_VARIABLE FESOM_EXEC_TEST_RESULT
)

if(NOT FESOM_EXEC_TEST_RESULT EQUAL 0)
    message(FATAL_ERROR "FESOM executable is not executable")
endif()

execute_process(
    COMMAND test -x "${TEST_DIR}/build/bin/fesom_meshpart"
    RESULT_VARIABLE MESHPART_EXEC_TEST_RESULT
)

if(NOT MESHPART_EXEC_TEST_RESULT EQUAL 0)
    message(FATAL_ERROR "meshpartitioner executable is not executable")
endif()

# Step 4: Test that the meshpartitioner executable can run
message(STATUS "Step 4: Testing meshpartitioner executable...")
execute_process(
    COMMAND ${TEST_DIR}/build/bin/fesom_meshpart --help
    RESULT_VARIABLE MESHPART_HELP_RESULT
    OUTPUT_VARIABLE MESHPART_HELP_OUTPUT
    ERROR_VARIABLE MESHPART_HELP_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

# Note: We don't fail here if --help doesn't work, as meshpartitioner might not support it
# But we log the result for debugging
if(MESHPART_HELP_RESULT EQUAL 0)
    message(STATUS "meshpartitioner --help command successful")
    message(STATUS "Help output: ${MESHPART_HELP_OUTPUT}")
else()
    message(STATUS "meshpartitioner --help command failed (this might be expected)")
    message(STATUS "Help error: ${MESHPART_HELP_ERROR}")
endif()

# Step 5: Test that the FESOM executable can run
message(STATUS "Step 5: Testing FESOM executable...")
execute_process(
    COMMAND ${TEST_DIR}/build/bin/fesom.x --help
    RESULT_VARIABLE FESOM_HELP_RESULT
    OUTPUT_VARIABLE FESOM_HELP_OUTPUT
    ERROR_VARIABLE FESOM_HELP_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

# Note: We don't fail here if --help doesn't work, as FESOM does not support it yet
# But we log the result for debugging
if(FESOM_HELP_RESULT EQUAL 0)
    message(STATUS "FESOM --help command successful")
    message(STATUS "Help output: ${FESOM_HELP_OUTPUT}")
else()
    message(STATUS "FESOM --help command failed (this might be expected)")
    message(STATUS "Help error: ${FESOM_HELP_ERROR}")
endif()

message(STATUS "ecbundle build test with meshpartitioner completed successfully!") 
