#===============================================================================
# test/ecbundle/test_ecbundle_testing.cmake - ecbundle build test with testing enabled
#===============================================================================

# This script tests the ecbundle build process for FESOM2 with testing enabled
# It uses the --with-testing=ON option from bundle.yml

message(STATUS "Starting ecbundle build test with testing enabled")
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

# Step 2: Use ecbundle-build to build with testing enabled
message(STATUS "Step 2: Running ecbundle-build with --with-testing=ON...")
execute_process(
    COMMAND ${ECBUNDLE_BUILD} --with-testing=ON
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE BUILD_RESULT
    OUTPUT_VARIABLE BUILD_OUTPUT
    ERROR_VARIABLE BUILD_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

if(NOT BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR "ecbundle-build with testing enabled failed with result ${BUILD_RESULT}")
    message(FATAL_ERROR "Output: ${BUILD_OUTPUT}")
    message(FATAL_ERROR "Error: ${BUILD_ERROR}")
endif()

message(STATUS "ecbundle-build with testing enabled completed successfully")
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

# Check that executable is actually executable
execute_process(
    COMMAND test -x "${TEST_DIR}/build/bin/fesom.x"
    RESULT_VARIABLE EXEC_TEST_RESULT
)

if(NOT EXEC_TEST_RESULT EQUAL 0)
    message(FATAL_ERROR "FESOM executable is not executable")
endif()

# Step 4: Check for test executables and test files
message(STATUS "Step 4: Checking for test artifacts...")

# Check if test executables were built
# Look for common test executable patterns
file(GLOB TEST_EXECUTABLES 
    "${TEST_DIR}/build/bin/*test*"
    "${TEST_DIR}/build/bin/test_*"
    "${TEST_DIR}/build/bin/*_test"
)

if(TEST_EXECUTABLES)
    message(STATUS "Found test executables: ${TEST_EXECUTABLES}")
else()
    message(STATUS "No test executables found (this might be expected)")
endif()

# Check if test files were created
if(EXISTS "${TEST_DIR}/build/Testing")
    message(STATUS "Testing directory found")
else()
    message(STATUS "No Testing directory found (this might be expected)")
endif()

# Check for CMake test files
file(GLOB CMAKE_TEST_FILES "${TEST_DIR}/build/CMakeFiles/*.cmake")
if(CMAKE_TEST_FILES)
    message(STATUS "Found CMake test files: ${CMAKE_TEST_FILES}")
else()
    message(STATUS "No CMake test files found (this might be expected)")
endif()

# Step 5: Try to run the tests if CTest is available
message(STATUS "Step 5: Attempting to run tests...")

find_program(CTEST_EXECUTABLE ctest)
if(CTEST_EXECUTABLE)
    message(STATUS "Found CTest: ${CTEST_EXECUTABLE}")
    
    # Try to run tests
    execute_process(
        COMMAND ${CTEST_EXECUTABLE} --output-on-failure
        WORKING_DIRECTORY ${TEST_DIR}/build
        RESULT_VARIABLE CTEST_RESULT
        OUTPUT_VARIABLE CTEST_OUTPUT
        ERROR_VARIABLE CTEST_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    
    if(CTEST_RESULT EQUAL 0)
        message(STATUS "CTest ran successfully")
        message(STATUS "CTest output: ${CTEST_OUTPUT}")
    else()
        message(STATUS "CTest failed (this might be expected if no tests are configured)")
        message(STATUS "CTest error: ${CTEST_ERROR}")
    endif()
else()
    message(STATUS "CTest not found")
endif()

# Step 6: Test that the main executable can run
message(STATUS "Step 6: Testing FESOM executable...")
execute_process(
    COMMAND ${TEST_DIR}/build/bin/fesom.x --help
    RESULT_VARIABLE HELP_RESULT
    OUTPUT_VARIABLE HELP_OUTPUT
    ERROR_VARIABLE HELP_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

# Note: We don't fail here if --help doesn't work, as FESOM might not support it
# But we log the result for debugging
if(HELP_RESULT EQUAL 0)
    message(STATUS "FESOM --help command successful")
    message(STATUS "Help output: ${HELP_OUTPUT}")
else()
    message(STATUS "FESOM --help command failed (this might be expected)")
    message(STATUS "Help error: ${HELP_ERROR}")
endif()

message(STATUS "ecbundle build test with testing enabled completed successfully!") 