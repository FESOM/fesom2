#===============================================================================
# tests/ecbundle/CMakeLists.txt - ecbundle build tests for FESOM2
#===============================================================================

# Check if ecbundle is available
find_program(ECBUNDLE_CREATE ecbundle-create)
find_program(ECBUNDLE_BUILD ecbundle-build)

if(NOT ECBUNDLE_CREATE OR NOT ECBUNDLE_BUILD)
    message(WARNING "ecbundle tools not found. Skipping ecbundle tests.")
    message(STATUS "To install ecbundle, see: https://github.com/ecmwf/ecbundle")
    return()
endif()

message(STATUS "Found ecbundle tools:")
message(STATUS "  ecbundle-create: ${ECBUNDLE_CREATE}")
message(STATUS "  ecbundle-build: ${ECBUNDLE_BUILD}")

# Create test directory
set(ECBUNDLE_TEST_DIR "${CMAKE_BINARY_DIR}/test_ecbundle")
file(MAKE_DIRECTORY "${ECBUNDLE_TEST_DIR}")

# Test 1: Basic ecbundle build (default options)
add_test(
    NAME ecbundle_basic_build
    COMMAND ${CMAKE_COMMAND}
        -DECBUNDLE_CREATE=${ECBUNDLE_CREATE}
        -DECBUNDLE_BUILD=${ECBUNDLE_BUILD}
        -DBUNDLE_FILE=${CMAKE_SOURCE_DIR}/bundle.yml
        -DTEST_DIR=${ECBUNDLE_TEST_DIR}/basic
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_ecbundle_basic.cmake
)

set_tests_properties(ecbundle_basic_build PROPERTIES
    TIMEOUT 1800  # 30 minutes timeout
    LABELS "ecbundle"
)

# Test 2: ecbundle build with meshpartitioner enabled
add_test(
    NAME ecbundle_with_meshpart
    COMMAND ${CMAKE_COMMAND}
        -DECBUNDLE_CREATE=${ECBUNDLE_CREATE}
        -DECBUNDLE_BUILD=${ECBUNDLE_BUILD}
        -DBUNDLE_FILE=${CMAKE_SOURCE_DIR}/bundle.yml
        -DTEST_DIR=${ECBUNDLE_TEST_DIR}/meshpart
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_ecbundle_meshpart.cmake
)

set_tests_properties(ecbundle_with_meshpart PROPERTIES
    TIMEOUT 1800  # 30 minutes timeout
    LABELS "ecbundle"
)

# Test 3: ecbundle build with OpenMP support
add_test(
    NAME ecbundle_with_omp
    COMMAND ${CMAKE_COMMAND}
        -DECBUNDLE_CREATE=${ECBUNDLE_CREATE}
        -DECBUNDLE_BUILD=${ECBUNDLE_BUILD}
        -DBUNDLE_FILE=${CMAKE_SOURCE_DIR}/bundle.yml
        -DTEST_DIR=${ECBUNDLE_TEST_DIR}/omp
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_ecbundle_omp.cmake
)

set_tests_properties(ecbundle_with_omp PROPERTIES
    TIMEOUT 1800  # 30 minutes timeout
    LABELS "ecbundle"
)

# Test 4: ecbundle build with testing enabled
add_test(
    NAME ecbundle_with_testing
    COMMAND ${CMAKE_COMMAND}
        -DECBUNDLE_CREATE=${ECBUNDLE_CREATE}
        -DECBUNDLE_BUILD=${ECBUNDLE_BUILD}
        -DBUNDLE_FILE=${CMAKE_SOURCE_DIR}/bundle.yml
        -DTEST_DIR=${ECBUNDLE_TEST_DIR}/testing
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_ecbundle_testing.cmake
)

set_tests_properties(ecbundle_with_testing PROPERTIES
    TIMEOUT 1800  # 30 minutes timeout
    LABELS "ecbundle"
)

# Create a convenience target to run all ecbundle tests
add_custom_target(run_ecbundle_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L ecbundle
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running FESOM2 ecbundle tests"
    VERBATIM
) 