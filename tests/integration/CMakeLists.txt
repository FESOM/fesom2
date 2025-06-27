#===============================================================================
# tests/integration/CMakeLists.txt - Integration tests for FESOM2
#===============================================================================

# MPI tests (if enabled)
if(ENABLE_MPI_TESTS)
    # Basic integration test with 2 processes (minimum required for mesh partitioning)
    add_fesom_test(integration_basic_run_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
    )
    
    # Short duration test to check initialization with 2 processes
    add_fesom_test(integration_init_test_mpi2
        MPI_TEST
        NP 2
        TIMEOUT 120  # 2 minutes for quick init test
    )
    
    # MPI test with 8 processes (as requested by user)
    add_fesom_test(integration_basic_run_mpi8
        MPI_TEST  
        NP 8
        TIMEOUT ${TEST_TIMEOUT}
    )
    
    # Additional MPI test configurations can be added here
    # For example, testing different domain decompositions
else()
    message(WARNING "MPI tests are disabled. Set ENABLE_MPI_TESTS=ON to enable them.")
endif()

# Create a test that verifies namelist configuration works
add_test(
    NAME integration_namelist_config
    COMMAND ${CMAKE_COMMAND}
        -DTEST_DATA_DIR=${TEST_DATA_DIR}
        -DCONFIG_DIR=${CMAKE_SOURCE_DIR}/config
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_namelist_config.cmake
)

set_tests_properties(integration_namelist_config PROPERTIES
    TIMEOUT 30
)

# Test that checks if test data is properly set up
add_test(
    NAME integration_test_data_check
    COMMAND ${CMAKE_COMMAND}
        -DTEST_DATA_DIR=${TEST_DATA_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test_data_check.cmake
)

set_tests_properties(integration_test_data_check PROPERTIES
    TIMEOUT 10
)
