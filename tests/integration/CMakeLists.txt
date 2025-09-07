#===============================================================================
# tests/integration/CMakeLists.txt - Integration tests for FESOM2
#===============================================================================

# MPI tests (if enabled)
if(ENABLE_MPI_TESTS)
    # Basic integration test with 2 processes (minimum required for mesh partitioning)
    add_fesom_test(integration_full_run_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
    )
    
    # MPI test with 8 processes (as requested by user)
    add_fesom_test(integration_full_run_mpi8
        MPI_TEST  
        NP 8
        TIMEOUT ${TEST_TIMEOUT}
    )
    
    # Pi cavity test with 2 processes (cavity enabled) - only dist_2 available
    add_fesom_test_with_options(integration_pi_cavity_mpi2
        "pi_cavity" "96" "1" "d" "1" "d" "10" ".true." ".true."
        MPI_TEST
        NP 2
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

# Meshdiag tests (only when BUILD_MESHDIAG is enabled)
if(BUILD_MESHDIAG AND TARGET fesom_meshdiag)
    message(STATUS "Adding meshdiag tests...")
    
    # Basic meshdiag test with 2 processes (matches integration_full_run_mpi2)
    add_fesom_meshdiag_test(meshdiag_pi_mpi2
        NP 2
        TIMEOUT 120
    )
    
    # Meshdiag test with 8 processes (matches integration_full_run_mpi8)
    add_fesom_meshdiag_test_with_options(meshdiag_pi_mpi8
        "pi" "mesh8"
        NP 8
        TIMEOUT 120
    )
    
else()
    message(STATUS "Meshdiag tests disabled. Set BUILD_MESHDIAG=ON to enable them.")
endif()

#===============================================================================
# Generic Mesh Pipeline Tests (New System)
#===============================================================================

# Include the new mesh pipeline tests
include(${CMAKE_CURRENT_SOURCE_DIR}/mesh_pipeline_tests.cmake)

#===============================================================================
# Legacy Core2 Tests (for backward compatibility)
#===============================================================================

# Legacy core2-specific tests removed - replaced by generic mesh pipeline system
