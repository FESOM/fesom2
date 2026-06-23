#===============================================================================
# tests/integration/CMakeLists.txt - Integration tests for FESOM2
#===============================================================================

# MPI tests (if enabled)
#
# These run a plain standalone fesom.x. A coupled build (FESOM_COUPLED or
# OIFS_COUPLED) links OASIS and aborts at startup without a coupling setup
# (namcouple), so the standalone pi run is not applicable there - skip
# registering them for coupled configurations (build-only verification).
if(ENABLE_MPI_TESTS AND NOT FESOM_COUPLED AND NOT OIFS_COUPLED)
    # Basic integration test with PI mesh and 2 processes (minimum required for mesh partitioning)
    add_fesom_test(integration_pi_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
    )
    
    # PI mesh test with 8 processes
    add_fesom_test(integration_pi_mpi8
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
elseif(NOT ENABLE_MPI_TESTS)
    message(WARNING "MPI tests are disabled. Set ENABLE_MPI_TESTS=ON to enable them.")
else()
    message(STATUS "Standalone integration tests skipped for coupled build (FESOM_COUPLED/OIFS_COUPLED): a coupled fesom.x requires an OASIS coupling setup and cannot run the standalone pi test.")
endif()


# Meshdiag tests (only when BUILD_MESHDIAG is enabled)
if(BUILD_MESHDIAG AND TARGET fesom_meshdiag)
    message(STATUS "Adding meshdiag tests...")
    
    # Basic meshdiag test with 2 processes (matches integration_full_run_mpi2)
    add_fesom_meshdiag_test_with_options(integration_meshdiag_pi_mpi2
        "pi" "fesom"
        NP 2
        TIMEOUT 120
    )
    
    # Meshdiag test with 8 processes (matches integration_full_run_mpi8)
    add_fesom_meshdiag_test_with_options(integration_meshdiag_pi_mpi8
        "pi" "fesom"
        NP 8
        TIMEOUT 120
    )
    
else()
    message(STATUS "Meshdiag tests disabled. Set BUILD_MESHDIAG=ON to enable them.")
endif()

#===============================================================================
# Local Integration Tests Only
#===============================================================================

# Remote mesh pipeline tests have been moved to tests/remote/
