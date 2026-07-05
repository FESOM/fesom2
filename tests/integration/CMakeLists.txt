#===============================================================================
# tests/integration/CMakeLists.txt - Integration tests for FESOM2
#===============================================================================

# MPI tests (if enabled)
#
# These run a plain standalone fesom.x. A build that activates an external
# coupler at runtime cannot run the standalone pi test and aborts at startup:
#   * FESOM_COUPLED / OIFS_COUPLED -> OASIS, aborts without 'namcouple'
#   * USE_YAC                      -> YAC,   aborts without 'coupling.yaml'
# Skip registering the standalone integration tests for those configurations
# (build-only verification).
if(ENABLE_MPI_TESTS AND NOT FESOM_COUPLED AND NOT OIFS_COUPLED AND NOT USE_YAC)
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
    
    # CVMix vertical mixing: run the pi mesh with CVMix schemes so the CVMix code
    # path is exercised end to end. NOTE: the stock namelist.oce ships
    # mix_scheme='KPP' -- FESOM's native (non-CVMix) KPP -- so a CVMix-enabled
    # build (CVMIX=ON) never runs CVMix by default; the scheme must be selected
    # explicitly. These tests do that via MIX_SCHEME.
    #
    # CVMix is a fixed double-precision library (kind cvmix_r8 = real64) built
    # into FESOM's tree; the FESOM-side WP<->cvmix_r8 conversion is what differs
    # between builds. So each test covers the previously-untested FESOM-DP+CVMix
    # path in a DP build and the FESOM-SP+CVMix path (the conversion shims) in an
    # SP build. Registered only when CVMix is compiled in (CVMIX=ON, the default).
    # Pass/fail is the clean run + artifacts: a wrong scheme string or missing
    # CVMix aborts at setup.
    if(CVMIX)
        add_fesom_test(integration_pi_cvmix_tke_mpi2
            MPI_TEST
            NP 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL cvmix
            MIX_SCHEME "cvmix_TKE"
        )
        add_fesom_test(integration_pi_cvmix_kpp_mpi2
            MPI_TEST
            NP 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL cvmix
            MIX_SCHEME "cvmix_KPP"
        )
    endif()

    # Note: there is no separate "SP" integration test. Precision is a property of
    # the build, not the test -- in a USE_SINGLE_PRECISION build every test above runs
    # single precision, and each fesom.x test auto-asserts the build's precision banner
    # (see add_fesom_test_with_options in cmake/FesomTesting.cmake). Run the suite in an
    # SP build dir to get SP coverage of everything, including the CVMix schemes.

    # Additional MPI test configurations can be added here
    # For example, testing different domain decompositions
elseif(NOT ENABLE_MPI_TESTS)
    message(WARNING "MPI tests are disabled. Set ENABLE_MPI_TESTS=ON to enable them.")
else()
    message(STATUS "Standalone integration tests skipped for coupled build (FESOM_COUPLED/OIFS_COUPLED/USE_YAC): a coupled fesom.x requires an external coupler setup (namcouple/coupling.yaml) and cannot run the standalone pi test.")
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
