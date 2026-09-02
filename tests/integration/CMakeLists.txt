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
    #
    # The EXTRA_SUCCESS_MARKERS pin the always-on SSH-CG solver diagnostics. Both the
    # end-of-run summary and the per-step line are plain writes with no artifact behind
    # them, so if either stops being emitted nothing else in the suite would notice --
    # which is exactly how the sprod(2) defect in this solver survived unfixed. Matching
    # the labels (not the numbers) keeps this stable across meshes and precisions.
    add_fesom_test(integration_pi_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
        EXTRA_SUCCESS_MARKERS
            "___SSH CG SOLVER"            # end-of-run summary block
            "non-convergences (maxiter)"  # the stall counter is reported
            "breakdowns (p.Ap<=0)"        # the breakdown guard is wired up
            "ssh_cg iters="               # per-step line, beside cfl_z
            "ssh_cg rms(r)="              # per-step residual vs its target
            "ssh CG soltol"               # the namelist knobs are read and echoed
            "ssh CG maxiter"
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
    
    # Regression test for ghost output records on time-step change (issue #123)
    add_fesom_output_restamp_test(integration_pi_output_restamp
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
    )

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
