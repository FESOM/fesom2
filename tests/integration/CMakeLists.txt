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

    # JRA55 atmospheric forcing: run the pi mesh for one day driven by the bundled
    # JRA55-do dataset (tests/data/FORCING/JRA55, year 1958, 3-hourly) instead of the
    # default CORE2 forcing. This exercises the JRA55 namelist wiring end to end:
    #   * FORCING JRA        -> namelist.forcing from config/namelist.forcing.JRA,
    #                           pool paths rewired to the bundled tree (CMOR var
    #                           names uas/tas/prra/..., relative to ClimateDataPath)
    #   * FORCING_YEAR 1958  -> fesom.clock start year, per-year forcing file, and
    #                           the required output artifact sst.fesom.1958.nc
    #   * LEAPYEAR .true.    -> include_fleapyear=.true., mandatory for JRA55's
    #                           gregorian time axis (fesom.x aborts otherwise; see
    #                           the calendar-consistency check in gen_surface_forcing.F90)
    # NOTE: passes in BOTH double- and single-precision builds (verified on pi/1-day).
    # The JRA55 SP-NaN concern is specific to large meshes (e.g. ng5), not the finer
    # sub-daily cadence itself, so it does not surface on the small pi mesh here.
    add_fesom_test(integration_pi_jra_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
        LABEL jra
        FORCING "JRA"
        FORCING_YEAR "1958"
        LEAPYEAR ".true."
    )

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
