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

    # JRA55 atmospheric forcing: run the pi mesh for one day driven by a JRA55-do
    # dataset (tests/data/FORCING/JRA55, year 1958, 3-hourly) instead of the default
    # CORE2 forcing. This exercises the JRA55 namelist wiring end to end:
    #   * FORCING JRA        -> namelist.forcing from config/namelist.forcing.JRA,
    #                           pool paths rewired to that tree (CMOR var names
    #                           uas/tas/prra/..., relative to ClimateDataPath)
    #   * FORCING_YEAR 1958  -> fesom.clock start year, per-year forcing file, and
    #                           the required output artifact sst.fesom.1958.nc
    #   * LEAPYEAR .true.    -> include_fleapyear=.true., mandatory for JRA55's
    #                           gregorian time axis (fesom.x aborts otherwise; see
    #                           the calendar-consistency check in gen_surface_forcing.F90)
    # NOTE: passes in BOTH double- and single-precision builds (verified on pi/1-day).
    # The JRA55 SP-NaN concern is specific to large meshes (e.g. ng5), not the finer
    # sub-daily cadence itself, so it does not surface on the small pi mesh here.
    #
    # The dataset itself is NOT bundled: 13 files / ~120 MB is too much to carry in the
    # repository, and it would land in git history permanently. It will be provided as a
    # downloadable dataset instead, like the meshes the remote_ tests fetch.
    #
    # Until then this test is simply not registered. All of the harness wiring for it --
    # FORCING/FORCING_YEAR/LEAPYEAR in add_fesom_test, configure_fesom_namelists_with_options
    # selecting config/namelist.forcing.JRA, and the pool-path rewriting in
    # cmake/FesomTesting.cmake -- is kept and unchanged, so dropping the data into
    # tests/data/FORCING/JRA55/ (or having a fetch step put it there) re-enables the test
    # with no further edits.
    if(EXISTS "${TEST_DATA_DIR}/FORCING/JRA55")
        add_fesom_test(integration_pi_jra_mpi2
            MPI_TEST
            NP 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL jra
            FORCING "JRA"
            FORCING_YEAR "1958"
            LEAPYEAR ".true."
        )
    else()
        message(STATUS "JRA55 forcing data not present (${TEST_DATA_DIR}/FORCING/JRA55) -- skipping integration_pi_jra_mpi2. Harness support for FORCING=JRA is unchanged; add the dataset to enable it.")
    endif()

    # Note: there is no separate "SP" integration test. Precision is a property of
    # the build, not the test -- in a USE_SINGLE_PRECISION build every test above runs
    # single precision, and each fesom.x test auto-asserts the build's precision banner
    # (see add_fesom_test_with_options in cmake/FesomTesting.cmake). Run the suite in an
    # SP build dir to get SP coverage of everything, including the CVMix schemes.

    # Regression test for ghost output records on time-step change
    add_fesom_output_restamp_test(integration_pi_output_restamp
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
    )

    # Idealized (toy_ocean) configurations. They supply their own forcing and
    # read only a mesh, so they need none of the bundled forcing or climatology
    # data; the meshes are used in place from test/meshes/ via MESH_DIR. Values
    # follow setups/test_souf and setups/test_neverworld2, i.e. what CI runs.
    #
    # Pass/fail is the clean run plus the sst artifact, so these show the
    # configurations still start and complete. They do not compare answers: a
    # numerical drift would go unnoticed. Adding that is designed in
    # docs/pr/PLAN-ctest-coverage.md.
    include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FesomToyTest.cmake)

    add_fesom_toy_test(integration_toy_soufflet_mpi8
        WHICH_TOY "soufflet"
        MESH_DIR "${FESOM_TESTING_ROOT}/test/meshes/soufflet"
        CYCLIC_LENGTH "4.5"
        NP 8
        TIMEOUT ${TEST_TIMEOUT}
    )

    add_fesom_toy_test(integration_toy_neverworld2_mpi8
        WHICH_TOY "neverworld2"
        MESH_DIR "${FESOM_TESTING_ROOT}/test/meshes/neverworld2"
        CYCLIC_LENGTH "60.0"
        NP 8
        TIMEOUT ${TEST_TIMEOUT}
    )

    # Ocean-only run: the pi mesh with the sea-ice model switched off, otherwise
    # identical to integration_pi_mpi2 (CORE2 forcing, WOA18 initial state, 1 day).
    # This is the only ice-free configuration in the suite that is not a toy_ocean
    # setup, so it is the only one exercising the ocean-side setup and forcing paths
    # with use_ice=.false.
    #
    # Pass/fail is the clean run plus the sst artifact. That detects a crash or an
    # abort in the ice-free path; it does NOT detect an ice-free run that completes
    # while the ocean receives no surface forcing at all. Checking for that needs a
    # non-zero assertion on the forcing diagnostics (tx_sur/fh/fw) and an ice-on
    # control run to keep the ice-output absence checks meaningful; see
    # docs/pr/PLAN-ctest-coverage.md.
    add_fesom_test(integration_pi_ice_free_mpi2
        MPI_TEST
        NP 2
        TIMEOUT ${TEST_TIMEOUT}
        LABEL ice_free
        USE_ICE ".false."
    )

    if(BUILD_MESHPARTITIONER)
        set(FRESH_DATA_DIR "${CMAKE_BINARY_DIR}/tests/fresh_data")
        file(MAKE_DIRECTORY "${FRESH_DATA_DIR}/MESHES")
        file(GLOB DATA_ENTRIES LIST_DIRECTORIES true "${TEST_DATA_DIR}/*")
        foreach(DATA_ENTRY IN LISTS DATA_ENTRIES)
            get_filename_component(DATA_NAME "${DATA_ENTRY}" NAME)
            if(NOT DATA_NAME STREQUAL "MESHES" AND IS_DIRECTORY "${DATA_ENTRY}" AND NOT EXISTS "${FRESH_DATA_DIR}/${DATA_NAME}")
                file(CREATE_LINK "${DATA_ENTRY}" "${FRESH_DATA_DIR}/${DATA_NAME}" SYMBOLIC)
            endif()
        endforeach()
        foreach(PARTITION_NAME flat_1 hier_1x1 hier_1x2)
            set(PARTITION_MESH_DIR "${CMAKE_BINARY_DIR}/tests/meshpartitioner/meshes/pi_${PARTITION_NAME}")
            if(NOT EXISTS "${PARTITION_MESH_DIR}")
                file(MAKE_DIRECTORY "${PARTITION_MESH_DIR}")
            endif()
            if(NOT EXISTS "${FRESH_DATA_DIR}/MESHES/pi_${PARTITION_NAME}")
                file(CREATE_LINK "${PARTITION_MESH_DIR}" "${FRESH_DATA_DIR}/MESHES/pi_${PARTITION_NAME}" SYMBOLIC)
            endif()
        endforeach()

        add_fesom_test_with_options(integration_pi_fresh_flat_1_mpi1_omp2
            "pi_flat_1" "96" "1" "s" "1" "s" "10" ".true." ".false."
            NP 1
            OMP_THREADS 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL fresh_partition
            DATA_DIR "${FRESH_DATA_DIR}"
        )
        set_tests_properties(integration_pi_fresh_flat_1_mpi1_omp2 PROPERTIES FIXTURES_REQUIRED "partition_pi_flat_1")

        add_fesom_test_with_options(integration_pi_fresh_hier_1x1_mpi1_omp2
            "pi_hier_1x1" "96" "1" "s" "1" "s" "10" ".true." ".false."
            NP 1
            OMP_THREADS 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL fresh_partition
            DATA_DIR "${FRESH_DATA_DIR}"
        )
        set_tests_properties(integration_pi_fresh_hier_1x1_mpi1_omp2 PROPERTIES FIXTURES_REQUIRED "partition_pi_hier_1x1")

        add_fesom_test_with_options(integration_pi_fresh_hier_1x2_mpi2_omp2
            "pi_hier_1x2" "96" "1" "s" "1" "s" "10" ".true." ".false."
            MPI_TEST
            NP 2
            OMP_THREADS 2
            TIMEOUT ${TEST_TIMEOUT}
            LABEL fresh_partition
            DATA_DIR "${FRESH_DATA_DIR}"
        )
        set_tests_properties(integration_pi_fresh_hier_1x2_mpi2_omp2 PROPERTIES FIXTURES_REQUIRED "partition_pi_hier_1x2")
    endif()

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
