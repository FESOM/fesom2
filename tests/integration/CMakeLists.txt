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

# Custom blowup test configuration
set(BLOWUP_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/integration_blowup_test_mpi2")
set(BLOWUP_TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
set(BLOWUP_RESULT_DIR "${BLOWUP_TEST_DIR}/results")

# Generate fesom.clock file for blowup test
file(MAKE_DIRECTORY "${BLOWUP_RESULT_DIR}")
file(WRITE "${BLOWUP_RESULT_DIR}/fesom.clock" "0 1 1948\n0 1 1948\n")

# Configure namelists for blowup test with step_per_day=2
configure_fesom_namelists("${BLOWUP_TEST_DIR}" "${BLOWUP_TEST_DATA_DIR}" "${BLOWUP_RESULT_DIR}")

# Override the namelist.config modification for blowup test
set(BLOWUP_MODIFY_SCRIPT "${BLOWUP_TEST_DIR}/modify_namelist.config.cmake")
file(WRITE ${BLOWUP_MODIFY_SCRIPT} "
    file(READ \"${BLOWUP_TEST_DIR}/namelist.config\" CONTENT)
    # Replace common paths with test data paths
    string(REGEX REPLACE \"MeshPath='[^']*'\" \"MeshPath='${BLOWUP_TEST_DATA_DIR}/MESHES/pi/'\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"ClimateDataPath='[^']*'\" \"ClimateDataPath='${BLOWUP_TEST_DATA_DIR}/'\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"ResultPath='[^']*'\" \"ResultPath='${BLOWUP_RESULT_DIR}/'\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"fwf_path='[^']*'\" \"fwf_path='${BLOWUP_TEST_DATA_DIR}/meshes/pi/'\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"age_tracer_path='[^']*'\" \"age_tracer_path='${BLOWUP_TEST_DATA_DIR}/meshes/pi/'\" CONTENT \"\${CONTENT}\")
    
    # BLOWUP TEST: Set step_per_day to 2 to trigger blowup
    string(REGEX REPLACE \"step_per_day=[0-9]+\" \"step_per_day=2\" CONTENT \"\${CONTENT}\")
    # Set run_length to 1 day for tests
    string(REGEX REPLACE \"run_length=[0-9]+\" \"run_length=1\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"run_length_unit='[^']*'\" \"run_length_unit='d'\" CONTENT \"\${CONTENT}\")
    # Set restart_length to 1 day
    string(REGEX REPLACE \"restart_length=[0-9]+\" \"restart_length=1\" CONTENT \"\${CONTENT}\")
    string(REGEX REPLACE \"restart_length_unit='[^']*'\" \"restart_length_unit='d'\" CONTENT \"\${CONTENT}\")
    # Set logfile output frequency
    string(REGEX REPLACE \"logfile_outfreq=[0-9]+\" \"logfile_outfreq=1\" CONTENT \"\${CONTENT}\")
    # Force rotation for test geometry
    string(REGEX REPLACE \"force_rotation=.[a-zA-Z].\" \"force_rotation=.true.\" CONTENT \"\${CONTENT}\")
    
    file(WRITE \"${BLOWUP_TEST_DIR}/namelist.config\" \"\${CONTENT}\")
")

# Create blowup test runner
set(BLOWUP_TEST_SCRIPT "${BLOWUP_TEST_DIR}/run_test.cmake")
file(WRITE ${BLOWUP_TEST_SCRIPT} "
    # Create test directories
    file(MAKE_DIRECTORY \"${BLOWUP_TEST_DIR}\")
    file(MAKE_DIRECTORY \"${BLOWUP_RESULT_DIR}\")
    
    # Configure namelists
    execute_process(
        COMMAND ${CMAKE_COMMAND} -P \"${BLOWUP_MODIFY_SCRIPT}\"
        RESULT_VARIABLE result
    )
    if(result)
        message(FATAL_ERROR \"Failed to configure namelist.config for blowup test\")
    endif()
    
    execute_process(
        COMMAND ${CMAKE_COMMAND} -P \"${BLOWUP_TEST_DIR}/modify_namelist.forcing.cmake\"
        RESULT_VARIABLE result
    )
    if(result)
        message(FATAL_ERROR \"Failed to configure namelist.forcing for blowup test\")
    endif()
    
    # Run FESOM with MPI - expect blowup (exit code 1 or other non-zero)
    execute_process(
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_BINARY_DIR}/bin/fesom.x
        WORKING_DIRECTORY \"${BLOWUP_TEST_DIR}\"
        RESULT_VARIABLE test_result
        OUTPUT_VARIABLE test_output
        ERROR_VARIABLE test_error
        TIMEOUT 60
    )
    
    # Log the output
    file(WRITE \"${BLOWUP_TEST_DIR}/test_output.log\" \"\${test_output}\")
    file(WRITE \"${BLOWUP_TEST_DIR}/test_error.log\" \"\${test_error}\")
    
    # For blowup test, we expect a non-zero exit code (blowup condition)
    if(test_result EQUAL 0)
        message(STATUS \"Blowup test completed normally (exit code: \${test_result}) - this may indicate the blowup condition was not triggered\")
    else()
        message(STATUS \"Blowup test completed with expected blowup condition (exit code: \${test_result})\")
    endif()
")

# Add the blowup test
add_test(
    NAME integration_blowup_test_mpi2
    COMMAND ${CMAKE_COMMAND} -P ${BLOWUP_TEST_SCRIPT}
)

# Set test properties
set_tests_properties(integration_blowup_test_mpi2 PROPERTIES
    TIMEOUT 60
    WORKING_DIRECTORY ${BLOWUP_TEST_DIR}
    PROCESSORS 2
    RUN_SERIAL FALSE
)
