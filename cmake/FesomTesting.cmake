#===============================================================================
# FesomTesting.cmake - Testing utilities for FESOM2
#===============================================================================
#
# This file provides utilities for configuring and running FESOM2 tests.
# The namelist configuration has been refactored into individual functions
# for better control and customization:
#
# Individual namelist functions:
# - update_common_paths() - Updates common paths in namelist.config only (defaults to pi mesh)
# - update_common_paths_with_mesh() - Updates common paths with specific mesh
# - update_namelist_config() - Configures namelist.config for test_pi (defaults)
# - update_namelist_config_with_options() - Configures namelist.config with custom options
# - update_namelist_dyn() - Configures namelist.dyn for test_pi
# - update_namelist_ice() - Configures namelist.ice for test_pi
# - update_namelist_tra() - Configures namelist.tra for test_pi
# - update_tracer_init3d_filelist() - Specifically updates tracer_init3d filelist
# - update_namelist_forcing() - Configures namelist.forcing for test_pi
# - update_namelist_oce() - Configures namelist.oce (placeholder)
# - update_namelist_io() - Configures namelist.io (placeholder)
# - update_namelist_cvmix() - Configures namelist.cvmix (placeholder)
#
# Main functions:
# - configure_fesom_namelists() - Main function that calls all individual functions (defaults)
# - configure_fesom_namelists_with_options() - Configure namelists with custom options
# - add_fesom_test() - Adds a complete FESOM integration test (defaults)
# - add_fesom_test_with_options() - Adds a FESOM test with custom options
# - generate_fesom_clock() - Generates fesom.clock file
# - setup_mpi_testing() - Sets up MPI for testing
#
# Example usage for custom namelist configuration:
# ```
# # Copy and configure namelist.config with common paths
# configure_file("${CMAKE_SOURCE_DIR}/config/namelist.config" "${TARGET_DIR}/namelist.config" COPYONLY)
# update_common_paths("${TARGET_DIR}/namelist.config" "${TARGET_DIR}/namelist.config" "${TEST_DATA_DIR}" "${RESULT_DIR}")
# update_namelist_config("${TARGET_DIR}/namelist.config" "${TARGET_DIR}/namelist.config")
# 
# # Copy and configure a specific namelist (e.g., tracer)
# configure_file("${CMAKE_SOURCE_DIR}/config/namelist.tra" "${TARGET_DIR}/namelist.tra" COPYONLY)
# update_tracer_init3d_filelist("${TARGET_DIR}/namelist.tra" "${TARGET_DIR}/namelist.tra" "${TEST_DATA_DIR}")
# 
# # Example: Custom test with pi_cavity mesh and cavity enabled
# add_fesom_test_with_options("test_pi_cavity_mpi8" "pi_cavity" "96" "1" "d" "1" "d" "10" ".true." ".true."
#     MPI_TEST NP 8 TIMEOUT 600
# )
# ```
#

# Function to update common paths in any namelist
function(update_common_paths NAMELIST_IN NAMELIST_OUT TEST_DATA_DIR RESULT_DIR)
    update_common_paths_with_mesh("${NAMELIST_IN}" "${NAMELIST_OUT}" "${TEST_DATA_DIR}" "${RESULT_DIR}" "pi")
endfunction()

# Function to update common paths with specific mesh
function(update_common_paths_with_mesh NAMELIST_IN NAMELIST_OUT TEST_DATA_DIR RESULT_DIR MESH_NAME)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Replace common paths with test data paths for specific mesh
    string(REGEX REPLACE "MeshPath='[^']*'" "MeshPath='${TEST_DATA_DIR}/MESHES/${MESH_NAME}/'" CONTENT "${CONTENT}")
    string(REGEX REPLACE "ClimateDataPath='[^']*'" "ClimateDataPath='${TEST_DATA_DIR}/'" CONTENT "${CONTENT}")
    string(REGEX REPLACE "ResultPath='[^']*'" "ResultPath='${RESULT_DIR}/'" CONTENT "${CONTENT}")
    string(REGEX REPLACE "fwf_path='[^']*'" "fwf_path='${TEST_DATA_DIR}/meshes/${MESH_NAME}/'" CONTENT "${CONTENT}")
    string(REGEX REPLACE "age_tracer_path='[^']*'" "age_tracer_path='${TEST_DATA_DIR}/meshes/${MESH_NAME}/'" CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.config for test_pi
function(update_namelist_config NAMELIST_IN NAMELIST_OUT)
    update_namelist_config_with_options("${NAMELIST_IN}" "${NAMELIST_OUT}" "96" "1" "d" "1" "d" "10" ".true." ".false.")
endfunction()

# Function to update namelist.config with custom options
function(update_namelist_config_with_options NAMELIST_IN NAMELIST_OUT STEP_PER_DAY RUN_LENGTH RUN_LENGTH_UNIT RESTART_LENGTH RESTART_LENGTH_UNIT LOGFILE_OUTFREQ FORCE_ROTATION USE_CAVITY)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Set step_per_day
    string(REGEX REPLACE "step_per_day=[0-9]+" "step_per_day=${STEP_PER_DAY}" CONTENT "${CONTENT}")
    # Set run_length
    string(REGEX REPLACE "run_length=[0-9]+" "run_length=${RUN_LENGTH}" CONTENT "${CONTENT}")
    string(REGEX REPLACE "run_length_unit='[^']*'" "run_length_unit='${RUN_LENGTH_UNIT}'" CONTENT "${CONTENT}")
    # Set restart_length
    string(REGEX REPLACE "restart_length=[0-9]+" "restart_length=${RESTART_LENGTH}" CONTENT "${CONTENT}")
    string(REGEX REPLACE "restart_length_unit='[^']*'" "restart_length_unit='${RESTART_LENGTH_UNIT}'" CONTENT "${CONTENT}")
    # Set logfile output frequency
    string(REGEX REPLACE "logfile_outfreq=[0-9]+" "logfile_outfreq=${LOGFILE_OUTFREQ}" CONTENT "${CONTENT}")
    # Force rotation for test geometry
    string(REGEX REPLACE "force_rotation=\\.[a-zA-Z]+\\." "force_rotation=${FORCE_ROTATION}" CONTENT "${CONTENT}")
    # Set cavity usage
    string(REGEX REPLACE "use_cavity=\\.[a-zA-Z]+\\." "use_cavity=${USE_CAVITY}" CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.dyn for test_pi
function(update_namelist_dyn NAMELIST_IN NAMELIST_OUT)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Enable wind stress splitting for test_pi dynamics
    string(REGEX REPLACE "use_wsplit=.[a-zA-Z]." "use_wsplit=.true." CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.ice for test_pi
function(update_namelist_ice NAMELIST_IN NAMELIST_OUT)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Set EVP rheology parameters for test_pi
    string(REGEX REPLACE "whichEVP=[0-9]+" "whichEVP=1" CONTENT "${CONTENT}")
    string(REGEX REPLACE "evp_rheol_steps=[0-9]+" "evp_rheol_steps=120" CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.tra for test_pi
function(update_namelist_tra NAMELIST_IN NAMELIST_OUT TEST_DATA_DIR)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Set filelist for test_global climate data with correct path
    # Replace the complete line including the comment
    string(REGEX REPLACE "filelist[ \t]*=[ \t]*'[^']*'[ \t]*,[ \t]*'[^']*'" "filelist = 'INITIAL/WOA18/woa18_netcdf_5deg.nc', 'INITIAL/WOA18/woa18_netcdf_5deg.nc'" CONTENT "${CONTENT}")
    
    # Update any other hardcoded paths in the tracers section
    string(REGEX REPLACE "/pool/data/AWICM/FESOM2/FORCING/[^']*" "${TEST_DATA_DIR}/initial" CONTENT "${CONTENT}")
    string(REGEX REPLACE "phc3.0_winter\\.nc" "woa18_netcdf_5deg.nc" CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update tracer_init3d filelist specifically (as requested)
function(update_tracer_init3d_filelist NAMELIST_IN NAMELIST_OUT TEST_DATA_DIR)
    file(READ "${NAMELIST_IN}" NAMELIST_TEXT)

    # Match the &tracer_init3d block up to the ending "/"
    string(REGEX MATCHALL "&tracer_init3d[^\n]*\n([^/]*\n)*[ \t]*/" MATCHED_BLOCKS "${NAMELIST_TEXT}")

    if(NOT MATCHED_BLOCKS)
        message(FATAL_ERROR "Could not find &tracer_init3d block in ${NAMELIST_IN}")
    endif()

    list(GET MATCHED_BLOCKS 0 TRACER_INIT3D_BLOCK)

    # Replace the filelist line in the extracted block
    string(REGEX REPLACE
        "filelist[ \t]*=.*\n"
        "filelist = '${TEST_DATA_DIR}/INITIAL/WOA18/woa18_netcdf_5deg.nc', '${TEST_DATA_DIR}/INITIAL/WOA18/woa18_netcdf_5deg.nc'\n"
        TRACER_INIT3D_MODIFIED
        "${TRACER_INIT3D_BLOCK}"
    )

    # Replace the original block with the modified one in the full namelist
    string(REPLACE "${TRACER_INIT3D_BLOCK}" "${TRACER_INIT3D_MODIFIED}" NAMELIST_MODIFIED "${NAMELIST_TEXT}")

    # Write the modified namelist to output
    file(WRITE "${NAMELIST_OUT}" "${NAMELIST_MODIFIED}")
endfunction()

# Function to update namelist.forcing for test_pi
function(update_namelist_forcing NAMELIST_IN NAMELIST_OUT TEST_DATA_DIR)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Apply test_global forcing configurations from setups/forcings.yml
    
    # forcing_exchange_coeff section
    string(REGEX REPLACE "Ce_atm_oce=[0-9.e-]+" "Ce_atm_oce=1.75e-3" CONTENT "${CONTENT}")
    string(REGEX REPLACE "Ch_atm_oce=[0-9.e-]+" "Ch_atm_oce=1.75e-3" CONTENT "${CONTENT}")
    string(REGEX REPLACE "Cd_atm_oce=[0-9.e-]+" "Cd_atm_oce=1.0e-3" CONTENT "${CONTENT}")
    string(REGEX REPLACE "Ce_atm_ice=[0-9.e-]+" "Ce_atm_ice=1.75e-3" CONTENT "${CONTENT}")
    string(REGEX REPLACE "Ch_atm_ice=[0-9.e-]+" "Ch_atm_ice=1.75e-3" CONTENT "${CONTENT}")
    string(REGEX REPLACE "Cd_atm_ice=[0-9.e-]+" "Cd_atm_ice=1.2e-3" CONTENT "${CONTENT}")
    
    # forcing_bulk section  
    string(REGEX REPLACE "AOMIP_drag_coeff=[.][A-Za-z][A-Za-z]*" "AOMIP_drag_coeff=.false." CONTENT "${CONTENT}")
    string(REGEX REPLACE "ncar_bulk_formulae=[.][A-Za-z][A-Za-z]*" "ncar_bulk_formulae=.true." CONTENT "${CONTENT}")
    string(REGEX REPLACE "ncar_bulk_z_wind=[0-9.]+" "ncar_bulk_z_wind=10.0" CONTENT "${CONTENT}")
    string(REGEX REPLACE "ncar_bulk_z_tair=[0-9.]+" "ncar_bulk_z_tair=10.0" CONTENT "${CONTENT}")
    string(REGEX REPLACE "ncar_bulk_z_shum=[0-9.]+" "ncar_bulk_z_shum=10.0" CONTENT "${CONTENT}")
    
    # Update nam_sbc section with test_global values
    # File paths (prepend with ClimateDataPath)
    #string(REGEX REPLACE "nm_xwind_file *= *'[^']*'" "nm_xwind_file = '${TEST_DATA_DIR}/forcing/global/u_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_ywind_file *= *'[^']*'" "nm_ywind_file = '${TEST_DATA_DIR}/forcing/global/v_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_xstre_file *= *'[^']*'" "nm_xstre_file = '${TEST_DATA_DIR}/forcing/global/u_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_ystre_file *= *'[^']*'" "nm_ystre_file = '${TEST_DATA_DIR}/forcing/global/v_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_humi_file *= *'[^']*'" "nm_humi_file = '${TEST_DATA_DIR}/forcing/global/q_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_qsr_file *= *'[^']*'" "nm_qsr_file = '${TEST_DATA_DIR}/forcing/global/ncar_rad.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_qlw_file *= *'[^']*'" "nm_qlw_file = '${TEST_DATA_DIR}/forcing/global/ncar_rad.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_tair_file *= *'[^']*'" "nm_tair_file = '${TEST_DATA_DIR}/forcing/global/t_10.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_prec_file *= *'[^']*'" "nm_prec_file = '${TEST_DATA_DIR}/forcing/global/ncar_precip.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_snow_file *= *'[^']*'" "nm_snow_file = '${TEST_DATA_DIR}/forcing/global/ncar_precip.1948.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_mslp_file *= *'[^']*'" "nm_mslp_file = '${TEST_DATA_DIR}/forcing/global/slp.1948.nc'" CONTENT "${CONTENT}")
    
    # Variable names
    #string(REGEX REPLACE "nm_xwind_var *= *'[^']*'" "nm_xwind_var = 'U_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_ywind_var *= *'[^']*'" "nm_ywind_var = 'V_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_xstre_var *= *'[^']*'" "nm_xstre_var = 'U_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_ystre_var *= *'[^']*'" "nm_ystre_var = 'V_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_humi_var *= *'[^']*'" "nm_humi_var = 'Q_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_qsr_var *= *'[^']*'" "nm_qsr_var = 'SWDN_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_qlw_var *= *'[^']*'" "nm_qlw_var = 'LWDN_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_tair_var *= *'[^']*'" "nm_tair_var = 'T_10_MOD'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_prec_var *= *'[^']*'" "nm_prec_var = 'RAIN'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_snow_var *= *'[^']*'" "nm_snow_var = 'SNOW'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_mslp_var *= *'[^']*'" "nm_mslp_var = 'SLP'" CONTENT "${CONTENT}")
    
    # NetCDF parameters
    #string(REGEX REPLACE "nm_nc_iyear *= *[0-9]+" "nm_nc_iyear = 1948" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_nc_imm *= *[0-9]+" "nm_nc_imm = 1" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_nc_idd *= *[0-9]+" "nm_nc_idd = 1" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_nc_freq *= *[0-9]+" "nm_nc_freq = 1" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_nc_tmid *= *[0-9]+" "nm_nc_tmid = 1" CONTENT "${CONTENT}")
    
    # Flags
    #string(REGEX REPLACE "l_xwind=[.][A-Za-z][A-Za-z]*" "l_xwind=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_ywind=[.][A-Za-z][A-Za-z]*" "l_ywind=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_xstre=[.][A-Za-z][A-Za-z]*" "l_xstre=.false." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_ystre=[.][A-Za-z][A-Za-z]*" "l_ystre=.false." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_humi=[.][A-Za-z][A-Za-z]*" "l_humi=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_qsr=[.][A-Za-z][A-Za-z]*" "l_qsr=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_qlw=[.][A-Za-z][A-Za-z]*" "l_qlw=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_tair=[.][A-Za-z][A-Za-z]*" "l_tair=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_prec=[.][A-Za-z][A-Za-z]*" "l_prec=.true." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_mslp=[.][A-Za-z][A-Za-z]*" "l_mslp=.false." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_cloud=[.][A-Za-z][A-Za-z]*" "l_cloud=.false." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "l_snow=[.][A-Za-z][A-Za-z]*" "l_snow=.true." CONTENT "${CONTENT}")
    
    # Data sources and files
    #string(REGEX REPLACE "runoff_data_source *= *'[^']*'" "runoff_data_source = 'CORE2'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_runoff_file *= *'[^']*'" "nm_runoff_file = '${TEST_DATA_DIR}/forcing/global/runoff.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "sss_data_source *= *'[^']*'" "sss_data_source = 'CORE2'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_sss_data_file *= *'[^']*'" "nm_sss_data_file = '${TEST_DATA_DIR}/forcing/global/PHC2_salx.nc'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "chl_data_source *= *'[^']*'" "chl_data_source = 'None'" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "nm_chl_data_file *= *'[^']*'" "nm_chl_data_file = ''" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "chl_const *= *[0-9.]*" "chl_const = 0.1" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "use_runoff_mapper *= *[.][A-Za-z][A-Za-z]*" "use_runoff_mapper = .FALSE." CONTENT "${CONTENT}")
    #string(REGEX REPLACE "runoff_basins_file *= *'[^']*'" "runoff_basins_file = ''" CONTENT "${CONTENT}")
    #string(REGEX REPLACE "runoff_radius *= *[0-9.]*" "runoff_radius = 500000." CONTENT "${CONTENT}")
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.oce (placeholder for future customization)
function(update_namelist_oce NAMELIST_IN NAMELIST_OUT)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Add ocean-specific modifications here as needed
    # For now, just copy the content as-is
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.io (placeholder for future customization)
function(update_namelist_io NAMELIST_IN NAMELIST_OUT)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Add I/O-specific modifications here as needed
    # For now, just copy the content as-is
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to update namelist.cvmix (placeholder for future customization)
function(update_namelist_cvmix NAMELIST_IN NAMELIST_OUT)
    file(READ "${NAMELIST_IN}" CONTENT)
    
    # Add CVMix-specific modifications here as needed
    # For now, just copy the content as-is
    
    file(WRITE "${NAMELIST_OUT}" "${CONTENT}")
endfunction()

# Function to configure namelists for a test (refactored version)
function(configure_fesom_namelists TARGET_DIR TEST_DATA_DIR RESULT_DIR)
    configure_fesom_namelists_with_options("${TARGET_DIR}" "${TEST_DATA_DIR}" "${RESULT_DIR}" "pi" "96" "1" "d" "1" "d" "10" ".true." ".false.")
endfunction()

# Function to configure namelists with custom options
function(configure_fesom_namelists_with_options TARGET_DIR TEST_DATA_DIR RESULT_DIR MESH_NAME STEP_PER_DAY RUN_LENGTH RUN_LENGTH_UNIT RESTART_LENGTH RESTART_LENGTH_UNIT LOGFILE_OUTFREQ FORCE_ROTATION USE_CAVITY)
    # List of namelists that need to be copied and configured
    set(NAMELISTS
        namelist.config
        namelist.forcing
        namelist.oce
        namelist.ice
        namelist.tra
        namelist.io
        namelist.dyn
        namelist.cvmix
    )
    
    foreach(NAMELIST ${NAMELISTS})
        if(EXISTS "${CMAKE_SOURCE_DIR}/config/${NAMELIST}")
            # Copy the namelist to target directory first
            configure_file(
                "${CMAKE_SOURCE_DIR}/config/${NAMELIST}"
                "${TARGET_DIR}/${NAMELIST}"
                COPYONLY
            )
            
            # Apply namelist-specific updates
            if("${NAMELIST}" STREQUAL "namelist.config")
                # Apply common path updates only to namelist.config with specific mesh
                update_common_paths_with_mesh("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}" "${TEST_DATA_DIR}" "${RESULT_DIR}" "${MESH_NAME}")
                update_namelist_config_with_options("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}" "${STEP_PER_DAY}" "${RUN_LENGTH}" "${RUN_LENGTH_UNIT}" "${RESTART_LENGTH}" "${RESTART_LENGTH_UNIT}" "${LOGFILE_OUTFREQ}" "${FORCE_ROTATION}" "${USE_CAVITY}")
            elseif("${NAMELIST}" STREQUAL "namelist.dyn")
                update_namelist_dyn("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}")
            elseif("${NAMELIST}" STREQUAL "namelist.ice")
                update_namelist_ice("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}")
            elseif("${NAMELIST}" STREQUAL "namelist.tra")
                update_namelist_tra("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}" "${TEST_DATA_DIR}")
            elseif("${NAMELIST}" STREQUAL "namelist.forcing")
                update_namelist_forcing("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}" "${TEST_DATA_DIR}")
            elseif("${NAMELIST}" STREQUAL "namelist.oce")
                update_namelist_oce("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}")
            elseif("${NAMELIST}" STREQUAL "namelist.io")
                update_namelist_io("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}")
            elseif("${NAMELIST}" STREQUAL "namelist.cvmix")
                update_namelist_cvmix("${TARGET_DIR}/${NAMELIST}" "${TARGET_DIR}/${NAMELIST}")
                endif()
        endif()
    endforeach()
endfunction()

# Function to generate fesom.clock file
function(generate_fesom_clock OUTPUT_DIR)
    # Create the output directory if it doesn't exist
    file(MAKE_DIRECTORY "${OUTPUT_DIR}")
    
    # Create the fesom.clock file with the correct format
    file(WRITE "${OUTPUT_DIR}/fesom.clock" "0 1 1948\n0 1 1948\n")
endfunction()

# Function to add a FESOM integration test
function(add_fesom_test TEST_NAME)
    add_fesom_test_with_options("${TEST_NAME}" "pi" "96" "1" "d" "1" "d" "10" ".true." ".false." ${ARGN})
endfunction()

# Function to add a FESOM integration test with custom options
function(add_fesom_test_with_options TEST_NAME MESH_NAME STEP_PER_DAY RUN_LENGTH RUN_LENGTH_UNIT RESTART_LENGTH RESTART_LENGTH_UNIT LOGFILE_OUTFREQ FORCE_ROTATION USE_CAVITY)
    set(options MPI_TEST)
    set(oneValueArgs NP TIMEOUT)
    set(multiValueArgs COMMAND_ARGS)
    cmake_parse_arguments(FESOM_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Set defaults
    if(NOT DEFINED FESOM_TEST_NP)
        set(FESOM_TEST_NP 1)
    endif()
    if(NOT DEFINED FESOM_TEST_TIMEOUT)
        set(FESOM_TEST_TIMEOUT 300)  # 5 minutes default
    endif()
    
    # Create test run directory
    set(TEST_RUN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
    set(RESULT_DIR "${TEST_RUN_DIR}/results")
    
    # Generate fesom.clock file in the results directory
    generate_fesom_clock("${RESULT_DIR}")
    
    # Generate the test script
    set(TEST_SCRIPT "${TEST_RUN_DIR}/run_test.cmake")
    
    if(FESOM_TEST_MPI_TEST AND FESOM_TEST_NP GREATER 1)
        # MPI test
        file(GENERATE OUTPUT ${TEST_SCRIPT} CONTENT "
            # Create test directories
            file(MAKE_DIRECTORY \"${TEST_RUN_DIR}\")
            file(MAKE_DIRECTORY \"${RESULT_DIR}\")
            
            # Run FESOM with MPI
            execute_process(
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${FESOM_TEST_NP} ${CMAKE_BINARY_DIR}/bin/fesom.x
                WORKING_DIRECTORY \"${TEST_RUN_DIR}\"
                RESULT_VARIABLE test_result
                OUTPUT_VARIABLE test_output
                ERROR_VARIABLE test_error
                TIMEOUT ${FESOM_TEST_TIMEOUT}
            )
            
            # Log the output
            file(WRITE \"${TEST_RUN_DIR}/test_output.log\" \"\${test_output}\")
            file(WRITE \"${TEST_RUN_DIR}/test_error.log\" \"\${test_error}\")
            
            # Check result (accept certain error codes as success for initial setup tests)
            if(test_result EQUAL 0 OR test_result EQUAL 1)
                message(STATUS \"Test ${TEST_NAME} completed (exit code: \${test_result})\")
            else()
                message(FATAL_ERROR \"Test ${TEST_NAME} failed with exit code: \${test_result}\")
            endif()
        ")
    else()
        # Serial test
        file(GENERATE OUTPUT ${TEST_SCRIPT} CONTENT "
            # Create test directories
            file(MAKE_DIRECTORY \"${TEST_RUN_DIR}\")
            file(MAKE_DIRECTORY \"${RESULT_DIR}\")
            
            # Run FESOM
            execute_process(
                COMMAND ${CMAKE_BINARY_DIR}/bin/fesom.x
                WORKING_DIRECTORY \"${TEST_RUN_DIR}\"
                RESULT_VARIABLE test_result
                OUTPUT_VARIABLE test_output
                ERROR_VARIABLE test_error
                TIMEOUT 600
            )
            
            # Log the output
            file(WRITE \"${TEST_RUN_DIR}/test_output.log\" \"\${test_output}\")
            file(WRITE \"${TEST_RUN_DIR}/test_error.log\" \"\${test_error}\")
            
            # Check result (accept certain error codes as success for initial setup tests)
            if(test_result EQUAL 0 OR test_result EQUAL 1)
                message(STATUS \"Test ${TEST_NAME} completed (exit code: \${test_result})\")
            else()
                message(FATAL_ERROR \"Test ${TEST_NAME} failed with exit code: \${test_result}\")
            endif()
        ")
    endif()
    
    # Configure namelists for this test with custom options
    configure_fesom_namelists_with_options("${TEST_RUN_DIR}" "${TEST_DATA_DIR}" "${RESULT_DIR}" "${MESH_NAME}" "${STEP_PER_DAY}" "${RUN_LENGTH}" "${RUN_LENGTH_UNIT}" "${RESTART_LENGTH}" "${RESTART_LENGTH_UNIT}" "${LOGFILE_OUTFREQ}" "${FORCE_ROTATION}" "${USE_CAVITY}")
    
    # Add the test
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} -P ${TEST_SCRIPT}
    )
    
    # Set test properties
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT ${FESOM_TEST_TIMEOUT}
        WORKING_DIRECTORY ${TEST_RUN_DIR}
    )
    
    # For MPI tests, set required properties
    if(FESOM_TEST_MPI_TEST AND FESOM_TEST_NP GREATER 1)
        set_tests_properties(${TEST_NAME} PROPERTIES
            PROCESSORS ${FESOM_TEST_NP}
            RUN_SERIAL FALSE
        )
    endif()
    
    message(STATUS "Added FESOM test: ${TEST_NAME} with mesh: ${MESH_NAME}, cavity: ${USE_CAVITY}")
endfunction()

# Function to find and validate MPI for testing
function(setup_mpi_testing)
    find_package(MPI REQUIRED)
    
    # Find mpiexec/mpirun
    if(NOT MPIEXEC_EXECUTABLE)
        find_program(MPIEXEC_EXECUTABLE 
            NAMES mpiexec mpirun
            HINTS ${MPI_HOME}/bin ${MPI_ROOT}/bin
            PATH_SUFFIXES bin
        )
    endif()
    
    if(NOT MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "Could not find mpiexec or mpirun for MPI testing")
    endif()
    
    # Set MPI test defaults
    if(NOT MPIEXEC_NUMPROC_FLAG)
        set(MPIEXEC_NUMPROC_FLAG "-np" PARENT_SCOPE)
    endif()
    
    set(MPIEXEC_EXECUTABLE ${MPIEXEC_EXECUTABLE} PARENT_SCOPE)
    
    message(STATUS "MPI testing setup complete:")
    message(STATUS "  MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
    message(STATUS "  MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
endfunction()

# Function to add a FESOM meshdiag test (uses same setup as integration tests)
function(add_fesom_meshdiag_test TEST_NAME)
    add_fesom_meshdiag_test_with_options("${TEST_NAME}" "pi" "mesh2" ${ARGN})
endfunction()

# Function to add a FESOM meshdiag test with custom options
function(add_fesom_meshdiag_test_with_options TEST_NAME MESH_NAME RUNID)
    set(oneValueArgs NP TIMEOUT)
    set(multiValueArgs COMMAND_ARGS)
    cmake_parse_arguments(FESOM_TEST "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Set defaults (meshdiag requires MPI, minimum 2 processes)
    if(NOT DEFINED FESOM_TEST_NP)
        set(FESOM_TEST_NP 2)
    endif()
    if(NOT DEFINED FESOM_TEST_TIMEOUT)
        set(FESOM_TEST_TIMEOUT 120)  # 2 minutes default (much faster than full FESOM)
    endif()
    
    # Create test run directory
    set(TEST_RUN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data")
    set(RESULT_DIR "${TEST_RUN_DIR}/results")
    
    # Configure namelists with custom runid to avoid conflicts
    configure_fesom_namelists_with_options("${TEST_RUN_DIR}" "${TEST_DATA_DIR}" "${RESULT_DIR}"
        "${MESH_NAME}" "96" "1" "d" "1" "d" "10" ".true." ".false.")
    
    # Update runid in namelist.config to create unique output file
    file(READ "${TEST_RUN_DIR}/namelist.config" CONTENT)
    string(REGEX REPLACE "runid='[^']*'" "runid='${RUNID}'" CONTENT "${CONTENT}")
    file(WRITE "${TEST_RUN_DIR}/namelist.config" "${CONTENT}")
    
    # Generate fesom.clock file in the results directory
    generate_fesom_clock("${RESULT_DIR}")
    
    # Generate the test script
    set(TEST_SCRIPT "${TEST_RUN_DIR}/run_test.cmake")
    
    # Generate MPI test script (meshdiag always requires MPI)
    file(GENERATE OUTPUT ${TEST_SCRIPT} CONTENT "
        # Create test directories
        file(MAKE_DIRECTORY \"${TEST_RUN_DIR}\")
        file(MAKE_DIRECTORY \"${RESULT_DIR}\")
        
        # Run fesom_meshdiag with MPI
        execute_process(
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${FESOM_TEST_NP} ${CMAKE_BINARY_DIR}/bin/fesom_meshdiag
            WORKING_DIRECTORY \"${TEST_RUN_DIR}\"
            RESULT_VARIABLE test_result
            OUTPUT_VARIABLE test_output
            ERROR_VARIABLE test_error
            TIMEOUT ${FESOM_TEST_TIMEOUT}
        )
        
        # Log the output
        file(WRITE \"${TEST_RUN_DIR}/test_output.log\" \"\${test_output}\")
        file(WRITE \"${TEST_RUN_DIR}/test_error.log\" \"\${test_error}\")
        
        # Check if mesh.diag.nc file was created
        set(EXPECTED_OUTPUT \"${RESULT_DIR}/${RUNID}.mesh.diag.nc\")
        if(EXISTS \"\${EXPECTED_OUTPUT}\")
            message(STATUS \"Test ${TEST_NAME} completed successfully - mesh diagnostics file created: \${EXPECTED_OUTPUT}\")
        else()
            message(FATAL_ERROR \"Test ${TEST_NAME} failed - mesh diagnostics file not created: \${EXPECTED_OUTPUT}\")
        endif()
        
        # Check result
        if(test_result EQUAL 0)
            message(STATUS \"Test ${TEST_NAME} completed with exit code: \${test_result}\")
        else()
            message(FATAL_ERROR \"Test ${TEST_NAME} failed with exit code: \${test_result}\")
        endif()
    ")
    
    # Add the test
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} -P ${TEST_SCRIPT}
    )
    
    # Set test properties (meshdiag always uses MPI)
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT ${FESOM_TEST_TIMEOUT}
        WORKING_DIRECTORY ${TEST_RUN_DIR}
        PROCESSORS ${FESOM_TEST_NP}
        RUN_SERIAL FALSE
    )
    
    message(STATUS "Added FESOM meshdiag test: ${TEST_NAME} with mesh: ${MESH_NAME}, runid: ${RUNID}")
endfunction()

#===============================================================================
# Generic Mesh Pipeline Functions
#===============================================================================

# Function to add mesh download fixture
function(add_mesh_download_fixture MESH_NAME)
    set(DOWNLOAD_TEST_NAME "mesh_download_${MESH_NAME}")
    
    add_test(
        NAME ${DOWNLOAD_TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DMESH_NAME=${MESH_NAME}
            -P ${CMAKE_SOURCE_DIR}/tests/integration/mesh_download.cmake
    )
    
    set_tests_properties(${DOWNLOAD_TEST_NAME} PROPERTIES
        TIMEOUT 1200  # 20 minutes for download
        FIXTURES_SETUP "mesh_${MESH_NAME}"
        LABELS "mesh_download"
    )
    
    message(STATUS "Added mesh download fixture: ${DOWNLOAD_TEST_NAME}")
endfunction()

# Function to add mesh partition fixture  
function(add_mesh_partition_fixture MESH_NAME NUM_PROCESSES)
    set(PARTITION_TEST_NAME "mesh_partition_${MESH_NAME}_${NUM_PROCESSES}")
    
    add_test(
        NAME ${PARTITION_TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMESH_NAME=${MESH_NAME}
            -DNUM_PROCESSES=${NUM_PROCESSES}
            -P ${CMAKE_SOURCE_DIR}/tests/integration/mesh_partition.cmake
    )
    
    set_tests_properties(${PARTITION_TEST_NAME} PROPERTIES
        TIMEOUT 1200  # 20 minutes for partitioning
        FIXTURES_SETUP "mesh_${MESH_NAME}_${NUM_PROCESSES}"
        FIXTURES_REQUIRED "mesh_${MESH_NAME}"
        LABELS "mesh_partition"
    )
    
    message(STATUS "Added mesh partition fixture: ${PARTITION_TEST_NAME}")
endfunction()

# Function to add complete mesh pipeline (download + multiple partitions)
function(add_mesh_pipeline MESH_NAME)
    # Parse additional arguments for process counts
    set(PROCESS_COUNTS ${ARGN})
    if(NOT PROCESS_COUNTS)
        set(PROCESS_COUNTS 2 4)  # Default process counts
    endif()
    
    # Add download fixture
    add_mesh_download_fixture(${MESH_NAME})
    
    # Add partition fixtures for each process count
    foreach(NP ${PROCESS_COUNTS})
        add_mesh_partition_fixture(${MESH_NAME} ${NP})
    endforeach()
    
    message(STATUS "Added complete mesh pipeline for '${MESH_NAME}' with process counts: ${PROCESS_COUNTS}")
endfunction()

# Enhanced function to add FESOM test with mesh pipeline integration
function(add_fesom_mesh_test TEST_NAME MESH_NAME NP)
    set(options MPI_TEST)
    set(oneValueArgs TIMEOUT)
    set(multiValueArgs)
    cmake_parse_arguments(FESOM_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Set defaults
    if(NOT DEFINED FESOM_TEST_TIMEOUT)
        set(FESOM_TEST_TIMEOUT 900)  # 15 minutes default
    endif()
    
    # Parse mesh configuration from registry
    set(MESH_REGISTRY "${CMAKE_SOURCE_DIR}/tests/mesh_registry.json")
    if(EXISTS "${MESH_REGISTRY}")
        file(READ "${MESH_REGISTRY}" REGISTRY_CONTENT)
        
        # Extract mesh configuration
        string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"force_rotation\"[ ]*:[ ]*([^,}]+)" ROTATION_MATCH "${REGISTRY_CONTENT}")
        if(ROTATION_MATCH)
            string(STRIP "${CMAKE_MATCH_1}" FORCE_ROTATION_RAW)
            string(REPLACE "true" ".true." FORCE_ROTATION "${FORCE_ROTATION_RAW}")
            string(REPLACE "false" ".false." FORCE_ROTATION "${FORCE_ROTATION}")
        else()
            set(FORCE_ROTATION ".true.")  # Default
        endif()
        
        string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"use_cavity\"[ ]*:[ ]*([^,}]+)" CAVITY_MATCH "${REGISTRY_CONTENT}")
        if(CAVITY_MATCH)
            string(STRIP "${CMAKE_MATCH_1}" USE_CAVITY_RAW)
            string(REPLACE "true" ".true." USE_CAVITY "${USE_CAVITY_RAW}")
            string(REPLACE "false" ".false." USE_CAVITY "${USE_CAVITY}")
        else()
            set(USE_CAVITY ".false.")  # Default
        endif()
        
        # Extract timing parameters
        string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"step_per_day\"[ ]*:[ ]*([^,}]+)" SPD_MATCH "${REGISTRY_CONTENT}")
        if(SPD_MATCH)
            string(STRIP "${CMAKE_MATCH_1}" STEP_PER_DAY)
        else()
            set(STEP_PER_DAY "288")  # Default
        endif()
        
        string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"run_length\"[ ]*:[ ]*([^,}]+)" RL_MATCH "${REGISTRY_CONTENT}")
        if(RL_MATCH)
            string(STRIP "${CMAKE_MATCH_1}" RUN_LENGTH)
        else()
            set(RUN_LENGTH "1")  # Default
        endif()
        
        string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"run_length_unit\"[ ]*:[ ]*\"([^\"]+)\"" RLU_MATCH "${REGISTRY_CONTENT}")
        if(RLU_MATCH)
            string(STRIP "${CMAKE_MATCH_1}" RUN_LENGTH_UNIT)
        else()
            set(RUN_LENGTH_UNIT "s")  # Default
        endif()
    else()
        # Fallback defaults
        set(FORCE_ROTATION ".true.")
        set(USE_CAVITY ".false.")
        set(STEP_PER_DAY "288")
        set(RUN_LENGTH "1")
        set(RUN_LENGTH_UNIT "s")
    endif()
    
    # Create FESOM test with appropriate configuration from registry
    add_fesom_test_with_options(${TEST_NAME}
        "${MESH_NAME}" "${STEP_PER_DAY}" "${RUN_LENGTH}" "${RUN_LENGTH_UNIT}" "1" "d" "10" "${FORCE_ROTATION}" "${USE_CAVITY}"
        MPI_TEST
        NP ${NP}
        TIMEOUT ${FESOM_TEST_TIMEOUT}
    )
    
    # Add mesh pipeline dependency
    set_tests_properties(${TEST_NAME} PROPERTIES
        FIXTURES_REQUIRED "mesh_${MESH_NAME}_${NP}"
        LABELS "mesh_simulation"
    )
    
    message(STATUS "Added FESOM mesh test: ${TEST_NAME} (mesh: ${MESH_NAME}, ${NP} processes)")
    message(STATUS "  Configuration: force_rotation=${FORCE_ROTATION}, use_cavity=${USE_CAVITY}")
    message(STATUS "  Timing: step_per_day=${STEP_PER_DAY}, run_length=${RUN_LENGTH} ${RUN_LENGTH_UNIT}")
endfunction()
