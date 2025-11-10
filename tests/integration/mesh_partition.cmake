#===============================================================================
# mesh_partition.cmake - Generic mesh partitioning utility for FESOM2
#===============================================================================

cmake_minimum_required(VERSION 3.16)

# Required parameters
if(NOT DEFINED MESH_NAME)
    message(FATAL_ERROR "MESH_NAME parameter is required")
endif()

if(NOT DEFINED NUM_PROCESSES)
    message(FATAL_ERROR "NUM_PROCESSES parameter is required")
endif()

if(NOT DEFINED SOURCE_DIR)
    message(FATAL_ERROR "SOURCE_DIR parameter is required")
endif()

if(NOT DEFINED BUILD_DIR)
    message(FATAL_ERROR "BUILD_DIR parameter is required")
endif()

# Set paths
set(MESH_REGISTRY "${SOURCE_DIR}/tests/mesh_registry.json")
set(TEST_DATA_DIR "${SOURCE_DIR}/tests/data")
set(MESH_DIR "${TEST_DATA_DIR}/MESHES")
set(TARGET_MESH_DIR "${MESH_DIR}/${MESH_NAME}")
set(CONFIG_DIR "${SOURCE_DIR}/config")
set(FESOM_MESHPART_EXECUTABLE "${BUILD_DIR}/bin/fesom_meshpart")

# Function to parse mesh configuration from registry
function(parse_mesh_config MESH_NAME)
    if(NOT EXISTS "${MESH_REGISTRY}")
        message(FATAL_ERROR "Mesh registry not found: ${MESH_REGISTRY}")
    endif()
    
    file(READ "${MESH_REGISTRY}" REGISTRY_CONTENT)
    
    # Extract force_rotation setting
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"force_rotation\"[ ]*:[ ]*([^,}]+)" ROTATION_MATCH "${REGISTRY_CONTENT}")
    if(ROTATION_MATCH)
        string(STRIP "${CMAKE_MATCH_1}" FORCE_ROTATION)
        set(MESH_FORCE_ROTATION "${FORCE_ROTATION}" PARENT_SCOPE)
    else()
        set(MESH_FORCE_ROTATION "true" PARENT_SCOPE)  # Default
    endif()
    
    # Extract use_cavity setting
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"use_cavity\"[ ]*:[ ]*([^,}]+)" CAVITY_MATCH "${REGISTRY_CONTENT}")
    if(CAVITY_MATCH)
        string(STRIP "${CMAKE_MATCH_1}" USE_CAVITY)
        set(MESH_USE_CAVITY "${USE_CAVITY}" PARENT_SCOPE)
    else()
        set(MESH_USE_CAVITY "false" PARENT_SCOPE)  # Default
    endif()
    
    # Extract timing parameters
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"step_per_day\"[ ]*:[ ]*([^,}]+)" SPD_MATCH "${REGISTRY_CONTENT}")
    if(SPD_MATCH)
        string(STRIP "${CMAKE_MATCH_1}" STEP_PER_DAY)
        set(MESH_STEP_PER_DAY "${STEP_PER_DAY}" PARENT_SCOPE)
    else()
        set(MESH_STEP_PER_DAY "288" PARENT_SCOPE)  # Default
    endif()
    
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"run_length\"[ ]*:[ ]*([^,}]+)" RL_MATCH "${REGISTRY_CONTENT}")
    if(RL_MATCH)
        string(STRIP "${CMAKE_MATCH_1}" RUN_LENGTH)
        set(MESH_RUN_LENGTH "${RUN_LENGTH}" PARENT_SCOPE)
    else()
        set(MESH_RUN_LENGTH "1" PARENT_SCOPE)  # Default
    endif()
    
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"run_length_unit\"[ ]*:[ ]*\"([^\"]+)\"" RLU_MATCH "${REGISTRY_CONTENT}")
    if(RLU_MATCH)
        string(STRIP "${CMAKE_MATCH_1}" RUN_LENGTH_UNIT)
        set(MESH_RUN_LENGTH_UNIT "${RUN_LENGTH_UNIT}" PARENT_SCOPE)
    else()
        set(MESH_RUN_LENGTH_UNIT "s" PARENT_SCOPE)  # Default
    endif()
endfunction()

# Function to partition mesh for specified number of processes
function(partition_mesh MESH_NAME NUM_PROC)
    message(STATUS "=== FESOM2 Mesh Partitioning: ${MESH_NAME} (${NUM_PROC} processes) ===")
    
    # Check if mesh exists
    if(NOT EXISTS "${TARGET_MESH_DIR}/nod2d.out")
        message(FATAL_ERROR "Mesh '${MESH_NAME}' not found at ${TARGET_MESH_DIR}. Run mesh download first.")
    endif()
    
    # Check if mesh partitioner exists
    if(NOT EXISTS "${FESOM_MESHPART_EXECUTABLE}")
        message(FATAL_ERROR "Mesh partitioner executable not found at ${FESOM_MESHPART_EXECUTABLE}")
    endif()
    
    # Create dist_N directory
    set(DIST_DIR "${TARGET_MESH_DIR}/dist_${NUM_PROC}")
    
    # Check if partition already exists and is complete
    if(EXISTS "${DIST_DIR}/rpart.out")
        message(STATUS "Partition for ${NUM_PROC} processes already exists")
        message(STATUS "Skipping partitioning")
        return()
    endif()
    
    file(MAKE_DIRECTORY "${DIST_DIR}")
    
    # Parse mesh configuration
    parse_mesh_config("${MESH_NAME}")
    
    message(STATUS "Mesh configuration:")
    message(STATUS "  - force_rotation: ${MESH_FORCE_ROTATION}")
    message(STATUS "  - use_cavity: ${MESH_USE_CAVITY}")
    message(STATUS "  - step_per_day: ${MESH_STEP_PER_DAY}")
    message(STATUS "  - run_length: ${MESH_RUN_LENGTH} ${MESH_RUN_LENGTH_UNIT}")
    
    # Create temporary namelist.config for partitioning
    set(TEMP_NAMELIST "${DIST_DIR}/namelist.config")
    if(NOT EXISTS "${CONFIG_DIR}/namelist.config")
        message(FATAL_ERROR "Base namelist.config not found at ${CONFIG_DIR}/namelist.config")
    endif()
    
    file(COPY "${CONFIG_DIR}/namelist.config" DESTINATION "${DIST_DIR}")
    
    # Update the namelist for this mesh and partition count
    file(READ "${TEMP_NAMELIST}" NAMELIST_CONTENT)
    
    # Set mesh path
    string(REGEX REPLACE "MeshPath='[^']*'" "MeshPath='${TARGET_MESH_DIR}/'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "ResultPath='[^']*'" "ResultPath='${DIST_DIR}/'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    
    # Set partitioning parameters - simple 1-level partitioning
    string(REGEX REPLACE "n_levels=[0-9]+" "n_levels=1" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "n_part=[ 0-9,]+" "n_part=${NUM_PROC}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    
    # Apply mesh-specific configuration
    string(REGEX REPLACE "force_rotation=\\.[a-zA-Z]+\\." "force_rotation=.${MESH_FORCE_ROTATION}." NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "use_cavity=\\.[a-zA-Z]+\\." "use_cavity=.${MESH_USE_CAVITY}." NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    
    # Apply mesh-specific timing parameters
    string(REGEX REPLACE "step_per_day=[0-9]+" "step_per_day=${MESH_STEP_PER_DAY}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "run_length=[0-9]+" "run_length=${MESH_RUN_LENGTH}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "run_length_unit='[^']*'" "run_length_unit='${MESH_RUN_LENGTH_UNIT}'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    
    file(WRITE "${TEMP_NAMELIST}" "${NAMELIST_CONTENT}")
    
    # Run the mesh partitioner
    message(STATUS "Running mesh partitioner...")
    message(STATUS "Command: ${FESOM_MESHPART_EXECUTABLE}")
    message(STATUS "Working directory: ${DIST_DIR}")
    
    execute_process(
        COMMAND "${FESOM_MESHPART_EXECUTABLE}"
        WORKING_DIRECTORY "${DIST_DIR}"
        RESULT_VARIABLE PARTITION_RESULT
        OUTPUT_VARIABLE PARTITION_OUTPUT
        ERROR_VARIABLE PARTITION_ERROR
        TIMEOUT 1200  # 20 minutes timeout
    )
    
    # Log the output
    file(WRITE "${DIST_DIR}/partition_output.log" "${PARTITION_OUTPUT}")
    file(WRITE "${DIST_DIR}/partition_error.log" "${PARTITION_ERROR}")
    
    if(NOT PARTITION_RESULT EQUAL 0)
        message(STATUS "Partition output: ${PARTITION_OUTPUT}")
        message(STATUS "Partition error: ${PARTITION_ERROR}")
        message(FATAL_ERROR "Mesh partitioning failed with exit code: ${PARTITION_RESULT}")
    endif()
    
    # Verify that partitioning files were created
    set(REQUIRED_PARTITION_FILES "rpart.out")
    foreach(REQUIRED_FILE ${REQUIRED_PARTITION_FILES})
        if(NOT EXISTS "${DIST_DIR}/${REQUIRED_FILE}")
            message(FATAL_ERROR "Partitioning failed - ${REQUIRED_FILE} not created")
        endif()
    endforeach()
    
    message(STATUS "✓ ${MESH_NAME} mesh successfully partitioned for ${NUM_PROC} processes")
    message(STATUS "Partition files location: ${DIST_DIR}")
endfunction()

# Execute partitioning if this script is run directly
if(CMAKE_SCRIPT_MODE_FILE)
    partition_mesh("${MESH_NAME}" "${NUM_PROCESSES}")
endif()