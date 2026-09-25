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
set(SRC_MESH_DIR "${MESH_DIR}/${MESH_NAME}")
set(CONFIG_DIR "${SOURCE_DIR}/config")
set(FESOM_MESHPART_EXECUTABLE "${BUILD_DIR}/bin/fesom_meshpart")

# Where partitioning runs and writes its output (dist_N/ plus the regenerated
# mesh-level files elvls.out/nlvls.out/elvls_raw.out, which the partitioner
# always rewrites in the mesh dir).
#
# - If MESH_WORKDIR is given, partition in that isolated working copy: base mesh
#   files are copied from the read-only SRC_MESH_DIR and all writes stay under
#   MESH_WORKDIR. This keeps the committed source tree clean and prevents one
#   test from mutating another test's mesh inputs. Used by the local tests.
# - Otherwise partition in place at SRC_MESH_DIR (used by the remote tests,
#   which download large meshes into the data dir and partition them there).
if(DEFINED MESH_WORKDIR AND NOT MESH_WORKDIR STREQUAL "")
    set(TARGET_MESH_DIR "${MESH_WORKDIR}")
else()
    set(TARGET_MESH_DIR "${SRC_MESH_DIR}")
endif()

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

    # When using an isolated working copy, seed it from the read-only source mesh
    # (base *.out files only; dist_*/ subdirs are regenerated, not copied).
    if(NOT TARGET_MESH_DIR STREQUAL SRC_MESH_DIR AND NOT EXISTS "${TARGET_MESH_DIR}/nod2d.out")
        if(NOT EXISTS "${SRC_MESH_DIR}/nod2d.out")
            message(FATAL_ERROR "Source mesh '${MESH_NAME}' not found at ${SRC_MESH_DIR}")
        endif()
        message(STATUS "Seeding isolated mesh working copy: ${SRC_MESH_DIR} -> ${TARGET_MESH_DIR}")
        file(MAKE_DIRECTORY "${TARGET_MESH_DIR}")
        file(GLOB SRC_MESH_FILES "${SRC_MESH_DIR}/*.out")
        file(COPY ${SRC_MESH_FILES} DESTINATION "${TARGET_MESH_DIR}")
    endif()

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
    
    # Set mesh path. NOTE: the namelist uses aligned assignments with spaces
    # around '=' (e.g. "MeshPath         = '...'"), so the regexes must tolerate
    # optional whitespace; the leading "([^A-Za-z0-9_])" anchors the key on a
    # non-word boundary (preserved via "\\1") to avoid substring matches.
    string(REGEX REPLACE "([^A-Za-z0-9_])MeshPath[ \t]*=[ \t]*'[^']*'" "\\1MeshPath='${TARGET_MESH_DIR}/'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "([^A-Za-z0-9_])ResultPath[ \t]*=[ \t]*'[^']*'" "\\1ResultPath='${DIST_DIR}/'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")

    # Set partitioning parameters - simple 1-level partitioning
    string(REGEX REPLACE "([^A-Za-z0-9_])n_levels[ \t]*=[ \t]*[0-9]+" "\\1n_levels=1" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "([^A-Za-z0-9_])n_part[ \t]*=[ \t]*[0-9, \t]+" "\\1n_part=${NUM_PROC}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")

    # Apply mesh-specific configuration
    string(REGEX REPLACE "([^A-Za-z0-9_])force_rotation[ \t]*=[ \t]*\\.[a-zA-Z]+\\." "\\1force_rotation=.${MESH_FORCE_ROTATION}." NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "([^A-Za-z0-9_])use_cavity[ \t]*=[ \t]*\\.[a-zA-Z]+\\." "\\1use_cavity=.${MESH_USE_CAVITY}." NAMELIST_CONTENT "${NAMELIST_CONTENT}")

    # Apply mesh-specific timing parameters
    string(REGEX REPLACE "([^A-Za-z0-9_])step_per_day[ \t]*=[ \t]*[0-9]+" "\\1step_per_day=${MESH_STEP_PER_DAY}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "([^A-Za-z0-9_])run_length[ \t]*=[ \t]*[0-9]+" "\\1run_length=${MESH_RUN_LENGTH}" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    string(REGEX REPLACE "([^A-Za-z0-9_])run_length_unit[ \t]*=[ \t]*'[^']*'" "\\1run_length_unit='${MESH_RUN_LENGTH_UNIT}'" NAMELIST_CONTENT "${NAMELIST_CONTENT}")
    
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
    
    # Robust verification: exit code + failure signatures + all expected
    # partition artifacts (rpart.out plus per-rank my_list/com_info files).
    # This catches the case where the partitioner prints an error but exits 0
    # (bare Fortran `stop`) or aborts non-zero, and where rpart.out is written
    # but per-rank distribution files are missing.
    set(REQUIRED_PARTITION_FILES "${DIST_DIR}/rpart.out")
    math(EXPR LAST_RANK "${NUM_PROC} - 1")
    foreach(RANK RANGE 0 ${LAST_RANK})
        # File names are zero-padded to 5 digits (e.g. my_list00000.out).
        string(LENGTH "${RANK}" RANK_LEN)
        math(EXPR PAD "5 - ${RANK_LEN}")
        set(RANK_STR "${RANK}")
        if(PAD GREATER 0)
            foreach(I RANGE 1 ${PAD})
                set(RANK_STR "0${RANK_STR}")
            endforeach()
        endif()
        list(APPEND REQUIRED_PARTITION_FILES
            "${DIST_DIR}/my_list${RANK_STR}.out"
            "${DIST_DIR}/com_info${RANK_STR}.out")
    endforeach()

    include("${SOURCE_DIR}/cmake/CheckFesomRun.cmake")
    check_fesom_run(
        NAME "meshpartitioner_${MESH_NAME}_${NUM_PROC}"
        RESULT "${PARTITION_RESULT}"
        OUTPUT_LOG "${DIST_DIR}/partition_output.log"
        ERROR_LOG "${DIST_DIR}/partition_error.log"
        REQUIRED_ARTIFACTS ${REQUIRED_PARTITION_FILES}
    )

    message(STATUS "✓ ${MESH_NAME} mesh successfully partitioned for ${NUM_PROC} processes")
    message(STATUS "Partition files location: ${DIST_DIR}")
endfunction()

# Execute partitioning if this script is run directly
if(CMAKE_SCRIPT_MODE_FILE)
    partition_mesh("${MESH_NAME}" "${NUM_PROCESSES}")
endif()