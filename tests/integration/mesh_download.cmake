#===============================================================================
# mesh_download.cmake - Generic mesh download utility for FESOM2
#===============================================================================

cmake_minimum_required(VERSION 3.16)

# Required parameters
if(NOT DEFINED MESH_NAME)
    message(FATAL_ERROR "MESH_NAME parameter is required")
endif()

if(NOT DEFINED SOURCE_DIR)
    message(FATAL_ERROR "SOURCE_DIR parameter is required")
endif()

# Set paths
set(MESH_REGISTRY "${SOURCE_DIR}/tests/mesh_registry.json")
set(TEST_DATA_DIR "${SOURCE_DIR}/tests/data")
set(MESH_DIR "${TEST_DATA_DIR}/MESHES")
set(TARGET_MESH_DIR "${MESH_DIR}/${MESH_NAME}")

# Function to parse JSON (simplified - could use more robust parser)
function(parse_mesh_registry MESH_NAME)
    if(NOT EXISTS "${MESH_REGISTRY}")
        message(FATAL_ERROR "Mesh registry not found: ${MESH_REGISTRY}")
    endif()
    
    file(READ "${MESH_REGISTRY}" REGISTRY_CONTENT)
    
    # Simple JSON parsing - look for mesh entry
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"url\"[ ]*:[ ]*\"([^\"]+)\"" URL_MATCH "${REGISTRY_CONTENT}")
    if(URL_MATCH)
        set(MESH_URL "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Mesh '${MESH_NAME}' not found in registry")
    endif()
    
    # Extract archive type
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"archive_type\"[ ]*:[ ]*\"([^\"]+)\"" ARCHIVE_MATCH "${REGISTRY_CONTENT}")
    if(ARCHIVE_MATCH)
        set(ARCHIVE_TYPE "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        set(ARCHIVE_TYPE "tar.gz" PARENT_SCOPE)
    endif()
    
    # Extract verification files
    string(REGEX MATCH "\"${MESH_NAME}\"[^}]*\"required_files\"[ ]*:[ ]*\\[([^\\]]+)\\]" FILES_MATCH "${REGISTRY_CONTENT}")
    if(FILES_MATCH)
        set(VERIFICATION_FILES "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        set(VERIFICATION_FILES "\"nod2d.out\", \"elem2d.out\"" PARENT_SCOPE)
    endif()
endfunction()

# Function to download and extract mesh
function(download_mesh MESH_NAME)
    message(STATUS "=== FESOM2 Mesh Download: ${MESH_NAME} ===")
    
    # Parse mesh configuration
    parse_mesh_registry("${MESH_NAME}")
    
    message(STATUS "Mesh URL: ${MESH_URL}")
    message(STATUS "Archive type: ${ARCHIVE_TYPE}")
    message(STATUS "Target directory: ${TARGET_MESH_DIR}")
    
    # Check if mesh already exists and is complete
    set(MESH_EXISTS TRUE)
    string(REPLACE "\"" "" CLEAN_FILES "${VERIFICATION_FILES}")
    string(REPLACE " " "" CLEAN_FILES "${CLEAN_FILES}")
    string(REPLACE "," ";" FILE_LIST "${CLEAN_FILES}")
    
    foreach(REQUIRED_FILE ${FILE_LIST})
        if(NOT EXISTS "${TARGET_MESH_DIR}/${REQUIRED_FILE}")
            set(MESH_EXISTS FALSE)
            break()
        endif()
    endforeach()
    
    if(MESH_EXISTS)
        message(STATUS "Mesh '${MESH_NAME}' already exists and is complete")
        message(STATUS "Skipping download")
        return()
    endif()
    
    # Create mesh directory
    file(MAKE_DIRECTORY "${TARGET_MESH_DIR}")
    
    # Download mesh
    message(STATUS "Downloading ${MESH_NAME} mesh...")
    message(STATUS "This may take a few minutes...")
    
    set(ARCHIVE_FILE "${TARGET_MESH_DIR}/${MESH_NAME}.${ARCHIVE_TYPE}")
    file(DOWNLOAD "${MESH_URL}" "${ARCHIVE_FILE}"
        SHOW_PROGRESS
        STATUS download_status
        TIMEOUT 1200  # 20 minutes timeout
    )
    
    # Check download status
    list(GET download_status 0 status_code)
    if(NOT status_code EQUAL 0)
        list(GET download_status 1 error_message)
        message(FATAL_ERROR "Failed to download ${MESH_NAME} mesh: ${error_message}")
    endif()
    
    message(STATUS "Extracting ${MESH_NAME} mesh...")
    
    # Extract based on archive type
    if(ARCHIVE_TYPE STREQUAL "tar.gz")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf "${ARCHIVE_FILE}"
            WORKING_DIRECTORY "${TARGET_MESH_DIR}"
            RESULT_VARIABLE extract_result
            OUTPUT_VARIABLE extract_output
            ERROR_VARIABLE extract_error
        )
    else()
        message(FATAL_ERROR "Unsupported archive type: ${ARCHIVE_TYPE}")
    endif()
    
    if(NOT extract_result EQUAL 0)
        message(STATUS "Extract output: ${extract_output}")
        message(STATUS "Extract error: ${extract_error}")
        message(FATAL_ERROR "Failed to extract ${MESH_NAME} mesh")
    endif()
    
    # Clean up archive file
    file(REMOVE "${ARCHIVE_FILE}")
    
    # Handle nested directory structure (most archives contain subdirectory)
    # Check for common subdirectory patterns: core2/, dars/, pi/, etc.
    set(POSSIBLE_SUBDIRS "core2" "dars" "pi" "glob" "arctic" "${MESH_NAME}")
    set(FOUND_SUBDIR "")
    
    foreach(SUBDIR ${POSSIBLE_SUBDIRS})
        if(EXISTS "${TARGET_MESH_DIR}/${SUBDIR}" AND IS_DIRECTORY "${TARGET_MESH_DIR}/${SUBDIR}")
            # Check if this subdirectory contains mesh files
            if(EXISTS "${TARGET_MESH_DIR}/${SUBDIR}/nod2d.out" OR EXISTS "${TARGET_MESH_DIR}/${SUBDIR}/elem2d.out")
                set(FOUND_SUBDIR "${SUBDIR}")
                break()
            endif()
        endif()
    endforeach()
    
    # If we found a subdirectory with mesh files, move them up
    if(NOT FOUND_SUBDIR STREQUAL "")
        message(STATUS "Found mesh files in subdirectory: ${FOUND_SUBDIR}")
        message(STATUS "Moving files to target directory...")
        
        # Get list of files in subdirectory
        file(GLOB MESH_FILES "${TARGET_MESH_DIR}/${FOUND_SUBDIR}/*")
        
        # Move each file
        foreach(MESH_FILE ${MESH_FILES})
            get_filename_component(FILENAME "${MESH_FILE}" NAME)
            file(RENAME "${MESH_FILE}" "${TARGET_MESH_DIR}/${FILENAME}")
        endforeach()
        
        # Remove empty subdirectory
        file(REMOVE_RECURSE "${TARGET_MESH_DIR}/${FOUND_SUBDIR}")
        message(STATUS "Files moved successfully")
    endif()
    
    # Verify extraction
    set(VERIFICATION_FAILED FALSE)
    foreach(REQUIRED_FILE ${FILE_LIST})
        if(NOT EXISTS "${TARGET_MESH_DIR}/${REQUIRED_FILE}")
            message(WARNING "Required file missing after extraction: ${REQUIRED_FILE}")
            set(VERIFICATION_FAILED TRUE)
        endif()
    endforeach()
    
    if(VERIFICATION_FAILED)
        message(FATAL_ERROR "Mesh extraction verification failed for ${MESH_NAME}")
    endif()
    
    message(STATUS "✓ ${MESH_NAME} mesh successfully downloaded and extracted")
    message(STATUS "Location: ${TARGET_MESH_DIR}")
endfunction()

# Execute download if this script is run directly
if(CMAKE_SCRIPT_MODE_FILE)
    download_mesh("${MESH_NAME}")
endif()