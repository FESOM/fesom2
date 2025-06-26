#===============================================================================
# test_data_check.cmake - Check test data setup for FESOM2
#===============================================================================

# Check that required variables are set
if(NOT DEFINED TEST_DATA_DIR)
    message(FATAL_ERROR "TEST_DATA_DIR not defined")
endif()

message(STATUS "Checking test data setup in: ${TEST_DATA_DIR}")

# Check that test data directory exists
if(NOT EXISTS "${TEST_DATA_DIR}")
    message(FATAL_ERROR "Test data directory does not exist: ${TEST_DATA_DIR}")
endif()

# Check that required subdirectories exist
set(REQUIRED_SUBDIRS
    MESHES
    FORCING 
    INITIAL 
)

foreach(SUBDIR ${REQUIRED_SUBDIRS})
    set(FULL_PATH "${TEST_DATA_DIR}/${SUBDIR}")
    if(NOT EXISTS "${FULL_PATH}")
        message(WARNING "Test data subdirectory missing: ${FULL_PATH}")
        message(STATUS "Creating directory: ${FULL_PATH}")
        file(MAKE_DIRECTORY "${FULL_PATH}")
    else()
        message(STATUS "Found test data directory: ${FULL_PATH}")
    endif()
endforeach()

# Check for any mesh files
file(GLOB MESH_FILES "${TEST_DATA_DIR}/MESHES/*")
if(NOT MESH_FILES)
    message(WARNING "No mesh files found in ${TEST_DATA_DIR}/MESHES/")
    message(STATUS "Note: You may need to add test mesh data for integration tests to work")

else()
    list(LENGTH MESH_FILES NUM_MESH_FILES)
    message(STATUS "Found ${NUM_MESH_FILES} items in meshes directory")
endif()

# Check for any forcing files
file(GLOB FORCING_FILES "${TEST_DATA_DIR}/FORCING/*")
if(NOT FORCING_FILES)
    message(WARNING "No forcing files found in ${TEST_DATA_DIR}/FORCING/")
    message(STATUS "Note: You may need to add minimal forcing data for integration tests")
else()
    list(LENGTH FORCING_FILES NUM_FORCING_FILES)
    message(STATUS "Found ${NUM_FORCING_FILES} items in forcing directory")
endif()

# Check for any initial files
file(GLOB INITIAL_FILES "${TEST_DATA_DIR}/INITIAL/*")
if(NOT INITIAL_FILES)
    message(WARNING "No initial files found in ${TEST_DATA_DIR}/INITIAL/")
    message(STATUS "Note: You may need to add minimal initial data for integration tests")
else()
    list(LENGTH INITIAL_FILES NUM_INITIAL_FILES)
    message(STATUS "Found ${NUM_INITIAL_FILES} items in initial directory")
endif()

message(STATUS "Test data check completed")
