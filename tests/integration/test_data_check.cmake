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
file(GLOB MESH_FILES "${TEST_DATA_DIR}/meshes/*")
if(NOT MESH_FILES)
    message(WARNING "No mesh files found in ${TEST_DATA_DIR}/meshes/")
    message(STATUS "Note: You may need to add test mesh data for integration tests to work")
    
    # Create a placeholder README
    file(WRITE "${TEST_DATA_DIR}/meshes/README.md"
"# Test Mesh Data

Place your test mesh files here. You can:

1. Copy a small mesh from your existing FESOM data
2. Create symbolic links to existing test meshes:
   ```bash
   ln -s /path/to/test/mesh ./test_mesh
   ```
3. Use the toy mesh configurations if available

The integration tests will look for mesh data in this directory.
")
else()
    list(LENGTH MESH_FILES NUM_MESH_FILES)
    message(STATUS "Found ${NUM_MESH_FILES} items in meshes directory")
endif()

# Check for any forcing files
file(GLOB FORCING_FILES "${TEST_DATA_DIR}/forcing/*")
if(NOT FORCING_FILES)
    message(WARNING "No forcing files found in ${TEST_DATA_DIR}/forcing/")
    message(STATUS "Note: You may need to add minimal forcing data for integration tests")
    
    # Create a placeholder README
    file(WRITE "${TEST_DATA_DIR}/forcing/README.md"
"# Test Forcing Data

Place your test forcing/climate files here. You can:

1. Copy minimal forcing data from existing FESOM datasets
2. Create symbolic links to existing test data:
   ```bash
   ln -s /path/to/minimal/forcing ./forcing_data
   ```
3. Use simplified/idealized forcing for testing

The integration tests will configure namelists to use this directory.
")
else()
    list(LENGTH FORCING_FILES NUM_FORCING_FILES)
    message(STATUS "Found ${NUM_FORCING_FILES} items in forcing directory")
endif()

message(STATUS "Test data check completed")
