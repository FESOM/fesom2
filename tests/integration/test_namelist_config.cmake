#===============================================================================
# test_namelist_config.cmake - Test namelist configuration for FESOM2
#===============================================================================

# Check that required variables are set
if(NOT DEFINED TEST_DATA_DIR)
    message(FATAL_ERROR "TEST_DATA_DIR not defined")
endif()

if(NOT DEFINED CONFIG_DIR)
    message(FATAL_ERROR "CONFIG_DIR not defined")
endif()

# Check that config directory exists
if(NOT EXISTS "${CONFIG_DIR}")
    message(FATAL_ERROR "Config directory does not exist: ${CONFIG_DIR}")
endif()

# Check that main namelists exist
set(REQUIRED_NAMELISTS
    namelist.config
    namelist.forcing
    namelist.oce
)

foreach(NAMELIST ${REQUIRED_NAMELISTS})
    if(NOT EXISTS "${CONFIG_DIR}/${NAMELIST}")
        message(FATAL_ERROR "Required namelist not found: ${CONFIG_DIR}/${NAMELIST}")
    endif()
endforeach()

# Test path replacement functionality
set(TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/namelist_test")
file(MAKE_DIRECTORY "${TEST_DIR}")

# Copy a test namelist
file(READ "${CONFIG_DIR}/namelist.config" ORIGINAL_CONTENT)
file(WRITE "${TEST_DIR}/namelist.config" "${ORIGINAL_CONTENT}")

# Test path replacements
set(TEST_MESH_PATH "${TEST_DATA_DIR}/meshes/")
set(TEST_FORCING_PATH "${TEST_DATA_DIR}/forcing/")
set(TEST_RESULT_PATH "${TEST_DIR}/results/")

# Apply the same replacements as in FesomTesting.cmake
file(READ "${TEST_DIR}/namelist.config" CONTENT)
string(REGEX REPLACE "MeshPath='[^']*'" "MeshPath='${TEST_MESH_PATH}'" CONTENT "${CONTENT}")
string(REGEX REPLACE "ClimateDataPath='[^']*'" "ClimateDataPath='${TEST_FORCING_PATH}'" CONTENT "${CONTENT}")
string(REGEX REPLACE "ResultPath='[^']*'" "ResultPath='${TEST_RESULT_PATH}'" CONTENT "${CONTENT}")
file(WRITE "${TEST_DIR}/namelist.config" "${CONTENT}")

# Verify the replacements worked
file(READ "${TEST_DIR}/namelist.config" MODIFIED_CONTENT)

# Check that paths were actually changed
if(NOT MODIFIED_CONTENT MATCHES "MeshPath='${TEST_MESH_PATH}'")
    message(FATAL_ERROR "MeshPath replacement failed")
endif()

if(NOT MODIFIED_CONTENT MATCHES "ClimateDataPath='${TEST_FORCING_PATH}'")
    message(FATAL_ERROR "ClimateDataPath replacement failed")
endif()

if(NOT MODIFIED_CONTENT MATCHES "ResultPath='${TEST_RESULT_PATH}'")
    message(FATAL_ERROR "ResultPath replacement failed")
endif()

message(STATUS "Namelist configuration test passed")

# Clean up
file(REMOVE_RECURSE "${TEST_DIR}")
