#===============================================================================
# tests/meshpartitioner/CMakeLists.txt - Mesh partitioner tests for FESOM2
#===============================================================================

# Mesh partitioner tests (only when BUILD_MESHPARTITIONER is enabled)
if(BUILD_MESHPARTITIONER)
    message(STATUS "Adding mesh partitioner tests...")
    
    # Test partitioning local PI mesh for 4 processes (not currently present)
    add_test(
        NAME meshpartitioner_partition_local_pi_mesh_4
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMESH_NAME=pi
            -DNUM_PROCESSES=4
            -P ${CMAKE_SOURCE_DIR}/tests/integration/mesh_partition.cmake
    )
    
    set_tests_properties(meshpartitioner_partition_local_pi_mesh_4 PROPERTIES
        TIMEOUT 300  # 5 minutes
        LABELS "meshpartitioner"
    )
    
    # Test partitioning local PI cavity mesh for 4 processes
    add_test(
        NAME meshpartitioner_partition_local_pi_cavity_mesh_4
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMESH_NAME=pi_cavity
            -DNUM_PROCESSES=4
            -P ${CMAKE_SOURCE_DIR}/tests/integration/mesh_partition.cmake
    )
    
    set_tests_properties(meshpartitioner_partition_local_pi_cavity_mesh_4 PROPERTIES
        TIMEOUT 300  # 5 minutes
        LABELS "meshpartitioner"
    )
    
    message(STATUS "Added mesh partitioner tests for local meshes")
    
else()
    message(STATUS "Mesh partitioner tests disabled. Set BUILD_MESHPARTITIONER=ON to enable them.")
endif()