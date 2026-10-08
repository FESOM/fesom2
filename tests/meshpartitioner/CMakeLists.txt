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
            -DSOURCE_DIR=${FESOM_TESTING_ROOT}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMESH_NAME=pi
            -DNUM_PROCESSES=4
            -DMESH_WORKDIR=${CMAKE_BINARY_DIR}/tests/meshpartitioner/meshes/pi
            -P ${FESOM_TESTING_ROOT}/tests/integration/mesh_partition.cmake
    )

    set_tests_properties(meshpartitioner_partition_local_pi_mesh_4 PROPERTIES
        TIMEOUT 300  # 5 minutes
        LABELS "meshpartitioner"
    )

    function(add_pi_partition_fixture PARTITION_NAME NUM_PROCESSES N_PART)
        add_test(
            NAME meshpartitioner_partition_local_pi_${PARTITION_NAME}
            COMMAND ${CMAKE_COMMAND}
                -DSOURCE_DIR=${FESOM_TESTING_ROOT}
                -DBUILD_DIR=${CMAKE_BINARY_DIR}
                -DMESH_NAME=pi
                -DNUM_PROCESSES=${NUM_PROCESSES}
                -DN_PART=${N_PART}
                -DMESH_WORKDIR=${CMAKE_BINARY_DIR}/tests/meshpartitioner/meshes/pi_${PARTITION_NAME}
                -P ${FESOM_TESTING_ROOT}/tests/integration/mesh_partition.cmake
        )
        set_tests_properties(meshpartitioner_partition_local_pi_${PARTITION_NAME} PROPERTIES
            TIMEOUT 300
            LABELS "meshpartitioner"
            FIXTURES_SETUP "partition_pi_${PARTITION_NAME}"
        )
    endfunction()

    add_pi_partition_fixture(flat_1 1 "")
    add_pi_partition_fixture(hier_1x1 1 "1,1")
    add_pi_partition_fixture(hier_1x2 2 "1,2")
    
    # Test partitioning local PI cavity mesh for 4 processes
    add_test(
        NAME meshpartitioner_partition_local_pi_cavity_mesh_4
        COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${FESOM_TESTING_ROOT}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMESH_NAME=pi_cavity
            -DNUM_PROCESSES=4
            -DMESH_WORKDIR=${CMAKE_BINARY_DIR}/tests/meshpartitioner/meshes/pi_cavity
            -P ${FESOM_TESTING_ROOT}/tests/integration/mesh_partition.cmake
    )
    
    set_tests_properties(meshpartitioner_partition_local_pi_cavity_mesh_4 PROPERTIES
        TIMEOUT 300  # 5 minutes
        LABELS "meshpartitioner"
    )
    
    message(STATUS "Added mesh partitioner tests for local meshes")
    
else()
    message(STATUS "Mesh partitioner tests disabled. Set BUILD_MESHPARTITIONER=ON to enable them.")
endif()
