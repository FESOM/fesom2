#===============================================================================
# mesh_pipeline_tests.cmake - Mesh pipeline tests using the new system
#===============================================================================

# This file demonstrates how to use the new mesh pipeline system.
# It replaces the manual core2-specific tests with the generic system.

# =============================================================================
# CORE2 Mesh Tests - Using the new pipeline system
# =============================================================================

# Add mesh pipeline for core2 (download + partitioning for 16 processes)
add_mesh_pipeline("core2" 16)

# Add FESOM simulation test that uses the mesh pipeline fixtures
add_fesom_mesh_test(integration_core2_mpi16_pipeline "core2" 16
    MPI_TEST
    TIMEOUT 900
)

# =============================================================================
# DARS Mesh Tests
# =============================================================================

# Add mesh pipeline for DARS (download + partitioning for 16 processes - min 16)
add_mesh_pipeline("dars" 16)

# Add FESOM simulation test for DARS mesh
add_fesom_mesh_test(integration_dars_mpi16 "dars" 16
    MPI_TEST
    TIMEOUT 1200
)

# =============================================================================
# PI Mesh Tests (Remote Version)
# =============================================================================

# Add mesh pipeline for PI remote (download + partitioning for 4, 8 processes - min 4)
add_mesh_pipeline("pi_remote" 4 8)

# Add FESOM simulation tests for PI remote mesh
add_fesom_mesh_test(integration_pi_remote_mpi4 "pi_remote" 4
    MPI_TEST
    TIMEOUT 600
)

add_fesom_mesh_test(integration_pi_remote_mpi8 "pi_remote" 8
    MPI_TEST
    TIMEOUT 600
)

# =============================================================================
# Example: Adding a new mesh (glob) - Just uncomment when ready
# =============================================================================

# # Add mesh pipeline for glob mesh (download + partitioning for 4, 8 processes)
# add_mesh_pipeline("glob" 4 8)
# 
# # Add FESOM simulation test for glob mesh
# add_fesom_mesh_test(integration_glob_mpi4 "glob" 4
#     MPI_TEST
#     TIMEOUT 1200  # Larger mesh may take longer
# )

# =============================================================================
# Example: Adding arctic mesh with cavity support
# =============================================================================

# # Add mesh pipeline for arctic mesh  
# add_mesh_pipeline("arctic" 8 16)
# 
# # Add FESOM simulation test for arctic mesh (with cavity)
# add_fesom_mesh_test(integration_arctic_mpi8 "arctic" 8
#     MPI_TEST
#     TIMEOUT 1500  # Arctic mesh with cavity may take longer
# )