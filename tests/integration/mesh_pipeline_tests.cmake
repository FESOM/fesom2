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

# Add mesh pipeline for PI remote (using recommended processor count from registry)
add_mesh_pipeline("pi_remote" 8)

# Add FESOM simulation test for PI remote mesh
add_fesom_mesh_test(integration_pi_remote_mpi8 "pi_remote" 8
    MPI_TEST
    TIMEOUT 600
)

# =============================================================================
# NG5 Mesh Tests (High Resolution ~5km)
# =============================================================================

# Add mesh pipeline for ng5 (using recommended processor count from registry)
add_mesh_pipeline("ng5" 16)

# Add FESOM simulation test for ng5 mesh
add_fesom_mesh_test(integration_ng5_mpi16 "ng5" 16
    MPI_TEST
    TIMEOUT 1800  # High resolution mesh may take longer
)

# =============================================================================
# ORCA25 Mesh Tests (Global ~0.25-degree)
# =============================================================================

# Add mesh pipeline for orca25 (using recommended processor count from registry)
add_mesh_pipeline("orca25" 16)

# Add FESOM simulation test for orca25 mesh
add_fesom_mesh_test(integration_orca25_mpi16 "orca25" 16
    MPI_TEST
    TIMEOUT 1500  # Larger global mesh may take longer
)