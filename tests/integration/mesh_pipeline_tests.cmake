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

# Add meshdiag test that uses the mesh pipeline fixtures
add_fesom_meshdiag_test_with_options(meshdiag_core2_mpi16 "core2" "core2"
    NP 16
    TIMEOUT 300
)

# Make meshdiag test depend on mesh pipeline
set_tests_properties(meshdiag_core2_mpi16 PROPERTIES
    FIXTURES_REQUIRED "mesh_core2_16"
    LABELS "mesh_validation"
)

# =============================================================================
# DARS Mesh Tests
# =============================================================================

# Add mesh pipeline for DARS (download + partitioning for 128 processes for large mesh)
add_mesh_pipeline("dars" 128)

# Add meshdiag test for DARS mesh (large mesh ~3M nodes, needs longer timeout)
add_fesom_meshdiag_test_with_options(meshdiag_dars_mpi128 "dars" "dars"
    NP 128
    TIMEOUT 1800  # 30 minutes for large DARS mesh
)

# Make meshdiag test depend on mesh pipeline
set_tests_properties(meshdiag_dars_mpi128 PROPERTIES
    FIXTURES_REQUIRED "mesh_dars_128"
    LABELS "mesh_validation"
)

# =============================================================================
# PI Mesh Tests (Remote Version)
# =============================================================================

# Add mesh pipeline for PI remote (using recommended processor count from registry)
add_mesh_pipeline("pi_remote" 8)

# Add meshdiag test for PI remote mesh
add_fesom_meshdiag_test_with_options(meshdiag_pi_remote_mpi8 "pi_remote" "pi_remote"
    NP 8
    TIMEOUT 300
)

# Make meshdiag test depend on mesh pipeline
set_tests_properties(meshdiag_pi_remote_mpi8 PROPERTIES
    FIXTURES_REQUIRED "mesh_pi_remote_8"
    LABELS "mesh_validation"
)

# =============================================================================
# NG5 Mesh Tests (High Resolution ~5km)
# =============================================================================

# Add mesh pipeline for ng5 (128 processes for very large mesh ~7.4M nodes)
add_mesh_pipeline("ng5" 128)

# Add meshdiag test for ng5 mesh (very large mesh ~7.4M nodes, needs long timeout)
add_fesom_meshdiag_test_with_options(meshdiag_ng5_mpi128 "ng5" "ng5"
    NP 128
    TIMEOUT 3600  # 60 minutes for very large NG5 mesh
)

# Make meshdiag test depend on mesh pipeline
set_tests_properties(meshdiag_ng5_mpi128 PROPERTIES
    FIXTURES_REQUIRED "mesh_ng5_128"
    LABELS "mesh_validation"
)

# =============================================================================
# ORCA25 Mesh Tests (Global ~0.25-degree)
# =============================================================================

# Add mesh pipeline for orca25 (128 processes for large global mesh)
add_mesh_pipeline("orca25" 128)

# Add meshdiag test for orca25 mesh (large global mesh, needs longer timeout)
add_fesom_meshdiag_test_with_options(meshdiag_orca25_mpi128 "orca25" "orca25"
    NP 128
    TIMEOUT 1800  # 30 minutes for large ORCA25 mesh
)

# Make meshdiag test depend on mesh pipeline
set_tests_properties(meshdiag_orca25_mpi128 PROPERTIES
    FIXTURES_REQUIRED "mesh_orca25_128"
    LABELS "mesh_validation"
)