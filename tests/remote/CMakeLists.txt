#===============================================================================
# tests/remote/CMakeLists.txt - Remote mesh pipeline tests for FESOM2
#===============================================================================

# Check for download capability (wget, curl, or similar)
find_program(WGET_EXECUTABLE wget)
find_program(CURL_EXECUTABLE curl)
set(HAS_DOWNLOAD_CAPABILITY FALSE)
if(WGET_EXECUTABLE OR CURL_EXECUTABLE)
    set(HAS_DOWNLOAD_CAPABILITY TRUE)
endif()

# Report component availability
message(STATUS "Remote mesh pipeline component availability:")
message(STATUS "  Download capability: ${HAS_DOWNLOAD_CAPABILITY}")
if(HAS_DOWNLOAD_CAPABILITY)
    if(WGET_EXECUTABLE)
        message(STATUS "    - wget: ${WGET_EXECUTABLE}")
    endif()
    if(CURL_EXECUTABLE)
        message(STATUS "    - curl: ${CURL_EXECUTABLE}")
    endif()
else()
    message(STATUS "    - No download tools found (wget or curl required)")
endif()

message(STATUS "  Mesh partitioner: ${BUILD_MESHPARTITIONER}")
if(NOT BUILD_MESHPARTITIONER)
    message(STATUS "    - Enable with -DBUILD_MESHPARTITIONER=ON")
endif()

message(STATUS "  Mesh diagnostics: ${BUILD_MESHDIAG}")
if(NOT BUILD_MESHDIAG)
    message(STATUS "    - Enable with -DBUILD_MESHDIAG=ON")
endif()

# =============================================================================
# CORE2 Remote Mesh Tests
# =============================================================================

if(HAS_DOWNLOAD_CAPABILITY)
    message(STATUS "Adding CORE2 remote mesh download test...")
    
    # Add mesh pipeline for core2 (download + partitioning for 16 processes)
    add_mesh_pipeline("core2" PREFIX remote PROCESS_COUNTS 16)
    
    if(BUILD_MESHPARTITIONER)
        message(STATUS "Adding CORE2 mesh partitioning test...")
        # Partitioning is included in add_mesh_pipeline above
    else()
        message(STATUS "CORE2 mesh partitioning test skipped - BUILD_MESHPARTITIONER=OFF")
    endif()
    
    if(BUILD_MESHDIAG)
        message(STATUS "Adding CORE2 mesh diagnostics test...")
        
        # Add meshdiag test that generates mesh diagnostics file
        add_fesom_meshdiag_test_with_options(generate_meshdiag_core2_mpi16 "core2" "fesom"
            NP 16
            TIMEOUT 300
            PREFIX remote
        )
        
        # Make meshdiag test depend on mesh pipeline
        set_tests_properties(remote_generate_meshdiag_core2_mpi16 PROPERTIES
            FIXTURES_REQUIRED "mesh_core2_16"
            LABELS "generate_meshdiag"
        )
    else()
        message(STATUS "CORE2 mesh diagnostics test skipped - BUILD_MESHDIAG=OFF")
    endif()
    
else()
    message(STATUS "CORE2 remote mesh tests skipped - no download capability (wget or curl required)")
endif()

# =============================================================================
# DARS Remote Mesh Tests
# =============================================================================

if(HAS_DOWNLOAD_CAPABILITY)
    message(STATUS "Adding DARS remote mesh download test...")
    
    # Add mesh pipeline for DARS (download + partitioning for 128 processes for large mesh)
    add_mesh_pipeline("dars" PREFIX remote PROCESS_COUNTS 128)
    
    if(BUILD_MESHPARTITIONER)
        message(STATUS "Adding DARS mesh partitioning test...")
        # Partitioning is included in add_mesh_pipeline above
    else()
        message(STATUS "DARS mesh partitioning test skipped - BUILD_MESHPARTITIONER=OFF")
    endif()
    
    if(BUILD_MESHDIAG)
        message(STATUS "Adding DARS mesh diagnostics test...")
        
        # Add meshdiag test that generates mesh diagnostics file (large mesh ~3M nodes, needs longer timeout)
        add_fesom_meshdiag_test_with_options(generate_meshdiag_dars_mpi128 "dars" "fesom"
            NP 128
            TIMEOUT 1800  # 30 minutes for large DARS mesh
            PREFIX remote
        )
        
        # Make meshdiag test depend on mesh pipeline
        set_tests_properties(remote_generate_meshdiag_dars_mpi128 PROPERTIES
            FIXTURES_REQUIRED "mesh_dars_128"
            LABELS "generate_meshdiag"
        )
    else()
        message(STATUS "DARS mesh diagnostics test skipped - BUILD_MESHDIAG=OFF")
    endif()
    
else()
    message(STATUS "DARS remote mesh tests skipped - no download capability (wget or curl required)")
endif()

# =============================================================================
# PI Remote Mesh Tests
# =============================================================================

if(HAS_DOWNLOAD_CAPABILITY)
    message(STATUS "Adding PI remote mesh download test...")
    
    # Add mesh pipeline for PI remote (using recommended processor count from registry)
    add_mesh_pipeline("pi_remote" PREFIX remote PROCESS_COUNTS 8)
    
    if(BUILD_MESHPARTITIONER)
        message(STATUS "Adding PI remote mesh partitioning test...")
        # Partitioning is included in add_mesh_pipeline above
    else()
        message(STATUS "PI remote mesh partitioning test skipped - BUILD_MESHPARTITIONER=OFF")
    endif()
    
    if(BUILD_MESHDIAG)
        message(STATUS "Adding PI remote mesh diagnostics test...")
        
        # Add meshdiag test that generates mesh diagnostics file
        add_fesom_meshdiag_test_with_options(generate_meshdiag_pi_remote_mpi8 "pi_remote" "fesom"
            NP 8
            TIMEOUT 300
            PREFIX remote
        )
        
        # Make meshdiag test depend on mesh pipeline
        set_tests_properties(remote_generate_meshdiag_pi_remote_mpi8 PROPERTIES
            FIXTURES_REQUIRED "mesh_pi_remote_8"
            LABELS "generate_meshdiag"
        )
    else()
        message(STATUS "PI remote mesh diagnostics test skipped - BUILD_MESHDIAG=OFF")
    endif()
    
else()
    message(STATUS "PI remote mesh tests skipped - no download capability (wget or curl required)")
endif()

# =============================================================================
# NG5 Remote Mesh Tests (High Resolution ~5km)
# =============================================================================

if(HAS_DOWNLOAD_CAPABILITY)
    message(STATUS "Adding NG5 remote mesh download test...")
    
    # Add mesh pipeline for ng5 (128 processes for very large mesh ~7.4M nodes)
    add_mesh_pipeline("ng5" PREFIX remote PROCESS_COUNTS 128)
    
    if(BUILD_MESHPARTITIONER)
        message(STATUS "Adding NG5 mesh partitioning test...")
        # Partitioning is included in add_mesh_pipeline above
    else()
        message(STATUS "NG5 mesh partitioning test skipped - BUILD_MESHPARTITIONER=OFF")
    endif()
    
    if(BUILD_MESHDIAG)
        message(STATUS "Adding NG5 mesh diagnostics test...")
        
        # Add meshdiag test that generates mesh diagnostics file (very large mesh ~7.4M nodes, needs long timeout)
        add_fesom_meshdiag_test_with_options(generate_meshdiag_ng5_mpi128 "ng5" "fesom"
            NP 128
            TIMEOUT 3600  # 60 minutes for very large NG5 mesh
            PREFIX remote
        )
        
        # Make meshdiag test depend on mesh pipeline
        set_tests_properties(remote_generate_meshdiag_ng5_mpi128 PROPERTIES
            FIXTURES_REQUIRED "mesh_ng5_128"
            LABELS "generate_meshdiag"
        )
    else()
        message(STATUS "NG5 mesh diagnostics test skipped - BUILD_MESHDIAG=OFF")
    endif()
    
else()
    message(STATUS "NG5 remote mesh tests skipped - no download capability (wget or curl required)")
endif()

# =============================================================================
# ORCA25 Remote Mesh Tests (Global ~0.25-degree)
# =============================================================================

if(HAS_DOWNLOAD_CAPABILITY)
    message(STATUS "Adding ORCA25 remote mesh download test...")
    
    # Add mesh pipeline for orca25 (128 processes for large global mesh)
    add_mesh_pipeline("orca25" PREFIX remote PROCESS_COUNTS 128)
    
    if(BUILD_MESHPARTITIONER)
        message(STATUS "Adding ORCA25 mesh partitioning test...")
        # Partitioning is included in add_mesh_pipeline above
    else()
        message(STATUS "ORCA25 mesh partitioning test skipped - BUILD_MESHPARTITIONER=OFF")
    endif()
    
    if(BUILD_MESHDIAG)
        message(STATUS "Adding ORCA25 mesh diagnostics test...")
        
        # Add meshdiag test that generates mesh diagnostics file (large global mesh, needs longer timeout)
        add_fesom_meshdiag_test_with_options(generate_meshdiag_orca25_mpi128 "orca25" "fesom"
            NP 128
            TIMEOUT 1800  # 30 minutes for large ORCA25 mesh
            PREFIX remote
        )
        
        # Make meshdiag test depend on mesh pipeline
        set_tests_properties(remote_generate_meshdiag_orca25_mpi128 PROPERTIES
            FIXTURES_REQUIRED "mesh_orca25_128"
            LABELS "generate_meshdiag"
        )
    else()
        message(STATUS "ORCA25 mesh diagnostics test skipped - BUILD_MESHDIAG=OFF")
    endif()
    
else()
    message(STATUS "ORCA25 remote mesh tests skipped - no download capability (wget or curl required)")
endif()