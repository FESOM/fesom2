# Remote Mesh Pipeline Tests

Comprehensive tests for the complete remote mesh pipeline: download → partition → validate.

## Overview

These tests validate the entire workflow for using remote meshes in FESOM2, from initial download through partitioning to final validation using mesh diagnostics.

## Test Pipeline Structure

Each remote mesh follows a **3-stage pipeline**:

1. **Download** (`remote_download_mesh_{name}`) - Download and extract mesh from remote repository
2. **Partition** (`remote_partition_mesh_{name}_{procs}`) - Partition mesh for specified MPI process count  
3. **Validate** (`remote_generate_meshdiag_{name}_mpi{procs}`) - Generate mesh diagnostics to verify setup

## Supported Meshes

### CORE2 (Medium Global Mesh)
- **Size**: ~127k nodes, global ocean mesh
- **Process Count**: 16
- **Timeout**: 5 minutes for meshdiag
- **Tests**: `remote_*_core2_*`

### DARS (Large Regional Mesh)  
- **Size**: ~3.2M nodes, Data Assimilation Research System mesh
- **Process Count**: 128
- **Timeout**: 30 minutes for meshdiag
- **Tests**: `remote_*_dars_*`

### PI Remote (Small Test Mesh)
- **Size**: ~3k nodes, same as local PI but downloaded
- **Process Count**: 8  
- **Timeout**: 5 minutes for meshdiag
- **Tests**: `remote_*_pi_remote_*`

### NG5 (Very Large High-Resolution Mesh)
- **Size**: ~7.4M nodes, ~5km resolution global mesh
- **Process Count**: 128
- **Timeout**: 60 minutes for meshdiag  
- **Tests**: `remote_*_ng5_*`

### ORCA25 (Large Global Mesh)
- **Size**: TBD, ~0.25 degree global mesh
- **Process Count**: 128
- **Timeout**: 30 minutes for meshdiag
- **Tests**: `remote_*_orca25_*`

## Requirements

### System Dependencies
```bash
# Download capability (at least one required)
wget  # Preferred
curl  # Alternative

# Build options
cmake .. -DBUILD_TESTING=ON \
         -DBUILD_MESHPARTITIONER=ON \  # For partition step
         -DBUILD_MESHDIAG=ON           # For validation step
```

### Runtime Requirements
- **Internet connection** for mesh downloads
- **Sufficient disk space** (meshes range from MB to GB)
- **Sufficient RAM** for large mesh processing (especially NG5)

## Component Detection

The system automatically detects available components and reports status:

```
Remote mesh pipeline component availability:
  Download capability: TRUE
    - wget: /usr/bin/wget
    - curl: /usr/bin/curl  
  Mesh partitioner: ON
  Mesh diagnostics: ON
```

### Missing Component Messages
- `Download capability: FALSE` → Install wget or curl
- `Mesh partitioner: OFF` → Add `-DBUILD_MESHPARTITIONER=ON`
- `Mesh diagnostics: OFF` → Add `-DBUILD_MESHDIAG=ON`

## Running Tests

### Complete Pipeline for Single Mesh
```bash
# Run complete CORE2 pipeline (download → partition → validate)
ctest -R "remote.*core2" -V

# Run complete DARS pipeline  
ctest -R "remote.*dars" -V
```

### By Pipeline Stage
```bash
# All download tests
ctest -L "download_mesh" -V

# All partition tests  
ctest -L "partition_mesh" -V

# All mesh diagnostics tests
ctest -L "generate_meshdiag" -V
```

### Individual Tests
```bash
# Download specific mesh
ctest -R remote_download_mesh_core2 -V

# Partition specific mesh
ctest -R remote_partition_mesh_core2_16 -V

# Generate diagnostics for specific mesh
ctest -R remote_generate_meshdiag_core2_mpi16 -V
```

## Test Dependencies

Tests use **CMake fixtures** to ensure proper ordering:

```
Download → Partition → Validate
    ↓         ↓         ↓
(no deps) → mesh_X → mesh_X_N
```

- **Download tests** create fixture `mesh_{name}`
- **Partition tests** require `mesh_{name}`, create `mesh_{name}_{procs}`
- **Meshdiag tests** require `mesh_{name}_{procs}`

## Output Files

### Downloaded Meshes
```
tests/data/MESHES/
├── core2/          # Downloaded and extracted
├── dars/           # Downloaded and extracted  
├── pi_remote/      # Downloaded and extracted
├── ng5/            # Downloaded and extracted
└── orca25/         # Downloaded and extracted
```

### Generated Partitions
```
tests/data/MESHES/{mesh}/
└── dist_{procs}/
    ├── rpart.out           # Main partition file
    ├── part.out            # Partition mapping
    ├── my_list*.out        # Process-specific node lists
    └── com_info*.out       # Communication information
```

### Mesh Diagnostics
```
build/tests/remote/remote_generate_meshdiag_{mesh}_mpi{procs}/
└── results/
    └── fesom.mesh.diag.nc  # Consistent filename across all tests
```

## Performance Expectations

### Download Times (depends on connection)
- **Small meshes** (pi_remote): ~1 minute
- **Medium meshes** (core2): ~5 minutes
- **Large meshes** (dars, ng5, orca25): ~10-30 minutes

### Partition Times
- **Small meshes**: ~30 seconds
- **Medium meshes**: ~2-5 minutes  
- **Large meshes**: ~10-30 minutes

### Meshdiag Times
- **Small meshes**: ~2 minutes
- **Medium meshes**: ~5 minutes
- **Large meshes**: ~30-60 minutes (hence the longer timeouts)

## Troubleshooting

### Download Issues
```bash
# Test download manually
wget https://swift.dkrz.de/v1/dkrz_035d8f6ff058403bb42f8302e6badfbc/FESOM_MESHES2.1/core2.tar.gz

# Check network connectivity
ping swift.dkrz.de
```

### Partition Issues  
```bash
# Check if partitioner executable exists
ls -la build/bin/fesom_meshpart

# Test partitioner manually
cd tests/data/MESHES/pi
mkdir -p dist_test
cd dist_test
../../../../build/bin/fesom_meshpart .. 4
```

### Meshdiag Issues
```bash
# Check if meshdiag executable exists
ls -la build/bin/fesom_meshdiag

# Test meshdiag manually
cd build/tests/remote/remote_generate_meshdiag_core2_mpi16
mpirun -np 16 ../../../bin/fesom_meshdiag
```

## Development

### Adding New Remote Meshes

1. **Add to mesh registry** (`tests/mesh_registry.json`)
2. **Add pipeline in `CMakeLists.txt`**:
   ```cmake
   if(HAS_DOWNLOAD_CAPABILITY)
       add_mesh_pipeline("new_mesh" PREFIX remote PROCESS_COUNTS 16)
       
       if(BUILD_MESHDIAG)
           add_fesom_meshdiag_test_with_options(generate_meshdiag_new_mesh_mpi16 "new_mesh" "fesom"
               NP 16
               TIMEOUT 600
               PREFIX remote
           )
       endif()
   endif()
   ```

### Best Practices
- **Use appropriate process counts** based on mesh size
- **Set realistic timeouts** based on mesh complexity  
- **Include conditional logic** for missing components
- **Use consistent naming** with `remote_` prefix
- **Generate clean output** (`fesom.mesh.diag.nc`)