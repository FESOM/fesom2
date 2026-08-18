# Mesh Partitioner Tests

Tests for the FESOM2 mesh partitioner tool (`fesom_meshpart`) using local mesh files.

## Overview

These tests validate the mesh partitioner functionality by creating missing partitions for local meshes and verifying the partitioner tool works correctly.

## Available Tests

### `meshpartitioner_partition_local_pi_mesh_4`
- **Purpose**: Create 4-process partition for local PI mesh
- **Mesh**: PI test mesh (~3k nodes) 
- **Duration**: ~30 seconds
- **Output**: `tests/data/MESHES/pi/dist_4/`
- **Status**: **Creates new partition** (was missing)

### `meshpartitioner_partition_local_pi_cavity_mesh_4`  
- **Purpose**: Create 4-process partition for local PI cavity mesh
- **Mesh**: PI cavity test mesh (~7k nodes)
- **Duration**: ~30 seconds  
- **Output**: `tests/data/MESHES/pi_cavity/dist_4/`
- **Status**: **Creates new partition** (was missing)

## Requirements

### Build Configuration
```bash
# Required for mesh partitioner tests
cmake .. -DBUILD_TESTING=ON \
         -DBUILD_MESHPARTITIONER=ON
```

### System Requirements
- **Local mesh files** in `tests/data/MESHES/`
- **METIS library** (automatically built with FESOM2)
- **Sufficient disk space** for partition files

## Test Execution

### Run All Mesh Partitioner Tests
```bash
make run_meshpartitioner_tests
# Or: ctest -L "meshpartitioner"
```

### Run Individual Tests
```bash
# Create PI mesh partition for 4 processes
ctest -R meshpartitioner_partition_local_pi_mesh_4 -V

# Create PI cavity mesh partition for 4 processes  
ctest -R meshpartitioner_partition_local_pi_cavity_mesh_4 -V
```

## Output Structure

### Created Partition Files
```
tests/data/MESHES/pi/dist_4/
├── rpart.out           # Main partition file (required by FESOM)
├── part.out            # Partition mapping
├── my_list00000.out    # Node list for process 0
├── my_list00001.out    # Node list for process 1
├── my_list00002.out    # Node list for process 2
├── my_list00003.out    # Node list for process 3
├── com_info00000.out   # Communication info for process 0
├── com_info00001.out   # Communication info for process 1
├── com_info00002.out   # Communication info for process 2
├── com_info00003.out   # Communication info for process 3
├── namelist.config     # Namelist used during partitioning
├── partition_output.log # Partitioner stdout
└── partition_error.log  # Partitioner stderr
```

## Validation

### Automatic Validation
The tests automatically verify:
- ✅ **Required files created** (`rpart.out`, `part.out`)
- ✅ **File sizes are reasonable** (not zero-length)
- ✅ **Partitioner exit code** (0 = success)
- ✅ **Process-specific files** created for each MPI rank

### Manual Verification
```bash
# Check partition file contents
head tests/data/MESHES/pi/dist_4/rpart.out

# Verify partition can be used by FESOM
cd build/tests/integration/integration_pi_mpi2
# Edit namelist.config to use 4 processes and dist_4
mpirun -np 4 ../../../bin/fesom.x
```

## Purpose and Benefits

### Why These Tests Matter
1. **Fill gaps** - Creates missing 4-process partitions for local development
2. **Test partitioner tool** - Validates `fesom_meshpart` executable works correctly
3. **Enable 4-process testing** - Allows local development with 4 MPI processes
4. **Verify partition quality** - Ensures partitioner produces valid output

### Integration with Other Tests
- **Enables new integration tests** that can use 4 processes
- **Supports mesh partitioner development** - changes can be tested locally
- **Validates partitioner before remote use** - same tool used for remote mesh pipeline

## Troubleshooting

### Common Issues

1. **"BUILD_MESHPARTITIONER=OFF"**:
   ```bash
   cmake .. -DBUILD_MESHPARTITIONER=ON
   make fesom_meshpart
   ```

2. **Missing mesh files**:
   ```bash
   # Check required files exist
   ls tests/data/MESHES/pi/nod2d.out
   ls tests/data/MESHES/pi/elem2d.out  
   ls tests/data/MESHES/pi/aux3d.out
   ```

3. **Partitioner fails**:
   ```bash
   # Check partitioner executable
   ./build/bin/fesom_meshpart tests/data/MESHES/pi 2
   ```

4. **Permission issues**:
   ```bash
   # Ensure write permissions to mesh directory
   chmod +w tests/data/MESHES/pi
   ```

## Development

### Adding New Mesh Partitioner Tests

To test partitioning for a new mesh or process count:

```cmake
# Add to CMakeLists.txt
add_test(
    NAME meshpartitioner_partition_local_{mesh}_mesh_{procs}
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DBUILD_DIR=${CMAKE_BINARY_DIR}
        -DMESH_NAME={mesh}
        -DNUM_PROCESSES={procs}
        -P ${CMAKE_SOURCE_DIR}/tests/integration/mesh_partition.cmake
)

set_tests_properties(meshpartitioner_partition_local_{mesh}_mesh_{procs} PROPERTIES
    TIMEOUT 300
    LABELS "meshpartitioner"
)
```

### Integration with Remote Tests
- **Same partitioner tool** used for both local and remote mesh partitioning
- **Same partition script** (`tests/integration/mesh_partition.cmake`)
- **Validates tool** before using it in remote pipeline

## Notes

- **Tests only run when** `BUILD_MESHPARTITIONER=ON`
- **Uses existing local meshes** - no downloads required
- **Creates real partition files** that can be used by FESOM
- **Tests are fast** - small local meshes partition quickly
- **Enables 4-process development** - fills partition gaps for local testing