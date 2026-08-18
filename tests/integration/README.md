# Integration Tests

Local FESOM2 integration tests that validate core functionality using local mesh files.

## Test Overview

These tests verify that FESOM2 can successfully run with local meshes and generate expected outputs.

### Available Tests

1. **`integration_pi_mpi2`** - PI mesh with 2 MPI processes
   - **Mesh**: PI test mesh (~3k nodes)
   - **Duration**: ~5-10 minutes
   - **Purpose**: Basic FESOM simulation test, minimum MPI configuration

2. **`integration_pi_mpi8`** - PI mesh with 8 MPI processes  
   - **Mesh**: PI test mesh (~3k nodes)
   - **Duration**: ~5-10 minutes
   - **Purpose**: Multi-process FESOM simulation test

3. **`integration_pi_cavity_mpi2`** - PI cavity mesh with 2 MPI processes
   - **Mesh**: PI cavity test mesh (~7k nodes) 
   - **Duration**: ~5-10 minutes
   - **Purpose**: Test cavity functionality

4. **`integration_meshdiag_pi_mpi2`** - PI mesh diagnostics with 2 processes
   - **Mesh**: PI test mesh (~3k nodes)
   - **Duration**: ~2 minutes
   - **Output**: `fesom.mesh.diag.nc`
   - **Purpose**: Fast mesh validation and diagnostics

5. **`integration_meshdiag_pi_mpi8`** - PI mesh diagnostics with 8 processes
   - **Mesh**: PI test mesh (~3k nodes)
   - **Duration**: ~2 minutes  
   - **Output**: `fesom.mesh.diag.nc`
   - **Purpose**: Multi-process mesh validation

## Requirements

### Local Mesh Data
These tests require local mesh files in `tests/data/MESHES/`:
```
tests/data/MESHES/
├── pi/
│   ├── nod2d.out      # Node coordinates
│   ├── elem2d.out     # Element connectivity
│   ├── aux3d.out      # 3D mesh information
│   ├── dist_2/        # 2-process partition
│   ├── dist_4/        # 4-process partition (created by meshpartitioner tests)
│   └── dist_8/        # 8-process partition
└── pi_cavity/
    ├── nod2d.out
    ├── elem2d.out
    ├── aux3d.out
    └── dist_2/        # 2-process partition
```

### Build Configuration
```bash
# Required
-DBUILD_TESTING=ON

# For meshdiag tests
-DBUILD_MESHDIAG=ON  

# For MPI tests  
-DENABLE_MPI_TESTS=ON
```

## Running Tests

### All Integration Tests
```bash
make run_integration_tests
# Or: ctest -R "integration_"
```

### Individual Tests
```bash
# Full FESOM simulation tests
ctest -R integration_pi_mpi2 -V
ctest -R integration_pi_mpi8 -V
ctest -R integration_pi_cavity_mpi2 -V

# Fast mesh diagnostics tests
ctest -R integration_meshdiag_pi_mpi2 -V
ctest -R integration_meshdiag_pi_mpi8 -V
```

### Quick Validation
```bash
# Fast validation of mesh and basic functionality
make run_meshdiag_tests
```

## Test Output

Each integration test creates its own directory under `build/tests/integration/`:

```
build/tests/integration/
├── integration_pi_mpi2/
│   ├── namelist.config    # Configured for test
│   ├── namelist.forcing   # Forcing configuration
│   ├── namelist.ice       # Ice configuration
│   ├── namelist.tra       # Tracer configuration
│   ├── test_output.log    # FESOM stdout
│   ├── test_error.log     # FESOM stderr
│   └── results/           # FESOM output files
└── integration_meshdiag_pi_mpi2/
    ├── namelist.config    # Configured for test
    ├── test_output.log    # meshdiag stdout
    ├── test_error.log     # meshdiag stderr (usually empty)
    └── results/
        └── fesom.mesh.diag.nc  # Mesh diagnostics file
```

## Debugging

### Test Failure Investigation
```bash
# Check test logs
cat build/tests/integration/integration_pi_mpi2/test_output.log
cat build/tests/integration/integration_pi_mpi2/test_error.log

# Run test manually
cd build/tests/integration/integration_pi_mpi2
mpirun -np 2 ../../../../bin/fesom.x
```

### Common Issues

1. **Missing mesh files**: Ensure local mesh data is present in `tests/data/MESHES/`
2. **MPI errors**: Check MPI installation with `mpirun -np 2 echo "test"`
3. **Timeout**: Integration tests may take 5-10 minutes on slower systems

## Adding New Integration Tests

To add a new integration test:

1. **Edit `CMakeLists.txt`** in this directory
2. **Use the `add_fesom_test*` functions** from `cmake/FesomTesting.cmake`
3. **Follow naming convention**: `integration_{mesh}_mpi{n}` or `integration_meshdiag_{mesh}_mpi{n}`

Example:
```cmake
# Add new PI mesh test with 4 processes
add_fesom_test(integration_pi_mpi4
    MPI_TEST
    NP 4
    TIMEOUT 600
)
```

## Notes

- **Integration tests use local meshes only** - no downloading required
- **Mesh partitions must exist** for the specified process count
- **All meshdiag tests create `fesom.mesh.diag.nc`** for consistency
- **Test directories are unique** so runid can be simple ("fesom")
- **Tests validate core FESOM functionality** without external dependencies