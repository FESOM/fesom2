# FESOM2 Testing Framework

This directory contains the comprehensive testing framework for FESOM2, organized into focused test groups for different aspects of the model.

## Table of Contents
- [Quick Start](#quick-start)
- [Test Enablement](#test-enablement)
- [Test Organization](#test-organization)
- [Listing and Filtering Tests](#listing-and-filtering-tests)
- [Running Tests](#running-tests)
- [Custom Test Targets](#custom-test-targets)
- [Test Configuration](#test-configuration)
- [Test Data](#test-data)
- [Troubleshooting](#troubleshooting)

## Quick Start

```bash
# 1. Build with all testing components
./configure.sh ubuntu -DBUILD_TESTING=ON -DBUILD_MESHDIAG=ON -DBUILD_MESHPARTITIONER=ON

# 2. List available tests (see what's enabled)
ctest -N

# 3. Run all tests
ctest --output-on-failure
```

## Test Enablement

**Tests are conditionally enabled based on build configuration.** The test system automatically detects available components and includes only relevant tests.

### Core Testing
```bash
-DBUILD_TESTING=ON              # Required - enables basic test framework
```
**Enables**: Integration tests, basic infrastructure

### Component-Specific Testing
```bash
-DBUILD_MESHDIAG=ON             # Enables mesh diagnostics tests
-DBUILD_MESHPARTITIONER=ON      # Enables mesh partitioner tests
```

### System Dependencies
The system **automatically detects** and reports:
- **Download capability**: wget or curl (for remote tests)
- **MPI availability**: mpiexec/mpirun (for MPI tests)

### Detection Report Example
```
Remote mesh pipeline component availability:
  Download capability: TRUE (wget: /usr/bin/wget)
  Mesh partitioner: ON
  Mesh diagnostics: ON
```

**Missing components result in clear skip messages** with instructions to enable them.

## Test Organization

The test framework is organized into **4 focused directories**:

### 📁 `integration/` - Local Integration Tests
**Purpose**: Core FESOM functionality using local mesh files
- `integration_pi_mpi2` - PI mesh with 2 processes
- `integration_pi_mpi8` - PI mesh with 8 processes  
- `integration_pi_cavity_mpi2` - PI cavity mesh with 2 processes
- `integration_meshdiag_pi_mpi2/8` - Generate mesh diagnostics (`fesom.mesh.diag.nc`)

**Requirements**: Local mesh files in `tests/data/MESHES/`

### 📁 `remote/` - Remote Mesh Pipeline Tests
**Purpose**: Test complete pipeline for downloading, partitioning, and validating remote meshes
- `remote_download_mesh_{name}` - Download mesh from remote repository
- `remote_partition_mesh_{name}_{procs}` - Partition mesh for specified process count
- `remote_generate_meshdiag_{name}_mpi{procs}` - Generate mesh diagnostics (`fesom.mesh.diag.nc`)

**Supported Meshes**: core2, dars, pi_remote, ng5, orca25
**Requirements**: wget/curl, internet connection

### 📁 `meshpartitioner/` - Local Mesh Partitioning Tests  
**Purpose**: Test mesh partitioner tool on local meshes
- `meshpartitioner_partition_local_pi_mesh_4` - Create 4-process partition for PI mesh
- `meshpartitioner_partition_local_pi_cavity_mesh_4` - Create 4-process partition for PI cavity

**Requirements**: `BUILD_MESHPARTITIONER=ON`

### 📁 `ecbundle/` - External Tool Integration Tests
**Purpose**: Test ECBundle integration and build system
- `ecbundle_basic_build` - Basic ECBundle functionality
- `ecbundle_with_meshpart` - ECBundle with mesh partitioner
- `ecbundle_with_omp` - ECBundle with OpenMP
- `ecbundle_with_testing` - ECBundle with testing enabled

## Listing and Filtering Tests

**Before running tests, it's helpful to see what's available** based on your build configuration.

### List All Available Tests
```bash
ctest -N                    # List all tests (typically 20-26 depending on build options)
ctest -N -R "integration_"  # List integration tests only
ctest -N -R "remote_"       # List remote pipeline tests  
ctest -N -R "meshpartitioner_" # List mesh partitioner tests
```

### Filter by Test Group
```bash
# Show tests by prefix
ctest -N -R "integration_"     # All integration tests
ctest -N -R "remote_"          # All remote mesh pipeline tests  
ctest -N -R "meshpartitioner_" # All mesh partitioner tests
ctest -N -R "ecbundle_"        # All ecbundle tests
```

### Filter by Functionality
```bash
# Show tests by label
ctest -N -L "download_mesh"     # All mesh download tests
ctest -N -L "partition_mesh"    # All mesh partition tests  
ctest -N -L "generate_meshdiag" # All meshdiag generation tests
ctest -N -R "_mpi"              # All MPI tests
```

### Example Output
```bash
$ ctest -N
Test project /path/to/build
  Test  #1: integration_pi_mpi2
  Test  #2: integration_pi_mpi8
  Test  #3: integration_pi_cavity_mpi2
  Test  #4: integration_meshdiag_pi_mpi2
  Test  #5: integration_meshdiag_pi_mpi8
  Test  #6: meshpartitioner_partition_local_pi_mesh_4
  Test  #7: meshpartitioner_partition_local_pi_cavity_mesh_4
  Test  #8: remote_download_mesh_core2
  ...
Total Tests: 26
```

## Running Tests

### Run by Test Group
```bash
ctest -R "integration_"     # All integration tests
ctest -R "remote_"          # All remote mesh pipeline tests  
ctest -R "meshpartitioner_" # All mesh partitioner tests
ctest -R "ecbundle_"        # All ecbundle tests
```

### Run by Functionality
```bash
ctest -L "download_mesh"     # All mesh download tests
ctest -L "partition_mesh"    # All mesh partition tests  
ctest -L "generate_meshdiag" # All meshdiag generation tests
ctest -R "_mpi"              # All MPI tests
```

### Run Specific Tests
```bash
ctest -R integration_pi_mpi2 -V                    # Single integration test
ctest -R remote_generate_meshdiag_core2_mpi16 -V   # Single meshdiag test
ctest -R meshpartitioner_partition_local_pi_mesh_4 # Single partition test
```

### Run All Tests
```bash
ctest --output-on-failure   # Run all available tests
```

## Custom Test Targets

Convenient make targets for common test workflows:

```bash
make run_tests                  # Run all tests
make run_integration_tests      # Integration tests only
make run_mpi_tests             # All MPI tests
make run_meshdiag_tests        # All mesh diagnostics tests
make run_meshpartitioner_tests # Mesh partitioner tests only
```

## Test Configuration

### CMake Build Options
```bash
# Core testing
-DBUILD_TESTING=ON              # Enable all tests (required)

# Component-specific testing  
-DBUILD_MESHDIAG=ON             # Enable mesh diagnostics tests
-DBUILD_MESHPARTITIONER=ON      # Enable mesh partitioner tests

# MPI configuration
-DENABLE_MPI_TESTS=ON           # Enable MPI tests (default: ON)
-DMPIEXEC_MAX_NUMPROCS=128      # Max MPI processes available

# Test timing
-DTEST_TIMEOUT=600              # Default test timeout (seconds)
```

### Component Dependencies
The test system automatically detects and reports available components:

```
Remote mesh pipeline component availability:
  Download capability: TRUE (wget: /usr/bin/wget, curl: /usr/bin/curl)
  Mesh partitioner: ON
  Mesh diagnostics: ON
```

**Missing components are clearly reported** with instructions to enable them.

## Test Data

### Local Test Data (`tests/data/`)
```
tests/data/
├── MESHES/
│   ├── pi/          # PI test mesh (3k nodes)
│   ├── pi_cavity/   # PI cavity test mesh (7k nodes)
│   └── soufflet/    # Soufflet test mesh (3k nodes)
├── FORCING/         # Minimal forcing data
└── INITIAL/         # Initial condition data
```

### Remote Test Data
Remote meshes are **automatically downloaded** by the test pipeline:
- **CORE2**: ~127k nodes, global ocean mesh
- **DARS**: ~3.2M nodes, large regional mesh  
- **NG5**: ~7.4M nodes, high-resolution global mesh
- **ORCA25**: Large global mesh (~0.25 degree)
- **PI Remote**: Same as local PI but downloaded

## Process Count Guidelines

Tests use **appropriate process counts** based on mesh size:
- **Small meshes** (pi, pi_cavity): 2, 4, 8 processes
- **Medium meshes** (core2): 16 processes
- **Large meshes** (dars, orca25): 128 processes  
- **Very large meshes** (ng5): 128 processes

## Output Files

### Integration Tests
All mesh diagnostics tests create: **`fesom.mesh.diag.nc`**

Directory structure tells the complete story:
```
integration_meshdiag_pi_mpi2/results/fesom.mesh.diag.nc    # PI mesh, 2 processes
integration_meshdiag_pi_mpi8/results/fesom.mesh.diag.nc    # PI mesh, 8 processes
remote_generate_meshdiag_core2_mpi16/results/fesom.mesh.diag.nc    # CORE2 mesh, 16 processes
```

### Test Logs
Each test creates detailed logs:
- `test_output.log` - Standard output from FESOM/meshdiag
- `test_error.log` - Error output (usually empty for successful tests)
- `namelist.*` - Configured namelists used by the test

## Troubleshooting

### Common Issues

1. **"No download capability"**:
   ```bash
   sudo apt install wget curl
   ```

2. **"Mesh partitioner tests skipped"**:
   ```bash
   cmake .. -DBUILD_MESHPARTITIONER=ON
   ```

3. **"Mesh diagnostics test skipped"**:
   ```bash
   cmake .. -DBUILD_MESHDIAG=ON
   ```

4. **Large mesh tests timeout**:
   - DARS: 30-minute timeout (3.2M nodes)
   - NG5: 60-minute timeout (7.4M nodes)
   - These are appropriate for mesh complexity

5. **MPI test failures**:
   ```bash
   # Check MPI installation
   which mpiexec mpirun
   mpirun -np 2 echo "MPI works"
   ```

### Debug Test Execution
```bash
# Run with maximum verbosity
ctest -R integration_pi_mpi2 -V --output-on-failure

# Check specific test logs
cat build/tests/integration/integration_pi_mpi2/test_output.log

# Run test manually
cd build/tests/integration/integration_pi_mpi2
mpirun -np 2 ../../../../bin/fesom.x
```

## Test Development

### Adding New Tests

1. **Integration tests** - Add to `integration/CMakeLists.txt`
2. **Remote meshes** - Add to mesh registry JSON and `remote/CMakeLists.txt`  
3. **Mesh partitioner tests** - Add to `meshpartitioner/CMakeLists.txt`

### Best Practices
- **Use group prefixes** (`integration_`, `remote_`, `meshpartitioner_`)
- **Include mesh name** in test names for clarity
- **Set appropriate timeouts** based on mesh complexity
- **Use conditional logic** with clear messages for missing dependencies
- **Create clean output files** (`fesom.mesh.diag.nc`)

## Test Infrastructure

The test system uses `cmake/FesomTesting.cmake` (757 lines) which provides:
- Automatic namelist configuration
- Mesh pipeline management  
- Test script generation
- MPI setup and execution

**Note**: A simplification roadmap is documented in `TODO_TEST_INFRASTRUCTURE_SIMPLIFICATION.md` for future maintenance.
