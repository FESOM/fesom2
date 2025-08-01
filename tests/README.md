# FESOM2 Testing Framework

This directory contains the testing framework for FESOM2, including unit tests, integration tests, and system tests.

## Table of Contents
- [Quick Start](#quick-start)
- [CMake Configuration](#cmake-configuration)
- [Running Tests](#running-tests)
- [Test Organization](#test-organization)
- [Test Data Setup](#test-data-setup)
- [Debugging Tests](#debugging-tests)
- [Adding New Tests](#adding-new-tests)

## Quick Start

1. **Build with testing enabled** (default):
   ```bash
   # Create build directory and configure
   mkdir -p build && cd build
   cmake -DBUILD_TESTING=ON ..
   
   # Build the code and tests
   make -j$(nproc)
   
   # Build specific test targets
   make fesom_test  # Build all tests
   make fesom.x     # Build main executable
   ```

2. **Run tests**:
   ```bash
   # Run all tests with verbose output
   ctest --output-on-failure -V
   
   # Run tests matching a pattern
   ctest -R "integration"       # Run integration tests
   ctest -R "mpi"              # Run all MPI tests
   ctest -R "basic"            # Run basic tests
   
   # Run specific test by name
   ctest -R integration_basic_run_mpi8
   
   # Run tests with more output
   ctest --output-on-failure --verbose
   
   # Run tests in parallel (e.g., 4 jobs)
   ctest -j4 --output-on-failure
   ```

3. **Use convenience targets**:
   ```bash
   make test                    # Run all tests (same as 'ctest')
   make run_tests              # Alias for 'make test'
   make run_integration_tests  # Run only integration tests
   make run_mpi_tests          # Run only MPI tests
   make run_unit_tests         # Run only unit tests
   ```

## CMake Configuration

### Build Configuration Options

When configuring the build with CMake, you can use these options:

```bash
cmake \
  -DBUILD_TESTING=ON \
  -DENABLE_MPI_TESTS=ON \
  -DMPIEXEC_MAX_NUMPROCS=8 \
  -DCMAKE_BUILD_TYPE=Debug \
  ..
```

Key CMake options:
- `BUILD_TESTING`: Enable testing (default: ON)
- `ENABLE_MPI_TESTS`: Enable MPI tests (default: ON)
- `MPIEXEC_MAX_NUMPROCS`: Maximum number of MPI processes for tests
- `CMAKE_BUILD_TYPE`: Build type (Debug, Release, RelWithDebInfo, MinSizeRel)
- `TEST_TIMEOUT`: Global test timeout in seconds (default: 600)

### Building Specific Test Targets

```bash
# Build all tests
make fesom_test

# Build specific test
make integration_basic_run_mpi8

# Build and run a specific test
make integration_basic_run_mpi8 && ctest -R integration_basic_run_mpi8 -V
```

## Running Tests

### Listing Available Tests

```bash
# List all available tests
ctest -N

# List tests matching a pattern
ctest -N -R "integration"

# Show detailed information about tests
ctest --show-only=json-v1
```

### Running Specific Tests

```bash
# Run a single test
ctest -R integration_basic_run_mpi8

# Run tests matching a pattern
ctest -R "mpi"

# Run tests excluding certain patterns
ctest -E "long_running|performance"

# Run tests with labels
ctest -L "mpi"
```

### Test Output and Debugging

```bash
# Show test output on failure
ctest --output-on-failure

# Run with verbose output
ctest -V

# Run with debug output
ctest --debug

# Run with extra output
ctest --extra-verbose
```

### Running Tests in Parallel

```bash
# Run tests in parallel (4 jobs)
ctest -j4

# Run tests with timeout (seconds)
ctest --timeout 30

# Run tests repeatedly (good for flaky tests)
ctest --repeat until-pass:3  # Run up to 3 times until pass
ctest --repeat until-fail:3  # Run up to 3 times until failure
ctest --repeat after-timeout:3  # Repeat only after timeout
```

## Test Organization

Tests are organized into different categories:

1. **Unit Tests**: Tests for individual components/functions
2. **Integration Tests**: Tests that verify components work together
3. **MPI Tests**: Tests that require MPI parallel execution
4. **System Tests**: End-to-end tests of the complete system

### Test Naming Convention

- `unit_*`: Unit tests
- `integration_*`: Integration tests
- `*_mpi*`: MPI parallel tests
- `*_basic_*`: Basic functionality tests
- `*_init_*`: Initialization tests

## Test Data Setup

For the integration tests to run successfully, you need minimal test data:

1. **Test mesh**:
   ```bash
   # Link to existing test mesh
   mkdir -p tests/data/meshes
   ln -s /path/to/fesom2/test/meshes/pi tests/data/meshes/
   ```

2. **Test forcing data**:
   ```bash
   # Link to minimal forcing data
   mkdir -p tests/data/forcing
   ln -s /path/to/forcing/global tests/data/forcing/
   ```

## Debugging Tests

### Running Tests in Debugger

```bash
# Run a test under gdb
export TEST_EXECUTABLE="$(pwd)/bin/fesom.x"
cd tests/integration/integration_basic_run_mpi8
gdb --args mpirun -np 2 $TEST_EXECUTABLE

# Or use the test script directly
cd build
gdb --args ctest -V -R integration_basic_run_mpi8
```

### Enabling Verbose Output

Add these to your `namelist.config`:
```
&io_nml
    verbose = .true.
    debug_level = 2
/
```

## Adding New Tests

1. **Unit Tests**:
   - Add new test files in `tests/unit/`
   - Follow the naming convention `test_*.f90`

2. **Integration Tests**:
   - Add new test configurations in `tests/integration/`
   - Update `tests/integration/CMakeLists.txt`

3. **MPI Tests**:
   - Use the `add_fesom_test` function with `MPI_TEST` option
   - Specify the number of processes with `NP`

Example of adding a new test:

```cmake
# In tests/integration/CMakeLists.txt
add_fesom_test(integration_my_new_test
    MPI_TEST
    NP 4
    TIMEOUT 300  # 5 minutes
    COMMAND_ARGS "--my-arg=value"
)
```

## Troubleshooting

### Common Issues

1. **Test Fails with MPI Errors**:
   - Ensure MPI is properly installed and in your PATH
   - Check that the number of processes matches the mesh partitioning

2. **Missing Test Data**:
   - Verify that test data is correctly linked in `tests/data/`
   - Check the test output for missing files

3. **Test Timeout**:
   - Increase the test timeout with `ctest --timeout 600`
   - Or modify the test's timeout in the CMake configuration

4. **Debugging Hanging Tests**:
   - Run with `ctest -V` for verbose output
   - Check system logs for OOM killer or other system events

### Getting Help

For issues with the test framework:
1. Check the test output in `build/Testing/Temporary/LastTest.log`
2. Consult the FESOM2 documentation
3. Open an issue on the FESOM2 GitHub repository

## Test Data Setup

For the integration tests to run successfully, you need minimal test data:

1. **Add test mesh**:
   ```bash
   # Link to existing test mesh
   ln -s /path/to/fesom2/test/meshes/test_mesh ./data/meshes/
   
   # Or copy minimal mesh files
   cp -r /path/to/minimal/mesh ./data/meshes/test_mesh
   ```

2. **Add test forcing data**:
   ```bash
   # Link to minimal forcing data
   ln -s /path/to/minimal/forcing ./data/forcing/
   ```

## Available Tests

- `integration_basic_run`: Serial FESOM execution test
- `integration_init_test`: Quick initialization test (2 min timeout)
- `integration_basic_run_mpi2`: MPI test with 2 processes
- `integration_basic_run_mpi4`: MPI test with 4 processes
- `integration_namelist_config`: Verify namelist configuration
- `integration_test_data_check`: Validate test data setup

## Test Configuration

Configure tests with CMake options:

```bash
# Disable integration tests
cmake -DENABLE_INTEGRATION_TESTS=OFF ..

# Disable MPI tests
cmake -DENABLE_MPI_TESTS=OFF ..

# Set custom timeout (default: 600s)
cmake -DTEST_TIMEOUT=300 ..
```

## Troubleshooting

1. **Missing test data**: Tests will warn about missing data but still pass for initial setup
2. **MPI not found**: Ensure MPI is installed and mpiexec/mpirun is in PATH
3. **Test timeouts**: Increase TEST_TIMEOUT for slower systems
4. **Build failures**: Ensure FESOM builds successfully before running tests

## Test Output

Each integration test creates its own directory with:
- Modified namelists pointing to test data
- `test_output.log`: Standard output from FESOM
- `test_error.log`: Error output from FESOM
- `results/`: Directory for FESOM output files

The test data directory should contain the following subdirectories:
- `meshes/`: Test mesh files
- `forcing/`: Test forcing/climate data files
