# FESOM2 ecbundle Tests

This directory contains CTest tests for building FESOM2 using ecbundle.

## Overview

ecbundle is a tool from ECMWF that helps manage and build software bundles. It provides two main scripts:
- `ecbundle-create`: Downloads and checks out sources based on a bundle configuration
- `ecbundle-build`: Builds the software using CMake with various options

## Prerequisites

1. **Install ecbundle**: Follow the installation instructions at https://github.com/ecmwf/ecbundle
2. **Ensure ecbundle tools are in PATH**: The tests look for `ecbundle-create` and `ecbundle-build` executables

## Available Tests

The following tests are available:

### 1. `ecbundle_basic_build`
Tests the basic ecbundle build process using default options from `bundle.yml`.

**What it tests:**
- `ecbundle-create` successfully downloads FESOM sources
- `ecbundle-build` successfully builds FESOM with default options
- FESOM executable is created and is executable
- Basic functionality test (--help command)

### 2. `ecbundle_with_meshpart`
Tests building FESOM with meshpartitioner enabled using `--with-meshpart=ON`.

**What it tests:**
- All basic build functionality
- Meshpartitioner executable is created and is executable
- Both FESOM and meshpartitioner executables are functional

### 3. `ecbundle_with_omp`
Tests building FESOM with OpenMP support using `--with-omp-parallel`.

**What it tests:**
- All basic build functionality
- Checks for OpenMP library linkage using `ldd` and `nm`
- Verifies OpenMP support was properly configured

### 4. `ecbundle_with_testing`
Tests building FESOM with testing enabled using `--with-testing=ON`.

**What it tests:**
- All basic build functionality
- Checks for test executables and test files
- Attempts to run tests using CTest if available

## Running the Tests

### From CMake build directory:
```bash
# Run all ecbundle tests
ctest -L ecbundle

# Run a specific test
ctest -R ecbundle_basic_build

# Run with verbose output
ctest -L ecbundle --output-on-failure
```

### Using the convenience target:
```bash
# Run all ecbundle tests
make run_ecbundle_tests
```

### From the test directory:
```bash
# Run all ecbundle tests
cd test/ecbundle
ctest --output-on-failure
```

## Test Configuration

### Timeout
Each test has a 30-minute timeout (1800 seconds) to account for download and build times.

### Test Directory Structure
Tests create the following structure in `${CMAKE_BINARY_DIR}/test_ecbundle/`:
```
test_ecbundle/
├── basic/           # Basic build test
├── meshpart/        # Meshpartitioner test
├── omp/            # OpenMP test
└── testing/        # Testing enabled test
```

Each subdirectory contains:
- `sources/` - Downloaded source code
- `build/` - Build artifacts and executables

## Bundle Configuration

The tests use the `bundle.yml` file from the FESOM2 root directory, which defines:
- FESOM2 source repository and version
- Default CMake options
- Available build options (meshpart, OpenMP, testing)

## Troubleshooting

### ecbundle tools not found
If you see a warning about ecbundle tools not being found:
1. Ensure ecbundle is properly installed
2. Check that `ecbundle-create` and `ecbundle-build` are in your PATH
3. Try running `which ecbundle-create` and `which ecbundle-build`

### Build failures
- Check that all required dependencies are installed
- Verify network connectivity for downloading sources
- Check build logs in the test directories for specific error messages

### Test timeouts
- Increase the timeout in the CMakeLists.txt if builds take longer on your system
- Check system resources (CPU, memory, disk space)

## Adding New Tests

To add a new ecbundle test:

1. Create a new test script (e.g., `test_ecbundle_newfeature.cmake`)
2. Add the test to `CMakeLists.txt` following the existing pattern
3. Update this README with documentation for the new test

## Integration with CI/CD

These tests can be integrated into continuous integration pipelines to ensure:
- ecbundle builds work correctly
- Different build configurations are tested
- Build artifacts are properly created
- No regressions in the build process 