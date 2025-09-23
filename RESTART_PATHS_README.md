# FESOM2 Separate Restart Input/Output Paths

## Overview

This feature implements separate input and output paths for restart files in FESOM2, solving issue #629. Previously, both restart reading and writing used the same directory (`ResultPath`), causing input restart files to be overwritten and breaking reproducibility.

## Configuration

### New Namelist Parameters

Add these parameters to your `namelist.config` file in the `/paths/` section:

```fortran
&paths
  MeshPath='./mesh/'
  ClimateDataPath='./hydrography/'  
  TideForcingPath='./tide_forcing/'
  ResultPath='./result/'
  MeshId='NONE'
  RestartInPath='./input_restarts/'     ! Optional: where to read restart files from
  RestartOutPath='./output_restarts/'   ! Optional: where to write restart files to
/
```

### Behavior

- **If both parameters are empty** (default): Uses `ResultPath` for both input and output (legacy behavior)
- **If only `RestartOutPath` is set**: Reads from `ResultPath`, writes to `RestartOutPath`
- **If only `RestartInPath` is set**: Reads from `RestartInPath`, writes to `ResultPath`  
- **If both are set**: Reads from `RestartInPath`, writes to `RestartOutPath`

### Backward Compatibility

Existing namelists without these parameters work exactly as before - no changes required.

## File Types Supported

All FESOM2 restart file types are supported:

1. **NetCDF Restart Files**: `.restart.nc` files for ocean, ice, bio, and icepack components
2. **Raw/Core Dump Files**: Binary dump files in `*_raw_restart/` directories
3. **Binary Restart Files**: Derived type binary files in `*_bin_restart/` directories

## Usage Examples

### Example 1: Reproducible Runs
```fortran
&paths
  ResultPath='./run1_output/'
  RestartInPath='./initial_conditions/'
  RestartOutPath='./run1_output/'
/
```
- Reads initial conditions from `./initial_conditions/`
- Writes all output including restarts to `./run1_output/`
- Input files remain untouched for future reproducible runs

### Example 2: Restart Chain
```fortran
&paths  
  ResultPath='./run2_output/'
  RestartInPath='./run1_output/'
  RestartOutPath='./run2_output/'
/
```
- Continues from restart files in `./run1_output/`
- Writes new output to `./run2_output/`

### Example 3: Legacy Mode (Default)
```fortran
&paths
  ResultPath='./result/'
  ! RestartInPath and RestartOutPath not specified
/
```
- Reads and writes from same directory (original behavior)

## Implementation Details

### Code Structure

The original monolithic `restart()` subroutine has been split into:

1. **`read_initial_conditions()`**: Handles reading restart files from `RestartInPath`
2. **`write_initial_conditions()`**: Handles writing restart files to `RestartOutPath`  
3. **`restart()`**: Wrapper maintaining backward compatibility

### Path Construction

New centralized helper functions build consistent file paths:
- `nc_restart_path()` - For NetCDF restart files
- `build_raw_restart_dirpath()` - For raw restart directories
- `build_bin_restart_dirpath()` - For binary restart directories

### File Path Examples

For a run with ID `test1` and year `2000`:

**NetCDF files:**
- Input: `{RestartInPath}/test1.2000.oce.restart.nc`
- Output: `{RestartOutPath}/test1.2000.oce.restart.nc`

**Raw restart files:**
- Input: `{RestartInPath}/test1_raw_restart/np{nproc}/`
- Output: `{RestartOutPath}/test1_raw_restart/np{nproc}/`

## Testing

### Backward Compatibility Tests
- Existing namelists produce identical results
- All restart file types work as before

### New Functionality Tests  
- Separate input/output directories work correctly
- Input files remain unmodified
- All file formats (NetCDF, raw, binary) respect new paths

## Error Handling

The implementation includes:
- Automatic path normalization (trailing slashes)
- User feedback on configuration
- Graceful fallback to legacy behavior
- Clear error messages for path issues

## Migration Guide

### For Existing Users
No action required - existing configurations work unchanged.

### For New Reproducible Workflows
1. Set `RestartInPath` to your archive of initial conditions
2. Set `RestartOutPath` to your run-specific output directory  
3. Your initial conditions will never be overwritten

## Performance Impact

Minimal - only adds a few string operations during initialization. No impact on computational performance or I/O patterns.

## Technical Notes

### Clock Files
Clock files (`.clock`) continue to be written to `ResultPath` as they are operational metadata, not restart data.

### Directory Creation
Output directories are created automatically if they don't exist.

### MPI Safety
All path operations are coordinated across MPI ranks to ensure consistency.

## Files Modified

- `src/gen_modules_config.F90` - Configuration variables
- `src/gen_model_setup.F90` - Default handling and validation
- `src/io_restart.F90` - Core restart logic refactoring
- `src/icepack_drivers/icedrv_io.F90` - Icepack integration

## Benefits

1. **Reproducibility**: Input restart files never overwritten
2. **Workflow Flexibility**: Easy to organize input/output separately
3. **Debugging Support**: Compare input vs output restart files
4. **Backward Compatibility**: Zero disruption to existing users
5. **Code Quality**: Cleaner separation of read/write operations

## Future Enhancements

Potential future improvements could include:
- Automatic restart file discovery and validation
- Restart file compression options
- Enhanced metadata tracking
- Integration with workflow management systems
