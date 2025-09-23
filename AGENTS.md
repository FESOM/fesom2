# FESOM2 Agent Information

This document contains important information for AI agents working with this FESOM2 codebase.

## Recent Changes

### Separate Input/Output Restart Paths (Issue #629)
- **Status**: IMPLEMENTED  
- **Branch**: Current branch
- **Description**: Added separate RestartInPath and RestartOutPath configuration parameters to prevent overwriting of input restart files, improving reproducibility.

#### Key Changes:
1. **Configuration**: Added `RestartInPath` and `RestartOutPath` to `namelist.config` `/paths/` section
2. **Backward Compatibility**: Empty values default to `ResultPath` (legacy behavior)
3. **Code Refactoring**: Split monolithic `restart()` subroutine into:
   - `read_initial_conditions()` - reads from RestartInPath
   - `write_initial_conditions()` - writes to RestartOutPath  
   - `restart()` - wrapper for backward compatibility

#### Files Modified:
- `src/gen_modules_config.F90` - Added new configuration variables
- `src/gen_model_setup.F90` - Added default logic and path validation
- `src/io_restart.F90` - Major refactoring with separate read/write subroutines
- `src/icepack_drivers/icedrv_io.F90` - Updated icepack restart paths

#### Testing:
- Backward compatibility: Existing namelists work without changes
- New functionality: Can specify separate input/output restart directories
- All restart types supported: NetCDF, raw, and binary

## Commands

- **Build:** `cmake -S . -B build --preset gcc && cmake --build build`
