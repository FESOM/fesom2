# FESOM2 Toy Ocean Model

This directory contains a minimal, simplified version of FESOM2 designed for idealized ocean-only test cases.

## Overview

The toy ocean model provides a streamlined implementation suitable for:
- Idealized channel configurations (e.g., Soufflet channel test case)
- Double gyre experiments
- Educational purposes
- Algorithm testing and development
- Performance benchmarking of core ocean physics

## What's Included

### Source Files

- **`fesom_toy_module.F90`** - Simplified module with three main subroutines:
  - `toy_fesom_init` - Initialize model (mesh, arrays, toy-specific initial conditions)
  - `toy_fesom_runloop` - Time-stepping loop (ocean physics only)
  - `toy_fesom_finalize` - Cleanup and statistics output

- **`fesom_toy_main.F90`** - Minimal driver program that calls the three subroutines

- **`CMakeLists.txt`** - Independent build configuration for the toy model
  - Creates `libfesom_toy.so` shared library
  - Builds `fesom_toy` executable
  - Handles installation and export

- **`fesom_toy-config.cmake.in`** - CMake package configuration template

### Key Simplifications

The toy model **removes** the following components compared to the full FESOM2:
- ❌ Sea ice model (replaced with dummy ice structure)
- ❌ Icebergs
- ❌ OASIS coupling
- ❌ OpenIFS coupling
- ❌ RECOM biogeochemistry
- ❌ Transient tracers (SF6, CFC-11, CFC-12, 14C)
- ❌ Global tides
- ❌ Age tracers
- ❌ Land ice coupling
- ❌ Icepack thermodynamics

The toy model **retains** essential ocean components:
- ✅ Full ocean dynamics (ALE formulation)
- ✅ Tracer advection-diffusion (T and S)
- ✅ Vertical mixing schemes (PP or KPP)
- ✅ MPI parallelization
- ✅ I/O and diagnostics
- ✅ Restart capability
- ✅ Toy-specific forcings (e.g., Soufflet zonal restoring)

### Size Comparison

- **Original `fesom_module.F90`**: 965 lines
- **Simplified `fesom_toy_module.F90`**: ~420 lines
- **Reduction**: ~56% smaller

## Building the Toy Model

The toy model can be built in two ways:

### Method 1: Build from FESOM2 Root Directory (Integrated Build)

#### Prerequisites

Same as the full FESOM2 model:
- CMake ≥ 3.16
- Fortran compiler (Intel, GNU, Cray, or NVHPC)
- MPI library
- NetCDF (C and Fortran)
- METIS (for mesh partitioning)

#### Build Instructions

1. **Configure with CMake** (from the FESOM2 root directory):

   ```bash
   mkdir -p build
   cd build
   cmake .. -DBUILD_TOY=ON
   ```

   Optional: combine with other build options:
   ```bash
   cmake .. -DBUILD_TOY=ON -DCMAKE_BUILD_TYPE=Release -DCVMIX=ON
   ```

2. **Build the executable**:

   ```bash
   make fesom_toy
   ```

   Or build everything including the toy model:
   ```bash
   make
   ```

3. **Install** (optional):

   ```bash
   make install
   ```

   The `fesom_toy` executable will be installed to `bin/` alongside `fesom.x`.

### Method 2: Standalone Build from dwarf/toy_ocean (Recommended for Development)

This method builds the toy model independently, following the FESOM2 dwarf pattern.

#### Build Instructions

1. **Setup symlinks to FESOM2 sources** (from `dwarf/toy_ocean/`):

   ```bash
   cd dwarf/toy_ocean
   ./dwarf_linkfiles.sh
   ```

   This creates symlinks to:
   - FESOM2 source files (97 minimal files)
   - Environment configuration (`env.sh`, `configure.sh`)
   - CMake build scripts

2. **Build the toy model**:

   ```bash
   ./configure.sh <platform> [options]
   ```

   Example for Ubuntu:
   ```bash
   ./configure.sh ubuntu -DCVMIX=OFF
   ```

   Example for other platforms:
   ```bash
   ./configure.sh levante.dkrz.de
   ```

3. **Build output**:

   After a successful build, you should have:
   - `build/src/fesom_toy` - The toy ocean executable
   - `build/src/libfesom_toy.so` - The toy model shared library (5.3 MB)

### Library Architecture

The toy model has an independent build structure located in `dwarf/toy_ocean/`:

```
dwarf/toy_ocean/
├── CMakeLists.txt              # Top-level build configuration
├── README.md                   # This file
├── dwarf_linkfiles.sh          # Script to create source symlinks
├── dwarf_ini/                  # Toy model core files
│   ├── CMakeLists.txt          # Build configuration (97 source files)
│   ├── fesom_toy_module.F90    # Simplified toy model module
│   ├── fesom_toy_main.F90      # Main program
│   └── fesom_toy-config.cmake.in
├── src/                        # Symlinks to FESOM2 sources (created by dwarf_linkfiles.sh)
├── work/                       # Working directory for simulations
└── build/                      # Build directory (created during build)
```

**Build dependency chain (standalone):**
```
fesom_toy (executable)
    └─→ libfesom_toy.so (independent toy model library)
            └─→ MPI, NetCDF, etc. (no dependency on main libfesom.so)
```

This architecture allows:
- **Fully independent compilation**: No dependency on main FESOM2 build
- **Minimal source footprint**: Only 97 essential FESOM2 source files
- **Library reuse**: Other programs can link against `libfesom_toy.so`
- **Clean separation**: Toy model has its own CMake namespace and install targets
- **Easy maintenance**: Changes to toy model don't affect main FESOM2 build
- **Development flexibility**: Ideal for algorithm testing and debugging

## Running the Toy Model

### Configuration Requirements

The toy model uses the same namelist structure as the full FESOM2 but with simplified settings.

#### Example: Soufflet Channel Configuration

Use a configuration like `config/namelist.config_toy_soufflet` with:

```fortran
&modelname
toy_ocean=.true.
which_toy='soufflet'
/

&paths
mesh_path='path/to/mesh/'
/

&timestep
step_per_day=72
run_length=365
run_length_unit='d'
/

&oce_dyn
use_ice=.false.
cartesian=.true.
fplane=.false.
cyclic_length=4.5
/
```

**Key settings:**
- `toy_ocean=.true.` - Enables toy ocean mode
- `which_toy='soufflet'` - Selects the Soufflet test case (or 'dbgyre' for double gyre)
- `use_ice=.false.` - Disables ice model
- Additional toy-specific parameters in the ocean namelist

### Execution

#### Setup Work Directory

1. **Copy configuration files** to the work directory:

   ```bash
   cd dwarf/toy_ocean/work
   # Copy namelist files from FESOM2 config directory
   cp ../../../config/namelist.config_toy_soufflet namelist.config
   cp ../../../config/namelist.oce .
   cp ../../../config/namelist.ice .
   # Adjust paths in namelist.config as needed
   ```

2. **Run the toy model**:

   For standalone build from `dwarf/toy_ocean/`:
   ```bash
   cd work
   mpirun -np 8 ../build/src/fesom_toy
   ```

   For integrated build from FESOM2 root:
   ```bash
   cd work_directory
   mpirun -n 4 /path/to/fesom_toy
   ```

   Or with a job scheduler:
   ```bash
   srun -n 96 /path/to/fesom_toy
   ```

#### Verified Test Case

✅ **Soufflet Channel Test** (Successfully tested with 8 MPI ranks):
```bash
cd dwarf/toy_ocean/work
# Copy Soufflet configuration from config/
mpirun -np 8 ../build/src/fesom_toy
```

The model runs successfully with the Soufflet channel configuration from `config/namelist.config_toy_soufflet`.

### Expected Behavior

The toy model will:
1. Initialize with analytical initial conditions (defined in `toy_channel_soufflet.F90` or similar)
2. Run ocean-only time-stepping
3. Apply toy-specific forcing (e.g., zonal restoring for Soufflet)
4. Output diagnostics and fields (same format as full FESOM2)
5. Report performance statistics at the end

## Code Structure

### Initialization (`toy_fesom_init`)

1. MPI initialization
2. Configuration reading (namelists)
3. Clock setup
4. Mesh loading and partitioning
5. Dynamics and tracer array allocation
6. **Ocean setup** - calls toy-specific initial conditions:
   - For Soufflet: `initial_state_soufflet()` sets idealized stratification
   - Computes geostrophically balanced initial velocity
7. Dummy ice structure initialization
8. Restart handling (if applicable)
9. Zonal mean computation (for Soufflet restoring)

### Time-Stepping (`toy_fesom_runloop`)

For each timestep:
1. Update clock
2. Compute velocity at nodes
3. Pre-ocean step preparation
4. **Core ocean timestep** (`oce_timestep_ale`) - THE MAIN COMPUTATION
   - Mixing coefficients
   - Pressure gradient
   - Velocity RHS (Coriolis, advection, diffusion)
   - SSH solution
   - Tracer advection-diffusion
   - Toy-specific restoring (Soufflet: zonal relaxation)
5. Diagnostics
6. Output writing
7. Restart writing

### Finalization (`toy_fesom_finalize`)

1. Finalize I/O
2. Gather runtime statistics across MPI ranks
3. Print performance report
4. MPI finalization

## Performance Considerations

The toy model is expected to run **faster** than the full FESOM2 for the same configuration because:
- No ice model time-stepping
- No ice-ocean coupling
- No atmospheric forcing updates
- Simpler forcing (analytical or zonal restoring)
- No biogeochemistry or additional tracers

Typical speedup: **10-30%** depending on the original ice fraction and forcing complexity.

## Development and Testing

### Adding New Toy Cases

To add a new idealized test case:

1. Create a new module (e.g., `toy_channel_newcase.F90`) with:
   - `initial_state_newcase()` subroutine
   - Any case-specific forcing or restoring functions

2. Add the module to the appropriate USE statements in `fesom_toy_module.F90`

3. Add a case in `ocean_setup()` (in `oce_setup_step.F90`):
   ```fortran
   CASE("newcase")
       call initial_state_newcase(dynamics, tracers, partit, mesh)
   ```

4. Create corresponding namelist files

### Debugging

The toy model supports the same debugging flags as the full FESOM2:
- Set `flag_debug=.true.` in namelist for verbose output
- Use `-DCMAKE_BUILD_TYPE=Debug` for compiler debug symbols
- Enable profiling with `-DFESOM_PROFILING=ON`

## Limitations

The toy model is **not suitable** for:
- Realistic simulations requiring ice-ocean interactions
- Coupled climate simulations
- Biogeochemical studies
- Studies requiring atmospheric forcing from reanalysis
- Tidal dynamics
- Ice shelf cavity processes

For these applications, use the full FESOM2 model.

## References

- **Soufflet et al. (2016)**: "On effective resolution in ocean models"
  *Ocean Modelling*, 98, 36-50. doi:10.1016/j.ocemod.2015.12.004

- **FESOM2 Documentation**: https://fesom2.readthedocs.io/

## Support

For questions or issues specific to the toy model:
1. Check that your configuration has `toy_ocean=.true.` and `use_ice=.false.`
2. Verify the toy case module exists (e.g., `toy_channel_soufflet.F90`)
3. Consult the main FESOM2 documentation for general issues
4. Report bugs via the FESOM2 GitHub repository

---

**Last Updated**: 2025-11-04
**FESOM2 Version**: 2.0
**Toy Model Version**: 1.0
