# CVMix Build and Integration Flow in FESOM2

This document explains how the CVMix (Community Vertical Mixing) library is integrated into FESOM2, covering the complete build process, file organization, and runtime usage.

## Overview: The Two-Piece CVMix Architecture

CVMix integration in FESOM2 consists of two distinct pieces:

1. **External CVMix Library** - Generic ocean mixing library from GitHub
2. **FESOM2 CVMix Drivers** - FESOM2-specific implementations and wrappers

---

## 1. External CVMix Library (from GitHub)

### What is it?

A generic ocean mixing library providing standard parameterizations for vertical mixing in ocean models.

- **Repository**: https://github.com/CVMix/CVMix-src.git
- **Version**: v1.0.0 (configurable in `cmake/BuildCVMix.cmake`)
- **Download location**: `build/external/cvmix-src/`

### CVMix Library Modules

The external CVMix library provides 10 core modules from `src/shared/`:

| Module | Purpose |
|--------|---------|
| `cvmix_kinds_and_types.F90` | Type definitions and kind parameters |
| `cvmix_kpp.F90` | KPP (K-Profile Parameterization) boundary layer scheme |
| `cvmix_shear.F90` | Shear instability (Pacanowski-Philander scheme) |
| `cvmix_tidal.F90` | Tidal mixing (Simmons et al.) |
| `cvmix_background.F90` | Background diffusivity |
| `cvmix_convection.F90` | Convective adjustment |
| `cvmix_ddiff.F90` | Double diffusion |
| `cvmix_math.F90` | Mathematical utilities |
| `cvmix_put_get.F90` | Setter/getter utilities |
| `cvmix_utils.F90` | General utilities |

### Build Output

**Location**: `build/external/cvmix/`

- **Static library**: `lib/libcvmix.a` (~318 KB)
- **Shared library**: `lib/libcvmix.so` (~227 KB)
- **Module files**: `include/cvmix/*.mod` (10 Fortran modules)

### Build Process (ExternalProject)

**Configured in**: `cmake/BuildCVMix.cmake`

**Key steps**:
1. Creates directory structure: `build/external/cvmix/{lib,include}`
2. Uses CMake's `ExternalProject_Add` to:
   - Clone from GitHub (shallow clone)
   - Configure with **same compiler flags as FESOM2** (critical for compatibility)
   - Build CVMix library
   - Install to `build/external/cvmix/`

**Compiler flag compatibility** :

| Compiler | Flags Added | Purpose |
|----------|-------------|---------|
| Intel | `-r8 -i4` | 8-byte reals, 4-byte integers |
| GNU | `-fdefault-real-8 -fdefault-double-8` | 8-byte reals and doubles |
| Cray | `-s real64` | 8-byte reals |
| NVHPC | `-r8` | 8-byte reals |

**Why this matters**: CVMix must use the same floating-point precision as FESOM2 to avoid ABI incompatibilities and ensure correct numerical results. 
**CURRENTLY WE KEEP THIS above duplication till we fix fesom's compilation to be more modular and respect CMAKE standards**
---

## 2. FESOM2 CVMix Driver Sources

### Location

`src/cvmix_driver/` - 9 Fortran files

### File Categories

#### A) FESOM2-Specific Physics Implementations 

These implement mixing schemes **not** available in the standard CVMix library:

**`cvmix_idemix.F90`** 
- **Purpose**: Internal wave energy mixing (IDEMIX parameterization)
- **Physics**: Olbers & Eden (2013) parameterization
- **What it does**:
  - Solves prognostic equation for internal wave energy `E_iw`
  - Computes vertical/horizontal group velocities (`c0`, `v0`)
  - Computes dissipation parameter `alpha_c`
  - Calculates wave breaking dissipation `iw_diss`
- **Uses from CVMix library**: `cvmix_kinds_and_types`
- **Uses from FESOM2**: `cvmix_kinds_and_types_addon`, `cvmix_utils_addon`
- **Key subroutines**:
  - `init_idemix()` - Initialize parameters
  - `calc_idemix_v0()` - Calculate horizontal group velocity
  - `integrate_idemix()` - Main computation (tridiagonal solver)
  - `cvmix_coeffs_idemix()` - Public interface

**`cvmix_tke.F90`** 
- **Purpose**: Turbulent kinetic energy mixing
- **Physics**: Gaspar et al. (1990) TKE scheme
- **What it does**:
  - Solves prognostic equation for TKE
  - Computes mixing length (Blanke & Delecluse 1993)
  - Calculates momentum (`KappaM`) and tracer (`KappaH`) diffusivities
  - Handles Langmuir turbulence
  - Couples with IDEMIX for wave breaking dissipation
- **Uses from CVMix library**: `cvmix_kinds_and_types`
- **Uses from FESOM2**: `cvmix_kinds_and_types_addon`, `cvmix_utils_addon`
- **Key subroutines**:
  - `init_tke()` - Initialize parameters
  - `integrate_tke()` - Main computation (mixing length, diffusivities)
  - `cvmix_coeffs_tke()` - Public interface

#### B) Addon/Extension Files 

**`cvmix_kinds_and_types_addon.F90`**
- **Purpose**: Extended data types for TKE/IDEMIX
- **What it does**: Extends CVMix's standard `cvmix_data_type` with FESOM2-specific fields:
  - IDEMIX: `E_iw`, `iw_diss`, `alpha_c`
  - TKE: `tke`, `tke_diss`, `KappaM_iface`, `KappaH_iface`
  - Shear/buoyancy: `Ssqr_iface`, `Nsqr_iface`
  - Forcing: `forc_tke_surf`, `forc_iw_surface`, `forc_iw_bottom`
- **Dependencies**: Standalone type definition (no CVMix use)

**`cvmix_utils_addon.F90`**
- **Purpose**: Utility functions for TKE/IDEMIX
- **What it does**:
  - `solve_tridiag()` - Tridiagonal matrix solver (used by TKE/IDEMIX integration)
  - `cvmix_update_tke()` - Update diffusivity values
- **Uses from CVMix library**: `cvmix_kinds_and_types`

#### C) FESOM2-CVMix Interface Wrappers 

These bridge between FESOM2's mesh/data structures and CVMix algorithms:

**`gen_modules_cvmix_kpp.F90`** 
- **Purpose**: FESOM2 interface to CVMix KPP boundary layer scheme
- **Uses from CVMix library**: `cvmix_kpp` (all KPP functions)
- **What it does**: Wraps CVMix KPP for FESOM2, handles FESOM mesh/data structures

**`gen_modules_cvmix_idemix.F90`** 
- **Purpose**: FESOM2 interface to IDEMIX implementation
- **Uses**: `cvmix_idemix` (FESOM2's implementation), `cvmix_put_get`, `cvmix_kinds_and_types`
- **What it does**: Manages IDEMIX forcing, horizontal propagation, I/O

**`gen_modules_cvmix_tke.F90`**
- **Purpose**: FESOM2 interface to TKE implementation
- **Uses**: `cvmix_tke` (FESOM2's implementation), coupling with IDEMIX
- **What it does**: Manages TKE forcing, IDEMIX-TKE coupling

**`gen_modules_cvmix_tidal.F90`**
- **Purpose**: FESOM2 interface to CVMix tidal mixing
- **Uses from CVMix library**: `cvmix_tidal` (Simmons et al. scheme)
- **What it does**: Wraps CVMix tidal mixing for FESOM2

**`gen_modules_cvmix_pp.F90`**
- **Purpose**: FESOM2 interface to CVMix PP scheme
- **Uses from CVMix library**: `cvmix_shear` (Pacanowski-Philander)
- **What it does**: Wraps CVMix shear mixing for FESOM2

---

### Include Directory Setup

**Critical for module resolution**:

- CVMix .mod files location: `build/external/cvmix/include/cvmix/`
- Added via line : `target_include_directories(fesom PRIVATE ${CVMIX_INCLUDE_DIR})`
- This allows FESOM2's cvmix_driver files to compile with statements like:

```fortran
use cvmix_kinds_and_types  ! Finds cvmix_kinds_and_types.mod
use cvmix_kpp              ! Finds cvmix_kpp.mod
use cvmix_idemix           ! Finds FESOM2's cvmix_idemix.mod
```

---

## 4. Linking and Dependencies

### Dependency Chain

```
fesom.x executable
    └─→ fesom library (target)
         ├─→ sources_Fortran (all .F90 including cvmix_driver/*.F90)
         ├─→ CVMix::cvmix (imported target)
         │    └─→ libcvmix.a (actual library file)
         └─→ DEPENDS: cvmix-build (ExternalProject target)
```

### How Linking Works

**From src/CMakeLists.txt**:
- Line : `target_link_libraries(fesom PRIVATE CVMix::cvmix)`
  - `CVMix::cvmix` is an IMPORTED target
  - Points to `build/external/cvmix/lib/libcvmix.a`
  - Provides include directories automatically via `INTERFACE_INCLUDE_DIRECTORIES`

**Symbol Resolution**:
- FESOM2's cvmix_driver files call CVMix library functions
- At link time, linker resolves symbols from `libcvmix.a`
- Example flows:
  - `cvmix_idemix.F90` calls `solve_tridiag()` from `cvmix_utils_addon.F90`
  - `gen_modules_cvmix_kpp.F90` calls `cvmix_init_kpp()` from CVMix library's `cvmix_kpp` module

---

## 5. Runtime Usage Examples

### Example 1: KPP Mixing (Uses External CVMix Library)

```
┌─────────────────────────────────────────────────────────────┐
│ 1. FESOM2 Ocean Code (e.g., oce_mixing.F90)                │
│    call compute_kpp_mixing(...)                              │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. gen_modules_cvmix_kpp.F90 (FESOM2 wrapper)               │
│    - Converts FESOM2 mesh data (T, S, depth) to CVMix format│
│    - Calls external CVMix library:                           │
│      use cvmix_kpp                                           │
│      call cvmix_coeffs_kpp(CVMixParams, ...)                 │
└────────────────────┬────────────────────────────────────────┘
                     │ (calls function from libcvmix.a)
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. CVMix Library (libcvmix.a - cvmix_kpp.F90)               │
│    - Performs KPP boundary layer calculations                │
│    - Returns diffusivity coefficients (KappaM, KappaH)       │
└────────────────────┬────────────────────────────────────────┘
                     │ (returns)
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. gen_modules_cvmix_kpp.F90 (wrapper)                      │
│    - Converts CVMix output back to FESOM2 structures        │
│    - Returns to caller                                       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. FESOM2 Ocean Code                                        │
│    - Uses mixing coefficients in momentum/tracer equations  │
└─────────────────────────────────────────────────────────────┘
```

### Example 2: TKE + IDEMIX Coupling (FESOM2-Specific)

```
┌─────────────────────────────────────────────────────────────┐
│ FESOM2 Ocean Code                                            │
│    call compute_tke_idemix_mixing(...)                       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ gen_modules_cvmix_idemix.F90 (wrapper)                      │
│    use cvmix_idemix  ! FESOM2's implementation, NOT CVMix   │
│    call cvmix_coeffs_idemix(...)                             │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ cvmix_idemix.F90 (FESOM2-specific, compiled with FESOM2)   │
│    - Computes internal wave energy E_iw                      │
│    - Computes dissipation iw_diss                            │
│    - Uses cvmix_utils_addon::solve_tridiag()                │
│    - Returns dissipation parameter alpha_c                   │
└────────────────────┬────────────────────────────────────────┘
                     │ iw_diss passed to TKE
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ gen_modules_cvmix_tke.F90 (wrapper)                         │
│    use cvmix_tke  ! FESOM2's implementation                 │
│    call cvmix_coeffs_tke(..., iw_diss_in=iw_diss)           │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ cvmix_tke.F90 (FESOM2-specific, compiled with FESOM2)      │
│    - Uses iw_diss from IDEMIX as forcing term               │
│    - Computes TKE via prognostic equation                   │
│    - Computes mixing length                                 │
│    - Computes KappaM and KappaH diffusivities              │
│    - Returns to wrapper → FESOM2                            │
└─────────────────────────────────────────────────────────────┘
```

---

## Visual Summary

### Two-Piece Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                  EXTERNAL CVMix LIBRARY                     │
│  (GitHub: CVMix/CVMix-src.git v1.0.0)                      │
│                                                             │
│  Built via ExternalProject BEFORE FESOM2                   │
│  Location: build/external/cvmix-src/                       │
│                                                             │
│  SOURCE FILES (src/shared/):                               │
│    • cvmix_kinds_and_types.F90  - Type definitions        │
│    • cvmix_kpp.F90              - KPP boundary layer      │
│    • cvmix_shear.F90            - PP scheme               │
│    • cvmix_tidal.F90            - Tidal mixing            │
│    • cvmix_background.F90       - Background diffusion    │
│    • cvmix_convection.F90       - Convective adjustment   │
│    • cvmix_ddiff.F90            - Double diffusion        │
│    • cvmix_put_get.F90          - Utilities               │
│    • cvmix_math.F90             - Math utilities          │
│    • cvmix_utils.F90            - General utilities       │
│                                                             │
│  BUILD OUTPUT:                                              │
│    • libcvmix.a (318 KB)     → build/external/cvmix/lib/   │
│    • *.mod files (10 modules)→ build/external/cvmix/include│
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │ Links against
                            │ Includes modules
┌───────────────────────────┼─────────────────────────────────┐
│                  FESOM2 CVMix DRIVERS                       │
│  (FESOM2: src/cvmix_driver/)                                │
│                                                             │
│  Built WITH FESOM2 (part of fesom target)                  │
│  Compiled with same flags as CVMix library                 │
│                                                             │
│  FESOM2-SPECIFIC IMPLEMENTATIONS (1886 lines):             │
│    • cvmix_idemix.F90 (788)   - Internal wave energy      │
│    • cvmix_tke.F90 (1098)     - Turbulent kinetic energy  │
│                                                             │
│  ADDON/EXTENSIONS (440 lines):                             │
│    • cvmix_kinds_and_types_addon.F90 (245)                │
│    • cvmix_utils_addon.F90 (195)                          │
│                                                             │
│  FESOM2-CVMix INTERFACE WRAPPERS (3259 lines):            │
│    • gen_modules_cvmix_kpp.F90 (1358)                     │
│    • gen_modules_cvmix_idemix.F90 (773)                   │
│    • gen_modules_cvmix_tke.F90 (508)                      │
│    • gen_modules_cvmix_tidal.F90 (341)                    │
│    • gen_modules_cvmix_pp.F90 (279)                       │
│                                                             │
│  ROLE: Bridge FESOM2 mesh/data ↔ CVMix algorithms         │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │ Called by
┌───────────────────────────┼─────────────────────────────────┐
│                    FESOM2 MAIN CODE                         │
│  (src/oce_*.F90, etc.)                                      │
│                                                             │
│  Calls wrapper modules to get mixing coefficients          │
└─────────────────────────────────────────────────────────────┘
```

### Build and Link Flow

```
┌──────────────────┐
│ cmake ../        │
│ (configure)      │
└────────┬─────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│ BuildCVMix.cmake                           │
│  1. Download CVMix-src from GitHub         │
│  2. CMake configure with FESOM2 flags      │
│  3. Build libcvmix.a                       │
│  4. Install to build/external/cvmix/       │
│  5. Create CVMix::cvmix imported target    │
└────────┬───────────────────────────────────┘
         │ (cvmix-build target completes)
         ▼
┌────────────────────────────────────────────┐
│ Compile FESOM2 sources                     │
│  1. src/*.F90 (main FESOM2 code)           │
│  2. src/cvmix_driver/*.F90 (9 files)       │
│     - Can find CVMix .mod files            │
│     - Include dir: build/external/cvmix/   │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│ Link fesom library                         │
│  - Combines all .o files                   │
│  - Links against CVMix::cvmix              │
│    (resolves symbols from libcvmix.a)      │
└────────┬───────────────────────────────────┘
         │
         ▼
┌────────────────────────────────────────────┐
│ Create fesom.x executable                  │
│  - Links fesom_main.F90 with fesom library │
└────────────────────────────────────────────┘
```
