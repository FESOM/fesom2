# Neural Network Feature Extraction - Implementation Status

## Completed

### ✅ Core Infrastructure
- `extract_nn_features()` - Full feature extraction with all 50 inputs in correct order
- `get_neural_net()` - Access to global NN module
- `apply_nn_corrections_tracer()` - Integration point (loops over nodes/levels)
- `nn_inference_for_node_layer()` - Orchestrates NN call for (node, level)
- Feature order specification (curl_u, ld_baroc1, N2, slopes, temp, velocities, all with neighbors except ld_baroc1)

### ✅ Normalization Framework
- `load_nn_normalization_params()` - Loads means/stds (cached with SAVE)
- Feature-by-feature normalization: `(raw - mean) / std`

### ✅ Call Chain Integration
- `oce_ale_tracer.F90` updated with call to `apply_nn_corrections_tracer()`
- Optional pre-computed fields support (curl_u_nn, ld_baroc, N2_nn, slope_x_nn, slope_y_nn)
- Template showing where to pass these fields

---

## Critical TODOs (Blocking Compile/Run)

### 1. 🔴 **Load Normalization Parameters**
**File**: `MOD_NEURALNET.F90` (lines ~475-520)  
**Status**: Currently placeholder (means=0, stds=1)  
**Required Action**:
```fortran
SUBROUTINE load_nn_normalization_params(means, stds)
    ! Load 44 normalization constants from:
    ! /albedo/home/rjuhrban/fesom2/src/neuralnet_params/mean/*.bin
    ! /albedo/home/rjuhrban/fesom2/src/neuralnet_params/std/*.bin
    
    ! Expected files (or concatenated single file):
    ! - mean/curl_u_nn.bin, std/curl_u_nn.bin
    ! - mean/ld_baroc1.bin, std/ld_baroc1.bin
    ! - mean/N2.bin, std/N2.bin
    ! - mean/slope_x.bin, std/slope_x.bin
    ! - mean/slope_y.bin, std/slope_y.bin
    ! - mean/temp.bin, std/temp.bin
    ! - mean/unod.bin, std/unod.bin
    ! - mean/vnod.bin, std/vnod.bin
    !
    ! If neighbors are separate files:
    ! - mean/curl_u_nn_nb.bin, etc.
    !
    ! Total: 44 values each (means and stds)
    ! Read and store in output arrays
END SUBROUTINE
```

### 2. 🔴 **Implement Mesh Neighbor Lookup**
**File**: `MOD_NEURALNET.F90` (lines ~521-535)  
**Status**: Returns zeros (placeholder)  
**Required Action**:
```fortran
SUBROUTINE get_node_neighbors(node, mesh, neighbors)
    ! Use mesh connectivity to get 6 neighbors of a node
    ! mesh%nod_in_elem2d or similar structure provides connectivity
    ! Return neighbors(1:6) with node indices, or 0 if not available
    
    ! Example pseudo-code:
    ! Loop over elements containing 'node'
    ! Extract the 3 vertices (triangle) → get 2 other neighbors
    ! Do this for all adjacent elements → collect all unique neighbors
    ! Take first 6 (or pad with 0 if fewer than 6)
    
    ! Note: Different mesh types have different connectivity structures
    ! Your double_gyre is triangular, so mesh traversal is needed
END SUBROUTINE
```

### 3. 🔴 **Output Integration**
**File**: `MOD_NEURALNET.F90` (lines ~370-390 in `apply_nn_corrections_tracer`)  
**Status**: TODO comment  
**Required Action**:
```fortran
! After CALL nn_inference_for_node_layer(..., nn_output):
! Example (MUST confirm with your NN training):

! Option A: Direct divergence
del_ttf(nz, n) = del_ttf(nz, n) + nn_output(1)

! Option B: Horizontal and vertical flux components
! flux_x = nn_output(1)
! flux_y = nn_output(2)
! dflux_x_dn = [gradient computation]
! dflux_y_dn = [gradient computation]
! del_ttf(nz, n) = del_ttf(nz, n) + (dflux_x_dn + dflux_y_dn)

! Option C: Scale if needed
! del_ttf(nz, n) = del_ttf(nz, n) + nn_output(1) * scaling_factor
```

---

## Important NOTEs (Non-Blocking But Needed)

### 4. ⚠️ **Optional Fields Availability**
**Where to pass from** (in `oce_ale_tracer.F90` solve_tracers_ale):

```fortran
! Need to extract from FESOM memory before calling apply_nn_corrections_tracer:

! From gen_modules_diag (if computed in this timestep):
! curl_u_nn => curl_vel3(:,:)

! From oce_fer_gm (only computed if Fer_GM is true):
! ld_baroc => rosb(:) - might be in dynamics% or o_ARRAYS

! From oce_ale_pressure_bv (before smoothing):
! N2_nn => bvfreq_unsmoothed(:,:)  ! Need to store unsmoothed version
! slope_x_nn => neutral_slope(1,:,:)
! slope_y_nn => neutral_slope(2,:,:)

! Then call:
if (use_GM_NN) then
    call apply_nn_corrections_tracer(dynamics, tracers, partit, mesh, tr_num, &
                                      curl_u_nn=curl_u_nn, &
                                      ld_baroc=ld_baroc, &
                                      N2_nn=N2_nn, &
                                      slope_x_nn=slope_x_nn, &
                                      slope_y_nn=slope_y_nn)
end if
```

### 5. ⚠️ **Handle use_GM_NN Conditioning**

These fields should only be computed if `use_GM_NN` is true:
- curl_u_nn: Currently computed in gen_modules_diag only if certain conditions. Add: `IF (use_GM_NN) THEN ...`
- ld_baroc (rosb): Only computed in oce_fer_gm if `Fer_GM` is true. Need to add: `IF (use_GM_NN .OR. Fer_GM)`
- N2: Already computed, but need unsmoothed version stored separately

### 6. ⚠️ **Neighbor Padding Strategy**

Currently: `raw_value = means(idx)` for boundary nodes (neighbors < 1)  
This normalizes to 0.0 after `(means - means) / stds`

**Alternative**: Could use node's own value repeated, or zero-pad differently. Verify this matches your training!

---

## Testing Checklist

- [ ] Compile successfully (check for undefined variables)
- [ ] Run with small test case
- [ ] Print NN input features to verify they're reasonable (bounds, magnitudes match training data)
- [ ] Verify neighbors are fetched correctly (sanity check on neighbor indices)
- [ ] Confirm NN forward pass produces outputs in expected range
- [ ] Test both with/without NN corrections to verify no crash on baseline

---

## File Locations Reference

```
NN Parameters:
  /albedo/home/rjuhrban/fesom2/src/neuralnet_params/
  └─ mean/ (probability .bin files for each feature)
  └─ std/ (probability .bin files for each feature)

FESOM Source:
  /albedo/home/rjuhrban/fesom2/src/
  ├─ MOD_NEURALNET.F90 (main NN module - edited)
  ├─ oce_ale_tracer.F90 (tracer solver - edited)
  ├─ gen_modules_diag.F90 (curl_u_nn computation)
  ├─ oce_fer_gm.F90 (ld_baroc/rosb computation)
  └─ oce_ale_pressure_bv.F90 (N2/bvfreq and neutral_slope)

Training Metadata:
  /albedo/home/rjuhrban/double_gyre/ (training setup)
  └─ various .ipynb for checking NN specs
```

---

## Feature Order Summary (50 total)

```
Indices  Variable         Count  Notes
1-7      curl_u_nn(node+nb)  7      node + 6 neighbors
8        ld_baroc1        1      2D field, no neighbors
9-15     N2(node+nb)      7      node + 6 neighbors
16-22    slope_x(node+nb) 7      node + 6 neighbors
23-29    slope_y(node+nb) 7      node + 6 neighbors
30-36    temp(node+nb)    7      node + 6 neighbors
37-43    unod(node+nb)    7      node + 6 neighbors
44-50    vnod(node+nb)    7      node + 6 neighbors
                          ──
                          50 total
```

**Verification**: 7 + 1 + 7 + 7 + 7 + 7 + 7 + 7 = 50 ✓
