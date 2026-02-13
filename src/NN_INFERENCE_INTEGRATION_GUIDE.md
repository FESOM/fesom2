# Neural Network Inference Integration Guide

## Overview

The neural network inference has been integrated into the FESOM2 tracer solver to add sub-grid (unresolved) flux corrections at every timestep. Here's the complete architecture:

### Call Chain

```
solve_tracers_ale (oce_ale_tracer.F90, main tracer loop per timestep)
  ├─ Loop: tr_num = 1, tracers%num_tracers
  │   ├─ Advection calculations
  │   ├─ [NEW] apply_nn_corrections_tracer()  ← Called if use_GM_NN=.TRUE.
  │   │   ├─ Loop: node = 1, myDim_nod2D (OpenMP parallelized)
  │   │   │   └─ Loop: layer = nzmin, nzmax
  │   │   │       ├─ nn_inference_for_node_layer()
  │   │   │       │   ├─ extract_nn_features()  ← [TODO: Implement based on your inputs]
  │   │   │       │   └─ forward_pass_full()
  │   │   │       │       ├─ forward_pass_single() layer 1
  │   │   │       │       ├─ forward_pass_single() layer 2
  │   │   │       │       └─ ... (loop through all layers)
  │   │   │       └─ Add nn_output to del_ttf  ← [TODO: Implement output integration]
  │   │
  │   ├─ diff_tracers_ale() (diffusion step)
  │   └─ Other tracer operations
```

## What Has Been Completed

1. ✅ **Neural Network Initialization** (`MOD_NEURALNET.F90`)
   - Type definition: `t_neural_net` with weights, biases, activations, layer sizes
   - Module-level storage: `nn_gm_module` (initialized once, shared across all ranks)
   - Initialization function: `neuralnet_init(mype, mpi_comm)` 
   - MPI broadcasting: All ranks receive initialized weights/biases
   - Access function: `get_neural_net()` returns pointer to module-level NN

2. ✅ **Forward Pass Subroutines**
   - `forward_pass_single()`: Single layer forward pass with activation
     - Computes: `output = TRANSPOSE(weights) · input + biases`
     - Applies activation (relu or id/identity)
   - `forward_pass_full()`: Full network forward pass
     - Chains all layers sequentially
     - No file I/O (data from initialized NN module)

3. ✅ **Integration Point** (`oce_ale_tracer.F90`)
   - Call location: After advection, before diffusion in tracer solver
   - Conditional on `use_GM_NN` flag
   - Position within timestep: Between advection and diffusion tendencies

## What Needs Completion

### 1. Feature Extraction: `extract_nn_features()`

**Location to implement**: Inside `nn_inference_for_node_layer()` in `MOD_NEURALNET.F90`

You need to extract and normalize inputs for your NN from FESOM data structures. Your NN expects inputs like:
- Temperature at (node, layer)
- Velocity components
- Slope magnitudes (neutral slope)
- Buoyancy frequency (stratification)
- Curl of velocity
- Neighbor information (6 neighbors in double_gyre setup)

**Example skeleton** (to be filled based on training data specs):

```fortran
SUBROUTINE extract_nn_features(node, layer, dynamics, tracers, mesh, nn_input)
    ! Extract and normalize NN input features from FESOM data
    INTEGER, INTENT(IN) :: node, layer
    TYPE(t_dyn), INTENT(IN), TARGET :: dynamics
    TYPE(t_tracer), INTENT(IN), TARGET :: tracers
    TYPE(t_mesh), INTENT(IN), TARGET :: mesh
    REAL(kind=WP), INTENT(INOUT) :: nn_input(:)
    
    ! Access FESOM data:
    ! - tracers%data(1)%values(layer, node)  [Temperature]
    ! - tracers%data(2)%values(layer, node)  [Salinity]
    ! - dynamics%uvnode(1, layer, node)      [U velocity]
    ! - dynamics%uvnode(2, layer, node)      [V velocity]
    ! - [...other fields...]
    
    ! Load normalization constants (once per simulation or cache them)
    ! - mean_{i} and std_{i} from /albedo/home/rjuhrban/fesom2/src/neuralnet_params/
    
    ! Extract local features
    ! Normalize: (raw_value - mean) / std
    ! Store in nn_input array in correct order
    
END SUBROUTINE extract_nn_features
```

**Note**: The `assign_input_values` subroutine in `MOD_NEURALNET.F90` has a partial implementation that could serve as a starting point, but it needs:
- Access to actual FESOM data pointers (dynamics, tracers)
- Proper handling of neighbors (mesh connectivity)
- Cached normalization parameters for efficiency

### 2. Output Integration: Add NN Predictions to Flux Divergence

**Location to implement**: Inside `apply_nn_corrections_tracer()` and `nn_inference_for_node_layer()` in `MOD_NEURALNET.F90`

After `forward_pass_full()` returns `nn_output`, you need to:
1. Interpret the output (what do the dimensions represent?)
2. Convert to flux divergence contribution
3. Add to `del_ttf` or other tendency arrays

**Example skeleton**:

```fortran
! Inside apply_nn_corrections_tracer, after getting nn_output:

! Option 1: If NN directly outputs flux divergence
del_ttf(nz, n) = del_ttf(nz, n) + nn_output(1)

! Option 2: If NN outputs flux components (need to compute divergence)
! vertical_flux = nn_output(1)
! horiz_flux_x = nn_output(2)
! horiz_flux_y = nn_output(3)
! divergence = gradient computation...
! del_ttf(nz, n) = del_ttf(nz, n) + divergence
```

## Parallelization Strategy

### MPI Level (No Additional Communication Needed)
- ✅ Weights/biases already broadcasted to all ranks during `neuralnet_init()`
- Each MPI rank evaluates NN independently on its local nodes
- **No inter-rank communication during inference** (efficient scaling!)

### OpenMP Level (Node-Level Parallelism)
```fortran
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmax, nzmin, nn_output)
DO n = 1, myDim_nod2D
    ! Loop over vertical layers
    DO nz = nzmin, nzmax
        ! NN inference per node/layer (embarrassingly parallel)
        CALL nn_inference_for_node_layer(n, nz, dynamics, tracers, nn, nn_output)
        ! Add output to tracer tendency
        ...
    END DO
END DO
!$OMP END PARALLEL DO
```

**Benefits**:
- Each OpenMP thread processes different nodes independently
- No data race conditions (each thread writes to different nodes)
- Near-perfect scaling up to number of nodes per MPI rank
- With 12 MPI ranks, NN evaluation is distributed: ~(total_nodes/12) per rank

### Performance Considerations

1. **Memory**: NN weights/biases ~10-100 MB per rank (fully replicated via broadcast)
2. **Compute per node/layer**: ~O(layer_sizes²) operations
   - With typical NN: 10→50→50→10 neurons → ~5000 ops per inference
   - 1000 nodes × 10 levels × 5000 ops = 50M ops per timestep (negligible)
3. **I/O**: Already optimized to one-time read at initialization
4. **Bottleneck**: Likely feature extraction from FESOM data structures (cache misses)

## Integration Checklist

- [x] Forward pass subroutines uncommented and fixed
- [x] Integration point added to solver
- [x] Module exports updated
- [ ] Feature extraction subroutine implemented
- [ ] Output integration logic implemented
- [ ] Testing on small problem
- [ ] Verification of NN predictions (sanity checks)
- [ ] Performance profiling (scaling with number of nodes)

## Key Files Modified

1. **MOD_NEURALNET.F90**
   - Uncommented `forward_pass_single()` and `forward_pass_full()`
   - Added `nn_inference_for_node_layer()` (template)
   - Added `apply_nn_corrections_tracer()` (template)
   - PUBLIC interface updated

2. **oce_ale_tracer.F90**
   - Added call to `apply_nn_corrections_tracer()` after advection

## Next Steps

1. **Implement `extract_nn_features()`**
   - Define exact input specification from training data
   - Access normalization constants from disk (once per simulation)
   - Extract 6-neighbor average where needed
   - Store in nn_input array in correct order

2. **Implement output integration**
   - Understand NN output structure (dimenisions, units, interpretation)
   - Compute flux divergence from raw output if needed
   - Apply scaling factors or corrections if needed
   - Add to `del_ttf` at correct location

3. **Test and validate**
   - Compile with `-O0 -g` for debugging
   - Check that NN produces reasonable predictions (check bounds, magnitudes)
   - Compare results with/without NN corrections
   - Profile scaling with more nodes/levels

4. **Optimize if needed**
   - Cache normalization constants instead of reading per timestep
   - Consider batching NN inference (multiple levels together) if beneficial
   - Profile to find bottlenecks

---

**Questions for refinement**:
- What is the exact structure of your NN inputs? (order, normalization, neighbor handling?)
- What are the outputs? (direct flux divergence, or components to be combined?)
- Are there specific scaling factors or corrections needed?
