 REcoM Biogeochemistry CI Test

This test validates FESOM2 with REcoM (Regulated Ecosystem Model) biogeochemistry module.

## Test Configuration

- **Mesh**: PI mesh (test_global) - ~3k nodes
- **Duration**: 1 day
- **Time steps**: 32 per day
- **Processes**: 2 MPI tasks
- **Features**: 
  - REcoM biogeochemistry (30 tracers)
  - 3 zooplankton, 2 detritus groups
  - Iron chemistry with dust input
  - CO2 system (without calcification for this test)

## Required Input Files

All files are located in `test/input/recom/`:

### REcoM-specific files:
- `DustClimMonthlyAlbani_pimesh.nc` - Dust deposition climatology
- `MonthlyAtmCO2_gcb2024.nc` - Atmospheric CO2 concentration
- `GLODAPv2.2016b.PI_TCO2_fesom2_mmol_fix_z_Fillvalue.nc` - Pre-industrial DIC (optional)

### Initial conditions for biogeochemical tracers:
- `fe5deg.nc` - Iron (Fe)
- `oxy5deg.nc` - Oxygen
- `si5deg.nc` - Silicate
- `talk5deg.nc` - Total alkalinity
- `tco2_5deg.nc` - Total dissolved inorganic carbon
- `din5deg.nc` - Dissolved inorganic nitrogen
- `woa18_netcdf_5deg.nc` - Temperature and salinity

## Build Instructions

Compile FESOM2 with REcoM support:

```bash
./configure.sh ubuntu -DRECOM_COUPLED=ON
```

## Running the Test

### Local Docker Test:
```bash
mkrun recom test_pi_recom -m docker
cd work_recom
./job_docker_new
```

### CI Test (GitHub Actions):
The test runs automatically via `.github/workflows/fesom2_recom.yml` on pull requests to main branch.

## Expected Outputs

The test outputs daily averages for:
- Physical variables: SST, temperature, salinity, velocities, sea ice
- Biogeochemical variables: DIN (nitrogen), DIC (carbon)

## Validation

The `fcheck` values in `setup.yml` validate that:
1. Physical fields match baseline test_pi values
2. Biogeochemical tracers are within expected ranges:
   - DIN: ~15 mmol/m³ (global mean)
   - DIC: ~2100 mmol/m³ (global mean)

## Notes

- This test uses `REcoM_restart=False` for initialization from climatology
- Cavity is disabled (`use_cavity=False`)
- Uses standard PI mesh without cavity features
- Test follows the same pattern as icebergs and cavity tests
