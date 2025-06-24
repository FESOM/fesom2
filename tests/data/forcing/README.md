# Test Forcing Data

Place your test forcing/climate files here. You can:

1. Copy minimal forcing data from existing FESOM datasets
2. Create symbolic links to existing test data:
   ```bash
   ln -s /path/to/minimal/forcing ./forcing_data
   ```
3. Use simplified/idealized forcing for testing

The integration tests will configure namelists to use this directory.
