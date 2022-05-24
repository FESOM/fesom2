.. _chap_forcing_configuration:

Forcing configuration (namelist.forcing)
****************************************

Sections of the namelist
========================

Section &forcing_exchange_coeff
"""""""""""""""""""""""""""""""

- **Ce_atm_oce=1.75e-3** Exchange coeff. of latent heat over open water.
- **Ch_atm_oce=1.75e-3** Exchange coeff. of sensible heat over open water.
- **Cd_atm_oce=1.0e-3** Drag coefficient between atmosphere and water.
- **Ce_atm_ice=1.75e-3** Exchange coeff. of latent heat over ice.
- **Ch_atm_ice=1.75e-3** Exchange coeff. of sensible heat over ice.
- **Cd_atm_ice=1.2e-3** Drag coefficient between atmosphere and ice.

Section &forcing_bulk
"""""""""""""""""""""

- **AOMIP_drag_coeff=.false.**
- **ncar_bulk_formulae=.true.**


Section &land_ice
"""""""""""""""""

**use_landice_water=.false.**
**landice_start_mon=5**
**landice_end_mon=10**

Section &nam_sbc
""""""""""""""""

Forcing file names should be in the form of ``variable.year.nc```. In the namelist you provide a full path to the file and ``variable.`` name in the form of::

    nm_xwind_file = '/path/to/forcing/CORE2/u_10.'

- **nm_xwind_file=''** Name of the file with winds.
- **nm_ywind_file=''** Name of the file with winds.
- **nm_humi_file=''** Name of the file with humidity.
- **nm_qsr_file=''** Name of the file with solar heat.
- **nm_qlw_file=''** Name of the file with Long wave.
- **nm_tair_file=''** Name of the file with 2m air temperature.
- **nm_prec_file=''** Name of the file with total precipitation.
- **nm_snow_file=''** Name of the file with snow  precipitation.
- **nm_mslp_file=''** Name of the file with air pressure at sea level.
- **nm_xwind_var=''** Name of the variable in file with wind.
- **nm_ywind_var=''** Name of the variable in file with wind.
- **nm_humi_var=''** Name of the variable in file with humidity.
- **nm_qsr_var=''** Name of the variable in file with solar heat.
- **nm_qlw_var=''** Name of the variable in file with Long wave.
- **nm_tair_var=''** Name of the variable in file with 2m air temperature.
- **nm_prec_var=''** Name of the variable in file with total precipitation.
- **nm_snow_var=''** Name of the variable in file with total precipitation.
- **nm_mslp_var=''** Name of the variable in file with air pressure at sea level.
- **nm_nc_iyear=1948** First year of the forcing.
- **nm_nc_imm=1** Initial month of time axis in netCDF.
- **nm_nc_idd=1** Initial day of time axis in netCDF.
- **nm_nc_freq=1** Data points per day (i.e. 86400 if the time axis is in seconds)
- **nm_nc_tmid=1** It's 1 if the time stamps are given at the mid points of the netcdf file, 0 otherwise (i.e. 1 in CORE1, CORE2; 0 in JRA55).

The following options control if the forcing files for particular variables are used or not.

- **l_xwind=.true.**
- **l_ywind=.true.**
- **l_humi=.true.**
- **l_qsr=.true.**
- **l_qlw=.true.**
- **l_tair=.true.**
- **l_prec=.true.**
- **l_mslp=.false.**
- **l_cloud=.false.**
- **l_snow=.true.**


- **nm_runoff_file=''** Name of the runoff file.
- **runoff_data_source ='CORE2'** Other options are ``Dai09``, ``CORE2``
- **nm_sss_data_file=''** Name of the sea surface salinity restoring data file.
- **sss_data_source='CORE2'**