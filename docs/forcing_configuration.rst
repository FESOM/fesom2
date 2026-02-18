.. _chap_forcing_configuration:

Forcing configuration (namelist.forcing)
****************************************

``namelist.forcing`` describes how atmospheric boundary conditions are read and processed. It controls the bulk formulae, optional land-ice freshwater, the age tracer mask, and all file names/variables for surface forcing, runoff and SSS restoring.

Sections of the namelist
========================

Section &forcing_exchange_coeff
"""""""""""""""""""""""""""""""

- **Ce_atm_oce=1.75e-3**, **Ch_atm_oce=1.75e-3**, **Cd_atm_oce=1.0e-3** turbulent exchange coefficients over open water for latent heat, sensible heat and momentum.
- **Ce_atm_ice=1.75e-3**, **Ch_atm_ice=1.75e-3**, **Cd_atm_ice=1.2e-3** analogous coefficients over sea ice.
- **Swind=0.0** couples ocean surface currents back to the wind stress following Renault et al. (2019); a value of ``1`` removes ocean-current influence entirely, values between 0 and 1 partially reduce the applied stress.

Section &forcing_bulk
"""""""""""""""""""""

- **AOMIP_drag_coeff=.false.** use AOMIP drag coefficients instead of the default NCAR set.
- **ncar_bulk_formulae=.true.** enable the Large & Yeager (2004/2009) NCAR bulk formulae. Set to ``.false.`` to bypass bulk flux calculations when direct fluxes are supplied.
- **ncar_bulk_z_wind=10.0**, **ncar_bulk_z_tair=10.0**, **ncar_bulk_z_shum=10.0** reference heights (m) for wind, air temperature and humidity that must match the forcing dataset (e.g. 10 m for CORE2/JRA55-do, 2 m for NCEP/JRA55).

Section &land_ice
"""""""""""""""""

- **use_landice_water=.false.** add freshwater from land ice (routed to the ocean surface and available for output when ``use_landice_water`` is true).
- **landice_start_mon=5**, **landice_end_mon=10** active months for land-ice freshwater input.
- **fwf_path='<add path>'** directory containing the land-ice runoff files used when the feature is enabled.

Section &age_tracer
"""""""""""""""""""

- **use_age_tracer=.false.** add a passive water-mass age tracer to the model state.
- **use_age_mask=.false.**, **age_tracer_path='<add path>'** optional mask file for age tracer initialization (e.g. to restrict the tracer to certain basins).
- **age_start_year=2000** reference year when the age tracer is set to zero; useful when replaying historical forcing.

Section &nam_sbc
""""""""""""""""

Forcing file names are provided without the year suffix; the code appends ``.<year>.nc`` when opening files. All paths are resolved relative to the working directory unless absolute paths are given.

- **nm_xwind_file**, **nm_ywind_file** 10 m wind speed components; **nm_xstre_file**, **nm_ystre_file** optionally provide wind stress directly instead of wind speed.
- **nm_humi_file**, **nm_qsr_file**, **nm_qlw_file**, **nm_tair_file**, **nm_prec_file**, **nm_snow_file**, **nm_mslp_file**, **nm_cloud_file** humidity, shortwave, longwave, 2 m air temperature, total precipitation, snowfall, mean sea level pressure and cloud cover forcing files.
- **nm_*_var** variable names inside each file (``nm_xwind_var``, ``nm_ywind_var``, ``nm_xstre_var``, ``nm_ystre_var``, ``nm_humi_var``, ``nm_qsr_var``, ``nm_qlw_var``, ``nm_tair_var``, ``nm_prec_var``, ``nm_snow_var``, **nm_mslp_var**, **nm_cloud_var**).
- **nm_nc_iyear=1900**, **nm_nc_imm=1**, **nm_nc_idd=1**, **nm_nc_freq=1**, **nm_nc_tmid=0** describe the time axis in the forcing files: start year/month/day, number of records per day (1 for daily, 4 for 6-hourly, etc.), and whether timestamps represent mid-interval values (CORE: 1, JRA55: 0).
- **y_perpetual=.false.** loop a single forcing year indefinitely (set ``nm_nc_iyear`` and provide the matching annual files).
- **l_xwind=.true.**, **l_ywind=.true.**, **l_xstre=.false.**, **l_ystre=.false.**, **l_humi=.true.**, **l_qsr=.true.**, **l_qlw=.true.**, **l_tair=.true.**, **l_prec=.true.**, **l_mslp=.false.**, **l_cloud=.false.**, **l_snow=.true.** switches for each forcing component. When ``l_xstre/l_ystre`` are ``.true.`` the stress fields are used instead of computing stress from wind speed.
- **nm_runoff_file='<add path>'**, **runoff_data_source='JRA55'**, **runoff_climatology=.false.** runoff source and file. Supported sources are ``'CORE1'``/``'CORE2'`` (single-field climatology), ``'JRA55'``/``'Dai09'`` (monthly climatology). Set ``runoff_climatology=.true.`` to treat the file as a repeating 12-month climatology; ``.false.`` expects a transient monthly time series.
- **nm_sss_data_file='<add path>'**, **sss_data_source='CORE2'** sea-surface salinity restoring data and source (``'CORE2'``, ``'WOA'``, ``'NONE'``). Only used when salinity restoring is enabled in ``namelist.oce``.
- **chl_data_source='Sweeney'**, **nm_chl_data_file='<add path>'**, **chl_const=0.1** chlorophyll data for shortwave penetration (active when ``use_sw_pene`` in ``namelist.config`` is ``.true.``). Choose ``'Sweeney'`` for the supplied climatology file or ``'None'`` to use the constant value ``chl_const``.
- **use_runoff_mapper=.false.**, **runoff_basins_file='<add path>'**, **runoff_radius=500000.** enable and configure the runoff mapper that spreads river discharge over neighbouring coastal cells using the provided basin map and radius.
