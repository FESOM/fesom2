! ============================================================================
! ============ Namelist file for FESOM2 output configuration =================
! ============================================================================
! This file contains configuration for model output and diagnostics:
! - Diagnostic flags for optional output fields
! - General output settings (compression, rotation)
! - Output variable list with frequency and precision
! - Complete catalog of all available output fields
!
! See the output catalog at the end of this file for all possible variables.
! Some outputs require specific flags in &diag_list or other namelists.
! ============================================================================

! ============================================================================
! DIAGNOSTIC FLAGS
! ============================================================================
! Enable/disable optional diagnostic computations and outputs.
! Setting these to .true. enables additional output fields (see catalog below).
! ============================================================================
&diag_list
ldiag_solver      = .false.  ! enables solver diagnostics (convergence, iterations)
lcurt_stress_surf = .false.  ! enables 'curl_surf' output (vorticity of surface stress)
ldiag_curl_vel3   = .false.  ! enables 'curl_u' output (relative vorticity from 3D velocity)
ldiag_Ri          = .false.  ! enables Richardson number diagnostics ('shear', 'Ri')
ldiag_turbflux    = .false.  ! enables turbulent flux diagnostics ('KvdTdz', 'KvdSdz')
ldiag_salt3D      = .false.  ! enables 3D salinity diagnostics
ldiag_dMOC        = .false.  ! enables 'dMOC' output (density MOC diagnostics)
ldiag_DVD         = .false.  ! enables 'DVD' output (Discrete Variance Decay diagnostics)
ldiag_forc        = .false.  ! enables 'FORC' output (comprehensive forcing diagnostics)
ldiag_extflds     = .false.  ! enables extended field diagnostics
ldiag_destine     = .false.  ! enables heat content computation ('hc300m', 'hc700m', 'hc')
ldiag_trflx       = .false.  ! enables tracer flux diagnostics ('utemp', 'vtemp', 'usalt', 'vsalt')
ldiag_uvw_sqr     = .false.  ! enables 'UVW_SQR' output (squared velocities: u2, v2, w2)
ldiag_trgrd_xyz   = .false.  ! enables 'TRGRD_XYZ' output (horizontal & vertical tracer gradients)
/

! ============================================================================
! GENERAL OUTPUT SETTINGS
! ============================================================================
&nml_general
io_listsize       = 120      ! total number of streams to allocate. Shall be larger or equal to the number of streams in &nml_list (max. 150)
vec_autorotate    = .false.  ! unrotate vector fields (velocities, winds) before writing to output files
compression_level = 1        ! compression level for netCDF output (1=fastest, 9=smallest)
/

! ============================================================================
! OUTPUT VARIABLE LIST
! ============================================================================
! Format: 'variable_id', frequency, unit, precision
!   frequency  = output frequency (integer)
!   unit       = 'y' (yearly), 'm' (monthly), 'd' (daily), 'h' (hourly), 's' (steps)
!   precision  = 4 (single precision) or 8 (double precision)
! ============================================================================
&nml_list
io_list = 'sst       ', 1, 'm', 4,
          'sss       ', 1, 'm', 4,
          'ssh       ', 1, 'm', 4,
          'uice      ', 1, 'd', 4,
          'vice      ', 1, 'd', 4,
          'a_ice     ', 1, 'm', 4,
          'm_ice     ', 1, 'm', 4,
          'm_snow    ', 1, 'm', 4,
          'MLD1      ', 1, 'm', 4,
          'MLD2      ', 1, 'm', 4,
          'MLD3      ', 1, 'm', 4,
          'tx_sur    ', 1, 'm', 4,
          'ty_sur    ', 1, 'm', 4,
          'temp      ', 1, 'm', 4,
          'salt      ', 1, 'm', 8,
          'N2        ', 1, 'y', 4,
          'Kv        ', 1, 'y', 4,
          'u         ', 1, 'y', 4,
          'v         ', 1, 'y', 4,
          'unod      ', 1, 'y', 4,
          'vnod      ', 1, 'y', 4,
          'w         ', 1, 'y', 4,
          'Av        ', 1, 'y', 4,
          'bolus_u   ', 1, 'y', 4,
          'bolus_v   ', 1, 'y', 4,
          'bolus_w   ', 1, 'y', 4,
          'fw        ', 1, 'm', 4,
          'fh        ', 1, 'm', 4,
          'otracers  ', 1, 'y', 4,
/

! ============================================================================
! COMPLETE CATALOG OF ALL POSSIBLE OUTPUT FIELDS
! ============================================================================
! Below is a comprehensive list of all valid io_list IDs available in FESOM2.
! To enable any field, copy the line to the &nml_list section above.
! NOTE: Some fields require specific flags to be enabled (see comments).
! ============================================================================

! --- 2D OCEAN SURFACE FIELDS ---
! 'sst       ',1, 'm', 4,  ! sea surface temperature [C]
! 'sss       ',1, 'm', 4,  ! sea surface salinity [psu]
! 'ssh       ',1, 'm', 4,  ! sea surface elevation [m]
! 'vve_5     ',1, 'm', 4,  ! vertical velocity at 5th level [m/s]
! 't_star    ',1, 'm', 4,  ! air temperature [C]
! 'qsr       ',1, 'm', 4,  ! solar radiation [W/s^2]

! --- 3D OCEAN FIELDS ---
! 'temp      ',1, 'm', 4,  ! temperature [C]
! 'salt      ',1, 'm', 8,  ! salinity [psu]
! 'sigma0    ',1, 'm', 4,  ! potential density [kg/m3]
! 'u         ',1, 'y', 4,  ! zonal velocity [m/s]
! 'v         ',1, 'y', 4,  ! meridional velocity [m/s]
! 'unod      ',1, 'y', 4,  ! zonal velocity at nodes [m/s]
! 'vnod      ',1, 'y', 4,  ! meridional velocity at nodes [m/s]
! 'w         ',1, 'y', 4,  ! vertical velocity [m/s]
! 'otracers  ',1, 'y', 4,  ! all other tracers if applicable
! 'age       ',1, 'm', 4,  ! water age tracer [year] (require use_age_tracer=.true.)

! --- 2D SSH DIAGNOSTIC VARIABLES ---
! 'ssh_rhs   ',1, 'm', 4,  ! ssh rhs [m/s]
! 'ssh_rhs_old',1, 'm', 4, ! ssh rhs old [m/s]
! 'd_eta     ',1, 'm', 4,  ! dssh from solver [m]
! 'hbar      ',1, 'm', 4,  ! ssh n+0.5 tstep [m]
! 'hbar_old  ',1, 'm', 4,  ! ssh n-0.5 tstep [m]
! 'dhe       ',1, 'm', 4,  ! dhbar @ elem [m]

! --- SEA ICE FIELDS (require use_ice=.true.) ---
! 'uice      ',1, 'd', 4,  ! ice velocity x [m/s]
! 'vice      ',1, 'd', 4,  ! ice velocity y [m/s]
! 'a_ice     ',1, 'm', 4,  ! ice concentration [%]
! 'm_ice     ',1, 'm', 4,  ! ice height per unit area [m]
! 'thdgr     ',1, 'm', 4,  ! thermodynamic growth rate ice [m/s]
! 'thdgrsn   ',1, 'm', 4,  ! thermodynamic growth rate snow [m/s]
! 'flice     ',1, 'm', 4,  ! flooding growth rate ice [m/s]
! 'm_snow    ',1, 'm', 4,  ! snow height per unit area [m]
! 'h_ice     ',1, 'm', 4,  ! ice thickness over ice-covered fraction [m]
! 'h_snow    ',1, 'm', 4,  ! snow thickness over ice-covered fraction [m]

! --- SEA ICE DEBUG VARIABLES (require use_ice=.true.) ---
! 'strength_ice',1, 'm', 4,  ! ice strength [?]
! 'inv_areamass',1, 'm', 4,  ! inv_areamass [?]
! 'rhs_a     ',1, 'm', 4,  ! rhs_a [?]
! 'rhs_m     ',1, 'm', 4,  ! rhs_m [?]
! 'sgm11     ',1, 'm', 4,  ! sgm11 [?]
! 'sgm12     ',1, 'm', 4,  ! sgm12 [?]
! 'sgm22     ',1, 'm', 4,  ! sgm22 [?]
! 'eps11     ',1, 'm', 4,  ! eps11 [?]
! 'eps12     ',1, 'm', 4,  ! eps12 [?]
! 'eps22     ',1, 'm', 4,  ! eps22 [?]
! 'u_rhs_ice ',1, 'm', 4,  ! u_rhs_ice [?]
! 'v_rhs_ice ',1, 'm', 4,  ! v_rhs_ice [?]
! 'metric_fac',1, 'm', 4,  ! metric_fac [?]
! 'elevat_ice',1, 'm', 4,  ! elevat_ice [?]
! 'uwice     ',1, 'm', 4,  ! uwice [?]
! 'vwice     ',1, 'm', 4,  ! vwice [?]
! 'twice     ',1, 'm', 4,  ! twice [?]
! 'swice     ',1, 'm', 4,  ! swice [?]

! --- MIXED LAYER DEPTH ---
! 'MLD1      ',1, 'm', 4,  ! Mixed Layer Depth [m]
! 'MLD2      ',1, 'm', 4,  ! Mixed Layer Depth [m]
! 'MLD3      ',1, 'm', 4,  ! Mixed Layer Depth [m]

! --- HEAT CONTENT (require ldiag_destine=.true.) ---
! 'hc300m    ',1, 'm', 4,  ! Vertically integrated heat content upper 300m [J m**-2]
! 'hc700m    ',1, 'm', 4,  ! Vertically integrated heat content upper 700m [J m**-2]
! 'hc        ',1, 'm', 4,  ! Vertically integrated heat content total column [J m**-2]

! --- WATER ISOTOPES IN SEA ICE (require lwiso=.true.) ---
! 'h2o18_ice ',1, 'm', 4,  ! h2o18 concentration in sea ice [kmol/m**3]
! 'hDo16_ice ',1, 'm', 4,  ! hDo16 concentration in sea ice [kmol/m**3]
! 'h2o16_ice ',1, 'm', 4,  ! h2o16 concentration in sea ice [kmol/m**3]

! --- FRESHWATER FLUX (require use_landice_water=.true.) ---
! 'landice   ',1, 'm', 4,  ! freshwater flux [m/s]

! --- SURFACE FORCING ---
! 'tx_sur    ',1, 'm', 4,  ! zonal wind str. to ocean [N/m2]
! 'ty_sur    ',1, 'm', 4,  ! meridional wind str. to ocean [N/m2]
! 'curl_surf ',1, 'm', 4,  ! vorticity of the surface stress [none] (require lcurt_stress_surf=.true.)
! 'fh        ',1, 'm', 4,  ! heat flux [W/m2]
! 'fw        ',1, 'm', 4,  ! fresh water flux [m/s]
! 'atmice_x  ',1, 'm', 4,  ! stress atmice x [N/m2]
! 'atmice_y  ',1, 'm', 4,  ! stress atmice y [N/m2]
! 'atmoce_x  ',1, 'm', 4,  ! stress atmoce x [N/m2]
! 'atmoce_y  ',1, 'm', 4,  ! stress atmoce y [N/m2]
! 'iceoce_x  ',1, 'm', 4,  ! stress iceoce x [N/m2]
! 'iceoce_y  ',1, 'm', 4,  ! stress iceoce y [N/m2]
! 'alpha     ',1, 'm', 4,  ! thermal expansion [none]
! 'beta      ',1, 'm', 4,  ! saline contraction [none]
! 'dens_flux ',1, 'm', 4,  ! density flux [kg/(m3*s)]
! 'runoff    ',1, 'm', 4,  ! river runoff [m/s]
! 'evap      ',1, 'm', 4,  ! evaporation [m/s]
! 'prec      ',1, 'm', 4,  ! precipitation rain [m/s]
! 'snow      ',1, 'm', 4,  ! precipitation snow [m/s]
! 'tair      ',1, 'm', 4,  ! surface air temperature [°C]
! 'shum      ',1, 'm', 4,  ! specific humidity []
! 'swr       ',1, 'm', 4,  ! short wave radiation [W/m^2]
! 'lwr       ',1, 'm', 4,  ! long wave radiation [W/m^2]
! 'uwind     ',1, 'm', 4,  ! 10m zonal surface wind velocity [m/s]
! 'vwind     ',1, 'm', 4,  ! 10m merid. surface wind velocity [m/s]
! 'virtsalt  ',1, 'm', 4,  ! virtual salt flux [m/s*psu]
! 'relaxsalt ',1, 'm', 4,  ! relaxation salt flux [m/s*psu]
! 'realsalt  ',1, 'm', 4,  ! real salt flux from sea ice [m/s*psu]

! --- KPP VERTICAL MIXING (require mix_scheme_nmb==1,17,3,37) ---
! 'kpp_obldepth',1, 'm', 4, ! KPP ocean boundary layer depth [m]
! 'kpp_sbuoyflx',1, 'm', 4, ! surface buoyancy flux [m2/s3]

! --- RECOM 2D BIOGEOCHEMISTRY (require use_REcoM=.true. and __recom) ---
! 'dpCO2s    ',1, 'm', 4,  ! Difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
! 'pCO2s     ',1, 'm', 4,  ! Partial pressure of oceanic CO2 [uatm]
! 'CO2f      ',1, 'm', 4,  ! CO2-flux into the surface water [mmolC/m2/d]
! 'O2f       ',1, 'm', 4,  ! O2-flux into the surface water [mmolO/m2/d]
! 'Hp        ',1, 'm', 4,  ! Mean of H-plus ions in the surface water [mol/kg]
! 'aFe       ',1, 'm', 4,  ! Atmospheric iron input [umolFe/m2/s]
! 'aN        ',1, 'm', 4,  ! Atmospheric DIN input [mmolN/m2/s]
! 'benN      ',1, 'm', 4,  ! Benthos Nitrogen [mmol]
! 'benC      ',1, 'm', 4,  ! Benthos Carbon [mmol]
! 'benSi     ',1, 'm', 4,  ! Benthos silicon [mmol]
! 'benCalc   ',1, 'm', 4,  ! Benthos calcite [mmol]
! 'NPPn      ',1, 'm', 4,  ! Mean NPP nanophytoplankton [mmolC/m2/d]
! 'NPPd      ',1, 'm', 4,  ! Mean NPP diatoms [mmolC/m2/d]
! 'GPPn      ',1, 'm', 4,  ! Mean GPP nanophytoplankton [mmolC/m2/d]
! 'GPPd      ',1, 'm', 4,  ! Mean GPP diatoms [mmolC/m2/d]
! 'NNAn      ',1, 'm', 4,  ! Net N-assimilation nanophytoplankton [mmolN/m2/d]
! 'NNAd      ',1, 'm', 4,  ! Net N-assimilation diatoms [mmolN/m2/d]
! 'Chldegn   ',1, 'm', 4,  ! Chlorophyll degradation nanophytoplankton [1/d]
! 'Chldegd   ',1, 'm', 4,  ! Chlorophyll degradation diatoms [1/d]
! 'NPPc      ',1, 'm', 4,  ! Mean NPP coccolithophores [mmolC/(m2*d)]
! 'GPPc      ',1, 'm', 4,  ! Mean GPP coccolithophores [mmolC/m2/d]
! 'NNAc      ',1, 'm', 4,  ! Net N-assimilation coccolithophores [mmolN/(m2*d)]
! 'Chldegc   ',1, 'm', 4,  ! Chlorophyll degradation coccolithophores [1/d]

! --- RECOM 3D BIOGEOCHEMISTRY (require use_REcoM=.true. and __recom) ---
! 'PAR       ',1, 'm', 4,  ! PAR [W/m2]
! 'respmeso  ',1, 'm', 4,  ! Respiration rate of mesozooplankton [mmolC/m2/d]
! 'respmacro ',1, 'm', 4,  ! Respiration rate of macrozooplankton [mmolC/m2/d]
! 'respmicro ',1, 'm', 4,  ! Respiration rate of microzooplankton [mmolC/m2/d]
! 'calcdiss  ',1, 'm', 4,  ! Calcite dissolution [mmolC/m2/d]
! 'calcif    ',1, 'm', 4,  ! Calcification [mmolC/m2/d]
! 'aggn      ',1, 'm', 4,  ! Aggregation of small phytoplankton [mmolC/m2/d]
! 'aggd      ',1, 'm', 4,  ! Aggregation of diatoms [mmolC/m2/d]
! 'aggc      ',1, 'm', 4,  ! Aggregation of coccolithophores [mmolC/m2/d]
! 'docexn    ',1, 'm', 4,  ! DOC excretion by small phytoplankton [mmolC/m2/d]
! 'docexd    ',1, 'm', 4,  ! DOC excretion by diatoms [mmolC/m2/d]
! 'docexc    ',1, 'm', 4,  ! DOC excretion by coccolithophores [mmolC/m2/d]
! 'respn     ',1, 'm', 4,  ! Respiration by small phytoplankton [mmolC/m2/d]
! 'respd     ',1, 'm', 4,  ! Respiration by diatoms [mmolC/m2/d]
! 'respc     ',1, 'm', 4,  ! Respiration by coccolithophores [mmolC/(m2*d)]
! 'NPPn3D    ',1, 'm', 4,  ! Net primary production of small phytoplankton [mmolC/m2/d]
! 'NPPd3D    ',1, 'm', 4,  ! Net primary production of diatoms [mmolC/m2/d]
! 'NPPc3D    ',1, 'm', 4,  ! Net primary production of coccolithophores [mmolC/m2/d]

! --- WATER ISOTOPES IN OCEAN (require lwiso=.true.) ---
! 'h2o18     ',1, 'm', 4,  ! h2o18 concentration [kmol/m**3]
! 'hDo16     ',1, 'm', 4,  ! hDo16 concentration [kmol/m**3]
! 'h2o16     ',1, 'm', 4,  ! h2o16 concentration [kmol/m**3]

! --- NEUTRAL SLOPES ---
! 'slopetap_x',1, 'y', 4,  ! neutral slope tapered X [none]
! 'slopetap_y',1, 'y', 4,  ! neutral slope tapered Y [none]
! 'slopetap_z',1, 'y', 4,  ! neutral slope tapered Z [none]
! 'slope_x   ',1, 'y', 4,  ! neutral slope X [none]
! 'slope_y   ',1, 'y', 4,  ! neutral slope Y [none]
! 'slope_z   ',1, 'y', 4,  ! neutral slope Z [none]

! --- MIXING AND DYNAMICS ---
! 'N2        ',1, 'y', 4,  ! brunt väisälä [1/s2]
! 'Kv        ',1, 'y', 4,  ! vertical diffusivity Kv [m2/s]
! 'Av        ',1, 'y', 4,  ! vertical viscosity Av [m2/s]

! --- VISCOSITY TENDENCIES (require dynamics%opt_visc==8) ---
! 'u_dis_tend',1, 'y', 4,  ! horizontal velocity viscosity tendency [m/s]
! 'v_dis_tend',1, 'y', 4,  ! meridional velocity viscosity tendency [m/s]
! 'u_back_tend',1, 'y', 4, ! horizontal velocity backscatter tendency [m2/s2]
! 'v_back_tend',1, 'y', 4, ! meridional velocity backscatter tendency [m2/s2]
! 'u_total_tend',1, 'y', 4,! horizontal velocity total viscosity tendency [m/s]
! 'v_total_tend',1, 'y', 4,! meridional velocity total viscosity tendency [m/s]

! --- FERRARI/GM PARAMETERISATION (require Fer_GM=.true.) ---
! 'bolus_u   ',1, 'y', 4,  ! GM bolus velocity U [m/s]
! 'bolus_v   ',1, 'y', 4,  ! GM bolus velocity V [m/s]
! 'bolus_w   ',1, 'y', 4,  ! GM bolus velocity W [m/s]
! 'fer_K     ',1, 'y', 4,  ! GM, stirring diff. [m2/s]
! 'fer_scal  ',1, 'y', 4,  ! GM surface scaling []
! 'fer_C     ',1, 'y', 4,  ! GM, depth independent speed [m/s]
! 'cfl_z     ',1, 'y', 4,  ! vertical CFL criteria [?]

! --- DENSITY MOC DIAGNOSTICS (require ldiag_dMOC=.true.) ---
! 'dMOC      ',1, 'y', 4,  ! fluxes for density MOC (multiple variables)

! --- PRESSURE GRADIENT FORCE ---
! 'pgf_x     ',1, 'y', 4,  ! zonal pressure gradient force [m/s^2]
! 'pgf_y     ',1, 'y', 4,  ! meridional pressure gradient force [m/s^2]

! --- ALE LAYER THICKNESS ---
! 'hnode     ',1, 'y', 4,  ! vertice layer thickness [m]
! 'hnode_new ',1, 'y', 4,  ! hnode_new [m]
! 'helem     ',1, 'y', 4,  ! elemental layer thickness [m]

! --- OIFS/IFS INTERFACE (require __oifs or __ifsinterface) ---
! 'alb       ',1, 'm', 4,  ! ice albedo [none]
! 'ist       ',1, 'm', 4,  ! ice surface temperature [K]
! 'qsi       ',1, 'm', 4,  ! ice heat flux [W/m^2]
! 'qso       ',1, 'm', 4,  ! oce heat flux [W/m^2]
! 'enthalpy  ',1, 'm', 4,  ! enthalpy of fusion [W/m^2]
! 'qcon      ',1, 'm', 4,  ! conductive heat flux [W/m^2]
! 'qres      ',1, 'm', 4,  ! residual heat flux [W/m^2]

! --- ICEBERG OUTPUTS (require use_icebergs=.true.) ---
! 'icb       ',1, 'm', 4,  ! iceberg outputs (multiple variables)

! --- TKE MIXING DIAGNOSTICS (require mix_scheme_nmb==5 or 56) ---
! 'TKE       ',1, 'm', 4,  ! TKE diagnostics (multiple variables)

! --- IDEMIX MIXING DIAGNOSTICS (require mod(mix_scheme_nmb,10)==6) ---
! 'IDEMIX    ',1, 'm', 4,  ! IDEMIX diagnostics (multiple variables)

! --- TIDAL MIXING DIAGNOSTICS (require mod(mix_scheme_nmb,10)==7) ---
! 'TIDAL     ',1, 'm', 4,  ! TIDAL diagnostics (multiple variables)

! --- FORCING DIAGNOSTICS (require ldiag_forc=.true.) ---
! 'FORC      ',1, 'm', 4,  ! forcing diagnostics (multiple variables)

! --- DISCRETE VARIANCE DECAY (require ldiag_DVD=.true.) ---
! 'DVD       ',1, 'y', 4,  ! DVD diagnostics (multiple variables)

! --- SPLIT-EXPLICIT SUBCYCLING (require dynamics%use_ssh_se_subcycl=.true.) ---
! 'SPLIT-EXPL',1, 'm', 4,  ! split-explicit diagnostics (multiple variables)

! --- SQUARED VELOCITIES (require ldiag_uvw_sqr=.true.) ---
! 'UVW_SQR   ',1, 'm', 4,  ! squared velocities (u2, v2, w2)

! --- TRACER GRADIENTS (require ldiag_trgrd_xyz=.true.) ---
! 'TRGRD_XYZ ',1, 'm', 4,  ! horizontal and vertical tracer gradients

! ============================================================================
! END OF CATALOG
! ============================================================================