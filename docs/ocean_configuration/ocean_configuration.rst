.. _chap_ocean_configuration:

Ocean configuration (namelist.oce, namelist.dyn, namelist.tra, namelist.cvmix, namelist.transit)
************************************************************************************************

Ocean settings are split across several namelists. ``namelist.oce`` controls large-scale mixing choices and Gent–McWilliams/Redi parameters, ``namelist.dyn`` configures the momentum solver, ``namelist.tra`` defines tracer advection/diffusion schemes and initial conditions, ``namelist.cvmix`` exposes CVMix vertical mixing parameters, and ``namelist.transit`` toggles transient tracer packages.

Sections of the namelists
=========================

Section &oce_dyn (namelist.oce)
"""""""""""""""""""""""""""""""

- **state_equation=1** selects the equation of state: ``1`` uses the full nonlinear formulation, ``0`` forces a linear equation of state in the pressure-gradient routines.
- **C_d=0.0025**, **A_ver=1.e-4**, **scale_area=5.8e9** control background drag/viscosity. ``scale_area`` is the reference element area used when scaling viscosity/diffusivity with resolution.
- **SPP=.false.**, **SPP_dep_N=-80.0**, **SPP_dep_S=-80.0**, **SPP_drhodz_cr_N=0.01**, **SPP_drhodz_cr_S=0.01**, **SPP_expon=5** brine rejection (“salt plume”) parameterization switches and shape parameters. When enabled, rejected salt is distributed with the given maximum depths and density-gradient thresholds per hemisphere.
- **N2smth_v=.false.**, **N2smth_h=.true.**, **N2smth_hidx=1** smoothing of the buoyancy frequency used by Redi/GM; vertical smoothing is off by default, one horizontal smoothing pass is applied when enabled.
- **visc_sh_limit=5.0e-3** ceiling for shear-induced viscosity in KPP; prevents runaway mixing when shear is strong.
- **mix_scheme='KPP'**, **Ricr=0.3**, **concv=1.6** choose the vertical mixing scheme (``'KPP'``, ``'PP'``, or CVMix variants). ``Ricr`` and ``concv`` set the bulk Richardson number threshold and convection constant used by the schemes.
- **which_pgf='shchepetkin'**, **alpha=1.0**, **theta=1.0**, **use_density_ref=.false.** control the pressure-gradient formulation: Shchepetkin/Jacobian by default, with implicitness parameters ``alpha`` and ``theta``. ``use_density_ref`` switches between a constant reference density and the profile in ``density_ref_T/S`` (see ``o_param``).
- **Fer_GM=.true.**, **K_GM_max=1000.0**, **K_GM_min=2.0**, **K_GM_bvref=1**, **K_GM_resscalorder=2**, **K_GM_rampmax=-1.0**, **K_GM_rampmin=-1.0**, **K_GM_cm=3.0**, **K_GM_cmin=0.1**, **K_GM_Ktaper=.false.** configure Gent–McWilliams thickness diffusivity. ``K_GM_resscalorder`` selects resolution vs area scaling, ``*_ramp*`` taper GM in coarse regions, ``K_GM_cm/K_GM_cmin`` bound the baroclinic wave speed used in scaling, and ``K_GM_Ktaper`` applies neutral-slope tapering to the coefficient itself.
- **scaling_Ferreira=.false.**, **scaling_Rossby=.false.**, **scaling_resolution=.true.**, **scaling_FESOM14=.false.** optional spatial scalings for GM: vertical Ferreira scaling, Rossby-radius switch-off, resolution-based scaling (default), and a FESOM1.4-style near-surface taper in the Northern Hemisphere.
- **scaling_GMzexp=.true.**, **GMzexp_zref=500.0**, **GMzexp_smin=0.6** apply an exponential decay of GM/Redi coefficients with depth; ``GMzexp_smin`` caps the minimum scaling factor.
- **scaling_GINsea=.false.**, **GINsea_fac=2.0** optionally up-scale GM/Redi inside the Greenland–Iceland–Norwegian seas by a constant factor.
- **Redi=.true.**, **Redi_Ktaper=.true.**, **Redi_Kmax=0.0**, **Redi_Kmin=100.0** enable isoneutral (Redi) diffusion. ``Redi_Kmax<=0`` ties Redi diffusivity to the GM coefficient; ``Redi_Ktaper`` applies slope tapering and ``Redi_Kmin`` sets the residual diffusivity retained when tapering is active.
- **scaling_ODM95=.true.**, **ODM95_Scr=0.2e-2**, **ODM95_Sd=1.0e-3** Danabasoglu & McWilliams (1995) critical-slope tapering parameters.
- **scaling_LDD97=.false.**, **LDD97_c=2.0e0**, **LDD97_rmin=15e3**, **LDD97_rmax=100e3** surface tapering following Large et al. (1997) with Rossby-radius limits.
- **which_ALE/use_partial_cell** remain defined in ``namelist.config``; GM/Redi operate on the resulting ALE grid.
- **use_global_tides=.false.** adds equilibrium tidal potential in the surface-pressure term when enabled.

Section &tracer_phys (namelist.oce)
"""""""""""""""""""""""""""""""""""

- **diff_sh_limit=5.0e-3** shear-instability diffusivity ceiling for KPP.
- **Kv0_const=.true.**, **K_ver=1.0e-5**, **K_hor=0.0** background vertical and horizontal diffusivities; set ``Kv0_const=.false.`` to use the latitude/depth-dependent background from FESOM1.4.
- **double_diffusion=.false.** enables double-diffusive mixing within the KPP framework.
- **surf_relax_T=0.0**, **surf_relax_S=1.929e-06**, **balance_salt_water=.true.** surface restoring coefficients (m/s). When salinity restoring is active, ``balance_salt_water`` converts virtual salt to freshwater fluxes to maintain mass conservation.
- **clim_relax=0.0** 3-D tracer restoring rate (1/s). A non-zero value expects spatial masks/weights and restores the full column.
- **ref_sss_local=.true.**, **ref_sss=34.** reference salinity used in virtual-salt to freshwater conversions. With ``ref_sss_local`` the local surface salinity is used instead of the constant value.
- **use_momix=.false.**, **momix_lat=-50.0**, **momix_kv=0.01** Monin–Obukhov mixed-layer enhancement (Timmermann & Beckmann 2004). Enabled south of ``momix_lat`` with mixing coefficient ``momix_kv``.
- **use_instabmix=.true.**, **instabmix_kv=0.1** extra vertical diffusion applied under static instability.
- **use_windmix=.false.**, **windmix_kv=1.e-3**, **windmix_nl=2** wind-driven near-surface diffusion added on top of PP to stabilise strong wind events; acts on ``windmix_nl`` top layers.
- **use_kpp_nonlclflx=.false.** toggles KPP non-local transport term.

Momentum solver (namelist.dyn)
""""""""""""""""""""""""""""""

Section &dynamics_visc
----------------------

- **opt_visc=5**, **check_opt_visc=.true.** selects the horizontal viscosity/backscatter scheme. Options: ``5`` kinematic (easy) backscatter, ``6/7`` biharmonic flow-aware flavours, ``8`` dynamic backscatter. ``check_opt_visc`` aborts if option 5 is used outside its intended resolution/Rossby-radius range.
- **visc_gamma0=0.003**, **visc_gamma1=0.1**, **visc_gamma2=0.285** baseline, flow-aware and auxiliary viscosity coefficients. ``visc_gamma2`` is only used by easy/dynamic backscatter.
- **visc_easybsreturn=1.5** fraction of subgrid energy returned to the flow when ``opt_visc=5``.
- **uke_scaling=.true.**, **uke_scaling_factor=1.0**, **uke_advection=.false.**, **rosb_dis=1.0** controls for the dynamic backscatter option. ``uke_scaling`` scales dissipation with the local Rossby radius vs grid spacing (factor in the denominator); ``rosb_dis`` shifts the Rossby-radius cutoff; ``uke_advection`` adds an explicit advection step of the backscatter energy.
- **smooth_back=2**, **smooth_dis=2**, **smooth_back_tend=4** Laplacian smoothing passes applied to the backscatter fields and their tendencies (opt. 8).
- **K_back=600.0**, **c_back=0.1** magnitude and damping coefficient for the dynamic backscatter reservoir.
- **use_ivertvisc=.true.** use implicit vertical viscosity to improve stability.

Section &dynamics_general
-------------------------

- **momadv_opt=2** momentum advection choice (currently only option 2 is supported).
- **use_freeslip=.false.** lateral boundary condition; ``.false.`` uses no-slip.
- **use_wsplit=.false.**, **wsplit_maxcfl=1.0** toggles implicit/explicit splitting of vertical velocity; ``wsplit_maxcfl`` is the allowed explicit CFL number.
- **ldiag_KE=.false.** enable kinetic-energy diagnostics (extra output/compute cost).
- **AB_order=2** Adams–Bashforth time-stepping order for momentum (2 or 3).
- **use_ssh_se_subcycl=.false.**, **se_BTsteps=50**, **se_BTtheta=0.14**, **se_bottdrag=.true.**, **se_bdrag_si=.true.**, **se_visc=.true.**, **se_visc_gamma0=10**, **se_visc_gamma1=19500**, **se_visc_gamma2=0** controls for split-explicit barotropic subcycling. When enabled, ``se_*`` parameters define the implicitness, viscosity and bottom drag applied in the fast barotropic solver.

Tracer transport (namelist.tra)
"""""""""""""""""""""""""""""""

Section &tracer_listsize
------------------------

- **num_tracers=100** allocates space for tracer definitions. This must be ≥ the total number of tracers after adding optional packages (age tracers, isotopes, transient tracers).

Section &tracer_list
--------------------

- **nml_tracer_list** defines advection/diffusion schemes per tracer using tuples ``id, hor_adv, vert_adv, limiter, hor_order, vert_order``. Default entries (IDs 1 and 2) configure temperature and salinity with multidimensional FCT and QUICR. Additional tracers can be appended; IDs 0/1 are reserved for T/S internally, and some packages add IDs automatically (e.g. water isotopes, transient tracers, age tracer).

Section &tracer_init3d
----------------------

- **n_ic3d=2**, **idlist=2,1**, **filelist='phc3.0_winter.nc','phc3.0_winter.nc'**, **varlist='salt','temp'**, **t_insitu=.true.** configure initial 3-D tracer fields. ``idlist`` follows the tracer IDs above; ``filelist`` entries are resolved relative to ``ClimateDataPath``. ``t_insitu`` converts in-situ temperature to potential on the fly.

Section &tracer_init2d
----------------------

- **n_ic2d=3**, **idlist=1,2,3**, **filelist='a_ice.nc','m_ice.nc','m_snow.nc'**, **varlist='a_ice','m_ice','m_snow'**, **ini_ice_from_file=.false.** initial conditions for 2-D sea-ice fields (concentration, ice thickness, snow thickness). When ``ini_ice_from_file`` is ``.false.`` the model starts from defaults instead of reading files.

Section &tracer_general
-----------------------

- **smooth_bh_tra=.false.**, **gamma0_tra=0.0005**, **gamma1_tra=0.0125**, **gamma2_tra=0.0** biharmonic diffusion settings for tracers (filter implementation). Recommended only at very high resolution; ``gamma2_tra`` is rarely used.
- **i_vert_diff=.true.** implicit vertical diffusion for tracers.
- **AB_order=2** Adams–Bashforth order for tracer advection (2 or 3).

Vertical mixing (namelist.cvmix)
"""""""""""""""""""""""""""""""""

Section &param_tke
------------------

- **tke_c_k=0.1**, **tke_c_eps=0.7**, **tke_alpha=30.0** core TKE closure constants.
- **tke_mxl_min=1.0e-8**, **tke_mxl_choice=2** minimum and formulation for the mixing length (option 2 = Blanke & Delecluse).
- **tke_kappaM_min=0.0**, **tke_kappaM_max=100.0** viscosity bounds.
- **tke_cd=3.75** surface boundary condition constant (3.75 for Dirichlet, 1.0 for Neumann).
- **tke_surf_min=1.0e-4**, **tke_min=1.0e-6** lower limits on surface/interior TKE.
- **tke_dolangmuir=.false.** include Langmuir turbulence effects in the TKE scheme.

Section &param_idemix
---------------------

- **idemix_tau_v=172800.0**, **idemix_tau_h=1296000.0** vertical and horizontal symmetrization timescales (s) for IDEMIX internal wave energy.
- **idemix_gamma=1.570**, **idemix_jstar=5.0**, **idemix_mu0=0.33333333** spectral shape/dissipation constants.
- **idemix_sforcusage=0.2**, **idemix_n_hor_iwe_prop_iter=5** fraction of surface forcing used and number of horizontal propagation iterations.
- **idemix_surforc_file='...nc'**, **idemix_surforc_vname='var706'** wind-generated internal-wave energy source file and variable.
- **idemix_botforc_file='...nc'**, **idemix_botforc_vname='stormt_M2_plus_nycand_CnoM2'** tidal internal-wave dissipation climatology (alternative commented options are provided in the namelist file).

Section &param_pp
-----------------

- **pp_use_fesompp=.true.** choose the FESOM flavour of Pacanowski–Philander vs the original formulation.
- **pp_Av0=0.01**, **pp_alpha=5.0**, **pp_exp=2.0** PP stability curve parameters (see Pacanowski & Philander 1981).
- **pp_Avbckg=1.0e-4**, **pp_Kvbckg=1.0e-5**, **pp_use_nonconstKvb=.true.** background viscosity/diffusivity and whether to use the depth/latitude-dependent background.

Section &param_kpp
------------------

- **kpp_use_fesomkpp=.false.** choose CVMix MOM5-like (``.true.``) vs MOM6-like (``.false.``) KPP.
- **kpp_use_enhanceKv=.true.**, **kpp_use_compEkman=.true.**, **kpp_use_monob=.true.** options that cap the boundary-layer depth using Ekman/Monin–Obukhov limits and add enhanced diffusivity at the base.
- **kpp_interptype_ri='linear'**, **kpp_interptype_atobl='LMD94'**, **kpp_matchtechc='ParabolicNonLocal'** interpolation/shape choices used in boundary-layer depth determination and profile matching.
- **kpp_internalmix='KPP'**, **kpp_pp_Av0=0.01** set the below-OBL mixing method (KPP or PP) and its coefficient.
- **kpp_Av0=5.0e-3**, **kpp_Kv0=5.0e-3**, **kpp_Ri0=0.7** shear-driven mixing parameters.
- **kpp_use_nonconstKvb=.true.**, **kpp_Avbckg=1.0e-4**, **kpp_Kvbckg=1.0e-5** background viscosity/diffusivity and whether to use the non-constant profile.
- **kpp_reduce_tauuice=.false.** optionally reduces wind stress under sea ice.
- **kpp_use_StokesMOST=.false.**, **kpp_A_stokes=0.005**, **kpp_langmuir_mixing='NONE'**, **kpp_langmuir_entrainment='NONE'** wave/Stokes-drift related enhancements (set to the documented options to activate the Langmuir packages).

Section &param_tidal
--------------------

- **tidal_mixscheme='Simmons'** selects the Simmons et al. (2004) tidal mixing scheme.
- **tidal_efficiency=0.2**, **tidal_lcl_mixfrac=0.33** set the mixing efficiency and local-vs-radiated dissipation fraction.
- **tidal_vert_decayscale=500.0** e-folding scale (m) for vertical decay of tidal mixing.

Transient tracers (namelist.transit)
""""""""""""""""""""""""""""""""""""

- **l_r14c=.false.**, **l_r39ar=.false.**, **l_f11=.false.**, **l_f12=.false.**, **l_sf6=.false.** enable individual transient tracers (radiocarbon, argon-39, CFC-11/12, SF6). Enabling any tracer automatically extends the tracer list and output diagnostics.
- **anthro_transit=.false.**, **paleo_transit=.false.** choose whether the forcing time series represents the industrial era or a paleo reconstruction. Only one of them should be true.
- **length_transit=1**, **ti_start_transit=1** length of the atmospheric forcing time series and the index within the file to start from (e.g. ``length_transit=166`` for 1765–2020 industrial runs).
- **ifile_transit='Table_CO2_isoC_CFCs1112_SF6.txt'** path to the atmospheric boundary-condition file; the code reads global mean mole fractions/ratios from here.
- **r14c_a=1.0**, **r39ar_a=1.0**, **xarg_a=9.34e-3**, **xco2_a=284.32e-6** mean atmospheric ratios/mole fractions used when constant forcing is requested.
- **dic_0=2.0**, **arg_0=0.01** initial mixed-layer concentrations for DIC and argon when transient tracers are enabled.
- **decay14=3.8561e-12**, **decay39=8.1708e-11** radioactive decay constants (1/s) for 14C and 39Ar; used to age tracers during integration.
