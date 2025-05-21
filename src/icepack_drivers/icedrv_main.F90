!=======================================================================
!
! Module that contains the whole icepack implementation in fesom2
!
! Author: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
! Adapted for Icepack 1.4.1 by F. Kauker (frank.kauker@awi.de)
!
!=======================================================================

module icedrv_main

  use icedrv_kinds
  use icedrv_constants
  use mod_partit
  
  implicit none
  
  !=======================================================================    
  !--------- list here all public variables and 
  !--------- subroutines to be seen outside of icepack 
  !=======================================================================

  public :: &
       ! Variables
       ncat, rdg_conv_elem, rdg_shear_elem, strength,          & 
       ! Subroutines
       set_icepack, alloc_icepack, init_icepack, step_icepack, &
       icepack_to_fesom, icepack_to_fesom_single_point,        &
       init_flux_atm_ocn, ini_icepack_io , ini_mean_icepack_io
    
  !=======================================================================
  !--------- Everything else is private
  !=======================================================================

  private 

  !=======================================================================    
  !--------- Declare all variables used or required by Icepack
  !=======================================================================

  !=======================================================================
  ! 1. Setting variables used Icepack
  !=======================================================================

  integer (kind=int_kind), save  :: nx         ! number of nodes and gost nodes for each mesh partition
  integer (kind=int_kind), save  :: nx_elem    ! number of elements and gost elements for each mesh partition
  integer (kind=int_kind), save  :: nx_nh      ! number of nodes for each mesh partition (NO GOST CELLS)
  integer (kind=int_kind), save  :: nx_elem_nh ! number of elements for each mesh partition (NO GOST CELLS)
  integer (kind=int_kind), save  :: ncat       ! number of categories in use
  integer (kind=int_kind), save  :: nfsd       ! number of floe size categories in use
  integer (kind=int_kind), save  :: nilyr      ! number of ice layers per category in use
  integer (kind=int_kind), save  :: nslyr      ! number of snow layers per category in use
  integer (kind=int_kind), save  :: n_iso      ! number of water isotopes (snow, sea ice)
  integer (kind=int_kind), save  :: n_aero     ! number of aerosols in use
  integer (kind=int_kind), save  :: n_zaero    ! number of z aerosols in use
  integer (kind=int_kind), save  :: n_algae    ! number of algae in use
  integer (kind=int_kind), save  :: n_doc      ! number of DOC pools in use
  integer (kind=int_kind), save  :: n_dic      ! number of DIC pools in use
  integer (kind=int_kind), save  :: n_don      ! number of DON pools in use
  integer (kind=int_kind), save  :: n_fed      ! number of Fe pools in use dissolved Fe
  integer (kind=int_kind), save  :: n_fep      ! number of Fe pools in use particulate Fe
  integer (kind=int_kind), save  :: nblyr      ! number of bio/brine layers per category
  integer (kind=int_kind), save  :: n_bgc      ! nit, am, sil, dmspp, dmspd, dms, pon, humic
  integer (kind=int_kind), save  :: nltrcr     ! number of zbgc (includes zaero) and zsalinity tracers
  integer (kind=int_kind), save  :: max_nsw    ! number of tracers active in shortwave calculation
  integer (kind=int_kind), save  :: max_ntrcr  ! number of tracers in total
  integer (kind=int_kind), save  :: nfreq      ! number of wave frequencies ! HARDWIRED FOR NOW
  integer (kind=int_kind), save  :: ndtd       ! dynamic time steps per thermodynamic time step
  type(t_partit), pointer, save  :: p_partit   ! a pointer to the mesh partitioning (has been accessed via "use g_parsup" in previous versions)
  integer (kind=int_kind), save  :: mype       ! a copy of a mype (has been accessed via "use g_parsup" in previous versions)

  !=======================================================================
  ! 2. State variabels of Icepack
  !=======================================================================
  
  integer (kind=int_kind), allocatable, save :: & ! DIM max_ntrcr
       trcr_depend(:),                          & ! = 0 for ice area tracers, = 1 for ice volume tracers, = 2 for snow volume tracers
       n_trcr_strata(:)                           ! number of underlying tracer layers
  
  integer (kind=int_kind), allocatable, save :: & ! DIM max_ntrcr,2
       nt_strata(:,:)                             ! indices of underlying tracer layers
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx           
       aice(:),                                 & ! concentration of ice
       vice(:),                                 & ! volume per unit area of ice (m)
       vsno(:),                                 & ! volume per unit area of snow (m)
       aice0(:),                                & ! concentration of open water
       uvel(:),                                 & ! x-component of velocity (m/s) on nodes
       vvel(:),                                 & ! y-component of velocity (m/s) on nodes
       divu(:),                                 & ! strain rate I component, velocity divergence (1/s)
       shear(:),                                & ! strain rate II component (1/s)
       strength(:),                             & ! ice strength (N/m)
       aice_init(:)                               ! initial concentration of ice, diagnostics
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM max_ntrcr,3
       trcr_base(:,:)                             ! = 0 or 1 depending on tracer dependency, argument 2:  (1) aice, (2) vice, (3) vsno

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,max_ntrcr
       trcr(:,:)                                  ! ice tracers
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat
       aicen(:,:),                              & ! concentration of ice
       vicen(:,:),                              & ! volume per unit area of ice (m)
       vsnon(:,:),                              & ! volume per unit area of snow (m)
       aicen_init(:,:),                         & ! initial ice concentration, for linear ITD
       vicen_init(:,:),                         & ! initial ice volume (m), for linear ITD
       vsnon_init(:,:)                            ! initial snow volume (m), for aerosol
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,max_ntrcr,ncat
       trcrn(:,:,:)                               ! tracers

  !=======================================================================
  ! 3. Flux variabels
  !=======================================================================
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx    
       strax(:), stray(:),                      & ! wind stress (N/m^2), IN FROM ATMOS (IF .NOT. CALC_STRAIR)
       uocn(:), vocn(:),                        & ! ocean current (m/s), IN FROM OCEAN
       strairxT(:), strairyT(:),                & ! air-ice stress, T-cell (kg/m s^2), out to atmosphere
       strocnxT(:),  strocnyT(:)                  ! ice-ocean stress, T-cell (kg/m s^2), out to ocean
  
  ! diagnostic
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       strairx(:),  strairy(:),                 & ! stress on ice by air
       daidtd(:),                               & ! ice area tendency due to transport   (1/s)
       dvidtd(:),                               & ! ice volume tendency due to transport (m/s)
       dagedtd(:),                              & ! ice age tendency due to transport (s/s)
       dardg1dt(:),                             & ! rate of area loss by ridging ice (1/s)
       dardg2dt(:),                             & ! rate of area gain by new ridges (1/s)
       dvirdgdt(:),                             & ! rate of ice volume ridged (m/s)
       closing(:), opening(:),                  & ! rate of closing/obening due to divergence/shear (1/s)
       dhi_dt(:), dhs_dt(:),                    & ! ice/snow volume tendency due to thermodynamics (m/s)
       dhi_t_dt(:), dhs_t_dt(:),                & ! ice/snow volume tendency due to thermodynamics (m/s),
                                                  ! WHAT IS THE DIFFERENCE TO ABOVE (FRANK.KAUKER@AWI.DE)???
       dhi_r_dt(:), dhs_r_dt(:)                   ! ice/snow volume tendency due to ridging (m/s)

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; RIDGING DIAGNOSTICS PER CATEGORY
       dardg1ndt(:,:),                          & ! rate of area loss by ridging ice (1/s)
       dardg2ndt(:,:),                          & ! rate of area gain by new ridges (1/s)
       dvirdgndt(:,:),                          & ! rate of ice volume ridged (m/s)
       aparticn(:,:),                           & ! participation function
       krdgn(:,:),                              & ! mean ridge thickness/thickness of ridging ice
       ardgn(:,:),                              & ! fractional area of ridged ice
       vrdgn(:,:),                              & ! volume of ridged ice
       araftn(:,:),                             & ! rafting ice area
       vraftn(:,:),                             & ! rafting ice volume
       aredistn(:,:),                           & ! redistribution function: fraction of new ridge area
       vredistn(:,:)                              ! redistribution function: fraction of new ridge volume
  
  real (kind=dbl_kind), save ::                 & ! IN FROM ATMOSPHERE (IF CALC_TSFC)
       zlvl_s,                                  & ! atm level height for scalars (temperature, humidity) (m)
       zlvl_v                                     ! atm level height for vectors (wind) (m)

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       uatm(:), vatm(:),                        & ! wind velocity (m/s)
       wind(:),                                 & ! wind speed (m/s)
       potT(:),                                 & ! air potential temperature (K)
       T_air(:),                                & ! air temperature (K)
       Qa(:),                                   & ! specific humidity (kg/kg)
       rhoa(:),                                 & ! air density (kg/m^3)
       swvdr(:), swvdf(:), swidr(:), swidf(:),  & ! sw down: visible direct/visible diffuse/near IR direct/near IR diffuse (W/m^2)
       flw(:),                                  & ! incoming longwave radiation (W/m^2)
       fsw(:)                                     ! incoming shortwave radiation (W/m^2) (internal use)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; IN FROM ATMOSPHERE (IF .NOT. CALC_TSFC), THESE ARE PER ICE AREA
       fsurfn_f(:,:),                           & ! net flux to top surface, excluding fcondtop
       fcondtopn_f(:,:),                        & ! downward cond flux at top surface (W m-2)
       fsensn_f(:,:),                           & ! sensible heat flux (W m-2)
       flatn_f(:,:)                               ! latent heat flux (W m-2)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; IN FROM ATMOSPHERE
       frain(:),                                & ! rainfall rate (kg/m^2 s)
       fsnow(:),                                & ! snowfall rate (kg/m^2 s)
       fsloss(:)                                  ! rate of snow loss to leads (kg/m^2/s)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       sss(:),                                  & ! sea surface salinity (ppt)
       sst(:),                                  & ! sea surface temperature (C)
       sstdat(:),                               & ! sea surface temperature (C) saved for restoring
       frzmlt(:),                               & ! freezing/melting potential (W/m^2)
       frzmlt_init(:),                          & ! frzmlt used in current time step (W/m^2)
       Tf(:),                                   & ! freezing temperature (C)
       qdp(:),                                  & ! deep ocean heat flux (W/m^2), negative upward
       hmix(:),                                 & ! mixed layer depth (m)
       ! water isotopes
       HDO_ocn(:),                              & ! seawater concentration of HDO (kg/kg)
       H2_16O_ocn(:),                           & ! seawater concentration of H2_16O (kg/kg)
       H2_18O_ocn(:)                              ! seawater concentration of H2_18O (kg/kg)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; OUT TO ATMOSPHERE (IF CALC_TSFC), NOTE THAT TSFC IS A TRACER
       fsens(:),                                & ! sensible heat flux (W/m^2)
       flat(:),                                 & ! latent heat flux   (W/m^2)
       fswabs(:),                               & ! shortwave flux absorbed in ice and ocean (W/m^2)
       fswint_ai(:),                            & ! SW absorbed in ice interior below surface (W/m^2)
       flwout(:),                               & ! outgoing longwave radiation (W/m^2)
       Tref(:),                                 & ! 2m atm reference temperature (K)
       Qref(:),                                 & ! 2m atm reference spec humidity (kg/kg)
       Uref(:),                                 & ! 10m atm reference wind speed (m/s)
       evap(:),                                 & ! evaporative water flux (kg/m^2/s)
       evaps(:),                                & ! evaporative water flux over snow (kg/m^2/s)
       evapi(:)                                   ! evaporative water flux over ice (kg/m^2/s)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; ALBEDOS AGGREGATED OVER CATEGORIES (IF CALC_TSFC)
       alvdr(:), alidr(:), alvdf(:), alidf(:),  & ! visible direct/near IR direct/visible diffuse/near IR diffuse (fraction)
                                                  ! grid-box-mean versions:
       alvdr_ai(:), alidr_ai(:), alvdf_ai(:), alidf_ai(:), & ! visible direct/near IR direct/visible diffuse/near IR diffuse (fraction)
                                                  ! for history:
       albice(:),                               & ! bare ice albedo
       albsno(:),                               & ! snow albedo
       albpnd(:),                               & ! melt pond albedo
       apeff_ai(:),                             & ! effective pond area used for radiation calculation
       snowfrac(:),                             & ! snow fraction used in radiation
                                                  ! for diagnostic:
       alvdr_init(:), alidr_init(:), alvdf_init(:), alidf_init(:) ! spectral albedos (fraction)

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; OUT TO OCEAN
       fpond(:),                                & ! fresh water flux to ponds (kg/m^2/s)
       fresh(:),                                & ! fresh water flux to ocean (kg/m^2/s)
       fsalt(:),                                & ! salt flux to ocean (kg/m^2/s)
       fhocn(:),                                & ! net heat flux to ocean (W/m^2)
       fswthru(:),                              & ! shortwave penetrating to ocean (W/m^2)
       fswthru_vdr(:),                          & ! vis dir shortwave penetrating to ocean (W/m^2)
       fswthru_vdf(:),                          & ! vis dif shortwave penetrating to ocean (W/m^2)
       fswthru_idr(:),                          & ! nir dir shortwave penetrating to ocean (W/m^2)
       fswthru_idf(:),                          & ! nir dif shortwave penetrating to ocean (W/m^2)   
       fresh_tot(:),                            & ! total fresh water flux to ocean (kg/m^2/s)
       fhocn_tot(:)                               ! total salt flux to ocean (kg/m^2/s)
  
  real (kind=dbl_kind), allocatable, public ::  & ! DIM nx; INTERNAL
       fswfac(:),                               & ! for history
       scale_factor(:)                            ! scaling factor for shortwave components
    
  logical (kind=log_kind), public ::            &
       update_ocn_f,                            & ! if true, update fresh water and salt fluxes
       l_mpond_fresh                              ! if true, include freshwater feedback from meltponds

  real (kind=dbl_kind), allocatable, save ::    & ! DIM NX,NCAT; WHEN RUNNING IN ICE-OCEAN OR COUPLED CONFIGURATION
       meltsn(:,:),                             & ! snow melt in category n (m)
       melttn(:,:),                             & ! top melt in category n (m)
       meltbn(:,:),                             & ! bottom melt in category n (m)
       congeln(:,:),                            & ! congelation ice formation in category n (m)
       snoicen(:,:),                            & ! snow-ice formation in category n (m)
       keffn_top(:,:)                             ! effective thermal conductivity of the top ice layer on categories (W/m^2/K)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM NX; QUANTITIES PASSED FROM OCEAN MIXED LAYER TO ATMOSPHERE
       strairx_ocn(:), strairy_ocn(:),          & ! stress on ocean by air
       fsens_ocn(:),                            & ! sensible heat flux (W/m^2)
       flat_ocn(:),                             & ! latent heat flux (W/m^2)
       flwout_ocn(:),                           & ! outgoing longwave radiation (W/m^2)
       evap_ocn(:),                             & ! evaporative water flux (kg/m^2/s)
       alvdr_ocn(:), alidr_ocn(:), alvdf_ocn(:), alidf_ocn(:), & ! spectral albedos (fraction)
       Tref_ocn(:),                             & ! 2m atm reference temperature (K)
       Qref_ocn(:)                                ! 2m atm reference spec humidity (kg/kg)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx, DIAGNOSTIC
       fsurf(:),                                & ! net surface heat flux (excluding fcondtop)(W/m^2)
       fcondtop(:),                             & ! top surface conductive flux (W/m^2)
       fcondbot(:),                             & ! bottom surface conductive flux (W/m^2)
       fbot(:),                                 & ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
       Tbot(:),                                 & ! Temperature at bottom surface of ice (deg C)
       Tsnice(:),                               & ! Temperature at snow ice interface (deg C)
       congel(:),                               & ! basal ice growth         (m/step-->cm/day)
       frazil(:),                               & ! frazil ice growth        (m/step-->cm/day)
       snoice(:),                               & ! snow-ice formation       (m/step-->cm/day)
       meltt(:),                                & ! top ice melt             (m/step-->cm/day)
       melts(:),                                & ! snow melt                (m/step-->cm/day)
       meltb(:),                                & ! basal ice melt           (m/step-->cm/day)
       meltl(:),                                & ! lateral ice melt         (m/step-->cm/day)
       dsnow(:),                                & ! change in snow thickness (m/step-->cm/day)
       daidtt(:),                               & ! ice area tendency thermo.   (s^-1)
       dvidtt(:),                               & ! ice volume tendency thermo. (m/s)
       dagedtt(:),                              & ! ice age tendency thermo.    (s/s)
       mlt_onset(:),                            & ! day of year that sfc melting begins
       frz_onset(:),                            & ! day of year that freezing begins (congel or frazil)
       frazil_diag(:)                             ! frazil ice growth diagnostic (m/step-->cm/day)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat
       fsurfn(:,:),                             & ! category fsurf
       fcondtopn(:,:),                          & ! category fcondtop
       fcondbotn(:,:),                          & ! category fcondbot
       fsensn(:,:),                             & ! category sensible heat flux
       flatn(:,:)                                 ! category latent heat flux
     
  ! As above but these remain grid box mean values i.e. they are not divided by aice at end of ice_dynamics.
  ! These are used for generating ice diagnostics as these are more accurate. (The others suffer from problem
  ! of incorrect values at grid boxes that change from an ice free state to an icy state.)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       fresh_ai(:),                             & ! fresh water flux to ocean (kg/m^2/s)
       fsalt_ai(:),                             & ! salt flux to ocean (kg/m^2/s)
       fhocn_ai(:),                             & ! net heat flux to ocean (W/m^2)
       fswthru_ai(:),                           & ! shortwave penetrating to ocean (W/m^2)
       rside(:),                                & ! fraction of ice that melts laterally
       fside(:),                                & ! lateral heat flux (W/m^2)
       wlat(:),                                 & ! lateral melt rate (m/s)
       cos_zen(:),                              & ! cosine solar zenith angle, < 0 for sun below horizon
       rdg_conv_elem(:),                        & ! convergence term for ridging on elements (1/s)
       rdg_shear_elem(:),                       & ! shear term for ridging on elements (1/s)
       rdg_conv(:),                             & ! convergence term for ridging on nodes (1/s)
       rdg_shear(:)                               ! shear term for ridging on nodes (1/s)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nilyr+1
       salinz(:,:),                             & ! initial salinity profile (ppt)
       Tmltz(:,:)                                 ! initial melting temperature (C)

  !=======================================================================
  ! 4. Flux variables for biogeochemistry
  !=======================================================================
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_aero; COUPLING VARIABLE FOR BOTH TR_AERO AND TR_ZAERO, IN FROM ATMOSPHERE
       faero_atm(:,:)                             ! aerosol deposition rate (kg/m^2 s)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_nbtrcr; IN FROM ATMOSPHERE
       flux_bio_atm(:,:)                          ! all bio fluxes to ice from atmosphere
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_iso; IN FROM ATMOSPHERE
       fiso_atm(:,:),                           & ! isotope deposition rate (kg/m^2 s)
       fiso_evap(:,:),                          & ! isotope evaporation rate (kg/m^2 s)         
       Qa_iso(:,:),                             & ! isotope specific humidity (kg/kg)
       Qref_iso(:,:)                              ! 2m atm reference isotope spec humidity (kg/kg)

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_aero; IN FROM OCEAN
       fiso_ocn(:,:),                           & ! water isotope flux to the ocean (kg/m^2/s)
       faero_ocn(:,:)                             ! aerosol flux to ocean (kg/m^2/s)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_nbtrcr; OUT TO OCEAN
       flux_bio(:,:),                           & ! all bio fluxes to ocean
       flux_bio_ai(:,:)                           ! all bio fluxes to ocean, averaged over grid cell
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; OUT TO OCEAN
       fzsal_ai(:), & ! salt flux to ocean from zsalinity (kg/m^2/s)
       fzsal_g_ai(:)  ! gravity drainage salt flux to ocean (kg/m^2/s)
  
  logical (kind=log_kind), save ::              & ! INTERNAL
       cpl_bgc                                    ! switch to couple BGC via drivers
        
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; INTERNAL
       hin_old(:,:),                            & ! old ice thickness
       dsnown(:,:)                                ! change in snow thickness in category n (m)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; INTERNAL
       nit(:),                                  & ! ocean nitrate (mmol/m^3)
       amm(:),                                  & ! ammonia/um (mmol/m^3)
       sil(:),                                  & ! silicate (mmol/m^3)
       dmsp(:),                                 & ! dmsp (mmol/m^3)
       dms(:),                                  & ! dms (mmol/m^3)
       hum(:),                                  & ! humic material carbon (mmol/m^3)
       fnit(:),                                 & ! ice-ocean nitrate flux (mmol/m^2/s), positive to ocean
       famm(:),                                 & ! ice-ocean ammonia/um flux (mmol/m^2/s), positive to ocean
       fsil(:),                                 & ! ice-ocean silicate flux (mmol/m^2/s), positive to ocean
       fdmsp(:),                                & ! ice-ocean dmsp (mmol/m^2/s), positive to ocean
       fdms(:),                                 & ! ice-ocean dms (mmol/m^2/s), positive to ocean
       fhum(:),                                 & ! ice-ocean humic material carbon (mmol/m^2/s), positive to ocean
       fdust(:)                                   ! ice-ocean dust flux (kg/m^2/s), positive to ocean
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_algae; INTERNAL
       algalN(:,:),                             & ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeo)
       falgalN(:,:)                               ! ice-ocean algal nitrogen flux (mmol/m^2/s) (diatoms, pico, phaeo)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_doc; INTERNAL
       doc(:,:),                                & ! ocean doc (mmol/m^3) (saccharids, lipids, tbd )
       fdoc(:,:)                                  ! ice-ocean doc flux (mmol/m^2/s) (saccharids, lipids, tbd)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_don; INTERNAL
       don(:,:),                                & ! ocean don (mmol/m^3) (proteins and amino acids)
       fdon(:,:)                                  ! ice-ocean don flux (mmol/m^2/s) (proteins and amino acids)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_dic; INTERNAL
       dic(:,:),                                & ! ocean dic (mmol/m^3)
       fdic(:,:)                                  ! ice-ocean dic flux (mmol/m^2/s)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_fe; INTERNAL
       fed(:,:), fep(:,:),                      & ! ocean dissolved and particulate fe (nM)
       ffed(:,:), ffep(:,:)                       ! ice-ocean dissolved and particulate fe flux (umol/m^2/s)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_aero; INTERNAL
       zaeros(:,:)                                ! ocean aerosols (mmol/m^3)
  
  !=======================================================================
  ! 5. Column variables
  !=======================================================================

  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       Cdn_atm(:),                              & ! atm drag coefficient
       Cdn_ocn(:),                              & ! ocn drag coefficient
                                                  ! form drag:
       hfreebd(:),                              & ! freeboard (m)
       hdraft(:),                               & ! draft of ice + snow column (Stoessel 1993)
       hridge(:),                               & ! ridge height
       distrdg(:),                              & ! distance between ridges
       hkeel(:),                                & ! keel depth
       dkeel(:),                                & ! distance between keels
       lfloe(:),                                & ! floe length
       dfloe(:),                                & ! distance between floes
       meltsliq(:),                             & ! snow melt mass (kg/m^2/step-->kg/m^2/day) - icepack_snow.F90
       Cdn_atm_skin(:),                         & ! neutral skin drag coefficient
       Cdn_atm_floe(:),                         & ! neutral floe edge drag coefficient
       Cdn_atm_pond(:),                         & ! neutral pond edge drag coefficient
       Cdn_atm_rdg(:),                          & ! neutral ridge drag coefficient
       Cdn_ocn_skin(:),                         & ! skin drag coefficient
       Cdn_ocn_floe(:),                         & ! floe edge drag coefficient
       Cdn_ocn_keel(:),                         & ! keel drag coefficient
       Cdn_atm_ratio(:)                           ! ratio drag atm / neutral drag atm
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM 0:ncat ?; SEE ICEPACK_ITD.F90
       hin_max(:)                                 ! category limits (m)
  
  character (len=35), allocatable, save :: c_hi_range(:) ! DIM ncat
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; SEE ICEPACK_MELTPOND_LVL.F90
       dhsn(:,:),                               & ! depth difference for snow on sea ice and pond ice
       ffracn(:,:),                             & ! fraction of fsurfn used to melt ipond
       meltsliqn(:,:)                             ! snow melt mass in category n (kg/m^2)         
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat;  category albedos, SEE ICEPACK_SHORTWAVE.F90
       alvdrn(:,:), alidrn(:,:), alvdfn(:,:),  alidfn(:,:) ! visible direct/near IR direct/visible diffuse/near IR diffuse albedo (fraction)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; ALBEDO COMPONENTS FOR HISTORY
       albicen(:,:),                            & ! bare ice
       albsnon(:,:),                            & ! snow
       albpndn(:,:),                            & ! pond
       apeffn(:,:)                                ! effective pond area used for radiation calculation
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat; SHORTWAVE COMPONENTS
       snowfracn(:,:),                          & ! Category snow fraction used in radiation
       fswsfcn(:,:),                            & ! SW absorbed at ice/snow surface (W/m^2)
       fswthrun(:,:),                           & ! SW through ice to ocean (W/m^2)
       fswthrun_vdr(:,:), fswthrun_vdf(:,:),    & ! visible direct/diffusive SW through ice to ocean (W/m^2)
       fswthrun_idr(:,:), fswthrun_idf(:,:),    & ! near IR dirrct/diffusive SW through ice to ocean (W/m^2)
       fswintn(:,:)                               ! SW absorbed in ice interior, below surface (W/m^2)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nilyr,ncat; SHORTWAVE COMPONENTS
       Iswabsn(:,:,:)                             ! SW radiation absorbed in ice layers (W/m^2)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nslyr,ncat; SHORTWAVE COMPONENTS
       Sswabsn(:,:,:)                             ! SW radiation absorbed in snow layers (W/m^2)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nilyr+1,ncat; SHORTWAVE COMPONENTS
       fswpenln(:,:,:)                            ! visible SW entering ice layers (W/m^2)
  
  ! aerosol optical properties   -> band  |
  !                                       v aerosol
  ! for combined dust category, use category 4 properties
  ! dEdd 3-band and 5-band data
  !!! WHY SHOULD THAT BE PUBLIC FOR FESOM2???? FRANK.KAUKER@AWI.DE
  ! real (kind=dbl_kind), allocatable, save ::   & ! DIM icepack_nspint,icepack_max_aero
  !    gaer_3bd(:,:),                            & ! gaer_3bd, aerosolAsymmetryParameter3band
  !    kaer_3bd(:,:),                            & ! kaer_3bd, aerosolMassExtinctionCrossSection3band
  !    waer_3bd(:,:),                            & ! waer_3bd, aerosolSingleScatterAlbedo3band
  !    gaer_5bd(:,:),                            & ! gaer_5bd, aerosolAsymmetryParameter5band
  !    kaer_5bd(:,:),                            & ! kaer_5bd, aerosolMassExtinctionCrossSection5band
  !    waer_5bd(:,:)                               ! waer_5bd, aerosolSingleScatterAlbedo5band

  ! real (kind=dbl_kind), allocatable, public :: &
  !    gaer_bc_3bd(:,:),                         & ! gaer_bc_3bd, modalAsymmetryParameter3band
  !    kaer_bc_3bd(:,:),                         & ! kaer_bc_3bd, modalMassExtinctionCrossSection3band
  !    waer_bc_3bd(:,:),                         & ! waer_bc_3bd, modalSingleScatterAlbedo3band
  !    gaer_bc_5bd(:,:),                         & ! gaer_bc_5bd, modalAsymmetryParameter5band
  !    kaer_bc_5bd(:,:),                         & ! kaer_bc_5bd, modalMassExtinctionCrossSection5band
  !    waer_bc_5bd(:,:)                            ! waer_bc_5bd, modalSingleScatterAlbedo5band

  ! real (kind=dbl_kind), allocatable, save ::   & ! DIM icepack_nspint,icepack_nmodal1,icepack_nmodal2
  !    bcenh_3bd(:,:,:),                         & ! bcenh_3bd, modalBCabsorptionParameter3band
  !    bcenh_5bd(:,:,:)                            ! bcenh_5bd, modalBCabsorptionParameter5band
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nblyr+2; BIOGEOCHEMISTRY VARIABLES
       bgrid(:)                                   ! biology nondimensional vertical grid points
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nblyr+1; BIOGEOCHEMISTRY VARIABLES
       igrid(:)                                   ! biology vertical interface points
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nilyr+1; BIOGEOCHEMISTRY VARIABLES
       cgrid(:),                                & ! CICE vertical coordinate
       icgrid(:),                               & ! interface grid for CICE (shortwave variable)
       swgrid(:)                                  ! grid for ice tracers used in dEdd scheme
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat
       first_ice_real(:,:)                        ! .true. = c1, .false. = c0
  
  logical (kind=log_kind), allocatable, save :: & ! DIM nx,ncat
       first_ice(:,:)                             ! distinguishes ice that disappears (e.g. melts) and reappears (e.g. transport) in a grid cell
                                                  ! during a single time step from ice that was there the entire time step (true until ice forms)
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_nbtrcr
       ocean_bio(:,:)                             ! contains all the ocean bgc tracer concentrations
  
  real (kind=dbl_kind), allocatable, save ::    & !DIM nx,icepack_max_nbtrcr; DIAGNOSTIC FLUXES
       fbio_snoice(:,:),                        & ! fluxes from snow to ice
       fbio_atmice(:,:)                           ! fluxes from atm to ice
  
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_nbtrcr
       ocean_bio_all(:,:)        ! fixed order, all values even for tracers false
                                 ! N(1:max_algae) = 1:max_algae
                                 ! Nit = max_algae + 1
                                 ! DOC(1:max_doc) = max_algae + 2 : max_algae +
                                 ! max_doc + 1
                                 ! DIC(1:max_dic) = max_algae + max_doc + 2 :
                                 ! max_algae + max_doc + 1 + max_dic
                                 ! chl(1:max_algae) =  max_algae + max_doc + 2 +
                                 ! max_dic :
                                 !                   2*max_algae + max_doc + 1 +
                                 !                   max_dic
                                 ! Am =  2*max_algae + max_doc + 2 + max_dic
                                 ! Sil=  2*max_algae + max_doc + 3 + max_dic
                                 ! DMSPp=  2*max_algae + max_doc + 4 + max_dic
                                 ! DMSPd=  2*max_algae + max_doc + 5 + max_dic
                                 ! DMS  =  2*max_algae + max_doc + 6 + max_dic
                                 ! PON  =  2*max_algae + max_doc + 7 + max_dic
                                 ! DON(1:max_don)  =  2*max_algae + max_doc + 8 +
                                 ! max_dic :
                                 !                    2*max_algae + max_doc + 7 +
                                 !                    max_dic + max_don
                                 ! Fed(1:max_fe) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don :
                                 !                 2*max_algae + max_doc + 7 +
                                 !                 max_dic + max_don + max_fe
                                 ! Fep(1:max_fe) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don + max_fe :
                                 !                 2*max_algae + max_doc + 7 +
                                 !                 max_dic + max_don + 2*max_fe
                                 ! zaero(1:max_aero) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don + 2*max_fe :
                                 !                     2*max_algae + max_doc + 7 +
                                 !                     max_dic + max_don + 2*max_fe
                                 !                     + max_aero
                                 ! humic =  2*max_algae + max_doc + 8 + max_dic +
                                 ! max_don + 2*max_fe + max_aero
    
  integer (kind=int_kind), allocatable, save :: & ! DIM nx,icepack_max_algae
       algal_peak(:,:)                            ! vertical location of algal maximum, 0 if no maximum
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nblyr+1,ncat
       zoo(:,:,:)                                 ! N losses accumulated in timestep (ie. zooplankton/bacteria) (mmol/m^3)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat
       dhbr_top(:,:),                           & ! brine top change
       dhbr_bot(:,:)                              ! brine bottom change
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       grow_net(:),                             & ! Specific growth rate (/s) per grid cell
       PP_net(:),                               & ! Total production (mg C/m^2/s) per grid cell
       hbri(:)                                    ! brine height, area-averaged for comparison with hi(m)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nblyr+2,ncat
       bphi(:,:,:),                             & ! porosity of layers
       bTiz(:,:,:)                                ! layer temperatures interpolated on bio grid (C)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,ncat
       darcy_V(:,:)                               ! darcy velocity positive up (m/s)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       zsal_tot(:),                             & ! Total ice salinity in per grid cell (g/m^2)
       chl_net(:),                              & ! Total chla (mg chla/m^2) per grid cell
       NO_net(:)                                  ! Total nitrate per grid cell
    
  logical (kind=log_kind), allocatable, save :: & ! DIM nx
       Rayleigh_criteria(:)                       ! .true. means Ra_c was reached
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       Rayleigh_real(:)                           ! .true. = c1, .false. = c0
    
  real (kind=dbl_kind), allocatable, public ::  & ! DIM nx,ncat
       sice_rho(:,:)                              ! avg sea ice density (kg/m^3): ech: diagnostic only?
    
  real (kind=dbl_kind), allocatable, public ::  & ! DIM nx,ncat
       fzsaln(:,:),                             & ! category fzsal(kg/m^2/s)
       fzsaln_g(:,:)                              ! salt flux from gravity drainage only
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       fzsal(:),                                & ! total flux of salt to ocean at time step for conservation
       fzsal_g(:)                                 ! total gravity drainage flux
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nblyr+1,ncat
       zfswin(:,:,:)                              ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nblyr+1,ncat
       iDi(:,:,:),                              & ! igrid Diffusivity (m^2/s)
       iki(:,:,:)                                 ! ice permeability (m^2)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx
       upNO(:),                                 & ! nitrate uptake rate (mmol/m^2/d) times aice
       upNH(:)                                    ! ammonium uptake rate (mmol/m^2/d) times aice
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,max_ntrcr,ncat
       trcrn_sw(:,:,:)                            ! bgc tracers active in the delta-Eddington shortwave calculation on the shortwave grid (swgrid)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,icepack_max_nbtrcr
       ice_bio_net(:,:),                        & ! depth integrated tracer (mmol/m^2)
       snow_bio_net(:,:)                          ! depth integrated snow tracer (mmol/m^2)

  real(kind=dbl_kind), allocatable, save ::     & ! DIM nfsd; FLOE SIZE DISTRIBUTION
       floe_rad_l(:),                           & ! fsd size lower bound in m (radius)
       floe_rad_c(:),                           & ! fsd size bin centre in m (radius)
       floe_binwidth(:)                           ! fsd size bin width in m (radius)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx; FLOE SIZE DISTRIBUTION
       wave_sig_ht(:)                             ! significant height of waves (m)
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nfreq; FLOE SIZE DISTRIBUTION
       wavefreq(:),                             & ! wave frequencies
       dwavefreq(:)                               ! wave frequency bin widths
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nfreq; FLOE SIZE DISTRIBUTION
       wave_spectrum(:,:)                         ! wave spectrum
    
  real (kind=dbl_kind), allocatable, save ::    & ! DIM nx,nfsd; FLOE SIZE DISTRIBUTION
       d_afsd_newi(:,:), d_afsd_latg(:,:), d_afsd_latm(:,:), d_afsd_wave(:,:), d_afsd_weld(:,:) ! change in floe size distribution due to processes
    
  character (len=35), allocatable, save ::      & ! DIM nfsd
       c_fsd_range(:)                             ! fsd floe_rad bounds (m)

  !=======================================================================
  ! 6. Grid variables
  !=======================================================================

  logical (kind=log_kind), allocatable, save ::  & ! DIM
       lmask_n(:),                               & ! northern hemisphere mask  
       lmask_s(:)                                  ! northern hemisphere mask        

  real (kind=dbl_kind), allocatable, save ::     & ! DIM nx
       lon_val(:),                               & ! mesh nodes longitude
       lat_val(:)                                  ! mesh nodes latitude

  logical (kind=log_kind) :: l_print_point         ! flag for printing debugging information (pole: mype=7, i=5)
  
  !=======================================================================
  ! 7. Clock variables
  !=======================================================================

  ! The following variables should be sufficient to run icepack coupled to fesom 2. Restart and output will be handeled by fesom

  integer (kind=int_kind), save  ::              &
       days_per_year,                            & ! number of days in one year
       daymo(12),                                & ! number of days in each month
       daycal(13),                               & ! day number at end of month
       daycal365(13),                            &
       daycal366(13)
    
  data daycal365 /0,31,59,90,120,151,181,212,243,273,304,334,365/
  data daycal366 /0,31,60,91,121,152,182,213,244,274,305,335,366/
  
  integer (kind=int_kind), save ::               &
       istep1,                                   & ! counter, number of steps at current timestep
       mday,                                     & ! day of the month
       month_i,                                  & ! month number, 1 to 12
       nyr,                                      & ! year number
       sec                                         ! elapsed seconds into date
    
  real (kind=dbl_kind), save  ::                 &
       time,                                     & ! total elapsed time (s)
       yday,                                     & ! day of the year
       dayyr,                                    & ! number of days per year
       nextsw_cday,                              & ! julian day of next shortwave calculation
       secday,                                   & ! seconds per day
       dt_dyn                                      ! dynamics/transport/ridging timestep (s)
    
  character (len=char_len), save :: calendar_type  ! differentiates Gregorian from other calendars; default = ' '

  !=======================================================================
  !--------- Define the interface for submodules
  !=======================================================================

  interface

     ! Read icepack namelists, setup the model parameter and write diagnostics               
     module subroutine set_icepack(ice, partit)
       use mod_partit
       use MOD_ICE
       implicit none 
       type(t_partit), intent(inout), target :: partit
       type(t_ice), intent(inout), target :: ice
     end subroutine set_icepack

     ! Set up hemispheric masks 
     module subroutine set_grid_icepack(mesh)
       use MOD_MESH
       implicit none
       type(t_mesh), intent(in), target :: mesh
     end subroutine set_grid_icepack

     ! Allocate all
     module subroutine alloc_icepack()
       implicit none
     end subroutine alloc_icepack
     
     ! Initialize fluxes to atmosphere and ocean
     module subroutine init_flux_atm_ocn()
       implicit none
     end subroutine init_flux_atm_ocn

     ! Initialize thermodynamic history
     module subroutine init_history_therm()
       implicit none
     end subroutine init_history_therm
     
     ! Initialize dynamic hystory
     module subroutine init_history_dyn()
       implicit none
     end subroutine init_history_dyn
     
     ! Initialize bgc hystory
     module subroutine init_history_bgc()
       implicit none
     end subroutine init_history_bgc
     
     ! Initialize all
     module subroutine init_icepack(flag_debug, ice, tracer, mesh)
       use MOD_MESH
       use mod_tracer
       use MOD_ICE
       implicit none
       logical (kind=log_kind), intent(in) :: flag_debug
       type(t_mesh), intent(in), target :: mesh
       type(t_tracer_data), intent(in), target :: tracer
       type(t_ice), intent(inout), target :: ice
     end subroutine init_icepack
     
     ! Copy variables from fesom to icepack
     module subroutine fesom_to_icepack(ice, mesh)
       use MOD_MESH
       use MOD_ICE
       implicit none
       type(t_mesh), intent(in), target :: mesh
       type(t_ice), intent(inout), target :: ice
     end subroutine fesom_to_icepack

     ! Copy variables from fesom parameter to icepack
     module subroutine fesom_to_icepack_para()
       implicit none
     end subroutine fesom_to_icepack_para

     ! Copy variables from icepack to fesom
     module subroutine icepack_to_fesom( &
          nx_in, aice_out, vice_out, vsno_out, fhocn_tot_out, fresh_tot_out, &
          strocnxT_out, strocnyT_out, dhs_dt_out, dhi_dt_out, fsalt_out,     &
          evap_ocn_out, evap_out)
       use MOD_MESH
       implicit none        
       integer (kind=int_kind), intent(in) :: nx_in ! block dimensions        
       real (kind=dbl_kind), dimension(nx_in), intent(out), optional :: &
            aice_out, vice_out, vsno_out, fhocn_tot_out, fresh_tot_out, &
            strocnxT_out, strocnyT_out, fsalt_out, dhs_dt_out,          &
            dhi_dt_out, evap_ocn_out, evap_out
     end subroutine icepack_to_fesom

     ! Copy variables from icepack to fesom (single node or element)
     module subroutine icepack_to_fesom_single_point(nx_in, strength_out)
       use MOD_MESH
       implicit none        
       integer (kind=int_kind), intent(in) :: nx_in ! block dimensions        
       real (kind=dbl_kind), intent(out), optional :: strength_out
     end subroutine icepack_to_fesom_single_point

     ! Trancers advection 
     module subroutine tracer_advection_icepack(ice, mesh)
       use MOD_MESH
       use MOD_ICE
       implicit none
       type(t_mesh), intent(in), target :: mesh
       type(t_ice), intent(in), target :: ice
     end subroutine tracer_advection_icepack

     ! Advection initialization
     module subroutine init_advection_icepack(mesh)
       use MOD_MESH
       implicit none
       type(t_mesh), intent(in), target :: mesh
     end subroutine init_advection_icepack

     ! Driving subroutine for column physics
     module subroutine step_icepack(flag_debug, ice, mesh, time_evp, time_advec, time_therm)
       use MOD_MESH
       use MOD_ICE
       use g_config, only: dt
       use icepack_intfc, only: icepack_ice_strength
       implicit none
       logical (kind=log_kind), intent(in) :: flag_debug
       real (kind=dbl_kind), intent(out) :: time_therm, time_advec, time_evp
       type(t_mesh), intent(in), target  :: mesh
       type(t_ice), intent(inout), target  :: ice
     end subroutine step_icepack

     ! Initialize output
     module subroutine ini_mean_icepack_io(mesh)
       use MOD_MESH
       implicit none
       type(t_mesh), intent(in), target :: mesh
     end subroutine ini_mean_icepack_io

     ! Initialize restart
     module subroutine ini_icepack_io(year, partit, mesh)
       use MOD_MESH
       use mod_partit
       use mod_parsup
       implicit none
       type(t_mesh), intent(in) , target :: mesh
       type(t_partit), intent(inout), target :: partit
       integer(kind=int_kind), intent(in) :: year
     end subroutine ini_icepack_io
              
     ! Cut off Icepack
     module subroutine cut_off_icepack
       use icepack_intfc, only: icepack_compute_tracers, icepack_aggregate, icepack_init_trcr, icepack_sea_freezing_temperature
       use icepack_therm_shared, only: calculate_Tin_from_qin
       use icepack_mushy_physics, only: icepack_mushy_temperature_mush
       implicit none
     end subroutine cut_off_icepack

  end interface
end module icedrv_main
