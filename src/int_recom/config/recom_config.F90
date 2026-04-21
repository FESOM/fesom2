!>
!! @par Copyright
!! This code is subject to the FESOM-REcoM - License - Agreement in it's most recent form.
!! Please see URL xxx
!! 
!! @brief Module for defining variables used in REcoM, ex constant sinking velocity and 
!! local time step dt
!!
!! @remarks This module contains namelist for recom
!! @author xxx, FESOM-REcoM, Bremerhaven (2019)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   2014: original implementation (V. Schourup-Kristensen)
!
module recom_config
  implicit none
  save

!! *** General constants ***

! *******************
! CASE 2phy 1zoo 1det
! *******************
  Integer :: idin  =  1, idic   =  2, ialk    =  3, iphyn   =  4, iphyc = 5, &
             ipchl =  6, idetn  =  7, idetc   =  8, ihetn   =  9,            &
             ihetc = 10, idon   = 11, idoc    = 12, idian   = 13,            &
             idiac = 14, idchl  = 15, idiasi  = 16, idetsi  = 17,            &
             isi   = 18, ife    = 19, iphycal = 20, idetcal = 21,            &
             ioxy  = 22

  Integer :: izoo2n  = 23, izoo2c   = 24, idetz2n    = 25,                   &
             idetz2c = 26, idetz2si = 27, idetz2calc = 28

  Integer :: idicremin = 0  ! added by Sina

  ! Microzooplankton (third zooplankton group)
  integer :: imiczoon = 0   ! Microzooplankton Nitrogen (set below)
  integer :: imiczooc = 0   ! Microzooplankton Carbon (set below)

  ! ---------------------------------------------------------------------------
  ! PHYTOPLANKTON GROUPS (coccos configuration)
  ! ---------------------------------------------------------------------------
  ! Coccolithophores and Phaeocystis when enable_coccos = .true.

  integer :: icocn = 0      ! Coccolithophore Nitrogen (set below)
  integer :: icocc = 0      ! Coccolithophore Carbon (set below)
  integer :: icchl = 0      ! Coccolithophore Chlorophyll (set below)

  integer :: iphan   = 0    ! Phaeocystis Nitrogen (set below)
  integer :: iphac   = 0    ! Phaeocystis Carbon (set below)
  integer :: iphachl = 0    ! Phaeocystis Chlorophyll (set below)

  ! Terrestrial DOC input (when enable_R2OMIP is enabled)
  integer :: idoct = 0 

!=============================================================================

  Integer :: ivphy = 1, ivdia = 2, ivdet = 3, ivdetsc = 4, ivcoc = 5, ivpha = 6

!=============================================================================

  integer, dimension(8)  :: recom_remin_tracer_id   = (/1001, 1002, 1003, 1018, 1019, 1022, 1302, 1402/)

! The static declaration integer, dimension(32) :: recom_sinking_tracer_id 
! must remain size 32 (the full-model case uses all 32 slots), and the = 0 
! reset before partial fills ensures unused slots are inert.
  integer, dimension(32) :: recom_sinking_tracer_id = (/1007, 1008, 1017, 1021, 1004, 1005, 1020, 1006, &
                                                        1013, 1014, 1016, 1015, 1025, 1026, 1027, 1028, &
                                                        1029, 1030, 1031, &
                                                        1032, 1033, 1034, &
                                                        1308, 1321, 1305, 1320, & 
                                                        1314, 1408, 1421, 1405, 1420, 1414/)
  integer, dimension(8)  :: recom_det_tracer_id     = (/1007, 1008, 1017, 1021, 1308, 1321, 1408, 1421/)
  integer, dimension(8)  :: recom_phy_tracer_id     = (/1004, 1005, 1020, 1305, 1320, 1405, 1420, 1006/)
  integer, dimension(6)  :: recom_dia_tracer_id     = (/1013, 1014, 1314, 1414, 1016, 1015/)

  ! Configuration-dependent tracer arrays (allocated during initialization)
  integer, dimension(3)  :: recom_cocco_tracer_id
  integer, dimension(3)  :: recom_phaeo_tracer_id
  integer, dimension(4)  :: recom_det2_tracer_id

!=============================================================================

  Real(kind=8)                 :: zero           = 0.d0
  Integer                      :: one            = 1
  Real(kind=8)                 :: tiny           = 2.23D-16
  Real(kind=8)                 :: tiny_chl       = 0.00001 
  Real(kind=8)                 :: SecondsPerDay  = 86400.d0     ! [s/day]
  Real(kind=8)                 :: Pa2atm         = 101325.d0    ! [Pa/atm] 
  Real(kind=8)                 :: redO2C         = 1.453        ! O2:C ratio Anderson and Sarmiento, 1994

!! *** REcoM setup ***
  Logical                :: enable_3zoo2det = .false.   ! Control extended zooplankton variables
  Logical                :: enable_coccos = .false.     ! Control coccolithophore variables
  Logical                :: enable_AWICM = .false.      ! Control AWICM
  Logical                :: enable_R2OMIP = .false.     ! Control R2OMIP
  namelist /parecomsetup/ enable_3zoo2det, enable_coccos, enable_AWICM, enable_R2OMIP

!! *** General configuration ***

  Logical                :: use_REcoM            = .true.
  Logical                :: REcoM_restart        = .false.

  Integer                :: bgc_num               = 36      ! NEW increased the number from 28 to 34 (added coccos and respiration) ! NEW 3Zoo changed from 31 to 33 ! added phaeocystis: changed from 33 to 36
  integer                :: bgc_base_num          = 22      ! standard tracers
  Integer                :: diags3d_num           = 31      ! Number of diagnostic 3d tracers to be saved
  Real(kind=8)           :: VDet                  = 20.d0   ! Sinking velocity, constant through the water column and positive downwards
  Real(kind=8)           :: VDet_zoo2             = 200.d0  ! Sinking velocity, constant through the water column 
  Real(kind=8)           :: VPhy                  = 0.d0    !!! If the number of sinking velocities are different from 3, code needs to be changed !!!
  Real(kind=8)           :: VDia                  = 0.d0 
  Real(kind=8)           :: VCocco                = 0.d0    ! NEW 
  Real(kind=8)           :: VPhaeo                = 0.d0    ! Phaeocystis
  Logical                :: allow_var_sinking     = .true.   
  Integer                :: biostep               = 1          ! Number of times biology should be stepped forward for each time step		 
  Logical                :: REcoM_Geider_limiter  = .false.    ! Decides what routine should be used to calculate limiters in sms
  Logical                :: REcoM_Grazing_Variable_Preference = .true. ! Allow grazing preference to vary with food availability
  Logical                :: Grazing_detritus      = .false.    ! Decides grazing on detritus                            
  Logical                :: het_resp_noredfield   = .true.     ! Decides respiratation of copepods              
  Logical                :: diatom_mucus          = .true.     ! Effect of nutrient limitation on the aggregation
  Logical                :: O2dep_remin           = .true.     ! NEW O2remin Add option for O2 dependency of organic matter remineralization
  Logical                :: use_ballasting        = .true.     ! NEW BALL
  Logical                :: use_density_scaling   = .true.     ! NEW BALL
  Logical                :: use_viscosity_scaling = .true.     ! NEW BALL
  Logical                :: OmegaC_diss           = .true.     ! NEW DISS Use mocsy calcite omega to compute calcite dissolution
  Logical                :: CO2lim                = .true.     ! NEW Use CO2 dependence of growth and calcification
  !Logical                :: inter_CT_CL           = .true.    ! NEW inter use interaction between CO2 and both, temperature and light
  Logical                :: Diags                 = .true.     !!!!!!!!!!!!!!!!!!!!!!Change in recom.F90 Diagnostics -> Diags
  Logical                :: constant_CO2          = .true.
  Logical                :: UseFeDust             = .true.     ! Turns dust input of iron off when set to.false.
  Logical                :: UseDustClim           = .true.
  Logical                :: UseDustClimAlbani     = .true.     ! Use Albani dustclim field (If it is false Mahowald will be used)
  Logical                :: use_photodamage       = .false.    ! use Alvarez et al (2018) for chlorophyll degradation
  logical                :: HetRespFlux_plus      = .true.     ! MB More stable computation of zooplankton respiration fluxes adding a small number to HetN
  character(100)         :: REcoMDataPath         = '/albedo/work/projects/MarESys/ogurses/input/mesh_CORE2_finaltopo_mean/'
  logical                :: restore_alkalinity    = .true.
  logical                :: useRivers             = .false.
  logical                :: constant_PI_Rivers    = .true.
  logical                :: useRivFe              = .false.    ! river input of Fe
  logical                :: useErosion            = .false.
  logical                :: NitrogenSS            = .false.    ! This one only activates rivers! And in principle denitrification, but denitrification is commented out. When set to true, external sources and sinks of nitrogen are activated (Riverine, aeolian and denitrification)
  logical                :: useAeolianN           = .false.    ! When set to true, aeolian nitrogen deposition is activated
  integer                :: firstyearoffesomcycle = 1958       ! The first year of the actual physical forcing (e.g. JRA-55) used
  integer                :: lastyearoffesomcycle  = 2023       ! Last year of the actual physical forcing used
  integer                :: numofCO2cycles        = 1          ! Number of cycles of the forcing planned 
  integer                :: currentCO2cycle       = 1          ! Which CO2 cycle we are currently running
  Logical                :: DIC_PI                = .true.
  integer                :: Nmocsy                = 1          ! Length of the vector that is passed to mocsy (always one for recom)
  logical                :: recom_debug           = .false.
  logical                :: ciso                  = .false.    !MB main switch to enable/disable carbon isotopes (13|14C)
  integer                :: benthos_num           = 4          !MB number of sediment tracers = 8 if ciso = .true.
  Logical                :: use_MEDUSA            = .false.    ! main switch for sediment model
  integer                :: sedflx_num            = 0         ! number of sedimentary fluxs from MEDUSA, = 7 if ciso
  Logical                :: add_loopback          = .false.
  real(kind=8)           :: lb_tscale             = 1.d0      ! time scale to balance the burial loss
  integer                :: bottflx_num           = 4         ! number of stored sinking fluxes from the bottom layer, = 6 if C13 and = 8 if C14
  Logical                :: use_atbox             = .false.   ! switch for atmospheric box model for CO2

  Logical                :: fe_2ligands           = .false.    ! consider Fe-ligand binding with two ligands
  Logical                :: fe_compl_nica         = .false.    ! use Fe-ligand parameterisation dependent on DOC and pH (Ye2020)

  namelist /pavariables/ use_REcoM,                       REcoM_restart,                                  &
                       bgc_num,                           diags3d_num,           bgc_base_num,            &
                       VDet,                              VDet_zoo2,                                      &
                       VPhy,                              VDia,                  VCocco,                  &
                                                                                 VPhaeo,                  &
                       allow_var_sinking,                 biostep,               REcoM_Geider_limiter,    &
                       REcoM_Grazing_Variable_Preference,                                                 &
                       Grazing_detritus,        &
                                                                                 het_resp_noredfield,     &
                       diatom_mucus,                                                                      &
                       O2dep_remin,                       use_ballasting,        use_density_scaling,     & ! O2remin, NEW BALL
                       use_viscosity_scaling,             OmegaC_diss,           CO2lim,                  & ! BALL, DISS added OmegaC_diss, added CO2lim
                       Diags      ,                       constant_CO2,                                   &
                       UseFeDust,                         UseDustClim,           UseDustClimAlbani,       &
                       use_photodamage,                   HetRespFlux_plus,      REcoMDataPath,           &
                       restore_alkalinity,                useRivers,             constant_PI_Rivers,      &
                       useRivFe,                &
                       useErosion,                        NitrogenSS,            useAeolianN,             &
                       firstyearoffesomcycle,             lastyearoffesomcycle,  numofCO2cycles,          &
                       currentCO2cycle,                   DIC_PI,                Nmocsy,                  &
                       recom_debug,                       ciso,                  benthos_num,             &
                       use_MEDUSA,                        sedflx_num,            bottflx_num,             &
                       add_loopback,                      lb_tscale,             use_atbox,               &
                       fe_2ligands,                       fe_compl_nica

!!------------------------------------------------------------------------------
!! *** Sinking ***
  Real(kind=8)                 :: Vdet_a         = 0.0288       ! [1/day]
  Real(kind=8)                 :: Vcalc          = 0.0216       ! [1/day] depth dependence of calc_diss               
                                                                                                                           
  namelist /pasinking/ Vdet_a, Vcalc
!!------------------------------------------------------------------------------
!! *** Initialization ***
  Real(kind=8)                 :: cPhyN          = 0.2d0
  Real(kind=8)                 :: cHetN          = 0.2d0
  Real(kind=8)                 :: cZoo2N         = 0.2d0

  namelist /painitialization_N/ cPhyN, cHetN, cZoo2N
!!------------------------------------------------------------------------------
!! *** Temperature and Arrhenius functions ***
  Real(kind=8)                 :: recom_Tref     = 288.15d0       ! [K]
  Real(kind=8)                 :: C2K            = 273.15d0       !     Conversion from degrees C to K
  Real(kind=8)                 :: Ae             = 4500.d0        ! [K] Slope of the linear part of the Arrhenius function

!! *** Temperature variables for Blanchard function ***
  Real(kind=8)                 :: Tmax_phaeo     = 16d0           ! [degC] For Blanchard temp fxn: maximum temperature
  Real(kind=8)                 :: Topt_phaeo     = 7.5272d0       ! [degC] For Blanchard temp fxn: optimum temperature
  Real(kind=8)                 :: uopt_phaeo     = 0.7328d0       ! [1/day] For Blanchard function: optimum growth date
  Real(kind=8)                 :: beta_phaeo     = 0.7829d0       ! [unitless] For Blanchard function

! NEW MODIFIED parameters
  Real(kind=8)                 :: ord_d          = -0.2216d0 ! parameters for diatom temperature function
  Real(kind=8)                 :: expon_d        = 0.0406d0 ! diatom exponent
  Real(kind=8)                 :: ord_phy        = -1.2154d0 ! small phyto ordonnee
  Real(kind=8)                 :: expon_phy      = 0.0599d0 ! small phyto exponent
  Real(kind=8)                 :: ord_cocco      = -0.2310d0 ! coccolith ordonnee
  Real(kind=8)                 :: expon_cocco    = 0.0327d0 ! small phyto ordonnee
  Real(kind=8)                 :: ord_phaeo      = -0.2310d0 ! phaeocystis ordonnee
  Real(kind=8)                 :: expon_phaeo    = 0.0327d0 ! phaeocystis ordonnee

  Real(kind=8)                 :: reminSi        = 0.02d0
  Real(kind=8)                 :: k_o2_remin     = 15.d0          ! NEW O2remin mmol m-3; Table 1 in Cram 2018 cites DeVries & Weber 2017 for a range of 0-30 mmol m-3
  namelist /paArrhenius/ recom_Tref, C2K, Ae, Tmax_phaeo, Topt_phaeo, uopt_phaeo, beta_phaeo, ord_d, expon_d, ord_phy, expon_phy, ord_cocco, expon_cocco, ord_phaeo, expon_phaeo, reminSi, k_o2_remin

!!------------------------------------------------------------------------------
!! *** For limiter function ***
  Real(kind=8)                 :: NMinSlope      = 50.d0 
  Real(kind=8)                 :: SiMinSlope     = 1000.d0
  Real(kind=8)                 :: NCmin          = 0.04d0
  Real(kind=8)                 :: NCmin_d        = 0.04d0
  Real(kind=8)                 :: NCmin_c        = 0.04d0         ! NEW
  Real(kind=8)                 :: NCmin_p        = 0.04d0         ! Phaeocystis
  Real(kind=8)                 :: SiCmin         = 0.04d0
  Real(kind=8)                 :: k_Fe           = 0.04d0
  Real(kind=8)                 :: k_Fe_d         = 0.12d0
  Real(kind=8)                 :: k_Fe_c         = 0.04           ! NEW
  Real(kind=8)                 :: k_Fe_p         = 0.09           ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: k_si           = 4.d0
  Real(kind=8)                 :: P_cm           = 3.0d0          ! [1/day]   the rate of C-specific photosynthesis
  Real(kind=8)                 :: P_cm_d         = 3.5d0
  Real(kind=8)                 :: P_cm_c         = 3.3d0          ! NEW
  Real(kind=8)                 :: P_cm_p         = 3.4d0          ! NEW for Phaeocystis ( to be tuned)
  namelist /palimiter_function/ NMinSlope, SiMinSlope, NCmin, NCmin_d, NCmin_c, NCmin_p, SiCmin, k_Fe, k_Fe_d, k_Fe_c, k_Fe_p, k_si, P_cm, P_cm_d, P_cm_c, P_cm_p

!!------------------------------------------------------------------------------
!! *** For light calculations ***
  Real(kind=8)                 :: k_w            = 0.04d0         ! [1/m]              Light attenuation coefficient
  Real(kind=8)                 :: a_chl          = 0.03d0         ! [1/m * 1/(mg Chl)] Chlorophyll specific attenuation coefficients
  namelist /palight_calculations/ k_w, a_chl   
!!------------------------------------------------------------------------------
!! *** Photosynthesis ***
  Real(kind=8)                 :: alfa           = 0.14d0         ! [(mmol C*m2)/(mg Chl*W*day)]
  Real(kind=8)                 :: alfa_d         = 0.19d0         ! An initial slope of the P-I curve
  Real(kind=8)                 :: alfa_c         = 0.10d0         ! NEW
  Real(kind=8)                 :: alfa_p         = 0.10d0         ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: parFrac        = 0.43d0
  namelist /paphotosynthesis/ alfa, alfa_d, alfa_c, alfa_p, parFrac
!!------------------------------------------------------------------------------
!! *** Assimilation ***
  Real(kind=8)                 :: V_cm_fact      = 0.7d0          ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)                 :: V_cm_fact_d    = 0.7d0  
  Real(kind=8)                 :: V_cm_fact_c    = 0.7d0          ! NEW
  Real(kind=8)                 :: V_cm_fact_p    = 0.7d0          ! Phaeocystis
  Real(kind=8)                 :: NMaxSlope      = 1000.d0        ! Max slope for limiting function
  Real(kind=8)                 :: SiMaxSlope     = 1000.d0
  Real(kind=8)                 :: NCmax          = 0.2d0          ! [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
  Real(kind=8)                 :: NCmax_d        = 0.2d0
  Real(kind=8)                 :: NCmax_c        = 0.15d0         ! NEW
  Real(kind=8)                 :: NCmax_p        = 0.1d0          ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: SiCmax         = 0.8d0
  Real(kind=8)                 :: NCuptakeRatio  = 0.2d0          ! [mmol N/mmol C] Maximum uptake ratio of N:C
  Real(kind=8)                 :: NCUptakeRatio_d = 0.2d0
  Real(kind=8)                 :: NCUptakeRatio_c = 0.2d0         ! NEW
  Real(kind=8)                 :: NCUptakeRatio_p = 0.2d0         ! Phaeocystis
  Real(kind=8)                 :: SiCUptakeRatio = 0.2d0
  Real(kind=8)                 :: k_din          = 0.55d0         ! [mmol N/m3] Half-saturation constant for nitrate uptake
  Real(kind=8)                 :: k_din_d        = 1.0d0
  Real(kind=8)                 :: k_din_c        = 0.55d0         ! NEW
  Real(kind=8)                 :: k_din_p        = 0.55d0         ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: Chl2N_max      = 3.15d0         ! [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
  Real(kind=8)                 :: Chl2N_max_d    = 4.2d0
  Real(kind=8)                 :: Chl2N_max_c    = 3.5d0          ! NEW
  Real(kind=8)                 :: Chl2N_max_p    = 3.5d0          ! Phaeocystis (to be tuned (?))
  Real(kind=8)                 :: res_phy        = 0.01d0         ! [1/day] Maintenance respiration rate constant
  Real(kind=8)                 :: res_phy_d      = 0.01d0
  Real(kind=8)                 :: res_phy_c      = 0.0075d0       ! NEW
  Real(kind=8)                 :: res_phy_p      = 0.008d0        ! Phaeocystis (to be tuned (?))
  Real(kind=8)                 :: biosynth       = 2.33d0         ! [mmol C/mmol N] Cost of biosynthesis
  Real(kind=8)                 :: biosynthSi     = 0.d0
  namelist /paassimilation/ V_cm_fact, V_cm_fact_d, V_cm_fact_c, V_cm_fact_p, NMaxSlope, SiMaxSlope, NCmax, NCmax_d, NCmax_c, NCmax_p, SiCmax, &
                       NCuptakeRatio, NCUptakeRatio_d, NCUptakeRatio_c, NCUptakeRatio_p, SiCUptakeRatio, k_din, k_din_d, k_din_c, k_din_p, &
                       Chl2N_max, Chl2N_max_d, Chl2N_max_c, Chl2N_max_p, res_phy, res_phy_d, res_phy_c, res_phy_p, biosynth, biosynthSi
!!------------------------------------------------------------------------------
!! *** Iron chemistry ***
  Real(kind=8)                 :: totalligand     = 1.d0        ! [mumol/m3] order 1. Total free ligand
  Real(kind=8)                 :: ligandStabConst = 100.d0      ! [m3/mumol] order 100. Ligand-free iron stability constant
  namelist /pairon_chem/ totalligand, ligandStabConst
!!------------------------------------------------------------------------------
!! *** Zooplankton ***
  Real(kind=8)                 :: graz_max      = 2.4d0           ! [mmol N/(m3 * day)] Maximum grazing loss parameter 
  Real(kind=8)                 :: epsilonr      = 0.35d0          ! [(mmol N)2 /m6] Half saturation constant for grazing loss 
  Real(kind=8)                 :: res_het       = 0.01d0          ! [1/day] Respiration by heterotrophs (loss to DIC)
  Real(kind=8)                 :: Redfield      = 6.625           ! [mmol C/mmol N] Redfield ratio of C:N = 106:16
  Real(kind=8)                 :: loss_het      = 0.05d0          ! [m3/(mmolN*day)] The quadratic mortality rate (loss to detritus)
  Real(kind=8)                 :: pzDia         = 0.5d0           ! Maximum diatom preference
  Real(kind=8)                 :: sDiaNsq       = 0.d0
  Real(kind=8)                 :: pzPhy         = 1.0d0           ! Maximum nano-phytoplankton preference 
  Real(kind=8)                 :: sPhyNsq       = 0.d0
  Real(kind=8)                 :: pzCocco       = 0.5d0           ! NEW (value is just a guess)
  Real(kind=8)                 :: sCoccoNsq     = 0.d0            ! NEW
  Real(kind=8)                 :: pzPhaeo       = 1.0d0           ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: sPhaeoNsq     = 0.d0            ! Phaeocystis 
  Real(kind=8)                 :: pzMicZoo      = 1.0d0           ! NEW 3Zoo Maximum nano-phytoplankton preference
  Real(kind=8)                 :: sMicZooNsq    = 0.d0            ! NEW 3Zoo
  real(kind=8)                 :: tiny_het      = 1.d-5           ! for more stable computation of HetRespFlux (_plus). Value can be > tiny because HetRespFlux ~ hetC**2.
  namelist /pazooplankton/ graz_max, epsilonr, res_het, Redfield, loss_het, pzDia, sDiaNsq, pzPhy, sPhyNsq, pzCocco, sCoccoNsq, pzPhaeo, sPhaeoNsq, pzMicZoo, sMicZooNsq, tiny_het
!!-------------------------------------------------------------------------------                                                                                                                                                     
!! *** SecondZooplankton (Macrozooplankton) ***
  Real(kind=8)                 :: graz_max2      = 0.1d0          ! [mmol N/(m3 * day)] Maximum grazing loss parameter
  Real(kind=8)                 :: epsilon2       = 0.0144d0       ! [(mmol N)2 /m6] Half saturation constant for grazing loss
  Real(kind=8)                 :: res_zoo2       = 0.0107d0       ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)                                                            
  Real(kind=8)                 :: loss_zoo2      = 0.003d0        ! [1/day] Temperature dependent N degradation of extracellular organic N
  Real(kind=8)                 :: fecal_rate_n   = 0.13d0
  Real(kind=8)                 :: fecal_rate_c   = 0.295d0
  Real(kind=8)                 :: fecal_rate_n_mes = 0.25d0       ! NEW 3Zoo
  Real(kind=8)                 :: fecal_rate_c_mes = 0.32d0       ! NEW 3Zoo
  Real(kind=8)                 :: pzDia2         = 1.d0           ! Maximum diatom preference
  Real(kind=8)                 :: sDiaNsq2       = 0.d0
  Real(kind=8)                 :: pzPhy2         = 0.5d0          ! Maximum diatom preference
  Real(kind=8)                 :: sPhyNsq2       = 0.d0
  Real(kind=8)                 :: pzCocco2       = 0.2d0          ! NEW (value is just a guess)
  Real(kind=8)                 :: sCoccoNsq2     = 0.d0           ! NEW
  Real(kind=8)                 :: pzPhaeo2       = 0.5d0          ! Phaeocystis (to be tuned)
  Real(kind=8)                 :: sPhaeoNsq2     = 0.d0           ! Phaeocystis
  Real(kind=8)                 :: pzHet          = 0.8d0          ! Maximum diatom preference
  Real(kind=8)                 :: sHetNsq        = 0.d0
  Real(kind=8)                 :: pzMicZoo2      = 0.8d0          ! NEW Zoo3 Maximum microzooplankton preference
  Real(kind=8)                 :: sMicZooNsq2    = 0.d0           ! NEW Zoo3
  Real(kind=8)                 :: t1_zoo2        = 28145.d0       ! Krill temp. function constant1                                                                                                       
  Real(kind=8)                 :: t2_zoo2        = 272.5d0        ! Krill temp. function constant2
  Real(kind=8)                 :: t3_zoo2        = 105234.d0      ! Krill temp. function constant3                                                                                                      
  Real(kind=8)                 :: t4_zoo2        = 274.15d0       ! Krill temp. function constant3                                                                                              
  namelist /pasecondzooplankton/ graz_max2, epsilon2, res_zoo2, & 
                                 loss_zoo2, fecal_rate_n, fecal_rate_c, fecal_rate_n_mes, fecal_rate_c_mes,     &    ! NEW 3Zoo
                                 pzDia2, sDiaNsq2, pzPhy2, sPhyNsq2, pzCocco2, sCoccoNsq2, pzPhaeo2, sPhaeoNsq2, pzHet, sHetNsq,      &
                                 pzMicZoo2, sMicZooNsq2, &
                                 t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2
!-------------------------------------------------------------------------------
!! *** Third Zooplankton (Microzooplankton) ***
  Real(kind=8)                 :: graz_max3     = 0.46d0          ! NEW 3Zoo [mmol N/(m3 * day)] Maximum grazing loss parameter
  Real(kind=8)                 :: epsilon3      = 0.64d0          ! NEW 3Zoo [(mmol N)2 /m6] Half saturation constant for grazing loss
  Real(kind=8)                 :: loss_miczoo   = 0.01d0          ! NEW 3Zoo [1/day] Temperature dependent N degradation of extracellular organic N (EON)
  Real(kind=8)                 :: res_miczoo    = 0.01d0          ! NEW 3Zoo [1/day] Respiration by heterotrophs and mortality (loss to detritus)
  Real(kind=8)                 :: pzDia3        = 0.5d0           ! NEW 3Zoo Maximum diatom preference
  Real(kind=8)                 :: sDiaNsq3      = 0.d0            ! NEW 3Zoo
  Real(kind=8)                 :: pzPhy3        = 1.0d0           ! NEW 3Zoo Maximum nano-phytoplankton preference
  Real(kind=8)                 :: sPhyNsq3      = 0.d0            ! NEW 3Zoo
  Real(kind=8)                 :: pzCocco3      = 0.d0            ! NEW 3Zoo Maximum coccolithophore preference ! ATTENTION: This value needs to be tuned; I start with zero preference!
  Real(kind=8)                 :: sCoccoNsq3    = 0.d0            ! NEW 3Zoo
  Real(kind=8)                 :: pzPhaeo3      = 1.0d0           ! Phaeocystis 3Zoo Maximum phaeocystis preference (to be tuned (?))
  Real(kind=8)                 :: sPhaeoNsq3    = 0.d0            ! Phaeocystis 3Zoo
  namelist /pathirdzooplankton/ graz_max3, epsilon3, loss_miczoo, res_miczoo, pzDia3, sDiaNsq3, pzPhy3, sPhyNsq3, pzCocco3, sCoccoNsq3, pzPhaeo3, sPhaeoNsq3

!-------------------------------------------------------------------------------                                                                                                                          
!! *** Detritus Grazing Params ***                                                                                                                                                                        
  Real(kind=8)                 :: pzDet         = 1.d0           ! Maximum small detritus prefence by first zooplankton 
  Real(kind=8)                 :: sDetNsq       = 0.d0  
  Real(kind=8)                 :: pzDetZ2       = 1.d0           ! Maximum large detritus preference by first zooplankton
  Real(kind=8)                 :: sDetZ2Nsq     = 0.d0
  Real(kind=8)                 :: pzDet2        = 1.d0           ! Maximum small detritus prefence by second zooplankton
  Real(kind=8)                 :: sDetNsq2      = 0.d0  
  Real(kind=8)                 :: pzDetZ22      = 1.d0           ! Maximum large detritus preference by second zooplankton
  Real(kind=8)                 :: sDetZ2Nsq2    = 0.d0
  namelist /pagrazingdetritus/ pzDet, sDetNsq, pzDetZ2, sDetZ2Nsq, &
                                 pzDet2, sDetNsq2, pzDetZ22, sDetZ2Nsq2
!!------------------------------------------------------------------------------
!! *** Aggregation ***
  Real(kind=8)                 :: agg_PD        = 0.165d0         ! [m3/(mmol N * day)] Maximum aggregation loss parameter for DetN
  Real(kind=8)                 :: agg_PP        = 0.015d0         ! [m3/(mmol N * day)] Maximum aggregation loss parameter for PhyN and DiaN (plankton)
  namelist /paaggregation/ agg_PD, agg_PP
!!------------------------------------------------------------------------------
!! *** DIN ***
  Real(kind=8)                 :: rho_N         = 0.11d0          ! [1/day] Temperature dependent N degradation of extracellular organic N (EON) (Remineralization of DON)
  Real(kind=8)                 :: rho_Nt        = 0.67d0          ! [1/day] terrestrial N degradation of extracellular organic N (DONt) (Remineralization of terrestrial DON) ~ scale of 1.5 years ! R2OMIP
  namelist /padin_rho_N/ rho_N, rho_Nt
!!------------------------------------------------------------------------------
!! *** DIC ***
  Real(kind=8)                 :: rho_C1        = 0.1d0           ! [1/day] Temperature dependent C degradation of extracellular organic C (EOC)
  Real(kind=8)                 :: rho_C1t       = 0.0018d0        ! [1/day] terrestrial C degradation of extracellular organic C (EOC)
  namelist /padic_rho_C1/ rho_C1, rho_C1t
!!------------------------------------------------------------------------------
!! *** Phytoplankton N ***
  Real(kind=8)                 :: lossN         = 0.05d0          ! [1/day] Phytoplankton loss of organic N compounds
  Real(kind=8)                 :: lossN_d       = 0.05d0
  Real(kind=8)                 :: lossN_c       = 0.05d0
  Real(kind=8)                 :: lossN_p       = 0.05d0          ! Phaeocystis
  namelist /paphytoplankton_N/ lossN, lossN_d, lossN_c, lossN_p
!!------------------------------------------------------------------------------
!! *** Phytoplankton C ***
  Real(kind=8)                 :: lossC         = 0.10d0          ! [1/day] Phytoplankton loss of carbon 
  Real(kind=8)                 :: lossC_d       = 0.10d0
  Real(kind=8)                 :: lossC_c       = 0.10d0
  Real(kind=8)                 :: lossC_p       = 0.10d0          ! Phaeocystis
  namelist /paphytoplankton_C/ lossC, lossC_d, lossC_c, lossC_p
!!------------------------------------------------------------------------------
!! *** Phytoplankton ChlA ***
  Real(8)                      :: deg_Chl       = 0.25d0          ! [1/day]
  Real(kind=8)                 :: deg_Chl_d     = 0.25d0
  Real(kind=8)                 :: deg_Chl_c     = 0.20d0          ! (value is just a guess)
  Real(kind=8)                 :: deg_Chl_p     = 0.25d0          ! Phaeocystis
  namelist /paphytoplankton_ChlA/ deg_Chl, deg_Chl_d, deg_Chl_c, deg_Chl_p
!!------------------------------------------------------------------------------
!! *** Detritus N ***
  Real(kind=8)                 :: gfin          = 0.3d0         ! 3Zoo [] Grazing efficiency (fraction of grazing flux into zooplankton pool) 
  Real(kind=8)                 :: grazEff2      = 0.8d0         ! 3Zoo [] Grazing efficiency (fraction of grazing flux into second zooplankton pool) 
  Real(kind=8)                 :: grazEff3      = 0.8d0         ! 3Zoo [] Grazing efficiency (fraction of grazing flux into microzooplankton pool)
  Real(kind=8)                 :: reminN        = 0.165d0       ! 3Zoo [1/day] Temperature dependent remineralisation rate of detritus 
  namelist /padetritus_N/ gfin, grazEff2, grazEff3, reminN      
!!------------------------------------------------------------------------------
!! *** Detritus C ***
  Real(kind=8)                 :: reminC        = 0.15d0        ! [1/day] Temperature dependent remineralisation rate of detritus
  Real(kind=8)                 :: rho_c2        = 0.1d0         ! [1/day] Temperature dependent C degradation of TEP-C
  namelist /padetritus_C/ reminC, rho_c2
!!------------------------------------------------------------------------------
!! *** Heterotrophs ***
  Real(kind=8)                 :: lossN_z       = 0.15d0
  Real(kind=8)                 :: lossC_z       = 0.15d0
  namelist /paheterotrophs/ lossN_z, lossC_z
!!------------------------------------------------------------------------------
!! *** Second Zooplankton ***             
  Real(kind=8)                 :: lossN_z2       = 0.02d0
  Real(kind=8)                 :: lossC_z2       = 0.02d0
  namelist /paseczooloss/ lossN_z2, lossC_z2
!!-----------------------------------------------------------------------------
!! *** Third Zooplankton ***
  Real(kind=8)                 :: lossN_z3      = 0.05d0         ! 3Zoo
  Real(kind=8)                 :: lossC_z3      = 0.05d0         ! 3Zoo
  namelist /pathirdzooloss/ lossN_z3, lossC_z3
!!------------------------------------------------------------------------------                                                                                                                                                                                               
!! *** Parameters for CO2 limitation ***                         

  Real(kind=8)                 :: Cunits         = 976.5625      ! Conversion factor between [mol/m3] (model) and [umol/kg] (function): (1000 * 1000) / 1024
  Real(kind=8)                 :: a_co2_phy      = 1.162e+00     ! [unitless]
  Real(kind=8)                 :: a_co2_dia      = 1.040e+00     ! [unitless]
  Real(kind=8)                 :: a_co2_cocco    = 1.109e+00     ! [unitless]
  Real(kind=8)                 :: a_co2_phaeo    = 1.162e+00     ! [unitless]
  Real(kind=8)                 :: a_co2_calc     = 1.102e+00     ! [unitless]
  Real(kind=8)                 :: b_co2_phy      = 4.888e+01     ! [mol/kg]
  Real(kind=8)                 :: b_co2_dia      = 2.890e+01     ! [mol/kg]
  Real(kind=8)                 :: b_co2_cocco    = 3.767e+01     ! [mol/kg]
  Real(kind=8)                 :: b_co2_phaeo    = 4.888e+01     ! [mol/kg]
  Real(kind=8)                 :: b_co2_calc     = 4.238e+01     ! [mol/kg]
  Real(kind=8)                 :: c_co2_phy      = 2.255e-01     ! [kg/mol]
  Real(kind=8)                 :: c_co2_dia      = 8.778e-01     ! [kg/mol]
  Real(kind=8)                 :: c_co2_cocco    = 3.912e-01     ! [kg/mol]
  Real(kind=8)                 :: c_co2_phaeo    = 2.255e-01     ! [kg/mol]
  Real(kind=8)                 :: c_co2_calc     = 7.079e-01     ! [kg/mol]
  Real(kind=8)                 :: d_co2_phy      = 1.023e+07     ! [kg/mol]
  Real(kind=8)                 :: d_co2_dia      = 2.640e+06     ! [kg/mol]
  Real(kind=8)                 :: d_co2_cocco    = 9.450e+06     ! [kg/mol]
  Real(kind=8)                 :: d_co2_phaeo    = 1.023e+07     ! [kg/mol]
  Real(kind=8)                 :: d_co2_calc     = 1.343e+07     ! [kg/mol]
  namelist /paco2lim/ Cunits, a_co2_phy, a_co2_dia, a_co2_cocco, a_co2_phaeo, a_co2_calc, &
                      b_co2_phy, b_co2_dia, b_co2_cocco, b_co2_phaeo, b_co2_calc, &
                      c_co2_phy, c_co2_dia, c_co2_cocco, c_co2_phaeo, c_co2_calc, &
                      d_co2_phy, d_co2_dia, d_co2_cocco, d_co2_phaeo, d_co2_calc
!!------------------------------------------------------------------------------
!! *** Iron ***
  Real(kind=8)                 :: Fe2N           = 0.033d0       ! Fe2C * 6.625 (Fe2C = 0.005d0)
  Real(kind=8)                 :: Fe2N_benthos   = 0.15d0        ! default was 0.14 Fe2C_benthos (=0.02125=0.68d0/32.d0) * 6.625 - will have to be tuned. [umol/m2/day]
  Real(kind=8)                 :: kScavFe        = 0.07d0
  Real(kind=8)                 :: dust_sol       = 0.02d0        !Dissolution of Dust for bioavaliable
  Real(kind=8)                 :: RiverFeConc   = 100d0        ! mean DFe concentration in rivers
  namelist /pairon/ Fe2N, Fe2N_benthos, kScavFe, dust_sol, RiverFeConc
!!------------------------------------------------------------------------------
!! *** Calcification ***
  Real(kind=8)                 :: calc_prod_ratio = 0.02d0
  Real(kind=8)                 :: calc_diss_guts  = 0.0d0
  Real(kind=8)                 :: calc_diss_rate  = 0.005714d0    !20.d0/3500.d0
  Real(kind=8)                 :: calc_diss_rate2 = 0.005714d0
  Real(kind=8)                 :: calc_diss_omegac = 0.197d0      ! NEW DISS value from Aumont et al. 2015, will be used with OmegaC_diss flag
  Real(kind=8)                 :: calc_diss_exp   = 1.d0          ! NEW DISS exponent in the dissolution rate of calcite, will be used with OmegaC_diss flag
  namelist /pacalc/ calc_prod_ratio, calc_diss_guts, calc_diss_rate, calc_diss_rate2, calc_diss_omegac, calc_diss_exp  ! NEW DISS added calc_diss_omegac, calc_diss_exp
!!------------------------------------------------------------------------------
!! *** Benthos ***
  Real(kind=8)                 :: decayRateBenN   = 0.005d0
  Real(kind=8)                 :: decayRateBenC   = 0.005d0
  Real(kind=8)                 :: decayRateBenSi  = 0.005d0
  Real(kind=8)                 :: q_NC_Denit      = 0.86d0         ! N:C quota of the denitrification process
  namelist /pabenthos_decay_rate/ decayRateBenN, decayRateBenC, decayRateBenSi, q_NC_Denit
!!------------------------------------------------------------------------------
!! *** CO2-flux ***
  Real(kind=8)                 :: permil          = 0.000000976    ! 1.e-3/1024.5d0 ! Converting DIC from [mmol/m3] to [mol/kg]
  Real(kind=8)                 :: permeg          = 1.e-6          ! [atm/uatm] Changes units from uatm to atm
  Real(kind=8)                 :: Xacc            = 1.e-12         ! Accuracy for ph-iteration (phacc)
!  Real(kind=8)                :: pCO2a           = 380.d0         ! [uatm] Atmospheric partial pressure of CO2
  Real(kind=8)                 :: CO2_for_spinup  = 278.d0         !  
  namelist /paco2_flux_param/ permil, permeg, Xacc, CO2_for_spinup
!!------------------------------------------------------------------------------
!! *** Alkalinity restoring ***
  Real(kind=8)                 :: surf_relax_Alk = 3.2e-07 !10.d0/31536000.d0
  namelist /paalkalinity_restoring/ surf_relax_Alk
!!-----------------------------------------------------------------------------
!! *** Ballasting ***                                              ! NEW BALL
  Real(kind=8)                 :: rho_POC              = 1033.d0   ! kg m-3; density of POC (see Table 1 in Cram et al., 2018)
  Real(kind=8)                 :: rho_PON              = 1033.d0   ! kg m-3; density of PON (see Table 1 in Cram et al., 2018)
  Real(kind=8)                 :: rho_CaCO3            = 2830.d0   ! kg m-3; density of CaCO3 (see Table 1 in Cram et al., 2018) 
  Real(kind=8)                 :: rho_opal             = 2090.d0   ! kg m-3; density of Opal (see Table 1 in Cram et al., 2018)
  Real(kind=8)                 :: rho_ref_part         = 1230.d0   ! kg m-3; reference particle density (see Cram et al., 2018)
  Real(kind=8)                 :: rho_ref_water        = 1027.d0   ! kg m-3; reference seawater density (see Cram et al., 2018)
  Real(kind=8)                 :: visc_ref_water       = 0.d00158  ! kg m-1 s-1; reference seawater viscosity, at Temp=4 degC (see Cram et al., 2018)
  Real(kind=8)                 :: w_ref1               = 10.d0     ! m s-1; reference sinking velocity of small detritus
  Real(kind=8)                 :: w_ref2               = 200.d0    ! m s-1; reference sinking velocity of large detritus
  Real(kind=8)                 :: depth_scaling1       = 0.d015    ! s-1; factor to increase sinking speed of det1 with depth, set to 0 if not wanted
  Real(kind=8)                 :: depth_scaling2       = 0.d0      ! s-1; factor to increase sinking speed of det2 with depth, set to 0 if not wanted
  Real(kind=8)                 :: max_sinking_velocity = 250.d0    ! d-1; for numerical stability, set a maximum possible sinking velocity here (applies to both detritus classes)
  namelist /paballasting/ rho_POC, rho_PON, rho_CaCO3, rho_opal, rho_ref_part, &
                          rho_ref_water, visc_ref_water, w_ref1, w_ref2, depth_scaling1,   &
                          depth_scaling2, max_sinking_velocity

contains

  ! ---------------------------------------------------------------------------
  ! SUBROUTINE: initialize_tracer_indices
  ! ---------------------------------------------------------------------------
  ! Purpose: Set up tracer indices based on model configuration
  ! ---------------------------------------------------------------------------
  subroutine initialize_tracer_indices()
    implicit none

    ! ---------------------------------------------------------------------------
    ! Base sinking tracers (always present in all configurations):
    !   Det1:  DetN(1007), DetC(1008), DetSi(1017), DetCalc(1021)
    !   Phy:   PhyN(1004), PhyC(1005), PhyCalc(1020), PhyChl(1006)
    !   Dia:   DiaN(1013), DiaC(1014), DiaSi(1016), DiaChl(1015)
    ! ---------------------------------------------------------------------------

    if (enable_3zoo2det .and. enable_coccos) then
        ! =======================================================================
        ! CASE: 4 phytoplankton + 3 zooplankton + 2 detritus
        ! =======================================================================
        ! Phytoplankton: small phyto, diatoms, coccolithophores, phaeocystis
        ! Zooplankton: mesozoo, macrozoo, microzoo
        ! Detritus: det1, det2

        icocn    = 29
        icocc    = 30
        icchl    = 31
        iphan    = 32
        iphac    = 33
        iphachl  = 34
        imiczoon = 35
        imiczooc = 36

        ! Terrestrial DOC input (when enable_R2OMIP is enabled)
        if (enable_R2OMIP) then
            idoct = 37
            idicremin = 38
        else
            idicremin = 37
        end if

        recom_cocco_tracer_id = (/1029, 1030, 1031/)
        recom_phaeo_tracer_id = (/1032, 1033, 1034/)
        recom_det2_tracer_id  = (/1025, 1026, 1027, 1028/)

        ! Base(12) + Det2(4) + Zoo2Det2(4) + Cocco(3) + Phaeo(3) + CoccoCalc(2)
        ! Det1 sinking: 1007,1008,1017,1021
        ! Phy  sinking: 1004,1005,1020,1006
        ! Dia  sinking: 1013,1014,1016,1015
        ! Det2 sinking: 1025,1026,1027,1028  (DetZ2N,DetZ2C,DetZ2Si,DetZ2Calc)
        ! Zoo2 det:     1308,1321            (via zoo2 sinking tracer IDs)
        ! Phy  cal:     1305,1320
        ! Dia  cal:     1314
        ! Cocco sinking:1029,1030,1031 -> mapped as 1308..? 
        ! 
        ! Reconstructing from original array for full model:
        ! Original: 1007,1008,1017,1021, 1004,1005,1020,1006,
        !           1013,1014,1016,1015, 1025,1026,1027,1028,
        !           1029,1030,1031,
        !           1032,1033,1034,
        !           1308,1321,1305,1320,
        !           1314,1408,1421,1405,1420,1414

        recom_sinking_tracer_id      = 0
        recom_sinking_tracer_id(1:22) = &
            (/1007, 1008, 1017, 1021, 1004, 1005, 1020, 1006, &
              1013, 1014, 1016, 1015, 1025, 1026, 1027, 1028, &
              1029, 1030, 1031, 1032, 1033, 1034/)
        if (ciso) then
            recom_sinking_tracer_id(23:32) = &
                (/1308, 1321, 1305, 1320, &
                  1314, 1408, 1421, 1405, 1420, 1414/)
        end if

    else if (enable_coccos .and. .not. enable_3zoo2det) then
        ! =======================================================================
        ! CASE: 4 phytoplankton + 1 zooplankton + 1 detritus
        ! =======================================================================
        ! Phytoplankton: small phyto, diatoms, coccolithophores, phaeocystis
        ! Zooplankton: mesozoo only
        ! Detritus: det1 only

        icocn   = 23
        icocc   = 24
        icchl   = 25
        iphan   = 26
        iphac   = 27
        iphachl = 28

        ! Terrestrial DOC input (when enable_R2OMIP is enabled)
        if (enable_R2OMIP) then
            idoct = 29
            idicremin = 30
        else
            idicremin = 29
        end if

        recom_cocco_tracer_id = (/1023, 1024, 1025/)
        recom_phaeo_tracer_id = (/1026, 1027, 1028/)

        ! Det1 sinking: 1007,1008,1017,1021
        ! Phy  sinking: 1004,1005,1020,1006
        ! Dia  sinking: 1013,1014,1016,1015
        ! Cocco sinking:1023,1024,1025
        ! Phaeo sinking:1026,1027,1028
        ! Ciso-style:   1308,1321,1305,1320,1314,1408,1421,1405,1420,1414
        ! No Det2, no Zoo2 sinking
        recom_sinking_tracer_id      = 0
        recom_sinking_tracer_id(1:18) = &
            (/1007, 1008, 1017, 1021, 1004, 1005, 1020, 1006, &
              1013, 1014, 1016, 1015, 1023, 1024, 1025,        &
              1026, 1027, 1028/)
        if (ciso) then
            recom_sinking_tracer_id(19:28) = &
                (/1308, 1321, 1305, 1320, &
                  1314, 1408, 1421, 1405, 1420, 1414/)
        end if

    else if (enable_3zoo2det .and. .not. enable_coccos) then
        ! =======================================================================
        ! CASE: 2 phytoplankton + 3 zooplankton + 2 detritus
        ! =======================================================================
        ! Phytoplankton: small phyto, diatoms only
        ! Zooplankton: mesozoo, macrozoo, microzoo
        ! Detritus: det1, det2

        imiczoon = 29
        imiczooc = 30

        ! Terrestrial DOC input (when enable_R2OMIP is enabled)
        if (enable_R2OMIP) then
            idoct = 31
            idicremin = 32
        else
            idicremin = 31
        end if

        recom_det2_tracer_id = (/1025, 1026, 1027, 1028/)

        ! Det1 sinking: 1007,1008,1017,1021
        ! Phy  sinking: 1004,1005,1020,1006
        ! Dia  sinking: 1013,1014,1016,1015
        ! Det2 sinking: 1025,1026,1027,1028
        ! Ciso-style:   1308,1321,1305,1320,1314,1408,1421,1405,1420,1414
        ! No Cocco/Phaeo sinking
        recom_sinking_tracer_id      = 0
        recom_sinking_tracer_id(1:16) = &
            (/1007, 1008, 1017, 1021, 1004, 1005, 1020, 1006, &
              1013, 1014, 1016, 1015, 1025, 1026, 1027, 1028/)
        if (ciso) then
            recom_sinking_tracer_id(17:26) = &
                (/1308, 1321, 1305, 1320, &
                  1314, 1408, 1421, 1405, 1420, 1414/)
        end if

    else
        ! =======================================================================
        ! CASE: 2 phytoplankton + 1 zooplankton + 1 detritus (BASE CONFIGURATION)
        ! =======================================================================
        ! Phytoplankton: small phyto, diatoms only
        ! Zooplankton: mesozoo only
        ! Detritus: det1 only
        ! (All indices already set to default values)

        ! Terrestrial DOC input (when enable_R2OMIP is enabled)
        if (enable_R2OMIP) then
            idoct = 23
            idicremin = 24
        else
            idicremin = 23
        end if

        ! Det1 sinking: 1007,1008,1017,1021
        ! Phy  sinking: 1004,1005,1020,1006
        ! Dia  sinking: 1013,1014,1016,1015
        ! Ciso-style:   1308,1321,1305,1320,1314,1408,1421,1405,1420,1414
        ! No Det2, no Cocco/Phaeo sinking
        recom_sinking_tracer_id      = 0
        recom_sinking_tracer_id(1:12) = &
            (/1007, 1008, 1017, 1021, 1004, 1005, 1020, 1006, &
              1013, 1014, 1016, 1015/)
        if (ciso) then
            recom_sinking_tracer_id(13:22) = &
                (/1308, 1321, 1305, 1320, &
                  1314, 1408, 1421, 1405, 1420, 1414/)
        end if

    endif
  end subroutine initialize_tracer_indices

! ==============================================================================
! SUBROUTINE: validate_recom_tracers
! ==============================================================================
! Purpose: Validate consistency between namelist tracer configuration and
!          biogeochemical model setup (enable_3zoo2det, enable_coccos)
! ==============================================================================
subroutine validate_recom_tracers(num_tracers, mype)
  
  use g_forcing_param, only: use_age_tracer !---age-code
  
  implicit none

  ! Arguments
  integer, intent(in) :: num_tracers  ! Total number of tracers from namelist
  integer, intent(in) :: mype         ! MPI rank

  ! Local variables
  integer :: expected_bgc_num
  integer :: actual_bgc_num
  integer :: expected_total_tracers
  integer :: num_physical_tracers
  integer :: doc_tracers
  logical :: config_error
  character(len=200) :: error_msg

  ! For tracer ID validation
  integer :: i, tracer_id
  integer, dimension(:), allocatable :: expected_tracer_ids
  logical, dimension(:), allocatable :: tracer_found
  integer :: num_expected_tracers
  logical :: id_error

  ! Physical tracers (temperature, salinity, etc.) - typically first 2
  num_physical_tracers = 2

  ! DOC tracers count
  if (enable_R2OMIP) then
    doc_tracers = 1
  else
    doc_tracers = 0
  end if

  ! Calculate actual BGC tracer count from namelist
  actual_bgc_num = num_tracers - num_physical_tracers

  ! ===========================================================================
  ! Determine expected BGC tracer count based on configuration
  ! ===========================================================================
  config_error = .false.

  if (enable_3zoo2det .and. enable_coccos) then
    ! ---------------------------------------------------------------------------
    ! Configuration 4: Full model (4 phyto + 3 zoo + 2 detritus)
    ! ---------------------------------------------------------------------------
    ! Base: 22 tracers (1001-1022)
    ! Additional 3zoo2det: 4 tracers for det2 (1025-1028)
    ! Additional coccos: 6 tracers for coccos (1029-1031)
    ! Additional phaeocystis: 3 tracers (1032-1034)
    ! Additional microzoo: 2 tracers (1035-1036)
    ! Total: 22 + 4 + 6 + 3 + 2 = 36 (actually 22 + 14 = 36)
    ! Additional DICremin: 1 tracer (1037) (added by Sina)
    ! Plus DOC: +1 if enable_R2OMIP
    expected_bgc_num = 37 + doc_tracers ! Sina: changed from 36 to 37

  else if (enable_coccos .and. .not. enable_3zoo2det) then
    ! ---------------------------------------------------------------------------
    ! Configuration 3: Coccos only (4 phyto + 1 zoo + 1 detritus)
    ! ---------------------------------------------------------------------------
    ! Base: 22 tracers (1001-1022)
    ! Additional coccos: 3 tracers (1023-1025)
    ! Additional phaeocystis: 3 tracers (1026-1028)
    ! Total: 22 + 6 = 28
    ! Additional DICremin: 1 tracer (1037) (added by Sina)
    ! Plus DOC: +1 if enable_R2OMIP
    expected_bgc_num = 29 + doc_tracers ! Sina: changed from 28 to 29

  else if (enable_3zoo2det .and. .not. enable_coccos) then
    ! ---------------------------------------------------------------------------
    ! Configuration 2: 3Zoo2Det only (2 phyto + 3 zoo + 2 detritus)
    ! ---------------------------------------------------------------------------
    ! Base: 22 tracers (1001-1022)
    ! Additional zoo2: 2 tracers (1023-1024)
    ! Additional det2: 4 tracers (1025-1028)
    ! Additional microzoo: 2 tracers (1029-1030)
    ! Total: 22 + 8 = 30
    ! Additional DICremin: 1 tracer (1037) (added by Sina)
    ! Plus DOC: +1 if enable_R2OMIP
    expected_bgc_num = 31 + doc_tracers ! Sina: changed from 30 to 31

  else
    ! ---------------------------------------------------------------------------
    ! Configuration 1: Base model (2 phyto + 1 zoo + 1 detritus)
    ! ---------------------------------------------------------------------------
    ! Base: 22 tracers (1001-1022)
    ! Additional DICremin: 1 tracer (1037) (added by Sina)
    ! Plus DOC: +1 if enable_R2OMIP
    expected_bgc_num = 23 + doc_tracers ! Sina: changed from 22 to 23

  end if

  if (use_age_tracer) then 
    expected_bgc_num = expected_bgc_num + 1 
  end if

  expected_total_tracers = num_physical_tracers + expected_bgc_num

  ! ===========================================================================
  ! Build expected tracer ID list for current configuration
  ! ===========================================================================

  ! Determine total expected tracers
  num_expected_tracers = expected_total_tracers
  allocate(expected_tracer_ids(num_expected_tracers))
  allocate(tracer_found(num_expected_tracers))
  tracer_found = .false.

  ! Physical tracers (always present)
  expected_tracer_ids(1) = 1    ! Temperature
  expected_tracer_ids(2) = 2    ! Salinity

  ! Base BGC tracers (always present for all configurations)
  do i = 1, 22
    expected_tracer_ids(num_physical_tracers + i) = 1000 + i
  end do

  ! Configuration-specific tracers
  if (enable_3zoo2det .and. enable_coccos) then
    ! Full model: 1001-1022 (base) + 1023-1024 (zoo2) + 1025-1028 (det2) + 1029-1036 (coccos+phaeo+zoo3)
    expected_tracer_ids(25) = 1023  ! Zoo2N
    expected_tracer_ids(26) = 1024  ! Zoo2C
    expected_tracer_ids(27) = 1025  ! DetZ2N
    expected_tracer_ids(28) = 1026  ! DetZ2C
    expected_tracer_ids(29) = 1027  ! DetZ2Si
    expected_tracer_ids(30) = 1028  ! DetZ2Calc
    expected_tracer_ids(31) = 1029  ! CoccoN
    expected_tracer_ids(32) = 1030  ! CoccoC
    expected_tracer_ids(33) = 1031  ! CoccoChl
    expected_tracer_ids(34) = 1032  ! PhaeoN
    expected_tracer_ids(35) = 1033  ! PhaeoC
    expected_tracer_ids(36) = 1034  ! PhaeoChl
    expected_tracer_ids(37) = 1035  ! Zoo3N
    expected_tracer_ids(38) = 1036  ! Zoo3C
    expected_tracer_ids(39) = 1037  ! DIC remin (added by Sina)
    if (enable_R2OMIP) then
      expected_tracer_ids(40) = 1038  ! DOCt (terrestrial DOC) ! Sina: each number increased by 1
    end if

  else if (enable_coccos .and. .not. enable_3zoo2det) then
    ! Coccos only: 1001-1022 (base) + 1023-1028 (coccos+phaeo)
    expected_tracer_ids(25) = 1023  ! CoccoN
    expected_tracer_ids(26) = 1024  ! CoccoC
    expected_tracer_ids(27) = 1025  ! CoccoChl
    expected_tracer_ids(28) = 1026  ! PhaeoN
    expected_tracer_ids(29) = 1027  ! PhaeoC
    expected_tracer_ids(30) = 1028  ! PhaeoChl
    expected_tracer_ids(31) = 1037  ! DIC remin (added by Sina)
    if (enable_R2OMIP) then
      expected_tracer_ids(32) = 1029  ! DOCt (terrestrial DOC) ! Sina: changed 31 to 32
    end if

  else if (enable_3zoo2det .and. .not. enable_coccos) then
    ! 3Zoo2Det only: 1001-1022 (base) + 1023-1030 (zoo2+det2+zoo3)
    expected_tracer_ids(25) = 1023  ! Zoo2N
    expected_tracer_ids(26) = 1024  ! Zoo2C
    expected_tracer_ids(27) = 1025  ! DetZ2N
    expected_tracer_ids(28) = 1026  ! DetZ2C
    expected_tracer_ids(29) = 1027  ! DetZ2Si
    expected_tracer_ids(30) = 1028  ! DetZ2Calc
    expected_tracer_ids(31) = 1029  ! Zoo3N
    expected_tracer_ids(32) = 1030  ! Zoo3C
    expected_tracer_ids(33) = 1037  ! DIC remin (added by Sina)
    if (enable_R2OMIP) then
      expected_tracer_ids(34) = 1031  ! DOCt (terrestrial DOC) ! Sina: changed 33 to 34 
    end if

  else
    expected_tracer_ids(25) = 1037 ! add DIC remin tracer to base BGC tracers (added by Sina)
    ! Base configuration: only tracers 1, 2, 1001-1022
    if (enable_R2OMIP) then
      expected_tracer_ids(26) = 1023  ! DOCt (terrestrial DOC) ! Sina: changed 25 to 26
    end if
  end if

  if (use_age_tracer) then 
    expected_tracer_ids(num_expected_tracers) = 100
  end if

  ! else: base configuration only needs tracers 1, 2, 1001-1022

  ! ===========================================================================
  ! Perform validation checks
  ! ===========================================================================

  if (mype == 0) then
    write(*,*) ''
    write(*,*) '=========================================================================='
    write(*,*) 'REcoM TRACER CONFIGURATION VALIDATION'
    write(*,*) '=========================================================================='
    write(*,*) 'Model configuration:'
    write(*,*) '  enable_3zoo2det = ', enable_3zoo2det
    write(*,*) '  enable_coccos   = ', enable_coccos
    write(*,*) '  enable_AWICM    = ', enable_AWICM
    write(*,*) '  enable_R2OMIP   = ', enable_R2OMIP
    write(*,*) '  useRivers       = ', useRivers
    write(*,*) ''
    write(*,*) 'Tracer counts:'
    write(*,*) '  Physical tracers (T, S, ...)      = ', num_physical_tracers
    write(*,*) '  DOC tracers (terrestrial input)   = ', doc_tracers
    write(*,*) '  Expected BGC tracers              = ', expected_bgc_num
    write(*,*) '  Expected TOTAL tracers            = ', expected_total_tracers
    write(*,*) '  Actual tracers from namelist      = ', num_tracers
    write(*,*) '  Actual BGC tracers from namelist  = ', actual_bgc_num
    write(*,*) ''
  end if

  ! Check for inconsistencies
  if (actual_bgc_num /= expected_bgc_num) then
    config_error = .true.
    if (mype == 0) then
      write(*,*) '=========================================================================='
      write(*,*) 'ERROR: TRACER COUNT MISMATCH!'
      write(*,*) '=========================================================================='
      write(*,*) 'The number of BGC tracers in the namelist does not match'
      write(*,*) 'the expected count for the current configuration.'
      write(*,*) ''
      write(*,*) '  Expected BGC tracers: ', expected_bgc_num
      write(*,*) '  Actual BGC tracers:   ', actual_bgc_num
      write(*,*) '  Difference:           ', actual_bgc_num - expected_bgc_num
      write(*,*) ''
      write(*,*) 'Required tracer IDs for current configuration:'
      write(*,*) '  Base tracers (always):  1001-1022 (22 tracers) + 1037 (DICremin)'

      if (enable_3zoo2det .and. .not. enable_coccos) then
        write(*,*) '  3Zoo2Det extension:     1023-1030 (8 tracers)'
        write(*,*) '    - Zoo2N, Zoo2C:       1023-1024'
        write(*,*) '    - DetZ2 pool:         1025-1028'
        write(*,*) '    - MicZooN, MicZooC:   1029-1030'
        write(*,*) '    - DICremin:           1037     ' ! added by Sina
        if (enable_R2OMIP) then
          write(*,*) '  Terrestrial DOC:        1031 (1 tracer)'
        end if
      else if (enable_coccos .and. .not. enable_3zoo2det) then
        write(*,*) '  Coccos extension:       1023-1028 (6 tracers)'
        write(*,*) '    - CoccoN, C, Chl:     1023-1025'
        write(*,*) '    - PhaeoN, C, Chl:     1026-1028'
        write(*,*) '    - DICremin:           1037     ' ! added by Sina
        if (enable_R2OMIP) then
          write(*,*) '  Terrestrial DOC:        1029 (1 tracer)'
        end if
      else if (enable_3zoo2det .and. enable_coccos) then
        write(*,*) '    - Zoo2N, Zoo2C:       1023-1024'
        write(*,*) '  3Zoo2Det extension:     1025-1028 (4 tracers for det2)'
        write(*,*) '  Coccos extension:       1029-1034 (6 tracers)'
        write(*,*) '    - CoccoN, C, Chl:     1029-1031'
        write(*,*) '    - PhaeoN, C, Chl:     1032-1034'
        write(*,*) '  MicroZoo extension:     1035-1036 (2 tracers)'
        write(*,*) '    - DICremin:           1037     ' ! added by Sina
        if (enable_R2OMIP) then
          write(*,*) '  Terrestrial DOC:        1038 (1 tracer)' ! Sina: increased from 1037 to 1038
        end if
      else
        if (enable_R2OMIP) then
          write(*,*) '  Terrestrial DOC:        1023 (1 tracer)'
        end if
      end if

      write(*,*) ''
      write(*,*) 'ACTION REQUIRED:'
      write(*,*) '  1. Check your namelist.config tracer_list section'
      write(*,*) '  2. Ensure enable_3zoo2det, enable_coccos, enable_R2OMIP, and useRivers match your setup'
      write(*,*) '  3. Add/remove tracers to match the expected configuration'
      write(*,*) '=========================================================================='
      write(*,*) ''
    end if
  else
    ! Validation passed
    if (mype == 0) then
      write(*,*) '=========================================================================='
      write(*,*) 'VALIDATION PASSED: Tracer configuration is consistent!'
      write(*,*) '=========================================================================='
      write(*,*) ''
    end if
  end if

  ! ===========================================================================
  ! Additional sanity check: verify bgc_num variable matches
  ! ===========================================================================
  if (bgc_num /= expected_bgc_num) then
    if (mype == 0) then
      write(*,*) '=========================================================================='
      write(*,*) 'WARNING: bgc_num variable inconsistency!'
      write(*,*) '=========================================================================='
      write(*,*) 'The bgc_num parameter does not match the expected value.'
      write(*,*) '  Current bgc_num value: ', bgc_num
      write(*,*) '  Expected value:        ', expected_bgc_num
      write(*,*) ''
      write(*,*) 'This may indicate that bgc_num was not updated after changing'
      write(*,*) 'enable_3zoo2det, enable_coccos, enable_R2OMIP, or useRivers flag.'
      write(*,*) '=========================================================================='
      write(*,*) ''
    end if
    config_error = .true.
  end if

  ! ===========================================================================
  ! Validate tracer IDs: Check for correct IDs and detect clashes
  ! ===========================================================================
  id_error = .false.

  ! This check requires access to the actual tracer IDs from the namelist
  ! We'll validate against the expected list
  if (mype == 0) then
    write(*,*) '=========================================================================='
    write(*,*) 'VALIDATING TRACER IDs'
    write(*,*) '=========================================================================='
    write(*,*) 'Expected tracer ID sequence:'
    write(*,*) ''

    ! Display expected IDs in a readable format
    write(*,*) 'Physical tracers:'
    write(*,*) '  ', expected_tracer_ids(1:num_physical_tracers)
    write(*,*) ''
    write(*,*) 'Base BGC tracers (1001-1022):'
    write(*,*) '  ', expected_tracer_ids(3:24)
    write(*,*) ''

    if (expected_bgc_num > 22) then
      write(*,*) 'Extended configuration tracers:'
      write(*,*) '  ', expected_tracer_ids(25:num_expected_tracers)
      write(*,*) ''
    end if

    write(*,*) 'CRITICAL: The tracer IDs in your namelist MUST match this sequence'
    write(*,*) '          exactly, in the same order!'
    write(*,*) ''
    write(*,*) 'Common errors to avoid:'
    write(*,*) '  - Using wrong tracer ID numbers (e.g., 1023 instead of 1025)'
    write(*,*) '  - Tracer ID clashes between configurations'
    write(*,*) '  - Incorrect order of tracer IDs in namelist'
    write(*,*) '  - Missing or duplicate tracer IDs'
    write(*,*) '  - Forgetting to add DOCt tracer when enable_R2OMIP is enabled'
    write(*,*) ''

    ! Configuration-specific warnings
    if (enable_3zoo2det .and. enable_coccos) then
      write(*,*) 'IMPORTANT for FULL MODEL (3zoo2det + coccos):'
    !  write(*,*) '  - Tracers 1023-1024 are NOT used (reserved for other configs)'
      write(*,*) '  - Zoo2 uses:         1023-1024'
      write(*,*) '  - Det2 pool uses:    1025-1028'
      write(*,*) '  - Coccos uses:       1029-1031'
      write(*,*) '  - Phaeocystis uses:  1032-1034'
      write(*,*) '  - Microzooplankton:  1035-1036'
      write(*,*) '  - DIC remin:         1037     ' ! added by Sina
      if (enable_R2OMIP) then
        write(*,*) '  - Terrestrial DOC:   1038'    ! Sina: increased from 1037 to 1038
      end if
      write(*,*) ''
    else if (enable_coccos .and. .not. enable_3zoo2det) then
      write(*,*) 'IMPORTANT for COCCOS-ONLY configuration:'
      write(*,*) '  - Coccos uses:       1023-1025 (NOT 1029-1031)'
      write(*,*) '  - Phaeocystis uses:  1026-1028 (NOT 1032-1034)'
      write(*,*) '  - DIC remin:         1037     ' ! added by Sina
      if (enable_R2OMIP) then
        write(*,*) '  - Terrestrial DOC:   1029'
        write(*,*) '  - Tracers 1030+ are NOT used in this configuration'
      else
        write(*,*) '  - Tracers 1029+ are NOT used in this configuration'
      end if
      write(*,*) ''
    else if (enable_3zoo2det .and. .not. enable_coccos) then
      write(*,*) 'IMPORTANT for 3ZOO2DET-ONLY configuration:'
      write(*,*) '  - Zoo2 uses:         1023-1024'
      write(*,*) '  - Det2 pool uses:    1025-1028'
      write(*,*) '  - Microzoo uses:     1029-1030 (NOT 1035-1036)'
      write(*,*) '  - DIC remin:         1037     ' ! added by Sina
      if (enable_R2OMIP) then
        write(*,*) '  - Terrestrial DOC:   1031'
        write(*,*) '  - Tracers 1032+ are NOT used in this configuration'
      else
        write(*,*) '  - Tracers 1031+ are NOT used in this configuration'
      end if
      write(*,*) ''
    else
      write(*,*) 'IMPORTANT for BASE configuration:'
      write(*,*) '  - Only tracers 1-2, 1001-1022 and 1037 should be present' ! 1037 added by Sina
      if (enable_R2OMIP) then
        write(*,*) '  - Plus terrestrial DOC: 1023'
        write(*,*) '  - Tracers 1024+ are NOT used in base configuration'
      else
        write(*,*) '  - Tracers 1023+ are NOT used in base configuration'
      endif
    end if
    write(*,*) ''
    write(*,*) '=========================================================================='
    write(*,*) ''
  end if

  ! ===========================================================================
  ! Check for tracer ID clashes based on configuration
  ! ===========================================================================
  if (mype == 0) then
    write(*,*) '=========================================================================='
    write(*,*) 'CHECKING FOR TRACER ID CONFLICTS'
    write(*,*) '=========================================================================='

    ! Warn about potential clashes between configurations
    if (enable_3zoo2det .and. enable_coccos) then
      write(*,*) 'Full model configuration active.'
    !  write(*,*) 'Ensure you are NOT using tracer IDs 1023-1024 in your namelist!'
    !  write(*,*) 'These are reserved for configurations WITHOUT full model.'
    else if (enable_coccos) then
      write(*,*) 'Coccos-only configuration active.'
      write(*,*) 'Coccos MUST use IDs 1023-1025 (NOT 1029-1031).'
      write(*,*) 'Phaeocystis MUST use IDs 1026-1028 (NOT 1032-1034).'
    else if (enable_3zoo2det) then
      write(*,*) '3Zoo2Det-only configuration active.'
      write(*,*) 'Microzoo MUST use IDs 1029-1030 (NOT 1035-1036).'
    end if

    if (useRivers) then
      write(*,*) ''
      write(*,*) 'Rivers enabled (useRivers = .true.)'
    end if

    if (enable_R2OMIP) then
      write(*,*) ''
      write(*,*) 'enable_R2OMIP enabled (enable_R2OMIP = .true.)'
      write(*,*) 'Terrestrial DOC tracer is REQUIRED as the last BGC tracer.'
    end if

    write(*,*) ''
   ! write(*,*) 'No automated clash detection available without tracer array access.'
    write(*,*) 'Please manually verify your namelist tracer_list against the'
    write(*,*) 'expected sequence shown above.'
    write(*,*) '=========================================================================='
    write(*,*) ''
  end if

  ! ===========================================================================
  ! Stop execution if configuration error detected
  ! ===========================================================================
  if (config_error) then
    if (mype == 0) then
      write(*,*) ''
      write(*,*) '******************************************************************'
      write(*,*) '***  FATAL ERROR: MODEL CONFIGURATION INCONSISTENCY DETECTED   ***'
      write(*,*) '***  MODEL EXECUTION STOPPED                                   ***'
      write(*,*) '******************************************************************'
      write(*,*) ''
    end if
    deallocate(expected_tracer_ids, tracer_found)
    call par_ex(0)  ! Stop execution (use appropriate stop routine for your model)
    stop
  end if

  ! Clean up
  deallocate(expected_tracer_ids, tracer_found)

end subroutine validate_recom_tracers

! ==============================================================================
! SUBROUTINE: validate_tracer_id_sequence
! ==============================================================================
! Purpose: Validate that actual tracer IDs from namelist match expected sequence
!          Call this after reading the tracer namelist
! ==============================================================================
subroutine validate_tracer_id_sequence(tracer_ids, num_tracers, mype)
  
  use g_forcing_param, only: use_age_tracer !---age-code

  implicit none

  ! Arguments
  integer, dimension(:), intent(in) :: tracer_ids   ! Actual IDs from namelist
  integer, intent(in) :: num_tracers                ! Number of tracers
  integer, intent(in) :: mype                        ! MPI rank

  ! Local variables
  integer :: i, j
  integer, dimension(:), allocatable :: expected_ids
  integer :: num_expected
  logical :: error_found
  logical :: duplicate_found
  integer :: num_physical_tracers
  integer :: doc_tracers

  error_found = .false.
  duplicate_found = .false.
  num_physical_tracers = 2

  ! DOC tracers count
  if (enable_R2OMIP) then
    doc_tracers = 1
  else
    doc_tracers = 0
  end if

  ! Allocate expected IDs array
  allocate(expected_ids(num_tracers))

  ! Build expected ID sequence
  expected_ids(1) = 1
  expected_ids(2) = 2

  do i = 1, 22
    expected_ids(num_physical_tracers + i) = 1000 + i
  end do

  if (enable_3zoo2det .and. enable_coccos) then
    ! Full model configuration
    expected_ids(25:30) = (/1023, 1024, 1025, 1026, 1027, 1028/)
    expected_ids(31:36) = (/1029, 1030, 1031, 1032, 1033, 1034/)
    expected_ids(37:38) = (/1035, 1036/)
    expected_ids(39)    = 1037 ! DICremin, added by Sina
    if (enable_R2OMIP) then
      expected_ids(40) = 1038  ! DOCt ! Sina: each number increased by 1
    end if

  else if (enable_coccos .and. .not. enable_3zoo2det) then
    expected_ids(25:30) = (/1023, 1024, 1025, 1026, 1027, 1028/)
    expected_ids(31)    = 1037 ! DICremin, added by Sina
    if (enable_R2OMIP) then
      expected_ids(32) = 1029  ! DOCt ! Sina: increased 31 to 32
    end if

  else if (enable_3zoo2det .and. .not. enable_coccos) then
    expected_ids(25:32) = (/1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030/)
    expected_ids(33)    = 1037 ! DICremin, added by Sina
    if (enable_R2OMIP) then
      expected_ids(34) = 1031  ! DOCt ! Sina: increased 33 to 34
    end if
    
  else
    ! Base configuration
    expected_ids(25) = 1037    ! DICremin, added by Sina
    if (enable_R2OMIP) then
      expected_ids(26) = 1023  ! DOCt ! Sina: increased 25 to 26
    end if
  end if

  if (use_age_tracer) then 
    expected_ids(num_tracers) = 100
  end if 

  ! ===========================================================================
  ! Check 1: Compare actual vs expected tracer IDs
  ! ===========================================================================
  if (mype == 0) then
    write(*,*) ''
    write(*,*) '=========================================================================='
    write(*,*) 'VALIDATING TRACER ID SEQUENCE FROM NAMELIST'
    write(*,*) '=========================================================================='
  end if

  do i = 1, num_tracers
    if (tracer_ids(i) /= expected_ids(i)) then
      error_found = .true.
      if (mype == 0) then
        write(*,*) 'ERROR at position ', i, ':'
        write(*,*) '  Expected tracer ID: ', expected_ids(i)
        write(*,*) '  Found tracer ID:    ', tracer_ids(i)
        write(*,*) ''
      end if
    end if
  end do

  ! ===========================================================================
  ! Check 2: Detect duplicate tracer IDs
  ! ===========================================================================
  do i = 1, num_tracers - 1
    do j = i + 1, num_tracers
      if (tracer_ids(i) == tracer_ids(j)) then
        duplicate_found = .true.
        if (mype == 0) then
          write(*,*) 'ERROR: Duplicate tracer ID detected!'
          write(*,*) '  Tracer ID ', tracer_ids(i), ' appears at positions ', i, ' and ', j
          write(*,*) ''
        end if
      end if
    end do
  end do

  ! ===========================================================================
  ! Check 3: Detect forbidden tracer IDs for current configuration
  ! ===========================================================================
  !if (enable_3zoo2det .and. enable_coccos) then
    ! Check for forbidden IDs 1023-1024 in full model
    !do i = 1, num_tracers
      !if (tracer_ids(i) == 1023 .or. tracer_ids(i) == 1024) then
        !error_found = .true.
        !if (mype == 0) then
          !write(*,*) 'ERROR: Forbidden tracer ID in full model configuration!'
          !write(*,*) '  Tracer ID ', tracer_ids(i), ' at position ', i
          !write(*,*) '  IDs 1023-1024 are NOT used when both flags are enabled'
          !write(*,*) ''
        !end if
      !end if
    !end do
  !end if

  ! ===========================================================================
  ! Report results
  ! ===========================================================================
  if (error_found .or. duplicate_found) then
    if (mype == 0) then
      write(*,*) '=========================================================================='
      write(*,*) 'TRACER ID VALIDATION FAILED!'
      write(*,*) '=========================================================================='
      write(*,*) ''
      write(*,*) 'Expected tracer ID sequence for current configuration:'
      write(*,*) expected_ids
      write(*,*) ''
      write(*,*) 'Actual tracer ID sequence from namelist:'
      write(*,*) tracer_ids
      write(*,*) ''
      write(*,*) 'ACTION REQUIRED:'
      write(*,*) '  Correct the tracer IDs in your namelist.config file'
      write(*,*) '  Ensure the sequence matches exactly as expected'
      write(*,*) '=========================================================================='
      write(*,*) ''
      write(*,*) '******************************************************************'
      write(*,*) '***  FATAL ERROR: INVALID TRACER ID SEQUENCE                   ***'
      write(*,*) '***  MODEL EXECUTION STOPPED                                   ***'
      write(*,*) '******************************************************************'
      write(*,*) ''
    end if
    deallocate(expected_ids)
    call par_ex(0)
    stop
  else
    if (mype == 0) then
      write(*,*) '=========================================================================='
      write(*,*) 'TRACER ID VALIDATION PASSED!'
      write(*,*) 'All tracer IDs match expected sequence - no clashes detected.'
      write(*,*) '=========================================================================='
      write(*,*) ''
    end if
  end if

  deallocate(expected_ids)

end subroutine validate_tracer_id_sequence

! ==============================================================================
! SUBROUTINE: print_sinking_config
! Purpose: Print the sinking tracer configuration for the current model setup
! ==============================================================================
subroutine print_sinking_config(mype)
    implicit none

    integer, intent(in) :: mype
    integer :: i, n_sinking

    if (mype /= 0) return

    ! Count active sinking tracers (non-zero entries)
    n_sinking = count(recom_sinking_tracer_id /= 0)

    write(*,*)
    write(*,*) '=========================================================================='
    write(*,*) 'REcoM SINKING TRACER CONFIGURATION'
    write(*,*) '=========================================================================='
    write(*,*) 'Model flags:'
    write(*,'(A,L1)') '  enable_3zoo2det : ', enable_3zoo2det
    write(*,'(A,L1)') '  enable_coccos   : ', enable_coccos
    write(*,'(A,L1)') '  useRivers       : ', useRivers
    write(*,*)
    write(*,'(A,I3)') '  Active sinking tracers: ', n_sinking
    write(*,*)

    write(*,*) '----------------------------------------------------------'
    write(*,*) '  Slot | Tracer ID | Description'
    write(*,*) '----------------------------------------------------------'
    do i = 1, 32
        if (recom_sinking_tracer_id(i) == 0) cycle
        write(*,'(A,I4,A,I6,A,A)') '  ', i, '   | ', &
            recom_sinking_tracer_id(i), '    | ', &
            trim(sinking_tracer_name(recom_sinking_tracer_id(i)))
    end do
    write(*,*) '----------------------------------------------------------'

    write(*,*)
    write(*,*) 'Sinking velocity assignments:'
    write(*,'(A,F8.3,A)') '  VDet  (Det1)  : ', VDet,      ' m/day'
    write(*,'(A,F8.3,A)') '  VDet_zoo2     : ', VDet_zoo2, ' m/day'
    write(*,'(A,F8.3,A)') '  VPhy          : ', VPhy,      ' m/day'
    write(*,'(A,F8.3,A)') '  VDia          : ', VDia,      ' m/day'
    if (enable_coccos) then
        write(*,'(A,F8.3,A)') '  VCocco        : ', VCocco,   ' m/day'
        write(*,'(A,F8.3,A)') '  VPhaeo        : ', VPhaeo,   ' m/day'
    end if
    write(*,'(A,L1)')     '  allow_var_sinking : ', allow_var_sinking
    write(*,*)
    write(*,*) '=========================================================================='
    write(*,*)

end subroutine print_sinking_config

! ==============================================================================
! FUNCTION: sinking_tracer_name
! Purpose: Return a human-readable name for a given tracer ID
! ==============================================================================
function sinking_tracer_name(id) result(name)
    implicit none

    integer, intent(in)  :: id
    character(len=24)    :: name

    select case (id)
        ! ------------------------------------------------------------------
        ! Detritus pool 1
        ! ------------------------------------------------------------------
        case (1007); name = 'DetN'
        case (1008); name = 'DetC'
        case (1017); name = 'DetSi'
        case (1021); name = 'DetCalc'
        ! ------------------------------------------------------------------
        ! Small phytoplankton
        ! ------------------------------------------------------------------
        case (1004); name = 'PhyN'
        case (1005); name = 'PhyC'
        case (1020); name = 'PhyCalc'
        case (1006); name = 'PhyChl'
        ! ------------------------------------------------------------------
        ! Diatoms
        ! ------------------------------------------------------------------
        case (1013); name = 'DiaN'
        case (1014); name = 'DiaC'
        case (1016); name = 'DiaSi'
        case (1015); name = 'DiaChl'
        ! ------------------------------------------------------------------
        ! Detritus pool 2 (3zoo2det)
        ! ------------------------------------------------------------------
        case (1025); name = 'DetZ2N'
        case (1026); name = 'DetZ2C'
        case (1027); name = 'DetZ2Si'
        case (1028); name = 'DetZ2Calc'
        ! ------------------------------------------------------------------
        ! Coccolithophores (coccos)
        ! ------------------------------------------------------------------
        case (1029); name = 'CoccoN'
        case (1030); name = 'CoccoC'
        case (1031); name = 'CoccoChl'
        ! ------------------------------------------------------------------
        ! Phaeocystis (coccos)
        ! ------------------------------------------------------------------
        case (1032); name = 'PhaeoN'
        case (1033); name = 'PhaeoC'
        case (1034); name = 'PhaeoChl'
        ! ------------------------------------------------------------------
        ! C-isotope (13C) counterparts
        ! ------------------------------------------------------------------
        case (1305); name = 'PhyCalc_13C'
        case (1308); name = 'DetN_13C'
        case (1314); name = 'DiaSi_13C'
        case (1320); name = 'PhyChl_13C'
        case (1321); name = 'DetC_13C'
        ! ------------------------------------------------------------------
        ! C-isotope (14C) counterparts
        ! ------------------------------------------------------------------
        case (1405); name = 'PhyCalc_14C'
        case (1408); name = 'DetN_14C'
        case (1414); name = 'DiaSi_14C'
        case (1420); name = 'PhyChl_14C'
        case (1421); name = 'DetC_14C'
        ! ------------------------------------------------------------------
        case default
            write(name,'(A,I6)') 'Unknown ID: ', id
    end select

end function sinking_tracer_name

end module recom_config

