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
!!   2019: ongoing
!
module recom_config
  implicit none
  save

!! *** General constants ***

  Integer :: idin = 1, idic = 2, ialk = 3, iphyn = 4, iphyc = 5, &
             ipchl = 6, idetn = 7, idetc = 8, ihetn = 9,         &
             ihetc = 10, idon = 11, idoc = 12, idian = 13,       &
             idiac = 14, idchl = 15, idiasi = 16, idetsi = 17,   &
             isi = 18, ife = 19, iphycal = 20, idetcal = 21,     &
             ioxy= 22, izoo2n = 23, izoo2c =24, idetz2n = 25,   &
             idetz2c = 26, idetz2si = 27, idetz2calc = 28    

  Integer :: ivphy = 1, ivdia = 2, ivdet = 3, ivdetsc = 4

!!MB TEST: tracer ids for revised remineralization and sinking in oce_ale_tracer.F90
  integer, dimension(8)  :: recom_remin_tracer_id   = (/1001, 1002, 1003, 1018, 1019, 1022, 1302, 1402/)
  integer, dimension(26) :: recom_sinking_tracer_id = (/1004, 1005, 1006, 1007, 1008, 1013, 1014, 1015, &
                                                        1016, 1017, 1020, 1021, 1025, 1026, 1027, 1028, &
                                                        1308, 1321, 1305, 1320, & 
                                                        1314, 1408, 1421, 1405, 1420, 1414/)
  integer, dimension(12)  :: recom_det_tracer_id     = (/1007, 1008, 1017, 1021, 1025, 1026, 1027, 1028, &
                                                        1308, 1321, 1408, 1421/)
  integer, dimension(8)  :: recom_phy_tracer_id     = (/1004, 1005, 1020, 1305, 1320, 1405, 1420, 1006/)
  integer, dimension(6)  :: recom_dia_tracer_id     = (/1013, 1014, 1314, 1414, 1016, 1015/)
!!MB

  Real(kind=8)                 :: zero           = 0.d0
  Integer                      :: one            = 1
  Real(kind=8)                 :: tiny           = 2.23D-16
  Real(kind=8)                 :: tiny_chl       = 0.00001 
  Real(kind=8)                 :: SecondsPerDay  = 86400.d0     ! [s/day]
  Real(kind=8)                 :: Pa2atm         = 101325.d0    ! [Pa/atm] 
  Real(kind=8)                 :: redO2C         = 1.453        ! O2:C ratio Anderson and Sarmiento, 1994

!! *** General configuration ***

  Logical                :: use_REcoM            = .true.
  Logical                :: REcoM_restart        = .false.
  logical                :: recom_binary_write   = .false.  ! Determines if tracervalue snapshots are saved. For fine grids it may crash the model to set this to true
  logical                :: recom_binary_init    = .false.  ! Restart from binary
  Integer                :: bgc_num               = 22
  integer                :: bgc_base_num          = 22      ! standard tracers
  Integer                :: diags3d_num           = 2       ! Number of diagnostic 3d tracers to be saved
  Real(kind=8)           :: VDet                  = 20.d0   ! Sinking velocity, constant through the water column and positive downwards
  Real(kind=8)           :: VDet_zoo2             = 200.d0   ! Sinking velocity, constant through the water column 
  Real(kind=8)           :: VPhy                  = 0.d0    !!! If the number of sinking velocities are different from 3, code needs to be changed !!!
  Real(kind=8)           :: VDia                  = 0.d0 
  Logical                :: allow_var_sinking     = .true.   
  Integer                :: biostep               = 1                    ! Number of times biology should be stepped forward for each time step		 
  Logical                :: REcoM_Geider_limiter  = .false.              ! Decides what routine should be used to calculate limiters in sms
  Logical                :: REcoM_Grazing_Variable_Preference = .true.  ! Decides if grazing should have preference for phyN or DiaN
  Logical                :: REcoM_Second_Zoo      = .false.    ! Decides whether having macrozooplankton and second detritus or not
  Logical                :: Grazing_detritus     = .false.    ! Decides grazing on detritus                            
  Logical                :: zoo2_fecal_loss     = .false.    ! Decides fecalloss for the second zooplankton            
  Logical                :: zoo2_initial_field     = .false.    ! Decides initialization of secondzoo        ! NOT CODED YET OG
  Logical                :: het_resp_noredfield    = .true.    ! Decides respiratation of copepods              
  Logical                :: diatom_mucus    = .true.    ! Effect of nutrient limitation on the aggregation
  Logical                :: Graz_pref_new    = .true.    ! If it is true Fasham et 1990, otherwise original recom variable preference
  Logical                :: Diags                 = .true.    !!!!!!!!!!!!!!!!!!!!!!Change in recom.F90 Diagnostics -> Diags
  Logical                :: constant_CO2          = .true.
  Logical                :: UseFeDust             = .true.               ! Turns dust input of iron off when set to.false.
  Logical                :: UseDustClim           = .true.
  Logical                :: UseDustClimAlbani     = .false.               ! Use Albani dustclim field (If it is false Mahowald will be used)
  Logical                :: use_Fe2N              = .true.               ! use Fe2N instead of Fe2C, as in MITgcm version
  Logical                :: use_photodamage       = .false.             ! use Alvarez et al (2018) for chlorophyll degradation
  logical                :: HetRespFlux_plus      = .true.     ! More stable computation of zooplankton respiration fluxes adding a small number to HetN
  character(100)         :: REcoMDataPath         = '/work/ollie/jhauck/forcing/core_new_REcoMforcing/'
  logical                :: restore_alkalinity    = .true.
  logical                :: useRivers             = .false.
  logical                :: useRivFe              = .false.    ! river input of Fe
  logical                :: useErosion            = .false.
  logical                :: NitrogenSS            = .false.   ! This one only activates rivers! And in principle denitrification, but denitrification is commented out. When set to true, external sources and sinks of nitrogen are activated (Riverine, aeolian and denitrification)
  logical                :: useAeolianN           = .false.   ! When set to true, aeolian nitrogen deposition is activated
  integer                :: firstyearoffesomcycle = 1948      ! The first year of the actual physical forcing (e.g. JRA-55) used
  integer                :: lastyearoffesomcycle  = 2009      ! Last year of the actual physical forcing used
  integer                :: numofCO2cycles        = 1         ! Number of cycles of the forcing planned 
  integer                :: currentCO2cycle       = 1         ! Which CO2 cycle we are currently running
  Logical                :: DIC_PI                = .true.
  integer                :: Nmocsy                = 1         ! Length of the vector that is passed to mocsy (always one for recom)
  logical                :: recom_debug           =.true.
  logical                :: ciso                  = .false.   ! main switch to enable/disable carbon isotopes (13|14C)
  integer                :: benthos_num           = 4         ! number of sediment tracers = 8 if ciso = .true.
  Logical                :: use_MEDUSA            = .true.    ! main switch for sediment model
  integer                :: sedflx_num            = 0         ! number of sedimentary fluxs from MEDUSA, = 7 if ciso
  Logical                :: add_loopback          = .false.
  real(kind=8)           :: lb_tscale             = 1.d0      ! time scale to balance the burial loss
  integer                :: bottflx_num           = 4         ! number of stored sinking fluxes from the bottom layer, = 6 if C13 and = 8 if C14
  Logical                :: use_atbox             = .false.   ! switch for atmospheric box model for CO2

  namelist /pavariables/ use_REcoM,                       REcoM_restart,         recom_binary_write,      &
                       recom_binary_init,                 bgc_num,               bgc_base_num,            &
                       diags3d_num,                                                                       &
                       VDet,                              VDet_zoo2,                                      &
                       VPhy,                              VDia,                                           &
                       allow_var_sinking,                 biostep,               REcoM_Geider_limiter,    &
                       REcoM_Grazing_Variable_Preference, REcoM_Second_Zoo,      Grazing_detritus,        &
                       zoo2_fecal_loss,                   zoo2_initial_field,    het_resp_noredfield,     &
                       diatom_mucus,                      Graz_pref_new,                                  &
                       Diags      ,                       constant_CO2,                                   &
                       UseFeDust,                         UseDustClim,           UseDustClimAlbani,       &
                       use_Fe2N,                          use_photodamage,       HetRespFlux_plus,        &
                       REcoMDataPath,                     restore_alkalinity,                             &
                       useRivers,                         useErosion,            useRivFe,                &
                       NitrogenSS,                        useAeolianN,           firstyearoffesomcycle,   &
                       lastyearoffesomcycle,              numofCO2cycles,        currentCO2cycle,         &
                       DIC_PI,                            Nmocsy,                recom_debug,             &
                       ciso,                              benthos_num,                                    &
                       use_MEDUSA,                        sedflx_num,            bottflx_num,             &
                       add_loopback,                      lb_tscale,             use_atbox
!!------------------------------------------------------------------------------
!! *** Sinking ***
  Real(kind=8)                 :: Vdet_a         = 0.0288       ! [1/day]
  Real(kind=8)                 :: Vcalc          = 0.0216       ! [1/day] depth dependence of calc_diss                                                                                                                           
  namelist /pasinking/ Vdet_a,  Vcalc
!!------------------------------------------------------------------------------
!! *** Initialization ***
  Real(kind=8)                 :: cPhyN          = 0.2d0
  Real(kind=8)                 :: cHetN          = 0.2d0
  Real(kind=8)                 :: cZoo2N          = 0.2d0

  namelist /painitialization_N/ cPhyN, cHetN, cZoo2N
!!------------------------------------------------------------------------------
!! *** Arrhenius function ***
  Real(kind=8)                 :: recom_Tref     = 288.15d0       ! [K]
  Real(kind=8)                 :: C2K            = 273.15d0       !     Conversion from degrees C to K
  Real(kind=8)                 :: Ae             = 4500.d0        ! [K] Slope of the linear part of the Arrhenius function
  Real(kind=8)                 :: reminSi        = 0.02d0
  namelist /paArrhenius/ recom_Tref, C2K, Ae, reminSi  
!!------------------------------------------------------------------------------
!! *** For limiter function ***
  Real(kind=8)                 :: NMinSlope      = 50.d0 
  Real(kind=8)                 :: SiMinSlope     = 1000.d0
  Real(kind=8)                 :: NCmin          = 0.04d0
  Real(kind=8)                 :: NCmin_d        = 0.04d0
  Real(kind=8)                 :: SiCmin         = 0.04d0
  Real(kind=8)                 :: k_Fe           = 0.04d0
  Real(kind=8)                 :: k_Fe_d         = 0.12d0
  Real(kind=8)                 :: k_si           = 4.d0
  Real(kind=8)                 :: P_cm           = 3.0d0          ! [1/day]   Rate of C-specific photosynthesis
  Real(kind=8)                 :: P_cm_d         = 3.5d0
  namelist /palimiter_function/ NMinSlope, SiMinSlope, NCmin, NCmin_d, SiCmin, k_Fe, k_Fe_d, k_si, P_cm, P_cm_d   
!!------------------------------------------------------------------------------
!! *** For light calculations ***
  Real(kind=8)                 :: k_w            = 0.04d0         ! [1/m]              Light attenuation coefficient
  Real(kind=8)                 :: a_chl          = 0.03d0         ! [1/m * 1/(mg Chl)] Chlorophyll specific attenuation coefficients
  namelist /palight_calculations/ k_w, a_chl   
!!------------------------------------------------------------------------------
!! *** Photosynthesis ***
  Real(kind=8)                 :: alfa           = 0.14d0	  ! [(mmol C*m2)/(mg Chl*W*day)] 
  Real(kind=8)                 :: alfa_d         = 0.19d0         ! An initial slope of the P-I curve
  Real(kind=8)                 :: parFrac        = 0.43d0
  namelist /paphotosynthesis/ alfa, alfa_d, parFrac   
!!------------------------------------------------------------------------------
!! *** Assimilation ***
  Real(kind=8)                 :: V_cm_fact     = 0.7d0          ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)                 :: V_cm_fact_d   = 0.7d0  
  Real(kind=8)                 :: NMaxSlope     = 1000.d0       ! Max slope for limiting function
  Real(kind=8)                 :: SiMaxSlope    = 1000.d0
  Real(kind=8)                 :: NCmax         = 0.2d0         ! [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
  Real(kind=8)                 :: NCmax_d       = 0.2d0
  Real(kind=8)                 :: SiCmax        = 0.8d0
  Real(kind=8)                 :: NCuptakeRatio = 0.2d0         ! [mmol N/mmol C] Maximum uptake ratio of N:C
  Real(kind=8)                 :: NCUptakeRatio_d = 0.2d0
  Real(kind=8)                 :: SiCUptakeRatio= 0.2d0
  Real(kind=8)                 :: k_din         = 0.55d0          ! [mmol N/m3] Half-saturation constant for nitrate uptake
  Real(kind=8)                 :: k_din_d       = 1.0d0
  Real(kind=8)                 :: Chl2N_max     = 3.15d0           ! [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
  Real(kind=8)                 :: Chl2N_max_d   = 4.2d0
  Real(kind=8)                 :: res_phy       = 0.01d0          ! [1/day] Maintenance respiration rate constant
  Real(kind=8)                 :: res_phy_d     = 0.01d0
  Real(kind=8)                 :: biosynth      = 2.33d0          ! [mmol C/mmol N] Cost of biosynthesis
  Real(kind=8)                 :: biosynthSi    = 0.d0
  namelist /paassimilation/ V_cm_fact, V_cm_fact_d, NMaxSlope, SiMaxSlope, NCmax, NCmax_d, SiCmax,        &
                       NCuptakeRatio, NCUptakeRatio_d, SiCUptakeRatio, k_din, k_din_d,                    &
                       Chl2N_max, Chl2N_max_d, res_phy, res_phy_d, biosynth, biosynthSi
!!------------------------------------------------------------------------------
!! *** Iron chemistry ***
  Real(kind=8)                 :: totalligand     = 1.d0        ! [mumol/m3] order 1. Total free ligand
  Real(kind=8)                 :: ligandStabConst = 100.d0      ! [m3/mumol] order 100. Ligand-free iron stability constant
  namelist /pairon_chem/ totalligand, ligandStabConst
!!------------------------------------------------------------------------------
!! *** Zooplankton ***
  Real(kind=8)                 :: graz_max      = 2.4d0           ! [mmol N/(m3 * day)] Maximum grazing loss parameter 
  Real(kind=8)                 :: epsilonr       = 0.35d0          ! [(mmol N)2 /m6] Half saturation constant for grazing loss 
  Real(kind=8)                 :: res_het       = 0.01d0          ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)
  Real(kind=8)                 :: Redfield      = 6.625           ! [mmol C/mmol N] Redfield ratio of C:N = 106:16
  Real(kind=8)                 :: loss_het      = 0.05d0          ! [1/day] Temperature dependent N degradation of extracellular organic N (EON)
  Real(kind=8)                 :: pzDia         = 0.5d0           ! Maximum diatom preference
  Real(kind=8)                 :: sDiaNsq       = 0.d0
  Real(kind=8)                 :: pzPhy         = 1.0d0           ! Maximum nano-phytoplankton preference                                                                                           
  Real(kind=8)                 :: sPhyNsq       = 0.d0
  real(kind=8)                 :: tiny_het       = 1.d-5          ! for more stable computation of HetRespFlux (_plus). Value can be > tiny because HetRespFlux ~ hetC**2.
  namelist /pazooplankton/ graz_max, epsilonr, res_het, Redfield, loss_het, pzDia, sDiaNsq, pzPhy, sPhyNsq, tiny_het
!-------------------------------------------------------------------------------                                   
                                                                                                                         
! Second Zooplankton                                                                                                                                                                                    
                             
  Real(kind=8)                 :: graz_max2      = 0.1d0          ! [mmol N/(m3 * day)] Maximum grazing loss parameter                                                                                        
  Real(kind=8)                 :: epsilon2       = 0.0144d0          ! [(mmol N)2 /m6] Half saturation constant for grazing loss                                                                              
  Real(kind=8)                 :: res_zoo2       = 0.0107d0          ! [1/day] Respiration by heterotrophs and mortality (loss to detritus)                                                            
! Real(kind=8)                 :: recip_res_zoo2 = 1./res_zoo2                                                                                                                                               
  Real(kind=8)                 :: loss_zoo2      = 0.003d0          ! [1/day] Temperature dependent N degradation of extracellular organic N
  Real(kind=8)                 :: fecal_rate_n      = 0.13d0                                                         
  Real(kind=8)                 :: fecal_rate_c      = 0.295d0                                                          
  Real(kind=8)                 :: pzDia2         = 1.d0           ! Maximum diatom preference                                                                                                                 
  Real(kind=8)                 :: sDiaNsq2       = 0.d0
  Real(kind=8)                 :: pzPhy2         = 0.5d0           ! Maximum diatom preference                                                                                                                
  Real(kind=8)                 :: sPhyNsq2       = 0.d0
  Real(kind=8)                 :: pzHet          = 0.8d0           ! Maximum diatom preference                                                                                                               
  Real(kind=8)                 :: sHetNsq        = 0.d0
  Real(kind=8)                 :: t1_zoo2        = 28145.d0           ! Krill temp. function constant1                                                                                                       
  Real(kind=8)                 :: t2_zoo2        = 272.5d0          ! Krill temp. function constant2                                                                                                         
  Real(kind=8)                 :: t3_zoo2        = 105234.d0          ! Krill temp. function constant3                                                                                                      
  Real(kind=8)                 :: t4_zoo2        = 274.15d0          ! Krill temp. function constant3                                                                                              
  namelist /pasecondzooplankton/ graz_max2, epsilon2, res_zoo2, & !recip_res_zoo2 &  
                                 loss_zoo2, fecal_rate_n, fecal_rate_c, pzDia2, pzPhy2, pzHet,        &
                                 sDiaNsq2, sPhyNsq2, sHetNsq,            &
                                 t1_zoo2, t2_zoo2, t3_zoo2, t4_zoo2
!-------------------------------------------------------------------------------                                                                                                                          
!! *** Detritus Grazing Params ***                                                                                                                                                                        
  Real(kind=8)                 :: pzDet         = 1.d0           ! Maximum small detritus prefence by first zooplankton                                                                                   
  Real(kind=8)                 :: sDetNsq       = 0.d0  
  Real(kind=8)                 :: pzDetZ2       = 1.d0           ! Maximum large detritus preference by first zooplankton                                                                                 
  Real(kind=8)                 :: sDetZ2Nsq     = 0.d0
  Real(kind=8)                 :: pzDet2         = 1.d0           ! Maximum small detritus prefence by second zooplankton                                                                                
  Real(kind=8)                 :: sDetNsq2       = 0.d0  
  Real(kind=8)                 :: pzDetZ22       = 1.d0           ! Maximum large detritus preference by second zooplankton                                                                              
  Real(kind=8)                 :: sDetZ2Nsq2     = 0.d0
  namelist /pagrazingdetritus/ pzDet, sDetNsq, pzDetZ2, sDetZ2Nsq, &
                                 pzDet2, sDetNsq2, pzDetZ22, sDetZ2Nsq2
!!------------------------------------------------------------------------------
!! *** Aggregation ***
  Real(kind=8)                 :: agg_PD        = 0.165d0          ! [m3/(mmol N * day)] Maximum aggregation loss parameter for DetN
  Real(kind=8)                 :: agg_PP        = 0.015d0          ! [m3/(mmol N * day)] Maximum aggregation loss parameter for PhyN and DiaN (plankton)
  namelist /paaggregation/ agg_PD, agg_PP
!!------------------------------------------------------------------------------
!! *** DIN ***
  Real(kind=8)                 :: rho_N         = 0.11d0           ! [1/day] Temperature dependent N degradation of extracellular organic N (EON) (Remineralization of DON)
  namelist /padin_rho_N/ rho_N
!!------------------------------------------------------------------------------
!! *** DIC ***
  Real(kind=8)                 :: rho_C1         = 0.1d0           ! [1/day] Temperature dependent C degradation of extracellular organic C (EOC)
  namelist /padic_rho_C1/ rho_C1
!!------------------------------------------------------------------------------
!! *** Phytoplankton N ***
  Real(kind=8)                 :: lossN         = 0.05d0          ! [1/day] Phytoplankton loss of organic N compounds
  Real(kind=8)                 :: lossN_d       = 0.05d0
  namelist /paphytoplankton_N/ lossN, lossN_d
!!------------------------------------------------------------------------------
!! *** Phytoplankton C ***
  Real(kind=8)                 :: lossC         = 0.10d0          ! [1/day] Phytoplankton loss of carbon 
  Real(kind=8)                 :: lossC_d       = 0.10d0
  namelist /paphytoplankton_C/ lossC, lossC_d
!!------------------------------------------------------------------------------
!! *** Phytoplankton ChlA ***
  Real(8)                      :: deg_Chl       = 0.25d0        ! [1/day]
  Real(kind=8)                 :: deg_Chl_d     = 0.25d0
  namelist /paphytoplankton_ChlA/ deg_Chl, deg_Chl_d
!!------------------------------------------------------------------------------
!! *** Detritus N ***
  Real(kind=8)                 :: grazEff       = 0.4d0         ! [] Grazing efficiency (fraction of grazing flux into zooplankton pool) 
  Real(kind=8)                 :: grazEff2      = 0.8d0         ! [] Grazing efficiency (fraction of grazing flux into second zooplankton pool)
  Real(kind=8)                 :: reminN        = 0.165d0        ! [1/day] Temperature dependent remineralisation rate of detritus	
  namelist /padetritus_N/ grazEff, grazEff2, reminN
!!------------------------------------------------------------------------------
!! *** Detritus C ***
  Real(kind=8)                 :: reminC        = 0.15d0        ! [1/day] Temperature dependent remineralisation rate of detritus
  Real(kind=8)                 :: rho_c2        = 0.1d0        ! [1/day] Temperature dependent C degradation of TEP-C
  namelist /padetritus_C/ reminC, rho_c2
!!------------------------------------------------------------------------------
!! *** Heterotrophs ***
  Real(kind=8)                 :: lossN_z       = 0.15d0
  Real(kind=8)                 :: lossC_z       = 0.15d0
  namelist /paheterotrophs/ lossN_z, lossC_z
!!------------------------------------------------------------------------------
! Second Zooplankton                                                                                                                                                                                        
                          
  Real(kind=8)                 :: lossN_z2       = 0.02d0
  Real(kind=8)                 :: lossC_z2       = 0.02d0
  namelist /paseczooloss/ lossN_z2, lossC_z2
!!------------------------------------------------------------------------------
!! *** Iron ***
!! only Fe2C or Fe2N is used, but I am not allowed to introduce an if-statement here
  Real(kind=8)                 :: Fe2N          = 0.033d0       ! Fe2C * 6.625
  Real(kind=8)                 :: Fe2N_benthos  = 0.15d0        ! test, default was 0.14 Fe2C_benthos * 6.625 - will have to be tuned. [umol/m2/day]
  Real(kind=8)                 :: Fe2C          = 0.005d0
  Real(kind=8)                 :: Fe2C_benthos  = 0.02125       !0.68d0/32.d0       ! [umol/m2/day]
  Real(kind=8)                 :: kScavFe       = 0.07d0
  Real(kind=8)                 :: dust_sol      = 0.02d0        !Dissolution of Dust for bioavaliable
  Real(kind=8)                 :: RiverFeConc   = 1000d0        ! mean DFe concentration in rivers   
  namelist /pairon/ Fe2N, Fe2N_benthos, Fe2C, Fe2C_benthos, kScavFe,    &
                    dust_sol, RiverFeConc
!!------------------------------------------------------------------------------
!! *** Calcification ***
  Real(kind=8)                 :: calc_prod_ratio = 0.02d0
  Real(kind=8)                 :: calc_diss_guts  = 0.0d0
  Real(kind=8)                 :: calc_diss_rate  = 0.005714d0   !20.d0/3500.d0
  Real(kind=8)                 :: calc_diss_rate2 = 0.005714d0
  namelist /pacalc/ calc_prod_ratio, calc_diss_guts, calc_diss_rate, calc_diss_rate2
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

end module recom_config
!
!===============================================================================
!
Module REcoM_declarations
  implicit none
  save

  Integer       :: save_count_recom
  Real(kind=8)  :: tiny_N                 ! Min PhyN
  Real(kind=8)  :: tiny_N_d               ! Min DiaN
  Real(kind=8)  :: tiny_C                 ! Min PhyC
  Real(kind=8)  :: tiny_C_d               ! Min DiaC
  Real(kind=8)  :: tiny_Si                ! Min DiaSi
!-------------------------------------------------------------------------------
! Temperature dependence of rates
  Real(kind=8)  :: rTref                  ! [1/K] Reciproque value of reference temp for Arrhenius function
  Real(kind=8)  :: rTloc                  ! [1/K] Reciproque of local ocean temp
  Real(kind=8)  :: arrFunc                ! []    Temp dependence of rates 
!  if (REcoM_Second_Zoo) then
  Real(kind=8)  :: arrFuncZoo2           ! []    Temperature function for krill     
!  endif
  Real(kind=8)  :: reminSiT
!-------------------------------------------------------------------------------
! Quotas
  Real(kind=8)  :: quota, quota_dia           ! [mmol N/mmol C]  Quota between phytoplankton N and C
  Real(kind=8)  :: recipQuota, recipQuota_dia ! [mmol C/mmol N]  Reciproque of 'quota'
  Real(kind=8)  :: Chl2C, Chl2C_dia           ! [mg ChlA/mmol C] Quota between phytoplankton ChlA and C
  Real(kind=8)  :: Chl2C_plast, Chl2C_plast_dia   ! [mg ChlA/mmol C] needed for photodamage
  Real(kind=8)  :: Chl2N, Chl2N_dia           ! [mg ChlA/mmol N] Quota between phytoplankton ChlA and N
  Real(kind=8)  :: qSiC
  Real(kind=8)  :: qSiN
  Real(kind=8)  :: recipQZoo                  ! [mmol C/mmol N]  Quota between heterotrophic C and N 
!  if (REcoM_Second_Zoo) then
  Real(kind=8)  :: recipQZoo2                  ! [mmol C/mmol N]  Quota between second zoo  C and N                                                                                          
!  endif
!!! Grazing detritus Quotas for converting                                                                                                                                                                
  Real(kind=8)  :: recipDet                   ! [mmol C/mmol N]  Quota between second zoo  C and N               
  Real(kind=8)  :: recipDet2                  ! [mmol C/mmol N]  Quota between second zoo  C and N 
!-------------------------------------------------------------------------------
! For limiter function
  Real(kind=8)          :: qlimitFac, qlimitFacTmp !            Factor that regulates photosynthesis
  Real(kind=8),external :: recom_limiter           !            Function calculating qlimitFac
  Real(kind=8)          :: FeLimitFac              ! [Mumol/m3] Half sat constant for iron  
  Real(kind=8)          :: pMax, pMax_dia          ! [1/day]    Maximum rate of C-specific photosynthesis 
!-------------------------------------------------------------------------------
! Light
  Real(kind=8)  :: kappar                  ! [1/m]  Light attenuation coefficient modified by chla
  Real(kind=8)  :: kappastar              ! []
  Real(kind=8)  :: kdzUpper,kdzLower      ! []     light attenuation * deltaZ at lower and upper control volume border
  Real(kind=8)  :: chl_upper,chl_lower    ! [mg/m3]     chl  at lower and upper control volume border
  Real(kind=8)  :: Chlave                 ! [mg/m3]     vertical average chl between two nodes
  Real(kind=8)  :: Upperlight, Lowerlight ! [?]    light at upper and lower border of control volume
  Real(kind=8)  :: PARave                 ! [?]    Average light in the control volumes
!-------------------------------------------------------------------------------
! Photosynthesis
  Real(kind=8)  :: Cphot, Cphot_dia       ! [1/day] C-specific rate of photosynthesis
!-------------------------------------------------------------------------------
! Assimilation
  Real(kind=8)  :: V_cm                   ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)  :: limitFacN,limitFacN_dia! Factor that regulates N-assimilation. Calc from function recom_limiter
  Real(kind=8)  :: limitFacSi
  Real(kind=8)  :: N_assim, N_assim_dia   ! [mmol N/(mmol C * day)] C specific N utilization rate
  Real(kind=8)  :: Si_assim
!-------------------------------------------------------------------------------
! Chlorophyll
  Real(kind=8)  :: ChlSynth, ChlSynth_dia           ! [mg CHL/ mmol N] CHL a synthesis regulation term 
  Real(kind=8)  :: phyRespRate, phyRespRate_dia     ! [1/day] Phytoplankton respiration rate
  Real(kind=8)  :: KOchl, KOchl_dia                 ! coefficient for damage to the photosynthetic apparatus 
!-------------------------------------------------------------------------------
! Iron chemistry
  Real(kind=8),external :: iron_chemistry 
!-------------------------------------------------------------------------------
! Zooplankton
  Real(kind=8)  :: DiaNsq  
  Real(kind=8)  :: varpzdia, fDiaN                  ! Part of Diatoms available for food
  Real(kind=8)  :: PhyNsq  
  Real(kind=8)  :: varpzPhy, fPhyN                  ! Part of Nano available for food
  Real(kind=8)  :: food, foodsq                     ! [(mmol N)2/m6]
  Real(kind=8)  :: grazingFlux_phy, grazingFlux_Dia ! [mmol N / (m3 * day)]
  Real(kind=8)  :: grazingFlux
  Real(kind=8)  :: HetRespFlux                      ! Zooplankton respiration
  Real(kind=8)  :: HetLossFlux                      ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
!-------------------------------------------------------------------------------
!  if (REcoM_Second_Zoo) then
! Second Zooplankton                                                                                           
       
  Real(kind=8)  :: DiaNsq2, PhyNsq2, HetNsq   
  Real(kind=8)  :: varpzDia2, fDiaN2, varpzPhy2, fPhyN2, varpzHet, fHetN ! Part of Diatoms available for food   
  Real(kind=8)  :: food2, foodsq2                     ! [(mmol N)2/m6]                                   
  Real(kind=8)  :: grazingFlux_phy2, grazingFlux_Dia2, grazingFlux_het2 ! [mmol N / (m3 * day)]        
  Real(kind=8)  :: grazingFlux2
  Real(kind=8)  :: Zoo2RespFlux     ! Zooplankton respiration                   
  Real(kind=8)  :: Zoo2LossFlux    ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)  
  Real(kind=8)  :: Zoo2fecalloss_n    ! [(mmol N)/(m3*day)] Second zoo fecal pellet                        
  Real(kind=8)  :: Zoo2fecalloss_c    ! [(mmol N)/(m3*day)] Second zoo fecal pellet                               
  Real(kind=8)  :: recip_res_zoo22 
!  endif
! Grazing Detritus                                                                                                                                                                                       \
  Real(kind=8)  :: DetNsq, DetZ2Nsq, DetNsq2, DetZ2Nsq2  
  Real(kind=8)  :: varpzDet, varpzDetZ2, varpzDet2, varpzDetZ22         ! Part of Diatoms available for food                                                                                              
  Real(kind=8)  :: fDetN, fDetZ2N, fDetN2, fDetZ2N2                                                                                                                         
  Real(kind=8)  :: grazingFlux_Det, grazingFlux_DetZ2 ! [mmol N / (m3 * day)]                                                                                                                            
  Real(kind=8)  :: grazingFlux_Det2, grazingFlux_DetZ22 ! [mmol N / (m3 * day)]                                                                                                                          \
!  endif
!---------------------------------------------------------------------------------
! Aggregation
  Real(kind=8)  :: AggregationRate                  ! [1/day] AggregationRate (of nitrogen)
!-------------------------------------------------------------------------------
! Calcification
  Real(kind=8)  :: calcification
  Real(kind=8)  :: calc_loss_agg
  Real(kind=8)  :: calc_loss_gra
  Real(kind=8)  :: calc_diss
  Real(kind=8)  :: calc_loss_gra2 !zoo2 detritus
  Real(kind=8)  :: calc_diss2     !zoo2 detritus
!-------------------------------------------------------------------------------
! Diagnostics
  Real(kind=8)  :: recipbiostep                         ! 1/number of steps per recom cycle
  Real(kind=8),allocatable,dimension(:,:) :: Diags3Dloc
!-------------------------------------------------------------------------------
! Benthos
  Real(kind=8),allocatable,dimension(:) :: decayBenthos ! [1/day] Decay rate of detritus in the benthic layer
  Real(kind=8),allocatable,dimension(:) :: wFluxDet     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
  Real(kind=8),allocatable,dimension(:) :: wFluxPhy     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
  Real(kind=8),allocatable,dimension(:) :: wFluxDia     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms 	
  Real(kind=8)              :: Vben_det     ! [m/day] speed of sinking into benthos from water column
  Real(kind=8)              :: Vben_det_seczoo !second zooplankton sinking benthos  
  Real(kind=8)              :: Vben_phy
  Real(kind=8)              :: Vben_dia 
  Real(kind=8)              :: Ironflux     ! [umol Fe/(m2*day)] Flux of Fe from sediment to water
!_______________________________________________________________________________
! Arrays added for RECOM implementation:
!!---- PAR
!real(kind=8),allocatable,dimension(:)     :: PAR

! --> multiplication factor for surface boundary condition in 
!     bc_surface for river and erosion
!     river on/off -->=1.0/0.0 
!     erosion on/off -->=1.0/0.0 

real(kind=8)                               :: is_riverinput
real(kind=8)                               :: is_erosioninput

end module REcoM_declarations

!===============================================================================
! For arrays needed for the whole 2D or 3D domain, but only needed in REcoM
!-------------------------------------------------------------------------------
Module REcoM_GloVar
  implicit none
  save

  Real(kind=8),allocatable,dimension(:,:)   :: Benthos        ! 4 types of benthos-tracers with size [4 n2d]
  Real(kind=8),allocatable,dimension(:,:,:) :: Benthos_tr     ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel

  Real(kind=8),allocatable,dimension(:)   :: GloFeDust        ! [umol/m2/s] Monthly 2D field of iron soluted in surface water from dust
  Real(kind=8),allocatable,dimension(:)   :: GloNDust         ! [mmol/m2/s] 10-year mean 2D fields of nitrogen soluted in surface water from dust
  Real(kind=8),dimension(12)              :: AtmCO2           ! [uatm] Atmospheric CO2 partial pressure. One value for the whole planet for each month

  Real(kind=8),allocatable,dimension(:)   :: AtmFeInput       ! [umol/m2/s] Includes ice, but is, other than that identlical to GloFeDust
  Real(kind=8),allocatable,dimension(:)   :: AtmNInput        ! [umol/m2/s] Includes ice, but is, other than that identlical to GloNDust
  Real(kind=8),allocatable,dimension(:)   :: GloPCO2surf      ! [uatm] Surface ocean CO2 partial pressure
  Real(kind=8),allocatable,dimension(:)   :: GloCO2flux       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloCO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloHplus         ! [mol/kg] Concentrations of H-plus ions in the surface ocean
  Real(kind=8),allocatable,dimension(:)   :: GlodPCO2surf       ! [mmol/m2/day] ocean-atmosphere  
  Real(kind=8),allocatable,dimension(:,:) :: GlodecayBenthos  ! [1/day] Decay rate of detritus in the benthic layer saved for oce_ale_tracer.F90

  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxDet    ! 
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxPhy    ! 
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxDia    ! 

  Real(kind=8),allocatable,dimension(:,:)   :: addtiny
  Real(kind=8),allocatable,dimension(:,:,:) :: Gloaddtiny
  Real(kind=8),allocatable,dimension(:,:,:) :: auxy 

!  Real(kind=8),allocatable,dimension(:,:)   :: GlowFlux         ! 

  Real(kind=8),allocatable,dimension(:,:) :: diags2D          ! Diagnostics in 2D [8 n2d]
  Real(kind=8),allocatable,dimension(:,:,:) :: diags3D          ! Diagnostics in 3D [2 nl-1 n2d]
  Real(kind=8),allocatable,dimension(:)   :: DenitBen         ! Benthic denitrification Field in 2D [n2d 1]

!  for using MEDUSA
  Real(kind=8),allocatable,dimension(:,:)   :: SinkFlx        ! Diagnostics in 2D [4 n2d] or [6 n2d] with ciso
  Real(kind=8),allocatable,dimension(:,:,:) :: SinkFlx_tr     ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel

  Real(kind=8),allocatable,dimension(:,:) :: GloSed           ! Yearly input into bottom water from sediments [n2d 5] or [n2d 7] with ciso
  Real(kind=8),allocatable,dimension(:,:) :: lb_flux            ! Yearly burial from medusa: [n2d 5] or [n2d 9] with ciso_14 

! atmospheric box model:
  real(kind=8),allocatable,dimension(:)   :: x_co2atm         ! atmospheric CO2 mixing ratio (mole fraction)

  Real(kind=8), allocatable,dimension(:)  :: Alk_surf         ! Surface alkalinity field used for restoring
  Real(kind=8), allocatable,dimension(:)  :: relax_alk
  Real(kind=8), allocatable,dimension(:)  :: virtual_alk
  real(kind=8), allocatable,dimension(:,:)  :: PAR3D            ! Light in the water column [nl-1 n2d]

! riverine input
  real(kind=8), allocatable,dimension(:)  :: RiverineLonOrig, RiverineLatOrig, RiverineDINOrig, RiverineDONOrig, RiverineDOCOrig, RiverineDSiOrig ! Variables to save original values for riverine nutrients
  real(kind=8), allocatable,dimension(:)  :: RiverDIN2D, RiverDON2D, RiverDOC2D, RiverDSi2D, RiverAlk2D, RiverDIC2D, RiverFe
  real(kind=8), allocatable,dimension(:)  :: ErosionTSi2D, ErosionTON2D, ErosionTOC2D

!! Cobeta, Cos(Angle of incidence)
  Real(kind=8), allocatable,dimension(:)  ::  cosAI
end module REcoM_GloVar

!===============================================================================
! For variables saved locally for each column and then used in REcoM
!-------------------------------------------------------------------------------
Module REcoM_locVar

  Real(kind=8),allocatable,dimension(:) :: LocBenthos ! Storing the values for benthos in current watercolumn: N,C,Si and Calc
  Real(kind=8) :: Hplus                     ! [mol/kg] Concentrations of H-plus ions in the surface node
  Real(kind=8) :: pCO2surf(1)               ! [uatm] Partial pressure of CO2 in surface layer at current 2D node	
  Real(kind=8) :: dflux(1)                  ! [mmol/m2/day] Flux of CO2 into the ocean
  Real(kind=8) :: o2ex(1)                   ! [mmol/m2/s] Flux of O2 into the ocean
  Real(kind=8) :: ULoc(1)                   ! Wind strength above current 2D node, change array size if used with mocsy input vector longer than one
  Real(kind=8) :: dpCO2surf(1)              ! [uatm] difference of oceanic pCO2 minus atmospheric pCO2

! mocsy output -----------------------------------------------------------------------------------------------------------------------------
  Real(kind=8) :: co2flux(1)                   ! air-to-sea flux of CO2 [mol/(m^2 * s)]
  Real(kind=8) :: co2ex(1)                     ! time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
  Real(kind=8) :: dpco2(1)                     ! difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
  Real(kind=8) :: ph(1)                        ! pH on total scale
  Real(kind=8) :: pco2(1)                      ! oceanic partial pressure of CO2 (uatm)
  Real(kind=8) :: fco2(1)                      ! oceanic fugacity of CO2 (uatm)
  Real(kind=8) :: co2(1)                       ! aqueous CO2 concentration [mol/m^3]
  Real(kind=8) :: hco3(1)                      ! bicarbonate (HCO3-) concentration [mol/m^3]
  Real(kind=8) :: co3(1)                       ! carbonate (CO3--) concentration [mol/m^3]
  Real(kind=8) :: OmegaA(1)                    ! Omega for aragonite, i.e., the aragonite saturation state
  Real(kind=8) :: OmegaC(1)                    ! Omega for calcite, i.e., the   calcite saturation state
  Real(kind=8) :: BetaD(1)                     ! BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  Real(kind=8) :: rhoSW(1)                     ! rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  Real(kind=8) :: p(1)                         ! pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  Real(kind=8) :: tempis(1)                    ! in-situ temperature [degrees C]
  Real(kind=8) :: kw660(1)                     ! gas transfer velocity (piston velocity) for CO2 [m/s] 
  Real(kind=8) :: co2flux_seaicemask(1)        ! air-to-sea flux of CO2 [mmol/m2/s]
  Real(kind=8) :: o2flux_seaicemask(1)        ! air-to-sea flux of CO2 [mmol/m2/s]
!-------------------------------------------------------------------------------

  Real(kind=8) :: bt, dic_molal, talk_molal ! Common block: Species
  Real(kind=8) :: k1, k2, kw, kb, ff        ! Common block: Equilibrium_constants
  Real(kind=8) :: FeDust                    ! [umol/m2/s]
  Real(kind=8) :: NDust                     ! [mmol/m2/s]
  Real(kind=8) :: Loc_ice_conc(1)           ! Used to calculate flux of DIC in REcoM 0 -> 1
  Real(kind=8) :: LocAtmCO2(1)              ! [uatm]
  Real(kind=8) :: LocDiags2D(8)
!  Real(kind=8) :: LocDenit
  Real(kind=8) :: LocRiverDIN, LocRiverDON, LocRiverDOC, LocRiverDSi, LocRiverDIC, LocRiverAlk
  Real(kind=8),allocatable,dimension(:) :: Loclbflux ! riverine input from loopback flux
!  if (REcoM_Second_Zoo) then
  Real(kind=8) :: res_zoo2_a, res_zoo2_f
  Real(kind=8) :: grazingFluxcarbonzoo2                      ! grazingfluxcarbon 
!  endif
  Integer      :: currentCO2year
	
end module REcoM_LocVar
!===============================================================================
! Specific declarations related to carbon isotope simulations
!-------------------------------------------------------------------------------
module REcoM_ciso
  implicit none
  save


! Options for carbon isotope simulations (see namelist.recom)
  logical :: ciso_init             = .false.  ! Initial fractionation of bulk organic matter
  logical :: ciso_14               = .false.  ! Include radiocarbon (-> 31 or 38 tracers)
  logical :: ciso_organic_14       = .false.  ! Include organic radiocarbon (-> 38 tracers)
  real(kind=8) :: delta_co2_13        = -6.61
  real(kind=8) :: big_delta_co2_14(3) = (/0., 0., 0./) 
  real(kind=8) :: lambda_14 = 3.8561e-12      ! Decay constant of carbon-14
! for revised atbox 14CO2 implementation
  logical      :: atbox_spinup    = .true.
  real(kind=8) :: cosmic_14_init  = 2.0       ! Initial 14C production flux (atoms / s / cm**2)


  namelist / paciso / ciso_init, ciso_14, ciso_organic_14, &
                      lambda_14, delta_co2_13, big_delta_co2_14, &
                      atbox_spinup, cosmic_14_init

! Extensions of other modules or subroutines
! Module REcoM_constants: ciso tracer indices
  integer :: idic_13, iphyc_13, idetc_13, ihetc_13, idoc_13, idiac_13, iphycal_13, idetcal_13, &
             idic_14, iphyc_14, idetc_14, ihetc_14, idoc_14, idiac_14, iphycal_14, idetcal_14

! Module REcoM_declarations:
  real(kind=8)  :: quota_13, quota_14, quota_dia_13, quota_dia_14,                     &  ! quotas
                   recipQuota_13, recipQuota_14, recipQuota_dia_13, recipQuota_dia_14, &  ! reciprocal quotas
                   recipQZoo_13, recipQZoo_14
  real(kind=8)  :: HetRespFlux_13, HetRespFlux_14           ! zooplankton respiration fluxes
  real(kind=8)  :: calcification_13, calcification_14, &    ! calcification
                   calc_loss_agg_13, calc_loss_agg_14, &
                   calc_loss_gra_13, calc_loss_gra_14, &
                   calc_diss_13, calc_diss_14

! Module REcoM_GloVar:
  real(kind=8),dimension(12)              :: AtmCO2_13                          ! [uatm] Atmospheric 13CO2 partial pressure. One value for the whole planet for each month
  real(kind=8),dimension(3,12)            :: AtmCO2_14                          ! [uatm] Atmospheric 14CO2 partial pressure. Three latitude zones for each month
  real(kind=8),allocatable,dimension(:)   :: GloPCO2surf_13, GloPCO2surf_14     ! [uatm] Surface ocean 13|14CO2 partial pressure
  real(kind=8),allocatable,dimension(:)   :: GloCO2flux_13, GloCO2flux_14       ! [mmol/m2/day] Positive downwards
  real(kind=8),allocatable,dimension(:)   :: GloCO2flux_seaicemask_13, GloCO2flux_seaicemask_14
  real(kind=8),allocatable,dimension(:)   :: RiverineDOCOrig_13, RiverineDOCOrig_14, RiverDOC2D_13, RiverDOC2D_14 

! Module REcoM_LocVar:
  real(kind=8) :: pCO2surf_13(1), pCO2surf_14(1), & ! [uatm] Partial pressure of 13|14CO2 in surface layer at current 2D node
                  co2flux_13(1), co2flux_14(1), &      ! mocsy output: air-to-sea flux of 13|14CO2 [mol/(m^2 * s)]
                  co2flux_seaicemask_13(1), co2flux_seaicemask_14(1) ! air-to-sea flux of CO2 [mmol/m2/s]
  real(kind=8) :: LocAtmCO2_13(1), LocAtmCO2_14(1)  ! [uatm]
  real(kind=8) :: LocRiverDOC_13, LocRiverDOC_14    ! CHECK

! Subroutines REcoM_main & REcoM_extra:
  real(kind=8) :: lat_val                           ! nodal latitude (of atmospheric input)

! Subroutine REcoM_extra:
  real(kind=8) :: delta_co2_14                      ! atmospheric Delta14CO2

! Subroutine REcoM_forcing:
  real(kind=8) :: co2sat,                         & ! dissolved CO2 at saturation (CO2*air) [mol / m**3]
                  kwco2                             ! piston velocity of CO2

! Subroutine REcoM_sms:
  real(kind=8) :: DIC_13, DIC_14,                 & ! [mmol/m3] Dissolved Inorganic 13|14Carbon
                  PhyC_13, PhyC_14,               & ! [mmol/m3] Intracellular conc of 13|14Carbon in small phytoplankton
                  DetC_13, DetC_14,               & ! [mmol/m3] Conc of 13|14C in Detritus
                  HetC_13, HetC_14,               & ! [mmol/m3] Conc of 13|14C in heterotrophs
                  EOC_13, EOC_14,                 & ! [mmol/m3] Extracellular Organic 13|14C conc
                  DiaC_13, DiaC_14,               & ! [mmol/m3] Intracellular conc of 13|14Carbon in diatoms
                  PhyCalc_13, PhyCalc_14,         & ! [mmol/m3] Conc of 13|14C in calcite of phytoplankton
                  DetCalc_13, DetCalc_14            ! [mmol/m3] Conc of 13|14C in calcite of detritus
  real(kind=8),allocatable,dimension(:) :: Cphot_z, Cphot_dia_z ! Vertical profiles of photosynthesis rates

! Subroutine REcoM_init:
  real(kind=8),allocatable,dimension(:,:) :: delta_dic_13_init, &     ! auxiliary initial
                                             delta_dic_14_init, &     ! d|Delta13|14C
                                             big_delta_dic_14_init    ! fields

! Atmospheric box model (global variables):
  real(kind=8),allocatable,dimension(:) ::   x_co2atm_13, x_co2atm_14, & ! atmospheric CO2 mixing ratio (mole fraction)
                                             cosmic_14                   ! cosmogenic 14 production (mol / s)
  real(kind=8) :: production_rate_to_flux_14, &                          ! conversion factor
                  r_atm_spinup_13, r_atm_spinup_14                       ! 13|14CO2 / 12CO2 spinup ratios

! Specific factors related the carbon-isotopic composition
! Isotopic ratios
  real(kind=8) :: r_atm_13, r_atm_14,             & ! atmospheric CO2
                  r_co2s_13, r_co2s_14,           & ! dissolved CO2
                  r_dic_13, r_dic_14,             & ! DIC in seawater
                  r_phyc_13, r_phyc_14,           & ! nanophytoplankton
                  r_diac_13, r_diac_14,           & ! diatoms
                  r_iorg_13 = 0.975,              & ! initial ratios of organic matter
                  r_iorg_14 = 0.950
                  
! Fractionation factors
  real(kind=8) :: alpha_k_13 = 0.99912,           & ! gas transfer (kinetic fractionation,
                  alpha_k_14 = 0.99824,           & ! mean values for 5-21C by Zhang et al., 1995)
                  alpha_aq_13, alpha_aq_14,       & ! dissolution of CO2 in sewater (equilibrium fractionation)
                  alpha_dic_13, alpha_dic_14,     & ! hydrolysis / dissociation of CO2 <-> DIC (equilibrium fract.)
                  alpha_p_13, alpha_p_14,         & ! photosynthesis of nanophytoplankton
                  alpha_p_dia_13, alpha_p_dia_14, & ! photosynthesis of diatoms
                  alpha_calc_13 = 1.000,          & ! calcification (Romanek et al., 1992: 1.001, 1.002)
                  alpha_calc_14 = 1.000,          &
                  alpha_dcal_13 = 1.000,          & ! dissolution of calcite (Romanek et al., 1992: 0.999, 0.998)
                  alpha_dcal_14 = 1.000
! Radioactive decay constant of carbon-14
! t1/2 = 5700 years (Bé et al., 2013; recommended by Orr et al., 2017, for OMIP-BGC)
! if 1 year := 365.25 days:  lambda_14 = 3.8534e-12 / second
! if 1 year := 365.00 days:  lambda_14 = 3.8561e-12 / second
! if 1 year := 360    days:  lambda_14 = 3.9096e-12 / second

! Tracer IDs to be considered in decay calculations (oce_ale_tracer.F90)
  integer, dimension(8) :: c14_tracer_id = (/1402, 1405, 1408, 1410, 1412, 1414, 1420, 1421/)

  contains


    subroutine recom_ciso_airsea(tempc, co3, dic)
!   ----------------------------------------------------------------------------------
!     Subroutine to calculate carbon-isotopic fractionation during air-sea exchange
!   ----------------------------------------------------------------------------------
!
!     Input variables:
!     tempc              lokal temperature in C
!     co3                carbonate ion concentration
!     dic                total carbon concentration
!
!     Output variables, defined in module REcoM_ciso:
!     alpha_k_13,14      kinetic fract. factors for gas transfer
!     alpha_aq_13,14     equilib. fract. factors for dissolution
!     alpha_dic_13,14    equilib. fract. factors for DIC <-> CO2
!
!     Internal variables:
!     epsilon_aq_13,14   equilib. fractionation for dissolution
!     epsilon_dic_13,14  equilib. fractionation for DIC <-> CO2
!     fco3               total carbon fraction
!
!     mbutzin, 2016 - 2019.


!     Declarations
      implicit none

      real(kind=8), intent(in) ::  tempc, co3, dic
      real(kind=8) ::  epsilon_aq_13, epsilon_dic_13, fco3

!     Calculation of carbon-isotopic fractionation factors, where
!
!     alpha_xy   = Rx / Ry               = fractionation factor
!     epsilon_xy = (alpha_xy - 1) * 1000 = fractionation (in per mill)
!     epsilon_14 = 2 * epsilon_13 => alpha_14 = 2 * alpha_13 - 1.

!     We use parametrisations and numerical values determined for carbon-13
!     by Zhang et al. (1995).

!     Kinetic fractionation during gas transfer, mean values between 5 and 21C
!     (values are defined in module REcoM_ciso)
!     epsilon_k_13 = -0.86 => alpha_k_13 =  0.99914, alpha_k_14 =  0.99828

!     Equilibrium fractionation during gas dissolution
!
      epsilon_aq_13 = 0.0049 * tempc - 1.31
      alpha_aq_13   = 1. + 0.001 * epsilon_aq_13

!     Equilibrium fractionation between DIC and CO2
!
!     The equilibrium fractionation between DIC and CO2 cannot be simply
!     calculated from the fractionation factors for HCO3, CO3 and CO2star.
!     Here, we employ an empirical function involving fCO3 = [CO3] / DIC
!     assuming that fCO3 is the same for all carbon isotopes
      fco3 = co3 / dic
      epsilon_dic_13 = (0.014 * fco3 - 0.107 ) * tempc + 10.53
      alpha_dic_13 = 1. + 0.001 * epsilon_dic_13

!     Fractionation of radiocarbon
      if (ciso_organic_14) then
        alpha_aq_14  = 2. * alpha_aq_13  - 1.
        alpha_dic_14 = 2. * alpha_dic_13 - 1.
      else
!       no fractionation in the inorganic approximation
        alpha_aq_14  = 1.
        alpha_dic_14 = 1.
      end if

    return
    end subroutine recom_ciso_airsea
!   ----------------------------------------------------------------------------------

    subroutine recom_ciso_photo (co2st)
!   ----------------------------------------------------------------------------------
!        Subroutine calculating carbon-isotopic fractionation during photosynthesis
!   ----------------------------------------------------------------------------------
!     Input:
!     dissolved CO2 (co2st) in mol / m**3
!
!     Output:
!     isotopic fractionation factors for phytoplankton and diatoms due to
!     photosynthesis (alpha_p_13|14, declared at the head of the module)
!
!     Note that we are interested in effective values (implictly including the 
!     fractionation of dissolved CO2) which are actually derived in field studies
!     or lab experiments. Young et al. 2013, eq. (5) with values from paragraph [35]
!
!     Here, we follow Young et al. 2013, eq. (5) with values from paragraph [35]
!     eps_p = eps_pm * (1. - rho / co2aq) = 17.6 * (1 - 2.02 / co2aq)
!     where co2aq is in umol / L
!
!     mbutzin, 2017 - 2021.

      implicit none
      real(kind=8), intent(in):: co2st
      real(kind=8)            :: co2aq

!     Convert dissolved CO2 from mol / m**3 to umol / L and prevent from division by zero
      co2aq = max(1.d-8, co2st * 1000.)

!     Fractionation wrt carbon-13
      alpha_p_13     = max(1., 1. + 0.001 * (17.6 * (1 - 2.02 / co2aq)))
      alpha_p_dia_13 = alpha_p_13 

!     Fractionation wrt carbon-14
      alpha_p_14     = 2. * alpha_p_13 - 1.
      alpha_p_dia_14 = 2. * alpha_p_dia_13 - 1.

    return
    end subroutine recom_ciso_photo
!   ----------------------------------------------------------------------------------  
 
 
    function lat_zone(lat_n)
!   ----------------------------------------------------------------------------------
!   Assign latitude zones from nodal latitude values 
!   ----------------------------------------------------------------------------------
    
      implicit none
      integer                  :: lat_zone
    
!     Input: Latitude value corresponding to node n
      real(kind=8), intent(in) :: lat_n

!     Binning of latitudes to three zones
      if (lat_n > 30.)  then       ! Northern Hemisphere polewards of 30°N
        lat_zone = 1
      else if (lat_n <- 30.) then  ! Southern Hemisphere polewards of 30°S
        lat_zone = 3
      else                         ! (Sub-) Tropical zone
        lat_zone = 2
      end if
    
      return
    end function lat_zone


    function wind_10(windstr_x, windstr_y)
!   ----------------------------------------------------------------------------------
!    computes wind speed at 10 m height "wind10" from wind stress fields tau_x, tau_y
!    as long as wind10 is not properly passed from ECHAM in coupled simulations.
!    We follow Peixoto & Oort (1992, Eq. (10.28), (10,29)) and Charnock (1955); 
!    also see MPI report 349 (2003), Eq. (5.7).
!   ----------------------------------------------------------------------------------
      implicit none
     
      real(kind=8) :: wind_10

!     Input
      real(kind=8), intent(in) :: windstr_x, windstr_y

!     Internal variables and parameters
!     Zonal and meridional velocities at 10 m height
      real(kind=8) :: u_10, v_10
!     Zonal and meridional friction velocities
      real(kind=8) :: u_fric, v_fric
!     Zonal and meridional roughness lengths
      real(kind=8) :: l_rough_x, l_rough_y
!     Inverse von-Karman constant (0.4), Charnock constant (0.018) divided by g, inverse density of air (1.3), log(10)
      real(kind=8), parameter :: inv_karm = 2.5, charn_g = 0.00173, inv_dens_air = 0.76923, log_10 = 2.30258
     
!     Calculate friction velocities (Peixoto & Oort, 1992, Eq. (10.28))
      u_fric = sqrt(abs(windstr_x) * inv_dens_air)
      v_fric = sqrt(abs(windstr_y) * inv_dens_air)

!     Calculate roughness lengths (MPI report 349, 2003, Eq. (5.7), quoting Charnock, 1955)
      l_rough_x = max((charn_g * u_fric**2), 1.5e-5)
      l_rough_y = max((charn_g * v_fric**2), 1.5e-5)

!     Calculate wind speed at 10 m (Peixoto & Oort, 1992, Eq. (10.29))
      u_10 = inv_karm * u_fric * (log_10 - log(l_rough_x))
      v_10 = inv_karm * v_fric * (log_10 - log(l_rough_y))
     
      wind_10 = sqrt(u_10**2 + v_10**2)
     
      return
    end function wind_10
!   ----------------------------------------------------------------------------------
   
   
end module REcoM_ciso

