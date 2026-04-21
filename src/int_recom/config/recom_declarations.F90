Module REcoM_declarations
  implicit none
  save

  Integer       :: save_count_recom
  Real(kind=8)  :: tiny_N                 ! Min PhyN
  Real(kind=8)  :: tiny_N_d               ! Min DiaN
  Real(kind=8)  :: tiny_N_c               ! Min CocN                 ! NEW
  Real(kind=8)  :: tiny_N_p               ! Min PhaN                 ! Phaeocystis
  Real(kind=8)  :: tiny_C                 ! Min PhyC
  Real(kind=8)  :: tiny_C_d               ! Min DiaC
  Real(kind=8)  :: tiny_C_c               ! Min CocC                 ! NEW
  Real(kind=8)  :: tiny_C_p               ! Min PhaC                 ! Phaeocystis
  Real(kind=8)  :: tiny_Si                ! Min DiaSi
!!------------------------------------------------------------------------------
!! *** Temperature dependence of rates ***
  Real(kind=8)  :: rTref                  ! [1/K] Reciproque value of reference temp for Arrhenius function
  Real(kind=8)  :: rTloc                  ! [1/K] Reciproque of local ocean temp
  Real(kind=8)  :: arrFunc                ! []    Temp dependence of rates (also for Phaeocystis)
  Real(kind=8)  :: CoccoTFunc             ! []    Temp dependence of coccolithophores
  Real(kind=8)  :: Temp_diatoms           ! []    Temp dependence of diatoms
  Real(kind=8)  :: Temp_phyto             ! []    Temp dependence of small phyto
  Real(kind=8)  :: Temp_cocco             ! []    Temp dependence of coccolithophores
  Real(kind=8)  :: Temp_phaeo             ! []    Temp dependence of phaeocystis
  Real(kind=8)  :: arrFuncZoo2            ! []    Temperature function for krill
  Real(kind=8)  :: q10_mic                ! 3Zoo
  Real(kind=8)  :: q10_mic_res            ! 3Zoo
  Real(kind=8)  :: q10_mes                ! 3Zoo
  Real(kind=8)  :: q10_mes_res            ! 3Zoo
  Real(kind=8)  :: reminSiT
  Real(kind=8)  :: O2Func                 ! O2remin
!!------------------------------------------------------------------------------
!! *** CO2 dependence of rates ! NEW CO2 ***
  Real(kind=8)  :: h_depth(1)             ! pH from mocsy is converted to proton concentration
  !Real(kind=8)  :: d_CT_CL_phy            ! NEW inter For the interaction term between CO2 and both temperature and light
  !Real(kind=8)  :: d_CT_CL_dia
  !Real(kind=8)  :: d_CT_CL_coc
  Real(kind=8)  :: CoccoCO2
  Real(kind=8)  :: DiaCO2
  Real(kind=8)  :: PhyCO2
  Real(kind=8)  :: PhaeoCO2

!!------------------------------------------------------------------------------
!! *** Quotas ***
  Real(kind=8)  :: quota, quota_dia, quota_cocco, quota_phaeo                         ! [mmol N/mmol C]  Quota between phytoplankton N and C (NEW changed term)
  Real(kind=8)  :: recipQuota, recipQuota_dia, recipQuota_cocco, recipQuota_phaeo     ! [mmol C/mmol N]  Reciproque of 'quota' (NEW changed term)
  Real(kind=8)  :: Chl2C, Chl2C_dia, Chl2C_cocco, Chl2C_phaeo                         ! [mg ChlA/mmol C] Quota between phytoplankton ChlA and C (NEW changed term)
  Real(kind=8)  :: Chl2C_plast, Chl2C_plast_dia, CHL2C_plast_cocco, CHL2C_plast_phaeo ! [mg ChlA/mmol C] needed for photodamage (NEW changed term)
  Real(kind=8)  :: Chl2N, Chl2N_dia, Chl2N_cocco, Chl2N_phaeo                         ! [mg ChlA/mmol N] Quota between phytoplankton ChlA and N (NEW changed term)
  Real(kind=8)  :: qSiC
  Real(kind=8)  :: qSiN
  Real(kind=8)  :: recipQZoo                                         ! [mmol C/mmol N]  Quota between heterotrophic C and N 
  Real(kind=8)  :: recipQZoo2                                        ! [mmol C/mmol N]  Quota between second zoo  C and N
  Real(kind=8)  :: recipQZoo3                                        ! Zoo3 [mmol C/mmol N] Quota between third zoo C and N
!!! Grazing detritus Quotas for converting                                         
  Real(kind=8)  :: recipDet                                          ! [mmol C/mmol N]  Quota between second zoo  C and N
  Real(kind=8)  :: recipDet2                                         ! [mmol C/mmol N]  Quota between second zoo  C and N 

!!------------------------------------------------------------------------------
!! *** For limiter function ***
  Real(kind=8)          :: qlimitFac, qlimitFacTmp                   ! Factor that regulates photosynthesis
!  Real(kind=8),external :: recom_limiter                             ! Function calculating qlimitFac
! recom_limiter is now a module procedure in recom_iron — use recom_iron in callers
  Real(kind=8)          :: FeLimitFac                                ! [Mumol/m3] Half sat constant for iron
  Real(kind=8)          :: pMax, pMax_dia, pMax_cocco, pMax_phaeo    ! [1/day]    Maximum rate of C-specific photosynthesis
!!------------------------------------------------------------------------------
!! *** Light ***
  Real(kind=8)  :: kappar                                             ! [1/m]  Light attenuation coefficient modified by chla
  Real(kind=8)  :: kappastar                                          ! []
  Real(kind=8)  :: kdzUpper,kdzLower                                  ! []     light attenuation * deltaZ at lower and upper control volume border
  Real(kind=8)  :: chl_upper,chl_lower                                ! [mg/m3]     chl  at lower and upper control volume border
  Real(kind=8)  :: Chlave                                             ! [mg/m3]     vertical average chl between two nodes
  Real(kind=8)  :: Upperlight, Lowerlight                             ! [?]    light at upper and lower border of control volume
  Real(kind=8)  :: PARave                                             ! [?]    Average light in the control volumes
!!------------------------------------------------------------------------------
!! *** Photosynthesis ***
  Real(kind=8)  :: Cphot, Cphot_dia, Cphot_cocco, Cphot_phaeo         ! [1/day] C-specific rate of photosynthesis
!!------------------------------------------------------------------------------
!! *** Assimilation ***
  Real(kind=8)  :: V_cm                                               ! scaling factor for temperature dependent maximum of C-specific N-uptake
  Real(kind=8)  :: limitFacN,limitFacN_dia,limitFacN_cocco, limitFacN_phaeo ! Factor that regulates N-assimilation. Calc from function recom_limiter
  Real(kind=8)  :: limitFacSi
  Real(kind=8)  :: N_assim, N_assim_dia, N_assim_Cocco, N_assim_phaeo       ! [mmol N/(mmol C * day)] C specific N utilization rate
  Real(kind=8)  :: Si_assim
!!------------------------------------------------------------------------------
!! *** Chlorophyll ***
  Real(kind=8)  :: ChlSynth, ChlSynth_dia, ChlSynth_cocco, ChlSynth_phaeo             ! [mg CHL/ mmol N] CHL a synthesis regulation term
  Real(kind=8)  :: phyRespRate, phyRespRate_dia, phyRespRate_cocco, phyRespRate_phaeo ! [1/day] Phytoplankton respiration rate
  Real(kind=8)  :: KOchl, KOchl_dia, KOchl_cocco, KOchl_phaeo                         ! coefficient for damage to the photosynthetic apparatus
!!------------------------------------------------------------------------------
!! *** Vertical only Decomposition of phytoplankton growth components ***
  Real(kind=8),allocatable,dimension(:)  :: VTTemp_diatoms, VTTemp_phyto, VTTemp_cocco, VTTemp_phaeo            ! Vertical 1D  temperature effect on phytoplankton photosynthesis
  Real(kind=8),allocatable,dimension(:)  :: VTPhyCO2, VTDiaCO2, VTCoccoCO2, VTPhaeoCO2                        ! CO2 effect
  Real(kind=8),allocatable,dimension(:)  :: VTqlimitFac_phyto, VTqlimitFac_diatoms, VTqlimitFac_cocco, VTqlimitFac_phaeo  ! nutrient effect
  Real(kind=8),allocatable,dimension(:)  :: VTCphotLigLim_phyto, VTCphotLigLim_diatoms, VTCphotLigLim_cocco, VTCphotLigLim_phaeo ! light limitation
  Real(kind=8),allocatable,dimension(:)  :: VTCphot_phyto, VTCphot_diatoms, VTCphot_cocco, VTCphot_phaeo
  Real(kind=8),allocatable,dimension(:)  :: VTSi_assimDia

!!------------------------------------------------------------------------------
!! *** Iron chemistry ***
!  Real(kind=8),external :: iron_chemistry, iron_chemistry_2ligands
! iron_chemistry and iron_chemistry_2ligands are now module procedures in recom_iron
  Real(kind=8)  :: logK1, logK2, Klig1, Klig2
!!------------------------------------------------------------------------------
!! *** Zooplankton ***
  Real(kind=8)  :: DiaNsq  
  Real(kind=8)  :: varpzdia, fDiaN                       ! Part of Diatoms available for food
  Real(kind=8)  :: PhyNsq
  Real(kind=8)  :: varpzPhy, fPhyN                       ! Part of Nano available for food
  Real(kind=8)  :: CoccoNsq
  Real(kind=8)  :: varpzCocco,fCoccoN
  Real(kind=8)  :: PhaeoNsq
  Real(kind=8)  :: varpzPhaeo,fPhaeoN
  Real(kind=8)  :: MicZooNsq                             ! NEW 3Zoo
  Real(kind=8)  :: varpzMicZoo, fMicZooN                 ! NEW 3Zoo Part of microzooplankton available for food 
  Real(kind=8)  :: food, foodsq                          ! [(mmol N)2/m6]
  Real(kind=8)  :: grazingFlux_phy, grazingFlux_Dia, grazingFlux_Cocco, grazingFlux_Phaeo ! [mmol N / (m3 * day)] (NEW changed term)
  Real(kind=8)  :: grazingFlux_miczoo                    ! NEW 3Zoo
  Real(kind=8)  :: grazingFlux
  Real(kind=8)  :: grazEff                               ! NEW 3Zoo
  Real(kind=8)  :: HetRespFlux                           ! Zooplankton respiration
  Real(kind=8)  :: HetLossFlux                           ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
!!------------------------------------------------------------------------------
!! *** Second Zooplankton  ***                                                                                          
     Real(kind=8)  :: DiaNsq2, PhyNsq2, CoccoNsq2, PhaeoNsq2, HetNsq                   ! NEW (changed term)
     Real(kind=8)  :: varpzDia2, fDiaN2, varpzPhy2, fPhyN2, varpzCocco2, fCoccoN2, varpzPhaeo2, fPhaeoN2, varpzHet, fHetN ! Part of Diatoms available for food
     Real(kind=8)  :: MicZooNsq2                         ! NEW Zoo3
     Real(kind=8)  :: varpzMicZoo2, fMicZooN2            ! NEW Zoo3
     Real(kind=8)  :: food2, foodsq2                     ! [(mmol N)2/m6]
     Real(kind=8)  :: grazingFlux_phy2, grazingFlux_Dia2, grazingFlux_Cocco2, grazingFlux_Phaeo2, grazingFlux_het2 ! [mmol N / (m3 * day)  (NEW changed term)
     Real(kind=8)  :: grazingFlux_miczoo2                ! NEW Zoo3
     Real(kind=8)  :: grazingFlux2
     Real(kind=8)  :: Zoo2RespFlux                       ! Zooplankton respiration                   
     Real(kind=8)  :: Zoo2LossFlux                       ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)  
     Real(kind=8)  :: Zoo2fecalloss_n                    ! [(mmol N)/(m3*day)] Second zoo fecal pellet                        
     Real(kind=8)  :: Zoo2fecalloss_c                    ! [(mmol N)/(m3*day)] Second zoo fecal pellet     
     Real(kind=8)  :: Mesfecalloss_n                     ! NEW Zoo3
     Real(kind=8)  :: Mesfecalloss_c                     ! NEW Zoo3                          
     Real(kind=8)  :: recip_res_zoo22 
!!------------------------------------------------------------------------------
!! *** Grazing Detritus  *** 
  Real(kind=8)  :: DetNsq, DetZ2Nsq, DetNsq2, DetZ2Nsq2  
  Real(kind=8)  :: varpzDet, varpzDetZ2, varpzDet2, varpzDetZ22         ! Part of Diatoms available for food
  Real(kind=8)  :: fDetN, fDetZ2N, fDetN2, fDetZ2N2
  Real(kind=8)  :: grazingFlux_Det, grazingFlux_DetZ2                   ! [mmol N / (m3 * day)]
  Real(kind=8)  :: grazingFlux_Det2, grazingFlux_DetZ22                 ! [mmol N / (m3 * day)]
!!------------------------------------------------------------------------------
!! *** Third zooplankton  ***       ! NEW 3Zoo
  Real(kind=8)  :: DiaNsq3
  Real(kind=8)  :: varpzDia3, fDiaN3                 ! Part of diatoms available for food
  Real(kind=8)  :: loss_hetfd
  Real(kind=8)  :: PhyNsq3
  Real(kind=8)  :: varpzPhy3, fPhyN3                 ! Part of small phytoplankton available for food
  Real(kind=8)  :: CoccoNsq3
  Real(kind=8)  :: varpzCocco3, fCoccoN3             ! Part of coccolithophores available for food
  Real(kind=8)  :: PhaeoNsq3
  Real(kind=8)  :: varpzPhaeo3, fPhaeoN3             ! Part of phaeocystis available for food
  Real(kind=8)  :: food3, foodsq3                    ! [(mmol N)2/m6]
  Real(kind=8)  :: grazingFlux_phy3, grazingFlux_Dia3, grazingFlux_Cocco3, grazingFlux_Phaeo3 ! [mmol N / (m3 * day)]
  Real(kind=8)  :: grazingFlux3
  Real(kind=8)  :: MicZooRespFlux                    ! Zooplankton respiration
  Real(kind=8)  :: MicZooLossFlux                    ! [(mmol N)2/(m6 * day)] Zooplankton mortality (quadratic loss)
!!------------------------------------------------------------------------------                                                                                
!! *** Aggregation  ***
  Real(kind=8)  :: AggregationRate                   ! [1/day] AggregationRate (of nitrogen)
!!------------------------------------------------------------------------------                                                                                
!! *** Calcification  ***
  Real(kind=8)  :: calc_prod_ratio_cocco             ! NEW (before it was defined as a fixed value, but now dependent on cocco and T)
  Real(kind=8)  :: calcification
  Real(kind=8)  :: calc_loss_agg
  Real(kind=8)  :: calc_loss_gra
  Real(kind=8)  :: calc_diss
  Real(kind=8)  :: calc_diss_ben                     ! NEW DISS
  Real(kind=8)  :: calc_loss_gra2                    ! zoo2 detritus
  Real(kind=8)  :: calc_diss2                        ! zoo2 detritus
  Real(kind=8)  :: calc_loss_gra3                    ! NEW Zoo3 detritus
  Real(kind=8)  :: Ca                                ! NEW DISS (calcium ion concentration)
  Real(kind=8)  :: CO3_sat                           ! NEW DISS (saturated CO3 concentration, calculated from kspc and Ca)
!!------------------------------------------------------------------------------                                                                                
!! *** Diagnostics  ***
  Real(kind=8)  :: recipbiostep                         ! 1/number of steps per recom cycle
  Real(kind=8),allocatable,dimension(:,:) :: Diags3Dloc

  ! ==================================================================
  ! SMALL PHYTOPLANKTON (n suffix)
  ! ==================================================================
  Real(kind=8)                          :: locNPPn, locGPPn, locNNAn, locChldegn
  Real(kind=8),allocatable,dimension(:) :: vertNPPn, vertGPPn, vertNNAn, vertChldegn
  Real(kind=8),allocatable,dimension(:) :: vertrespn
  Real(kind=8),allocatable,dimension(:) :: vertdocexn
  Real(kind=8),allocatable,dimension(:) :: vertaggn

  ! ==================================================================
  ! DIATOMS (d suffix)
  ! ==================================================================
  Real(kind=8)                          :: locNPPd, locGPPd, locNNAd, locChldegd
  Real(kind=8),allocatable,dimension(:) :: vertNPPd, vertGPPd, vertNNAd, vertChldegd
  Real(kind=8),allocatable,dimension(:) :: vertrespd
  Real(kind=8),allocatable,dimension(:) :: vertdocexd
  Real(kind=8),allocatable,dimension(:) :: vertaggd

  ! ==================================================================
  ! COCCOLITHOPHORES (c suffix)
  ! ==================================================================
  Real(kind=8)                          :: locNPPc, locGPPc, locNNAc, locChldegc
  Real(kind=8),allocatable,dimension(:) :: vertNPPc, vertGPPc, vertNNAc, vertChldegc
  Real(kind=8),allocatable,dimension(:) :: vertrespc
  Real(kind=8),allocatable,dimension(:) :: vertdocexc
  Real(kind=8),allocatable,dimension(:) :: vertaggc
  Real(kind=8),allocatable,dimension(:) :: vertcalcdiss, vertcalcif

  ! ==================================================================
  ! PHAEOCYSTIS (p suffix)
  ! ==================================================================
  Real(kind=8)                          :: locNPPp, locGPPp, locNNAp, locChldegp
  Real(kind=8),allocatable,dimension(:) :: vertNPPp, vertGPPp, vertNNAp, vertChldegp
  Real(kind=8),allocatable,dimension(:) :: vertrespp
  Real(kind=8),allocatable,dimension(:) :: vertdocexp
  Real(kind=8),allocatable,dimension(:) :: vertaggp

  ! ==================================================================
  ! MICROZOOPLANKTON
  ! ==================================================================
  Real(kind=8)                          :: locgrazmicro_tot, locgrazmicro_n, locgrazmicro_d, locgrazmicro_c, locgrazmicro_p
  Real(kind=8),allocatable,dimension(:) :: vertgrazmicro_tot, vertgrazmicro_n, vertgrazmicro_d, vertgrazmicro_c, vertgrazmicro_p
  Real(kind=8),allocatable,dimension(:) :: vertrespmicro

  ! ==================================================================
  ! MESOZOOPLANKTON
  ! ==================================================================
  Real(kind=8)                          :: locgrazmeso_tot, locgrazmeso_n, locgrazmeso_d, locgrazmeso_c, locgrazmeso_p
  Real(kind=8)                          :: locgrazmeso_det, locgrazmeso_mic, locgrazmeso_det2
  Real(kind=8),allocatable,dimension(:) :: vertgrazmeso_tot, vertgrazmeso_n, vertgrazmeso_d, vertgrazmeso_c, vertgrazmeso_p
  Real(kind=8),allocatable,dimension(:) :: vertgrazmeso_det, vertgrazmeso_mic, vertgrazmeso_det2
  Real(kind=8),allocatable,dimension(:) :: vertrespmeso

  ! ==================================================================
  ! MACROZOOPLANKTON
  ! ==================================================================
  Real(kind=8)                          :: locgrazmacro_tot, locgrazmacro_n, locgrazmacro_d, locgrazmacro_c, locgrazmacro_p
  Real(kind=8)                          :: locgrazmacro_mes, locgrazmacro_det, locgrazmacro_mic, locgrazmacro_det2
  Real(kind=8),allocatable,dimension(:) :: vertgrazmacro_tot, vertgrazmacro_n, vertgrazmacro_d, vertgrazmacro_c, vertgrazmacro_p
  Real(kind=8),allocatable,dimension(:) :: vertgrazmacro_mes, vertgrazmacro_det, vertgrazmacro_mic, vertgrazmacro_det2
  Real(kind=8),allocatable,dimension(:) :: vertrespmacro

  Real(kind=8),allocatable,dimension(:) :: vertfastcdis
  Real(kind=8),allocatable,dimension(:) :: vertphotn, vertphotd, vertphotc, vertphotp          !Phytoplankton photosynthesis terms
  Real(kind=8),allocatable,dimension(:) :: vertmesocdis, vertmicrocdis, vertmacrocdis          !Additional zooplankton calcium dissolution terms 
  Real(kind=8),allocatable,dimension(:) :: vertNassimn, vertNassimd, vertNassimc, vertNassimp  !N assimilation by phytoplanktons
  Real(kind=8),allocatable,dimension(:) :: vertDONremin      !DON remineralization term
  Real(kind=8),allocatable,dimension(:) :: vertDOCremin      !DOC remineralization term 
  Real(kind=8),allocatable,dimension(:,:,:)   :: dtr_bf_din  ! Diagnostics for DIN bottom flux RP on 26.09.2025
  Real(kind=8),allocatable,dimension(:,:,:)   :: dtr_bf_dic  ! Diagnostics for DIC bottom flux
  Real(kind=8),allocatable,dimension(:,:,:)   :: dtr_bf_alk  ! Diagnostics for Alk bottom flux
  Real(kind=8),allocatable,dimension(:,:,:)   :: dtr_bf_dsi  ! Diagnostics for DSi bottom flux

  !===================================================================
  ! DISSOLUTION AND REMINERALIZATION ! R2OMIP
  !===================================================================
  Real(kind=8)                          :: locDISSOC, locDISSON, locDISSOSi, locREMOC, locREMOCt, locREMON
  Real(kind=8),allocatable,dimension(:) :: vertDISSOC, vertDISSON, vertDISSOSi, vertREMOC, vertREMOCt, vertREMON
!!------------------------------------------------------------------------------                                                                                
!! *** Benthos  ***
  Real(kind=8),allocatable,dimension(:) :: decayBenthos ! [1/day] Decay rate of detritus in the benthic layer
  Real(kind=8),allocatable,dimension(:) :: wFluxDet     ! [mmol/(m2 * day)] Flux of N,C,Si and calc through sinking of detritus
  Real(kind=8),allocatable,dimension(:) :: wFluxPhy     ! [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of phytoplankton
  Real(kind=8),allocatable,dimension(:) :: wFluxDia     ! [mmol/(m2 * day)] Flux of N,C, Si and chl through sinking of diatoms
  Real(kind=8),allocatable,dimension(:) :: wFluxCocco   ! NEW [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of coccos
  Real(kind=8),allocatable,dimension(:) :: wFluxPhaeo   ! NEW [mmol/(m2 * day)] Flux of N,C, calc and chl through sinking of Phaeocystis
  Real(kind=8)              :: Vben_det     ! [m/day] speed of sinking into benthos from water column
  Real(kind=8)              :: Vben_det_seczoo !second zooplankton sinking benthos  
  Real(kind=8)              :: Vben_phy
  Real(kind=8)              :: Vben_dia
  Real(kind=8)              :: Vben_coc     
  Real(kind=8)              :: Vben_pha     ! Phaeocystis
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

real(kind=8)                               :: is_3zoo2det
real(kind=8)                               :: is_coccos
real(kind=8)                               :: is_R2OMIP

end module REcoM_declarations


