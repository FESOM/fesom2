!===============================================================================
! For arrays needed for the whole 2D or 3D domain, but only needed in REcoM
!-------------------------------------------------------------------------------
Module REcoM_GloVar
  implicit none
  save

  Real(kind=8),allocatable,dimension(:,:) :: Benthos          ! 4 types of benthos-tracers with size [4 n2d]
  Real(kind=8),allocatable,dimension(:,:,:) :: Benthos_tr     ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel

  Real(kind=8),allocatable,dimension(:)   :: GloFeDust        ! [umol/m2/s] Monthly 2D field of iron soluted in surface water from dust
  Real(kind=8),allocatable,dimension(:)   :: GloNDust         ! [mmol/m2/s] 10-year mean 2D fields of nitrogen soluted in surface water from dust
  Real(kind=8),dimension(12)              :: AtmCO2           ! [uatm] Atmospheric CO2 partial pressure. One value for the whole planet for each month
  Real(kind=8),allocatable,dimension(:)   :: OmegaC_bottom    !< calcite saturation state at the ocean bottom !R2OMIP

  Real(kind=8),allocatable,dimension(:)   :: AtmFeInput       ! [umol/m2/s] Includes ice, but is, other than that identlical to GloFeDust
  Real(kind=8),allocatable,dimension(:)   :: AtmNInput        ! [umol/m2/s] Includes ice, but is, other than that identlical to GloNDust
  Real(kind=8),allocatable,dimension(:)   :: GloPCO2surf      ! [uatm] Surface ocean CO2 partial pressure
  Real(kind=8),allocatable,dimension(:)   :: GloCO2flux       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloO2flux        ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloCO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloO2flux_seaicemask       ! [mmol/m2/day] Positive downwards
  Real(kind=8),allocatable,dimension(:)   :: GloHplus         ! [mol/kg] Concentrations of H-plus ions in the surface ocean

  Real(kind=8),allocatable,dimension(:,:)   :: CO23D            ! MOCSY: [mol/m3] Aqueous CO2 concentration for all depths
  Real(kind=8),allocatable,dimension(:,:)   :: pH3D             ! MOCSY: total scale
  Real(kind=8),allocatable,dimension(:,:)   :: pCO23D           ! MOCSY: [uatm] CO2 partial pressure
  Real(kind=8),allocatable,dimension(:,:)   :: HCO33D           ! MOCSY: [mol/m3] Bicarbonate ion concentration
  Real(kind=8),allocatable,dimension(:,:)   :: CO33D            ! DISS: [mol/m3] Carbonate ion concentration
  Real(kind=8),allocatable,dimension(:,:)   :: OmegaC3D         ! DISS: calcite saturation state
  Real(kind=8),allocatable,dimension(:,:)   :: kspc3D           ! DISS: [mol^2/kg^2] stoichiometric solubility product of calcite
  Real(kind=8),allocatable,dimension(:,:)   :: rhoSW3D          ! DISS: [mol/m3] in-situ density of seawater

  Real(kind=8),allocatable,dimension(:,:)   :: rho_particle1       ! BALL: density of particle class 1
  Real(kind=8),allocatable,dimension(:,:)   :: rho_particle2       ! BALL: density of particle class 2
  Real(kind=8),allocatable,dimension(:,:)   :: scaling_density1_3D ! BALL: scaling factor
  Real(kind=8),allocatable,dimension(:,:)   :: scaling_density2_3D ! BALL: scaling factor
  Real(kind=8),allocatable,dimension(:,:)   :: scaling_visc_3D     ! BALL: scaling factor
  Real(kind=8),allocatable,dimension(:,:)   :: seawater_visc_3D    ! BALL: scaling factor

  Real(kind=8),allocatable,dimension(:)     :: GlodPCO2surf       ! [mmol/m2/day] ocean-atmosphere  
  Real(kind=8),allocatable,dimension(:,:)   :: GlodecayBenthos  ! [1/day] Decay rate of detritus in the benthic layer saved for oce_ale_tracer.F90
  Real(kind=8),allocatable,dimension(:,:)   :: Sed_2_Ocean_Flux ! Flux from the sediment back to the bottom ocean ! R2OMIP
  Real(kind=8),allocatable,dimension(:,:)   :: Ocean_2_Sed_Flux ! Flux from the bottom ocean to the sediment ! R2OMIP
  Real(kind=8),allocatable,dimension(:,:)   :: Burial         ! Benthic permanent burial Field in 2D [n2d 1] ! R2OMIP
  Real(kind=8),allocatable,dimension(:)     :: PistonVelocity   ! [m s-1]
  Real(kind=8),allocatable,dimension(:)     :: alphaCO2         ! [mol L-1 atm-1]

  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxDet    ! 
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxPhy    ! 
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxDia    ! 
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxCocco  !
  Real(kind=8),allocatable,dimension(:,:)   :: GlowFluxPhaeo  !

  Real(kind=8),allocatable,dimension(:,:)   :: diags2D          ! Diagnostics in 2D [8 n2d]
  Real(kind=8),allocatable,dimension(:)     :: NPPn
  Real(kind=8),allocatable,dimension(:)     :: NPPd
  Real(kind=8),allocatable,dimension(:)     :: GPPn
  Real(kind=8),allocatable,dimension(:)     :: GPPd
  Real(kind=8),allocatable,dimension(:)     :: NNAn
  Real(kind=8),allocatable,dimension(:)     :: NNAd
  Real(kind=8),allocatable,dimension(:)     :: Chldegn
  Real(kind=8),allocatable,dimension(:)     :: Chldegd
  Real(kind=8),allocatable,dimension(:)     :: NPPc
  Real(kind=8),allocatable,dimension(:)     :: GPPc
  Real(kind=8),allocatable,dimension(:)     :: NNAc
  Real(kind=8),allocatable,dimension(:)     :: Chldegc
  Real(kind=8),allocatable,dimension(:)     :: NPPp            ! Phaeocystis
  Real(kind=8),allocatable,dimension(:)     :: GPPp
  Real(kind=8),allocatable,dimension(:)     :: NNAp
  Real(kind=8),allocatable,dimension(:)     :: Chldegp
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_tot
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_n
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_d
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_c
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_p
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_det
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_mic
  Real(kind=8),allocatable,dimension(:)     :: grazmeso_det2
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_tot
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_n
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_d
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_c
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_p
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_mes
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_det
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_mic
  Real(kind=8),allocatable,dimension(:)     :: grazmacro_det2
  Real(kind=8),allocatable,dimension(:)     :: grazmicro_tot
  Real(kind=8),allocatable,dimension(:)     :: grazmicro_n
  Real(kind=8),allocatable,dimension(:)     :: grazmicro_d
  Real(kind=8),allocatable,dimension(:)     :: grazmicro_c
  Real(kind=8),allocatable,dimension(:)     :: grazmicro_p
  Real(kind=8),allocatable,dimension(:)     :: DISSOC
  Real(kind=8),allocatable,dimension(:)     :: DISSON
  Real(kind=8),allocatable,dimension(:)     :: DISSOSi
  Real(kind=8),allocatable,dimension(:)     :: REMOC
  Real(kind=8),allocatable,dimension(:)     :: REMOCt
  Real(kind=8),allocatable,dimension(:)     :: REMON
  Real(kind=8),allocatable,dimension(:,:)   :: respmeso
  Real(kind=8),allocatable,dimension(:,:)   :: respmacro
  Real(kind=8),allocatable,dimension(:,:)   :: respmicro
  Real(kind=8),allocatable,dimension(:,:)   :: calcdiss
  Real(kind=8),allocatable,dimension(:,:)   :: calcif
  Real(kind=8),allocatable,dimension(:,:)   :: aggn
  Real(kind=8),allocatable,dimension(:,:)   :: aggd
  Real(kind=8),allocatable,dimension(:,:)   :: aggc
  Real(kind=8),allocatable,dimension(:,:)   :: aggp             ! Phaeocystis
  Real(kind=8),allocatable,dimension(:,:)   :: docexn
  Real(kind=8),allocatable,dimension(:,:)   :: docexd
  Real(kind=8),allocatable,dimension(:,:)   :: docexc
  Real(kind=8),allocatable,dimension(:,:)   :: docexp           ! Phaeocystis
  Real(kind=8),allocatable,dimension(:,:)   :: respn
  Real(kind=8),allocatable,dimension(:,:)   :: respd
  Real(kind=8),allocatable,dimension(:,:)   :: respc
  Real(kind=8),allocatable,dimension(:,:)   :: respp            ! Phaeocystis
  Real(kind=8),allocatable,dimension(:,:)   :: NPPn3D
  Real(kind=8),allocatable,dimension(:,:)   :: NPPd3D
  Real(kind=8),allocatable,dimension(:,:)   :: NPPc3D
  Real(kind=8),allocatable,dimension(:,:)   :: NPPp3D           ! Phaeocystis

  Real(kind=8),allocatable,dimension(:,:)   :: photn
  Real(kind=8),allocatable,dimension(:,:)   :: photd
  Real(kind=8),allocatable,dimension(:,:)   :: photc
  Real(kind=8),allocatable,dimension(:,:)   :: photp
  Real(kind=8),allocatable,dimension(:,:)   :: DOCremin
  Real(kind=8),allocatable,dimension(:,:)   :: Nassimn
  Real(kind=8),allocatable,dimension(:,:)   :: Nassimd
  Real(kind=8),allocatable,dimension(:,:)   :: Nassimc
  Real(kind=8),allocatable,dimension(:,:)   :: Nassimp 
  Real(kind=8),allocatable,dimension(:,:)   :: DONremin
  Real(kind=8),allocatable,dimension(:,:)   :: fastcdis
  Real(kind=8),allocatable,dimension(:,:)   :: mesocdis
  Real(kind=8),allocatable,dimension(:,:)   :: microcdis
  Real(kind=8),allocatable,dimension(:,:)   :: macrocdis

  Real(kind=8),allocatable,dimension(:,:)   :: TTemp_diatoms ! my new variables to track
  Real(kind=8),allocatable,dimension(:,:)   :: TTemp_phyto ! new Temperature effect
  Real(kind=8),allocatable,dimension(:,:)   :: TTemp_cocco ! new
  Real(kind=8),allocatable,dimension(:,:)   :: TTemp_phaeo ! new
  Real(kind=8),allocatable,dimension(:,:)   :: TPhyCO2 ! new CO2 effect
  Real(kind=8),allocatable,dimension(:,:)   :: TDiaCO2
  Real(kind=8),allocatable,dimension(:,:)   :: TCoccoCO2
  Real(kind=8),allocatable,dimension(:,:)   :: TPhaeoCO2
  Real(kind=8),allocatable,dimension(:,:)   :: TqlimitFac_phyto ! new nutrient limitation
  Real(kind=8),allocatable,dimension(:,:)   :: TqlimitFac_diatoms
  Real(kind=8),allocatable,dimension(:,:)   :: TqlimitFac_cocco
  Real(kind=8),allocatable,dimension(:,:)   :: TqlimitFac_phaeo
  Real(kind=8),allocatable,dimension(:,:)   :: TCphotLigLim_phyto ! new light limitation
  Real(kind=8),allocatable,dimension(:,:)   :: TCphot_phyto       ! new
  Real(kind=8),allocatable,dimension(:,:)   :: TCphotLigLim_diatoms ! new light limitation
  Real(kind=8),allocatable,dimension(:,:)   :: TCphot_diatoms
  Real(kind=8),allocatable,dimension(:,:)   :: TCphotLigLim_cocco ! new light limitation
  Real(kind=8),allocatable,dimension(:,:)   :: TCphot_cocco
  Real(kind=8),allocatable,dimension(:,:)   :: TCphotLigLim_phaeo ! new light limitation
  Real(kind=8),allocatable,dimension(:,:)   :: TCphot_phaeo
  Real(kind=8),allocatable,dimension(:,:)   :: TSi_assimDia ! tracking the assimilation of Si by Diatoms

  Real(kind=8),allocatable,dimension(:)     :: DenitBen         ! Benthic denitrification Field in 2D [n2d 1]

!  for using MEDUSA
  Real(kind=8),allocatable,dimension(:,:)   :: SinkFlx         ! Diagnostics in 2D [4 n2d] or [6 n2d] with ciso
  Real(kind=8),allocatable,dimension(:,:,:) :: SinkFlx_tr      ! kh 25.03.22 buffer sums per tracer index to avoid non bit identical results regarding global sums when running the tracer loop in parallel
  Real(kind=8),allocatable,dimension(:,:)   :: Sinkingvel1     ! Diagnostics for vertical sinking
  Real(kind=8),allocatable,dimension(:,:)   :: Sinkingvel2     ! Diagnostics for vertical sinking  
  Real(kind=8),allocatable,dimension(:,:,:) :: Sinkvel1_tr     ! Sinking speed of particle class 1 OG 16.03.23 
  Real(kind=8),allocatable,dimension(:,:,:) :: Sinkvel2_tr     ! Sinking speed of particle class 2 OG 16.03.23
  Real(kind=8),allocatable,dimension(:,:)   :: dtr_bflux_din   ! Diagnostics for DIN bottom flux RP on 30.09.2025 
  Real(kind=8),allocatable,dimension(:,:)   :: dtr_bflux_dic   ! Diagnostics for DIC bottom flux
  Real(kind=8),allocatable,dimension(:,:)   :: dtr_bflux_alk   ! Diagnostics for Alk bottom flux
  Real(kind=8),allocatable,dimension(:,:)   :: dtr_bflux_dsi   ! Diagnostics for DSi bottom flux

  Real(kind=8),allocatable,dimension(:,:)   :: GloSed           ! Yearly input into bottom water from sediments [n2d 5] or [n2d 7] with ciso
  Real(kind=8),allocatable,dimension(:,:)   :: lb_flux          ! Yearly burial from medusa: [n2d 5] or [n2d 9] with ciso_14 

! atmospheric box model:
  Real(kind=8),allocatable,dimension(:)     :: x_co2atm         ! atmospheric CO2 mixing ratio (mole fraction)

  Real(kind=8), allocatable,dimension(:)    :: Alk_surf         ! Surface alkalinity field used for restoring
  Real(kind=8), allocatable,dimension(:)    :: relax_alk
  Real(kind=8), allocatable,dimension(:)    :: virtual_alk

  real(kind=8), allocatable,dimension(:,:)  :: PAR3D            ! Light in the water column [nl-1 n2d]
  real(kind=8), allocatable,dimension(:)    :: RiverineLonOrig, RiverineLatOrig, RiverineDINOrig, RiverineDONOrig, RiverineDOCOrig, RiverineDSiOrig ! Variables to save original values for riverine nutrients
  real(kind=8), allocatable,dimension(:)    :: RiverDIC2D, RiverDIN2D, RiverDOCl2D, RiverDOCsl2D, RiverPOC2D, RiverFe
  real(kind=8), allocatable,dimension(:)    :: RiverDON2D, RiverDOC2D, RiverDSi2D, RiverAlk2D
  Real(kind=8),allocatable,dimension(:)     :: LocDenit
  Real(kind=8),allocatable,dimension(:,:)   :: LocBurial ! R2OMIP
  Real(kind=8),allocatable,dimension(:)     :: BurialBen ! R2OMIP

  real(kind=8), allocatable,dimension(:)    :: ErosionTSi2D, ErosionTON2D, ErosionTOC2D
!! Cobeta, Cos(Angle of incidence)
  Real(kind=8), allocatable,dimension(:)    ::  cosAI
end module REcoM_GloVar

