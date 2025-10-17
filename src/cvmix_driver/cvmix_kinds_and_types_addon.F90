module cvmix_kinds_and_types_addon

!BOP
! !MODULE:  cvmix_kinds_and_types_addon
!
! !AUTHOR:
!  Michael Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains the declarations for all required vertical mixing
!  data types. It also contains several global parameters used by the cvmix
!  package, such as kind numbers and string lengths.
!  \\
!  \\

! !USES:
!  uses no other modules
!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! The cvmix package uses double precision for floating point computations.
  integer, parameter, public :: cvmix_r8       = selected_real_kind(15, 307), &
                                cvmix_log_kind = kind(.true.),                &
                                cvmix_strlen   = 256

  ! Parameters to allow CVMix to store integers instead of strings
  integer, parameter, public :: CVMIX_OVERWRITE_OLD_VAL    = 1
  integer, parameter, public :: CVMIX_SUM_OLD_AND_NEW_VALS = 2
  integer, parameter, public :: CVMIX_MAX_OLD_AND_NEW_VALS = 3

  ! Global parameters:
  ! The constant 1 is used repeatedly in PP and double-diff mixing.
  ! The value for pi is needed for Bryan-Lewis mixing.
  real(cvmix_r8), parameter, public :: cvmix_zero = real(0,cvmix_r8),         &
                                       cvmix_one  = real(1,cvmix_r8)
  real(cvmix_r8), parameter, public :: cvmix_PI   = &
                                       3.14159265358979323846_cvmix_r8

! !PUBLIC TYPES:

  ! cvmix_data_type_addon contains variables for time-dependent and column-specific
  ! mixing. Time-independent physical parameters should be stored in
  ! cvmix_global_params_type and *-mixing specific parameters should be
  ! stored in cvmix_*_params_type (found in the cvmix_* module).
  type, public :: cvmix_data_type_addon
    integer        :: nlev = -1      ! Number of active levels in column
    integer        :: max_nlev = -1  ! Number of levels in column
                                     ! Setting defaults to -1 might be F95...

    ! Scalar quantities
    ! -----------------
    ! distance from sea level to ocean bottom (positive => below sea level)
    real(cvmix_r8) :: OceanDepth
                    ! units: m
    ! distance from sea level to OBL bottom (positive => below sea level)
    real(cvmix_r8) :: BoundaryLayerDepth
                    ! units: m
    ! sea surface height (positive => above sea level)
    real(cvmix_r8) :: SeaSurfaceHeight
                    ! units: m
    ! turbulent friction velocity at surface
    real(cvmix_r8) :: SurfaceFriction
                    ! units: m/s
    ! buoyancy forcing at surface
    real(cvmix_r8) :: SurfaceBuoyancyForcing
                    ! units: m^2 s^-3
    ! latitude of column
    real(cvmix_r8) :: lat
                    ! units: degrees
    ! longitude of column
    real(cvmix_r8) :: lon
                    ! units: degrees
    ! Coriolis parameter
    real(cvmix_r8) :: Coriolis
                    ! units: s^-1
    ! Index of cell containing OBL (fraction > .5 => below cell center)
    real(cvmix_r8) :: kOBL_depth
                    ! units: unitless
    ! Langmuir mixing induced enhancement factor to turbulent velocity scale
    real(cvmix_r8) :: LangmuirEnhancementFactor
                    ! units: unitless
    ! Langmuir number
    real(cvmix_r8) :: LangmuirNumber
                    ! units: unitless
    ! Stokes Similarity Parameter
    real(cvmix_r8) :: StokesMostXi
                    ! units: unitless
    ! Numerical limit of Ocean Boundary Layer Depth
    real(cvmix_r8) :: zBottomOceanNumerics
                    ! units: m
    ! A time-invariant coefficient needed for Simmons, et al. tidal mixing
    real(cvmix_r8) :: SimmonsCoeff

    !IDEMIX
    real(cvmix_r8) :: forc_iw_bottom              ! energy flux into IW at bottom, [m^3/s^3]
    real(cvmix_r8) :: forc_iw_surface             ! energy flux into IW at surface, [m^3/s^3]           

    !Tidal (Simmons et al. (2004)
    !real(cvmix_r8) :: forc_tidal_bottom           ! energy flux into IW at bottom, [m^3/s^3]

    !TKE
    ! wind forcing for TKE, input should be sqrt(taux^2+tauy^2)^(3/2) for flux
    ! boundary condition; 
    !with taux, tauy the x and y-components of the surface wind stress; for
    !Dirichlet BC assignment of uppermost values in tridiag-matrix need to be
    !adjusted
    real(cvmix_r8) :: forc_tke_surf
                    ! units: m^3/s^3
    !bottom friction 
    real(cvmix_r8) :: bottom_fric
                    !units: m^2/s^3
    !surface buoyancy forcing for TKE 
    real(cvmix_r8) :: forc_rho_surf
                    !units: kg/m^2/s
    !needed
    real(cvmix_r8) :: dtime 
    real(cvmix_r8) :: rho_ref 


    ! Values on interfaces (dimsize = nlev+1)
    ! --------------------
    ! height of interfaces in column (positive up => most are negative)
    real(cvmix_r8), dimension(:), pointer :: zw_iface => NULL()
                                           ! units: m

    ! distance between neighboring cell centers (first value is top of ocean to
    ! middle of first cell, last value is middle of last cell to ocean bottom
    real(cvmix_r8), dimension(:), pointer :: dzw                   => NULL()
                                           ! units: m

    ! diffusivity coefficients at interfaces
    ! different coefficients for momentum (Mdiff), temperature (Tdiff), and
    ! salinity / non-temp tracers (Sdiff)
    real(cvmix_r8), dimension(:), pointer :: Mdiff_iface => NULL()
    real(cvmix_r8), dimension(:), pointer :: Tdiff_iface => NULL()
    real(cvmix_r8), dimension(:), pointer :: Sdiff_iface => NULL()
                                           ! units: m^2/s

    ! shear Richardson number at column interfaces
    real(cvmix_r8), dimension(:), pointer :: ShearRichardson_iface => NULL()
                                           ! units: unitless

    ! For tidal mixing, we need the squared buoyancy frequency and vertical
    ! deposition function
    real(cvmix_r8), dimension(:), pointer :: SqrBuoyancyFreq_iface => NULL()
                                           ! units: s^-2
    real(cvmix_r8), dimension(:), pointer :: VertDep_iface => NULL()
                                           ! units: unitless

    ! A time-dependent coefficient needed for Schmittner 2014
    real(cvmix_r8), dimension(:), pointer   :: SchmittnerCoeff => NULL()

    ! A time-invariant coefficient needed in Schmittner tidal mixing
    real(cvmix_r8), dimension(:), pointer   :: SchmittnerSouthernOcean => NULL()

    ! Another time-invariant coefficient needed in Schmittner tidal mixing
    real(cvmix_r8), dimension(:,:), pointer :: exp_hab_zetar => NULL()



    ! For KPP, need to store non-local transport term
    real(cvmix_r8), dimension(:), pointer :: kpp_Tnonlocal_iface => NULL()
    real(cvmix_r8), dimension(:), pointer :: kpp_Snonlocal_iface => NULL()
                                           ! units: unitless (see note below)
    ! Note that kpp_transport_iface is the value of K_x*gamma_x/flux_x: in
    ! other words, the user must multiply this value by either the freshwater
    ! flux or the penetrative shortwave heat flux to come the values in Eqs.
    ! (7.128) and (7.129) of the CVMix manual.
    ! Currently only provide nonlocal term for temperature tracer and salinity
    ! (non-temperature) tracers. Eventually may add support for momentum terms
    ! (would be 2D for x- and y-, respectively) but current implementation
    ! assumes momentum term is 0 everywhere.

    
    ! TKE / IDEMIX  
    ! diffusivity coefficients at interfaces
    ! different coefficients for momentum (KappaM), temperature and
    ! salinity  (KappaH) units: m^2/s
    real(cvmix_r8), dimension(:), pointer :: KappaM_iface => NULL()
    real(cvmix_r8), dimension(:), pointer :: KappaH_iface => NULL()
    !currently not in use; necessary for energy consistent linking of parameterizations 
    real(cvmix_r8), dimension(:), pointer :: Kappa_GM => NULL()
    ! total energy of internal waves (m^2/s^2) 
    real(cvmix_r8), dimension(:), pointer ::  E_iw    => NULL()
    ! total dissipation of IW energy (m^2/s^3)
    real(cvmix_r8), dimension(:), pointer :: iw_diss   => NULL()
    ! dissipation parameter (s/m^2)
    real(cvmix_r8), dimension(:), pointer :: alpha_c   => NULL()
    ! Turbulent kinetic energy (m^2/s^2)
    real(cvmix_r8), dimension(:), pointer :: tke   => NULL()
    ! dissipation of TKE (m^2/s^3)
    real(cvmix_r8), dimension(:), pointer :: tke_diss   => NULL()
    ! square of vertical shear of horizontal velocity at column interfaces (1/s^2)
    real(cvmix_r8), dimension(:), pointer :: Ssqr_iface => NULL()
    ! squared buoyancy frequency (s^-2) 
    real(cvmix_r8), dimension(:), pointer :: Nsqr_iface => NULL()



    ! Values at tracer points (dimsize = nlev)
    ! -----------------------
    ! height of cell centers in column (positive up => most are negative)
    real(cvmix_r8), dimension(:), pointer :: zt_cntr => NULL()
                                           ! units: m

    ! level thicknesses (positive semi-definite)
    real(cvmix_r8), dimension(:), pointer :: dzt => NULL()
                                           ! units: m

    ! Two density values are stored: the actual density of water at a given
    ! level and the the density of water after adiabatic displacement to the
    ! level below where the water actually is
    real(cvmix_r8), dimension(:), pointer :: WaterDensity_cntr      => NULL()
    real(cvmix_r8), dimension(:), pointer :: AdiabWaterDensity_cntr => NULL()
                                           ! units: kg m^-3

    ! bulk Richardson number
    real(cvmix_r8), dimension(:), pointer :: BulkRichardson_cntr => NULL()
                                           ! units: unitless

    ! For double diffusion mixing, we need to calculate the stratification
    ! parameter R_rho. Since the denominator of this ratio may be zero, we
    ! store the numerator and denominator separately and make sure the
    ! denominator is non-zero before performing the division.
    real(cvmix_r8), dimension(:), pointer :: strat_param_num   => NULL()
    real(cvmix_r8), dimension(:), pointer :: strat_param_denom => NULL()
                                           ! units: unitless

    ! For KPP we need velocity (in both x direction and y direction)
    real(cvmix_r8), dimension(:), pointer :: Vx_cntr => NULL()
    real(cvmix_r8), dimension(:), pointer :: Vy_cntr => NULL()
                                           ! units: m/s
  end type cvmix_data_type_addon


end module cvmix_kinds_and_types_addon

