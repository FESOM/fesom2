! =====================
!  Sea ice
!  Finite-volume implementation
!  Modules for coupled version 
!  Only EVP solver is available in this distrib. memory setup
! ======================
!  Ice velocity is defined at nodes
!===========================================================================
!
MODULE i_PARAM
  !
  ! Ice specific parameters
  !
  USE o_PARAM
  IMPLICIT NONE
  SAVE
  ! ice model parameters:
  ! RHEOLOGY
  REAL(kind=WP)             :: Pstar = 30000.0_WP        ![N/m^2]
  REAL(kind=WP)             :: ellipse =2.0_WP           !
  REAL(kind=WP)             :: c_pressure =20.0_WP       !
  REAL(kind=WP)             :: delta_min=1.0e-11         ! [s^(-1)]
  REAL(kind=WP)             :: Clim_evp=615              ! kg/m^2
  REAL(kind=WP)             :: zeta_min=4.0e+8           ! kg/s
  INTEGER                   :: evp_rheol_steps=120       ! EVP rheology
                                                         ! cybcycling steps
  REAL(kind=WP)             :: ice_gamma_fct=0.25_WP     ! smoothing parameter
                                                         ! in ice fct advection
  REAL(kind=WP)             :: ice_diff=10.0_WP          ! diffusion to stabilize
                                                         ! ice advection
  REAL(kind=WP)             :: Tevp_inv                  
  real(kind=WP)             :: theta_io=0.0_WP           ! rotation angle
                                                         ! (ice-ocean), available
						         ! in EVP
  real(kind=WP)             :: alpha_evp=250, beta_evp=250
  real(kind=WP)             :: c_aevp=0.15 ! 0.1--0.2, but should be adjusted experimentally   
  ! Ice forcing averaging
  integer		    :: ice_ave_steps=1 !ice step=ice_ave_steps*oce_step
  real(kind=WP)             :: cd_oce_ice = 5.5e-3       ! drag coef. oce - ice      

  logical                   :: ice_free_slip=.false.
  integer                   :: whichEVP=0 !0=standart; 1=mEVP; 2=aEVP
  real(kind=WP)             :: ice_dt !ice step=ice_ave_steps*oce_step

NAMELIST /ice_dyn/ whichEVP, Pstar, ellipse, c_pressure, delta_min, evp_rheol_steps, Cd_oce_ice, &
ice_gamma_fct, ice_diff, theta_io, ice_ave_steps, alpha_evp, beta_evp, c_aevp

END MODULE i_PARAM
!
!=============================================================================
!
MODULE i_ARRAYS
!
! Arrays used to store ice variables and organize coupling
!
USE o_PARAM
implicit none
save
  logical                   :: ice_update = .true. !
  integer                   :: ice_steps_since_upd = 0 !
  real(kind=WP),allocatable,dimension(:,:)         :: ice_grad_vel
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: U_ice, V_ice, m_ice, a_ice  
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: U_ice_old, V_ice_old, m_ice_old, a_ice_old, m_snow_old,thdgr_old !PS
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: U_rhs_ice, V_rhs_ice, m_snow
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: rhs_m, rhs_a, rhs_ms, ths_temp
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: U_w, V_w
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: u_ice_aux, v_ice_aux  ! of the size of u_ice, v_ice
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: rhs_mdiv, rhs_adiv, rhs_msdiv
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: elevation
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: sigma11, sigma12, sigma22
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: eps11, eps12, eps22
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: fresh_wa_flux
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: net_heat_flux
#if defined (__oasis)
  real(kind=WP),target, allocatable, dimension(:)  :: ice_alb, ice_temp ! new fields for OIFS coupling
  real(kind=WP),target, allocatable, dimension(:)  :: oce_heat_flux, ice_heat_flux
  real(kind=WP),target, allocatable, dimension(:)  :: tmp_oce_heat_flux, tmp_ice_heat_flux
							!temporary flux fields
							!(for flux correction)
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: rhs_temp, m_templ, dm_temp, rhs_tempdiv
#if defined (__oifs)
  real(kind=WP),target, allocatable, dimension(:)  :: enthalpyoffuse
#endif
#endif /* (__oasis) */

  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: S_oc_array, T_oc_array
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_x         
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_y
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_x         
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_y
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: t_skin
 ! FCT implementation
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: m_icel, a_icel, m_snowl
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: dm_ice, da_ice, dm_snow
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:,:)        :: icefluxes
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: icepplus, icepminus
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: mass_matrix  
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: alpha_evp_array(:)   ! of myDim_elem2D
 REAL(kind=WP), ALLOCATABLE, DIMENSION(:)          :: beta_evp_array(:)    ! of myDim_node2D+eDim_node2D

! Mean arrays
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: U_ice_mean, V_ice_mean
  REAL(kind=WP), ALLOCATABLE, DIMENSION(:)         :: m_ice_mean, a_ice_mean, m_snow_mean
  END MODULE i_ARRAYS
!=====================================================================
module i_therm_param
USE o_PARAM
  implicit none
REAL(kind=WP), parameter  :: rhoair=  1.3            ! Air density,  LY2004 !1.3 AOMIP
REAL(kind=WP), parameter  :: inv_rhoair=  1./1.3     ! Air density,  LY2004 !1.3 AOMIP
REAL(kind=WP), parameter  :: rhowat= 1025.            ! Water density
REAL(kind=WP), parameter  :: inv_rhowat= 1./1025.     ! Inverse Water density
REAL(kind=WP), parameter  :: rhoice=  910.            ! Ice density, AOMIP
REAL(kind=WP), parameter  :: inv_rhoice=  1./910.     ! Ice density, AOMIP
REAL(kind=WP), parameter  :: rhosno=  290.            ! Snow density, AOMIP
REAL(kind=WP), parameter  :: inv_rhosno=  1./290.     ! Snow density, AOMIP

REAL(kind=WP), parameter  :: cpair=1005.       ! Specific heat of air [J/(kg * K)] 
REAL(kind=WP), parameter  :: cpice=2106.       ! Specific heat of ice [J/(kg * K)] 
REAL(kind=WP), parameter  :: cpsno=2090.       ! Specific heat of snow [J/(kg * K)] 
REAL(kind=WP), parameter  :: cc=rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
REAL(kind=WP), parameter  :: cl=rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf) 
REAL(kind=WP), parameter  :: clhw=2.501e6      ! Specific latent heat [J/kg]: water	-> water vapor
REAL(kind=WP), parameter  :: clhi=2.835e6      !                              sea ice-> water vapor
 
REAL(kind=WP), parameter  :: tmelt=273.15      ! 0 deg C expressed in K 
REAL(kind=WP), parameter  :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity

REAL(kind=WP)    :: con   = 2.1656    ! Thermal conductivities: ice; W/m/K
REAL(kind=WP)    :: consn = 0.31      !                         snow

REAL(kind=WP)    :: Sice = 4.0        ! Ice salinity 3.2--5.0 ppt.

integer          :: iclasses=7        ! Number of ice thickness gradations for ice growth calcs.
REAL(kind=WP)    :: h0=1.0	      ! Lead closing parameter [m] ! 0.5

REAL(kind=WP)    :: hmin= 0.01        ! Cut-off ice thickness     !!
REAL(kind=WP)    :: Armin=0.01        ! Minimum ice concentration !!

REAL(kind=WP)    :: emiss_ice=0.97        ! Emissivity of Snow/Ice, 
REAL(kind=WP)    :: emiss_wat=0.97        ! Emissivity of open water

REAL(kind=WP)    :: albsn=   0.81     ! Albedo: frozen snow
REAL(kind=WP)    :: albsnm=  0.77     !         melting snow
REAL(kind=WP)    :: albi=    0.70     !         frozen ice
REAL(kind=WP)    :: albim=   0.68     !         melting ice
REAL(kind=WP)    :: albw=    0.066    !         open water, LY2004

  NAMELIST /ice_therm/ Sice, h0, emiss_ice, &
  emiss_wat, albsn, albsnm, albi, albim, albw, con, consn

end module i_therm_param


!==============================================================================
