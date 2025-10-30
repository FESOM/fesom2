! Adapted from FESOM code of Q. Wang
!
! 1) modules for forcing parameters and arrays
! 2) modules for updating forcing in time and interpolation in space.
!------------------------------------------------------------------------

module g_forcing_param
use o_param
implicit none
save   

  ! *** exchange coefficients ***
real(kind=WP)  :: Ce_atm_oce=1.75e-3 ! exch.coeff. of lat. heat over open water
real(kind=WP)  :: Ch_atm_oce=1.75e-3 ! exch.coeff. of sens. heat over open water
real(kind=WP)  :: Cd_atm_oce=1.0e-3  ! drag coeff. between atmosphere and water

real(kind=WP)  :: Ce_atm_ice=1.75e-3 ! exch. coeff. of latent heat over ice
real(kind=WP)  :: Ch_atm_ice=1.75e-3 ! exch. coeff. of sensible heat over ice
real(kind=WP)  :: Cd_atm_ice=1.32e-3 ! drag coeff. between atmosphere and ice 
real(kind=WP)  :: Swind     =0.0_WP  ! parameterization for coupled current feedback after Renault et al. 2019; varies from 0 (not parameterized) to 1 (no ocean contribution to wind stress)

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice, Swind


  ! *** forcing source and type ***
  logical                       :: use_virt_salt ! will be set TRUE in case of which_ALE='linfs', otherwise FALSE

  ! *** coefficients in bulk formulae ***
  logical                       :: AOMIP_drag_coeff=.false.
  logical                       :: ncar_bulk_formulae=.false.
  real(kind=WP)                 :: ncar_bulk_z_wind=10.0_WP
  real(kind=WP)                 :: ncar_bulk_z_tair=10.0_WP
  real(kind=WP)                 :: ncar_bulk_z_shum=10.0_WP

  namelist /forcing_bulk/ AOMIP_drag_coeff, ncar_bulk_formulae, ncar_bulk_z_wind, ncar_bulk_z_tair, ncar_bulk_z_shum

  ! *** add land ice melt water ***
  logical                       :: use_landice_water=.false.
  integer                       :: landice_start_mon=1
  integer                       :: landice_end_mon=12
  !---fwf-code-begin
  character(MAX_PATH)           :: fwf_path='./mesh/'
  !---fwf-code-end

  namelist /land_ice/ use_landice_water, landice_start_mon, landice_end_mon, fwf_path !---fwf-code, add fwf_path

  !---age-code-begin
  logical                       :: use_age_tracer=.false.
  logical                       :: use_age_mask=.false.
  character(MAX_PATH)           :: age_tracer_path='./mesh/'
  integer                       :: age_start_year=2000

  namelist /age_tracer/ use_age_tracer, use_age_mask, age_tracer_path, age_start_year
  !---age-code-end

end module g_forcing_param
! ====================================================================
module g_forcing_arrays
use o_param
  implicit none
  save    

  ! forcing arrays
  real(kind=WP), allocatable, dimension(:)         :: u_wind, v_wind 
  real(kind=WP), allocatable, dimension(:)         :: u_wind_ib, v_wind_ib ! kh 19.02.21 additional arrays for asynchronous iceberg computations

  real(kind=WP), allocatable, dimension(:)         :: Tair, shum
  real(kind=WP), allocatable, dimension(:,:)       :: u_wind_t, v_wind_t 
  real(kind=WP), allocatable, dimension(:,:)       :: Tair_t, shum_t
  real(kind=WP), allocatable, dimension(:)         :: shortwave, longwave
  real(kind=WP), allocatable, dimension(:)         :: prec_rain, prec_snow
  real(kind=WP), allocatable, dimension(:)         :: runoff, evaporation, ice_sublimation
  real(kind=WP), allocatable, dimension(:)         :: cloudiness, press_air
  !---wiso-code
  real(kind=WP), allocatable, dimension(:)         :: www1,www2,www3,iii1,iii2,iii3
  real(kind=WP), allocatable, dimension(:)         :: tmp_iii1,tmp_iii2,tmp_iii3
  !---wiso-code-end
  !---age-code-begin
  integer, allocatable, dimension(:)               :: age_tracer_loc_index
  !---age-code-end

#if defined (__oasis) || defined (__ifsinterface) || defined (__yac) /* todo:use a single shared definition  */
  real(kind=WP), target, allocatable, dimension(:) :: sublimation, evap_no_ifrac
#endif
#if defined (__oasis) || defined (__yac)
  real(kind=WP), target, allocatable, dimension(:) :: tmp_sublimation, tmp_evap_no_ifrac !temporary flux fields
  real(kind=WP), target, allocatable, dimension(:) :: tmp_shortwave 			!(for flux correction) 
  real(kind=WP), allocatable, dimension(:)         :: atm_net_fluxes_north, atm_net_fluxes_south
  real(kind=WP), allocatable, dimension(:)         :: oce_net_fluxes_north, oce_net_fluxes_south
  real(kind=WP), allocatable, dimension(:)         :: flux_correction_north, flux_correction_south, flux_correction_total
#endif

#if defined (__oasis) || defined (__ifsinterface)
  real(kind=WP), allocatable, dimension(:)         :: residualifwflx
#endif

  real(kind=WP), allocatable, dimension(:)         :: runoff_landice
  real(kind=WP)                                    :: landice_season(12)

  ! shortwave penetration
  real(kind=WP), allocatable, dimension(:)         :: chl
  real(kind=WP), allocatable, dimension(:,:)       :: sw_3d

!   real(kind=WP), allocatable, dimension(:)         :: thdgr, thdgrsn
  real(kind=WP), allocatable, dimension(:)         :: flice
  real(kind=WP), allocatable, dimension(:)         :: hf_Qlat, hf_Qsen, hf_Qradtot, hf_Qswr, hf_Qlwr, hf_Qlwrout
  real(kind=WP), allocatable, dimension(:)         :: real_salt_flux !PS

  ! drag coefficient Cd_atm_oce and transfer coefficients for evaporation
  ! Ce_atm_oce and sensible heat Ch_atm_oce between atmosphere and ocean
  real(kind=WP), allocatable, dimension(:)	  :: Cd_atm_oce_arr
  real(kind=WP), allocatable, dimension(:)	  :: Ch_atm_oce_arr
  real(kind=WP), allocatable, dimension(:)	  :: Ce_atm_oce_arr

  ! drag coefficient Cd_atm_oce between atmosphere and ice
  real(kind=WP), allocatable, dimension(:)	  :: Cd_atm_ice_arr
  
  real(kind=WP),allocatable,dimension(:)          ::  aver_temp
end module g_forcing_arrays
