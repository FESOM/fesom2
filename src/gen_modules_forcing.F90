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

  namelist /land_ice/ use_landice_water, landice_start_mon, landice_end_mon

end module g_forcing_param
! ====================================================================
module g_forcing_arrays
use o_param
  implicit none
  save    

  ! forcing arrays
  real(kind=WP), allocatable, dimension(:)         :: u_wind, v_wind 
  real(kind=WP), allocatable, dimension(:)         :: Tair, shum
  real(kind=WP), allocatable, dimension(:,:)       :: u_wind_t, v_wind_t 
  real(kind=WP), allocatable, dimension(:,:)       :: Tair_t, shum_t
  real(kind=WP), allocatable, dimension(:)         :: shortwave, longwave
  real(kind=WP), allocatable, dimension(:)         :: prec_rain, prec_snow
  real(kind=WP), allocatable, dimension(:)         :: runoff, evaporation, ice_sublimation
  real(kind=WP), allocatable, dimension(:)         :: cloudiness, press_air

#if defined (__oasis)
  real(kind=WP), target, allocatable, dimension(:) :: sublimation, evap_no_ifrac
  real(kind=WP), target, allocatable, dimension(:) :: tmp_sublimation, tmp_evap_no_ifrac !temporary flux fields
  real(kind=WP), target, allocatable, dimension(:) :: tmp_shortwave 			!(for flux correction) 
  real(kind=WP), allocatable, dimension(:)         :: atm_net_fluxes_north, atm_net_fluxes_south
  real(kind=WP), allocatable, dimension(:)         :: oce_net_fluxes_north, oce_net_fluxes_south
  real(kind=WP), allocatable, dimension(:)         :: flux_correction_north, flux_correction_south, flux_correction_total
#endif
  
  real(kind=WP), allocatable, dimension(:)         :: runoff_landice
  real(kind=WP)                                    :: landice_season(12)

  ! shortwave penetration
  real(kind=WP), allocatable, dimension(:)         :: chl
  real(kind=WP), allocatable, dimension(:,:)       :: sw_3d

  real(kind=WP), allocatable, dimension(:)         :: thdgr, thdgrsn, flice
  real(kind=WP), allocatable, dimension(:)         :: olat_heat, osen_heat, olwout
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
