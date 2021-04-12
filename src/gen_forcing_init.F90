! Adapted from FESOM code by Q. Wang. 
! Added the driving routine forcing_setup.
! S.D 05.04.12
! ==========================================================
subroutine forcing_setup
use g_parsup
use g_CONFIG
use g_sbf, only: sbc_ini
implicit none

  if (mype==0) write(*,*) '****************************************************'
  if (use_ice) then
     call forcing_array_setup
#ifndef __oasis
     call sbc_ini         ! initialize forcing fields
#endif
  endif 
end subroutine forcing_setup
! ==========================================================
subroutine forcing_array_setup
  !inializing forcing fields 
  use o_param
  use o_mesh
  use i_arrays
  use g_forcing_arrays
  use g_forcing_param
  use g_parsup
  use g_config
  use g_sbf, only: l_mslp, l_cloud
#if defined (__oasis)
  use cpl_driver, only : nrecv
#endif   
  implicit none

  integer    :: n2

! kh 19.02.21  
  integer    :: i
   
  n2=myDim_nod2D+eDim_nod2D
  ! Allocate memory for atmospheric forcing 
  allocate(shortwave(n2), longwave(n2))
  allocate(prec_rain(n2), prec_snow(n2))

! kh 19.02.21
  if (ib_async_mode == 0) then
      allocate(u_wind(n2), v_wind(n2))
      allocate(u_wind_ib(n2), v_wind_ib(n2))
  else
! kh 19.02.21 support "first touch" idea
!$omp parallel sections num_threads(2)
!$omp section
      allocate(u_wind(n2), v_wind(n2))
      do i = 1, n2
          u_wind(i) = 0._WP
          v_wind(i) = 0._WP
      end do
!$omp section
      allocate(u_wind_ib(n2), v_wind_ib(n2))
      do i = 1, n2
          u_wind_ib(i) = 0._WP
          v_wind_ib(i) = 0._WP
      end do
!$omp end parallel sections
  end if

  allocate(Tair(n2), shum(n2))
  allocate(runoff(n2), evaporation(n2))

#if defined (__oasis)
  allocate(sublimation(n2), evap_no_ifrac(n2))
  allocate(tmp_sublimation(n2),tmp_evap_no_ifrac(n2), tmp_shortwave(n2))
  allocate(atm_net_fluxes_north(nrecv), atm_net_fluxes_south(nrecv))
  allocate(oce_net_fluxes_north(nrecv), oce_net_fluxes_south(nrecv))
  allocate(flux_correction_north(nrecv), flux_correction_south(nrecv))
  allocate(flux_correction_total(nrecv))
  sublimation=0.
  evap_no_ifrac=0.
  tmp_sublimation = 0.
  tmp_evap_no_ifrac = 0.
  tmp_shortwave = 0.
  atm_net_fluxes_north=0.
  atm_net_fluxes_south=0.
  oce_net_fluxes_north=0.
  oce_net_fluxes_south=0.
  flux_correction_north=0.
  flux_correction_south=0.
  flux_correction_total=0.  
#endif 


! Temp storage for averaging
  allocate(aver_temp(n2))
  shortwave=0.
  longwave=0.
  prec_rain=0.
  prec_snow=0.
  u_wind=0.
  v_wind=0.
  Tair=0.
  shum=0.
  runoff=0.

  if (l_cloud) then
     allocate(cloudiness(n2))
     cloudiness=0.
  end if
  if (l_mslp) then
     allocate(Pair(n2))
     Pair=0.
  end if
 
  allocate(u_wind_t(2,n2),v_wind_t(2,n2))
  allocate(Tair_t(2,n2), shum_t(2,n2))
  u_wind_t=0.
  v_wind_t=0.
  Tair_t=0.
  shum_t=0.

  if(use_landice_water) then
    allocate(runoff_landice(n2))
    runoff_landice=0.0
  end if
 
  ! shortwave penetration
  if(use_sw_pene) then
    allocate(chl(n2))
    allocate(sw_3d(nl,n2))
    chl=.1
  endif

  !for ice diagnose
  if(use_ice) then
    allocate(thdgr(n2), thdgrsn(n2), flice(n2))
    allocate(olat_heat(n2), osen_heat(n2), olwout(n2))
    thdgr=0.
    thdgrsn=0.
    flice=0.
    olat_heat=0.
    osen_heat=0.
    olwout=0.
  endif 

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  allocate(Ce_atm_oce_arr(n2))
  allocate(Ch_atm_oce_arr(n2))
  Cd_atm_oce_arr=Cd_atm_oce
  Ce_atm_oce_arr=Ce_atm_oce 
  Ch_atm_oce_arr=Ch_atm_oce
  if(use_ice) then
    allocate(Cd_atm_ice_arr(n2)) 
    Cd_atm_ice_arr=Cd_atm_ice   
  endif

  if(mype==0) write(*,*) 'forcing arrays have been set up'  
end subroutine forcing_array_setup
!
!----------------------------------------------------------------------
!
