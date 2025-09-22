module forcing_array_setup_interfaces
  interface
    subroutine forcing_array_setup(partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
    end subroutine forcing_array_setup
  end interface
end module forcing_array_setup_interfaces


module forcing_array_setup_dbgyre_interfaces
  interface
    subroutine forcing_array_setup_dbgyre(partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_mesh),   intent(in),    target :: mesh
      type(t_partit), intent(inout), target :: partit
    end subroutine forcing_array_setup_dbgyre
  end interface
end module forcing_array_setup_dbgyre_interfaces

! Adapted from FESOM code by Q. Wang. 
! Added the driving routine forcing_setup.
! S.D 05.04.12
! ==========================================================
subroutine forcing_setup(partit, mesh)
use g_CONFIG
use g_sbf, only: sbc_ini
use mod_mesh
USE MOD_PARTIT
USE MOD_PARSUP
use forcing_array_setup_interfaces
use forcing_array_setup_dbgyre_interfaces
implicit none
type(t_mesh),   intent(in),    target :: mesh
type(t_partit), intent(inout), target :: partit

  if (partit%mype==0) write(*,*) '****************************************************'
  if (use_ice) then
     call forcing_array_setup(partit, mesh)
#if !defined(__oasis) && !defined(__yac)
     call sbc_ini(partit, mesh)         ! initialize forcing fields
#endif
  endif 
  if ((toy_ocean) .AND. TRIM(which_toy)=="dbgyre" .AND. (use_sw_pene)) then
     call forcing_array_setup_dbgyre(partit, mesh)
  endif
end subroutine forcing_setup

! ==========================================================
subroutine forcing_array_setup_dbgyre(partit, mesh)
  use g_forcing_arrays
  use mod_mesh
  USE MOD_PARTIT
  implicit none
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer    :: n2

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  n2=myDim_nod2D+eDim_nod2D 

  allocate(chl(n2))
  allocate(sw_3d(nl,n2))
  chl=0.1_WP

end subroutine forcing_array_setup_dbgyre

! ==========================================================
subroutine forcing_array_setup(partit, mesh)
  !inializing forcing fields 
  use o_param
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_forcing_arrays
  use g_forcing_param
  use g_config
  use g_sbf, only: l_mslp, l_cloud
#if defined (__oasis)
  use cpl_driver, only : nrecv
#elif defined(__yac)
  use cpl_yac_driver, only : nrecv
#endif
  implicit none
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer    :: n2

! kh 19.02.21  
  integer    :: i
   
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  n2=myDim_nod2D+eDim_nod2D      
  ! Allocate memory for atmospheric forcing 
  allocate(shortwave(n2), longwave(n2))
  shortwave=0.0_WP
  longwave=0.0_WP
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

  prec_rain=0.0_WP
  prec_snow=0.0_WP
  !allocate(u_wind(n2), v_wind(n2))
  u_wind=0.0_WP
  v_wind=0.0_WP
  allocate(Tair(n2), shum(n2))
  Tair=0.0_WP
  shum=0.0_WP
  allocate(runoff(n2), evaporation(n2),ice_sublimation(n2))
  runoff=0.0_WP
  evaporation = 0.0_WP
  ice_sublimation = 0.0_WP
#if defined (__oasis) || defined (__ifsinterface) || defined (__yac)
  allocate(sublimation(n2), evap_no_ifrac(n2))
  sublimation=0.0_WP
  evap_no_ifrac=0.0_WP
#endif
#if defined (__oasis) || defined (__yac)
  allocate(tmp_sublimation(n2),tmp_evap_no_ifrac(n2), tmp_shortwave(n2))
  allocate(atm_net_fluxes_north(nrecv), atm_net_fluxes_south(nrecv))
  allocate(oce_net_fluxes_north(nrecv), oce_net_fluxes_south(nrecv))
  allocate(flux_correction_north(nrecv), flux_correction_south(nrecv))
  allocate(flux_correction_total(nrecv))
  tmp_sublimation = 0.0_WP
  tmp_evap_no_ifrac = 0.0_WP
  tmp_shortwave = 0.0_WP
  atm_net_fluxes_north=0.0_WP
  atm_net_fluxes_south=0.0_WP
  oce_net_fluxes_north=0.0_WP
  oce_net_fluxes_south=0.0_WP
  flux_correction_north=0.0_WP
  flux_correction_south=0.0_WP
  flux_correction_total=0.0_WP  
#endif 

#if defined (__oasis) || defined (__ifsinterface)
  allocate(residualifwflx(n2))
  residualifwflx = 0.0_WP

  !---wiso-code
  IF (lwiso) THEN
    allocate(tmp_iii1(n2), tmp_iii2(n2), tmp_iii3(n2))
    allocate(www1(n2), www2(n2), www3(n2), iii1(n2), iii2(n2), iii3(n2))
    tmp_iii1=0.0_WP
    tmp_iii2=0.0_WP
    tmp_iii3=0.0_WP
    www1=0.0_WP
    www2=0.0_WP
    www3=0.0_WP
    iii1=0.0_WP
    iii2=0.0_WP
    iii3=0.0_WP
  END IF
  !---wiso-code-end
#endif 

! Temp storage for averaging
!!PS   allocate(aver_temp(n2))

!!PS   allocate(Tair_mo(n2),shum_mo(n2))
!!PS   Tair_mo=0.0_WP !!PS
!!PS   shum_mo=0.0_WP !!PS
!!PS   
!!PS   allocate(mo_index(n2),dv10_mo(n2),dv10(n2),auxt1(n2),auxt2(n2),auxt3(n2))  
!!PS   mo_index = 0.0_WP
!!PS   dv10_mo  =0.0_WP
!!PS   dv10     =0.0_WP
!!PS   auxt1    =0.0_WP
!!PS   auxt2    =0.0_WP
!!PS   auxt3    =0.0_WP

  if (l_cloud) then
     allocate(cloudiness(n2))
     cloudiness=0.0_WP
  end if
!  if (l_mslp) then
     allocate(press_air(n2))
     press_air=0.0_WP
!  end if
 
  allocate(u_wind_t(2,n2),v_wind_t(2,n2))
  allocate(Tair_t(2,n2), shum_t(2,n2))
  u_wind_t=0.0_WP
  v_wind_t=0.0_WP
  Tair_t=0.0_WP
  shum_t=0.0_WP

  if(use_landice_water) then
    allocate(runoff_landice(n2))
    runoff_landice=0.0_WP
  end if
  !---age-code-begin
  if(use_age_tracer) then
    allocate(age_tracer_loc_index(n2))
    age_tracer_loc_index=0
  end if
  !---age-code-end
 
  ! shortwave penetration
  if(use_sw_pene) then
    allocate(chl(n2))
    allocate(sw_3d(nl,n2))
    chl=0.1_WP
  endif

  !for ice diagnose
  if(use_ice) then
!     allocate(thdgr(n2), thdgrsn(n2))
    allocate(flice(n2))
    allocate(olat_heat(n2), osen_heat(n2), olwout(n2))
!     thdgr=0.0_WP
!     thdgrsn=0.0_WP
    flice=0.0_WP
    olat_heat=0.0_WP
    osen_heat=0.0_WP
    olwout=0.0_WP
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
