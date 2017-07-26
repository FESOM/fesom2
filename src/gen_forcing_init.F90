! Adapted from FESOM code by Q. Wang. 
! Added the driving routine forcing_setup.
! S.D 05.04.12
! ==========================================================
subroutine forcing_setup
! forcing: arrays, initialization, interpolation preparation  
use g_forcing_index
use g_parsup
use g_CONFIG
use g_forcing_interp
implicit none

  if (mype==0) write(*,*) '****************************************************'
  call forcing_index
  if (use_ice) then
     call forcing_array_setup
     call init_forcing_interp      ! calculates the forcing interpolation weights
     call init_atm_forcing         ! initialize forcing fields
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
  implicit none

  integer    :: n2
   
  n2=myDim_nod2D+eDim_nod2D      
  ! Allocate memory for atmospheric forcing 
  allocate(shortwave(n2), longwave(n2))
  allocate(prec_rain(n2), prec_snow(n2))
  allocate(u_wind(n2), v_wind(n2))
  allocate(Tair(n2), shum(n2))
  allocate(runoff(n2), evaporation(n2))
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
  if(rad_data_source=='NCEP') then
    allocate(cloudiness(n2), Pair(n2))
    cloudiness=0.
    Pair=0.
  end if

  if(wind_ttp_ind==1) then
    allocate(u_wind_t(2,n2),v_wind_t(2,n2))
    allocate(Tair_t(2,n2), shum_t(2,n2))
    u_wind_t=0.
    v_wind_t=0.
    Tair_t=0.
    shum_t=0.
  end if

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
