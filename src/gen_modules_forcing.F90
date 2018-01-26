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

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice


  ! *** forcing source and type ***
  character(10)                 :: wind_data_source='CORE2'
  character(10)                 :: rad_data_source='CORE2'
  character(10)                 :: precip_data_source='CORE2'
  character(10)                 :: runoff_data_source='CORE2'
  character(10)                 :: sss_data_source='CORE2'
  integer                       :: wind_ttp_ind=1
  integer                       :: rad_ttp_ind=2
  integer                       :: precip_ttp_ind=3
  integer                       :: runoff_ttp_ind=0
  integer                       :: sss_ttp_ind=4
  logical                       :: use_virt_salt ! will be set TRUE in case of which_ALE='linfs', otherwise FALSE

  namelist /forcing_source/ wind_data_source, rad_data_source, &
       precip_data_source, runoff_data_source, sss_data_source, &
       wind_ttp_ind, rad_ttp_ind, precip_ttp_ind, &
       runoff_ttp_ind, sss_ttp_ind

  ! *** coefficients in bulk formulae ***
  logical                       :: AOMIP_drag_coeff=.false.
  logical                       :: ncar_bulk_formulae=.false.

  namelist /forcing_bulk/ AOMIP_drag_coeff, ncar_bulk_formulae

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
  real(kind=WP), allocatable, dimension(:)         :: runoff, evaporation
  real(kind=WP), allocatable, dimension(:)         :: cloudiness, Pair

#if defined (__oasis)
  real(kind=8), target, allocatable, dimension(:) :: sublimation, evap_no_ifrac
  real(kind=8), target, allocatable, dimension(:) :: tmp_sublimation, tmp_evap_no_ifrac !temporary flux fields
  real(kind=8), target, allocatable, dimension(:) :: tmp_shortwave 			!(for flux correction) 
  real(kind=8), allocatable, dimension(:)         :: atm_net_fluxes_north, atm_net_fluxes_south
  real(kind=8), allocatable, dimension(:)         :: oce_net_fluxes_north, oce_net_fluxes_south
  real(kind=8), allocatable, dimension(:)         :: flux_correction_north, flux_correction_south, flux_correction_total
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
  
  real*8,allocatable,dimension(:) ::  aver_temp
end module g_forcing_arrays
! ====================================================================
module g_forcing_index
  ! Control parameters for updating forcing in time 	
  !
  ! Coded by Qiang Wang
  ! Uncomment the part required for interpolation in time by Qiang Wang, 12,Oct. 2010
  ! Recover the interpolation for 6hourly data, Qiang, 07.02.2012
  ! Reviewed by ??
  ! --------------------------------------------------------
  use o_param
  implicit none
  save

  ! arrays for temporal interpolation
  integer                                         :: update_forcing_flag(4)
  integer                                         :: forcing_rec(4)
  real(kind=WP)                                   :: interp_coef(4)

 contains
! ===================================================================
  subroutine forcing_index
    use o_param
    use g_clock
    use g_PARSUP, only : par_ex
    implicit none

    real(kind=8)          :: sixhour_sec
    real(kind=8)          :: oneday_sec
    real(kind=8)          :: modtimeold

    data sixhour_sec /21600.0/, oneday_sec /86400.0/

    modtimeold=mod(timeold,oneday_sec)

    ! if update forcing or not
    update_forcing_flag=0
    if(mod(timeold, sixhour_sec)==0.0) update_forcing_flag(1)=1
    if(modtimeold==0.0) update_forcing_flag(2)=1
    if(day_in_month==1 .and. modtimeold==0.0) update_forcing_flag(3:4)=1

    ! which record should be used as the first one in interpolation
    forcing_rec(1) = 1+int(modtimeold/sixhour_sec)+4*(daynew-1)
    forcing_rec(2) = daynew
    forcing_rec(3) = month
    forcing_rec(4) = month

  end subroutine forcing_index

end module g_forcing_index
!
!----------------------------------------------------------------------------
!
module g_forcing_interp
  ! This module initializes the weights for interpolating 
  ! forcing data Tair, humidity, wind velocities,
  ! precipitation, radiations, etc.
  ! Based on the assumption that forcing data are on the T62 NCEP/NCAR grid
  !
  ! Follows the old version, but uses different 
  ! way of weighting (Currently using bilinear interpolation).
  !
  ! Modified by Qiang Wang for bilinear interpolation
  ! Reviewed by ??
  !------------------------------------------------------------------	
  implicit none
  integer, allocatable      :: lint_ind(:,:,:)
  real(kind=8), allocatable :: lint_weight(:,:)

  contains

  !------------------------------------------------------------------
  subroutine  init_forcing_interp

    ! routine to calculate neighbours and weights for linear interpolation
    !
    ! required information
    !    xmod(nmp)  longitudes of model point on geographical grid in degree (no order assumed)
    !    ymod(nmp)  latitudes of model points on geographical grid in degree (no order assumed)
    !         nmp   number of model points where data are needed
    !    cx(ni)     longitudes of data points on regular geographical grid
    !               by now must be in range[0,..360] in ascending order
    !    cy(nj)     latitudes of data points on regular geographical grid 
    !               by now must be in range [-90,..90] in ascending order
    !
    ! OUTPUT
    !    lint_ind(4,2,nmp)   /i,j/ adress for four neighbors of each model node
    !    lint_weight(4,nmp)  interpolation weights

    use o_mesh
    use o_param
    use g_config
    use g_parsup
    implicit none

    integer, parameter 	:: ni=192, nj=94  ! NCEP and CORE are on the same grid.
    integer     	:: i, ii, j, n, row, n2
    real(kind=8)      	:: rlon, rlat, aux
    real(kind=8)      	:: cx(ni), cy(nj)
    real(kind=8)      	:: xmod(myDim_nod2D+eDim_nod2D), ymod(myDim_nod2D+eDim_nod2D)
    real(kind=8)        :: wt(4)

    n2=myDim_nod2D+eDim_nod2D     

    allocate(lint_ind(4,2,n2))   
    allocate(lint_weight(4,n2))  

    ! NCEP/CORE latitude
    data cy /-88.542, -86.6531, -84.7532, -82.8508,  -80.9473,  -79.0435, &  
            -77.1394, -75.2351, -73.3307, -71.4262,  -69.5217,  -67.6171, &  
            -65.7125, -63.8079, -61.9033, -59.9986,  -58.0939,  -56.1893, &  
            -54.2846, -52.3799, -50.4752, -48.5705,  -46.6658,  -44.7611, &  
            -42.8564, -40.9517, -39.0470, -37.1422,  -35.2375,  -33.3328, &  
            -31.4281, -29.5234, -27.6186, -25.7139,  -23.8092,  -21.9044, &  
            -19.9997, -18.0950, -16.1902, -14.2855,  -12.3808, -10.47604, &  
            -8.57131, -6.66657, -4.76184,  -2.8571, -0.952368,  0.952368, &  
              2.8571,  4.76184,  6.66657,  8.57131,  10.47604,   12.3808, &  
             14.2855,  16.1902,   18.095,  19.9997,   21.9044,   23.8092, &  
             25.7139,  27.6186,  29.5234,  31.4281,   33.3328,   35.2375, &  
             37.1422,   39.047,  40.9517,  42.8564,   44.7611,   46.6658, &  
             48.5705,  50.4752,  52.3799,  54.2846,   56.1893,   58.0939, &  
             59.9986,  61.9033,  63.8079,  65.7125,   67.6171,   69.5217, &  
             71.4262,  73.3307,  75.2351,  77.1394,   79.0435,   80.9473, &  
             82.8508,  84.7532,  86.6531,  88.542 /

    ! NCEP/CORE longitude
    cx(1)=0.0
    do i=2,ni
       cx(i)=cx(i-1)+1.875
    enddo

    ! in the following we need cx and cy in unit radian
    cx=cx*rad
    cy=cy*rad

    ! model grid geographical coordinates (in radian, between 0 and 2*pi)
    do row=1,n2 
       xmod(row)=geo_coord_nod2D(1,row)
       ymod(row)=geo_coord_nod2d(2,row)
       if (xmod(row)<0.0) xmod(row)=2*pi+xmod(row)	
    enddo

    ! linear interpolation: nodes and weight
    do row=1,n2
       if(xmod(row)<cx(ni)) then
          do i=1,ni
             if(xmod(row)<cx(i)) then
                lint_ind(1,1,row)=i-1
                lint_ind(2,1,row)=i-1
                lint_ind(3,1,row)=i
                lint_ind(4,1,row)=i
                aux=(cx(i)-xmod(row))/(cx(i)-cx(i-1))
                exit
             end if
          end do
       else
          lint_ind(1,1,row)=ni
          lint_ind(2,1,row)=ni
          lint_ind(3,1,row)=1
          lint_ind(4,1,row)=1
          aux=(360.0*rad-xmod(row))/(360.0*rad-cx(ni))
       end if
       wt(1)=aux
       wt(2)=aux
       wt(3)=1.0_8-aux
       wt(4)=1.0_8-aux

       if(ymod(row)<cy(nj)) then
          do j=1,nj
             if(ymod(row)<cy(j)) then
                lint_ind(1,2,row)=j-1
                lint_ind(2,2,row)=j
                lint_ind(3,2,row)=j
                lint_ind(4,2,row)=j-1
                aux=(cy(j)-ymod(row))/(cy(j)-cy(j-1))
                exit
             end if
          end do
       else
          lint_ind(1,2,row)=nj
          lint_ind(2,2,row)=nj
          lint_ind(3,2,row)=nj
          lint_ind(4,2,row)=nj
          aux=1.0_8
       end if
       lint_weight(1,row)=wt(1)*aux
       lint_weight(2,row)=wt(2)*(1.0_8-aux)
       lint_weight(3,row)=wt(3)*(1.0_8-aux)
       lint_weight(4,row)=wt(4)*aux
    end do

    if(mype==0)  write(*,*) 'weights for interpolating forcing / 2D fields computed'
    return     
  end subroutine init_forcing_interp


  !------------------------------------------------------------------------
  subroutine forcing_linear_ip(zd,idim,jdim,ind,weights,zi,nmpt)
    !  this subroutine interpolates data using prepared weights- 
    !  see subroutine init_forcing_interp
    !
    !  INPUT
    !        zd(idim,jdim)         available data set
    !        nmpt                  number of model positions, where data are wanted
    !        indx(4,2,nmpt)        i,j index of four neighbors of each node
    !        weights(4,nmpt)       interpolation weights
    !        
    !  OUTPUT
    !        zi(nmpt)              array of interpolated values for model points

    use g_parsup
    implicit none                                             

    integer      :: idim, jdim, nmpt                          
    integer      :: ind(4,2,nmpt)
    integer      :: i, n                      
    real(kind=8) :: zd(idim,jdim), zi(nmpt)
    real(kind=8) :: weights(4,nmpt)
    real(kind=8) :: fx

    do n=1,nmpt
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
    enddo

    return
  end subroutine forcing_linear_ip

end module g_forcing_interp
