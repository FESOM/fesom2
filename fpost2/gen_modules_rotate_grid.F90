module g_rotate_grid
  use g_config
  implicit none
  save
  real(kind=8)        :: rotate_matrix(3,3), rotate_matrix_inv(3,3)

contains


  !----------------------------------------------------------------
  !
  subroutine calculate_rotate_matrix
    ! A, B, G [radian] are Euler angles.
    ! The convention (A, B, G) used here is: the first rotation is by an
    ! angle A around z-axis, the second is by an angle B about the new 
    ! x-axis, and the third is by an angle G about the new z-axis.   
    implicit none

    real(kind=8)      :: al, be, ga, det

    al=alphaEuler*rad
    be=betaEuler*rad
    ga=gammaEuler*rad

    ! rotation matrix
    rotate_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
    rotate_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
    rotate_matrix(1,3)=sin(ga)*sin(be)
    rotate_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
    rotate_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
    rotate_matrix(2,3)=cos(ga)*sin(be)
    rotate_matrix(3,1)=sin(be)*sin(al) 
    rotate_matrix(3,2)=-sin(be)*cos(al)  
    rotate_matrix(3,3)=cos(be)
    ! inverse for vector_r2g
    call matrix_inverse_3x3(rotate_matrix, rotate_matrix_inv, det)
  end subroutine calculate_rotate_matrix
  !
  !----------------------------------------------------------------
  !
  subroutine r2g(lon, lat, rlon, rlat)
    ! Convert the rotated coordinates to geographical coordinates  
    ! lon, lat		:: [radian] geographical coordinates
    ! rlon, rlat	:: [radian] rotated coordinates
    !
    implicit none
    real(kind=8), intent(out)      :: lon, lat
    real(kind=8), intent(in)       :: rlon, rlat
    real(kind=8)                   :: xr, yr, zr, xg, yg, zg
    !
    ! Rotated Cartesian coordinates:
    xr=cos(rlat)*cos(rlon)
    yr=cos(rlat)*sin(rlon)
    zr=sin(rlat)

    ! Geographical Cartesian coordinates:
    xg=rotate_matrix(1,1)*xr + rotate_matrix(2,1)*yr + rotate_matrix(3,1)*zr
    yg=rotate_matrix(1,2)*xr + rotate_matrix(2,2)*yr + rotate_matrix(3,2)*zr  
    zg=rotate_matrix(1,3)*xr + rotate_matrix(2,3)*yr + rotate_matrix(3,3)*zr  

    ! Geographical coordinates:
    lat=asin(zg)
    if(yg==0. .and. xg==0.) then
       lon=0.0     ! exactly at the poles
    else
       lon=atan2(yg,xg)
    end if
  end subroutine r2g
  !
  !----------------------------------------------------------------
  !
  subroutine g2r(lon, lat, rlon, rlat)
    ! Convert the geographical coordinates to rotated coordinates  
    ! lon, lat		:: [radian] geographical coordinates
    ! rlon, rlat	:: [radian] rotated coordinates
    !
    implicit none
    real(kind=8), intent(in)       :: lon, lat
    real(kind=8), intent(out)      :: rlon, rlat
    real(kind=8)                   :: xr, yr, zr, xg, yg, zg
    !
    ! geographical Cartesian coordinates:
    xg=cos(lat)*cos(lon)
    yg=cos(lat)*sin(lon)
    zg=sin(lat)

    ! rotated Cartesian coordinates:
    xr=rotate_matrix(1,1)*xg + rotate_matrix(1,2)*yg + rotate_matrix(1,3)*zg
    yr=rotate_matrix(2,1)*xg + rotate_matrix(2,2)*yg + rotate_matrix(2,3)*zg  
    zr=rotate_matrix(3,1)*xg + rotate_matrix(3,2)*yg + rotate_matrix(3,3)*zg  

    ! rotated coordinates:
    rlat=asin(zr)
    if(yr==0. .and. xr==0.) then
       rlon=0.0     ! exactly at the poles
    else
       rlon=atan2(yr,xr)
    end if
  end subroutine g2r
  !
  !--------------------------------------------------------------------
  !
  subroutine vector_g2r(tlon, tlat, lon, lat, flag_coord)
    ! rotate 2d vector (tlon, tlat) to be in the rotated coordinates
    ! tlon, tlat (in)	:: lon & lat components of a vector in geographical coordinates 
    !            (out)	:: lon & lat components of the vector in rotated coordinates              
    ! lon, lat	        :: [radian] coordinates
    ! flag_coord        :: 1, (lon,lat) is the geographical coord.; else, rotated coord.
    !
    implicit none
    integer, intent(in)           :: flag_coord
    real(kind=8), intent(inout)   :: tlon, tlat
    real(kind=8), intent(in)      :: lon, lat
    real(kind=8)                  :: rlon, rlat, glon, glat
    real(kind=8)		  :: txg, tyg, tzg, txr, tyr, tzr
    !
    ! geographical coordinate
    if(flag_coord==1) then  ! input is geographical coordinates
       glon=lon
       glat=lat
       call g2r(glon,glat,rlon,rlat)
    else                    ! input is rotated coordinates 
       rlon=lon
       rlat=lat
       call r2g(glon,glat,rlon,rlat)
    end if
    !
    ! vector in Cartesian
    txg=-tlat*sin(glat)*cos(glon)-tlon*sin(glon)
    tyg=-tlat*sin(glat)*sin(glon)+tlon*cos(glon)
    tzg=tlat*cos(glat)
    !
    ! vector in rotated Cartesian
    txr=rotate_matrix(1,1)*txg + rotate_matrix(1,2)*tyg + rotate_matrix(1,3)*tzg 
    tyr=rotate_matrix(2,1)*txg + rotate_matrix(2,2)*tyg + rotate_matrix(2,3)*tzg 
    tzr=rotate_matrix(3,1)*txg + rotate_matrix(3,2)*tyg + rotate_matrix(3,3)*tzg 
    !
    ! vector in rotated coordinate
    tlat=-sin(rlat)*cos(rlon)*txr - sin(rlat)*sin(rlon)*tyr + cos(rlat)*tzr
    tlon=-sin(rlon)*txr + cos(rlon)*tyr

  end subroutine vector_g2r
  !
  subroutine vector_r2g(tlon, tlat, lon, lat, flag_coord)
   ! rotate 2d vector (tlon, tlat) to be in geo. coordinates
   ! tlon, tlat (in)    :: lon & lat components of a vector in rotated coordinates
   !            (out)   :: lon & lat components of the vector in geo. coordinates
   ! lon, lat           :: [radian] coordinates
   ! flag_coord         :: 1, (lon,lat) is the geographical coord.; else, rotated coord.
   !
   implicit none
   integer, intent(in)           :: flag_coord
   real(kind=8), intent(inout)   :: tlon, tlat
   real(kind=8), intent(in)      :: lon, lat
   real(kind=8)                  :: rlon, rlat, glon, glat
   real(kind=8)                :: txg, tyg, tzg, txr, tyr, tzr
   !
   ! geographical coordinate
   if(flag_coord==1) then  ! input is geographical coordinates
      glon=lon
      glat=lat
      call g2r(glon,glat,rlon,rlat)
   else                    ! input is rotated coordinates
      rlon=lon
      rlat=lat
      call r2g(glon,glat,rlon,rlat)
   end if
   !
   ! vector rotated Cartesian
   txg=-tlat*sin(rlat)*cos(rlon)-tlon*sin(rlon)
   tyg=-tlat*sin(rlat)*sin(rlon)+tlon*cos(rlon)
   tzg=tlat*cos(rlat)
   !
   ! vector in geo Cartesian
   txr=rotate_matrix_inv(1,1)*txg + rotate_matrix_inv(1,2)*tyg + rotate_matrix_inv(1,3)*tzg
   tyr=rotate_matrix_inv(2,1)*txg + rotate_matrix_inv(2,2)*tyg + rotate_matrix_inv(2,3)*tzg
   tzr=rotate_matrix_inv(3,1)*txg + rotate_matrix_inv(3,2)*tyg + rotate_matrix_inv(3,3)*tzg
   !
   ! vector in geo coordinate
   tlat=-sin(glat)*cos(glon)*txr - sin(glat)*sin(glon)*tyr + cos(glat)*tzr
   tlon=-sin(glon)*txr + cos(glon)*tyr
  end subroutine vector_r2g
!  
end module g_rotate_grid
