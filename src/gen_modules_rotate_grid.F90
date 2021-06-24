! Routines needed to support displaced poles: 
! The new pole position is set with 
! alphaEuler, betaEuler and gammaEuler. The default values 
! alphaEuler=50.   [degree] Euler angles, convention:
! betaEuler=15.	   [degree] first around z, then around new x,
! gammaEuler=-90.  [degree] then around new z.
!
! A helpful animation may be found online here:
! https://en.wikipedia.org/wiki/Euler_angles
!
! The first two define the new pole position 
! as phi_p=alphaEuler-90, theta_p=90-betaEuler.
! The third, gammaEuler, is in reality irrelevant and just
! sets the position of new meridian. For -90 it is just the previous
! phi_p. We can live without gammaEuler, but we keep it to be similar
! to FESOM
! Q. Wang

module g_rotate_grid
  use g_config
  implicit none
  save
  real(kind=WP)        :: r2g_matrix(3,3)

 contains
  
  !
  !----------------------------------------------------------------
  !
  subroutine set_mesh_transform_matrix
    ! A, B, G [radian] are Euler angles.
    ! The convention (A, B, G) used here is: the first rotation is by an
    ! angle A around z-axis, the second is by an angle B about the new 
    ! x-axis, and the third is by an angle G about the new z-axis.   
    use o_PARAM
    use g_PARSUP, only : mype
    implicit none
    real(kind=WP)      :: al, be, ga

    al=alphaEuler
    be=betaEuler
    ga=gammaEuler

    ! rotation matrix
    r2g_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
    r2g_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
    r2g_matrix(1,3)=sin(ga)*sin(be)
    r2g_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
    r2g_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
    r2g_matrix(2,3)=cos(ga)*sin(be)
    r2g_matrix(3,1)=sin(be)*sin(al) 
    r2g_matrix(3,2)=-sin(be)*cos(al)  
    r2g_matrix(3,3)=cos(be)
    if(mype==0) write(*,*) 'rotation matrix for rotated model grids prepared'
  end subroutine set_mesh_transform_matrix
!
!----------------------------------------------------------------
!
  subroutine r2g(glon, glat, rlon, rlat)
    ! Transform from the mesh (rotated) coordinates 
    !           to geographical coordinates  
    ! glon, glat        :: [radian] geographical coordinates
    ! rlon, rlat        :: [radian] rotated coordinates
    !
    implicit none
    real(kind=WP), intent(out)      :: glon, glat
    real(kind=WP), intent(in)       :: rlon, rlat
    real(kind=WP)                   :: xr, yr, zr, xg, yg, zg
    !
    ! Rotated Cartesian coordinates:
    xr=cos(rlat)*cos(rlon)
    yr=cos(rlat)*sin(rlon)
    zr=sin(rlat)

    ! Geographical Cartesian coordinates:
    xg=r2g_matrix(1,1)*xr + r2g_matrix(2,1)*yr + r2g_matrix(3,1)*zr
    yg=r2g_matrix(1,2)*xr + r2g_matrix(2,2)*yr + r2g_matrix(3,2)*zr  
    zg=r2g_matrix(1,3)*xr + r2g_matrix(2,3)*yr + r2g_matrix(3,3)*zr  

    ! Geographical lon-lat coordinates:
    glat=asin(zg)
    if(yg==0. .and. xg==0._WP) then
       glon=0.0_WP     ! exactly at the poles
    else
       glon=atan2(yg,xg)
    end if
  end subroutine r2g
  !
  !----------------------------------------------------------------
  !
  subroutine g2r(glon, glat, rlon, rlat)
    ! Transform the geographical coordinates to rotated coordinates  
    ! glon, glat        :: [radian] geographical coordinates
    ! rlon, rlat        :: [radian] rotated coordinates
    !
    implicit none
    real(kind=WP), intent(in)       :: glon, glat
    real(kind=WP), intent(out)      :: rlon, rlat
    real(kind=WP)                   :: xr, yr, zr, xg, yg, zg
    !
    ! geographical Cartesian coordinates:
    xg=cos(glat)*cos(glon)
    yg=cos(glat)*sin(glon)
    zg=sin(glat)

    ! rotated Cartesian coordinates:
    xr=r2g_matrix(1,1)*xg + r2g_matrix(1,2)*yg + r2g_matrix(1,3)*zg
    yr=r2g_matrix(2,1)*xg + r2g_matrix(2,2)*yg + r2g_matrix(2,3)*zg  
    zr=r2g_matrix(3,1)*xg + r2g_matrix(3,2)*yg + r2g_matrix(3,3)*zg  

    ! rotated coordinates:
    rlat=asin(zr)
    if(yr==0._WP .and. xr==0._WP) then
       rlon=0.0_WP     ! exactly at the poles
    else
       rlon=atan2(yr,xr)
    end if
  end subroutine g2r
  !
  !--------------------------------------------------------------------
  !
  subroutine vector_g2r(tlon, tlat, lon, lat, flag_coord)
    ! Transform a 2d vector with components (tlon, tlat) in
    ! geographical coordinates to the rotated mesh coordinates
    ! tlon, tlat (in)	:: lon & lat comp. of a vector in geo. coord. 
    !            (out)	:: lon & lat comp. in rot. coord.              
    ! lon, lat	        :: [radian] coordinates of vector position
    ! flag_coord        :: 1, (lon,lat)= geo. coord., else, rot. coord.
    !
    implicit none
    integer, intent(in)           :: flag_coord
    real(kind=WP), intent(inout)  :: tlon, tlat
    real(kind=WP), intent(in)     :: lon, lat
    real(kind=WP)                 :: rlon, rlat, glon, glat
    real(kind=WP)		  :: txg, tyg, tzg, txr, tyr, tzr
    !
    ! geographical coordinates
    if(flag_coord==1) then  ! input is in geographical coordinates
       glon=lon
       glat=lat
       call g2r(glon,glat,rlon,rlat)
    else                    ! input is in rotated coordinates 
       rlon=lon
       rlat=lat
       call r2g(glon,glat,rlon,rlat)
    end if
    !
    ! vector in Cartesian geo. coordinates
    txg=-tlat*sin(glat)*cos(glon)-tlon*sin(glon)
    tyg=-tlat*sin(glat)*sin(glon)+tlon*cos(glon)
    tzg=tlat*cos(glat)
    !
    ! vector in rot. Cartesian coordinates
    txr=r2g_matrix(1,1)*txg + r2g_matrix(1,2)*tyg + r2g_matrix(1,3)*tzg 
    tyr=r2g_matrix(2,1)*txg + r2g_matrix(2,2)*tyg + r2g_matrix(2,3)*tzg 
    tzr=r2g_matrix(3,1)*txg + r2g_matrix(3,2)*tyg + r2g_matrix(3,3)*tzg 
    !
    ! vector in rotated coordinates
    tlat=-sin(rlat)*cos(rlon)*txr - sin(rlat)*sin(rlon)*tyr + cos(rlat)*tzr
    tlon=-sin(rlon)*txr + cos(rlon)*tyr

  end subroutine vector_g2r
  !
  !----------------------------------------------------------------------------
  !
  subroutine vector_r2g(tlon, tlat, lon, lat, flag_coord)
   ! rotate 2d vector (tlon, tlat) to be in geo. coordinates
   ! tlon, tlat (in)    :: lon & lat components of a vector in rotated coordinates
   !            (out)   :: lon & lat components of the vector in geo. coordinates
   ! lon, lat           :: [radian] coordinates
   ! flag_coord         :: 1, (lon,lat) is the geographical coord.; else, rotated coord.
   !
   implicit none
   integer, intent(in)            :: flag_coord
   real(kind=WP), intent(inout)   :: tlon, tlat
   real(kind=WP), intent(in)      :: lon, lat
   real(kind=WP)                  :: rlon, rlat, glon, glat
   real(kind=WP)                  :: txg, tyg, tzg, txr, tyr, tzr
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
   txr=r2g_matrix(1,1)*txg + r2g_matrix(2,1)*tyg + r2g_matrix(3,1)*tzg
   tyr=r2g_matrix(1,2)*txg + r2g_matrix(2,2)*tyg + r2g_matrix(3,2)*tzg
   tzr=r2g_matrix(1,3)*txg + r2g_matrix(2,3)*tyg + r2g_matrix(3,3)*tzg
   !
   ! vector in geo coordinate
   tlat=-sin(glat)*cos(glon)*txr - sin(glat)*sin(glon)*tyr + cos(glat)*tzr
   tlon=-sin(glon)*txr + cos(glon)*tyr
  end subroutine vector_r2g
!===================================================================
  subroutine trim_cyclic(b)
    use o_PARAM
    use g_config
    implicit none
    real(kind=WP) :: b
    if(b> cyclic_length/2.0_WP) b=b-cyclic_length
    if(b<-cyclic_length/2.0_WP) b=b+cyclic_length
  end subroutine trim_cyclic
end module g_rotate_grid
!===================================================================

