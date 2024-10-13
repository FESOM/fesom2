MODULE Toy_Channel_Nemo
  use mod_mesh
  USE o_ARRAYS
  USE o_PARAM
  USE MOD_PARSUP
  USE MOD_PARTIT
  USE g_config
  !USE i_therm_param, only: albw

  implicit none
  SAVE 
  private
  public :: initial_state_nemo
  real(kind=WP)     ::   ysize       =3180000.0            ! the meridional lenght of the channel [m]
  real(kind=WP)     ::   xsize       =2120000.0   ! the zonal lenght of the channel [m] !4.5*pi*r_earth=90018410.49779853
  real(kind=WP)     ::   zsize       =4000.0      ! m  The depth
  real(kind=WP)     ::   lat0        =17.5  
  
!-----------------------------------------------------------------------------------------
!
  contains
!
!--------------------------------------------------------------------------------------------

subroutine initial_state_nemo(mesh)
 ! Profiles NEMO 2010 
  implicit none
  type(t_mesh), intent(in) , target  :: mesh
  integer                            :: elem, n, nz, elnodes(3),m,k,nzmax
  real(kind=8)                       :: dst!, yn, Fy, Lx
  real(kind=WP)                      :: lon, lat
  real(kind=WP)                      :: ztau, zTstar, ztrp ! wind intensity, atmospheric intensity, heat flux intensity
  REAL(kind=WP)                      :: zrhoa ! Air density kg/m3
  REAL(kind=WP)                      :: a_constant, b_constant ! constants
  REAL(kind=WP)                      :: t_star_top, t_star_bottom
  integer                            :: wind, tpert, tprofile, sprofile, tflux, elevation
  real(kind=WP):: c, c2, c3, c4, c5
  real(kind=WP):: v1, v2, sc1, sc2
  real(kind=WP):: swsurf, aux

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

 ! default values
 wind=1
 tpert=0
 tprofile=1
 sprofile=0
 tflux=1
 elevation=0
 stress_surf   = 0.0_8 !make sure how it works (y-component is 0)
 heat_flux     = 0.0_8

 t_star        = 0.0_8
 qsr_c       = 0.0_8
 relax2clim    = 0.0_8
 tr_arr(:,:,2) = 35.0_8
 !tr_arr(:,:,1) = 20.0_8
 
 Tsurf         = tr_arr(1,:,1)
 Ssurf         = tr_arr(1,:,2)
 water_flux    = 0.0_8
 !clim_relax   = 0.0_8
 
    

! ========
!  wind forcing (momentum flux):
! ========
  ! take the mean intensity as "spring" (and we don't need rotation coefficient as we don't have i,j coordinates)
  ztau = -0.12  !0.105
  zrhoa  = 1.22         
  
  
  if (wind==1) then
      DO elem=1, myDim_elem2D
         elnodes=elem2d_nodes(:,elem)
         lat=sum(coord_nod2D(2,elnodes))/3.0_WP
         lon=sum(coord_nod2D(1,elnodes))/3.0_WP

         stress_surf(1,elem)= ztau * SIN( pi * (lat/rad - 15.0) / (29.0-15.0) )
         
         ! write to the file (in the first time step both x and y components)
     
       END DO
  else 
  stress_surf   = 0.0
  endif 
  
  
! ========
!  elevation
! ========  
  if (elevation==1) then 
      eta_n=0.01*(coord_nod2D(2,:)-30.0*rad)/(15.0*rad)
  end if
  
! ========
!  Initial temperature profile (use Philander analytic profile)
! ======== 
  !IF (tprofile==1) THEN
  !    DO nz=1, nlevels_nod2D(n)-1
  !        tr_arr(nz,:,1) = (1+Z(nz)/4000)*(coord_nod2D(2,:)/rad)
  !    END DO
  !ENDIF
  if (tprofile==1) then  
     DO n=1, myDim_nod2D+eDim_nod2D
         !lat = coord_nod2D(2,n)/rad
         DO nz=1, nlevels_nod2D(n)-1
             
             !tr_arr(nz,n,1) = (1+Z(nz)/4000)*(coord_nod2D(2,n)/rad)
             !tr_arr(nz,n,1)=tr_arr(nz,n,1)-0.95_WP*20*tanh(abs(Z(nz))/300)-abs(Z(nz))/2400.0_WP
             tr_arr(nz,n,1) = 7.5 * ( 1.-TANH((ABS(Z(nz))-80.)/30.) )   &
                &               + 10. * ( 5000. - ABS(Z(nz)) ) /5000.   
             !tr_arr(nz,n,1)=25.0 + .05*tanh((Z(nz) - 300.0)/5.0) - Z(nz)*0.008_WP
         
         END DO
     END DO
  ENDIF
  !write(*,*) mype, 'Temperature', maxval(tr_arr(:,:,1)), minval(tr_arr(:,:,1))
  
  
  Tsurf=tr_arr(1,:,1)


  if (sprofile==1) then  
   DO n=1, myDim_nod2D+eDim_nod2D

       DO nz=1, nlevels_nod2D(n)-1
             
           tr_arr(nz,n,2) = 36. - ABS(Z(nz)) /4000.

         
       END DO
   END DO
  ENDIF
  !write(*,*) mype, 'Salinity', maxval(tr_arr(:,:,2)), minval(tr_arr(:,:,2))
  
  !Ssurf=tr_arr(1,:,2)
  
! ========
!  temperature flux from atmosphere
! ======== 
  zTstar    = 28.3      ! intensity from 28.3 a -5 deg
  ztrp = 4.0        ! retroaction term on heat fluxes (W/m2/K)
  
  IF (tflux==1) then 
      DO n=1, myDim_nod2D+eDim_nod2D 
          lat = coord_nod2D(2,n)
          !the function varies from 27C at 50N to 7C at 15N, again the analogy with spring
          t_star (n) = zTstar * cos(pi*((lat / rad - 5.0 ) / 107.0))
          !solar radiation component
          qsr_c(n) = 230 * COS( pi * ((lat / rad) / 162.0 ))
         
          heat_flux(n) = ztrp * (tr_arr(1,n,1) - t_star(n)) !+ qsr_c(n)
          
          !write(*,*) mype, 't_star', t_star
          
          ! shortwave rad.
          !swsurf=(1.0_WP-albw)*qsr_c(n)
          ! the visible part (300nm-750nm)
          !swsurf=swsurf*0.54_WP
          ! subtract visible sw rad. from heat_flux, which is '+' for upward
          !heat_flux(n)=heat_flux(n)+swsurf
          
          ! limit chl from below
          !if (chl(n) < 0.02_WP) chl(n)=0.02_WP
     
     
          !c=log10(chl(n))
          !c2=c*c
          !c3=c2*c
          !c4=c3*c
          !c5=c4*c
          ! --> coefficients come from Sweeney et al. 2005, "Impacts of shortwave
          ! penetration depthon large scale ocean circulation and heat transport" see 
          ! Appendix A
     

     
          !v1=0.008_WP*c+0.132_WP*c2+0.038_WP*c3-0.017_WP*c4-0.007_WP*c5
          !v2=0.679_WP-v1
          !v1=0.321_WP+v1
          !sc1=1.54_WP-0.197_WP*c+0.166_WP*c2-0.252_WP*c3-0.055_WP*c4+0.042_WP*c5
          !sc2=7.925_WP-6.644_WP*c+3.662_WP*c2-1.815_WP*c3-0.218_WP*c4+0.502_WP*c5
          
               ! convert from heat flux [W/m2] to temperature flux [K m/s]
          !swsurf=swsurf/vcpw
          ! vis. sw. rad. in the colume
          !nzmax=(nlevels(n))
          !sw_3d(1, n)=swsurf
          !do k=2, nzmax
          !   aux=(v1*exp(zbar_3d_n(k,n)/sc1)+v2*exp(zbar_3d_n(k,n)/sc2))
          !   sw_3d(k, n)=swsurf*aux
          !   if (aux < 1.e-5_WP .OR. k==nzmax) then 
          !      sw_3d(k, n)=0.0_WP
          !     exit
          !   end if
          !end do
          
          
      END DO
  ENDIF
  

 
  !if (mype==0) write(*,*) "heat_flux_channel", heat_flux(100)
  !if (mype==0) write(*,*) "radiation_channel", qsr_c(100)


 ! =========
 ! Coriolis calculation
 ! =========

  DO n=1,myDim_elem2D
      elnodes=elem2D_nodes(:,n)
      !dst - change in latitude
      dst=(sum(coord_nod2D(2, elnodes))/3.0-lat0*rad)*r_earth-ysize/2
      coriolis(n)=1.0e-4+dst*1.8e-11	 
  END DO
  write(*,*) mype, 'COR', maxval(coriolis*10000.0), minval(coriolis*10000.0) 
  
  DO n=1,myDim_nod2D+eDim_nod2D
      dst=(coord_nod2D(2, n)-lat0*rad)*r_earth-ysize/2
      coriolis_node(n)=1.0e-4+dst*1.8e-11	 
  END DO
  write(*,*) mype, 'COR_n', maxval(coriolis_node*10000.0), minval(coriolis_node*10000.0) 

 
 
 ! =========
 ! Temperature perturbation to check
 ! =========
  
if (tpert==1) then
  do n=1, myDim_nod2D+eDim_nod2D
     lat=coord_nod2D(2,n)-15.0*rad
     lon=coord_nod2D(1,n)
     if (lon>cyclic_length/2) lon=lon-cyclic_length
     if (lon<-cyclic_length/2) lon=lon+cyclic_length
     dst=sqrt((lat)**2+(lon)**2)
     if (dst>1.5*rad) cycle
     do nz=1, nlevels_nod2D(n)-1 
        !if(abs(Z(nz)+500)<300) then
        tr_arr(nz,n,1)=tr_arr(nz,n,1)+1.0*cos(pi*dst/2.0/1.5/rad)       !exp(-(dst/(1.5*rad))**2)
	!end if
     end do
  end do 
end if
!write(*,*) mype, 'dst', dst
!write(*,*) mype, 'rad', rad

end subroutine initial_state_nemo
END MODULE Toy_Channel_Nemo
