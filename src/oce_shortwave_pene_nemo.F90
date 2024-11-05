subroutine cal_shortwave_rad_nemo(ice, tracers, partit, mesh)
  ! This routine is inherited from FESOM 1.4 and adopted appropreately. It calculates 
  ! shortwave penetration into the ocean assuming the constant chlorophyll concentration.
  ! No penetration under the ice is applied. A decent way for ice region is to be discussed.
  ! This routine should be called after ice2oce coupling done if ice model is used.
  ! Ref.: Morel and Antoine 1994, Sweeney et al. 2005
  USE MOD_MESH
  USE MOD_ICE
  USE o_PARAM
  USE o_ARRAYS
  USE MOD_PARSUP
  USE MOD_PARTIT
  USE MOD_TRACER
  USE g_CONFIG
  use g_forcing_arrays!, only: chl, sw_3d
  use g_comm_auto
 !  use i_param
 !  use i_arrays
 ! use i_therm_param
 ! use Toy_Channel_Nemo!, only: qsr_c, t_star
  IMPLICIT NONE
  type(t_ice) ,intent(in), target :: ice
  type(t_tracer), intent(inout), target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh), intent(in) , target :: mesh


  integer      :: m, n2, n3, k, nzmax, nzmin
  real(kind=WP):: swsurf, aux, zTstar, ztrp
  real(kind=WP):: c, c2, c3, c4, c5
  real(kind=WP):: v1, v2, sc1, sc2
  real(kind=WP), pointer :: albw


#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  sw_3d=0.0_WP
  albw => ice%thermo%albw  

  
  do n2=1, myDim_nod2D+eDim_nod2D     
     
     !calculate heat flux
     zTstar    = 28.3      ! intensity from 28.3 a -5 deg
     ztrp = 4.0        ! retroaction term on heat fluxes (W/m2/K)
     
     !lat = coord_nod2D(2,n2)
     
     !t_star (n2) = zTstar * cos(pi*((lat / rad - 5.0 ) / 107.0))
     heat_flux(n2) = ztrp * (tracers%data(1)%values(1,n2) - t_star(n2)) 
     
     
     
     ! shortwave rad.
     swsurf=(1.0_WP-albw)*qsr_c(n2)
     ! the visible part (300nm-750nm)
     swsurf=swsurf*0.54_WP
     ! subtract visible sw rad. from heat_flux, which is '+' for upward
     
     !if (mype==0) write(*,*) "heat_flux_routine", heat_flux(100)
     !if (mype==0) write(*,*) "radiation_routine", qsr_c(100)
     
     heat_flux(n2)=heat_flux(n2)+swsurf
 
     ! attenuation func. for vis. sw rad. according to Morel/Antoine param.
     ! the four parameters in the func.
     


     ! limit chl from below
     if (chl(n2) < 0.02_WP) chl(n2)=0.02_WP
     
     
     c=log10(chl(n2))
     c2=c*c
     c3=c2*c
     c4=c3*c
     c5=c4*c
     ! --> coefficients come from Sweeney et al. 2005, "Impacts of shortwave
     ! penetration depthon large scale ocean circulation and heat transport" see 
     ! Appendix A
     

     
     v1=0.008_WP*c+0.132_WP*c2+0.038_WP*c3-0.017_WP*c4-0.007_WP*c5
     v2=0.679_WP-v1
     v1=0.321_WP+v1
     sc1=1.54_WP-0.197_WP*c+0.166_WP*c2-0.252_WP*c3-0.055_WP*c4+0.042_WP*c5
     sc2=7.925_WP-6.644_WP*c+3.662_WP*c2-1.815_WP*c3-0.218_WP*c4+0.502_WP*c5
    

    
     ! convert from heat flux [W/m2] to temperature flux [K m/s]
     swsurf=swsurf/vcpw
     ! vis. sw. rad. in the colume
     nzmax=(mesh%nlevels(n2))
     nzmin=(mesh%ulevels(n2))
     sw_3d(nzmin, n2)=swsurf
     do k=nzmin+1, nzmax
        aux=(v1*exp(mesh%zbar_3d_n(k,n2)/sc1)+v2*exp(mesh%zbar_3d_n(k,n2)/sc2))
        sw_3d(k, n2)=swsurf*aux
        if (aux < 1.e-5_WP .OR. k==nzmax) then 
           sw_3d(k, n2)=0.0_WP
           exit
        end if
     end do
     
     ! sw_3d --> TEMPERATURE FLUX through full depth level interfaces into/out off 
     ! the tracer volume 
     ! sum(sw_3d(1:nlevels(n2)-1,n2)-sw_3d(2:nlevels(n2),n2)) = swsurf !!!
     
!for testing the subroutine
!if (mype==30 .and. n2==100) then
!write(*,*) 'heat_flux=', heat_flux(n2)
!write(*,*) 'short/longwave=', shortwave(n2), longwave(n2), swsurf*vcpw
!do k=1, nzmax
!   write(*,*) sw_3d(k, n2)*vcpw
!end do
!end if

  end do
!call par_ex
!stop
end subroutine cal_shortwave_rad_nemo
