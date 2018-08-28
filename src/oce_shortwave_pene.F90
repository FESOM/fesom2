subroutine cal_shortwave_rad
  ! This routine is inherited from FESOM 1.4 and adopted appropreately. It calculates 
  ! shortwave penetration into the ocean assuming the constant chlorophyll concentration.
  ! No penetration under the ice is applied. A decent way for ice region is to be discussed.
  ! This routine should be called after ice2oce coupling done if ice model is used.
  ! Ref.: Morel and Antoine 1994, Sweeney et al. 2005
  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS
  USE g_PARSUP
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto
  use i_param
  use i_arrays
  use i_therm_param
  IMPLICIT NONE

  integer      :: m, n2, n3, k, nzmax
  real(kind=WP):: swsurf, aux
  real(kind=WP):: c, c2, c3, c4, c5
  real(kind=WP):: v1, v2, sc1, sc2

  sw_3d=0.0

  do n2=1, myDim_nod2D+eDim_nod2D     
     if (use_ice .and. a_ice(n2)> 0._WP) cycle !assume in ice region no penetration
     ! shortwave rad.
     swsurf=(1.0-albw)*shortwave(n2)
     ! the visible part (300nm-750nm)
     swsurf=swsurf*0.54
     ! subtract visible sw rad. from heat_flux, which is '+' for upward
     heat_flux(n2)=heat_flux(n2)+swsurf
 
     ! attenuation func. for vis. sw rad. according to Morel/Antoine param.
     ! the four parameters in the func.

     ! limit chl from below
     if (chl(n2) < 0.02) chl(n2)=0.02
     c=log10(chl(n2))
     c2=c*c
     c3=c2*c
     c4=c3*c
     c5=c4*c
     v1=0.008*c+0.132*c2+0.038*c3-0.017*c4-0.007*c5
     v2=0.679-v1
     v1=v1+0.321
     sc1=1.54-0.197*c+0.166*c2-0.252*c3-0.055*c4+0.042*c5
     sc2=7.925-6.644*c+3.662*c2-1.815*c3-0.218*c4+0.502*c5

     ! convert from heat flux [W/m2] to temperature flux [K m/s]
     swsurf=swsurf/vcpw
     ! vis. sw. rad. in the colume
     nzmax=(nlevels(n2))
     sw_3d(1, n2)=swsurf
     do k=2, nzmax
        aux=(v1*exp(zbar_3d_n(k,n2)/sc1)+v2*exp(zbar_3d_n(k,n2)/sc2))
        sw_3d(k, n2)=swsurf*aux
        if (aux < 1.e-5 .OR. k==nzmax) then 
           sw_3d(k, n2)=0._wp
           exit
        end if
     end do
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
end subroutine cal_shortwave_rad
