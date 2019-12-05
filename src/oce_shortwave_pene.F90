subroutine cal_shortwave_rad(mesh)
  ! This routine is inherited from FESOM 1.4 and adopted appropreately. It calculates 
  ! shortwave penetration into the ocean assuming the constant chlorophyll concentration.
  ! No penetration under the ice is applied. A decent way for ice region is to be discussed.
  ! This routine should be called after ice2oce coupling done if ice model is used.
  ! Ref.: Morel and Antoine 1994, Sweeney et al. 2005
  USE MOD_MESH
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
  type(t_mesh), intent(in) :: mesh
  associate(nod2D=>mesh%nod2D, elem2D=>mesh%elem2D, edge2D=>mesh%edge2D, elem2D_nodes=>mesh%elem2D_nodes, elem_neighbors=>mesh%elem_neighbors, nod_in_elem2D_num=>mesh%nod_in_elem2D_num, &
            nod_in_elem2D=>mesh%nod_in_elem2D, elem_area=>mesh%elem_area, depth=>mesh%depth, nl=>mesh%nl, zbar=>mesh%zbar, z=>mesh%z, nlevels_nod2D=>mesh%nlevels_nod2D, elem_cos=>mesh%elem_cos, &
            coord_nod2D=>mesh%coord_nod2D, geo_coord_nod2D=>mesh%geo_coord_nod2D, metric_factor=>mesh%metric_factor, edges=>mesh%edges, edge_dxdy=>mesh%edge_dxdy, edge_tri=>mesh%edge_tri, &
            edge_cross_dxdy=>mesh%edge_cross_dxdy, gradient_sca=>mesh%gradient_sca, gradient_vec=>mesh%gradient_vec, elem_edges=>mesh%elem_edges, bc_index_nod2D=>mesh%bc_index_nod2D, &
            edge2D_in=>mesh%edge2D_in, area=>mesh%area, ssh_stiff=>mesh%ssh_stiff, nlevels=>mesh%nlevels)

  sw_3d=0.0_WP

  do n2=1, myDim_nod2D+eDim_nod2D     
     if (use_ice .and. a_ice(n2)> 0._WP) cycle !assume in ice region no penetration
     ! shortwave rad.
     swsurf=(1.0_WP-albw)*shortwave(n2)
     ! the visible part (300nm-750nm)
     swsurf=swsurf*0.54_WP
     ! subtract visible sw rad. from heat_flux, which is '+' for upward
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
     nzmax=(nlevels(n2))
     sw_3d(1, n2)=swsurf
     do k=2, nzmax
        aux=(v1*exp(zbar_3d_n(k,n2)/sc1)+v2*exp(zbar_3d_n(k,n2)/sc2))
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
  end associate
end subroutine cal_shortwave_rad
