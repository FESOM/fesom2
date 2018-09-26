subroutine mo_convect
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
use i_arrays
IMPLICIT NONE

integer               :: node, elem, nz, elnodes(3)
real(kind=WP)         :: kv_conv=0.1_WP, av_conv=0.1_WP

if (mo_on) then
do node=1, myDim_nod2D+eDim_nod2D
!
! calcualte monin obukhov length
   call mo_length(water_flux(node),heat_flux(node), &         
        stress_atmoce_x(node),stress_atmoce_y(node), &    
        u_ice(node),v_ice(node),a_ice(node), &                             
        dt, mixlength(node))
!_______________________________________________________________
! increase vertical diffusion within monin obukov length to namelist
! parameter modiff. modiff in moment set to 0.01 --> that means 
! very strong vertical mixing within mixlength
   do nz = 2,nlevels_nod2d(node)-1
      mo(nz,node) = 0._WP
      if(abs(zbar_3d_n(nz,node)) <= mixlength(node)) mo(nz,node)=modiff    ! Potentialy bad place 
   end do 
end do
end if
!
!___________________________________________________________________________
DO node=1, myDim_nod2D+eDim_nod2D
   DO nz=2, nlevels_nod2d(node)-1
      if (mo_on) Kv(nz,node)=Kv(nz,node)+mo(nz,node)
      if (bvfreq(nz, node) < 0.) Kv(nz,node)=max(kv_conv, Kv(nz,node))
!     Kv(nz,node)=min(Kv(nz,node), kv_conv)
   END DO
END DO
	
DO elem=1, myDim_elem2D+eDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   DO nz=2,nlevels(elem)-1
      if (mo_on) Av(nz,elem)=Av(nz,elem)+sum(mo(nz,elnodes))/3.0 
      if (any(bvfreq(nz, elnodes) < 0.)) Av(nz,elem)=max(av_conv, Av(nz,elem))
!     Av(nz,elem)=min(Av(nz,elem), av_conv)
   END DO
END DO

end subroutine mo_convect

