!============================================================================================
MODULE o_tracers
IMPLICIT NONE
CONTAINS
!=======================================================================
SUBROUTINE tracer_gradient_elements(ttf)
!computes elemental gradient of tracer 

USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
IMPLICIT NONE
real(kind=WP)     :: ttf(nl-1,myDim_nod2D+eDim_nod2D)
integer           :: elem,  elnodes(3)
integer           :: n, nz

DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   DO nz=1, nlevels(elem)-1   
      tr_xy(1,nz, elem)=sum(gradient_sca(1:3,elem)*ttf(nz,elnodes))
      tr_xy(2,nz, elem)=sum(gradient_sca(4:6,elem)*ttf(nz,elnodes))
   END DO
 END DO     
END SUBROUTINE tracer_gradient_elements
!========================================================================================
SUBROUTINE init_tracers_AB(tr_num)
use g_parsup
use o_PARAM, only: tracer_adv
use o_arrays
use o_mesh
use g_comm_auto
IMPLICIT NONE
integer :: tr_num,n,nz 


!filling work arrays
del_ttf=0d0
!AB interpolation
tr_arr_old(:,:,tr_num)=-(0.5+epsilon)*tr_arr_old(:,:,tr_num)+(1.5+epsilon)*tr_arr(:,:,tr_num)
call tracer_gradient_elements(tr_arr(:,:,tr_num))
call exchange_elem_begin(tr_xy)

call tracer_gradient_z(tr_arr(:,:,tr_num))

call exchange_elem_end()! tr_xy used in fill_up_dn_grad
call exchange_nod_begin(tr_z) ! not used in fill_up_dn_grad 

call fill_up_dn_grad

call exchange_nod_end() ! tr_z halos should have arrived by now.

END SUBROUTINE init_tracers_AB
!========================================================================================
SUBROUTINE relax_to_clim(tr_num)
use o_mesh
use g_config,only: dt
USE g_PARSUP
use o_arrays
use o_PARAM, only: tracer_adv
IMPLICIT NONE
integer :: tr_num,n,nz

if ((clim_relax>1.0e-8).and.(tr_num==1)) then
   DO n=1, myDim_nod2D
      tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+&
            relax2clim(n)*dt*(Tclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
  END DO
END if
if ((clim_relax>1.0e-8).and.(tr_num==2)) then
   DO n=1, myDim_nod2D
      tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+&
            relax2clim(n)*dt*(Sclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
   END DO
END IF 
END SUBROUTINE relax_to_clim
END MODULE o_tracers
!========================================================================================
SUBROUTINE tracer_gradient_z(ttf)
!computes vertical gradient of tracer
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
real(kind=WP)     :: ttf(nl-1,myDim_nod2D+eDim_nod2D)
real(kind=WP)     :: dz
integer           :: n, nz, nlev

if (use_ALE) then
DO n=1, myDim_nod2D+eDim_nod2D
   nlev=nlevels_nod2D(n)
   DO nz=2,  nlev-1
      dz=0.5_WP*(hnode_new(nz-1,n)+hnode_new(nz,n))
      tr_z(nz, n)=(ttf(nz-1,n)-ttf(nz,n))/dz
   END DO
   tr_z(1,    n)=0.0_WP
   tr_z(nlev, n)=0.0_WP
END DO
else
DO n=1, myDim_nod2D+eDim_nod2D
     nlev=nlevels_nod2D(n)
     DO nz=2, nlev-1
        dz=0.5_WP*(Z(nz-1)-Z(nz))
        tr_z(nz, n)=(ttf(nz-1,n)-ttf(nz,n))/dz
     END DO
     tr_z(1,    n)=0.0_WP
     tr_z(nlev, n)=0.0_WP
END DO
end if

END SUBROUTINE tracer_gradient_z
