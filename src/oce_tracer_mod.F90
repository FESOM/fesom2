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
call tracer_gradient_elements(tr_arr_old(:,:,tr_num))
call exchange_elem(tr_xy)
call fill_up_dn_grad
END SUBROUTINE init_tracers_AB
!========================================================================================
SUBROUTINE adv_tracers(tr_num)
use g_parsup
use o_mesh,only: nl,nod2D
use o_PARAM, only: tracer_adv
use o_arrays
IMPLICIT NONE
integer :: tr_num 

select case (tracer_adv)
   case(1) !MUSCL
      call adv_tracer_muscl(tr_arr(:,:,tr_num), del_ttf, tr_arr_old(:,:,tr_num))
   case(2) !MUSCL+FCT(3D)
      call adv_tracer_fct(tr_arr(:,:,tr_num),del_ttf,tr_arr_old(:,:,tr_num), 0.75_WP)
   CASE DEFAULT !unknown
      IF (mype==0) write(*,*) 'Unknown advection type. Check your namelists.'
      call par_ex(1)
END SELECT
END SUBROUTINE adv_tracers
!========================================================================================
SUBROUTINE diff_tracers(tr_num)
use o_mesh
USE g_PARSUP
use o_arrays
IMPLICIT NONE
integer, intent(in) :: tr_num
integer             :: n, nl1

tr_arr_old(:,:,tr_num)=tr_arr(:,:,tr_num) !DS: check that this is the right place!
call diff_part_hor
if (not(i_vert_diff)) call diff_ver_part_expl(tr_num)
!Update tracer berofe implicit operator is applied
DO n=1, myDim_nod2D
   nl1=nlevels_nod2D(n)-1
   tr_arr(1:nl1,n,tr_num)=tr_arr(1:nl1,n,tr_num)+del_ttf(1:nl1,n)
END DO
!We DO not set del_ttf to zero because it will not be used in this timestep anymore
!init_tracers will set it to zero for the next timestep
if (i_vert_diff) call diff_ver_part_impl(tr_num)
END SUBROUTINE diff_tracers
!========================================================================================
SUBROUTINE diff_part_hor
use o_ARRAYS
use g_PARSUP
use o_MESH
USE o_param
use g_config
IMPLICIT NONE
real(kind=WP)   :: deltaX1,deltaY1,deltaX2,deltaY2
integer         :: edge
integer         :: nl1,nz,el(2),elnodes(3),n,nl2,enodes(2)
real(kind=WP)   :: temp(3),d,c1,Tx,Ty,Kh

ttrhs=0d0
DO edge=1, myDim_edge2D
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   el=edge_tri(:,edge)
   enodes=edges(:,edge)
   nl1=nlevels(el(1))-1
   elnodes=elem2d_nodes(:,el(1))
   Kh=elem_area(el(1))
   IF (el(2)>0) then 
      Kh=0.5*(Kh+elem_area(el(2)))
      nl2=nlevels(el(2))-1
   ENDif
    Kh=K_hor*Kh/scale_area

   DO nz=1,nl1
      Tx=0d0
      Ty=0d0
      IF ((el(2)>0).and.(nz<=nl2)) then
         Tx=Tx+0.5*Kh*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
         Ty=Ty+0.5*Kh*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
      ELSE
         Tx=Tx+Kh*tr_xy(1,nz,el(1))
         Ty=Ty+Kh*tr_xy(2,nz,el(1))
      END IF  
      ttrhs(nz,enodes(1)) = ttrhs(nz,enodes(1)) - Ty*deltaX1 + Tx*deltaY1
      ttrhs(nz,enodes(2)) = ttrhs(nz,enodes(2)) + Ty*deltaX1 - Tx*deltaY1
   ENDDO

   IF (el(2)>0) then
      deltaX2=edge_cross_dxdy(3,edge)
      deltaY2=edge_cross_dxdy(4,edge)
      DO nz=1,nl2
         Tx=0d0
         Ty=0d0

      IF (nz<=nl1) then
         Tx=Tx+0.5*Kh*(tr_xy(1,nz,el(1))+tr_xy(1,nz,el(2)))
         Ty=Ty+0.5*Kh*(tr_xy(2,nz,el(1))+tr_xy(2,nz,el(2)))
      ELSE
         Tx=Tx+Kh*tr_xy(1,nz,el(2))
         Ty=Ty+Kh*tr_xy(2,nz,el(2))
      END if

      ttrhs(nz,enodes(1)) = ttrhs(nz,enodes(1)) + Ty*deltaX2 - Tx*deltaY2
      ttrhs(nz,enodes(2)) = ttrhs(nz,enodes(2)) - Ty*deltaX2 + Tx*deltaY2
    ENDDO
   ENDif
ENDDO

DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      del_ttf(nz,n)=del_ttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
   END DO
END DO
END SUBROUTINE diff_part_hor
!========================================================================================
SUBROUTINE diff_ver_part_expl(tr_num)
!Vertical diffusive flux(explicit scheme):                                                                            
use o_ARRAYS
use o_MESH
use g_PARSUP
use g_config,only: dt
IMPLICIT NONE

real(kind=WP) :: vd_flux(nl-1)
real(kind=WP) :: rdata,flux,rlx
integer       :: nz,nl1,tr_num,n
real(kind=WP) :: zinv1,Ty

DO n=1, myDim_nod2D
   nl1=nlevels_nod2D(n)-1
   vd_flux=0d0
   IF (tr_num==1) then
      flux  = -heat_flux(n)/density_0/4200.0
      rdata =  Tsurf(n)
      rlx   =  surf_relax_T
   ELSEIF (tr_num==2) then
      flux  =  water_flux(n)*tr_arr(1,n,2)
      rdata =  Ssurf(n)
      rlx   =  surf_relax_S
   ELSE
      flux  = 0d0
      rdata = 0d0
      rlx=0d0
   ENDif

!Surface forcing
   vd_flux(1)= flux + rlx*(rdata-tr_arr(1,n,tr_num))

   DO nz=2,nl1
      zinv1=1.0_WP/(Z(nz-1)-Z(nz))
      Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
              Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2

      vd_flux(nz) = (Kv(nz,n)+Ty)*(tr_arr(nz-1,n,tr_num)-tr_arr(nz,n,tr_num))*zinv1*area(nz,n)
ENDDO
DO nz=1,nl1-1
 del_ttf(nz,n) =del_ttf(nz,n) + (vd_flux(nz) - vd_flux(nz+1))/(zbar(nz)-zbar(nz+1))*dt/area(nz,n)
ENDDO
 del_ttf(nl1,n) = del_ttf(nl1,n) + (vd_flux(nl1)/(zbar(nl1)-zbar(nl1+1)))*dt/area(nl1,n)
ENDDO

END SUBROUTINE diff_ver_part_expl
!========================================================================================
SUBROUTINE diff_ver_part_impl(tr_num)
USE o_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer, intent(in)   :: tr_num
real(kind=WP)         :: a(nl), b(nl), c(nl), tr(nl)
real(kind=WP)         ::  cp(nl), tp(nl)
integer               ::  nz, n, nzmax
real(kind=WP)         ::  m, zinv, dt_inv
real(kind=WP)         ::  rsss, Ty,Ty1,c1,zinv1,zinv2, v_adv

dt_inv=1.0_WP/dt
DO n=1, myDim_nod2D
   nzmax=nlevels_nod2D(n)
! the first row
   nz=1
   zinv2=1.0_WP/(Z(nz)-Z(nz+1))
   zinv=1.0_WP*dt/(zbar(1)-zbar(2))
   Ty1= Kd(4,nz,n)*(Z(nz)-zbar(nz+1))*zinv2 *neutral_slope(3,nz,n)**2 + \
        Kd(4,nz+1,n)*(zbar(nz+1)-Z(nz+1))*zinv2 *neutral_slope(3,nz+1,n)**2
   c(1)=-(Kv(2,n)+Ty1)*zinv2*zinv*area(2,n)/area(1,n)
   a(1)=0.0_WP
   b(1)=-c(1)+1.0_WP
! update from the vertical advection
   v_adv=zinv*area(2,n)/area(1,n)
   b(1)=b(1)+Wvel_i(1, n)*zinv-min(0._WP, Wvel_i(2, n))*v_adv
   c(1)=c(1)-max(0._WP, Wvel_i(2, n))*v_adv

   zinv1=zinv2
! regular part of coefficients:
   DO nz=2, nzmax-2
      zinv2=1.0_WP/(Z(nz)-Z(nz+1))
      Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
          Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2
      Ty1= Kd(4,nz,n)*(Z(nz)-zbar(nz+1))*zinv2 *neutral_slope(3,nz,n)**2 + \
           Kd(4,nz+1,n)*(zbar(nz+1)-Z(nz+1))*zinv2 *neutral_slope(3,nz+1,n)**2
        
      zinv=1.0_WP*dt/(zbar(nz)-zbar(nz+1))
      a(nz)=-(Kv(nz,n)+Ty)*zinv1*zinv
      c(nz)=-(Kv(nz+1,n)+Ty1)*zinv2*zinv*area(nz+1,n)/area(nz,n)
      b(nz)=-a(nz)-c(nz)+1.0_WP
      zinv1=zinv2

! update from the vertical advection
      v_adv=zinv
      a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv
      b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv

      v_adv=v_adv*area(nz+1,n)/area(nz,n)
      b(nz)=b(nz)-min(0._WP, Wvel_i(nz+1, n))*v_adv
      c(nz)=c(nz)-max(0._WP, Wvel_i(nz+1, n))*v_adv
   END DO
! the last row
   nz=nzmax-1
   zinv=1.0_WP*dt/(zbar(nzmax-1)-zbar(nzmax))
   Ty= Kd(4,nz-1,n)*(Z(nz-1)-zbar(nz))*zinv1 *neutral_slope(3,nz-1,n)**2 + \
       Kd(4,nz,n)*(zbar(nz)-Z(nz))*zinv1 *neutral_slope(3,nz,n)**2
   a(nzmax-1)=-(Kv(nzmax-1,n)+Ty)*zinv1*zinv
   b(nzmax-1)=-a(nzmax-1)+1.0_WP
   c(nzmax-1)=0.0_WP

! update from the vertical advection
   v_adv=zinv
   a(nz)=a(nz)+min(0._WP, Wvel_i(nz, n))*v_adv       
   b(nz)=b(nz)+max(0._WP, Wvel_i(nz, n))*v_adv

! the rhs:
   DO nz=2,nzmax-2
      tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-1.0_WP)*tr_arr(nz,n,tr_num)-c(nz)*tr_arr(nz+1,n,tr_num)
   END DO
   nz=nzmax-1
   tr(nz)=-a(nz)*tr_arr(nz-1,n,tr_num)-(b(nz)-1.0_WP)*tr_arr(nz,n,tr_num)
! the first row CONTAINS also surface forcing
   tr(1)=-(b(1)-1.0_WP)*tr_arr(1,n,tr_num)-c(1)*tr_arr(2,n,tr_num)

   zinv=1.0_WP*dt/(zbar(1)-zbar(2))
   IF (tr_num==1) then
      tr(1)= tr(1)  -  &
             zinv*(1d-3*heat_flux(n)/(4.2_WP)/density_0 - surf_relax_T*(Tsurf(n)-tr_arr(1,n,1)))
   ELSEIF (tr_num==2) then
      rsss=ref_sss
      IF (ref_sss_local) rsss = tr_arr(1,n,2)
      tr(1)= tr(1)  +  &
             zinv*(rsss*water_flux(n) + surf_relax_S*(Ssurf(n)-tr_arr(1,n,2)))
   END IF
! the sweep algorithm
! initialize c-prime and s,t-prime
      cp(1) = c(1)/b(1)
      tp(1) = tr(1)/b(1)
! solve for vectors c-prime and t, s-prime
      DO nz = 2,nzmax-1
         m = b(nz)-cp(nz-1)*a(nz)
         cp(nz) = c(nz)/m
         tp(nz) = (tr(nz)-tp(nz-1)*a(nz))/m
      END DO
! initialize x
      tr(nzmax-1) = tp(nzmax-1)
! solve for x from the vectors c-prime and d-prime
      DO nz = nzmax-2, 1, -1
         tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
      END DO

      DO nz=1,nzmax-1
         tr_arr(nz,n,tr_num)=tr_arr(nz,n,tr_num)+tr(nz)
      END DO
        
END DO   !!! cycle over nodes
END SUBROUTINE diff_ver_part_impl
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
      tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+\
            relax2clim(n)*dt*(Tclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
  END DO
END if
if ((clim_relax>1.0e-8).and.(tr_num==2)) then
   DO n=1, myDim_nod2D
      tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)=tr_arr(1:nlevels_nod2D(n)-1,n,tr_num)+\
            relax2clim(n)*dt*(Sclim(1:nlevels_nod2D(n)-1,n)-tr_arr(1:nlevels_nod2D(n)-1,n,tr_num))
   END DO
END IF 
END SUBROUTINE relax_to_clim
END MODULE o_tracers
!========================================================================================
SUBROUTINE solve_tracers
use g_parsup
use o_PARAM, only: tracer_adv,num_tracers
use o_arrays
use o_mesh
use g_comm_auto
use o_tracers
IMPLICIT NONE
integer :: tr_num

!only MUSCL(+FCT) supported

DO tr_num=1,num_tracers
   call init_tracers_AB(tr_num)
   call adv_tracers(tr_num)
   call diff_tracers(tr_num)
   call relax_to_clim(tr_num)
   call exchange_nod(tr_arr(:,:,tr_num))
END DO
END SUBROUTINE solve_tracers
