! ===================================================================
! Contains routines needed for computations of dynamics.
! includes: update_vel, compute_vel_nodes, viscosity_filt2x
! ===================================================================
SUBROUTINE update_vel(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer       :: elem, elnodes(3), nz, m
real(kind=WP) :: eta(3) 
real(kind=WP) :: Fx, Fy
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 DO elem=1, myDim_elem2D
    elnodes=elem2D_nodes(:,elem)
    eta=-g*theta*dt*d_eta(elnodes)
    Fx=sum(gradient_sca(1:3,elem)*eta)
    Fy=sum(gradient_sca(4:6,elem)*eta)
    DO nz=1, nlevels(elem)-1
     UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
     UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
    END DO
 END DO
eta_n=eta_n+d_eta
call exchange_elem(UV)
end subroutine update_vel
!==========================================================================
subroutine compute_vel_nodes(mesh)
USE MOD_MESH
USE o_PARAM
USE o_ARRAYS
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
integer            :: n, nz, k, elem
real(kind=WP)      :: tx, ty, tvol
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

DO n=1, myDim_nod2D 
   DO nz=1, nlevels_nod2D(n)-1
      tvol=0.0_WP
      tx=0.0_WP
      ty=0.0_WP
      DO k=1, nod_in_elem2D_num(n)
         elem=nod_in_elem2D(k,n)
         if (nlevels(elem)-1<nz) cycle
            tvol=tvol+elem_area(elem)
            tx=tx+UV(1,nz,elem)*elem_area(elem)
            ty=ty+UV(2,nz,elem)*elem_area(elem)
      END DO
      Unode(1,nz,n)=tx/tvol
      Unode(2,nz,n)=ty/tvol
   END DO
END DO
call exchange_nod(Unode)
end subroutine compute_vel_nodes
!===========================================================================
SUBROUTINE viscosity_filt2x(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_config
USE g_comm_auto
IMPLICIT NONE
type(t_mesh), intent(in) , target :: mesh
real(kind=WP) :: u1, v1, tau_inv, s, factor, factor1, factor2
integer       :: ed, el(2), nz, elem
real(kind=WP) :: UV_c(2,mesh%nl-1,myDim_elem2D+eDim_elem2D), UV_f(2,mesh%nl-1,myDim_elem2D+eDim_elem2D)

#include "associate_mesh.h"

! Filter is applied twice. It should be approximately 
! equivalent to biharmonic operator with the coefficient
! (tau_c/day)a^3/9. Scaling inside is found to help 
! with smoothness in places of mesh transition. *(it makes a^3 from a^4) 


UV_c=0.0_WP
!NR UV_f=0.0_WP
tau_inv=dt*tau_c/3600.0_WP/24.0_WP     ! SET IT experimentally 
  
DO ed=1, myDim_edge2D+eDim_edge2D
   if(myList_edge2D(ed)>edge2D_in) cycle
   el=edge_tri(:,ed)

   factor1 = -tau_inv*sqrt(scale_area/elem_area(el(1)))
   factor2 = -tau_inv*sqrt(scale_area/elem_area(el(2)))

   DO nz=1,min(nlevels(el(1)),nlevels(el(2)))-1
      u1 = UV(1,nz,el(1)) - UV(1,nz,el(2))
      v1 = UV(2,nz,el(1)) - UV(2,nz,el(2))

      UV_c(1,nz,el(1)) = UV_c(1,nz,el(1)) - u1*factor1
      UV_c(2,nz,el(1)) = UV_c(2,nz,el(1)) - v1*factor1
      UV_c(1,nz,el(2)) = UV_c(1,nz,el(2)) + u1*factor2
      UV_c(2,nz,el(2)) = UV_c(2,nz,el(2)) + v1*factor2
   END DO 
END DO
 
! ============ 
! Contribution from boundary edges (Dirichlet boundary conditions)
! ============
!NR If it should be used again, please integrate into the loop above.
!NR if (inner edge) then [...] else [Dirichlet] endif
! DO ed=1, myDim_edge2D+eDim_edge2D
!    if(myList_edge2D(ed)<=edge2D_in) cycle
!    el=edge_tri(:, ed)
!    DO  nz=1, nlevels(el(1))-1
!        UV_c(1,nz,el(1))=UV_c(1,nz,el(1))-2.0_WP*UV(1,nz,el(1))
!        UV_c(2,nz,el(1))=UV_c(2,nz,el(1))-2.0_WP*UV(2,nz,el(1))
!    END DO
! END DO

!NR  moved this factor to the edge loop. If the commented Dirichlet
!NR boundary condition shall be reactivated, do not forget the factor.
!NR do elem=1,myDim_elem2D
!NR   factor = tau_inv*sqrt(scale_area/elem_area(elem))
!NR   Do nz=1,nlevels(elem)-1
!NR       UV_c(1,nz,elem) = -UV_c(1,nz,elem)*factor
!NR       UV_c(2,nz,elem) = -UV_c(2,nz,elem)*factor
!NR   END DO
!NR end do

call exchange_elem(UV_c)


DO ed=1, myDim_edge2D+eDim_edge2D
   if(myList_edge2D(ed)>edge2D_in) cycle
   el=edge_tri(:,ed)
   DO nz=1,min(nlevels(el(1)),nlevels(el(2)))-1 
      u1 = UV_c(1,nz,el(1)) - UV_c(1,nz,el(2))
      v1 = UV_c(2,nz,el(1)) - UV_c(2,nz,el(2))

     !NR UV_f(1,nz,el(1)) = UV_f(1,nz,el(1)) - u1
     !NR UV_f(2,nz,el(1)) = UV_f(2,nz,el(1)) - v1
     !NR UV_f(1,nz,el(2)) = UV_f(1,nz,el(2)) + u1
     !NR UV_f(2,nz,el(2)) = UV_f(2,nz,el(2)) + v1
     !NR As UV_f is simply added to UV_rhs, we can skip this intermediate field
     !NR and add to UV_rhs directly.
     !NR However, if the limit to the contribution of the filter is to be applied,
     !NR this does not work any more. At the time of my optimizations, the limit
     !NR was commented.
     !NR If UV_f is to be used again, please remember to initialize it with UV_f=0. 
     UV_rhs(1,nz,el(1)) = UV_rhs(1,nz,el(1)) - u1
     UV_rhs(2,nz,el(1)) = UV_rhs(2,nz,el(1)) - v1
     UV_rhs(1,nz,el(2)) = UV_rhs(1,nz,el(2)) + u1
     UV_rhs(2,nz,el(2)) = UV_rhs(2,nz,el(2)) + v1

   END DO 
END DO
 
!NR DO elem=1, myDim_elem2D
!NR   DO nz=1, nlevels(elem)-1 
!NR !      u1 = sqrt(UV_f(1,nz,elem)**2+UV_f(2,nz,elem)**2)+tiny(u1)
!NR !      v1 = sqrt(UV(1,nz,elem)**2+UV(2,nz,elem)**2)
!NR    ! we limit the maximum contribution from the filter such, that the update is less than the N (N=2 currently) times velocity
!NR    ! this is done to force the CFL, which is otherwise exceeded in some points
!NR    ! some other criteria is welcome (i.e. like computing the eigenvalues from filtering)
!NR       UV_rhs(1,nz,elem) = UV_rhs(1,nz,elem)+UV_f(1,nz,elem)!*min(1.0_WP, 2.0_WP*v1/u1)
!NR       UV_rhs(2,nz,elem) = UV_rhs(2,nz,elem)+UV_f(2,nz,elem)!*min(1.0_WP, 2.0_WP*v1/u1)
!NR    END DO 
!NR END DO
end subroutine viscosity_filt2x
!===========================================================================
subroutine viscosity_filter(option, mesh)
use o_PARAM
use g_PARSUP
use MOD_MESH
IMPLICIT NONE 
integer                  :: option
type(t_mesh), intent(in) , target :: mesh
! Driving routine 
! Background viscosity is selected in terms of Vl, where V is 
! background velocity scale and l is the resolution. V is 0.005 
! or 0.01, perhaps it would be better to pass it as a parameter.

! h_viscosity_leiht needs vorticity, so vorticity array should be 
! allocated. At present, there are two rounds of smoothing in 
! h_viscosity. 

SELECT CASE (option)
CASE (1)
     ! ====
     ! Laplacian+Leith parameterization + harmonic background
     ! ====
     call h_viscosity_leith(mesh)
     call viscosity_filtxx(mesh)
CASE (2)
     ! ===
     ! Laplacian+Leith+biharmonic background
     ! ===
     call h_viscosity_leith(mesh)
     call viscosity_filtxxx(mesh)
CASE (3)
     ! ===
     ! Biharmonic+Leith+ background 
     ! ===
     call h_viscosity_leith(mesh)
     call viscosity_filt2xx(2, mesh)
CASE (4)
     ! ===
     ! Biharmonic+upwind-type+ background 
     ! ===
     call viscosity_filt2xx(1, mesh)
CASE (5)
     call viscosity_filt_h_backscatter(mesh)
CASE DEFAULT
     if (mype==0) write(*,*) 'mixing scheme with option ' , option, 'has not yet been implemented'
     call par_ex
     stop
END SELECT
end subroutine viscosity_filter  
! ===================================================================
SUBROUTINE viscosity_filtxx(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE

real(kind=WP) :: u1, v1, le(2), len, crosslen, vi 
integer       :: nz, ed, el(2)
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! An analog of harmonic viscosity operator.  
 ! It adds to the rhs(0) Visc*(u1+u2+u3-3*u0)/area
 ! on triangles, which is Visc*Laplacian/4 on equilateral triangles. 
 ! The contribution from boundary edges is neglected (free slip). 
 
 
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    le=edge_dxdy(:,ed)
    le(1)=le(1)*sum(elem_cos(el))*0.25_WP
    len=sqrt(le(1)**2+le(2)**2)*r_earth
    le(1)=edge_cross_dxdy(1,ed)-edge_cross_dxdy(3,ed)
    le(2)=edge_cross_dxdy(2,ed)-edge_cross_dxdy(4,ed)
    crosslen=sqrt(le(1)**2+le(2)**2) 
    DO  nz=1,minval(nlevels(el))-1
     vi=dt*len*(Visc(nz,el(1))+Visc(nz,el(2)))/crosslen
     vi=max(vi,0.005_WP*len*dt) ! This helps to reduce noise in places where 
                              ! Visc is small and decoupling might happen 
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))*vi
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))*vi
     
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO
 END DO
end subroutine viscosity_filtxx
! ===================================================================
SUBROUTINE viscosity_filt2xx(option, mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
! An energy conserving version
! Also, we use the Leith viscosity
!

real(kind=WP) :: u1, v1, vi, len
integer       :: ed, el(2), nz, option
real(kind=WP), allocatable  :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! Filter is applied twice. 
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_WP
 V_c=0.0_WP
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     U_c(nz,el(1))=U_c(nz,el(1))-u1
     U_c(nz,el(2))=U_c(nz,el(2))+u1
     V_c(nz,el(1))=V_c(nz,el(1))-v1
     V_c(nz,el(2))=V_c(nz,el(2))+v1
    END DO 
 END DO
 if(option==1) then
 Do ed=1,myDim_elem2D
    len=sqrt(elem_area(ed))                     
    len=dt*len/30.0_WP
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coefficient because of 
     ! the division by area in the end 
     ! ====
     ! Case 1 -- an analog to the third-order upwind (vi=|u|l/12)
     ! ====
     vi=max(0.2_WP,sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2))*len 
     U_c(nz,ed)=-U_c(nz,ed)*vi                             
     V_c(nz,ed)=-V_c(nz,ed)*vi
    END DO
 end do
 end if
if(option==2) then
 Do ed=1,myDim_elem2D
    len=sqrt(elem_area(ed))                     
    len=dt*len/30.0_WP
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coefficient because of 
     ! the division by area in the end 
     ! ===   
     ! Case 2 -- Leith +background (call h_viscosity_leith)
     ! ===
     vi=max(Visc(nz,ed),0.15_WP*len)
     !    
     U_c(nz,ed)=-U_c(nz,ed)*vi                             
     V_c(nz,ed)=-V_c(nz,ed)*vi
    END DO
 end do
 end if

 call exchange_elem(U_c)
 call exchange_elem(V_c)
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
     el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
     u1=(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=(V_c(nz,el(1))-V_c(nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
     
 deallocate(V_c,U_c)
end subroutine viscosity_filt2xx
! ===================================================================
SUBROUTINE viscosity_filtxxx(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
! An energy and momentum conserving version.
! We use the harmonic Leith viscosity + 
! biharmonic background viscosity
!

real(kind=WP) :: u1, v1, vi, len, crosslen, le(2)
integer       :: ed, el(2), nz
real(kind=WP), allocatable  :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! Filter is applied twice. 
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_WP
 V_c=0.0_WP
  DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    le=edge_dxdy(:,ed)
    le(1)=le(1)*sum(elem_cos(el))*0.25_WP
    len=sqrt(le(1)**2+le(2)**2)*r_earth
    le(1)=edge_cross_dxdy(1,ed)-edge_cross_dxdy(3,ed)
    le(2)=edge_cross_dxdy(2,ed)-edge_cross_dxdy(4,ed)
    crosslen=sqrt(le(1)**2+le(2)**2) 
    DO  nz=1,minval(nlevels(el))-1
     vi=dt*len*(Visc(nz,el(1))+Visc(nz,el(2)))/crosslen
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     U_c(nz,el(1))=U_c(nz,el(1))-u1
     U_c(nz,el(2))=U_c(nz,el(2))+u1
     V_c(nz,el(1))=V_c(nz,el(1))-v1
     V_c(nz,el(2))=V_c(nz,el(2))+v1
     u1=u1*vi
     v1=v1*vi
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
 
 Do ed=1,myDim_elem2D
    len=0.01*sqrt(elem_area(ed))                     
    Do nz=1,nlevels(ed)-1
     vi=dt*abs(min(Visc(nz,ed)-len,0.0_WP))    
     U_c(nz,ed)=-U_c(nz,ed)*vi                             
     V_c(nz,ed)=-V_c(nz,ed)*vi
    END DO
 end do
 call exchange_elem(U_c)
 call exchange_elem(V_c)
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
     el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
     u1=(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=(V_c(nz,el(1))-V_c(nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
     
 deallocate(V_c,U_c)
end subroutine viscosity_filtxxx

! ===================================================================

SUBROUTINE h_viscosity_leith(mesh)
!
! Coefficient of horizontal viscosity is a combination of 
! the Leith and modified Leith (with Div_c)
! and background (with A_hor)
! Background is scaled as sqrt(A/A_0), others are scaled
! in natural way by construction
! Can only be used with the  momentum invariant 
! form
!
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
real(kind=WP)  ::  dz, div_elem(3), xe, ye, vi
integer        :: elem, nl1, nz, elnodes(3),n,k, nt
real(kind=WP)  :: leithx, leithy
real(kind=WP), allocatable :: aux(:,:) 
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

	!  
	if(mom_adv<4) call relative_vorticity(mesh)  !!! vorticity array should be allocated
	! Fill in viscosity:
	DO  elem=1, myDim_elem2D    	!! m=1, myDim_elem2D
									!! elem=myList_elem2D(m)
		!_______________________________________________________________________
		! Here can not exchange zbar_n & Z_n with zbar_3d_n & Z_3d_n because 
		! they run over elements here 
		nl1 =nlevels(elem)-1
		zbar_n=0.0_WP
		! in case of partial cells zbar_n(nzmax) is not any more at zbar(nzmax), 
		! zbar_n(nzmax) is now zbar_e_bot(elem), 
		zbar_n(nl1+1)=zbar_e_bot(elem)
		do nz=nl1,2,-1
			zbar_n(nz) = zbar_n(nz+1) + helem(nz,elem)
		end do
		zbar_n(1) = zbar_n(2) + helem(1,elem)
		
		!_______________________________________________________________________
		elnodes=elem2D_nodes(:,elem)
		do nz=1,nl1
			dz=zbar_n(nz)-zbar_n(nz+1)
			div_elem=(Wvel(nz,elnodes)-Wvel(nz+1,elnodes))/dz
			xe=sum(gradient_sca(1:3,elem)*div_elem)
			ye=sum(gradient_sca(4:6,elem)*div_elem)
			div_elem=vorticity(nz,elnodes)
			leithx=sum(gradient_sca(1:3,elem)*div_elem)
			leithy=sum(gradient_sca(4:6,elem)*div_elem)
			vi=A_hor*sqrt(elem_area(elem)/scale_area)
			Visc(nz,elem)=0.2_WP*min(elem_area(elem)*sqrt((Div_c*(xe**2+ye**2) &
				+ Leith_c*(leithx**2+leithy**2))*elem_area(elem)) &
				+ vi, elem_area(elem)/dt)
		end do                        !! 0.1 here comes from (2S)^{3/2}/pi^3
		do nz=nl1+1, nl-1
			Visc(nz, elem)=0.0_WP
		end do		    
	END DO 
    
	allocate(aux(nl-1,myDim_nod2D+eDim_nod2D))
	DO nt=1,2 
		DO n=1, myDim_nod2D 
			DO nz=1, nlevels_nod2D(n)-1
				dz=0.0_WP
				vi=0.0_WP
				DO k=1, nod_in_elem2D_num(n)
					elem=nod_in_elem2D(k,n)
					dz=dz+elem_area(elem)
					vi=vi+Visc(nz,elem)*elem_area(elem)
				END DO
				aux(nz,n)=vi/dz
			END DO
		END DO 
		call exchange_nod(aux)
		do elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:,elem)
			nl1=nlevels(elem)-1
			Do nz=1, nl1 
				Visc(nz,elem)=sum(aux(nz,elnodes))/3.0_WP
			END DO
			DO nz=nl1+1, nl-1
				Visc(nz,elem)=0.0_WP
			END Do
		end do
	end do
	call exchange_elem(Visc)
	deallocate(aux)
END subroutine h_viscosity_leith
! =======================================================================
!A collection of viscosity routines
!(i)          visc_filt_bh == biharmonic
!(ii)         visc_filt_bh2== biharmonic, with less scale selective viscosity
!(iii)        visc_filt_h  == harmonic 
!(iv)         visc_filt_h_backscatter == harmonic, with backscatter (the least dissipative)
!In qualitative terms, viscosity in (i) is as the Leith, and in the rest it is as 
!the Smagorinsky. 
!(a) In all cases except (iii) auxliary arrays are allocated/deallocated inside. 
!May be they should be made static instead.  
!(b) There are tuning parameters inside. They should go to namelist.
!Poor-man-Backscatter does not bring a lot if added to the biharmonic operator,
!presumably because we need more filtering.
!
!In test cases, using viscosity_h_backscatter I am able to return back with a factor 1.5,
!but in the real world we have to start from 1.0 and see what we get.
!
!Important:  We do not need to call routines computing the Leith viscosity and velocity gradient.
!This makes code faster. Do not forget this!
! ===================================================================
SUBROUTINE viscosity_filt_bh(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity Laplacian, i.e., on an analog of
! the Leith viscosity (Lapl==second derivatives)
! \nu=|3u_c-u_n1-u_n2-u_n3|*sqrt(S_c)/100. There is an additional term
! in viscosity that is proportional to the velocity amplitude squared.
! The coefficient has to be selected experimentally.
real(kind=8)  :: u1, v1, vi, len
integer       :: ed, el(2), nz
real(kind=8), allocatable  :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"
!
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_WP
 V_c=0.0_WP
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     U_c(nz,el(1))=U_c(nz,el(1))-u1
     U_c(nz,el(2))=U_c(nz,el(2))+u1
     V_c(nz,el(1))=V_c(nz,el(1))-v1
     V_c(nz,el(2))=V_c(nz,el(2))+v1
    END DO 
 END DO
 
 Do ed=1,myDim_elem2D
    len=sqrt(elem_area(ed))                     
    len=len/100.0_WP
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coef. because of 
     ! division by area in the end 
     vi=dt*max(sqrt(U_c(nz,ed)**2+V_c(nz,ed)**2), 5.0_WP*(U_c(nz,ed)**2+V_c(nz,ed)**2))*len   
     U_c(nz,ed)=-U_c(nz,ed)*vi                             
     V_c(nz,ed)=-V_c(nz,ed)*vi
    END DO
 end do
 call exchange_elem(U_c)
 call exchange_elem(V_c)
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
     el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
     u1=(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=(V_c(nz,el(1))-V_c(nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
     
 deallocate(V_c,U_c)
end subroutine viscosity_filt_bh
! ===================================================================
SUBROUTINE viscosity_filt_bh2(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity differences, and is introduced symmetrically 
! into both stages of biharmonic operator
! On each edge, \nu=sqrt(|u_c1-u_c2|*sqrt(S_c1+S_c2)/100)
! The effect is \nu^2
! Quadratic in velocity term can be introduced if needed.
real(kind=8)  :: u1, v1, vi, len
integer       :: ed, el(2), nz
real(kind=8), allocatable  :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"
!
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_WP
 V_c=0.0_WP
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    len=sqrt(sum(elem_area(el)))/100.0_WP
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     vi=sqrt(len*sqrt(u1*u1+v1*v1))
     u1=u1*vi
     v1=v1*vi    
     U_c(nz,el(1))=U_c(nz,el(1))-u1
     U_c(nz,el(2))=U_c(nz,el(2))+u1
     V_c(nz,el(1))=V_c(nz,el(1))-v1
     V_c(nz,el(2))=V_c(nz,el(2))+v1
    END DO 
 END DO
 
 call exchange_elem(U_c)
 call exchange_elem(V_c)
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
     el=edge_tri(:,ed)
     len=sqrt(sum(elem_area(el)))/100.0_WP
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     vi=-dt*sqrt(len*sqrt(u1*u1+v1*v1))
     u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
     
 deallocate(V_c,U_c)
end subroutine viscosity_filt_bh2
! ===================================================================
SUBROUTINE viscosity_filt_h(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE

real(kind=8)  :: u1, v1, le, vi 
integer       :: nz, ed, el(2)
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! An analog of harmonic viscosity operator.  
 ! The contribution from boundary edges is neglected (free slip).
 !
 ! The  viscosity coefficient \nu=|u_c1-u_c2|*le/30
 !
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    le=sqrt(sum(elem_area(el)))/30.0_WP
    DO  nz=1,minval(nlevels(el))-1
      u1=UV(1,nz,el(1))-UV(1,nz,el(2))
      v1=UV(2,nz,el(1))-UV(2,nz,el(2))
      vi=dt*sqrt(u1*u1+v1*v1)*le
      u1=u1*vi
      v1=v1*vi
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
end subroutine viscosity_filt_h
! ===================================================================
SUBROUTINE viscosity_filt_h_backscatter(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE

real(kind=8)  :: u1, v1, le, vi 
integer       :: nz, ed, el(2), nelem(3),k, elem
real(kind=8), allocatable  ::  U_b(:,:), V_b(:,:), U_c(:,:), V_c(:,:)  
type(t_mesh), intent(in) , target :: mesh

#include "associate_mesh.h"

 ! An analog of harmonic viscosity operator.
 ! Same as visc_filt_h, but with the backscatter. 
 ! Here the contribution from squared velocities is added to the viscosity.    
 ! The contribution from boundary edges is neglected (free slip). 

 ed=myDim_elem2D+eDim_elem2D
 allocate(U_b(nl-1,ed), V_b(nl-1, ed))
 ed=myDim_nod2D+eDim_nod2D
 allocate(U_c(nl-1,ed), V_c(nl-1,ed))
 U_b=0.0_WP
 V_b=0.0_WP
 U_c=0.0_WP
 V_c=0.0_WP
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    le=sqrt(sum(elem_area(el)))/easy_bs_scale
    DO  nz=1,minval(nlevels(el))-1
      u1=UV(1,nz,el(1))-UV(1,nz,el(2))
      v1=UV(2,nz,el(1))-UV(2,nz,el(2))
      vi=dt*max(sqrt(u1*u1+v1*v1),10*(u1*u1+v1*v1))*le
      u1=u1*vi
      v1=v1*vi
     U_b(nz,el(1))=U_b(nz,el(1))-u1/elem_area(el(1))
     U_b(nz,el(2))=U_b(nz,el(2))+u1/elem_area(el(2))
     V_b(nz,el(1))=V_b(nz,el(1))-v1/elem_area(el(1))
     V_b(nz,el(2))=V_b(nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
 call exchange_elem(U_b)
 call exchange_elem(V_b)
 ! ===========
 ! Compute smoothed viscous term: 
 ! ===========
   DO ed=1, myDim_nod2D 
    DO nz=1, nlevels_nod2D(ed)-1
       vi=0.0_WP
       u1=0.0_WP
       v1=0.0_WP
       DO k=1, nod_in_elem2D_num(ed)
           elem=nod_in_elem2D(k,ed)
           vi=vi+elem_area(elem)
           u1=u1+U_b(nz,elem)*elem_area(elem)
           v1=v1+V_b(nz,elem)*elem_area(elem)
        END DO
	U_c(nz,ed)=u1/vi
        V_c(nz,ed)=v1/vi
    END DO
    END DO
  call exchange_nod(U_c)
  call exchange_nod(V_c)
  do ed=1, myDim_elem2D
         nelem=elem2D_nodes(:,ed)
         Do nz=1, nlevels(ed)-1
         UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -easy_bs_return*sum(U_c(nz,nelem))/3.0_WP
         UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -easy_bs_return*sum(V_c(nz,nelem))/3.0_WP
         END DO
  end do
 deallocate(V_c,U_c,V_b,U_b)
end subroutine viscosity_filt_h_backscatter

! ===================================================================

