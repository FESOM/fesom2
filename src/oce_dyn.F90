! A set of routines for computing the horizonlal viscosity
! the control parameters (their default values) are:
! gamma0 (0.01 [m/s]), gamma1 (0.1 [no dim.]), gamma2 (10.[s/m]), Div_c [1.], Leith_c[1.?]
! 1. gamma0 has the dimension of velocity. It should be as small as possible, but in any case smaller than 0.01 m/s. 
!    All major ocean circulation models are stable with harmonic viscosity 0.01*len.
! 2. gamma1 is nondimensional. In commonly used Leith or Smagorinsky parameterizations it is C/pi^2=0.1 (C is about 1). 
!    We therefore try to follow this, allowing some adjustments (because our mesh is triangular, our resolution is different, etc.). 
!    We however, try to keep gamma1<0.1
! 3. gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on. It is only used in: 
!    (5) visc_filt_bcksct, (6) visc_filt_bilapl, (7) visc_filt_bidiff
! 4. Div_c  =1.    should be default
! 5. Leith_c=?    (need to be adjusted)
module h_viscosity_leith_interface
  interface
    subroutine h_viscosity_leith(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_harmon_interface
  interface
    subroutine visc_filt_harmon(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_hbhmix_interface
  interface
    subroutine visc_filt_hbhmix(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_biharm_interface
  interface
    subroutine visc_filt_biharm(option, mesh)
      use mod_mesh
      integer       :: option
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_bcksct_interface
  interface
    subroutine visc_filt_bcksct(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_bilapl_interface
  interface
    subroutine visc_filt_bilapl(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module visc_filt_bidiff_interface
  interface
    subroutine visc_filt_bidiff(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
! ===================================================================
! Contains routines needed for computations of dynamics.
! includes: update_vel, compute_vel_nodes
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
subroutine viscosity_filter(option, mesh)
use o_PARAM
use g_PARSUP
use MOD_MESH
use h_viscosity_leith_interface
use visc_filt_harmon_interface
use visc_filt_hbhmix_interface
use visc_filt_biharm_interface
use visc_filt_bcksct_interface
use visc_filt_bilapl_interface
use visc_filt_bidiff_interface
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
     ! Harmonic Leith parameterization
     ! ====
     call h_viscosity_leith(mesh)
     call visc_filt_harmon(mesh)
CASE (2)
     ! ===
     ! Laplacian+Leith+biharmonic background
     ! ===
     call h_viscosity_leith(mesh)
     call visc_filt_hbhmix(mesh)
CASE (3)
     ! ===
     ! Biharmonic Leith parameterization
     ! ===
     call h_viscosity_leith(mesh)
     call visc_filt_biharm(2, mesh)
CASE (4)
     ! ===
     ! Biharmonic+upwind-type
     ! ===
     call visc_filt_biharm(1, mesh)
CASE (5)
     call visc_filt_bcksct(mesh)
CASE (6)
     call visc_filt_bilapl(mesh)
CASE (7)
     call visc_filt_bidiff(mesh)
CASE DEFAULT
     if (mype==0) write(*,*) 'mixing scheme with option ' , option, 'has not yet been implemented'
     call par_ex
     stop
END SELECT
end subroutine viscosity_filter  
! ===================================================================
SUBROUTINE visc_filt_harmon(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE

real(kind=WP) :: u1, v1, le(2), len, vi 
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
    len=sqrt(sum(elem_area(el(1:2))))
    DO  nz=1,minval(nlevels(el))-1
        vi=0.5_WP*(Visc(nz,el(1))+Visc(nz,el(2)))
        vi=max(vi, gamma0*len)*dt ! limited from below by backgroung
        u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))*vi
        v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))*vi

        UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
        UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
        UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
        UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO
 END DO
end subroutine visc_filt_harmon
! ===================================================================
SUBROUTINE visc_filt_biharm(option, mesh)
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
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coefficient because of 
     ! the division by area in the end 
     ! ====
     ! Case 1 -- an analog to the third-order upwind (vi=gamma1 * |u| * l)
     ! ====
     vi=max(gamma0, gamma1*sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2))*len*dt
     U_c(nz,ed)=-U_c(nz,ed)*vi                             
     V_c(nz,ed)=-V_c(nz,ed)*vi
    END DO
 end do
 end if
if(option==2) then
 Do ed=1,myDim_elem2D
    len=sqrt(elem_area(ed))                     
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coefficient because of 
     ! the division by area in the end 
     ! ===   
     ! Case 2 -- Leith +background (do not forget to call h_viscosity_leith before using this option)
     ! ===
     vi=max(Visc(nz,ed), gamma0*len)*dt ! limited from below by backgroung
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
end subroutine visc_filt_biharm
! ===================================================================
SUBROUTINE visc_filt_hbhmix(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE

! An energy and momentum conserving version.
! We use the harmonic Leith viscosity + biharmonic background viscosity
!

real(kind=WP)                      :: u1, v1, vi, len, crosslen, le(2)
integer                            :: ed, el(2), nz
real(kind=WP), allocatable         :: U_c(:,:), V_c(:,:) 
type(t_mesh),  intent(in) , target :: mesh

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
     vi=dt*0.5_WP*(Visc(nz,el(1))+Visc(nz,el(2)))
     ! backgroung is added later (biharmonically)
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
    len=sqrt(elem_area(ed))
    Do nz=1,nlevels(ed)-1
     vi=dt*gamma0*len ! add biharmonic backgroung
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
end subroutine visc_filt_hbhmix

! ===================================================================

SUBROUTINE h_viscosity_leith(mesh)
!
! Coefficient of horizontal viscosity is a combination of the Leith (with Leith_c) and modified Leith (with Div_c)
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
			Visc(nz,elem)=min(gamma1*elem_area(elem)*sqrt((Div_c*(xe**2+ye**2) &
				+ Leith_c*(leithx**2+leithy**2))*elem_area(elem)), elem_area(elem)/dt)
		end do            !! 0.1 here comes from (2S)^{3/2}/pi^3
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
SUBROUTINE visc_filt_bcksct(mesh)
    USE MOD_MESH
    USE o_ARRAYS
    USE o_PARAM
    USE g_PARSUP
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE

    real(kind=8)  :: u1, v1, len, vi 
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
        len=sqrt(sum(elem_area(el)))
        DO  nz=1,minval(nlevels(el))-1
            u1=UV(1,nz,el(1))-UV(1,nz,el(2))
            v1=UV(2,nz,el(1))-UV(2,nz,el(2))
            vi=dt*max(gamma0, max(gamma1*sqrt(u1*u1+v1*v1), gamma2*(u1*u1+v1*v1)))*len
!            vi=dt*max(gamma0, gamma1*max(sqrt(u1*u1+v1*v1), gamma2*(u1*u1+v1*v1)))*len 
            !here gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on.
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
end subroutine visc_filt_bcksct

! ===================================================================
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity Laplacian, i.e., on an analog of
! the Leith viscosity (Lapl==second derivatives)
! \nu=|3u_c-u_n1-u_n2-u_n3|*sqrt(S_c)/100. There is an additional term
! in viscosity that is proportional to the velocity amplitude squared.
! The coefficient has to be selected experimentally.
SUBROUTINE visc_filt_bilapl(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE
real(kind=8)  :: u1, v1, vi, len
integer       :: ed, el(2), nz
real(kind=8), allocatable         :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
!
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_8
 V_c=0.0_8
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
    Do nz=1,nlevels(ed)-1
     ! vi has the sense of harmonic viscosity coef. because of 
     ! division by area in the end 
     u1=U_c(nz,ed)**2+V_c(nz,ed)**2
     vi=max(gamma0, max(gamma1*sqrt(u1), gamma2*u1))*len*dt
!     vi=max(gamma0, gamma1*max(sqrt(u1), gamma2*u1))*len*dt
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
end subroutine visc_filt_bilapl
! ===================================================================
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity differences, and is introduced symmetrically 
! into both stages of biharmonic operator
! On each edge, \nu=sqrt(|u_c1-u_c2|*sqrt(S_c1+S_c2)/100)
! The effect is \nu^2
! Quadratic in velocity term can be introduced if needed.
SUBROUTINE visc_filt_bidiff(mesh)
USE MOD_MESH
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE
real(kind=8)  :: u1, v1, vi, len
integer       :: ed, el(2), nz
real(kind=8), allocatable         :: U_c(:,:), V_c(:,:) 
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
!
ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
 U_c=0.0_8
 V_c=0.0_8
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    len=sqrt(sum(elem_area(el)))
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     vi=u1*u1+v1*v1
     vi=sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
!     vi=sqrt(max(gamma0, gamma1*max(sqrt(vi), gamma2*vi))*len)
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
     len=sqrt(sum(elem_area(el)))
    DO  nz=1,minval(nlevels(el))-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     vi=u1*u1+v1*v1
     vi=-dt*sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
!     vi=-dt*sqrt(max(gamma0, gamma1*max(sqrt(vi), gamma2*vi))*len)
     u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1/elem_area(el(1))
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1/elem_area(el(2))
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1/elem_area(el(1))
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
 deallocate(V_c, U_c)
end subroutine visc_filt_bidiff
! ===================================================================
