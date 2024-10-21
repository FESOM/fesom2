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
module visc_filt_dbcksc_interface
  interface
    subroutine visc_filt_dbcksc(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module
module backscatter_coef_interface
  interface
    subroutine backscatter_coef(mesh)
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

!=========================================================================


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
use visc_filt_dbcksc_interface
use backscatter_coef_interface
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
CASE (8)
     call backscatter_coef(mesh)
     call visc_filt_dbcksc(mesh) 
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
 UV_total_tend=0.0_8
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
        UV_total_tend(1,nz,el(1))=UV_total_tend(1,nz,el(1))-u1/elem_area(el(1))
        UV_total_tend(1,nz,el(2))=UV_total_tend(1,nz,el(2))+u1/elem_area(el(2))
        UV_total_tend(2,nz,el(1))=UV_total_tend(2,nz,el(1))-v1/elem_area(el(1))
        UV_total_tend(2,nz,el(2))=UV_total_tend(2,nz,el(2))+v1/elem_area(el(2))
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
 UV_total_tend=0.0_8
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
     UV_total_tend(1,nz,el(1))=UV_total_tend(1,nz,el(1))-u1/elem_area(el(1))
     UV_total_tend(1,nz,el(2))=UV_total_tend(1,nz,el(2))+u1/elem_area(el(2))
     UV_total_tend(2,nz,el(1))=UV_total_tend(2,nz,el(1))-v1/elem_area(el(1))
     UV_total_tend(2,nz,el(2))=UV_total_tend(2,nz,el(2))+v1/elem_area(el(2))
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
 UV_total_tend=0.0_8
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
        UV_total_tend(1,nz,el(1))=UV_total_tend(1,nz,el(1))-u1/elem_area(el(1))
        UV_total_tend(1,nz,el(2))=UV_total_tend(1,nz,el(2))+u1/elem_area(el(2))
        UV_total_tend(2,nz,el(1))=UV_total_tend(2,nz,el(1))-v1/elem_area(el(1))
        UV_total_tend(2,nz,el(2))=UV_total_tend(2,nz,el(2))+v1/elem_area(el(2))        
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
    INTEGER, PARAMETER    :: R=selected_int_kind(9)
    !integer(R), save :: seed
    !real(kind=WP) ::ran
    !real                  :: t10  !for random initialization using the computer clock
    integer :: istep

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
    UV_total_tend=0.0_8
    UV_back_tend=0.0_8
    UV_dis_tend=0.0_8
    !call seed randomizer and subroutine myrandom
    !IF(istep==1) THEN
    !t10 = MPI_Wtime()
    !seed = 345 !int(t10)
    !END IF
    !call myrandom(seed, ran)
    do ed=1, myDim_elem2D
        nelem=elem2D_nodes(:,ed)
        Do nz=1, nlevels(ed)-1
            UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -easy_bs_return*sum(U_c(nz,nelem))/3.0_WP
            UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -easy_bs_return*sum(V_c(nz,nelem))/3.0_WP
            UV_total_tend(1,nz,ed)=UV_total_tend(1,nz,ed)+U_b(nz,ed) -easy_bs_return*sum(U_c(nz,nelem))/3.0_WP
            UV_total_tend(2,nz,ed)=UV_total_tend(2,nz,ed)+V_b(nz,ed) -easy_bs_return*sum(V_c(nz,nelem))/3.0_WP
            UV_back_tend(1,nz,ed)=UV_back_tend(1,nz,ed)-easy_bs_return*sum(U_c(nz,nelem))/3.0_WP
            UV_back_tend(2,nz,ed)=UV_back_tend(2,nz,ed)-easy_bs_return*sum(V_c(nz,nelem))/3.0_WP
            UV_dis_tend(1,nz,ed)=UV_dis_tend(1,nz,ed)+U_b(nz,ed)
            UV_dis_tend(2,nz,ed)=UV_dis_tend(2,nz,ed)+V_b(nz,ed)
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
 UV_total_tend=0.0_8
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
     UV_total_tend(1,nz,el(1))=UV_total_tend(1,nz,el(1))-u1/elem_area(el(1))
     UV_total_tend(1,nz,el(2))=UV_total_tend(1,nz,el(2))+u1/elem_area(el(2))
     UV_total_tend(2,nz,el(1))=UV_total_tend(2,nz,el(1))-v1/elem_area(el(1))
     UV_total_tend(2,nz,el(2))=UV_total_tend(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
 deallocate(V_c,U_c)
end subroutine visc_filt_bilapl
! ===================================================================
! Viscosity depends on velocity differences, and is introduced symmetrically 
! into both stages of biharmonic operator
! On each edge, \nu=sqrt(|u_c1-u_c2|*sqrt(S_c1+S_c2)/100)
! The effect is \nu^2
! Quadratic in velocity term can be introduced if needed.
SUBROUTINE visc_filt_bidiff(mesh)
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
 UV_total_tend=0.0_8
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
     UV_total_tend(1,nz,el(1))=UV_total_tend(1,nz,el(1))-u1/elem_area(el(1))
     UV_total_tend(1,nz,el(2))=UV_total_tend(1,nz,el(2))+u1/elem_area(el(2))
     UV_total_tend(2,nz,el(1))=UV_total_tend(2,nz,el(1))-v1/elem_area(el(1))
     UV_total_tend(2,nz,el(2))=UV_total_tend(2,nz,el(2))+v1/elem_area(el(2))
    END DO 
 END DO
 deallocate(V_c, U_c)
end subroutine visc_filt_bidiff
! ===================================================================


! ===================================================================
SUBROUTINE visc_filt_dbcksc(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
USE g_support
IMPLICIT NONE

real(kind=8)  :: u1, v1, le(2), len, len2, crosslen, vi, uke1 
integer       :: nz, ed, el(2), nzmin, nzmax
real(kind=8), allocatable  :: U_c(:,:), V_c(:,:), UV_back(:,:,:), UV_dis(:,:,:), uke_d(:,:) 
real(kind=8), allocatable  :: uuu(:,:)
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
 ! An analog of harmonic viscosity operator.  
 ! It adds to the rhs(0) Visc*(u1+u2+u3-3*u0)/area
 ! on triangles, which is Visc*Laplacian/4 on equilateral triangles. 
 ! The contribution from boundary edges is neglected (free slip). 
 ! Filter is applied twice. 

ed=myDim_elem2D+eDim_elem2D
allocate(U_c(nl-1,ed), V_c(nl-1, ed)) 
allocate(UV_back(2,nl-1,ed), UV_dis(2,nl-1, ed)) 
allocate(uke_d(nl-1,ed)) 
allocate(uuu(nl-1,ed))
 
 U_c=0.0_8
 V_c=0.0_8
 UV_back=0.0_8
 UV_dis=0.0_8
 uke_d=0.0_8
!---------------------- 
 ! New linew instead:
!---------------------- 

  DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    ! New viscosity lines
    len=sqrt(sum(elem_area(el)))
    nzmin = 1
    nzmax = minval(nlevels(el))
    DO  nz=nzmin,nzmax-1
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     ! New viscosity lines
     vi=u1*u1+v1*v1
     vi=sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len)
     u1=u1*vi
     v1=v1*vi
     
     U_c(nz,el(1))=U_c(nz,el(1))-u1
     U_c(nz,el(2))=U_c(nz,el(2))+u1
     V_c(nz,el(1))=V_c(nz,el(1))-v1
     V_c(nz,el(2))=V_c(nz,el(2))+v1
    END DO 
 END DO




 !DO ed=1, myDim_edge2D+eDim_edge2D
 !   if(myList_edge2D(ed)>edge2D_in) cycle
 !   el=edge_tri(:,ed)
 !   DO  nz=1,minval(nlevels(el))-1
 !    u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
 !    v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
 !    !if(mype==0 .AND. ed==1)  write(*,*) 'u1', u1
 !    
 !    U_c(nz,el(1))=U_c(nz,el(1))-u1
 !    U_c(nz,el(2))=U_c(nz,el(2))+u1
 !    V_c(nz,el(1))=V_c(nz,el(1))-v1
 !    V_c(nz,el(2))=V_c(nz,el(2))+v1
 !   END DO 
 !END DO
 
!if(mype==0) write(*,*) 'U_c_1', maxval(U_c)
 

! Do ed=1,myDim_elem2D
!    len=sqrt(elem_area(ed))                     
!    len=dt*len/30.0_8
!    Do nz=1,nlevels(ed)-1
!     ! vi has the sense of harmonic viscosity coefficient because of 
!     ! the division by area in the end 
!     ! ====
!     ! Case 1 -- an analog to the third-order upwind (vi=|u|l/12)
!     ! ====
!     vi=max(0.2_8,sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2))*len 
!     U_c(nz,ed)=-U_c(nz,ed)*vi                             
!     V_c(nz,ed)=-V_c(nz,ed)*vi
!    END DO
! end do



 call exchange_elem(U_c)
 call exchange_elem(V_c) 

 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    le=edge_dxdy(:,ed)
    le(1)=le(1)*sum(elem_cos(el))*0.25_8
    len=sqrt(le(1)**2+le(2)**2)*r_earth
    le(1)=edge_cross_dxdy(1,ed)-edge_cross_dxdy(3,ed)
    le(2)=edge_cross_dxdy(2,ed)-edge_cross_dxdy(4,ed)
    crosslen=sqrt(le(1)**2+le(2)**2) 
    ! Weighting new viscosity
    len2=sqrt(sum(elem_area(el)))
    nzmin = 1
    nzmax = minval(nlevels(el))
    Do nz=nzmin,nzmax-1
     vi=dt*len*(v_back(nz,el(1))+v_back(nz,el(2)))/crosslen

     !Backscatter contribution
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))*vi
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))*vi
     
     !UKE diffusion
     vi=dt*len*(K_back*sqrt(elem_area(el(1))/scale_area)+K_back*sqrt(elem_area(el(2))/scale_area))/crosslen

     uke1=(uke(nz,el(1))-uke(nz,el(2)))*vi

     
     UV_back(1,nz,el(1))=UV_back(1,nz,el(1))-u1/elem_area(el(1))
     UV_back(1,nz,el(2))=UV_back(1,nz,el(2))+u1/elem_area(el(2))
     UV_back(2,nz,el(1))=UV_back(2,nz,el(1))-v1/elem_area(el(1))
     UV_back(2,nz,el(2))=UV_back(2,nz,el(2))+v1/elem_area(el(2))  
     
     !Correct scaling for the diffusion?
     uke_d(nz,el(1))=uke_d(nz,el(1))-uke1/elem_area(el(1))
     uke_d(nz,el(2))=uke_d(nz,el(2))+uke1/elem_area(el(2))
     
     
     
     !Biharmonic contribution
     !u1=(U_c(nz,el(1))-U_c(nz,el(2)))
     !v1=(V_c(nz,el(1))-V_c(nz,el(2)))
     
     ! New viscosity param part
     u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
     vi=u1*u1+v1*v1
     vi=-dt*sqrt(max(gamma0, max(gamma1*sqrt(vi), gamma2*vi))*len2) 
     u1=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
     v1=vi*(V_c(nz,el(1))-V_c(nz,el(2)))

     UV_dis(1,nz,el(1))=UV_dis(1,nz,el(1))-u1/elem_area(el(1))
     UV_dis(1,nz,el(2))=UV_dis(1,nz,el(2))+u1/elem_area(el(2))
     UV_dis(2,nz,el(1))=UV_dis(2,nz,el(1))-v1/elem_area(el(1))
     UV_dis(2,nz,el(2))=UV_dis(2,nz,el(2))+v1/elem_area(el(2))
     
    END DO 
 END DO



call exchange_elem(UV_back)
uuu=UV_back(1,:,:)
call smooth_elem(uuu,smooth_back_tend, mesh)
UV_back(1,:,:)=uuu
uuu=UV_back(2,:,:)
call smooth_elem(uuu,smooth_back_tend, mesh)
UV_back(2,:,:)=uuu 

!DO  nz=1, nl-1
!    uuu=0.0_8
!    uuu=UV_back(1,nz,:)
!    call smooth_elem(uuu,smooth_back_tend, mesh)
!    UV_back(1,nz,:)=uuu
!    uuu=0.0_8
!    uuu=UV_back(2,nz,:)
!    call smooth_elem(uuu,smooth_back_tend, mesh)
!    UV_back(2,nz,:)=uuu 
!END DO

! DO ed=1, myDim_elem2D
!     DO  nz=1,nlevels(ed)-1
!     UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_dis(1,nz,ed)+UV_back(1,nz,ed)
!     UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_dis(2,nz,ed)+UV_back(2,nz,ed)               
!     END DO 
! END DO



 DO ed=1, myDim_elem2D
     nzmin = 1
     nzmax = nlevels(ed)
     Do nz=nzmin,nzmax-1
     UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_dis(1,nz,ed)+UV_back(1,nz,ed)
     UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_dis(2,nz,ed)+UV_back(2,nz,ed)               
     END DO 
 END DO
 
 UV_dis_tend=UV_dis!+UV_back
 UV_total_tend=UV_dis+UV_back
 UV_back_tend=UV_back
 uke_dif=uke_d

 call uke_update(mesh)

 deallocate(V_c,U_c)
 deallocate(UV_dis,UV_back) 
 deallocate(uke_d)
 deallocate(uuu)



end subroutine visc_filt_dbcksc
!===========================================================================
SUBROUTINE backscatter_adv(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
type(t_mesh), intent(in) , target :: mesh

integer                  :: n, nz, el1, el2
integer                  :: nl1, nl2, nod(2), el, ed, k, nle
real(kind=WP)            :: contr1(1:mesh%nl-1), contr2(1:mesh%nl-1)
real(kind=WP)            :: advnode(mesh%nl-1,myDim_nod2d+eDim_nod2D)
real(kind=WP)            :: wuke(1:mesh%nl)


#include "associate_mesh.h"

advnode=0.0_8
 
!vertical advection
do n=1,myDim_nod2d
   nl1 = nlevels_nod2D(n)-1
   wuke(1:nl1+1) = 0._WP
   
   ! here we find an average of uke between vertical layers
   do k=1,nod_in_elem2D_num(n)
       el = nod_in_elem2D(k,n)
       nle = nlevels(el)-1
   
       !wuke(1) = wuke(1) + uke(1,n)*elem_area(el)
       !wuke(2:nl1) = wuke(2:nl1) + 0.5_WP*(uke(2:nl1,n)+uke(1:nl1-1,n))*elem_area(el)
       wuke(1) = wuke(1) + uke(1,el)*elem_area(el)
       wuke(2:nle) = wuke(2:nle) + 0.5_WP*(uke(2:nle,el)+uke(1:nle-1,el))*elem_area(el)      
       


   end do
   ! here we find an average horizontal velocity between vertical layers 
   
   wuke(1:nl) = wuke(1:nl1)*Wvel_e(1:nl1,n)

   do nz=1,nl1
      ! Here 1/3 because 1/3 of the area is related to the node
      advnode(nz,n) = - (wuke(nz) - wuke(nz+1) ) / (3._WP*hnode(nz,n)) 
   enddo
   
   advnode(nl1+1:nl-1,n) = 0._WP

end do

!horizontal advection
DO ed=1, myDim_edge2D
   nod = edges(:,ed) 
   el1  = edge_tri(1,ed)   
   el2  = edge_tri(2,ed)
   nl1 = nlevels(el1)-1

   ! here we calculate fluxes throught the each edge 
   !for the edges that are on the boundary and have only one element 
   contr1(1:nl1) = UV(2,1:nl1,el1)*edge_cross_dxdy(1,ed) - UV(1,1:nl1,el1)*edge_cross_dxdy(2,ed)  
   
   ! if we have the other element as well
   if (el2>0) then
      nl2 = nlevels(el2)-1
      
      contr2(1:nl2) = - UV(2,1:nl2,el2)*edge_cross_dxdy(3,ed) + UV(1,1:nl2,el2)*edge_cross_dxdy(4,ed)
      
      contr1(nl1+1:max(nl1,nl2)) = 0._WP
      contr2(nl2+1:max(nl1,nl2)) = 0._WP
      
      if (nod(1) <= myDim_nod2d) then
         do nz=1, max(nl1,nl2)
            advnode(nz,nod(1)) = advnode(nz,nod(1)) + contr1(nz)*uke(nz,el1) + contr2(nz)*uke(nz,el2)             
         end do
      endif
      
      if (nod(2) <= myDim_nod2d) then
         do nz=1, max(nl1,nl2)
            advnode(nz,nod(2)) = advnode(nz,nod(2)) - contr1(nz)*uke(nz,el1) - contr2(nz)*uke(nz,el2)
         end do
      endif
   
   else  ! ed is a boundary edge, there is only the contribution from el1
      if (nod(1) <= myDim_nod2d) then
         do nz=1, nl1
            advnode(nz,nod(1)) = advnode(nz,nod(1)) + contr1(nz)*uke(nz,el1)
         end do
      endif
      ! second edge node
      if  (nod(2) <= myDim_nod2d) then
         do nz=1, nl1
            advnode(nz,nod(2)) = advnode(nz,nod(2)) - contr1(nz)*uke(nz,el1)
         end do
      endif
   endif

end do

! Now we need to multiply by the segment area (here we already make a loop by nodes)

do n=1,myDim_nod2d
   nl1 = nlevels_nod2D(n)-1
   
   advnode(1:nl1,n) = advnode(1:nl1,n) * area_inv(1:nl1,n)
   
end do

call exchange_nod(advnode)

! Now we calculate advection as an average across the nodes

do n=1,myDim_nod2d
   ul1 = ulevels_nod2D(n)-1
   advnode(ul1:nl1,n) = advnode(1:nl1,n) *areasvol_inv(1:nl1,n)
end do !-->do n=1,myDim_nod2d



do el=1, myDim_elem2D
   nl1 = nlevels(el)-1
   uke_adv(1:nl1,el) = uke_adv(1:nl1,el) &
        + elem_area(el)*(advnode(1:nl1,elem2D_nodes(1,el)) &
        + advnode(1:nl1,elem2D_nodes(2,el)) & 
        + advnode(1:nl1,elem2D_nodes(3,el)))/ 3.0_WP
   
end do

do el=1, myDim_elem2D
    do nz=1,nlevels(el)-1   
        
        uke_adv(nz,el)=dt*uke_adv(nz,el)/elem_area(el)

    end do
end do


!if(mype==0) write(*,*) 'elem_area', elem_area(1)

end subroutine backscatter_adv


!===========================================================================
SUBROUTINE backscatter_coef(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
type(t_mesh), intent(in) , target :: mesh
integer                 :: elem, nz
#include "associate_mesh.h"

!Potentially add the Rossby number scaling to the script...
!check if sign is right! Different in the Jansen paper
!Also check with the normalization by area; as before we use element length sqrt(2*elem_area(ed))

v_back=0.0_8
DO  elem=1, myDim_elem2D 
     DO  nz=1,nlevels(elem)-1
!v_back(1,ed)=c_back*sqrt(2.0_WP*elem_area(ed))*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))*(3600.0_WP*24.0_WP/tau_c)*4.0_WP/sqrt(2.0_WP*elem_area(ed))**2 !*sqrt(max(2.0_WP*uke(1,ed),0.0_WP))
!v_back(nz,elem)=-c_back*sqrt(4._8/sqrt(3.0_8)*elem_area(elem))*sqrt(max(2.0_8*uke(nz,elem),0.0_8)) !Is the scaling correct
    v_back(nz,elem)=min(-c_back*sqrt(elem_area(elem))*sqrt(max(2.0_8*uke(nz,elem),0.0_8)),0.2*elem_area(elem)/dt) !Is the scaling correct
    !v_back(nz,elem)=-c_back*sqrt(elem_area(elem))*sqrt(max(2.0_WP*uke(nz,elem),0.0_WP))
     END DO
END DO

call exchange_elem(v_back)

end subroutine backscatter_coef
!===========================================================================



SUBROUTINE uke_update(mesh)
USE MOD_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
USE g_support
USE g_rotate_grid
IMPLICIT NONE

!I had to change uke(:) to uke(:,:) to make output and restart work!!

!Why is it necessary to implement the length of the array? It doesn't work without!
!integer, intent(in)        :: t_levels
type(t_mesh), intent(in) , target :: mesh
!real(kind=8), dimension(:,:,:), intent(in)   :: UV_dis, UV_back
!real(kind=8), dimension(:,:), intent(in)   :: UV_dif
!real(kind=8), intent(in)   :: UV_dis(nl-1,myDim_elem2D+eDim_elem2D), UV_back(nl-1, myDim_elem2D+eDim_elem2D)
!real(kind=8), intent(in)   :: UV_dif(nl-1,myDim_elem2D+eDim_elem2D)
real(kind=8)		   :: hall, h1_eta, hnz, vol
integer			   :: elnodes(3), nz, ed, edi, node, j, elem, q
real(kind=8), allocatable  :: uuu(:), work_array(:), U_work(:,:), V_work(:,:), rosb_array(:,:), work_uv(:)
integer			   :: kk, nzmax, el, clock
real(kind=8)       :: c1, rosb, vel_u, vel_v, vel_uv, scaling, reso
real*8             :: c_min=0.5, f_min=1.e-6, r_max=200000., ex, ey, a1, a2, len_reg, dist_reg(2), t1, t2 ! Are those values still correct?
#include "associate_mesh.h"
!rosb_dis=1._8 !Should be variable to control how much of the dissipated energy is backscattered
!rossby_num=2

t1=MPI_Wtime()

ed=myDim_elem2D+eDim_elem2D
allocate(uuu(ed)) 

uke_back=0.0_8
uke_dis=0.0_8

!call PCA

DO ed=1, myDim_elem2D
DO nz=1, nlevels(ed)-1  
   uke_dis(nz,ed)=(UV(1,nz,ed)*UV_dis_tend(1,nz,ed)+UV(2,nz,ed)*UV_dis_tend(2,nz,ed))   
   uke_back(nz,ed)=(UV(1,nz,ed)*UV_back_tend(1,nz,ed)+UV(2,nz,ed)*UV_back_tend(2,nz,ed))

END DO
END DO

DO  nz=1,nl-1
    uuu=0.0_8
    uuu=uke_back(nz,:)
    call smooth_elem(uuu,smooth_back, mesh) !3) ?
    uke_back(nz,:)=uuu
END DO




!Timestepping use simple backward timestepping; all components should have dt in it, unless they need it twice
!Amplitudes should be right given the correction of the viscosities; check for all, also for biharmonic
!uke(1,ed)=uke(1,ed)-uke_dis(1,ed)-uke_back(1,ed)+uke_dif(1,ed)
ed=myDim_elem2D+eDim_elem2D
allocate(U_work(nl-1,myDim_nod2D+eDim_nod2D),V_work(nl-1,myDim_nod2D+eDim_nod2D))
allocate(work_uv(myDim_nod2D+eDim_nod2D))
allocate(rosb_array(nl-1,ed))
call exchange_elem(UV)
rosb_array=0._8
DO nz=1, nl-1
    work_uv=0._WP
    DO node=1, myDim_nod2D
       vol=0._WP
       U_work(nz,node)=0._WP 
       V_work(nz,node)=0._WP 
       DO j=1, nod_in_elem2D_num(node)
           elem=nod_in_elem2D(j, node)
           U_work(nz,node)=U_work(nz,node)+UV(1,nz,elem)*elem_area(elem)
           V_work(nz,node)=V_work(nz,node)+UV(2,nz,elem)*elem_area(elem)
           vol=vol+elem_area(elem)
       END DO
       U_work(nz,node)=U_work(nz,node)/vol
       V_work(nz,node)=U_work(nz,node)/vol
    END DO
    work_uv=U_work(nz,:)
    call exchange_nod(work_uv)
    U_work(nz,:)=work_uv
    work_uv=V_work(nz,:)
    call exchange_nod(work_uv)
    V_work(nz,:)=work_uv    
END DO

    DO el=1,myDim_elem2D
     DO nz=1, nlevels(el)-1     
        rosb_array(nz,el)=sqrt((sum(gradient_sca(1:3,el)*U_work(nz,elem2D_nodes(1:3,el)))-&
              sum(gradient_sca(4:6, el)*V_work(nz,elem2D_nodes(1:3,el))))**2+&
              (sum(gradient_sca(4:6, el)*U_work(nz,elem2D_nodes(1:3,el)))+&
              +sum(gradient_sca(1:3, el)*V_work(nz,elem2D_nodes(1:3,el))))**2)
!        hall=hall+hnz
     END DO
!     rosb_array(el)=rosb_array(el)/hall
    END DO
DO ed=1, myDim_elem2D
    scaling=1._WP
    IF(uke_scaling) then
      reso=sqrt(elem_area(ed)*4._wp/sqrt(3._wp))
      rosb=0._wp
      elnodes=elem2D_nodes(:, ed)   
      DO kk=1,3
        c1=0._wp
        nzmax=minval(nlevels(nod_in_elem2D(1:nod_in_elem2D_num(elnodes(kk)), elnodes(kk))), 1)
        !Vertical average; same scaling in the vertical
        DO nz=1, nzmax-1
         c1=c1+hnode_new(nz,elnodes(kk))*(sqrt(max(bvfreq(nz,elnodes(kk)), 0._WP))+sqrt(max(bvfreq(nz+1,elnodes(kk)), 0._WP)))/2.
        END DO
        c1=max(c_min, c1/pi) !ca. first baroclinic gravity wave speed limited from below by c_min
        !Cutoff K_GM depending on (Resolution/Rossby radius) ratio
        rosb=rosb+min(c1/max(abs(coriolis_node(elnodes(kk))), f_min), r_max)
      END DO
      rosb=rosb/3._8
      scaling=1._WP/(1._WP+(uke_scaling_factor*reso/rosb))!(4._wp*reso/rosb))
    END IF
    
    DO nz=1, nlevels(ed)-1  
    elnodes=elem2D_nodes(:,ed)
    
    !Taking out that one place where it is always weird (Pacific Southern Ocean)
    !Should not really be used later on, once we fix the issue with the 1/4 degree grid
    if(.not. (TRIM(which_toy)=="soufflet") .AND. .not. (TRIM(which_toy)=="nemo") ) then
      call elem_center(ed, ex, ey)
      !a1=-104.*rad
      !a2=-49.*rad
      call g2r(-104.*rad, -49.*rad, a1, a2)
      dist_reg(1)=ex-a1
      dist_reg(2)=ey-a2
      call trim_cyclic(dist_reg(1))
      dist_reg(1)=dist_reg(1)*elem_cos(ed)
      dist_reg=dist_reg*r_earth
      len_reg=sqrt(dist_reg(1)**2+dist_reg(2)**2)
    
    
      !if(mype==0) write(*,*) 'len_reg ', len_reg , ' and dist_reg' , dist_reg, ' and ex, ey', ex, ey, ' and a ', a1, a2
      rosb_array(nz,ed)=rosb_array(nz,ed)/max(abs(sum(coriolis_node(elnodes(:)))), f_min)
      !uke_dif(nz, ed)=scaling*(1-exp(-len_reg/300000))*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)!UV_dif(1,ed)
      uke_dis(nz,ed)=scaling*(1-exp(-len_reg/300000))*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)*uke_dis(nz,ed)
    else
      rosb_array(nz,ed)=rosb_array(nz,ed)/max(abs(sum(coriolis_node(elnodes(:)))), f_min)
      !uke_dif(nz, ed)=scaling*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)!UV_dif(1,ed)
      uke_dis(nz,ed)=scaling*1._8/(1._8+rosb_array(nz,ed)/rosb_dis)*uke_dis(nz,ed)
    end if
    
    END DO
END DO
deallocate(U_work, V_work)
deallocate(rosb_array)
deallocate(work_uv)
call exchange_elem(uke_dis)
call exchange_elem(uke_dif)
!call smooth_elem(uke_dis,smooth_dis, mesh)
DO nz=1, nl-1  
    uuu=uke_dis(nz,:)
    call smooth_elem(uuu,smooth_dis, mesh)
    uke_dis(nz,:)=uuu
END DO




call PCA

!IF(mype==0) write(*,*) 'OK '

call backscatter_adv(mesh)
call exchange_elem(uke_adv)


!IF(mype==0) write(*,*) "SB_PCs", SB_PC1(1), SB_PC2(1), SB_PC3(1), SB_PC4(1), SB_PC5(1), SB_PC6(1), SB_PC7(1) 
!IF(mype==0) write(*,*) "uke_adv", maxval(uke_adv(1,:))

DO ed=1, myDim_elem2D
    DO  nz=1,nlevels(ed)-1
    uke_rhs_old(nz,ed)=uke_rhs(nz,ed)
    uke_rhs(nz,ed)=-uke_dis(nz,ed)-uke_back(nz,ed)+uke_dif(nz,ed)-uke_adv(nz,ed)-uke(nz,ed)*0.001*(SB_EOF1(ed) &
    * SB_PC1(1)+SB_EOF2(ed)*SB_PC2(1)&
    + SB_EOF3(ed)*SB_PC3(1)+SB_EOF4(ed)*SB_PC4(1)+SB_EOF5(ed)*SB_PC5(1)+SB_EOF6(ed)*SB_PC6(1) &
    + SB_EOF7(ed)*SB_PC7(1))!+SB_EOF8(ed)*SB_PC8(1)+SB_EOF9(ed)*SB_PC9(1)+SB_EOF10(ed)*SB_PC10(1) &
    !+ SB_EOF11(ed)*SB_PC11(1)+SB_EOF12(ed)*SB_PC12(1)+SB_EOF13(ed)*SB_PC13(1)+SB_EOF14(ed)*SB_PC14(1) &
    !+ SB_EOF15(ed)*SB_PC15(1)+SB_EOF16(ed)*SB_PC16(1)+SB_EOF17(ed)*SB_PC17(1))
    
    EOF_check(nz,ed) = uke(nz,ed)*0.001*(SB_EOF1(ed)*SB_PC1(1)+SB_EOF2(ed)*SB_PC2(1) &
    + SB_EOF3(ed)*SB_PC3(1)+SB_EOF4(ed)*SB_PC4(1)+SB_EOF5(ed)*SB_PC5(1)+SB_EOF6(ed)*SB_PC6(1) &
    + SB_EOF7(ed)*SB_PC7(1))!+SB_EOF8(ed)*SB_PC8(1)+SB_EOF9(ed)*SB_PC9(1)+SB_EOF10(ed)*SB_PC10(1) &
    !+ SB_EOF11(ed)*SB_PC11(1)+SB_EOF12(ed)*SB_PC12(1)+SB_EOF13(ed)*SB_PC13(1)+SB_EOF14(ed)*SB_PC14(1) &
    !+ SB_EOF15(ed)*SB_PC15(1)+SB_EOF16(ed)*SB_PC16(1)+SB_EOF17(ed)*SB_PC17(1))
    
    uke(nz,ed)=uke(nz,ed)+1.5_8*uke_rhs(nz,ed)-0.5_8*uke_rhs_old(nz,ed)
    
    !abs(UV_rhs(1,nz,ed))
    !(abs(UV_rhs(1,nz,ed))+abs(UV_rhs(2,nz,ed)))
    
    END DO
END DO
call exchange_elem(uke)

!IF(mype==0) write(*,*) "SB_PCs", SB_PC1(1), SB_PC2(1), SB_PC3(1), SB_PC4(1), SB_PC5(1), SB_PC6(1), SB_PC7(1) 
!IF(mype==0) write(*,*) "uke_adv", maxval(uke_adv(1,:))
!IF(mype==0) write(*,*) "uke_back", maxval(uke_back(1,:))
!IF(mype==0) write(*,*) "uke_dis", minval(uke_dis(1,:))
!IF(mype==0) write(*,*) "EOF_check", maxval(EOF_check(1,:)), minval(EOF_check(1,:))


deallocate(uuu)
!t2=MPI_Wtime()
!IF(mype==0) write(*,*) 'Step took ', t2-t1
end subroutine uke_update

! ===================================================================