! VP rheology routines
!
!  Distributed memory version :::  It is not tested in the 
!  distributed memory version, but it was for the global memory version.
!  
!  Add to i_param: integer :: ice_VP_iter=0   (0 just standard solve, 1-.. Picard iterations)        
!  FESOM1.4 icestiff is available and it is used here.
!  FESOM2.0 ssh_stiff has the same structure and can be used for colptr and rowind.
!  icestiff_values have to be allocated separately. Its size is 
!  the same as ssh_stiff%values, i.e., ssh_stiff%rowptr(myDim_nod2d+1)-ssh_stiff%rowptr(1)   
! ===============================================
! The driving routine is VPdynamics; it calls  other routines.
! They include vp_matrix_rhs, vpbc and solve_vp.
! vp_matrix_rhs assembles matrix and rhs in a single cycle
! over elements.
! vpbc sets boundary conditions of no slip. Free-slip 
! boundary conditions are not implemented yet.
! solve_vp calls PETSc solver. Inside the routine there 
! are several variants for Pmode which defines how solver is working.
! It is recommended to try several choices.

!
! icestiff is used to store the velocity stiffness matrix.
! If implicit advection is selected icestiff is also used to store
! the advection stiffness matrix. Since the 
! matrices are updated on each time step there is no conflict.
 
 
! !!! ATTENTION Metrics terms are introduced, but not tested
!   
! The algorithm is basically of Zhang and Hibler 1997, but the
! operator part is used as recommended in Hutchings et al (2003)

module ice_VP_interfaces
  interface
  subroutine VPmatrix_rhs(mesh)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target  :: mesh
  end subroutine

  subroutine VPsolve(mesh, new_values)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target :: mesh
      integer,      intent(in)         :: new_values
  end subroutine

  subroutine VPdynamics(mesh)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target :: mesh
  end subroutine

  subroutine VPmake_implicit(mesh)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target :: mesh
  end subroutine

  subroutine VPbc(mesh)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target :: mesh
  end subroutine

  subroutine VP_init(mesh)
      use g_parsup
      use mod_mesh
      type(t_mesh), intent(in), target :: mesh
  end subroutine
  end interface
end module
! =====================================================================
! ============================================================================
subroutine VP_init(mesh)
use g_parsup
use mod_mesh
use i_ARRAYS
use ice_VP_interfaces, except_this_one => VP_init
IMPLICIT NONE 
type(t_mesh), intent(in), target  :: mesh
#include "associate_mesh.h"

allocate(ice_stiff_values(ssh_stiff%rowptr_loc(myDim_nod2D+1)-1))
allocate(rhs_u(myDim_nod2D+eDim_nod2D), rhs_v(myDim_nod2D+eDim_nod2D))
ice_stiff_values=0.
rhs_u=0.0_WP
rhs_v=0.0_WP
end subroutine VP_init

SUBROUTINE VPmatrix_rhs(mesh)
! Fills in the entries of the vp rheology stiffness matrix
! and simultaneously computes the rhs vector
use g_parsup
use mod_mesh
use o_MESH
use i_PARAM
use i_ARRAYS
USE o_ARRAYS
use i_therm_param
use ice_VP_interfaces, except_this_one => VPmatrix_rhs
IMPLICIT NONE 

REAL(kind=8)                  :: m_e, a_e, drag
REAL(kind=8)                  :: dx(3),dy(3), pressure, xi, eta
REAL(kind=8)                  :: epsu, epsv, delta
REAL(kind=8)                  :: inv_mass, elex, eley, vale, entries(3)
INTEGER                       :: col, elem, elnodes(3), n2, offset
INTEGER                       :: ipos, row, q, n, i,nini, nend
real(kind=8)                  :: ac, as, ue, ve
integer, allocatable          :: n_num(:)
integer                       :: is, ie 

type(t_mesh), intent(in), target  :: mesh

#include "associate_mesh.h"

allocate(n_num(myDim_nod2D+eDim_nod2D))
ac=cos(theta_io)
as=sin(theta_io)
! The cos part is always implicit, the sin part --explicit

 do row=1,myDim_nod2D       ! Put the mass matrix part+ocean/ice friction part
                            ! + the Coriolis part. This is done in a lumped way,
			    ! so the cycle is over the owned nodes   
                            ! PP  (Distrib. memory)           
                            ! PP row=myList_nod2D(i)   
   
    inv_mass=(rhoice*m_ice(row)+rhosno*m_snow(row))
    inv_mass=max(inv_mass, 9.0_8)      ! Limit the weighted mass 
                                       ! if it is too small
    inv_mass=a_ice(row)/inv_mass
    drag=sqrt((u_ice(row)-u_w(row))**2+(v_ice(row)-v_w(row))**2)* &
         Cd_oce_ice*density_0*inv_mass
                    ! PP DO q=nini,nend 

!    offset=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)
!    ice_stiff_values(offset+1)=area(1,row)*(1.0/ice_dt+drag*ac)  ! Diagonal value
!    DO q=2,ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(row)
!       ice_stiff_values(offset+q)=0.0_8                      ! Initialize with 0.0
!    END DO

    !NR By construction, the first entry for each row is the diagonal one.
    !NR We can replace the loop with the if and colind_loc by a simple
    !NR direct statement.by a simple direct statement.
!!$    DO q=ssh_stiff%rowptr_loc(row), ssh_stiff%rowptr_loc(row+1)-1
!!$       if (ssh_stiff%colind_loc(q)==row) then
!!$          ice_stiff_values(q)=area(1, row)*(1.0/ice_dt+drag*ac)
!!$       else
!!$          ice_stiff_values(q)=0.0_8
!!$       end if
!!$    END DO
    is = ssh_stiff%rowptr_loc(row)
    ie = ssh_stiff%rowptr_loc(row+1)-1
    ice_stiff_values(is)      = area(1, row)*(1.0/ice_dt+drag*ac)
    ice_stiff_values(is+1:ie) = 0. 

    ! =========
     rhs_u(row)=area(1, row)*(rhs_m(row)/ice_dt+drag*ac*u_w(row) &
                +coriolis_node(row)*v_ice(row)+&
                stress_atmice_x(row)*inv_mass-&
		drag*as*(v_w(row)-v_ice(row)))
     rhs_v(row)=area(1, row)*(rhs_a(row)/ice_dt+drag*ac*v_w(row) &
                -coriolis_node(row)*u_ice(row)+&
                stress_atmice_y(row)*inv_mass+&
		drag*as*(u_w(row)-u_ice(row)))
 end do

! Cycle over the elements to assemble the rheology part. 
! Only the owned part of the stiffness matrix is filled 
  vale=1.0_8/(ellipse**2) 
  DO elem=1, myDim_elem2D      ! PP n=1, myDim_elem2D     
                               ! PP elem=myList_elem2D(n) 
     elnodes=elem2D_nodes(:,elem)
     
     delta=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (delta==0) CYCLE       ! Skip if there is no ice
       
             dx=gradient_sca(1:3,elem)    !! FESOM2
             dy=gradient_sca(4:6,elem)    !! FESOM2 

	  ! Deformation rate tensor on element elem 
	  ! (metric contribution is below):
     eps11(elem)=sum(dx*u_ice(elnodes))
     eps22(elem)=sum(dy*v_ice(elnodes))
     epsu=sum(dy*u_ice(elnodes))
     epsv=sum(dx*v_ice(elnodes))
     eps12(elem)=0.5_8*(epsu+epsv)
         ! delta and mean ice properties on element:
     delta=(eps11(elem)**2+eps22(elem)**2)*(1.0_8+vale)+4.0_8*vale*eps12(elem)**2 + &
            2.0_8*eps11(elem)*eps22(elem)*(1.0_8-vale)
     delta=sqrt(delta)
     m_e=sum(m_ice(elnodes))/3.0
     a_e=sum(a_ice(elnodes))/3.0
         ! pressure and viscosities 
     pressure=0.5_8*pstar*m_e*exp(-c_pressure*(1.0_8-a_e)) 
     !xi=pressure/max(delta,delta_min)
     xi = pressure/(delta+delta_min)   ! new
     pressure=xi*delta                 ! new
     eta=xi*vale
         ! contribution from ssh
     elex=-g*sum(dx*elevation(elnodes))*elem_area(elem)/3.0
     eley=-g*sum(dy*elevation(elnodes))*elem_area(elem)/3.0
     ue=sum(u_ice(elnodes))/3.0
     ve=sum(v_ice(elnodes))/3.0
     DO i=1,3  ! Cycle over rows elem contributes to
	row=elnodes(i)
        if(row>myDim_nod2D) cycle       !! PP if(part(row).ne.mype) cycle

        DO q=SSH_stiff%rowptr_loc(row), SSH_stiff%rowptr_loc(row+1)-1
           n2=SSH_stiff%colind_loc(q)
           n_num(n2)=q
        END DO

        inv_mass=(rhoice*m_ice(row)+rhosno*m_snow(row)) !/a_ice(row)
        inv_mass=max(inv_mass, 9.0_8)        ! Limit the weighted mass 
                                             ! if it is too small
        inv_mass=1.0/inv_mass
	
	!====== Rheology part
	                 ! u,v martices: contains phi_x(eta+xi)du/dx +
		         ! phi_y (eta+xi) du/dy; the rest is delegated to  
			 ! the rhs			 
	DO q=1,3         ! all columns 
	   col=elnodes(q)
           ipos = n_num(col)
	  
           entries(q)=(eta+xi)*(dx(i)*dx(q)+dy(i)*dy(q))* &
	              elem_area(elem)*inv_mass  
	   ice_stiff_values(ipos)=ice_stiff_values(ipos)+entries(q)  
	END DO 	
	                 ! division by inv_mass is to be after computing
			 ! contribution from viscosities (otherwise 
			 ! information on m drops out)
			 
	                 ! u,v rhs: 
			 ! u-row: contains -eta(phi_x du/dx+
	                 ! phi_y dv/dx) + (xi/(xi+eta))*entries*u
			 ! -psi_x(xi-eta)div u +psi_x p/2
			 ! v-row: -eta (phi_x du/dy+
			 ! phi_y dv/dy)+
			 ! (xi/(xi+eta))*entries*v 
			 ! -psi_y(xi-eta)div u +psi_y p/2

	rhs_u(row)=rhs_u(row)+(-eta*(dx(i)*eps11(elem)+dy(i)*epsv) &
	   -(xi-eta)*dx(i)*(eps11(elem)+eps22(elem))+dx(i)*pressure)*elem_area(elem) &
	   *inv_mass &
	   +0.8*sum(entries*u_ice(elnodes))+elex
	rhs_v(row)=rhs_v(row)+(-eta*(dx(i)*epsu+dy(i)*eps22(elem)) &
	   -(xi-eta)*dy(i)*(eps11(elem)+eps22(elem))+dy(i)*pressure)*elem_area(elem) &
	   *inv_mass &
	   +0.8*sum(entries*v_ice(elnodes))+eley
	                 ! 0.8=xi/(xi+eta)
	!!!!! Metric terms: They all should go to the rhs because we
	!!!!! want to keep the operator (matrix) same for u and v.

        !rhs_u(row)=rhs_u(row)+metric_factor(elem)*inv_mass*(-(eta+xi)*(dx(i)*ve-epsu/3.0) &
        !   -eta*dy(i)*ue-eta*(epsv+ue*metric_factor(elem))/3.0+xi*epsu/3.0)*elem_area(elem)	 
        
        !rhs_v(row)=rhs_v(row)+metric_factor(elem)*inv_mass*(-(eta+xi)*(-dx(i)*ue-eps11(elem)/3.0 &
        !   -metric_factor(elem)*ve/3.0)-xi*dx(i)*ue+(xi-eta)*(dy(i)*ve+epsv/3.0))*elem_area(elem)	 
        rhs_u(row)=rhs_u(row)+metric_factor(elem)*inv_mass*((eta+xi)*dx(i)*ve &
           -eta*dy(i)*ue-eta*(epsv+epsu+ue*metric_factor(elem))/3.0)*elem_area(elem)
       
        rhs_v(row)=rhs_v(row)+metric_factor(elem)*inv_mass*(-eta*ue*dx(i)+((eta+xi)*(eps11(elem) &
           -metric_factor(elem)*ve)+(xi-eta)*eps22(elem)-pressure/2.0)/3.0)*elem_area(elem)
     
     END DO	       
  END DO 
 
END SUBROUTINE VPmatrix_rhs
!=========================================================================
!
subroutine VPsolve(mesh, new_values)

! Calls PETSc to solve parts of the VP problem. 
! Only owned rows of rhs and matrices are used
! Experimenting with Pmode1 according to the lines below may 
! speed up the performance. The comments there are obtained with 
! matrices of relatively small size and are not necessarily 
! universally valid 

use g_parsup
use mod_mesh
use i_ARRAYS
use ice_VP_interfaces, except_this_one => VPsolve
use g_comm_auto
use iso_c_binding, only: C_INT, C_DOUBLE
implicit none
#include "fparms.h"

logical, save        :: lfirst=.true.
integer(kind=C_INT)  :: ident
integer(kind=C_INT)  :: n3, reuse
integer(kind=C_INT)  :: maxiter, restart, lutype, fillin
real(kind=C_DOUBLE)  :: droptol, soltol
integer :: n
type(t_mesh), intent(in) , target :: mesh
integer,      intent(in)          :: new_values

interface
   subroutine psolver_init(ident, SOL, PCGLOB, PCLOC, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part, rowptr, colind, values, reuse, MPI_COMM) bind(C)
     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, SOL, PCGLOB, PCLOC, lutype, &
                            fillin,  maxiter, restart, &
                            part(*), rowptr(*), colind(*), reuse, MPI_COMM
     real(kind=C_DOUBLE) :: droptol,  soltol, values(*)
   end subroutine psolver_init

   subroutine psolve(ident, ssh_rhs, values, d_eta, newvalues) bind(C)

     use iso_c_binding, only: C_INT, C_DOUBLE
     integer(kind=C_INT) :: ident, newvalues
     real(kind=C_DOUBLE) :: values(*), ssh_rhs(*), d_eta(*)

   end subroutine psolve
end interface

#include "associate_mesh.h"

ident=2
maxiter=2000
restart=15
fillin=3
lutype=2
droptol=1.e-8
soltol=1.e-12

!NR  As the matrix changes so much between the steps, the preconditioner
!NR  must be recomputed. However, also the nonzero-pattern of the matrix
!NR  changes, as non-ice-nodes have zeros on the off-diagonal.
!NR  Best choice is therefor no preconditioner at all: SOLBICGS_NOPC
!NR  The preconditioner settings are (basically) ignored, i.e., the
!NR  preconditioner is set up, but never used.
!NR  Yes, the solution is a bit hacky, but it is only intermediate, hopefully.
!NR
!NR  Remark: The wrapper to parms does employ diagonal scaling to the matrix
!NR  and rhs, which already helps a lot!

if (lfirst) then
    call psolver_init(ident, SOLBICGS_NOPC, PCRAS, PCILUK, lutype, &
        fillin, droptol, maxiter, restart, soltol, &
        part-1, ssh_stiff%rowptr(:)-ssh_stiff%rowptr(1), &
        ssh_stiff%colind-1, ice_stiff_values, 1, MPI_COMM_FESOM)
    lfirst=.false.
 end if
 
 call psolve(ident, rhs_u, ice_stiff_values, u_ice, new_values)
 call psolve(ident, rhs_v, ice_stiff_values, v_ice, 0)  ! values remain constant
   
! Sets boundary conditions to the solution - might leave parms as /= 0
DO n=1,myDim_nod2D
   if((a_ice(n)<0.01).or.(bc_index_nod2D(n)<0.5)) then    
     u_ice(n)=0._8
     v_ice(n)=0._8
   end if   
END DO

end subroutine VPsolve 
! ============================================================================
subroutine VPdynamics(mesh)
use g_parsup
use mod_mesh
use i_PARAM
use i_ARRAYS
use ice_VP_interfaces, except_this_one => VPdynamics
use g_comm_auto
! Driving routine of VP rheology
IMPLICIT NONE
integer        :: n, row, is, ie
logical, save  :: lfirst=.true.
type(t_mesh), intent(in), target  :: mesh
REAL(kind=WP), ALLOCATABLE, DIMENSION(:), SAVE         :: rhs_diag_ice !!! just for testing, remove it afterwards
#include "associate_mesh.h"

 if (lfirst) then
    call VP_init(mesh)
    lfirst=.false.

    allocate(rhs_diag_ice(myDim_nod2D)) !!! just for testing, remove it afterwards
    rhs_diag_ice=0.0                    !!! just for testing, remove it afterwards
 end if


 DO row=1,myDim_nod2D+eDim_nod2D     !! PP n=1,myDim_nod2D+eDim_nod2D       
                                     !! PP row=myList_nod2D(n)
   rhs_m(row)=u_ice(row)
   rhs_a(row)=v_ice(row)
 END DO   
 call VPmatrix_rhs(mesh)
 call VPbc(mesh)
 call exchange_nod(rhs_u,rhs_v)
 call VPsolve(mesh, 1)

 DO row=1,myDim_nod2D+eDim_nod2D     !! PP n=1,myDim_nod2D+eDim_nod2D       
                                     !! PP row=myList_nod2D(n)
   u_ice(row)=0.5*(u_ice(row)+rhs_m(row))
   v_ice(row)=0.5*(v_ice(row)+rhs_a(row))
 END DO  
 ! Previous and predicted states are mixed as in Zhang and Hibler 
 ! 1997
 DO n=1, ice_VP_iter   !(Picard iterations)
  call VPmatrix_rhs(mesh)
  call VPbc(mesh)
  call exchange_nod(rhs_u,rhs_v)
  call VPsolve(mesh, 1)
 END DO  
 DO row=1,myDim_nod2D+eDim_nod2D     !! PP  n=1,myDim_nod2D+eDim_nod2D       
                                     !! PP row=myList_nod2D(n)
   rhs_m(row)=u_ice(row)  ! It has been used to estimate drag coefficient
   rhs_a(row)=v_ice(row)  ! and of-diagonal terms (Coriolis and rotated stress)
 END DO  
 call VPmatrix_rhs(mesh)
 call VPbc(mesh)
 call exchange_nod(rhs_u,rhs_v)
 call VPsolve(mesh, 1)
 call VPmake_implicit(mesh)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  do n=1, myDim_nod2D
!     is=ssh_stiff%rowptr_loc(n)
!     ie=ssh_stiff%rowptr_loc(n+1)-1     
!     rhs_diag_ice(n)=sum(ice_stiff_values(is:ie)*u_ice(ssh_stiff%colind_loc(is:ie)))
!  end do
!
!if (mype==61) then
!     do n=1, myDim_nod2D
!        is=ssh_stiff%rowptr_loc(n)
!        ie=ssh_stiff%rowptr_loc(n+1)-1  
!        write(*,*) rhs_diag_ice(n), rhs_u(n), u_ice(n)
!        write(*,*) n, '***********values*******'
!        write(*,*) ice_stiff_values(is:ie)
!        write(*,*) n, '***********u_ice********'       
!        write(*,*) u_ice(ssh_stiff%colind_loc(is:ie))
!        write(*,*) n, '***********RHS**********'       
!        write(*,*) rhs_u(n)
!        write(*,*) '***************************'       
!     end do
!end if
!call par_ex
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine VPdynamics
! ============================================================================
subroutine VPmake_implicit(mesh)
use g_parsup
use mod_mesh
USE o_ARRAYS
use i_PARAM
use i_ARRAYS
use i_therm_param
use ice_VP_interfaces, except_this_one => VPmake_implicit
IMPLICIT NONE 
real(kind=8)     ::  as, ac, drag, inv_mass, c1, c2, r1, r2, d
integer          ::  n, row


type(t_mesh), intent(in), target  :: mesh

#include "associate_mesh.h"


! updates Coriolis and drag terms to full implicit formulation
! as adviced in Zhang and Hibler 
 ac=cos(theta_io)
 as=sin(theta_io)
 DO row=1,myDim_nod2D+eDim_nod2D    !! PP n=1,myDim_nod2D+eDim_nod2D       
                                    !! row=myList_nod2D(n) 
   if(a_ice(row)<0.01) cycle
   inv_mass=(rhoice*m_ice(row)+rhosno*m_snow(row))/a_ice(row)
   inv_mass=max(inv_mass, 9.0)        ! Limit the weighted mass 
                                      ! if it is too small
   inv_mass=1.0/inv_mass
   drag=sqrt((rhs_m(row)-u_w(row))**2+(rhs_a(row)-v_w(row))**2)* &
         Cd_oce_ice*density_0*inv_mass
   c1=1.0_8/ice_dt+drag*ac
   c2=coriolis_node(row)+drag*as
   
   r1=c1*c1*u_ice(row)+c1*c2*(v_ice(row)-rhs_a(row))+c2*c2*rhs_m(row)
   r2=c1*c1*v_ice(row)-c1*c2*(u_ice(row)-rhs_m(row))+c2*c2*rhs_a(row)  
   d=c1*c1+c2*c2
   u_ice(row)=r1/d
   v_ice(row)=r2/d
  END DO 
end subroutine VPmake_implicit
! ============================================================================
subroutine VPbc(mesh)
use g_parsup
use mod_mesh
use i_PARAM
use i_ARRAYS
use ice_VP_interfaces, except_this_one => VPbc
IMPLICIT NONE 

integer                          :: row, j, is, ie
type(t_mesh), intent(in), target :: mesh

#include "associate_mesh.h"

! Sets boundary conditions to the rhs and operator part.
DO row=1,myDim_nod2D
   if((a_ice(row)<0.01).or.(bc_index_nod2D(row)<0.5)) then

     is=ssh_stiff%rowptr_loc(row)
     ie=ssh_stiff%rowptr_loc(row+1)-1

     !NR By construction, the first entry for each row is the diagonal one.
     !NR We can replace the loop with the if and colind_loc by a simple
     !NR direct statement.
!!$     where (ssh_stiff%colind_loc(is:ie)==row)
!!$           ice_stiff_values(is:ie)=1.0_8
!!$     elsewhere
!!$           ice_stiff_values(is:ie)=0.0_8
!!$     end where
     is = ssh_stiff%rowptr_loc(row)
     ie = ssh_stiff%rowptr_loc(row+1)-1
     ice_stiff_values(is)      = 1.
     ice_stiff_values(is+1:ie) = 0.
    
     rhs_u(row)=0._8
     rhs_v(row)=0._8
   end if   
END DO  

end subroutine VPbc
!=========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  do n=1, myDim_nod2D
!     is=ssh_stiff%rowptr_loc(n)
!     ie=ssh_stiff%rowptr_loc(n+1)-1     
!     rhs_diag_ice(n)=sum(ice_stiff_values(is:ie)*u_ice(ssh_stiff%colind_loc(is:ie)))
!
!  end do
!
!if (mype==147) then
!     do n=1, myDim_nod2D
!        is=ssh_stiff%rowptr_loc(n)
!        ie=ssh_stiff%rowptr_loc(n+1)-1  
!        write(*,*) rhs_diag_ice(n), rhs_u(n), u_ice(n)
!        write(*,*) n, '***********values*******'
!        write(*,*) ice_stiff_values(is:ie)
!        write(*,*) n, '***********u_ice********'       
!        write(*,*) u_ice(ssh_stiff%colind_loc(is:ie))
!        write(*,*) n, '***********RHS**********'       
!        write(*,*) rhs_u(n)
!        write(*,*) '***************************'       
!     end do
!end if
!call par_ex
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

