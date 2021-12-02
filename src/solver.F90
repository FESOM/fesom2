module ssh_solve_preconditioner_interface
    interface
        subroutine ssh_solve_preconditioner(partit, mesh)
        use MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        end subroutine
    end interface
end module

module ssh_solve_cg_interface
    interface
        subroutine ssh_solve_cg(x, rhs, tol, maxiter, partit, mesh)
        use MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(inout), target :: mesh
        real(kind=WP),  intent(inout) :: x(partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(in)    :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),  intent(in)    :: tol
        integer,        intent(in)    :: maxiter
        end subroutine
    end interface
end module
!=========================================================================
subroutine ssh_solve_preconditioner(partit, mesh)
  ! Preconditioner follows MITgcm (JGR, 102,5753-5766, 1997)
  ! If the row r of the ssh equation is a_r eta_r +\sum a_i\eta_i=rhs_row_r
  ! where summation is over all nodes neighboring node r,
  ! the inverse of the preconditioner matrix has the coefficients
  ! 1/a_r, .... -2*a_i/a_r/(a_r+(a_diag)_i) ....
  ! Here (a_diag)_i is the diagonal value in row i of the ssh matrix.

  ! The  inverse of preconditioner matrix (M^{-1} in general notation and K in the
  ! paper cited) is, in reality, one iteration of the
  ! Jacobi method, with symmetrization. We need symmetrization to be able to use
  ! the conjugate gradient method.    
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE g_comm_auto
    IMPLICIT NONE
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    integer                      :: nend, row, node, n, offset
    real(kind=WP), allocatable   :: diag_values(:)
    real(kind=WP), pointer       :: pr_values(:)
    integer,       pointer       :: rptr(:), cind(:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    nend=ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)
    allocate(mesh%ssh_stiff%pr_values(nend))     ! Will store the values of inverse preconditioner matrix 
    pr_values=>mesh%ssh_stiff%pr_values
    cind     =>mesh%ssh_stiff%colind_loc
    rptr     =>mesh%ssh_stiff%rowptr_loc
    allocate(diag_values(myDim_nod2D+eDim_nod2D)) ! Temporary, will be thrown away

    DO row=1, myDim_nod2D
       offset=ssh_stiff%rowptr(row)- ssh_stiff%rowptr(1)+1
       diag_values(row)=ssh_stiff%values(offset)
    END DO
    call exchange_nod(diag_values, partit)              ! We have diagonal values
    ! ==========
    ! Fill in the preconditioner
    ! ==========
    DO row=1, myDim_nod2D
       offset=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)
       nend=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(row)
       pr_values(offset+1)=1.0_WP/ssh_stiff%values(offset+1)
       DO n=2, nend   
          node=cind(offset+n)    ! Will be ssh_stiff$colind(offset+n) 
          pr_values(n+offset)=-0.5_WP*(ssh_stiff%values(n+offset)/ssh_stiff%values(1+offset))/  &
                               (ssh_stiff%values(1+offset)+ diag_values(node)) 
       END DO
    END DO
   deallocate(diag_values)
end subroutine ssh_solve_preconditioner

! ========================================================
subroutine ssh_solve_cg(x,rhs, tol, maxiter, partit, mesh)
  ! Conjugate gradient solver
  ! Our ssh matrix is symmetric,  because we  compute divergencethe contributions as
  ! integrated over area of scalar control volume.
  ! 
  ! I tried first to follow the MITgcm paper, but I have doubts about
  ! their computations of beta. The variant below -- see Wikipedia.
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE g_comm_auto
  IMPLICIT NONE
  type(t_partit),intent(inout), target :: partit
  type(t_mesh),  intent(inout), target :: mesh
  real(kind=WP), intent(inout) :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), intent(in)    :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP), intent(in)    :: tol
  integer, intent(in)          :: maxiter
  integer                      :: row, nini, nend, iter 
  real(kind=WP)                :: sprod(2), s_old, s_aux, al, be, rtol
  integer                      :: req
  real(kind=WP), allocatable   :: rr(:), zz(:), pp(:), App(:)
  real(kind=WP), pointer       :: pr_values(:)
  integer,       pointer       :: rptr(:), cind(:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  pr_values=>mesh%ssh_stiff%pr_values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc

  row=myDim_nod2D+eDim_nod2D
  allocate(rr(row), zz(row), pp(row), App(row))
  rr=0.0_WP
  zz=0.0_WP
  pp=0.0_WP
  App=0.0_WP
  
  ! ============== 
  ! Initialization. We solve AX=b, r_0=b-AX_0
  ! ============== 
  ! Define working tolerance: 
  ! ============== 
  s_old=sum(rhs(1:myDim_nod2D)*rhs(1:myDim_nod2D))
  call MPI_Iallreduce(MPI_IN_PLACE, s_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, req, MPIerr)
  call MPI_Wait(req, MPI_STATUS_IGNORE, MPIerr)
  rtol=tol*sqrt(s_old/real(nod2D,WP))
  ! ==============
  ! Compute r0
  ! ==============
  DO row=1, myDim_nod2D
     !nini=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
     !nend=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
     !rr(row)=rhs(row)-sum(ssh_stiff%values(nini:nend)*X(cind(nini:nend)))
     rr(row)=rhs(row)-sum(ssh_stiff%values(rptr(row):rptr(row+1)-1)* &
                      X(cind(rptr(row):rptr(row+1)-1)))
  END DO
  call exchange_nod(rr, partit)
  ! =============
  ! z_0=M^{-1} r_0  (M^{-1} is the precondit. matrix)
  ! pp is the search direction
  ! =============
  DO row=1, myDim_nod2D
     zz(row)= sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
     pp(row)=zz(row)
  END DO
  ! ===============
  ! Scalar product of r*z
  ! ===============
  s_old=sum(rr(1:myDim_nod2D)*zz(1:myDim_nod2D))
  call  MPI_Iallreduce(MPI_IN_PLACE, s_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, req, MPIerr)
  call MPI_Wait(req, MPI_STATUS_IGNORE, MPIerr)
  
  ! ===============
  ! Iterations
  ! ===============
  Do iter=1, maxiter
     ! ============
     ! Compute Ap
     ! ============
  call exchange_nod(pp, partit)     !  Update before matrix-vector multiplications

  DO row=1, myDim_nod2D
     App(row)=sum(ssh_stiff%values(rptr(row):rptr(row+1)-1)*pp(cind(rptr(row):rptr(row+1)-1)))
  END DO
     ! ============
     ! Scalar products for alpha
     ! ============
 
  s_aux=sum(pp(1:myDim_nod2D)*App(1:myDim_nod2D))
 
  call  MPI_Iallreduce(MPI_IN_PLACE, s_aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, req, MPIerr)
  call MPI_Wait(req, MPI_STATUS_IGNORE, MPIerr)
  al=s_old/s_aux
     ! ===========
     ! New X and residual r
     ! ===========
  DO row=1, myDim_nod2D
     X(row)=X(row)+al*pp(row)
     rr(row)=rr(row)-al*App(row)
  END DO
     ! =========== 
     ! New z
     ! ===========
  call exchange_nod(rr, partit)     ! Update before matrix-vector multiplications
 
  DO row=1, myDim_nod2D
     zz(row)= sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
  END DO

     ! ===========
     ! Scalar products for beta
     ! ===========
  sprod(1)=sum(rr(1:myDim_nod2D)*zz(1:myDim_nod2D))
  sprod(2)=sum(rr(1:myDim_nod2D)*rr(1:myDim_nod2D))
  
  call MPI_Iallreduce(MPI_IN_PLACE, sprod, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, req, MPIerr)
  call MPI_Wait(req, MPI_STATUS_IGNORE, MPIerr)
     ! ===========
     ! Exit if tolerance is achieved
     ! ===========
  if(sqrt(sprod(2)/nod2D)< rtol) then
  !write(*,*) mype, 'exit', iter
  exit
  endif
  be=sprod(1)/s_old
  s_old=sprod(1)
     ! ===========
     ! New p
     ! ===========
  DO row=1,myDim_nod2D
  pp(row)=zz(row)+be*pp(row)
  end do
  END DO
  deallocate(App,zz,pp,rr)

 ! At the end: The result is in X, but it needs a halo exchange.
  call exchange_nod(x, partit)
end subroutine ssh_solve_cg 

! ===================================================================

