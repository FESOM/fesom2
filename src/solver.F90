module ssh_solve_preconditioner_interface
    interface
        subroutine ssh_solve_preconditioner(solverinfo, partit, mesh)
        use MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_solverinfo),  intent(inout), target :: solverinfo
        type(t_partit),      intent(inout), target :: partit
        type(t_mesh),        intent(inout), target :: mesh
        end subroutine ssh_solve_preconditioner
    end interface
end module ssh_solve_preconditioner_interface

module ssh_solve_cg_interface
    interface
        subroutine ssh_solve_cg(x, rhs, solverinfo, partit, mesh)
        use MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_DYN
        type(t_solverinfo),  intent(inout), target :: solverinfo
        type(t_partit),      intent(inout), target :: partit
        type(t_mesh),        intent(inout), target :: mesh
        real(kind=WP),       intent(inout) :: x(partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP),       intent(in)    :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
        end subroutine ssh_solve_cg
    end interface
end module ssh_solve_cg_interface
!=========================================================================
subroutine ssh_solve_preconditioner(solverinfo, partit, mesh)
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
    type(t_solverinfo),  intent(inout), target :: solverinfo
    type(t_partit),      intent(inout), target :: partit
    type(t_mesh),        intent(inout), target :: mesh
    integer                      :: nend, row, node, n, offset
    real(kind=WP), allocatable   :: diag_values(:)
    real(kind=WP), pointer       :: pr_values(:)
    integer,       pointer       :: rptr(:), cind(:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    nend=ssh_stiff%rowptr(myDim_nod2D+1)-ssh_stiff%rowptr(1)
    ! Guarded: building the preconditioner twice (a restart path can do this)
    ! otherwise aborts with "allocatable array is already allocated".
    if (.not. allocated(mesh%ssh_stiff%pr_values)) then
        allocate(mesh%ssh_stiff%pr_values(nend))  ! values of the inverse preconditioner matrix
    end if
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
          if (solverinfo%precond_variant == 0) then
             ! As coded since 60b46bdc. Note two things about it. The comment at
             ! the top of this routine specifies -2*a_i/a_r/(a_r+(a_diag)_i),
             ! which is this expression times 4. And entry (r,i) here is
             ! -c*a_ri/(a_rr*(a_rr+a_ii)) while its transpose (i,r) is
             ! -c*a_ir/(a_ii*(a_ii+a_rr)); with A symmetric those agree only
             ! where a_rr == a_ii, i.e. on a uniform mesh. CG needs an SPD M^-1,
             ! so this is worth measuring rather than assuming -- see the r.z<0
             ! sensor in ssh_solve_cg.
             pr_values(n+offset)=-0.5_WP*(ssh_stiff%values(n+offset)/ssh_stiff%values(1+offset))/  &
                                  (ssh_stiff%values(1+offset)+ diag_values(node)) 
          else
             ! Textbook symmetrised Jacobi, M^-1 = 2D^-1 - D^-1 A D^-1. The
             ! off-diagonal -a_ri/(a_rr*a_ii) is symmetric by construction, and
             ! its diagonal reduces to 1/a_rr, which is what is coded above.
             pr_values(n+offset)=-ssh_stiff%values(n+offset)/                                     &
                                  (ssh_stiff%values(1+offset)*diag_values(node))
          end if
       END DO
    END DO
   deallocate(diag_values)

   n=myDim_nod2D+eDim_nod2D
   if (.not. allocated(solverinfo%rr)) then
      allocate(solverinfo%rr(n), solverinfo%zz(n), solverinfo%pp(n), solverinfo%App(n))
   end if
   solverinfo%rr =0.0_WP
   solverinfo%zz =0.0_WP
   solverinfo%pp =0.0_WP
   solverinfo%App=0.0_WP
end subroutine ssh_solve_preconditioner

! ========================================================
subroutine ssh_solve_cg(x, rhs, solverinfo, partit, mesh)
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
#if defined(FESOM_PROFILING)
  USE fesom_rtdiagnostics, only: diag_count
#endif
  IMPLICIT NONE
  type(t_solverinfo),  intent(inout), target :: solverinfo
  type(t_partit),      intent(inout), target :: partit
  type(t_mesh),        intent(inout), target :: mesh
  real(kind=WP),       intent(inout)         :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),       intent(in)            :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  integer                      :: row, nini, nend, iter
  real(kind=WP)                :: sprod(2), s_old, s_aux, al, be, rtol
  integer                      :: req
  logical                      :: converged
  real(kind=WP), pointer       :: pr_values(:), rr(:), zz(:), pp(:), App(:)
  integer,       pointer       :: rptr(:), cind(:)


#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  pr_values=>mesh%ssh_stiff%pr_values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc

  rr =>solverinfo%rr
  zz =>solverinfo%zz
  pp =>solverinfo%pp
  App=>solverinfo%App
 
  ! ============== 
  ! Initialization. We solve AX=b, r_0=b-AX_0
  ! ============== 
  ! Define working tolerance: 
  ! ==============
#if !defined(__openmp_reproducible)
  s_old=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) REDUCTION(+:s_old)
  DO row=1, myDim_nod2D
     s_old=s_old+rhs(row)*rhs(row)
  END DO
!$OMP END PARALLEL DO
#else
 s_old = sum(rhs(1:myDim_nod2D) * rhs(1:myDim_nod2D))
#endif

  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)
  rtol=solverinfo%soltol*sqrt(s_old/real(nod2D,WP))
  ! ==============
  ! Compute r0
  ! ==============
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, myDim_nod2D
     rr(row)=rhs(row)-sum(ssh_stiff%values(rptr(row):rptr(row+1)-1)* &
                      X(cind(rptr(row):rptr(row+1)-1)))
  END DO
!$OMP END PARALLEL DO 
  call exchange_nod(rr, partit)
!$OMP BARRIER
  ! =============
  ! z_0=M^{-1} r_0  (M^{-1} is the precondit. matrix)
  ! pp is the search direction
  ! =============
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, myDim_nod2D
     zz(row)= sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
     pp(row)=zz(row)
  END DO
!$OMP END PARALLEL DO 
  ! ===============
  ! Scalar product of r*z
  ! ===============

#if !defined(__openmp_reproducible)
  s_old=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) REDUCTION(+:s_old)
  DO row=1, myDim_nod2D
     s_old=s_old+rr(row)*zz(row)
  END DO
!$OMP END PARALLEL DO
#else
  s_old = sum(rr(1:myDim_nod2D) * zz(1:myDim_nod2D))
#endif

  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)
  
  ! ===============
  ! Iterations
  ! ===============
  converged = .false.
  Do iter=1, solverinfo%maxiter
     ! ============
     ! Compute Ap
     ! ============
  call exchange_nod(pp, partit)     !  Update before matrix-vector multiplications
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, myDim_nod2D
     App(row)=sum(ssh_stiff%values(rptr(row):rptr(row+1)-1)*pp(cind(rptr(row):rptr(row+1)-1)))
  END DO
!$OMP END PARALLEL DO 
     ! ============
     ! Scalar products for alpha
     ! ============
 
#if !defined(__openmp_reproducible)
  s_aux=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) REDUCTION(+:s_aux)
  DO row=1, myDim_nod2D
     s_aux=s_aux+pp(row)*App(row)
  END DO
!$OMP END PARALLEL DO
#else
 s_aux = sum(pp(1:myDim_nod2D) * App(1:myDim_nod2D))
#endif

  call MPI_Allreduce(MPI_IN_PLACE, s_aux, 1, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)

     ! ===========
     ! Breakdown guard. An equivalent check existed until 4eb2f21d ("Removed
     ! debug output in solver."). p.Ap must be positive for an SPD system; if it
     ! reaches zero or NaN then al = s_old/s_aux poisons the entire d_eta vector
     ! with Inf/NaN, and that surfaces several stages downstream -- typically as
     ! a NaN in Wvel -- with nothing pointing back here. Keep the current X and
     ! leave the loop; do not abort, and do not mark the solve converged.
     ! Evaluated on the already-reduced value, so every rank branches alike.
     ! ===========
  if (s_aux <= 0.0_WP .or. s_aux /= s_aux) then
     solverinfo%nbreakdown = solverinfo%nbreakdown + 1
     if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
        write(*,*) 'WARNING: ssh CG breakdown at iter ', iter, ': p.Ap= ', s_aux, &
                   ' is not positive; keeping current solution and leaving the loop'
     endif
     exit
  endif

  al=s_old/s_aux
     ! ===========
     ! New X and residual r
     ! ===========
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  DO row=1, myDim_nod2D
     X(row) =X(row) +al* pp(row)
     rr(row)=rr(row)-al*App(row)
  END DO
!$OMP END PARALLEL DO
     ! =========== 
     ! New z
     ! ===========
  call exchange_nod(rr, partit)     ! Update before matrix-vector multiplications
!$OMP BARRIER
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  DO row=1, myDim_nod2D
     zz(row)= sum(pr_values(rptr(row):rptr(row+1)-1)*rr(cind(rptr(row):rptr(row+1)-1)))
  END DO
!$OMP END PARALLEL DO
     ! ===========
     ! Scalar products for beta
     ! ===========
#if !defined(__openmp_reproducible)
sprod(1:2)=0.0_WP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) REDUCTION(+:sprod)
  DO row=1, myDim_nod2D
     sprod(1)=sprod(1)+rr(row)*zz(row)
     sprod(2)=sprod(2)+rr(row)*rr(row)
  END DO
!$OMP END PARALLEL DO
#else
    sprod(1) = sum(rr(1:myDim_nod2D) * zz(1:myDim_nod2D))
    sprod(2) = sum(rr(1:myDim_nod2D) * rr(1:myDim_nod2D))
#endif
  
  call MPI_Allreduce(MPI_IN_PLACE, sprod, 2, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, MPIerr)

!$OMP BARRIER
     ! ===========
     ! Exit if tolerance is achieved
     ! ===========
  if (sqrt(sprod(2)/nod2D)< rtol) then
     converged = .true.
     exit
  endif
     ! ===========
     ! SPD sensor. r.z cannot be negative for an SPD preconditioner, so any count
     ! here is direct evidence that M^-1 is not SPD. Cheap, and evaluated on the
     ! already-reduced value so all ranks agree.
     ! ===========
  if (sprod(1) < 0.0_WP) then
     solverinfo%nnegrz = solverinfo%nnegrz + 1
     if (partit%mype==0 .and. solverinfo%nnegrz == 1) then
        write(*,*) 'WARNING: ssh CG r.z < 0 at iter ', iter, ': r.z= ', sprod(1), &
                   ' -- preconditioner is not symmetric positive definite'
     endif
  endif

  if (s_old == 0.0_WP) then
     solverinfo%nbreakdown = solverinfo%nbreakdown + 1
     if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
        write(*,*) 'WARNING: ssh CG breakdown at iter ', iter, ': previous r.z is zero'
     endif
     exit
  endif
  be=sprod(1)/s_old
  s_old=sprod(1)
     ! ===========
     ! New p
     ! ===========
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  DO row=1,myDim_nod2D
     pp(row)=zz(row)+be*pp(row)
  END DO
!$OMP END PARALLEL DO
  END DO

  ! Record SSH-CG iteration count (and non-convergence) for fesom.stats.
#if defined(FESOM_PROFILING)
  call diag_count("ssh_cg_iters", min(iter, solverinfo%maxiter))
  if (.not. converged) call diag_count("ssh_cg_nonconv", 1)
#endif

  !_____________________________________________________________________________
  ! Always-on counters. FESOM_PROFILING defaults OFF, so the block above is
  ! invisible in a normal run; these feed the per-step line in write_step_info
  ! and the end-of-run summary, which are not gated on anything.
  ! No MPI needed: every rank exits on the same iteration, because the exit test
  ! above is on the globally reduced sprod.
  solverinfo%iters_last = min(iter, solverinfo%maxiter)
  solverinfo%iters_max  = max(solverinfo%iters_max, solverinfo%iters_last)
  solverinfo%iters_sum  = solverinfo%iters_sum + solverinfo%iters_last
  solverinfo%nsolves    = solverinfo%nsolves + 1
  solverinfo%resid_last = sqrt(sprod(2)/nod2D)
  solverinfo%rtol_last  = rtol
  if (.not. converged) solverinfo%nonconv = solverinfo%nonconv + 1

  !_____________________________________________________________________________
  ! Safety-net log line: a stall at maxiter used to be completely silent. Unlike
  ! cfl_z, which has check_blowup as a hard backstop, non-convergence has none --
  ! we deliberately do not abort, because a stall on step 1 is legitimate in some
  ! configurations, so this line is the only thing between a silently
  ! under-converged run and nobody noticing. Rate-limited because a stalled
  ! solver usually stalls every step, and one line per step would bury the
  ! surrounding diagnostics.
  if (.not. converged .and. partit%mype==0) then
     if (solverinfo%nonconv <= 5 .or. mod(solverinfo%nonconv, 100) == 0) then
        write(*,*) 'WARNING: ssh CG did not converge in ', solverinfo%maxiter, &
                   ' iters; rms(resid)=', solverinfo%resid_last,               &
                   ' target rtol=', rtol, ' (occurrence ', solverinfo%nonconv, ')'
     endif
  endif

 ! At the end: The result is in X, but it needs a halo exchange.
  call exchange_nod(x, partit)
!$OMP BARRIER
end subroutine ssh_solve_cg

! ===================================================================

