module solver_module
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_DYN
    USE g_comm_auto
    USE halo_timing_mod, only: ht_on, ht_sum_mark, ht_sum_done
#if defined(FESOM_PROFILING)
  USE fesom_rtdiagnostics, only: diag_count
#endif

    implicit none

    private
    public :: ssh_solve_preconditioner, ssh_solve_cg

    ! Chebyshev preconditioner (precond_variant = 2): polynomial in D^-1 A on [cheb_a, cheb_b]
    integer,            parameter        :: cheb_powit = 80
    real(kind=WP_full), parameter        :: cheb_margin = 1.1_WP_full, cheb_ratio = 30.0_WP_full
    integer,            save             :: cheb_deg = 5
    real(kind=WP_full), save             :: cheb_a = 0.0_WP_full, cheb_b = 0.0_WP_full
    real(kind=WP), allocatable, save     :: cheb_dinv(:), cheb_d(:), cheb_res(:)
    ! one-reduction CG (cg_variant = 1): s = A p
    real(kind=WP), allocatable, save     :: cg1_s(:)

contains

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
    IMPLICIT NONE
    type(t_solverinfo),  intent(inout), target :: solverinfo
    type(t_partit),      intent(inout), target :: partit
    type(t_mesh),        intent(inout), target :: mesh
    integer                      :: nend, row, node, n, offset
    real(kind=WP), allocatable   :: diag_values(:)
    real(kind=WP), pointer       :: pr_values(:)
    integer,       pointer       :: rptr(:), cind(:)


    nend=mesh%ssh_stiff%rowptr(partit%myDim_nod2D+1)-mesh%ssh_stiff%rowptr(1)
    ! Guarded: building the preconditioner twice (a restart path can do this)
    ! otherwise aborts with "allocatable array is already allocated".
    if (.not. allocated(mesh%ssh_stiff%pr_values)) then
        allocate(mesh%ssh_stiff%pr_values(nend))  ! values of the inverse preconditioner matrix
    end if
    pr_values=>mesh%ssh_stiff%pr_values
    cind     =>mesh%ssh_stiff%colind_loc
    rptr     =>mesh%ssh_stiff%rowptr_loc
    allocate(diag_values(partit%myDim_nod2D+partit%eDim_nod2D)) ! Temporary, will be thrown away

    DO row=1, partit%myDim_nod2D
       offset=mesh%ssh_stiff%rowptr(row)- mesh%ssh_stiff%rowptr(1)+1
       diag_values(row)=mesh%ssh_stiff%values(offset)
    END DO
    call exchange_nod(diag_values, partit)              ! We have diagonal values
    ! ==========
    ! Fill in the preconditioner
    ! ==========
    DO row=1, partit%myDim_nod2D
       offset=mesh%ssh_stiff%rowptr(row)-mesh%ssh_stiff%rowptr(1)
       nend=mesh%ssh_stiff%rowptr(row+1)-mesh%ssh_stiff%rowptr(row)
       pr_values(offset+1)=1.0_WP/mesh%ssh_stiff%values(offset+1)
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
             pr_values(n+offset)=-0.5_WP*(mesh%ssh_stiff%values(n+offset)/mesh%ssh_stiff%values(1+offset))/  &
                                  (mesh%ssh_stiff%values(1+offset)+ diag_values(node)) 
          else
             ! Textbook symmetrised Jacobi, M^-1 = 2D^-1 - D^-1 A D^-1. The
             ! off-diagonal -a_ri/(a_rr*a_ii) is symmetric by construction, and
             ! its diagonal reduces to 1/a_rr, which is what is coded above.
             pr_values(n+offset)=-mesh%ssh_stiff%values(n+offset)/                                     &
                                  (mesh%ssh_stiff%values(1+offset)*diag_values(node))
          end if
       END DO
    END DO
   deallocate(diag_values)

   n=partit%myDim_nod2D+partit%eDim_nod2D
   if (.not. allocated(solverinfo%rr)) then
      allocate(solverinfo%rr(n), solverinfo%zz(n), solverinfo%pp(n), solverinfo%App(n))
   end if
   solverinfo%rr =0.0_WP
   solverinfo%zz =0.0_WP
   solverinfo%pp =0.0_WP
   solverinfo%App=0.0_WP
   if (solverinfo%precond_variant == 2) call cheb_setup(solverinfo, partit, mesh)
end subroutine ssh_solve_preconditioner

! ========================================================
subroutine ssh_solve_cg(x, rhs, solverinfo, partit, mesh)
  ! Conjugate gradient solver
  ! Our ssh matrix is symmetric,  because we  compute divergencethe contributions as
  ! integrated over area of scalar control volume.
  ! 
  ! I tried first to follow the MITgcm paper, but I have doubts about
  ! their computations of beta. The variant below -- see Wikipedia.
  IMPLICIT NONE
  type(t_solverinfo),  intent(inout), target :: solverinfo
  type(t_partit),      intent(inout), target :: partit
  type(t_mesh),        intent(inout), target :: mesh
  real(kind=WP),       intent(inout)         :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),       intent(in)            :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  integer                      :: row, nini, nend, iter
  ! CG scalars in WP_full. Every one of these is either a global reduction of a
  ! sum over all nod2D, or derived from one, and they set the search direction and
  ! the stopping test -- so rounding here steers the iteration itself rather than
  ! just reporting on it. The vectors (rr/zz/pp/App) and the matrix stay WP: the
  ! bandwidth is theirs, the accuracy is these.
  ! The MPI_Allreduce calls below pass MPI_WP_FULL to match; declaring these
  ! WP_full makes the buffer type and the declared type agree by construction
  ! rather than by WP happening to be real64.
  real(kind=WP_full)           :: sprod(2), s_old, s_aux, al, be, rtol
  integer                      :: req
  logical                      :: converged
  real(kind=WP), pointer       :: pr_values(:), rr(:), zz(:), pp(:), App(:)
  integer,       pointer       :: rptr(:), cind(:)


  
  pr_values=>mesh%ssh_stiff%pr_values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc

  rr =>solverinfo%rr
  zz =>solverinfo%zz
  pp =>solverinfo%pp
  App=>solverinfo%App
#ifdef ENABLE_OPENACC
  ! the same iteration on the device (ssh_cg_iterate_acc); the bookkeeping below is shared
  if (solverinfo%cg_variant == 1) then
     call ssh_cg1_iterate_acc(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  else
     call ssh_cg_iterate_acc(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  end if
#else
  if (solverinfo%cg_variant == 1) then
     call ssh_cg1_iterate(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  else
 
  ! ============== 
  ! Initialization. We solve AX=b, r_0=b-AX_0
  ! ============== 
  ! Define working tolerance: 
  ! ==============
 s_old = sum(real(rhs(1:partit%myDim_nod2D), WP_full) * real(rhs(1:partit%myDim_nod2D), WP_full))

  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()
  rtol=solverinfo%soltol*sqrt(s_old/real(mesh%nod2D,WP_full))
  ! ==============
  ! Compute r0
  ! ==============
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, partit%myDim_nod2D
     rr(row)=real(real(rhs(row),WP_full)                                        &
              - sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)*     &
                    real(X(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  END DO
!$OMP END PARALLEL DO 
  ! =============
  ! z_0=M^{-1} r_0  (M^{-1} is the precondit. matrix)
  ! pp is the search direction
  ! =============
  if (solverinfo%precond_variant == 2) then
     call cheb_apply(rr, zz, partit, mesh)
  else
  call exchange_nod(rr, partit)
!$OMP BARRIER
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, partit%myDim_nod2D
     zz(row)= real(sum(real(pr_values(rptr(row):rptr(row+1)-1),WP_full)*       &
                       real(rr(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  END DO
!$OMP END PARALLEL DO 
  end if
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  DO row=1, partit%myDim_nod2D
     pp(row)=zz(row)
  END DO
!$OMP END PARALLEL DO
  ! ===============
  ! Scalar product of r*z
  ! ===============

  s_old = sum(real(rr(1:partit%myDim_nod2D), WP_full) * real(zz(1:partit%myDim_nod2D), WP_full))

  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()
  
  ! ===============
  ! Iterations
  ! ===============
  converged = .false.
  sprod     = 0.0_WP_full ! defined also if the loop is left before the first residual
  Do iter=1, solverinfo%maxiter
     ! ============
     ! Compute Ap
     ! ============
  call exchange_nod(pp, partit)     !  Update before matrix-vector multiplications
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row) 
  DO row=1, partit%myDim_nod2D
     App(row)=real(sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                        real(pp(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  END DO
!$OMP END PARALLEL DO 
     ! ============
     ! Scalar products for alpha
     ! ============
 
 s_aux = sum(real(pp(1:partit%myDim_nod2D), WP_full) * real(App(1:partit%myDim_nod2D), WP_full))

  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, s_aux, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()

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
  DO row=1, partit%myDim_nod2D
     X(row) =X(row) +al* pp(row)
     rr(row)=rr(row)-al*App(row)
  END DO
!$OMP END PARALLEL DO
     ! =========== 
     ! New z
     ! ===========
  if (solverinfo%precond_variant == 2) then
     call cheb_apply(rr, zz, partit, mesh)
  else
  call exchange_nod(rr, partit)     ! Update before matrix-vector multiplications
!$OMP BARRIER
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  DO row=1, partit%myDim_nod2D
     zz(row)= real(sum(real(pr_values(rptr(row):rptr(row+1)-1),WP_full)*       &
                       real(rr(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  END DO
!$OMP END PARALLEL DO
  end if
     ! ===========
     ! Scalar products for beta
     ! ===========
    sprod(1) = sum(real(rr(1:partit%myDim_nod2D), WP_full) * real(zz(1:partit%myDim_nod2D), WP_full))
    sprod(2) = sum(real(rr(1:partit%myDim_nod2D), WP_full) * real(rr(1:partit%myDim_nod2D), WP_full))
  
  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, sprod, 2, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()

!$OMP BARRIER
     ! ===========
     ! Exit if tolerance is achieved
     ! ===========
  if (sqrt(sprod(2)/mesh%nod2D)< rtol) then
     converged = .true.
     exit
  endif
     ! ===========
     ! SPD sensor. r.z cannot be negative for an SPD preconditioner, so any count
     ! here is direct evidence that M^-1 is not SPD. Cheap, and evaluated on the
     ! already-reduced value so all ranks agree.
     ! ===========
  if (sprod(1) < 0.0_WP_full) then
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
  DO row=1,partit%myDim_nod2D
     pp(row)=zz(row)+be*pp(row)
  END DO
!$OMP END PARALLEL DO
  END DO
  end if
#endif

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
  solverinfo%resid_last = sqrt(sprod(2)/mesh%nod2D)
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

#ifdef ENABLE_OPENACC
!=========================================================================
! Initialisation and iterations of ssh_solve_cg on the device: explicit row products, scalar products
! as device reductions in WP_full (summation order differs from the host); ssh_stiff%values is copied every step.
subroutine ssh_cg_iterate_acc(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  IMPLICIT NONE
  type(t_solverinfo),  intent(inout), target :: solverinfo
  type(t_partit),      intent(inout), target :: partit
  type(t_mesh),        intent(inout), target :: mesh
  real(kind=WP),       intent(inout)         :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),       intent(in)            :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  integer,             intent(out)           :: iter
  logical,             intent(out)           :: converged
  real(kind=WP_full),  intent(out)           :: sprod(2), rtol
  real(kind=WP_full)           :: s_old, s_aux, al, be, acc, sp1, sp2
  integer                      :: row, j, nrow
  logical                      :: cheb
  logical, save                :: on_device = .false.
  real(kind=WP), pointer       :: pr_values(:), values(:), rr(:), zz(:), pp(:), App(:)
  integer,       pointer       :: rptr(:), cind(:)
  pr_values=>mesh%ssh_stiff%pr_values
  values   =>mesh%ssh_stiff%values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc
  rr =>solverinfo%rr
  zz =>solverinfo%zz
  pp =>solverinfo%pp
  App=>solverinfo%App
  nrow=partit%myDim_nod2D
  cheb=(solverinfo%precond_variant == 2)
  if (.not. on_device) then
     !$ACC ENTER DATA COPYIN(cind, rptr, pr_values) CREATE(values, rr, zz, pp, App, x, rhs)
     on_device = .true.
  end if
  !$ACC UPDATE DEVICE(values, x, rhs)

  ! working tolerance
  s_old = 0.0_WP_full
  !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+:s_old)
  do row=1, nrow
     s_old = s_old + real(rhs(row),WP_full)*real(rhs(row),WP_full)
  end do
  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()
  rtol=solverinfo%soltol*sqrt(s_old/real(mesh%nod2D,WP_full))

  ! r0 = b - A x0 (one thread per row: the rows have a few entries only)
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j)
  do row=1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j=rptr(row), rptr(row+1)-1
        acc = acc + real(values(j),WP_full)*real(x(cind(j)),WP_full)
     end do
     rr(row)=real(real(rhs(row),WP_full)-acc, WP)
  end do

  ! z0 = M^-1 r0, p0 = z0, r.z
  s_old = 0.0_WP_full
  if (cheb) then
     call cheb_apply_acc(rr, zz, values, rptr, cind, nrow, partit)
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+:s_old)
     do row=1, nrow
        pp(row)=zz(row)
        s_old = s_old + real(rr(row),WP_full)*real(zz(row),WP_full)
     end do
  else
  call exchange_nod(rr, partit, luse_g2g=.true.)
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j) REDUCTION(+:s_old)
  do row=1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j=rptr(row), rptr(row+1)-1
        acc = acc + real(pr_values(j),WP_full)*real(rr(cind(j)),WP_full)
     end do
     zz(row)=real(acc, WP)
     pp(row)=zz(row)
     s_old = s_old + real(rr(row),WP_full)*real(zz(row),WP_full)
  end do
  end if
  if (ht_on) call ht_sum_mark()
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  if (ht_on) call ht_sum_done()

  converged = .false.
  sprod = 0.0_WP_full
  do iter=1, solverinfo%maxiter
     ! A p and p.Ap
     call exchange_nod(pp, partit, luse_g2g=.true.)
     s_aux = 0.0_WP_full
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j) REDUCTION(+:s_aux)
     do row=1, nrow
        acc = 0.0_WP_full
        !$ACC LOOP SEQ
        do j=rptr(row), rptr(row+1)-1
           acc = acc + real(values(j),WP_full)*real(pp(cind(j)),WP_full)
        end do
        App(row)=real(acc, WP)
        s_aux = s_aux + real(pp(row),WP_full)*real(App(row),WP_full)
     end do
     if (ht_on) call ht_sum_mark()
     call MPI_Allreduce(MPI_IN_PLACE, s_aux, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
     if (ht_on) call ht_sum_done()
     if (s_aux <= 0.0_WP .or. s_aux /= s_aux) then
        solverinfo%nbreakdown = solverinfo%nbreakdown + 1
        if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
           write(*,*) 'WARNING: ssh CG breakdown at iter ', iter, ': p.Ap= ', s_aux, &
                      ' is not positive; keeping current solution and leaving the loop'
        endif
        exit
     endif
     al=s_old/s_aux
     !$ACC PARALLEL LOOP DEFAULT(PRESENT)
     do row=1, nrow
        x(row) =x(row) +al* pp(row)
        rr(row)=rr(row)-al*App(row)
     end do
     ! z = M^-1 r, r.z and r.r
     sp1 = 0.0_WP_full
     sp2 = 0.0_WP_full
     if (cheb) then
        call cheb_apply_acc(rr, zz, values, rptr, cind, nrow, partit)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+:sp1, sp2)
        do row=1, nrow
           sp1 = sp1 + real(rr(row),WP_full)*real(zz(row),WP_full)
           sp2 = sp2 + real(rr(row),WP_full)*real(rr(row),WP_full)
        end do
     else
     call exchange_nod(rr, partit, luse_g2g=.true.)
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j) REDUCTION(+:sp1, sp2)
     do row=1, nrow
        acc = 0.0_WP_full
        !$ACC LOOP SEQ
        do j=rptr(row), rptr(row+1)-1
           acc = acc + real(pr_values(j),WP_full)*real(rr(cind(j)),WP_full)
        end do
        zz(row)=real(acc, WP)
        sp1 = sp1 + real(rr(row),WP_full)*real(zz(row),WP_full)
        sp2 = sp2 + real(rr(row),WP_full)*real(rr(row),WP_full)
     end do
     end if
     sprod(1)=sp1
     sprod(2)=sp2
     if (ht_on) call ht_sum_mark()
     call MPI_Allreduce(MPI_IN_PLACE, sprod, 2, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
     if (ht_on) call ht_sum_done()
     if (sqrt(sprod(2)/mesh%nod2D)< rtol) then
        converged = .true.
        exit
     endif
     if (sprod(1) < 0.0_WP_full) then
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
     !$ACC PARALLEL LOOP DEFAULT(PRESENT)
     do row=1, nrow
        pp(row)=zz(row)+be*pp(row)
     end do
  end do
  ! the halo of x is exchanged by the caller (on the host)
  !$ACC UPDATE HOST(x)
end subroutine ssh_cg_iterate_acc
#endif

!=========================================================================
! Chebyshev preconditioner setup: D^-1 and lmax(D^-1 A) by power iteration.
! The interval is [b/cheb_ratio, b] with b = cheb_margin*lmax.
subroutine cheb_setup(solverinfo, partit, mesh)
  IMPLICIT NONE
  type(t_solverinfo), intent(inout), target :: solverinfo
  type(t_partit),     intent(inout), target :: partit
  type(t_mesh),       intent(inout), target :: mesh
  integer                    :: n, nt, row, it
  integer, pointer           :: rptr(:), cind(:)
  real(kind=WP_full)         :: sbuf(2), lam
  real(kind=WP), allocatable :: v(:), w(:)

  cheb_deg = max(1, solverinfo%precond_cheb_degree)
  rptr => mesh%ssh_stiff%rowptr_loc
  cind => mesh%ssh_stiff%colind_loc
  n  = partit%myDim_nod2D
  nt = partit%myDim_nod2D + partit%eDim_nod2D
  if (allocated(cheb_dinv)) then
     !$ACC EXIT DATA DELETE(cheb_dinv, cheb_d, cheb_res)
     deallocate(cheb_dinv, cheb_d, cheb_res)
  end if
  allocate(cheb_dinv(n), cheb_d(nt), cheb_res(n), v(nt), w(n))
  do row = 1, n
     cheb_dinv(row) = 1.0_WP / mesh%ssh_stiff%values(rptr(row))
  end do
  ! deterministic start vector that does not depend on the partitioning
  do row = 1, nt
     v(row) = 1.0_WP + 0.1_WP*real(mod(partit%myList_nod2D(row), 13), WP)
  end do
  lam = 0.0_WP_full
  do it = 1, cheb_powit
     call exchange_nod(v, partit)
     do row = 1, n
        w(row) = real(sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                          real(v(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
     end do
     ! Rayleigh quotient v.Av / v.Dv
     sbuf(1) = sum(real(v(1:n),WP_full)*real(w(1:n),WP_full))
     sbuf(2) = sum(real(v(1:n),WP_full)**2 / real(cheb_dinv(1:n),WP_full))
     call MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
     lam = sbuf(1)/sbuf(2)
     v(1:n) = real(real(cheb_dinv(1:n),WP_full)*real(w(1:n),WP_full)/lam, WP)
  end do
  deallocate(v, w)
  cheb_b = cheb_margin*lam
  cheb_a = cheb_b/cheb_ratio
  cheb_d   = 0.0_WP
  cheb_res = 0.0_WP
  !$ACC ENTER DATA COPYIN(cheb_dinv, cheb_d, cheb_res)
  if (partit%mype == 0) write(*,'(a,i3,a,es12.5,a,es12.5,a)') '     ssh CG Chebyshev: degree=', cheb_deg, &
       '  interval=[', cheb_a, ',', cheb_b, ']'
end subroutine cheb_setup

!=========================================================================
! z = p(D^-1 A) D^-1 r by cheb_deg steps of the Chebyshev iteration on D^-1 A z = D^-1 r from z=0.
! Uses r on the owned rows only, so r needs no halo exchange; z is valid on the owned rows.
subroutine cheb_apply(r, z, partit, mesh)
  IMPLICIT NONE
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(inout), target :: mesh
  real(kind=WP),  intent(in)            :: r(:)
  real(kind=WP),  intent(inout)         :: z(:)
  integer            :: n, row, i
  integer, pointer   :: rptr(:), cind(:)
  real(kind=WP_full) :: th, de, s1, rho, rhon, c1, c2

  rptr => mesh%ssh_stiff%rowptr_loc
  cind => mesh%ssh_stiff%colind_loc
  n  = partit%myDim_nod2D
  th = 0.5_WP_full*(cheb_b + cheb_a)
  de = 0.5_WP_full*(cheb_b - cheb_a)
  s1 = th/de
  rho = 1.0_WP_full/s1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  do row = 1, n
     z(row) = real(real(cheb_dinv(row),WP_full)*real(r(row),WP_full)/th, WP)
     cheb_d(row) = z(row)
  end do
!$OMP END PARALLEL DO
  call exchange_nod(z, partit)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  do row = 1, n
     cheb_res(row) = real(real(cheb_dinv(row),WP_full)*(real(r(row),WP_full) - &
                     sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                         real(z(cind(rptr(row):rptr(row+1)-1)),WP_full))), WP)
  end do
!$OMP END PARALLEL DO
  do i = 1, cheb_deg
     rhon = 1.0_WP_full/(2.0_WP_full*s1 - rho)
     c1 = rhon*rho
     c2 = 2.0_WP_full*rhon/de
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
     do row = 1, n
        cheb_d(row) = real(c1*real(cheb_d(row),WP_full) + c2*real(cheb_res(row),WP_full), WP)
        z(row) = z(row) + cheb_d(row)
     end do
!$OMP END PARALLEL DO
     rho = rhon
     if (i < cheb_deg) then
        call exchange_nod(cheb_d, partit)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
        do row = 1, n
           cheb_res(row) = cheb_res(row) - real(real(cheb_dinv(row),WP_full)* &
                sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                    real(cheb_d(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
        end do
!$OMP END PARALLEL DO
     end if
  end do
end subroutine cheb_apply

#ifdef ENABLE_OPENACC
!=========================================================================
! cheb_apply on the device; all arrays are present on the device.
subroutine cheb_apply_acc(r, z, values, rptr, cind, nrow, partit)
  IMPLICIT NONE
  real(kind=WP),  intent(inout)         :: r(:), z(:)
  real(kind=WP),  intent(in)            :: values(:)
  integer,        intent(in)            :: rptr(:), cind(:)
  integer,        intent(in)            :: nrow
  type(t_partit), intent(inout), target :: partit
  integer            :: row, j, i
  real(kind=WP_full) :: acc, th, de, s1, rho, rhon, c1, c2

  th = 0.5_WP_full*(cheb_b + cheb_a)
  de = 0.5_WP_full*(cheb_b - cheb_a)
  s1 = th/de
  rho = 1.0_WP_full/s1
  !$ACC PARALLEL LOOP DEFAULT(PRESENT)
  do row = 1, nrow
     z(row) = real(real(cheb_dinv(row),WP_full)*real(r(row),WP_full)/th, WP)
     cheb_d(row) = z(row)
  end do
  call exchange_nod(z, partit, luse_g2g=.true.)
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j)
  do row = 1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j = rptr(row), rptr(row+1)-1
        acc = acc + real(values(j),WP_full)*real(z(cind(j)),WP_full)
     end do
     cheb_res(row) = real(real(cheb_dinv(row),WP_full)*(real(r(row),WP_full) - acc), WP)
  end do
  do i = 1, cheb_deg
     rhon = 1.0_WP_full/(2.0_WP_full*s1 - rho)
     c1 = rhon*rho
     c2 = 2.0_WP_full*rhon/de
     !$ACC PARALLEL LOOP DEFAULT(PRESENT)
     do row = 1, nrow
        cheb_d(row) = real(c1*real(cheb_d(row),WP_full) + c2*real(cheb_res(row),WP_full), WP)
        z(row) = z(row) + cheb_d(row)
     end do
     rho = rhon
     if (i < cheb_deg) then
        call exchange_nod(cheb_d, partit, luse_g2g=.true.)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j)
        do row = 1, nrow
           acc = 0.0_WP_full
           !$ACC LOOP SEQ
           do j = rptr(row), rptr(row+1)-1
              acc = acc + real(values(j),WP_full)*real(cheb_d(cind(j)),WP_full)
           end do
           cheb_res(row) = cheb_res(row) - real(real(cheb_dinv(row),WP_full)*acc, WP)
        end do
     end if
  end do
end subroutine cheb_apply_acc
#endif

!=========================================================================
! z = M^-1 r for the selected preconditioner; r on the owned rows, z valid on the owned rows.
subroutine ssh_precond_apply(r, z, solverinfo, partit, mesh)
  IMPLICIT NONE
  type(t_solverinfo), intent(inout), target :: solverinfo
  type(t_partit),     intent(inout), target :: partit
  type(t_mesh),       intent(inout), target :: mesh
  real(kind=WP),      intent(inout)         :: r(:), z(:)
  integer                :: row
  real(kind=WP), pointer :: pr_values(:)
  integer,       pointer :: rptr(:), cind(:)

  if (solverinfo%precond_variant == 2) then
     call cheb_apply(r, z, partit, mesh)
     return
  end if
  pr_values=>mesh%ssh_stiff%pr_values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc
  call exchange_nod(r, partit)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  do row=1, partit%myDim_nod2D
     z(row)= real(sum(real(pr_values(rptr(row):rptr(row+1)-1),WP_full)*       &
                      real(r(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  end do
!$OMP END PARALLEL DO
end subroutine ssh_precond_apply

!=========================================================================
! Preconditioned CG in the Chronopoulos-Gear form: the three scalar products of an iteration
! (r.z, Az.z, r.r) go into one Allreduce; s = A p is carried by a recurrence.
subroutine ssh_cg1_iterate(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  IMPLICIT NONE
  type(t_solverinfo),  intent(inout), target :: solverinfo
  type(t_partit),      intent(inout), target :: partit
  type(t_mesh),        intent(inout), target :: mesh
  real(kind=WP),       intent(inout)         :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),       intent(in)            :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  integer,             intent(out)           :: iter
  logical,             intent(out)           :: converged
  real(kind=WP_full),  intent(out)           :: sprod(2), rtol
  real(kind=WP_full)           :: s_old, dlt, gam, al, be, den, s3(3)
  integer                      :: row, n
  real(kind=WP), pointer       :: rr(:), zz(:), pp(:), App(:)
  integer,       pointer       :: rptr(:), cind(:)

  cind=>mesh%ssh_stiff%colind_loc
  rptr=>mesh%ssh_stiff%rowptr_loc
  rr  =>solverinfo%rr
  zz  =>solverinfo%zz
  pp  =>solverinfo%pp
  App =>solverinfo%App
  n=partit%myDim_nod2D
  if (.not. allocated(cg1_s)) allocate(cg1_s(n))

  s_old = sum(real(rhs(1:n), WP_full) * real(rhs(1:n), WP_full))
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  rtol=solverinfo%soltol*sqrt(s_old/real(mesh%nod2D,WP_full))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  do row=1, n
     rr(row)=real(real(rhs(row),WP_full)                                        &
              - sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)*     &
                    real(x(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
  end do
!$OMP END PARALLEL DO
  call ssh_precond_apply(rr, zz, solverinfo, partit, mesh)
  call exchange_nod(zz, partit)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
  do row=1, n
     pp(row) =zz(row)
     App(row)=real(sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                        real(zz(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
     cg1_s(row)=App(row)
  end do
!$OMP END PARALLEL DO
  s3(1) = sum(real(rr(1:n), WP_full) * real(zz(1:n), WP_full))
  s3(2) = sum(real(zz(1:n), WP_full) * real(App(1:n), WP_full))
  call MPI_Allreduce(MPI_IN_PLACE, s3, 2, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  s_old = s3(1)
  dlt   = s3(2)

  converged = .false.
  sprod     = 0.0_WP_full
  iter      = 0
  if (dlt <= 0.0_WP_full .or. dlt /= dlt) then
     solverinfo%nbreakdown = solverinfo%nbreakdown + 1
     if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
        write(*,*) 'WARNING: ssh CG breakdown at start: z.Az= ', dlt, '; keeping the initial guess'
     endif
     return
  end if
  al = s_old/dlt
  do iter=1, solverinfo%maxiter
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
     do row=1, n
        x(row) =x(row) +al*pp(row)
        rr(row)=rr(row)-al*cg1_s(row)
     end do
!$OMP END PARALLEL DO
     call ssh_precond_apply(rr, zz, solverinfo, partit, mesh)
     call exchange_nod(zz, partit)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
     do row=1, n
        App(row)=real(sum(real(mesh%ssh_stiff%values(rptr(row):rptr(row+1)-1),WP_full)* &
                           real(zz(cind(rptr(row):rptr(row+1)-1)),WP_full)), WP)
     end do
!$OMP END PARALLEL DO
     s3(1) = sum(real(rr(1:n), WP_full) * real(zz(1:n), WP_full))
     s3(2) = sum(real(App(1:n), WP_full) * real(zz(1:n), WP_full))
     s3(3) = sum(real(rr(1:n), WP_full) * real(rr(1:n), WP_full))
     call MPI_Allreduce(MPI_IN_PLACE, s3, 3, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
     gam      = s3(1)
     dlt      = s3(2)
     sprod(1) = s3(1)
     sprod(2) = s3(3)
     if (sqrt(sprod(2)/mesh%nod2D) < rtol) then
        converged = .true.
        exit
     endif
     if (gam < 0.0_WP_full) then
        solverinfo%nnegrz = solverinfo%nnegrz + 1
        if (partit%mype==0 .and. solverinfo%nnegrz == 1) then
           write(*,*) 'WARNING: ssh CG r.z < 0 at iter ', iter, ': r.z= ', gam, &
                      ' -- preconditioner is not symmetric positive definite'
        endif
     endif
     be  = gam/s_old
     den = dlt - be*gam/al
     if (den <= 0.0_WP_full .or. den /= den) then
        solverinfo%nbreakdown = solverinfo%nbreakdown + 1
        if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
           write(*,*) 'WARNING: ssh CG breakdown at iter ', iter, ': p.Ap= ', den, &
                      ' is not positive; keeping current solution and leaving the loop'
        endif
        exit
     endif
     al    = gam/den
     s_old = gam
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
     do row=1, n
        pp(row)   =zz(row) +be*pp(row)
        cg1_s(row)=App(row)+be*cg1_s(row)
     end do
!$OMP END PARALLEL DO
  end do
end subroutine ssh_cg1_iterate

#ifdef ENABLE_OPENACC
!=========================================================================
! ssh_precond_apply on the device
subroutine ssh_precond_apply_acc(r, z, values, pr_values, rptr, cind, nrow, pv, partit)
  IMPLICIT NONE
  real(kind=WP),  intent(inout)         :: r(:), z(:)
  real(kind=WP),  intent(in)            :: values(:), pr_values(:)
  integer,        intent(in)            :: rptr(:), cind(:)
  integer,        intent(in)            :: nrow, pv
  type(t_partit), intent(inout), target :: partit
  integer            :: row, j
  real(kind=WP_full) :: acc

  if (pv == 2) then
     call cheb_apply_acc(r, z, values, rptr, cind, nrow, partit)
     return
  end if
  call exchange_nod(r, partit, luse_g2g=.true.)
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j)
  do row=1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j=rptr(row), rptr(row+1)-1
        acc = acc + real(pr_values(j),WP_full)*real(r(cind(j)),WP_full)
     end do
     z(row)=real(acc, WP)
  end do
end subroutine ssh_precond_apply_acc

!=========================================================================
! ssh_cg1_iterate on the device; x leaves on the host, its halo is exchanged by the caller.
subroutine ssh_cg1_iterate_acc(x, rhs, solverinfo, partit, mesh, iter, converged, sprod, rtol)
  IMPLICIT NONE
  type(t_solverinfo),  intent(inout), target :: solverinfo
  type(t_partit),      intent(inout), target :: partit
  type(t_mesh),        intent(inout), target :: mesh
  real(kind=WP),       intent(inout)         :: x(partit%myDim_nod2D+partit%eDim_nod2D)
  real(kind=WP),       intent(in)            :: rhs(partit%myDim_nod2D+partit%eDim_nod2D)
  integer,             intent(out)           :: iter
  logical,             intent(out)           :: converged
  real(kind=WP_full),  intent(out)           :: sprod(2), rtol
  real(kind=WP_full)           :: s_old, dlt, gam, al, be, den, acc, sa, sb, sc, s3(3)
  integer                      :: row, j, nrow, pv
  logical, save                :: on_device = .false.
  real(kind=WP), pointer       :: pr_values(:), values(:), rr(:), zz(:), pp(:), App(:)
  integer,       pointer       :: rptr(:), cind(:)

  pr_values=>mesh%ssh_stiff%pr_values
  values   =>mesh%ssh_stiff%values
  cind     =>mesh%ssh_stiff%colind_loc
  rptr     =>mesh%ssh_stiff%rowptr_loc
  rr =>solverinfo%rr
  zz =>solverinfo%zz
  pp =>solverinfo%pp
  App=>solverinfo%App
  nrow=partit%myDim_nod2D
  pv  =solverinfo%precond_variant
  if (.not. on_device) then
     allocate(cg1_s(nrow))
     cg1_s = 0.0_WP
     !$ACC ENTER DATA COPYIN(cind, rptr, pr_values, cg1_s) CREATE(values, rr, zz, pp, App, x, rhs)
     on_device = .true.
  end if
  !$ACC UPDATE DEVICE(values, x, rhs)

  s_old = 0.0_WP_full
  !$ACC PARALLEL LOOP DEFAULT(PRESENT) REDUCTION(+:s_old)
  do row=1, nrow
     s_old = s_old + real(rhs(row),WP_full)*real(rhs(row),WP_full)
  end do
  call MPI_Allreduce(MPI_IN_PLACE, s_old, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  rtol=solverinfo%soltol*sqrt(s_old/real(mesh%nod2D,WP_full))

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j)
  do row=1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j=rptr(row), rptr(row+1)-1
        acc = acc + real(values(j),WP_full)*real(x(cind(j)),WP_full)
     end do
     rr(row)=real(real(rhs(row),WP_full)-acc, WP)
  end do
  call ssh_precond_apply_acc(rr, zz, values, pr_values, rptr, cind, nrow, pv, partit)
  call exchange_nod(zz, partit, luse_g2g=.true.)
  sa = 0.0_WP_full
  sb = 0.0_WP_full
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j) REDUCTION(+:sa, sb)
  do row=1, nrow
     acc = 0.0_WP_full
     !$ACC LOOP SEQ
     do j=rptr(row), rptr(row+1)-1
        acc = acc + real(values(j),WP_full)*real(zz(cind(j)),WP_full)
     end do
     pp(row)   = zz(row)
     App(row)  = real(acc, WP)
     cg1_s(row)= App(row)
     sa = sa + real(rr(row),WP_full)*real(zz(row),WP_full)
     sb = sb + real(zz(row),WP_full)*acc
  end do
  sprod(1)=sa
  sprod(2)=sb
  call MPI_Allreduce(MPI_IN_PLACE, sprod, 2, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
  s_old=sprod(1)
  dlt  =sprod(2)

  converged = .false.
  sprod     = 0.0_WP_full
  iter      = 0
  if (dlt <= 0.0_WP_full .or. dlt /= dlt) then
     solverinfo%nbreakdown = solverinfo%nbreakdown + 1
     if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
        write(*,*) 'WARNING: ssh CG breakdown at start: z.Az= ', dlt, '; keeping the initial guess'
     endif
     !$ACC UPDATE HOST(x)
     return
  end if
  al = s_old/dlt
  do iter=1, solverinfo%maxiter
     !$ACC PARALLEL LOOP DEFAULT(PRESENT)
     do row=1, nrow
        x(row) = x(row) + al*pp(row)
        rr(row)= rr(row)- al*cg1_s(row)
     end do
     call ssh_precond_apply_acc(rr, zz, values, pr_values, rptr, cind, nrow, pv, partit)
     call exchange_nod(zz, partit, luse_g2g=.true.)
     sa = 0.0_WP_full
     sb = 0.0_WP_full
     sc = 0.0_WP_full
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(acc, j) REDUCTION(+:sa, sb, sc)
     do row=1, nrow
        acc = 0.0_WP_full
        !$ACC LOOP SEQ
        do j=rptr(row), rptr(row+1)-1
           acc = acc + real(values(j),WP_full)*real(zz(cind(j)),WP_full)
        end do
        App(row) = real(acc, WP)
        sa = sa + real(rr(row),WP_full)*real(zz(row),WP_full)
        sb = sb + real(zz(row),WP_full)*acc
        sc = sc + real(rr(row),WP_full)*real(rr(row),WP_full)
     end do
     s3 = (/ sa, sb, sc /)
     call MPI_Allreduce(MPI_IN_PLACE, s3, 3, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
     gam      = s3(1)
     dlt      = s3(2)
     sprod(1) = s3(1)
     sprod(2) = s3(3)
     if (sqrt(sprod(2)/mesh%nod2D) < rtol) then
        converged = .true.
        exit
     endif
     if (gam < 0.0_WP_full) then
        solverinfo%nnegrz = solverinfo%nnegrz + 1
        if (partit%mype==0 .and. solverinfo%nnegrz == 1) then
           write(*,*) 'WARNING: ssh CG r.z < 0 at iter ', iter, ': r.z= ', gam, &
                      ' -- preconditioner is not symmetric positive definite'
        endif
     endif
     be  = gam/s_old
     den = dlt - be*gam/al
     if (den <= 0.0_WP_full .or. den /= den) then
        solverinfo%nbreakdown = solverinfo%nbreakdown + 1
        if (partit%mype==0 .and. solverinfo%nbreakdown <= 5) then
           write(*,*) 'WARNING: ssh CG breakdown at iter ', iter, ': p.Ap= ', den, &
                      ' is not positive; keeping current solution and leaving the loop'
        endif
        exit
     endif
     al    = gam/den
     s_old = gam
     !$ACC PARALLEL LOOP DEFAULT(PRESENT)
     do row=1, nrow
        pp(row)   = zz(row) + be*pp(row)
        cg1_s(row)= App(row)+ be*cg1_s(row)
     end do
  end do
  !$ACC UPDATE HOST(x)
end subroutine ssh_cg1_iterate_acc
#endif

! ===================================================================

end module solver_module
