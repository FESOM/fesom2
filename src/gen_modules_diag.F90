module diagnostics

  use g_config
  use o_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use i_ARRAYS
  use o_mixing_KPP_mod
  implicit none

  private
  public :: compute_diagnostics, rhs_diag, ldiag_solver

  ! Arrays used for diagnostics, some shall be accessible to the I/O
  ! 1. solver diagnostics: A*x=rhs? 
  ! A=ssh_stiff, x=d_eta, rhs=ssh_rhs; rhs_diag=A*x;
  real(kind=8),  save, allocatable, target      :: rhs_diag(:)
  logical                                       :: ldiag_solver=.false.
  contains

! ==============================================================
!rhs_diag=ssh_rhs?
subroutine diag_solver(mode)
  implicit none
  integer, intent(in)           :: mode
  integer                       :: n, is, ie
  logical, save                 :: firstcall=.true.

!=====================

  if (firstcall) then !allocate the stuff at the first call
     allocate(rhs_diag(mydim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  do n=1, myDim_nod2D
     is=ssh_stiff%rowptr_loc(n)
     ie=ssh_stiff%rowptr_loc(n+1)-1
     rhs_diag(n)=sum(ssh_stiff%values(is:ie)*d_eta(ssh_stiff%colind_loc(is:ie)))
  end do
end subroutine diag_solver
! ==============================================================
subroutine compute_diagnostics(mode)
  implicit none
  integer, intent(in)           :: mode !constructor mode (0=only allocation; any other=do diagnostic)

  !1. solver diagnostic
  if (ldiag_solver) call diag_solver(mode)

end subroutine compute_diagnostics
end module diagnostics
