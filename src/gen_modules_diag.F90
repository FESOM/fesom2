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
  public :: compute_diagnostics, ldiag_solver, rhs_diag, lcurt_stress_surf, curl_stress_surf, ldiag_curl_vel3, curl_vel3, ldiag_energy, wzmid, wzmidrho, rho

  ! Arrays used for diagnostics, some shall be accessible to the I/O
  ! 1. solver diagnostics: A*x=rhs? 
  ! A=ssh_stiff, x=d_eta, rhs=ssh_rhs; rhs_diag=A*x;
  real(kind=8),  save, allocatable, target      :: rhs_diag(:)
  real(kind=8),  save, allocatable, target      :: curl_stress_surf(:)
  real(kind=8),  save, allocatable, target      :: curl_vel3(:,:)
  real(kind=8),  save, allocatable, target      :: wzmid(:,:), wzmidrho(:,:), rho(:,:)

  logical                                       :: ldiag_solver     =.false.
  logical                                       :: lcurt_stress_surf=.false.
  logical                                       :: ldiag_curl_vel3  =.false.
  logical                                       :: ldiag_energy     =.false.
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
!curt(stress_surf)
subroutine diag_curl_stress_surf(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer        :: enodes(2), el(2), ed, n
  real(kind=8)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

!=====================

  if (firstcall) then  !allocate the stuff at the first call
     allocate(curl_stress_surf(myDim_nod2D+eDim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  curl_stress_surf=0.

  DO ed=1, myDim_edge2D
     enodes=edges(:,ed)
     el=edge_tri(:,ed)
     deltaX1=edge_cross_dxdy(1,ed)
     deltaY1=edge_cross_dxdy(2,ed)
     if (el(2)>0) then
        deltaX2=edge_cross_dxdy(3,ed)
        deltaY2=edge_cross_dxdy(4,ed)
     end if
     c1=deltaX1*stress_surf(1,el(1))+deltaY1*stress_surf(2,el(1))
     curl_stress_surf(enodes(1))=curl_stress_surf(enodes(1))+c1
     curl_stress_surf(enodes(2))=curl_stress_surf(enodes(2))-c1
     if (el(2)>0) then
        c1= -deltaX2*stress_surf(1,el(2))-deltaY2*stress_surf(2,el(2))
        curl_stress_surf(enodes(1))=curl_stress_surf(enodes(1))+c1
        curl_stress_surf(enodes(2))=curl_stress_surf(enodes(2))-c1
     end if
  END DO
  DO n=1, myDim_nod2D+eDim_nod2D
     curl_stress_surf(n)=curl_stress_surf(n)/area(1,n)
  END DO
end subroutine diag_curl_stress_surf
! ==============================================================
!3D curl(velocity)
subroutine diag_curl_vel3(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer        :: enodes(2), el(2), ed, n, nz, nl1, nl2
  real(kind=8)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(curl_vel3(nl-1, myDim_nod2D+eDim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  curl_vel3=0.

  DO ed=1,myDim_edge2D
     enodes=edges(:,ed)
     el=edge_tri(:,ed)
     nl1=nlevels(el(1))-1
     deltaX1=edge_cross_dxdy(1,ed)
     deltaY1=edge_cross_dxdy(2,ed)
     nl2=0
     if (el(2)>0) then
        deltaX2=edge_cross_dxdy(3,ed)
        deltaY2=edge_cross_dxdy(4,ed)
        nl2=nlevels(el(2))-1
      end if     
      DO nz=1,min(nl1,nl2)
         c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
         deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
      DO nz=min(nl1,nl2)+1,nl1
         c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
      DO nz=min(nl1,nl2)+1,nl2
         c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
         curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
         curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
      END DO
   END DO

   DO n=1, myDim_nod2D
      DO nz=1, nlevels_nod2D(n)-1
         curl_vel3(nz,n)=curl_vel3(nz,n)/area(nz,n)
      END DO
   END DO    
end subroutine diag_curl_vel3
! ==============================================================
!energy budget
subroutine diag_energy(mode)
  implicit none
  integer, intent(in)           :: mode
  logical, save                 :: firstcall=.true.
  integer        :: enodes(2), el(2), ed, n, nz, nl1, nl2
  real(kind=8)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(wzmid(nl-1, myDim_nod2D), wzmidrho(nl-1, myDim_nod2D), rho(nl-1, myDim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if
  wzmid=0.
  rho  =0.
  wzmidrho=0.

  wzmid=0.5*(Wvel(1:nl-1, 1:myDim_nod2D)+Wvel(2:nl, 1:myDim_nod2D))
  rho  =density_m_rho0(1:nl-1, 1:myDim_nod2D)
  where(abs(rho)>0.)
     rho=rho+density_0
  end where
  wzmidrho=rho(1:nl-1, 1:myDim_nod2D)*wzmid(1:nl-1, 1:myDim_nod2D)
end subroutine diag_energy
! ==============================================================
subroutine compute_diagnostics(mode)
  implicit none
  integer, intent(in)           :: mode !constructor mode (0=only allocation; any other=do diagnostic)

  !1. solver diagnostic
  if (ldiag_solver)      call diag_solver(mode)
  !2. compute curl(stress_surf)
  if (lcurt_stress_surf) call diag_curl_stress_surf(mode)
  !3. compute curl(velocity)
  if (ldiag_curl_vel3)   call diag_curl_vel3(mode)
  !4. compute energy budget
  if (ldiag_energy)      call diag_energy(mode)
end subroutine compute_diagnostics
end module diagnostics
