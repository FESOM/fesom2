module diagnostics

  use g_config
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use MOD_TRACER
  use MOD_DYN
  use MOD_ICE
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use o_mixing_KPP_mod
  use g_rotate_grid
  use g_support
  implicit none

  private
  public :: ldiag_solver, lcurt_stress_surf, ldiag_Ri, ldiag_TurbFlux, ldiag_dMOC, ldiag_DVD,        &
            ldiag_forc, ldiag_salt3D, ldiag_curl_vel3, diag_list, ldiag_vorticity, ldiag_extflds, ldiag_ice,   &
            compute_diagnostics, rhs_diag, curl_stress_surf, curl_vel3, shear, Ri, KvdTdZ, KvdSdZ,   & 
            std_dens_min, std_dens_max, std_dens_N, std_dens, ldiag_trflx,                           &
            std_dens_UVDZ, std_dens_DIV, std_dens_DIV_fer, std_dens_Z, std_dens_H, std_dens_dVdT, std_dens_flux,       &
            dens_flux_e, vorticity, zisotherm, tempzavg, saltzavg, compute_diag_dvd_2ndmoment_klingbeil_etal_2014,       &
            compute_diag_dvd_2ndmoment_burchard_etal_2008, compute_diag_dvd, vol_ice, vol_snow, compute_ice_diag, thetao, tuv, suv
  ! Arrays used for diagnostics, some shall be accessible to the I/O
  ! 1. solver diagnostics: A*x=rhs? 
  ! A=ssh_stiff, x=d_eta, rhs=ssh_rhs; rhs_diag=A*x;
  real(kind=WP),  save, allocatable, target      :: rhs_diag(:)
  real(kind=WP),  save, allocatable, target      :: curl_stress_surf(:)
  real(kind=WP),  save, allocatable, target      :: curl_vel3(:,:)

  real(kind=WP),  save, allocatable, target      :: shear(:,:), Ri(:,:), KvdTdZ(:,:), KvdSdZ(:,:)
  real(kind=WP),  save, allocatable, target      :: stress_bott(:,:), u_bott(:), v_bott(:), u_surf(:), v_surf(:)
  real(kind=WP),  save, allocatable, target      :: vorticity(:,:)
  real(kind=WP),  save, allocatable, target      :: zisotherm(:)              !target temperature is specified as whichtemp in compute_extflds
  real(kind=WP),  save, allocatable, target      :: tempzavg(:), saltzavg(:)  !target depth for averaging is specified as whichdepth in compute_extflds
  real(kind=WP),  save, allocatable, target      :: vol_ice(:),  vol_snow(:)
  ! defining a set of standard density bins which will be used for computing densMOC
! integer,        parameter                      :: std_dens_N  = 100
! real(kind=WP),  save, target                   :: std_dens(std_dens_N)
  integer,        parameter                      :: std_dens_N  =89
  real(kind=WP),  save, target                   :: std_dens(std_dens_N)=(/ &
                                                            0.0000,   30.00000, 30.55556, 31.11111, 31.36000, 31.66667, 31.91000, 32.22222, 32.46000, &
                                                            32.77778, 33.01000, 33.33333, 33.56000, 33.88889, 34.11000, 34.44444, 34.62000, 35.00000, &
                                                            35.05000, 35.10622, 35.20319, 35.29239, 35.37498, 35.41300, 35.45187, 35.52380, 35.59136, &
                                                            35.65506, 35.71531, 35.77247, 35.82685, 35.87869, 35.92823, 35.97566, 35.98000, 36.02115, &
                                                            36.06487, 36.10692, 36.14746, 36.18656, 36.22434, 36.26089, 36.29626, 36.33056, 36.36383, &
                                                            36.39613, 36.42753, 36.45806, 36.48778, 36.51674, 36.54495, 36.57246, 36.59500, 36.59932, &
                                                            36.62555, 36.65117, 36.67621, 36.68000, 36.70071, 36.72467, 36.74813, 36.75200, 36.77111, &
                                                            36.79363, 36.81570, 36.83733, 36.85857, 36.87500, 36.87940, 36.89985, 36.91993, 36.93965, &
                                                            36.95904, 36.97808, 36.99682, 37.01524, 37.03336, 37.05119, 37.06874, 37.08602, 37.10303, &
                                                            37.11979, 37.13630, 37.15257, 37.16861, 37.18441, 37.50000, 37.75000, 40.00000/)
  real(kind=WP),  save, target                   :: std_dd(std_dens_N-1)
  real(kind=WP),  save, target                   :: std_dens_min=1030., std_dens_max=1040.
  real(kind=WP),  save, allocatable, target      :: std_dens_UVDZ(:,:,:), std_dens_flux(:,:,:), std_dens_dVdT(:,:), std_dens_DIV(:,:), std_dens_DIV_fer(:,:), std_dens_Z(:,:), std_dens_H(:,:)
  real(kind=WP),  save, allocatable, target      :: dens_flux_e(:)
  real(kind=WP),  save, allocatable, target      :: thetao(:) ! sst in K
  real(kind=WP),  save, allocatable, target      :: tuv(:,:,:), suv(:,:,:)

  logical                                       :: ldiag_solver     =.false.
  logical                                       :: lcurt_stress_surf=.false.
  logical                                       :: ldiag_curl_vel3  =.false.
  logical                                       :: ldiag_Ri         =.false.
  logical                                       :: ldiag_TurbFlux   =.false.
  logical                                       :: ldiag_KE         =.false.
  logical                                       :: ldiag_salt3D     =.false.
  ! this option activates writing the horizintal velocity transports within the density bins (U_rho_x_DZ and V_rho_x_DZ)
  ! an additional field (RHO_Z) will be computed which allows for diagnosing the numerical diapycnal mixing after A. Megann 2018
  logical                                       :: ldiag_dMOC       =.false.
  
  ! flag for calculating the Discrete Variance Decay --> estimator for numerical/
  ! spurious mixing in the advection schemes
  logical                                       :: ldiag_DVD        =.false.
  
  logical                                       :: ldiag_forc       =.false.
  
  logical                                       :: ldiag_vorticity  =.false.
  logical                                       :: ldiag_extflds    =.false.
  logical                                       :: ldiag_ice        =.false.
  logical                                       :: ldiag_trflx      =.false.

  namelist /diag_list/ ldiag_solver, lcurt_stress_surf, ldiag_curl_vel3, ldiag_Ri, ldiag_TurbFlux, ldiag_dMOC, ldiag_DVD, &
                       ldiag_salt3D, ldiag_forc, ldiag_vorticity, ldiag_extflds, ldiag_trflx
  
  contains

! ==============================================================
!rhs_diag=ssh_rhs?
subroutine diag_solver(mode, dynamics, partit, mesh)
  implicit none
  type(t_mesh)  , intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_dyn)   , intent(inout), target :: dynamics
  integer,        intent(in)            :: mode
  integer                               :: n, is, ie
  logical, save                         :: firstcall=.true.
  real(kind=WP), dimension(:)    , pointer :: d_eta
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  d_eta    =>dynamics%d_eta(:)
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
subroutine diag_curl_stress_surf(mode, partit, mesh)
  implicit none
  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer,        intent(in)            :: mode
  logical,        save                  :: firstcall=.true.
  integer                               :: enodes(2), el(2), ed, n
  real(kind=WP)                         :: deltaX1, deltaY1, deltaX2, deltaY2, c1
!=====================
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

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
     !!PS curl_stress_surf(n)=curl_stress_surf(n)/area(1,n)
     curl_stress_surf(n)=curl_stress_surf(n)/areasvol(ulevels_nod2D(n),n)
  END DO
end subroutine diag_curl_stress_surf
! ==============================================================
!3D curl(velocity)
subroutine diag_curl_vel3(mode, dynamics, partit, mesh)
    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    integer,        intent(in)            :: mode
    logical,        save                  :: firstcall=.true.
    integer                               :: enodes(2), el(2), ed, n, nz, nl1, nl2, nl12, nu1, nu2, nu12
    real(kind=WP)                         :: deltaX1, deltaY1, deltaX2, deltaY2, c1
    real(kind=WP), dimension(:,:,:), pointer :: UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV => dynamics%uv(:,:,:)

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
        nu1=ulevels(el(1))
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        nl2=0
        nu2=0
        if (el(2)>0) then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
        end if   
        
        nl12 = min(nl1,nl2)
        nu12 = min(nu1,nu2)
        DO nz=nu1,nu12-1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
            curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
        END DO
        if (nu2>0) then
            DO nz=nu2,nu12-1
                c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            END DO
        end if
        DO nz=nu12,nl12
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
            deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
            curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
        END DO
        DO nz=nl12+1,nl1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
            curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
        END DO
        DO nz=nl12+1,nl2
            c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
            curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
        END DO
    END DO

    DO n=1, myDim_nod2D
        !!PS DO nz=1, nlevels_nod2D(n)-1
        DO nz=ulevels_nod2D(n), nlevels_nod2D(n)-1
            curl_vel3(nz,n)=curl_vel3(nz,n)/areasvol(nz,n)
        END DO
    END DO
end subroutine diag_curl_vel3
! ==============================================================
! 
subroutine diag_turbflux(mode, dynamics, tracers, partit, mesh)
  implicit none
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_tracer), intent(in)   , target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh)  , intent(in)   , target :: mesh
  integer,        intent(in)            :: mode
  logical,        save                     :: firstcall=.true.
  integer                                  :: n, nz, nzmax, nzmin
  real(kind=WP), dimension(:,:,:), pointer :: UVnode
  real(kind=WP), dimension(:,:),   pointer :: temp, salt
  real(kind=WP)                            :: dz_inv

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
UVnode=>dynamics%uvnode(:,:,:)
temp   => tracers%data(1)%values(:,:)
salt   => tracers%data(2)%values(:,:)
!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(KvdTdZ(nl, myDim_nod2D+eDim_nod2D), KvdSdZ(nl, myDim_nod2D+eDim_nod2D))
     KvdTdZ =0.0_WP
     KvdSdZ =0.0_WP
     firstcall=.false.
     if (mode==0) return
  end if  

  do n=1, myDim_nod2D+eDim_nod2D
     nzmin = ulevels_nod2d(n)
     nzmax = nlevels_nod2d(n)
     do nz=nzmin+1,nzmax-1
        dz_inv=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))
        KvdTdZ(nz,n) = -Kv(nz,n)*(temp(nz-1,n)-temp(nz,n))*dz_inv
        KvdSdZ(nz,n) = -Kv(nz,n)*(salt(nz-1,n)-salt(nz,n))*dz_inv
     end do
  end do
end subroutine diag_turbflux
! ==============================================================
!
subroutine diag_trflx(mode, dynamics, tracers, partit, mesh)
  implicit none
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_tracer), intent(in)   , target :: tracers
  type(t_partit), intent(inout), target :: partit
  type(t_mesh)  , intent(in)   , target :: mesh
  integer,        intent(in)            :: mode
  logical,        save                     :: firstcall=.true.
  integer                                  :: elem, nz, nzu, nzl, elnodes(3)
  real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV
  real(kind=WP), dimension(:,:),   pointer :: temp, salt

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
UV     => dynamics%uv(:,:,:)
temp   => tracers%data(1)%values(:,:)
salt   => tracers%data(2)%values(:,:)
fer_UV => dynamics%fer_uv(:,:,:)
!=====================
  if (firstcall) then !allocate the stuff at the first call
      allocate(tuv(2,nl-1,myDim_elem2D+eDim_elem2D))
      allocate(suv(2,nl-1,myDim_elem2D+eDim_elem2D))
      tuv = 0.0_WP
      suv = 0.0_WP
      firstcall=.false.
      if (mode==0) return
  end if

  !___________________________________________________________________________
  ! compute tracer fluxes
  do elem=1,myDim_elem2D
      elnodes = elem2D_nodes(:,elem)
      nzu     = ulevels(elem)
      nzl     = nlevels(elem)-1
      if (Fer_GM) then
          do nz=nzu, nzl
              tuv(1,nz,elem) = (UV(1,nz,elem) + fer_UV(1,nz, elem)) * sum(temp(nz,elnodes))/3._WP
              tuv(2,nz,elem) = (UV(2,nz,elem) + fer_UV(2,nz, elem)) * sum(temp(nz,elnodes))/3._WP
              suv(1,nz,elem) = (UV(1,nz,elem) + fer_UV(1,nz, elem)) * sum(salt(nz,elnodes))/3._WP
              suv(2,nz,elem) = (UV(2,nz,elem) + fer_UV(2,nz, elem)) * sum(salt(nz,elnodes))/3._WP
          end do
      else
          do nz=nzu, nzl
              tuv(1,nz,elem) = UV(1,nz,elem) * sum(temp(nz,elnodes))/3._WP
              tuv(2,nz,elem) = UV(2,nz,elem) * sum(temp(nz,elnodes))/3._WP
              suv(1,nz,elem) = UV(1,nz,elem) * sum(salt(nz,elnodes))/3._WP
              suv(2,nz,elem) = UV(2,nz,elem) * sum(salt(nz,elnodes))/3._WP
          end do
      end if
  end do
end subroutine diag_trflx
! ==============================================================
! 
subroutine diag_Ri(mode, dynamics, partit, mesh)
  implicit none
  type(t_dyn)   , intent(inout), target :: dynamics
  type(t_partit), intent(inout), target :: partit
  type(t_mesh)  , intent(in)   , target :: mesh
  integer,        intent(in)            :: mode
  logical,        save                     :: firstcall=.true.
  integer                                  :: n, nz, nzmax, nzmin
  real(kind=WP), dimension(:,:,:), pointer :: UVnode
  real(kind=WP)                            :: val, dz_inv

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
UVnode=>dynamics%uvnode(:,:,:)
!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(shear(nl, myDim_nod2D+eDim_nod2D), Ri(nl, myDim_nod2D+eDim_nod2D))
     shear=0.0_WP
     Ri   =0.0_WP
     firstcall=.false.
     if (mode==0) return
  end if  

  do n=1, myDim_nod2D+eDim_nod2D
     nzmin = ulevels_nod2d(n)
     nzmax = nlevels_nod2d(n)
     do nz=nzmin+1,nzmax-1
        dz_inv=1.0_WP/(Z_3d_n(nz-1,n)-Z_3d_n(nz,n))
        val =   (UVnode(1,nz-1,n)-UVnode(1,nz,n))**2 +&
                (UVnode(2,nz-1,n)-UVnode(2,nz,n))**2
        shear(nz,n) = val*dz_inv*dz_inv
        Ri(nz,n)    = bvfreq(nz,n)/max(shear(nz,n), 1.e-12)
     end do
  end do
end subroutine diag_Ri
! ==============================================================
subroutine diag_densMOC(mode, dynamics, tracers, partit, mesh)
  implicit none
  integer, intent(in)                     :: mode
  type(t_mesh)  , intent(in)   , target   :: mesh
  type(t_partit), intent(inout), target   :: partit
  type(t_tracer), intent(in)   , target   :: tracers
  type(t_dyn)   , intent(in)   , target   :: dynamics
  integer                                 :: nz, snz, elem, nzmax, nzmin, elnodes(3), is, ie, pos
  integer                                 :: e, edge, enodes(2), eelems(2)
  real(kind=WP)                           :: div, deltaX, deltaY, locz
  integer                                 :: jj
  real(kind=WP), save                     :: dd
  real(kind=WP)                           :: uvdz_el(2), rhoz_el, vol_el, dz, weight, dmin, dmax, ddiff, test, test1, test2, test3
  real(kind=WP), save, allocatable        :: dens(:), aux(:), el_depth(:)
  real(kind=WP), save, allocatable        :: std_dens_w(:,:), std_dens_VOL1(:,:), std_dens_VOL2(:,:)
  logical, save                           :: firstcall_s=.true., firstcall_e=.true.
  real(kind=WP), dimension(:,:), pointer  :: temp, salt
  real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
  UV     => dynamics%uv(:,:,:)
  temp   => tracers%data(1)%values(:,:)
  salt   => tracers%data(2)%values(:,:)
  fer_UV => dynamics%fer_uv(:,:,:)

  if (firstcall_s) then !allocate the stuff at the first call
     allocate(std_dens_UVDZ(2,std_dens_N, myDim_elem2D))
     allocate(std_dens_w   (  std_dens_N, myDim_elem2D))
     allocate(std_dens_dVdT(  std_dens_N, myDim_elem2D))
     allocate(std_dens_DIV (  std_dens_N, myDim_nod2D+eDim_nod2D))
     if (Fer_GM) allocate(std_dens_DIV_fer(  std_dens_N, myDim_nod2D+eDim_nod2D))
     allocate(std_dens_VOL1(  std_dens_N, myDim_elem2D))
     allocate(std_dens_VOL2(  std_dens_N, myDim_elem2D))
     allocate(std_dens_flux(3,std_dens_N, myDim_elem2D))
     allocate(std_dens_Z   (  std_dens_N, myDim_elem2D))
     allocate(std_dens_H   (  std_dens_N, myDim_elem2D))
     allocate(dens_flux_e(elem2D))
     allocate(aux  (nl-1))
     allocate(dens (nl))
     allocate(el_depth(nl))
!
!std_dens(1)=0.
!std_dens(2)=30.
!do nz=3, std_dens_N-1
!std_dens(nz)=std_dens(nz-1)+10.5/real(std_dens_N-2)
!end do
!std_dens(std_dens_N)=40.
!
     std_dd(:)=std_dens(2:)-std_dens(:std_dens_N-1)
     dens         =0.
     std_dens_UVDZ=0. !will be U & V transports within the density class
     std_dens_dVdT=0. !rate of change of a bin volume (for estimating the 'model drift')
     std_dens_DIV =0. !meridional divergence within a density bin (for reconstruction of the diapycnal velocity) !TOP PRIORITY
     if (Fer_GM) std_dens_DIV_fer =0. !meridional divergence of bolus velocity within a density bin (for reconstruction of the diapycnal velocity) !TOP PRIORITY
     std_dens_VOL1=0. !temporal arrays for computing std_dens_dVdT
     std_dens_VOL2=0.
     std_dens_flux=0. !bouyancy flux for computation of surface bouyancy transformations
     std_dens_Z   =0. !will be the vertical position of the density class (for convertion between dAMOC <-> zMOC)
     std_dens_H   =0. !will be the vertical layerthickness of the density class (for convertion between dAMOC <-> zMOC)
     depth        =0.
     el_depth     =0.
     firstcall_s=.false.
     if (mode==0) return
  end if

  std_dens_UVDZ=0.
  std_dens_w   =0.! temporat thing for wiighting (ageraging) mean fields within a bin
  std_dens_flux=0.
  dens_flux_e    =0.
  std_dens_VOL2=0.
  std_dens_DIV =0.
  if (Fer_GM) std_dens_DIV_fer =0. !meridional divergence of bolus velocity within a density bin (for reconstruction of the diapycnal velocity) !TOP PRIORITY
  std_dens_Z   =0.
  std_dens_H   =0.
  
  ! proceed with fields at elements...
  do elem=1, myDim_elem2D
     elnodes=elem2D_nodes(:,elem)     
     nzmax = nlevels(elem)
     nzmin = ulevels(elem)
     
     ! density flux on elements (although not related to binning it might be usefull for diagnostic and to verify the consistency)
     do jj=1,3
       dens_flux_e(elem)=dens_flux_e(elem) + (sw_alpha(ulevels_nod2D(elnodes(jj)),elnodes(jj)) * heat_flux_in(elnodes(jj))  / vcpw + &
                                            sw_beta(ulevels_nod2D(elnodes(jj)),elnodes(jj)) * (relax_salt  (elnodes(jj)) + water_flux(elnodes(jj)) * & 
                                            salt(ulevels_nod2D(elnodes(jj)),elnodes(jj))))
     end do 
     dens_flux_e(elem) =dens_flux_e(elem)/3.0_WP
     ! density_dmoc is the sigma_2 density given at nodes. it is computed in oce_ale_pressure_bv
     do nz=nzmin, nzmax-1
        aux(nz)=sum(density_dmoc(nz, elnodes))/3.-1000.
     end do

     ! dens will be the density within the column at nodes     
     el_depth(nzmax)=zbar_e_bot(elem)
     do nz=nzmax-1,nzmin+1,-1
        dens(nz)       = (aux(nz)     * helem(nz-1,elem)+&
                          aux(nz-1)   * helem(nz,  elem))/sum(helem(nz-1:nz,elem))
        el_depth(nz)   = el_depth(nz+1) + helem(nz, elem)
     end do
     dens(nzmax)=dens(nzmax-1)+(dens(nzmax-1)-dens(nzmax-2))*helem(nzmax-1,elem)/helem(nzmax-2,elem)
     dens(nzmin)    =dens(nzmin+1)      +(dens(nzmin+1)-dens(nzmin+2))            *helem(nzmin, elem)/helem(nzmin+1,elem)
     el_depth(1)=0.

     ! heat, freshwater and restoring at density classes
     is=minloc(abs(std_dens-dens(1)),1)
     std_dens_flux(1, is,elem)=std_dens_flux(1, is,elem)+elem_area(elem)*sum(sw_alpha(1,elnodes) * heat_flux_in(elnodes   ))/3./vcpw
     std_dens_flux(2, is,elem)=std_dens_flux(2, is,elem)+elem_area(elem)*sum(sw_beta (1,elnodes) * relax_salt  (elnodes   ))/3.
     
     dd = 0.0_WP
     do jj=1,3
        dd = dd + (sw_beta (1,elnodes(jj)) * water_flux(elnodes(jj)) * salt(ulevels_nod2D(elnodes(jj)),  elnodes(jj)))
     end do
     std_dens_flux(3, is,elem)=std_dens_flux(3, is,elem)+elem_area(elem)*dd/3.
     
     do nz=nzmax-1,nzmin,-1
        dmin=minval(dens(nz:nz+1))
        dmax=maxval(dens(nz:nz+1))
        ddiff=abs(dens(nz)-dens(nz+1))
        ! do vertical  binning onto prescribed density classes
        is=std_dens_N
        do jj = 1, std_dens_N
           if (std_dens(jj) > dmin) then
              is = jj
              exit
           endif
        end do

        ie=1
        do jj = std_dens_N,1,-1
           if (std_dens(jj) < dmax) then
              ie = jj
              exit
           endif
        end do

        if (std_dens(is)>=dmax) is=ie
        if (std_dens(ie)<=dmin) ie=is
        if (Fer_GM) then
           uvdz_el=(UV(:,nz,elem)+fer_uv(:,nz,elem))*helem(nz,elem)
        else
           uvdz_el=UV(:,nz,elem)*helem(nz,elem)
        end if
        rhoz_el=(dens(nz)-dens(nz+1))/helem(nz,elem)
        vol_el =helem(nz,elem)*elem_area(elem)
        if (ie-is > 0) then
           weight=(std_dens(is)-dmin)+std_dd(is)/2.
           weight=max(weight, 0.)/ddiff
           std_dens_UVDZ(:, is, elem)=std_dens_UVDZ(:, is, elem)+weight*uvdz_el
           std_dens_VOL2(   is, elem)=std_dens_VOL2(   is, elem)+weight*vol_el
           locz=el_depth(nz+1)+weight*helem(nz,elem)
           std_dens_Z   (   is, elem)=std_dens_Z   (   is, elem)+locz*weight
           std_dens_w(      is, elem)=std_dens_w   (   is, elem)+weight
           do snz=is+1, ie-1
              weight=(sum(std_dd(snz-1:snz))/2.)/ddiff
              std_dens_UVDZ(:, snz, elem)=std_dens_UVDZ(:, snz, elem)+weight*uvdz_el
              std_dens_VOL2(   snz, elem)=std_dens_VOL2(   snz, elem)+weight*vol_el
              locz=locz+weight*helem(nz,elem)
              std_dens_Z   (   snz, elem)=std_dens_Z   (   snz, elem)+locz*weight
              std_dens_w   (   snz, elem)=std_dens_w   (   snz, elem)+weight
           end do
           weight=(dmax-std_dens(ie))+std_dd(ie-1)/2.
           weight=max(weight, 0.)/ddiff
           std_dens_UVDZ(:, ie, elem)=std_dens_UVDZ(:, ie, elem)+weight*uvdz_el
           std_dens_VOL2(   ie, elem)=std_dens_VOL2(   ie, elem)+weight*vol_el
           locz=locz+weight*helem(nz,elem)
           std_dens_Z   (   ie, elem)=std_dens_Z   (   ie, elem)+locz*weight
           std_dens_w   (   ie, elem)=std_dens_w   (   ie, elem)+weight
        else
           std_dens_UVDZ(:, is, elem)=std_dens_UVDZ(:, is, elem)+uvdz_el
           std_dens_VOL2(   is, elem)=std_dens_VOL2(   is, elem)+vol_el
           std_dens_Z   (   is, elem)=std_dens_Z   (   is, elem)+el_depth(nz+1)+helem(nz,elem)/2.
           std_dens_w   (   is, elem)=std_dens_w   (   is, elem)+1._wp
        end if
     end do
  end do

    !___________________________________________________________________________
    ! proceed with fields at nodes (cycle over edges to compute the divergence)...
    do edge=1, myDim_edge2D
        if (myList_edge2D(edge) > edge2D_in) cycle
        enodes=edges(:,edge)
        eelems=edge_tri(:,edge)
        nzmax =nlevels(eelems(1))
        nzmin =ulevels(eelems(1))
        if (eelems(2)>0) nzmax=max(nzmax, nlevels(eelems(2)))
        do nz=nzmin, nzmax-1
            aux(nz)=sum(density_dmoc(nz, enodes))/2.-1000.
        end do
        
        !_______________________________________________________________________
        do e=1,2
            elem=eelems(e)
            if (elem<=0) CYCLE
            deltaX=edge_cross_dxdy(1+(e-1)*2,edge) 
            deltaY=edge_cross_dxdy(2+(e-1)*2,edge)
            nzmax =nlevels(elem)
            nzmin =ulevels(elem)
            
            !___________________________________________________________________
            do nz=nzmax-1,nzmin+1,-1
                dens(nz)   = (aux(nz)     * helem(nz-1,elem)+&
                            aux(nz-1)   * helem(nz,  elem))/sum(helem(nz-1:nz,elem))
            end do
            dens(nzmax)=dens(nzmax-1)+(dens(nzmax-1)-dens(nzmax-2))*helem(nzmax-1,elem)/helem(nzmax-2,elem)
            dens(nzmin)    =dens(nzmin+1)      +(dens(nzmin+1)-dens(nzmin+2))            *helem(nzmin, elem)     /helem(nzmin+1,elem)       
            is=minloc(abs(std_dens-dens(nzmin)),1)
            
            !___________________________________________________________________
            do nz=nzmax-1,nzmin,-1
                div=(UV(2,nz,elem)*deltaX-UV(1,nz,elem)*deltaY)*helem(nz,elem)
                if (e==2) div=-div
                dmin =minval(dens(nz:nz+1))
                dmax =maxval(dens(nz:nz+1))
                ddiff=abs(dens(nz)-dens(nz+1))
                
                ! do vertical  binning onto prescribed density classes
                is=std_dens_N
                do jj = 1, std_dens_N
                    if (std_dens(jj) > dmin) then
                        is = jj
                        exit
                    endif
                end do
                ie=1
                do jj = std_dens_N,1,-1
                    if (std_dens(jj) < dmax) then
                        ie = jj
                        exit
                    endif
                end do
                
                if (std_dens(is)>=dmax) is=ie
                if (std_dens(ie)<=dmin) ie=is
                if (ie-is > 0) then
                    weight=(std_dens(is)-dmin)+std_dd(is)/2.
                    weight=max(weight, 0.)/ddiff
                    std_dens_DIV(is, enodes(1))=std_dens_DIV(is, enodes(1))+weight*div
                    std_dens_DIV(is, enodes(2))=std_dens_DIV(is, enodes(2))-weight*div
                    do snz=is+1, ie-1
                        weight=(sum(std_dd(snz-1:snz))/2.)/ddiff
                        std_dens_DIV(snz, enodes(1))=std_dens_DIV(snz, enodes(1))+weight*div
                        std_dens_DIV(snz, enodes(2))=std_dens_DIV(snz, enodes(2))-weight*div
                    end do
                    weight=(dmax-std_dens(ie))+std_dd(ie-1)/2.
                    weight=max(weight, 0.)/ddiff
                    std_dens_DIV(ie, enodes(1))=std_dens_DIV(ie, enodes(1))+weight*div
                    std_dens_DIV(ie, enodes(2))=std_dens_DIV(ie, enodes(2))-weight*div
                else
                    std_dens_DIV(is, enodes(1))=std_dens_DIV(is, enodes(1))+div
                    std_dens_DIV(is, enodes(2))=std_dens_DIV(is, enodes(2))-div
                end if ! --> if (ie-is > 0) then
            end do ! --> do nz=nzmax-1,nzmin,-1
            
            !___________________________________________________________________
            ! compute density class divergence from GM Bolus velocity
            if (Fer_GM) then
                do nz=nzmax-1,nzmin,-1
                    div=(fer_uv(2,nz,elem)*deltaX-fer_uv(1,nz,elem)*deltaY)*helem(nz,elem)
                    if (e==2) div=-div
                    dmin =minval(dens(nz:nz+1))
                    dmax =maxval(dens(nz:nz+1))
                    ddiff=abs(dens(nz)-dens(nz+1))
                    
                    ! do vertical  binning onto prescribed density classes
                    is=std_dens_N
                    do jj = 1, std_dens_N
                        if (std_dens(jj) > dmin) then
                            is = jj
                            exit
                        endif
                    end do
                    ie=1
                    do jj = std_dens_N,1,-1
                        if (std_dens(jj) < dmax) then
                            ie = jj
                            exit
                        endif
                    end do
                    
                    if (std_dens(is)>=dmax) is=ie
                    if (std_dens(ie)<=dmin) ie=is
                    if (ie-is > 0) then
                        weight=(std_dens(is)-dmin)+std_dd(is)/2.
                        weight=max(weight, 0.)/ddiff
                        std_dens_DIV_fer(is, enodes(1))=std_dens_DIV_fer(is, enodes(1))+weight*div
                        std_dens_DIV_fer(is, enodes(2))=std_dens_DIV_fer(is, enodes(2))-weight*div
                        do snz=is+1, ie-1
                            weight=(sum(std_dd(snz-1:snz))/2.)/ddiff
                            std_dens_DIV_fer(snz, enodes(1))=std_dens_DIV_fer(snz, enodes(1))+weight*div
                            std_dens_DIV_fer(snz, enodes(2))=std_dens_DIV_fer(snz, enodes(2))-weight*div
                        end do
                        weight=(dmax-std_dens(ie))+std_dd(ie-1)/2.
                        weight=max(weight, 0.)/ddiff
                        std_dens_DIV_fer(ie, enodes(1))=std_dens_DIV_fer(ie, enodes(1))+weight*div
                        std_dens_DIV_fer(ie, enodes(2))=std_dens_DIV_fer(ie, enodes(2))-weight*div
                    else
                        std_dens_DIV_fer(is, enodes(1))=std_dens_DIV_fer(is, enodes(1))+div
                        std_dens_DIV_fer(is, enodes(2))=std_dens_DIV_fer(is, enodes(2))-div
                    end if ! --> if (ie-is > 0) then
                end do ! --> do nz=nzmax-1,nzmin,-1
            end if ! --> if (Fer_GM) then
        end do ! --> do e=1,2
    end do ! --> do edge=1, myDim_edge2D
  
  !_____________________________________________________________________________
  where (std_dens_w > 0.)
        std_dens_Z   =std_dens_Z / std_dens_w
  end where
  
  !_____________________________________________________________________________
  ! compute density class volume change over time 
  if (.not. firstcall_e) then
     std_dens_dVdT=(std_dens_VOL2-std_dens_VOL1)/dt
  end if
  std_dens_VOL1=std_dens_VOL2
  
  !_____________________________________________________________________________
  ! compute mean thickness of density class, try to extract better vertical position
  ! when do projection into zcoord. 
  std_dens_H = std_dens_VOL2
  do jj = 1, std_dens_N
        std_dens_H(jj,1:myDim_elem2D) = std_dens_H(jj,1:myDim_elem2D)/elem_area(1:myDim_elem2D)
  end do
  
  firstcall_e=.false.
end subroutine diag_densMOC
!
!
!_______________________________________________________________________________
subroutine relative_vorticity(mode, dynamics, partit, mesh)
    IMPLICIT NONE
    integer        :: n, nz, el(2), enodes(2), nl1, nl2, edge, ul1, ul2, nl12, ul12
    real(kind=WP)  :: deltaX1, deltaY1, deltaX2, deltaY2, c1
    integer,        intent(in)            :: mode
    logical,        save                  :: firstcall=.true.
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    real(kind=WP), dimension(:,:,:), pointer :: UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    UV => dynamics%uv(:,:,:)
    
    !___________________________________________________________________________
    if (firstcall) then  !allocate the stuff at the first call
        allocate(vorticity(nl-1, myDim_nod2D+eDim_nod2D))
        firstcall=.false.
        if (mode==0) return
    end if
    !!PS DO n=1,myDim_nod2D
    !!PS    nl1 = nlevels_nod2D(n)-1
    !!PS    ul1 = ulevels_nod2D(n)
    !!PS    vorticity(ul1:nl1,n)=0.0_WP
    !!PS    !!PS DO nz=1, nlevels_nod2D(n)-1
    !!PS    !!PS    vorticity(nz,n)=0.0_WP
    !!PS    !!PS END DO
    !!PS END DO      
    vorticity = 0.0_WP
    DO edge=1,myDim_edge2D
                                    !! edge=myList_edge2D(m)
        enodes=edges(:,edge)
        el=edge_tri(:,edge)
        nl1=nlevels(el(1))-1
        ul1=ulevels(el(1))
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        nl2=0
        ul2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            nl2=nlevels(el(2))-1
            ul2=ulevels(el(2))
        end if  
        nl12 = min(nl1,nl2)
        ul12 = max(ul1,ul2)
        
        DO nz=ul1,ul12-1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        if (ul2>0) then
            DO nz=ul2,ul12-1
                c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
                vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
            END DO
        endif 
        !!PS DO nz=1,min(nl1,nl2)
        DO nz=ul12,nl12
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
            deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        !!PS DO nz=min(nl1,nl2)+1,nl1
        DO nz=nl12+1,nl1
            c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
        !!PS DO nz=min(nl1,nl2)+1,nl2
        DO nz=nl12+1,nl2
            c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
            vorticity(nz,enodes(1))=vorticity(nz,enodes(1))+c1
            vorticity(nz,enodes(2))=vorticity(nz,enodes(2))-c1
        END DO
    END DO
    
    ! vorticity = vorticity*area at this stage
    ! It is correct only on myDim nodes
    DO n=1,myDim_nod2D
                                !! n=myList_nod2D(m)
        ul1 = ulevels_nod2D(n)
        nl1 = nlevels_nod2D(n)
        !!PS DO nz=1,nlevels_nod2D(n)-1
        DO nz=ul1,nl1-1
            vorticity(nz,n)=vorticity(nz,n)/areasvol(nz,n)
        END DO
    END DO      
    
    call exchange_nod(vorticity, partit)
    
! Now it the relative vorticity known on neighbors too
end subroutine relative_vorticity
!
!
!_______________________________________________________________________________
subroutine compute_extflds(mode, dynamics, tracers, partit, mesh)
    IMPLICIT NONE
    integer,        intent(in)              :: mode
    logical,        save                    :: firstcall=.true.
    type(t_dyn)   , intent(in),     target  :: dynamics
    type(t_tracer), intent(in)   ,  target  :: tracers
    type(t_partit), intent(inout),  target  :: partit
    type(t_mesh)  , intent(in)   ,  target  :: mesh
    real(kind=WP),  dimension(:,:), pointer :: temp, salt
    real(kind=WP)                           :: zn, zint, tup, tlo
    integer                                 :: n, nz, nzmin, nzmax
    real(kind=WP)                           :: whichtemp=  20.0_WP ! which isotherm to compute (set 20 per default)
    real(kind=WP)                           :: whichdepth=300.0_WP ! for which tepth to average for tempzavg & saltzavg

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    if (firstcall) then  !allocate the stuff at the first call
        allocate(zisotherm(myDim_nod2D+eDim_nod2D))
        allocate(tempzavg(myDim_nod2D+eDim_nod2D), saltzavg(myDim_nod2D+eDim_nod2D))
        zisotherm=0.0_WP
        tempzavg =0.0_WP
        saltzavg =0.0_WP
        firstcall=.false.
        if (mode==0) return
    end if  
    temp   => tracers%data(1)%values(:,:)
    salt   => tracers%data(2)%values(:,:)    

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, zn, tup, tlo)
    DO n=1, myDim_nod2D
       saltzavg(n) =0.0_WP
       nzmax=nlevels_nod2D(n)
       nzmin=ulevels_nod2D(n)
       zn  =0.0_WP
       do nz=nzmin+1, nzmax-1
          tup=temp(nz-1, n)
          if (tup < whichtemp) exit
          tlo=temp(nz,   n)                 
          if ((tup-whichtemp)*(tlo-whichtemp)<0) then
             zn=zn+0.5_WP*(hnode(nz-1, n)+(whichtemp-tup)*sum(hnode(nz-1:nz, n))/(tlo-tup))
             exit
          end if
          zn=zn+hnode(nz-1, n) 
       end do
       zisotherm(n)=zn
    END DO
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, nzmin, nzmax, zint)
    DO n=1, myDim_nod2D
       tempzavg(n) =0.0_WP
       saltzavg(n) =0.0_WP
       nzmax=nlevels_nod2D(n)
       nzmin=ulevels_nod2D(n)
       zint=0.0_WP
       do nz=nzmin, nzmax-1
          zint=zint+hnode(nz, n) 
          tempzavg(n)=tempzavg(n)+temp(nz,   n)*hnode(nz, n)
          saltzavg(n)=saltzavg(n)+salt(nz,   n)*hnode(nz, n)
          if (zint>=whichdepth) exit
       end do
       tempzavg(n)=tempzavg(n)/zint
       saltzavg(n)=saltzavg(n)/zint
    END DO
!$OMP END PARALLEL DO 

  call exchange_nod(zisotherm, partit)
  call exchange_nod(tempzavg, partit)
  call exchange_nod(saltzavg, partit)
end subroutine compute_extflds
!_______________________________________________________________________________
subroutine compute_ice_diag(mode, ice, partit, mesh)
    IMPLICIT NONE
    integer,        intent(in)              :: mode
    logical,        save                    :: firstcall=.true.
    type(t_ice)   , intent(in),     target  :: ice
    type(t_partit), intent(inout),  target  :: partit
    type(t_mesh)  , intent(in)   ,  target  :: mesh
    integer                                 :: n

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    if (firstcall) then  !allocate the stuff at the first call
        allocate(vol_ice(myDim_nod2D+eDim_nod2D), vol_snow(myDim_nod2D+eDim_nod2D))
        vol_ice =0.0_WP
        vol_snow=0.0_WP
        firstcall=.false.
        if (mode==0) return
    end if
!$OMP PARALLEL DO
DO n=1, myDim_nod2D
   vol_ice(n) =ice%data(1)%values(n)*ice%data(2)%values(n)
   vol_snow(n)=ice%data(1)%values(n)*ice%data(3)%values(n)
END DO
!$OMP END PARALLEL DO

end subroutine compute_ice_diag

! SST in K
subroutine compute_thetao(mode, tracers, partit, mesh)
  implicit none
  integer,        intent(in)            :: mode
  type(t_tracer), intent(in) ,  target  :: tracers
  type(t_mesh)  , intent(in) ,  target  :: mesh
  type(t_partit), intent(in), target :: partit
  logical, save                         :: firstcall=.true.
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (firstcall) then !allocate the stuff at the first call
     allocate(thetao(mydim_nod2D))
     firstcall=.false.
     if (mode==0) return
  end if

  !skipping loop 
  thetao(:) = tracers%data(1)%values(1,1:myDim_nod2D)+273.15_WP 
end subroutine compute_thetao

! ==============================================================
subroutine compute_diagnostics(mode, dynamics, tracers, ice, partit, mesh)
  implicit none
  type(t_mesh)  , intent(in)   , target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(inout), target :: tracers
  type(t_ice),    intent(inout), target :: ice
  type(t_dyn)   , intent(inout), target :: dynamics
  integer, intent(in)                   :: mode !constructor mode (0=only allocation; any other=do diagnostic)
  real(kind=WP)                         :: val  !1. solver diagnostic
  if (ldiag_solver)      call diag_solver(mode, dynamics, partit, mesh)
  !2. compute curl(stress_surf)
  if (lcurt_stress_surf) call diag_curl_stress_surf(mode, partit, mesh)
  !3. compute curl(velocity)
  if (ldiag_curl_vel3)   call diag_curl_vel3(mode, dynamics, partit, mesh)
  !4. compute energy budget
  if (ldiag_Ri)          call diag_Ri(mode, dynamics, partit, mesh)
  !5. print integrated temperature 
  if (ldiag_salt3d) then
     if (mod(mstep,logfile_outfreq)==0) then
        call integrate_nod(tracers%data(2)%values(:,:), val, partit, mesh)
        if (partit%mype==0) then
           write(*,*) 'total integral of salinity at timestep :', mstep, val
        end if
     end if
  end if
  !6. MOC in density coordinate
  if (ldiag_dMOC)        call diag_densMOC(mode, dynamics, tracers, partit, mesh)
  !7. compute turbulent fluxes
  if (ldiag_turbflux)    call diag_turbflux(mode, dynamics, tracers, partit, mesh)
  !8. compute tracers fluxes
  if (ldiag_trflx)       call diag_trflx(mode, dynamics, tracers, partit, mesh)
  ! compute relative vorticity
  if (ldiag_vorticity)   call relative_vorticity(mode, dynamics, partit, mesh)
  ! soe exchanged fields requested by IFS/FESOM in NextGEMS.
  if (ldiag_extflds)     call compute_extflds(mode, dynamics, tracers, partit, mesh)
  !fields required for for destinE
  if (ldiag_ice)         call compute_ice_diag(mode, ice, partit, mesh)
  call compute_thetao(mode, tracers, partit, mesh) 
end subroutine compute_diagnostics

!
!
!_______________________________________________________________________________
! calculate horizintal and vertical advection for squared tracer (2nd moments)
! see: 
! Burchard and Rennau, 2008, Comparative quantification of physically and 
!                      numerically induced mixing in ocean models ...
! Rennau, and Burchard, 2009, Quantitative analysis of numerically induced mixing 
!                      in a coastal model application ...
! Klingbeil et al., 2014, Quantification of spurious dissipation and mixing – 
!                      Discrete variance decay in a Finite-Volume framework ...
subroutine compute_diag_dvd_2ndmoment_burchard_etal_2008(tr_num, tracers, partit, mesh)
    use o_arrays
    use oce_adv_tra_driver_interfaces    
    implicit none
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(inout), target :: tracers
    integer,        intent(in)            :: tr_num
    integer                  :: node, nz, nzmin, nzmax
    real(kind=WP)            :: tr_sqr(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D), trAB_sqr(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
    !___________________________________________________________________________
    ! square up fields for actual tracers and Adams Bashfort tracer
    ! --> dont forget to square also the halo !!! --> that why do node = ...+eDim_nod2D
    tr_sqr   = 0.0_WP
    trAB_sqr = 0.0_WP
    do node = 1, myDim_nod2D+eDim_nod2D
        !!PS do nz = 1, nlevels_nod2D(node)-1
        nzmax = nlevels_nod2D(node)-1
        nzmin = ulevels_nod2D(node)
        do nz = nzmin, nzmax
            tr_sqr(nz,node)   = tracers%data(tr_num)%values  (nz,node)**2
            trAB_sqr(nz,node) = tracers%data(tr_num)%valuesAB(nz,node)**2
        end do
    end do
        
    !___________________________________________________________________________
    ! calculate horizintal and vertical advection for squared tracer (2nd moments)
    ! see Burchard and Rennau, 2008, Comparative quantification of physically and 
    ! numerically induced mixing in ocean models ...
    tracers%work%del_ttf_advhoriz = 0.0_WP
    tracers%work%del_ttf_advvert  = 0.0_WP
!   maybe just to introduce an another tharer of t_tracer type with **do_Xmoment?
!   call do_oce_adv_tra(dt, UV, wvel, wvel_i, wvel_e, tr_sqr, trAB_sqr, 1, tracers%work%del_ttf_advhoriz, tracers%work%del_ttf_advvert, tra_adv_ph, tra_adv_pv, partit, mesh)   
    !___________________________________________________________________________
    ! add target second moment to DVD
    do node = 1,mydim_nod2D
        !!PS do nz = 1,nlevels_nod2D(node)-1
        nzmax = nlevels_nod2D(node)-1
        nzmin = ulevels_nod2D(node)
        do nz = nzmin, nzmax
            ! eq 16 & 17 Klingbeil et al. 2014
            !
            ! (phi^2)^(n+1) = 1/V^(n+1)*[ V^(n)*(phi^2)^(n) + dt*ADV[phi^2]  ]
            !
            !  DVD = 1/dt * [ (phi^2)^(n+1) -  ( phi^(n+1) )^2 ]
            !                     |
            !                     v
            !                first this part
            ! --> split it up in DVD contribution from horizontal and vertical 
            ! advection since for the horizontal advection Adams Bashfort tracer 
            ! are used and for the vertical the normal tracer values.
            tracers%work%tr_dvd_horiz(nz,node,tr_num) = hnode(nz,node)/hnode_new(nz,node)*trAB_sqr(nz,node) - tracers%work%del_ttf_advhoriz(nz,node)/hnode_new(nz,node)
            tracers%work%tr_dvd_vert (nz,node,tr_num) = hnode(nz,node)/hnode_new(nz,node)*tr_sqr(  nz,node) - tracers%work%del_ttf_advvert( nz,node)/hnode_new(nz,node)
        end do
    end do
end subroutine compute_diag_dvd_2ndmoment_burchard_etal_2008
!
!
!
!
!_______________________________________________________________________________
! calculate horizintal and vertical advection for squared tracer (2nd moments)
! see: 
! Klingbeil et al., 2014, Quantification of spurious dissipation and mixing – 
!                      Discrete variance decay in a Finite-Volume framework ...
subroutine compute_diag_dvd_2ndmoment_klingbeil_etal_2014(tr_num, tracers, partit, mesh)
    use o_arrays
    use oce_adv_tra_driver_interfaces
    implicit none
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(inout), target :: tracers
    integer                               :: node, nz, nzmin, nzmax
    integer,        intent(in)            :: tr_num

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    !___________________________________________________________________________
    ! calculate horizintal and vertical advection for squared tracer (2nd moments)
    ! see Burchard and Rennau, 2008, Comparative quantification of physically and 
    ! numerically induced mixing in ocean models ...
    tracers%work%del_ttf_advhoriz = 0.0_WP
    tracers%work%del_ttf_advvert  = 0.0_WP
!   maybe just to introduce an another tharer of t_tracer type with **do_Xmoment?
!   call do_oce_adv_tra(dt, UV, wvel, wvel_i, wvel_e, tracers%data(tr_num)%values, tracers%data(tr_num)%valuesAB(:,:), 2, tracers%work%del_ttf_advhoriz, tracers%work%del_ttf_advvert, tra_adv_ph, tra_adv_pv, partit, mesh)   
    !___________________________________________________________________________
    ! add target second moment to DVD
    do node = 1,mydim_nod2D
        !!PS do nz = 1,nlevels_nod2D(node)-1
        nzmax = nlevels_nod2D(node)-1
        nzmin = ulevels_nod2D(node)
        do nz = nzmin, nzmax
            ! eq 23 Klingbeil et al. 2014
            !
            ! phi^(n+1) = 1/V^(n+1)*[ V^(n)*phi^(n) + dt*ADV[phi]  ]
            !
            !  DVD = -1/dt * [ (phi^(n+1))^2 -  V^n/V^n+1*(phi^n)^2 + dt/V^(n+1)*ADV[phi^2]]
            !
            !      = 1/dt * [-(V^n/V^n+1*phi^n - dt/V^(n+1)*ADV[phi])^2 
            !                + V^n/V^n+1*(phi^n)^2 - dt/V^(n+1)*ADV[phi^2] ]
            !                  \_________________________________________/
            !                                      |
            !                                      v
            !                                first this part
            ! --> split it up in DVD contribution from horizontal and vertical 
            ! advection since for the horizontal advection Adams Bashfort tracer 
            ! are used and for the vertical the normal tracer values.
            tracers%work%tr_dvd_horiz(nz,node,tr_num) = hnode(nz,node)/hnode_new(nz,node)*(tracers%data(tr_num)%valuesAB(nz,node)**2) &
                                           - tracers%work%del_ttf_advhoriz(nz,node)/hnode_new(nz,node)
            tracers%work%tr_dvd_vert(nz,node,tr_num)  = hnode(nz,node)/hnode_new(nz,node)*(tracers%data(tr_num)%values  (nz,node)**2) &
                                           - tracers%work%del_ttf_advvert( nz,node)/hnode_new(nz,node)
        end do
    end do
end subroutine compute_diag_dvd_2ndmoment_klingbeil_etal_2014
!
!
!_______________________________________________________________________________
! calculate horizintal and vertical advection for squared tracer (2nd moments)
! see: 
! Burchard and Rennau, 2008, Comparative quantification of physically and 
!                      numerically induced mixing in ocean models ...
! Rennau, and Burchard, 2009, Quantitative analysis of numerically induced mixing 
!                      in a coastal model application ...
! Klingbeil et al., 2014, Quantification of spurious dissipation and mixing – 
!                      Discrete variance decay in a Finite-Volume framework ...
subroutine compute_diag_dvd(tr_num, tracers, partit, mesh)
    use g_config, only: dt
    use o_arrays   
    implicit none
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(inout), target :: tracers
    integer                               :: node, nz, nzmin, nzmax
    integer,        intent(in)            :: tr_num


#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    !___________________________________________________________________________
    ! add discret second moment to DVD
    do node = 1,mydim_nod2D
        !!PS do nz = 1,nlevels_nod2D(node)-1
        nzmax = nlevels_nod2D(node)-1
        nzmin = ulevels_nod2D(node)
        do nz = nzmin, nzmax
            ! eq 16 & 17  and eq 23. Klingbeil et al. 2014
            !
            ! (phi^2)^(n+1) = 1/V^(n+1)*[ V^(n)*(phi^2)^(n) + dt*ADV[phi^2]  ]
            !
            !  DVD = 1/dt * [ (phi^2)^(n+1) -  ( phi^(n+1) )^2 ]
            !                                         |
            !                                         v
            !                                    now add this part
            ! --> tracers%work%tr_dvd_horiz contains already the expected target second moments
            ! from subroutine compute_diag_dvd_2ndmoment
            tracers%work%tr_dvd_horiz(nz,node,tr_num) = (tracers%work%tr_dvd_horiz(nz,node,tr_num)                                     &
                                            -( hnode(nz,node)/hnode_new(nz,node)*tracers%data(tr_num)%valuesAB(nz,node)          &
                                              -tracers%work%del_ttf_advhoriz(nz,node)/hnode_new(nz,node)                        &
                                              )**2                                                                 &
                                            )/dt
            tracers%work%tr_dvd_vert(nz,node,tr_num)  = (tracers%work%tr_dvd_vert(nz,node,tr_num)                                      &
                                            -( hnode(nz,node)/hnode_new(nz,node)*tracers%data(tr_num)%values  (nz,node)          &
                                              -tracers%work%del_ttf_advvert( nz,node)/hnode_new(nz,node)                        &
                                              )**2                                                                 &
                                            )/dt
        end do
    end do
end subroutine compute_diag_dvd

end module diagnostics
