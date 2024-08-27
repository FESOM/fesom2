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
  use o_tracers
  use g_forcing_arrays
  use o_mixing_KPP_mod
  use g_rotate_grid
  use g_support
  use Toy_Channel_Soufflet
  implicit none

  private
  public :: ldiag_solver, lcurt_stress_surf, ldiag_Ri, ldiag_TurbFlux, ldiag_dMOC, ldiag_DVD,        &
            ldiag_forc, ldiag_salt3D, ldiag_curl_vel3, diag_list, ldiag_vorticity, ldiag_extflds, ldiag_ice,   &
            compute_diagnostics, rhs_diag, curl_stress_surf, curl_vel3, shear, Ri, KvdTdZ, KvdSdZ,   & 
            std_dens_min, std_dens_max, std_dens_N, std_dens, ldiag_trflx,                           &
            std_dens_UVDZ, std_dens_DIV, std_dens_DIV_fer, std_dens_Z, std_dens_H, std_dens_dVdT, std_dens_flux,       &
            dens_flux_e, vorticity, zisotherm, tempzavg, saltzavg, vol_ice, vol_snow, compute_ice_diag, thetao,   &
            tuv, suv, compute_dvd, dvd_KK_tot, dvd_SD_tot, dvd_SD_chi_adv_h, dvd_SD_chi_adv_v, dvd_SD_chi_dif_he, &
            dvd_SD_chi_dif_heR, dvd_SD_chi_dif_hbh, dvd_SD_chi_dif_veR, dvd_SD_chi_dif_viR, dvd_SD_chi_dif_vi,    &
            dvd_SD_chi_dif_ve, dvd_SD_chi_dif_sbc, dvd_xdfac

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

  !_____________________________________________________________________________
  ! DVD diagnostics
  real(kind=WP),  save, allocatable, target      :: dvd_KK_tot(:,:,:), dvd_SD_tot(:,:,:), dvd_SD_chi_adv_h(:,:,:), &
                                                    dvd_SD_chi_adv_v( :,:,:), dvd_SD_chi_dif_heR(:,:,:), dvd_SD_chi_dif_veR(:,:,:), &
                                                    dvd_SD_chi_dif_viR(:,:,:), dvd_SD_chi_dif_vi(:,:,:), dvd_SD_chi_dif_hbh(:,:,:), &
                                                    dvd_SD_chi_dif_ve(:,:,:), dvd_SD_chi_dif_he(:,:,:), dvd_SD_chi_dif_sbc(:,:,:), trstar(:,:)
  real(kind=WP),  parameter                      :: dvd_xdfac=0.5_WP  ! Xchi distribution factor, default distribute 
                                                                      ! equal amount (50:50) of xchi on both side of face
  !_____________________________________________________________________________
  ! Define diagnostic flags + with corresponding namelist
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
  
  namelist /diag_list/ ldiag_solver, lcurt_stress_surf, ldiag_curl_vel3, ldiag_Ri, & 
                       ldiag_TurbFlux, ldiag_dMOC, ldiag_DVD, ldiag_salt3D, ldiag_forc, &
                       ldiag_vorticity, ldiag_extflds, ldiag_trflx, ldiag_ice
  
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

    !___________________________________________________________________________
    if (firstcall) then  !allocate the stuff at the first call
        allocate(curl_vel3(nl-1, myDim_nod2D+eDim_nod2D))
        firstcall=.false.
        if (mode==0) return
    end if

    !___________________________________________________________________________
    curl_vel3=0.
    do ed=1,myDim_edge2D
        enodes=edges(:,ed)
        el=edge_tri(:,ed)
        !_______________________________________________________________________
        nl1=nlevels(el(1))-1
        nu1=ulevels(el(1))
        deltaX1=edge_cross_dxdy(1,ed)
        deltaY1=edge_cross_dxdy(2,ed)
        
        !_______________________________________________________________________
        if (el(2)>0) then
            deltaX2=edge_cross_dxdy(3,ed)
            deltaY2=edge_cross_dxdy(4,ed)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            nl12 = min(nl1,nl2)
            nu12 = max(nu1,nu2)
            !___________________________________________________________________
            do nz=nu1,nu12-1
                c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
            do nz=nu2,nu12-1
                c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
            !___________________________________________________________________
            do nz=nu12,nl12
                c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))- &
                deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
            !___________________________________________________________________
            do nz=nl12+1,nl1
                c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
            do nz=nl12+1,nl2
                c1= -deltaX2*UV(1,nz,el(2))-deltaY2*UV(2,nz,el(2))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
        !_______________________________________________________________________    
        else
            do nz=nu1,nl1
                c1=deltaX1*UV(1,nz,el(1))+deltaY1*UV(2,nz,el(1))
                curl_vel3(nz,enodes(1))=curl_vel3(nz,enodes(1))+c1
                curl_vel3(nz,enodes(2))=curl_vel3(nz,enodes(2))-c1
            end do
        end if
    end do
    
    !___________________________________________________________________________
    do n=1, myDim_nod2D
        do nz=ulevels_nod2D(n), nlevels_nod2D(n)-1
            curl_vel3(nz,n)=curl_vel3(nz,n)/areasvol(nz,n)
        end do
    end do
    
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
  ! 1. compute solver diagnostic
  if (ldiag_solver)      call diag_solver(mode, dynamics, partit, mesh)
  
  ! 2. compute curl(stress_surf)
  if (lcurt_stress_surf) call diag_curl_stress_surf(mode, partit, mesh)
  
  ! 3. compute curl(velocity)
  if (ldiag_curl_vel3)   call diag_curl_vel3(mode, dynamics, partit, mesh)
  
  ! 4. compute energy budget
  if (ldiag_Ri)          call diag_Ri(mode, dynamics, partit, mesh)
  
  ! 5. print integrated temperature 
  if (ldiag_salt3d) then
     if (mod(mstep,logfile_outfreq)==0) then
        call integrate_nod(tracers%data(2)%values(:,:), val, partit, mesh)
        if (partit%mype==0) then
           write(*,*) 'total integral of salinity at timestep :', mstep, val
        end if
     end if
  end if
  
  ! 6. MOC in density coordinate
  if (ldiag_dMOC)        call diag_densMOC(mode, dynamics, tracers, partit, mesh)
  
  ! 7. compute turbulent fluxes
  if (ldiag_turbflux)    call diag_turbflux(mode, dynamics, tracers, partit, mesh)
  
  ! 8. compute tracers fluxes
  if (ldiag_trflx)       call diag_trflx(mode, dynamics, tracers, partit, mesh)

  ! 9. compute relative vorticity
  if (ldiag_vorticity)   call relative_vorticity(mode, dynamics, partit, mesh)
  
  ! 10. compute some exchanged fields requested by IFS/FESOM in NextGEMS.
  if (ldiag_extflds)     call compute_extflds(mode, dynamics, tracers, partit, mesh)

  ! 11. compute Discrete Variance Decay (DVD) diagnostic of: 
  ! K. Klingbeil et al., 2014, Quantification of spurious dissipation and mixing – 
  ! Discrete variance decay in Finite-Volume framework ...
  ! T. Banerjee, S. Danilov, K. Klingbeil, 2023, Discrete variance decay analysis of 
  ! spurious mixing
  if (ldiag_DVD)     call compute_dvd(mode, dynamics, tracers, partit, mesh)
 
  !fields required for for destinE
  if (ldiag_ice)         call compute_ice_diag(mode, ice, partit, mesh)
  
  call compute_thetao(mode, tracers, partit, mesh) 

end subroutine compute_diagnostics


!
!
!_______________________________________________________________________________
! compute spurious mixing diagnostic dvd 
subroutine compute_dvd(mode, dynamics, tracers, partit, mesh)
    IMPLICIT NONE
    integer,        intent(in)                :: mode
    logical,        save                      :: firstcall=.true.
    type(t_dyn)   , intent(in),     target    :: dynamics
    type(t_tracer), intent(in)   ,  target    :: tracers
    type(t_partit), intent(inout),  target    :: partit
    type(t_mesh)  , intent(in)   ,  target    :: mesh
    integer                                   :: tr_num, node, elem
    real(kind=WP),  dimension(:,:)  , pointer :: trflx_h, trflx_v, tr, trold, Wvel, dump, fer_Wvel
    real(kind=WP),  dimension(:,:,:), pointer :: UV, fer_UV
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    if (firstcall) then  !allocate the stuff at the first call
        allocate(trstar(                nl-1, myDim_nod2D+eDim_nod2D))
        allocate(dvd_KK_tot(            nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_tot(            nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_chi_adv_h(      nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_chi_adv_v(      nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_chi_dif_he(     nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_chi_dif_vi(     nl-1, myDim_nod2D+eDim_nod2D, 2))
        allocate(dvd_SD_chi_dif_sbc(    nl-1, myDim_nod2D+eDim_nod2D, 2))
        if (tracers%data(1)%smooth_bh_tra) then
            allocate(dvd_SD_chi_dif_hbh(nl-1, myDim_nod2D+eDim_nod2D, 2))
        end if 
        if (Redi) then
            allocate(dvd_SD_chi_dif_heR(nl-1, myDim_nod2D+eDim_nod2D, 2))
            allocate(dvd_SD_chi_dif_veR(nl-1, myDim_nod2D+eDim_nod2D, 2))
            allocate(dvd_SD_chi_dif_viR(nl-1, myDim_nod2D+eDim_nod2D, 2))
        end if 
        if (.not. tracers%data(1)%i_vert_diff) then
            allocate(dvd_SD_chi_dif_ve( nl-1, myDim_nod2D+eDim_nod2D, 2))
        end if 
        firstcall=.false.
        if (mode==0) return
    end if  
    
    !___________________________________________________________________________
    ! initialise each time diagnostic is computed
    trstar            = 0.0_WP
    dvd_KK_tot        = 0.0_WP ! --> DVD diagnostic after Klingbeil et al. 2014
    dvd_SD_tot        = 0.0_WP ! --> DVD diagnostic after Banjerjee et al. 2023 (Sergeys way!!!)
    dvd_SD_chi_adv_h  = 0.0_WP
    dvd_SD_chi_adv_v  = 0.0_WP
    dvd_SD_chi_dif_he = 0.0_WP 
    dvd_SD_chi_dif_vi = 0.0_WP ! implicite part
    dvd_SD_chi_dif_sbc= 0.0_WP ! implicite part
    if (tracers%data(1)%smooth_bh_tra) then
        dvd_SD_chi_dif_hbh= 0.0_WP 
    end if 
    if (Redi) then
        dvd_SD_chi_dif_heR= 0.0_WP 
        dvd_SD_chi_dif_veR= 0.0_WP ! explicite part (in case of Redi==.true.)
        dvd_SD_chi_dif_viR= 0.0_WP 
    end if 
    if (.not. tracers%data(1)%i_vert_diff) then
        dvd_SD_chi_dif_ve = 0.0_WP 
    end if 
    
    !___________________________________________________________________________
    UV     => dynamics%uv(:,:,:)
    Wvel   => dynamics%w(:,:)
    if (Fer_GM) then
        fer_UV     => dynamics%fer_uv(:,:,:)
        fer_Wvel   => dynamics%fer_w(:,:)
    end if
    
    !___________________________________________________________________________
    ! update 3D velocities with the bolus velocities:
    ! 1. bolus velocities are computed according to GM implementation after R. Ferrari et al., 2010
    ! 2. bolus velocities are used only for advecting tracers and shall be subtracted back afterwards
    if (Fer_GM) then
!$OMP PARALLEL DO
        do elem=1, myDim_elem2D+eDim_elem2D
           UV(:, :, elem) = UV(:, :, elem) + fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
            Wvel(:, node) = Wvel(:, node) + fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    ! loop over temp and salt tracer
    do tr_num=1,2 
        !_______________________________________________________________________
        ! create pointer to reconstructed tracer fluxes through mid edge face, and 
        ! upper lower cell prism face and for tracer^n and tracer^(n+1)
        trflx_h   => tracers%work%dvd_trflx_hor( :,:,tr_num) ! horizontal advectiv tracer flux through mid edge face using AB tracer
        trflx_v   => tracers%work%dvd_trflx_ver( :,:,tr_num) ! vertical advectiv tracer flux 
        trold     => tracers%data(tr_num)%valuesold(1, :, :) ! tracer^n
        tr        => tracers%data(tr_num)%values(:,:)        ! tracer^(n+1)
        dump      => tracers%work%del_ttf
        
        !_______________________________________________________________________
        ! need to recompute tracer gradients 
        call tracer_gradient_elements(trold, partit, mesh)
        call exchange_elem(tr_xy, partit)
        call tracer_gradient_z(trold, partit, mesh)
        call exchange_nod(tr_z, partit)
        
        !=== DVD Knut Klingbeil et al. 2014 ====================================
        ! add time derivativ of 2nd. moment tracer
        ! --> at this point trstar corresponds to tr_old
        call dvd_add_time_deriv(tr_num, dvd_KK_tot, trold, tr, partit, mesh)            
        
        ! compute horizontal 2nd moment tracer flux from advective tracer fluxes at mid 
        ! edge face
        call dvd_add_advflux_hor( .false., tr_num, dvd_KK_tot, trflx_h, UV, trold, dump, partit, mesh)
        
        ! add vertical 2nd moment tracer flux at upper/lower scalar cell prism interface
        call dvd_add_advflux_ver( .false., tr_num, dvd_KK_tot, trflx_v, Wvel, trold, partit, mesh)
        
        ! add contribution from horizontal diffusion flux (after klingbeil et al. 2014)
        ! --> keep in mind here trstar corresponds to tr_old
        call dvd_add_difflux_horexpl( .false., tr_num, dvd_KK_tot, trold, Ki, tr_xy, dump, partit, mesh)
        
        ! add contribution from vertical diffusion flux (after klingbeil et al. 2014)
        if (.not. tracers%data(tr_num)%i_vert_diff) then
            call dvd_add_difflux_vertexpl(.false., tr_num, dvd_KK_tot, tr, trold, Kv, partit, mesh)
        end if     
        if (Redi) then
            call dvd_add_difflux_horexplredi( .false., tr_num, dvd_KK_tot, trold, Ki, slope_tapered, tr_z, dump, partit, mesh)
            call dvd_add_difflux_vertexplredi(.false., tr_num, dvd_KK_tot, trold, Ki, slope_tapered, tr_xy, partit, mesh)
            call dvd_add_difflux_vertimplredi(.false., tr_num, dvd_KK_tot, tr, trold, Ki, slope_tapered, partit, mesh)
        end if     
        
        call dvd_add_difflux_vertimpl(.false., tr_num, dvd_KK_tot, tr, trold, Kv, partit, mesh)
        
        call dvd_add_difflux_sbc(.false., tr_num, dvd_KK_tot, tr, trold, partit, mesh)
        
        ! add contribution from horizontal biharmonic diffusion flux if applied
        if (tracers%data(tr_num)%smooth_bh_tra) then
            call dvd_add_difflux_bhvisc(.false., tr_num, dvd_KK_tot, tr, trold,& 
            tracers%data(tr_num)%gamma0_tra, tracers%data(tr_num)%gamma1_tra, &
            tracers%data(tr_num)%gamma2_tra, dump, partit, mesh)
        end if 
        
        ! add contribution from climatological 3d restoring 
        if ((toy_ocean) .AND. ((tr_num==1) .AND. (TRIM(which_toy)=="soufflet"))) then
            call dvd_add_clim_relax_channel(.false., tr_num, dvd_KK_tot, trold, partit, mesh)
        elseif (clim_relax>1.0e-8_WP) then 
            call dvd_add_clim_relax(        .false., tr_num, dvd_KK_tot, trold, partit, mesh)
        end if 
        
        
        !
        !=== DVD Sergey Danilov after T. Banerjee et al. 2023 =======================
        ! from here on compute Tstar = ( T^(n+1) + T^(n) )*0.5
        trstar  = (trold + tr)*0.5_WP
        
        ! add contribution from horizontal advection
        call dvd_add_advflux_hor( .true., tr_num, dvd_SD_chi_adv_h, trflx_h, UV, trstar, dump, partit, mesh)
        
        ! add contribution from vertical advection
        call dvd_add_advflux_ver( .true., tr_num, dvd_SD_chi_adv_v, trflx_v, Wvel, trstar,  partit, mesh)
        
        ! add contribution from horizontal diffusion
        call dvd_add_difflux_horexpl( .true., tr_num, dvd_SD_chi_dif_he, trstar, Ki, tr_xy, dump, partit, mesh)
        
        ! add contribution from vertical diffusion
        if (.not. tracers%data(tr_num)%i_vert_diff) then 
            call dvd_add_difflux_vertexpl(.true., tr_num, dvd_SD_chi_dif_ve, tr, trstar, Kv, partit, mesh)
        end if     
        if (Redi) then
            call dvd_add_difflux_horexplredi( .true., tr_num, dvd_SD_chi_dif_heR, trstar, Ki, slope_tapered, tr_z, dump, partit, mesh)
            call dvd_add_difflux_vertexplredi(.true., tr_num, dvd_SD_chi_dif_veR, trstar, Ki, slope_tapered, tr_xy, partit, mesh)
            call dvd_add_difflux_vertimplredi(.true., tr_num, dvd_SD_chi_dif_viR, tr, trstar, Ki, slope_tapered, partit, mesh)
        end if     
        
        call dvd_add_difflux_vertimpl(.true., tr_num, dvd_SD_chi_dif_vi, tr, trstar, Kv, partit, mesh)
        
        call dvd_add_difflux_sbc(.true., tr_num, dvd_SD_chi_dif_sbc, tr, trstar, partit, mesh)
        
        ! add contribution from horizontal biharmonic diffusion flux if applied
        if (tracers%data(tr_num)%smooth_bh_tra) then
            call dvd_add_difflux_bhvisc(.true., tr_num, dvd_SD_chi_dif_hbh, tr, trstar,& 
            tracers%data(tr_num)%gamma0_tra, tracers%data(tr_num)%gamma1_tra, &
            tracers%data(tr_num)%gamma2_tra, dump, partit, mesh)
        end if 
        
        ! compute total Xchi 
        dvd_SD_tot(:,:,tr_num) = dvd_SD_chi_adv_h(  :,:,tr_num) + dvd_SD_chi_adv_v( :,:,tr_num) + &
                                 dvd_SD_chi_dif_vi( :,:,tr_num) + dvd_SD_chi_dif_he(:,:,tr_num) + &
                                 dvd_SD_chi_dif_sbc(:,:,tr_num)
                                 
        if (.not. tracers%data(tr_num)%i_vert_diff) then 
            dvd_SD_tot(:,:,tr_num) = dvd_SD_tot(:,:,tr_num) + dvd_SD_chi_dif_ve(:,:,tr_num)
        end if 
        if (Redi) then
            dvd_SD_tot(:,:,tr_num) = dvd_SD_tot(:,:,tr_num) + dvd_SD_chi_dif_heR(:,:,tr_num) &
                                                            + dvd_SD_chi_dif_veR(:,:,tr_num) &
                                                            + dvd_SD_chi_dif_viR(:,:,tr_num)
            
        end if 
        if (tracers%data(tr_num)%smooth_bh_tra) then
            dvd_SD_tot(:,:,tr_num) = dvd_SD_tot(:,:,tr_num) + dvd_SD_chi_dif_hbh(:,:,tr_num)
        end if 
        
    end do !--> do tr_num=1,2 
    
    !___________________________________________________________________________
    ! subtract the the bolus velocities back from 3D velocities:
    if (Fer_GM) then
!$OMP PARALLEL DO
        do elem=1, myDim_elem2D+eDim_elem2D
           UV(:, :, elem) = UV(:, :, elem) - fer_UV(:, :, elem)
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
        do node=1, myDim_nod2D+eDim_nod2D
           Wvel(:, node) = Wvel(:, node) - fer_Wvel(:, node)
        end do
!$OMP END PARALLEL DO
    end if
    
end subroutine compute_dvd
!
!
!_______________________________________________________________________________
! add time derivativ of 2nd. moment tracer 
! Xchi^(n+1) =  ...+ ( V^(n+1)*(Tr^(n+1))^2 + V^(n)*(Tr^(n))^2 )/ V^(n+1) +...
subroutine dvd_add_time_deriv(tr_num, dvd_tot, tr_old, tr, partit, mesh)
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    integer       , intent(in)             :: tr_num
    real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
    real(kind=WP) , intent(in)             :: tr_old( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer                                :: node, nz, nu1, nl1
    real(kind=WP)                          :: tr_n, tr_np1
    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"  
   
    ! add time derivative of 2nd moment tracer´
    do node=1, myDim_nod2D
        nu1 = ulevels_nod2D(node) 
        nl1 = nlevels_nod2D(node)
        do nz=nu1, nl1-1
            ! (T^(n+1))^2
            tr_np1 = hnode(nz,node) * tr(nz, node)*tr(nz, node)
            ! (T^(n))^2
            ! !!! ATTENTION: !!!
            ! hnode^(n+1) is stored in the variable hnode and hnode^(n) is 
            ! stored here in the variable hnode_new, since hnode^(n) was other
            ! wise overwritten by hnode^(n+1) in the variable hnode
            tr_n = hnode_new(nz,node) * tr_old(nz, node)*tr_old(nz, node)
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) - (tr_np1-tr_n)/dt/hnode(nz, node)
            !                                                     |
            !         +---------this is the minus sign------------+
            !         |
            !         v
            !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
            !   |
            !   v
            !  Kv*(nabla*Tr)^2
        end do !--> do nz=nu1, nl1-1
    end do !-->do node=1, myDim_nod2D
    
end subroutine dvd_add_time_deriv
!
!
!_______________________________________________________________________________
! add vertical 2nd. moment tracer flux
! Xchi^(n+1) =  ...+ ( Tr^(n+1) * Advflx_ver[Tr^(n+1)] )/ V^(n+1) +...
subroutine dvd_add_advflux_hor(do_SDdvd, tr_num, dvd_tot, trflx_h, UV, trstar, dump, partit, mesh)
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    integer       , intent(in)             :: tr_num
    logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
    real(kind=WP) , intent(inout)          :: dvd_tot( mesh%nl-1, partit%myDim_nod2D +partit%eDim_nod2D, 2)
    real(kind=WP) , intent(in)             :: trflx_h( mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP) , intent(in)             :: UV(   2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP) , intent(in)             :: trstar(  mesh%nl-1, partit%myDim_nod2D +partit%eDim_nod2D)
    real(kind=WP) , intent(inout)          :: dump(    mesh%nl-1, partit%myDim_nod2D +partit%eDim_nod2D)
    integer                                :: node, edge, nz, nu12, nl12, nl1, nl2, nu1, nu2 
    integer                                :: ednodes(2), edelem(2), elnodes(3)
    real(kind=WP)                          :: deltaX1, deltaY1, deltaX2, deltaY2, vflx, tr2flx, xchi
    real(kind=WP)                          :: eps_vflx=1.0e-16
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 
    !___________________________________________________________________________
    ! reuse here an already allocated working array --> initialise first 
    do node=1, myDim_nod2D+eDim_nod2D
        dump(:, node)=0.0_WP
    end do

    ! compute horizontal 2nd moment tracer flux from advective tracer fluxes at mid 
    ! edge face (tr^2) and add it to dvd_tot
    do edge = 1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        ednodes = edges(:,edge)
        ! local index of element that contribute to edge
        edelem  = edge_tri(:,edge)
        !nl1 ... num of layers -1 at elem el(1), nu1...idx off surf layer in case of cavity !=1
        nl1     = nlevels(edelem(1))-1
        nu1     = ulevels(edelem(1))
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1 = edge_cross_dxdy(1,edge)
        deltaY1 = edge_cross_dxdy(2,edge)
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2     = 0
        nu2     = 0
        if(edelem(2)>0) then
            deltaX2 = edge_cross_dxdy(3,edge)
            deltaY2 = edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2     = nlevels(edelem(2))-1
            nu2     = ulevels(edelem(2))
        end if !-->if(el(2)>0) then
        !_______________________________________________________________________
        ! nl12 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed, nu12 ... upper index of layers between element 
        ! el(1) & el(2) that contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)
        !_______________________________________________________________________
        ! (A) this loop only when edge has 1 facing element el(1) -> bnd edge
        do nz=nu1, nu12-1
            ! volume and 2nd moment tracer flux across across mid-edge faces
            ! --> (tr^2)*vflx/V = (tr*vflx)^2/vflx/(A*h)
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during advection, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new ! 
            elnodes= elem2D_nodes(:, edelem(1))
            vflx   = (-UV(2, nz, edelem(1))*deltaX1 + UV(1, nz, edelem(1))*deltaY1)*sum(hnode_new(nz, elnodes))/3.0_WP
            
            ! knut way to compute 2nd moment horizontal advection flux
            tr2flx = (trflx_h(nz, edge)*trflx_h(nz, edge))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx)
            
            ! sergey way to compute advective component of DVD eq. 26
            !  U = u*h
            !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
            !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]
            xchi   = 2.0_WP * -trflx_h(nz, edge)*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                 |-> this minus because trflx_v contains the 
            !                      negative tracer flx we need it positive
            xchi   = xchi - vflx*(trstar(nz,ednodes(1))+trstar(nz,ednodes(2)))*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                |-> 2.0_WP*0.5_WP
            
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), tr2flx, do_SDdvd)
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), tr2flx, do_SDdvd)
        end do !--> do nz=nu1, nu12-1
        
        !_______________________________________________________________________
        ! (B) this loop only when edge has 1 facing element el(2) -> bnd edge --> in cavity case
        if (nu2 > 0) then
            do nz=nu2, nu12-1
                ! volume and 2nd moment tracer flux across across mid-edge faces
                ! --> (tr^2)*vflx/V = (tr*vflx)^2/vflx/(A*h)
                ! !!! ATTENTION: !!!
                ! --> here helem must be the elemental thickness that was used 
                ! during advection, which was helem^(n), but helem^(n) got overwritten
                ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
                ! in the variable hnode_new (Although this sounds confusing!!!), so
                ! hnode_new contains here in reality hnode_old. So we recompute
                ! helem from hnode_new ! 
                elnodes= elem2D_nodes(:, edelem(2))
                vflx   = (UV(2, nz, edelem(2))*deltaX2 - UV(1, nz, edelem(2))*deltaY2)*sum(hnode_new(nz, elnodes))/3.0_WP
                
                ! knut way to compute 2nd moment horizontal advection flux
                tr2flx = (trflx_h(nz, edge)*trflx_h(nz, edge))
                tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx)
                
                ! sergey way to compute advective component of DVD eq. 26
                !  U = u*h
                !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
                !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]
                xchi   = 2.0_WP * -trflx_h(nz, edge)*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
                !                 |-> this minus because trflx_v contains the 
                !                      negative tracer flx we need it positive                         
                xchi   = xchi - vflx*(trstar(nz,ednodes(1))+trstar(nz,ednodes(2)))*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
                !                |-> 2.0_WP*0.5_WP
                
                dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), tr2flx, do_SDdvd)
                dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), tr2flx, do_SDdvd)
            end do !--> do nz=nu2, nu12-1
        end if !--> if (nu2 > 0) then
        
        !_______________________________________________________________________
        ! (C) this loop if edge has both triangles valid 
        do nz=nu12, nl12
            ! volume and 2nd moment tracer flux across across mid-edge faces
            ! --> (tr^2)*vflx/V = (tr*vflx)^2/vflx/(A*h)
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during advection, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new !
            elnodes= elem2D_nodes(:, edelem(1))
            vflx   =  ( -UV(2, nz, edelem(1))*deltaX1 + UV(1, nz, edelem(1))*deltaY1 )*sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(2))
            vflx   = vflx + (  UV(2, nz, edelem(2))*deltaX2 - UV(1, nz, edelem(2))*deltaY2 )*sum(hnode_new(nz, elnodes))/3.0_WP
            
            ! knut way to compute 2nd moment horizontal advection flux
            tr2flx = (trflx_h(nz, edge)*trflx_h(nz, edge))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx)
            
            ! sergey way to compute advective component of DVD eq. 26
            !  U = u*h
            !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
            !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]
            xchi   = 2.0_WP * -trflx_h(nz, edge)*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                 |-> this minus because trflx_v contains the 
            !                      negative tracer flx we need it positive
            xchi   = xchi - vflx*(trstar(nz,ednodes(1))+trstar(nz,ednodes(2)))*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                |-> 2.0_WP*0.5_WP
            
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), tr2flx, do_SDdvd)
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), tr2flx, do_SDdvd) 
        end do !--> do nz=nu12, nl12
        
        !_______________________________________________________________________
        ! (D) this loop only when edge has 1 facing element el(1) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl1
            ! volume and 2nd moment tracer flux across across mid-edge faces
            ! --> (tr^2)*vflx/V = (tr*vflx)^2/vflx/(A*h)
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during advection, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new ! 
            elnodes= elem2D_nodes(:, edelem(1))
            vflx   = (-UV(2, nz, edelem(1))*deltaX1 + UV(1, nz, edelem(1))*deltaY1)*sum(hnode_new(nz, elnodes))/3.0_WP
            
            ! knut way to compute 2nd moment horizontal advection flux
            tr2flx = (trflx_h(nz, edge)*trflx_h(nz, edge))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx)
            
            ! sergey way to compute advective component of DVD eq. 26
            !  U = u*h
            !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
            !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]
            xchi   = 2.0_WP * -trflx_h(nz, edge)*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                 |-> this minus because trflx_v contains the 
            !                      negative tracer flx we need it positive
            xchi   = xchi - vflx*(trstar(nz,ednodes(1))+trstar(nz,ednodes(2)))*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                |-> 2.0_WP*0.5_WP
            
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), tr2flx, do_SDdvd)
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), tr2flx, do_SDdvd)
        end do !--> do nz=nl12+1, nl1
        
        !_______________________________________________________________________
        ! (E) this loop only when edge has 1 facing element el(2) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl2
            ! volume and 2nd moment tracer flux across across mid-edge faces
            ! --> (tr^2)*vflx/V = (tr*vflx)^2/vflx/(A*h)
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during advection, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new ! 
            elnodes= elem2D_nodes(:, edelem(2))
            vflx   = (UV(2, nz, edelem(2))*deltaX2 - UV(1, nz, edelem(2))*deltaY2)*sum(hnode_new(nz, elnodes))/3.0_WP
            
            ! knut way to compute 2nd moment horizontal advection flux
            tr2flx = (trflx_h(nz, edge)*trflx_h(nz, edge))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx)
            
            ! sergey way to compute advective component of DVD eq. 26 (approximation small dt)
            !  U = u*h
            !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
            !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]
            xchi   = 2.0_WP * -trflx_h(nz, edge)*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                 |-> this minus because trflx_v contains the 
            !                      negative tracer flx we need it positive
            xchi   = xchi - vflx*(trstar(nz,ednodes(1))+trstar(nz,ednodes(2)))*(trstar(nz,ednodes(1))-trstar(nz,ednodes(2)))
            !                |-> 2.0_WP*0.5_WP         
            
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), tr2flx, do_SDdvd)
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), tr2flx, do_SDdvd)
        end do !--> do nz=nl12+1, nl2
    end do !--> do edge=1, myDim_edge2D
    call exchange_nod(dump, partit)
    
    !___________________________________________________________________________
    ! switch from Xchi_tilde --> Xchi
    if (do_SDdvd) then 
        ! In Sergeys case normalization should be done with volume at moment of 
        ! diffusion V^(n) 
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) + dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode_new(nu1:nl1, node))
        end do
    else
        ! In Knuts case volume normalization is done with V^(n+1)
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) - dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode(nu1:nl1, node))
            !                                                               |
            !                     +--------this is the minus sign-----------+
            !                     |
            !                     v
            !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
            !   |
            !   v
            !  Kv*(nabla*Tr)^2
        end do
    end if     
end subroutine dvd_add_advflux_hor
!
!
!_______________________________________________________________________________
! add vertical 2nd. moment tracer flux
! Xchi^(n+1) =  ...+ ( Tr^(n+1) * Advflx_ver[Tr^(n+1)] )/ V^(n+1) +...
subroutine dvd_add_advflux_ver(do_SDdvd, tr_num, dvd_tot, trflx_v, Wvel, trstar, partit, mesh)
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    integer       , intent(in)             :: tr_num
    logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
    real(kind=WP) , intent(inout)          :: dvd_tot( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
    real(kind=WP) , intent(in)             :: trflx_v( mesh%nl  , partit%myDim_nod2D)
    real(kind=WP) , intent(in)             :: Wvel(    mesh%nl  , partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: trstar(  mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer                                :: node, nz, nu1, nl1
    real(kind=WP)                          :: tr2flx, vflx, xchi(mesh%nl), trstar_zlev
    real(kind=WP)                          :: eps_vflx=1.0e-16
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"  
    if (do_SDdvd) then 
        !_______________________________________________________________________
        !  compute DVD contribution from vertical advection after Sergeys method
        !   --> see eq. 26 in T. Banerjee et al. 2023 --> U = u*h, T^tilde face reconstructed
        !   value
        !  Xchi_(i+0.5) = 2*[U_(i+0.5)* T^tilde_(i+0.5) * ( Tstar_i-Tstar_(i-1) ) -
        !                    U_(i+0.5)*0.5*(Tstar_i+Tstar_(i-1))*(Tstar_i-Tstar_(i-1))]     
        do node=1,myDim_nod2D
            nu1 = ulevels_nod2D(node) 
            nl1 = nlevels_nod2D(node)
            !___________________________________________________________________
            xchi = 0.0_WP
            
            !___________________________________________________________________
            ! surface xchi
            nz=nu1
            trstar_zlev = trstar(nz, node)
            xchi(nz) = 2.0_WP*-trflx_v(nz, node)*(trstar(nz, node)-trstar(nz+1, node))
            xchi(nz) = xchi(nz) - 2.0_WP*Wvel(nz, node)*area(nz, node)*trstar_zlev*(trstar(nz, node)-trstar(nz+1, node))
            
            ! surface and bulk dvd, bulk xchi 
            do nz=nu1+1,nl1-1
                
                ! compute trstar on zlev interface
                !PStrstar_zlev = 0.5_WP*(trstar(nz-1, node)+trstar(nz, node))
                trstar_zlev = (trstar(nz-1, node)*hnode(nz-1, node)+trstar(nz, node)*hnode(nz, node))/(hnode(nz-1,node)+hnode(nz,node))
                
                ! --> here we are on full depth levels eq. 26 for small dt
                xchi(nz) = 2.0_WP*-trflx_v(nz, node)*(trstar(nz-1, node)-trstar(nz, node))
                !                 |-> this minus because trflx_v contains the 
                !                     negative tracer flx we need it positive
                xchi(nz) = xchi(nz) - 2.0_WP*Wvel(nz, node)*area(nz, node)*trstar_zlev*(trstar(nz-1, node)-trstar(nz, node))
                
                ! !!! ATTENTION: !!!
                ! --> here hnode must be the vertice thickness that was used 
                ! during advection, which was hnode^(n), but hnode^(n) got overwritten
                ! with hnode^(n+1) in call update_thickness(...). But we saved hnode^(n)
                ! in the variable hnode_new (Although this sounds confusing!!!)                
                ! --> here we are on mid-depth levels
                dvd_tot(nz-1, node, tr_num) = dvd_tot(nz-1, node, tr_num) - (dvd_xdfac*xchi(nz-1)+(1.0_WP-dvd_xdfac)*xchi(nz))/(areasvol(nz-1, node)*hnode_new(nz-1, node))
                
            end do !--> do nz=nu1+1,nl1-1
            !___________________________________________________________________
            ! same for bottom 
            nz=nl1-1
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) - (dvd_xdfac*xchi(nz)+(1.0_WP-dvd_xdfac)*xchi(nz+1))/(areasvol(nz, node)*hnode_new(nz, node))
                
        end do !--> do node=1,myDim_nod2D
    
    else
        !_______________________________________________________________________
        ! compute vertical 2nd moment tracer flux from advective tracer fluxes 
        ! at scalar cell prism interface after Klingbeil etal 2014
        do node=1,myDim_nod2D
            nu1 = ulevels_nod2D(node) 
            nl1 = nlevels_nod2D(node)
            !___________________________________________________________________
            ! surface
            nz     = nu1
            ! in case of cavity nu1 is also a hard boundary with boundary condition 
            ! Wvel(nu1)==0 --> 1/(Wvel(nz,n)*area(nz,n))--> Inf
            ! !!! ATTENTION: !!!
            ! --> here hnode must be the vertice thickness that was used 
            ! during advection, which was hnode^(n), but hnode^(n) got overwritten
            ! with hnode^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!)
            ! volume and tracer flux across upper cell prism face
            ! --> Keep in mind that the advective tracer fluxes are saved as negative fluxes
            !     thats why 
            vflx   = Wvel(nz, node)*area(nz, node)
            tr2flx = (trflx_v(nz, node)*trflx_v(nz, node))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx) ! prevent division by zero vflx
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) - tr2flx/(areasvol(nz, node)*hnode(nz, node))
            
            ! volume and tracer flux across lower cell prism face
            vflx   = Wvel(nz+1, node)*area(nz+1, node)
            tr2flx = (trflx_v(nz+1, node)*trflx_v(nz+1, node))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx) ! prevent division by zero vflx
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + tr2flx/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk 
            do nz=nu1+1,nl1-2 ! stop two layers over botrtom since Wvel(nl1)==0 --> 1/Wvel(nl1) --> Inf
                ! volume and tracer flux across upper cell prism face
                vflx   = Wvel(nz, node)*area(nz, node)
                tr2flx = (trflx_v(nz, node)*trflx_v(nz, node))
                tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx) ! prevent division by zero vflx
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) - tr2flx/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |
                ! volume and tracer flux across lower cell prism face |
                vflx   = Wvel(nz+1, node)*area(nz+1, node)!           |
                tr2flx = (trflx_v(nz+1, node)*trflx_v(nz+1, node))!  |
                tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx) ! prevent division by zero vflx
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + tr2flx/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |
                !                     +----this is the minus sign-----+
                !                     |
                !                     v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
                !
                !
            end do !--> do nz=nu1,nl1
            
            !___________________________________________________________________
            ! bottom
            nz     = nl1-1
            vflx   = Wvel(nz,node)*area(nz,node)
            tr2flx = (trflx_v(nz, node)*trflx_v(nz, node))
            tr2flx = merge(tr2flx/vflx, 0.0_WP, abs(vflx)>eps_vflx) ! prevent division by zero vflx
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) - tr2flx/(areasvol(nz, node)*hnode(nz, node))
            ! --> do here not nz+1 since  Wvel(nz+1)==0 (bottom bnd condition) --> 1/(Wvel(nz+1,n)*area(nz+1,n)) --> Inf
        end do !--> do node=1,myDim_nod2D
    end if 
end subroutine dvd_add_advflux_ver
!
!
!_______________________________________________________________________________
! add horizontal diffusive flux after Klingbeil et al. 2014, take Redi into account
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx_hor[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx_hor[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_horexpl(do_SDdvd, tr_num, dvd_tot, tr, Ki, tr_xy, dump, partit, mesh)
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    integer       , intent(in)             :: tr_num
    logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
    real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
    real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: Ki(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
   real(kind=WP) , intent(in)              :: tr_xy(2,mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D)
    real(kind=WP) , intent(inout)          :: dump(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    
    integer                                :: node, edge, nz, nu12, nl12, nl1, nl2, nu1, nu2 
    integer                                :: ednodes(2), edelem(2), elnodes(3)
    real(kind=WP)                          :: deltaX1, deltaY1, deltaX2, deltaY2, dZ
    real(kind=WP)                          :: Kh, Trx, Try, Trz(2), SxTz, SyTz, Dflx, Trc, xchi
    real(kind=WP)                          :: isredi=0.0_WP
    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

    !___________________________________________________________________________
    ! reuse here an already allocated working array --> initialise first 
    do node=1, myDim_nod2D+eDim_nod2D
        dump(:, node)=0.0_WP
    end do
    
    ! add contribution of redi diffusivity when activated
    if (Redi) isredi=1.0_WP
    
    ! compute horizontal 2nd moment tracer flux from advective tracer fluxes at mid 
    ! edge face (tr^2) and add it to dvd_tot
    do edge = 1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        ednodes  = edges(:,edge)
        ! local index of element that contribute to edge
        edelem  = edge_tri(:,edge)
        !nl1 ... num of layers -1 at elem el(1), nu1...idx off surf layer in case of cavity !=1
        nl1     = nlevels(edelem(1))-1
        nu1     = ulevels(edelem(1))
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1 = edge_cross_dxdy(1,edge)
        deltaY1 = edge_cross_dxdy(2,edge)
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2     = 0
        nu2     = 0
        if(edelem(2)>0) then
            deltaX2 = edge_cross_dxdy(3,edge)
            deltaY2 = edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2     = nlevels(edelem(2))-1
            nu2     = ulevels(edelem(2))
        end if !-->if(el(2)>0) then
        !_______________________________________________________________________
        ! nl12 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed, nu12 ... upper index of layers between element 
        ! el(1) & el(2) that contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)
        
        !_______________________________________________________________________
        ! (A) this loop only when edge has 1 facing element el(1) -> bnd edge
        do nz=nu1, nu12-1
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trx  = tr_xy(1, nz, edelem(1)) 
            Try  = tr_xy(2, nz, edelem(1))
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(-deltaX1*Try + deltaY1*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value --> knut variant
            ! --> here variable tr corresponds to tr_old = Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = -2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar dfiference --> sergey variant
            ! --> here variable tr corresponds to trstar = (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !--> do nz=nu1, nu12-1
        
        !_______________________________________________________________________
        ! (B) this loop only when edge has 1 facing element el(2) -> bnd edge --> in cavity case
        if (nu2 > 0) then
            do nz=nu2, nu12-1
                !_______________________________________________________________
                ! !!! ATTENTION: !!!
                ! --> here helem must be the elemental thickness that was used 
                ! during diffusion, which was helem^(n), but helem^(n) got overwritten
                ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
                ! in the variable hnode_new (Although this sounds confusing!!!), so
                ! hnode_new contains here in reality hnode_old. So we recompute
                ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
                elnodes= elem2D_nodes(:, edelem(2))
                dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
                ! also Ki is defined at mid depth levels
                Kh   = sum(Ki(nz, ednodes))*0.5_WP
                ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
                ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
                ! tracer gradient back to mid depth levels
                Trx  = tr_xy(1, nz, edelem(2))
                Try  = tr_xy(2, nz, edelem(2))
                ! Dflx is here the horizontal diffusive volume flux 
                Dflx = Kh*(deltaX2*Try - deltaY2*Trx)*dZ
                !_______________________________________________________________
                ! compute edge centered tracer value
                ! --> here variable tr corresponds to tr_old == Tr^(n)
                Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
                Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
                !_______________________________________________________________
                ! compute Tstar diference --> sergey variant
                ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
                xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
                !_______________________________________________________________
                dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
                dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
            end do !--> do nz=nu2, nu12-1
        end if !--> if (nu2 > 0) then
        
        !_______________________________________________________________________
        ! (C) this loop if edge has both triangles valid 
        do nz=nu12, nl12
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(2))
            dZ   = (dZ + sum(hnode_new(nz, elnodes))/3.0_WP)*0.5_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trx  = 0.5_WP*( tr_xy(1, nz, edelem(1))+tr_xy(1, nz, edelem(2)) )
            Try  = 0.5_WP*( tr_xy(2, nz, edelem(1))+tr_xy(2, nz, edelem(2)) )
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*((deltaX2-deltaX1)*Try-(deltaY2-deltaY1)*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !--> do nz=nu12, nl12
        
        !_______________________________________________________________________
        ! (D) this loop only when edge has 1 facing element el(1) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl1
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trx  = tr_xy(1, nz, edelem(1))
            Try  = tr_xy(2, nz, edelem(1))
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(-deltaX1*Try + deltaY1*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !-->do nz=nl12+1, nl1
        
        !_______________________________________________________________________
        ! (E) this loop only when edge has 1 facing element el(2) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl2
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(2))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trx  = tr_xy(1, nz, edelem(2))
            Try  = tr_xy(2, nz, edelem(2))
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(deltaX2*Try - deltaY2*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1))*hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2))*hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
            
        end do !--> do nz=nl12+1, nl2
    end do !--> do edge=1, myDim_edge2D
    call exchange_nod(dump, partit)
    
    !___________________________________________________________________________
    ! switch to volume normalized Xchi
    if (do_SDdvd) then 
        ! In Sergeys case normalization should be done with volume at moment of 
        ! diffusion V^(n) 
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) + dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode_new(nu1:nl1, node))
        end do
    else
        ! In Knuts case volume normalization is done with V^(n+1)
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) + dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode(nu1:nl1, node))
            !                                                               |
            !                                          this is the          |
            !                                   +--------plus sign----------+
            !                                   |
            !                                   v
            !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
            !   |
            !   v
            !  Kv*(nabla*Tr)^2
        end do
    end if 

end subroutine dvd_add_difflux_horexpl
!
!
!_______________________________________________________________________________
! add horizontal diffusive flux after Klingbeil et al. 2014, take Redi into account
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx_hor[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx_hor[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_horexplredi(do_SDdvd, tr_num, dvd_tot, tr, Ki, slope, tr_z, dump, partit, mesh)
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh)  , intent(in)   , target  :: mesh
    integer       , intent(in)             :: tr_num
    logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
    real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
    real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: Ki(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: slope(3,mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(in)             :: tr_z(   mesh%nl  , partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP) , intent(inout)          :: dump(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    
    integer                                :: node, edge, nz, nu12, nl12, nl1, nl2, nu1, nu2 
    integer                                :: ednodes(2), edelem(2), elnodes(3)
    real(kind=WP)                          :: deltaX1, deltaY1, deltaX2, deltaY2, dZ
    real(kind=WP)                          :: Kh, Trx, Try, Trz(2), SxTz, SyTz, Dflx, Trc, xchi
    real(kind=WP)                          :: isredi=0.0_WP
    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h" 

    !___________________________________________________________________________
    ! reuse here an already allocated working array --> initialise first 
    do node=1, myDim_nod2D+eDim_nod2D
        dump(:, node)=0.0_WP
    end do
    
    ! add contribution of redi diffusivity when activated
    if (Redi) isredi=1.0_WP
    
    ! compute horizontal 2nd moment tracer flux from advective tracer fluxes at mid 
    ! edge face (tr^2) and add it to dvd_tot
    do edge = 1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        ednodes  = edges(:,edge)
        ! local index of element that contribute to edge
        edelem  = edge_tri(:,edge)
        !nl1 ... num of layers -1 at elem el(1), nu1...idx off surf layer in case of cavity !=1
        nl1     = nlevels(edelem(1))-1
        nu1     = ulevels(edelem(1))
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1 = edge_cross_dxdy(1,edge)
        deltaY1 = edge_cross_dxdy(2,edge)
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2     = 0
        nu2     = 0
        if(edelem(2)>0) then
            deltaX2 = edge_cross_dxdy(3,edge)
            deltaY2 = edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2     = nlevels(edelem(2))-1
            nu2     = ulevels(edelem(2))
        end if !-->if(el(2)>0) then
        !_______________________________________________________________________
        ! nl12 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed, nu12 ... upper index of layers between element 
        ! el(1) & el(2) that contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)
        
        !_______________________________________________________________________
        ! (A) this loop only when edge has 1 facing element el(1) -> bnd edge
        do nz=nu1, nu12-1
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trz  = 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes))
            SxTz = sum(Trz*slope(1, nz, ednodes))*0.5_WP
            SyTz = sum(Trz*slope(2, nz, ednodes))*0.5_WP
            Trx  = SxTz*isredi
            Try  = SyTz*isredi
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(-deltaX1*Try + deltaY1*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value --> knut variant
            ! --> here variable tr corresponds to tr_old = Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = -2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar dfiference --> sergey variant
            ! --> here variable tr corresponds to trstar = (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !--> do nz=nu1, nu12-1
        
        !_______________________________________________________________________
        ! (B) this loop only when edge has 1 facing element el(2) -> bnd edge --> in cavity case
        if (nu2 > 0) then
            do nz=nu2, nu12-1
                !_______________________________________________________________
                ! !!! ATTENTION: !!!
                ! --> here helem must be the elemental thickness that was used 
                ! during diffusion, which was helem^(n), but helem^(n) got overwritten
                ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
                ! in the variable hnode_new (Although this sounds confusing!!!), so
                ! hnode_new contains here in reality hnode_old. So we recompute
                ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
                elnodes= elem2D_nodes(:, edelem(2))
                dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
                ! also Ki is defined at mid depth levels
                Kh   = sum(Ki(nz, ednodes))*0.5_WP
                ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
                ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
                ! tracer gradient back to mid depth levels
                Trz  = 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes))
                SxTz = sum(Trz*slope(1, nz, ednodes))*0.5_WP
                SyTz = sum(Trz*slope(2, nz, ednodes))*0.5_WP
                Trx  = SxTz*isredi
                Try  = SyTz*isredi
                ! Dflx is here the horizontal diffusive volume flux 
                Dflx = Kh*(deltaX2*Try - deltaY2*Trx)*dZ
                !_______________________________________________________________
                ! compute edge centered tracer value
                ! --> here variable tr corresponds to tr_old == Tr^(n)
                Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
                Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
                !_______________________________________________________________
                ! compute Tstar diference --> sergey variant
                ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
                xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
                !_______________________________________________________________
                dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
                dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
            end do !--> do nz=nu2, nu12-1
        end if !--> if (nu2 > 0) then
        
        !_______________________________________________________________________
        ! (C) this loop if edge has both triangles valid 
        do nz=nu12, nl12
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(2))
            dZ   = (dZ + sum(hnode_new(nz, elnodes))/3.0_WP)*0.5_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trz  = 0.5_WP*(tr_z(nz, ednodes)+tr_z(nz+1, ednodes))
            SxTz = sum(Trz*slope(1, nz, ednodes))*0.5_WP
            SyTz = sum(Trz*slope(2, nz, ednodes))*0.5_WP
            Trx  = SxTz*isredi
            Try  = SyTz*isredi
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*((deltaX2-deltaX1)*Try-(deltaY2-deltaY1)*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !--> do nz=nu12, nl12
        
        !_______________________________________________________________________
        ! (D) this loop only when edge has 1 facing element el(1) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl1
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(1))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trz  = 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes))
            SxTz = sum(Trz*slope(1, nz, ednodes))*0.5_WP
            SyTz = sum(Trz*slope(2, nz, ednodes))*0.5_WP
            Trx  = SxTz*isredi
            Try  = SyTz*isredi
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(-deltaX1*Try + deltaY1*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1)) * hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
        end do !-->do nz=nl12+1, nl1
        
        !_______________________________________________________________________
        ! (E) this loop only when edge has 1 facing element el(2) -> bnd edge 
        ! at bottom topography
        do nz=nl12+1, nl2
            !___________________________________________________________________
            ! !!! ATTENTION: !!!
            ! --> here helem must be the elemental thickness that was used 
            ! during diffusion, which was helem^(n), but helem^(n) got overwritten
            ! with helem^(n+1) in call update_thickness(...). But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            elnodes= elem2D_nodes(:, edelem(2))
            dZ   = sum(hnode_new(nz, elnodes))/3.0_WP
            ! also Ki is defined at mid depth levels
            Kh   = sum(Ki(nz, ednodes))*0.5_WP
            ! reminder: tr_z ... vertical tracer gradient sits on full depth levels
            ! with 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes)) we bring the vertical 
            ! tracer gradient back to mid depth levels
            Trz  = 0.5_WP*(tr_z(nz,ednodes)+tr_z(nz+1,ednodes))
            SxTz = sum(Trz*slope(1, nz, ednodes))*0.5_WP
            SyTz = sum(Trz*slope(2, nz, ednodes))*0.5_WP
            Trx  = SxTz*isredi
            Try  = SyTz*isredi
            ! Dflx is here the horizontal diffusive volume flux 
            Dflx = Kh*(deltaX2*Try - deltaY2*Trx)*dZ
            !___________________________________________________________________
            ! compute edge centered tracer value
            ! --> here variable tr corresponds to tr_old == Tr^(n)
            Trc  = tr(nz, ednodes(1))*hnode_new(nz, ednodes(1)) + tr(nz, ednodes(2))*hnode_new(nz, ednodes(2))
            Trc  = 2.0_WP * Trc/(hnode_new(nz, ednodes(1))+hnode_new(nz, ednodes(2)))
            !___________________________________________________________________
            ! compute Tstar diference --> sergey variant
            ! --> here variable tr corresponds to trstar == (Tr^(n+1) + Tr^(n))*0.5
            xchi = -2.0_WP * (tr(nz, ednodes(1))-tr(nz, ednodes(2)))
            !___________________________________________________________________
            dump(nz, ednodes(1)) = dump(nz, ednodes(1)) + merge( xchi*(       dvd_xdfac), Trc, do_SDdvd)*Dflx
            dump(nz, ednodes(2)) = dump(nz, ednodes(2)) - merge(-xchi*(1.0_WP-dvd_xdfac), Trc, do_SDdvd)*Dflx
            
        end do !--> do nz=nl12+1, nl2
    end do !--> do edge=1, myDim_edge2D
    call exchange_nod(dump, partit)
    
    !___________________________________________________________________________
    ! switch to volume normalized Xchi
    if (do_SDdvd) then 
        ! In Sergeys case normalization should be done with volume at moment of 
        ! diffusion V^(n) 
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) + dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode_new(nu1:nl1, node))
        end do
    else
        ! In Knuts case volume normalization is done with V^(n+1)
        do node=1, myDim_nod2D+eDim_nod2D
            nu1 = ulevels_nod2D(node)
            nl1 = nlevels_nod2D(node)-1
            dvd_tot(nu1:nl1, node, tr_num) = dvd_tot(nu1:nl1, node, tr_num) + dump(nu1:nl1, node)/(areasvol(nu1:nl1, node)*hnode(nu1:nl1, node))
            !                                                               |
            !                                          this is the          |
            !                                   +--------plus sign----------+
            !                                   |
            !                                   v
            !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
            !   |
            !   v
            !  Kv*(nabla*Tr)^2
        end do
    end if 

end subroutine dvd_add_difflux_horexplredi
!
!
!_______________________________________________________________________________
! add explicite vertical REdi diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! only in case of Redi!!!
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_vertexpl(do_SDdvd, tr_num, dvd_tot, tr, trstar, Kv, partit, mesh)
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        integer       , intent(in)             :: tr_num
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: trstar( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: Kv(     mesh%nl  , partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, elem, nz, nu1, nl1, k
        real(kind=WP)                          :: tr_up, tr_dwn, Dflx(mesh%nl), dzup, dzdwn, trnx, trny, sbc
        real(kind=WP)                          :: tr_xynodes(2,mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP)                          :: trN(mesh%nl-1)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    ! add contribution of implicite vertical diffusion (after Klingbeil et al.2014)
    do node=1, myDim_nod2D
        nu1  = ulevels_nod2D(node) 
        nl1  = nlevels_nod2D(node)
        Dflx = 0.0_WP
        
        !_______________________________________________________________________
        if (do_SDdvd) then
            ! here: trstar = (T^(n+1) + T^(n))*0.5, tr = T^(n+1)
            trN = 2.0_WP*trstar(:, node) - tr(:, node)
        else
            ! here: trstar is T^(n)
            trN = trstar(:, node)
        end if
        
        !_______________________________________________________________________
        ! compute diffusive flux at full surface levels: ---> Dflx[Tr^(n+1)]
        ! Dflx = Kv * dT/dz
        !   ---------------- zbar_1  -+
        !                             |
        !          o tr_1 ------------|- hnode_1 -+
        !                             |           |
        !   --Dflx-x-Tr----- zbar_2 --+           |- (hnode_1+hnod_2)/2
        !                             |           |
        !          o tr_2 ------------|- hnode_2 -+
        !                             |
        !   ---------------- zbar_3  -+
        !_______________________________________________________________________
        ! !!! ATTENTION: !!!
        ! Implicite vertical Diffusion is done with hnode^(n+1) which is stored 
        ! in the variable hnode (keep in mind at this point hnode_new was used
        ! to rescue the varaible hnode^(n) ).
        do nz=nu1+1, nl1-1
            ! But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            dzup     = hnode_new(nz-1, node)*0.5_WP
            dzdwn    = hnode_new(nz  , node)*0.5_WP
            
            !PS Dflx(nz) = Kv(nz, node) * (tr(nz-1, node)-tr(nz, node)) /(dzup+dzdwn) * area(nz, node)
            Dflx(nz) = Kv(nz, node) * (trN(nz-1)-trN(nz)) /(dzup+dzdwn) * area(nz, node)
            
            ! extent Dflx after Sergey to compute the diffusive component of DVD 
            ! see T. Banerjee et al. 2023 eq. 7
            ! Xchi_(i+0.5) = 2*Kv * (T^(n+1)_(i+1)-T^(n+!)_i)/dz * (Tstar_(i+1)-Tstar_(i)/dz
            ! Dflx(nz) = Dflx * 2*(Tstar_(i+1)-Tstar_(i)/dz
            ! --> if Sergey diagonsitc is used Dflx correspond from here on to 
            !     Xchi at full depth level interface
            ! --> at this points tr must be trstar
            Dflx(nz) = Dflx(nz) * merge(-2.0_WP*( trstar(nz, node)-trstar(nz-1, node) ), 1.0_WP, do_SDdvd)
            !                       |
            !                       +-> merge: if do_SDdvd==True use first argument, if
            !                           do_SDdvd==False use second argument of merge!
        end do !--> do nz=nu1+1, nl1-1
        
        !_______________________________________________________________________
        ! decide between Sergeys and Knuts DVD method
        if (do_SDdvd) then 
            !___________________________________________________________________
            ! In Sergeys DVD diagnostic we do not take into account the surface 
            ! and interior fluxes of surface boundary condition and short wave 
            ! penetration, since akthough they are source terms for variance, 
            ! they not contribute to the decay variance !!!
            ! If they arte considered than only in a separate term but not added
            ! to the variance decay term of diffusivity. However in Knuts DVD |
            ! diagnostic this is different. Here the source term have to be explicitly 
            ! taken into account since they already contributed to (T^(n+1))^2-(T^(n))^2
            
            ! now compute flx of Xchi into control volume, see T. Banerjee et al. 2023
            ! eq. 9. Each control volume gets one half of the Xchi interface value 
            ! --> keep in mind we are now again on mid depth levels and in this case
            !     Dflx contains the interface value of Xchi
            ! --> compute volume normlized Xchi contribution
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + (dvd_xdfac*Dflx(nz)+(1.0_WP-dvd_xdfac)*Dflx(nz+1))/(areasvol(nz, node)*hnode_new(nz, node))
            end do !--> do nz=nu1+1, nl1-1
        else
            !___________________________________________________________________
            ! apply surface boundary condition --> have to do this here from hand -->
            ! i was somehow not able to bind in the function module for the surface 
            ! boundary condition
            ! d/dt(h*T) + div_h*(u*h*T) + h*d/dz(w*T) - div_3(h*Kv*div_3*T) = -(Q_hf + W*T_w)_k=1
            !     |
            ! need eq. for 2nd. moment    
            ! tracer thickness variance
            !     |
            !     V
            ! d/dt(h*T^2) = T*d/dt(h*T) + h*T*d/dt(T)+
            !                                      |
            !     +-------eq. for d/dt(T)----------+
            !     |
            ! d/dt(T) + div_h*(u*T) + d/dz(w*T) - div_3(Kv*div_3*T) = -(Q_hf + W*(T_w-T_k))_k=1 / h_k
            !
            ! d/dt(h*T^2) + div*(u*h*T^2) + h*d/dz(w*T^2) - ... = (-2*T_k*Q_hf - W*T_w*T_k - W*T_k*(T_w-T_k) )_k=1
            !                                                   = 2*(-T_k*Q_hf - 0.5*W*T_w*T_k - 0.5*W*T_k*(T_w-T_k) )_k=1
            ! --> temp: T_w=T_k --> 2*(-T_k*Q_hf - 0.5*W*T_w*T_k)
            ! --> salt: T_w=0   --> 2*(-T_k*Q_hf + 0.5*W*T_w*T_k)
            nz = nu1
            if     (tr_num==1) then 
                sbc = -(heat_flux(node)/vcpw + trN(nz)*water_flux(node)*is_nonlinfs*0.5_WP)
            elseif (tr_num==2) then 
                sbc =  (virtual_salt(node) + relax_salt(node) - real_salt_flux(node)*is_nonlinfs*0.5_WP)
            end if 
            Dflx(nz) = Dflx(nz) + sbc*area(nz, node)
            
            !___________________________________________________________________
            ! apply shortwave penetration fluxes
            if (use_sw_pene .and. tr_num==1) then
                nz=nu1
                Dflx(nz) = Dflx(nz) + sw_3d(nz, node)*area(nz,node)
                do nz=nu1+1, nl1-1
                    Dflx(nz) = Dflx(nz) + sw_3d(nz, node)*area(nz,node)
                end do
            end if
            
            !___________________________________________________________________
            ! compute diffusive flux flxdiff(nz), after Klingbeil etal 2014
            ! flxdiff(nz)= ([Kv*dT/dz]_1*T_1*A_1-[Kv*dT/dz]_1*T_1*A_2)/V_1 into the volume 
            ! surface
            nz     = nu1
            tr_up  = Dflx(nz  )*trN(nz) 
            tr_dwn = (trN(nz)*hnode_new(nz, node)+trN(nz+1)*hnode_new(nz+1, node))/(hnode_new(nz, node)+hnode_new(nz+1, node))
            tr_dwn = Dflx(nz+1)*tr_dwn
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk
            do nz=nu1+1, nl1-2
                ! volume (hnode) weighted temperature reconstruction at upper/lower
                ! scalar cell interface
                tr_up  = (trN(nz-1)*hnode_new(nz-1, node)+trN(nz  )*hnode_new(nz  , node))/(hnode_new(nz-1, node)+hnode_new(nz  , node))
                tr_dwn = (trN(nz  )*hnode_new(nz  , node)+trN(nz+1)*hnode_new(nz+1, node))/(hnode_new(nz  , node)+hnode_new(nz+1, node))
                ! compute  T^(n+1)_(i+0.5) * Dflx_(i+0.5) through upper and lower face
                tr_up  =  Dflx(nz  )*tr_up 
                tr_dwn =  Dflx(nz+1)*tr_dwn
                ! copute dvd contribution normalized with volume 
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |    |-> factor 2 comes here from Klingbeil et al.2014
                !                                      this is the    |         (2*Tr^(n+1) * Dflx[Tr^(n+1)] 
                !                                   +----plus sign----+
                !                                   |
                !                                   v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
            end do !--> do nz=nu1+1, nl1-1
            
            !___________________________________________________________________
            ! bottom
            nz    = nl1-1
            tr_up = (trN(nz-1)*hnode_new(nz-1, node)+trN(nz)*hnode_new(nz  , node))/(hnode_new(nz-1, node)+hnode_new(nz  , node))
            tr_up = Dflx(nz)*tr_up
            tr_dwn= 0.0_WP
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
        end if 
    end do !--> do node=1, myDim_nod2D
end subroutine dvd_add_difflux_vertexpl
!
!
!_______________________________________________________________________________
! add explicite vertical REdi diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! only in case of Redi!!!
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_vertexplredi(do_SDdvd, tr_num, dvd_tot, tr, Ki, slope, tr_xy, partit, mesh)
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        integer       , intent(in)             :: tr_num
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: Ki(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: slope(3,mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: tr_xy(2,mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D+partit%eXDim_elem2D)
        integer                                :: node, elem, nz, nu1, nl1, k
        real(kind=WP)                          :: tr_up, tr_dwn, Dflx(mesh%nl), dzup, dzdwn, trnx, trny
        real(kind=WP)                          :: tr_xynodes(2,mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    ! compute horizontal tracer gradient on nodes
    do node=1, myDim_nod2D
        nu1=ulevels_nod2D(node)
        nl1=nlevels_nod2D(node)
        do nz=nu1, nl1-1
            trnx=0.0_WP
            trny=0.0_WP
            do k=1, nod_in_elem2D_num(node)
                elem=nod_in_elem2D(k,node)
                if( nz.LE.(nlevels(elem)-1) .and. nz.GE.(ulevels(elem))) then
                    trnx=trnx+tr_xy(1, nz, elem)*elem_area(elem)
                    trny=trny+tr_xy(2, nz, elem)*elem_area(elem)
                endif
            end do
            tr_xynodes(1, nz, node)=trnx/3.0_WP/areasvol(nz, node)
            tr_xynodes(2, nz, node)=trny/3.0_WP/areasvol(nz, node)
        end do
    end do
    
    ! add contribution of implicite vertical diffusion (after Klingbeil et al.2014)
    do node=1, myDim_nod2D
        nu1  = ulevels_nod2D(node) 
        nl1  = nlevels_nod2D(node)
        Dflx = 0.0_WP
        
        !_______________________________________________________________________
        ! compute diffusive flux at full surface levels: ---> Dflx[Tr^(n+1)]
        ! Dflx = Kv * dT/dz
        !   ---------------- zbar_1  -+
        !                             |
        !          o tr_1 ------------|- hnode_1 -+
        !                             |           |
        !   --Dflx-x-Tr----- zbar_2 --+           |- (hnode_1+hnod_2)/2
        !                             |           |
        !          o tr_2 ------------|- hnode_2 -+
        !                             |
        !   ---------------- zbar_3  -+
        !_______________________________________________________________________
        ! !!! ATTENTION: !!!
        ! Implicite vertical Diffusion is done with hnode^(n+1) which is stored 
        ! in the variable hnode (keep in mind at this point hnode_new was used
        ! to rescue the varaible hnode^(n) ).
        do nz=nu1+1, nl1-1
            ! But we saved hnode^(n)
            ! in the variable hnode_new (Although this sounds confusing!!!), so
            ! hnode_new contains here in reality hnode_old. So we recompute
            ! helem from hnode_new --> sum(hnode_new(nz, elnodes))/3.0_WP
            dzup     = hnode_new(nz-1, node)*0.5_WP
            dzdwn    = hnode_new(nz  , node)*0.5_WP
            
            Dflx(nz) = dzup  * (slope(1, nz-1, node)*tr_xynodes(1, nz-1, node) + slope(2, nz-1, node)*tr_xynodes(2, nz-1, node)) * Ki(nz-1, node) + &
                       dzdwn * (slope(1, nz  , node)*tr_xynodes(1, nz  , node) + slope(2, nz  , node)*tr_xynodes(2, nz  , node)) * Ki(nz  , node)
            Dflx(nz) = Dflx(nz) / (dzup+dzdwn) * area(nz, node)
            
            ! extent Dflx after Sergey to compute the diffusive component of DVD 
            ! see T. Banerjee et al. 2023 eq. 7
            ! Xchi_(i+0.5) = 2*Kv * (T^(n+1)_(i+1)-T^(n+!)_i)/dz * (Tstar_(i+1)-Tstar_(i)/dz
            ! Dflx(nz) = -Dflx * 2*(Tstar_(i+1)-Tstar_(i)/dz
            ! --> if Sergey diagonsitc is used Dflx correspond from here on to 
            !     Xchi at full depth level interface
            ! --> at this points tr must be trstar
            Dflx(nz) = Dflx(nz) * merge(-2.0_WP*( tr(nz, node)-tr(nz-1, node) ), 1.0_WP, do_SDdvd)
            !                       |
            !                       +-> merge: if do_SDdvd==True use first argument, if
            !                           do_SDdvd==False use second argument of merge!
        end do !--> do nz=nu1+1, nl1-1
       
        !_______________________________________________________________________
        ! decide between Sergeys and Knuts DVD method
        if (do_SDdvd) then 
            !___________________________________________________________________
            ! now compute flx of Xchi into control volume, see T. Banerjee et al. 2023
            ! eq. 9. Each control volume gets one half of the Xchi interface value 
            ! --> keep in mind we are now again on mid depth levels and in this case
            !     Dflx contains the interface value of Xchi
            ! --> compute volume normlized Xchi contribution
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + (dvd_xdfac*Dflx(nz)+(1.0_WP-dvd_xdfac)*Dflx(nz+1))/(areasvol(nz, node)*hnode_new(nz, node))
            end do !--> do nz=nu1+1, nl1-1
            
        else
            !___________________________________________________________________
            ! compute diffusive flux flxdiff(nz), after Klingbeil etal 2014
            ! flxdiff(nz)= ([Kv*dT/dz]_1*T_1*A_1-[Kv*dT/dz]_1*T_1*A_2)/V_1 into the volume 
            ! surface
            nz     = nu1
            tr_up  = Dflx(nz  )*tr(nz, node)
            tr_dwn = (tr(nz, node)*hnode_new(nz, node)+tr(nz+1, node)*hnode_new(nz+1, node))/(hnode_new(nz, node)+hnode_new(nz+1, node))
            tr_dwn = Dflx(nz+1)*tr_dwn
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk
            do nz=nu1+1, nl1-2
                ! volume (hnode) weighted temperature reconstruction at upper/lower
                ! scalar cell interface
                tr_up  = (tr(nz-1, node)*hnode_new(nz-1, node)+tr(nz  , node)*hnode_new(nz  , node))/(hnode_new(nz-1, node)+hnode_new(nz  , node))
                tr_dwn = (tr(nz  , node)*hnode_new(nz  , node)+tr(nz+1, node)*hnode_new(nz+1, node))/(hnode_new(nz  , node)+hnode_new(nz+1, node))
                ! compute  T^(n+1)_(i+0.5) * Dflx_(i+0.5) through upper and lower face
                tr_up  =  Dflx(nz  )*tr_up 
                tr_dwn =  Dflx(nz+1)*tr_dwn
                ! copute dvd contribution normalized with volume 
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |    |-> factor 2 comes here from Klingbeil et al.2014
                !                                      this is thes   |         (2*Tr^(n+1) * Dflx[Tr^(n+1)] 
                !                                   +----plus sign----+
                !                                   |
                !                                   v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
            end do !--> do nz=nu1+1, nl1-1
            
            !___________________________________________________________________
            ! bottom
            nz    = nl1-1
            tr_up = (tr(nz-1, node)*hnode_new(nz-1, node)+tr(nz  , node)*hnode_new(nz  , node))/(hnode_new(nz-1, node)+hnode_new(nz  , node))
            tr_up = Dflx(nz)*tr_up
            tr_dwn= 0.0_WP
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
        end if 
    end do !--> do node=1, myDim_nod2D
end subroutine dvd_add_difflux_vertexplredi
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_vertimplredi(do_SDdvd, tr_num, dvd_tot, tr, trstar, Ki, slope, partit, mesh)
    use g_cvmix_kpp, only: kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        integer       , intent(in)             :: tr_num
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: trstar( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: Ki(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: slope(3,mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, nz, nu1, nl1
        real(kind=WP)                          :: tr_up, tr_dwn, Dflx(mesh%nl), dzup, dzdwn, dz, Kv_R
        real(kind=WP)                          :: isredi=0.0_WP, sbc, rsss
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    ! add contribution of redi diffusivity when activated
    if (Redi) isredi=1.0_WP
    
    ! add contribution of implicite vertical diffusion (after Klingbeil et al.2014)
    do node=1, myDim_nod2D
        nu1  = ulevels_nod2D(node) 
        nl1  = nlevels_nod2D(node)
        Dflx = 0.0_WP
        
        !_______________________________________________________________________
        ! compute diffusive flux at full surface levels: ---> Dflx[Tr^(n+1)]
        ! Dflx = Kv * dT/dz
        !   ---------------- zbar_1  -+
        !                             |
        !          o tr_1 ------------|- hnode_1 -+
        !                             |           |
        !   --Dflx-x-Tr----- zbar_2 --+           |- (hnode_1+hnod_2)/2
        !                             |           |
        !          o tr_2 ------------|- hnode_2 -+
        !                             |
        !   ---------------- zbar_3  -+
        !_______________________________________________________________________
        ! !!! ATTENTION: !!!
        ! Implicite vertical Diffusion is done with hnode^(n+1) which is stored 
        ! in the variable hnode (keep in mind at this point hnode_new was used
        ! to rescue the varaible hnode^(n) ).
        
        ! keep in mind: nz = nu1 --> there is no diffusive flux through the 
        ! surface zlevel therfor we start with nu1+1
        do nz=nu1+1, nl1-1
            ! take into account isoneutral Redi diffusivity Kd*s^2 --> K_33 = Kv + Kd*s^2
            dzup     = hnode(nz-1, node)*0.5_WP
            dzdwn    = hnode(nz  , node)*0.5_WP
            Kv_R     = dzup /(dzup + dzdwn) * slope(3, nz-1, node)**2 * Ki(nz-1,node) + &
                       dzdwn/(dzup + dzdwn) * slope(3, nz  , node)**2 * Ki(nz  ,node)
            Kv_R     = Kv_R*isredi
            
            ! fdiff = Kv*dT/dz = Kv (T(nz-1) - T(nz))/ ((hnode(nz-1)+hnode(nz))/2)
            ! --> Dflx is until here the diffusive volume flux on full depth levels
            Dflx(nz) = -(Kv_R) * (tr(nz, node)-tr(nz-1, node)) / (dzup+dzdwn) * area(nz, node)
            
            ! extent Dflx after Sergey to compute the diffusive component of DVD 
            ! see T. Banerjee et al. 2023 eq. 7
            ! Xchi_(i+0.5) = 2*Kv * (T^(n+1)_(i+1)-T^(n+!)_i)/dz * (Tstar_(i+1)-Tstar_(i)/dz
            ! Dflx(nz) = -2*Dflx*(Tstar_(i+1)-Tstar_(i)/dz
            ! --> if Sergey diagonsitc is used Dflx correspond from here on to 
            
            !     Xchi at full depth level interface
            Dflx(nz) = Dflx(nz) * merge(-2.0_WP*( trstar(nz, node)-trstar(nz-1, node) ), 1.0_WP, do_SDdvd)
            !                       |
            !                       +-> merge: if do_SDdvd==True use first argument, if
            !                           do_SDdvd==False use second argument of merge!
        end do !--> do nz=nu1+1, nl1-1
        
        !_______________________________________________________________________
        ! decide between Sergeys and Knuts DVD method
        if (do_SDdvd) then 
            !___________________________________________________________________
            ! In Sergeys DVD diagnostic we do not take into account the surface 
            ! and interior fluxes of surface boundary condition and short wave 
            ! penetration, since akthough they are source terms for variance, 
            ! they not contribute to the decay variance !!!
            ! If they arte considered than only in a separate term but not added
            ! to the variance decay term of diffusivity. However in Knuts DVD |
            ! diagnostic this is different. Here the source term have to be explicitly 
            ! taken into account since they already contributed to (T^(n+1))^2-(T^(n))^2
            
            ! now compute flx of Xchi into control volume, see T. Banerjee et al. 2023
            ! eq. 9. Each control volume gets one half of the Xchi interface value 
            ! --> keep in mind we are now again on mid depth levels and in this case
            !     Dflx contains the interface value of Xchi
            ! --> compute volume normlized Xchi contribution
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + (dvd_xdfac*Dflx(nz)+(1.0_WP-dvd_xdfac)*Dflx(nz+1))/(areasvol(nz, node)*hnode(nz, node))
            end do !--> do nz=nu1+1, nl1-1
        else
            !___________________________________________________________________
            ! compute diffusive flux flxdiff(nz), after Klingbeil etal 2014
            ! flxdiff(nz)= ([Kv*dT/dz]_1*T_1*A_1-[Kv*dT/dz]_1*T_1*A_2)/V_1 into the volume 
            ! surface
            nz     = nu1
            tr_up  = Dflx(nz  )*tr(nz, node)
            tr_dwn = (tr(nz, node)*hnode(nz, node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz, node)+hnode(nz+1, node))
            tr_dwn = Dflx(nz+1)*tr_dwn
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk
            do nz=nu1+1, nl1-2
                ! volume (hnode) weighted temperature reconstruction at upper/lower
                ! scalar cell interface
                tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
                tr_dwn = (tr(nz  , node)*hnode(nz  , node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
                ! compute  T^(n+1)_(i+0.5) * Dflx_(i+0.5) through upper and lower face
                tr_up  =  Dflx(nz  )*tr_up 
                tr_dwn =  Dflx(nz+1)*tr_dwn
                ! copute dvd contribution normalized with volume 
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |    |-> factor 2 comes here from Klingbeil et al.2014
                !                                      this is the    |         (2*Tr^(n+1) * Dflx[Tr^(n+1)] 
                !                                   +----plus sign----+
                !                                   |
                !                                   v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
            end do !--> do nz=nu1+1, nl1-1
            
            !___________________________________________________________________
            ! bottom
            nz     = nl1-1
            tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
            tr_dwn = tr(nz, node)
            tr_up  = Dflx(nz  )*tr_up
            tr_dwn = 0.0_WP
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
        end if 
    end do !--> do node=1, myDim_nod2D
end subroutine dvd_add_difflux_vertimplredi
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_vertimpl(do_SDdvd, tr_num, dvd_tot, tr, trstar, Kv, partit, mesh)
    use g_cvmix_kpp, only: kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        integer       , intent(in)             :: tr_num
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: trstar( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: Kv(     mesh%nl  , partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, nz, nu1, nl1
        real(kind=WP)                          :: tr_up, tr_dwn, Dflx(mesh%nl), dzup, dzdwn, dz, Kv_R
        real(kind=WP)                          :: sbc, rsss
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    ! add contribution of implicite vertical diffusion (after Klingbeil et al.2014)
    do node=1, myDim_nod2D
        nu1  = ulevels_nod2D(node) 
        nl1  = nlevels_nod2D(node)
        Dflx = 0.0_WP
        
        !_______________________________________________________________________
        ! compute diffusive flux at full surface levels: ---> Dflx[Tr^(n+1)]
        ! Dflx = Kv * dT/dz
        !   ---------------- zbar_1  -+
        !                             |
        !          o tr_1 ------------|- hnode_1 -+
        !                             |           |
        !   --Dflx-x-Tr----- zbar_2 --+           |- (hnode_1+hnod_2)/2
        !                             |           |
        !          o tr_2 ------------|- hnode_2 -+
        !                             |
        !   ---------------- zbar_3  -+
        !_______________________________________________________________________
        ! !!! ATTENTION: !!!
        ! Implicite vertical Diffusion is done with hnode^(n+1) which is stored 
        ! in the variable hnode (keep in mind at this point hnode_new was used
        ! to rescue the varaible hnode^(n) ).
        
        ! keep in mind: nz = nu1 --> there is no diffusive flux through the 
        ! surface zlevel therfor we start with nu1+1
        do nz=nu1+1, nl1-1
            ! take into account isoneutral Redi diffusivity Kd*s^2 --> K_33 = Kv + Kd*s^2
            dzup     = hnode(nz-1, node)*0.5_WP
            dzdwn    = hnode(nz  , node)*0.5_WP
           
            ! fdiff = Kv*dT/dz = Kv (T(nz-1) - T(nz))/ ((hnode(nz-1)+hnode(nz))/2)
            ! --> Dflx is until here the diffusive volume flux on full depth levels
            Dflx(nz) = -(Kv(nz, node)) * (tr(nz, node)-tr(nz-1, node)) / (dzup+dzdwn) * area(nz, node)
            
            ! extent Dflx after Sergey to compute the diffusive component of DVD 
            ! see T. Banerjee et al. 2023 eq. 7
            ! Xchi_(i+0.5) = 2*Kv * (T^(n+1)_(i+1)-T^(n+!)_i)/dz * (Tstar_(i+1)-Tstar_(i)/dz
            ! Dflx(nz) = -2*Dflx*(Tstar_(i+1)-Tstar_(i)/dz
            ! --> if Sergey diagonsitc is used Dflx correspond from here on to 
            
            !     Xchi at full depth level interface
            Dflx(nz) = Dflx(nz) * merge(-2.0_WP*( trstar(nz, node)-trstar(nz-1, node) ), 1.0_WP, do_SDdvd)
            !                       |
            !                       +-> merge: if do_SDdvd==True use first argument, if
            !                           do_SDdvd==False use second argument of merge!
        end do !--> do nz=nu1+1, nl1-1
        
        !_______________________________________________________________________
        ! decide between Sergeys and Knuts DVD method
        if (do_SDdvd) then 
            !___________________________________________________________________
            ! In Sergeys DVD diagnostic we do not take into account the surface 
            ! and interior fluxes of surface boundary condition and short wave 
            ! penetration, since akthough they are source terms for variance, 
            ! they not contribute to the decay variance !!!
            ! If they arte considered than only in a separate term but not added
            ! to the variance decay term of diffusivity. However in Knuts DVD |
            ! diagnostic this is different. Here the source term have to be explicitly 
            ! taken into account since they already contributed to (T^(n+1))^2-(T^(n))^2
            
            ! now compute flx of Xchi into control volume, see T. Banerjee et al. 2023
            ! eq. 9. Each control volume gets one half of the Xchi interface value 
            ! --> keep in mind we are now again on mid depth levels and in this case
            !     Dflx contains the interface value of Xchi
            ! --> compute volume normlized Xchi contribution
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + (dvd_xdfac*Dflx(nz)+(1.0_WP-dvd_xdfac)*Dflx(nz+1))/(areasvol(nz, node)*hnode(nz, node))
            end do !--> do nz=nu1+1, nl1-1
        else

            !___________________________________________________________________
            ! compute diffusive flux flxdiff(nz), after Klingbeil etal 2014
            ! flxdiff(nz)= ([Kv*dT/dz]_1*T_1*A_1-[Kv*dT/dz]_1*T_1*A_2)/V_1 into the volume 
            ! surface
            nz     = nu1
            tr_up  = Dflx(nz  )*tr(nz, node)
            tr_dwn = (tr(nz, node)*hnode(nz, node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz, node)+hnode(nz+1, node))
            tr_dwn = Dflx(nz+1)*tr_dwn
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk
            do nz=nu1+1, nl1-2
                ! volume (hnode) weighted temperature reconstruction at upper/lower
                ! scalar cell interface
                tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
                tr_dwn = (tr(nz  , node)*hnode(nz  , node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
                ! compute  T^(n+1)_(i+0.5) * Dflx_(i+0.5) through upper and lower face
                tr_up  =  Dflx(nz  )*tr_up 
                tr_dwn =  Dflx(nz+1)*tr_dwn
                ! copute dvd contribution normalized with volume 
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |    |-> factor 2 comes here from Klingbeil et al.2014
                !                                      this is the    |         (2*Tr^(n+1) * Dflx[Tr^(n+1)] 
                !                                   +----plus sign----+
                !                                   |
                !                                   v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
            end do !--> do nz=nu1+1, nl1-1
            
            !___________________________________________________________________
            ! bottom
            nz     = nl1-1
            tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
            tr_dwn = tr(nz, node)
            tr_up  = Dflx(nz  )*tr_up
            tr_dwn = 0.0_WP
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
        end if 
    end do !--> do node=1, myDim_nod2D
end subroutine dvd_add_difflux_vertimpl
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux for surface boundary condition after 
! Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_sbc(do_SDdvd, tr_num, dvd_tot, tr, trstar, partit, mesh)
    use g_cvmix_kpp, only: kpp_nonlcltranspT, kpp_nonlcltranspS, kpp_oblmixc
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        integer       , intent(in)             :: tr_num
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: trstar( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, nz, nu1, nl1
        real(kind=WP)                          :: tr_up, tr_dwn, Dflx(mesh%nl), dzup, dzdwn, dz, Kv_R
        real(kind=WP)                          :: sbc, rsss
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    ! add contribution of implicite vertical diffusion (after Klingbeil et al.2014)
    do node=1, myDim_nod2D
        nu1  = ulevels_nod2D(node) 
        nl1  = nlevels_nod2D(node)
        Dflx = 0.0_WP
        
        !_______________________________________________________________________
        ! apply apply surface boundarie condition
        nz = nu1
        if     (tr_num==1) then 
            sbc = -(heat_flux(node)/vcpw + tr(nz, node)*water_flux(node)*is_nonlinfs*0.5_WP)
        elseif (tr_num==2) then 
                sbc =  (virtual_salt(node) + relax_salt(node) - (real_salt_flux(node) - tr(nz, node)*water_flux(node)*0.5_WP)*is_nonlinfs)
        end if 
        Dflx(nz) = Dflx(nz) + sbc*area(nz, node)
            
        !_______________________________________________________________________
        ! apply shortwave penetration fluxes
        if (use_sw_pene .and. tr_num==1) then
            nz=nu1
            Dflx(nz) = Dflx(nz) + sw_3d(nz, node)*area(nz,node)
            do nz=nu1+1, nl1-1
                Dflx(nz) = Dflx(nz) + sw_3d(nz, node)*area(nz,node) 
            end do
        end if
            
        !_______________________________________________________________________
        ! Add KPP nonlocal fluxes to the rhs (only T and S currently)
        ! use here blmc or kpp_oblmixc instead of Kv, since Kv already contains
        ! at this point the mixing enhancments from momix, instable
        ! mixing or windmixing which are to much for nonlocal
        ! transports and lead to instability of the model
        if (use_kpp_nonlclflx) then
            if (tr_num==2) then
                rsss=ref_sss
                if (ref_sss_local) rsss=tr(1,node)
            end if
            !___________________________________________________________________
            ! use fesom1.4 KPP
            if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then
                if     (tr_num==1) then ! temperature
                    do nz=nu1+1, nl1-1
                        Dflx(nz) = Dflx(nz) + MIN(ghats(nz, node)*blmc(nz, node, 2), 1.0_WP)*heat_flux(node)/vcpw*area(nz, node) 
                    end do
                elseif (tr_num==2) then ! salinity
                    do nz=nu1+1, nl1-1
                        Dflx(nz) = Dflx(nz) + MIN(ghats(nz, node)*blmc(nz, node, 3), 1.0_WP)*rsss*water_flux(node)*area(nz, node) 
                    end do
                end if
            !___________________________________________________________________
            ! use cvmix KPP
            elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then
                if     (tr_num==1) then ! temperature
                    do nz=nu1+1, nl1-1
                        Dflx(nz) = Dflx(nz) + MIN(kpp_nonlcltranspT(nz, node)*kpp_oblmixc(nz, node, 2), 1.0_WP)*heat_flux(node)/vcpw*area(nz, node) 
                    end do
                elseif (tr_num==2) then ! salinity
                    do nz=nu1+1, nl1-1
                        Dflx(nz) = Dflx(nz) + MIN(kpp_nonlcltranspT(nz, node)*kpp_oblmixc(nz, node, 3), 1.0_WP)*rsss*water_flux(node)*area(nz, node) 
                    end do
                end if    
            end if
        end if ! --> if (use_kpp_nonlclflx) then
        
        !_______________________________________________________________________
        if (do_SDdvd) then 
            ! --> if Sergey diagonsitc is used Dflx correspond from here on to 
            !     Xchi at full depth level interface
            nz = 1 
            Dflx(nz) = Dflx(nz) * -2.0_WP*( trstar(nz, node) )
            do nz=nu1+1, nl1-1
                Dflx(nz) = Dflx(nz) * -2.0_WP*( trstar(nz, node)-trstar(nz-1, node) )
            end do ! --> do nz=nu1+1, nl1-1
        end if 
        
        !_______________________________________________________________________
        ! decide between Sergeys and Knuts DVD method
        if (do_SDdvd) then 
            !___________________________________________________________________
            ! In Sergeys DVD diagnostic we do not take into account the surface 
            ! and interior fluxes of surface boundary condition and short wave 
            ! penetration, since akthough they are source terms for variance, 
            ! they not contribute to the decay variance !!!
            ! If they arte considered than only in a separate term but not added
            ! to the variance decay term of diffusivity. However in Knuts DVD |
            ! diagnostic this is different. Here the source term have to be explicitly 
            ! taken into account since they already contributed to (T^(n+1))^2-(T^(n))^2
            
            ! now compute flx of Xchi into control volume, see T. Banerjee et al. 2023
            ! eq. 9. Each control volume gets one half of the Xchi interface value 
            ! --> keep in mind we are now again on mid depth levels and in this case
            !     Dflx contains the interface value of Xchi
            ! --> compute volume normlized Xchi contribution
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + (dvd_xdfac*Dflx(nz)+(1.0_WP-dvd_xdfac)*Dflx(nz+1))/(areasvol(nz, node)*hnode(nz, node))
            end do !--> do nz=nu1+1, nl1-1
        else
            !___________________________________________________________________
            ! compute diffusive flux flxdiff(nz), after Klingbeil etal 2014
            ! flxdiff(nz)= ([Kv*dT/dz]_1*T_1*A_1-[Kv*dT/dz]_1*T_1*A_2)/V_1 into the volume 
            ! surface
            nz     = nu1
            tr_up  = Dflx(nz  )*tr(nz, node)
            tr_dwn = (tr(nz, node)*hnode(nz, node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz, node)+hnode(nz+1, node))
            tr_dwn = Dflx(nz+1)*tr_dwn
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
            
            !___________________________________________________________________
            ! bulk
            do nz=nu1+1, nl1-2
                ! volume (hnode) weighted temperature reconstruction at upper/lower
                ! scalar cell interface
                tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
                tr_dwn = (tr(nz  , node)*hnode(nz  , node)+tr(nz+1, node)*hnode(nz+1, node))/(hnode(nz  , node)+hnode(nz+1, node))
                ! compute  T^(n+1)_(i+0.5) * Dflx_(i+0.5) through upper and lower face
                tr_up  =  Dflx(nz  )*tr_up 
                tr_dwn =  Dflx(nz+1)*tr_dwn
                ! copute dvd contribution normalized with volume 
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
                !                                                     |    |-> factor 2 comes here from Klingbeil et al.2014
                !                                      this is the    |         (2*Tr^(n+1) * Dflx[Tr^(n+1)] 
                !                                   +----plus sign----+
                !                                   |
                !                                   v
                !  Xchi = -d/dt(Tr^2) - div(v*Tr^2) + div(Kv*nabla*Tr^2)
                !   |
                !   v
                !  Kv*(nabla*Tr)^2
            end do !--> do nz=nu1+1, nl1-1
            
            !___________________________________________________________________
            ! bottom
            nz     = nl1-1
            tr_up  = (tr(nz-1, node)*hnode(nz-1, node)+tr(nz  , node)*hnode(nz  , node))/(hnode(nz-1, node)+hnode(nz  , node))
            tr_dwn = tr(nz, node)
            tr_up  = Dflx(nz  )*tr_up
            tr_dwn = 0.0_WP
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*(tr_up-tr_dwn)/(areasvol(nz, node)*hnode(nz, node))
        end if 
    end do !--> do node=1, myDim_nod2D
end subroutine dvd_add_difflux_sbc
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx_ver[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx_ver[...] are reconstructed values at the interfase
subroutine dvd_add_difflux_bhvisc(do_SDdvd, tr_num, dvd_tot, tr, trstar, gamma0_tra, gamma1_tra, gamma2_tra, dump, partit, mesh)
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        integer       , intent(in)             :: tr_num
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: trstar( mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        real(kind=WP) , intent(in)             :: gamma0_tra
        real(kind=WP) , intent(in)             :: gamma1_tra
        real(kind=WP) , intent(in)             :: gamma2_tra
        real(kind=WP) , intent(inout)          :: dump(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, edge, nz, nu1, nl1
        integer                                :: ednodes(2), edelem(2), elnodes_l(3), elnodes_r(3)
        real(kind=WP)                          :: len, du, dv, vi, dtr(mesh%nl-1), trc(mesh%nl-1)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"    
    
    !___________________________________________________________________________
    ! reuse here an already allocated working array --> initialise first 
    do node=1, myDim_nod2D+eDim_nod2D
        dump(:, node)=0.0_WP
    end do

    !___________________________________________________________________________
    ! first round 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, nz, ednodes, edelem, elnodes_l, elnodes_r, &
!$OMP nu1, nl1, du, dv, dt, len, vi)
!$OMP DO    
    do edge=1, myDim_edge2D!+eDim_edge2D
        ! skip boundary edges only consider inner edges 
        if (myList_edge2D(edge) > edge2D_in) cycle
        edelem   = edge_tri(:,edge)
        ednodes  = edges(:,edge)
        len      = sqrt(sum(elem_area(edelem)))
        nl1      = minval(nlevels(edelem))
        nu1      = maxval(ulevels(edelem))
        elnodes_l= elem2d_nodes(:, edelem(1))
        elnodes_r= elem2d_nodes(:, edelem(2))
        do nz=nu1, nl1-1
            du     = maxval(tr(nz, elnodes_l))-minval(tr(nz, elnodes_r))
            dv     = minval(tr(nz, elnodes_l))-maxval(tr(nz, elnodes_r))
            vi     = du*du+dv*dv
            dtr(nz)= tr(nz, ednodes(1))-tr(nz, ednodes(2))
            vi     = sqrt(max(gamma0_tra,           &
                          max(gamma1_tra*sqrt(vi),   &
                              gamma2_tra*     vi)    &
                            )*len)
            dtr(nz)=dtr(nz)*vi
        END DO
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_set_lock  (partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif
       dump(nu1:nl1-1, ednodes(1)) = dump(nu1:nl1-1, ednodes(1))-dtr(nu1:nl1-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(ednodes(1)))
       call omp_set_lock  (partit%plock(ednodes(2)))
#endif
       dump(nu1:nl1-1, ednodes(2)) = dump(nu1:nl1-1, ednodes(2))+dtr(nu1:nl1-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
       call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif
    end do !--> do edge=1, myDim_edge2D!+eDim_edge2D
!$OMP END DO
!$OMP MASTER
    call exchange_nod(dump, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! second round:
!$OMP DO
    do edge=1, myDim_edge2D!+eDim_edge2D
        ! skip boundary edges only consider inner edges 
        if (myList_edge2D(edge) > edge2D_in) cycle
        edelem   = edge_tri(:,edge)
        ednodes  = edges(:,edge)
        len      = sqrt(sum(elem_area(edelem)))
        nl1      = minval(nlevels(edelem))
        nu1      = maxval(ulevels(edelem))
        elnodes_l= elem2d_nodes(:, edelem(1))
        elnodes_r= elem2d_nodes(:, edelem(2))
        do nz=nu1, nl1-1
            du    = maxval(tr(nz, elnodes_l))-minval(tr(nz, elnodes_r))
            dv    = minval(tr(nz, elnodes_l))-maxval(tr(nz, elnodes_r))
            vi    = du*du+dv*dv
            dtr(nz)= dump(nz, ednodes(1))-dump(nz, ednodes(2))
            vi    = sqrt(max(gamma0_tra,            &
                         max(gamma1_tra*sqrt(vi),   &
                             gamma2_tra*     vi)    &
                             )*len)
            dtr(nz)=-dtr(nz)*vi
        end do !-->do nz=nu1, nl1-1  
        
        !_______________________________________________________________________
        ! make the difference between Sergeys and Knuts method
        if (do_SDdvd) then
            do nz=nu1, nl1-1
                !_______________________________________________________________
                ! compute Tstar as in T. Banerjee et al. 2023 --> Tstar = (Tr^(n+1) + Tr^n)/2
                ! compute net edge Tracer contribution with Tstar
                trc(nz) = 2.0_WP*( trstar(nz, ednodes(1)) - trstar(nz, ednodes(2)) )
            end do !-->do nz=nu1, nl1-1
            !___________________________________________________________________
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif
            dvd_tot(nu1:nl1-1, ednodes(1), tr_num) = dvd_tot(nu1:nl1-1, ednodes(1), tr_num) + & 
                                    trc(nu1:nl1-1)*dtr(nu1:nl1-1)/(areasvol(nu1:nl1-1,ednodes(1))) 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(1)))
            call omp_set_lock  (partit%plock(ednodes(2)))
#endif
            dvd_tot(nu1:nl1-1, ednodes(2), tr_num) = dvd_tot(nu1:nl1-1, ednodes(2), tr_num) - &
                                    trc(nu1:nl1-1)*dtr(nu1:nl1-1)/(areasvol(nu1:nl1-1,ednodes(2))) 
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif
        else
            do nz=nu1, nl1-1
                !_______________________________________________________________
                ! compute edge centered tracer value see Klingbeil et al 2014
                trc(nz) = tr(nz, ednodes(1)) * hnode(nz, ednodes(1)) + tr(nz, ednodes(2)) * hnode(nz, ednodes(2))
                trc(nz) = -2.0_WP*trc(nz)/(hnode(nz, ednodes(1))+hnode(nz, ednodes(2)))
                !          |
                !          +-> factor 2.0 comes here from Klingbeil et al. 2014
                !              Dflux[ tr^2] = 2*tr*Dflux[tr]          
            end do !-->do nz=nu1, nl1-1
            
            !___________________________________________________________________
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_set_lock  (partit%plock(ednodes(1)))
#else
!$OMP ORDERED
#endif
            dvd_tot(nu1:nl1-1, ednodes(1), tr_num) = dvd_tot(nu1:nl1-1, ednodes(1), tr_num) + & 
                                        trc(nu1:nl1-1)*dtr(nu1:nl1-1)  / ( areasvol(nu1:nl1-1,ednodes(1)) )                                     
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(1)))
            call omp_set_lock  (partit%plock(ednodes(2)))
#endif
            dvd_tot(nu1:nl1-1, ednodes(2), tr_num) = dvd_tot(nu1:nl1-1, ednodes(2), tr_num) - &
                                        trc(nu1:nl1-1)*dtr(nu1:nl1-1)  / ( areasvol(nu1:nl1-1,ednodes(2)) )
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
            call omp_unset_lock(partit%plock(ednodes(2)))
#else
!$OMP END ORDERED
#endif
        end if 
        
    end do !--> do edge=1, myDim_edge2D!+eDim_edge2D
!$OMP END DO
!$OMP END PARALLEL

end subroutine dvd_add_difflux_bhvisc    
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx_ver[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx_ver[...] are reconstructed values at the interfase
subroutine dvd_add_clim_relax_channel(do_SDdvd, tr_num, dvd_tot, tr, partit, mesh)
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        integer       , intent(in)             :: tr_num
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, nz, nu1, nl1, nn, nn1
        real(kind=WP)                          :: yy, a, Tzon
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"  
    
    !___________________________________________________________________________
    do node=1, myDim_nod2D+eDim_nod2D
        nu1      = ulevels_nod2D(node)
        nl1      = nlevels_nod2D(node)
        
        !_______________________________________________________________________
        yy= coord_nod2D(2, node)-lat0
        a = 0 
        if (yy<dy/2) then  ! southward from the center of the first bin
            nn=1
            nn1=1
        else 
            nn=floor(yy/dy-0.5)+1
            nn1=nn+1
            if (nn1>100) nn1=nn  ! northward of the center of the last bin 
                                ! Linear interpolation (nearest if close to boundary)
                a=yy/dy+0.5-nn
        end if
        
        !_______________________________________________________________________
        do nz=nu1, nl1-1
            Tzon=(1.0-a)*ztem(nz, nn)+a*ztem(nz, nn1)
            dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*tr(nz, node)*tau_inv*(Tclim(nz, node)-Tzon)
        end do
    end do
end subroutine dvd_add_clim_relax_channel
!
!
!_______________________________________________________________________________
! add implicite vertical diffusive flux after Klingbeil et al. 2014, see between eq. 18 
! ans eq. 20 --> or eq. 24
! Xchi^(n+1) =  ...+ (2*Tr^(n+1) * Dflx_ver[Tr^(n+1)] )/ V^(n+1) +...
! --> here Tr^(n+1) und Dflx_ver[...] are reconstructed values at the interfase
subroutine dvd_add_clim_relax(do_SDdvd, tr_num, dvd_tot, tr, partit, mesh)
    implicit none
        type(t_partit), intent(inout), target  :: partit
        type(t_mesh)  , intent(in)   , target  :: mesh
        logical       , intent(in)             :: do_SDdvd ! if Sergey DVD==1, if Knut DvD==0
        integer       , intent(in)             :: tr_num
        real(kind=WP) , intent(inout)          :: dvd_tot(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D, 2)
        real(kind=WP) , intent(in)             :: tr(     mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
        integer                                :: node, nz, nu1, nl1
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"  
    
    !___________________________________________________________________________
    do node=1, myDim_nod2D+eDim_nod2D
        nu1      = ulevels_nod2D(node)
        nl1      = nlevels_nod2D(node)
        
        !_______________________________________________________________________
        if     (tr_num == 1 ) then
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*tr(nz, node)*relax2clim(node)*(Tclim(nz, node)-tr(nz, node))
            end do
        !_______________________________________________________________________    
        elseif (tr_num == 2 ) then
            do nz=nu1, nl1-1
                dvd_tot(nz, node, tr_num) = dvd_tot(nz, node, tr_num) + 2.0_WP*tr(nz, node)*relax2clim(node)*(Sclim(nz, node)-tr(nz, node))
            end do
        end if    
    end do
end subroutine dvd_add_clim_relax

end module diagnostics
