module diagnostics

  use g_config
  use mod_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use i_ARRAYS
  use o_mixing_KPP_mod
  use g_rotate_grid
  use g_support
  implicit none

  private
!!PS   
  public :: ldiag_solver, lcurt_stress_surf, ldiag_energy, ldiag_dMOC, ldiag_DVD, ldiag_forc, ldiag_salt3D, ldiag_curl_vel3, diag_list, &
            compute_diagnostics, rhs_diag, curl_stress_surf, curl_vel3, wrhof, rhof, &
            u_x_u, u_x_v, v_x_v, v_x_w, u_x_w, dudx, dudy, dvdx, dvdy, dudz, dvdz, utau_surf, utau_bott, av_dudz_sq, av_dudz, av_dvdz, stress_bott, u_surf, v_surf, u_bott, v_bott, &
            std_dens_min, std_dens_max, std_dens_N, std_dens, std_dens_UVDZ, std_dens_RHOZ, &
            compute_diag_dvd_2ndmoment_klingbeil_etal_2014, compute_diag_dvd_2ndmoment_burchard_etal_2008, compute_diag_dvd
  ! Arrays used for diagnostics, some shall be accessible to the I/O
  ! 1. solver diagnostics: A*x=rhs? 
  ! A=ssh_stiff, x=d_eta, rhs=ssh_rhs; rhs_diag=A*x;
  real(kind=WP),  save, allocatable, target      :: rhs_diag(:)
  real(kind=WP),  save, allocatable, target      :: curl_stress_surf(:)
  real(kind=WP),  save, allocatable, target      :: curl_vel3(:,:)
  real(kind=WP),  save, allocatable, target      :: wrhof(:,:), rhof(:,:)
  real(kind=WP),  save, allocatable, target      :: u_x_u(:,:), u_x_v(:,:), v_x_v(:,:), v_x_w(:,:), u_x_w(:,:)
  real(kind=WP),  save, allocatable, target      :: dudx(:,:), dudy(:,:), dvdx(:,:), dvdy(:,:), dudz(:,:), dvdz(:,:), av_dudz(:,:), av_dvdz(:,:), av_dudz_sq(:,:)
  real(kind=WP),  save, allocatable, target      :: utau_surf(:), utau_bott(:)
  real(kind=WP),  save, allocatable, target      :: stress_bott(:,:), u_bott(:), v_bott(:), u_surf(:), v_surf(:)

! defining a set of standard density bins which will be used for computing densMOC
! integer,        parameter                      :: std_dens_N  = 100
! real(kind=WP),  save, target                   :: std_dens(std_dens_N)
  integer,        parameter                      :: std_dens_N  =72
  real(kind=WP),  save, target                   :: std_dens(std_dens_N)=(/ &
                                                            0.0000,   30.00000, 30.55556, 31.11111, 31.66667, 32.22222, &
                                                            32.77778, 33.33333, 33.88889, 34.44444, 35.00000, 35.10622, &
                                                            35.20319, 35.29239, 35.37498, 35.45187, 35.52380, 35.59136, &
                                                            35.65506, 35.71531, 35.77247, 35.82685, 35.87869, 35.92823, &
                                                            35.97566, 36.02115, 36.06487, 36.10692, 36.14746, 36.18656, &
                                                            36.22434, 36.26089, 36.29626, 36.33056, 36.36383, 36.39613, &
                                                            36.42753, 36.45806, 36.48778, 36.51674, 36.54495, 36.57246, &
                                                            36.59932, 36.62555, 36.65117, 36.67621, 36.70071, 36.72467, &
                                                            36.74813, 36.77111, 36.79363, 36.81570, 36.83733, 36.85857, &
                                                            36.87940, 36.89985, 36.91993, 36.93965, 36.95904, 36.97808, &
                                                            36.99682, 37.01524, 37.03336, 37.05119, 37.06874, 37.08602, &
                                                            37.10303, 37.11979, 37.13630, 37.15257, 37.16861, 37.18441/)
  real(kind=WP),  save, target                   :: std_dd(std_dens_N-1)
  real(kind=WP),  save, target                   :: std_dens_min=1012., std_dens_max=1039.
  real(kind=WP),  save, allocatable, target      :: std_dens_UVDZ(:,:,:), std_dens_RHOZ(:,:)

  logical                                       :: ldiag_solver     =.false.
  logical                                       :: lcurt_stress_surf=.false.
  logical                                       :: ldiag_curl_vel3  =.false.
  logical                                       :: ldiag_energy     =.false.
  logical                                       :: ldiag_salt3D     =.false.
  ! this option activates writing the horizintal velocity transports within the density bins (U_rho_x_DZ and V_rho_x_DZ)
  ! an additional field (RHO_Z) will be computed which allows for diagnosing the numerical diapycnal mixing after A. Megann 2018
  logical                                       :: ldiag_dMOC       =.false.
  
  ! flag for calculating the Discrete Variance Decay --> estimator for numerical/
  ! spurious mixing in the advection schemes
  logical                                       :: ldiag_DVD        =.false.
  
  logical                                       :: ldiag_forc       =.false.
  
  namelist /diag_list/ ldiag_solver, lcurt_stress_surf, ldiag_curl_vel3, ldiag_energy, &
                       ldiag_dMOC, ldiag_DVD, ldiag_salt3D, ldiag_forc
  
  contains

! ==============================================================
!rhs_diag=ssh_rhs?
subroutine diag_solver(mode, mesh)
  implicit none
  integer, intent(in)           :: mode
  integer                       :: n, is, ie
  logical, save                 :: firstcall=.true.
  type(t_mesh), intent(in)     , target :: mesh
#include "associate_mesh.h"
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
subroutine diag_curl_stress_surf(mode, mesh)
  implicit none
  integer, intent(in)        :: mode
  logical, save              :: firstcall=.true.
  integer                    :: enodes(2), el(2), ed, n
  real(kind=WP)              :: deltaX1, deltaY1, deltaX2, deltaY2, c1
  type(t_mesh), intent(in)  , target :: mesh
!=====================
#include "associate_mesh.h"

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
subroutine diag_curl_vel3(mode, mesh)
  implicit none
  integer, intent(in)        :: mode
  logical, save              :: firstcall=.true.
  integer                    :: enodes(2), el(2), ed, n, nz, nl1, nl2
  real(kind=WP)              :: deltaX1, deltaY1, deltaX2, deltaY2, c1
  type(t_mesh), intent(in)  , target :: mesh

#include "associate_mesh.h"

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
subroutine diag_energy(mode, mesh)
  implicit none
  integer, intent(in)        :: mode
  type(t_mesh), intent(in)  , target :: mesh
  logical, save              :: firstcall=.true.
  integer                    :: n, nz, k, i, elem, nzmax, elnodes(3)
  integer                    :: iup, ilo
  real(kind=WP)              :: ux, vx, uy, vy, tvol, rval(2)
  real(kind=WP)              :: geo_grad_x(3), geo_grad_y(3), geo_u(3), geo_v(3)

#include "associate_mesh.h"
!=====================
  if (firstcall) then  !allocate the stuff at the first call
     allocate(wrhof(nl, myDim_nod2D), rhof(nl, myDim_nod2D))
     allocate(u_x_u(nl-1, myDim_nod2D), u_x_v(nl-1, myDim_nod2D), v_x_v(nl-1, myDim_nod2D), v_x_w(nl-1, myDim_elem2D), u_x_w(nl-1, myDim_elem2D))
     allocate(dudx(nl-1, myDim_nod2D), dudy(nl-1, myDim_nod2D), dvdx(nl-1, myDim_nod2D), dvdy(nl-1, myDim_nod2D), dudz(nl, myDim_elem2D), dvdz(nl, myDim_elem2D))
     allocate(utau_surf(myDim_elem2D), utau_bott(myDim_elem2D), av_dudz_sq(nl, myDim_elem2D), av_dudz(nl, myDim_elem2D), av_dvdz(nl, myDim_elem2D))
     allocate(u_surf(myDim_elem2D), v_surf(myDim_elem2D), u_bott(myDim_elem2D), v_bott(myDim_elem2D), stress_bott(2, myDim_elem2D))
     rhof  =0.
     wrhof=0.
     u_x_u=0.
     u_x_v=0.
     v_x_v=0.
     u_x_w=0.
     v_x_w=0.
     dudx=0.
     dudy=0.
     dvdx=0.
     dvdy=0.
     dudz=0.
     dvdz=0.
     av_dudz_sq=0.
     av_dudz=0.
     av_dvdz=0.

     u_surf=0.
     v_surf=0.

     u_bott=0.
     v_bott=0.

     stress_bott=0.

     utau_surf=0.
     utau_bott=0.

     firstcall=.false.
     if (mode==0) return
  end if
  
  u_x_u=Unode(1,1:nl-1,1:myDim_nod2D)*Unode(1,1:nl-1,1:myDim_nod2D)
  u_x_v=Unode(1,1:nl-1,1:myDim_nod2D)*Unode(2,1:nl-1,1:myDim_nod2D)
  v_x_v=Unode(2,1:nl-1,1:myDim_nod2D)*Unode(2,1:nl-1,1:myDim_nod2D)
  ! this loop might be very expensive
  DO n=1, myDim_elem2D
     nzmax =nlevels(n)
     zbar_n=0.0_WP
     Z_n   =0.0_WP
     ! in case of partial cells zbar_n(nzmax) is not any more at zbar(nzmax), 
     ! zbar_n(nzmax) is now zbar_e_bot(n), 
     zbar_n(nzmax)=zbar_e_bot(n)
     Z_n(nzmax-1)=zbar_n(nzmax) + helem(nzmax-1,n)/2.0_WP
     do nz=nzmax-1,2,-1
        zbar_n(nz) = zbar_n(nz+1) + helem(nz,n)
        Z_n(nz-1) = zbar_n(nz) + helem(nz-1,n)/2.0_WP
     end do
     zbar_n(1) = zbar_n(2) + helem(1,n)
     !compute du/dz & dv/dz
     dudz(2:nzmax-1, n)=(UV(1, 1:nzmax-2, n)-UV(1, 2:nzmax-1, n))/(Z_n(1:nzmax-2)-Z_n(2:nzmax-1)) !central differences in vertical
     dvdz(2:nzmax-1, n)=(UV(2, 1:nzmax-2, n)-UV(2, 2:nzmax-1, n))/(Z_n(1:nzmax-2)-Z_n(2:nzmax-1))


     rval=-C_d*sqrt(UV(1, nzmax-1,n)**2+ UV(2, nzmax-1, n)**2)*UV(:,nzmax-1,n)
     
     dudz(1, n)  =0.!stress_surf(1,n)/density_0/Av(2, n)
     dvdz(1, n)  =0.!stress_surf(2,n)/density_0/Av(2, n)
     dudz(nzmax, n) =0.!rval(1)/Av(nzmax-1, n)
     dvdz(nzmax, n) =0.!rval(2)/Av(nzmax-1, n)

     !compute int(Av * (du/dz)^2)
     av_dudz_sq(1:nzmax, n)=(dudz(1:nzmax, n)**2+dvdz(1:nzmax, n)**2)*Av(1:nzmax, n)
     av_dudz   (1:nzmax, n)= dudz(1:nzmax, n)*Av(1:nzmax, n)
     av_dvdz   (1:nzmax, n)= dvdz(1:nzmax, n)*Av(1:nzmax, n)

     utau_surf(n)=sum(stress_surf(:,n)/density_0*UV(:,1,n))
     utau_bott(n)=sum(rval*UV(:,nzmax-1,n))    !a scalar product tau_bottom times u 

     stress_bott(:,n)=rval

     u_surf(n)=UV(1,1,n)
     v_surf(n)=UV(2,1,n)

     u_bott(n)=UV(1,nzmax-1,n)
     v_bott(n)=UV(2,nzmax-1,n)
     
     elnodes=elem2D_nodes(:,n)
     DO nz=1, nzmax-1
        iup=max(nz-1, 1)
        ilo=min(nz, nzmax-1)
        u_x_w(nz,n)=sum(Wvel(nz, elnodes))/3.*(UV(1, iup, n)*helem(iup ,n)+UV(1, ilo, n)*helem(ilo,n))/(helem(iup ,n)+helem(ilo ,n))
        v_x_w(nz,n)=sum(Wvel(nz, elnodes))/3.*(UV(2, iup, n)*helem(iup ,n)+UV(2, ilo, n)*helem(ilo,n))/(helem(iup ,n)+helem(ilo ,n))
     END DO
  END DO
  ! this loop might be very expensive
  DO n=1, myDim_nod2D
     nzmax=nlevels_nod2D(n)
     ! compute Z_n which is the mid depth of prisms (ALE supports changing layer thicknesses)
     zbar_n=0.0_WP
     Z_n   =0.0_WP
     zbar_n(nzmax) =zbar_n_bot(n)
     rhof(nzmax,n) =density_m_rho0(nzmax-1, n)
     rhof(1,n)     =density_m_rho0(1, n)

     Z_n(nzmax-1) =zbar_n(nzmax)  + hnode_new(nzmax-1,n)/2.0_WP
     do nz=nzmax-1,2,-1
        zbar_n(nz) = zbar_n(nz+1) + hnode_new(nz,n)
        Z_n(nz-1)  = zbar_n(nz)   + hnode_new(nz-1,n)/2.0_WP
        rhof(nz,n) = (hnode_new(nz,n)*density_m_rho0(nz, n)+hnode_new(nz-1,n)*density_m_rho0(nz-1, n))/(hnode_new(nz,n)+hnode_new(nz-1,n))
     end do
     zbar_n(1)         = zbar_n(2) + hnode_new(1,n)
     wrhof(1:nzmax, n) = rhof(1:nzmax, n)*Wvel(1:nzmax, n)

     DO nz=1, nzmax-1
        tvol=0.0_WP
        ux  =0.0_WP
        uy  =0.0_WP
        vx  =0.0_WP
        vy  =0.0_WP
        DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           if (nlevels(elem)-1 < nz) cycle
	   elnodes=elem2D_nodes(:, elem)
           tvol=tvol+elem_area(elem)
           ux=ux+sum(gradient_sca(1:3,elem)*Unode(1,nz,elnodes))*elem_area(elem)         !accumulate tensor of velocity derivatives
           vx=vx+sum(gradient_sca(1:3,elem)*Unode(2,nz,elnodes))*elem_area(elem)
           uy=uy+sum(gradient_sca(4:6,elem)*Unode(1,nz,elnodes))*elem_area(elem)
           vy=vy+sum(gradient_sca(4:6,elem)*Unode(2,nz,elnodes))*elem_area(elem)
        END DO
        dudx(nz,n)=ux/tvol!/area(nz, n)/3.
        dvdx(nz,n)=vx/tvol
        dudy(nz,n)=uy/tvol
        dvdy(nz,n)=vy/tvol
     END DO
  END DO
end subroutine diag_energy
! ==============================================================
subroutine diag_densMOC(mode, mesh)
  implicit none
  integer, intent(in)                 :: mode
  integer                             :: nz, snz, elem, nzmax, elnodes(3), is, ie, pos
  integer                             :: jj
  real(kind=WP), save                 :: dd
  real(kind=WP)                       :: uvdz_el(2), rhoz_el, dz, weight, dmin, dmax, ddiff, test, test1, test2, test3
  real(kind=WP), save, allocatable    :: dens(:), aux(:)
  real(kind=WP), save, allocatable    :: std_dens_w(:,:)
  logical, save                       :: firstcall=.true.
  type(t_mesh), intent(in)           , target :: mesh

#include "associate_mesh.h"
!=====================


  if (firstcall) then !allocate the stuff at the first call
     allocate(std_dens_UVDZ(2,std_dens_N, myDim_elem2D))
     allocate(std_dens_RHOZ(  std_dens_N, myDim_elem2D))
     allocate(std_dens_w   (  std_dens_N, myDim_elem2D))
     allocate(aux(nl-1))
     allocate(dens(nl))
!
!std_dens(1)=27.5
!do nz=2, std_dens_N
!std_dens(nz)=std_dens(nz-1)+10.5/real(std_dens_N)
!end do
!
     std_dd(:)=std_dens(2:)-std_dens(:std_dens_N-1)
     dens         =0.
     std_dens_UVDZ=0.
     std_dens_RHOZ=0.
     firstcall=.false.
     if (mode==0) return
  end if

  std_dens_UVDZ=0.
  std_dens_RHOZ=0.
  std_dens_w   =0.
  do elem=1, myDim_elem2D
     elnodes=elem2D_nodes(:,elem)    
     nzmax =nlevels(elem)
     do nz=1, nzmax-1
       aux(nz)=sum(density_dmoc(nz, elnodes))/3.-1000.
     end do
     do nz=nzmax-1,2,-1
        dens(nz)   = (aux(nz)     * helem(nz-1,elem)+&
                      aux(nz-1)   * helem(nz,  elem))/sum(helem(nz-1:nz,elem))
     end do
     dens(nzmax)=dens(nzmax-1)+(dens(nzmax-1)-dens(nzmax-2))*helem(nzmax-1,elem)/helem(nzmax-2,elem)
     dens(1)    =dens(2)      +(dens(2)-dens(3))            *helem(1, elem)/helem(2,elem)

     do nz=nzmax-1,1,-1
        dmin=minval(dens(nz:nz+1))
        dmax=maxval(dens(nz:nz+1))
!       is=findloc(std_dens > dmin, value=.true., dim=1)
	is=1
        do jj = 1, std_dens_N
           if (std_dens(jj) > dmin) then
              is = jj
              exit
           endif
        end do

!       ie=findloc(std_dens < dmax, value=.true., back=.true., dim=1)
	ie=std_dens_N
        do jj = std_dens_N,1,-1
           if (std_dens(jj) < dmin) then
              ie = jj
              exit
           endif
        end do
        if (std_dens(is)>=dmax) is=ie
        if (std_dens(ie)<=dmin) ie=is
        uvdz_el=(UV(:,nz,elem)+fer_uv(:,nz,elem))*helem(nz,elem)
        rhoz_el=(dens(nz)-dens(nz+1))/helem(nz,elem)
        ddiff=abs(dens(nz)-dens(nz+1))
        if (ie-is > 0) then
           weight=(std_dens(is)-dmin)+std_dd(is)/2.
           weight=max(weight, 0.)/ddiff
           std_dens_UVDZ(:, is, elem)=std_dens_UVDZ(:, is, elem)+weight*uvdz_el
           std_dens_RHOZ(   is, elem)=std_dens_RHOZ(   is, elem)+weight*rhoz_el
           std_dens_w(   is, elem)   =std_dens_w(   is, elem)   +weight
           do snz=is+1, ie-1
              weight=(sum(std_dd(snz-1:snz))/2.)/ddiff
              std_dens_UVDZ(:, snz, elem)=std_dens_UVDZ(:, snz, elem)+weight*uvdz_el
              std_dens_RHOZ(   snz, elem)=std_dens_RHOZ(   snz, elem)+weight*rhoz_el
              std_dens_w   (   snz, elem)=std_dens_w   (   snz, elem)+weight
           end do
           weight=(dmax-std_dens(ie))+std_dd(ie-1)/2.
           weight=max(weight, 0.)/ddiff
           std_dens_UVDZ(:, ie, elem)=std_dens_UVDZ(:, ie, elem)+weight*uvdz_el
           std_dens_RHOZ(   ie, elem)=std_dens_RHOZ(   ie, elem)+weight*rhoz_el
           std_dens_w   (   ie, elem)=std_dens_w   (   ie, elem)+weight
        else
           std_dens_UVDZ(:, is, elem)=std_dens_UVDZ(:, is, elem)+uvdz_el
           std_dens_RHOZ(   is, elem)=std_dens_RHOZ(   is, elem)+rhoz_el
           std_dens_w   (   is, elem)=std_dens_w   (   is, elem)+1._wp
        end if
     end do
  end do
  where (std_dens_w > 0.)
        std_dens_RHOZ=std_dens_RHOZ/std_dens_w
  end where
end subroutine diag_densMOC
! ==============================================================
subroutine compute_diagnostics(mode, mesh)
  implicit none
  integer, intent(in)           :: mode !constructor mode (0=only allocation; any other=do diagnostic)
  real(kind=WP)                 :: val
  type(t_mesh), intent(in)  , target :: mesh
  !1. solver diagnostic
  if (ldiag_solver)      call diag_solver(mode, mesh)
  !2. compute curl(stress_surf)
  if (lcurt_stress_surf) call diag_curl_stress_surf(mode, mesh)
  !3. compute curl(velocity)
  if (ldiag_curl_vel3)   call diag_curl_vel3(mode, mesh)
  !4. compute energy budget
  if (ldiag_energy)      call diag_energy(mode, mesh)
  !5. print integrated temperature 
  if (ldiag_salt3d) then
     if (mod(mstep,logfile_outfreq)==0) then
        call integrate_nod(tr_arr(:,:,2), val, mesh)
        if (mype==0) then
           write(*,*) 'total integral of salinity at timestep :', mstep, val
        end if
     end if
  end if
  !6. MOC in density coordinate
  if (ldiag_dMOC)        call diag_densMOC(mode, mesh)

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
subroutine compute_diag_dvd_2ndmoment_burchard_etal_2008(tr_num, mesh)
    use o_arrays
    use o_PARAM, only: tracer_adv
    use g_PARSUP
    
    implicit none
    type(t_mesh), intent(in), target :: mesh
    integer, intent(in)      :: tr_num 
    integer                  :: node, nz
    real(kind=WP)            :: tr_sqr(mesh%nl-1,myDim_nod2D+eDim_nod2D), trAB_sqr(mesh%nl-1,myDim_nod2D+eDim_nod2D)

#include "associate_mesh.h"
  
    !___________________________________________________________________________
    ! square up fields for actual tracers and Adams Bashfort tracer
    ! --> dont forget to square also the halo !!! --> that why do node = ...+eDim_nod2D
    tr_sqr   = 0.0_WP
    trAB_sqr = 0.0_WP
    do node = 1, myDim_nod2D+eDim_nod2D
        do nz = 1, nlevels_nod2D(node)-1
            tr_sqr(nz,node)   = tr_arr(nz,node,tr_num)**2
            trAB_sqr(nz,node) = tr_arr_old(nz,node,tr_num)**2
        end do
    end do
        
    !___________________________________________________________________________
    ! calculate horizintal and vertical advection for squared tracer (2nd moments)
    ! see Burchard and Rennau, 2008, Comparative quantification of physically and 
    ! numerically induced mixing in ocean models ...
    del_ttf_advhoriz = 0.0_WP
    del_ttf_advvert  = 0.0_WP
    select case (tracer_adv)
        case(1) !MUSCL
            ! --> tr_arr_old ... AB interpolated tracer from call init_tracers_AB(tr_num)
            call adv_tracers_muscle_ale(trAB_sqr, .25_WP, 1)
            !                2nd tracer moment is used <--'
            call adv_tracers_vert_ppm_ale(tr_sqr, 1)
            !        2nd tracer moment is used <--'
        case(2) !MUSCL+FCT(3D)
            call adv_tracer_fct_ale(trAB_sqr,tr_sqr, 1.0_WP, 1)
            !                   2nd tracer moment is used <--'
        case default !unknown
            if (mype==0) write(*,*) 'Unknown ALE advection type. Check your namelists.'
            call par_ex(1)
    end select

    !___________________________________________________________________________
    ! add target second moment to DVD
    do node = 1,mydim_nod2D
        do nz = 1,nlevels_nod2D(node)-1
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
            tr_dvd_horiz(nz,node,tr_num) = hnode(nz,node)/hnode_new(nz,node)*trAB_sqr(nz,node) - del_ttf_advhoriz(nz,node)/hnode_new(nz,node)
            tr_dvd_vert(nz,node,tr_num)  = hnode(nz,node)/hnode_new(nz,node)*tr_sqr(  nz,node) - del_ttf_advvert( nz,node)/hnode_new(nz,node)
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
subroutine compute_diag_dvd_2ndmoment_klingbeil_etal_2014(tr_num, mesh)
    use o_arrays
    use o_PARAM, only: tracer_adv
    use g_PARSUP
    
    implicit none
    integer, intent(in)      :: tr_num 
    integer                  :: node, nz
    type(t_mesh), intent(in), target :: mesh

#include "associate_mesh.h"
    !___________________________________________________________________________
    ! calculate horizintal and vertical advection for squared tracer (2nd moments)
    ! see Burchard and Rennau, 2008, Comparative quantification of physically and 
    ! numerically induced mixing in ocean models ...
    del_ttf_advhoriz = 0.0_WP
    del_ttf_advvert  = 0.0_WP
    select case (tracer_adv)
        case(1) !MUSCL
            ! --> tr_arr_old ... AB interpolated tracer from call init_tracers_AB(tr_num)
            call adv_tracers_muscle_ale(tr_arr_old(:,:,tr_num), .25_WP, 2)
            !                                                      |    | 
            !             fraction of fourth-order contribution <--'    |
            !                              2nd tracer moment is used <--'
            call adv_tracers_vert_ppm_ale(tr_arr(:,:,tr_num), 2)
            !                                                 | 
            !                    2nd tracer moment is used <--'
        case(2) !MUSCL+FCT(3D)
            call adv_tracer_fct_ale(tr_arr_old(:,:,tr_num),tr_arr(:,:,tr_num), 1.0_WP, 2)
            !                                                                     |    | 
            !                            fraction of fourth-order contribution <--'    | 
            !                                             2nd tracer moment is used <--'
        case default !unknown
            if (mype==0) write(*,*) 'Unknown ALE advection type. Check your namelists.'
            call par_ex(1)
    end select

    !___________________________________________________________________________
    ! add target second moment to DVD
    do node = 1,mydim_nod2D
        do nz = 1,nlevels_nod2D(node)-1
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
            tr_dvd_horiz(nz,node,tr_num) = hnode(nz,node)/hnode_new(nz,node)*(tr_arr_old(nz,node,tr_num)**2) &
                                           - del_ttf_advhoriz(nz,node)/hnode_new(nz,node)
            tr_dvd_vert(nz,node,tr_num)  = hnode(nz,node)/hnode_new(nz,node)*(tr_arr(    nz,node,tr_num)**2) &
                                           - del_ttf_advvert( nz,node)/hnode_new(nz,node)
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
subroutine compute_diag_dvd(tr_num, mesh)
    use g_config, only: dt
    use o_arrays
    use o_PARAM, only: tracer_adv
    use g_PARSUP
    
    implicit none
    integer, intent(in)      :: tr_num 
    integer                  :: node, nz
    type(t_mesh), intent(in), target :: mesh

#include "associate_mesh.h"
    !___________________________________________________________________________
    ! add discret second moment to DVD
    do node = 1,mydim_nod2D
        do nz = 1,nlevels_nod2D(node)-1
            ! eq 16 & 17  and eq 23. Klingbeil et al. 2014
            !
            ! (phi^2)^(n+1) = 1/V^(n+1)*[ V^(n)*(phi^2)^(n) + dt*ADV[phi^2]  ]
            !
            !  DVD = 1/dt * [ (phi^2)^(n+1) -  ( phi^(n+1) )^2 ]
            !                                         |
            !                                         v
            !                                    now add this part
            ! --> tr_dvd_horiz contains already the expected target second moments
            ! from subroutine compute_diag_dvd_2ndmoment
            tr_dvd_horiz(nz,node,tr_num) = (tr_dvd_horiz(nz,node,tr_num)                                    & 
                                            -( hnode(nz,node)/hnode_new(nz,node)*tr_arr_old(nz,node,tr_num) &
                                              -del_ttf_advhoriz(nz,node)/hnode_new(nz,node)                 &
                                              )**2                                                          &
                                            )/dt
            tr_dvd_vert(nz,node,tr_num)  = (tr_dvd_vert(nz,node,tr_num)                                     &
                                            -( hnode(nz,node)/hnode_new(nz,node)*tr_arr(    nz,node,tr_num) &
                                              -del_ttf_advvert( nz,node)/hnode_new(nz,node)                 &
                                              )**2                                                          &
                                            )/dt
        end do
    end do
end subroutine compute_diag_dvd

end module diagnostics
