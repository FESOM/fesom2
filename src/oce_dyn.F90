
! A set of routines for computing the horizonlal viscosity
! the control parameters (their default values) are:
! dynamics%visc_gamma0 (0.01 [m/s]), dynamics%visc_gamma1 (0.1 [no dim.]), dynamics%visc_gamma2 (10.[s/m]), Div_c [1.], Leith_c[1.?]
! 1. dynamics%visc_gamma0 has the dimension of velocity. It should be as small as possible, but in any case smaller than 0.01 m/s. 
!    All major ocean circulation models are stable with harmonic viscosity 0.01*len.
! 2. dynamics%visc_gamma1 is nondimensional. In commonly used Leith or Smagorinsky parameterizations it is C/pi^2=0.1 (C is about 1). 
!    We therefore try to follow this, allowing some adjustments (because our mesh is triangular, our resolution is different, etc.). 
!    We however, try to keep dynamics%visc_gamma1<0.1
! 3. dynamics%visc_gamma2 is dimensional (1/velocity). If it is 10, then the respective term dominates starting from |u|=0.1 m/s an so on. It is only used in: 
!    (5) visc_filt_bcksct, (6) visc_filt_bilapl, (7) visc_filt_bidiff
module visc_filt_bcksct_interface
  interface
    subroutine visc_filt_bcksct(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module

module visc_filt_bilapl_interface
  interface
    subroutine visc_filt_bilapl(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module

module visc_filt_bidiff_interface
  interface
    subroutine visc_filt_bidiff(dynamics, partit, mesh)
      use mod_mesh
      USE MOD_PARTIT
      USE MOD_PARSUP
      USE MOD_DYN
      type(t_dyn)   , intent(inout), target :: dynamics
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      
    end subroutine
  end interface
end module

module check_validviscopt_interface
    interface
        subroutine check_validviscopt_5(partit, mesh)
            USE MOD_MESH
            USE MOD_PARTIT
            USE MOD_PARSUP
            type(t_partit), intent(inout), target :: partit
            type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine    
    end interface
end module 

! 
! Contains routines needed for computations of dynamics.
! includes: update_vel, compute_vel_nodes
!_______________________________________________________________________________
SUBROUTINE update_vel(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    use g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer       :: n, elem, elnodes(3), nz, nzmin, nzmax
    real(kind=WP) :: eta(3) 
    real(kind=WP) :: Fx, Fy
    real(kind=WP) :: usum(2), udiff(2)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:)    , pointer :: eta_n, d_eta
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    eta_n  => dynamics%eta_n(:)
    d_eta  => dynamics%d_eta(:)
    
    !___________________________________________________________________________    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax, eta, Fx, Fy, usum, udiff)
    DO elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        eta=-g*theta*dt*d_eta(elnodes)
        Fx=sum(gradient_sca(1:3,elem)*eta)
        Fy=sum(gradient_sca(4:6,elem)*eta)
        nzmin = ulevels(elem)
        nzmax = nlevels(elem)

        if (dynamics%ldiag_ke) then
            DO nz=nzmin, nzmax-1
                dynamics%ke_pre(1,nz,elem)  =dynamics%ke_pre(1,nz,elem)   + Fx
                dynamics%ke_pre(2,nz,elem)  =dynamics%ke_pre(2,nz,elem)   + Fy
                
                
                ! U_(n+1)   - U_n   =  Urhs_n |* (U_(n+1)+U_n)
                ! U_(n+1)^2 - U_n^2 =  Urhs_n * (U_(n+1)+U_n) 
                !                                  | 
                !                                  +-> U_(n+1) = U_n+Urhs_n
                ! U_(n+1)^2 - U_n^2 =  Urhs_n * (2*U_n + Urhs) 
                !                          |            |
                !                          v            v 
                !                        udiff         usum 
                usum(1)  = 2.0_WP*UV(1,nz,elem)+(UV_rhs(1,nz,elem) + Fx)
                usum(2)  = 2.0_WP*UV(2,nz,elem)+(UV_rhs(2,nz,elem) + Fy)
                udiff(1) = UV_rhs(1,nz,elem) + Fx
                udiff(2) = UV_rhs(2,nz,elem) + Fy
                
                ! (U_(n+1)^2 - U_n^2)/2 = usum*udiff/2
                dynamics%ke_du2 (:,nz,elem)       = usum*udiff/2.0_WP
                
                dynamics%ke_pre_xVEL (:,nz,elem)  = usum*dynamics%ke_pre (:,nz,elem)/2.0_WP
                dynamics%ke_adv_xVEL (:,nz,elem)  = usum*dynamics%ke_adv (:,nz,elem)/2.0_WP
                dynamics%ke_cor_xVEL (:,nz,elem)  = usum*dynamics%ke_cor (:,nz,elem)/2.0_WP
                dynamics%ke_hvis_xVEL(:,nz,elem)  = usum*dynamics%ke_hvis(:,nz,elem)/2.0_WP
                dynamics%ke_vvis_xVEL(:,nz,elem)  = usum*dynamics%ke_vvis(:,nz,elem)/2.0_WP
                
                ! U_(n+0.5)    = U_n + 0.5*Urhs
                dynamics%ke_umean(    :,nz,elem)  = usum/2.0_WP
                ! U_(n+0.5)^2
                dynamics%ke_u2mean(   :,nz,elem)  = (usum*usum)/4.0_WP
                
                if (nz==nzmin) then
                    dynamics%ke_wind_xVEL(:,elem) = usum*dynamics%ke_wind(:,elem)/2.0_WP
                end if
                
                if (nz==nzmax-1) then
                    dynamics%ke_drag_xVEL(:,elem) = usum*dynamics%ke_drag(:,elem)/2.0_WP
                end if
            END DO
        end if

        DO nz=nzmin, nzmax-1
            UV(1,nz,elem)= UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
            UV(2,nz,elem)= UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
        END DO
    END DO
!$OMP END PARALLEL DO
!$OMP BARRIER
!!PS Why we do this here eta_n is anyway overwriten through ...
!!PS do node=1, myDim_nod2D+eDim_nod2D
!!PS     if (ulevels_nod2D(node)==1) eta_n(node)=alpha*hbar(node)+(1.0_WP-alpha)*hbar_old(node)
!!PS end do
!!PS !$OMP PARALLEL DO
!!PS     DO n=1, myDim_nod2D+eDim_nod2D
!!PS        eta_n(n)=eta_n(n)+d_eta(n)
!!PS     END DO
!!PS !$OMP END PARALLEL DO
    call exchange_elem(UV, partit)
end subroutine update_vel
!
!
!_______________________________________________________________________________
subroutine compute_vel_nodes(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_DYN
    USE o_PARAM
    use g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer            :: n, nz, k, elem, nln, uln, nle, ule
    real(kind=WP)      :: tx, ty, tvol
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UVnode
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV=>dynamics%uv(:,:,:)
    UVnode=>dynamics%uvnode(:,:,:)
    !___________________________________________________________________________
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, k, elem, nln, uln, nle, ule, tx, ty, tvol)
    DO n=1, myDim_nod2D 
        uln = ulevels_nod2D(n)
        nln = nlevels_nod2D(n) 
        DO nz=uln, nln-1
            tvol=0.0_WP
            tx  =0.0_WP
            ty  =0.0_WP
            DO k=1, nod_in_elem2D_num(n)
                elem=nod_in_elem2D(k,n)
                ule = ulevels(elem)
                nle = nlevels(elem)
                !!PS if (nlevels(elem)-1<nz) cycle
                if (nle-1<nz .or. nz<ule) cycle
                tvol=tvol+elem_area(elem)
                tx=tx+UV(1,nz,elem)*elem_area(elem)
                ty=ty+UV(2,nz,elem)*elem_area(elem)
            END DO
            UVnode(1,nz,n)=tx/tvol
            UVnode(2,nz,n)=ty/tvol
        END DO
    END DO
!$OMP END PARALLEL DO
!$OMP BARRIER
    call exchange_nod(UVnode, partit)
end subroutine compute_vel_nodes
!
!
!_______________________________________________________________________________
subroutine viscosity_filter(option, dynamics, partit, mesh)
    use o_PARAM
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    use visc_filt_bcksct_interface
    use visc_filt_bilapl_interface
    use visc_filt_bidiff_interface
    use g_backscatter
    IMPLICIT NONE 
    integer                               :: option
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    
    ! Driving routine 
    ! Background viscosity is selected in terms of Vl, where V is 
    ! background velocity scale and l is the resolution. V is 0.005 
    ! or 0.01, perhaps it would be better to pass it as a parameter.

    ! h_viscosity_leiht needs vorticity, so vorticity array should be 
    ! allocated. At present, there are two rounds of smoothing in 
    ! h_viscosity. 
    SELECT CASE (option)
        CASE (5)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bcksct'//achar(27)//'[0m'
            call visc_filt_bcksct(dynamics, partit, mesh)
        CASE (6)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bilapl'//achar(27)//'[0m'
            call visc_filt_bilapl(dynamics, partit, mesh)
        CASE (7)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_bidiff'//achar(27)//'[0m'
            call visc_filt_bidiff(dynamics, partit, mesh)
        CASE (8)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call backscatter_coef'//achar(27)//'[0m'
            call backscatter_coef(dynamics, partit, mesh)
            if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[37m'//'         --> call visc_filt_dbcksc'//achar(27)//'[0m'
            call visc_filt_dbcksc(dynamics, partit, mesh) 
        CASE DEFAULT
            if (partit%mype==0) write(*,*) 'mixing scheme with option ' , option, 'has not yet been implemented'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            stop
    END SELECT
end subroutine viscosity_filter  
!
!
!_______________________________________________________________________________
SUBROUTINE visc_filt_bcksct(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, len, vi
    integer       :: nz, ed, el(2), nelem(3),k, elem, nzmin, nzmax
    ! still to be understood but if you allocate these arrays statically the results will be different:
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c, U_b, V_b
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     => dynamics%uv(    :,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    U_b    => dynamics%work%u_b(:,:)
    V_b    => dynamics%work%v_b(:,:)

    !___________________________________________________________________________
    ! An analog of harmonic viscosity operator.
    ! Same as visc_filt_h, but with the backscatter. 
    ! Here the contribution from squared velocities is added to the viscosity.    
    ! The contribution from boundary edges is neglected (free slip).
!$OMP PARALLEL DO
    DO elem=1, myDim_elem2D+eDim_elem2D
       U_b(:, elem) = 0.0_WP
       V_b(:, elem) = 0.0_WP
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, nz, ed, el, nelem, k, elem, nzmin, nzmax, update_u, update_v)
!$OMP DO
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmax = minval(nlevels(el))
        nzmin = maxval(ulevels(el))
        DO  nz=nzmin,nzmax-1
            u1=UV(1,nz,el(1))-UV(1,nz,el(2))
            v1=UV(2,nz,el(1))-UV(2,nz,el(2))
            vi=dt*max(dynamics%visc_gamma0,                         &
                      max(dynamics%visc_gamma1*sqrt(u1*u1+v1*v1),   &
                      dynamics%visc_gamma2*(u1*u1+v1*v1))           &
                    )*len
            update_u(nz)=u1*vi
            update_v(nz)=v1*vi
        END DO 
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(el(1)))
#else
!$OMP ORDERED
#endif
        U_b(nzmin:nzmax-1, el(1))=U_b(nzmin:nzmax-1, el(1))-update_u(nzmin:nzmax-1)/elem_area(el(1))
        V_b(nzmin:nzmax-1, el(1))=V_b(nzmin:nzmax-1, el(1))-update_v(nzmin:nzmax-1)/elem_area(el(1))
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(1)))
        call omp_set_lock  (partit%plock(el(2)))
#endif
        U_b(nzmin:nzmax-1, el(2))=U_b(nzmin:nzmax-1, el(2))+update_u(nzmin:nzmax-1)/elem_area(el(2))
        V_b(nzmin:nzmax-1, el(2))=V_b(nzmin:nzmax-1, el(2))+update_v(nzmin:nzmax-1)/elem_area(el(2))
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(2)))
#else
!$OMP END ORDERED 
#endif
    END DO
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER
    call exchange_elem(U_b, partit)
    call exchange_elem(V_b, partit)
!$OMP END MASTER
!$OMP BARRIER
    
    !___________________________________________________________________________
    ! Compute smoothed viscous term, smoth from elements to nodes
!$OMP DO
    DO ed=1, myDim_nod2D 
        nzmin = ulevels_nod2D(ed)
        nzmax = nlevels_nod2D(ed)
        DO nz=nzmin, nzmax-1
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
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER
    call exchange_nod(U_c, V_c, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! smooth from nodes back to elements, take unsmooth viscos term minus 
    ! smoothed viscos term multiplied backscatter parameter
!$OMP DO
    do ed=1, myDim_elem2D
        nelem=elem2D_nodes(:,ed)
        nzmin = ulevels(ed)
        nzmax = nlevels(ed)
        !_______________________________________________________________________
        if (dynamics%use_ssh_se_subcycl) then
            !SD Approximate update for transports. We do not care about accuracy 
            !SD here. --> of course, helem will be better.
            Do nz=nzmin, nzmax-1
                !PS UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+(U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP)*(zbar(nz)-zbar(nz+1)) 
                !PS UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+(V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP)*(zbar(nz)-zbar(nz+1))
                UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+(U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP)*helem(nz, ed)
                UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+(V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP)*helem(nz, ed)
            END DO
        else
            Do nz=nzmin, nzmax-1
                UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP
                UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP
            END DO
        end if
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine visc_filt_bcksct
!
!
!_______________________________________________________________________________
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity Laplacian, i.e., on an analog of
! the Leith viscosity (Lapl==second derivatives)
! \nu=|3u_c-u_n1-u_n2-u_n3|*sqrt(S_c)/100. There is an additional term
! in viscosity that is proportional to the velocity amplitude squared.
! The coefficient has to be selected experimentally.
SUBROUTINE visc_filt_bilapl(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, vi, len
    integer       :: ed, el(2), elem, nz, nzmin, nzmax
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    
    !___________________________________________________________________________
!$OMP PARALLEL DO
    DO elem=1, myDim_elem2D+eDim_elem2D
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! Sum up velocity differences over edge with respect to elemtnal index
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, ed, el, nz, nzmin, nzmax)
!$OMP DO
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el    = edge_tri(:,ed)
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        DO  nz=nzmin, nzmax-1
            update_u(nz)=(UV(1,nz,el(1))-UV(1,nz,el(2)))
            update_v(nz)=(UV(2,nz,el(1))-UV(2,nz,el(2)))
        END DO
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(el(1)))
#else
!$OMP ORDERED
#endif
        U_c(nzmin:nzmax-1, el(1))=U_c(nzmin:nzmax-1, el(1))-update_u(nzmin:nzmax-1)
        V_c(nzmin:nzmax-1, el(1))=V_c(nzmin:nzmax-1, el(1))-update_v(nzmin:nzmax-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(1)))
        call omp_set_lock  (partit%plock(el(2)))
#endif
        U_c(nzmin:nzmax-1, el(2))=U_c(nzmin:nzmax-1, el(2))+update_u(nzmin:nzmax-1)
        V_c(nzmin:nzmax-1, el(2))=V_c(nzmin:nzmax-1, el(2))+update_v(nzmin:nzmax-1)
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO

    !___________________________________________________________________________
    ! compute viscosity on element
!$OMP DO
    DO ed=1,myDim_elem2D
        len   = sqrt(elem_area(ed))
        nzmin = ulevels(ed)
        nzmax = nlevels(ed)
        Do nz=nzmin, nzmax-1
            ! vi has the sense of harmonic viscosity coef. because of 
            ! division by area in the end 
            u1=U_c(nz,ed)**2+V_c(nz,ed)**2
            vi=max(dynamics%visc_gamma0,                &
                   max(dynamics%visc_gamma1*sqrt(u1),   &
                   dynamics%visc_gamma2*u1)             &
                  )*len*dt
            U_c(nz,ed)=-U_c(nz,ed)*vi                             
            V_c(nz,ed)=-V_c(nz,ed)*vi
        END DO
    END DO
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER
    call exchange_elem(U_c, partit)
    call exchange_elem(V_c, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
!$OMP DO
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !_______________________________________________________________________
        if (dynamics%use_ssh_se_subcycl) then
            !SD Approximate update for transports. We do not care about accuracy 
            !SD here. --> of course, helem will be better.
            do  nz=nzmin,nzmax-1
                !PS update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                !PS update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2))) * (helem(nz, el(1))+helem(nz, el(2)))*0.5_WP
                update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2))) * (helem(nz, el(1))+helem(nz, el(2)))*0.5_WP
            end do
        else
            do  nz=nzmin,nzmax-1
                update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2)))
                update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2)))
            end do
        end if
        
        !_______________________________________________________________________
#if defined(_OPENMP) && ! defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(el(1)))
#else
!$OMP ORDERED
#endif
        UV_rhs(1, nzmin:nzmax-1, el(1))=UV_rhs(1, nzmin:nzmax-1, el(1))-update_u(nzmin:nzmax-1)/elem_area(el(1))
        UV_rhs(2, nzmin:nzmax-1, el(1))=UV_rhs(2, nzmin:nzmax-1, el(1))-update_v(nzmin:nzmax-1)/elem_area(el(1))
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(1)))
        call omp_set_lock  (partit%plock(el(2)))
#endif
        UV_rhs(1, nzmin:nzmax-1, el(2))=UV_rhs(1, nzmin:nzmax-1, el(2))+update_u(nzmin:nzmax-1)/elem_area(el(2))
        UV_rhs(2, nzmin:nzmax-1, el(2))=UV_rhs(2, nzmin:nzmax-1, el(2))+update_v(nzmin:nzmax-1)/elem_area(el(2))
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO
!$OMP END PARALLEL
end subroutine visc_filt_bilapl
!
!
!_______________________________________________________________________________
! Strictly energy dissipative and momentum conserving version
! Viscosity depends on velocity differences, and is introduced symmetrically 
! into both stages of biharmonic operator
! On each edge, \nu=sqrt(|u_c1-u_c2|*sqrt(S_c1+S_c2)/100)
! The effect is \nu^2
! Quadratic in velocity term can be introduced if needed.
SUBROUTINE visc_filt_bidiff(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, len, vi
    integer       :: ed, el(2), nz, nzmin, nzmax, elem
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    
    !___________________________________________________________________________
!$OMP PARALLEL DO
    DO elem=1, myDim_elem2D+eDim_elem2D
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, ed, el, nz, nzmin, nzmax, update_u, update_v)
!$OMP DO
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        DO  nz=nzmin,nzmax-1
            u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
            v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
            vi=u1*u1+v1*v1
            vi=sqrt(max(dynamics%visc_gamma0,           &
                    max(dynamics%visc_gamma1*sqrt(vi),  &
                    dynamics%visc_gamma2*vi)            &
                   )*len)
            ! vi=sqrt(max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(vi), dynamics%visc_gamma2*vi))*len)
            update_u(nz)=u1*vi
            update_v(nz)=v1*vi
        END DO
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(el(1)))
#else
!$OMP ORDERED
#endif
        U_c(nzmin:nzmax-1, el(1))=U_c(nzmin:nzmax-1, el(1))-update_u(nzmin:nzmax-1)
        V_c(nzmin:nzmax-1, el(1))=V_c(nzmin:nzmax-1, el(1))-update_v(nzmin:nzmax-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(1)))
        call omp_set_lock  (partit%plock(el(2)))
#endif
        U_c(nzmin:nzmax-1, el(2))=U_c(nzmin:nzmax-1, el(2))+update_u(nzmin:nzmax-1)
        V_c(nzmin:nzmax-1, el(2))=V_c(nzmin:nzmax-1, el(2))+update_v(nzmin:nzmax-1)
#if defined(_OPENMP) && !defined(__openmp_reproducible) 
        call omp_unset_lock(partit%plock(el(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER
    call exchange_elem(U_c, partit)
    call exchange_elem(V_c, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
!$OMP DO
    DO ed=1, myDim_edge2D+eDim_edge2D
        if(myList_edge2D(ed)>edge2D_in) cycle
        el=edge_tri(:,ed)
        len=sqrt(sum(elem_area(el)))
        nzmin = maxval(ulevels(el))
        nzmax = minval(nlevels(el))
        !_______________________________________________________________________
        if (dynamics%use_ssh_se_subcycl) then
            !SD Approximate update for transports. We do not care about accuracy 
            !SD here. --> of course, helem will be better.
            do nz=nzmin,nzmax-1
                u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
                v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
                vi=u1*u1+v1*v1
                vi=-dt*sqrt(max(dynamics%visc_gamma0,           &
                            max(dynamics%visc_gamma1*sqrt(vi),  &
                            dynamics%visc_gamma2*vi)            &
                        )*len)
                ! vi=-dt*sqrt(max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(vi), dynamics%visc_gamma2*vi))*len)
                !PS update_u(nz)=vi*(U_c(nz,el(1))-U_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                !PS update_v(nz)=vi*(V_c(nz,el(1))-V_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                update_u(nz)=vi*(U_c(nz,el(1))-U_c(nz,el(2)))*(helem(nz, el(1))+helem(nz, el(2)))*0.5_WP
                update_v(nz)=vi*(V_c(nz,el(1))-V_c(nz,el(2)))*(helem(nz, el(1))+helem(nz, el(2)))*0.5_WP
            end do
        else
            do nz=nzmin,nzmax-1
                u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
                v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))
                vi=u1*u1+v1*v1
                vi=-dt*sqrt(max(dynamics%visc_gamma0,           &
                            max(dynamics%visc_gamma1*sqrt(vi),  &
                            dynamics%visc_gamma2*vi)            &
                        )*len)
                ! vi=-dt*sqrt(max(dynamics%visc_gamma0, dynamics%visc_gamma1*max(sqrt(vi), dynamics%visc_gamma2*vi))*len)
                update_u(nz)=vi*(U_c(nz,el(1))-U_c(nz,el(2)))
                update_v(nz)=vi*(V_c(nz,el(1))-V_c(nz,el(2)))
            end do
        end if 
        
        !_______________________________________________________________________
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_set_lock(partit%plock(el(1)))
#else
!$OMP ORDERED
#endif
        UV_rhs(1, nzmin:nzmax-1, el(1))=UV_rhs(1, nzmin:nzmax-1, el(1))-update_u(nzmin:nzmax-1)/elem_area(el(1))
        UV_rhs(2, nzmin:nzmax-1, el(1))=UV_rhs(2, nzmin:nzmax-1, el(1))-update_v(nzmin:nzmax-1)/elem_area(el(1))
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(1)))
        call omp_set_lock  (partit%plock(el(2)))
#endif
        UV_rhs(1, nzmin:nzmax-1, el(2))=UV_rhs(1, nzmin:nzmax-1, el(2))+update_u(nzmin:nzmax-1)/elem_area(el(2))
        UV_rhs(2, nzmin:nzmax-1, el(2))=UV_rhs(2, nzmin:nzmax-1, el(2))+update_v(nzmin:nzmax-1)/elem_area(el(2))
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(el(2)))
#else
!$OMP END ORDERED
#endif
    END DO
!$OMP END DO
!$OMP END PARALLEL
end subroutine visc_filt_bidiff
!
!
!_______________________________________________________________________________
SUBROUTINE compute_ke_wrho(dynamics, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    USE o_ARRAYS
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    real(kind=WP) , pointer               :: inv_rhowat
    !___________________________________________________________________________
    integer        :: n, nz, nzmin, nzmax
    real(kind=WP)  :: wu, wl
    real(kind=WP)  :: dW, P
    !___________________________________________________________________________
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  ! we use pressure to compute this term as it appeares much easier (P*dW) instead of (dP*w)
  dynamics%ke_wrho=0.0_WP
  DO n=1, myDim_nod2D
     nzmin=ulevels_nod2D(n)
     nzmax=nlevels_nod2D(n)
     do nz=nzmin, nzmax-1
        wu=0.5*(dynamics%w(nz,   n)+dynamics%w_old(nz,   n))
        wl=0.5*(dynamics%w(nz+1, n)+dynamics%w_old(nz+1, n))
        P =(g*dynamics%eta_n(n)+hpressure(nz, n)/density_0)
        dW=(wu*area(nz, n)-wl*area(nz+1, n))
        dW=dW/area(nz, n)/hnode_new(nz,n)
!       r=(g*dynamics%eta_n(n)+hpressure(nz, n)/density_0)*(wu*area(nz, n)-wl*area(nz+1, n)+0.5_WP*(hnode_new(nz, n)-mesh%hnode_old(nz, n))/dt*area(nz, n) )
        dynamics%ke_Pfull(nz, n) = P
        dynamics%ke_dW   (nz, n) = dW
        dynamics%ke_wrho (nz, n) = P*dW*dt
     end do
  END DO
  call exchange_nod(dynamics%ke_wrho, partit)
END SUBROUTINE compute_ke_wrho
!
!
!_______________________________________________________________________________
! APE generation stuff
SUBROUTINE compute_apegen(dynamics, tracers, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_TRACER
    USE MOD_PARSUP
    use MOD_DYN
    USE o_PARAM
    USE g_comm_auto
    USE o_ARRAYS
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target   :: dynamics
    type(t_tracer), intent(in)   , target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh)  , intent(in)   , target   :: mesh
    real(kind=WP), dimension(:,:), pointer  :: salt
    !___________________________________________________________________________
    integer        :: n, nzmin
    real(kind=WP)  :: JS, GS, D
    !___________________________________________________________________________
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

salt   => tracers%data(2)%values(:,:)

  DO n=1, myDim_nod2D
     nzmin=ulevels_nod2D(n)
     ! heat and salinity fluxes at node=n
     JS= heat_flux_in(n) / vcpw
     GS=(relax_salt(n) + water_flux(n) * salt(nzmin,n))
     D=density_m_rho0(nzmin,n)
     dynamics%ke_J   (n)=JS
     dynamics%ke_G   (n)=GS
     dynamics%ke_D   (n)=D
     dynamics%ke_D2  (n)=D**2
     dynamics%ke_n0  (n)=-bvfreq (nzmin,n)*density_0/g
     dynamics%ke_JD(n)  =JS*D
     dynamics%ke_GD(n)  =GS*D
     dynamics%ke_D (n)  =D
     dynamics%ke_swA (n)=sw_alpha(nzmin, n)
     dynamics%ke_swB (n)=sw_beta (nzmin, n)
  END DO
  call exchange_nod(dynamics%ke_J, partit)
  call exchange_nod(dynamics%ke_G, partit)
  call exchange_nod(dynamics%ke_D, partit)

  call exchange_nod(dynamics%ke_D2, partit)
  call exchange_nod(dynamics%ke_n0, partit)

  call exchange_nod(dynamics%ke_JD, partit)
  call exchange_nod(dynamics%ke_GD, partit)

  call exchange_nod(dynamics%ke_swA, partit)
  call exchange_nod(dynamics%ke_swB, partit)
END SUBROUTINE compute_apegen


!
!
!_______________________________________________________________________________
! check validity of visc_opt=5 selection on basis of ratio betweem resolution and 
! 1st baroclinic rossby radius. visc_opt=5 ("easy backscatter") in eddy 
! resolving/partly eddy permitting setups can lead to problems in form of 
! exedingly strong near boundary currents that can lead e.g to very weak 
! Drake Passage throughflow (<80Sv). In this case better use visc_opt=7, 
! which is the flow aware viscosity option 
subroutine check_viscopt(dynamics, partit, mesh)
    USE MOD_DYN
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE check_validviscopt_interface
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    
    if (dynamics%check_opt_visc) then
        select case (dynamics%opt_visc)
            case(5) 
                ! check validity of visc_opt=5 especially in higher resolved setups
                ! --> there it can lead to problems
                call check_validviscopt_5(partit, mesh)
        end select   
    end if 
end subroutine check_viscopt

!
!
!_______________________________________________________________________________
! check if viscopt=5 is a valid and recommended option for the used configuration
subroutine check_validviscopt_5(partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE o_PARAM , ONLY: rad
    USE o_ARRAYS, ONLY: bvfreq
    USE g_CONFIG
    USE g_comm_auto
    IMPLICIT NONE
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer                               :: node, nz, nzmax, nzmin
    real(kind=WP)                         :: f_min=1.e-6_WP, z_min=100.0_WP, r_max=200000._WP, r_min=2000.0_WP
    real(kind=WP)                         :: c_min=0.5_WP, c1
    real(kind=WP)                         :: excl_EQ=30.0, thresh_ratio=1.5
    real(kind=WP)                         :: loc_R, loc_A
    real(kind=WP)                         :: glb_R, glb_A
    real(kind=WP)                         :: fac_ResR1barocl, rossbyr_1barocl
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    loc_R = 0.0_WP
    loc_A = 0.0_WP
    
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, nz, nzmax, nzmin, c1, rossbyr_1barocl, &
!$OMP                                  fac_ResR1barocl) REDUCTION(+:loc_R, loc_A)
!$OMP DO
    do node=1, myDim_nod2D
        !_______________________________________________________________________
        ! Exlude the equator |lat|<30° from checking the ration between resolution
        ! and first baroclinic rossby radius  
        if (abs(mesh%geo_coord_nod2D(2,node)/rad)<excl_EQ) cycle
        
        !_______________________________________________________________________
        !ca. first baroclinic gravity wave speed limited from below by c_min
        nzmax=mesh%nlevels_nod2D_min(node)
        nzmin=mesh%ulevels_nod2D_max(node)
        c1=0.0_WP
        do nz=nzmin, nzmax-1
            c1=c1+hnode_new(nz,node)*(                                         &
                                    sqrt(abs(max(bvfreq(nz  ,node), 0._WP)))+  &
                                    sqrt(abs(max(bvfreq(nz+1,node), 0._WP)))   &
                                  )/2._WP ! add abs() for -0 case, cray
        end do
        c1=max(c_min, c1/pi)
        
        !_______________________________________________________________________
        ! compute 1st. baroclinic rossby radius
        rossbyr_1barocl = max( min( c1/max( abs(mesh%coriolis_node(node)), f_min), r_max), r_min)
        
        !_______________________________________________________________________
        ! compute ratio between rossby and resolution, if ratio >=1 setup is not eddy 
        ! resolving in best case eddy permitting --> GM/Redi will still be neccessary
        fac_ResR1barocl =min(mesh_resolution(node)/rossbyr_1barocl, 5._WP)
        
        !_______________________________________________________________________
        ! compute local mean ratio but exclude equator |lat|<30°, since Rossby radius 
        ! at equator becomes very large
        loc_R   = loc_R + fac_ResR1barocl
        loc_A   = loc_A + 1.0_WP
    end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP BARRIER

    !___________________________________________________________________________
    ! compute global mean ratio --> core2 Ratio=4.26 (eddy parameterizted), 
    ! dart Ratio=0.97 (eddy resolving/permitting)
    call MPI_AllREDUCE(loc_R, glb_R, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(loc_A, glb_A, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    glb_R  = glb_R/glb_A
    
    !___________________________________________________________________________
    ! create warning message when ratio<thresh_ratio
    if (glb_R<thresh_ratio) then 
        if (mype==0) then 
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' --check opt_visc--> Mean Ratio Resol/Rrossby = ',  glb_R 
        write(*,*) '____________________________________________________________________'
        write(*,*) ' WARNING: You want to use opt_visc=5 (easy backscatter viscosity) in'
        write(*,*) '          an eddy resolving to eddy permitting (Resol/Rrossby<1.5)  '
        write(*,*) '          configuration. It revealed to us that this can lead to    '
        write(*,*) '          problems in form of unrealistically strong near boundary  '
        write(*,*) '          currents that can ultimatively lead to an e.g very weak   '
        write(*,*) '          Drake Passage throughflow (<80Sv in spinup up dart        '
        write(*,*) '          configuration). For these kind of setups we recommend to  '
        write(*,*) '          use opt_visc=7, which is the flow aware viscosity         '
        write(*,*) '          parameterisation.                                         '
        write(*,*) '          The easy backscatter viscosity was designed to energetise '
        write(*,*) '          coarse resolved configurations!!!                         '
        write(*,*) '          --> So please change the viscosity option in namelist.dyn '
        write(*,*) '              to opt_visc=7 .                                       '
        write(*,*) '          --> If you insist in using opt_visc=5 in your simulation, '
        write(*,*) '              you can switch off this check in namelist.dyn by      '
        write(*,*) '              setting check_opt_visc=.false.                        '
        write(*,*) '                                                                    '
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
        end if
    else
        if (mype==0) then 
        write(*,*) ' --check opt_visc--> Mean Ratio Resol/Rrossby = ',  glb_R 
        end if
    end if     
    
end subroutine check_validviscopt_5
