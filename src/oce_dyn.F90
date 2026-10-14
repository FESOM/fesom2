module oce_dyn_module
    USE MOD_MESH
    USE MOD_PARTIT
    use par_support_module, only: par_abort
    USE MOD_DYN
    USE o_PARAM
    USE g_CONFIG
    USE g_comm_auto
    USE g_backscatter
    USE o_ARRAYS
    USE MOD_TRACER

    implicit none

    private
    public :: update_vel, compute_vel_nodes, viscosity_filter, &
              visc_filt_bcksct, visc_filt_bilapl, visc_filt_bidiff, &
              compute_ke_wrho, compute_apegen, compute_PePm, &
              check_viscopt, check_validviscopt_5
    public :: mom_on_device
    ! OpenACC: the momentum rhs of this step is computed on the device (set in
    ! compute_vel_rhs; .false. for configurations the device path does not cover)
    logical :: mom_on_device = .false.

contains

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
! 
! Contains routines needed for computations of dynamics.
! includes: update_vel, compute_vel_nodes
!_______________________________________________________________________________
SUBROUTINE update_vel(dynamics, partit, mesh)
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
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    eta_n  => dynamics%eta_n(:)
    d_eta  => dynamics%d_eta(:)
    
    !___________________________________________________________________________    
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax, eta, Fx, Fy, usum, udiff)
    DO elem=1, partit%myDim_elem2D
        elnodes=mesh%elem2D_nodes(:,elem)
        eta=-g*theta*dt*d_eta(elnodes)
        Fx=sum(mesh%gradient_sca(1:3,elem)*eta)
        Fy=sum(mesh%gradient_sca(4:6,elem)*eta)
        nzmin = mesh%ulevels(elem)
        nzmax = mesh%nlevels(elem)

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
    UV=>dynamics%uv(:,:,:)
    UVnode=>dynamics%uvnode(:,:,:)
    !___________________________________________________________________________
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz, k, elem, nln, uln, nle, ule, tx, ty, tvol)
    DO n=1, partit%myDim_nod2D 
        uln = mesh%ulevels_nod2D(n)
        nln = mesh%nlevels_nod2D(n) 
        DO nz=uln, nln-1
            tvol=0.0_WP
            tx  =0.0_WP
            ty  =0.0_WP
            DO k=1, mesh%nod_in_elem2D_num(n)
                elem=mesh%nod_in_elem2D(k,n)
                ule = mesh%ulevels(elem)
                nle = mesh%nlevels(elem)
                !!PS if (nlevels(elem)-1<nz) cycle
                if (nle-1<nz .or. nz<ule) cycle
                tvol=tvol+mesh%elem_area(elem)
                tx=tx+UV(1,nz,elem)*mesh%elem_area(elem)
                ty=ty+UV(2,nz,elem)*mesh%elem_area(elem)
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
    ! on the device path UV_rhs comes from compute_vel_rhs on the device; only option 5
    ! has a device version
    if (mom_on_device .and. option/=5) then
        !$ACC UPDATE HOST(dynamics%uv_rhs)
    end if
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
            call par_abort(partit%MPI_COMM_FESOM, partit%mype)
            stop
    END SELECT
    ! UV_rhs stays on the device for the implicit vertical viscosity
    if (mom_on_device .and. option/=5) then
        !$ACC UPDATE DEVICE(dynamics%uv_rhs)
    end if
end subroutine viscosity_filter  
!
!
!_______________________________________________________________________________
SUBROUTINE visc_filt_bcksct(dynamics, partit, mesh)
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, len, vi
    integer       :: nz, ed, el(2), nelem(3),k, elem, nzmin, nzmax, eledges(3), q
    ! still to be understood but if you allocate these arrays statically the results will be different:
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c, U_b, V_b
    UV     => dynamics%uv(    :,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    U_b    => dynamics%work%u_b(:,:)
    V_b    => dynamics%work%v_b(:,:)
#ifdef ENABLE_OPENACC
    if (mom_on_device) then
        call visc_filt_bcksct_acc(dynamics, partit, mesh)
        return
    end if
#endif

    !___________________________________________________________________________
    ! An analog of harmonic viscosity operator.
    ! Same as visc_filt_h, but with the backscatter. 
    ! Here the contribution from squared velocities is added to the viscosity.    
    ! The contribution from boundary edges is neglected (free slip).
!$OMP PARALLEL DO
    DO elem=1, partit%myDim_elem2D+partit%eDim_elem2D
       U_b(:, elem) = 0.0_WP
       V_b(:, elem) = 0.0_WP
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, nz, ed, el, nelem, k, elem, nzmin, nzmax, update_u, update_v, eledges, q)
    ! Each owned element gathers the contributions of its three edges in ascending
    ! edge order, with the sign of the edge orientation. No thread writes an element
    ! it does not own, and the sum does not depend on the number of threads.
!$OMP DO
    DO elem=1, partit%myDim_elem2D
        eledges=mesh%elem_edges(:,elem)
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        if (eledges(2)>eledges(3)) then; q=eledges(2); eledges(2)=eledges(3); eledges(3)=q; end if
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        DO q=1,3
            ed=eledges(q)
            if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
            el=mesh%edge_tri(:,ed)
            len=sqrt(sum(mesh%elem_area(el)))
            nzmax = minval(mesh%nlevels(el))
            nzmin = maxval(mesh%ulevels(el))
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
            if (el(1)==elem) then
                U_b(nzmin:nzmax-1, elem)=U_b(nzmin:nzmax-1, elem)-update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                V_b(nzmin:nzmax-1, elem)=V_b(nzmin:nzmax-1, elem)-update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            else
                U_b(nzmin:nzmax-1, elem)=U_b(nzmin:nzmax-1, elem)+update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                V_b(nzmin:nzmax-1, elem)=V_b(nzmin:nzmax-1, elem)+update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            end if
        END DO ! --> DO q=1,3
    END DO ! --> DO elem=1, myDim_elem2D
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
    DO ed=1, partit%myDim_nod2D 
        nzmin = mesh%ulevels_nod2D(ed)
        nzmax = mesh%nlevels_nod2D(ed)
        DO nz=nzmin, nzmax-1
            vi=0.0_WP
            u1=0.0_WP
            v1=0.0_WP
            DO k=1, mesh%nod_in_elem2D_num(ed)
                elem=mesh%nod_in_elem2D(k,ed)
                vi=vi+mesh%elem_area(elem)
                u1=u1+U_b(nz,elem)*mesh%elem_area(elem)
                v1=v1+V_b(nz,elem)*mesh%elem_area(elem)
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
    do ed=1, partit%myDim_elem2D
        nelem=mesh%elem2D_nodes(:,ed)
        nzmin = mesh%ulevels(ed)
        nzmax = mesh%nlevels(ed)
        !_______________________________________________________________________
        if (dynamics%use_ssh_se_subcycl) then
            !SD Approximate update for transports. We do not care about accuracy 
            !SD here. --> of course, helem will be better.
            Do nz=nzmin, nzmax-1
                !PS UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+(U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP)*(zbar(nz)-zbar(nz+1)) 
                !PS UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+(V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP)*(zbar(nz)-zbar(nz+1))
                UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+(U_b(nz,ed) -dynamics%visc_easybsreturn*sum(U_c(nz,nelem))/3.0_WP)*mesh%helem(nz, ed)
                UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+(V_b(nz,ed) -dynamics%visc_easybsreturn*sum(V_c(nz,nelem))/3.0_WP)*mesh%helem(nz, ed)
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
#ifdef ENABLE_OPENACC
!
!
!_______________________________________________________________________________
! visc_filt_bcksct on the device (UV and UV_rhs present there), same arithmetic per element and level;
! update_u/update_v are replaced by direct accumulation (no private arrays of run-time size).
SUBROUTINE visc_filt_bcksct_acc(dynamics, partit, mesh)
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, len, vi, gamma0, gamma1, gamma2, easybsreturn
    integer       :: nz, nl1, ed, el1, el2, k, elem, nzmin, nzmax, e1, e2, e3, q
    !___________________________________________________________________________
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c, U_b, V_b
    UV     => dynamics%uv(    :,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    U_b    => dynamics%work%u_b(:,:)
    V_b    => dynamics%work%v_b(:,:)
    gamma0       = dynamics%visc_gamma0
    gamma1       = dynamics%visc_gamma1
    gamma2       = dynamics%visc_gamma2
    easybsreturn = dynamics%visc_easybsreturn
    nl1          = size(U_b, 1)

!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(32)
    DO elem=1, partit%myDim_elem2D+partit%eDim_elem2D
!$ACC LOOP VECTOR
        DO nz=1, nl1
            U_b(nz, elem) = 0.0_WP
            V_b(nz, elem) = 0.0_WP
            U_c(nz, elem) = 0.0_WP
            V_c(nz, elem) = 0.0_WP
        END DO
    END DO
!$ACC END PARALLEL LOOP

    !___________________________________________________________________________
    ! each owned element gathers the contributions of its three edges in ascending
    ! edge order, as on the host
    ! one gang per element, the levels over the vector lanes; each level adds its three edges in edge order
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE(e1, e2, e3, q) VECTOR_LENGTH(32)
    DO elem=1, partit%myDim_elem2D
        e1=mesh%elem_edges(1,elem)
        e2=mesh%elem_edges(2,elem)
        e3=mesh%elem_edges(3,elem)
        if (e1>e2) then; q=e1; e1=e2; e2=q; end if
        if (e2>e3) then; q=e2; e2=e3; e3=q; end if
        if (e1>e2) then; q=e1; e1=e2; e2=q; end if
!$ACC LOOP VECTOR PRIVATE(q, ed, el1, el2, len, nzmin, nzmax, u1, v1, vi)
        DO nz=1, nl1
!$ACC LOOP SEQ
            DO q=1,3
                if (q==1) ed=e1
                if (q==2) ed=e2
                if (q==3) ed=e3
                if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
                el1=mesh%edge_tri(1,ed)
                el2=mesh%edge_tri(2,ed)
                nzmax = min(mesh%nlevels(el1), mesh%nlevels(el2))
                nzmin = max(mesh%ulevels(el1), mesh%ulevels(el2))
                if (nz<nzmin .or. nz>nzmax-1) cycle
                len=sqrt(mesh%elem_area(el1)+mesh%elem_area(el2))
                u1=UV(1,nz,el1)-UV(1,nz,el2)
                v1=UV(2,nz,el1)-UV(2,nz,el2)
                vi=dt*max(gamma0, max(gamma1*sqrt(u1*u1+v1*v1), gamma2*(u1*u1+v1*v1)))*len
                if (el1==elem) then
                    U_b(nz, elem)=U_b(nz, elem)-(u1*vi)/mesh%elem_area(elem)
                    V_b(nz, elem)=V_b(nz, elem)-(v1*vi)/mesh%elem_area(elem)
                else
                    U_b(nz, elem)=U_b(nz, elem)+(u1*vi)/mesh%elem_area(elem)
                    V_b(nz, elem)=V_b(nz, elem)+(v1*vi)/mesh%elem_area(elem)
                end if
            END DO
        END DO
    END DO
!$ACC END PARALLEL LOOP
    call exchange_elem(U_b, partit, luse_g2g=.true.)
    call exchange_elem(V_b, partit, luse_g2g=.true.)

    !___________________________________________________________________________
    ! smooth from elements to nodes
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE(nzmin, nzmax) VECTOR_LENGTH(32)
    DO ed=1, partit%myDim_nod2D
        nzmin = mesh%ulevels_nod2D(ed)
        nzmax = mesh%nlevels_nod2D(ed)
!$ACC LOOP VECTOR PRIVATE(k, elem, vi, u1, v1)
        DO nz=nzmin, nzmax-1
            vi=0.0_WP
            u1=0.0_WP
            v1=0.0_WP
!$ACC LOOP SEQ
            DO k=1, mesh%nod_in_elem2D_num(ed)
                elem=mesh%nod_in_elem2D(k,ed)
                vi=vi+mesh%elem_area(elem)
                u1=u1+U_b(nz,elem)*mesh%elem_area(elem)
                v1=v1+V_b(nz,elem)*mesh%elem_area(elem)
            END DO
            U_c(nz,ed)=u1/vi
            V_c(nz,ed)=v1/vi
        END DO
    END DO
!$ACC END PARALLEL LOOP
    call exchange_nod(U_c, V_c, partit, luse_g2g=.true.)

    !___________________________________________________________________________
    ! back to elements: unsmoothed term minus the backscatter part of the smoothed one
!$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) PRIVATE(nzmin, nzmax, e1, e2, e3) VECTOR_LENGTH(32)
    DO ed=1, partit%myDim_elem2D
        e1 = mesh%elem2D_nodes(1,ed)
        e2 = mesh%elem2D_nodes(2,ed)
        e3 = mesh%elem2D_nodes(3,ed)
        nzmin = mesh%ulevels(ed)
        nzmax = mesh%nlevels(ed)
!$ACC LOOP VECTOR
        DO nz=nzmin, nzmax-1
            UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+U_b(nz,ed) -easybsreturn*(U_c(nz,e1)+U_c(nz,e2)+U_c(nz,e3))/3.0_WP
            UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+V_b(nz,ed) -easybsreturn*(V_c(nz,e1)+V_c(nz,e2)+V_c(nz,e3))/3.0_WP
        END DO
    END DO
!$ACC END PARALLEL LOOP
end subroutine visc_filt_bcksct_acc
#endif
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
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, vi, len
    integer       :: ed, el(2), elem, nz, nzmin, nzmax, eledges(3), q
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    
    !___________________________________________________________________________
!$OMP PARALLEL DO
    DO elem=1, partit%myDim_elem2D+partit%eDim_elem2D
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! Sum up velocity differences over edge with respect to elemtnal index
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, ed, el, nz, nzmin, nzmax, update_u, update_v, elem, eledges, q)
    ! Each owned element gathers the contributions of its three edges in ascending
    ! edge order, with the sign of the edge orientation. No thread writes an element
    ! it does not own, and the sum does not depend on the number of threads.
!$OMP DO
    DO elem=1, partit%myDim_elem2D
        eledges=mesh%elem_edges(:,elem)
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        if (eledges(2)>eledges(3)) then; q=eledges(2); eledges(2)=eledges(3); eledges(3)=q; end if
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        DO q=1,3
            ed=eledges(q)
            if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
            el=mesh%edge_tri(:,ed)
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
            DO  nz=nzmin, nzmax-1
                update_u(nz)=(UV(1,nz,el(1))-UV(1,nz,el(2)))
                update_v(nz)=(UV(2,nz,el(1))-UV(2,nz,el(2)))
            END DO
            if (el(1)==elem) then
                U_c(nzmin:nzmax-1, elem)=U_c(nzmin:nzmax-1, elem)-update_u(nzmin:nzmax-1)
                V_c(nzmin:nzmax-1, elem)=V_c(nzmin:nzmax-1, elem)-update_v(nzmin:nzmax-1)
            else
                U_c(nzmin:nzmax-1, elem)=U_c(nzmin:nzmax-1, elem)+update_u(nzmin:nzmax-1)
                V_c(nzmin:nzmax-1, elem)=V_c(nzmin:nzmax-1, elem)+update_v(nzmin:nzmax-1)
            end if
        END DO ! --> DO q=1,3
    END DO ! --> DO elem=1, myDim_elem2D
!$OMP END DO

    !___________________________________________________________________________
    ! compute viscosity on element
!$OMP DO
    DO ed=1,partit%myDim_elem2D
        len   = sqrt(mesh%elem_area(ed))
        nzmin = mesh%ulevels(ed)
        nzmax = mesh%nlevels(ed)
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
    ! Each owned element gathers the contributions of its three edges in ascending
    ! edge order, with the sign of the edge orientation. No thread writes an element
    ! it does not own, and the sum does not depend on the number of threads.
!$OMP DO
    DO elem=1, partit%myDim_elem2D
        eledges=mesh%elem_edges(:,elem)
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        if (eledges(2)>eledges(3)) then; q=eledges(2); eledges(2)=eledges(3); eledges(3)=q; end if
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        DO q=1,3
            ed=eledges(q)
            if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
            el=mesh%edge_tri(:,ed)
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
            !_______________________________________________________________________
            if (dynamics%use_ssh_se_subcycl) then
                !SD Approximate update for transports. We do not care about accuracy 
                !SD here. --> of course, helem will be better.
                do  nz=nzmin,nzmax-1
                    !PS update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                    !PS update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                    update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2))) * (mesh%helem(nz, el(1))+mesh%helem(nz, el(2)))*0.5_WP
                    update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2))) * (mesh%helem(nz, el(1))+mesh%helem(nz, el(2)))*0.5_WP
                end do
            else
                do  nz=nzmin,nzmax-1
                    update_u(nz)=(U_c(nz,el(1))-U_c(nz,el(2)))
                    update_v(nz)=(V_c(nz,el(1))-V_c(nz,el(2)))
                end do
            end if
        
            !_______________________________________________________________________
            if (el(1)==elem) then
                UV_rhs(1, nzmin:nzmax-1, elem)=UV_rhs(1, nzmin:nzmax-1, elem)-update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                UV_rhs(2, nzmin:nzmax-1, elem)=UV_rhs(2, nzmin:nzmax-1, elem)-update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            else
                UV_rhs(1, nzmin:nzmax-1, elem)=UV_rhs(1, nzmin:nzmax-1, elem)+update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                UV_rhs(2, nzmin:nzmax-1, elem)=UV_rhs(2, nzmin:nzmax-1, elem)+update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            end if
        END DO ! --> DO q=1,3
    END DO ! --> DO elem=1, myDim_elem2D
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
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=8)  :: u1, v1, len, vi, viLapl
    integer       :: ed, el(2), nz, nzmin, nzmax, elem, eledges(3), q
    real(kind=8)  :: update_u(mesh%nl-1), update_v(mesh%nl-1)
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
    real(kind=WP), dimension(:,:)  , pointer :: U_c, V_c
    UV     => dynamics%uv(:,:,:)
    UV_rhs => dynamics%uv_rhs(:,:,:)
    U_c    => dynamics%work%u_c(:,:)
    V_c    => dynamics%work%v_c(:,:)
    
    !___________________________________________________________________________
!$OMP PARALLEL DO
    DO elem=1, partit%myDim_elem2D+partit%eDim_elem2D
       U_c(:, elem) = 0.0_WP
       V_c(:, elem) = 0.0_WP
    END DO
!$OMP END PARALLEL DO

    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u1, v1, len, vi, viLapl, ed, el, nz, nzmin, nzmax, update_u, update_v, elem, eledges, q)
    ! Each owned element gathers the contributions of its three edges in ascending
    ! edge order, with the sign of the edge orientation. No thread writes an element
    ! it does not own, and the sum does not depend on the number of threads.
!$OMP DO
    DO elem=1, partit%myDim_elem2D
        eledges=mesh%elem_edges(:,elem)
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        if (eledges(2)>eledges(3)) then; q=eledges(2); eledges(2)=eledges(3); eledges(3)=q; end if
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        DO q=1,3
            ed=eledges(q)
            if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
            el=mesh%edge_tri(:,ed)
            len=sqrt(sum(mesh%elem_area(el)))
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
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
            if (el(1)==elem) then
                U_c(nzmin:nzmax-1, elem)=U_c(nzmin:nzmax-1, elem)-update_u(nzmin:nzmax-1)
                V_c(nzmin:nzmax-1, elem)=V_c(nzmin:nzmax-1, elem)-update_v(nzmin:nzmax-1)
            else
                U_c(nzmin:nzmax-1, elem)=U_c(nzmin:nzmax-1, elem)+update_u(nzmin:nzmax-1)
                V_c(nzmin:nzmax-1, elem)=V_c(nzmin:nzmax-1, elem)+update_v(nzmin:nzmax-1)
            end if
        END DO ! --> DO q=1,3
    END DO ! --> DO elem=1, myDim_elem2D
!$OMP END DO

    !___________________________________________________________________________
!$OMP MASTER
    call exchange_elem(U_c, partit)
    call exchange_elem(V_c, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
    ! Each owned element gathers the contributions of its three edges in ascending
    ! edge order, with the sign of the edge orientation. No thread writes an element
    ! it does not own, and the sum does not depend on the number of threads.
!$OMP DO
    DO elem=1, partit%myDim_elem2D
        eledges=mesh%elem_edges(:,elem)
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        if (eledges(2)>eledges(3)) then; q=eledges(2); eledges(2)=eledges(3); eledges(3)=q; end if
        if (eledges(1)>eledges(2)) then; q=eledges(1); eledges(1)=eledges(2); eledges(2)=q; end if
        DO q=1,3
            ed=eledges(q)
            if(partit%myList_edge2D(ed)>mesh%edge2D_in) cycle
            el=mesh%edge_tri(:,ed)
            len=sqrt(sum(mesh%elem_area(el)))
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
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
                    ! Optional Laplacian viscosity (active when visc_gamma0_h or visc_gamma1_h > 0)
                    viLapl=dt*max(dynamics%visc_gamma0_h, dynamics%visc_gamma1_h*sqrt(u1*u1+v1*v1))*len
                    !PS update_u(nz)=vi*(U_c(nz,el(1))-U_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                    !PS update_v(nz)=vi*(V_c(nz,el(1))-V_c(nz,el(2)))*(zbar(nz)-zbar(nz+1)) 
                    update_u(nz)=(vi*(U_c(nz,el(1))-U_c(nz,el(2)))+viLapl*u1)*(mesh%helem(nz, el(1))+mesh%helem(nz, el(2)))*0.5_WP
                    update_v(nz)=(vi*(V_c(nz,el(1))-V_c(nz,el(2)))+viLapl*v1)*(mesh%helem(nz, el(1))+mesh%helem(nz, el(2)))*0.5_WP
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
                    ! Optional Laplacian viscosity (active when visc_gamma0_h or visc_gamma1_h > 0)
                    viLapl=dt*max(dynamics%visc_gamma0_h, dynamics%visc_gamma1_h*sqrt(u1*u1+v1*v1))*len
                    update_u(nz)=vi*(U_c(nz,el(1))-U_c(nz,el(2)))+viLapl*u1
                    update_v(nz)=vi*(V_c(nz,el(1))-V_c(nz,el(2)))+viLapl*v1
                end do
            end if 
        
            !_______________________________________________________________________
            if (el(1)==elem) then
                UV_rhs(1, nzmin:nzmax-1, elem)=UV_rhs(1, nzmin:nzmax-1, elem)-update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                UV_rhs(2, nzmin:nzmax-1, elem)=UV_rhs(2, nzmin:nzmax-1, elem)-update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            else
                UV_rhs(1, nzmin:nzmax-1, elem)=UV_rhs(1, nzmin:nzmax-1, elem)+update_u(nzmin:nzmax-1)/mesh%elem_area(elem)
                UV_rhs(2, nzmin:nzmax-1, elem)=UV_rhs(2, nzmin:nzmax-1, elem)+update_v(nzmin:nzmax-1)/mesh%elem_area(elem)
            end if
        END DO ! --> DO q=1,3
    END DO ! --> DO elem=1, myDim_elem2D
!$OMP END DO
!$OMP END PARALLEL
end subroutine visc_filt_bidiff
!
!
!_______________________________________________________________________________
SUBROUTINE compute_ke_wrho(dynamics, partit, mesh)
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
  ! we use pressure to compute this term as it appeares much easier (P*dW) instead of (dP*w)
  dynamics%ke_wrho=0.0_WP
  DO n=1, partit%myDim_nod2D
     nzmin=mesh%ulevels_nod2D(n)
     nzmax=mesh%nlevels_nod2D(n)
     do nz=nzmin, nzmax-1
        wu=0.5*(dynamics%w(nz,   n)+dynamics%w_old(nz,   n))
        wl=0.5*(dynamics%w(nz+1, n)+dynamics%w_old(nz+1, n))
        P =(g*dynamics%eta_n(n)+hpressure(nz, n)/density_0)
        dW=(wu*mesh%area(nz, n)-wl*mesh%area(nz+1, n))
        dW=dW/mesh%area(nz, n)/mesh%hnode_new(nz,n)
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
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target   :: dynamics
    type(t_tracer), intent(in)   , target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh)  , intent(in)   , target   :: mesh
    real(kind=WP), dimension(:,:), pointer  :: salt
    !___________________________________________________________________________
    integer        :: n, k, nzmin, nzmax
    real(kind=WP)  :: JS, GS, D
    !___________________________________________________________________________

salt   => tracers%data(2)%values(:,:)

  DO n=1, partit%myDim_nod2D
     nzmin=mesh%ulevels_nod2D(n)
     nzmax=mesh%nlevels_nod2D(n)
     ! heat and salinity fluxes at node=n
     JS= heat_flux_in(n) / vcpw
     GS=(relax_salt(n) + water_flux(n) * salt(nzmin,n))
     D=density_m_rho0(nzmin,n)
     dynamics%ke_J   (n)=JS
     dynamics%ke_G   (n)=GS
     dynamics%ke_D   (n)=D
     dynamics%ke_D2  (n)=D**2
     !this we need in 3D
     !dynamics%ke_n0  (n)=-bvfreq (nzmin,n)*density_0/g
     DO k=nzmin, nzmax-1
        dynamics%ke_n0(k, n)=-bvfreq (k, n)*density_0/g
     END DO
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
! compute energy conversion (Pm<->Pe), as well as Pm & Pe
SUBROUTINE compute_PePm(dynamics, tracers, partit, mesh)
    IMPLICIT NONE
    type(t_dyn)   , intent(inout), target   :: dynamics
    type(t_tracer), intent(in)   , target   :: tracers
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh)  , intent(in)   , target   :: mesh
   !___________________________________________________________________________
    real(kind=WP), dimension(:,:,:), pointer :: UV
    integer        :: elem, k, nzmin, nzmax, elnodes(3)
    real(kind=WP)  :: dens
    !___________________________________________________________________________
  UV     => dynamics%uv(:,:,:)
  DO elem=1, partit%myDim_elem2D
     elnodes=mesh%elem2D_nodes(:,elem)
     nzmin = mesh%ulevels(elem)
     nzmax = mesh%nlevels(elem)
     DO k=nzmin, nzmax-1
        dynamics%ke_Dx(k,elem)=sum(mesh%gradient_sca(1:3,elem)*density_m_rho0(k,elnodes))
        dynamics%ke_Dy(k,elem)=sum(mesh%gradient_sca(4:6,elem)*density_m_rho0(k,elnodes))
        dens=sum(density_m_rho0(k,elnodes))/3.0_WP
        dynamics%ke_DU(k,elem)=dens*UV(1,k,elem)
        dynamics%ke_DV(k,elem)=dens*UV(2,k,elem)
        dynamics%ke_elemD (k,elem)=dens
        dynamics%ke_elemD2(k,elem)=dens**2
     END DO
  END DO
  call exchange_elem(dynamics%ke_Dx,   partit)
  call exchange_elem(dynamics%ke_Dy,   partit)
  call exchange_elem(dynamics%ke_DU,   partit)
  call exchange_elem(dynamics%ke_DV,   partit)
  call exchange_elem(dynamics%ke_elemD,  partit)
  call exchange_elem(dynamics%ke_elemD2, partit)
END SUBROUTINE compute_PePm
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
    
    !___________________________________________________________________________
    loc_R = 0.0_WP
    loc_A = 0.0_WP
    
#if !defined(__openmp_reproducible)
! A "+" reduction combines the per-thread partial sums in a thread-count
! dependent order, so it is not bit-reproducible: run the loop serially.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(node, nz, nzmax, nzmin, c1, rossbyr_1barocl, &
!$OMP                                  fac_ResR1barocl) REDUCTION(+:loc_R, loc_A)
!$OMP DO
#endif
    do node=1, partit%myDim_nod2D
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
            c1=c1+mesh%hnode_new(nz,node)*(                                         &
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
        fac_ResR1barocl =min(mesh%mesh_resolution(node)/rossbyr_1barocl, 5._WP)
        
        !_______________________________________________________________________
        ! compute local mean ratio but exclude equator |lat|<30°, since Rossby radius 
        ! at equator becomes very large
        loc_R   = loc_R + fac_ResR1barocl
        loc_A   = loc_A + 1.0_WP
    end do
#if !defined(__openmp_reproducible)
! A "+" reduction combines the per-thread partial sums in a thread-count
! dependent order, so it is not bit-reproducible: run the loop serially.
!$OMP END DO
!$OMP END PARALLEL
#endif
!$OMP BARRIER

    !___________________________________________________________________________
    ! compute global mean ratio --> core2 Ratio=4.26 (eddy parameterizted), 
    ! dart Ratio=0.97 (eddy resolving/permitting)
    call MPI_AllREDUCE(loc_R, glb_R, 1, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
    call MPI_AllREDUCE(loc_A, glb_A, 1, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, partit%MPIerr)
    glb_R  = glb_R/glb_A
    
    !___________________________________________________________________________
    ! create warning message when ratio<thresh_ratio
    if (glb_R<thresh_ratio) then 
        if (partit%mype==0) then 
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
        call par_abort(partit%MPI_COMM_FESOM, partit%mype)
        end if
    else
        if (partit%mype==0) then 
        write(*,*) ' --check opt_visc--> Mean Ratio Resol/Rrossby = ',  glb_R 
        end if
    end if     
    
end subroutine check_validviscopt_5

end module oce_dyn_module
