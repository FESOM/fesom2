module oce_dyn_visc
    ! Horizontal momentum viscosity (the SEPARATE operator run AFTER compute_vel_rhs in
    ! the timestep, FESOM2 oce_ale.F90:3822 — NOT inside compute_vel_rhs like momadv).
    ! Transcribed from FESOM2 v2.7.3 oce_dyn.F90 (viscosity_filter 230-274 +
    ! visc_filt_bidiff 591-744).
    !
    ! v1 supports only opt_visc=7 (`visc_filt_bidiff`) — the flow-aware BIHARMONIC
    ! viscosity used by the reduced-M2 oracle namelist. opt_visc 5/6/8 (backscatter,
    ! plain bi-Laplacian, dynamic backscatter) are later milestones.
    !
    ! visc_filt_bidiff is a biharmonic (Del^4) operator built as TWO edge-based
    ! Laplacian sweeps over the INTERIOR edges only (free-slip on the boundary edges):
    !   Pass 1: per interior edge, the across-edge velocity jump u1=UV(el1)-UV(el2)
    !           times a flow-aware coefficient vi=sqrt(max(g0,max(g1*|du|,g2*|du|^2))*len)
    !           is scattered (-/+) into the element field U_c/V_c -> the first Laplacian.
    !   (exchange_elem(U_c/V_c) — 1-rank no-op.)
    !   Pass 2: per interior edge, the across-edge jump of U_c (the SECOND Laplacian)
    !           times vi2=-dt*sqrt(...) (plus an optional Laplacian term viLapl, zero on
    !           pi) is scattered (-/+)/elem_area into UV_rhs.
    ! The net effect ~ -dt * Del^2(nu * Del^2 u), strictly dissipative + momentum
    ! conserving (the symmetric -/+ edge scatter conserves the area-weighted sum).
    !
    ! Pinned to the pi / reduced-M2 gated path (all v1-out / UNGATED branches dropped,
    ! re-added with their own gates later):
    !   - opt_visc==7                     (5/6/8 are M2.x)
    !   - use_ssh_se_subcycl=.false.      -> the non-subcycl update (no helem weighting)
    !   - visc_gamma0_h==0, visc_gamma1_h==0 (pi) -> viLapl==0, pure biharmonic. The
    !     viLapl term is transcribed (so the code is complete) but contributes 0 on pi;
    !     the harmonic-addition path (gamma_h>0) is a DEFERRED sub-gate (see LESSONS L17).
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern — every operand is already
    ! byte-pinned, so faithful transcription byte-matches like M1.1-M2.4):
    !  - The two per-element edge-scatters are FP-order-sensitive; the 1-rank global
    !    edge order matches FESOM2's exactly (the geometry gate pinned edges/edge_tri
    !    order; the area gate pinned elem_area). The FESOM2 oracle is built
    !    ENABLE_OPENMP=OFF, so its omp locks / !$OMP ORDERED compile out and the scatter
    !    runs serially in edge order 1..edge2D — exactly this serial loop (L16).
    !  - INTERIOR edges only: FESOM2 `if(myList_edge2D(ed)>edge2D_in) cycle`. At 1-rank
    !    myList_edge2D is identity, and the mesh orders interior edges first (1..edge2D_in,
    !    boundary edge2D_in+1..edge2D — fvom_init), so this is `if(ed>edge2D_in) cycle`.
    !    edge2D_in is byte-pinned by the geometry gate (same edgenum.out, same ordering).
    !  - `len=sqrt(sum(elem_area(el)))`, `*sqrt(...)`, `-dt*sqrt(...)` and the final
    !    `/elem_area(el)` are RUNTIME divisors/operands but byte-identical on both sides
    !    (geom-proven), so the -no-prec-div reciprocal matches (L7/L10/L14).
    !  - edge_tri(:,ed) is (2,edge2D) — the `:` is the size-2 left/right slot, NOT a
    !    MAX_NV dim, so no L15 shape trap. Interior edges: el(2) /= el(1).
    !  - 1-rank only: exchange_elem(U_c/V_c) is a no-op (lifted at M2.12). The two passes
    !    are separate serial loops, so pass 1 fully fills U_c before pass 2 reads it
    !    (FESOM2 uses an !$OMP BARRIER for the same effect).
    !  - Local intermediates u1/v1/vi/len/update_* are real(WP). FESOM2 hard-codes them
    !    real(kind=8); at the DP anchor WP==8 so byte-identical. (A future SP gate would
    !    need to replicate FESOM2's kind=8 promotion here — a deferred precision nuance.)
    use mod_precision, only: WP
    use mod_mesh,      only: t_mesh
    use mod_dyn,       only: t_dyn
    use mod_partit,    only: t_partit
    use mod_part_bounds, only: local_dims, is_multirank
    use mod_halo,      only: exchange_elem
    implicit none
    private
    public :: viscosity_filter

contains

    subroutine viscosity_filter(option, dynamics, mesh, dt, partit)
        ! Driving routine — dispatch on the horizontal-viscosity scheme. v1 supports
        ! only opt_visc=7 (the reduced-M2 / pi-gated biharmonic). dt is passed
        ! explicitly (FESOM2 reads it from g_config; FESOM3 has no global dt — cf.
        ! compute_vel_rhs).
        integer,       intent(in)            :: option
        type(t_dyn),   intent(inout), target :: dynamics
        type(t_mesh),  intent(in),    target :: mesh
        real(kind=WP), intent(in)            :: dt
        type(t_partit), intent(in), optional :: partit
        select case (option)
        case (7)
            call visc_filt_bidiff(dynamics, mesh, dt, partit)
        case default
            write(*,*) 'viscosity_filter: opt_visc=', option, &
                       ' not implemented in v1 (only opt_visc=7)'
            error stop 1
        end select
    end subroutine viscosity_filter

    !==========================================================================
    subroutine visc_filt_bidiff(dynamics, mesh, dt, partit)
        ! Strictly energy-dissipative, momentum-conserving biharmonic viscosity.
        ! Transcribed from FESOM2 v2.7.3 oce_dyn.F90:591-744 (the non-subcycl branch).
        ! M2.12c: optional partit -> U_c/V_c zeroed over owned+eDim elements (FESOM2
        ! :622 do elem=1,myDim_elem2D+eDim_elem2D); both edge passes run over owned+halo
        ! edges (FESOM2 :631/:679 do ed=1,myDim_edge2D+eDim_edge2D) with the interior
        ! test on the GLOBAL edge id (myList_edge2D(ed)>edge2D_in); exchange_elem(U_c/V_c)
        ! between the passes so pass 2 reads the halo-element Laplacian. UV is halo-valid
        ! (update_vel exchanges it / prescribed+exchanged); the pass-2 UV_rhs scatter
        ! reaches halo elements exactly as FESOM2, so the post-visc halo UV_rhs matches.
        type(t_dyn),   intent(inout), target :: dynamics
        type(t_mesh),  intent(in),    target :: mesh
        real(kind=WP), intent(in)            :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        real(kind=WP) :: u1, v1, len, vi, viLapl
        integer       :: ed, el(2), nz, nzmin, nzmax, elem
        integer       :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        logical       :: lmr
        real(kind=WP) :: update_u(mesh%nl-1), update_v(mesh%nl-1)
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:,:),   pointer :: U_c, V_c

        UV     => dynamics%uv
        UV_rhs => dynamics%uv_rhs
        U_c    => dynamics%work%u_c
        V_c    => dynamics%work%v_c
        call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
        lmr = is_multirank(partit)

        !______________________________________________________________________
        ! zero the first-stage Laplacian accumulator over owned+eDim elements
        do elem = 1, nElemL
            U_c(:, elem) = 0.0_WP
            V_c(:, elem) = 0.0_WP
        end do

        !______________________________________________________________________
        ! Pass 1: first Laplacian -> U_c/V_c (interior edges only; free slip on bnd).
        do ed = 1, nEdgeL
            if (lmr) then
                if (partit%myList_edge2D(ed) > mesh%edge2D_in) cycle   ! boundary edge
            else
                if (ed > mesh%edge2D_in) cycle        ! boundary edge -> free slip
            end if
            el    = mesh%edge_tri(:, ed)
            len   = sqrt(sum(mesh%elem_area(el)))
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
            do nz = nzmin, nzmax-1
                u1 = UV(1,nz,el(1)) - UV(1,nz,el(2))
                v1 = UV(2,nz,el(1)) - UV(2,nz,el(2))
                vi = u1*u1 + v1*v1
                vi = sqrt(max(dynamics%visc_gamma0,                  &
                          max(dynamics%visc_gamma1*sqrt(vi),         &
                              dynamics%visc_gamma2*vi)               &
                         )*len)
                update_u(nz) = u1*vi
                update_v(nz) = v1*vi
            end do
            ! symmetric edge scatter (OpenMP-off oracle -> serial, L16)
            U_c(nzmin:nzmax-1, el(1)) = U_c(nzmin:nzmax-1, el(1)) - update_u(nzmin:nzmax-1)
            V_c(nzmin:nzmax-1, el(1)) = V_c(nzmin:nzmax-1, el(1)) - update_v(nzmin:nzmax-1)
            U_c(nzmin:nzmax-1, el(2)) = U_c(nzmin:nzmax-1, el(2)) + update_u(nzmin:nzmax-1)
            V_c(nzmin:nzmax-1, el(2)) = V_c(nzmin:nzmax-1, el(2)) + update_v(nzmin:nzmax-1)
        end do

        if (lmr) then
            call exchange_elem(U_c, partit)       ! FESOM2 :672
            call exchange_elem(V_c, partit)       ! FESOM2 :673
        end if

        !______________________________________________________________________
        ! Pass 2: second Laplacian (across-edge jump of U_c) + optional Laplacian
        ! term -> scatter into UV_rhs / elem_area. use_ssh_se_subcycl=.false. branch.
        do ed = 1, nEdgeL
            if (lmr) then
                if (partit%myList_edge2D(ed) > mesh%edge2D_in) cycle   ! boundary edge
            else
                if (ed > mesh%edge2D_in) cycle        ! boundary edge -> free slip
            end if
            el    = mesh%edge_tri(:, ed)
            len   = sqrt(sum(mesh%elem_area(el)))
            nzmin = maxval(mesh%ulevels(el))
            nzmax = minval(mesh%nlevels(el))
            do nz = nzmin, nzmax-1
                u1 = UV(1,nz,el(1)) - UV(1,nz,el(2))
                v1 = UV(2,nz,el(1)) - UV(2,nz,el(2))
                vi = u1*u1 + v1*v1
                vi = -dt*sqrt(max(dynamics%visc_gamma0,             &
                              max(dynamics%visc_gamma1*sqrt(vi),    &
                                  dynamics%visc_gamma2*vi)          &
                             )*len)
                ! optional harmonic (Laplacian) viscosity — zero on pi (gamma_h=0)
                viLapl = dt*max(dynamics%visc_gamma0_h,                       &
                                dynamics%visc_gamma1_h*sqrt(u1*u1+v1*v1))*len
                update_u(nz) = vi*(U_c(nz,el(1))-U_c(nz,el(2))) + viLapl*u1
                update_v(nz) = vi*(V_c(nz,el(1))-V_c(nz,el(2))) + viLapl*v1
            end do
            ! symmetric edge scatter into UV_rhs (OpenMP-off oracle -> serial, L16)
            UV_rhs(1,nzmin:nzmax-1,el(1)) = UV_rhs(1,nzmin:nzmax-1,el(1)) - update_u(nzmin:nzmax-1)/mesh%elem_area(el(1))
            UV_rhs(2,nzmin:nzmax-1,el(1)) = UV_rhs(2,nzmin:nzmax-1,el(1)) - update_v(nzmin:nzmax-1)/mesh%elem_area(el(1))
            UV_rhs(1,nzmin:nzmax-1,el(2)) = UV_rhs(1,nzmin:nzmax-1,el(2)) + update_u(nzmin:nzmax-1)/mesh%elem_area(el(2))
            UV_rhs(2,nzmin:nzmax-1,el(2)) = UV_rhs(2,nzmin:nzmax-1,el(2)) + update_v(nzmin:nzmax-1)/mesh%elem_area(el(2))
        end do
    end subroutine visc_filt_bidiff

end module oce_dyn_visc
