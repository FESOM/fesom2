module oce_tracer_mod
    ! Per-tracer pre-advection setup, transcribed from FESOM2 v2.7.3
    ! oce_tracer_mod.F90:13-145 (SUBROUTINE init_tracers_AB).
    !
    ! Called once per tracer per step, before do_oce_adv_tra. It:
    !   1. zeroes the tendency accumulators del_ttf / del_ttf_advhoriz / del_ttf_advvert
    !   2. Adams-Bashforth-interpolates values -> valuesAB (the field the high-order
    !      advection scheme reconstructs from):
    !        AB2: valuesAB = -(0.5+eps)*valuesold(1) + (1.5+eps)*values
    !        AB3: valuesAB = (5*valuesold(2) - 16*valuesold(1) + 23*values)/12
    !   3. rolls the history valuesold
    !   4. rebuilds the ELEMENTAL gradient of `values` (NOT valuesAB — FESOM2 commented
    !      the AB variant out at oce_tracer_mod.F90:126-127; see LESSONS L11), from which
    !      the MUSCL kernels take the up/downwind reconstruction on the fly (FESOM2
    !      stores a per-edge array edge_up_dn_grad instead).
    !
    ! The AB2 offset `eps` is FESOM2's o_PARAM `epsilon = 0.1` (oce_modules.F90:92),
    ! a runtime module VARIABLE there (not a parameter). It is kept non-parameter here
    ! (mod_config) so the anchor compiler does not fold `(1.5_WP+eps)` to a literal at
    ! compile time, which could differ in the last bit from FESOM2's runtime sum under
    ! the fast-math anchor flags (the L7 literal-vs-runtime trap, with + not /).
    !
    ! SCOPE / clean-architecture deviations (D7, "USE-globals -> explicit args"):
    !  * 1-rank only (myDim_* == global). FESOM2's halo exchanges of tr_xy / tr_z and
    !    the begin/end overlap are no-ops at 1 rank and are dropped; multi-rank is M1.5.
    !  * tr_xy (FESOM2 o_ARRAYS global) is twork%tr_xy here (also reused by the
    !    horizontal diffusion of the same tracer). The redundant SECOND tracer_gradient_elements(values) call
    !    (FESOM2 line 142, "redefine to current timestep") is omitted: `values` is
    !    unchanged across init_tracers_AB, so it recomputes the identical tr_xy.
    !  * tracer_gradient_z -> tr_z (FESOM2 line 131) is omitted: tr_z is the vertical
    !    tracer gradient, consumed by vertical diffusion (M2), NOT by advection. M1.4
    !    is advection-only, so it has no effect on the gated del_ttf.
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_tracer,       only: t_tracer
    use mod_config,       only: ab_epsilon
    use mod_part_bounds,  only: owned_bounds, is_multirank
    use mod_halo,         only: exchange_elem_full
    use oce_tracer_grad,  only: tracer_gradient_elements
    implicit none
    private
    public :: init_tracers_AB

contains

    subroutine init_tracers_AB(tr_num, tracers, mesh, partit)
        ! M2.12b: optional partit. The del_ttf zero / AB interpolation / history roll
        ! run over OWNED+HALO nodes (values/valuesold are halo-consistent so valuesAB is
        ! valid at the halo the HO kernels read). tr_xy is built on OWNED elements then
        ! exchanged over the FULL element halo (FESOM2 init_tracers_AB:128-132): the MUSCL
        ! kernels read it at the up/downwind elements of owned edges (halo elements of a halo node's
        ! element list.
        integer,        intent(in)    :: tr_num
        type(t_mesh),   intent(in)    :: mesh
        type(t_tracer), intent(inout) :: tracers
        type(t_partit), intent(in), optional :: partit
        integer :: n, nz
        integer :: nNodO, nNodL, nEdgeO, nElemO, nElemA

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! del_ttf will accumulate all advection/diffusion contributions this step.
        do n = 1, nNodL
            do nz = 1, mesh%nl-1
                tracers%work%del_ttf         (nz, n) = 0.0_WP
                tracers%work%del_ttf_advhoriz(nz, n) = 0.0_WP
                tracers%work%del_ttf_advvert (nz, n) = 0.0_WP
            end do
        end do

        ! Adams-Bashforth interpolation -> valuesAB
        if (tracers%data(tr_num)%AB_order == 2) then
            do n = 1, nNodL
                tracers%data(tr_num)%valuesAB(:, n) = -(0.5_WP+ab_epsilon)*tracers%data(tr_num)%valuesold(1, :, n) &
                                                     + (1.5_WP+ab_epsilon)*tracers%data(tr_num)%values(:, n)
            end do
        elseif (tracers%data(tr_num)%AB_order == 3) then
            do n = 1, nNodL
                tracers%data(tr_num)%valuesAB(:, n) = 5.0_WP*tracers%data(tr_num)%valuesold(2, :, n) &
                                                    - 16.0_WP*tracers%data(tr_num)%valuesold(1, :, n) &
                                                    + 23.0_WP*tracers%data(tr_num)%values(:, n)
                tracers%data(tr_num)%valuesAB(:, n) = tracers%data(tr_num)%valuesAB(:, n)/12.0_WP
            end do
        else
            error stop 'init_tracers_AB: AB_order must be 2 or 3 (check namelist.tra)'
        end if

        ! roll the history
        if (tracers%data(tr_num)%AB_order == 2) then
            do n = 1, nNodL
                tracers%data(tr_num)%valuesold(1, :, n) = tracers%data(tr_num)%values(:, n)
            end do
        elseif (tracers%data(tr_num)%AB_order == 3) then
            do n = 1, nNodL
                tracers%data(tr_num)%valuesold(2, :, n) = tracers%data(tr_num)%valuesold(1, :, n)
                tracers%data(tr_num)%valuesold(1, :, n) = tracers%data(tr_num)%values(:, n)
            end do
        end if

        ! elemental gradient of values over the FULL element halo: the MUSCL kernels take
        ! the up/downwind reconstruction from it on the fly (the reference RK3 rule, no
        ! stored per-edge array), and the horizontal diffusion / Redi of this tracer reuses
        ! it (solve_tracers_ale). Persistent scratch, allocated on first use. Zeroed once:
        ! tracer_gradient_elements writes only wet levels, and the dry ones are never read.
        nElemA = mesh%elem2D
        if (is_multirank(partit)) &
            nElemA = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
        if (.not. allocated(tracers%work%tr_xy)) then
            allocate(tracers%work%tr_xy(2, mesh%nl-1, nElemA))
            tracers%work%tr_xy = 0.0_WP
        end if
        call tracer_gradient_elements(tracers%data(tr_num)%values, tracers%work%tr_xy, mesh, partit)
        if (is_multirank(partit)) call exchange_elem_full(tracers%work%tr_xy, partit)
    end subroutine init_tracers_AB

end module oce_tracer_mod
