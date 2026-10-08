module oce_adv_tra_driver
    ! Tracer-advection driver, transcribed from FESOM2 v2.7.3
    ! oce_adv_tra_driver.F90:46-419 (subroutine do_oce_adv_tra).
    !
    ! Assembles the byte-proven kernels (M1.1 horizontal, M1.2 vertical, M1.3 FCT
    ! limiter, the edge->node scatter) into the per-tracer advection operator. For a
    ! tracer with tra_adv_lim=='FCT' it runs the Zalesak path:
    !   LO = UPW1(values) horiz + vert            -> fct_LO  (low-order solution)
    !   HO = {hor scheme, ver scheme}(valuesAB)   -> antidiffusive flux (HO-LO)
    !   oce_tra_adv_fct clips the antidiffusive flux
    !   oce_tra_adv_flux2dtracer(use_lo) scatters -> del_ttf_advhoriz/advvert
    ! For a non-FCT tracer (do_zero_flux=.true.) the HO scheme is applied directly to
    ! valuesAB (o_init_zero=.true.) and scattered without the low-order reconstruction.
    !
    ! Dispatch is on the per-tracer scheme selectors (tracers%data(tr_num)%tra_adv_*):
    !   tra_adv_hor : UPW1 / MUSCL
    !   tra_adv_ver : UPW1 / QR4C            (CDIFF / PPM declared but not yet ported)
    !   tra_adv_lim : FCT or anything else (non-FCT)
    ! tra_adv_ph / tra_adv_pv are the horizontal / vertical 4th-order fractions
    ! (opth / optv) passed to MUSCL/QR4C.
    !
    ! SCOPE / clean-architecture deviations (D7, "USE-globals -> explicit args"):
    !  * 1-rank only (myDim_* == global). The two FESOM2 halo exchanges
    !    (exchange_nod(fct_LO) after the LO build, the implicit ones in the kernels)
    !    are no-ops at 1 rank and are dropped here; multi-rank is M1.5.
    !  * use_wsplit (w = w_e + w_i, oce_wsplit / compute_Wvel_split): the FCT low-order
    !    solution is advanced explicitly with w_e, then implicitly with w_i
    !    (adv_tra_vert_impl, upwind backward Euler on fct_LO), and the low-order upwind
    !    flux is RECOMPUTED with the full w so that the antidiffusive flux is
    !    HO(w) - LO(w) (FESOM2 :282-292; wiring pinned by test_wimpl_tra C6). With
    !    use_wsplit=.false. w_e == w and w_i == 0 (compute_Wvel_split), the branch is
    !    skipped and the LO/HO vertical fluxes use the same field, as in FESOM2.
    !  * the DVD (downgradient-variance-decomposition) diagnostic blocks are dropped (not v1).
    !  * work arrays live in tracers%work (MP); kernels take WP. At the DP/SP anchor
    !    MP==WP so the pointers/args bind by kind value; revisit for FP16 (D3).
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use mod_partit,       only: t_partit
    use mod_part_bounds,  only: owned_bounds, is_multirank
    use mod_halo,         only: exchange_nod
    use oce_adv_tra_hor,  only: adv_tra_hor_upw1, adv_tra_hor_muscl
    use oce_adv_tra_ver,  only: adv_tra_ver_upw1, adv_tra_ver_qr4c, adv_tra_vert_impl
    use oce_adv_tra_fct,  only: oce_tra_adv_fct
    use oce_adv_tra_flux, only: oce_tra_adv_flux2dtracer
    implicit none
    private
    public :: do_oce_adv_tra

contains

    subroutine do_oce_adv_tra(dt, vel, w, wi, we, tr_num, dynamics, tracers, mesh, partit)
        ! do_oce_adv_tra (oce_adv_tra_driver.F90:46-419). vel = horizontal velocity at
        ! elements; w / we / wi = full / explicit / implicit vertical velocity at nodes.
        ! M2.12b: optional partit -> owned/halo loop bounds + exchange_nod(fct_LO) after
        ! the low-order solution (FESOM2 :294), passed on to every kernel.
        real(kind=WP),  intent(in)            :: dt
        type(t_mesh),   intent(in)            :: mesh
        type(t_dyn),    intent(in)            :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        integer,        intent(in)            :: tr_num
        real(kind=WP),  intent(in)            :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP),  intent(in), target    :: w (mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(in)            :: wi(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(in), target    :: we(mesh%nl, mesh%nod2D)
        type(t_partit), intent(in), optional  :: partit

        ! pointers into the per-tracer data / shared work (mirrors FESOM2 lines 92-106)
        real(kind=WP), pointer :: ttf(:,:), ttfAB(:,:), fct_LO(:,:)
        real(kind=WP), pointer :: adv_flux_hor(:,:), adv_flux_ver(:,:), dttf_h(:,:), dttf_v(:,:)
        real(kind=WP), pointer :: fct_ttf_min(:,:), fct_ttf_max(:,:), fct_plus(:,:), fct_minus(:,:)
        real(kind=WP), pointer :: tr_xy(:,:,:)
        integer,       pointer :: edge_up_dn_tri(:,:)
        real(kind=WP), pointer :: pwvel(:,:)

        real(kind=WP) :: opth, optv
        logical       :: do_zero_flux
        integer       :: e, n, nz, enodes(2), el(2), nl1, nl2, nu1, nu2, nl12, nu12
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ttf             => tracers%data(tr_num)%values
        ttfAB           => tracers%data(tr_num)%valuesAB
        opth            =  tracers%data(tr_num)%tra_adv_ph
        optv            =  tracers%data(tr_num)%tra_adv_pv
        fct_LO          => tracers%work%fct_LO
        adv_flux_ver    => tracers%work%adv_flux_ver
        adv_flux_hor    => tracers%work%adv_flux_hor
        tr_xy           => tracers%work%tr_xy
        edge_up_dn_tri  => tracers%work%edge_up_dn_tri
        fct_ttf_min     => tracers%work%fct_ttf_min
        fct_ttf_max     => tracers%work%fct_ttf_max
        fct_plus        => tracers%work%fct_plus
        fct_minus       => tracers%work%fct_minus
        dttf_h          => tracers%work%del_ttf_advhoriz
        dttf_v          => tracers%work%del_ttf_advvert

        !_______________________________________________________________________
        ! FCT: low-order horizontal+vertical solution + low-order antidiffusive part
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            ! low-order upwind horizontal flux (zero the flux first)
            call adv_tra_hor_upw1(vel, ttf, mesh, adv_flux_hor, o_init_zero=.true., partit=partit)

            ! fct_LO = scatter of the LO horizontal flux over edges
            do n = 1, nNodL
                do nz = 1, mesh%nl-1
                    fct_LO(nz,n) = 0.0_WP
                end do
            end do
            do e = 1, nEdgeO
                enodes = mesh%edges(:, e)
                el     = mesh%edge_tri(:, e)
                nl1 = mesh%nlevels(el(1))-1
                nu1 = mesh%ulevels(el(1))
                nl2 = 0; nu2 = 0
                if (el(2) /= el(1)) then
                    nl2 = mesh%nlevels(el(2))-1
                    nu2 = mesh%ulevels(el(2))
                end if
                nl12 = max(nl1, nl2)
                nu12 = nu1
                if (nu2 > 0) nu12 = min(nu1, nu2)
                do nz = nu12, nl12
                    fct_LO(nz, enodes(1)) = fct_LO(nz, enodes(1)) + adv_flux_hor(nz, e)
                    fct_LO(nz, enodes(2)) = fct_LO(nz, enodes(2)) - adv_flux_hor(nz, e)
                end do
            end do

            ! low-order upwind vertical flux (explicit part: we), then finish the LO solution
            call adv_tra_ver_upw1(we, ttf, mesh, adv_flux_ver, o_init_zero=.true., partit=partit)
            do n = 1, nNodO
                nu1 = mesh%ulevels_nod2D(n)
                nl1 = mesh%nlevels_nod2D(n)
                do nz = nu1, nl1-1
                    fct_LO(nz,n) = (ttf(nz,n)*mesh%hnode(nz,n) &
                                  + (fct_LO(nz,n) + (adv_flux_ver(nz,n)-adv_flux_ver(nz+1,n)))*dt/mesh%areasvol(n)) &
                                  / mesh%hnode_new(nz,n)
                end do
            end do

            if (dynamics%use_wsplit) then
                ! w-split (FESOM2 :282-292): the implicit part w_i advects the low-order
                ! solution implicitly (upwind backward Euler, in place), then the low-order
                ! upwind vertical flux is recomputed with the FULL w, so that the HO call
                ! below (o_init_zero=.false.) leaves HO(w) - LO(w) as the antidiffusive flux.
                ! (LO(w_e) there would double-count the w_i transport; test_wimpl_tra C6.)
                call adv_tra_vert_impl(dt, wi, fct_LO, mesh, partit=partit)
                call adv_tra_ver_upw1(w, ttf, mesh, adv_flux_ver, o_init_zero=.true., partit=partit)
            end if
            ! M2.12b: share the low-order solution to the halo (FESOM2 :294) — the FCT
            ! limiter's a1 reads fct_LO at owned+halo nodes.
            if (is_multirank(partit)) call exchange_nod(fct_LO, partit)
        end if

        do_zero_flux = .true.
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') do_zero_flux = .false.

        !_______________________________________________________________________
        ! horizontal advection (FCT: high-order antidiffusive; else: full high-order)
        select case (trim(tracers%data(tr_num)%tra_adv_hor))
        case ('MUSCL')
            call adv_tra_hor_muscl(vel, ttfAB, mesh, opth, adv_flux_hor, tr_xy, edge_up_dn_tri, o_init_zero=do_zero_flux, partit=partit)
        case ('UPW1')
            call adv_tra_hor_upw1 (vel, ttfAB, mesh,       adv_flux_hor,                                 o_init_zero=do_zero_flux, partit=partit)
        case default
            error stop 'do_oce_adv_tra: unknown tra_adv_hor (expected UPW1/MUSCL)'
        end select

        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            pwvel => w
        else
            pwvel => we
        end if

        !_______________________________________________________________________
        ! vertical advection (FCT: high-order antidiffusive; else: full high-order)
        select case (trim(tracers%data(tr_num)%tra_adv_ver))
        case ('QR4C')
            call adv_tra_ver_qr4c(pwvel, ttfAB, mesh, optv, adv_flux_ver, o_init_zero=do_zero_flux, partit=partit)
        case ('UPW1')
            call adv_tra_ver_upw1(pwvel, ttfAB, mesh,       adv_flux_ver, o_init_zero=do_zero_flux, partit=partit)
        case ('CDIFF', 'PPM')
            error stop 'do_oce_adv_tra: tra_adv_ver CDIFF/PPM not yet ported (M2)'
        case default
            error stop 'do_oce_adv_tra: unknown tra_adv_ver (expected UPW1/QR4C)'
        end select

        !_______________________________________________________________________
        ! limit (FCT) and scatter to the per-node tendencies del_ttf_advhoriz/advvert
        if (trim(tracers%data(tr_num)%tra_adv_lim) == 'FCT') then
            call oce_tra_adv_fct(dt, ttf, fct_LO, adv_flux_hor, adv_flux_ver, &
                                 fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, mesh, partit=partit)
            call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh, &
                                          use_lo=.true., ttf=ttf, lo=fct_LO, partit=partit)
        else
            call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, adv_flux_hor, adv_flux_ver, mesh, partit=partit)
        end if
    end subroutine do_oce_adv_tra

end module oce_adv_tra_driver
