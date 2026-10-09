module oce_adv_tra_hor
    ! Horizontal tracer advection — edge finite-volume fluxes, transcribed from
    ! FESOM2 v2.7.3 oce_adv_tra_hor.F90:
    !   adv_tra_hor_upw1  (64-257)  low-order upwind
    !   adv_tra_hor_muscl (546-834) MUSCL (3rd/4th-order) = FESOM2's adv_tra_hor_mfct.
    !     FESOM2's clamped adv_tra_hor_muscl (261-542, nboundary_lay) is not ported: as in
    !     the reference RK3 code there is a single MUSCL kernel without a boundary clamp.
    ! The MUSCL kernel forms the up/downwind reconstruction ON THE FLY, as the reference RK3
    ! code (qq/oce_stepRK3.F90, t_hor_adv_muscl_RK3): the elemental gradient tr_xy of the
    ! upwind / downwind element (edge_up_dn_tri; triangle or quad) where it is wet at the level,
    ! a zero increment where it is missing or dry. FESOM2 instead stores a per-edge array
    ! edge_up_dn_grad with a Miura node-averaged gradient on those levels; FESOM3 dropped
    ! both (docs/plans/completed/2026-10-06-muscl-onthefly.md). Identical to FESOM2 on the
    ! levels where both edge nodes share both elements (test_muscl_onthefly part R).
    !
    ! Each returns an EDGE flux that contributes with +sign to the 1st edge node and
    ! -sign to the 2nd (the driver scatters it to nodes). o_init_zero=.true. zeroes
    ! the flux first; .false. SUBTRACTS the new contribution from the input flux
    ! (so a HO call after an LO call yields the antidiffusive flux). flux is NOT
    ! multiplied by dt.
    !
    ! 1-rank only (myDim_edge2D == edge2D; el(2)==el(1) marks a boundary edge — its
    ! el(2)-using depth segments have empty ranges so el(2) is never used).
    ! Multi-rank loop bounds are deferred to M1.5. The volume flux scales the edge
    ! velocity by helem (layer thickness at element); a = mean r_earth*elem_cos over
    ! the edge's two elements (used by the MUSCL reconstruction, dead in upw1 — kept
    ! to mirror FESOM2).
    use mod_precision,   only: WP
    use mod_constants,   only: r_earth
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: adv_tra_hor_upw1, adv_tra_hor_muscl

contains

    !===========================================================================
    subroutine adv_tra_hor_upw1(vel, ttf, mesh, flux, o_init_zero, partit)
        ! Low-order upwind horizontal flux (oce_adv_tra_hor.F90:64-257).
        ! M2.12b: optional partit -> loop over OWNED edges (myDim_edge2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, nEdgeO
                do nz = 1, mesh%nl-1
                    flux(nz, edge) = 0.0_WP
                end do
            end do
        end if

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            nl2 = 0; nu2 = 0
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A) el(1)-only top segments (cavity surface)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_lo(nz, vflux)
            end do
            ! (B) el(2)-only top segments (cavity surface)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_lo(nz, vflux)
                end do
            end if
            ! (C) both segments
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_lo(nz, vflux)
            end do
            ! (D) remaining el(1) segments
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_lo(nz, vflux)
            end do
            ! (E) remaining el(2) segments
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_lo(nz, vflux)
            end do
        end do

    contains
        subroutine flux_lo(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            flux(nz, edge) = -0.5_WP * ( ttf(nz, enodes(1))*(vflux+abs(vflux)) &
                                       + ttf(nz, enodes(2))*(vflux-abs(vflux)) ) - flux(nz, edge)
        end subroutine flux_lo
    end subroutine adv_tra_hor_upw1



    !===========================================================================
    subroutine adv_tra_hor_muscl(vel, ttf, mesh, num_ord, flux, tr_xy, edge_up_dn_tri, o_init_zero, partit)
        ! MUSCL horizontal flux (FESOM2 adv_tra_hor_mfct, oce_adv_tra_hor.F90:546-834).
        ! num_ord = fraction of 4th-order (centered) contribution; (1-num_ord) is
        ! 3rd-order upwind. No boundary clamp: where the up/downwind element is missing
        ! or dry the face value is 2/3*T_up + 1/3*T_dn. Not monotone; usually run under FCT.
        ! M2.12b: optional partit -> loop over OWNED edges (myDim_edge2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(in)    :: tr_xy(2, mesh%nl-1, mesh%elem2D)   ! elemental gradient, full halo
        integer,       intent(in)    :: edge_up_dn_tri(2, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        integer       :: tri1, tri2, u1, l1, u2, l2   ! up/dn elements and their wet level ranges

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, nEdgeO
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            ! up/downwind element of this edge and the levels on which it is wet; a missing
            ! element gets an empty range (u = 1, l = 0)
            tri1 = edge_up_dn_tri(1, edge); u1 = 1; l1 = 0
            tri2 = edge_up_dn_tri(2, edge); u2 = 1; l2 = 0
            if (tri1 > 0) then
                u1 = mesh%ulevels(tri1); l1 = mesh%nlevels(tri1) - 1
            end if
            if (tri2 > 0) then
                u2 = mesh%ulevels(tri2); l2 = mesh%nlevels(tri2) - 1
            end if
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            nl2 = 0; nu2 = 0
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
        end do

    contains
        subroutine flux_ho(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO, g1(2), g2(2)
            ! the reference RK3 rule (qq/oce_stepRK3.F90 t_hor_adv_muscl_RK3): the
            ! up/downwind element's gradient where that element is wet, else 0
            g1 = 0.0_WP; g2 = 0.0_WP
            if (nz >= u1 .and. nz <= l1) g1 = tr_xy(:, nz, tri1)
            if (nz >= u2 .and. nz <= l2) g2 = tr_xy(:, nz, tri2)
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*g2(1) + &
                      mesh%edge_dxdy(2,edge)*g2(2))/6.0_WP
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*g1(1) + &
                      mesh%edge_dxdy(2,edge)*g1(2))/6.0_WP
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine adv_tra_hor_muscl

end module oce_adv_tra_hor
