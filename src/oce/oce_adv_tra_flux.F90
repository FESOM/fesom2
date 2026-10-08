module oce_adv_tra_flux
    ! Scatter the edge/interface advective fluxes onto the per-node tracer-change
    ! arrays. Transcribed from FESOM2 v2.7.3 oce_adv_tra_driver.F90:423-575
    ! (subroutine oce_tra_adv_flux2dtracer).
    !
    !   dttf_h(nz,node) = sum over edges of  +/- flux_h(nz,edge)*dt/areasvol(nz,node)
    !   dttf_v(nz,node) =                     (flux_v(nz)-flux_v(nz+1))*dt/areasvol
    !
    ! The horizontal scatter accumulates over edges in edge order; with FESOM2's
    ! edge list (byte-identical here) the per-node summation order matches, so the
    ! result is byte-identical. dttf_h/dttf_v are NOT zeroed here (the caller does,
    ! as FESOM2 does in init_tracers_AB).
    !
    ! The use_lo / hnode branch (FCT low-order ALE reconstruction) is transcribed
    ! faithfully but is dead in M1.1 (no FCT, no ALE thickness yet): M1.1 calls this
    ! with flux_v=0 and without use_lo, so only the horizontal scatter is live and
    ! dttf_v stays 0. mesh%hnode/hnode_new are referenced only in that dead branch.
    ! 1-rank only (myDim_* == global); multi-rank halo of dttf is deferred to M1.5.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: oce_tra_adv_flux2dtracer

contains

    subroutine oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, flux_h, flux_v, mesh, use_lo, ttf, lo, partit)
        ! M2.12b: optional partit -> scatter over OWNED nodes/edges (myDim_*).
        real(kind=WP), intent(in)    :: dt
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(inout) :: dttf_h(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: dttf_v(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: flux_h(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(inout) :: flux_v(mesh%nl,   mesh%nod2D)
        logical,       optional      :: use_lo
        real(kind=WP), optional      :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), optional      :: lo (mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer :: n, nz, el(2), enodes(3), nu12, nl12, nu1, nu2, nl1, nl2, edge
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        ! Vertical (FCT low-order ALE reconstruct) — dead in M1.1
        if (present(use_lo)) then
            if (use_lo) then
                do n = 1, nNodO
                    nu1 = mesh%ulevels_nod2D(n)
                    nl1 = mesh%nlevels_nod2D(n)
                    do nz = nu1, nl1-1
                        dttf_v(nz,n) = dttf_v(nz,n) - ttf(nz,n)*mesh%hnode(nz,n) + lo(nz,n)*mesh%hnode_new(nz,n)
                    end do
                end do
            end if
        end if
        ! Vertical flux divergence
        do n = 1, nNodO
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            do nz = nu1, nl1-1
                dttf_v(nz,n) = dttf_v(nz,n) + (flux_v(nz,n)-flux_v(nz+1,n))*dt/mesh%areasvol(n)
            end do
        end do
        ! Horizontal edge-flux scatter -> del_ttf_advhoriz
        do edge = 1, nEdgeO
            enodes(1:2) = mesh%edges(:, edge)
            el = mesh%edge_tri(:, edge)
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
                dttf_h(nz,enodes(1)) = dttf_h(nz,enodes(1)) + flux_h(nz,edge)*dt/mesh%areasvol(enodes(1))
                dttf_h(nz,enodes(2)) = dttf_h(nz,enodes(2)) - flux_h(nz,edge)*dt/mesh%areasvol(enodes(2))
            end do
        end do
    end subroutine oce_tra_adv_flux2dtracer

end module oce_adv_tra_flux
