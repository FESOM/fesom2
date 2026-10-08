module oce_adv_tra_fct
    ! 3D Flux-Corrected-Transport (Zalesak) limiter, transcribed from FESOM2 v2.7.3
    ! oce_adv_tra_fct.F90 (subroutine oce_tra_adv_fct, lines 72-500).
    !
    ! Limits the antidiffusive fluxes (the difference HO-LO already stored in adf_h/
    ! adf_v by the o_init_zero=.false. HO kernel calls) so that adding them to the
    ! low-order solution `lo` does not create new extrema. Steps (FESOM2 labels):
    !   a1  per-node max/min of (lo, ttf)                         -> fct_ttf_max/min
    !   a2  per-element admissible bounds over its 4 vertex slots -> AUX(1:2,nz,elem)
    !       (layers at/below the element bottom set to -/+bignumber: unconstrained)
    !   a3  cluster (node + vertical nz-1:nz+1) bounds            -> tvert_max/min
    !       then admissible increment w.r.t. lo                   -> fct_ttf_max/min
    !   b1  split antidiffusive contributions into fct_plus/minus (vertical + horiz)
    !   b2  nodal limiting factors  min(1, bound / (flux*dt/areasvol/hnode_new))
    !   b3  clip adf_v and adf_h by the upwind-side limiting factor
    !
    ! On the gate adf_h/adf_v ENTER as the antidiffusive fluxes and LEAVE clipped;
    ! fct_ttf_max/min leave as the admissible increments (a3); fct_plus/minus leave
    ! as the limiting factors (b2). The caller then scatters the clipped fluxes with
    ! oce_tra_adv_flux2dtracer(... use_lo=.TRUE., ttf, lo).
    !
    ! Deviations from FESOM2, both byte-identical on pi (1-rank, no cavity):
    !  - AUX is a LOCAL allocatable here; FESOM2 reuses twork%edge_up_dn_grad as
    !    scratch. a2 writes AUX(1:2, ulevels(elem):nl-1, elem) for every element, and
    !    a3 only reads AUX(:,nz,elem) for nz in [ulevels_nod2D(n), nlevels_nod2D(n)-1];
    !    with ulevels==1 everywhere (no cavity) every read entry is written first, so
    !    the (uninitialised) initial AUX values never reach a result. On a CAVITY mesh
    !    a node could read AUX(:,nz,elem) at nz < ulevels(elem) (unwritten) -> FESOM2
    !    read stale edge_up_dn_grad there. FIXED (FESOM2 d8dcdfb): the a2 loop now pads
    !    ABOVE the cavity too, so every entry a3 can reach is written. Bit-unchanged
    !    without a cavity; the cavity path still wants a gate on a cavity mesh (M2+).
    !  - exchange_nod(fct_plus,fct_minus) (FESOM2 b2->b3) is a no-op at 1 rank, dropped.
    !  - dmax1/dmin1 -> generic max/min: identical IEEE result on real(WP=8) at the
    !    anchor, and portable to a single-precision build.
    !
    ! Runtime divisors areasvol / hnode_new (b2) match FESOM2 because both codes divide
    ! by byte-identical operands (geom-proven areasvol; gated hnode_new); cf. LESSONS L7.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    implicit none
    private
    public :: oce_tra_adv_fct

contains

    subroutine oce_tra_adv_fct(dt, ttf, lo, adf_h, adf_v, fct_ttf_min, fct_ttf_max, &
                               fct_plus, fct_minus, mesh, partit)
        ! M2.12b: optional partit. a1 runs over OWNED+HALO nodes (myDim+eDim, reads the
        ! exchanged lo / halo ttf); a2 over OWNED elements; a3/b1/b2/b3 over OWNED
        ! nodes/edges. exchange_nod(fct_plus,fct_minus) between b2 and b3 (FESOM2 :401)
        ! so the owned-edge b3 limiting reads correct halo-node limiting factors.
        real(kind=WP), intent(in)    :: dt
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: lo (mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: adf_h(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(inout) :: adf_v(mesh%nl,   mesh%nod2D)
        real(kind=WP), intent(inout) :: fct_ttf_min(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: fct_ttf_max(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: fct_plus (mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: fct_minus(mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit

        integer :: n, nz, elem, enodes(MAX_NV), el(2), nl1, nl2, nu1, nu2, nl12, nu12, edge
        real(kind=WP) :: flux, ae
        real(kind=WP), allocatable :: tvert_max(:,:), tvert_min(:,:), AUX(:,:,:)
        real(kind=WP) :: flux_eps=1e-16
        real(kind=WP) :: bignumber=1e3
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        ! local scratch: tvert over owned+halo nodes; AUX over OWNED elements (an owned
        ! node's element list is all owned, so a3 only reads AUX at owned elements).
        allocate(tvert_max(mesh%nl-1, nNodL), tvert_min(mesh%nl-1, nNodL))
        allocate(AUX(2, mesh%nl-1, nElemO))

        !_______________________________________________________________________
        ! a1. max, min between old solution and updated low-order solution per node
        do n=1, nNodL
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            do nz=nu1, nl1-1
                fct_ttf_max(nz,n)=max(lo(nz,n), ttf(nz,n))
                fct_ttf_min(nz,n)=min(lo(nz,n), ttf(nz,n))
            end do
        end do

        !_______________________________________________________________________
        ! a2. Admissible increments on elements (max/min bound per element). Layers
        !     at/below the element bottom (nz>=nlevels(elem)-1) are set to -/+bignumber
        !     so a shallow element does not constrain a deeper node's bounds.
        do elem=1, nElemO
            enodes = mesh%elem2D_nodes(:, elem)   ! a triangle's slot 4 repeats slot 1 (no effect on max/min)
            nu1 = mesh%ulevels(elem)
            nl1 = mesh%nlevels(elem)
            ! FESOM2 d8dcdfb: pad ABOVE a cavity element as well as below it. a3 reads
            ! AUX at nz in [ulevels_nod2D(n), nlevels_nod2D(n)-1], which on a cavity mesh
            ! reaches nz < ulevels(elem) -- never written, so an UNINITIALISED value used
            ! to enter the max/min reduction (first-timestep NaNs upstream; AUX is a fresh
            ! local allocatable here, so it is uninitialised rather than stale). Neutral
            ! fill, same convention as the below-bottom pad. Zero iterations without a
            ! cavity, so core2/pi are bit-unchanged.
            do nz=1, nu1-1
                AUX(1,nz,elem)=-bignumber
                AUX(2,nz,elem)= bignumber
            end do
            do nz=nu1, nl1-1
                AUX(1,nz,elem)=max(fct_ttf_max(nz,enodes(1)), fct_ttf_max(nz,enodes(2)), fct_ttf_max(nz,enodes(3)), &
                                   fct_ttf_max(nz,enodes(4)))
                AUX(2,nz,elem)=min(fct_ttf_min(nz,enodes(1)), fct_ttf_min(nz,enodes(2)), fct_ttf_min(nz,enodes(3)), &
                                   fct_ttf_min(nz,enodes(4)))
            end do
            if (nl1<=mesh%nl-1) then
                do nz=nl1, mesh%nl-1
                    AUX(1,nz,elem)=-bignumber
                    AUX(2,nz,elem)= bignumber
                end do
            endif
        end do

        !_______________________________________________________________________
        ! a3. Bounds on clusters (node neighbourhood) and admissible increments.
        do n=1, nNodO
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            do nz=nu1, nl1-1
                tvert_max(nz, n) = AUX(1,nz, mesh%nod_in_elem2D(1, n))
                tvert_min(nz, n) = AUX(2,nz, mesh%nod_in_elem2D(1, n))
                do elem=2, mesh%nod_in_elem2D_num(n)
                    tvert_max(nz, n) = max(tvert_max(nz, n), AUX(1,nz, mesh%nod_in_elem2D(elem,n)))
                    tvert_min(nz, n) = min(tvert_min(nz, n), AUX(2,nz, mesh%nod_in_elem2D(elem,n)))
                end do
            end do
        end do

        do n=1, nNodO
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            ! surface layer increment w.r.t. low-order solution
            fct_ttf_max(nu1,n)=tvert_max(nu1, n)-lo(nu1,n)
            fct_ttf_min(nu1,n)=tvert_min(nu1, n)-lo(nu1,n)
            ! interior: increment from nz-1:nz+1
            do nz=nu1+1, nl1-2
                fct_ttf_max(nz,n)=max(tvert_max(nz-1, n), tvert_max(nz, n), tvert_max(nz+1, n))-lo(nz,n)
                fct_ttf_min(nz,n)=min(tvert_min(nz-1, n), tvert_min(nz, n), tvert_min(nz+1, n))-lo(nz,n)
            end do
            ! bottom layer increment
            nz=nl1-1
            fct_ttf_max(nz,n)=tvert_max(nz, n)-lo(nz,n)
            fct_ttf_min(nz,n)=tvert_min(nz, n)-lo(nz,n)
        end do

        !_______________________________________________________________________
        ! b1. Split positive (fct_plus) and negative (fct_minus) antidiffusive
        !     contributions, accumulated per node from vertical and horizontal fluxes.
        do n=1, nNodO
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            do nz=nu1, nl1-1
                fct_plus(nz,n)=0._WP
                fct_minus(nz,n)=0._WP
            end do
        end do
        ! Vertical
        do n=1, nNodO
            nu1 = mesh%ulevels_nod2D(n)
            nl1 = mesh%nlevels_nod2D(n)
            do nz=nu1, nl1-1
                fct_plus(nz,n) =fct_plus(nz,n) +(max(0.0_WP,adf_v(nz,n))+max(0.0_WP,-adf_v(nz+1,n)))
                fct_minus(nz,n)=fct_minus(nz,n)+(min(0.0_WP,adf_v(nz,n))+min(0.0_WP,-adf_v(nz+1,n)))
            end do
        end do
        ! Horizontal
        do edge=1, nEdgeO
            enodes(1:2)=mesh%edges(:,edge)
            el=mesh%edge_tri(:,edge)
            nl1=mesh%nlevels(el(1))-1
            nu1=mesh%ulevels(el(1))
            nl2=0
            nu2=0
            if (el(2)/=el(1)) then
                nl2=mesh%nlevels(el(2))-1
                nu2=mesh%ulevels(el(2))
            end if
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
            do nz=nu12, nl12
                fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0_WP, adf_h(nz,edge))
                fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0_WP, adf_h(nz,edge))
                fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0_WP,-adf_h(nz,edge))
                fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0_WP,-adf_h(nz,edge))
            end do
        end do

        !_______________________________________________________________________
        ! b2. Limiting factors
        do n=1, nNodO
            nu1=mesh%ulevels_nod2D(n)
            nl1=mesh%nlevels_nod2D(n)
            do nz=nu1, nl1-1
                flux=fct_plus(nz,n)*dt/mesh%areasvol(n)/mesh%hnode_new(nz,n)+flux_eps
                fct_plus(nz,n)=min(1.0_WP,fct_ttf_max(nz,n)/flux)
                flux=fct_minus(nz,n)*dt/mesh%areasvol(n)/mesh%hnode_new(nz,n)-flux_eps
                fct_minus(nz,n)=min(1.0_WP,fct_ttf_min(nz,n)/flux)
            end do
        end do
        ! M2.12b: the owned-edge b3 limiting reads fct_plus/fct_minus at halo nodes
        ! (an owned edge can touch a halo node) -> exchange owner->halo (FESOM2 :401).
        if (is_multirank(partit)) then
            call exchange_nod(fct_plus,  partit)
            call exchange_nod(fct_minus, partit)
        end if

        !_______________________________________________________________________
        ! b3. Limiting
        ! Vertical
        do n=1, nNodO
            nu1=mesh%ulevels_nod2D(n)
            nl1=mesh%nlevels_nod2D(n)
            ! surface interface
            nz=nu1
            ae=1.0_WP
            flux=adf_v(nz,n)
            if(flux>=0.0_WP) then
                ae=min(ae,fct_plus(nz,n))
            else
                ae=min(ae,fct_minus(nz,n))
            end if
            adf_v(nz,n)=ae*adf_v(nz,n)
            ! interior interfaces
            do nz=nu1+1, nl1-1
                ae=1.0_WP
                flux=adf_v(nz,n)
                if(flux>=0._WP) then
                    ae=min(ae,fct_minus(nz-1,n))
                    ae=min(ae,fct_plus(nz,n))
                else
                    ae=min(ae,fct_plus(nz-1,n))
                    ae=min(ae,fct_minus(nz,n))
                end if
                adf_v(nz,n)=ae*adf_v(nz,n)
            end do
            ! the bottom flux is always zero
        end do
        ! Horizontal
        do edge=1, nEdgeO
            enodes(1:2)=mesh%edges(:,edge)
            el=mesh%edge_tri(:,edge)
            nu1=mesh%ulevels(el(1))
            nl1=mesh%nlevels(el(1))-1
            nl2=0
            nu2=0
            if(el(2)/=el(1)) then
                nu2=mesh%ulevels(el(2))
                nl2=mesh%nlevels(el(2))-1
            end if
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
            do nz=nu12, nl12
                ae=1.0_WP
                flux=adf_h(nz,edge)
                if(flux>=0._WP) then
                    ae=min(ae,fct_plus(nz,enodes(1)))
                    ae=min(ae,fct_minus(nz,enodes(2)))
                else
                    ae=min(ae,fct_minus(nz,enodes(1)))
                    ae=min(ae,fct_plus(nz,enodes(2)))
                endif
                adf_h(nz,edge)=ae*adf_h(nz,edge)
            end do
        end do

        deallocate(tvert_max, tvert_min, AUX)
    end subroutine oce_tra_adv_fct

end module oce_adv_tra_fct
