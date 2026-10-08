module oce_dyn_velrhs
    ! Momentum-RHS assembly: Coriolis + Adams-Bashforth(2) blend + hydrostatic
    ! pressure-gradient force + sea-surface-height (SSH) gradient + momentum advection.
    ! Transcribed from FESOM2 v2.7.3 oce_ale_vel_rhs.F90 (compute_vel_rhs 35-328 +
    ! momentum_adv_scalar 335-589).
    !
    ! M2.3 ported the Coriolis/AB2/PGF/SSH-gradient pieces; M2.4 adds momentum_adv_scalar
    ! (called from inside compute_vel_rhs at the momadv_opt==2 site, FESOM2 :271-273),
    ! which ADDS the advective momentum flux into the SAME UV_rhsAB(1,1:2,:) slot as the
    ! Coriolis term, BEFORE the AB blend. The gate is now the FULL UV_rhs (max|delta|=0).
    !
    ! Pinned to the pi gated path (all v1-out / UNGATED branches dropped, re-added
    ! with their own gates later):
    !   - AB_order==2                    (AB_order==3 is M2.x)
    !   - which_ale='linfs'  => use_pice=0  -> no floating-ice loading p_ice
    !   - l_mslp=.false.                  -> no sea-level air pressure p_air
    !   - use_global_tides=.false.        -> no ssh_gp tidal term
    !   - use_ssh_se_subcycl=.false.      -> no split-explicit transport velocities
    !   - ldiag_ke=.false.                -> no kinetic-energy diagnostics
    !
    ! ALGORITHM (per element, mid-depth levels nz = ulevels..nlevels-1):
    !   ab1 = -(0.5 + eps),  ab2 = (1.5 + eps)          [eps = ab_epsilon = 0.1]
    !   (1) AB-init :  UV_rhs(:,nz)  = ab1 * UV_rhsAB_prev(:,nz)
    !   (2) SSH grad:  pre = -g*eta_n(elnodes) ; Fx/Fy = gradient_sca(1:4/5:8) . pre
    !                  ff  = coriolis(elem) * elem_area(elem)
    !   (3) PGF+Cor :  UV_rhs(1,nz) += (Fx - pgf_x(nz)) * elem_area
    !                  UV_rhs(2,nz) += (Fy - pgf_y(nz)) * elem_area
    !                  UV_rhsAB(1,nz) =  UV(2,nz)*ff ;  UV_rhsAB(2,nz) = -UV(1,nz)*ff
    !   (3b) momadv :  if momadv_opt==2, momentum_adv_scalar ADDS w*du/dz + u*du/dx
    !                  into the SAME UV_rhsAB(1,1:2,nz,elem) slot (M2.4)
    !   (4) update  :  ff = ab2 (or 1.0 on the first Euler step) ;
    !                  UV_rhs(:,nz) = dt*(UV_rhs(:,nz) + UV_rhsAB(:,nz)*ff)/elem_area
    !
    ! BIT-IDENTITY NOTES:
    !  - epsilon -> ab_epsilon (mod_config; non-parameter so 0.5+0.1->0.6 / 1.5+0.1->
    !    1.6 are not compile-folded; L12 proved both sums fold-safe regardless).
    !  - coriolis is a geometry field (2*omega*sin(lat_geo) via r2g, mod_mesh_areas);
    !    elem_area and gradient_sca are geometry-gate-proven (max|delta|=0). pgf_x/pgf_y
    !    come from oce_pgf (M2.2, gated). So every operand is already byte-pinned.
    !  - the final /elem_area is a RUNTIME divisor but byte-identical on both sides
    !    (geom-proven operand), so the -no-prec-div reciprocal matches (L7/L10/L14).
    !  - UV_rhsAB layout matches FESOM2: (AB_order-1, 2, nl-1, elem2D) -> the index
    !    convention (ab_lvl=1, comp, nz, elem) for AB_order==2.
    !  - 1-rank only: the ported part of compute_vel_rhs has NO halo exchange (those
    !    live in the skipped momentum_adv_scalar); UV_rhs is element-local.
    use mod_precision, only: WP
    use mod_mesh,      only: t_mesh, MAX_NV
    use mod_dyn,       only: t_dyn
    use mod_constants, only: g
    use mod_config,    only: ab_epsilon
    use mod_partit,    only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,      only: exchange_nod
    use oce_dyn_vinv,  only: momentum_adv_vinv
    implicit none
    private
    public :: compute_vel_rhs

contains

    subroutine compute_vel_rhs(dynamics, mesh, dt, lfirst, partit)
        ! dynamics: uv/uv_rhs/uv_rhsAB/eta_n + work%pgf_x/pgf_y (intent inout).
        ! lfirst: first Euler timestep (ff=1.0). FESOM2's guard is
        !   lfirst .and. .not. r_restart ; v1 has no restart yet (M2.11), so r_restart
        !   is folded out into the lfirst argument supplied by the caller.
        type(t_dyn),   intent(inout), target :: dynamics
        type(t_mesh),  intent(in),    target :: mesh
        real(kind=WP), intent(in) :: dt
        logical,       intent(in) :: lfirst
        type(t_partit), intent(in), optional :: partit
        !_______________________________________________________________________
        integer       :: elem, elnodes(MAX_NV), nz, nzmax, nzmin
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: ff, Fx, Fy, pre(MAX_NV), p_eta(MAX_NV)
        real(kind=WP) :: ab1, ab2
        real(kind=WP), dimension(:,:,:),   pointer :: UV, UV_rhs
        real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB
        real(kind=WP), dimension(:),       pointer :: eta_n
        real(kind=WP), dimension(:,:),     pointer :: pgf_x, pgf_y

        UV       => dynamics%uv
        UV_rhs   => dynamics%uv_rhs
        UV_rhsAB => dynamics%uv_rhsAB
        eta_n    => dynamics%eta_n
        pgf_x    => dynamics%work%pgf_x
        pgf_y    => dynamics%work%pgf_y
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! 2nd-order Adams-Bashforth coefficients (FESOM2 :97-100). eps = ab_epsilon=0.1
        ab1 = -(0.5_WP + ab_epsilon)
        ab2 =  (1.5_WP + ab_epsilon)

        !_______________________________________________________________________
        ! Coriolis + AB2 + PGF (+ SSH gradient) assembly
        do elem = 1, nElemO
            nzmax = mesh%nlevels(elem)
            nzmin = mesh%ulevels(elem)

            ! (1) AB part: init UV_rhs from the previous-step AB array (AB2)
            do nz = nzmin, nzmax-1
                UV_rhs(1,nz,elem) = ab1*UV_rhsAB(1,1,nz,elem)
                UV_rhs(2,nz,elem) = ab1*UV_rhsAB(1,2,nz,elem)
            end do

            ! (2) SSH gradient contribution -g*grad(eta) + Coriolis factor.
            !     p_air (l_mslp) and p_ice (use_pice, linfs) are 0 on pi -> dropped.
            elnodes = mesh%elem2D_nodes(:,elem)     ! 4 slots (a triangle's slot 4 has gradient 0)
            p_eta   = g*eta_n(elnodes)
            ff      = mesh%coriolis(elem)*mesh%elem_area(elem)
            pre = -(p_eta)
            Fx  = sum(mesh%gradient_sca(1:4, elem)*pre)
            Fy  = sum(mesh%gradient_sca(5:8, elem)*pre)

            ! (3) add PGF + SSH-gradient to UV_rhs; init this step's AB array with
            !     the Coriolis term
            do nz = nzmin, nzmax-1
                UV_rhs(1,nz,elem) = UV_rhs(1,nz,elem) + (Fx-pgf_x(nz,elem))*mesh%elem_area(elem)
                UV_rhs(2,nz,elem) = UV_rhs(2,nz,elem) + (Fy-pgf_y(nz,elem))*mesh%elem_area(elem)
                UV_rhsAB(1,1,nz,elem) =  UV(2,nz,elem)*ff
                UV_rhsAB(1,2,nz,elem) = -UV(1,nz,elem)*ff
            end do
        end do

        !_______________________________________________________________________
        ! Momentum advection -> ADD into this-step UV_rhsAB(1,1:2,:) (FESOM2 :271-273),
        ! AFTER the Coriolis/PGF elem loop and BEFORE the AB blend. v1 has no
        ! split-explicit subcycling (use_ssh_se_subcycl=.false.), so momadv_opt==2
        ! routes to momentum_adv_scalar (the _transpv variant is M2.x). momadv_opt==1 is
        ! the VECTOR-INVARIANT form (oce_dyn_vinv): FESOM2 aborts on it (never adapted to
        ! ALE, oce_ale_vel_rhs.F90:268), so ours has no oracle -- see that module's header.
        ! Both branches only ADD into UV_rhsAB, so (4) below is common to them.
        if (dynamics%momadv_opt == 1) then
            call momentum_adv_vinv(dynamics, mesh, partit)
        else if (dynamics%momadv_opt == 2) then
            call momentum_adv_scalar(dynamics, mesh, partit)
        end if

        !_______________________________________________________________________
        ! (4) AB blend + scale by dt/elem_area. First Euler step -> ff=1.0 else ab2.
        ff = ab2
        if (lfirst) ff = 1.0_WP
        do elem = 1, nElemO
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                UV_rhs(1,nz,elem) = dt*(UV_rhs(1,nz,elem)+UV_rhsAB(1,1,nz,elem)*ff)/mesh%elem_area(elem)
                UV_rhs(2,nz,elem) = dt*(UV_rhs(2,nz,elem)+UV_rhsAB(1,2,nz,elem)*ff)/mesh%elem_area(elem)
            end do
        end do
    end subroutine compute_vel_rhs

    !==========================================================================
    subroutine momentum_adv_scalar(dynamics, mesh, partit)
        ! Momentum advection on scalar (nodal) control volumes with ALE adaption.
        ! Transcribed VERBATIM from FESOM2 v2.7.3 oce_ale_vel_rhs.F90:335-589. Three
        ! passes:
        !   1. vertical w*du/dz : per node, average elemental UV to full-depth prism
        !      faces (weighted by the node's kite share elem_area*nod_in_elem2D_w), multiply
        !      by w_e, then the layer flux divergence -(wu(nz)-wu(nz+1))/hnode(nz,n)
        !      -> UVnode_rhs.
        !   2. horizontal u*du/dx + v*du/dy : per edge, the edge-normal velocity
        !      un1/un2 from edge_cross_dxdy, scatter +-un*UV into UVnode_rhs at the two
        !      edge nodes (the boundary el2==el1 branch keeps only the el1 contribution).
        !   3. UVnode_rhs *= areasvol_inv ; exchange_nod (1-rank no-op) ; vertice->element
        !      kite-weighted mean (elem_w) ADD into UV_rhsAB(1,1:2,ul:nl-1,el).
        !
        ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern — every operand is already
        ! byte-pinned, so faithful transcription byte-matches like M1.1-M2.3):
        !  - The per-node accumulation over nod_in_elem2D(:,n) (pass 1) and the per-node
        !    edge-scatter (pass 2) are FP-order-sensitive; the 1-rank global element/edge
        !    order matches FESOM2's exactly (same L8/L9 argument as the area gate / FCT
        !    scatter: the `area` gate already pinned nod_in_elem2D order, the geometry
        !    gate pinned edges/edge_tri order).
        !  - Hybrid meshes: the FESOM2 `/(3._WP*hnode)` and vertice->element `/3.0_WP`
        !    are the triangle kite weights 1/3; they are now nod_in_elem2D_w / elem_w (unified
        !    4-slot formulas, triangles equal to roundoff, not bitwise). elem2D_nodes(k,el)
        !    is accessed PER-COMPONENT (k = 1..4; a triangle's slot 4 repeats slot 1 with
        !    weight 0).
        !  - 1-rank only: exchange_nod(UVnode_rhs) is a no-op (lifted at M2.12). The
        !    `nod(1:2) <= nod2D` guards are always true at 1-rank (kept for multi-rank
        !    fidelity). OpenMP locks/ordered are dropped (FESOM3 v1 is serial).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer :: n, nz, el1, el2, nl1, nl2, ul1, ul2, nod(2), el, ed, k, nle, ule
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: un1(1:mesh%nl-1), un2(1:mesh%nl-1)
        real(kind=WP) :: wu(1:mesh%nl),    wv(1:mesh%nl), wa
        real(kind=WP), dimension(:,:,:),   pointer :: UV, UVnode_rhs
        real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB
        real(kind=WP), dimension(:,:),     pointer :: Wvel_e

        UV         => dynamics%uv
        UV_rhsAB   => dynamics%uv_rhsAB
        UVnode_rhs => dynamics%work%uvnode_rhs
        Wvel_e     => dynamics%w_e
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !______________________________________________________________________
        ! 1st. vertical momentum advection component: w*du/dz, w*dv/dz
        do n = 1, nNodO
            nl1 = mesh%nlevels_nod2D(n) - 1
            ul1 = mesh%ulevels_nod2D(n)
            wu(1:nl1+1) = 0._WP
            wv(1:nl1+1) = 0._WP

            ! loop over adjacent elements of vertice n; wa = the node's kite share of el
            do k = 1, mesh%nod_in_elem2D_num(n)
                el  = mesh%nod_in_elem2D(k, n)
                nle = mesh%nlevels(el) - 1
                ule = mesh%ulevels(el)
                wa  = mesh%elem_area(el)*mesh%nod_in_elem2D_w(k, n)

                ! accumulate horizontal velocities at full-depth faces (top/bottom of
                ! prism). The surface face (ule==1) takes the surface-layer velocity;
                ! interior faces average the two adjacent mid-depth-layer velocities.
                if (ule == 1) then
                    wu(ule) = wu(ule) + UV(1,ule,el)*wa
                    wv(ule) = wv(ule) + UV(2,ule,el)*wa
                end if
                wu(ule+1:nle) = wu(ule+1:nle) + 0.5_WP*(UV(1,ule+1:nle,el)+UV(1,ule:nle-1,el))*wa
                wv(ule+1:nle) = wv(ule+1:nle) + 0.5_WP*(UV(2,ule+1:nle,el)+UV(2,ule:nle-1,el))*wa
            end do

            ! multiply face momentum by the (explicit) vertical velocity
            wu(ul1:nl1) = wu(ul1:nl1)*Wvel_e(ul1:nl1,n)
            wv(ul1:nl1) = wv(ul1:nl1)*Wvel_e(ul1:nl1,n)

            ! w*du/dz, w*dv/dz (the kite weights already split each element area)
            do nz = ul1, nl1
                UVnode_rhs(1,nz,n) = - (wu(nz) - wu(nz+1)) / mesh%hnode(nz,n)
                UVnode_rhs(2,nz,n) = - (wv(nz) - wv(nz+1)) / mesh%hnode(nz,n)
            end do

            ! clean checksum: zero the remaining (below-bottom / above-surface) entries
            UVnode_rhs(1:2,nl1+1:mesh%nl-1,n) = 0._WP
            UVnode_rhs(1:2,1:ul1-1        ,n) = 0._WP
        end do

        !______________________________________________________________________
        ! 2nd. horizontal advection component: u*du/dx, v*du/dx & u*dv/dy, v*dv/dy
        do ed = 1, nEdgeO
            nod = mesh%edges(:,ed)
            el1 = mesh%edge_tri(1,ed)
            el2 = mesh%edge_tri(2,ed)
            nl1 = mesh%nlevels(el1) - 1
            ul1 = mesh%ulevels(el1)

            ! edge-normal velocity from el1 centroid towards the edge mid-point
            un1(ul1:nl1) =   UV(2,ul1:nl1,el1)*mesh%edge_cross_dxdy(1,ed) &
                           - UV(1,ul1:nl1,el1)*mesh%edge_cross_dxdy(2,ed)

            if (el2 /= el1) then  ! interior edge: el2 is the other element
                nl2 = mesh%nlevels(el2) - 1
                ul2 = mesh%ulevels(el2)

                un2(ul2:nl2) = - UV(2,ul2:nl2,el2)*mesh%edge_cross_dxdy(3,ed) &
                               + UV(1,ul2:nl2,el2)*mesh%edge_cross_dxdy(4,ed)

                ! zero-fill to combine the two columns over the common level range
                un1(nl1+1:max(nl1,nl2)) = 0._WP
                un2(nl2+1:max(nl1,nl2)) = 0._WP
                un1(1:ul1-1)            = 0._WP
                un2(1:ul2-1)            = 0._WP

                ! first edge node (always owned at 1-rank)
                if (nod(1) <= nNodO) then
                    do nz = min(ul1,ul2), max(nl1,nl2)
                        UVnode_rhs(1,nz,nod(1)) = UVnode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1) + un2(nz)*UV(1,nz,el2)
                        UVnode_rhs(2,nz,nod(1)) = UVnode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1) + un2(nz)*UV(2,nz,el2)
                    end do
                end if
                ! second edge node
                if (nod(2) <= nNodO) then
                    do nz = min(ul1,ul2), max(nl1,nl2)
                        UVnode_rhs(1,nz,nod(2)) = UVnode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1) - un2(nz)*UV(1,nz,el2)
                        UVnode_rhs(2,nz,nod(2)) = UVnode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1) - un2(nz)*UV(2,nz,el2)
                    end do
                end if

            else  ! boundary edge: only el1 contributes
                ! first edge node
                if (nod(1) <= nNodO) then
                    do nz = ul1, nl1
                        UVnode_rhs(1,nz,nod(1)) = UVnode_rhs(1,nz,nod(1)) + un1(nz)*UV(1,nz,el1)
                        UVnode_rhs(2,nz,nod(1)) = UVnode_rhs(2,nz,nod(1)) + un1(nz)*UV(2,nz,el1)
                    end do
                end if
                ! second edge node
                if (nod(2) <= nNodO) then
                    do nz = ul1, nl1
                        UVnode_rhs(1,nz,nod(2)) = UVnode_rhs(1,nz,nod(2)) - un1(nz)*UV(1,nz,el1)
                        UVnode_rhs(2,nz,nod(2)) = UVnode_rhs(2,nz,nod(2)) - un1(nz)*UV(2,nz,el1)
                    end do
                end if
            end if
        end do

        !______________________________________________________________________
        ! 3rd. divide the total nodal advection by the scalar control-volume area.
        !
        ! FESOM3 bottom at vertices: this is a FLUX DIVERGENCE over the control volume
        ! (wu/wv above accumulate UV*elem_area, i.e. a transport), not an average, so the
        ! full prism area is the correct divisor and it is deliberately NOT renormalised
        ! by the wet element area the way tr_xynodes is (oce_ale_tracer.F90). Under the
        ! vertex bottom the CV is the full prism at every wet layer and the lateral walls
        ! run through the velocity points, so a dry adjacent element contributes no flux
        ! and needs no share of the area.
        do n = 1, nNodO
            nl1 = mesh%nlevels_nod2D(n) - 1
            ul1 = mesh%ulevels_nod2D(n)
            UVnode_rhs(1,ul1:nl1,n) = UVnode_rhs(1,ul1:nl1,n)*mesh%areasvol_inv(n)
            UVnode_rhs(2,ul1:nl1,n) = UVnode_rhs(2,ul1:nl1,n)*mesh%areasvol_inv(n)
        end do

        ! M2.12c: share the nodal advection to the halo (FESOM2 :559) — the vertice->
        ! element conversion below reads UVnode_rhs at an owned element's vertices (halo).
        if (is_multirank(partit)) call exchange_nod(UVnode_rhs, partit)

        ! convert nodal advection vertice -> element (kite-weighted mean) and ADD into UV_rhsAB
        do el = 1, nElemO
            nl1 = mesh%nlevels(el) - 1
            ul1 = mesh%ulevels(el)
            UV_rhsAB(1,1:2,ul1:nl1,el) = UV_rhsAB(1,1:2,ul1:nl1,el) &
                    + mesh%elem_area(el)*(mesh%elem_w(1,el)*UVnode_rhs(1:2,ul1:nl1,mesh%elem2D_nodes(1,el)) &
                    + mesh%elem_w(2,el)*UVnode_rhs(1:2,ul1:nl1,mesh%elem2D_nodes(2,el)) &
                    + mesh%elem_w(3,el)*UVnode_rhs(1:2,ul1:nl1,mesh%elem2D_nodes(3,el)) &
                    + mesh%elem_w(4,el)*UVnode_rhs(1:2,ul1:nl1,mesh%elem2D_nodes(4,el)))
        end do
    end subroutine momentum_adv_scalar

end module oce_dyn_velrhs
