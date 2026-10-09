module oce_dyn_vinv
    ! Vector-invariant momentum advection (momadv_opt==1):
    !
    !     (curl u + f) x u  +  grad(u^2/2)  +  w du/dz
    !
    ! PROVENANCE. Transcribed from /home/a/a270029/qq/oce_vinv_mom_adv.F90 (routines
    ! relative_vorticity + v_inv_mom_adv). There is NO FESOM2 oracle for this path: the
    ! v2.7.3 branch exists but ABORTS --
    !   oce_ale_vel_rhs.F90:268  'in moment not adapted mom_adv advection typ for ALE'
    ! so FESOM2 never adapted the vector-invariant form to ALE and it was never tested
    ! there. Correctness therefore rests on analytic unit tests (test/test_vinv.F90) plus
    ! the no-work invariant u.[(f+zeta) x u] == 0 checked in fesom_conserve, NOT on a
    ! reference dump. See docs/plans/2026-10-01-momadv-vector-invariant.md.
    !
    ! INTEGRATION. momentum_adv_vinv is a pure ADDITION into UV_rhsAB, exactly like
    ! momentum_adv_scalar. compute_vel_rhs step (3) already wrote f x u into that slot, so
    ! this routine adds  zeta x u - grad(KE) + w du/dz, which completes (f+zeta) x u.
    ! compute_vel_rhs steps (1)/(2)/(4) are untouched and momadv_opt==2 is unreachable
    ! from here, so that path stays bit-identical.
    !
    ! FORCED DEVIATIONS from qq (it is pre-ALE research code):
    !   qq                                   -> here                  why
    !   dz = zbar(1:nl-1)-zbar(2:nl) global  -> mesh%helem(nz,elem)   ALE: per-elem/step
    !   area(nz,n)                           -> mesh%area(n)          collapsed to 1-D
    !   w_cv(:,elem)                         -> mesh%elem_w(:,elem)   kite weights
    !   exchange_nod3D(x)                    -> exchange_nod(x,partit) halo API
    !   loops from 1                         -> ulevels..nlevels-1    cavity bounds
    !   Wvel                                 -> dynamics%w_e          match momadv_opt==2
    !   whole-RHS assembly                   -> advection blocks only rest is in velrhs
    ! qq's `else` branch (explicit vertical viscosity: wind stress + a hard-coded
    ! friction=0.005 "Soufflet" replacing the real bottom drag) is NOT transcribed; we
    ! take its i_vert_visc=.true. branch, whose own comment reads "Do only advection".
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_dyn,         only: t_dyn
    use mod_partit,      only: t_partit
    use mod_param_phys,  only: rvo_upwind
    use mod_constants,   only: r_earth
    use mod_mesh_rotate, only: trim_cyclic
    use mod_mesh_areas,  only: elem_center
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod, exchange_elem
    implicit none
    private
    public :: relative_vorticity, momentum_adv_vinv

contains

    subroutine relative_vorticity(dynamics, mesh, partit)
        ! Block A. zeta = curl(u) at nodes, as the circulation integral around the median-
        ! dual contour (qq oce_vinv_mom_adv.F90:5-61).
        !
        ! Per edge the contour crosses the two adjacent element centres, so the segment
        ! contribution is  d . u  with d = edge_cross_dxdy (METRES here, as qq assumes),
        ! taken with opposite signs for the two edge nodes. qq's THREE level ranges are
        ! kept verbatim: both elements wet, then el1 only, then el2 only. Those one-sided
        ! ranges are what makes the operator correct at a bathymetry step -- below the
        ! shallower element only the deeper one contributes, i.e. the dry element's velocity
        ! is treated as exactly zero while area(n) keeps its FULL weight. That is precisely
        ! the bottom-at-vertices rule (mod_mesh.F90), so area(nz,n) -> area(n) and the
        ! one-sided ranges are mutually consistent, not two independent choices.
        !
        ! Halo: every edge incident to an owned node is owned, and both elements of an
        ! owned edge are owned (LESSONS L31), so the loop over 1..nEdgeO is complete for
        ! owned rows and UV needs NO element exchange. One exchange_nod publishes zeta to
        ! the halo, which the per-element (f+zeta) average then reads.
        !
        ! Cavity: the TOP asymmetry (one adjacent element shallower at the surface) is NOT
        ! handled -- the ranges below start at the deeper ulevels. Untested without a
        ! cavity mesh, like the rest of the tree's cavity paths.
        type(t_mesh),   intent(in)              :: mesh
        type(t_dyn),    intent(inout), target   :: dynamics
        type(t_partit), intent(in),    optional :: partit
        integer       :: n, nz, edge, el(2), enodes(2), ul1, ul2, nl1, nl2, nboth, ultop
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dx1, dy1, dx2, dy2, c1
        real(kind=WP), dimension(:,:),   pointer :: vort
        real(kind=WP), dimension(:,:,:), pointer :: UV

        vort => dynamics%work%vorticity
        UV   => dynamics%uv
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! Zero the FULL local range -- owned AND halo, all levels. The scatter below writes
        ! into halo rows, and -init=zero does not cover allocatables, so accumulating onto
        ! fresh heap could trap under -fpe0 before the exchange made the value irrelevant.
        vort(:, 1:nNodL) = 0.0_WP

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            ul1    = mesh%ulevels(el(1))
            nl1    = mesh%nlevels(el(1)) - 1
            dx1    = real(mesh%edge_cross_dxdy(1, edge), WP)
            dy1    = real(mesh%edge_cross_dxdy(2, edge), WP)
            ul2    = ul1
            nl2    = 0
            dx2    = 0.0_WP
            dy2    = 0.0_WP
            if (el(2) /= el(1)) then
                dx2 = real(mesh%edge_cross_dxdy(3, edge), WP)
                dy2 = real(mesh%edge_cross_dxdy(4, edge), WP)
                ul2 = mesh%ulevels(el(2))
                nl2 = mesh%nlevels(el(2)) - 1
            end if
            nboth = min(nl1, nl2)
            ultop = max(ul1, ul2)

            ! both elements wet
            do nz = ultop, nboth
                c1 = dx1*UV(1,nz,el(1)) + dy1*UV(2,nz,el(1)) &
                   - dx2*UV(1,nz,el(2)) - dy2*UV(2,nz,el(2))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
            ! el1 only (el2 dry or absent): its velocity contributes zero
            do nz = max(nboth+1, ul1), nl1
                c1 = dx1*UV(1,nz,el(1)) + dy1*UV(2,nz,el(1))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
            ! el2 only
            do nz = max(nboth+1, ul2), nl2
                c1 = -dx2*UV(1,nz,el(2)) - dy2*UV(2,nz,el(2))
                vort(nz, enodes(1)) = vort(nz, enodes(1)) + c1
                vort(nz, enodes(2)) = vort(nz, enodes(2)) - c1
            end do
        end do

        ! circulation -> vorticity.
        !   area(n): the full median-dual area, sum over the node's elements of the
        !            kite shares elem_area*nod_in_elem2D_w (DEVIATION: qq divides by
        !            area(nz,n)).
        !   SIGN   : qq's combination c1 = d1.u(el1) - d2.u(el2), added to enodes(1) and
        !            subtracted from enodes(2), traverses the dual contour in the sense that
        !            matches OUR edge_tri left/right convention -- no flip needed.
        !            VERIFIED: test_vinv V2 gives zeta == -alpha for u = alpha*y to
        !            relative 5e-13 over 81729 node-levels, i.e. exactly
        !            curl(u) = dv/dx - du/dy. That same exactness also PROVES the
        !            area(nz,n) -> area(n) deviation is right: zeta = circulation/area
        !            only reproduces -alpha if area(n) equals the sum of the kite shares
        !            over the node's WET adjacent elements. On quads the exactness holds
        !            for parallelograms only (general quads: first order -- test_vinv header).
        do n = 1, nNodO
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
                vort(nz,n) = vort(nz,n)/real(mesh%area(n), WP)
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(vort, partit)
    end subroutine relative_vorticity

    subroutine momentum_adv_vinv(dynamics, mesh, partit)
        ! Blocks B and C. ADDS  zeta x u - grad(KE) - w du/dz  into UV_rhsAB(1,1:2,:).
        ! compute_vel_rhs step (3) already put f x u there, so the sum is (f+zeta) x u.
        type(t_mesh),   intent(in)              :: mesh
        type(t_dyn),    intent(inout), target   :: dynamics
        type(t_partit), intent(in),    optional :: partit
        integer       :: n, nz, elem, ed, elnodes(MAX_NV), ul, nl1, k, k2, nb, nElemF
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: ea, pre(MAX_NV), Fx, Fy, zb, w, umean, vmean, h1, h2, da, hinv, ke_e
        real(kind=WP) :: xv(MAX_NV), yv(MAX_NV), tx, ty, px, py, dref, nxk(MAX_NV), nyk(MAX_NV)
        real(kind=WP) :: d, wk, sumw, zup, xc, yc, dx
        real(kind=WP) :: we(MAX_NV)
        real(kind=WP) :: uvert(2, mesh%nl)
        real(kind=WP), allocatable :: KE(:,:), omega_e(:,:)
        real(kind=WP), dimension(:,:),     pointer :: vort, Wv
        real(kind=WP), dimension(:,:,:),   pointer :: UV
        real(kind=WP), dimension(:,:,:,:), pointer :: UV_rhsAB

        vort     => dynamics%work%vorticity
        UV       => dynamics%uv
        UV_rhsAB => dynamics%uv_rhsAB
        Wv       => dynamics%w_e
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        call relative_vorticity(dynamics, mesh, partit)

        !___________________________________________________________________
        ! rvo_upwind: element-centre vorticity for the upwind face reconstruction.
        ! omega_e(nz,e) = the SAME kite-weighted vertex mean Block B2 uses as its face value; one
        ! STANDARD exchange_elem (eDim ring) makes every halo neighbour id readable, so a
        ! neighbour's VALUE is available without ever needing its vertex list
        ! (elem2D_nodes is owned-only). The eDim ring suffices because an edge neighbour
        ! of an owned element is by construction in com_elem2D, never in eXDim (see
        ! build_elem_adjacency; test_vinv A5 guards it). A neighbour that is absent or
        ! DRY at the level is skipped in the blend, so the 0 initial value is never read.
        ! Local allocatable: nothing persists, nothing reaches the restart, and at
        ! rvo_upwind = 0 none of this executes -- the default path is bit-identical by
        ! construction, not by tolerance.
        if (rvo_upwind > 0.0_WP) then
            nElemF = size(UV, 3)
            allocate(omega_e(mesh%nl-1, nElemF))
            omega_e = 0.0_WP
            do elem = 1, nElemO
                elnodes = mesh%elem2D_nodes(:, elem)
                we      = real(mesh%elem_w(:, elem), WP)
                do nz = mesh%ulevels(elem), mesh%nlevels(elem)-1
                    omega_e(nz, elem) = sum(we*vort(nz,elnodes))
                end do
            end do
            if (is_multirank(partit)) call exchange_elem(omega_e, partit)
        end if

        !___________________________________________________________________
        ! Block B1: kinetic energy at nodes (qq :78-115), the node gather with the kite
        ! shares elem_w*elem_area. Scalar scatter over the 4 slots: a triangle's slot 4
        ! repeats slot 1 (with weight 0), so a vector-subscript update would be many-one.
        allocate(KE(mesh%nl-1, nNodL))
        KE = 0.0_WP                       ! full local range: the scatter writes halo rows
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            we      = real(mesh%elem_w(:, elem), WP)
            ea      = real(mesh%elem_area(elem), WP)
            do nz = mesh%ulevels(elem), mesh%nlevels(elem)-1
                ke_e = 0.5_WP*(UV(1,nz,elem)**2 + UV(2,nz,elem)**2)
                do k = 1, MAX_NV
                    KE(nz, elnodes(k)) = KE(nz, elnodes(k)) + ke_e*we(k)*ea
                end do
            end do
        end do
        do n = 1, nNodO
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
                KE(nz,n) = KE(nz,n)/real(mesh%area(n), WP)
            end do
        end do
        ! Lateral-wall BC (qq :106-112): KE = 0 at boundary nodes. edge_tri(2,ed) ==
        ! edge_tri(1,ed) is the partition-robust boundary test; qq's `myList_edge2D(n) <= edge2D_in` form
        ! needs the GLOBAL edge id and is silently wrong on local indices at npes>1.
        ! NOTE this makes grad(KE) nonzero along every coastline even for a uniform flow --
        ! it is qq's idealised-channel choice, transcribed as-is.
        do ed = 1, nEdgeO
            if (mesh%edge_tri(2, ed) /= mesh%edge_tri(1, ed)) cycle
            KE(:, mesh%edges(1, ed)) = 0.0_WP
            KE(:, mesh%edges(2, ed)) = 0.0_WP
        end do
        if (is_multirank(partit)) call exchange_nod(KE, partit)

        !___________________________________________________________________
        ! Block B2: zeta x u  and  -grad(KE)   (qq :157-179)
        ! DEVIATION: qq averages coriolis_node(elnodes) TOGETHER with vorticity using the
        ! same w_cv weights. Here f comes from compute_vel_rhs step (3) as the element-
        ! centre mesh%coriolis(elem), so f and zeta are weighted differently:
        ! coriolis(elem) + sum(w*zeta) rather than sum(w*(coriolis_node + zeta)). Deliberate:
        ! it keeps momadv_opt==2 bit-identical and uses FESOM2's standard Coriolis.
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            we      = real(mesh%elem_w(:, elem), WP)
            ea      = real(mesh%elem_area(elem), WP)

            if (rvo_upwind > 0.0_WP) then
                ! The INWARD side normals, rebuilt from the element's OWN vertices:
                ! per-edge geometry (edge_dxdy) exists for OWNED edges only and an owned
                ! element can carry one HALO edge, so edge arrays are unusable here. Side k
                ! connects elnodes(k) and elnodes(k+1) (slot 4 -> slot 1; a triangle's side
                ! 4 has zero length, so its normal is 0 and it never gets weight); of the two
                ! perpendiculars of t = (dlon*elem_cos*r_earth, dlat*r_earth), the INWARD one
                ! points toward the element centroid (mod_mesh_areas elem_center; its
                ! longitude wrapped to the side's start vertex) -- a local,
                ! convention-free orientation (no edge_tri left/right branch; proven
                ! geometrically by test_vinv V10: the slot-k neighbour always lies on the
                ! outward side).
                do k = 1, MAX_NV
                    xv(k) = real(mesh%coord_nod2D(1, elnodes(k)), WP)
                    yv(k) = real(mesh%coord_nod2D(2, elnodes(k)), WP)
                end do
                do k = 2, MAX_NV
                    tx = xv(k) - xv(1); call trim_cyclic(tx); xv(k) = xv(1) + tx
                end do
                call elem_center(mesh, elem, xc, yc)
                do k = 1, MAX_NV
                    k2 = mod(k,MAX_NV)+1
                    tx = (xv(k2)-xv(k))*real(mesh%elem_cos(elem),WP)*r_earth
                    ty = (yv(k2)-yv(k))*r_earth
                    px =  ty; py = -tx
                    dx = xc - xv(k); call trim_cyclic(dx)
                    dref = px*dx*real(mesh%elem_cos(elem),WP)*r_earth &
                         + py*(yc-yv(k))*r_earth
                    if (dref < 0.0_WP) then
                        px = -px; py = -py
                    end if
                    nxk(k) = px; nyk(k) = py
                end do
            end if

            do nz = mesh%ulevels(elem), mesh%nlevels(elem)-1
                pre = -KE(nz, elnodes)
                Fx  = sum(real(mesh%gradient_sca(1:4, elem), WP)*pre)
                Fy  = sum(real(mesh%gradient_sca(5:8, elem), WP)*pre)
                zb  = sum(we*vort(nz,elnodes))

                ! rvo_upwind face blend:  zb += rvo*(zb_upwind - zb).
                !  - d = u . n_in  (own-element velocity): d > 0 means flow ENTERS through
                !    edge k, so that neighbour is upwind; w = d+|d| zeroes outflow edges
                !    and the normalisation weights a stronger inflow neighbour more.
                !  - dry or boundary neighbours are EXCLUDED, not zero-valued: zb_upwind
                !    is a pointwise VALUE estimate, not a flux (estimator rule, in
                !    contrast to the zero-with-full-weight rule for flux-forming sums).
                !    Everything excluded, or no inflow (u=0): sumw = 0 -> blend inert.
                !  - convex combination -> max principle: no new vorticity extrema.
                !  - energy-neutral for ANY zb: u.[(f+zeta) x u] == 0 identically, so the
                !    blend acts on the vorticity/enstrophy dynamics only.
                if (rvo_upwind > 0.0_WP) then
                    sumw = 0.0_WP
                    zup  = 0.0_WP
                    do k = 1, MAX_NV
                        nb = mesh%elem_neighbors(k, elem)
                        if (nb == elem) cycle
                        if (nz < mesh%ulevels(nb) .or. nz > mesh%nlevels(nb)-1) cycle
                        d  = UV(1,nz,elem)*nxk(k) + UV(2,nz,elem)*nyk(k)
                        wk = d + abs(d)
                        sumw = sumw + wk
                        zup  = zup  + wk*omega_e(nz, nb)
                    end do
                    if (sumw > 0.0_WP) zb = zb + rvo_upwind*(zup/sumw - zb)
                end if

                UV_rhsAB(1,1,nz,elem) = UV_rhsAB(1,1,nz,elem) + ( UV(2,nz,elem)*zb + Fx)*ea
                UV_rhsAB(1,2,nz,elem) = UV_rhsAB(1,2,nz,elem) + (-UV(1,nz,elem)*zb + Fy)*ea
            end do
        end do
        deallocate(KE)
        if (allocated(omega_e)) deallocate(omega_e)

        !___________________________________________________________________
        ! Block C: -w du/dz, as d(wu)/dz - u dw/dz  (qq :180-246, i_vert_visc branch)
        ! NO Av term: FESOM3 solves vertical viscosity implicitly in impl_vert_visc_ale,
        ! so transcribing qq's Av*(du/dz) here would apply it twice.
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            we      = real(mesh%elem_w(:, elem), WP)
            ul      = mesh%ulevels(elem)
            nl1     = mesh%nlevels(elem)-1
            ea      = real(mesh%elem_area(elem), WP)
            if (nl1 < ul) cycle

            w = sum(we*Wv(ul, elnodes))
            uvert(1,ul) = -w*UV(1,ul,elem)
            uvert(2,ul) = -w*UV(2,ul,elem)

            do nz = ul+1, nl1
                w  = sum(we*Wv(nz, elnodes))
                h1 = real(mesh%helem(nz-1,elem), WP)      ! qq's dz(nz-1), now ALE
                h2 = real(mesh%helem(nz,  elem), WP)      ! qq's dz(nz)
                umean = (UV(1,nz-1,elem)*h2 + UV(1,nz,elem)*h1)/(h1+h2)
                vmean = (UV(2,nz-1,elem)*h2 + UV(2,nz,elem)*h1)/(h1+h2)
                uvert(1,nz) = -umean*w
                uvert(2,nz) = -vmean*w
            end do

            ! DEVIATION FROM qq (a FIX, not a forced change): qq sets uvert(nl1+1) = 0.
            ! Under bottom-at-vertices that breaks the telescoping, because da below reads
            ! w_e(nlevels(e)) which is NOT zero: vert_vel_ale zeroes Wvel only at
            ! nlevels_nod2D(n) (oce_ale.F90:373-392) while nlevels(e) = minval over the
            ! element's nodes. Every element beside a deeper vertex column would then gain
            ! a spurious -w_bot*U/helem in its bottom layer. Treating the bottom face like
            ! the surface restores exact cancellation for a uniform u and keeps the real
            ! bottom flux. test_vinv V6 is the net for this.
            w = sum(we*Wv(nl1+1, elnodes))
            uvert(1,nl1+1) = -w*UV(1,nl1,elem)
            uvert(2,nl1+1) = -w*UV(2,nl1,elem)

            do nz = ul, nl1
                da   = sum(we*(Wv(nz,elnodes) - Wv(nz+1,elnodes)))
                hinv = 1.0_WP/real(mesh%helem(nz,elem), WP)
                UV_rhsAB(1,1,nz,elem) = UV_rhsAB(1,1,nz,elem) &
                    + (uvert(1,nz) - uvert(1,nz+1) + da*UV(1,nz,elem))*ea*hinv
                UV_rhsAB(1,2,nz,elem) = UV_rhsAB(1,2,nz,elem) &
                    + (uvert(2,nz) - uvert(2,nz+1) + da*UV(2,nz,elem))*ea*hinv
            end do
        end do
    end subroutine momentum_adv_vinv

end module oce_dyn_vinv
