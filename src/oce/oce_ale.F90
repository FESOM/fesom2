module oce_ale
    ! M2.7 — ALE (linfs) velocity / SSH / thickness-W update: the post-CG tail of the
    ! FESOM2 timestep (oce_ale.F90:3946-4081 + the inline eta_n blend at :3973-3977).
    !
    ! After the free-surface CG solve (M2.6 -> d_eta) the timestep does, in order:
    !   1. update_vel        — u^(n+1) = u* + UV_rhs + [-g*theta*dt*grad(d_eta)]
    !                          (FESOM2 oce_dyn.F90:88-173)
    !   2. compute_hbar_ale  — hbar_old=hbar; hbar += dt/areasvol * div(UV); dhe
    !                          (FESOM2 oce_ale.F90:2165-2307)
    !   3. update_eta_n      — eta_n = alpha*hbar + (1-alpha)*hbar_old (inline :3973-3977)
    !   4. vert_vel_ale      — W = -cumsum(div(UV*h))/area; linfs leaves hnode_new=hnode;
    !                          then compute_CFLz + compute_Wvel_split
    !                          (FESOM2 oce_ale.F90:2323-2884, the which_ale='linfs' path)
    !
    ! Scope reductions (all the gated-config branches the pi/reduced-M2 namelist takes):
    !  - which_ale='linfs' -> the zlevel/zstar thickness redistribution in vert_vel_ale is
    !    NOT taken, so hnode_new stays = hnode (its init value) and there is no surface
    !    Wvel/hnode correction; the compute_hbar_ale water_flux term ((1-alpha) / non-linfs)
    !    vanishes. The full-free-surface thickness evolution gets its own later gate.
    !  - Fer_GM=.false. -> the fer_UV/fer_Wvel transport in vert_vel_ale is dropped.
    !  - ldiag_ke=.false. -> the update_vel ke_* energy diagnostics are dropped (no ke_*
    !    in t_dyn, as in compute_vel_rhs / impl_vert_visc_ale).
    !  - 1-rank: myDim_nod2D+eDim_nod2D == myDim_nod2D == nod2D, myDim_edge2D == edge2D;
    !    every exchange_nod/exchange_elem is a no-op (lifted at M2.12). The OpenMP-off
    !    oracle (L16/L17/L19) compiles the omp locks/ORDERED out -> serial edge order, which
    !    matches these serial loops, so the divergence/W scatters byte-match by L9.
    !
    ! Byte-gate (M2.7, tools/run_pressure_gate.sh): every operand is already pinned —
    ! d_eta (M2.6), UV/UV_rhs (the post-TDMA M2.5 value), gradient_sca/edge_cross_dxdy/
    ! helem/areasvol/area/levels (geometry + M2.5), g/density_0/dt/alpha/theta — so the
    ! whole chain is max|delta|=0 by transitivity (L9). hbar is prescribed identically.
    use mod_precision,  only: WP, MP
    use mod_constants,  only: g
    use mod_param_phys, only: alpha, theta, Fer_GM
    use mod_mesh,       only: t_mesh, MAX_NV
    use mod_dyn,        only: t_dyn
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,       only: exchange_nod, exchange_elem, exchange_elem_full, exchange_elem_shared
    use mod_timer,      only: timer_start, timer_stop, TMR_OCE_SHSYNC
    use mod_config,     only: which_ALE          ! M6 ALE: 'linfs'/'zlevel'/'zstar'
    use oce_wsplit,     only: wsplit_implicit_fraction   ! the C^1 w_e/w_i split function
    implicit none
    private
    public :: update_vel, compute_hbar_ale, update_eta_n, vert_vel_ale
    ! M2.9b (step assembly): compute_vel_nodes is the REAL dyn%uvnode source (prescribed
    ! at M2.8); update_thickness_ale is the linfs thickness commit (a no-op for which_ale=
    ! 'linfs'). Both join the faithful oce_timestep sequence assembled in mod_step_oce.
    public :: compute_vel_nodes, update_thickness_ale
    ! The vertical-CFL / w-split pair at the end of vert_vel_ale, exported for
    ! test/test_wsplit.F90 (parts S and X drive them on prescribed w) and
    ! test/test_wimpl_tra.F90 (C2-C6 build their split with them).
    public :: compute_CFLz, compute_Wvel_split

contains

    !===========================================================================
    subroutine compute_vel_nodes(dynamics, mesh, partit)
        ! FESOM2 oce_dyn.F90:177-226. Horizontal velocity averaged from elements to
        ! nodes (area-weighted over the node's element star) -> dyn%uvnode, the input
        ! PP mixing (M2.8 oce_mixing_pp) reads for the vertical shear. Called at the top
        ! of the ocean step (FESOM2 fesom_module.F90:654, before oce_timestep_ale).
        ! M2.12c: optional partit -> owned node loop (FESOM2 :202 do n=1,myDim_nod2D) +
        ! the trailing exchange_nod(UVnode). The per-node accumulation over nod_in_elem2D
        ! is order-sensitive but byte-pinned by the geometry gate (L9), so it matches
        ! FESOM2's local element order on the same partition.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: n, nz, k, elem, nln, uln, nle, ule
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: tx, ty, tvol
        real(kind=WP), dimension(:,:,:), pointer :: UV, UVnode

        UV     => dynamics%uv
        UVnode => dynamics%uvnode
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! FESOM3 BOTTOM AT VERTICES. An adjacent element that is dry at level nz is NOT
        ! dropped from the average: its velocity there is genuinely ZERO -- that is what
        ! "velocities touching topography are zero" means -- so it contributes zero while
        ! still carrying its FULL area weight. The denominator is therefore the full
        ! median-dual area area(n) = sum_k nod_in_elem2D_w(k,n)*elem_area (the kite share of
        ! every adjacent element), not the area of the elements that happen to be wet.
        ! This is the same statement that makes the scalar column's area depth-independent,
        ! and it is why area/areasvol are 1-D.
        !
        ! Consequence, and it is intended rather than a side effect: a uniform velocity
        ! field does NOT node-average to the same value at every level -- it decreases where
        ! the wet fraction does. That is the real damping of the cell-mean velocity by
        ! topography, not an artefact; UVnode feeds the vertical shear used by PP/KPP/TKE.
        do n = 1, nNodO
            uln = mesh%ulevels_nod2D(n)
            nln = mesh%nlevels_nod2D(n)
            tvol = real(mesh%area(n), WP)             ! full median-dual area at this vertex
            do nz = uln, nln-1
                tx = 0.0_WP; ty = 0.0_WP
                do k = 1, mesh%nod_in_elem2D_num(n)
                    elem = mesh%nod_in_elem2D(k, n)
                    ule  = mesh%ulevels(elem)
                    nle  = mesh%nlevels(elem)
                    if (nle-1 < nz .or. nz < ule) cycle   ! dry: contributes 0, weight kept
                    tx   = tx   + UV(1,nz,elem)*mesh%elem_area(elem)*mesh%nod_in_elem2D_w(k,n)
                    ty   = ty   + UV(2,nz,elem)*mesh%elem_area(elem)*mesh%nod_in_elem2D_w(k,n)
                end do
                UVnode(1,nz,n) = tx/tvol
                UVnode(2,nz,n) = ty/tvol
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(UVnode, partit)   ! FESOM2 :225
    end subroutine compute_vel_nodes

    !===========================================================================
    subroutine update_vel(dynamics, mesh, dt, partit)
        ! FESOM2 oce_dyn.F90:88-173. New horizontal velocity after the SSH solve:
        !   u^(n+1) = u* + UV_rhs + [-g*theta*dt*grad(d_eta)]
        ! The SSH-gradient correction is a gradient_sca contraction of -g*theta*dt*d_eta
        ! over the element's 4 vertex slots (same shape as the M2.3 vel_rhs SSH-gradient
        ! and the M2.2 PGF; a triangle's slot 4 has coefficient 0). UV_rhs here is the
        ! post-TDMA value (M2.5; the CG solve does not touch it).
        ! M2.12c-3: optional partit -> OWNED element loop (FESOM2 :105 do elem=1,
        ! myDim_elem2D), then exchange_elem(UV) (:172). d_eta is read at the owned element's
        ! vertices (within owned+halo, invariant iii) — valid from the CG's exchange_nod(x).
        ! The UV halo is refreshed owner->halo for the NEXT step (compute_vel_nodes/advection);
        ! the owned UV the post-update kernels read is computed directly here. The rank-3
        ! (2,nl-1,elem) UV exchanges over the eDim com_elem2D (exchange_elem_blk_r) — EXACTLY
        ! FESOM2's exchange_elem(UV), refreshing only the eDim halo the owned-edge viscosity
        ! reads. (The earlier com_elem2D_full superset corrupted the CORE2 multi-step halo —
        ! see the production-validation lessons; eDim is the faithful + correct choice.)
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: elem, nz, nzmin, nzmax, elnodes(MAX_NV)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: eta(MAX_NV), Fx, Fy
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:),     pointer :: d_eta

        UV     => dynamics%uv
        UV_rhs => dynamics%uv_rhs
        d_eta  => dynamics%d_eta
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            eta = -g*theta*dt*d_eta(elnodes)
            Fx  = sum(mesh%gradient_sca(1:4, elem)*eta)
            Fy  = sum(mesh%gradient_sca(5:8, elem)*eta)
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                UV(1,nz,elem) = UV(1,nz,elem) + UV_rhs(1,nz,elem) + Fx
                UV(2,nz,elem) = UV(2,nz,elem) + UV_rhs(2,nz,elem) + Fy
            end do
        end do
        if (is_multirank(partit)) then
            ! One owner per shared element (src/infra/mod_shared_elem.F90): a border element
            ! is owned by every rank that owns one of its vertices, and each owner has just
            ! updated its own copy. The copies differ by roundoff (order-dependent sums in the
            ! viscosity filter), which the rvo_upwind blend can turn into O(1); the ranks
            ! then compute different fluxes across their common edges and ssh_rhs_old no
            ! longer sums to zero (linfs mass leak). Every other owner takes the designated
            ! owner's copy -- what a restart read does -- BEFORE the halo exchange, so the
            ! eDim halo copies carry the synced value too. FESOM3_SHARED_SYNC=0 skips it.
            ! (Nested ifs: partit is optional and Fortran does not short-circuit.)
            if (partit%sync_shared_elem) then
                call timer_start(TMR_OCE_SHSYNC)
                call exchange_elem_shared(UV, partit)
                call timer_stop(TMR_OCE_SHSYNC)
            end if
            call exchange_elem(UV, partit)                         ! FESOM2 :172 (eDim com_elem2D)
        end if
    end subroutine update_vel

    !===========================================================================
    subroutine compute_hbar_ale(dynamics, mesh, dt, partit, water_flux)
        ! FESOM2 oce_ale.F90:2165-2307. Advance hbar (the elevation on semi-integer
        ! timesteps) by the depth-integrated horizontal divergence of the UPDATED UV:
        !   ssh_rhs_old = div_h( int(UV) dz )        (edge-scatter, NO alpha here)
        !   hbar_old    = hbar
        !   hbar        = hbar_old + ssh_rhs_old * dt / areasvol
        !   dhe         = sum(elem_w*(hbar - hbar_old)) (per element; feeds update_stiff,
        !                 which linfs never calls — gated anyway)
        ! Same edge-divergence structure as M2.6 compute_ssh_rhs_ale (same edge order,
        ! gated). linfs skips the water_flux term + its exchange_nod. Writes mesh%hbar/
        ! hbar_old/dhe and dynamics%ssh_rhs_old.
        ! M2.12c-3: optional partit -> FESOM2 :2197-2306 bounds: ssh_rhs_old zero owned+halo,
        ! the edge-divergence over OWNED edges (owned-node ssh_rhs_old complete, invariant i),
        ! hbar_old=hbar owned+halo, the hbar update OWNED, then exchange_nod(hbar) (:2293,
        ! ALWAYS — outside the linfs water_flux guard) so the dhe loop over OWNED elements can
        ! read hbar at the element's vertices (within owned+halo). The linfs water_flux term +
        ! its exchange_nod(ssh_rhs_old) stay skipped.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(inout), target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        real(kind=WP),  intent(in), optional :: water_flux(:)   ! M6a-3: non-linfs freshwater (zstar)
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), elem, elnodes(MAX_NV), n, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2
        real(kind=WP), dimension(:,:,:), pointer :: UV
        real(kind=WP), dimension(:),     pointer :: ssh_rhs_old

        UV          => dynamics%uv
        ssh_rhs_old => dynamics%ssh_rhs_old
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodL
            ssh_rhs_old(n) = 0.0_WP
        end do

        do ed = 1, nEdgeO
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            !__________________________________________________________________
            ! depth integral div(int(U_n)dz) for el(1)
            c1      = 0.0_WP
            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmin, nzmax
                c1 = c1 + (UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
            end do

            !__________________________________________________________________
            ! ... and for el(2), unless ed is a boundary edge
            c2 = 0.0_WP
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmin, nzmax
                    c2 = c2 - (UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
                end do
            end if

            ssh_rhs_old(enodes(1)) = ssh_rhs_old(enodes(1)) + (c1+c2)
            ssh_rhs_old(enodes(2)) = ssh_rhs_old(enodes(2)) - (c1+c2)
        end do

        ! M6a-3: non-linfs (zstar) freshwater flux enters the SSH balance — the real
        ! water_flux (positive out of the ocean) is subtracted from the divergence before
        ! hbar is advanced (FESOM2 oce_ale.F90:2262-2271). linfs skips it (virtual salt
        ! instead). water_flux*areasvol so it matches the edge-divergence units (un-divided
        ! by area until the hbar update below).
        if (trim(which_ALE)/='linfs' .and. present(water_flux)) then
            do n = 1, nNodO
                if (mesh%ulevels_nod2D(n) > 1) cycle
                ssh_rhs_old(n) = ssh_rhs_old(n) - water_flux(n)*mesh%areasvol(n)
            end do
            if (is_multirank(partit)) call exchange_nod(ssh_rhs_old, partit)   ! FESOM2 :2270
        end if

        do n = 1, nNodL
            mesh%hbar_old(n) = mesh%hbar(n)
        end do

        do n = 1, nNodO
            if (mesh%ulevels_nod2D(n) > 1) cycle          ! cavity node: hbar == hbar_old
            mesh%hbar(n) = mesh%hbar_old(n) + ssh_rhs_old(n)*dt/mesh%areasvol(n)
        end do
        if (is_multirank(partit)) call exchange_nod(mesh%hbar, partit)   ! FESOM2 :2293

        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            if (mesh%ulevels(elem) > 1) then
                mesh%dhe(elem) = 0.0_WP
            else
                mesh%dhe(elem) = sum(mesh%elem_w(:, elem)*(mesh%hbar(elnodes) - mesh%hbar_old(elnodes)))
            end if
        end do
    end subroutine compute_hbar_ale

    !===========================================================================
    subroutine update_eta_n(dynamics, mesh, partit)
        ! FESOM2 oce_ale.F90:3973-3977 (inline in the step). Current dynamic elevation:
        !   eta_n = alpha*hbar + (1-alpha)*hbar_old   (only where ulevels_nod2D==1, i.e.
        ! no rigid-lid cavity). On pi alpha=1 -> eta_n = hbar (the new elevation).
        ! M2.12c-3: optional partit -> the FESOM2 owned+HALO node loop (:3974 do node=1,
        ! myDim_nod2D+eDim_nod2D); both hbar (compute_hbar_ale's exchange_nod) and hbar_old
        ! (its owned+halo copy) are halo-valid, so eta_n is computed at owned+halo with no
        ! exchange (matches FESOM2 — no exchange_nod(eta_n) after this loop).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer :: node, nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do node = 1, nNodL
            if (mesh%ulevels_nod2D(node) == 1) &
                dynamics%eta_n(node) = alpha*mesh%hbar(node) + (1.0_WP-alpha)*mesh%hbar_old(node)
        end do
    end subroutine update_eta_n

    !===========================================================================
    subroutine vert_vel_ale(dynamics, mesh, dt, partit, water_flux)
        ! FESOM2 oce_ale.F90:2323-2884, the which_ale='linfs' path. Vertical velocity from
        ! the layer-thickness divergence:
        !   Wvel(:,n)  = 0
        !   edge loop  -> Wvel = div(UV*h) per layer (edge-scatter, both signs)
        !   cumsum up  -> Wvel(nz) = Wvel(nz) + Wvel(nz+1)        (W_k = W_{k+1} - div(h_k u_k))
        !   /area      -> physical m/s
        ! linfs: dh/dt=0, so NO zlevel/zstar surface correction and hnode_new stays = hnode.
        ! Then the CFL diagnostic + the explicit/implicit Wvel split.
        ! M2.12c-3: optional partit -> FESOM2 :2381-2537 bounds: Wvel zero owned+halo, the
        ! div(UV*h) edge-scatter over OWNED edges (owned-node Wvel complete, invariant i),
        ! the cumulative sum + the area divide over OWNED nodes; then exchange_nod(Wvel) +
        ! exchange_nod(hnode_new) (:2870-2871) so compute_CFLz/compute_Wvel_split (which loop
        ! owned+halo) read a complete halo. The owned Wvel the ALE dump emits is computed
        ! directly (no halo dependence). hnode_new stays = hnode (linfs); its exchange is
        ! value-neutral but kept for faithfulness.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(inout), target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        real(kind=WP),  intent(in), optional :: water_flux(:)   ! M6a-3: non-linfs freshwater (zstar)
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), n, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2
        real(kind=WP) :: dd, dd1, dddt          ! M6a-2 zstar surface-stretch scalars
        ! c1 is an array over levels (FESOM2 keeps it an array — a deadlock-avoidance note
        ! under OpenMP; here it just carries the per-level edge flux for the slice scatter).
        real(kind=WP) :: c1(mesh%nl-1), c2(mesh%nl-1)
        real(kind=WP), dimension(:,:,:), pointer :: UV, fer_UV
        real(kind=WP), dimension(:,:),   pointer :: Wvel, fer_Wvel

        UV   => dynamics%uv
        Wvel => dynamics%w
        ! M4b GM bolus vertical velocity: fer_w = vertical divergence of fer_uv, mirrors
        ! Wvel exactly under if(Fer_GM) (FESOM2 oce_ale.F90:2366-2524). Byte-neutral when
        ! Fer_GM=.false. (the M2/M3 config) — fer_uv/fer_w are then unallocated + untouched.
        if (Fer_GM) then
            fer_UV   => dynamics%fer_uv
            fer_Wvel => dynamics%fer_w
        end if
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !______________________________________________________________________
        ! zero the vertical velocity over the whole column (below-bottom stays 0)
        do n = 1, nNodL
            Wvel(:, n) = 0.0_WP
            if (Fer_GM) fer_Wvel(:, n) = 0.0_WP
        end do

        !______________________________________________________________________
        ! contributions from levels in the divergence (edge-scatter of div(UV*h))
        do ed = 1, nEdgeO
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmax, nzmin, -1
                c1(nz) = (UV(2,nz,el(1))*deltaX1 - UV(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
                if (Fer_GM) c2(nz) = (fer_UV(2,nz,el(1))*deltaX1 - fer_UV(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
            end do
            Wvel(nzmin:nzmax, enodes(1)) = Wvel(nzmin:nzmax, enodes(1)) + c1(nzmin:nzmax)
            Wvel(nzmin:nzmax, enodes(2)) = Wvel(nzmin:nzmax, enodes(2)) - c1(nzmin:nzmax)
            if (Fer_GM) then
                fer_Wvel(nzmin:nzmax, enodes(1)) = fer_Wvel(nzmin:nzmax, enodes(1)) + c2(nzmin:nzmax)
                fer_Wvel(nzmin:nzmax, enodes(2)) = fer_Wvel(nzmin:nzmax, enodes(2)) - c2(nzmin:nzmax)
            end if

            c1 = 0.0_WP
            if (Fer_GM) c2 = 0.0_WP
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmax, nzmin, -1
                    c1(nz) = -(UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
                    if (Fer_GM) c2(nz) = -(fer_UV(2,nz,el(2))*deltaX2 - fer_UV(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
                end do
                Wvel(nzmin:nzmax, enodes(1)) = Wvel(nzmin:nzmax, enodes(1)) + c1(nzmin:nzmax)
                Wvel(nzmin:nzmax, enodes(2)) = Wvel(nzmin:nzmax, enodes(2)) - c1(nzmin:nzmax)
                if (Fer_GM) then
                    fer_Wvel(nzmin:nzmax, enodes(1)) = fer_Wvel(nzmin:nzmax, enodes(1)) + c2(nzmin:nzmax)
                    fer_Wvel(nzmin:nzmax, enodes(2)) = fer_Wvel(nzmin:nzmax, enodes(2)) - c2(nzmin:nzmax)
                end if
            end if
        end do

        !______________________________________________________________________
        ! cumulative summation of div(UV*h) vertically: W_k = W_k + W_{k+1}
        do n = 1, nNodO
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = nzmax, nzmin, -1
                Wvel(nz, n) = Wvel(nz, n) + Wvel(nz+1, n)
                if (Fer_GM) fer_Wvel(nz, n) = fer_Wvel(nz, n) + fer_Wvel(nz+1, n)
            end do
        end do

        !______________________________________________________________________
        ! divide by depth-dependent cell area -> physical vertical velocity (m/s)
        do n = 1, nNodO
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do nz = nzmin, nzmax
                Wvel(nz, n) = Wvel(nz, n)/mesh%area(n)
                if (Fer_GM) fer_Wvel(nz, n) = fer_Wvel(nz, n)/mesh%area(n)
            end do
        end do

        !______________________________________________________________________
        ! M6a-2 zstar free-surface correction (FESOM2 oce_ale.F90:2755-2821): distribute
        ! the total ssh change (hbar-hbar_old) PROPORTIONALLY over all (full-prism) layers
        !   dd = (hbar-hbar_old)/H,  H = zbar_3d_n(1) - zbar_3d_n(nzmax)
        !   Wvel(nz)      -= (zbar_3d_n(nz)-zbar_3d_n(nzmax)) * dd/dt   (W from the bottom-up
        !                                                                integration of the stretch)
        !   hnode_new(nz)  = hnode(nz) + (zbar_3d_n(nz)-zbar_3d_n(nz+1)) * dd
        ! All divides scalar (no L29 trap). linfs leaves hnode_new=hnode. (Cavity nzmin>1
        ! DEAD; water_flux term is M6a-3 forced — omitted here since unforced water_flux=0.)
        !
        ! FESOM3 bottom at vertices: nzmax is the vertex's OWN column bottom, not
        ! nlevels_nod2D_min. The stretch is meant to span the part of the column where
        ! area(nz)=area(1) -- with the scalar cell now a straight prism that is the WHOLE
        ! column, so more of it is stretched than under FESOM2. (A real physics change,
        ! not a rename: nlevels_nod2D_min is several levels shallower than nlevels_nod2D
        ! over most of a real mesh.)
        if (trim(which_ALE)=='zstar') then
            do n = 1, nNodO
                nzmin = mesh%ulevels_nod2D(n)
                nzmax = mesh%nlevels_nod2D(n) - 1
                if (nzmin == 1) then
                    dd1  = mesh%zbar_3d_n(nzmax, n)
                    dd   = mesh%zbar_3d_n(nzmin, n) - dd1
                    dd   = (mesh%hbar(n) - mesh%hbar_old(n)) / dd
                    dddt = dd/dt
                    do nz = nzmin, nzmax-1
                        Wvel(nz, n)          = Wvel(nz, n) - (mesh%zbar_3d_n(nz,n) - dd1)*dddt
                        mesh%hnode_new(nz,n) = mesh%hnode(nz,n) + (mesh%zbar_3d_n(nz,n) - mesh%zbar_3d_n(nz+1,n))*dd
                    end do
                    ! M6a-3: surface freshwater flux as the upper continuity BC (FESOM2 :2809).
                    ! linfs uses virtual salt; unforced => water_flux=0 (no-op). Present only
                    ! when step_oce threads it (the pressure shim / unforced path omit it).
                    if (present(water_flux)) Wvel(nzmin, n) = Wvel(nzmin, n) - water_flux(n)
                end if
            end do
        end if

        ! linfs: no zlevel/zstar free-surface correction; hnode_new unchanged (= hnode).
        if (is_multirank(partit)) then
            call exchange_nod(Wvel, partit)             ! FESOM2 :2870
            call exchange_nod(mesh%hnode_new, partit)   ! FESOM2 :2871
            if (Fer_GM) call exchange_nod(fer_Wvel, partit)
        end if

        call compute_CFLz(dynamics, mesh, dt, partit)
        call compute_Wvel_split(dynamics, mesh, partit)
    end subroutine vert_vel_ale

    !===========================================================================
    subroutine compute_CFLz(dynamics, mesh, dt, partit)
        ! FESOM2 oce_ale.F90:3126-3213. Vertical CFL per layer (used by the Wvel split):
        !   CFL_z(nz)   += |Wvel(nz)  *dt/hnode_new(nz)|
        !   CFL_z(nz+1)  = |Wvel(nz+1)*dt/hnode_new(nz)|
        ! The two assignments are kept SEPARATE (not folded into one accumulate) — FESOM2
        ! notes this exact form is "for the sake of reproducibility ... (rounding error)"
        ! (L19: the reduction/accumulation form is part of the bits). cflmax + the warning
        ! prints are diagnostics, dropped.
        ! M2.12c-3: optional partit -> the FESOM2 owned+HALO node loops (:3154/:3161 do node=
        ! 1, myDim_nod2D+eDim_nod2D); Wvel/hnode_new are halo-valid from vert_vel_ale's two
        ! exchange_nod's, so CFL_z is built over owned+halo (feeds compute_Wvel_split's
        ! owned+halo Wvel_e/Wvel_i, which the NEXT step's momentum TDMA and vector-invariant
        ! advection read at the halo vertices of owned elements; test_wsplit S5).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        real(kind=WP), intent(in) :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: node, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: c1, c2
        real(kind=WP), dimension(:,:), pointer :: Wvel, CFL_z

        Wvel  => dynamics%w
        CFL_z => dynamics%cfl_z
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do node = 1, nNodL
            CFL_z(:, node) = 0._WP
        end do

        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node) - 1
            do nz = nzmin, nzmax
                c1 = abs(Wvel(nz,   node)*dt/mesh%hnode_new(nz, node))
                c2 = abs(Wvel(nz+1, node)*dt/mesh%hnode_new(nz, node))
                CFL_z(nz,   node) = CFL_z(nz, node) + c1
                CFL_z(nz+1, node) = c2
            end do
        end do
    end subroutine compute_CFLz

    !===========================================================================
    subroutine compute_Wvel_split(dynamics, mesh, partit)
        ! Split Wvel into the explicitly (Wvel_e) and implicitly (Wvel_i) advected parts,
        ! w = w_e + w_i, face by face: w_i = f*w with the implicit share f from
        ! oce_wsplit::wsplit_implicit_fraction(CFL_z, wsplit_mincfl, wsplit_maxcfl), the C^1
        ! Shchepetkin function that replaces FESOM2's hard switch (oce_ale.F90:3217-3265;
        ! WHY, the three branches and the no-hard-switch-mode decision are documented ONCE,
        ! in oce_wsplit). CFL_z is compute_CFLz's (the face flux counted against BOTH
        ! adjacent cells, FESOM2's definition, unchanged), so wsplit_maxcfl keeps its FESOM2
        ! meaning: the explicit CFL_z never exceeds it.
        !
        ! Consumers. Same step, OWNED columns: the FCT adv_tra_vert_impl and the non-FCT
        ! do_wimpl TDMA (both run in the tracer step after vert_vel_ale). NEXT step, at the
        ! halo vertices of owned elements: the momentum TDMA impl_vert_visc_ale (w_i) and
        ! the vector-invariant advection momentum_adv_vinv (w_e); momentum_adv_scalar reads
        ! w_e at owned nodes only. M2.12c-3: optional partit -> the FESOM2 owned+HALO node
        ! loop (:3251 do node=1, myDim_nod2D+eDim_nod2D); Wvel/CFL_z are halo-valid so
        ! Wvel_e/Wvel_i are produced at owned+halo with no trailing exchange (matches
        ! FESOM2; halo == owner bitwise, test_wsplit S5).
        !
        ! w_e = w - w_i (not a second product) makes w_e + w_i == w to 1 ulp, bitwise when
        ! |w_i| >= |w|/2 (Sterbenz); no consumer needs more. use_wsplit=.false. is literally
        ! the previous off path: w_e = w, w_i = +0 over the same face range, validated
        ! nothing. OVERWRITES dynamics%w_e/w_i (the same arrays prescribed as M2.4/M2.5
        ! inputs — fully consumed before this runs). Pinned by test/test_wsplit.F90 part S
        ! (the identities above on every owned+halo face) and part X (transition
        ! smoothness on the tracer path).
        !
        ! PRECONDITIONS, validated at SETUP and assumed here (stateless kernel, like the
        ! rest of src/oce): (wsplit_mincfl, wsplit_maxcfl) admissible -- the driver that
        ! enables use_wsplit called wsplit_assert_params once (the t_dyn default mincfl =
        ! 0.5 is an ABSOLUTE number, so a driver that lowers maxcfl must set mincfl too);
        ! every owned column has >= 2 layers, the row layout of the implicit consumers
        ! adv_tra_vert_impl / do_wimpl (read_mesh's assert_min_layers).
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer       :: node, nz, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: f
        real(kind=WP), dimension(:,:), pointer :: Wvel, Wvel_e, Wvel_i, CFL_z

        Wvel   => dynamics%w
        Wvel_e => dynamics%w_e
        Wvel_i => dynamics%w_i
        CFL_z  => dynamics%cfl_z
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        if (.not. dynamics%use_wsplit) then
            do node = 1, nNodL
                nzmin = mesh%ulevels_nod2D(node)
                nzmax = mesh%nlevels_nod2D(node)
                do nz = nzmin, nzmax
                    Wvel_e(nz, node) = Wvel(nz, node)
                    Wvel_i(nz, node) = 0.0_WP
                end do
            end do
            return
        end if

        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            do nz = nzmin, nzmax
                f = wsplit_implicit_fraction(CFL_z(nz, node), dynamics%wsplit_mincfl, dynamics%wsplit_maxcfl)
                Wvel_i(nz, node) = f*Wvel(nz, node)
                Wvel_e(nz, node) = Wvel(nz, node) - Wvel_i(nz, node)
            end do
        end do
    end subroutine compute_Wvel_split

    !===========================================================================
    subroutine update_thickness_ale(mesh, partit)
        ! FESOM2 oce_ale.F90:1226-1444. Commit the new layer thicknesses at the end of the
        ! step: hnode=hnode_new and recompute helem / zbar_3d_n / Z_3d_n.
        !  - linfs (the M2 path): dh/dt=0, hnode_new=hnode, helem already full-halo-valid
        !    -> NO redistribution; exchange_elem(helem) is value-neutral and SKIPPED (L33).
        !  - M6a-2 zstar (FESOM2 :1378-1436): the stretch from vert_vel_ale is committed.
        !    Node loop (owned+halo: hnode_new was exchange_nod'd in vert_vel_ale) rebuilds
        !    hnode + the depth levels BOTTOM-UP from the fixed anchor zbar_3d_n(nzmax+1)
        !    (nzmax=nlevels_nod2D-2; deepest layer not stretched). Element loop (owned)
        !    averages hnode -> helem, then exchange_elem(helem). ldiag_DVD off (no
        !    rescue_hnode_old bookkeeping); cavity off (nzmin>1 dead). zlevel not ported
        !    (not a target — see docs/plans/2026-06-24-m6-zstar.md).
        type(t_mesh),   intent(inout), target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: n, nz, elem, elnodes(MAX_NV), nzmin, nzmax
        integer :: nNodO, nNodL, nEdgeO, nElemO

        if (trim(which_ALE)=='linfs') return   ! no-op (helem unchanged + full-halo-valid)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        if (trim(which_ALE)=='zstar') then
            !__________________________________________________________________
            ! commit layer thickness + depth levels at node (owned+halo)
            do n = 1, nNodL
                nzmin = mesh%ulevels_nod2D(n)
                nzmax = mesh%nlevels_nod2D(n) - 2      ! vertex's own column (see vert_vel_ale)
                if (nzmin > 1) cycle                       ! cavity (dead)
                do nz = nzmax, nzmin, -1                   ! bottom-up: zbar_3d_n(nz+1) already updated
                    mesh%hnode(nz,n)     = mesh%hnode_new(nz,n)
                    mesh%zbar_3d_n(nz,n) = mesh%zbar_3d_n(nz+1,n) + mesh%hnode_new(nz,n)
                    mesh%Z_3d_n(nz,n)    = mesh%zbar_3d_n(nz+1,n) + mesh%hnode_new(nz,n)/2.0_WP
                end do
            end do
            !__________________________________________________________________
            ! mean layer thickness at element (owned): kite-weighted mean of hnode
            do elem = 1, nElemO
                nzmin = mesh%ulevels(elem)
                nzmax = mesh%nlevels(elem) - 1
                if (nzmin > 1) cycle                       ! cavity (dead)
                elnodes = mesh%elem2D_nodes(:, elem)
                ! FULL element range. FESOM2 stopped one layer short (nzmax-1) and was
                ! safe: nlevels_nod2D_min(n) <= nlevels(e) for every e adjacent to n, so
                ! no node ever stretched layer nlevels(e)-1 and the frozen helem there
                ! still matched. Under the vertex bottom that flips to nlevels(e) <=
                ! nlevels_nod2D(n) -- the element's deepest layer IS stretched at its
                ! deeper nodes, so a frozen helem would break helem == sum(elem_w*hnode), and
                ! helem is the thickness compute_hbar_ale, vert_vel_ale and adv_tra_hor
                ! all use: the SSH/continuity budget and the tracer volume would disagree.
                do nz = nzmin, nzmax
                    mesh%helem(nz,elem) = sum(mesh%elem_w(:,elem)*mesh%hnode(nz,elnodes))
                end do
            end do
        end if

        if (is_multirank(partit)) call exchange_elem(mesh%helem, partit)   ! FESOM2 :1440
    end subroutine update_thickness_ale

end module oce_ale
