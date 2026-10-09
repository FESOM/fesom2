module oce_ssh_rhs
    ! SSH (free-surface) implicit-solve assembly: the elevation stiffness matrix and
    ! the elevation right-hand side. Run in the timestep AFTER impl_vert_visc_ale and
    ! BEFORE the CG solve (FESOM2 oce_ale.F90:3920-3930, the use_ssh_se_subcycl=.false.
    ! branch). Transcribed from FESOM2 v2.7.3 oce_ale.F90:
    !   init_stiff_mat_ale   (:1584-1875) — builds the CSR stiffness operator once
    !   compute_ssh_rhs_ale  (:2012-2147) — assembles ssh_rhs each step
    !
    ! The split-implicit free-surface equation (eq. 18 of "FESOM2: from finite elements
    ! to finite volumes") is
    !       [ M/dt - alpha*theta*g*dt*div(H grad) ] d_eta = ssh_rhs
    ! a symmetric sparse linear system in the SSH increment d_eta, solved by the
    ! preconditioned CG in oce_ssh_solve. This module builds the left-hand-side matrix
    ! (init_stiff_mat_ale) and the right-hand side (compute_ssh_rhs_ale).
    !
    ! STIFFNESS MATRIX (init_stiff_mat_ale):
    !  - CSR sparsity from the node neighbourhood: self, the edge-connected nodes (in
    !    edge order, as before), then the remaining vertices of the elements around the
    !    node (the quad diagonals; none on a triangle mesh, whose CSR is unchanged). Row
    !    lengths are dynamic. The diagonal (self) is the FIRST entry of every row
    !    (n_pos(1,n)=n) -- the Jacobi preconditioner in oce_ssh_solve relies on it.
    !  - Stiffness part (loop over edges, both adjacent elements): the H*div operator
    !    factor*(zbar_e_bot-zbar_e_srf)*(gradient_sca_x*edge_cross_dy - gradient_sca_y*
    !    edge_cross_dx) over the 4 vertex slots (x 1:4, y 5:8; a triangle's slot 4 is 0
    !    and repeats vertex 1, so it adds an exact zero), scattered +/- into the two
    !    edge-node rows, factor=g*dt*alpha*theta. With the Green-Gauss coefficients the
    !    matrix stays symmetric on quads (dual-face normals ~ gradient coefficients).
    !  - Checkerboard: the Green-Gauss quad gradient vanishes for (+1,-1,+1,-1) on the
    !    vertices of ANY quad (coefficient j ~ x_{j+1}-x_{j-1}, so c_1=-c_3, c_2=-c_4),
    !    so on an all-quad mesh with a consistent 2-colouring the stiffness part has that
    !    null mode (asserted in test/test_ssh_stiff.F90); only the mass part sees it.
    !  - Mass part: + areasvol(surf)/dt added to each row diagonal.
    !  - For which_ale=='linfs' the matrix is built ONCE here with the unperturbed depth
    !    and never updated (FESOM2 oce_ale.F90:3921 calls update_stiff_mat_ale only for
    !    NON-linfs). So this single assembly IS the production matrix on pi/reduced-M2.
    !
    ! SSH RHS (compute_ssh_rhs_ale): the depth-integrated horizontal-divergence of the
    ! provisional transport alpha*(UV+UV_rhs) over the scalar control volume, as an
    ! edge flux scattered +/- into the two edge nodes; plus the linfs (1-alpha)*ssh_rhs_old
    ! term (which vanishes here since alpha=1, theta=1 on pi/reduced-M2). UV_rhs is the
    ! post-impl_vert_visc_ale velocity rhs (the M2.5-gated uv_rhs_ivv).
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern again — every operand already
    ! byte-pinned, so faithful transcription byte-matches like M1.1-M2.5):
    !  - The stiffness value scatter (+= fy*factor over edges) and the ssh_rhs scatter
    !    (+= (c1+c2) over edges) are FP-order-sensitive, but the 1-rank global edge order
    !    + edge_tri + edges + gradient_sca + edge_cross_dxdy + zbar_e_bot + areasvol +
    !    helem are all geometry/M2.5-proven (max|delta|=0), so the accumulation order
    !    matches FESOM2 exactly. The FESOM2 oracle is ENABLE_OPENMP=OFF (L16), so its
    !    !$OMP ORDERED / omp_set_lock edge-scatter compiles to a plain serial loop in
    !    edge order — identical to this serial transcription.
    !  - factor = g*dt*alpha*theta and areasvol/dt are runtime products/quotients with
    !    byte-identical operands on both sides; the Intel -no-prec-div reciprocal matches
    !    (L7/L10/L14). dt MUST be the pi namelist dt=86400/36 (the value FESOM2's
    !    init_stiff_mat_ale saw at ocean_setup), computed the SAME way (not a 2400.0
    !    literal — see the driver) so the bits agree.
    !  - zbar_e_srf = zbar(ulevels(elem)) = zbar(1) = 0 on pi (no cavity; FESOM2 default
    !    oce_ale.F90:525). Computed inline (no new mesh field): faithful to the non-cavity
    !    formula, and 0 on pi.
    !  - elem2D_nodes(1:4,el) and gradient_sca(1:4|5:8): unified 4-slot form (hybrid
    !    meshes); the scatter is a scalar loop because a triangle's npos(4)=npos(1).
    !
    ! 1-RANK SCOPE (multi-rank lifted at M2.12):
    !  - The CSR is built in LOCAL numbering (rowptr_loc/colind_loc, used by the CG); at
    !    1-rank local==global natural, so the FESOM2 global-contiguous remap (rpart.out
    !    mapping, the per-PE nza offset) is the identity and is dropped. rowptr/colind are
    !    set equal to rowptr_loc/colind_loc (their 1-rank values).
    !  - myDim_nod2D -> nod2D, myDim_edge2D -> edge2D, eDim_nod2D -> 0. The
    !    row<=myDim_nod2D / n1<=myDim_nod2D guards are always true at 1-rank (kept for
    !    multi-rank faithfulness). exchange_nod(ssh_rhs) at the tail is a no-op, dropped.
    use mod_precision,  only: WP, MP
    use mod_constants,  only: g
    use mod_param_phys, only: alpha, theta
    use mod_config,     only: which_ALE          ! M6a-3: 'linfs'/'zlevel'/'zstar'
    use mod_mesh,       only: t_mesh, MAX_NV
    use mod_dyn,        only: t_dyn
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,       only: exchange_nod
    implicit none
    private
    public :: init_stiff_mat_ale, compute_ssh_rhs_ale, update_stiff_mat_ale

contains

    !===========================================================================
    subroutine init_stiff_mat_ale(mesh, dt, partit)
        ! Build mesh%ssh_stiff (CSR sparsity + the linfs stiffness/mass values). Call
        ! ONCE; dt is the pi namelist timestep (see module/driver notes). Idempotent
        ! it is NOT (allocates ssh_stiff) — call exactly once per mesh.
        !
        ! M2.12c-2 (multi-rank, OPTIONAL partit): the owned rows assemble FULLY LOCALLY
        ! — NO halo exchange, NO mesh-infra extension. Verified partition invariants
        ! (pi dist_2/8): (i) every edge incident to an owned node is owned, so looping
        ! the owned edges (1..nEdgeO) visits every contribution to an owned row; (ii)
        ! both elements of an owned edge are OWNED (edge_tri(:,ed) <= nElemO), so the
        ! elem2D_nodes/gradient_sca/zbar_e_bot reads stay within the owned-only arrays
        ! (M2.12a) — the FESOM2 oracle relies on the SAME invariant (its elem2D_nodes/
        ! gradient_sca are owned-only too, oce_mesh.F90:497,2466); (iii) an owned
        ! element's nodes are within owned+halo (<= nNodL), so n_num(elnodes) is in
        ! bounds and the CSR can have HALO columns (colind_loc up to nNodL). partit
        ! absent OR npes==1 => nNodO=nNodL=mesh%nod2D etc. (the proven 1-rank path
        ! VERBATIM). FESOM2's global-contiguous rowptr/colind remap (for an external
        ! solver) is dropped — the CG uses rowptr_loc/colind_loc only.
        type(t_mesh),  intent(inout), target :: mesh
        real(kind=WP), intent(in)            :: dt
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        integer              :: n, n1, n2, i, k, j, row, ed, el1
        integer              :: elnodes(MAX_NV), el(2), npos(MAX_NV)
        integer              :: offset, nini, nend, maxrow
        integer              :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP)        :: zsrf
        integer, allocatable :: n_num(:), n_pos(:,:), mark(:)
        real(kind=WP), allocatable :: depth_e(:)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        associate(ssh_stiff => mesh%ssh_stiff)
        !__________________________________________________________________
        ! a) neighbourhood: n_num(n) = #neighbours of OWNED node n (incl. self),
        !    n_pos(:,n) = their local indices (a neighbour may be a HALO node),
        !    n_pos(1,n)=n (self -> diagonal first), then the edge neighbours in edge
        !    order, then the other vertices of the elements around n (quad diagonals).
        !    Row length <= 1 + #edges + #elements (a quad adds at most its diagonal
        !    vertex), counted first. n_num is sized owned+halo because it is reused
        !    below as the reverse-map indexed by a (possibly halo) colind.
        allocate(n_num(nNodL), mark(nNodL))
        n_num = 0
        do n = 1, nEdgeO
            n1 = mesh%edges(1,n)
            n2 = mesh%edges(2,n)
            if (n1 <= nNodO) n_num(n1) = n_num(n1)+1
            if (n2 <= nNodO) n_num(n2) = n_num(n2)+1
        end do
        maxrow = 0
        do n = 1, nNodO
            maxrow = max(maxrow, 1 + n_num(n) + mesh%nod_in_elem2D_num(n))
        end do
        allocate(n_pos(maxrow, nNodO))
        n_pos = 0
        do n = 1, nNodO
            n_num(n)   = 1
            n_pos(1,n) = n
        end do
        do n = 1, nEdgeO
            n1 = mesh%edges(1,n)
            n2 = mesh%edges(2,n)
            if (n1 <= nNodO) then
                n_pos(n_num(n1)+1, n1) = n2
                n_num(n1) = n_num(n1)+1
            end if
            if (n2 <= nNodO) then
                n_pos(n_num(n2)+1, n2) = n1
                n_num(n2) = n_num(n2)+1
            end if
        end do
        mark = 0
        do n = 1, nNodO
            mark(n_pos(1:n_num(n), n)) = n
            do k = 1, mesh%nod_in_elem2D_num(n)
                el1 = mesh%nod_in_elem2D(k, n)
                do j = 1, mesh%elem2D_nnodes(el1)
                    n1 = mesh%elem2D_nodes(j, el1)
                    if (mark(n1) == n) cycle
                    mark(n1) = n
                    n_num(n) = n_num(n)+1
                    n_pos(n_num(n), n) = n1
                end do
            end do
        end do

        !__________________________________________________________________
        ! b) CSR row pointers (OWNED rows) + nonzero count
        ssh_stiff%dim = mesh%nod2D   ! GLOBAL node count (FESOM2 ssh_stiff%dim=nod2D)
        allocate(ssh_stiff%rowptr(nNodO+1), ssh_stiff%rowptr_loc(nNodO+1))
        ssh_stiff%rowptr_loc(1) = 1
        do n = 1, nNodO
            ssh_stiff%rowptr_loc(n+1) = ssh_stiff%rowptr_loc(n) + n_num(n)
        end do
        ssh_stiff%nza = ssh_stiff%rowptr_loc(nNodO+1) - 1

        !__________________________________________________________________
        ! c) CSR column indices (local; may be HALO) + zero the values
        allocate(ssh_stiff%colind(ssh_stiff%nza), ssh_stiff%colind_loc(ssh_stiff%nza))
        allocate(ssh_stiff%values(ssh_stiff%nza))
        ssh_stiff%values = 0.0_WP
        do n = 1, nNodO
            nini = ssh_stiff%rowptr_loc(n)
            nend = ssh_stiff%rowptr_loc(n+1) - 1
            ssh_stiff%colind_loc(nini:nend) = n_pos(1:n_num(n), n)
        end do
        ! the CG uses only rowptr_loc/colind_loc; keep rowptr/colind as their copies
        ! (FESOM2's global-contiguous remap is for an external solver, dropped here).
        ssh_stiff%rowptr = ssh_stiff%rowptr_loc
        ssh_stiff%colind = ssh_stiff%colind_loc

        !__________________________________________________________________
        ! d) stiffness part: factor * H * div over OWNED edges, H = zbar_e_bot - zbar_e_srf
        !    (zbar_e_srf = zbar(ulevels) = 0, no cavity). n_num is reused as the reverse map.
        allocate(depth_e(nElemO))
        do el1 = 1, nElemO
            zsrf = mesh%zbar(mesh%ulevels(el1))
            depth_e(el1) = mesh%zbar_e_bot(el1) - zsrf
        end do
        call add_edge_stiffness(mesh, nNodO, nEdgeO, depth_e, g*dt*alpha*theta, n_num)
        deallocate(depth_e)

        !__________________________________________________________________
        ! e) mass part: + areasvol(surf)/dt on the row diagonal (the first CSR entry).
        !    Skip cavity nodes (ulevels_nod2D>1: rigid-lid, no eta time-derivative).
        do row = 1, nNodO
            if (mesh%ulevels_nod2D(row) > 1) cycle
            offset = ssh_stiff%rowptr_loc(row)
            ssh_stiff%values(offset) = ssh_stiff%values(offset) &
                                     + mesh%areasvol(row)/dt
        end do

        deallocate(n_pos, n_num, mark)
        end associate
    end subroutine init_stiff_mat_ale

    !===========================================================================
    subroutine update_stiff_mat_ale(mesh, dt, partit)
        ! Per-step SSH-stiffness 2nd-term update for non-linfs ALE (zlevel/zstar).
        ! Transcribed from FESOM2 v2.7.3 oce_ale.F90:1892-2001 (update_stiff_mat_ale),
        ! called every step (oce_ale.F90:3921, `if(.not. linfs)`) BEFORE compute_ssh_rhs_ale.
        !
        ! init_stiff_mat_ale built the stiffness 2nd term over the UNPERTURBED depth
        ! (zbar_e_bot-zsrf). This ADDS -dhe(elem) per step, where dhe = hbar(n+1/2)-
        ! hbar(n-1/2) is the element-interpolated SSH change from the previous step's
        ! compute_hbar_ale (LAGGED; step-1 dhe=0 -> stiffness unchanged step 1). The
        ! additions ACCUMULATE across steps so the matrix tracks the moving free surface.
        !
        ! The same edge assembly as init_stiff_mat_ale's stiffness part (d), shared in
        ! add_edge_stiffness — same owned-edge order, same CSR (rowptr_loc/colind_loc),
        ! same gradient_sca/edge_cross_dxdy — only (zbar_e_bot-zsrf) -> -dhe(el). The
        ! edges(1)=+fy / edges(2)=-fy split == the oracle's j-loop with the j==2 flip.
        ! M6a-2: optional partit -> owned-edge loop; mutates mesh%ssh_stiff%values.
        type(t_mesh),   intent(inout), target :: mesh
        real(kind=WP),  intent(in)            :: dt
        type(t_partit), intent(in), optional  :: partit
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        integer, allocatable :: rev(:)   ! reverse-map local-node -> CSR position (init's n_num reuse)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        allocate(rev(nNodL))
        call add_edge_stiffness(mesh, nNodO, nEdgeO, -mesh%dhe(1:nElemO), g*dt*alpha*theta, rev)
        deallocate(rev)
    end subroutine update_stiff_mat_ale

    !===========================================================================
    subroutine add_edge_stiffness(mesh, nNodO, nEdgeO, depth_e, factor, rev)
        ! The edge assembly of init_stiff_mat_ale (d) and update_stiff_mat_ale: for every
        ! OWNED edge and each of its elements (both OWNED, invariant ii; a boundary edge
        ! has one), fy(k) = depth_e(el)*(gradient_sca_x(k)*edge_cross_dy -
        ! gradient_sca_y(k)*edge_cross_dx) over the 4 slots, sign-flipped for the 2nd
        ! element, added as +fy*factor to the row of edges(1) and -fy*factor to the row of
        ! edges(2) when that row is owned. rev (size >= owned+halo) is the work reverse map
        ! local node -> CSR position in the current row.
        type(t_mesh),  intent(inout), target :: mesh
        integer,       intent(in)            :: nNodO, nEdgeO
        real(kind=WP), intent(in)            :: depth_e(:), factor
        integer,       intent(inout)         :: rev(:)
        integer       :: ed, el(2), elnodes(MAX_NV), i, j, k, n, row, npos(MAX_NV)
        real(kind=WP) :: fy(MAX_NV)
        associate(ssh_stiff => mesh%ssh_stiff)
        rev = 0
        do ed = 1, nEdgeO
            el = mesh%edge_tri(:, ed)
            do i = 1, 2   ! the two elements sharing edge ed (both OWNED, invariant ii)
                if (i == 2 .and. el(2) == el(1)) cycle   ! boundary edge has only one element
                elnodes = mesh%elem2D_nodes(1:MAX_NV, el(i))
                fy = depth_e(el(i)) * &
                     ( mesh%gradient_sca(1:4,el(i)) * mesh%edge_cross_dxdy(2*i  ,ed)  &
                      -mesh%gradient_sca(5:8,el(i)) * mesh%edge_cross_dxdy(2*i-1,ed) )
                if (i==2) fy = -fy
                do j = 1, 2   ! the rows of the edge's two nodes: +fy, -fy
                    row = mesh%edges(j, ed)
                    if (row > nNodO) cycle
                    do n = ssh_stiff%rowptr_loc(row), ssh_stiff%rowptr_loc(row+1)-1
                        rev(ssh_stiff%colind_loc(n)) = n
                    end do
                    npos = rev(elnodes)
                    if (j == 1) then
                        do k = 1, MAX_NV
                            ssh_stiff%values(npos(k)) = ssh_stiff%values(npos(k)) + fy(k)*factor
                        end do
                    else
                        do k = 1, MAX_NV
                            ssh_stiff%values(npos(k)) = ssh_stiff%values(npos(k)) - fy(k)*factor
                        end do
                    end if
                end do
            end do
        end do
        end associate
    end subroutine add_edge_stiffness

    !===========================================================================
    subroutine compute_ssh_rhs_ale(dynamics, mesh, partit, water_flux)
        ! Assemble dynamics%ssh_rhs = depth-integrated horizontal divergence of
        ! alpha*(UV+UV_rhs), scattered as an edge flux into the two edge nodes
        ! (+ the linfs (1-alpha)*ssh_rhs_old term, = 0 on pi since alpha=1).
        ! M6a-3: for non-linfs (zstar) the real freshwater flux enters here too
        ! (ssh_rhs -= alpha*water_flux*areasvol; FESOM2 oce_ale.F90:2122-2134).
        ! M2.12c: optional partit -> ssh_rhs zeroed at owned+halo (FESOM2 :2045
        ! do n=1,myDim_nod2D+eDim_nod2D), the flux scatter over OWNED edges (FESOM2 :2053
        ! do ed=1,myDim_edge2D — every edge incident to an owned node is owned, so the
        ! owned-node ssh_rhs is fully accumulated), the linfs term over owned nodes, then
        ! exchange_nod(ssh_rhs) (FESOM2 :2145) to fill the halo for the CG.
        type(t_dyn),  intent(inout), target :: dynamics
        type(t_mesh), intent(in),    target :: mesh
        type(t_partit), intent(in), optional :: partit
        real(kind=WP),  intent(in), optional :: water_flux(:)   ! M6a-3: non-linfs freshwater (zstar)
        !______________________________________________________________________
        integer       :: ed, el(2), enodes(2), nz, n, nzmin, nzmax
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: c1, c2, deltaX1, deltaX2, deltaY1, deltaY2
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:),     pointer :: ssh_rhs, ssh_rhs_old

        UV          => dynamics%uv
        UV_rhs      => dynamics%uv_rhs
        ssh_rhs     => dynamics%ssh_rhs
        ssh_rhs_old => dynamics%ssh_rhs_old
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodL
            ssh_rhs(n) = 0.0_WP
        end do

        do ed = 1, nEdgeO
            enodes = mesh%edges(:, ed)
            el     = mesh%edge_tri(:, ed)

            !__________________________________________________________________
            ! alpha*nabla*int(U_n + U_rhs) dz contribution from el(1)
            c1      = 0.0_WP
            deltaX1 = mesh%edge_cross_dxdy(1, ed)
            deltaY1 = mesh%edge_cross_dxdy(2, ed)
            nzmin   = mesh%ulevels(el(1))
            nzmax   = mesh%nlevels(el(1)) - 1
            do nz = nzmin, nzmax
                c1 = c1 + alpha*((UV(2,nz,el(1))+UV_rhs(2,nz,el(1)))*deltaX1 -  &
                                 (UV(1,nz,el(1))+UV_rhs(1,nz,el(1)))*deltaY1)*mesh%helem(nz,el(1))
            end do

            !__________________________________________________________________
            ! ... and from el(2), unless ed is a boundary edge
            c2 = 0.0_WP
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, ed)
                deltaY2 = mesh%edge_cross_dxdy(4, ed)
                nzmin   = mesh%ulevels(el(2))
                nzmax   = mesh%nlevels(el(2)) - 1
                do nz = nzmin, nzmax
                    c2 = c2 - alpha*((UV(2,nz,el(2))+UV_rhs(2,nz,el(2)))*deltaX2 -  &
                                     (UV(1,nz,el(2))+UV_rhs(1,nz,el(2)))*deltaY2)*mesh%helem(nz,el(2))
                end do
            end if

            !__________________________________________________________________
            ! net "flux" scattered into the two edge nodes
            ssh_rhs(enodes(1)) = ssh_rhs(enodes(1)) + (c1+c2)
            ssh_rhs(enodes(2)) = ssh_rhs(enodes(2)) - (c1+c2)
        end do

        ! water-flux term (FESOM2 oce_ale.F90:2122-2138). linfs: ssh_rhs += (1-alpha)*ssh_rhs_old
        ! (= 0 since alpha=1). M6a-3 non-linfs (zstar): the REAL freshwater flux enters the SSH rhs
        ! — ssh_rhs -= alpha*water_flux*areasvol (+ the (1-alpha)*ssh_rhs_old term). Open ocean
        ! only (cavity dead, use_cavity=.false. -> use_cavity_fw2press branch skipped).
        if (trim(which_ALE)/='linfs' .and. present(water_flux)) then
            do n = 1, nNodO
                if (mesh%ulevels_nod2D(n) > 1) cycle   ! cavity (dead)
                ssh_rhs(n) = ssh_rhs(n) - alpha*water_flux(n)*mesh%areasvol(n) &
                                        + (1.0_WP-alpha)*ssh_rhs_old(n)
            end do
        else
            do n = 1, nNodO
                ssh_rhs(n) = ssh_rhs(n) + (1.0_WP-alpha)*ssh_rhs_old(n)
            end do
        end if
        if (is_multirank(partit)) call exchange_nod(ssh_rhs, partit)   ! FESOM2 :2145
    end subroutine compute_ssh_rhs_ale

end module oce_ssh_rhs
