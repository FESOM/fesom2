module mod_mesh_read
    ! Read a partitioned FESOM mesh into t_mesh (read_mesh: 1 rank, local == global;
    ! read_mesh_local: multi-rank). File formats + topology transcribed from FESOM2
    ! v2.7.3 oce_mesh.F90 (file:line cites inline). Builds nod_in_elem2D, enforces
    ! clockwise element orientation (test_tri), sets up the vertical levels.
    ! Hybrid meshes (docs/HANDOFF.md "Hybrid (mixed triangle + quad) meshes"):
    !  - elem2d.out with 3 or 4 node ids per line, detected on the first element line
    !    (read_elem_line; a triangle repeats its 1st vertex as the 4th);
    !  - zbar and node depths from the NEMO depth.out when aux3d.out is absent
    !    (open_depth_file; node levels still come from nlvls.out only);
    !  - with dist_N/format.out the double-halo nodes (eXDim_nod2D) are read too, with
    !    their coordinates and levels (read_mesh_local);
    !  - an anticlockwise element is turned by swapping slots 2 and nv (a quad 2<->4,
    !    enforce_cw_orientation).
    use mod_precision,   only: WP, MP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh, MAX_NV, MAX_ADJACENT
    use mod_partit,      only: t_partit
    use mod_partitioning, only: read_mesh_dims
    use mod_halo,        only: exchange_nod, allreduce_sum
    use mod_shared_elem, only: build_com_elem2D_shared
    use mod_mesh_rotate, only: init_mesh_rotation, g2r, r2g, trim_cyclic
    implicit none
    private
    public :: read_mesh, derive_vertical_bounds, enforce_cw_orientation

contains

    subroutine read_mesh(mesh, partit, mesh_dir, alpha_deg, beta_deg, gamma_deg, &
                         cyclic_deg, force_rotation, n_cw_swaps)
        type(t_mesh),     intent(inout) :: mesh
        type(t_partit),   intent(inout) :: partit   ! inout: read_mesh_local adds com_elem2D_shared
        character(len=*), intent(in)    :: mesh_dir
        real(kind=WP),    intent(in)    :: alpha_deg, beta_deg, gamma_deg, cyclic_deg
        logical,          intent(in)    :: force_rotation
        integer,          intent(out)   :: n_cw_swaps
        integer :: u, ios, i, n, id, bnd, ncol
        real(kind=WP) :: lon, lat, gx, gy
        logical :: nemo

        if (partit%npes /= 1) then
            ! M2.12a: multi-rank local-mesh remap (needs set_partition to have run).
            call read_mesh_local(mesh, partit, mesh_dir, alpha_deg, beta_deg, gamma_deg, &
                                 cyclic_deg, force_rotation, n_cw_swaps)
            return
        end if
        call init_mesh_rotation(alpha_deg, beta_deg, gamma_deg, cyclic_deg)

        ! ---- nodes: nod2d.out (id lon_deg lat_deg bnd) ----
        open(newunit=u, file=trim(mesh_dir)//'/nod2d.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%nod2D
        allocate(mesh%coord_nod2D(2, mesh%nod2D), mesh%geo_coord_nod2D(2, mesh%nod2D))
        allocate(mesh%bc_index_nod2D(mesh%nod2D))
        do i = 1, mesh%nod2D
            read(u,*) id, lon, lat, bnd
            gx = lon * rad; gy = lat * rad           ! geographic, radians
            if (force_rotation) then
                call g2r(gx, gy, mesh%coord_nod2D(1, id), mesh%coord_nod2D(2, id))
            else
                mesh%coord_nod2D(1, id) = gx; mesh%coord_nod2D(2, id) = gy
            end if
            mesh%bc_index_nod2D(id) = bnd
        end do
        close(u)
        ! geo_coord_nod2D = r2g(coord_nod2D)  (oce_mesh.F90:2480-2492)
        do i = 1, mesh%nod2D
            call r2g(mesh%geo_coord_nod2D(1,i), mesh%geo_coord_nod2D(2,i), &
                     mesh%coord_nod2D(1,i), mesh%coord_nod2D(2,i))
        end do

        ! ---- elements: elem2d.out (3 or 4 node ids per line, see read_elem_line) ----
        open(newunit=u, file=trim(mesh_dir)//'/elem2d.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%elem2D
        allocate(mesh%elem2D_nodes(MAX_NV, mesh%elem2D), mesh%elem2D_nnodes(mesh%elem2D))
        mesh%elem2D_nodes = 0
        ncol = 0
        do i = 1, mesh%elem2D
            call read_elem_line(u, i, mesh%nod2D, ncol, mesh%elem2D_nodes(:, i), mesh%elem2D_nnodes(i))
        end do
        close(u)

        ! ---- edges: edgenum.out, edges.out, edge_tri.out ----
        open(newunit=u, file=trim(mesh_dir)//'/edgenum.out', status='old', action='read', iostat=ios)
        read(u,*) mesh%edge2D
        read(u,*, iostat=ios) mesh%edge2D_in
        if (ios /= 0) mesh%edge2D_in = mesh%edge2D
        close(u)
        allocate(mesh%edges(2, mesh%edge2D), mesh%edge_tri(2, mesh%edge2D))
        open(newunit=u, file=trim(mesh_dir)//'/edges.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%edge2D
            read(u,*) mesh%edges(1, i), mesh%edges(2, i)
        end do
        close(u)
        open(newunit=u, file=trim(mesh_dir)//'/edge_tri.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%edge2D
            read(u,*) mesh%edge_tri(1, i), mesh%edge_tri(2, i)
        end do
        close(u)
        ! Boundary edges store -999 (no 2nd element) on disk. FESOM3 convention (hybrid
        ! meshes): a missing 2nd element is the 1st one, so el(2) == el(1) marks a
        ! boundary edge and no index is ever zero (edge_cross_dxdy(3:4) = 0 there).
        do i = 1, mesh%edge2D
            if (mesh%edge_tri(2, i) <= 0) mesh%edge_tri(2, i) = mesh%edge_tri(1, i)
        end do

        ! ---- vertical level counts: nlvls.out (node) ONLY ----
        ! FESOM3 bottom at vertices: the vertex column is authoritative, so elvls.out is
        ! NOT read. The element bounds are derived in setup_vertical below.
        allocate(mesh%nlevels_nod2D(mesh%nod2D))
        open(newunit=u, file=trim(mesh_dir)//'/nlvls.out', status='old', action='read', iostat=ios)
        do i = 1, mesh%nod2D
            read(u,*) mesh%nlevels_nod2D(i)
        end do
        close(u)

        ! ---- aux3d.out (or depth.out): nl, zbar(nl), per-node depth ----
        call open_depth_file(mesh, mesh_dir, u, nemo)
        allocate(mesh%depth(mesh%nod2D))
        do i = 1, mesh%nod2D
            read(u,*, iostat=ios) mesh%depth(i)
            if (ios /= 0) then; mesh%depth(i) = 0.0_MP; end if
            if (nemo) mesh%depth(i) = -mesh%depth(i)
        end do
        close(u)

        ! ---- topology + orientation + vertical structure ----
        call build_nod_in_elem(mesh)
        call enforce_cw_orientation(mesh, n_cw_swaps)
        call setup_vertical(mesh)
    end subroutine read_mesh

    subroutine read_elem_line(u, iel, nod2D, ncol, nodes, nv)
        ! Element iel from its line of elem2d.out: 3 node ids (FESOM2) or 4 (hybrid
        ! meshes; a triangle repeats its 1st vertex as the 4th). The column count is
        ! detected on the first element line (ncol = 0 on entry) and every later line must
        ! have it; 3 columns -> 4th := 1st. nv = 3 if the 4th vertex is the 1st, else 4.
        ! A line with another number of values, a non-integer or a node id outside
        ! 1..nod2D stops with a message.
        integer, intent(in)    :: u, iel, nod2D
        integer, intent(inout) :: ncol
        integer, intent(out)   :: nodes(MAX_NV), nv
        character(len=256) :: line
        integer :: ios, ntok, k
        logical :: blank_prev
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) call elem_line_error(iel, 'is missing', '')
        ntok = 0; blank_prev = .true.
        do k = 1, len_trim(line)
            if (line(k:k) == ' ' .or. line(k:k) == char(9)) then
                blank_prev = .true.
            else
                if (blank_prev) ntok = ntok + 1
                blank_prev = .false.
            end if
        end do
        if (ncol == 0) then
            if (ntok /= 3 .and. ntok /= 4) &
                call elem_line_error(iel, 'has neither 3 nor 4 node ids:', line)
            ncol = ntok
        end if
        if (ntok /= ncol) call elem_line_error(iel, 'has a wrong number of node ids:', line)
        read(line, *, iostat=ios) nodes(1:ncol)
        if (ios /= 0) call elem_line_error(iel, 'is not a list of integers:', line)
        if (ncol == 3) nodes(4) = nodes(1)
        if (any(nodes < 1) .or. any(nodes > nod2D)) &
            call elem_line_error(iel, 'has a node id outside 1..nod2D:', line)
        nv = 4
        if (nodes(4) == nodes(1)) nv = 3
    end subroutine read_elem_line

    subroutine elem_line_error(iel, what, line)
        integer,          intent(in) :: iel
        character(len=*), intent(in) :: what, line
        write(*,'(a,i0,a)') 'read_mesh: elem2d.out: the line of element ', iel, &
            ' '//what//' "'//trim(line)//'"'
        error stop 1
    end subroutine elem_line_error

    subroutine open_depth_file(mesh, mesh_dir, u, nemo)
        ! zbar (nl levels, made negative if given positive) and Z from aux3d.out (nl,
        ! zbar(nl), node depths) or, if aux3d.out is absent, from the NEMO depth.out (nl,
        ! positive zbar(nl), positive node depths, node layer counts -- the counts are not
        ! read: node levels come from nlvls.out only). Returns the unit u positioned at
        ! the node depths; nemo = depth.out (its depths are positive: callers negate).
        type(t_mesh),     intent(inout) :: mesh
        character(len=*), intent(in)    :: mesh_dir
        integer,          intent(out)   :: u
        logical,          intent(out)   :: nemo
        integer :: ios
        logical :: has_aux3d
        inquire(file=trim(mesh_dir)//'/aux3d.out', exist=has_aux3d)
        nemo = .not. has_aux3d
        if (nemo) then
            open(newunit=u, file=trim(mesh_dir)//'/depth.out', status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,'(a)') 'read_mesh: neither aux3d.out nor depth.out in '//trim(mesh_dir)
                error stop 1
            end if
        else
            open(newunit=u, file=trim(mesh_dir)//'/aux3d.out', status='old', action='read', iostat=ios)
        end if
        read(u,*) mesh%nl
        allocate(mesh%zbar(mesh%nl), mesh%Z(mesh%nl-1))
        read(u,*) mesh%zbar
        if (mesh%zbar(2) > 0.0_MP) mesh%zbar = -mesh%zbar            ! oce_mesh.F90:578
        mesh%Z = 0.5_MP * (mesh%zbar(1:mesh%nl-1) + mesh%zbar(2:mesh%nl))  ! :580-581
    end subroutine open_depth_file

    subroutine build_nod_in_elem(mesh)
        ! node->element adjacency, dense (oce_mesh.F90:2010-2044), 2-pass.
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, node, kmax
        allocate(mesh%nod_in_elem2D_num(mesh%nod2D))
        mesh%nod_in_elem2D_num = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
            end do
        end do
        kmax = maxval(mesh%nod_in_elem2D_num)
        if (kmax > MAX_ADJACENT) then
            write(*,'(a,i0)') 'build_nod_in_elem: MAX_ADJACENT exceeded: ', kmax
            error stop 1
        end if
        allocate(mesh%nod_in_elem2D(MAX_ADJACENT, mesh%nod2D))
        mesh%nod_in_elem2D = 0
        mesh%nod_in_elem2D_num = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
                mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
            end do
        end do
    end subroutine build_nod_in_elem

    subroutine enforce_cw_orientation(mesh, n_swaps, nelem)
        ! test_tri (oce_mesh.F90:1679-1729): r = b1*c2 - b2*c1; if r>0 swap nodes 2,3.
        ! A quad is tested on its sub-triangle (1,2,3) and turned by swapping nodes 2,4
        ! (as fesom_ini does; it has already oriented and checked every quad).
        ! nelem (optional): number of elements to process (multi-rank: owned only,
        ! myDim_elem2D). Defaults to mesh%elem2D (the 1-rank global path, unchanged).
        type(t_mesh),      intent(inout) :: mesh
        integer,           intent(out)   :: n_swaps
        integer, optional, intent(in)    :: nelem
        integer :: n, k, tmp, ne_
        real(kind=WP) :: a1, a2, b1, b2, c1, c2, r
        ne_ = mesh%elem2D
        if (present(nelem)) ne_ = nelem
        n_swaps = 0
        do n = 1, ne_
            a1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(1, n))
            a2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(1, n))
            b1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(2, n)) - a1
            b2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(2, n)) - a2
            c1 = mesh%coord_nod2D(1, mesh%elem2D_nodes(3, n)) - a1
            c2 = mesh%coord_nod2D(2, mesh%elem2D_nodes(3, n)) - a2
            call trim_cyclic(b1)
            call trim_cyclic(c1)
            r = b1*c2 - b2*c1
            if (r > 0.0_WP) then
                ! reverse the cycle keeping vertex 1: swap slots 2 and nv (a triangle
                ! 2<->3, a quad 2<->4; a triangle's 4th slot stays == the 1st)
                k = mesh%elem2D_nnodes(n)
                tmp = mesh%elem2D_nodes(2, n)
                mesh%elem2D_nodes(2, n) = mesh%elem2D_nodes(k, n)
                mesh%elem2D_nodes(k, n) = tmp
                n_swaps = n_swaps + 1
            end if
        end do
    end subroutine enforce_cw_orientation

    subroutine setup_vertical(mesh)
        ! 1-rank vertical structure. ulevels_nod2D = 1 (no cavity); everything else is
        ! derived from the vertex columns by derive_vertical_bounds.
        type(t_mesh), intent(inout) :: mesh
        allocate(mesh%nlevels(mesh%elem2D))
        allocate(mesh%ulevels(mesh%elem2D), mesh%ulevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels_nod2D_max(mesh%nod2D), mesh%nlevels_nod2D_min(mesh%nod2D))
        allocate(mesh%elem_depth(mesh%elem2D))
        mesh%ulevels_nod2D = 1                          ! no cavity
        call derive_vertical_bounds(mesh, mesh%elem2D, mesh%nod2D, 'setup_vertical')
    end subroutine setup_vertical

    subroutine derive_vertical_bounds(mesh, nElem, nNod, where, derive_elems, partit)
        ! FESOM3 BOTTOM AT VERTICES. Shared by all three mesh paths (1-rank read,
        ! analytic, multi-rank). The element vertical bounds are DERIVED from the vertex
        ! columns -- elvls.out is never read:
        !     ulevels(e) = maxval(ulevels_nod2D(elnodes))     (all 1: no cavity)
        !     nlevels(e) = minval(nlevels_nod2D(elnodes))
        ! so an element's layer range ulevels(e)..nlevels(e)-1 is exactly its FULLY WET
        ! prisms. See the contract block in mod_mesh.F90.
        !
        ! nlevels_nod2D_min / ulevels_nod2D_max are then rebuilt from the DERIVED element
        ! bounds. They are NOT aliases of the vertex column: they bound work that reaches
        ! ADJACENT ELEMENTS from a node (oce_muscl_adv.F90:303, oce_fer_gm.F90:82/231).
        !
        ! derive_elems=.false. (multi-rank): nlevels/ulevels were already filled for the
        ! full element halo from the GLOBAL node levels during the mesh scatter, because
        ! elem2D_nodes is populated for OWNED elements only. Only elem_depth, the node
        ! ring bounds and the invariants are done here.
        !
        ! Two vertical invariants are asserted at the end, once, at setup (the kernels
        ! assume them and validate nothing): every owned column has >= 2 layers
        ! (assert_min_layers; partit, when given, makes the verdict collective) and the
        ! vertex-bottom invariant (assert_bottom_invariant).
        type(t_mesh),      intent(inout) :: mesh
        integer,           intent(in)    :: nElem, nNod
        character(len=*),  intent(in)    :: where
        logical, optional, intent(in)    :: derive_elems
        type(t_partit), optional, intent(in) :: partit
        integer :: n, k, nv
        logical :: do_elems
        do_elems = .true.
        if (present(derive_elems)) do_elems = derive_elems
        do n = 1, nElem
            if (do_elems) then
                nv = mesh%elem2D_nnodes(n)
                if (nv <= 0) cycle
                mesh%ulevels(n) = maxval(mesh%ulevels_nod2D(mesh%elem2D_nodes(1:nv, n)))
                mesh%nlevels(n) = minval(mesh%nlevels_nod2D(mesh%elem2D_nodes(1:nv, n)))
            end if
            if (mesh%nlevels(n) > 0) mesh%elem_depth(n) = mesh%zbar(mesh%nlevels(n))
        end do
        do n = 1, nNod
            k = mesh%nod_in_elem2D_num(n)
            if (k <= 0) cycle
            mesh%nlevels_nod2D_min(n) = minval(mesh%nlevels(mesh%nod_in_elem2D(1:k, n)))
            mesh%ulevels_nod2D_max(n) = maxval(mesh%ulevels(mesh%nod_in_elem2D(1:k, n)))
        end do
        call assert_min_layers(mesh, nNod, where, partit)
        call assert_bottom_invariant(mesh, nNod, where)
    end subroutine derive_vertical_bounds

    subroutine assert_min_layers(mesh, nNod, where, partit)
        ! REQUIRED INVARIANT of the vertical kernels: every column (nlevels_nod2D > 0)
        ! has at least 2 layers,
        !     nlevels_nod2D(n) - ulevels_nod2D(n) >= 2.
        ! The explicit vertical advection reads ttf(nzmin-1) in a 1-layer column
        ! (adv_tra_ver_upw1, the trap its module header notes) and the implicit solves
        ! (adv_tra_vert_impl, do_wimpl, the momentum TDMA) lay out a surface row nzmin
        ! and a bottom row nlevels-1 that coincide with one layer. pi/core2 minimum is 4.
        ! Counted over the nNod OWNED columns; with partit at npes > 1 the count is
        ! MPI-summed so every rank stops together (a rank-local stop would leave the
        ! others in the next collective). Pinned by test/test_minlayers.F90.
        type(t_mesh),     intent(in) :: mesh
        integer,          intent(in) :: nNod
        character(len=*), intent(in) :: where
        type(t_partit), optional, intent(in) :: partit
        integer :: n, nbad
        nbad = 0
        do n = 1, nNod
            if (mesh%nod_in_elem2D_num(n) <= 0) cycle
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            if (mesh%nlevels_nod2D(n) - mesh%ulevels_nod2D(n) < 2) then
                if (nbad == 0) write(*,'(a,i0,a,i0,a,i0)') trim(where)// &
                    ': a column with fewer than 2 layers at node ', n, ': ulevels_nod2D=', &
                    mesh%ulevels_nod2D(n), ' nlevels_nod2D=', mesh%nlevels_nod2D(n)
                nbad = nbad + 1
            end if
        end do
        if (present(partit)) then
            if (partit%npes > 1) call allreduce_sum(nbad, partit)
        end if
        if (nbad > 0) then
            write(*,'(a,i0,a)') trim(where)//': ', nbad, &
                ' owned column(s) with fewer than 2 layers; the vertical kernels need 2'
            error stop 1
        end if
    end subroutine assert_min_layers

    subroutine assert_bottom_invariant(mesh, nNod, where)
        ! REQUIRED INVARIANT of the vertex-bottom scheme:
        !     maxval over e in adj(n) of nlevels(e) == nlevels_nod2D(n)
        ! i.e. every vertex's deepest scalar cell has at least one wet adjacent element.
        ! It holds whenever the mesh's nlvls.out is the max over adjacent elvls.out
        ! (verified on pi and core2). This is NOT a quality metric: three UNGUARDED
        ! divides depend on it -- oce_ale.F90 compute_vel_nodes (tx/tvol), vert_vel_ale
        ! (Wvel/area) and oce_pressure_bv.F90 smooth_nod (1/vol) -- so a mesh that violates it produces
        ! NaN rather than a wrong-but-finite answer. Fail loudly at setup instead.
        type(t_mesh),     intent(in) :: mesh
        integer,          intent(in) :: nNod
        character(len=*), intent(in) :: where
        integer :: n, k, deepest
        do n = 1, nNod
            k = mesh%nod_in_elem2D_num(n)
            if (k <= 0) cycle
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            deepest = maxval(mesh%nlevels(mesh%nod_in_elem2D(1:k, n)))
            if (deepest /= mesh%nlevels_nod2D(n)) then
                write(*,'(a,i0,a,i0,a,i0)') trim(where)//': bottom invariant VIOLATED at node ', &
                    n, ': deepest adjacent nlevels=', deepest, ' but nlevels_nod2D=', &
                    mesh%nlevels_nod2D(n)
                write(*,'(a)') '  the vertex bottom cell would have no wet adjacent element'
                error stop 1
            end if
        end do
    end subroutine assert_bottom_invariant

    !==========================================================================
    ! M2.12a: multi-rank local-mesh remap. Builds the per-rank LOCAL mesh from the
    ! GLOBAL mesh files + the partition maps (partit%myList_*, populated by
    ! set_partition BEFORE this call). Transcribed from FESOM2 read_mesh
    ! (oce_mesh.F90:212): a deterministic global->local scatter. We use a full-size
    ! inverse map instead of FESOM2's chunked `mapping` — the result (local arrays) is
    ! identical (a pure scatter). File-read arrays (coords/levels/depth) fill owned+halo
    ! directly; elem2D_nodes is OWNED-only (its nodes are guaranteed local); edges/
    ! edge_tri fill all local edges. nod_in_elem2D is built for OWNED nodes from owned
    ! elements (the find_neighbors halo dance is M2.12b — NOT needed for the owned-node
    ! area gate: compute_node_areas exchanges area owner->halo). The mesh%nod2D/elem2D/
    ! edge2D fields hold the GLOBAL counts (FESOM2 semantics); compute_geometry derives
    ! local bounds from partit.
    subroutine read_mesh_local(mesh, partit, mesh_dir, alpha_deg, beta_deg, gamma_deg, &
                               cyclic_deg, force_rotation, n_cw_swaps)
        type(t_mesh),     intent(inout) :: mesh
        type(t_partit),   intent(inout) :: partit
        character(len=*), intent(in)    :: mesh_dir
        real(kind=WP),    intent(in)    :: alpha_deg, beta_deg, gamma_deg, cyclic_deg
        logical,          intent(in)    :: force_rotation
        integer,          intent(out)   :: n_cw_swaps
        integer :: u, ios, i, g, lid, bnd, gn1, gn2, ge1, ge2, lev, ncol, nv
        integer :: gnodes(MAX_NV)
        integer :: nNodG, nElemG, nEdgeG, edge2D_in_g
        integer :: nNodL, nNodX, nElemO, nElemF, nEdgeL
        real(kind=WP) :: lon, lat, gx, gy, dval
        integer, allocatable :: imap_nod(:), imap_elem(:), imap_edge(:), nlvls_g(:)
        logical :: nemo

        call init_mesh_rotation(alpha_deg, beta_deg, gamma_deg, cyclic_deg)

        ! Node arrays of the mesh FILES (coordinates, boundary index, depth, levels) also
        ! hold the double-halo nodes myDim+eDim+1..nNodX (eXDim_nod2D, dist_N/format.out);
        ! everything derived from the elements keeps nNodL.
        nNodL  = partit%myDim_nod2D  + partit%eDim_nod2D
        nNodX  = nNodL + partit%eXDim_nod2D
        nElemO = partit%myDim_elem2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
        nEdgeL = partit%myDim_edge2D + partit%eDim_edge2D

        ! ---- global counts (file headers) ----
        call read_mesh_dims(mesh_dir, nNodG, nElemG, nEdgeG, edge2D_in_g)
        mesh%nod2D = nNodG; mesh%elem2D = nElemG
        mesh%edge2D = nEdgeG; mesh%edge2D_in = edge2D_in_g

        ! ---- global->local inverse maps (0 = not local) ----
        allocate(imap_nod(nNodG), imap_elem(nElemG), imap_edge(nEdgeG))
        imap_nod = 0; imap_elem = 0; imap_edge = 0
        do i = 1, nNodX;  imap_nod(partit%myList_nod2D(i))   = i; end do
        do i = 1, nElemF; imap_elem(partit%myList_elem2D(i)) = i; end do
        do i = 1, nEdgeL; imap_edge(partit%myList_edge2D(i)) = i; end do

        ! ---- nodes: nod2d.out (id lon lat bnd); scatter owned+halo(+double-halo) coords ----
        allocate(mesh%coord_nod2D(2, nNodX), mesh%geo_coord_nod2D(2, nNodX))
        allocate(mesh%bc_index_nod2D(nNodX))
        open(newunit=u, file=trim(mesh_dir)//'/nod2d.out', status='old', action='read', iostat=ios)
        read(u,*) g                                  ! header == nNodG (line k => global node k)
        do g = 1, nNodG
            read(u,*) i, lon, lat, bnd
            lid = imap_nod(g)
            if (lid > 0) then
                gx = lon * rad; gy = lat * rad
                if (force_rotation) then
                    call g2r(gx, gy, mesh%coord_nod2D(1, lid), mesh%coord_nod2D(2, lid))
                else
                    mesh%coord_nod2D(1, lid) = gx; mesh%coord_nod2D(2, lid) = gy
                end if
                mesh%bc_index_nod2D(lid) = bnd
            end if
        end do
        close(u)
        do i = 1, nNodX
            call r2g(mesh%geo_coord_nod2D(1,i), mesh%geo_coord_nod2D(2,i), &
                     mesh%coord_nod2D(1,i), mesh%coord_nod2D(2,i))
        end do

        ! ---- vertex level counts: nlvls.out (node) ONLY ----
        ! FESOM3 bottom at vertices: the vertex column is authoritative and elvls.out is
        ! NOT read. Kept here in a GLOBAL temp because the element derivation below needs
        ! the levels of an element's vertices by their GLOBAL ids.
        allocate(mesh%nlevels(nElemF), mesh%nlevels_nod2D(nNodX), mesh%ulevels(nElemF))
        mesh%nlevels = 0; mesh%nlevels_nod2D = 0; mesh%ulevels = 0
        allocate(nlvls_g(nNodG))
        nlvls_g = 0
        open(newunit=u, file=trim(mesh_dir)//'/nlvls.out', status='old', action='read', iostat=ios)
        do g = 1, nNodG
            read(u,*) lev
            nlvls_g(g) = lev
            lid = imap_nod(g); if (lid > 0) mesh%nlevels_nod2D(lid) = lev
        end do
        close(u)

        ! ---- elements: elem2d.out (global node ids); store OWNED, localize nodes ----
        ! The SAME pass derives nlevels/ulevels for EVERY local element -- owned AND halo.
        ! It has to: elem2D_nodes is stored for owned elements only, but oce_adv_tra_hor
        ! and vert_vel_ale read nlevels(el(2)) where el(2) can be a halo element, so the
        ! halo entries must be valid. Deriving them here from the global vertex levels
        ! needs no MPI at all and is bit-identical on every rank; the alternative,
        ! exchange_elem_full_2D_i, would add a collective to a path that does not need one.
        ! ulevels is 1 everywhere (no cavity), i.e. maxval over the element's vertices.
        allocate(mesh%elem2D_nodes(MAX_NV, nElemO), mesh%elem2D_nnodes(nElemO))
        mesh%elem2D_nodes = 0
        open(newunit=u, file=trim(mesh_dir)//'/elem2d.out', status='old', action='read', iostat=ios)
        read(u,*) g                                  ! header
        ncol = 0
        do g = 1, nElemG
            call read_elem_line(u, g, nNodG, ncol, gnodes, nv)
            lid = imap_elem(g)
            if (lid > 0) then                         ! ANY local element (owned or halo)
                mesh%nlevels(lid) = minval(nlvls_g(gnodes(1:nv)))
                mesh%ulevels(lid) = 1
            end if
            if (lid > 0 .and. lid <= nElemO) then     ! owned element: also the node list
                mesh%elem2D_nodes(:, lid) = imap_nod(gnodes)
                mesh%elem2D_nnodes(lid)   = nv
            end if
        end do
        close(u)
        deallocate(nlvls_g)

        ! ---- edges + edge_tri: all local edges; localize node + element ids ----
        allocate(mesh%edges(2, nEdgeL), mesh%edge_tri(2, nEdgeL))
        mesh%edges = 0; mesh%edge_tri = 0
        open(newunit=u, file=trim(mesh_dir)//'/edges.out', status='old', action='read', iostat=ios)
        do g = 1, nEdgeG
            read(u,*) gn1, gn2
            lid = imap_edge(g)
            if (lid > 0) then
                mesh%edges(1, lid) = imap_nod(gn1); mesh%edges(2, lid) = imap_nod(gn2)
            end if
        end do
        close(u)
        open(newunit=u, file=trim(mesh_dir)//'/edge_tri.out', status='old', action='read', iostat=ios)
        do g = 1, nEdgeG
            read(u,*) ge1, ge2                        ! ge2 = -999 on a boundary edge
            lid = imap_edge(g)
            if (lid > 0) then
                mesh%edge_tri(1, lid) = imap_elem(ge1)
                if (ge2 > 0) then
                    mesh%edge_tri(2, lid) = imap_elem(ge2)
                else                                  ! boundary: missing 2nd = 1st
                    mesh%edge_tri(2, lid) = mesh%edge_tri(1, lid)
                end if
            end if
        end do
        close(u)

        ! ---- aux3d.out (or depth.out): nl, zbar(nl) [global], per-node depth [scatter] ----
        call open_depth_file(mesh, mesh_dir, u, nemo)
        allocate(mesh%depth(nNodX))
        mesh%depth = 0.0_MP
        do g = 1, nNodG
            read(u,*, iostat=ios) dval
            if (ios /= 0) exit
            if (nemo) dval = -dval
            lid = imap_nod(g); if (lid > 0) mesh%depth(lid) = dval
        end do
        close(u)

        deallocate(imap_nod, imap_edge)

        ! ---- topology + orientation + vertical structure (LOCAL) ----
        call build_nod_in_elem_local(mesh, partit)
        ! M2.12b: complete nod_in_elem2D for HALO nodes (the find_neighbors halo dance)
        ! — needed by MUSCL find_up_downwind_elements, which reads a
        ! halo node's full element list (reaching eXDim). imap_elem is the full-halo
        ! global->local inverse map (1..nElemF) used to re-localize.
        call complete_nod_in_elem_halo(mesh, partit, imap_elem)
        call enforce_cw_orientation(mesh, n_cw_swaps, nelem=nElemO)
        call setup_vertical_local(mesh, partit)
        deallocate(imap_elem)
        ! one designated owner per shared (border) element: the list update_vel syncs uv over
        call build_com_elem2D_shared(mesh, partit)
    end subroutine read_mesh_local

    subroutine build_nod_in_elem_local(mesh, partit)
        ! Local node->element adjacency for OWNED nodes (FESOM2 find_neighbors
        ! oce_mesh.F90:2021-2049): loop owned elements, count/store only for OWNED
        ! nodes (node <= myDim_nod2D). Halo nodes' lists stay empty (filled via the
        ! M2.12b find_neighbors dance; not needed for the owned-node area gate).
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        integer :: n, j, node, kmax, nNodL, nElemO
        nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
        nElemO = partit%myDim_elem2D
        allocate(mesh%nod_in_elem2D_num(nNodL))
        mesh%nod_in_elem2D_num = 0
        do n = 1, nElemO
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                if (node > partit%myDim_nod2D) cycle
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
            end do
        end do
        kmax = maxval(mesh%nod_in_elem2D_num)
        if (kmax > MAX_ADJACENT) then
            write(*,'(a,i0)') 'build_nod_in_elem_local: MAX_ADJACENT exceeded: ', kmax
            error stop 1
        end if
        allocate(mesh%nod_in_elem2D(MAX_ADJACENT, nNodL))
        mesh%nod_in_elem2D = 0
        mesh%nod_in_elem2D_num = 0
        do n = 1, nElemO
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                if (node > partit%myDim_nod2D) cycle
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
                mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
            end do
        end do
    end subroutine build_nod_in_elem_local

    subroutine complete_nod_in_elem_halo(mesh, partit, imap_elem)
        ! M2.12b find_neighbors halo dance (FESOM2 oce_mesh.F90:2051-2074). On entry
        ! nod_in_elem2D holds LOCAL element indices for OWNED nodes only (halo nodes
        ! empty, build_nod_in_elem_local). Completes it for HALO nodes:
        !   (1) exchange_nod the per-node count (owner -> halo);
        !   (2) per slot: convert OWNED entries local->GLOBAL, exchange owner->halo,
        !       store back (now GLOBAL on owned+halo);
        !   (3) re-localize every entry GLOBAL->LOCAL via the full-halo inverse map
        !       imap_elem (covers myDim+eDim+eXDim). The eXDim element halo guarantees
        !       a halo node's full element neighbourhood is local, so imap_elem(g)>0.
        ! After this every owned+halo node's element list is LOCAL and complete.
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        integer,        intent(in)    :: imap_elem(:)
        integer :: slot, n, kmax, g, nNodO, nNodL
        integer, allocatable :: temp_i(:)
        nNodO = partit%myDim_nod2D
        nNodL = partit%myDim_nod2D + partit%eDim_nod2D
        ! (1) per-node element count owner -> halo
        call exchange_nod(mesh%nod_in_elem2D_num, partit)
        ! (2) per slot: local id -> global on owned, exchange, store back (global on all)
        kmax = size(mesh%nod_in_elem2D, 1)
        allocate(temp_i(nNodL))
        do slot = 1, kmax
            temp_i = 0
            do n = 1, nNodO
                if (mesh%nod_in_elem2D(slot, n) > 0) &
                    temp_i(n) = partit%myList_elem2D(mesh%nod_in_elem2D(slot, n))
            end do
            call exchange_nod(temp_i, partit)
            mesh%nod_in_elem2D(slot, :) = temp_i
        end do
        deallocate(temp_i)
        ! (3) re-localize global -> local through the full-halo inverse map (eXDim)
        do n = 1, nNodL
            do slot = 1, mesh%nod_in_elem2D_num(n)
                g = mesh%nod_in_elem2D(slot, n)
                if (g > 0) mesh%nod_in_elem2D(slot, n) = imap_elem(g)
            end do
        end do
    end subroutine complete_nod_in_elem_halo

    subroutine setup_vertical_local(mesh, partit)
        ! Local vertical structure. nlevels/ulevels were already derived for the FULL
        ! element halo during the mesh scatter (from the global vertex levels -- see
        ! read_mesh_local), because elem2D_nodes is owned-only. Here we only add
        ! elem_depth, the node ring bounds and the invariant, hence derive_elems=.false.
        !
        ! nlevels_nod2D_min / ulevels_nod2D_max are computed on OWNED nodes and then
        ! halo-exchanged: an owned node's element neighbourhood is complete, a halo
        ! node's is not.
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        integer :: nNodL, nElemF, nNodO
        nNodO  = partit%myDim_nod2D
        nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
        allocate(mesh%ulevels_nod2D(nNodL))
        allocate(mesh%ulevels_nod2D_max(nNodL), mesh%nlevels_nod2D_min(nNodL))
        allocate(mesh%elem_depth(nElemF))
        mesh%ulevels_nod2D = 1; mesh%ulevels_nod2D_max = 1   ! no cavity
        mesh%nlevels_nod2D_min = 0
        mesh%elem_depth = 0.0_MP
        call derive_vertical_bounds(mesh, nElemF, nNodO, 'setup_vertical_local', &
                                    derive_elems=.false., partit=partit)
        call exchange_nod(mesh%nlevels_nod2D_min, partit)
        call exchange_nod(mesh%ulevels_nod2D_max, partit)
    end subroutine setup_vertical_local

end module mod_mesh_read
