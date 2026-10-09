module mod_mesh_analytic
    ! In-memory analytic mesh generator (no file I/O), for fast tests + the
    ! fesom_analytic driver. A regular nx*ny Cartesian grid, 2 triangles per cell,
    ! nl uniform levels. cartesian=.true. so the metric is off. Synthesizes a
    ! 1-rank identity partition too. (Doubly-periodic variant deferred; closed only.)
    ! generate_mixed_patch: the same grid with jittered interior nodes and a mix of
    ! quads and triangle pairs (hybrid meshes). generate_quad_patch: all quads, the
    ! grid sheared into parallelograms (checkerboard test of the SSH matrix). All use
    ! the 4-slot conventions of
    ! mod_mesh (triangle 4th vertex = 1st, boundary edge_tri(2) = edge_tri(1)) and
    ! FESOM's edge order (interior edges first, edge2D_in).
    use mod_precision, only: WP, MP
    use mod_constants, only: r_earth
    use mod_mesh,      only: t_mesh, MAX_NV, MAX_ADJACENT
    use mod_partit,    only: t_partit
    use mod_partitioning, only: set_empty_comms
    use mod_mesh_areas,  only: compute_geometry
    use mod_mesh_rotate, only: init_mesh_rotation
    use mod_mesh_read,   only: enforce_cw_orientation, derive_vertical_bounds
    implicit none
    private
    public :: generate_analytic_mesh, generate_mixed_patch, generate_quad_patch, build_edges
    integer, parameter :: PATCH_TRI = 0, PATCH_MIXED = 1, PATCH_QUAD = 2

contains

    subroutine generate_analytic_mesh(mesh, partit, nx, ny, nl, Lx, Ly, max_depth)
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(inout) :: partit
        integer,        intent(in)    :: nx, ny, nl
        real(kind=WP),  intent(in)    :: Lx, Ly, max_depth
        call build_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth, PATCH_TRI, 0.0_WP)
    end subroutine generate_analytic_mesh

    subroutine generate_mixed_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth)
        ! Mixed triangle/quad patch: cell (i,j) is a quad unless mod(i+2*j,3) == 0, then
        ! it is cut into two triangles (diagonal alternating with i); interior nodes are
        ! displaced by up to 0.15 of the spacing (deterministic), so the quads are general
        ! convex quadrilaterals, not parallelograms.
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(inout) :: partit
        integer,        intent(in)    :: nx, ny, nl
        real(kind=WP),  intent(in)    :: Lx, Ly, max_depth
        call build_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth, PATCH_MIXED, 0.0_WP)
    end subroutine generate_mixed_patch

    subroutine generate_quad_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth, skew)
        ! All-quad patch of parallelograms: node (i,j) at ((i-1)+skew*(j-1))*dx, (j-1)*dy,
        ! cell (i,j) = quad (sw, nw, ne, se). Node (i,j) has colour (-1)**(i+j), so the
        ! vertices of every quad carry (+1,-1,+1,-1).
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(inout) :: partit
        integer,        intent(in)    :: nx, ny, nl
        real(kind=WP),  intent(in)    :: Lx, Ly, max_depth, skew
        call build_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth, PATCH_QUAD, skew)
    end subroutine generate_quad_patch

    subroutine build_patch(mesh, partit, nx, ny, nl, Lx, Ly, max_depth, ptype, skew)
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(inout) :: partit
        integer,        intent(in)    :: nx, ny, nl
        real(kind=WP),  intent(in)    :: Lx, Ly, max_depth
        integer,        intent(in)    :: ptype    ! PATCH_TRI / PATCH_MIXED / PATCH_QUAD
        real(kind=WP),  intent(in)    :: skew     ! PATCH_QUAD: x shift per row, in dx
        integer :: i, j, nod, el, k, sw, se, nw, ne
        logical :: mixed
        real(kind=WP) :: dx, dy

        mixed = ptype == PATCH_MIXED
        mesh%nod2D  = nx * ny
        if (ptype == PATCH_QUAD) then
            mesh%elem2D = (nx-1) * (ny-1)
        else if (mixed) then
            mesh%elem2D = 0
            do j = 1, ny-1
                do i = 1, nx-1
                    mesh%elem2D = mesh%elem2D + merge(2, 1, is_split(i, j))
                end do
            end do
        else
            mesh%elem2D = 2 * (nx-1) * (ny-1)
        end if
        mesh%nl     = nl
        ! Coords stored radians-like (physical metres / r_earth) so the shared
        ! geometry pipeline's `*r_earth^2` (sphere convention) yields physical m^2.
        dx = Lx / real(nx-1, WP) / r_earth; dy = Ly / real(ny-1, WP) / r_earth

        ! ---- nodes (Cartesian) ----
        allocate(mesh%coord_nod2D(2, mesh%nod2D), mesh%geo_coord_nod2D(2, mesh%nod2D))
        allocate(mesh%bc_index_nod2D(mesh%nod2D))
        mesh%bc_index_nod2D = 0
        do j = 1, ny
            do i = 1, nx
                nod = (j-1)*nx + i
                mesh%coord_nod2D(1, nod) = real(i-1, WP) * dx
                mesh%coord_nod2D(2, nod) = real(j-1, WP) * dy
                if (ptype == PATCH_QUAD) mesh%coord_nod2D(1, nod) = &
                    (real(i-1, WP) + skew * real(j-1, WP)) * dx
                if (mixed .and. i > 1 .and. i < nx .and. j > 1 .and. j < ny) then
                    mesh%coord_nod2D(1, nod) = mesh%coord_nod2D(1, nod) &
                        + 0.15_WP * dx * sin(12.9898_WP*real(i, WP) + 78.233_WP*real(j, WP))
                    mesh%coord_nod2D(2, nod) = mesh%coord_nod2D(2, nod) &
                        + 0.15_WP * dy * cos(39.346_WP*real(i, WP) + 11.135_WP*real(j, WP))
                end if
            end do
        end do
        mesh%geo_coord_nod2D = mesh%coord_nod2D

        ! ---- elements (2 triangles per cell; mixed: quads and triangle pairs) ----
        ! A triangle stores its 1st vertex again as the 4th (mod_mesh convention).
        allocate(mesh%elem2D_nodes(MAX_NV, mesh%elem2D), mesh%elem2D_nnodes(mesh%elem2D))
        mesh%elem2D_nodes = 0; mesh%elem2D_nnodes = 3
        el = 0
        do j = 1, ny-1
            do i = 1, nx-1
                sw = (j-1)*nx + i; se = sw + 1; nw = j*nx + i; ne = nw + 1
                if (ptype == PATCH_QUAD) then
                    el = el + 1                                 ! parallelogram, clockwise
                    mesh%elem2D_nodes(:, el) = [sw, nw, ne, se]
                    mesh%elem2D_nnodes(el)   = 4
                else if (.not. mixed) then
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [sw, se, nw]
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [se, ne, nw]
                else if (.not. is_split(i, j)) then
                    el = el + 1                                 ! quad, clockwise
                    mesh%elem2D_nodes(:, el) = [sw, nw, ne, se]
                    mesh%elem2D_nnodes(el)   = 4
                else if (mod(i, 2) == 0) then
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [sw, nw, ne]
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [sw, ne, se]
                else
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [sw, nw, se]
                    el = el + 1
                    mesh%elem2D_nodes(1:3, el) = [se, nw, ne]
                end if
            end do
        end do
        do el = 1, mesh%elem2D
            if (mesh%elem2D_nnodes(el) == 3) mesh%elem2D_nodes(4, el) = mesh%elem2D_nodes(1, el)
        end do

        ! ---- vertical: uniform levels 0 .. -max_depth ----
        allocate(mesh%zbar(nl), mesh%Z(nl-1))
        do k = 1, nl
            mesh%zbar(k) = -max_depth * real(k-1, MP) / real(nl-1, MP)
        end do
        mesh%Z = 0.5_MP * (mesh%zbar(1:nl-1) + mesh%zbar(2:nl))
        allocate(mesh%nlevels(mesh%elem2D), mesh%nlevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels(mesh%elem2D), mesh%ulevels_nod2D(mesh%nod2D))
        allocate(mesh%ulevels_nod2D_max(mesh%nod2D), mesh%nlevels_nod2D_min(mesh%nod2D))
        allocate(mesh%elem_depth(mesh%elem2D), mesh%depth(mesh%nod2D))
        ! FESOM3 bottom at vertices: only the VERTEX columns are set here; the element
        ! bounds and the node ring bounds are derived (below, once nod_in_elem2D exists).
        ! This mesh is flat-bottomed -- every vertex at nl -- so the derived values equal
        ! the old hard-set ones and nothing moves; the point is that the contract lives in
        ! exactly one place.
        mesh%nlevels_nod2D = nl; mesh%ulevels_nod2D = 1
        mesh%nlevels = 0; mesh%ulevels = 1; mesh%ulevels_nod2D_max = 1
        mesh%nlevels_nod2D_min = nl; mesh%elem_depth = -max_depth; mesh%depth = -max_depth

        ! ---- orientation, topology, 1-rank partition, geometry ----
        ! Cartesian coords are in metres; disable cyclic wrapping (huge cyclic
        ! length) so trim_cyclic inside enforce/compute is a no-op here. CW
        ! enforcement must precede geometry (gradient_sca assumes CW node order).
        call init_mesh_rotation(0.0_WP, 0.0_WP, 0.0_WP, 1.0e12_WP)
        call enforce_cw_orientation(mesh, el)
        call build_nod_in_elem_local(mesh)
        call derive_vertical_bounds(mesh, mesh%elem2D, mesh%nod2D, 'generate_analytic_mesh')
        call build_edges(mesh)
        call synth_partit(mesh, partit)
        call compute_geometry(mesh, partit, cartesian=.true.)
    end subroutine build_patch

    pure logical function is_split(i, j)
        ! mixed patch: cell (i,j) is cut into two triangles
        integer, intent(in) :: i, j
        is_split = mod(i + 2*j, 3) == 0
    end function is_split

    subroutine build_nod_in_elem_local(mesh)
        type(t_mesh), intent(inout) :: mesh
        integer :: n, j, node
        allocate(mesh%nod_in_elem2D_num(mesh%nod2D))
        allocate(mesh%nod_in_elem2D(MAX_ADJACENT, mesh%nod2D))
        mesh%nod_in_elem2D_num = 0; mesh%nod_in_elem2D = 0
        do n = 1, mesh%elem2D
            do j = 1, mesh%elem2D_nnodes(n)
                node = mesh%elem2D_nodes(j, n)
                mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
                mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
            end do
        end do
    end subroutine build_nod_in_elem_local

    subroutine build_edges(mesh)
        ! Generic edge + edge_tri builder from elem2D_nodes: edge k of an element joins
        ! its slots k and k+1 over the nv vertices. edge_tri(1) = first element seen,
        ! edge_tri(2) = the second, or the first again on a boundary edge (mod_mesh
        ! convention). FESOM's orientation: edge_tri(1) lies LEFT of edges(1)->edges(2)
        ! (the elements are clockwise, so the first element's side a->b is stored b->a);
        ! the SSH stiffness/rhs dual-face fluxes rely on it. Interior edges first,
        ! boundary edges last (FESOM's order, edge2D_in = number of interior edges).
        type(t_mesh), intent(inout) :: mesh
        integer, allocatable :: head_hi(:,:), head_eid(:,:), head_n(:), perm(:)
        integer, allocatable :: edges(:,:), edge_tri(:,:)
        integer :: n, j, a, b, lo, hi, k, eid, nedge, nv, nin, nbd
        allocate(head_hi(MAX_ADJACENT, mesh%nod2D), head_eid(MAX_ADJACENT, mesh%nod2D))
        allocate(head_n(mesh%nod2D))
        head_n = 0
        allocate(edges(2, MAX_NV*mesh%elem2D), edge_tri(2, MAX_NV*mesh%elem2D))
        edge_tri = 0
        nedge = 0
        do n = 1, mesh%elem2D
            nv = mesh%elem2D_nnodes(n)
            do j = 1, nv
                a = mesh%elem2D_nodes(j, n); b = mesh%elem2D_nodes(mod(j,nv)+1, n)
                lo = min(a,b); hi = max(a,b)
                eid = 0
                do k = 1, head_n(lo)
                    if (head_hi(k, lo) == hi) then; eid = head_eid(k, lo); exit; end if
                end do
                if (eid == 0) then
                    nedge = nedge + 1
                    edges(1, nedge) = b; edges(2, nedge) = a
                    edge_tri(1, nedge) = n
                    head_n(lo) = head_n(lo) + 1
                    head_hi(head_n(lo), lo) = hi; head_eid(head_n(lo), lo) = nedge
                else
                    edge_tri(2, eid) = n
                end if
            end do
        end do
        ! interior edges first, boundary edges last (stable); a boundary edge stores its
        ! one element twice
        allocate(perm(nedge))
        nin = count(edge_tri(2, 1:nedge) > 0)
        k = 0; nbd = nin
        do n = 1, nedge
            if (edge_tri(2, n) > 0) then
                k = k + 1; perm(k) = n
            else
                nbd = nbd + 1; perm(nbd) = n
                edge_tri(2, n) = edge_tri(1, n)
            end if
        end do
        mesh%edge2D = nedge; mesh%edge2D_in = nin
        if (allocated(mesh%edges))    deallocate(mesh%edges)
        if (allocated(mesh%edge_tri)) deallocate(mesh%edge_tri)
        allocate(mesh%edges(2, nedge), mesh%edge_tri(2, nedge))
        mesh%edges    = edges(:, perm)
        mesh%edge_tri = edge_tri(:, perm)
        deallocate(head_hi, head_eid, head_n, perm, edges, edge_tri)
    end subroutine build_edges

    subroutine synth_partit(mesh, partit)
        type(t_mesh),   intent(in)    :: mesh
        type(t_partit), intent(inout) :: partit
        integer :: i
        partit%myDim_nod2D = mesh%nod2D; partit%eDim_nod2D = 0; partit%eXDim_nod2D = 0
        partit%myDim_elem2D = mesh%elem2D; partit%eDim_elem2D = 0; partit%eXDim_elem2D = 0
        partit%myDim_edge2D = mesh%edge2D; partit%eDim_edge2D = 0
        if (allocated(partit%myList_nod2D)) deallocate(partit%myList_nod2D)
        allocate(partit%myList_nod2D(mesh%nod2D)); partit%myList_nod2D = [(i, i=1, mesh%nod2D)]
        call set_empty_comms(partit)      ! as synthesize_1rank: the exchanges are no-ops
    end subroutine synth_partit

end module mod_mesh_analytic
