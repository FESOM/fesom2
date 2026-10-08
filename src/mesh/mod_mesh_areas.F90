module mod_mesh_areas
    ! Core 2D mesh geometry for mixed triangle/quad meshes (hybrid-mesh plan, Task 7),
    ! derived from FESOM2 v2.7.3 oce_mesh.F90 and the qq prototype (oce_mesh.F90
    ! elem_center / mesh_arrays). Every element has 4 slots; a triangle repeats its 1st
    ! vertex as the 4th (elem2D_nnodes == 3). Unified formulas (radians with elem_cos, then
    ! r_earth; cyclic wrap on all vertices):
    !   centroid       area-weighted centroid of the sub-triangles (1,3,4) and (1,2,3)
    !   elem_area      their summed lon/lat area * elem_cos (cos of the centroid latitude)
    !   elem_w         kite weights 1/2*[A(j,j-1,c) + A(j,j+1,c)] / A (triangle 1/3,1/3,1/3,0)
    !   area           sum over elements and slots of elem_area*elem_w(k)
    !   gradient_sca   Green-Gauss: vertex j ~ (x_{j+1}-x_{j-1})/(2A) over the nv vertices,
    !                  x in 1:4, y in 5:8 (triangle slots 4/8 = 0)
    ! For triangles all of them equal FESOM2's P1 formulas up to roundoff (not bitwise:
    ! the FESOM2-oracle geometry byte-gate is retired, docs/HANDOFF.md).
    !
    ! Provides elem_cos, metric_factor, elem_area, elem_w, nod_in_elem2D_w, gradient_sca,
    ! edge_dxdy, edge_cross_dxdy, elem_edges, elem_neighbors, area/areasvol(+inv),
    ! coriolis and mesh_resolution. FESOM3 bottom at vertices: area/areasvol(+inv) are 1-D
    ! (per vertex) -- a scalar column's horizontal area does not vary with depth. See the
    ! contract in mod_mesh.F90. gradient_vec is still deferred.
    !
    ! ORDER (kept from FESOM2 mesh_areas + mesh_auxiliary_arrays): elem_area and the node
    ! areas accumulate UNSCALED, then one *r_earth^2 on elem_area/area/areasvol together,
    ! then area_inv/areasvol_inv; gradient_sca uses the SCALED elem_area.
    !
    ! MULTI-RANK (M2.12a). The geometry is partition-agnostic in arithmetic; only the
    ! loop bounds + allocation sizes become LOCAL. local_bounds() returns mesh global
    ! counts at npes==1 (the proven 1-rank path is byte-for-byte unchanged) and the
    ! partit myDim/eDim dims at npes>1. Element CENTERS are precomputed for OWNED
    ! elements and halo-exchanged, because edge_cross_dxdy of an owned edge can read the
    ! center of a halo (eDim) neighbour element whose elem2D_nodes is not stored (FESOM2
    ! mesh_auxiliary_arrays:2526-2529 + 2552-2578). Owned-node areas are computed locally
    ! (the partition guarantees an owned node's full element-neighbourhood is owned), so
    ! the area/elem_area halo exchanges (FESOM2 mesh_areas:2220/2322) are NOT needed for
    ! the owned-entry gate and are deferred to M2.12b (where the dynamics consume halos).
    use mod_precision,   only: WP, MP
    use mod_constants,   only: r_earth, omega, pi
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_partit,      only: t_partit
    use mod_halo,        only: exchange_elem, exchange_elem_full, exchange_nod, allreduce_sum
    use mod_mesh_rotate, only: trim_cyclic, get_cyclic_length, r2g
    implicit none
    private
    public :: compute_geometry, elem_center

contains

    subroutine local_bounds(mesh, partit, nNodO, nNodL, nElemO, nElemF, nEdgeO, nEdgeL)
        ! Local array bounds. npes==1: the global mesh counts (proven path, unchanged).
        ! npes>1: owned (myDim) and local (myDim+eDim[+eXDim]) dims from the partition.
        type(t_mesh),   intent(in)  :: mesh
        type(t_partit), intent(in)  :: partit
        integer,        intent(out) :: nNodO, nNodL, nElemO, nElemF, nEdgeO, nEdgeL
        if (partit%npes == 1) then
            nNodO  = mesh%nod2D;  nNodL  = mesh%nod2D
            nElemO = mesh%elem2D; nElemF = mesh%elem2D
            nEdgeO = mesh%edge2D; nEdgeL = mesh%edge2D
        else
            nNodO  = partit%myDim_nod2D
            nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
            nElemO = partit%myDim_elem2D
            nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
            nEdgeO = partit%myDim_edge2D
            nEdgeL = partit%myDim_edge2D + partit%eDim_edge2D
        end if
    end subroutine local_bounds

    subroutine compute_geometry(mesh, partit, cartesian)
        ! Same data-dependency order as FESOM2 mesh_areas + mesh_auxiliary_arrays:
        ! elem_cos before edge/gradient; elem_area UNSCALED before area accumulation;
        ! elem_area scaled (in compute_node_areas) before gradient_sca.
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        logical,        intent(in)    :: cartesian
        integer :: nNodO, nNodL, nElemO, nElemF, nEdgeO, nEdgeL
        real(kind=WP), allocatable :: center_x(:), center_y(:), sarea(:)
        call local_bounds(mesh, partit, nNodO, nNodL, nElemO, nElemF, nEdgeO, nEdgeL)
        allocate(center_x(nElemF), center_y(nElemF), sarea(nElemO))
        center_x = 0.0_WP; center_y = 0.0_WP
        call compute_centers(mesh, nElemO, center_x, center_y, sarea)   ! owned centroids + lon/lat areas
        call compute_elem_metric(mesh, cartesian, nElemO, nElemF, center_y) ! elem_cos, metric_factor
        call compute_coriolis(mesh, cartesian, nNodL, nElemO)           ! coriolis, coriolis_node (M2.3)
        call compute_elem_area(mesh, nElemO, nElemF, sarea)             ! elem_area (UNSCALED radians^2)
        call compute_elem_weights(mesh, nElemO, center_x, center_y)     ! elem_w (kite weights)
        call compute_nod_weights(mesh, nNodO, nNodL, partit)            ! nod_in_elem2D_w
        if (partit%npes > 1) then
            ! Halo element centers for owned-edge edge_cross_dxdy (FESOM2:2528-2529).
            call exchange_elem(center_x, partit)
            call exchange_elem(center_y, partit)
            call exchange_elem_cos(mesh, partit, nElemF)
        end if
        call compute_node_areas(mesh, nNodO, nNodL, partit)             ! accumulate area, then scale
        call build_elem_adjacency(mesh, nElemO, nEdgeL)                 ! elem_edges, elem_neighbors
        if (partit%npes > 1) then
            ! M2.12b: halo elem_area (FULL halo; FESOM2 used it for the area weighting of
            ! fill_up_dn_grad at halo elements reached through a halo node's element list),
            ! and halo node areas owner->halo (FESOM2 mesh_areas:2220 + 2322-2323). The
            ! accumulation stays owned-only (an owned node's element neighbourhood is
            ! complete) so OWNED area is unchanged; the exchange fills only the halo
            ! entries (exchanging the SCALED arrays == FESOM2's scale-then-broadcast).
            call exchange_elem_full(mesh%elem_area, partit)
            call exchange_nod(mesh%area, partit)
            call exchange_nod(mesh%areasvol, partit)
            call exchange_nod(mesh%area_inv, partit)
            call exchange_nod(mesh%areasvol_inv, partit)
        end if
        call compute_mesh_resolution(mesh, nNodO, nNodL, partit)        ! M4 GM: scalar cell resolution
        call compute_edge_geometry(mesh, nEdgeO, center_x, center_y)    ! edge_dxdy, edge_cross_dxdy
        call compute_gradient_sca(mesh, nElemO)                         ! uses SCALED elem_area + elem_cos
        deallocate(center_x, center_y, sarea)
    end subroutine compute_geometry

    subroutine exchange_elem_cos(mesh, partit, nElemF)
        ! Halo-fill elem_cos via a WP scratch (elem_cos is MP; MP==WP at dp/sp).
        type(t_mesh),   intent(inout) :: mesh
        type(t_partit), intent(in)    :: partit
        integer,        intent(in)    :: nElemF
        real(kind=WP), allocatable :: tmp(:)
        allocate(tmp(nElemF))
        tmp = real(mesh%elem_cos(1:nElemF), WP)
        call exchange_elem(tmp, partit)
        mesh%elem_cos(1:nElemF) = real(tmp, MP)
        deallocate(tmp)
    end subroutine exchange_elem_cos

    subroutine compute_centers(mesh, nElemO, cx, cy, sa)
        ! Element centroids and lon/lat areas (radians^2, no cosine) for OWNED elements
        ! (elem_center needs elem2D_nodes, stored owned-only). The halo centroids are
        ! filled by exchange in compute_geometry.
        type(t_mesh),  intent(in)  :: mesh
        integer,       intent(in)  :: nElemO
        real(kind=WP), intent(out) :: cx(:), cy(:), sa(:)
        integer :: n
        real(kind=WP) :: ax, ay, s
        do n = 1, nElemO
            call elem_center(mesh, n, ax, ay, s)
            cx(n) = ax; cy(n) = ay; sa(n) = s
        end do
    end subroutine compute_centers

    !--------------------------------------------------------------------------
    subroutine compute_coriolis(mesh, cartesian, nNodL, nElemO)
        ! Coriolis parameter f = 2*omega*sin(lat_geo) at elements (coriolis) and nodes
        ! (coriolis_node). FESOM2 mesh_auxiliary_arrays (oce_mesh.F90:2476-2503): the
        ! geographical latitude is r2g applied to the ROTATED element centroid
        ! (elem_center) and to the rotated node coordinate. r2g uses the same rotation
        ! matrix that g2r built coord_nod2D with (geometry-gate-proven), so coriolis
        ! byte-matches FESOM2 by construction (L9 transitive-gate pattern). M2.3 needs
        ! coriolis (elements); coriolis_node (nodes) is faithful + cheap, used later.
        type(t_mesh), intent(inout) :: mesh
        logical,      intent(in)    :: cartesian
        integer,      intent(in)    :: nNodL, nElemO
        integer :: n
        real(kind=WP) :: ax, ay, lon, lat
        allocate(mesh%coriolis(nElemO), mesh%coriolis_node(nNodL))
        if (.not. cartesian) then
            do n = 1, nNodL
                call r2g(lon, lat, mesh%coord_nod2D(1, n), mesh%coord_nod2D(2, n))
                mesh%coriolis_node(n) = 2 * omega * sin(lat)
            end do
            do n = 1, nElemO
                call elem_center(mesh, n, ax, ay)
                call r2g(lon, lat, ax, ay)
                mesh%coriolis(n) = 2 * omega * sin(lat)
            end do
        else
            ! cartesian/analytic mesh (no rotated->geo transform): use the stored
            ! latitude directly. coriolis is NOT byte-gated on the analytic mesh; this
            ! is a benign finite fill that avoids r2g/asin on cartesian coords.
            do n = 1, nNodL
                mesh%coriolis_node(n) = 2 * omega * sin(mesh%coord_nod2D(2, n))
            end do
            do n = 1, nElemO
                call elem_center(mesh, n, ax, ay)
                mesh%coriolis(n) = 2 * omega * sin(ay)
            end do
        end if
    end subroutine compute_coriolis

    !--------------------------------------------------------------------------
    pure subroutine elem_center(mesh, n, cx, cy, s)
        ! Centroid of element n (qq oce_mesh.F90 elem_center): the area-weighted mean of
        ! the centroids of the sub-triangles (1,3,4) and (1,2,3), in (rotated) lon/lat.
        ! Each vertex longitude is first wrapped relative to the element's MINIMUM
        ! longitude (FESOM2 oce_mesh.F90:2135-2155: >= / < against cyclic_length/2, NOT
        ! the symmetric trim_cyclic). s (optional) = the element's lon/lat area (radians^2,
        ! no cosine) = the sum of the two sub-triangle areas. A triangle repeats vertex 1
        ! as vertex 4, so its (1,3,4) part has area exactly 0 and (cx,cy) is the vertex
        ! mean up to roundoff. Public: also the element coordinate of the output stores
        ! (mod_io_coords).
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: n
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP), intent(out), optional :: s
        integer :: j
        real(kind=WP) :: ax(MAX_NV), ay(MAX_NV), amin, cl, s1, s2
        do j = 1, MAX_NV
            ax(j) = mesh%coord_nod2D(1, mesh%elem2D_nodes(j, n))
            ay(j) = mesh%coord_nod2D(2, mesh%elem2D_nodes(j, n))
        end do
        amin = minval(ax)
        cl = get_cyclic_length()
        do j = 1, MAX_NV
            if (ax(j) - amin >=  cl/2.0_WP) ax(j) = ax(j) - cl
            if (ax(j) - amin <  -cl/2.0_WP) ax(j) = ax(j) + cl
        end do
        s1 = abs((ax(3)-ax(4))*(ay(1)-ay(4)) - (ay(3)-ay(4))*(ax(1)-ax(4)))   ! 2*A(1,3,4)
        s2 = abs((ax(3)-ax(2))*(ay(1)-ay(2)) - (ay(3)-ay(2))*(ax(1)-ax(2)))   ! 2*A(1,2,3)
        cx = ((ax(1)+ax(3)+ax(4))*s1 + (ax(1)+ax(2)+ax(3))*s2) / (3.0_WP*(s1+s2))
        cy = ((ay(1)+ay(3)+ay(4))*s1 + (ay(1)+ay(2)+ay(3))*s2) / (3.0_WP*(s1+s2))
        if (present(s)) s = 0.5_WP*(s1+s2)
    end subroutine elem_center

    pure subroutine edge_center(mesh, e, cx, cy)
        ! FESOM2 oce_mesh.F90:2115-2133: a = coord(n1), b = coord(n2); shift a(1) DOWN
        ! if a-b too positive, shift b(1) DOWN if a-b too negative; cx = 0.5*(a1+b1).
        type(t_mesh), intent(in)  :: mesh
        integer,      intent(in)  :: e
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP) :: a1, a2, b1, b2, cl
        a1 = mesh%coord_nod2D(1, mesh%edges(1, e)); a2 = mesh%coord_nod2D(2, mesh%edges(1, e))
        b1 = mesh%coord_nod2D(1, mesh%edges(2, e)); b2 = mesh%coord_nod2D(2, mesh%edges(2, e))
        cl = get_cyclic_length()
        if (a1 - b1 >  cl/2.0_WP) a1 = a1 - cl
        if (a1 - b1 < -cl/2.0_WP) b1 = b1 - cl
        cx = 0.5_WP * (a1 + b1)
        cy = 0.5_WP * (a2 + b2)
    end subroutine edge_center

    !--------------------------------------------------------------------------
    subroutine compute_elem_metric(mesh, cartesian, nElemO, nElemF, cy)
        ! elem_cos = cos(lat_center); metric_factor = tan(lat_center)/r_earth, using
        ! the precomputed owned centers cy. (oce_mesh.F90:2508-2528)
        type(t_mesh),  intent(inout) :: mesh
        logical,       intent(in)    :: cartesian
        integer,       intent(in)    :: nElemO, nElemF
        real(kind=WP), intent(in)    :: cy(:)
        integer :: n
        allocate(mesh%elem_cos(nElemF), mesh%metric_factor(nElemF))
        mesh%elem_cos = 0.0_MP; mesh%metric_factor = 0.0_MP
        do n = 1, nElemO
            mesh%elem_cos(n)      = cos(cy(n))
            mesh%metric_factor(n) = tan(cy(n)) / r_earth
        end do
        if (cartesian) then
            mesh%elem_cos = 1.0_MP; mesh%metric_factor = 0.0_MP
        end if
    end subroutine compute_elem_metric

    subroutine compute_elem_area(mesh, nElemO, nElemF, sa)
        ! elem_area = lon/lat area (elem_center) * elem_cos, i.e. the planar area in
        ! (lon*cos(lat_c), lat) with lat_c the centroid latitude (elem_cos = 1 for a
        ! cartesian mesh). This is the area the Green-Gauss gradient needs to be exact for
        ! linear fields. elem_area is left UNSCALED here (radians^2); the *r_earth^2
        ! happens in compute_node_areas (FESOM2's deferred scaling, mesh_areas).
        type(t_mesh),  intent(inout) :: mesh
        integer,       intent(in)    :: nElemO, nElemF
        real(kind=WP), intent(in)    :: sa(:)
        integer :: n
        allocate(mesh%elem_area(nElemF))
        mesh%elem_area = 0.0_MP
        do n = 1, nElemO
            mesh%elem_area(n) = sa(n) * mesh%elem_cos(n)
        end do
    end subroutine compute_elem_area

    subroutine compute_elem_weights(mesh, nElemO, cx, cy)
        ! Kite (median-dual) weights of OWNED elements: the share of vertex j is half the
        ! areas of the triangles (j, j-1, c) and (j, j+1, c), c = the centroid, i.e. the
        ! kite between vertex j, its two mid-edges and c (qq w_cv). Normalised by the sum
        ! of the four side triangles (= the element area for a convex element), so the
        ! weights sum to 1 to roundoff. Lon/lat areas: the cosine factor cancels. Triangles
        ! take the exact 1/3 (slot 4 = 0).
        type(t_mesh),  intent(inout) :: mesh
        integer,       intent(in)    :: nElemO
        real(kind=WP), intent(in)    :: cx(:), cy(:)
        integer :: e, j, jp
        real(kind=WP) :: dx(MAX_NV), dy(MAX_NV), t(MAX_NV)
        allocate(mesh%elem_w(MAX_NV, nElemO))
        mesh%elem_w = 0.0_MP
        do e = 1, nElemO
            if (mesh%elem2D_nnodes(e) == 3) then
                mesh%elem_w(1:3, e) = 1.0_WP/3.0_WP
                cycle
            end if
            do j = 1, MAX_NV
                dx(j) = mesh%coord_nod2D(1, mesh%elem2D_nodes(j, e)) - cx(e)
                call trim_cyclic(dx(j))
                dy(j) = mesh%coord_nod2D(2, mesh%elem2D_nodes(j, e)) - cy(e)
            end do
            do j = 1, MAX_NV                           ! t(j) = 2*A(j, j+1, c)
                jp = mod(j, MAX_NV) + 1
                t(j) = abs(dx(j)*dy(jp) - dy(j)*dx(jp))
            end do
            do j = 1, MAX_NV
                mesh%elem_w(j, e) = 0.5_WP*(t(mod(j+2, MAX_NV)+1) + t(j)) / sum(t)
            end do
        end do
    end subroutine compute_elem_weights

    subroutine compute_nod_weights(mesh, nNodO, nNodL, partit)
        ! nod_in_elem2D_w(k,n) = elem_w of node n in its k-th element nod_in_elem2D(k,n)
        ! (the slot among the element's nv vertices). Owned nodes from their (owned)
        ! elements; halo nodes by exchange -- complete_nod_in_elem_halo copied the owner's
        ! slot order, so the owner's weights are the halo node's.
        type(t_mesh),   intent(inout) :: mesh
        integer,        intent(in)    :: nNodO, nNodL
        type(t_partit), intent(in)    :: partit
        integer :: n, j, k, e
        allocate(mesh%nod_in_elem2D_w(size(mesh%nod_in_elem2D, 1), nNodL))
        mesh%nod_in_elem2D_w = 0.0_MP
        do n = 1, nNodO
            do j = 1, mesh%nod_in_elem2D_num(n)
                e = mesh%nod_in_elem2D(j, n)
                do k = 1, mesh%elem2D_nnodes(e)
                    if (mesh%elem2D_nodes(k, e) == n) then
                        mesh%nod_in_elem2D_w(j, n) = mesh%elem_w(k, e)
                        exit
                    end if
                end do
            end do
        end do
        if (partit%npes > 1) call exchange_nod(mesh%nod_in_elem2D_w, partit)
    end subroutine compute_nod_weights

    subroutine compute_gradient_sca(mesh, nElemO)
        ! Green-Gauss gradient coefficients over the element's nv vertices (qq
        ! mesh_arrays; FESOM2's P1 coefficients for a triangle, up to roundoff): with the
        ! vertices clockwise, vertex j gets (-(y_{j+1}-y_{j-1}), (x_{j+1}-x_{j-1})*elem_cos)
        ! * r_earth/(2*elem_area). x in slots 1:4, y in 5:8; a triangle's slots 4/8 stay 0.
        ! The longitudes are taken relative to vertex 1 (trim_cyclic, as FESOM2), so all
        ! differences come from one planar embedding of the element. Uses the SCALED
        ! elem_area. Owned elements.
        type(t_mesh), intent(inout) :: mesh
        integer,      intent(in)    :: nElemO
        integer :: e, j, nv, jp, jm
        real(kind=WP) :: x(MAX_NV), y(MAX_NV), dfac
        allocate(mesh%gradient_sca(2*MAX_NV, nElemO))   ! (2*MAX_NV, owned elems)
        mesh%gradient_sca = 0.0_MP
        do e = 1, nElemO
            nv   = mesh%elem2D_nnodes(e)
            dfac = 0.5_WP * r_earth / mesh%elem_area(e)
            do j = 1, nv
                x(j) = mesh%coord_nod2D(1, mesh%elem2D_nodes(j, e)) &
                     - mesh%coord_nod2D(1, mesh%elem2D_nodes(1, e))
                call trim_cyclic(x(j))
                y(j) = mesh%coord_nod2D(2, mesh%elem2D_nodes(j, e))
            end do
            do j = 1, nv
                jp = mod(j, nv) + 1
                jm = mod(j+nv-2, nv) + 1
                mesh%gradient_sca(j,          e) = -(y(jp) - y(jm)) * dfac
                mesh%gradient_sca(j + MAX_NV, e) =  mesh%elem_cos(e) * (x(jp) - x(jm)) * dfac
            end do
        end do
    end subroutine compute_gradient_sca

    subroutine build_elem_adjacency(mesh, nElemO, nEdgeL)
        ! mesh%elem_edges(k, e) / mesh%elem_neighbors(k, e): the edge k of OWNED element
        ! e, connecting slots k and k+1 (1-2, 2-3, 3-4, 4-1 over the nv vertices), and the
        ! element across it -- e ITSELF on a domain boundary (no zero index). A triangle's
        ! 4th entry repeats its 1st (as elem2D_nodes). Built ONCE from edges/edge_tri: the
        ! neighbour across an edge is "the other edge_tri entry", which on a boundary edge
        ! (edge_tri(2) == edge_tri(1)) is the element itself.
        ! Neighbour ids may be HALO element ids, and then ALWAYS from the eDim ring: by
        ! construction (FESOM2 gen_comm.F90) com_elem2D holds exactly the non-owned elements
        ! sharing an EDGE with an owned one, while eXDim (com_elem2D_full minus com_elem2D)
        ! holds node-only neighbours, which are never across an edge. So a neighbour's value
        ! is readable from any element field after the STANDARD exchange_elem (the
        ! rvo_upwind omega_e consumer does that); exchange_elem_full is NOT needed.
        ! test_vinv A5 asserts this with the production fill/exchange/read pattern.
        !
        ! The sweep covers ALL local edges (owned + halo): an owned element can have ONE
        ! halo edge (both endpoints halo nodes), carried only by the eDim part of the edge
        ! list. NOTE this is also why no per-edge GEOMETRY can be used for owned elements:
        ! edge_dxdy/edge_cross_dxdy exist for OWNED edges only (see their allocation), so
        ! consumers rebuild edge vectors from the element's own vertices instead.
        ! The -1 init + post-sweep check turns the partition assumption "an owned element
        ! sees all of its edges locally" into a hard invariant (same style as
        ! assert_bottom_invariant).
        type(t_mesh), intent(inout) :: mesh
        integer,      intent(in)    :: nElemO, nEdgeL
        integer :: ed, s, e, k, nv, a, b, v1, v2

        allocate(mesh%elem_neighbors(MAX_NV, nElemO), mesh%elem_edges(MAX_NV, nElemO))
        mesh%elem_neighbors = -1
        mesh%elem_edges     = -1
        do ed = 1, nEdgeL
            a = mesh%edges(1, ed)
            b = mesh%edges(2, ed)
            do s = 1, 2
                e = mesh%edge_tri(s, ed)
                if (s == 2 .and. e == mesh%edge_tri(1, ed)) cycle   ! boundary edge: one element
                if (e < 1 .or. e > nElemO) cycle                     ! halo element
                nv = mesh%elem2D_nnodes(e)
                do k = 1, nv
                    v1 = mesh%elem2D_nodes(k, e)
                    v2 = mesh%elem2D_nodes(mod(k,nv)+1, e)
                    if ((v1 == a .and. v2 == b) .or. (v1 == b .and. v2 == a)) then
                        mesh%elem_neighbors(k, e) = mesh%edge_tri(3-s, ed)   ! = e on a boundary
                        mesh%elem_edges(k, e)     = ed
                        exit
                    end if
                end do
            end do
        end do
        do e = 1, nElemO
            nv = mesh%elem2D_nnodes(e)
            if (any(mesh%elem_neighbors(1:nv, e) < 0)) then
                write(*,'(a,i0)') 'build_elem_adjacency: unresolved edge slot at owned element ', e
                error stop 'build_elem_adjacency: local edge list does not cover an owned element'
            end if
            if (nv == 3) then                                          ! triangle: 4th = 1st
                mesh%elem_neighbors(4, e) = mesh%elem_neighbors(1, e)
                mesh%elem_edges(4, e)     = mesh%elem_edges(1, e)
            end if
        end do
    end subroutine build_elem_adjacency

    subroutine compute_edge_geometry(mesh, nEdgeO, center_x, center_y)
        ! edge_dxdy (along-edge) + edge_len + edge_cross_dxdy (edge-center to elem
        ! centers). ALL THREE IN METRES. Uses the precomputed (and, at npes>1,
        ! halo-exchanged) element centers so an owned edge with a halo (eDim) neighbour
        ! element resolves its center without elem_center on a halo element.
        !
        ! R7 (FESOM3): FESOM2 stored edge_dxdy in RADIAN measure and multiplied by
        ! r_earth * mean(elem_cos over the edge's elements) at each point of use --
        ! oce_adv_tra_hor.F90 computed exactly that `a` inline. FESOM3 folds the factor in
        ! here, so edge_dxdy is a physical length and the consumers just use it. The mean
        ! cosine is over the two elements adjacent to the edge, or the one available at a
        ! boundary edge, which is why the two loops are merged: the radian pass had no
        ! el1/el2 in hand.
        !
        ! No `cartesian` special case is needed: that mode sets elem_cos = 1 while still
        ! scaling everything else by r_earth, so the same expression is consistent.
        type(t_mesh),  intent(inout) :: mesh
        integer,       intent(in)    :: nEdgeO
        real(kind=WP), intent(in)    :: center_x(:), center_y(:)
        integer :: n, el1, el2
        real(kind=WP) :: a1, a2, ecx, ecy, b1, b2, cosm
        allocate(mesh%edge_dxdy(2, nEdgeO), mesh%edge_cross_dxdy(4, nEdgeO))
        allocate(mesh%edge_len(nEdgeO))
        do n = 1, nEdgeO
            call edge_center(mesh, n, ecx, ecy)
            el1 = mesh%edge_tri(1, n); el2 = mesh%edge_tri(2, n)

            ! along-edge separation, in metres
            a1 = mesh%coord_nod2D(1, mesh%edges(2, n)) - mesh%coord_nod2D(1, mesh%edges(1, n))
            a2 = mesh%coord_nod2D(2, mesh%edges(2, n)) - mesh%coord_nod2D(2, mesh%edges(1, n))
            call trim_cyclic(a1)
            cosm = mesh%elem_cos(el1)
            if (el2 /= el1) cosm = 0.5_WP*(cosm + mesh%elem_cos(el2))
            mesh%edge_dxdy(1, n) = a1 * cosm * r_earth                ! [m]
            mesh%edge_dxdy(2, n) = a2 * r_earth                       ! [m]
            mesh%edge_len(n)     = sqrt(real(mesh%edge_dxdy(1,n),WP)**2 &
                                      + real(mesh%edge_dxdy(2,n),WP)**2)   ! [m]
            b1 = center_x(el1) - ecx; b2 = center_y(el1) - ecy; call trim_cyclic(b1)
            b1 = b1 * mesh%elem_cos(el1)
            mesh%edge_cross_dxdy(1, n) = b1 * r_earth
            mesh%edge_cross_dxdy(2, n) = b2 * r_earth
            if (el2 /= el1) then
                b1 = center_x(el2) - ecx; b2 = center_y(el2) - ecy; call trim_cyclic(b1)
                b1 = b1 * mesh%elem_cos(el2)
                mesh%edge_cross_dxdy(3, n) = b1 * r_earth
                mesh%edge_cross_dxdy(4, n) = b2 * r_earth
            else
                mesh%edge_cross_dxdy(3, n) = 0.0_WP; mesh%edge_cross_dxdy(4, n) = 0.0_WP
            end if
        end do
    end subroutine compute_edge_geometry

    subroutine compute_node_areas(mesh, nNodO, nNodL, partit)
        ! Control-volume area per level. The accumulation runs on UNSCALED elem_area;
        ! then elem_area, area and areasvol are all multiplied by r_earth^2 together (a
        ! single deferred scaling, as in FESOM2 mesh_areas:2313-2315) so the per-node
        ! sums round identically. Owned nodes are accumulated locally (the partition
        ! guarantees a complete owned element-neighbourhood); the halo exchange is done
        ! by the caller.
        !
        ! FESOM3 BOTTOM AT VERTICES: the scalar cell of vertex n is a STRAIGHT PRISM --
        ! one bottom level, and the FULL median-dual area at every wet layer. So area is
        ! depth-independent over the vertex's wet range, which is the design note's
        ! `area(1:myDim+eDim)`.
        !
        ! FESOM2 instead gathered elem_area/3 only from elements deep enough to reach
        ! level nz (oce_mesh.F90:2252-2351), so the cell narrowed with depth as its
        ! elements bottomed out. That is exactly what a vertex-defined bottom removes:
        ! with nlevels(e) = min over the element's nodes, every adjacent element is at
        ! most as deep as the node, and a depth-gathered area would shrink to the
        ! deepest element's share alone rather than the cell's true area.
        !
        ! The array stays 2-D. The entry at nz = nlevels_nod2D(n) is deliberately left
        ! ZERO: area(nz,n) doubles as the INTERFACE area at level nz, so a zero there is
        ! a closed bottom. Nothing depends on it today -- cal_shortwave_rad already
        ! forces sw_3d(nzmax,n)=0 (oce_shortwave_pene.F90:77-81) and the vertical
        ! advection hard-zeroes its own bottom flux (oce_adv_tra_ver.F90:68-69, :120-121)
        ! -- but it costs nothing and keeps the idiom available.
        type(t_mesh),   intent(inout) :: mesh
        integer,        intent(in)    :: nNodO, nNodL
        type(t_partit), intent(in)    :: partit
        integer       :: n, j
        real(kind=MP) :: acell
        allocate(mesh%area(nNodL), mesh%area_inv(nNodL))
        allocate(mesh%areasvol(nNodL), mesh%areasvol_inv(nNodL))
        mesh%area = 0.0_MP
        do n = 1, nNodO
            ! full median-dual area: EVERY adjacent element (its kite share), no depth test.
            acell = 0.0_MP
            do j = 1, mesh%nod_in_elem2D_num(n)
                acell = acell + mesh%elem_area(mesh%nod_in_elem2D(j, n)) * mesh%nod_in_elem2D_w(j, n)
            end do
            mesh%area(n) = acell
        end do
        ! non-cavity: the scalar-volume area IS the control-volume area
        mesh%areasvol = mesh%area
        ! deferred single scaling to physical m^2 (mesh_areas:2313-2315)
        mesh%elem_area = mesh%elem_area * r_earth * r_earth
        mesh%area      = mesh%area      * r_earth * r_earth
        mesh%areasvol  = mesh%areasvol  * r_earth * r_earth
        mesh%area_inv = 0.0_MP
        do n = 1, nNodO
            if (mesh%area(n) > 0.0_MP) mesh%area_inv(n) = 1.0_MP / mesh%area(n)
        end do
        mesh%areasvol_inv = mesh%area_inv
        ! ocean_area / ocean_areawithcav: faithful FESOM2 oce_mesh.F90:2385 sequential
        ! accumulation over areasvol (cavity-aware), NOT sum(area(1,:)). ocean_area is the
        ! divisor in the M3e oce_fluxes flux balancing (net/ocean_area), so its summation
        ! order must match FESOM2 bit-for-bit (-fp-model precise stops reassociation, but a
        ! sum() intrinsic and a do-loop are not guaranteed identical — the explicit loop is).
        ! At 1-rank this local partial sum IS the global value; at npes>1 (M3f-4) the owned
        ! partial sums are summed across ranks by allreduce_sum (FESOM2 oce_mesh.F90:2389
        ! MPI_AllREDUCE(vol/vol2, MPI_SUM, MPI_DOUBLE_PRECISION) — the reduction is over the
        ! SAME owned-node partial sums in the same comm, so it is byte-identical, L6/L33).
        block
            real(kind=MP) :: vol, vol2
            real(kind=WP) :: gvol, gvol2
            vol = 0.0_MP; vol2 = 0.0_MP
            do n = 1, nNodO
                vol2 = vol2 + mesh%areasvol(n)
                if (mesh%ulevels_nod2D(n) > 1) cycle
                vol  = vol + mesh%areasvol(n)
            end do
            if (partit%npes > 1) then
                gvol = real(vol, WP); gvol2 = real(vol2, WP)
                call allreduce_sum(gvol,  partit)     ! ocean_area      (surface, no cavity)
                call allreduce_sum(gvol2, partit)     ! ocean_areawithcav (incl. cavity)
                mesh%ocean_area        = real(gvol,  MP)
                mesh%ocean_areawithcav = real(gvol2, MP)
            else
                mesh%ocean_area        = vol
                mesh%ocean_areawithcav = vol2
            end if
        end block
    end subroutine compute_node_areas

    subroutine compute_mesh_resolution(mesh, nNodO, nNodL, partit)
        ! Scalar cell resolution (oce_mesh.F90:2358-2383). Raw resolution
        ! 2*sqrt(areasvol(ulevel)/pi) on OWNED+HALO nodes, then 3 mass-matrix smoothing
        ! sweeps: area-weighted neighbour average over nod_in_elem2D, each followed by
        ! exchange_nod so the next sweep's halo reads are valid (FESOM2:2381). Runs after
        ! compute_node_areas + its halo exchange (so areasvol is SCALED and halo-valid).
        ! Consumed by the M4 GM/Redi K scaling (init_Redi_GM, oce_fer_gm.F90:258). The
        ! element mean is the kite-weighted one over the 4 slots (FESOM2:2373 had /3).
        type(t_mesh),   intent(inout) :: mesh
        integer,        intent(in)    :: nNodO, nNodL
        type(t_partit), intent(in)    :: partit
        integer :: n, j, q, elem
        integer :: elnodes(MAX_NV)
        real(kind=WP) :: vol, acc
        real(kind=WP), allocatable :: work_array(:)
        allocate(mesh%mesh_resolution(nNodL))
        do n = 1, nNodL
            mesh%mesh_resolution(n) = sqrt(mesh%areasvol(n) / pi) * 2.0_WP
        end do
        allocate(work_array(nNodO))
        do q = 1, 3                                     ! apply mass matrix 3x to smooth
            do n = 1, nNodO
                vol = 0.0_WP
                acc = 0.0_WP
                do j = 1, mesh%nod_in_elem2D_num(n)
                    elem = mesh%nod_in_elem2D(j, n)
                    elnodes = mesh%elem2D_nodes(:, elem)
                    acc = acc + sum(mesh%elem_w(:, elem)*mesh%mesh_resolution(elnodes)) * mesh%elem_area(elem)
                    vol = vol + mesh%elem_area(elem)
                end do
                work_array(n) = acc / vol
            end do
            do n = 1, nNodO
                mesh%mesh_resolution(n) = work_array(n)
            end do
            if (partit%npes > 1) call exchange_nod(mesh%mesh_resolution, partit)
        end do
        deallocate(work_array)
    end subroutine compute_mesh_resolution

end module mod_mesh_areas
