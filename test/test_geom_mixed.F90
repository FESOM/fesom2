program test_geom_mixed
    ! Mixed-mesh geometry (plan 2026-10-08-hybrid-meshes, Task 7): compute_geometry on a
    ! mesh read from FESOM3_MESH_DIR (default the NEMO-mixed cut-out; pi is registered
    ! too, rotated) at 1 or N ranks and, at 1 rank, on the analytic mixed patch
    ! (generate_mixed_patch: jittered quads + triangle pairs). Identities:
    !   G1 Sum node area (owned, all ranks) == Sum element area (each element once)
    !   G2 Sum elem_w == 1 per element; slots 1..nv > 0, a triangle's slot 4 == 0
    !   G3 Sum_k nod_in_elem2D_w(k,n)*elem_area(e_k) == area(n) (owned and halo nodes)
    !   G4 gradient_sca exact for linear fields in the element metric (lon*cos, lat),
    !      triangles and quads, all 4 slots
    !   G5 cross-edge vectors close around each element: with d_k = edge_cross_dxdy of
    !      the element at its edge k (edge centre -> centroid), elem_area*gradient_sca of
    !      vertex j == rot(d_j - d_{j-1}) -- the Green-Gauss coefficients ARE the dual
    !      faces, which telescope to zero around the element
    !   G6 conventions: edge_tri has no zero, a boundary edge stores its element twice
    !      with edge_cross_dxdy(3:4) == 0; elem_edges(k) joins slots k, k+1;
    !      elem_neighbors(k) is the other element of that edge (the element itself on a
    !      boundary); a triangle's slot 4 repeats slot 1
    !   G7 triangles: centroid (elem_cos), elem_area, gradient_sca and the node area of
    !      all-triangle nodes equal FESOM2's P1 formulas (the previous FESOM3 code) to
    !      roundoff
    !   G8 quads: elem_w == the kite areas (vertex, its two mid-edges, centroid) by the
    !      shoelace formula / the shoelace polygon area, and elem_center == the shoelace
    !      polygon centroid (independent of the sub-triangle code); on the jittered patch
    !      the weights must differ from 1/4 (no vacuous pass)
    !   G9 compute_mesh_resolution (np 1, parallelogram patch): 2*sqrt(areasvol/pi) after
    !      the 3 kite-weighted smoothing sweeps equals the raw value of the uniform cells at
    !      the nodes 4+ rings from the rim (their 3-ring neighbourhood has no rim node), and
    !      the rim-near nodes differ (non-vacuity)
    !   G10 output element coordinates (mod_io_coords io_coords_compute, the lon/lat of every
    !      element store): the geographic image (r2g) of the shoelace polygon centroid in the
    !      rotated lon/lat plane, to 1e-9 degrees (on the rotated pi_w across the dateline too)
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth, rad
    use mod_mesh,          only: t_mesh, MAX_NV
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init
    use mod_mesh_areas,    only: elem_center
    use mod_mesh_analytic, only: generate_mixed_patch, generate_quad_patch
    use mod_mesh_rotate,   only: trim_cyclic, get_cyclic_length, r2g
    use mod_io_coords,     only: io_coords_compute
    use mod_io_decomp,     only: DECOMP_ELEM
    use mod_part_bounds,   only: owned_bounds
    implicit none

    type(t_partit) :: partit, apartit, qpartit
    type(t_mesh)   :: mesh, amesh, qmesh

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.false., partition_np1=.true.)
    call check_mesh(mesh, partit, 'mesh', .false.)
    call check_output_coords(mesh, partit)

    if (partit%npes == 1) then
        ! analytic mixed patch (1-rank identity partition of its own)
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=4, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check(count(amesh%elem2D_nnodes == 4) > 0 .and. count(amesh%elem2D_nnodes == 3) > 0, &
                   'patch: triangles and quads', 0)
        call check(patch_convex(amesh), 'patch: every element clockwise and convex', 0)
        call check(amesh%edge2D_in > 0 .and. &
                   all(amesh%edge_tri(2, 1:amesh%edge2D_in) /= amesh%edge_tri(1, 1:amesh%edge2D_in)) .and. &
                   all(amesh%edge_tri(2, amesh%edge2D_in+1:) == amesh%edge_tri(1, amesh%edge2D_in+1:)), &
                   'patch: interior edges first, boundary edges last', 0)
        call check_mesh(amesh, apartit, 'patch', .true.)
        ! parallelogram patch: mesh_resolution of the uniform cells
        call generate_quad_patch(qmesh, qpartit, nx=11, ny=11, nl=4, &
                                 Lx=10000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP, skew=0.4_WP)
        call check_resolution(qmesh, 11, 11)
    end if

    call finish(partit, 'test_geom_mixed')

contains

    subroutine check_mesh(m, p, tag, cartesian)
        type(t_mesh),     intent(in) :: m
        type(t_partit),   intent(in) :: p
        character(len=*), intent(in) :: tag
        logical,          intent(in) :: cartesian
        integer :: nNodO, nNodL, nElemO, nElemF, nEdgeO, nEdgeL, ntri, nquad
        real(kind=WP) :: e1, e2, e3, e4, e5, e6(4)
        if (p%npes == 1) then
            nNodO = m%nod2D; nNodL = m%nod2D; nElemO = m%elem2D; nElemF = m%elem2D
            nEdgeO = m%edge2D; nEdgeL = m%edge2D
        else
            nNodO  = p%myDim_nod2D;  nNodL  = nNodO + p%eDim_nod2D
            nElemO = p%myDim_elem2D; nElemF = nElemO + p%eDim_elem2D + p%eXDim_elem2D
            nEdgeO = p%myDim_edge2D; nEdgeL = nEdgeO + p%eDim_edge2D
        end if
        ntri = count(m%elem2D_nnodes(1:nElemO) == 3); nquad = count(m%elem2D_nnodes(1:nElemO) == 4)

        e1 = area_balance(m, p, nNodO, nElemO)
        e2 = weight_sum_error(m, nElemO)
        e3 = node_gather_error(m, nNodL)
        e4 = linear_gradient_error(m, nElemO)
        e5 = closure_error(m, nElemO, nEdgeO)
        call max_all(e2, p); call max_all(e3, p); call max_all(e4, p); call max_all(e5, p)
        if (p%mype == 0) write(*,'(a,a,i0,a,i0,a,5es10.2)') tag, ': owned triangles=', ntri, &
            ' quads=', nquad, '  G1..G5 errors:', e1, e2, e3, e4, e5
        call check(e1 < 1.0e-12_WP, tag//': G1 sum node area == sum element area', p%mype)
        call check(e2 < 1.0e-13_WP, tag//': G2 kite weights sum to 1, positive, triangle slot 4 = 0', p%mype)
        call check(e3 < 1.0e-13_WP, tag//': G3 node gather weights reproduce the node area', p%mype)
        call check(e4 < 1.0e-9_WP,  tag//': G4 gradient_sca exact for linear fields', p%mype)
        call check(e5 < 1.0e-9_WP,  tag//': G5 cross-edge vectors close around each element', p%mype)
        call check(conventions_ok(m, nElemO, nElemF, nEdgeO, nEdgeL), &
                   tag//': G6 edge_tri / elem_edges / elem_neighbors conventions', p%mype)
        if (nquad > 0) then
            call quad_kite_error(m, nElemO, e6(1:3))
            call max_all(e6(1), p); call max_all(e6(2), p); call max_all(e6(3), p)
            if (p%mype == 0) write(*,'(a,a,3es10.2)') tag, ': G8 quads vs shoelace kites ' // &
                '(elem_w, centroid, max|w-1/4|):', e6(1:3)
            call check(e6(1) < 1.0e-12_WP .and. e6(2) < 1.0e-12_WP, &
                       tag//': G8 quad kite weights and centroid == shoelace formulas', p%mype)
            if (cartesian) call check(e6(3) > 1.0e-3_WP, tag//': G8 weights differ from 1/4', p%mype)
        end if
        if (ntri > 0 .and. .not. cartesian) then
            call previous_triangle_formulas(m, nNodO, nElemO, e6)
            call max_all(e6(1), p); call max_all(e6(2), p); call max_all(e6(3), p); call max_all(e6(4), p)
            if (p%mype == 0) write(*,'(a,a,4es10.2)') tag, ': G7 triangles vs previous formulas ' // &
                '(elem_cos, elem_area, gradient_sca, node area):', e6
            call check(e6(1) < 1.0e-13_WP .and. e6(2) < 1.0e-12_WP .and. e6(3) < 1.0e-10_WP &
                       .and. e6(4) < 1.0e-12_WP, tag//': G7 triangles equal the previous formulas', p%mype)
        end if
    end subroutine

    real(kind=WP) function area_balance(m, p, nNodO, nElemO)
        ! G1: relative difference of the global sums; an element owned by several ranks
        ! is counted once (by its global id, all copies are identical)
        type(t_mesh),   intent(in) :: m
        type(t_partit), intent(in) :: p
        integer,        intent(in) :: nNodO, nElemO
        real(kind=WP), allocatable :: ea(:), eag(:)
        real(kind=WP) :: an, ang
        integer :: e, ierr2
        allocate(ea(m%elem2D), eag(m%elem2D))
        ea = 0.0_WP
        do e = 1, nElemO
            if (p%npes == 1) then
                ea(e) = m%elem_area(e)
            else
                ea(p%myList_elem2D(e)) = m%elem_area(e)
            end if
        end do
        an = sum(real(m%area(1:nNodO), WP))
        if (p%npes > 1) then
            call MPI_Allreduce(ea, eag, m%elem2D, MPI_DOUBLE_PRECISION, MPI_MAX, p%MPI_COMM_FESOM, ierr2)
            call MPI_Allreduce(an, ang, 1, MPI_DOUBLE_PRECISION, MPI_SUM, p%MPI_COMM_FESOM, ierr2)
        else
            eag = ea; ang = an
        end if
        area_balance = abs(ang - sum(eag)) / sum(eag)
    end function

    real(kind=WP) function weight_sum_error(m, nElemO)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElemO
        integer :: e, nv
        weight_sum_error = 0.0_WP
        do e = 1, nElemO
            nv = m%elem2D_nnodes(e)
            weight_sum_error = max(weight_sum_error, abs(sum(real(m%elem_w(:, e), WP)) - 1.0_WP))
            if (any(m%elem_w(1:nv, e) <= 0.0_MP)) weight_sum_error = huge(1.0_WP)
            if (nv == 3 .and. m%elem_w(4, e) /= 0.0_MP) weight_sum_error = huge(1.0_WP)
        end do
    end function

    real(kind=WP) function node_gather_error(m, nNodL)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nNodL
        integer :: n, k
        real(kind=WP) :: s
        node_gather_error = 0.0_WP
        do n = 1, nNodL
            s = 0.0_WP
            do k = 1, m%nod_in_elem2D_num(n)
                s = s + m%nod_in_elem2D_w(k, n) * m%elem_area(m%nod_in_elem2D(k, n))
            end do
            node_gather_error = max(node_gather_error, abs(s - m%area(n)) / m%area(n))
        end do
    end function

    real(kind=WP) function linear_gradient_error(m, nElemO)
        ! f = a*X + b*Y with X = dlon*elem_cos*r_earth, Y = dlat*r_earth relative to
        ! vertex 1 (the element's own metric); error relative to |(a,b)|
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElemO
        real(kind=WP), parameter :: a = 1.7_WP, b = -0.6_WP
        integer :: e, k, n1, nk
        real(kind=WP) :: f(MAX_NV), dx, gx, gy
        linear_gradient_error = 0.0_WP
        do e = 1, nElemO
            n1 = m%elem2D_nodes(1, e)
            do k = 1, MAX_NV
                nk = m%elem2D_nodes(k, e)
                dx = m%coord_nod2D(1, nk) - m%coord_nod2D(1, n1); call trim_cyclic(dx)
                f(k) = a*dx*m%elem_cos(e)*r_earth + b*(m%coord_nod2D(2, nk) - m%coord_nod2D(2, n1))*r_earth
            end do
            gx = sum(m%gradient_sca(1:MAX_NV, e) * f)
            gy = sum(m%gradient_sca(MAX_NV+1:2*MAX_NV, e) * f)
            linear_gradient_error = max(linear_gradient_error, sqrt((gx - a)**2 + (gy - b)**2) / sqrt(a*a + b*b))
        end do
    end function

    real(kind=WP) function closure_error(m, nElemO, nEdgeO)
        ! G5 over owned elements whose edges are all owned (edge geometry is owned-only)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElemO, nEdgeO
        integer :: e, j, jm, nv, ed
        real(kind=WP) :: d(2, MAX_NV), wx, wy, sc
        closure_error = 0.0_WP
        elem_loop: do e = 1, nElemO
            nv = m%elem2D_nnodes(e)
            do j = 1, nv
                ed = m%elem_edges(j, e)
                if (ed > nEdgeO) cycle elem_loop
                if (m%edge_tri(1, ed) == e) then
                    d(:, j) = m%edge_cross_dxdy(1:2, ed)
                else
                    d(:, j) = m%edge_cross_dxdy(3:4, ed)
                end if
            end do
            sc = maxval(abs(d(:, 1:nv)))
            do j = 1, nv
                jm = mod(j+nv-2, nv) + 1
                wx = d(1, j) - d(1, jm); wy = d(2, j) - d(2, jm)
                closure_error = max(closure_error, &
                    abs(m%elem_area(e)*m%gradient_sca(j, e)        - wy) / sc, &
                    abs(m%elem_area(e)*m%gradient_sca(j+MAX_NV, e) + wx) / sc)
            end do
        end do elem_loop
    end function

    logical function conventions_ok(m, nElemO, nElemF, nEdgeO, nEdgeL)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElemO, nElemF, nEdgeO, nEdgeL
        integer :: e, k, nv, ed, v1, v2, other
        conventions_ok = all(m%edge_tri(:, 1:nEdgeL) >= 1 .and. m%edge_tri(:, 1:nEdgeL) <= nElemF)
        do ed = 1, nEdgeO
            if (m%edge_tri(2, ed) == m%edge_tri(1, ed)) then
                if (any(m%edge_cross_dxdy(3:4, ed) /= 0.0_MP)) conventions_ok = .false.
            end if
        end do
        do e = 1, nElemO
            nv = m%elem2D_nnodes(e)
            if (nv == 3) then
                if (m%elem_edges(4, e) /= m%elem_edges(1, e) .or. &
                    m%elem_neighbors(4, e) /= m%elem_neighbors(1, e)) conventions_ok = .false.
            end if
            do k = 1, nv
                ed = m%elem_edges(k, e)
                v1 = m%elem2D_nodes(k, e); v2 = m%elem2D_nodes(mod(k, nv) + 1, e)
                if (.not. ((m%edges(1, ed) == v1 .and. m%edges(2, ed) == v2) .or. &
                           (m%edges(1, ed) == v2 .and. m%edges(2, ed) == v1))) conventions_ok = .false.
                if (m%edge_tri(1, ed) == e) then
                    other = m%edge_tri(2, ed)
                else if (m%edge_tri(2, ed) == e) then
                    other = m%edge_tri(1, ed)
                else
                    other = -1
                end if
                if (m%elem_neighbors(k, e) /= other) conventions_ok = .false.
            end do
        end do
    end function

    subroutine previous_triangle_formulas(m, nNodO, nElemO, err)
        ! G7: FESOM2 oce_mesh.F90 elem_center (vertex mean after the min-longitude wrap),
        ! mesh_areas (cos of the mean latitude) and the P1 gradient, as the previous
        ! FESOM3 code had them; err = (elem_cos, elem_area, gradient_sca, node area)
        type(t_mesh),  intent(in)  :: m
        integer,       intent(in)  :: nNodO, nElemO
        real(kind=WP), intent(out) :: err(4)
        integer :: e, j, n, n1, n2, n3
        real(kind=WP) :: ax(3), amin, cl, cy, ay, a1, a2, b1, b2, aold, g(6), gmax, s
        real(kind=WP) :: dX31, dX21, dY31, dY21, dfac
        real(kind=WP), allocatable :: aold_e(:)
        logical :: alltri
        allocate(aold_e(nElemO))
        err = 0.0_WP
        cl = get_cyclic_length()
        do e = 1, nElemO
            if (m%elem2D_nnodes(e) /= 3) cycle
            n1 = m%elem2D_nodes(1, e); n2 = m%elem2D_nodes(2, e); n3 = m%elem2D_nodes(3, e)
            ax = m%coord_nod2D(1, [n1, n2, n3])
            amin = minval(ax)
            do j = 1, 3
                if (ax(j) - amin >= cl/2.0_WP) ax(j) = ax(j) - cl
            end do
            cy = sum(m%coord_nod2D(2, [n1, n2, n3])) / 3.0_WP
            err(1) = max(err(1), abs(m%elem_cos(e) - cos(cy)))
            ay = cos(cy)
            a1 = m%coord_nod2D(1, n2) - m%coord_nod2D(1, n1); a2 = m%coord_nod2D(2, n2) - m%coord_nod2D(2, n1)
            b1 = m%coord_nod2D(1, n3) - m%coord_nod2D(1, n1); b2 = m%coord_nod2D(2, n3) - m%coord_nod2D(2, n1)
            call trim_cyclic(a1); call trim_cyclic(b1)
            aold = 0.5_WP * abs(a1*ay*b2 - b1*ay*a2) * r_earth * r_earth
            aold_e(e) = aold
            err(2) = max(err(2), abs(m%elem_area(e) - aold) / aold)
            dX31 = b1 * m%elem_cos(e); dX21 = a1 * m%elem_cos(e); dY31 = b2; dY21 = a2
            dfac = -0.5_WP * r_earth / aold
            g = [(-dY31 + dY21), dY31, -dY21, (dX31 - dX21), -dX31, dX21] * dfac
            gmax = maxval(abs(g))
            err(3) = max(err(3), maxval(abs(m%gradient_sca([1, 2, 3, 5, 6, 7], e) - g)) / gmax, &
                                 maxval(abs(m%gradient_sca([4, 8], e))) / gmax)
        end do
        do n = 1, nNodO
            alltri = .true.
            s = 0.0_WP
            do j = 1, m%nod_in_elem2D_num(n)
                e = m%nod_in_elem2D(j, n)
                if (m%elem2D_nnodes(e) /= 3) alltri = .false.
                if (alltri) s = s + aold_e(e) / 3.0_WP
            end do
            if (alltri) err(4) = max(err(4), abs(m%area(n) - s) / s)
        end do
    end subroutine

    subroutine check_output_coords(m, p)
        ! G10 over the owned elements
        type(t_mesh),   intent(in) :: m
        type(t_partit), intent(in) :: p
        integer :: nNodO, nNodL, nEdgeO, nElemO, e, j, jp, nv
        integer,  allocatable :: nlev(:)
        real(kind=WP), allocatable :: lon(:), lat(:), rlon(:), rlat(:)
        real(kind=WP) :: x(4), y(4), a, cx, cy, cr, glon, glat, dlon, err
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        allocate(lon(nElemO), lat(nElemO), rlon(nElemO), rlat(nElemO), nlev(nElemO))
        call io_coords_compute(DECOMP_ELEM, m, nElemO, lon, lat, rlon, rlat, nlev)
        err = 0.0_WP
        do e = 1, nElemO
            nv = m%elem2D_nnodes(e)
            do j = 1, nv
                x(j) = m%coord_nod2D(1, m%elem2D_nodes(j, e)) - m%coord_nod2D(1, m%elem2D_nodes(1, e))
                call trim_cyclic(x(j))
                y(j) = m%coord_nod2D(2, m%elem2D_nodes(j, e))
            end do
            a = 0.0_WP; cx = 0.0_WP; cy = 0.0_WP
            do j = 1, nv
                jp = mod(j, nv) + 1
                cr = x(j)*y(jp) - x(jp)*y(j)
                a = a + 0.5_WP*cr
                cx = cx + (x(j) + x(jp))*cr; cy = cy + (y(j) + y(jp))*cr
            end do
            cx = cx/(6.0_WP*a) + m%coord_nod2D(1, m%elem2D_nodes(1, e)); cy = cy/(6.0_WP*a)
            call r2g(glon, glat, cx, cy)
            dlon = modulo(lon(e) - glon/rad + 180.0_WP, 360.0_WP) - 180.0_WP
            err = max(err, abs(dlon), abs(lat(e) - glat/rad))
        end do
        call max_all(err, p)
        if (p%mype == 0) write(*,'(a,es10.2)') 'mesh: G10 output element lon/lat vs r2g(shoelace centroid), deg:', err
        call check(err < 1.0e-9_WP, 'mesh: G10 output element coordinates == the polygon centroid', p%mype)
    end subroutine check_output_coords

    subroutine check_resolution(m, nx, ny)
        ! G9: node (i,j) of the nx x ny parallelogram patch; ring = distance to the rim
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nx, ny
        integer :: n, i, j, ring, nin
        real(kind=WP) :: r0, err, dev
        n = (ny/2)*nx + nx/2 + 1                      ! the centre node
        r0 = 2.0_WP*sqrt(real(m%areasvol(n), WP)/acos(-1.0_WP))
        err = 0.0_WP; dev = 0.0_WP; nin = 0
        do n = 1, m%nod2D
            i = mod(n-1, nx) + 1; j = (n-1)/nx + 1
            ring = min(i-1, nx-i, j-1, ny-j)
            if (ring >= 4) then
                err = max(err, abs(real(m%mesh_resolution(n), WP) - r0)/r0); nin = nin + 1
            else if (ring == 1) then
                dev = max(dev, abs(real(m%mesh_resolution(n), WP) - r0)/r0)
            end if
        end do
        write(*,'(a,i0,a,2es10.2)') 'parallelogram patch: G9 mesh_resolution, ', nin, &
            ' deep nodes rel. err, rim-near deviation:', err, dev
        call check(nin > 0 .and. err < 1.0e-13_WP, 'parallelogram patch: G9 mesh_resolution of uniform cells', 0)
        call check(dev > 1.0e-3_WP, 'parallelogram patch: G9 rim-near nodes differ (non-vacuity)', 0)
    end subroutine

    subroutine quad_kite_error(m, nElemO, err)
        ! G8 for owned quads, in lon/lat relative to vertex 1 (trim_cyclic): polygon area
        ! A and centroid c by the shoelace formulas; kite j = (p_j, mid(j,j+1), c,
        ! mid(j-1,j)) by the shoelace formula; err = (max|elem_w - kite/A|, max distance
        ! elem_center - c relative to sqrt(A), max|elem_w - 1/4|)
        type(t_mesh),  intent(in)  :: m
        integer,       intent(in)  :: nElemO
        real(kind=WP), intent(out) :: err(3)
        integer :: e, j, jp, jm
        real(kind=WP) :: x(4), y(4), cr, a, cx, cy, ex, ey, kx(4), ky(4), kite
        err = 0.0_WP
        do e = 1, nElemO
            if (m%elem2D_nnodes(e) /= 4) cycle
            do j = 1, 4
                x(j) = m%coord_nod2D(1, m%elem2D_nodes(j, e)) - m%coord_nod2D(1, m%elem2D_nodes(1, e))
                call trim_cyclic(x(j))
                y(j) = m%coord_nod2D(2, m%elem2D_nodes(j, e)) - m%coord_nod2D(2, m%elem2D_nodes(1, e))
            end do
            a = 0.0_WP; cx = 0.0_WP; cy = 0.0_WP
            do j = 1, 4
                jp = mod(j, 4) + 1
                cr = x(j)*y(jp) - x(jp)*y(j)
                a = a + 0.5_WP*cr
                cx = cx + (x(j) + x(jp))*cr; cy = cy + (y(j) + y(jp))*cr
            end do
            cx = cx / (6.0_WP*a); cy = cy / (6.0_WP*a)
            do j = 1, 4
                jp = mod(j, 4) + 1; jm = mod(j+2, 4) + 1
                kx = [x(j), 0.5_WP*(x(j)+x(jp)), cx, 0.5_WP*(x(j)+x(jm))]
                ky = [y(j), 0.5_WP*(y(j)+y(jp)), cy, 0.5_WP*(y(j)+y(jm))]
                kite = 0.5_WP*(sum(kx*cshift(ky, 1)) - sum(cshift(kx, 1)*ky))
                err(1) = max(err(1), abs(m%elem_w(j, e) - kite/a))
                err(3) = max(err(3), abs(m%elem_w(j, e) - 0.25_WP))
            end do
            call elem_center(m, e, ex, ey)
            ex = ex - m%coord_nod2D(1, m%elem2D_nodes(1, e)); call trim_cyclic(ex)
            ey = ey - m%coord_nod2D(2, m%elem2D_nodes(1, e))
            err(2) = max(err(2), sqrt((ex - cx)**2 + (ey - cy)**2) / sqrt(abs(a)))
        end do
    end subroutine

    logical function patch_convex(m)
        ! every corner of every element turns clockwise: (next - v) x (previous - v) < 0
        type(t_mesh), intent(in) :: m
        integer :: e, j, nv, p0, p1, p2
        real(kind=WP) :: b1, b2, c1, c2
        patch_convex = .true.
        do e = 1, m%elem2D
            nv = m%elem2D_nnodes(e)
            do j = 1, nv
                p0 = m%elem2D_nodes(j, e)
                p1 = m%elem2D_nodes(mod(j, nv) + 1, e)
                p2 = m%elem2D_nodes(mod(j+nv-2, nv) + 1, e)
                b1 = m%coord_nod2D(1, p1) - m%coord_nod2D(1, p0); b2 = m%coord_nod2D(2, p1) - m%coord_nod2D(2, p0)
                c1 = m%coord_nod2D(1, p2) - m%coord_nod2D(1, p0); c2 = m%coord_nod2D(2, p2) - m%coord_nod2D(2, p0)
                if (b1*c2 - b2*c1 >= 0.0_WP) patch_convex = .false.
            end do
        end do
    end function

end program test_geom_mixed
