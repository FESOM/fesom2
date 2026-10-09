program test_updown_mixed
    ! MUSCL up/downwind elements on mixed triangle/quad meshes (plan 2026-10-08-hybrid-meshes,
    ! Task 11): find_up_downwind_elements (oce_muscl_adv) searches, at each edge node, the
    ! element whose vertex angle (between the sides to the next and the previous vertex)
    ! contains the edge direction. Checked on the mesh read from FESOM3_MESH_DIR (default the
    ! NEMO-mixed cut-out, unrotated; FESOM3_TEST_ROTATE=1 for pi/CORE2) at 1 or N ranks and,
    ! at 1 rank, on the analytic mixed patch and the skewed parallelogram patch:
    !   U1 every edge node that is INTERIOR (its vertex angles sum to 2*pi) has its
    !      up/downwind element (side 1 = upwind at node 1 for -edge, side 2 = downwind at
    !      node 2 for +edge)
    !   U2 every element found lies around its anchor node and the search direction lies in
    !      its vertex angle there -- checked independently with cross products (the search
    !      uses atan2), over the full element halo at N ranks
    !   L1 a tracer linear in (x, y) = r_earth*(lon, lat) (cartesian geometry): the elemental
    !      gradient is exact, and on every edge with both up- and downwind elements the
    !      MUSCL face values equal the exact edge-midpoint value, for both velocity signs
    !      and any 3rd/4th-order mix (edge flux = -vflux*(T1+T2)/2)
    !   L2 uniform flow advects the linear tracer exactly: the node tendency of the MUSCL
    !      fluxes = -u.grad(T) at interior nodes whose edges all have both elements --
    !      roundoff where the dual cell is made of triangles and parallelograms only (the
    !      median-dual face values are midpoint values, exact on those shapes); first
    !      order on general quads (reported, bounded)
    ! "MUSCL equals the reference rule" is test_muscl_onthefly (part R) on the same meshes.
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all, sum_all
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth
    use mod_mesh,          only: t_mesh, MAX_NV
    use mod_partit,        only: t_partit
    use mod_tracer,        only: t_tracer_work
    use mod_partitioning,  only: par_init
    use mod_mesh_analytic, only: generate_mixed_patch, generate_quad_patch
    use mod_mesh_rotate,   only: get_cyclic_length
    use mod_part_bounds,   only: owned_bounds, is_multirank
    use mod_halo,          only: exchange_elem_full
    use oce_tracer_grad,   only: tracer_gradient_elements
    use oce_muscl_adv,     only: muscl_adv_init
    use oce_adv_tra_hor,   only: adv_tra_hor_upw1, adv_tra_hor_muscl
    implicit none

    type(t_partit) :: partit, apartit, qpartit
    type(t_mesh)   :: mesh, amesh, qmesh
    real(kind=WP), parameter :: pi2 = 2.0_WP*3.141592653589793_WP

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.true., partition_np1=.false.)
    call check_mesh(mesh, partit, 'mesh')

    if (partit%npes == 1) then
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=4, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_mesh(amesh, apartit, 'mixed patch')
        call generate_quad_patch(qmesh, qpartit, nx=11, ny=9, nl=4, &
                                 Lx=10000.0_WP, Ly=8000.0_WP, max_depth=100.0_WP, skew=0.4_WP)
        call check_mesh(qmesh, qpartit, 'parallelogram patch')
    end if

    call finish(partit, 'test_updown_mixed')

contains

    pure function wrap(v, cl) result(w)
        real(kind=WP), intent(in) :: v(2), cl
        real(kind=WP) :: w(2)
        w = v
        if (w(1) >  cl/2.0_WP) w(1) = w(1) - cl
        if (w(1) < -cl/2.0_WP) w(1) = w(1) + cl
    end function wrap

    pure real(kind=WP) function cross(a, b)
        real(kind=WP), intent(in) :: a(2), b(2)
        cross = a(1)*b(2) - a(2)*b(1)
    end function cross

    subroutine check_mesh(m, p, tag)
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(inout) :: p
        character(len=*), intent(in)    :: tag
        type(t_tracer_work) :: twork
        integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, e, k, j, kk, nz, ed, isg
        integer :: anchor, el, n_int, n_miss, n_bad, n_found, n_lin, n_tri, n_par, n_gen, ntri, nquad, nk
        integer :: tri1, tri2, u1, l1, u2, l2, nodes(2)
        logical :: lmr, inside
        real(kind=WP) :: cl, x(2), q(2), pv(2), s, ang, tol, sc, err_g, err_f, err_tt, err_tp, err_tg, expd
        real(kind=WP), allocatable :: xe(:,:,:), tmp(:), angsum(:)
        real(kind=WP), allocatable :: ttf(:,:), one(:,:), tr_xy(:,:,:), vel(:,:,:)
        real(kind=WP), allocatable :: f_lo(:,:), f_ho(:,:), tend(:), xn(:), yn(:)
        integer,       allocatable :: ide(:,:), tmp_i(:), nvf(:), gid(:)
        logical,       allocatable :: full(:), seam(:)
        integer,       allocatable :: shape_e(:)
        real(kind=WP), parameter :: A = 2.0e-4_WP, B = -3.0e-4_WP     ! tracer gradient [1/m]
        real(kind=WP), parameter :: U = 0.3_WP, V = -0.2_WP           ! uniform flow [m/s]

        lmr = is_multirank(p)
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        nElemF = m%elem2D
        if (lmr) nElemF = p%myDim_elem2D + p%eDim_elem2D + p%eXDim_elem2D
        nl = m%nl
        cl = get_cyclic_length()
        ntri  = count(m%elem2D_nnodes(1:nElemO) == 3); call sum_all(ntri, p)
        nquad = count(m%elem2D_nnodes(1:nElemO) == 4); call sum_all(nquad, p)
        if (p%mype == 0) write(*,'(3a,i0,a,i0)') '  [', tag, '] owned triangles=', ntri, ' quads=', nquad

        if (allocated(twork%edge_up_dn_tri)) deallocate(twork%edge_up_dn_tri)
        call muscl_adv_init(twork, m, p)

        ! vertex coordinates and GLOBAL vertex ids over the full element halo
        allocate(xe(2, MAX_NV, nElemF), ide(MAX_NV, nElemF), nvf(nElemF))
        if (.not. lmr) then
            do e = 1, nElemF
                do j = 1, MAX_NV
                    xe(:, j, e) = real(m%coord_nod2D(:, m%elem2D_nodes(j, e)), WP)
                    ide(j, e)   = m%elem2D_nodes(j, e)
                end do
            end do
        else
            allocate(tmp(nElemF), tmp_i(nElemF))
            do j = 1, MAX_NV
                do kk = 1, 2
                    tmp = 0.0_WP
                    do e = 1, nElemO
                        tmp(e) = real(m%coord_nod2D(kk, m%elem2D_nodes(j, e)), WP)
                    end do
                    call exchange_elem_full(tmp, p)
                    xe(kk, j, :) = tmp
                end do
                tmp_i = 0
                do e = 1, nElemO
                    tmp_i(e) = p%myList_nod2D(m%elem2D_nodes(j, e))
                end do
                call exchange_elem_full(tmp_i, p)
                ide(j, :) = tmp_i
            end do
            deallocate(tmp, tmp_i)
        end if
        do e = 1, nElemF
            nvf(e) = MAX_NV
            if (ide(MAX_NV, e) == ide(1, e)) nvf(e) = 3
        end do
        allocate(gid(nNodL))
        do n = 1, nNodL
            gid(n) = n
            if (lmr) gid(n) = p%myList_nod2D(n)
        end do

        ! vertex angle sums of the local nodes: interior <=> 2*pi
        allocate(angsum(nNodL))
        do n = 1, nNodL
            angsum(n) = 0.0_WP
            do k = 1, m%nod_in_elem2D_num(n)
                call sides(m%nod_in_elem2D(k, n), gid(n), ide, xe, nvf, cl, q, pv)
                angsum(n) = angsum(n) + acos(max(-1.0_WP, min(1.0_WP, &
                            sum(q*pv)/sqrt(sum(q*q)*sum(pv*pv)))))
            end do
        end do

        !--- U1/U2 ---
        n_int = 0; n_miss = 0; n_bad = 0; n_found = 0
        tol = 1.0e-9_WP
        do ed = 1, nEdgeO
            nodes = m%edges(:, ed)
            do k = 1, 2
                ! k = 1: direction -edge at node 1; k = 2: +edge at node 2
                anchor = nodes(k)
                x = wrap(real(m%coord_nod2D(:, nodes(2)) - m%coord_nod2D(:, nodes(1)), WP), cl)
                if (k == 1) x = -x
                el = twork%edge_up_dn_tri(k, ed)
                if (abs(angsum(anchor) - pi2) < 1.0e-8_WP) n_int = n_int + 1
                if (el == 0) then
                    if (abs(angsum(anchor) - pi2) < 1.0e-8_WP) n_miss = n_miss + 1
                    cycle
                end if
                n_found = n_found + 1
                inside = any(m%nod_in_elem2D(1:m%nod_in_elem2D_num(anchor), anchor) == el)
                if (inside) then
                    call sides(el, gid(anchor), ide, xe, nvf, cl, q, pv)
                    s = sign(1.0_WP, cross(q, pv))
                    sc = sqrt(sum(x*x))
                    inside = s*cross(q, x)/(sqrt(sum(q*q))*sc)   >= -tol .and. &
                             s*cross(x, pv)/(sqrt(sum(pv*pv))*sc) >= -tol
                end if
                if (.not. inside) n_bad = n_bad + 1
            end do
        end do
        call sum_all(n_int, p); call sum_all(n_miss, p); call sum_all(n_bad, p)
        call sum_all(n_found, p)
        if (p%mype == 0) write(*,'(3a,i0,a,i0,a,i0,a,i0)') '  [', tag, '] edge sides at interior nodes=', &
            n_int, ' missing=', n_miss, '  found=', n_found, ' outside the vertex angle=', n_bad
        call check(n_int > 0,   tag//' U1: interior edge nodes exist', p%mype)
        call check(n_miss == 0, tag//' U1: up/downwind element found at every interior edge node', p%mype)
        call check(n_bad == 0,  tag//' U2: found element around the anchor, direction in its vertex angle', p%mype)

        !--- L1/L2: linear tracer, uniform flow ---
        ! node coordinates in metres relative to this rank's first node; elements spanning
        ! the cyclic seam (global meshes) are excluded
        allocate(xn(nNodL), yn(nNodL))
        do n = 1, nNodL
            xn(n) = real(m%coord_nod2D(1, n), WP)*r_earth   ! NOT wrapped: the linear field
            yn(n) = real(m%coord_nod2D(2, n), WP)*r_earth   ! jumps only at the cyclic seam
        end do
        allocate(seam(nElemF), shape_e(nElemF))
        do e = 1, nElemF
            seam(e) = maxval(xe(1, :, e)) - minval(xe(1, :, e)) > cl/2.0_WP
            ! triangle (1), parallelogram v1 + v3 = v2 + v4 (2), general quad (3)
            x = xe(:, 1, e) + xe(:, 3, e) - xe(:, 2, e) - xe(:, 4, e)
            shape_e(e) = 3                                      ! general quad
            if (sqrt(sum(x*x)) < 1.0e-9_WP*sqrt(sum((xe(:, 3, e) - xe(:, 1, e))**2))) shape_e(e) = 2
            if (nvf(e) == 3) shape_e(e) = 1                    ! triangle
        end do
        allocate(ttf(nl-1, nNodL), one(nl-1, nNodL), tr_xy(2, nl-1, nElemF), vel(2, nl-1, nElemF))
        allocate(f_lo(nl-1, nEdgeO), f_ho(nl-1, nEdgeO), tend(nNodO), full(nNodL))
        if (allocated(m%helem)) deallocate(m%helem)
        allocate(m%helem(nl-1, nElemF))
        m%helem = 1.0_MP
        do n = 1, nNodL
            ttf(:, n) = 1.0_WP + A*xn(n) + B*yn(n)
        end do
        one = 1.0_WP
        tr_xy = 0.0_WP
        call tracer_gradient_elements(ttf, tr_xy, m, p)
        if (lmr) call exchange_elem_full(tr_xy, p)
        err_g = 0.0_WP
        do e = 1, nElemO
            if (seam(e)) cycle
            do nz = m%ulevels(e), m%nlevels(e)-1
                err_g = max(err_g, abs(tr_xy(1, nz, e) - A), abs(tr_xy(2, nz, e) - B))
            end do
        end do
        call max_all(err_g, p); err_g = err_g/abs(B)

        err_f = 0.0_WP; err_tt = 0.0_WP; err_tp = 0.0_WP; err_tg = 0.0_WP
        n_lin = 0; n_tri = 0; n_par = 0; n_gen = 0
        do isg = 1, 2
            vel(1, :, :) = U; vel(2, :, :) = V
            if (isg == 2) vel = -vel
            ! -vflux from the upwind flux of T = 1; then the MUSCL flux of the linear T
            call adv_tra_hor_upw1(vel, one, m, f_lo, o_init_zero=.true., partit=p)
            call adv_tra_hor_muscl(vel, ttf, m, 0.3_WP, f_ho, tr_xy, twork%edge_up_dn_tri, &
                                   o_init_zero=.true., partit=p)
            full(:) = .true.
            tend(:) = 0.0_WP
            do ed = 1, nEdgeO
                nodes = m%edges(:, ed)
                tri1 = twork%edge_up_dn_tri(1, ed); tri2 = twork%edge_up_dn_tri(2, ed)
                if (tri1 == 0 .or. tri2 == 0) then
                    full(nodes) = .false.
                    cycle
                end if
                if (seam(tri1) .or. seam(tri2) .or. seam(m%edge_tri(1, ed)) .or. seam(m%edge_tri(2, ed))) then
                    full(nodes) = .false.
                    cycle
                end if
                u1 = m%ulevels(tri1); l1 = m%nlevels(tri1) - 1
                u2 = m%ulevels(tri2); l2 = m%nlevels(tri2) - 1
                do nz = max(u1, u2), min(l1, l2)
                    expd = f_lo(nz, ed)*0.5_WP*(ttf(nz, nodes(1)) + ttf(nz, nodes(2)))
                    err_f = max(err_f, abs(f_ho(nz, ed) - expd)/(abs(f_lo(nz, ed)) + tiny(1.0_WP)))
                    if (isg == 1) n_lin = n_lin + 1
                end do
                ! surface-level node tendency (helem = 1: the volume is the node area)
                if (nodes(1) <= nNodO) tend(nodes(1)) = tend(nodes(1)) + f_ho(1, ed)
                if (nodes(2) <= nNodO) tend(nodes(2)) = tend(nodes(2)) - f_ho(1, ed)
            end do
            expd = -(vel(1, 1, 1)*A + vel(2, 1, 1)*B)
            do n = 1, nNodO
                if (.not. full(n) .or. abs(angsum(n) - pi2) >= 1.0e-8_WP) cycle
                if (m%ulevels_nod2D(n) > 1) cycle
                ang = abs(tend(n)/real(m%area(n), WP) - expd)/abs(expd)
                nk = m%nod_in_elem2D_num(n)
                if (all(shape_e(m%nod_in_elem2D(1:nk, n)) == 1)) then
                    err_tt = max(err_tt, ang)
                    if (isg == 1) n_tri = n_tri + 1
                else if (all(shape_e(m%nod_in_elem2D(1:nk, n)) == 2)) then
                    err_tp = max(err_tp, ang)
                    if (isg == 1) n_par = n_par + 1
                else
                    err_tg = max(err_tg, ang)
                    if (isg == 1) n_gen = n_gen + 1
                end if
            end do
        end do
        call max_all(err_f, p); call max_all(err_tt, p); call max_all(err_tp, p); call max_all(err_tg, p)
        call sum_all(n_lin, p); call sum_all(n_tri, p); call sum_all(n_par, p); call sum_all(n_gen, p)
        if (p%mype == 0) then
            write(*,'(3a,es10.3)') '  [', tag, '] L1 linear tracer: gradient rel err ', err_g
            write(*,'(3a,i0,a,es10.3)') '  [', tag, '] L1 MUSCL face value == midpoint value: edge-levels ', &
                n_lin, '  max rel err ', err_f
            write(*,'(3a,2(i0,a,es10.3,a),i0,a,es10.3)') '  [', tag, '] L2 node tendency vs -u.grad(T): all-triangle nodes ', &
                n_tri, ' max rel err ', err_tt, ';  all-parallelogram nodes ', n_par, ' max rel err ', err_tp, &
                ';  mixed/general-quad nodes ', n_gen, ' max rel err ', err_tg
        end if
        call check(err_g < 1.0e-9_WP, tag//' L1: elemental gradient of a linear tracer exact', p%mype)
        call check(n_lin > 0,         tag//' L1: edges with both up- and downwind elements exist', p%mype)
        call check(err_f < 1.0e-9_WP, tag//' L1: MUSCL face values of a linear tracer exact', p%mype)
        call check(err_tt < 1.0e-9_WP, tag//' L2: exact advection at all-triangle nodes', p%mype)
        call check(err_tp < 1.0e-9_WP, tag//' L2: exact advection at all-parallelogram nodes', p%mype)
        call check(err_tg < 0.5_WP,    tag//' L2: bounded error at triangle-quad / general-quad nodes', p%mype)
        ! non-vacuity: every case has L2 nodes; all but the jittered mixed patch (general
        ! quads only, by design) have nodes with the exact expectation
        call check(n_tri + n_par + n_gen > 0, tag//' L2: nodes checked (non-vacuity)', p%mype)
        if (tag /= 'mixed patch') &
            call check(n_tri + n_par > 0, tag//' L2: all-triangle / all-parallelogram nodes checked', p%mype)

        deallocate(xe, ide, nvf, gid, angsum, xn, yn, seam, shape_e, ttf, one, tr_xy, vel, f_lo, f_ho, tend, full)
        deallocate(twork%edge_up_dn_tri)

    end subroutine check_mesh

    subroutine sides(eel, gidn, ide, xe, nvf, cl, qq, pp)
        ! sides of element eel from its vertex with GLOBAL id gidn to the next (qq) and the
        ! previous (pp) unique vertex
        integer,       intent(in)  :: eel, gidn, ide(:,:), nvf(:)
        real(kind=WP), intent(in)  :: xe(:,:,:), cl
        real(kind=WP), intent(out) :: qq(2), pp(2)
        integer :: jj, ja, nv
        nv = nvf(eel)
        ja = 0
        do jj = 1, nv
            if (ide(jj, eel) == gidn) ja = jj
        end do
        if (ja == 0) then
            qq = 0.0_WP; pp = 0.0_WP
            return
        end if
        qq = wrap(xe(:, mod(ja, nv) + 1, eel)          - xe(:, ja, eel), cl)
        pp = wrap(xe(:, mod(ja + nv - 2, nv) + 1, eel) - xe(:, ja, eel), cl)
    end subroutine sides

end program test_updown_mixed
