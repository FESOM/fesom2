program test_ssh_stiff
    ! SSH stiffness matrix on mixed meshes (plan 2026-10-08-hybrid-meshes, Task 9):
    ! init_stiff_mat_ale / update_stiff_mat_ale on a mesh read from FESOM3_MESH_DIR
    ! (default the NEMO quad cut-out; FESOM3_TEST_ROTATE=1 reads it rotated as pi/CORE2)
    ! at 1 or N ranks and, at 1 rank, on two analytic patches (the jittered mixed patch
    ! and a parallelogram all-quad patch). With A = M/dt + K (M = areasvol on the
    ! diagonal, K the g*dt*alpha*theta*H*div-grad stiffness part):
    !   S1 sparsity: row n = n first (Jacobi preconditioner), then exactly the other
    !      vertices of the elements around n (quad diagonals included), no duplicates
    !   S2 row sums = areasvol/dt (K annihilates constants)
    !   S3 every quad diagonal entry is present and nonzero
    !   S4 A symmetric (all rows gathered on rank 0 in global numbering)
    !   S5 triangle-only meshes: CSR and values equal the previous assembly (edge-only
    !      sparsity with 12 slots, 3-slot gradient) to roundoff
    !   S6 update_stiff_mat_ale with dhe = -H/2 adds exactly K/2
    !   S7 the CG (solve_ssh_ale) converges to x for b = A x, x smooth
    !   S8 checkerboard: f = (+1,-1,+1,-1) on the vertices of every quad is in the null
    !      space of K -- on the parallelogram patch (the documented property) and, when
    !      the mesh is all-quad and 2-colourable along its edges, on the mesh read (np 1;
    !      FESOM3_EXPECT_CHECKERBOARD=1: a mesh that is not is a failure, not a skip);
    !      a stripe pattern is NOT (non-vacuity). The Green-Gauss quad gradient vanishes
    !      for it on any quad, so only the mass part controls this mode.
    !   S9 values: on the parallelogram patch K f == -g*dt*alpha*theta*H*lap(f)*area(n) at
    !      the interior nodes for quadratic f (the scale of the quad contributions; S2/S4/
    !      S6/S8 are blind to a per-element factor)
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all, sum_all, set_depths
    use mod_precision,     only: WP, MP
    use mod_mesh,          only: t_mesh
    use mod_partit,        only: t_partit
    use mod_dyn,           only: t_dyn
    use mod_partitioning,  only: par_init
    use mod_mesh_analytic, only: generate_mixed_patch, generate_quad_patch
    use mod_part_bounds,   only: owned_bounds
    use mod_constants,     only: g, r_earth
    use mod_param_phys,    only: alpha, theta
    use oce_ssh_rhs,       only: init_stiff_mat_ale, update_stiff_mat_ale
    use oce_ssh_solve,     only: solve_ssh_ale
    implicit none

    real(kind=WP), parameter :: dt = 1200.0_WP
    type(t_partit) :: partit, ppartit, qpartit
    type(t_mesh)   :: mesh, pmesh, qmesh
    type(t_dyn)    :: dyn

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025q_cut', cartesian=.false., partition_np1=.true.)
    call set_depths(mesh, partit)
    call init_stiff_mat_ale(mesh, dt, partit)
    call check_matrix(mesh, partit, 'mesh')
    call check_previous_triangles(mesh, partit, 'mesh')
    call check_update(mesh, partit, 'mesh')
    call check_cg(mesh, partit, 'mesh')
    if (partit%npes == 1) call check_checkerboard_mesh(mesh, partit, 'mesh')

    if (partit%npes == 1) then
        ! analytic patches (1-rank identity partitions of their own)
        call generate_mixed_patch(pmesh, ppartit, nx=13, ny=11, nl=4, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call set_depths(pmesh, ppartit)
        call init_stiff_mat_ale(pmesh, dt, ppartit)
        call check_matrix(pmesh, ppartit, 'mixed patch')
        call check_update(pmesh, ppartit, 'mixed patch')

        call generate_quad_patch(qmesh, qpartit, nx=11, ny=9, nl=4, &
                                 Lx=10000.0_WP, Ly=8000.0_WP, max_depth=100.0_WP, skew=0.4_WP)
        call set_depths(qmesh, qpartit)
        call init_stiff_mat_ale(qmesh, dt, qpartit)
        call check_matrix(qmesh, qpartit, 'parallelogram patch')
        call check_checkerboard_patch(qmesh, 11, 'parallelogram patch')
        call check_quadratic_patch(qmesh, 11, 9, 'parallelogram patch')
    end if

    call finish(partit, 'test_ssh_stiff')

contains

    integer function gid(p, n)
        type(t_partit), intent(in) :: p
        integer,        intent(in) :: n
        if (p%npes == 1) then
            gid = n
        else
            gid = p%myList_nod2D(n)
        end if
    end function

    real(kind=WP) function mass(m, row)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: row
        mass = 0.0_WP
        if (m%ulevels_nod2D(row) == 1) mass = real(m%areasvol(row), WP)/dt
    end function

    subroutine check_matrix(m, p, tag)
        ! S1 sparsity, S2 row sums, S3 quad diagonals, S4 symmetry
        type(t_mesh),     intent(in) :: m
        type(t_partit),   intent(in) :: p
        character(len=*), intent(in) :: tag
        integer :: nNodO, nNodL, nEdgeO, nElemO, row, k, j, el, nv, c, nini, nend, nexp, nquad
        integer, allocatable :: mark(:)
        logical :: ok_sp
        real(kind=WP) :: s, sa, e_sum, e_sym, offmax, dmin, a
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        allocate(mark(nNodL)); mark = 0
        ok_sp = .true.; e_sum = 0.0_WP; nquad = 0; dmin = huge(1.0_WP)
        do row = 1, nNodO
            nini = m%ssh_stiff%rowptr_loc(row); nend = m%ssh_stiff%rowptr_loc(row+1) - 1
            ! S1: expected column set = vertices of the elements around row
            if (m%ssh_stiff%colind_loc(nini) /= row) ok_sp = .false.
            nexp = 0
            do k = 1, m%nod_in_elem2D_num(row)
                el = m%nod_in_elem2D(k, row)
                do j = 1, m%elem2D_nnodes(el)
                    c = m%elem2D_nodes(j, el)
                    if (mark(c) /= row) then
                        mark(c) = row; nexp = nexp + 1
                    end if
                end do
            end do
            if (nend - nini + 1 /= nexp) ok_sp = .false.
            do k = nini, nend
                c = m%ssh_stiff%colind_loc(k)
                if (mark(c) /= row) ok_sp = .false.   ! not a vertex, or a duplicate
                mark(c) = -row
            end do
            ! S2: row sum = areasvol/dt
            s = 0.0_WP; sa = 0.0_WP; offmax = 0.0_WP
            do k = nini, nend
                a = real(m%ssh_stiff%values(k), WP)
                s = s + a; sa = sa + abs(a)
                if (k > nini) offmax = max(offmax, abs(a))
            end do
            e_sum = max(e_sum, abs(s - mass(m, row)) / sa)
            ! S3: the diagonal vertex of every quad around row, relative to the row's
            ! largest off-diagonal entry
            do k = 1, m%nod_in_elem2D_num(row)
                el = m%nod_in_elem2D(k, row)
                if (m%elem2D_nnodes(el) /= 4) cycle
                do j = 1, 4
                    if (m%elem2D_nodes(j, el) == row) exit
                end do
                c = m%elem2D_nodes(mod(j+1, 4) + 1, el)
                nquad = nquad + 1
                a = 0.0_WP
                do nv = nini, nend
                    if (m%ssh_stiff%colind_loc(nv) == c) a = abs(real(m%ssh_stiff%values(nv), WP))
                end do
                dmin = min(dmin, a / offmax)
            end do
        end do
        e_sym = symmetry_error(m, p)
        call max_all(e_sum, p); dmin = -dmin; call max_all(dmin, p); dmin = -dmin
        call sum_all(nquad, p)
        if (nquad == 0) dmin = 0.0_WP
        k = merge(0, 1, ok_sp); call sum_all(k, p)
        if (p%mype == 0) write(*,'(a,a,i0,a,i0,a,3es10.2)') tag, ': nza(rank 0)=', &
            m%ssh_stiff%nza, ' quad corners=', nquad, &
            '  S2 row-sum err, S3 min |diag|/|offdiag|, S4 sym err:', e_sum, dmin, e_sym
        call check(k == 0, tag//': S1 rows = self first + all element vertices, no duplicates', p%mype)
        call check(e_sum < 1.0e-12_WP, tag//': S2 row sums = areasvol/dt', p%mype)
        if (nquad > 0) call check(dmin > 1.0e-3_WP, tag//': S3 quad diagonal entries present and nonzero', p%mype)
        call check(e_sym < 1.0e-12_WP, tag//': S4 matrix symmetric', p%mype)
    end subroutine

    real(kind=WP) function symmetry_error(m, p)
        ! max |a_ij - a_ji| / max |a_ij| (i /= j) over the global matrix, gathered on rank 0
        type(t_mesh),   intent(in) :: m
        type(t_partit), intent(in) :: p
        integer :: nNodO, nNodL, nEdgeO, nElemO, row, k, nz, ntot, ierr2, i, j, q
        integer, allocatable :: ri(:), ci(:), cnt(:), dsp(:), gri(:), gci(:), gptr(:), gcol(:), fill(:)
        real(kind=WP), allocatable :: v(:), gv(:), gval(:)
        real(kind=WP) :: emax, amax, aji
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        nz = m%ssh_stiff%nza
        allocate(ri(nz), ci(nz), v(nz))
        do row = 1, nNodO
            do k = m%ssh_stiff%rowptr_loc(row), m%ssh_stiff%rowptr_loc(row+1) - 1
                ri(k) = gid(p, row); ci(k) = gid(p, m%ssh_stiff%colind_loc(k))
                v(k)  = real(m%ssh_stiff%values(k), WP)
            end do
        end do
        allocate(cnt(p%npes), dsp(p%npes))
        if (p%npes > 1) then
            call MPI_Gather(nz, 1, MPI_INTEGER, cnt, 1, MPI_INTEGER, 0, p%MPI_COMM_FESOM, ierr2)
        else
            cnt(1) = nz
        end if
        dsp(1) = 0
        do i = 2, p%npes
            dsp(i) = dsp(i-1) + cnt(i-1)
        end do
        ntot = sum(cnt)
        if (p%mype == 0) then
            allocate(gri(ntot), gci(ntot), gv(ntot))
        else
            allocate(gri(1), gci(1), gv(1))
        end if
        if (p%npes > 1) then
            call MPI_Gatherv(ri, nz, MPI_INTEGER, gri, cnt, dsp, MPI_INTEGER, 0, p%MPI_COMM_FESOM, ierr2)
            call MPI_Gatherv(ci, nz, MPI_INTEGER, gci, cnt, dsp, MPI_INTEGER, 0, p%MPI_COMM_FESOM, ierr2)
            call MPI_Gatherv(v, nz, MPI_DOUBLE_PRECISION, gv, cnt, dsp, MPI_DOUBLE_PRECISION, 0, &
                             p%MPI_COMM_FESOM, ierr2)
        else
            gri = ri; gci = ci; gv = v
        end if
        emax = 0.0_WP
        if (p%mype == 0) then
            ! global CSR (rows in any order); every row must appear exactly once
            allocate(gptr(m%nod2D+1), fill(m%nod2D), gcol(ntot), gval(ntot))
            gptr = 0
            do k = 1, ntot
                gptr(gri(k)+1) = gptr(gri(k)+1) + 1
            end do
            gptr(1) = 1
            do i = 1, m%nod2D
                gptr(i+1) = gptr(i+1) + gptr(i)
            end do
            fill = gptr(1:m%nod2D)
            do k = 1, ntot
                i = gri(k)
                gcol(fill(i)) = gci(k); gval(fill(i)) = gv(k); fill(i) = fill(i) + 1
            end do
            amax = 0.0_WP
            do k = 1, ntot
                if (gri(k) /= gci(k)) amax = max(amax, abs(gv(k)))
            end do
            if (count(gri == gci) /= m%nod2D) emax = huge(1.0_WP)   ! a row missing or twice
            do k = 1, ntot
                i = gri(k); j = gci(k)
                if (i == j) cycle
                aji = huge(1.0_WP)
                do q = gptr(j), gptr(j+1) - 1
                    if (gcol(q) == i) aji = gval(q)
                end do
                emax = max(emax, abs(gv(k) - aji) / amax)
            end do
        end if
        if (p%npes > 1) call MPI_Bcast(emax, 1, MPI_DOUBLE_PRECISION, 0, p%MPI_COMM_FESOM, ierr2)
        symmetry_error = emax
    end function

    subroutine check_previous_triangles(m, p, tag)
        ! S5: the assembly before Task 9 (edge-only sparsity, fixed 12 slots, 3-slot
        ! gradient_sca(1:3|5:7)), transcribed; only on all-triangle meshes
        type(t_mesh),     intent(in) :: m
        type(t_partit),   intent(in) :: p
        character(len=*), intent(in) :: tag
        integer :: nNodO, nNodL, nEdgeO, nElemO, n, n1, n2, ed, i, row, k
        integer :: elnodes(3), el(2), npos(3), nquad
        integer, allocatable :: n_num(:), n_pos(:,:), rowptr(:), colind(:)
        real(kind=WP), allocatable :: values(:)
        real(kind=WP) :: factor, fy(3), zsrf, e_val, amax
        logical :: same_csr
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        nquad = count(m%elem2D_nnodes(1:nElemO) == 4); call sum_all(nquad, p)
        if (nquad > 0) return
        allocate(n_num(nNodL), n_pos(12, nNodO))
        n_pos = 0
        do n = 1, nNodO
            n_num(n) = 1; n_pos(1,n) = n
        end do
        do ed = 1, nEdgeO
            n1 = m%edges(1,ed); n2 = m%edges(2,ed)
            if (n1 <= nNodO) then
                n_pos(n_num(n1)+1, n1) = n2; n_num(n1) = n_num(n1)+1
            end if
            if (n2 <= nNodO) then
                n_pos(n_num(n2)+1, n2) = n1; n_num(n2) = n_num(n2)+1
            end if
        end do
        allocate(rowptr(nNodO+1))
        rowptr(1) = 1
        do n = 1, nNodO
            rowptr(n+1) = rowptr(n) + n_num(n)
        end do
        allocate(colind(rowptr(nNodO+1)-1), values(rowptr(nNodO+1)-1))
        values = 0.0_WP
        do n = 1, nNodO
            colind(rowptr(n):rowptr(n+1)-1) = n_pos(1:n_num(n), n)
        end do
        n_num  = 0
        factor = g*dt*alpha*theta
        do ed = 1, nEdgeO
            el = m%edge_tri(:, ed)
            do i = 1, 2
                if (i == 2 .and. el(2) == el(1)) cycle
                elnodes = m%elem2D_nodes(1:3, el(i))
                zsrf    = m%zbar(m%ulevels(el(i)))
                fy(1:3) = (m%zbar_e_bot(el(i)) - zsrf) * &
                          ( m%gradient_sca(1:3,el(i)) * m%edge_cross_dxdy(2*i  ,ed)  &
                           -m%gradient_sca(5:7,el(i)) * m%edge_cross_dxdy(2*i-1,ed) )
                if (i==2) fy = -fy
                row = m%edges(1, ed)
                if (row <= nNodO) then
                    do k = rowptr(row), rowptr(row+1)-1
                        n_num(colind(k)) = k
                    end do
                    npos = n_num(elnodes)
                    values(npos) = values(npos) + fy*factor
                end if
                row = m%edges(2, ed)
                if (row <= nNodO) then
                    do k = rowptr(row), rowptr(row+1)-1
                        n_num(colind(k)) = k
                    end do
                    npos = n_num(elnodes)
                    values(npos) = values(npos) - fy*factor
                end if
            end do
        end do
        do row = 1, nNodO
            values(rowptr(row)) = values(rowptr(row)) + mass(m, row)
        end do
        same_csr = size(colind) == m%ssh_stiff%nza
        if (same_csr) same_csr = all(rowptr == m%ssh_stiff%rowptr_loc) .and. &
                                 all(colind == m%ssh_stiff%colind_loc)
        e_val = huge(1.0_WP)
        if (same_csr) then
            amax  = maxval(abs(values))
            e_val = maxval(abs(real(m%ssh_stiff%values, WP) - values)) / amax
        end if
        call max_all(e_val, p)
        k = merge(0, 1, same_csr); call sum_all(k, p)
        if (p%mype == 0) write(*,'(a,a,es10.2)') tag, &
            ': S5 triangles vs previous assembly, max rel value diff:', e_val
        call check(k == 0, tag//': S5 triangle CSR equals the previous one', p%mype)
        call check(e_val < 1.0e-13_WP, tag//': S5 triangle values equal the previous ones', p%mype)
    end subroutine

    subroutine check_update(m, p, tag)
        ! S6: dhe = -H/2 adds K/2 (K = A - mass diagonal); values restored afterwards
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(in)    :: p
        character(len=*), intent(in)    :: tag
        integer :: nNodO, nNodL, nEdgeO, nElemO, e, row, k
        real(kind=MP), allocatable :: v0(:)
        real(kind=WP), allocatable :: kk(:)
        real(kind=WP) :: err
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        allocate(v0(m%ssh_stiff%nza), kk(m%ssh_stiff%nza))
        v0 = m%ssh_stiff%values
        kk = real(v0, WP)
        do row = 1, nNodO
            k = m%ssh_stiff%rowptr_loc(row)
            kk(k) = kk(k) - mass(m, row)
        end do
        do e = 1, nElemO
            m%dhe(e) = -0.5_MP * (m%zbar_e_bot(e) - m%zbar(m%ulevels(e)))
        end do
        call update_stiff_mat_ale(m, dt, p)
        err = maxval(abs(real(m%ssh_stiff%values, WP) - (real(v0, WP) + 0.5_WP*kk))) / maxval(abs(kk))
        call max_all(err, p)
        if (p%mype == 0) write(*,'(a,a,es10.2)') tag, ': S6 update_stiff_mat_ale(dhe=-H/2) - K/2:', err
        call check(err < 1.0e-12_WP, tag//': S6 update_stiff_mat_ale adds the 4-slot stiffness', p%mype)
        m%ssh_stiff%values = v0
        m%dhe = 0.0_MP
    end subroutine

    subroutine check_cg(m, p, tag)
        ! S7: b = A x for a smooth x (owned+halo from coordinates, no exchange), CG from 0
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(in)    :: p
        character(len=*), intent(in)    :: tag
        integer :: nNodO, nNodL, nEdgeO, nElemO, row, k, n_iter
        real(kind=WP), allocatable :: x(:)
        real(kind=WP) :: err, xmax
        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        allocate(x(nNodL), dyn%d_eta(nNodL), dyn%ssh_rhs(nNodL))
        do row = 1, nNodL
            x(row) = cos(3.0_WP*m%coord_nod2D(2, row)) * sin(5.0_WP*m%coord_nod2D(1, row)) &
                   + 0.3_WP*sin(11.0_WP*m%coord_nod2D(2, row))
        end do
        dyn%ssh_rhs = 0.0_WP; dyn%d_eta = 0.0_WP
        do row = 1, nNodO
            do k = m%ssh_stiff%rowptr_loc(row), m%ssh_stiff%rowptr_loc(row+1) - 1
                dyn%ssh_rhs(row) = dyn%ssh_rhs(row) + m%ssh_stiff%values(k)*x(m%ssh_stiff%colind_loc(k))
            end do
        end do
        dyn%solverinfo%soltol = 1.0e-12_WP
        call solve_ssh_ale(dyn, m, n_iter, p)
        err  = maxval(abs(dyn%d_eta(1:nNodO) - x(1:nNodO)))
        xmax = maxval(abs(x(1:nNodO)))
        call max_all(err, p); call max_all(xmax, p)
        err = err / xmax
        if (p%mype == 0) write(*,'(a,a,i0,a,es10.2)') tag, ': S7 CG iterations=', n_iter, &
            ' max rel error:', err
        call check(n_iter < dyn%solverinfo%maxiter, tag//': S7 CG converges', p%mype)
        call check(err < 1.0e-8_WP, tag//': S7 CG solution', p%mype)
    end subroutine

    subroutine stiffness_apply(m, f, nNodO, kf, kmax)
        ! kf = K f over the owned rows; kmax = max row sum of |K|
        type(t_mesh),  intent(in)  :: m
        real(kind=WP), intent(in)  :: f(:)
        integer,       intent(in)  :: nNodO
        real(kind=WP), intent(out) :: kf(:), kmax
        integer :: row, k
        real(kind=WP) :: a, sa
        kmax = 0.0_WP
        do row = 1, nNodO
            kf(row) = -mass(m, row) * f(row); sa = mass(m, row)
            do k = m%ssh_stiff%rowptr_loc(row), m%ssh_stiff%rowptr_loc(row+1) - 1
                a = real(m%ssh_stiff%values(k), WP)
                kf(row) = kf(row) + a * f(m%ssh_stiff%colind_loc(k))
                sa = sa + abs(a - merge(mass(m, row), 0.0_WP, k == m%ssh_stiff%rowptr_loc(row)))
            end do
            kmax = max(kmax, sa)
        end do
    end subroutine

    subroutine check_checkerboard_patch(m, nx, tag)
        ! S8 on the parallelogram patch: node (i,j) has colour (-1)**(i+j)
        type(t_mesh),     intent(in) :: m
        integer,          intent(in) :: nx
        character(len=*), intent(in) :: tag
        integer :: n, i, j
        real(kind=WP), allocatable :: f(:), s(:), kf(:)
        real(kind=WP) :: kmax, e_cb, e_st
        allocate(f(m%nod2D), s(m%nod2D), kf(m%nod2D))
        do n = 1, m%nod2D
            i = mod(n-1, nx) + 1; j = (n-1)/nx + 1
            f(n) = merge(1.0_WP, -1.0_WP, mod(i+j, 2) == 0)
            s(n) = merge(1.0_WP, -1.0_WP, mod(i, 2) == 0)
        end do
        call stiffness_apply(m, f, m%nod2D, kf, kmax)
        e_cb = maxval(abs(kf)) / kmax
        call stiffness_apply(m, s, m%nod2D, kf, kmax)
        e_st = maxval(abs(kf)) / kmax
        write(*,'(a,a,2es10.2)') tag, ': S8 |K checkerboard|/|K|, |K stripes|/|K|:', e_cb, e_st
        call check(e_cb < 1.0e-12_WP, tag//': S8 checkerboard is in the null space of K', 0)
        call check(e_st > 1.0e-2_WP,  tag//': S8 stripes are not (non-vacuity)', 0)
    end subroutine

    subroutine check_quadratic_patch(m, nx, ny, tag)
        ! S9: f = a*x^2 + b*x*y + c*y^2 (x, y in metres); interior node (i,j), 1 < i < nx,
        ! 1 < j < ny; H = the (flat) element depth -zbar_e_bot
        type(t_mesh),     intent(in) :: m
        integer,          intent(in) :: nx, ny
        character(len=*), intent(in) :: tag
        real(kind=WP), parameter :: a = 1.3e-6_WP, b = -0.7e-6_WP, c = 0.4e-6_WP
        integer :: n, i, j
        real(kind=WP), allocatable :: f(:), kf(:)
        real(kind=WP) :: x, y, kmax, ex, err, h
        allocate(f(m%nod2D), kf(m%nod2D))
        do n = 1, m%nod2D
            x = m%coord_nod2D(1, n)*r_earth; y = m%coord_nod2D(2, n)*r_earth
            f(n) = a*x*x + b*x*y + c*y*y
        end do
        h = -real(m%zbar_e_bot(1), WP)
        call check(all(m%zbar_e_bot == m%zbar_e_bot(1)), tag//': S9 flat bottom', 0)
        call stiffness_apply(m, f, m%nod2D, kf, kmax)
        err = 0.0_WP
        do n = 1, m%nod2D
            i = mod(n-1, nx) + 1; j = (n-1)/nx + 1
            if (i == 1 .or. i == nx .or. j == 1 .or. j == ny) cycle
            ex = -g*dt*alpha*theta*h*(2.0_WP*a + 2.0_WP*c)*real(m%area(n), WP)
            err = max(err, abs(kf(n) - ex)/abs(ex))
        end do
        write(*,'(a,a,es10.2)') tag, ': S9 interior K f vs -g*dt*alpha*theta*H*lap(f)*area, rel:', err
        call check(err < 1.0e-10_WP, tag//': S9 K f == -g*dt*alpha*theta*H*lap(f)*area (quadratic f)', 0)
    end subroutine

    subroutine check_checkerboard_mesh(m, p, tag)
        ! S8 on the mesh read (np 1): 2-colour the nodes along the edges (breadth first);
        ! when every element is a quad and the colouring is consistent, K f = 0
        type(t_mesh),     intent(in) :: m
        type(t_partit),   intent(in) :: p
        character(len=*), intent(in) :: tag
        integer, allocatable :: col(:), queue(:), nptr(:), nadj(:), fillp(:)
        real(kind=WP), allocatable :: f(:), kf(:)
        integer :: n, ed, a, b, head, tail, k, start
        logical :: bip, expect
        character(len=8) :: ebuf
        real(kind=WP) :: kmax, e_cb
        ebuf = ''
        call get_environment_variable('FESOM3_EXPECT_CHECKERBOARD', ebuf)
        expect = trim(ebuf) == '1'
        if (any(m%elem2D_nnodes /= 4)) then
            if (expect) call check(.false., tag//': S8 expected an all-quad mesh', p%mype)
            return
        end if
        allocate(nptr(m%nod2D+1), fillp(m%nod2D), nadj(2*m%edge2D))
        nptr = 0
        do ed = 1, m%edge2D
            nptr(m%edges(1,ed)+1) = nptr(m%edges(1,ed)+1) + 1
            nptr(m%edges(2,ed)+1) = nptr(m%edges(2,ed)+1) + 1
        end do
        nptr(1) = 1
        do n = 1, m%nod2D
            nptr(n+1) = nptr(n+1) + nptr(n)
        end do
        fillp = nptr(1:m%nod2D)
        do ed = 1, m%edge2D
            a = m%edges(1,ed); b = m%edges(2,ed)
            nadj(fillp(a)) = b; fillp(a) = fillp(a) + 1
            nadj(fillp(b)) = a; fillp(b) = fillp(b) + 1
        end do
        allocate(col(m%nod2D), queue(m%nod2D))
        col = 0; bip = .true.
        do start = 1, m%nod2D
            if (col(start) /= 0) cycle
            col(start) = 1; head = 1; tail = 1; queue(1) = start
            do while (head <= tail)
                a = queue(head); head = head + 1
                do k = nptr(a), nptr(a+1) - 1
                    b = nadj(k)
                    if (col(b) == 0) then
                        col(b) = -col(a); tail = tail + 1; queue(tail) = b
                    else if (col(b) == col(a)) then
                        bip = .false.
                    end if
                end do
            end do
        end do
        if (.not. bip) then
            write(*,'(a,a)') tag, ': S8 all-quad mesh is not 2-colourable along its edges (skipped)'
            if (expect) call check(.false., tag//': S8 expected a 2-colourable all-quad mesh', p%mype)
            return
        end if
        allocate(f(m%nod2D), kf(m%nod2D))
        f = real(col, WP)
        call stiffness_apply(m, f, m%nod2D, kf, kmax)
        e_cb = maxval(abs(kf)) / kmax
        write(*,'(a,a,es10.2)') tag, ': S8 all-quad mesh, |K checkerboard|/|K|:', e_cb
        call check(e_cb < 1.0e-12_WP, tag//': S8 checkerboard is in the null space of K', p%mype)
    end subroutine

end program test_ssh_stiff
