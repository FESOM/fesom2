program test_ice_mass
    ! Sea-ice FCT mass matrix and FCT advection on mixed triangle/quad meshes (plan
    ! 2026-10-08-hybrid-meshes, Task 16, the user's specification): consistent mass matrix
    ! from the element contributions (triangle S/6 diagonal, S/12 off-diagonal; quad row j
    ! (S/36)*(4,2,1,2) rotated), lumped mass = its row sums replacing the node area in every
    ! FCT step, antidiffusive operator consistent - lumped, TG rhs with the element mass
    ! contributions and the kite-weighted mean velocity. Mesh from FESOM3_MESH_DIR (default
    ! the NEMO-mixed cut-out, unrotated; FESOM3_TEST_ROTATE=1 for pi/CORE2) at 1 or N ranks
    ! and, at 1 rank, the analytic triangle, mixed and parallelogram patches (planar). Node
    ! fields are hashes of the global id (identical on owned and halo copies):
    !   F1 fct_massmatrix == an independent assembly of the element contributions
    !      (coefficient by the cyclic vertex distance)
    !   F2 row sums == fct_lumpedmass == sum_e S_e/nv_e (owned and, after the exchange,
    !      halo nodes); sum of the lumped mass == the total (wet) element area; triangle
    !      meshes: lumped mass == the node area
    !   F3 a constant ice field at rest stays constant (with ice_diff /= 0): TG rhs 0,
    !      low order (valuesl) == c, high order (dvalues) 0, final == c; planar patches with
    !      a uniform velocity: constant at the nodes 3+ rings away from the boundary
    !   F4 sum(lumped*f) is conserved by TG rhs + FCT for the three ice tracers with a
    !      non-uniform velocity and fields with zeros (limiter active)
    !   F5 triangle meshes: TG rhs, high order, low order and the final fields == the
    !      former FCT (FESOM2 formulas: S/12 mass, node area, (-2,1,1)/12, sum(u)/3)
    !   F6 boundedness of the F4 run (all meshes): the final field of every owned node
    !      lies within the min/max of the low-order and the old field over its CSR row
    !      (the node and its element neighbours, the limiter's cluster)
    ! FESOM3_TEST_MISSING_COL=1: a column removed from the CSR must stop
    ! ice_mass_matrix_fill with a message (ctest regex; reaching the end = 'did NOT stop').
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all, sum_all, hval, set_depths
    use mod_precision,     only: WP, MP
    use mod_mesh,          only: t_mesh, MAX_NV
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex
    use mod_mesh_analytic, only: generate_analytic_mesh, generate_mixed_patch, generate_quad_patch
    use mod_part_bounds,   only: owned_bounds, is_multirank
    use mod_halo,          only: exchange_nod
    use mod_param_phys,    only: scale_area
    use mod_ice,           only: t_ice
    use mod_ice_setup,     only: ice_allocate, ice_mass_matrix_fill
    use mod_ice_fct,       only: ice_TG_rhs, ice_fct_solve, ice_mcoef, ice_acoef
    use oce_ssh_rhs,       only: init_stiff_mat_ale
    implicit none

    real(kind=WP), parameter :: dt = 1800.0_WP
    character(len=512) :: cbuf
    type(t_partit) :: partit, tpartit, apartit, qpartit
    type(t_mesh)   :: mesh, tmesh, amesh, qmesh

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.false., partition_np1=.false.)

    cbuf = ''
    call get_environment_variable('FESOM3_TEST_MISSING_COL', cbuf)
    if (trim(cbuf) == '1') call missing_column_case(mesh, partit)

    call check_tables()
    call check_fct(mesh, partit, 'mesh', .false.)
    if (partit%npes == 1) then
        call generate_analytic_mesh(tmesh, tpartit, nx=13, ny=11, nl=4, &
                                    Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_fct(tmesh, tpartit, 'patch_tri', .true.)
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=4, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_fct(amesh, apartit, 'patch_mixed', .true.)
        call generate_quad_patch(qmesh, qpartit, nx=13, ny=11, nl=4, &
                                 Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP, skew=0.3_WP)
        call check_fct(qmesh, qpartit, 'patch_quad', .true.)
    end if

    call finish(partit, 'test_ice_mass')

contains

    ! element mass coefficient (per unit area) of the vertex pair (j,k), by the cyclic
    ! vertex distance -- independent of the mod_ice_fct table
    real(kind=WP) function mcoef(j, k, nv)
        integer, intent(in) :: j, k, nv
        if (nv == 3) then
            mcoef = 1.0_WP/12.0_WP
            if (j == k) mcoef = 1.0_WP/6.0_WP
        else
            select case (mod(k - j + 4, 4))
            case (0);    mcoef = 4.0_WP/36.0_WP
            case (2);    mcoef = 1.0_WP/36.0_WP
            case default; mcoef = 2.0_WP/36.0_WP
            end select
        end if
    end function mcoef

    subroutine check_tables()
        ! F0 the mod_ice_fct tables: ice_mcoef == the cyclic-distance coefficients over the
        ! nv unique vertices (triangle slot 4 = 0), row sums 1/3 and 1/4, ice_acoef ==
        ! ice_mcoef - its row sums on the diagonal
        integer :: nv, j, k
        real(kind=WP) :: e0(3), ex
        e0 = 0.0_WP
        do nv = 3, 4
            do j = 1, MAX_NV
                do k = 1, MAX_NV
                    ex = 0.0_WP
                    if (j <= nv .and. k <= nv) ex = mcoef(j, k, nv)
                    e0(1) = max(e0(1), abs(ice_mcoef(j, k, nv) - ex))
                    if (j == k) ex = ex - sum(ice_mcoef(j, :, nv))
                    e0(3) = max(e0(3), abs(ice_acoef(j, k, nv) - ex))
                end do
                if (j <= nv) e0(2) = max(e0(2), abs(sum(ice_mcoef(j, :, nv)) - 1.0_WP/real(nv, WP)))
            end do
        end do
        call check(e0(1) <= 1.0e-16_WP, 'F0 ice_mcoef == element contributions', partit%mype)
        call check(e0(2) <= 1.0e-16_WP, 'F0 ice_mcoef row sums == 1/3, 1/4', partit%mype)
        call check(e0(3) <= 1.0e-16_WP, 'F0 ice_acoef == consistent - lumped', partit%mype)
        if (partit%mype == 0) write(*,'(a,3es9.2)') '  F0 tables mcoef/rowsum/acoef ', e0
    end subroutine check_tables

    subroutine setup_ice(ice, m, p)
        type(t_ice),    intent(inout) :: ice
        type(t_mesh),   intent(inout) :: m
        type(t_partit), intent(in)    :: p
        call ice_allocate(ice, m, p)
        ice%ice_dt        = dt
        ice%ice_diff      = 10.0_WP     ! the t_ice default (CORE2 runs 0): exercises diffusion
        ice%ice_gamma_fct = 0.5_WP
        call ice_mass_matrix_fill(ice, m, p)
    end subroutine

    subroutine check_fct(m, p, tag, planar)
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(inout) :: p
        character(len=*), intent(in)    :: tag
        logical,          intent(in)    :: planar
        type(t_ice) :: ice
        integer :: nNodO, nNodL, nEdgeO, nElemO, n, e, j, k, nv, row, ipos, nini, nend, ed, r
        integer :: nquad, ntri, it
        integer, allocatable :: ngid(:), dist(:)
        logical :: found, alltri
        real(kind=WP), allocatable :: ex(:), lex(:), f0(:,:), rhs(:,:), fl(:,:), df(:,:), f1(:,:)
        real(kind=WP), allocatable :: rrhs(:), rdf(:), rfl(:), rf1(:)
        real(kind=WP) :: f1e, f2(4), f3(5), f4(3), f4m, f5(4), s(2), cmax, sc, tmp(1), f6, bmin, bmax
        real(kind=WP), parameter :: cst(3) = [0.7_WP, 1.3_WP, 0.2_WP]

        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        call set_depths(m, p)
        call init_stiff_mat_ale(m, dt, p)
        call setup_ice(ice, m, p)
        allocate(ngid(nNodL))
        do n = 1, nNodL
            ngid(n) = n
            if (p%npes > 1) ngid(n) = p%myList_nod2D(n)
        end do
        ntri  = count(m%elem2D_nnodes(1:nElemO) == 3 .and. m%elem2D_nodes(1, 1:nElemO) <= nNodO)
        nquad = count(m%elem2D_nnodes(1:nElemO) == 4 .and. m%elem2D_nodes(1, 1:nElemO) <= nNodO)
        call sum_all(nquad, p); call sum_all(ntri, p)
        alltri = (nquad == 0)

        !--- F1: matrix entries vs an independent assembly ---
        associate(st => m%ssh_stiff)
        allocate(ex(st%nza)); ex = 0.0_WP
        found = .true.
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            nv = m%elem2D_nnodes(e)
            do j = 1, nv
                row = m%elem2D_nodes(j, e)
                if (row > nNodO) cycle
                do k = 1, nv
                    ipos = findloc(st%colind_loc(st%rowptr_loc(row):st%rowptr_loc(row+1)-1), &
                                   m%elem2D_nodes(k, e), 1)
                    if (ipos == 0) then
                        found = .false.; cycle
                    end if
                    ipos = ipos + st%rowptr_loc(row) - 1
                    ex(ipos) = ex(ipos) + real(m%elem_area(e), WP)*mcoef(j, k, nv)
                end do
            end do
        end do
        call check(found, tag//' F1 all element vertex pairs in the CSR', p%mype)
        f1e = maxval(abs(ice%work%fct_massmatrix(1:st%nza) - ex))/maxval(abs(ex))
        call max_all(f1e, p)

        !--- F2: lumped mass ---
        allocate(lex(nNodL)); lex = 0.0_WP
        f2 = 0.0_WP; s = 0.0_WP
        do n = 1, nNodO
            nini = st%rowptr_loc(n); nend = st%rowptr_loc(n+1) - 1
            f2(1) = max(f2(1), abs(sum(ice%work%fct_massmatrix(nini:nend)) - ice%work%fct_lumpedmass(n)) &
                               /ice%work%fct_lumpedmass(n))
            do j = 1, m%nod_in_elem2D_num(n)
                e = m%nod_in_elem2D(j, n)
                if (m%ulevels(e) > 1) cycle
                if (m%elem2D_nnodes(e) == 3) then
                    lex(n) = lex(n) + real(m%elem_area(e), WP)/3.0_WP
                else
                    lex(n) = lex(n) + real(m%elem_area(e), WP)/4.0_WP
                end if
            end do
            s(1) = s(1) + ice%work%fct_lumpedmass(n)
        end do
        do e = 1, nElemO     ! each element once: on the rank that owns its 1st vertex
            if (m%ulevels(e) == 1 .and. m%elem2D_nodes(1, e) <= nNodO) s(2) = s(2) + real(m%elem_area(e), WP)
        end do
        if (is_multirank(p)) call exchange_nod(lex, p)
        do n = 1, nNodL
            f2(2) = max(f2(2), abs(ice%work%fct_lumpedmass(n) - lex(n))/lex(n))
            if (alltri) f2(4) = max(f2(4), abs(ice%work%fct_lumpedmass(n) - real(m%area(n), WP))/lex(n))
        end do
        call sum_all(s, p)
        f2(3) = abs(s(1) - s(2))/s(2)
        do k = 1, 4
            call max_all(f2(k), p)
        end do
        end associate
        call check(f1e <= 1.0e-14_WP, tag//' F1 fct_massmatrix == element contributions', p%mype)
        call check(f2(1) <= 1.0e-14_WP, tag//' F2 row sums == fct_lumpedmass', p%mype)
        call check(f2(2) <= 1.0e-14_WP, tag//' F2 fct_lumpedmass == sum S/nv (owned + halo)', p%mype)
        call check(f2(3) <= 1.0e-13_WP, tag//' F2 total lumped mass == total area', p%mype)
        call check(f2(4) <= 1.0e-13_WP, tag//' F2 triangles: lumped mass == node area', p%mype)

        !--- F3: constant field at rest (and uniform drift on the patches) ---
        allocate(f0(nNodL,3), rhs(nNodL,3), fl(nNodL,3), df(nNodL,3), f1(nNodL,3))
        f3 = 0.0_WP
        do it = 1, merge(2, 1, planar)
            do k = 1, 3
                ice%data(k)%values = cst(k)
            end do
            ice%uice = 0.0_WP; ice%vice = 0.0_WP
            if (it == 2) then
                ice%uice = 0.08_WP; ice%vice = -0.05_WP
            end if
            call ice_TG_rhs(ice, m, p)
            call ice_fct_solve(ice, m, p)
            if (it == 1) then
                do k = 1, 3
                    do n = 1, nNodO
                        if (m%ulevels_nod2D(n) > 1) cycle
                        f3(1) = max(f3(1), abs(ice%data(k)%values_rhs(n))/(cst(k)*ice%work%fct_lumpedmass(n)))
                        f3(2) = max(f3(2), abs(ice%data(k)%valuesl(n) - cst(k))/cst(k))
                        f3(3) = max(f3(3), abs(ice%data(k)%dvalues(n))/cst(k))
                        f3(4) = max(f3(4), abs(ice%data(k)%values(n) - cst(k))/cst(k))
                    end do
                end do
            else
                ! ring distance from the boundary over the CSR neighbours (1 rank)
                allocate(dist(nNodL)); dist = huge(1)
                do ed = m%edge2D_in + 1, nEdgeO
                    dist(m%edges(1:2, ed)) = 0
                end do
                do r = 0, 3
                    do n = 1, nNodO
                        if (dist(n) /= r) cycle
                        do ipos = m%ssh_stiff%rowptr_loc(n), m%ssh_stiff%rowptr_loc(n+1) - 1
                            j = m%ssh_stiff%colind_loc(ipos)
                            dist(j) = min(dist(j), r + 1)
                        end do
                    end do
                end do
                cmax = 0.0_WP
                do k = 1, 3
                    do n = 1, nNodO
                        if (dist(n) >= 3) f3(5) = max(f3(5), abs(ice%data(k)%values(n) - cst(k))/cst(k))
                        cmax = max(cmax, abs(ice%data(k)%values(n) - cst(k))/cst(k))
                    end do
                end do
                call check(cmax > 1.0e-6_WP, tag//' F3 uniform drift changes the boundary nodes (non-vacuity)', p%mype)
                deallocate(dist)
            end if
        end do
        do k = 1, 4
            call max_all(f3(k), p)
        end do
        call check(f3(1) <= 1.0e-13_WP, tag//' F3 constant at rest: TG rhs 0', p%mype)
        call check(f3(2) <= 1.0e-13_WP, tag//' F3 constant at rest: low order == c', p%mype)
        call check(f3(3) <= 1.0e-13_WP, tag//' F3 constant at rest: high order 0', p%mype)
        call check(f3(4) <= 1.0e-13_WP, tag//' F3 constant at rest: final == c', p%mype)
        if (planar) call check(f3(5) <= 1.0e-12_WP, tag//' F3 uniform drift: interior (ring >= 3) == c', p%mype)

        !--- F4: conservation of sum(lumped*f) ---
        do n = 1, nNodL
            ice%data(1)%values(n) = 0.5_WP + 0.4_WP*hval(ngid(n), 1)
            ice%data(2)%values(n) = max(0.0_WP, 1.5_WP*hval(ngid(n), 2))
            ice%data(3)%values(n) = 0.3_WP*(1.0_WP + hval(ngid(n), 3))
            ice%uice(n) = 0.1_WP*hval(ngid(n), 4)
            ice%vice(n) = 0.1_WP*hval(ngid(n), 5)
        end do
        do k = 1, 3
            f0(:,k) = ice%data(k)%values
        end do
        call ice_TG_rhs(ice, m, p)
        call ice_fct_solve(ice, m, p)
        f4 = 0.0_WP; f4m = 0.0_WP
        do k = 1, 3
            rhs(:,k) = ice%data(k)%values_rhs; fl(:,k) = ice%data(k)%valuesl
            df(:,k)  = ice%data(k)%dvalues;    f1(:,k) = ice%data(k)%values
            s = 0.0_WP
            do n = 1, nNodO
                if (m%ulevels_nod2D(n) > 1) cycle
                s(1) = s(1) + ice%work%fct_lumpedmass(n)*f0(n,k)
                s(2) = s(2) + ice%work%fct_lumpedmass(n)*f1(n,k)
                f4m  = max(f4m, abs(f1(n,k) - f0(n,k)))
            end do
            call sum_all(s, p)
            f4(k) = abs(s(2) - s(1))/s(1)
        end do
        call max_all(f4m, p)
        call check(maxval(f4) <= 1.0e-13_WP, tag//' F4 sum(lumped*f) conserved (a, m, m_snow)', p%mype)
        call check(f4m > 1.0e-6_WP, tag//' F4 fields advected (non-vacuity)', p%mype)

        !--- F6: the FCT result within the low-order / old-field bounds of its cluster ---
        f6 = 0.0_WP
        associate(st => m%ssh_stiff)
        do k = 1, 3
            sc = maxval(abs(f0(1:nNodO,k))); tmp(1) = sc; call max_all(tmp(1), p); sc = tmp(1)
            do n = 1, nNodO
                if (m%ulevels_nod2D(n) > 1) cycle
                nini = st%rowptr_loc(n) - st%rowptr_loc(1) + 1
                nend = st%rowptr_loc(n+1) - st%rowptr_loc(1)
                bmax = max(maxval(fl(st%colind_loc(nini:nend),k)), maxval(f0(st%colind_loc(nini:nend),k)))
                bmin = min(minval(fl(st%colind_loc(nini:nend),k)), minval(f0(st%colind_loc(nini:nend),k)))
                f6 = max(f6, (f1(n,k) - bmax)/sc, (bmin - f1(n,k))/sc)
            end do
        end do
        end associate
        call max_all(f6, p)
        call check(f6 <= 1.0e-13_WP, tag//' F6 FCT result within the low-order/old bounds (a, m, m_snow)', p%mype)

        !--- F5: triangle meshes == the former FCT ---
        f5 = 0.0_WP
        if (alltri) then
            allocate(rrhs(nNodL), rdf(nNodL), rfl(nNodL), rf1(nNodL))
            do k = 1, 3
                call former_fct(m, p, ice, f0(:,k), rrhs, rdf, rfl, rf1)
                sc = maxval(abs(rrhs(1:nNodO))); tmp(1) = sc; call max_all(tmp(1), p); sc = tmp(1)
                f5(1) = max(f5(1), maxval(abs(rhs(1:nNodO,k) - rrhs(1:nNodO)))/sc)
                sc = maxval(abs(f0(1:nNodO,k))); tmp(1) = sc; call max_all(tmp(1), p); sc = tmp(1)
                f5(2) = max(f5(2), maxval(abs(df(1:nNodO,k) - rdf(1:nNodO)))/sc)
                f5(3) = max(f5(3), maxval(abs(fl(1:nNodO,k) - rfl(1:nNodO)))/sc)
                f5(4) = max(f5(4), maxval(abs(f1(1:nNodO,k) - rf1(1:nNodO)))/sc)
            end do
            do k = 1, 4
                call max_all(f5(k), p)
            end do
            call check(f5(1) <= 1.0e-12_WP, tag//' F5 triangles: TG rhs == former', p%mype)
            call check(f5(2) <= 1.0e-12_WP, tag//' F5 triangles: high order == former', p%mype)
            call check(f5(3) <= 1.0e-12_WP, tag//' F5 triangles: low order == former', p%mype)
            call check(f5(4) <= 1.0e-12_WP, tag//' F5 triangles: final == former', p%mype)
        end if

        if (p%mype == 0) then
            write(*,'(a,i0,a,i0,a,i0)') '  '//tag//': triangles ', ntri, ' quads ', nquad, ' npes ', p%npes
            write(*,'(a,es9.2,a,4es9.2)') '    F1 entries ', f1e, '  F2 rowsum/lumped/total/area ', f2
            write(*,'(a,5es9.2)') '    F3 rest rhs/low/high/final, drift interior ', f3
            write(*,'(a,3es9.2,a,es9.2)') '    F4 conservation a/m/s ', f4, '  max|change| ', f4m
            write(*,'(a,es9.2)') '    F6 max bound excess (<= 0 when bounded) ', f6
            if (alltri) write(*,'(a,4es9.2)') '    F5 vs former rhs/high/low/final ', f5
        end if
    end subroutine check_fct

    subroutine former_fct(m, p, ice, f, rhs, df, fl, fnew)
        ! the ice FCT for one tracer before Task 16 (FESOM2 formulas, triangles only): mass
        ! matrix S/12 (+S/12 on the diagonal), node area as the lumped mass, TG rhs with
        ! sum(u)/12 and the /9 mean, antidiffusive operator (-2,1,1)/12. Same velocity,
        ! ice_dt, ice_diff and gamma as ice%; results on owned nodes (exchanged like FESOM2).
        type(t_mesh),   intent(in)  :: m
        type(t_partit), intent(in)  :: p
        type(t_ice),    intent(in)  :: ice
        real(kind=WP),  intent(in)  :: f(:)
        real(kind=WP),  intent(out) :: rhs(:), df(:), fl(:), fnew(:)
        integer :: nNodO, nNodL, nEdgeO, nElemO, e, n, q, row, ipos, nini, nend, el(3), it, icoef(3,3)
        logical :: lmr
        real(kind=WP) :: entries(3), dx(3), dy(3), um, vm, vol, diff, gamma, ae, flux
        real(kind=WP), allocatable :: mm(:), tmp(:), fx(:,:), tmax(:), tmin(:), pp(:), pm(:), ar(:)

        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        lmr = is_multirank(p)
        gamma = ice%ice_gamma_fct
        allocate(mm(m%ssh_stiff%nza), tmp(nNodL), fx(nElemO,3), ar(nNodL))
        allocate(tmax(nNodL), tmin(nNodL), pp(nNodL), pm(nNodL))
        ar = real(m%area(1:nNodL), WP)
        associate(rp => m%ssh_stiff%rowptr_loc, ci => m%ssh_stiff%colind_loc)
        mm = 0.0_WP
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el = m%elem2D_nodes(1:3, e)
            do n = 1, 3
                row = el(n)
                if (row > nNodO) cycle
                do q = 1, 3
                    do ipos = rp(row), rp(row+1) - 1
                        if (ci(ipos) == el(q)) exit
                    end do
                    mm(ipos) = mm(ipos) + m%elem_area(e)/12.0_WP
                    if (q == n) mm(ipos) = mm(ipos) + m%elem_area(e)/12.0_WP
                end do
            end do
        end do
        ! TG rhs
        rhs = 0.0_WP
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el  = m%elem2D_nodes(1:3, e)
            dx  = m%gradient_sca(1:3, e)
            dy  = m%gradient_sca(5:7, e)
            vol = m%elem_area(e)
            um  = sum(ice%uice(el))
            vm  = sum(ice%vice(el))
            diff = ice%ice_diff*sqrt(m%elem_area(e)/scale_area)
            do n = 1, 3
                row = el(n)
                do q = 1, 3
                    entries(q) = vol*ice%ice_dt*((dx(n)*(um+ice%uice(el(q))) + &
                                 dy(n)*(vm+ice%vice(el(q))))/12.0_WP - &
                                 diff*(dx(n)*dx(q) + dy(n)*dy(q)) - &
                                 0.5_WP*ice%ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)
                end do
                rhs(row) = rhs(row) + sum(entries*f(el))
            end do
        end do
        ! high order (3 Jacobi-style iterations)
        df = 0.0_WP
        do row = 1, nNodO
            df(row) = rhs(row)/ar(row)
        end do
        if (lmr) call exchange_nod(df, p)
        do it = 1, 2
            do row = 1, nNodO
                nini = rp(row); nend = rp(row+1) - 1
                tmp(row) = df(row) + (rhs(row) - sum(mm(nini:nend)*df(ci(nini:nend))))/ar(row)
            end do
            df(1:nNodO) = tmp(1:nNodO)
            if (lmr) call exchange_nod(df, p)
        end do
        ! low order
        fl = 0.0_WP
        do row = 1, nNodO
            nini = rp(row); nend = rp(row+1) - 1
            fl(row) = (rhs(row) + gamma*sum(mm(nini:nend)*f(ci(nini:nend))))/ar(row) + (1.0_WP-gamma)*f(row)
        end do
        if (lmr) call exchange_nod(fl, p)
        ! Zalesak limiter
        icoef = 1
        do n = 1, 3
            icoef(n,n) = -2
        end do
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el  = m%elem2D_nodes(1:3, e)
            vol = m%elem_area(e)
            do q = 1, 3
                fx(e,q) = -sum(icoef(:,q)*(gamma*f(el) + df(el)))*(vol/ar(el(q)))/12.0_WP
            end do
        end do
        tmax = 0.0_WP; tmin = 0.0_WP
        do row = 1, nNodO
            nini = rp(row); nend = rp(row+1) - 1
            tmax(row) = max(maxval(fl(ci(nini:nend))), maxval(f(ci(nini:nend)))) - fl(row)
            tmin(row) = min(minval(fl(ci(nini:nend))), minval(f(ci(nini:nend)))) - fl(row)
        end do
        pp = 0.0_WP; pm = 0.0_WP
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el = m%elem2D_nodes(1:3, e)
            do q = 1, 3
                if (fx(e,q) > 0) then
                    pp(el(q)) = pp(el(q)) + fx(e,q)
                else
                    pm(el(q)) = pm(el(q)) + fx(e,q)
                end if
            end do
        end do
        do n = 1, nNodO
            flux = pp(n)
            pp(n) = 0.0_WP
            if (abs(flux) > 0) pp(n) = min(1.0_WP, tmax(n)/max(flux, 1.e-12))
            flux = pm(n)
            pm(n) = 0.0_WP
            if (abs(flux) > 0) pm(n) = min(1.0_WP, tmin(n)/min(flux, -1.e-12))
        end do
        if (lmr) then
            call exchange_nod(pm, p)
            call exchange_nod(pp, p)
        end if
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el = m%elem2D_nodes(1:3, e)
            ae = 1.0_WP
            do q = 1, 3
                if (fx(e,q) >= 0.0_WP) ae = min(ae, pp(el(q)))
                if (fx(e,q) <  0.0_WP) ae = min(ae, pm(el(q)))
            end do
            fx(e,:) = ae*fx(e,:)
        end do
        fnew = f
        fnew(1:nNodO) = fl(1:nNodO)
        do e = 1, nElemO
            if (m%ulevels(e) > 1) cycle
            el = m%elem2D_nodes(1:3, e)
            do q = 1, 3
                fnew(el(q)) = fnew(el(q)) + fx(e,q)
            end do
        end do
        if (lmr) call exchange_nod(fnew, p)
        end associate
    end subroutine former_fct

    subroutine missing_column_case(m, p)
        ! remove the first off-diagonal column of the first owned row from the CSR: the
        ! mass-matrix fill must stop with a message
        type(t_mesh),   intent(inout) :: m
        type(t_partit), intent(inout) :: p
        type(t_ice) :: ice
        call set_depths(m, p)
        call init_stiff_mat_ale(m, dt, p)
        m%ssh_stiff%colind_loc(m%ssh_stiff%rowptr_loc(1) + 1) = -1
        call ice_allocate(ice, m, p)
        call ice_mass_matrix_fill(ice, m, p)
        write(*,'(a)') 'test_ice_mass: ice_mass_matrix_fill did NOT stop on a missing CSR column'
        call par_ex(p%MPI_COMM_FESOM, p%mype)
        error stop 1
    end subroutine missing_column_case

end program test_ice_mass
