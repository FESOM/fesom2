program test_pgf_mixed
    ! Ocean dynamics kernels on mixed triangle/quad meshes (plan 2026-10-08-hybrid-meshes,
    ! Task 10): the unified 4-slot formulas (gradient_sca 1:4 / 5:8, element means with the
    ! kite weights elem_w, node gathers with nod_in_elem2D_w) on the mesh read from
    ! FESOM3_MESH_DIR (default the NEMO-mixed cut-out, unrotated) at 1 or N ranks and, at 1
    ! rank, on the analytic mixed patch (generate_mixed_patch: jittered quads + triangle
    ! pairs). Geometry with cartesian=.true., so fields linear in (x, y) = r_earth*(lon, lat)
    ! are linear in the element metric and every assertion below is exact to roundoff:
    !   P1 linfs full-cell PGF, resting ocean: hpressure depends on the level only
    !      -> pgf_x = pgf_y = 0
    !   P2 linfs full-cell PGF, hpressure = p(z) + a*x + b*y -> pgf = (a, b)/density_0
    !   P3 Shchepetkin PGF, resting stratified ocean on TILTED levels: Z_3d_n linear in
    !      (x, y) on every level, density linear in z -> pgf = 0 (needs the 4-slot
    !      gradients AND the vertex drho_dz averaged with weights summing to 1: a /3 over
    !      the 4 slots of a quad leaves (1/3)*rho_z*dz/dx)
    !   P4 Shchepetkin PGF, flat levels, density = c*x (no z dependence) -> pgf_x(nz) =
    !      g/density_0 * c * (depth of the mid-layer below the element top), pgf_y = 0
    !   N1 compute_vel_nodes: a uniform velocity gathers to the same velocity at every node
    !      level where all adjacent elements are wet (sum_k nod_in_elem2D_w*elem_area ==
    !      area; the former 3*area denominator gives 4/3 at an all-quad node)
    !   N2 update_vel: d_eta = a*x + b*y -> UV increment = -g*theta*dt*(a, b)
    !   R1 compute_sigma_xy (kite-share node gather of the 4-slot element gradients): T and
    !      S linear in (x, y) -> sigma_xy = (-alpha*grad T + beta*grad S)*density_0 at every
    !      owned node level (each wet element gradient is exact, the shares sum to vol)
    !   R2 smooth_nod (kite-weighted element means, kite-share node means): a constant
    !      field stays constant over 3 sweeps
    !   R3 fer_gamma2vel: Gamma linear in zbar, horizontally uniform -> fer_uv = the
    !      vertical derivative (elem_w sums to 1 over the 4 slots)
    !   E1 extrap_nod3D (analytic patch, np 1): one missing node next to a triangle gets
    !      the mean over its elements' DISTINCT vertices (a triangle's repeated slot 4
    !      counted once); the mean with slot 4 counted twice differs (non-vacuity)
    ! On a global mesh (pi_w, rotated) a field linear in longitude jumps across the cyclic
    ! seam, so elements whose vertices span more than pi in longitude are skipped.
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth, g, density_0
    use mod_mesh,          only: t_mesh
    use mod_dyn,           only: t_dyn
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init
    use mod_param_phys,    only: theta
    use mod_mesh_analytic, only: generate_mixed_patch
    use mod_part_bounds,   only: owned_bounds
    use oce_pgf,           only: pressure_force_4_linfs_fullcell, pressure_force_4_zxxxx_shchepetkin
    use oce_ale,           only: compute_vel_nodes, update_vel
    use oce_pressure_bv,   only: compute_sigma_xy, smooth_nod
    use oce_fer_gm,        only: fer_gamma2vel
    use oce_initial_state, only: t_ic3d_config, extrap_nod3D
    implicit none

    type(t_partit) :: partit, apartit
    type(t_mesh)   :: mesh, amesh
    logical, allocatable :: use_el(:)     ! owned element off the cyclic seam

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.true., partition_np1=.false.)
    call check_kernels(mesh, partit, 'mesh')

    if (partit%npes == 1) then
        ! analytic mixed patch (1-rank identity partition of its own)
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=6, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_kernels(amesh, apartit, 'patch')
    end if

    call finish(partit, 'test_pgf_mixed')

contains

    subroutine check_kernels(m, p, tag)
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(inout) :: p
        character(len=*), intent(in)    :: tag
        type(t_dyn) :: dyn
        integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, nz, e, ntri, nquad
        real(kind=WP), allocatable :: x(:), y(:), hp(:,:), rho(:,:), pgx(:,:), pgy(:,:)
        real(kind=WP) :: pz, err(6), scal, dep, expx, rerr(3)
        real(kind=WP), allocatable :: tt(:,:), ss(:,:), swa(:,:), swb(:,:), sxy(:,:,:), gam(:,:,:), fuv(:,:,:)
        real(kind=WP), parameter :: GT = 2.0e-5_WP, GS = -1.0e-5_WP   ! R1 dT/dx, dS/dy
        real(kind=WP), parameter :: SWA0 = 2.0e-4_WP, SWB0 = 7.5e-4_WP, GG = 3.0e-3_WP
        logical :: allel
        integer :: k
        real(kind=WP), parameter :: A = 2.0_WP, B = -3.0_WP      ! P2 pressure gradient [Pa/m]
        real(kind=WP), parameter :: TX = 1.0e-6_WP, TY = -0.5e-6_WP  ! P3 level tilt [1/m]
        real(kind=WP), parameter :: RZ = 2.0e-3_WP                ! P3 d(rho)/dz [kg/m^4]
        real(kind=WP), parameter :: CX = 1.0e-6_WP                ! P4 d(rho)/dx [kg/m^4]
        real(kind=WP), parameter :: EX = 1.0e-6_WP, EY = 2.0e-6_WP   ! N2 d_eta gradient
        real(kind=WP), parameter :: DT = 100.0_WP

        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        if (p%npes == 1) then
            nElemF = m%elem2D
        else
            nElemF = p%myDim_elem2D + p%eDim_elem2D + p%eXDim_elem2D
        end if
        nl = m%nl
        ntri = count(m%elem2D_nnodes(1:nElemO) == 3); nquad = count(m%elem2D_nnodes(1:nElemO) == 4)
        if (allocated(use_el)) deallocate(use_el)
        allocate(use_el(nElemO))
        do e = 1, nElemO
            use_el(e) = maxval(m%coord_nod2D(1, m%elem2D_nodes(:, e))) &
                      - minval(m%coord_nod2D(1, m%elem2D_nodes(:, e))) < 3.14159_MP
        end do

        ! node coordinates in metres, relative to this rank's first node (a per-rank
        ! constant offset leaves every gradient unchanged)
        allocate(x(nNodL), y(nNodL))
        do n = 1, nNodL
            x(n) = real(m%coord_nod2D(1,n) - m%coord_nod2D(1,1), WP)*r_earth
            y(n) = real(m%coord_nod2D(2,n) - m%coord_nod2D(2,1), WP)*r_earth
        end do

        ! kernel arrays sized as the kernels declare them (global nod2D/elem2D)
        allocate(hp(nl, m%nod2D), rho(nl-1, m%nod2D), pgx(nl-1, m%elem2D), pgy(nl-1, m%elem2D))
        hp = 0.0_WP; rho = 0.0_WP

        !--- P1: resting ocean, pressure depends on the level only ---
        do n = 1, nNodL
            do nz = 1, nl
                hp(nz, n) = 1.0e4_WP*real(nz, WP)
            end do
        end do
        pgx = 0.0_WP; pgy = 0.0_WP
        call pressure_force_4_linfs_fullcell(hp, m, pgx, pgy, p)
        err(1) = owned_max(abs(pgx) + abs(pgy), m, p, nElemO)
        scal = 1.0e4_WP*real(nl, WP)/density_0/sqrt(minval(real(m%elem_area(1:nElemO), WP)))
        err(1) = err(1)/scal

        !--- P2: linear horizontal pressure ---
        do n = 1, nNodL
            do nz = 1, nl
                hp(nz, n) = 1.0e4_WP*real(nz, WP) + A*x(n) + B*y(n)
            end do
        end do
        pgx = 0.0_WP; pgy = 0.0_WP
        call pressure_force_4_linfs_fullcell(hp, m, pgx, pgy, p)
        err(2) = max(owned_max(abs(pgx - A/density_0), m, p, nElemO), &
                     owned_max(abs(pgy - B/density_0), m, p, nElemO)) &
                 /(abs(B)/density_0)

        !--- P3 / P4: Shchepetkin ---
        if (allocated(m%Z_3d_n))     deallocate(m%Z_3d_n)
        if (allocated(m%helem))      deallocate(m%helem)
        if (allocated(m%zbar_e_bot)) deallocate(m%zbar_e_bot)
        allocate(m%Z_3d_n(nl-1, nNodL), m%helem(nl-1, nElemF), m%zbar_e_bot(nElemF))
        do e = 1, nElemF
            do nz = 1, nl-1
                m%helem(nz, e) = m%zbar(nz) - m%zbar(nz+1)
            end do
            m%zbar_e_bot(e) = m%zbar(m%nlevels(e))
        end do
        ! P3: tilted levels, density linear in z
        do n = 1, nNodL
            do nz = 1, nl-1
                m%Z_3d_n(nz, n) = m%Z(nz)*(1.0_MP + real(TX*x(n) + TY*y(n), MP))
                rho(nz, n) = 0.5_WP + RZ*real(m%Z_3d_n(nz, n), WP)
            end do
        end do
        pgx = 0.0_WP; pgy = 0.0_WP
        call pressure_force_4_zxxxx_shchepetkin(rho, m, pgx, pgy, p)
        err(3) = owned_max(abs(pgx) + abs(pgy), m, p, nElemO)
        scal = g/density_0*RZ*real(maxval(abs(m%zbar)), WP)**2*(abs(TX) + abs(TY))
        err(3) = err(3)/scal
        ! P4: flat levels, density = CX*x
        do n = 1, nNodL
            do nz = 1, nl-1
                m%Z_3d_n(nz, n) = m%Z(nz)
                rho(nz, n) = CX*x(n)
            end do
        end do
        pgx = 0.0_WP; pgy = 0.0_WP
        call pressure_force_4_zxxxx_shchepetkin(rho, m, pgx, pgy, p)
        err(4) = 0.0_WP
        do e = 1, nElemO
            if (.not. use_el(e)) cycle
            dep = 0.0_WP
            do nz = m%ulevels(e), m%nlevels(e)-1
                expx = g/density_0*CX*(dep + 0.5_WP*real(m%helem(nz, e), WP))
                dep  = dep + real(m%helem(nz, e), WP)
                err(4) = max(err(4), abs(pgx(nz, e) - expx), abs(pgy(nz, e)))
            end do
        end do
        call max_all(err(4), p)
        err(4) = err(4)/(g/density_0*CX*real(maxval(abs(m%zbar)), WP))

        !--- N1: node gather of a uniform velocity ---
        allocate(dyn%uv(2, nl-1, nElemF), dyn%uvnode(2, nl-1, nNodL))
        dyn%uv(1,:,:) = 0.3_WP; dyn%uv(2,:,:) = -0.2_WP
        dyn%uvnode = 0.0_WP
        call compute_vel_nodes(dyn, m, p)
        err(5) = 0.0_WP
        do n = 1, nNodO
            do nz = m%ulevels_nod2D(n), m%nlevels_nod2D_min(n)-1
                err(5) = max(err(5), abs(dyn%uvnode(1,nz,n) - 0.3_WP), abs(dyn%uvnode(2,nz,n) + 0.2_WP))
            end do
        end do
        call max_all(err(5), p)
        err(5) = err(5)/0.3_WP

        !--- N2: update_vel with a linear d_eta ---
        allocate(dyn%uv_rhs(2, nl-1, nElemF), dyn%d_eta(nNodL))
        dyn%uv = 0.0_WP; dyn%uv_rhs = 0.0_WP
        do n = 1, nNodL
            dyn%d_eta(n) = EX*x(n) + EY*y(n)
        end do
        call update_vel(dyn, m, DT, p)
        err(6) = 0.0_WP
        do e = 1, nElemO
            if (.not. use_el(e)) cycle
            do nz = m%ulevels(e), m%nlevels(e)-1
                err(6) = max(err(6), abs(dyn%uv(1,nz,e) + g*theta*DT*EX), abs(dyn%uv(2,nz,e) + g*theta*DT*EY))
            end do
        end do
        call max_all(err(6), p)
        err(6) = err(6)/(g*theta*DT*EY)

        !--- R1: sigma_xy of linear T, S ---
        allocate(tt(nl-1, m%nod2D), ss(nl-1, m%nod2D), swa(nl-1, m%nod2D), swb(nl-1, m%nod2D), &
                 sxy(2, nl-1, m%nod2D))
        tt = 0.0_WP; ss = 0.0_WP; sxy = 0.0_WP
        swa = SWA0; swb = SWB0
        do n = 1, nNodL
            tt(:, n) = 10.0_WP + GT*x(n)
            ss(:, n) = 35.0_WP + GS*y(n)
        end do
        call compute_sigma_xy(tt, ss, swa, swb, m, sxy, p)
        rerr = 0.0_WP
        do n = 1, nNodO
            allel = .true.
            do k = 1, m%nod_in_elem2D_num(n)
                if (.not. use_el(m%nod_in_elem2D(k, n))) allel = .false.
            end do
            if (.not. allel) cycle
            do nz = m%ulevels_nod2D(n), m%nlevels_nod2D(n)-1
                rerr(1) = max(rerr(1), abs(sxy(1, nz, n) + SWA0*GT*density_0), &
                                       abs(sxy(2, nz, n) - SWB0*GS*density_0))
            end do
        end do
        call max_all(rerr(1), p)
        rerr(1) = rerr(1)/(SWB0*abs(GS)*density_0)
        !--- R2: smooth_nod keeps a constant ---
        tt = 3.7_WP
        call smooth_nod(tt, 3, m, p)
        do n = 1, nNodO
            do nz = m%ulevels_nod2D(n), min(nl-1, m%nlevels_nod2D(n))
                rerr(2) = max(rerr(2), abs(tt(nz, n) - 3.7_WP)/3.7_WP)
            end do
        end do
        call max_all(rerr(2), p)
        !--- R3: fer_gamma2vel of a horizontally uniform Gamma linear in zbar ---
        allocate(gam(2, nl, m%nod2D), fuv(2, nl-1, m%elem2D))
        gam = 0.0_WP; fuv = 0.0_WP
        do n = 1, nNodL
            gam(1, :, n) =  GG*real(m%zbar(1:nl), WP)
            gam(2, :, n) = -GG*real(m%zbar(1:nl), WP)
        end do
        call fer_gamma2vel(m, gam, fuv, p)
        do e = 1, nElemO
            do nz = m%ulevels(e), m%nlevels(e)-1
                rerr(3) = max(rerr(3), abs(fuv(1, nz, e) - GG), abs(fuv(2, nz, e) + GG))
            end do
        end do
        call max_all(rerr(3), p)
        rerr(3) = rerr(3)/GG
        if (p%mype == 0) write(*,'(a,a,3es10.2)') tag, ': R1 sigma_xy, R2 smooth_nod, R3 fer_gamma2vel rel. errors:', rerr
        call check(rerr(1) < 1.0e-10_WP, tag//': R1 sigma_xy of linear T, S exact', p%mype)
        call check(rerr(2) < 1.0e-13_WP, tag//': R2 smooth_nod keeps a constant', p%mype)
        call check(rerr(3) < 1.0e-12_WP, tag//': R3 fer_gamma2vel of a uniform Gamma exact', p%mype)
        if (tag == 'patch') call check_extrap(m, x, y)
        deallocate(tt, ss, swa, swb, sxy, gam, fuv)

        if (p%mype == 0) write(*,'(a,a,i0,a,i0,a,i0,a,6es10.2)') tag, ': owned triangles=', ntri, &
            ' quads=', nquad, ' (seam skipped: ', count(.not. use_el), ')  P1 P2 P3 P4 N1 N2 rel. errors:', err
        call check(err(1) < 1.0e-12_WP, tag//': P1 linfs PGF of a resting ocean == 0', p%mype)
        call check(err(2) < 1.0e-10_WP, tag//': P2 linfs PGF of a linear pressure exact', p%mype)
        call check(err(3) < 1.0e-10_WP, tag//': P3 Shchepetkin PGF == 0 on tilted levels (resting, stratified)', p%mype)
        call check(err(4) < 1.0e-10_WP, tag//': P4 Shchepetkin PGF of a horizontal density gradient exact', p%mype)
        call check(err(5) < 1.0e-12_WP, tag//': N1 node gather of a uniform velocity exact', p%mype)
        call check(err(6) < 1.0e-10_WP, tag//': N2 update_vel SSH gradient exact', p%mype)
        deallocate(x, y, hp, rho, pgx, pgy)
    end subroutine check_kernels

    subroutine check_extrap(m, x, y)
        ! E1 on the analytic patch (1 rank): the first node off the patch rim (2 rings)
        ! with a triangle around it is set missing at the surface; expected = mean of f
        ! over the distinct vertices (slot k skipped if it repeats an earlier slot) of
        ! its elements, the node itself excluded
        type(t_mesh),  intent(in) :: m
        real(kind=WP), intent(in) :: x(:), y(:)
        type(t_ic3d_config) :: ic
        real(kind=WP), allocatable :: arr(:,:)
        real(kind=WP) :: want, dup, xmin, xmax, ymin, ymax, hx, hy
        integer :: n0, n, k, j, e, v, cnt, cntd
        logical :: tri
        xmin = minval(x); xmax = maxval(x); ymin = minval(y); ymax = maxval(y)
        hx = 0.2_WP*(xmax - xmin); hy = 0.2_WP*(ymax - ymin)
        n0 = 0
        do n = 1, m%nod2D
            if (x(n) < xmin + hx .or. x(n) > xmax - hx .or. y(n) < ymin + hy .or. y(n) > ymax - hy) cycle
            tri = .false.
            do k = 1, m%nod_in_elem2D_num(n)
                if (m%elem2D_nnodes(m%nod_in_elem2D(k, n)) == 3) tri = .true.
            end do
            if (tri) then
                n0 = n; exit
            end if
        end do
        call check(n0 > 0, 'patch: E1 interior node next to a triangle found', 0)
        if (n0 == 0) return
        allocate(arr(m%nl-1, m%nod2D))
        do n = 1, m%nod2D
            arr(:, n) = 1.0_WP + 1.0e-4_WP*x(n) - 2.0e-4_WP*y(n)
        end do
        arr(1, n0) = ic%dummy
        want = 0.0_WP; dup = 0.0_WP; cnt = 0; cntd = 0
        do k = 1, m%nod_in_elem2D_num(n0)
            e = m%nod_in_elem2D(k, n0)
            do j = 1, 4
                v = m%elem2D_nodes(j, e)
                if (v == n0) cycle
                dup = dup + arr(1, v); cntd = cntd + 1
                if (any(m%elem2D_nodes(1:j-1, e) == v)) cycle
                want = want + arr(1, v); cnt = cnt + 1
            end do
        end do
        want = want/real(cnt, WP); dup = dup/real(cntd, WP)
        call extrap_nod3D(arr, ic, m)
        write(*,'(a,i0,a,3es12.4)') 'patch: E1 extrap_nod3D at node ', n0, ': got, distinct-vertex mean, ' // &
            'slot-4-twice mean:', arr(1, n0), want, dup
        call check(abs(arr(1, n0) - want) < 1.0e-13_WP*abs(want), &
                   'patch: E1 extrap_nod3D = mean over the distinct element vertices', 0)
        call check(abs(dup - want) > 1.0e-8_WP*abs(want), 'patch: E1 slot-4-twice mean differs (non-vacuity)', 0)
    end subroutine check_extrap

    real(kind=WP) function owned_max(f, m, p, nElemO)
        ! max over the owned elements and their wet levels of a (nl-1, elem) field
        real(kind=WP),  intent(in) :: f(:,:)
        type(t_mesh),   intent(in) :: m
        type(t_partit), intent(in) :: p
        integer,        intent(in) :: nElemO
        integer :: ee, kk
        owned_max = 0.0_WP
        do ee = 1, nElemO
            if (.not. use_el(ee)) cycle
            do kk = m%ulevels(ee), m%nlevels(ee)-1
                owned_max = max(owned_max, f(kk, ee))
            end do
        end do
        call max_all(owned_max, p)
    end function owned_max

end program test_pgf_mixed
