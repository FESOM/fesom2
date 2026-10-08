program test_wimpl_tra
    ! Implicit vertical advection of the FCT low-order solution with the implicit part
    ! w_i of the split vertical velocity (adv_tra_vert_impl, src/oce/oce_adv_tra_ver.F90,
    ! the port of FESOM2 oce_adv_tra_ver.F90:90-240) and its assembly in do_oce_adv_tra
    ! (src/oce/oce_adv_tra_driver.F90, FESOM2 oce_adv_tra_driver.F90:282-292).
    ! docs/plans/completed/2026-10-02-wsplit-smooth.md, Task 3.
    !
    ! WHY THIS TEST EXISTS
    ! --------------------
    ! With use_wsplit the vertical velocity is w = w_e + w_i (oce_wsplit). On the FCT
    ! tracer path the low-order solution is advanced explicitly (1st-order upwind) with
    ! w_e, then implicitly (backward-Euler upwind TDMA) with w_i, after which the low-order
    ! upwind flux is RECOMPUTED with the full w so that the antidiffusive flux is
    ! HO(w) - LO(w). There is no byte oracle (the FESOM2 gates are retired since
    ! bottom-at-vertices and the split itself is FESOM3's smooth function), so the kernel
    ! and the wiring are pinned by exact identities (L54).
    !
    ! THE KERNEL. Per owned column, v = dt*area(n)/areasvol(n) (f3: 1-D areas), h' =
    ! hnode_new, rows
    !   surface  a = 0             b = h' + w(nz)*v - min(0,w(nz+1))*v        c = -max(0,w(nz+1))*v
    !   interior a = min(0,w(nz))*v b = h' + max(0,w(nz))*v - min(0,w(nz+1))*v c = -max(0,w(nz+1))*v
    !   bottom   a = min(0,w(nz))*v b = h' + max(0,w(nz))*v                   c = 0
    ! rhs = -a*T(nz-1) - (b-h')*T(nz) - c*T(nz+1), T += solve, i.e. M*T^{n+1} = h'*T*.
    ! Column sums of M equal h' except the surface column (h' + w(nzmin)*v, the UNSIGNED
    ! surface flux, the implicit twin of adv_tra_ver_upw1's -w*T*area): the solve
    ! conserves h'*T up to the surface transport w(nzmin)*v*T^{n+1}(nzmin).
    !
    ! WHY CONSTANTS ARE PRESERVED (the consistency proof, C2). hnode_new is advanced with
    ! the FULL w: h' = h - dt*(w(nz) - w(nz+1))*area/areasvol. For uniform T the explicit
    ! upwind step with w_e gives T* = T*(h' + dt*dw_i*a)/h' (dw_i = w_i(nz) - w_i(nz+1),
    ! a = area/areasvol), and the implicit step solves (h' + dt*dw_i*a)*T^{n+1} = h'*T*,
    ! so T^{n+1} = T for ANY split with w_e + w_i = w. C2 therefore pins the split
    ! identity, the flux-form coefficients and that BOTH steps use hnode_new.
    !
    !   C1  identity: w_i == 0 -> ttf unchanged bitwise (even with hnode_new /= hnode)
    !   C2  constancy: uniform T, divergent w (0 at the surface and bottom faces), both
    !       signs, hnode_new consistent, splits (0.5,1), (0,1), (0.9,1):
    !       max|T^{n+1} - T| <= 1e-13*|T|; every cfl_z class (below Cu_min / bend /
    !       capped) populated, so all three branches of the split are exercised
    !   C2b constancy with w(surface) /= 0 of both signs (hnode_new consistent): the same
    !       bound (the unsigned surface row)
    !   C3  conservation: non-uniform T, w = 0 at the surface and bottom: per column
    !       sum h'*T^{n+1}*areasvol == sum h*T*areasvol to 1e-13 relative (+ the owned total)
    !   C3b with w(surface) /= 0: the column content changes by EXACTLY
    !       -dt*area*(w_e(1)*T(1) + w_i(1)*T^{n+1}(1)) -- the explicit part carries the old,
    !       the implicit part the new surface value (for uniform T that is the plan's
    !       -dt*w(1)*T*area, asserted in C2b)
    !   C4  boundedness at CFL_z ~ 10: uniform w through a column with one thin cell
    !       (0.08*h0, so CFL_z = |w|dt*(1/h0 + 1/(0.08 h0)) = 10.8 at its faces), surface
    !       open, bottom closed, both signs. Precondition asserted: per-cell explicit
    !       outflow <= hnode (the cap gives w_e*dt <= 1/(1/h_above + 1/h_below) < h_thin).
    !       T^{n+1} finite and within [min T, max T] of the column. Positive control: the
    !       UNCAPPED explicit step (w_e = w) leaves the range by more than half of it
    !       (the thin cell loses 10x its content; measured 0.95 / 5.6 of the range for
    !       w > 0 / w < 0 against ~1e-15 with the cap).
    !   C5  fully implicit limit, w at a single interior face k (both signs), hnode_new
    !       consistent: (a) the kernel alone == its two-row closed form, donor
    !       T_d*h'_d/(h'_d + |w|v), receiver T_r + |w|v*T_d^{n+1}/h'_r; (b) the sequence
    !       explicit(w_e = 0) -> implicit(w_i = w) == implicit upwind in flux form with the
    !       hnode_new mass: donor unchanged, receiver (h_r*T_r + |w|v*T_d)/(h_r + |w|v),
    !       and h'_d*T_d + h'_r*T_r == h_d*T_d + h_r*T_r; errors < 1e-12 relative
    !   C6  driver-level assembly: do_oce_adv_tra (FCT; UPW1 horizontal with vel = 0; QR4C
    !       vertical, num_ord = 1) with use_wsplit (0.5,1), a linear T(z) (horizontally
    !       uniform, both signs of the slope), uniform w0 of BOTH signs through the surface
    !       and the interior faces down to the second-last cell, the last cell inert (w = 0
    !       at both its faces), Cu = 2|w0|dt/h0 = 1.2. For w0 > 0 the FCT limiter is
    !       provably inactive (below); for w0 < 0 it is active, so the test's own assembly
    !       carries the limiter (oce_tra_adv_fct) and the wiring is checked with it live.
    !       Asserted:
    !       (i)   [w > 0] the limiter leaves the antidiffusive flux unchanged (the driver's
    !             clipped adv_flux_ver == the unclipped HO(w) - LO(w) at every owned face,
    !             every limiting factor b3 applies to a nonzero face is exactly 1, and the
    !             number of nonzero antidiffusive faces is > 0 -- the non-vacuity guard of
    !             the factor check, which an identically-zero flux would pass);
    !             [w < 0] the limiter IS active (some applied factor < 1), as derived
    !       (ii)  [both signs] the driver's fct_LO and del_ttf_advvert == the test's own
    !             assembly from the public parts: adv_tra_ver_upw1(w_e) -> the LO update ->
    !             adv_tra_vert_impl(w_i) -> exchange_nod(lo) -> HO(w) - LO(w) ->
    !             oce_tra_adv_fct -> oce_tra_adv_flux2dtracer(use_lo)
    !       (iii) positive control: the WRONG assembly HO(w) - LO(w_e). [w > 0, unlimited]
    !             differs from the driver's by exactly the double-counted w_i transport,
    !             (-w_i(nz)*T(nz) + w_i(nz+1)*T(nz+1))*area*dt/areasvol, far above the
    !             tolerance (the wiring bug a conservation gate cannot see: it telescopes);
    !             [w < 0, limited like the driver] differs by more than the tolerance (the
    !             limiter does not mask the wiring bug)
    !       (iv)  [w > 0] the closed form of the whole step. s = T(nz) - T(nz+1), r_i(nz) =
    !             w_i(nz)*dt/h0 (the split's own w_i), cells nzmin..N (N = nzmax-1):
    !             LO: d_{N-1} = 0 (drained end cell: split invariance), d_N = 0 (inert),
    !             d_n = ((Cu/2)*s + r_i(n+1)*d_{n+1})/(1 + r_i(n)) upward, lo = T - d;
    !             HO - LO increments: +Cu*s/4 at the surface cell, -(Cu/4)*s/(1 - Cu/2) at
    !             cell N-1 (h' = h0*(1 - Cu/2) there), 0 elsewhere (uniform w and linear T:
    !             the centred and the 4th-order face values are T(z_face), the antidiffusive
    !             flux -w0*s*area/2 is uniform over the interior faces)
    !
    ! WHY THE LIMITER IS INACTIVE ONLY FOR UPWARD FLOW, OPEN SURFACE, INERT BOTTOM CELL.
    ! At the drained end of the flow (no inflow face) the low-order cell keeps its value
    ! (lo = T, the split invariance) while the HO flux through its outflow face is the
    ! centred (T_above + T)/2 instead of the upwind T: the HO wants to move that cell PAST
    ! its own value, which is the cluster extremum of a monotone profile -- unless an inert
    ! cell below it supplies the room (its value T - s; the HO move is (Cu/4)*s/(1 - Cu/2)
    ! < s for Cu < 4/3). For downward flow the drained end is the surface cell and nothing
    ! is above it (and with an open surface the unsigned surface flux carries T(1) both
    ! ways, so lo(1) = T(1) and the HO increase Cu*s/4 is clipped). Hence the closed form
    ! and the inactive-limiter premise hold for w > 0 only; the plan's "uniform interior w"
    ! with w = 0 at the bottom face is clipped at its last cell, which is why the last cell
    ! is made inert here. (i) asserts the premise for w > 0 and its negation for w < 0.
    !
    ! SYNTHESISED COLUMNS (synth_columns; by the GLOBAL node index so np 1/2 agree). pi has
    ! ulevels_nod2D = 1 everywhere and at least 4 layers, so the nzmin-based row layout of
    ! adv_tra_vert_impl and of compute_CFLz / compute_Wvel_split would otherwise be tested
    ! at nzmin = 1 only, and the 2-layer column the row layout admits as its minimum (empty
    ! interior loop; surface row nzmin, bottom row nzmin+1) never. Every 10th global node
    ! with >= 14 levels gets ulevels_nod2D = 3 (cavity), every 97th with >= 5 levels gets
    ! nlevels_nod2D = ulevels + 2 (2 layers; the selector constants are named at
    ! synth_columns); C1-C5 run on the mix (C5 skips columns with fewer than 3 layers: its
    ! single interior face is nzmin+2). The originals are restored before C6: the driver's
    ! FCT limiter builds its clusters from ulevels(elem), which the node-only synthesis
    ! would leave inconsistent. Columns with fewer than 2 layers are refused at SETUP by
    ! read_mesh (assert_min_layers, ctest test_minlayers), not by the kernel, which is
    ! stateless and validates nothing.
    use mpi
    use mod_precision,      only: WP, MP
    use mod_constants,      only: pi
    use mod_mesh,           only: t_mesh
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex, set_partition
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_part_bounds,    only: owned_bounds, is_multirank
    use mod_halo,           only: exchange_nod
    use oce_ale,            only: compute_CFLz, compute_Wvel_split
    use oce_adv_tra_ver,    only: adv_tra_ver_upw1, adv_tra_ver_qr4c, adv_tra_vert_impl
    use oce_adv_tra_flux,   only: oce_tra_adv_flux2dtracer
    use oce_adv_tra_fct,    only: oce_tra_adv_fct
    use oce_adv_tra_driver, only: do_oce_adv_tra
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    type(t_tracer) :: tracers
    integer :: nfail, nsw
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl
    integer, allocatable :: ulev0(:), nlev0(:)      ! the mesh's own ulevels/nlevels (restored for C6)
    real(kind=WP), parameter :: dt   = 1800.0_WP     ! s
    real(kind=WP), parameter :: h0   = 10.0_WP       ! reference layer thickness [m]
    real(kind=WP), parameter :: t0   = 10.0_WP
    real(kind=WP), parameter :: thin = 0.08_WP       ! C4: thin cell = thin*h0
    real(kind=WP), parameter :: ctl_margin = 1.0e3_WP ! C6 (iii): the wrong assembly must differ by > ctl_margin*tol
                                                     ! (measured 29.9/58.8 K m for w > 0, 3.9/1.4 for w < 0, tol 1e-11)
    integer,       parameter :: nsplit = 3
    real(kind=WP), parameter :: smin(nsplit) = [0.5_WP, 0.0_WP, 0.9_WP]
    real(kind=WP), parameter :: smax(nsplit) = [1.0_WP, 1.0_WP, 1.0_WP]
    real(kind=WP), allocatable :: tin(:,:), tlo(:,:), tnew(:,:), tother(:,:)
    real(kind=WP), allocatable :: flux_v(:,:), flux_h(:,:), dttf_h(:,:), dttf_v(:,:)
    real(kind=WP), allocatable :: vel(:,:,:)

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)

    !=========================================================================
    ! pi mesh (the scaffold of test_vinv / test_wsplit)
    !=========================================================================
    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.true.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    nElemF = mesh%elem2D
    if (partit%npes > 1) nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    nl = mesh%nl
    allocate(dyn%w(nl, nNodL), dyn%w_e(nl, nNodL), dyn%w_i(nl, nNodL), dyn%cfl_z(nl, nNodL))
    dyn%w = 0.0_WP; dyn%w_e = 0.0_WP; dyn%w_i = 0.0_WP; dyn%cfl_z = 0.0_WP
    allocate(mesh%hnode(nl-1, nNodL), mesh%hnode_new(nl-1, nNodL))
    allocate(mesh%zbar_3d_n(nl, nNodL), mesh%Z_3d_n(nl-1, nNodL))
    allocate(tin(nl-1, nNodL), tlo(nl-1, nNodL), tnew(nl-1, nNodL), tother(nl-1, nNodL))
    allocate(flux_v(nl, nNodL), flux_h(nl-1, nEdgeO), dttf_h(nl-1, nNodL), dttf_v(nl-1, nNodL))
    flux_v = 0.0_WP; flux_h = 0.0_WP; dttf_h = 0.0_WP; dttf_v = 0.0_WP

    call synth_columns()

    call part_c1()
    call part_c2_c3()
    call part_c4()
    call part_c5()
    mesh%ulevels_nod2D(1:nNodL) = ulev0          ! C6: the driver's FCT limiter needs the
    mesh%nlevels_nod2D(1:nNodL) = nlev0          ! mesh's own (element-consistent) levels
    call part_c6()

    if (partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_wimpl_tra: OK'
        else
            write(*,'(a,i0,a)') 'test_wimpl_tra: ', nfail, ' FAILURE(S)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)   ! NO third arg: abort is presence-based
    if (nfail /= 0) error stop 1

contains
    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            if (partit%mype == 0) write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check_true

    ! --- MPI reductions (identity at np 1) ------------------------------------------
    function gmax(x) result(m)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: m
        integer :: ierr
        m = x
        if (partit%npes > 1) call MPI_Allreduce(x, m, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gmax

    function gmin(x) result(m)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: m
        integer :: ierr
        m = x
        if (partit%npes > 1) call MPI_Allreduce(x, m, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gmin

    function gsumr(x) result(s)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: s
        integer :: ierr
        s = x
        if (partit%npes > 1) call MPI_Allreduce(x, s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsumr

    function gsum(i) result(s)
        integer, intent(in) :: i
        integer :: s, ierr
        s = i
        if (partit%npes > 1) call MPI_Allreduce(i, s, 1, MPI_INTEGER, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsum

    function gall(l) result(a)
        logical, intent(in) :: l
        logical :: a
        integer :: ierr
        a = l
        if (partit%npes > 1) call MPI_Allreduce(l, a, 1, MPI_LOGICAL, MPI_LAND, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gall

    !=========================================================================
    ! Synthesised columns (header)
    !=========================================================================
    subroutine synth_columns()
        ! Selectors (by the global node index g): the cavity set is test_wsplit's
        ! synth_cavity set (same constants), so the two tests exercise the same columns.
        integer, parameter :: cav_every = 10       ! every 10th global node ...
        integer, parameter :: cav_nlev_min = 14    ! ... with >= 14 levels keeps >= 11 layers below the cavity surface
        integer, parameter :: cav_ulev = 3         ! the cavity surface level (nzmin = 3)
        integer, parameter :: two_every = 97       ! every 97th global node (not a cavity one) ...
        integer, parameter :: two_nlev_min = 5     ! ... with >= 4 layers (pi's minimum), so the cut only ever shrinks
        integer :: n, g, ncav, n2
        allocate(ulev0(nNodL), nlev0(nNodL))
        ulev0 = mesh%ulevels_nod2D(1:nNodL)
        nlev0 = mesh%nlevels_nod2D(1:nNodL)
        ncav = 0; n2 = 0
        do n = 1, nNodL
            g = n
            if (partit%npes > 1) g = partit%myList_nod2D(n)
            if (mod(g, cav_every) == 0 .and. mesh%nlevels_nod2D(n) >= cav_nlev_min) then
                mesh%ulevels_nod2D(n) = cav_ulev
                if (n <= nNodO) ncav = ncav + 1
            else if (mod(g, two_every) == 0 .and. mesh%nlevels_nod2D(n) >= two_nlev_min) then
                mesh%nlevels_nod2D(n) = mesh%ulevels_nod2D(n) + 2
                if (n <= nNodO) n2 = n2 + 1
            end if
        end do
        ncav = gsum(ncav); n2 = gsum(n2)
        if (partit%mype == 0) write(*,'(a,i0,a,i0)') '  synthesised owned columns: cavity (ulevels = 3) ', ncav, &
                                                     '  2-layer ', n2
        call check_true('cavity columns synthesised (ulevels_nod2D = 3)', ncav > 0)
        call check_true('2-layer columns synthesised (nlevels_nod2D = ulevels + 2)', n2 > 0)
    end subroutine synth_columns

    !=========================================================================
    ! State builders
    !=========================================================================
    subroutine build_layers(layer_kind)
        ! Prescribed layer thickness hnode at owned+halo nodes, with zbar_3d_n/Z_3d_n
        ! consistent with it (QR4C reads them); hnode_new = hnode until hnew_from_w.
        !   layer_kind 1: non-uniform h = h0*(1 + 0.25*cos(nz)) in [0.75, 1.25]*h0
        !   layer_kind 2: uniform h0
        !   layer_kind 3: uniform h0 with one thin cell thin*h0 at nzmin+2 (C4)
        integer, intent(in) :: layer_kind
        integer :: n, nz, nzmin, nzmax
        real(kind=WP) :: h, zb
        mesh%hnode = 0.0_MP; mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
            zb = 0.0_WP
            mesh%zbar_3d_n(nzmin, n) = real(zb, MP)
            do nz = nzmin, nzmax - 1
                select case (layer_kind)
                case (1)
                    h = h0*(1.0_WP + 0.25_WP*cos(real(nz, WP)))
                case (2)
                    h = h0
                case default
                    h = h0
                    if (nz == nzmin + 2) h = thin*h0
                end select
                mesh%hnode(nz, n)  = real(h, MP)
                mesh%Z_3d_n(nz, n) = real(zb - 0.5_WP*h, MP)
                zb = zb - h
                mesh%zbar_3d_n(nz+1, n) = real(zb, MP)
            end do
        end do
        mesh%hnode_new = mesh%hnode
    end subroutine build_layers

    function hnew_from_w() result(rmin)
        ! hnode_new = hnode - dt*(w(nz) - w(nz+1))*area/areasvol: the thickness the model
        ! advances with the FULL w (a pure column, no horizontal transport). Returns the
        ! global min of hnode_new/hnode over owned+halo cells (the precondition > 0).
        real(kind=WP) :: rmin, a
        integer :: n, nz
        rmin = huge(1.0_WP)
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            a = real(mesh%area(n), WP)/real(mesh%areasvol(n), WP)
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                mesh%hnode_new(nz, n) = real(real(mesh%hnode(nz, n), WP) &
                                             - dt*(dyn%w(nz, n) - dyn%w(nz+1, n))*a, MP)
                rmin = min(rmin, real(mesh%hnode_new(nz, n), WP)/real(mesh%hnode(nz, n), WP))
            end do
        end do
        rmin = gmin(rmin)
    end function hnew_from_w

    subroutine w_divergent(wsign, surf_frac)
        ! w = A*sin(pi*(nz - nzmin)/N) at the faces nzmin..nzmax (0 at the surface and the
        ! bottom face), A = wsign*(h0/(2dt))*min(3, 0.4*N): CFL_z up to ~4 on deep columns
        ! while every cell keeps hnode_new > 0 (|dw|*dt <= h0*min(3,0.4N)*sin(pi/(2N)) <
        ! 0.62*h0 < 0.75*h0). surf_frac /= 0 opens the surface face: w(nzmin) = surf_frac*A.
        real(kind=WP), intent(in) :: wsign, surf_frac
        integer :: n, nz, nzmin, nzmax, nlay
        real(kind=WP) :: amp
        dyn%w = 0.0_WP
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
            nlay  = nzmax - nzmin
            amp   = wsign*(h0/(2.0_WP*dt))*min(3.0_WP, 0.4_WP*real(nlay, WP))
            do nz = nzmin, nzmax
                dyn%w(nz, n) = amp*sin(pi*real(nz - nzmin, WP)/real(nlay, WP))
            end do
            dyn%w(nzmax, n) = 0.0_WP        ! exactly (sin(pi) is not)
            dyn%w(nzmin, n) = surf_frac*amp
        end do
    end subroutine w_divergent

    subroutine split(cmin, cmax, n_lo, n_mid, n_hi)
        ! The model's split: cfl_z from hnode_new (compute_CFLz), then compute_Wvel_split.
        ! Counts the owned faces below Cu_min / on the bend / capped.
        real(kind=WP), intent(in)  :: cmin, cmax
        integer,       intent(out) :: n_lo, n_mid, n_hi
        integer :: n, nz
        real(kind=WP) :: ccut
        dyn%use_wsplit    = .true.
        dyn%wsplit_mincfl = cmin
        dyn%wsplit_maxcfl = cmax
        call compute_CFLz(dyn, mesh, dt, partit)
        call compute_Wvel_split(dyn, mesh, partit)
        ccut = 2.0_WP*cmax - cmin
        n_lo = 0; n_mid = 0; n_hi = 0
        do n = 1, nNodO
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)
                if (dyn%cfl_z(nz, n) <= cmin) then
                    n_lo = n_lo + 1
                else if (dyn%cfl_z(nz, n) >= ccut) then
                    n_hi = n_hi + 1
                else
                    n_mid = n_mid + 1
                end if
            end do
        end do
        n_lo = gsum(n_lo); n_mid = gsum(n_mid); n_hi = gsum(n_hi)
    end subroutine split

    subroutine lo_step(we, t_in, t_lo)
        ! The driver's explicit low-order vertical step (oce_adv_tra_driver.F90, FCT branch)
        ! REPLICATED: the upwind flux with we, then
        !     lo = (T*hnode + (flux(nz) - flux(nz+1))*dt/areasvol)/hnode_new
        ! (the horizontal low-order flux is zero: a pure column). Halo cells copy T.
        real(kind=WP), intent(in)  :: we(nl, nNodL), t_in(nl-1, nNodL)
        real(kind=WP), intent(out) :: t_lo(nl-1, nNodL)
        integer :: n, nz
        call adv_tra_ver_upw1(we, t_in, mesh, flux_v, o_init_zero=.true., partit=partit)
        t_lo = t_in
        do n = 1, nNodO
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                t_lo(nz, n) = (t_in(nz, n)*mesh%hnode(nz, n) &
                               + (flux_v(nz, n) - flux_v(nz+1, n))*dt/mesh%areasvol(n)) &
                              /mesh%hnode_new(nz, n)
            end do
        end do
    end subroutine lo_step

    subroutine profile_b(t)
        ! a non-uniform, non-monotone, strictly positive profile (varies per column too)
        real(kind=WP), intent(out) :: t(nl-1, nNodL)
        integer :: n, nz
        t = t0
        do n = 1, nNodL
            if (mesh%nlevels_nod2D(n) <= 0) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                t(nz, n) = t0 + 2.0_WP*sin(0.7_WP*real(nz, WP)) + 0.5_WP*cos(0.3_WP*real(n, WP) + real(nz, WP))
            end do
        end do
    end subroutine profile_b

    !=========================================================================
    ! C1 - identity
    !=========================================================================
    subroutine part_c1()
        real(kind=WP) :: rmin
        logical :: ok
        call build_layers(1)
        call w_divergent(1.0_WP, 0.3_WP)
        rmin = hnew_from_w()                  ! hnode_new /= hnode, w_i = 0
        call check_true('C1 precondition hnode_new > 0', rmin > 0.0_WP)
        call profile_b(tin)
        tnew  = tin
        dyn%w_i = 0.0_WP
        call adv_tra_vert_impl(dt, dyn%w_i, tnew, mesh, partit)
        ok = all(tnew(:, 1:nNodO) == tin(:, 1:nNodO))
        call check_true('C1 w_i == 0: ttf unchanged bitwise (hnode_new /= hnode)', gall(ok))
    end subroutine part_c1

    !=========================================================================
    ! C2 / C2b / C3 / C3b - constancy and conservation
    !=========================================================================
    subroutine part_c2_c3()
        integer :: isg, iopen, iset, n, nz, nzmin, nzmax, n_lo, n_mid, n_hi
        real(kind=WP) :: wsign, surf_frac, rmin, err2, err3, errb, c_old, c_new, c_abs, pred, cflmax
        real(kind=WP) :: tot_old, tot_new, tot_pred, tot_abs
        character(len=8)  :: stag
        character(len=24) :: tag

        call build_layers(1)
        do isg = 1, 2
            wsign = merge(1.0_WP, -1.0_WP, isg == 1)
            stag  = merge('(w > 0)', '(w < 0)', isg == 1)
            do iopen = 0, 1
                surf_frac = merge(0.0_WP, 0.3_WP, iopen == 0)
                call w_divergent(wsign, surf_frac)
                rmin = hnew_from_w()
                call check_true('C2 precondition hnode_new >= 0.1*hnode '//stag, rmin >= 0.1_WP)
                do iset = 1, nsplit
                    write(tag,'(a,f4.2,a,f4.2,a)') '(', smin(iset), ',', smax(iset), ')'
                    call split(smin(iset), smax(iset), n_lo, n_mid, n_hi)
                    cflmax = gmax(maxval(dyn%cfl_z(:, 1:nNodO)))

                    ! ---- C2 / C2b: uniform T stays uniform -------------------------------
                    tin = t0
                    call lo_step(dyn%w_e, tin, tlo)
                    tnew = tlo
                    call adv_tra_vert_impl(dt, dyn%w_i, tnew, mesh, partit)
                    err2 = 0.0_WP; errb = 0.0_WP
                    do n = 1, nNodO
                        if (mesh%nlevels_nod2D(n) <= 0) cycle
                        nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                        err2 = max(err2, maxval(abs(tnew(nzmin:nzmax-1, n) - t0))/t0)
                        if (iopen == 1) then
                            ! the plan's literal C3b form for uniform T: -dt*w(1)*T*area
                            c_old = sum(mesh%hnode(nzmin:nzmax-1, n)*t0)*mesh%areasvol(n)
                            c_new = sum(mesh%hnode_new(nzmin:nzmax-1, n)*tnew(nzmin:nzmax-1, n))*mesh%areasvol(n)
                            pred  = -dt*dyn%w(nzmin, n)*t0*mesh%area(n)
                            errb  = max(errb, abs(c_new - c_old - pred)/c_old)
                        end if
                    end do
                    err2 = gmax(err2); errb = gmax(errb)
                    if (partit%mype == 0) then
                        if (iopen == 0) then
                            write(*,'(5a,es10.2,a,f6.2,a,i0,a,i0,a,i0)') '  C2  ', stag, ' ', tag, &
                                ': max|T - T0|/T0 = ', err2, '  max cfl_z = ', cflmax, &
                                '  faces below/bend/capped = ', n_lo, '/', n_mid, '/', n_hi
                        else
                            write(*,'(5a,es10.2,a,es10.2)') '  C2b ', stag, ' ', tag, &
                                ': max|T - T0|/T0 = ', err2, '  uniform-T surface budget err = ', errb
                        end if
                    end if
                    if (iopen == 0) then
                        call check_true('C2 constancy, closed surface '//stag//' '//tag, err2 <= 1.0e-13_WP)
                        call check_true('C2 every cfl_z class populated '//stag//' '//tag, &
                                        n_lo > 0 .and. n_mid > 0 .and. n_hi > 0)
                    else
                        call check_true('C2b constancy, open surface '//stag//' '//tag, err2 <= 1.0e-13_WP)
                        call check_true('C2b uniform T: content change == -dt*w(1)*T*area '//stag//' '//tag, &
                                        errb <= 1.0e-13_WP)
                    end if

                    ! ---- C3 / C3b: conservation with a non-uniform T ----------------------
                    call profile_b(tin)
                    call lo_step(dyn%w_e, tin, tlo)
                    tnew = tlo
                    call adv_tra_vert_impl(dt, dyn%w_i, tnew, mesh, partit)
                    err3 = 0.0_WP; tot_old = 0.0_WP; tot_new = 0.0_WP; tot_pred = 0.0_WP; tot_abs = 0.0_WP
                    do n = 1, nNodO
                        if (mesh%nlevels_nod2D(n) <= 0) cycle
                        nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                        c_old = sum(mesh%hnode(nzmin:nzmax-1, n)*tin(nzmin:nzmax-1, n))*mesh%areasvol(n)
                        c_abs = sum(mesh%hnode(nzmin:nzmax-1, n)*abs(tin(nzmin:nzmax-1, n)))*mesh%areasvol(n)
                        c_new = sum(mesh%hnode_new(nzmin:nzmax-1, n)*tnew(nzmin:nzmax-1, n))*mesh%areasvol(n)
                        ! surface transport: explicit part with the old, implicit with the new value
                        pred  = -dt*mesh%area(n)*(dyn%w_e(nzmin, n)*tin(nzmin, n) + dyn%w_i(nzmin, n)*tnew(nzmin, n))
                        err3  = max(err3, abs(c_new - c_old - pred)/c_abs)
                        tot_old = tot_old + c_old; tot_new = tot_new + c_new
                        tot_pred = tot_pred + pred; tot_abs = tot_abs + c_abs
                    end do
                    err3 = gmax(err3)
                    tot_old = gsumr(tot_old); tot_new = gsumr(tot_new); tot_pred = gsumr(tot_pred); tot_abs = gsumr(tot_abs)
                    if (partit%mype == 0) then
                        if (iopen == 0) then
                            write(*,'(5a,es10.2,a,es10.2)') '  C3  ', stag, ' ', tag, &
                                ': max column budget err = ', err3, '  owned total err = ', abs(tot_new - tot_old)/tot_abs
                        else
                            write(*,'(5a,es10.2,a,es10.2)') '  C3b ', stag, ' ', tag, &
                                ': max column budget err = ', err3, '  owned total err = ', &
                                abs(tot_new - tot_old - tot_pred)/tot_abs
                        end if
                    end if
                    if (iopen == 0) then
                        call check_true('C3 conservation, closed surface '//stag//' '//tag, err3 <= 1.0e-13_WP)
                    else
                        call check_true('C3b content change == -dt*area*(w_e*T_old + w_i*T_new)(surface) '//stag//' '//tag, &
                                        err3 <= 1.0e-13_WP)
                    end if
                    call check_true('C3 owned total budget '//stag//' '//tag, &
                                    abs(tot_new - tot_old - tot_pred)/tot_abs <= 1.0e-13_WP)
                end do
            end do
        end do
    end subroutine part_c2_c3

    !=========================================================================
    ! C4 - boundedness at CFL_z ~ 10
    !=========================================================================
    subroutine part_c4()
        integer :: isg, n, nz, nzmin, nzmax, n_lo, n_mid, n_hi
        real(kind=WP) :: wsign, w0, rmin, cflmax, outmax, viol, viol_ctl, tmn, tmx, a
        logical :: finite
        character(len=8) :: stag

        call build_layers(3)
        do isg = 1, 2
            wsign = merge(1.0_WP, -1.0_WP, isg == 1)
            stag  = merge('(w > 0)', '(w < 0)', isg == 1)
            w0 = wsign*0.8_WP*h0/dt
            dyn%w = 0.0_WP
            do n = 1, nNodL
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                dyn%w(nzmin:nzmax-1, n) = w0          ! surface open, bottom face closed
            end do
            rmin = hnew_from_w()
            call check_true('C4 precondition hnode_new > 0 '//stag, rmin > 0.0_WP)
            call split(0.5_WP, 1.0_WP, n_lo, n_mid, n_hi)
            cflmax = gmax(maxval(dyn%cfl_z(:, 1:nNodO)))

            ! precondition: per-cell explicit outflow <= hnode (the cap)
            outmax = 0.0_WP
            do n = 1, nNodO
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                a = real(mesh%area(n), WP)/real(mesh%areasvol(n), WP)
                do nz = nzmin, nzmax - 1
                    outmax = max(outmax, (max(0.0_WP, dyn%w_e(nz, n)) + max(0.0_WP, -dyn%w_e(nz+1, n)))*dt*a &
                                         /real(mesh%hnode(nz, n), WP))
                end do
            end do
            outmax = gmax(outmax)

            ! the step on a non-monotone profile
            tin = t0
            do n = 1, nNodL
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                    tin(nz, n) = t0 + 3.0_WP*sin(1.3_WP*real(nz, WP))
                end do
            end do
            call lo_step(dyn%w_e, tin, tlo)
            tnew = tlo
            call adv_tra_vert_impl(dt, dyn%w_i, tnew, mesh, partit)
            ! positive control: the uncapped explicit step (w_e = w, w_i = 0)
            call lo_step(dyn%w, tin, tother)
            viol = 0.0_WP; viol_ctl = 0.0_WP; finite = .true.
            do n = 1, nNodO
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                tmn = minval(tin(nzmin:nzmax-1, n)); tmx = maxval(tin(nzmin:nzmax-1, n))
                do nz = nzmin, nzmax - 1
                    if (.not. abs(tnew(nz, n)) < huge(1.0_WP)) finite = .false.
                    viol     = max(viol,     (tmn - tnew(nz, n))/(tmx - tmn),   (tnew(nz, n) - tmx)/(tmx - tmn))
                    viol_ctl = max(viol_ctl, (tmn - tother(nz, n))/(tmx - tmn), (tother(nz, n) - tmx)/(tmx - tmn))
                end do
            end do
            viol = gmax(viol); viol_ctl = gmax(viol_ctl)
            if (partit%mype == 0) then
                write(*,'(a,a,a,f6.2,a,f6.3,a,i0,a,i0,a,i0)') '  C4  ', stag, ': max cfl_z = ', cflmax, &
                    '  max explicit outflow/hnode = ', outmax, '  faces below/bend/capped = ', n_lo, '/', n_mid, '/', n_hi
                write(*,'(a,a,a,es10.2,a,f8.2)') '  C4  ', stag, ': overshoot/range = ', viol, &
                    '   uncapped control overshoot/range = ', viol_ctl
            end if
            call check_true('C4 precondition: explicit outflow <= hnode in every cell '//stag, outmax <= 1.0_WP)
            call check_true('C4 the split reaches the capped branch '//stag, n_hi > 0)
            call check_true('C4 T^{n+1} finite '//stag, gall(finite))
            call check_true('C4 T^{n+1} within [min T, max T] '//stag, viol <= 1.0e-12_WP)
            call check_true('C4 positive control: the uncapped explicit step overshoots by > 0.5 of the range '//stag, &
                            viol_ctl > 0.5_WP)
        end do
    end subroutine part_c4

    !=========================================================================
    ! C5 - the fully implicit limit, one face
    !=========================================================================
    subroutine part_c5()
        integer :: isg, n, nz, k, d, r, nzmin, nzmax
        real(kind=WP) :: wsign, w0, rmin, v, wa, err_a, err_b, err_m, s_d, s_r, tmax
        character(len=8) :: stag

        call build_layers(2)
        do isg = 1, 2
            wsign = merge(1.0_WP, -1.0_WP, isg == 1)
            stag  = merge('(w > 0)', '(w < 0)', isg == 1)
            w0 = wsign*0.4_WP*h0/dt
            dyn%w = 0.0_WP
            do n = 1, nNodL
                if (mesh%nlevels_nod2D(n) - mesh%ulevels_nod2D(n) < 3) cycle   ! face nzmin+2 must be interior
                dyn%w(mesh%ulevels_nod2D(n) + 2, n) = w0
            end do
            rmin = hnew_from_w()
            call check_true('C5 precondition hnode_new > 0 '//stag, rmin > 0.0_WP)
            call profile_b(tin)
            tmax = gmax(maxval(abs(tin(:, 1:nNodO))))

            ! (a) the kernel alone on T
            tnew = tin
            call adv_tra_vert_impl(dt, dyn%w, tnew, mesh, partit)
            ! (b) the sequence: explicit with w_e = 0 (T*hnode/hnode_new), implicit with w_i = w
            dyn%w_e = 0.0_WP
            call lo_step(dyn%w_e, tin, tlo)
            tother = tlo
            call adv_tra_vert_impl(dt, dyn%w, tother, mesh, partit)

            err_a = 0.0_WP; err_b = 0.0_WP; err_m = 0.0_WP
            do n = 1, nNodO
                if (mesh%nlevels_nod2D(n) - mesh%ulevels_nod2D(n) < 3) cycle   ! (w = 0 there: untouched)
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                k  = nzmin + 2
                v  = dt*real(mesh%area(n), WP)/real(mesh%areasvol(n), WP)
                wa = abs(w0)
                if (w0 > 0.0_WP) then
                    d = k; r = k - 1          ! upward: the donor is the cell below the face
                else
                    d = k - 1; r = k
                end if
                ! (a) two-row closed form of the kernel
                s_d = tin(d, n)*mesh%hnode_new(d, n)/(mesh%hnode_new(d, n) + wa*v)
                s_r = tin(r, n) + wa*v*s_d/mesh%hnode_new(r, n)
                err_a = max(err_a, abs(tnew(d, n) - s_d), abs(tnew(r, n) - s_r))
                ! (b) implicit upwind in flux form with the hnode_new mass
                s_d = tin(d, n)
                s_r = (mesh%hnode(r, n)*tin(r, n) + wa*v*tin(d, n))/(mesh%hnode(r, n) + wa*v)
                err_b = max(err_b, abs(tother(d, n) - s_d), abs(tother(r, n) - s_r))
                ! every other cell untouched
                do nz = nzmin, nzmax - 1
                    if (nz == d .or. nz == r) cycle
                    err_a = max(err_a, abs(tnew(nz, n)   - tin(nz, n)))
                    err_b = max(err_b, abs(tother(nz, n) - tin(nz, n)))
                end do
                err_m = max(err_m, abs(mesh%hnode_new(d, n)*tother(d, n) + mesh%hnode_new(r, n)*tother(r, n) &
                                       - mesh%hnode(d, n)*tin(d, n) - mesh%hnode(r, n)*tin(r, n)) &
                                   /(mesh%hnode(d, n)*abs(tin(d, n)) + mesh%hnode(r, n)*abs(tin(r, n))))
            end do
            err_a = gmax(err_a)/tmax; err_b = gmax(err_b)/tmax; err_m = gmax(err_m)
            if (partit%mype == 0) write(*,'(a,a,a,es10.2,a,es10.2,a,es10.2)') '  C5  ', stag, &
                ': kernel vs two-row closed form = ', err_a, '  sequence vs flux form = ', err_b, &
                '  two-cell mass err = ', err_m
            call check_true('C5a kernel == two-row closed form (1e-12) '//stag, err_a <= 1.0e-12_WP)
            call check_true('C5b sequence == implicit upwind flux form with the hnode_new mass (1e-12) '//stag, &
                            err_b <= 1.0e-12_WP)
            call check_true('C5b two-cell mass conserved (1e-12) '//stag, err_m <= 1.0e-12_WP)
        end do
    end subroutine part_c5

    !=========================================================================
    ! C6 - driver-level assembly
    !=========================================================================
    subroutine part_c6()
        integer :: isg, iw, n, nz, nzmin, nzmax, n_lo, n_mid, n_hi, nnz
        real(kind=WP) :: bz, w0, wsign, cu, rmin, s, ri, rin, tol, err_lo, err_v, err_w, err_cf, dmax, pmax
        real(kind=WP) :: d(nl), inc(nl), tdrv, tcf, hmin_fp, hmin_fm
        logical :: ok_lim, upward
        character(len=8) :: stag, wtag
        real(kind=WP), allocatable :: adf_test(:,:), adf_unlim(:,:), dttf_wrong(:,:)
        real(kind=WP), allocatable :: fmin_t(:,:), fmax_t(:,:), fplus_t(:,:), fminus_t(:,:), adfh_t(:,:)

        allocate(adf_test(nl, nNodL), adf_unlim(nl, nNodL), dttf_wrong(nl-1, nNodL))
        allocate(fmin_t(nl-1, nNodL), fmax_t(nl-1, nNodL), fplus_t(nl-1, nNodL), fminus_t(nl-1, nNodL))
        allocate(adfh_t(nl-1, nEdgeO))
        ! the one-tracer FCT state the driver reads (shapes as fesom_conserve)
        tracers%num_tracers = 1
        allocate(tracers%data(1))
        allocate(tracers%data(1)%values(nl-1, nNodL), tracers%data(1)%valuesAB(nl-1, nNodL))
        tracers%data(1)%ID          = 1
        tracers%data(1)%tra_adv_hor = 'UPW1'
        tracers%data(1)%tra_adv_ver = 'QR4C'
        tracers%data(1)%tra_adv_lim = 'FCT'
        tracers%data(1)%tra_adv_ph  = 1.0_WP
        tracers%data(1)%tra_adv_pv  = 1.0_WP     ! num_ord = 1: pure 4th-order centred
        allocate(tracers%work%fct_LO(nl-1, nNodL), tracers%work%adv_flux_hor(nl-1, nEdgeO))
        allocate(tracers%work%adv_flux_ver(nl, nNodL))
        allocate(tracers%work%fct_ttf_max(nl-1, nNodL), tracers%work%fct_ttf_min(nl-1, nNodL))
        allocate(tracers%work%fct_plus(nl-1, nNodL), tracers%work%fct_minus(nl-1, nNodL))
        allocate(tracers%work%del_ttf_advhoriz(nl-1, nNodL), tracers%work%del_ttf_advvert(nl-1, nNodL))
        ! MUSCL gradients off: no up/downwind triangles make every on-the-fly up/dn
        ! increment 0 (the C6 closed forms assume the 2*dT-only reconstruction)
        allocate(tracers%work%edge_up_dn_tri(2, nEdgeO))
        allocate(tracers%work%tr_xy(2, nl-1, nElemF))
        tracers%work%edge_up_dn_tri = 0
        tracers%work%tr_xy = 0.0_WP
        allocate(vel(2, nl-1, nElemF)); vel = 0.0_WP
        allocate(mesh%helem(nl-1, nElemF)); mesh%helem = real(h0, MP)   ! read by adv_tra_hor_upw1 (x vel = 0)

        call build_layers(2)
        cu = 1.2_WP
        do iw = 1, 2
            upward = (iw == 1)
            wsign  = merge(1.0_WP, -1.0_WP, upward)
            wtag   = merge('(w > 0)', '(w < 0)', upward)
            w0 = wsign*cu*h0/(2.0_WP*dt)
            dyn%w = 0.0_WP
            do n = 1, nNodL
                if (mesh%nlevels_nod2D(n) <= 0) cycle
                nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                dyn%w(nzmin:nzmax-2, n) = w0        ! open surface ... cell nzmax-2 is the flow's end; cell nzmax-1 inert
            end do
            rmin = hnew_from_w()
            call check_true('C6 precondition hnode_new > 0 '//wtag, rmin > 0.0_WP)
            call split(0.5_WP, 1.0_WP, n_lo, n_mid, n_hi)
            call check_true('C6 the split is active (faces on the bend) '//wtag, n_mid > 0)

            do isg = 1, 2
                bz   = merge(0.02_WP, -0.02_WP, isg == 1)      ! K/m; s = T(nz) - T(nz+1) = bz*h0
                stag = merge('(s > 0)', '(s < 0)', isg == 1)
                s    = bz*h0
                do n = 1, nNodL
                    tin(:, n) = t0
                    if (mesh%nlevels_nod2D(n) <= 0) cycle
                    do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                        tin(nz, n) = t0 + bz*real(mesh%Z_3d_n(nz, n), WP)
                    end do
                end do
                tol = 1.0e-13_WP*t0*h0

                ! ---- the driver ----------------------------------------------------------------
                tracers%data(1)%values   = tin
                tracers%data(1)%valuesAB = tin
                tracers%work%del_ttf_advhoriz = 0.0_MP
                tracers%work%del_ttf_advvert  = 0.0_MP
                call do_oce_adv_tra(dt, vel, dyn%w, dyn%w_i, dyn%w_e, 1, dyn, tracers, mesh, partit)

                ! ---- the test's own assembly from the public parts, limiter included ---------
                call lo_step(dyn%w_e, tin, tlo)                                 ! explicit LO with w_e
                call adv_tra_vert_impl(dt, dyn%w_i, tlo, mesh, partit)          ! implicit with w_i
                if (is_multirank(partit)) call exchange_nod(tlo, partit)        ! the limiter's a1 reads lo at the halo
                call adv_tra_ver_upw1(dyn%w, tin, mesh, adf_test, o_init_zero=.true., partit=partit)  ! LO(w)
                call adv_tra_ver_qr4c(dyn%w, tin, mesh, tracers%data(1)%tra_adv_pv, adf_test, &
                                      o_init_zero=.false., partit=partit)                              ! HO(w) - LO(w)
                adf_unlim = adf_test                                            ! the unclipped flux, for (i)
                adfh_t = 0.0_WP; fmin_t = 0.0_WP; fmax_t = 0.0_WP; fplus_t = 0.0_WP; fminus_t = 0.0_WP
                call oce_tra_adv_fct(dt, tin, tlo, adfh_t, adf_test, fmin_t, fmax_t, fplus_t, fminus_t, &
                                     mesh, partit=partit)
                dttf_h = 0.0_WP; dttf_v = 0.0_WP; flux_h = 0.0_WP
                call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_v, flux_h, adf_test, mesh, &
                                              use_lo=.true., ttf=tin, lo=tlo, partit=partit)
                ! ---- the WRONG assembly: HO(w) - LO(w_e); unlimited for w > 0 (exact identity),
                !      limited like the driver for w < 0 (the limiter must not mask the bug) ----
                call adv_tra_ver_upw1(dyn%w_e, tin, mesh, flux_v, o_init_zero=.true., partit=partit)
                call adv_tra_ver_qr4c(dyn%w, tin, mesh, tracers%data(1)%tra_adv_pv, flux_v, &
                                      o_init_zero=.false., partit=partit)
                if (.not. upward) then
                    adfh_t = 0.0_WP; fmin_t = 0.0_WP; fmax_t = 0.0_WP; fplus_t = 0.0_WP; fminus_t = 0.0_WP
                    call oce_tra_adv_fct(dt, tin, tlo, adfh_t, flux_v, fmin_t, fmax_t, fplus_t, fminus_t, &
                                         mesh, partit=partit)
                end if
                dttf_h = 0.0_WP; dttf_wrong = 0.0_WP; flux_h = 0.0_WP
                call oce_tra_adv_flux2dtracer(dt, dttf_h, dttf_wrong, flux_h, flux_v, mesh, &
                                              use_lo=.true., ttf=tin, lo=tlo, partit=partit)

                ok_lim = .true.; err_lo = 0.0_WP; err_v = 0.0_WP; err_w = 0.0_WP; err_cf = 0.0_WP
                dmax = 0.0_WP; pmax = 0.0_WP; hmin_fp = 1.0_WP; hmin_fm = 1.0_WP; nnz = 0
                do n = 1, nNodO
                    if (mesh%nlevels_nod2D(n) <= 0) cycle
                    nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
                    ! (i) the limiter's effect: the driver's clipped flux against the unclipped
                    ! one, and the factors b3 APPLIES to nonzero faces (b3: a face flux >= 0 is
                    ! clipped by fct_minus(nz-1) and fct_plus(nz), < 0 by fct_plus(nz-1) and
                    ! fct_minus(nz); the surface face by its own cell's factor only). Unused
                    ! factors can be 0 -- e.g. the surface cell's one-sided factor on the side
                    ! with no antidiffusive contribution, whose a3 cluster is the surface level
                    ! alone -- and must not be counted. nnz counts the nonzero faces so the
                    ! factor check cannot pass on an identically-zero flux.
                    do nz = nzmin, nzmax
                        if (tracers%work%adv_flux_ver(nz, n) /= adf_unlim(nz, n)) ok_lim = .false.
                    end do
                    do nz = nzmin, nzmax - 1
                        if (adf_unlim(nz, n) == 0.0_WP) cycle
                        nnz = nnz + 1
                        if (adf_unlim(nz, n) >= 0.0_WP) then
                            hmin_fp = min(hmin_fp, real(tracers%work%fct_plus(nz, n), WP))
                            if (nz > nzmin) hmin_fm = min(hmin_fm, real(tracers%work%fct_minus(nz-1, n), WP))
                        else
                            hmin_fm = min(hmin_fm, real(tracers%work%fct_minus(nz, n), WP))
                            if (nz > nzmin) hmin_fp = min(hmin_fp, real(tracers%work%fct_plus(nz-1, n), WP))
                        end if
                    end do
                    ! (ii) driver == test assembly
                    err_lo = max(err_lo, maxval(abs(real(tracers%work%fct_LO(nzmin:nzmax-1, n), WP) - tlo(nzmin:nzmax-1, n))))
                    err_v  = max(err_v,  maxval(abs(real(tracers%work%del_ttf_advvert(nzmin:nzmax-1, n), WP) - dttf_v(nzmin:nzmax-1, n))))
                    ! (iii) wrong - driver (== the double-counted w_i transport for w > 0)
                    do nz = nzmin, nzmax - 1
                        tcf = -dyn%w_i(nz, n)*tin(nz, n)
                        if (nz + 1 <= nzmax - 1) tcf = tcf + dyn%w_i(nz+1, n)*tin(nz+1, n)
                        tcf  = tcf*mesh%area(n)*dt/mesh%areasvol(n)
                        tdrv = dttf_wrong(nz, n) - real(tracers%work%del_ttf_advvert(nz, n), WP)
                        dmax  = max(dmax, abs(tdrv))
                        pmax  = max(pmax, abs(tcf))
                        err_w = max(err_w, abs(tdrv - tcf))
                    end do
                    ! (iv) closed form of the whole step (w > 0): T_new = T - d + inc
                    d = 0.0_WP; inc = 0.0_WP
                    do nz = nzmax - 3, nzmin, -1
                        ri  = dyn%w_i(nz, n)*dt/h0
                        rin = dyn%w_i(nz+1, n)*dt/h0
                        d(nz) = (0.5_WP*cu*s + rin*d(nz+1))/(1.0_WP + ri)
                    end do
                    inc(nzmin)   = 0.25_WP*cu*s
                    inc(nzmax-2) = inc(nzmax-2) - 0.25_WP*cu*s/(1.0_WP - 0.5_WP*cu)
                    do nz = nzmin, nzmax - 1
                        ! the ALE reconstruct of the driver's tendency: T + (dttf_v + T*(h - h'))/h'
                        tdrv = tin(nz, n) + (real(tracers%work%del_ttf_advvert(nz, n), WP) &
                                             + tin(nz, n)*(mesh%hnode(nz, n) - mesh%hnode_new(nz, n)))/mesh%hnode_new(nz, n)
                        tcf  = tin(nz, n) - d(nz) + inc(nz)
                        err_cf = max(err_cf, abs(tdrv - tcf))
                    end do
                end do
                err_lo = gmax(err_lo); err_v = gmax(err_v); err_w = gmax(err_w); err_cf = gmax(err_cf)
                dmax = gmax(dmax); pmax = gmax(pmax); hmin_fp = gmin(hmin_fp); hmin_fm = gmin(hmin_fm); nnz = gsum(nnz)
                if (partit%mype == 0) then
                    write(*,'(a,a,a,a,a,f6.3,a,f6.3,a,i0,a)') '  C6  ', wtag, ' ', stag, &
                        ': min limiting factor b3 applies to a nonzero face: fct_plus = ', &
                        hmin_fp, '  fct_minus = ', hmin_fm, '  (', nnz, ' nonzero faces)'
                    write(*,'(a,a,a,a,a,es10.2,a,es10.2,a,es9.2,a)') '  C6  ', wtag, ' ', stag, ': |fct_LO - test LO| = ', err_lo, &
                        '  |del_ttf_advvert - test| = ', err_v, '  (tol ', tol, ' K m)'
                    write(*,'(a,a,a,a,a,es10.2,a,es10.2,a,es10.2)') '  C6  ', wtag, ' ', stag, ': wrong assembly |delta| = ', dmax, &
                        '  predicted = ', pmax, '  |delta - predicted| = ', err_w
                    if (upward) write(*,'(a,a,a,a,a,es10.2,a)') '  C6  ', wtag, ' ', stag, &
                        ': |T_new(driver) - closed form| = ', err_cf, ' K'
                end if
                if (upward) then
                    call check_true('C6 (i) the limiter leaves the antidiffusive flux unchanged '//wtag//' '//stag, gall(ok_lim))
                    call check_true('C6 (i) every limiting factor applied to a nonzero face == 1 '//wtag//' '//stag, &
                                    hmin_fp == 1.0_WP .and. hmin_fm == 1.0_WP)
                    call check_true('C6 (i) nonzero antidiffusive faces > 0 (non-vacuity of the factor check) '//wtag//' '//stag, &
                                    nnz > 0)
                else
                    call check_true('C6 (i) the limiter is active for downward flow, as derived '//wtag//' '//stag, &
                                    nnz > 0 .and. (hmin_fp < 1.0_WP .or. hmin_fm < 1.0_WP))
                end if
                call check_true('C6 (ii) driver fct_LO == test assembly LO '//wtag//' '//stag, err_lo <= tol)
                call check_true('C6 (ii) driver del_ttf_advvert == test assembly (limiter included) '//wtag//' '//stag, &
                                err_v <= tol)
                call check_true('C6 (iii) positive control: HO(w) - LO(w_e) differs by more than the tolerance '//wtag//' '//stag, &
                                dmax > ctl_margin*tol)
                if (upward) then
                    call check_true('C6 (iii) positive control: the difference == the double-counted w_i transport '//wtag//' '//stag, &
                                    err_w <= 1.0e-12_WP*pmax)
                    call check_true('C6 (iv) driver result == closed form of the FCT+wsplit step (1e-12) '//wtag//' '//stag, &
                                    err_cf <= 1.0e-12_WP*t0)
                end if
            end do
        end do
    end subroutine part_c6
end program test_wimpl_tra
