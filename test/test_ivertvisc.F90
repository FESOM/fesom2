program test_ivertvisc
    ! Implicit vertical viscosity + advection solver (impl_vert_visc_ale, src/oce/
    ! oce_dyn_ivertvisc.F90): the FORM of the implicit vertical advection must match the
    ! momentum advection scheme.
    !
    ! WHY THIS TEST EXISTS
    ! --------------------
    ! With use_wsplit the vertical velocity is split, w = w_e + w_i: only w_e is advected
    ! explicitly, w_i goes into the tridiagonal solve. FESOM2 -- and the scalar scheme
    ! momadv_opt==2, which is flux form in both directions -- puts it there in upwind FLUX
    ! form d(w_i u)/dz, whose row sums are (wu-wd)*zinv: a uniform u gets the tendency
    ! -u*dw_i/dz, balanced in the flux-form scheme by the horizontal u*div(u) through
    ! continuity. The vector-invariant scheme (momadv_opt==1) is ADVECTIVE form in both
    ! directions (zeta x u + grad KE; Block C does w_e du/dz as d(w_e u)/dz - u dw_e/dz),
    ! so the flux-form operator would leave a spurious u*dw_i/dz wherever the split is
    ! active. The solver therefore subtracts (wu-wd)*zinv from the diagonal when
    ! momadv_opt==1: upwind ADVECTIVE form, zero row sums, constants preserved implicitly.
    !
    ! WHAT IS ASSERTED. Av = 0, stress = 0, C_d = 0 and UV_rhs(in) = 0, so the output
    ! UV_rhs = u^{n+1} - u^n is the implicit advection increment alone.
    !   I1 momadv_opt=1, uniform u, divergent w_i (surface AND bottom faces nonzero)
    !        -> increment == 0 to round-off
    !   I2 momadv_opt=2, same input -> increment /= 0: the flux-form signature, i.e. the
    !        proof that I1 has teeth
    !   I3 momadv_opt=1, one face w0 > 0 at nz0, u = u0 + D in cell nz0, u0 elsewhere:
    !        d(nz0-1) = D*r1/(1+r1),  r1 = dt*w0/h(nz0-1);  all other cells 0
    !   I4 momadv_opt=2, same input: closed form of the ORIGINAL flux form -- the
    !        regression pin that the ==2 path is untouched:
    !        d(nz0)   = -(u0+D)*r2/(1+r2),  r2 = dt*w0/h(nz0)
    !        d(nz0-1) =  (u0+D)*r1/(1+r2);  all other cells 0
    !   I5 finiteness
    ! The closed forms hold per element column because w_i is set to w0 at EVERY node of
    ! face nz0 (the kite-weighted element mean sum(elem_w*w_i) == w0 to 1 ulp on triangles
    ! AND quads -- a 3-node /3 mean over the 4 slots of a quad would give 4/3*w0) and helem
    ! is what the solver rebuilds zinv from. Both velocity components are checked (v = -u/2
    ! on input, so d_v = -d_u/2).
    !
    ! MESHES (hybrid meshes, Task 10): FESOM3_MESH_DIR (default pi, rotated; with
    ! FESOM3_TEST_ROTATE=0 read unrotated, e.g. the NEMO-mixed cut-out) or, with
    ! FESOM3_TEST_PATCH=mixed at np 1, the analytic mixed triangle/quad patch.
    use mpi
    use mod_precision,     only: WP, MP
    use mod_mesh,          only: t_mesh
    use mod_dyn,           only: t_dyn
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_param_phys,    only: C_d
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_mesh_analytic, only: generate_mixed_patch
    use mod_part_bounds,   only: owned_bounds
    use oce_dyn_ivertvisc, only: impl_vert_visc_ale
    implicit none

    character(len=512) :: mesh_dir, cbuf
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    integer :: nfail, nsw
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, nz, e, nzmin, nzmax, nchk
    real(kind=WP), allocatable :: Av(:,:), stress(:,:), expct(:)
    real(kind=WP), parameter   :: dt  = 1800.0_WP
    real(kind=WP), parameter   :: u0  = 0.3_WP      ! uniform part
    real(kind=WP), parameter   :: dlt = 0.05_WP     ! perturbation in cell nz0 (I3/I4)
    real(kind=WP), parameter   :: w0  = 2.0e-3_WP   ! single-face w_i (I3/I4): r ~ 0.25
    integer,       parameter   :: nz0 = 3           ! the face; cells nz0-1 and nz0 see it
    real(kind=WP) :: zmax, zmax2, h1, h2, r1, r2, err

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)
    cbuf = ''
    call get_environment_variable('FESOM3_TEST_PATCH', cbuf)
    if (trim(cbuf) == 'mixed') then
        if (partit%npes /= 1) error stop 'test_ivertvisc: FESOM3_TEST_PATCH needs np 1'
        call generate_mixed_patch(mesh, partit, nx=13, ny=11, nl=6, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
    else
        if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
        cbuf = ''
        call get_environment_variable('FESOM3_TEST_ROTATE', cbuf)
        if (trim(cbuf) == '0') then
            call read_mesh(mesh, partit, trim(mesh_dir), 0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP, &
                           force_rotation=.false., n_cw_swaps=nsw)
        else
            call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                           force_rotation=.true., n_cw_swaps=nsw)
        end if
        call compute_geometry(mesh, partit, cartesian=.true.)
    end if
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    write(*,'(a,i0,a,i0)') '  owned triangles=', count(mesh%elem2D_nnodes(1:nElemO) == 3), &
                           ' quads=', count(mesh%elem2D_nnodes(1:nElemO) == 4)
    if (partit%npes == 1) then
        nElemF = mesh%elem2D
    else
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    nl = mesh%nl

    allocate(dyn%uv(2, nl-1, nElemF), dyn%uv_rhs(2, nl-1, nElemF), dyn%w_i(nl, nNodL))
    allocate(mesh%helem(nl-1, nElemF), mesh%zbar_e_bot(nElemF))
    allocate(Av(nl, nElemF), stress(2, nElemF), expct(nl-1))
    Av = 0.0_WP; stress = 0.0_WP; C_d = 0.0_WP          ! isolate the advection operator
    do nz = 1, nl-1
        mesh%helem(nz,:) = 10.0_MP + 2.0_MP*real(nz, MP)  ! any positive ALE thickness
    end do
    do e = 1, nElemF
        mesh%zbar_e_bot(e) = -sum(mesh%helem(mesh%ulevels(e):mesh%nlevels(e)-1, e))
    end do

    !=========================================================================
    ! I1 / I2 - uniform u, w_i with vertical divergence everywhere (incl. the surface
    ! face, which the surface row carries with its full sign, and the bottom face)
    !=========================================================================
    do n = 1, nNodL
        do nz = 1, nl
            dyn%w_i(nz,n) = 1.0e-3_WP*sin(0.7_WP*real(nz,WP))*(1.0_WP + 0.05_WP*real(mod(n,5),WP))
        end do
    end do
    dyn%uv(1,:,:) = u0
    dyn%uv(2,:,:) = -0.5_WP*u0

    dyn%momadv_opt = 1
    dyn%uv_rhs = 0.0_WP
    call impl_vert_visc_ale(dyn, mesh, dt, Av, stress, partit)
    zmax = owned_absmax(dyn%uv_rhs)
    write(*,'(a,es12.4)') '  I1 vinv, uniform u, divergent w_i: max|du| = ', zmax
    call check_true('I1 momadv_opt=1: uniform u unchanged by implicit w_i advection', zmax < 1.0e-13_WP)

    dyn%momadv_opt = 2
    dyn%uv_rhs = 0.0_WP
    call impl_vert_visc_ale(dyn, mesh, dt, Av, stress, partit)
    zmax2 = owned_absmax(dyn%uv_rhs)
    write(*,'(a,es12.4)') '  I2 scalar, same input           : max|du| = ', zmax2
    call check_true('I2 momadv_opt=2: flux form moves a uniform u (I1 has teeth)', zmax2 > 1.0e-6_WP)

    !=========================================================================
    ! I3 / I4 - single face w0 at nz0, closed forms
    !=========================================================================
    dyn%w_i = 0.0_WP
    dyn%w_i(nz0, 1:nNodL) = w0
    dyn%uv(1,:,:)   = u0
    dyn%uv(1,nz0,:) = u0 + dlt
    dyn%uv(2,:,:)   = -0.5_WP*dyn%uv(1,:,:)

    dyn%momadv_opt = 1
    dyn%uv_rhs = 0.0_WP
    call impl_vert_visc_ale(dyn, mesh, dt, Av, stress, partit)
    err = 0.0_WP; nchk = 0
    do e = 1, nElemO
        nzmin = mesh%ulevels(e); nzmax = mesh%nlevels(e)
        if (nzmin > nz0-2 .or. nzmax-1 < nz0) cycle      ! both cells must be regular rows
        nchk = nchk + 1
        h1 = real(mesh%helem(nz0-1,e), WP); r1 = dt*w0/h1
        expct = 0.0_WP
        expct(nz0-1) = dlt*r1/(1.0_WP + r1)
        err = max(err, column_err(e, nzmin, nzmax))
    end do
    err = gmax(err); nchk = gsum(nchk)
    write(*,'(a,es12.4,a,i0,a)') '  I3 vinv  single-face closed form  : max err = ', err, '  (', nchk, ' columns)'
    call check_true('I3 momadv_opt=1: upwind ADVECTIVE closed form', err < 1.0e-12_WP .and. nchk > 0)

    dyn%momadv_opt = 2
    dyn%uv_rhs = 0.0_WP
    call impl_vert_visc_ale(dyn, mesh, dt, Av, stress, partit)
    err = 0.0_WP
    do e = 1, nElemO
        nzmin = mesh%ulevels(e); nzmax = mesh%nlevels(e)
        if (nzmin > nz0-2 .or. nzmax-1 < nz0) cycle
        h1 = real(mesh%helem(nz0-1,e), WP); r1 = dt*w0/h1
        h2 = real(mesh%helem(nz0,  e), WP); r2 = dt*w0/h2
        expct = 0.0_WP
        expct(nz0)   = -(u0 + dlt)*r2/(1.0_WP + r2)
        expct(nz0-1) =  (u0 + dlt)*r1/(1.0_WP + r2)
        err = max(err, column_err(e, nzmin, nzmax))
    end do
    err = gmax(err)
    write(*,'(a,es12.4)') '  I4 scalar single-face closed form : max err = ', err
    call check_true('I4 momadv_opt=2: upwind FLUX closed form (path untouched)', err < 1.0e-12_WP)

    !=========================================================================
    ! I5 - finiteness of the last solve
    !=========================================================================
    call check_true('I5 finite', all(abs(dyn%uv_rhs(:,:,1:nElemO)) < huge(1.0_WP)))

    if (nfail == 0) then
        write(*,'(a)') 'test_ivertvisc: OK'
    else
        write(*,'(a,i0,a)') 'test_ivertvisc: ', nfail, ' FAILURE(S)'
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

    ! max over the column of |uv_rhs - expected| for both components (v = -u/2)
    function column_err(e, nzmin, nzmax) result(cerr)
        integer, intent(in) :: e, nzmin, nzmax
        real(kind=WP) :: cerr
        integer :: k
        cerr = 0.0_WP
        do k = nzmin, nzmax-1
            cerr = max(cerr, abs(dyn%uv_rhs(1,k,e) - expct(k)), &
                             abs(dyn%uv_rhs(2,k,e) + 0.5_WP*expct(k)))
        end do
    end function column_err

    function owned_absmax(x) result(m)
        real(kind=WP), intent(in) :: x(:,:,:)
        real(kind=WP) :: m
        integer :: k, ee
        m = 0.0_WP
        do ee = 1, nElemO
            do k = mesh%ulevels(ee), mesh%nlevels(ee)-1
                m = max(m, abs(x(1,k,ee)), abs(x(2,k,ee)))
            end do
        end do
        m = gmax(m)
    end function owned_absmax

    function gmax(x) result(m)
        real(kind=WP), intent(in) :: x
        real(kind=WP) :: m
        integer :: ierr
        m = x
        if (partit%npes > 1) call MPI_Allreduce(x, m, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gmax

    function gsum(i) result(s)
        integer, intent(in) :: i
        integer :: s, ierr
        s = i
        if (partit%npes > 1) call MPI_Allreduce(i, s, 1, MPI_INTEGER, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsum
end program test_ivertvisc
