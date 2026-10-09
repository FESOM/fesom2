program fesom_advhordump_mr
    ! M2.12b MULTI-RANK tracer-advection byte-gate driver. Runs on dist_<NP> (NP>=1),
    ! prescribes the SAME analytic velocity/tracer fields as the 1-rank M1 gate (and the
    ! FESOM2 oracle src/fesom_advhor_dump.F90 multi-rank branch), drives the REAL
    ! advection assembly (advect_tracer = init_tracers_AB + do_oce_adv_tra, with all the
    ! multi-rank halo exchanges), and dumps the per-rank OWNED advection fields. The gate
    ! tools/run_advhor_gate_multirank.sh compares FESOM3 dist_NP vs FESOM2 dist_NP per
    ! rank on owned entries (1..myDim_*) for max|delta|=0 (the L8 same-partition rule).
    !
    ! Inputs are prescribed at owned+halo from the rotated coordinates (byte-identical to
    ! the oracle on the same partition) — node fields directly, the element velocity at
    ! owned elements then exchange_elem_full'd (halo elements have no local elem2D_nodes).
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): the FESOM2-oracle
    ! byte gate is retired (docs/HANDOFF.md "Gate migration"). The driver has no triangle-
    ! specific code of its own (it runs the converted kernels) and stays a diagnostic dump.
    !
    !   FESOM3_MESH_DIR    mesh dir (must contain dist_<NP>/)
    !   FESOM3_ADVHOR_OUT  out path prefix (per rank: <prefix>.<mype5>)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_halo,         only: exchange_elem_full
    use mod_tracer,       only: t_tracer
    use mod_dyn,          only: t_dyn
    use oce_muscl_adv,    only: muscl_adv_init
    use oce_tracer_mod,   only: init_tracers_AB
    use oce_adv_tra_driver, only: do_oce_adv_tra
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r2, wr_r3
    implicit none

    real(kind=WP), parameter :: dt   = 1800.0_WP
    real(kind=WP), parameter :: opth = 0.0_WP    ! MUSCL 4th-order fraction (FCT config)
    real(kind=WP), parameter :: optv = 1.0_WP    ! QR4C 4th-order fraction (FCT config)

    character(len=512) :: mesh_dir, out_prefix, out_path
    character(len=6)   :: rsuf
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_tracer) :: tr
    type(t_dyn)    :: dyn
    integer :: nsw, e, n, nz, nl, u, nzmin, nzmax
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: ttf(:,:), ttfAB(:,:), wvel(:,:)
    ! captured FCT-config fields (overwritten by the non-FCT advect call)
    real(kind=MP), allocatable :: valuesAB_s(:,:), fctLO_s(:,:)
    real(kind=MP), allocatable :: fmax_s(:,:), fmin_s(:,:), fplus_s(:,:), fminus_s(:,:)
    real(kind=MP), allocatable :: dh_s(:,:), dv_s(:,:), dt_s(:,:)
    real(kind=MP), allocatable :: dh_n(:,:), dv_n(:,:), dt_n(:,:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_ADVHOR_OUT', out_prefix)
    if (len_trim(out_prefix) == 0) out_prefix = 'advhor_mr_f3.bin'

    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl

    if (partit%npes == 1) then
        nNodO = mesh%nod2D;  nNodL = mesh%nod2D
        nEdgeO = mesh%edge2D; nElemO = mesh%elem2D; nElemF = mesh%elem2D
    else
        nNodO  = partit%myDim_nod2D
        nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
        nEdgeO = partit%myDim_edge2D
        nElemO = partit%myDim_elem2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if

    ! --- thickness / ALE depth arrays (initial full-cell linfs; pure functions of the
    !     halo-filled nlevels, so built directly over the full local range) -----------
    allocate(mesh%helem(nl-1, nElemF)); mesh%helem = 0.0_MP
    do e = 1, nElemF
        if (mesh%nlevels(e) <= 0) cycle
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    allocate(mesh%hnode(nl-1, nNodL), mesh%hnode_new(nl-1, nNodL)); mesh%hnode = 0.0_MP
    do n = 1, nNodL
        if (mesh%nlevels_nod2D(n) <= 0) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    mesh%hnode_new = mesh%hnode
    allocate(mesh%zbar_3d_n(nl, nNodL), mesh%Z_3d_n(nl-1, nNodL))
    mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
    do n = 1, nNodL
        if (mesh%nlevels_nod2D(n) <= 0) cycle
        nzmin = mesh%ulevels_nod2D(n); nzmax = mesh%nlevels_nod2D(n)
        zbar_srf = mesh%zbar(nzmin); zbar_bot = mesh%zbar(nzmax)
        mesh%zbar_3d_n(1:nzmin-1, n)       = mesh%zbar(1:nzmin-1)
        mesh%zbar_3d_n(nzmin, n)           = zbar_srf
        mesh%zbar_3d_n(nzmin+1:nzmax-1, n) = mesh%zbar(nzmin+1:nzmax-1)
        mesh%zbar_3d_n(nzmax, n)           = zbar_bot
        mesh%Z_3d_n(1:nzmin-1, n)          = mesh%Z(1:nzmin-1)
        mesh%Z_3d_n(nzmin, n)              = mesh%zbar_3d_n(nzmin,n)   + (mesh%zbar_3d_n(nzmin+1,n)-zbar_srf)/2
        mesh%Z_3d_n(nzmin+1:nzmax-2, n)    = mesh%Z(nzmin+1:nzmax-2)
        mesh%Z_3d_n(nzmax-1, n)            = mesh%zbar_3d_n(nzmax-1,n) + (zbar_bot-mesh%zbar_3d_n(nzmax-1,n))/2
    end do

    ! --- prescribe analytic fields (MUST equal the FESOM2 oracle) -------------------
    ! node fields at owned+halo (coords are byte-identical on the same partition);
    ! element velocity at owned elements then exchanged owner->halo.
    allocate(ttf(nl-1, nNodL), ttfAB(nl-1, nNodL), wvel(nl, nNodL))
    do n = 1, nNodL
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            ttf(nz, n)   = sin(2.0_WP*lon)*cos(lat) + 0.02_WP*real(nz, WP)
            ttfAB(nz, n) = 20.0_WP*sin(8.0_WP*lon)*cos(6.0_WP*lat) &
                         + 10.0_WP*sin(5.0_WP*lon)*cos(1.5_WP*real(nz, WP))
        end do
        do nz = 1, nl
            wvel(nz, n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz, WP))
        end do
    end do

    ! --- t_tracer / t_dyn (local sizes) --------------------------------------------
    tr%num_tracers = 1
    allocate(tr%data(1))
    allocate(tr%data(1)%values   (nl-1, nNodL))
    allocate(tr%data(1)%valuesAB (nl-1, nNodL))
    allocate(tr%data(1)%valuesold(2, nl-1, nNodL))
    allocate(tr%work%fct_LO          (nl-1, nNodL))
    allocate(tr%work%adv_flux_hor    (nl-1, nEdgeO))
    allocate(tr%work%adv_flux_ver    (nl,   nNodL))
    allocate(tr%work%fct_ttf_max     (nl-1, nNodL))
    allocate(tr%work%fct_ttf_min     (nl-1, nNodL))
    allocate(tr%work%fct_plus        (nl-1, nNodL))
    allocate(tr%work%fct_minus       (nl-1, nNodL))
    allocate(tr%work%del_ttf         (nl-1, nNodL))
    allocate(tr%work%del_ttf_advhoriz(nl-1, nNodL))
    allocate(tr%work%del_ttf_advvert (nl-1, nNodL))

    allocate(dyn%uv(2, nl-1, nElemF)); dyn%uv = 0.0_WP
    allocate(dyn%w(nl, nNodL), dyn%w_e(nl, nNodL), dyn%w_i(nl, nNodL))
    dyn%use_wsplit = .false.
    ! velocity at owned elements (analytic), then fill the halo by exchange
    do e = 1, nElemO
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1, e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1, e))
        do nz = 1, nl-1
            dyn%uv(1, nz, e) = 0.20_WP*cos(lat)*sin(lon) + 0.001_WP*real(nz, WP)
            dyn%uv(2, nz, e) = 0.15_WP*cos(2.0_WP*lat)
        end do
    end do
    if (partit%npes > 1) call exchange_elem_full(dyn%uv, partit)
    dyn%w = wvel; dyn%w_e = wvel; dyn%w_i = wvel

    call muscl_adv_init(tr%work, mesh, partit)

    ! --- FCT config (MUSCL/QR4C/FCT, AB2) -------------------------------------------
    tr%data(1)%values(:,:)      = ttf
    tr%data(1)%valuesold(1,:,:) = ttfAB
    tr%data(1)%AB_order    = 2
    tr%data(1)%tra_adv_hor = 'MUSCL'
    tr%data(1)%tra_adv_ver = 'QR4C'
    tr%data(1)%tra_adv_lim = 'FCT'
    tr%data(1)%tra_adv_ph  = opth
    tr%data(1)%tra_adv_pv  = optv
    call init_tracers_AB(1, tr, mesh, partit)
    allocate(valuesAB_s(nl-1, nNodO), fctLO_s(nl-1, nNodO))
    allocate(fmax_s(nl-1, nNodO), fmin_s(nl-1, nNodO), fplus_s(nl-1, nNodO), fminus_s(nl-1, nNodO))
    allocate(dh_s(nl-1, nNodO), dv_s(nl-1, nNodO), dt_s(nl-1, nNodO))
    valuesAB_s = real(tr%data(1)%valuesAB(1:nl-1, 1:nNodO),   MP)
    call do_oce_adv_tra(dt, dyn%uv, dyn%w, dyn%w_i, dyn%w_e, 1, dyn, tr, mesh, partit)
    do n = 1, nNodO
        tr%work%del_ttf(:,n) = tr%work%del_ttf(:,n) &
                             + tr%work%del_ttf_advhoriz(:,n) + tr%work%del_ttf_advvert(:,n)
    end do
    fctLO_s    = real(tr%work%fct_LO(1:nl-1, 1:nNodO),        MP)
    fmax_s     = real(tr%work%fct_ttf_max(1:nl-1, 1:nNodO),   MP)
    fmin_s     = real(tr%work%fct_ttf_min(1:nl-1, 1:nNodO),   MP)
    fplus_s    = real(tr%work%fct_plus(1:nl-1, 1:nNodO),      MP)
    fminus_s   = real(tr%work%fct_minus(1:nl-1, 1:nNodO),     MP)
    dh_s       = real(tr%work%del_ttf_advhoriz(1:nl-1, 1:nNodO), MP)
    dv_s       = real(tr%work%del_ttf_advvert(1:nl-1, 1:nNodO),  MP)
    dt_s       = real(tr%work%del_ttf(1:nl-1, 1:nNodO),          MP)

    ! --- non-FCT config (MUSCL/QR4C/NON, ph=pv=0.75) -------------------------------
    tr%data(1)%values(:,:)      = ttf
    tr%data(1)%valuesold(1,:,:) = ttfAB
    tr%data(1)%tra_adv_hor = 'MUSCL'
    tr%data(1)%tra_adv_ver = 'QR4C'
    tr%data(1)%tra_adv_lim = 'NON'
    tr%data(1)%tra_adv_ph  = 0.75_WP
    tr%data(1)%tra_adv_pv  = 0.75_WP
    call init_tracers_AB(1, tr, mesh, partit)
    call do_oce_adv_tra(dt, dyn%uv, dyn%w, dyn%w_i, dyn%w_e, 1, dyn, tr, mesh, partit)
    do n = 1, nNodO
        tr%work%del_ttf(:,n) = tr%work%del_ttf(:,n) &
                             + tr%work%del_ttf_advhoriz(:,n) + tr%work%del_ttf_advvert(:,n)
    end do
    allocate(dh_n(nl-1, nNodO), dv_n(nl-1, nNodO), dt_n(nl-1, nNodO))
    dh_n = real(tr%work%del_ttf_advhoriz(1:nl-1, 1:nNodO), MP)
    dv_n = real(tr%work%del_ttf_advvert(1:nl-1, 1:nNodO),  MP)
    dt_n = real(tr%work%del_ttf(1:nl-1, 1:nNodO),          MP)

    ! --- per-rank dump (OWNED entries) ---------------------------------------------
    out_path = trim(out_prefix)
    if (partit%npes /= 1) then
        write(rsuf,'(i5.5)') partit%mype
        out_path = trim(out_prefix)//'.'//rsuf
    end if
    call advhor_dump_open(u, trim(out_path), nNodO, nElemO, nEdgeO, nl)
    call wr_r2(u, 'valuesAB',             valuesAB_s)
    call wr_r2(u, 'fct_LO',               fctLO_s)
    call wr_r2(u, 'fct_ttf_max',          fmax_s)
    call wr_r2(u, 'fct_ttf_min',          fmin_s)
    call wr_r2(u, 'fct_plus',             fplus_s)
    call wr_r2(u, 'fct_minus',            fminus_s)
    call wr_r2(u, 'del_ttf_advhoriz_fct', dh_s)
    call wr_r2(u, 'del_ttf_advvert_fct',  dv_s)
    call wr_r2(u, 'del_ttf_fct',          dt_s)
    call wr_r2(u, 'del_ttf_advhoriz_non', dh_n)
    call wr_r2(u, 'del_ttf_advvert_non',  dv_n)
    call wr_r2(u, 'del_ttf_non',          dt_n)
    call advhor_dump_close(u)
    write(*,'(a,i0,a)') 'fesom_advhordump_mr: rank ', partit%mype, ' wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_advhordump_mr
