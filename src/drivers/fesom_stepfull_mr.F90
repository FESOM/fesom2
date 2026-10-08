program fesom_stepfull_mr
    ! M2.12c-3 MULTI-RANK WHOLE-STEP byte-gate driver. Runs on dist_<NP> (NP>=1) and drives
    ! the ASSEMBLED ocean step (mod_step_oce::step_oce) ONCE through the optional partit, so
    ! EVERY kernel runs its multi-rank path: the c-1 pre-SSH dynamics chain, the c-2 SSH
    ! stiffness + free-surface CG, AND the c-3 post-SSH ALE update (update_vel /
    ! compute_hbar_ale / update_eta_n / vert_vel_ale(+CFLz/Wvel_split) / update_thickness_ale)
    ! + tracer SOLVE (advection [M2.12b] + horizontal/implicit-vertical diffusion). It dumps
    ! every per-substep node field per rank (mod_dump, gid-keyed); tools/run_step_gate_multirank.sh
    ! compares FESOM3 dist_NP vs FESOM2 dist_NP (the REAL oce_timestep_ale via the oracle's
    ! fesom_step_dump npes>1 branch + its built-in dump_shim) by GLOBAL probe id for
    ! max|delta|=0 — the L8 same-partition rule.
    !
    ! This is the multi-rank analog of the 1-rank fesom_stepdump (whole step via step_oce);
    ! the c-1/c-2 fesom_stepdump_mr (dynamics chain -> d_eta, direct kernel calls) stays as a
    ! faster regression probe. The PRESCRIBED state + forcing inputs are byte-identical to the
    ! 1-rank fesom_stepdump AND the FESOM2 oracle (same analytic formulas, same rotated coords
    ! on the same partition), prescribed at owned+halo (node fields directly; element velocity
    ! owned then exchange_elem_full — halo elements have no local elem2D_nodes).
    !
    !   FESOM3_MESH_DIR   mesh dir (must contain dist_<NP>/)
    !   FESOM_DUMP_FILE   node dump prefix (mod_dump; per rank <prefix>.<mype5>)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: density_0
    use mod_param_phys,   only: N2smth_h, alpha, theta
    use mod_param_phys,   only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,   only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex, set_partition
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_halo,         only: exchange_elem_full
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use oce_muscl_adv,    only: muscl_adv_init
    use oce_ssh_rhs,      only: init_stiff_mat_ale
    use mod_step_oce,     only: step_oce
    use mod_dump,         only: dump_init, dump_finalize
    implicit none

    ! pi namelist timestep: dt = 86400/36 = 2400 s (= the 1-rank step gate / FESOM2 oracle).
    real(kind=WP), parameter :: dt = 86400.0_WP / real(36, WP)

    character(len=512) :: mesh_dir
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    integer :: nsw, n, nz, nl, nzmin, nzmax, e, tr_num
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), heat_flux(:), water_flux(:), virtual_salt(:)
    real(kind=WP), allocatable :: relax_salt(:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    integer, allocatable :: idlist(:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'

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

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells; pure functions of the halo-filled
    ! nlevels, built over the full local range — as in fesom_stepdump_mr).
    allocate(mesh%hnode(nl-1, nNodL)); mesh%hnode = 0.0_MP
    do n = 1, nNodL
        if (mesh%nlevels_nod2D(n) <= 0) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    allocate(mesh%helem(nl-1, nElemF)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(nElemF)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, nElemF
        if (mesh%nlevels(e) <= 0) cycle
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
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
    allocate(mesh%hbar(nNodL), mesh%hbar_old(nNodL), mesh%dhe(nElemF))
    allocate(mesh%hnode_new(nl-1, nNodL))
    mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode

    !===========================================================================
    ! dynamics state (local sizes; elements full halo, nodes owned+halo).
    allocate(dyn%uv(2, nl-1, nElemF), dyn%uv_rhs(2, nl-1, nElemF))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, nElemF))     ! (AB_order-1, 2, nl-1, elem)
    allocate(dyn%uvnode(2, nl-1, nNodL))
    allocate(dyn%eta_n(nNodL), dyn%d_eta(nNodL))
    allocate(dyn%ssh_rhs(nNodL), dyn%ssh_rhs_old(nNodL))
    allocate(dyn%w(nl, nNodL), dyn%w_e(nl, nNodL), dyn%w_i(nl, nNodL))
    allocate(dyn%cfl_z(nl, nNodL))
    allocate(dyn%work%density_ref(nl-1, nNodL), dyn%work%density_m_rho0(nl-1, nNodL))
    allocate(dyn%work%hpressure(nl, nNodL), dyn%work%bvfreq(nl, nNodL), dyn%work%bvfreq_raw(nl, nNodL), dyn%work%shear2(nl, nNodL))
    allocate(dyn%work%pgf_x(nl-1, nElemF), dyn%work%pgf_y(nl-1, nElemF))
    allocate(dyn%work%u_c(nl-1, nElemF), dyn%work%v_c(nl-1, nElemF))
    allocate(dyn%work%uvnode_rhs(2, nl-1, nNodL))
    allocate(dyn%work%vorticity(nl-1, nNodL))   ! momadv_opt==1 (vector-invariant) scratch
    dyn%work%vorticity = 0.0_WP
    allocate(dyn%work%Kv(nl, nNodL), dyn%work%Av(nl, nElemF))
    dyn%uv = 0.0_WP; dyn%uv_rhs = 0.0_WP; dyn%uv_rhsAB = 0.0_WP; dyn%uvnode = 0.0_WP
    dyn%eta_n = 0.0_WP; dyn%d_eta = 0.0_WP; dyn%ssh_rhs = 0.0_WP; dyn%ssh_rhs_old = 0.0_WP
    dyn%w = 0.0_WP; dyn%w_e = 0.0_WP; dyn%w_i = 0.0_WP; dyn%cfl_z = 0.0_WP
    dyn%work%density_m_rho0 = 0.0_WP; dyn%work%hpressure = 0.0_WP; dyn%work%bvfreq = 0.0_WP; dyn%work%bvfreq_raw = 0.0_WP; dyn%work%shear2 = 0.0_WP
    dyn%work%pgf_x = 0.0_WP; dyn%work%pgf_y = 0.0_WP
    dyn%work%u_c = 0.0_WP; dyn%work%v_c = 0.0_WP; dyn%work%uvnode_rhs = 0.0_WP
    dyn%work%Kv = 0.0_WP; dyn%work%Av = 0.0_WP
    dyn%work%density_ref = density_0
    dyn%AB_order      = 2
    dyn%momadv_opt    = 2
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    dyn%use_wsplit    = .false.
    dyn%wsplit_maxcfl = 1.0_WP

    !===========================================================================
    ! tracers (T/S) + advection machinery (local sizes; adv_flux_hor owned edges).
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%values   (nl-1, nNodL))
        allocate(tracers%data(tr_num)%valuesAB (nl-1, nNodL))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, nNodL))
    end do
    allocate(tracers%work%fct_LO          (nl-1, nNodL))
    allocate(tracers%work%adv_flux_hor    (nl-1, nEdgeO))
    allocate(tracers%work%adv_flux_ver    (nl,   nNodL))
    allocate(tracers%work%fct_ttf_max     (nl-1, nNodL))
    allocate(tracers%work%fct_ttf_min     (nl-1, nNodL))
    allocate(tracers%work%fct_plus        (nl-1, nNodL))
    allocate(tracers%work%fct_minus       (nl-1, nNodL))
    allocate(tracers%work%del_ttf         (nl-1, nNodL))
    allocate(tracers%work%del_ttf_advhoriz(nl-1, nNodL))
    allocate(tracers%work%del_ttf_advvert (nl-1, nNodL))

    !===========================================================================
    ! prescribe the analytic state at owned+halo (= the FESOM2 oracle / 1-rank
    ! fesom_stepdump, same coords on the same partition) — node fields directly;
    ! element velocity owned + exchange_elem_full.
    do n = 1, nNodL
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            tracers%data(1)%values(nz,n) = 12.0_WP + 8.0_WP*cos(lat)*cos(lon) - 0.20_WP*real(nz,WP) &
                        + 8.0_WP*max(0.0_WP, cos(lat)*cos(lon))*min(real(nz,WP),8.0_WP)/8.0_WP
            tracers%data(2)%values(nz,n) = 34.5_WP + 0.5_WP*sin(2.0_WP*lon)*cos(lat) + 0.03_WP*real(nz,WP)
        end do
        dyn%eta_n(n) = 0.5_WP*cos(lat)*sin(lon) + 0.3_WP*sin(2.0_WP*lat)
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz,WP))
            dyn%w_i(nz,n) = 2.0e-4_WP*sin(lon)*cos(2.0_WP*lat)*cos(0.25_WP*real(nz,WP))
        end do
        mesh%hbar(n) = 0.4_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
    end do

    allocate(stress_surf(2, nElemF)); stress_surf = 0.0_WP
    do e = 1, nElemO
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        do nz = 1, nl-1
            dyn%uv(1,nz,e) =  0.50_WP*cos(lat)*sin(lon)        - 0.005_WP*real(nz,WP)
            dyn%uv(2,nz,e) = -0.40_WP*sin(lat)*cos(2.0_WP*lon) + 0.004_WP*real(nz,WP)
            dyn%uv_rhsAB(1,1,nz,e) =  1.0e6_WP*sin(lon)*cos(lat)        + 1.0e4_WP*real(nz,WP)
            dyn%uv_rhsAB(1,2,nz,e) = -1.0e6_WP*cos(lon)*sin(2.0_WP*lat) - 1.0e4_WP*real(nz,WP)
        end do
        stress_surf(1,e) =  0.10_WP*cos(lat)*sin(lon)
        stress_surf(2,e) = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
    end do
    if (partit%npes > 1) call exchange_elem_full(dyn%uv, partit)

    ! FCT config + cold-start valuesold (= oracle / 1-rank fesom_stepdump).
    do tr_num = 1, 2
        tracers%data(tr_num)%ID          = tr_num
        tracers%data(tr_num)%AB_order    = 2
        tracers%data(tr_num)%tra_adv_hor = 'MUSCL'
        tracers%data(tr_num)%tra_adv_ver = 'QR4C'
        tracers%data(tr_num)%tra_adv_lim = 'FCT'
        tracers%data(tr_num)%tra_adv_ph  = 0.0_WP    ! MUSCL
        tracers%data(tr_num)%tra_adv_pv  = 1.0_WP    ! QR4C
        tracers%data(tr_num)%i_vert_diff = .true.
        tracers%data(tr_num)%valuesAB        = 0.0_WP
        tracers%data(tr_num)%valuesold(1,:,:) = tracers%data(tr_num)%values   ! cold start
        tracers%data(tr_num)%valuesold(2,:,:) = tracers%data(tr_num)%values
    end do
    call muscl_adv_init(tracers%work, mesh, partit)   ! edge_up_dn_tri

    !===========================================================================
    ! forcing inputs (D7 explicit args) at owned+halo — Ki is read at owned edges' nodes
    ! (halo) by the horizontal diffusion; the surface fluxes are read at owned nodes only.
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, nNodL), heat_flux(nNodL), water_flux(nNodL))
    allocate(virtual_salt(nNodL), relax_salt(nNodL), real_salt_flux(nNodL))
    do n = 1, nNodL
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            Ki(nz,n) = 300.0_WP + 150.0_WP*cos(lat)*cos(lon) + 10.0_WP*real(nz,WP)
        end do
        heat_flux(n)      =  50.0_WP*cos(lat)*sin(lon)
        water_flux(n)     =  1.0e-6_WP*sin(2.0_WP*lon)*cos(lat)
        virtual_salt(n)   =  3.0e-5_WP*sin(lon)*cos(lat)
        relax_salt(n)     =  2.0e-5_WP*cos(2.0_WP*lon)*cos(lat)
        real_salt_flux(n) =  0.0_WP
    end do

    !===========================================================================
    ! reduced-M2 module config (= the FESOM2 oracle / pi namelist).
    alpha = 1.0_WP; theta = 1.0_WP
    N2smth_h     = .true.
    mix_scheme_nmb = 2           ! reduced-M2 mixing = PP
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.
    use_windmix   = .false.

    ! SSH stiffness (built ONCE; dt = the pi namelist timestep).
    call init_stiff_mat_ale(mesh, dt, partit)

    !===========================================================================
    ! per-rank dump (mod_dump: gid-keyed probes, the rank owning a probe writes it).
    call dump_init(partit%mype, nNodO, partit%myList_nod2D, nElemO, partit%myList_elem2D)

    !===========================================================================
    ! run the ASSEMBLED ocean step ONCE through the optional partit (lfirst=.true. ->
    ! Euler start, like FESOM2 n=1); step_oce emits every per-substep node dump.
    call step_oce(1, dt, .true., dyn, tracers, mesh, Ki, &
                  heat_flux, water_flux, virtual_salt, relax_salt, &
                  real_salt_flux, is_nonlinfs, stress_surf, partit)

    call dump_finalize()
    write(*,'(a,i0,a)') 'fesom_stepfull_mr: rank ', partit%mype, ' done (whole step + per-substep dump).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_stepfull_mr
