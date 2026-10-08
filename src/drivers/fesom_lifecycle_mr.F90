program fesom_lifecycle_mr
    ! Production-validation MULTI-RANK lifecycle driver (the multi-rank analog of the
    ! 1-rank fesom_lifecycle). Runs on dist_<NP> (NP>=1) the REAL free-running reduced-M2
    ! lifecycle: model_init (MR mesh remap + geometry + ALE depths + do_ic3d phc3.0 IC
    ! through the optional partit + cold-start state + SSH stiffness) -> runloop(N) calling
    ! step_oce(n, dt, lfirst=(n==1), ..., partit) with persistent dyn/tracer state (multi-
    ! step AB2). step_oce emits the per-substep gid-keyed node dumps (mod_dump); the byte-
    ! gate (tools/run_lifecycle_gate_multirank.sh) compares them PER RANK vs the REAL FESOM2
    ! oce_timestep_ale lifecycle on the SAME dist_<NP> partition (L8 same-partition rule).
    !
    ! This closes the one combination the M2.12 single-step prescribe-and-stop gates never
    ! exercised: multi-rank x FREE-RUNNING x many steps (the AB2 history + cross-step halo
    ! state w_e/w_i/eta carried forward). UNFORCED (oracle use_ice=.false. -> all surface
    ! fluxes 0; FESOM3 prescribes Ki/heat_flux/water_flux/virtual_salt/relax_salt/stress = 0).
    !
    !   FESOM3_MESH_DIR   mesh dir (must contain dist_<NP>/)  default: CORE2
    !   FESOM3_IC_FILE    IC netcdf                           default: pool phc3.0_winter.nc
    !   FESOM_DUMP_FILE   per-rank node dump prefix (mod_dump; <prefix>.<mype5>)
    !   FESOM_DUMP_MAXSTEPS  dump step cap (mod_dump)
    !   FESOM3_NSTEPS     number of steps (default: 3)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32
    use mod_precision,      only: WP, MP
    use mod_constants,      only: density_0
    use mod_param_phys,     only: N2smth_h, alpha, theta
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex, set_partition
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_muscl_adv,      only: muscl_adv_init
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_step_oce,       only: step_oce
    use mod_dump,           only: dump_init, dump_finalize
    implicit none

    ! CORE2 namelist timestep: dt = 86400/step_per_day, step_per_day=48 -> 1800 s. Computed
    ! the SAME way FESOM2 does (gen_model_setup.F90:92), NOT an 1800.0 literal.
    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)

    character(len=512) :: mesh_dir, ic_file, env
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ic3d_config):: ic
    integer :: n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), heat_flux(:), water_flux(:), virtual_salt(:)
    real(kind=WP), allocatable :: relax_salt(:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    integer, allocatable :: idlist(:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    nsteps = 3
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps

    !===========================================================================
    ! model_init: MR mesh remap + geometry (set_partition -> read_mesh dispatches to
    ! read_mesh_local at npes>1; npes==1 reads the global mesh) — same as fesom_stepfull_mr.
    call par_init(partit)
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=n)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    if (partit%npes == 1) then
        nNodO = mesh%nod2D;   nNodL  = mesh%nod2D
        nEdgeO = mesh%edge2D; nElemO = mesh%elem2D; nElemF = mesh%elem2D
    else
        nNodO  = partit%myDim_nod2D
        nNodL  = partit%myDim_nod2D + partit%eDim_nod2D
        nEdgeO = partit%myDim_edge2D
        nElemO = partit%myDim_elem2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    if (partit%mype == 0) &
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_lifecycle_mr: nod2D=', mesh%nod2D, &
            ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells; local sizes) — as fesom_stepfull_mr.
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
    ! ALE elevation state — COLD START (eta=hbar=0).
    allocate(mesh%hbar(nNodL), mesh%hbar_old(nNodL), mesh%dhe(nElemF))
    allocate(mesh%hnode_new(nl-1, nNodL))
    mesh%hbar = 0.0_MP; mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode

    !===========================================================================
    ! dynamics state — COLD START (UV/eta_n/w/uv_rhsAB/d_eta/ssh_rhs_old = 0; local sizes).
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
    dyn%work%density_ref = density_0              ! use_density_ref=.false.
    dyn%AB_order      = 2
    dyn%momadv_opt    = 2
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    dyn%use_wsplit    = .false.          ! CORE2 production (M1.4/M2.9b precedent)
    dyn%wsplit_maxcfl = 1.0_WP

    !===========================================================================
    ! 2-tracer state (data(1)=T ID 1, data(2)=S ID 2; local sizes) + do_ic3d phc3.0 IC.
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, nNodL), tracers%data(2)%values(nl-1, nNodL))
    tracers%data(1)%values = 0.0_WP;  tracers%data(1)%ID = 1
    tracers%data(2)%values = 0.0_WP;  tracers%data(2)%ID = 2
    ! IC config = work_core/namelist.tra &tracer_init3d (idlist=2,1 -> salt first, temp second)
    ic%n_ic3d      = 2
    ic%idlist(1:2) = [2, 1]
    ic%filelist(1) = trim(ic_file)
    ic%filelist(2) = trim(ic_file)
    ic%varlist(1)  = 'salt'
    ic%varlist(2)  = 'temp'
    ic%t_insitu    = .true.
    ic%ic_cyclic   = .true.
    ic%dummy       = 1.e10_WP
    call do_ic3d(tracers, ic, mesh, partit)
    if (partit%mype == 0) &
        write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle_mr: IC(rank0 owned) T=', &
            minval(tracers%data(1)%values(:,1:nNodO)), maxval(tracers%data(1)%values(:,1:nNodO)), &
            '  S=', minval(tracers%data(2)%values(:,1:nNodO)), maxval(tracers%data(2)%values(:,1:nNodO))

    !===========================================================================
    ! tracer advection machinery: valuesAB / valuesold (cold start = values) + FCT work
    ! arrays (local sizes; adv_flux_hor owned edges) + muscl_adv_init.
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, nNodL))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, nNodL))
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
    call muscl_adv_init(tracers%work, mesh, partit)

    !===========================================================================
    ! surface forcing — UNFORCED: all zero (oracle runs use_ice=.false.; local sizes).
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, nNodL), heat_flux(nNodL), water_flux(nNodL))
    allocate(virtual_salt(nNodL), relax_salt(nNodL), real_salt_flux(nNodL))
    allocate(stress_surf(2, nElemF))
    Ki = 0.0_WP; heat_flux = 0.0_WP; water_flux = 0.0_WP; virtual_salt = 0.0_WP
    relax_salt = 0.0_WP; real_salt_flux = 0.0_WP; stress_surf = 0.0_WP

    !===========================================================================
    ! reduced-M2 module config (= the FESOM2 oracle / CORE2 namelist).
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

    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt, partit)

    !===========================================================================
    ! per-rank dump (mod_dump: gid-keyed probes, the rank owning a probe writes it).
    allocate(idlist(max(nNodO, nElemO)))
    do n = 1, size(idlist); idlist(n) = n; end do
    if (partit%npes == 1) then
        call dump_init(partit%mype, nNodO, idlist(1:nNodO), nElemO, idlist(1:nElemO))
    else
        call dump_init(partit%mype, nNodO, partit%myList_nod2D, nElemO, partit%myList_elem2D)
    end if

    !===========================================================================
    ! runloop: N steps, state persists (multi-step AB2 evolution). lfirst=(n==1).
    do n = 1, nsteps
        call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                      heat_flux, water_flux, virtual_salt, relax_salt, &
                      real_salt_flux, is_nonlinfs, stress_surf, partit)
        if (partit%mype == 0) &
            write(*,'(a,i0,a,es12.4,a,es12.4)') 'fesom_lifecycle_mr: step ', n, &
                '  max|eta_n(owned)|=', maxval(abs(dyn%eta_n(1:nNodO))), &
                '  max|uv(owned)|=', maxval(abs(dyn%uv(:,:,1:nElemO)))
    end do

    call dump_finalize()
    if (partit%mype == 0) &
        write(*,'(a,i0,a)') 'fesom_lifecycle_mr: done (', nsteps, ' steps + per-substep node dump).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_lifecycle_mr
