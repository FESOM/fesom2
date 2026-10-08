program fesom_stepdump
    ! M2.9b step-assembly byte-gate driver. Loads pi 1-rank with the SAME rotation as the
    ! FESOM2 pi run (50/15/-90, cyclic 360), builds the ALE depths / thicknesses, PRESCRIBES
    ! the clean initial ocean state (analytic T/S + UV/eta_n/w_e/w_i/uv_rhsAB/hbar + the
    ! forcing inputs Ki/heat_flux/water_flux/virtual_salt/relax_salt/stress_surf — identical
    ! formulas to the FESOM2 oracle src/fesom_step_dump.F90, from the byte-identical coords),
    ! builds the SSH stiffness, then runs the ASSEMBLED ocean step (mod_step_oce::step_oce)
    ! ONCE with LIVE data flow and dumps every per-substep node field (density/pressure/bvfreq
    ! / Kv / ssh_rhs / d_eta / hbar / eta_n / hnode_new / w / T / S / hnode) in the oracle's
    ! mod_dump format. tools/run_step_gate.sh compares against the REAL FESOM2 oce_timestep_ale
    ! (its built-in dump_shim_record_node) for max|delta|=0 (tools/dump_diff.py).
    !
    ! Unlike fesom_pressuredump (which runs each kernel ISOLATED on prescribed intermediates),
    ! this driver prescribes ONLY the true step inputs and lets step_oce produce uvnode (from
    ! UV), Kv/Av (PP+convection), d_eta (CG) etc. LIVE — so the gate tests the ASSEMBLY (the
    ! data flow between kernels), not the kernels in isolation.
    !
    !   FESOM3_MESH_DIR   mesh dir         (default: pi)
    !   FESOM_DUMP_FILE   node dump prefix (mod_dump; set by run_stepdump_pi.sh)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: density_0
    use mod_param_phys,   only: N2smth_h, alpha, theta
    use mod_param_phys,   only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,   only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_dyn,          only: t_dyn
    use mod_tracer,       only: t_tracer
    use oce_muscl_adv,    only: muscl_adv_init
    use oce_ssh_rhs,      only: init_stiff_mat_ale
    use mod_step_oce,     only: step_oce
    use mod_dump,         only: dump_init, dump_finalize
    implicit none

    ! pi namelist timestep (the value FESOM2 uses uniformly in oce_timestep_ale):
    ! dt = 86400/step_per_day, step_per_day=36 -> 2400 s. Compute it the SAME way FESOM2
    ! does (gen_model_setup.F90:92), NOT a 2400.0 literal, so the -no-prec-div bits agree.
    real(kind=WP), parameter :: dt = 86400.0_WP / real(36, WP)

    character(len=512) :: mesh_dir
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    type(t_dyn)         :: dyn
    type(t_tracer)      :: tracers
    integer :: nsw, n, nz, nl, nzmin, nzmax, e, tr_num
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
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_stepdump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_stepdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells) — identical to fesom_pressuredump.
    allocate(mesh%hnode(nl-1, mesh%nod2D)); mesh%hnode = 0.0_MP
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
    allocate(mesh%zbar_3d_n(nl, mesh%nod2D), mesh%Z_3d_n(nl-1, mesh%nod2D))
    mesh%zbar_3d_n = 0.0_MP; mesh%Z_3d_n = 0.0_MP
    do n = 1, mesh%nod2D
        nzmin = mesh%ulevels_nod2D(n)
        nzmax = mesh%nlevels_nod2D(n)
        zbar_srf = mesh%zbar(nzmin)
        zbar_bot = mesh%zbar(nzmax)
        mesh%zbar_3d_n(1:nzmin-1, n)       = mesh%zbar(1:nzmin-1)
        mesh%zbar_3d_n(nzmin, n)           = zbar_srf
        mesh%zbar_3d_n(nzmin+1:nzmax-1, n) = mesh%zbar(nzmin+1:nzmax-1)
        mesh%zbar_3d_n(nzmax, n)           = zbar_bot
        mesh%Z_3d_n(1:nzmin-1, n)          = mesh%Z(1:nzmin-1)
        mesh%Z_3d_n(nzmin, n)              = mesh%zbar_3d_n(nzmin,n)   + (mesh%zbar_3d_n(nzmin+1,n)-zbar_srf)/2
        mesh%Z_3d_n(nzmin+1:nzmax-2, n)    = mesh%Z(nzmin+1:nzmax-2)
        mesh%Z_3d_n(nzmax-1, n)            = mesh%zbar_3d_n(nzmax-1,n) + (zbar_bot-mesh%zbar_3d_n(nzmax-1,n))/2
    end do
    ! ALE elevation state the update reads/writes (hbar prescribed below).
    allocate(mesh%hbar(mesh%nod2D), mesh%hbar_old(mesh%nod2D), mesh%dhe(mesh%elem2D))
    allocate(mesh%hnode_new(nl-1, mesh%nod2D))
    mesh%hbar_old  = 0.0_MP
    mesh%dhe       = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode

    !===========================================================================
    ! dynamics state + reduced-M2 config (= fesom_pressuredump + the FESOM2 namelist).
    allocate(dyn%uv(2, nl-1, mesh%elem2D), dyn%uv_rhs(2, nl-1, mesh%elem2D))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, mesh%elem2D))     ! (AB_order-1, 2, nl-1, elem2D)
    allocate(dyn%uvnode(2, nl-1, mesh%nod2D))
    allocate(dyn%eta_n(mesh%nod2D), dyn%d_eta(mesh%nod2D))
    allocate(dyn%ssh_rhs(mesh%nod2D), dyn%ssh_rhs_old(mesh%nod2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%w_e(nl, mesh%nod2D), dyn%w_i(nl, mesh%nod2D))
    allocate(dyn%cfl_z(nl, mesh%nod2D))
    ! work scratch (zeroed at setup, mirroring FESOM2 arrays_init; pressure_bv/PP overwrite)
    allocate(dyn%work%density_ref(nl-1, mesh%nod2D), dyn%work%density_m_rho0(nl-1, mesh%nod2D))
    allocate(dyn%work%hpressure(nl, mesh%nod2D), dyn%work%bvfreq(nl, mesh%nod2D), dyn%work%bvfreq_raw(nl, mesh%nod2D), dyn%work%shear2(nl, mesh%nod2D))
    allocate(dyn%work%pgf_x(nl-1, mesh%elem2D), dyn%work%pgf_y(nl-1, mesh%elem2D))
    allocate(dyn%work%u_c(nl-1, mesh%elem2D), dyn%work%v_c(nl-1, mesh%elem2D))
    allocate(dyn%work%uvnode_rhs(2, nl-1, mesh%nod2D))
    allocate(dyn%work%vorticity(nl-1, mesh%nod2D))   ! momadv_opt==1 (vector-invariant) scratch
    dyn%work%vorticity = 0.0_WP
    allocate(dyn%work%Kv(nl, mesh%nod2D), dyn%work%Av(nl, mesh%elem2D))
    dyn%work%density_m_rho0 = 0.0_WP; dyn%work%hpressure = 0.0_WP; dyn%work%bvfreq = 0.0_WP; dyn%work%bvfreq_raw = 0.0_WP; dyn%work%shear2 = 0.0_WP
    dyn%work%pgf_x = 0.0_WP; dyn%work%pgf_y = 0.0_WP
    dyn%work%u_c = 0.0_WP; dyn%work%v_c = 0.0_WP; dyn%work%uvnode_rhs = 0.0_WP
    dyn%work%Kv = 0.0_WP; dyn%work%Av = 0.0_WP    ! surface/bottom stay 0 (PP sets interior)
    dyn%work%density_ref = density_0              ! use_density_ref=.false. on pi
    dyn%uv_rhs = 0.0_WP
    dyn%d_eta = 0.0_WP; dyn%ssh_rhs = 0.0_WP; dyn%ssh_rhs_old = 0.0_WP
    dyn%w = 0.0_WP; dyn%cfl_z = 0.0_WP
    dyn%AB_order      = 2
    dyn%momadv_opt    = 2
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    ! use_wsplit=.false. for the step gate (M1.4 precedent): the FESOM2 oracle of this
    ! dump has the hard CFL switch, FESOM3 the smooth Shchepetkin split (oce_wsplit), so
    ! the use_wsplit=.true. path (incl. the FCT adv_tra_vert_impl of do_oce_adv_tra) is
    ! not byte-comparable here; its nets are test_wsplit / test_wimpl_tra and the
    ! fesom_conserve wsplit gate configs. With .false. the explicit/implicit split is
    ! trivial (w_e=w, w_i=0) AFTER vert_vel_ale; impl_vert_visc_ale still runs on the
    ! prescribed w_i (it precedes compute_Wvel_split), so vertical momentum advection is
    ! exercised.
    dyn%use_wsplit    = .false.
    dyn%wsplit_maxcfl = 1.0_WP

    ! prescribe T/S (with the M2.8b unstable band) — identical to fesom_pressuredump.
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            tracers%data(1)%values(nz,n) = 12.0_WP + 8.0_WP*cos(lat)*cos(lon) - 0.20_WP*real(nz,WP) &
                        + 8.0_WP*max(0.0_WP, cos(lat)*cos(lon))*min(real(nz,WP),8.0_WP)/8.0_WP
            tracers%data(2)%values(nz,n) = 34.5_WP + 0.5_WP*sin(2.0_WP*lon)*cos(lat) + 0.03_WP*real(nz,WP)
        end do
    end do

    ! prescribe UV (elements) / uv_rhsAB (previous step) + w_e/w_i + eta_n + stress_surf,
    ! identical to fesom_pressuredump (the strong-current 2.0/1.5 stress test).
    allocate(stress_surf(2, mesh%elem2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        dyn%eta_n(n) = 0.5_WP*cos(lat)*sin(lon) + 0.3_WP*sin(2.0_WP*lat)
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz,WP))
            dyn%w_i(nz,n) = 2.0e-4_WP*sin(lon)*cos(2.0_WP*lat)*cos(0.25_WP*real(nz,WP))
        end do
    end do
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        ! UV amplitude 0.50/0.40 m/s (a strong but PHYSICAL current): keeps the
        ! single-step elevation eta_n within the +/-10 m blowup guard so the REAL FESOM2
        ! oce_timestep_ale completes cleanly (the viscosity-gate's 2.0/1.5 m/s stress test
        ! drives eta_n ~13 -> check_blowup; that branch coverage is gated at M2.4). The
        ! depth terms (-0.005*nz/+0.004*nz) give vertical shear -> non-trivial uvnode/PP and
        ! momentum advection. MUST equal the FESOM2 oracle src/fesom_step_dump.F90.
        do nz = 1, nl-1
            dyn%uv(1,nz,e) =  0.50_WP*cos(lat)*sin(lon)        - 0.005_WP*real(nz,WP)
            dyn%uv(2,nz,e) = -0.40_WP*sin(lat)*cos(2.0_WP*lon) + 0.004_WP*real(nz,WP)
            dyn%uv_rhsAB(1,1,nz,e) =  1.0e6_WP*sin(lon)*cos(lat)        + 1.0e4_WP*real(nz,WP)
            dyn%uv_rhsAB(1,2,nz,e) = -1.0e6_WP*cos(lon)*sin(2.0_WP*lat) - 1.0e4_WP*real(nz,WP)
        end do
        stress_surf(1,e) =  0.10_WP*cos(lat)*sin(lon)
        stress_surf(2,e) = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    ! prescribe the previous-step elevation hbar (sign-varying ~0.6 m).
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        mesh%hbar(n) = 0.4_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
    end do

    !===========================================================================
    ! tracer advection machinery: valuesAB / valuesold (cold start = values) + the FCT
    ! work arrays + muscl_adv_init, and the pi FCT config (MUSCL/QR4C/FCT, opth=0/optv=1).
    tracers%num_tracers = 2
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, mesh%nod2D))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, mesh%nod2D))
        tracers%data(tr_num)%ID         = tr_num
        tracers%data(tr_num)%AB_order   = 2
        tracers%data(tr_num)%tra_adv_hor = 'MUSCL'
        tracers%data(tr_num)%tra_adv_ver = 'QR4C'
        tracers%data(tr_num)%tra_adv_lim = 'FCT'
        tracers%data(tr_num)%tra_adv_ph  = 0.0_WP    ! MUSCL
        tracers%data(tr_num)%tra_adv_pv  = 1.0_WP    ! QR4C
        tracers%data(tr_num)%i_vert_diff = .true.
        tracers%data(tr_num)%valuesAB        = 0.0_WP
        tracers%data(tr_num)%valuesold(1,:,:) = tracers%data(tr_num)%values   ! cold start: T^{n-1}=T^n
        tracers%data(tr_num)%valuesold(2,:,:) = tracers%data(tr_num)%values
    end do
    allocate(tracers%work%fct_LO          (nl-1, mesh%nod2D))
    allocate(tracers%work%adv_flux_hor    (nl-1, mesh%edge2D))
    allocate(tracers%work%adv_flux_ver    (nl,   mesh%nod2D))
    allocate(tracers%work%fct_ttf_max     (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_ttf_min     (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_plus        (nl-1, mesh%nod2D))
    allocate(tracers%work%fct_minus       (nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf         (nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf_advhoriz(nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf_advvert (nl-1, mesh%nod2D))
    call muscl_adv_init(tracers%work, mesh)     ! edge_up_dn_tri

    !===========================================================================
    ! forcing inputs (D7 explicit args; M2.10 forcing / M4 resolution will source them).
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, mesh%nod2D), heat_flux(mesh%nod2D), water_flux(mesh%nod2D))
    allocate(virtual_salt(mesh%nod2D), relax_salt(mesh%nod2D), real_salt_flux(mesh%nod2D))
    do n = 1, mesh%nod2D
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
    N2smth_h     = .true.          ! horizontal N^2 smoothing ON (the real-step path)
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
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! open the per-substep node dump (mod_dump; FESOM_DUMP_FILE set by the run script).
    allocate(idlist(max(mesh%nod2D, mesh%elem2D)))
    do n = 1, size(idlist); idlist(n) = n; end do   ! 1-rank: identity global ids
    call dump_init(partit%mype, mesh%nod2D, idlist(1:mesh%nod2D), &
                   mesh%elem2D, idlist(1:mesh%elem2D))

    !===========================================================================
    ! run the assembled ocean step ONCE (lfirst=.true. -> Euler start, like FESOM2 n=1).
    call step_oce(1, dt, .true., dyn, tracers, mesh, Ki, &
                  heat_flux, water_flux, virtual_salt, relax_salt, &
                  real_salt_flux, is_nonlinfs, stress_surf)

    call dump_finalize()
    write(*,'(a)') 'fesom_stepdump: done (assembled step + per-substep node dump).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_stepdump
