program fesom_icefctdump
    ! M3c sea-ice FCT-advection byte-gate driver. Mirrors fesom_evpdump (M3b) init up to
    ! and including ocean2ice + EVPdynamics_solve (the EVP ice velocity uice/vice is M3b-
    ! byte-proven), then runs the ice FCT advection ice_TG_rhs + ice_fct_solve (the M3c
    ! kernels) on the cold-start ice tracers (a_ice/m_ice/m_snow) with that velocity, and
    ! dumps the TG rhs (rhs_a/m/ms — localises a TG_rhs vs an fct_solve bug) AND the post-
    ! advection a_ice/m_ice/m_snow. tools/run_icefct_gate_core2.sh compares vs the FESOM2
    ! oracle (src/fesom_ice_dump.F90::ice_fct_dump_write) for max|delta|=0, CORE2 1-rank.
    !
    ! ice_diff / ice_gamma_fct are set to the CORE2 namelist.ice values (0.0 / 0.5) NOT the
    ! t_ice single->WP defaults (10.0 / 0.25) — the M3a precision note (the oracle reads the
    ! namelist). cd_oce_ice / delta_min are the namelist doubles (as in fesom_evpdump).
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): the FESOM2-oracle
    ! byte gate is retired (docs/HANDOFF.md "Gate migration"). The driver has no triangle-
    ! specific code of its own; the ice kernels it runs are converted in Tasks 15/16. It stays
    ! a diagnostic dump.
    !
    !   FESOM3_MESH_DIR  mesh dir   (default: CORE2)
    !   FESOM3_IC_FILE   IC netcdf  (default: pool phc3.0_winter.nc — same file the oracle reads)
    !   FESOM3_WHICHEVP  0=EVP / 1=mEVP (default 0)
    !   FESOM3_FCT_OUT   out path   (default: icefct_f3.bin)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_param_phys,     only: alpha, theta
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use mod_ice_dyn,        only: ocean2ice, EVPdynamics_solve
    use mod_ice_fct,        only: ice_TG_rhs, ice_fct_solve
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_advhor_dump,    only: advhor_dump_open, advhor_dump_close, wr_r1
    implicit none

    ! CORE2 namelist timestep: dt = 86400/48 = 1800 s (computed like FESOM2, not a literal).
    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)

    character(len=512) :: mesh_dir, out_path, ic_file, whichevp_str
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ice)        :: ice
    type(t_ic3d_config):: ic
    integer :: nsw, n, nl, nzmin, nzmax, e, u, whichevp, ios
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    call get_environment_variable('FESOM3_FCT_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'icefct_f3.bin'
    ! whichEVP: 0 = standard EVP, 1 = modified EVP (mEVP). Default 0.
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', whichevp_str)
    if (len_trim(whichevp_str) > 0) read(whichevp_str, *, iostat=ios) whichevp

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_icefctdump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! mesh + geometry (CORE2 rotation 50/15/-90, identical to the lifecycle/IC/ice gate).
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_icefctdump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depths Z_3d_n (do_ic3d) + zbar_e_bot (init_stiff_mat_ale). linfs full cells,
    ! eta=0 — identical to fesom_evpdump.
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
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do
    ! ALE elevation state ocean2ice copies into srfoce_ssh (prescribed below).
    allocate(mesh%hbar(mesh%nod2D)); mesh%hbar = 0.0_MP

    !===========================================================================
    ! 2-tracer ocean state + do_ic3d phc3.0 IC (gives the surface T the cold start reads,
    ! and the surface T/S ocean2ice copies into srfoce_temp/salt).
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
    tracers%data(1)%values = 0.0_WP;  tracers%data(1)%ID = 1
    tracers%data(2)%values = 0.0_WP;  tracers%data(2)%ID = 2
    ic%n_ic3d      = 2
    ic%idlist(1:2) = [2, 1]
    ic%filelist(1) = trim(ic_file)
    ic%filelist(2) = trim(ic_file)
    ic%varlist(1)  = 'salt'
    ic%varlist(2)  = 'temp'
    ic%t_insitu    = .true.
    ic%ic_cyclic   = .true.
    ic%dummy       = 1.e10_WP
    call do_ic3d(tracers, ic, mesh)
    write(*,'(a,2es12.4)') 'fesom_icefctdump: surface T range = ', &
        minval(tracers%data(1)%values(1,:)), maxval(tracers%data(1)%values(1,:))

    !===========================================================================
    ! ssh_stiff CSR (reduced-M2: alpha=theta=1) — the ice mass matrix rides its sparsity.
    alpha = 1.0_WP; theta = 1.0_WP
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! ice setup: allocate + ice_mass_matrix_fill + cold-start ice_initial_state (uice=
    ! vice=0; a_ice/m_ice/m_snow from the SST<0 sign test). Then override the EVP drag /
    ! deformation params AND the FCT advection params to the CORE2 namelist values.
    call ice_setup(ice, tracers, mesh, dt)
    ice%whichEVP      = whichevp        ! 0 = standard EVP, 1 = mEVP (FESOM3_WHICHEVP)
    ice%cd_oce_ice    = 0.0055_WP       ! namelist Cd_oce_ice (double), not the 5.5e-3 single->WP default
    ice%delta_min     = 1.0e-11_WP      ! namelist delta_min  (double)
    ice%ice_diff      = 0.0_WP          ! namelist ice_diff (NO artificial diffusion), not the 10.0 default
    ice%ice_gamma_fct = 0.5_WP          ! namelist ice_gamma_fct, not the 0.25 default
    write(*,'(a,i0,a,3es12.4)') 'fesom_icefctdump: whichEVP=', whichevp, &
        '  max a_ice/m_ice/m_snow (IC) = ', &
        maxval(ice%data(1)%values), maxval(ice%data(2)%values), maxval(ice%data(3)%values)

    !===========================================================================
    ! dynamics: only ocean2ice's dyn%uv(:,1,:) is read. Prescribe the analytic surface
    ! ocean velocity (the deeper levels are zeroed, never read). MUST equal the oracle.
    allocate(dyn%uv(2, nl-1, mesh%elem2D)); dyn%uv = 0.0_WP
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        dyn%uv(1,1,e) =  0.20_WP*cos(lat)*sin(lon)
        dyn%uv(2,1,e) = -0.15_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    !===========================================================================
    ! prescribe the elevation (hbar, ~0.5 m) + the wind-on-ice stress (~0.1 N/m^2).
    ! MUST equal the FESOM2 oracle src/fesom_ice_dump.F90::ice_fct_dump_write.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        mesh%hbar(n)              =  0.5_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
        ice%stress_atmice_x(n)    =  0.10_WP*cos(lat)*sin(lon)
        ice%stress_atmice_y(n)    = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    !===========================================================================
    ! ocean -> ice coupling + EVP dynamics (120 subcycles) -> uice/vice (M3b-proven), then
    ! ice FCT advection (ice_TG_rhs -> ice_fct_solve) -> rhs + post-advection tracers (M3c).
    call ocean2ice(ice, dyn, tracers, mesh)
    call EVPdynamics_solve(ice, mesh)
    write(*,'(a,2es12.4)') 'fesom_icefctdump: max|uice| max|vice| = ', &
        maxval(abs(ice%uice)), maxval(abs(ice%vice))
    call ice_TG_rhs(ice, mesh)
    call ice_fct_solve(ice, mesh)
    write(*,'(a,3es12.4)') 'fesom_icefctdump: post-adv max a_ice/m_ice/m_snow = ', &
        maxval(ice%data(1)%values(1:mesh%nod2D)), &
        maxval(ice%data(2)%values(1:mesh%nod2D)), &
        maxval(ice%data(3)%values(1:mesh%nod2D))

    !===========================================================================
    ! dump (FADVHDMP, 1-rank identity ids). TG rhs + post-advection tracers (all nod2D).
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r1(u, 'ice_rhs_a',  ice%data(1)%values_rhs(1:mesh%nod2D))
    call wr_r1(u, 'ice_rhs_m',  ice%data(2)%values_rhs(1:mesh%nod2D))
    call wr_r1(u, 'ice_rhs_ms', ice%data(3)%values_rhs(1:mesh%nod2D))
    call wr_r1(u, 'ice_a_ice',  ice%data(1)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_ice',  ice%data(2)%values(1:mesh%nod2D))
    call wr_r1(u, 'ice_m_snow', ice%data(3)%values(1:mesh%nod2D))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_icefctdump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_icefctdump
