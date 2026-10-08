program fesom_lifecycle_native
    ! M3f native-flux forced lifecycle driver — the multi-step coupled sea-ice + ocean run
    ! with NATIVELY computed air-sea fluxes (the M3 payoff assembled into a lifecycle). It is
    ! the analog of fesom_lifecycle (M2.11c) but REPLACES the prescribed FESOM3_FLUX_FILE dump
    ! (heat_flux/water_flux/virtual_salt/relax_salt/stress_surf) with the native chain:
    !
    !   ocean2ice          ocean surface -> srfoce_u/v/temp/salt/ssh             (M3b)
    !   ice_timestep       EVP -> FCT -> cut_off -> thermodynamics                (M3f-1)
    !   oce_fluxes_mom     ice+atm momentum stress -> stress_surf                 (M3e)
    !   oce_fluxes         heat/freshwater/salt budget -> the tracer/SSH BCs      (M3e)
    !   step_oce           the whole ocean step on the NATIVE fluxes              (M2.9b)
    !
    ! Per-step ATMOSPHERIC forcing (the netCDF read + NCAR bulk: shortwave/longwave/Tair/
    ! shum/prec/runoff/wind/Ch-Ce_atm_oce/stress_atmoce/stress_atmice/Ssurf) is PRESCRIBED
    ! from the oracle dump FESOM3_ATMFLUX_FILE (written by port2 fesom_atmflux_dump.F90) so
    ! the native sea-ice + air-sea coupling is isolated from the CORE2 forcing read (which
    ! gets its own gate in M3f-3). This mirrors the established prescribe-the-input discipline:
    ! everything UPSTREAM of the kernel under test is fed from the oracle, the kernel runs
    ! native, and the byte-gate compares the DOWNSTREAM ocean substeps (the 195-record set).
    !
    ! The ocean init + step_oce are byte-for-byte the proven M2.11c-2 forced lifecycle; the
    ! ONLY change is the flux source (oracle dump -> native ice/oce_fluxes). The prognostic
    ! ice state (sigma elastic memory, uice/vice, t_skin, values_old) persists across steps in
    ! the live `ice` object — so this is the FIRST test of the MULTI-STEP ice evolution (the
    ! M3a-e gates were all single-step from cold start).
    !
    !   FESOM3_MESH_DIR     mesh dir            (default: CORE2)
    !   FESOM3_IC_FILE      IC netcdf           (default: pool phc3.0_winter.nc)
    !   FESOM3_ATMFLUX_FILE per-step atm dump   (REQUIRED; from fesom_atmflux_dump)
    !   FESOM3_FLUX_FILE    oracle flux dump    (OPTIONAL; if set, per-step native-vs-oracle
    !                                            flux self-check is printed to localize drift)
    !   FESOM3_WHICHEVP     0=EVP / 1=mEVP       (default 0)
    !   FESOM_DUMP_FILE     node dump prefix    (mod_dump; set by run script)
    !   FESOM_DUMP_MAXSTEPS dump step cap
    !   FESOM3_NSTEPS       number of steps     (default: 3)
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP
    use mod_constants,      only: density_0
    use mod_param_phys,     only: N2smth_h, alpha, theta, scale_area
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    ! M4e: GM bolus (FESOM3_FER_GM) + Redi isopycnal diffusion (FESOM3_REDI) in the
    ! forced/native lifecycle — the work_core GM+Redi config (copied from fesom_lifecycle.F90).
    use mod_param_phys,     only: Fer_GM, Redi, K_GM_max, K_GM_min, K_GM_bvref, &
                                  K_GM_rampmax, K_GM_rampmin, K_GM_resscalorder, K_GM_cm, &
                                  K_GM_cmin, K_GM_Ktaper, scaling_Ferreira, scaling_Rossby, &
                                  scaling_resolution, scaling_FESOM14, scaling_GMzexp, &
                                  scaling_GINsea, GMzexp_zref, GMzexp_smin
    use mod_param_phys,     only: Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                                  scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    ! M5c: KPP vertical mixing (FESOM3_MIX_KPP) + shortwave penetration (FESOM3_SW_PENE) +
    ! the KPP nonlocal counter-gradient flux (FESOM3_KPP_NONLCL — dead in production).
    use mod_param_phys,     only: Ricr, concv, visc_sh_limit, diff_sh_limit, &
                                  use_kpp_nonlclflx, ref_sss, ref_sss_local
    use mod_param_phys,     only: tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, &
                                  tke_kappaM_min, tke_kappaM_max, tke_min, tke_surf_min, &
                                  tke_mxl_choice, tke_only, tke_use_ubound_dirichlet, &
                                  tke_use_lbound_dirichlet, tke_dolangmuir
    use mod_config,         only: use_sw_pene, which_ALE
    use oce_mixing_kpp,     only: oce_mixing_kpp_init
    use oce_mixing_tke,     only: tke_init
    use oce_shortwave_pene, only: cal_shortwave_rad
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex
    use mod_mesh_read,      only: read_mesh
    use mod_mesh_areas,     only: compute_geometry
    use mod_dyn,            only: t_dyn
    use mod_tracer,         only: t_tracer
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use mod_ice_dyn,        only: ocean2ice
    use mod_ice_step,       only: ice_timestep
    use mod_ice_thermo,     only: t_atmflux
    use mod_ice_oce_coupling, only: oce_fluxes_mom, oce_fluxes
    use mod_forcing_read    ! M3f-3: native CORE2 atmospheric forcing read (self-check)
    use mod_forcing_bulk,   only: forcing_bulk_ncar, forcing_wind_stress, forcing_ice_stress
    use mod_forcing_other,  only: read_other_NetCDF   ! M3f-3b: runoff + SSS climatology read
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_muscl_adv,      only: muscl_adv_init
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_step_oce,       only: step_oce
    use mod_dump,           only: dump_init, dump_finalize
    implicit none

    real(kind=WP), parameter :: dt = 86400.0_WP / real(48, WP)   ! CORE2 dt = 1800 s

    character(len=512) :: mesh_dir, ic_file, env, atm_file, flux_file, whichevp_str
    type(t_partit)     :: partit
    type(t_mesh)       :: mesh
    type(t_dyn)        :: dyn
    type(t_tracer)     :: tracers
    type(t_ice)        :: ice
    type(t_atmflux)    :: atm
    type(t_ic3d_config):: ic
    integer :: n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len, nsw, whichevp
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), heat_flux(:), water_flux(:), virtual_salt(:)
    real(kind=WP), allocatable :: relax_salt(:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    integer :: atm_unit, flux_unit
    integer(int32) :: fstep, fnn, fne
    logical :: selfcheck, prescribe_atm
    ! M4e: enable GM bolus (FESOM3_FER_GM) / Redi isopycnal diffusion (FESOM3_REDI).
    logical :: use_fer_gm, use_redi
    ! M5c: KPP vertical mixing (use_kpp) + shortwave penetration (do_swpene) + the KPP
    ! nonlocal counter-gradient flux (do_nonlcl — DEAD in the production config). chl is the
    ! constant chlorophyll (work_core chl_const=0.1) cal_shortwave_rad reads (floored in place).
    logical :: use_kpp, do_swpene, do_nonlcl
    logical :: use_tke      ! M7c: FESOM3_MIX_TKE -> cvmix_TKE producer (forced; pair FESOM3_SW_PENE=1)
    real(kind=WP), allocatable :: chl(:)
    ! oracle flux arrays for the optional self-check
    real(kind=WP), allocatable :: o_hf(:), o_wf(:), o_vs(:), o_rs(:), o_ss(:,:)
    ! M3f-3 native CORE2 forcing read (self-check): when FESOM3_FORCING_DIR is set, the 8 NCAR
    ! fields + the NCAR bulk (Ch/Ce/stress_atmoce) + the wind-on-ice stress are recomputed
    ! NATIVELY each step and compared to the prescribed atm (the step itself still USES the
    ! prescribed values — the proven M3f-2 path — so the 195-record gate cannot regress; the
    ! self-check independently proves the native forcing is byte-exact). runoff + Ssurf stay
    ! prescribed (the monthly-climatology read is M3f-3b).
    character(len=512)  :: forcing_dir, runoff_file, sss_file
    logical             :: native_forcing
    type(t_atm_forcing) :: frc
    integer             :: fld, fyear
    real(kind=WP)       :: rdate_cold, rdate, timenew
    real(kind=WP), allocatable :: nuw(:), nvw(:), nta(:), nsh(:), nswr(:), nlw(:), npr(:), nps(:)
    real(kind=WP), allocatable :: ncd(:), nch(:), nce(:), nsx(:), nsy(:), nix(:), niy(:)
    ! M3f-3b: native runoff + SSS (read once — CORE runoff is time-constant, the SSS climatology
    ! is read once at mstep==1, i=month=1; both stay constant across a short January run).
    real(kind=WP), allocatable :: nro(:), nss(:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    ! FESOM3_ATMFLUX_FILE present => PRESCRIBE the per-step atmosphere from the oracle dump
    ! (the M3f-2 isolate-the-coupling path). ABSENT => FULLY NATIVE (M3f-3 complete): the
    ! whole atmosphere (NCAR read + bulk + 2 stresses + runoff + Ssurf) is computed natively,
    ! which then REQUIRES FESOM3_FORCING_DIR (the native NCAR read).
    call get_environment_variable('FESOM3_ATMFLUX_FILE', atm_file)
    prescribe_atm = (len_trim(atm_file) > 0)
    nsteps = 3
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', whichevp_str)
    if (len_trim(whichevp_str) > 0) read(whichevp_str, *, iostat=ios) whichevp
    call get_environment_variable('FESOM3_FER_GM', env, length=env_len, status=ios)
    use_fer_gm = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_REDI', env, length=env_len, status=ios)
    use_redi = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_MIX_KPP', env, length=env_len, status=ios)
    use_kpp = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_MIX_TKE', env, length=env_len, status=ios)
    use_tke = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_SW_PENE', env, length=env_len, status=ios)
    do_swpene = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_KPP_NONLCL', env, length=env_len, status=ios)
    do_nonlcl = (ios == 0 .and. env_len > 0)
    ! M6a-3: ALE vertical coordinate (default linfs). 'zstar' -> the full free surface
    ! (Shchepetkin PGF + per-step stiffness + thickness stretch) AND the REAL freshwater
    ! flux path (use_virt_salt=.false., is_nonlinfs=1 -> the water_flux terms in compute_hbar_ale
    ! /vert_vel_ale + real_salt_flux in bc_surface; oce_fluxes takes its levitating-ice branch).
    call get_environment_variable('FESOM3_WHICH_ALE', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) which_ALE = trim(env)

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_lifecycle_native: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    !===========================================================================
    ! mesh + geometry (CORE2 rotation 50/15/-90, identical to the lifecycle/IC/ice gate)
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_lifecycle_native: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' nl=', nl, ' CW swaps=', nsw

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells) — identical to fesom_lifecycle.
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
    allocate(mesh%hbar(mesh%nod2D), mesh%hbar_old(mesh%nod2D), mesh%dhe(mesh%elem2D))
    allocate(mesh%hnode_new(nl-1, mesh%nod2D))
    mesh%hbar = 0.0_MP; mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode

    !===========================================================================
    ! dynamics state — COLD START (UV/eta_n/w/uv_rhsAB/d_eta/ssh_rhs_old = 0).
    allocate(dyn%uv(2, nl-1, mesh%elem2D), dyn%uv_rhs(2, nl-1, mesh%elem2D))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, mesh%elem2D))
    allocate(dyn%uvnode(2, nl-1, mesh%nod2D))
    allocate(dyn%eta_n(mesh%nod2D), dyn%d_eta(mesh%nod2D))
    allocate(dyn%ssh_rhs(mesh%nod2D), dyn%ssh_rhs_old(mesh%nod2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%w_e(nl, mesh%nod2D), dyn%w_i(nl, mesh%nod2D))
    allocate(dyn%cfl_z(nl, mesh%nod2D))
    allocate(dyn%work%density_ref(nl-1, mesh%nod2D), dyn%work%density_m_rho0(nl-1, mesh%nod2D))
    allocate(dyn%work%hpressure(nl, mesh%nod2D), dyn%work%bvfreq(nl, mesh%nod2D), dyn%work%bvfreq_raw(nl, mesh%nod2D), dyn%work%shear2(nl, mesh%nod2D))
    allocate(dyn%work%pgf_x(nl-1, mesh%elem2D), dyn%work%pgf_y(nl-1, mesh%elem2D))
    allocate(dyn%work%u_c(nl-1, mesh%elem2D), dyn%work%v_c(nl-1, mesh%elem2D))
    allocate(dyn%work%uvnode_rhs(2, nl-1, mesh%nod2D))
    allocate(dyn%work%vorticity(nl-1, mesh%nod2D))   ! momadv_opt==1 (vector-invariant) scratch
    dyn%work%vorticity = 0.0_WP
    allocate(dyn%work%Kv(nl, mesh%nod2D), dyn%work%Av(nl, mesh%elem2D))
    dyn%uv = 0.0_WP; dyn%uv_rhs = 0.0_WP; dyn%uv_rhsAB = 0.0_WP; dyn%uvnode = 0.0_WP
    dyn%eta_n = 0.0_WP; dyn%d_eta = 0.0_WP
    dyn%ssh_rhs = 0.0_WP; dyn%ssh_rhs_old = 0.0_WP
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
    ! 2-tracer ocean state (data(1)=T, data(2)=S) + do_ic3d phc3.0 IC.
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
    write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle_native: IC T=', &
        minval(tracers%data(1)%values), maxval(tracers%data(1)%values), &
        '  S=', minval(tracers%data(2)%values), maxval(tracers%data(2)%values)

    !===========================================================================
    ! tracer advection machinery (cold start = values; MUSCL/QR4C/FCT, opth=0/optv=1).
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, mesh%nod2D))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, mesh%nod2D))
        tracers%data(tr_num)%AB_order   = 2
        tracers%data(tr_num)%tra_adv_hor = 'MUSCL'
        tracers%data(tr_num)%tra_adv_ver = 'QR4C'
        tracers%data(tr_num)%tra_adv_lim = 'FCT'
        tracers%data(tr_num)%tra_adv_ph  = 0.0_WP
        tracers%data(tr_num)%tra_adv_pv  = 1.0_WP
        tracers%data(tr_num)%i_vert_diff = .true.
        tracers%data(tr_num)%valuesAB        = 0.0_WP
        tracers%data(tr_num)%valuesold(1,:,:) = tracers%data(tr_num)%values
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
    call muscl_adv_init(tracers%work, mesh)

    !===========================================================================
    ! reduced-M2 module config (= the FESOM2 oracle / CORE2 namelist).
    alpha = 1.0_WP; theta = 1.0_WP
    N2smth_h     = .true.
    ! &oce_dyn scale_area: the CORE2 namelist.oce value, NOT the o_PARAM default
    ! 2.0e8. read_param_phys() is never called here, so the default used to stand.
    ! The ice FCT is its only consumer (diff = ice_diff*sqrt(elem_area/scale_area),
    ! mod_ice_fct.F90); inert while ice_diff=0, but at ice_diff/=0 the default gives
    ! sqrt(5.8e9/2e8) = 5.4x the reference ice diffusion -- the cause of the C port's
    ! systematic bi-hemispheric ice loss and Southern Ocean warming drift over 63-yr
    ! hindcasts (port2 21da79e).
    scale_area   = 5.8e9_WP
    mix_scheme_nmb = 2           ! reduced-M2 mixing = PP (KPP is M5c)
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.
    use_windmix   = .false.

    !===========================================================================
    ! M4e GM/Redi (FESOM3_FER_GM / FESOM3_REDI): enable Gent-McWilliams bolus advection
    ! (+ Redi isopycnal diffusion) in the FORCED / fully-native lifecycle. Identical block to
    ! fesom_lifecycle.F90 (M4c/M4d) — force the work_core GM+Redi config, allocate + init the
    ! GM/Redi work arrays to the FESOM2 oce_setup_step.F90:962-967 values (fer_K=500/fer_c=1/
    ! fer_scal=0/fer_gamma=0, fer_uv/fer_w=0; Ki=0; fer_tapfac=1) so below-bottom byte-matches.
    ! The GM+Redi chain + bolus run inside step_oce / solve_tracers_ale (guarded by Fer_GM/Redi);
    ! the native fluxes feed the surface BC so this is the production tracer physics under forcing.
    if (use_fer_gm .or. use_redi) then
        Fer_GM = .true.        ! Redi rides the GM coupling (Ki = max(fer_scal*Redi_Kmax, K_GM_min))
        K_GM_max = 1000.0_WP; K_GM_min = 2.0_WP; K_GM_bvref = 1
        K_GM_rampmax = -1.0_WP; K_GM_rampmin = -1.0_WP; K_GM_resscalorder = 2.0_WP
        K_GM_cm = 3.0_WP; K_GM_cmin = 0.1_WP; K_GM_Ktaper = .false.
        scaling_Ferreira = .false.; scaling_Rossby = .false.; scaling_resolution = .true.
        scaling_FESOM14 = .false.; scaling_GMzexp = .true.; scaling_GINsea = .false.
        GMzexp_zref = 500.0_WP; GMzexp_smin = 0.6_WP
        allocate(dyn%fer_uv(2, nl-1, mesh%elem2D), dyn%fer_w(nl, mesh%nod2D))
        allocate(dyn%work%sw_alpha(nl-1, mesh%nod2D), dyn%work%sw_beta(nl-1, mesh%nod2D))
        allocate(dyn%work%sigma_xy(2, nl-1, mesh%nod2D))
        allocate(dyn%work%fer_K(nl, mesh%nod2D), dyn%work%fer_c(mesh%nod2D), dyn%work%fer_scal(mesh%nod2D))
        allocate(dyn%work%fer_gamma(2, nl, mesh%nod2D))
        dyn%fer_uv = 0.0_WP; dyn%fer_w = 0.0_WP
        dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP; dyn%work%sigma_xy = 0.0_WP
        dyn%work%fer_K = 500.0_WP; dyn%work%fer_c = 1.0_WP; dyn%work%fer_scal = 0.0_WP
        dyn%work%fer_gamma = 0.0_WP
        if (use_redi) then
            ! M4d Redi (work_core): Redi=T; Redi_Kmax<=0 -> synced to K_GM_max; Redi_Kmin=100;
            ! Redi_Ktaper sqrt-split; ODM95 slope taper (work_core ODM95_Scr=0.2e-2); K_hor=0.
            Redi = .true.; Redi_Ktaper = .true.; Redi_Kmax = 0.0_WP; Redi_Kmin = 100.0_WP
            scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP; scaling_LDD97 = .false.
            K_hor = 0.0_WP
            allocate(dyn%work%Ki(nl-1, mesh%nod2D), dyn%work%fer_tapfac(nl-1, mesh%nod2D))
            allocate(dyn%work%neutral_slope(3, nl-1, mesh%nod2D), dyn%work%slope_tapered(3, nl-1, mesh%nod2D))
            allocate(tracers%work%tr_z(nl, mesh%nod2D))
            dyn%work%Ki = 0.0_WP; dyn%work%fer_tapfac = 1.0_WP
            dyn%work%neutral_slope = 0.0_WP; dyn%work%slope_tapered = 0.0_WP
            tracers%work%tr_z = 0.0_WP
            write(*,'(a)') 'fesom_lifecycle_native: Fer_GM + Redi ENABLED (work_core GM+Redi config)'
        else
            Redi = .false.
            write(*,'(a)') 'fesom_lifecycle_native: Fer_GM ENABLED (work_core GM config; Redi off)'
        end if
    end if

    !===========================================================================
    ! M5c KPP vertical mixing (FESOM3_MIX_KPP) — the production CORE2 boundary-layer scheme,
    ! identical to the M5b unforced fesom_lifecycle block. mix_scheme_nmb=1 routes step_oce to
    ! oce_mixing_kpp_driver (Av element + Kv_double node), then Kv=Kv_double(:,:,1) + mo_convect;
    ! Av stays element-based so impl_vert_visc_ale is UNCHANGED. The KPP arrays live in dyn%work
    ! (allocated only here). ref_sss/ref_sss_local mirror work_core (used only by the gated-off
    ! ghats salinity branch). sw_alpha/sw_beta (Bo) are allocated by the GM/Redi block above.
    if (use_kpp) then
        mix_scheme_nmb = 1
        Ricr = 0.3_WP; concv = 1.6_WP
        visc_sh_limit = 5.0e-3_WP; diff_sh_limit = 5.0e-3_WP
        ref_sss = 34.0_WP; ref_sss_local = .true.
        if (.not. allocated(dyn%work%sw_alpha)) then
            allocate(dyn%work%sw_alpha(nl-1, mesh%nod2D), dyn%work%sw_beta(nl-1, mesh%nod2D))
            dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP
        end if
        allocate(dyn%work%Kv_double(nl, mesh%nod2D, tracers%num_tracers))
        allocate(dyn%work%viscA_kpp(nl, mesh%nod2D), dyn%work%blmc(nl, mesh%nod2D, 3))
        allocate(dyn%work%ghats(nl-1, mesh%nod2D), dyn%work%dkm1(mesh%nod2D, 3))
        allocate(dyn%work%dbsfc(nl, mesh%nod2D), dyn%work%dVsq(nl, mesh%nod2D))
        allocate(dyn%work%sw_3d(nl, mesh%nod2D))
        allocate(dyn%work%hbl(mesh%nod2D), dyn%work%bfsfc(mesh%nod2D))
        allocate(dyn%work%stable(mesh%nod2D), dyn%work%caseA(mesh%nod2D))
        allocate(dyn%work%ustar(mesh%nod2D), dyn%work%Bo(mesh%nod2D), dyn%work%kbl(mesh%nod2D))
        dyn%work%Kv_double = 0.0_WP; dyn%work%viscA_kpp = 0.0_WP; dyn%work%blmc = 0.0_WP
        dyn%work%ghats = 0.0_WP; dyn%work%dkm1 = 0.0_WP; dyn%work%dbsfc = 0.0_WP
        dyn%work%dVsq = 0.0_WP; dyn%work%sw_3d = 0.0_WP
        dyn%work%hbl = 0.0_WP; dyn%work%bfsfc = 0.0_WP; dyn%work%stable = 0.0_WP
        dyn%work%caseA = 0.0_WP; dyn%work%ustar = 0.0_WP; dyn%work%Bo = 0.0_WP; dyn%work%kbl = 0
        call oce_mixing_kpp_init(Ricr, concv)   ! wmt/wst lookup tables + Vtc/cg (once)
        write(*,'(a)') 'fesom_lifecycle_native: KPP vertical mixing ENABLED (work_core KPP)'
    end if

    !===========================================================================
    ! M7c TKE vertical mixing (FESOM3_MIX_TKE) — the prognostic cvmix_TKE Av/Kv producer, the
    ! FORCED/native analog of the M7b unforced block. mix_scheme_nmb=5 routes step_oce to
    ! calc_cvmix_tke (assemble vshear2/bvfreq2/dz_trr + the REAL forc_tke_surf=|stress_node_surf|/
    ! density_0 -> integrate_tke -> Kv=tke_Kv + Av=avg(tke_Av)) + mo_convect. stress_node_surf is
    ! the real node stress from oce_fluxes_mom (allocated by the forcing path). ⚠️ work_*_tke has
    ! use_sw_pene=.true. -> the runner MUST set FESOM3_SW_PENE=1 (the sw_3d tracer term is allocated
    ! by the do_swpene block below, NOT here). tke restart-write is M8 (tke=0 at cold start).
    if (use_tke) then
        mix_scheme_nmb = 5
        allocate(dyn%work%tke(nl, mesh%nod2D), dyn%work%tke_Av(nl, mesh%nod2D), &
                 dyn%work%tke_Kv(nl, mesh%nod2D))
        dyn%work%tke = 0.0_WP; dyn%work%tke_Av = 0.0_WP; dyn%work%tke_Kv = 0.0_WP
        call tke_init(tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, tke_kappaM_min, &
                      tke_kappaM_max, tke_min, tke_surf_min, tke_mxl_choice, tke_only, &
                      tke_use_ubound_dirichlet, tke_use_lbound_dirichlet, tke_dolangmuir)
        write(*,'(a)') 'fesom_lifecycle_native: TKE vertical mixing ENABLED (cvmix_TKE)'
    end if

    !===========================================================================
    ! M5c shortwave penetration (FESOM3_SW_PENE) — flip use_sw_pene + allocate the constant
    ! chlorophyll (work_core chl_data_source='None', chl_const=0.1). cal_shortwave_rad runs
    ! after oce_fluxes each step, filling dyn%work%sw_3d (consumed by the tracer TDMA) and
    ! adding the visible band back to atm%heat_flux. sw_3d is allocated by the KPP block.
    if (do_swpene) then
        use_sw_pene = .true.
        allocate(chl(mesh%nod2D)); chl = 0.1_WP
        if (.not. allocated(dyn%work%sw_3d)) then
            allocate(dyn%work%sw_3d(nl, mesh%nod2D)); dyn%work%sw_3d = 0.0_WP
        end if
        write(*,'(a)') 'fesom_lifecycle_native: shortwave penetration ENABLED (chl_const=0.1)'
    end if

    !===========================================================================
    ! M5c KPP nonlocal counter-gradient flux (FESOM3_KPP_NONLCL) — DEAD in production
    ! (work_core leaves use_kpp_nonlclflx=.false.); the KPP_NONLCL=1 gate variant turns it on
    ! BOTH sides to byte-verify the otherwise-untested ghats term in diff_ver_part_impl_ale.
    if (do_nonlcl) then
        use_kpp_nonlclflx = .true.
        write(*,'(a)') 'fesom_lifecycle_native: KPP nonlocal counter-gradient flux ENABLED (ghats)'
    end if

    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt)

    !===========================================================================
    ! ice setup (allocate + ice_mass_matrix_fill + cold-start ice_initial_state) + the CORE2
    ! &ice_dyn / &ice_therm config overrides (the M3c/M3d/M3e namelist doubles), identical to
    ! fesom_icefluxdump.
    call ice_setup(ice, tracers, mesh, dt)
    ice%whichEVP      = whichevp
    ice%cd_oce_ice    = 0.0055_WP
    ice%delta_min     = 1.0e-11_WP
    ice%ice_diff      = 0.0_WP
    ice%ice_gamma_fct = 0.5_WP
    ice%thermo%Sice              = 4.0_WP
    ice%thermo%iclasses          = 7
    ice%thermo%new_iclasses      = .false.
    ice%thermo%h_cutoff          = 3.0_WP
    ice%thermo%h0                = 0.5_WP
    ice%thermo%h0_s              = 0.5_WP
    ice%thermo%hmin              = 0.01_WP
    ice%thermo%Armin             = 0.01_WP
    ice%thermo%emiss_ice         = 0.97_WP
    ice%thermo%emiss_wat         = 0.97_WP
    ice%thermo%albsn             = 0.81_WP
    ice%thermo%albsnm            = 0.77_WP
    ice%thermo%albi              = 0.7_WP
    ice%thermo%albim             = 0.68_WP
    ice%thermo%albw              = 0.1_WP
    ice%thermo%open_water_albedo = 0
    ice%thermo%con               = 2.1656_WP
    ice%thermo%consn             = 0.31_WP
    ice%thermo%snowdist          = .true.
    ice%thermo%c_melt            = 0.5_WP

    !===========================================================================
    ! atmflux: allocate (nod2D) + scalar config (the CORE2 namelist scalars).
    call alloc_atm(atm, mesh%nod2D)
    atm%Ch_atm_ice    = 0.00175_WP
    atm%Ce_atm_ice    = 0.00175_WP
    atm%ref_sss       = 34.0_WP
    atm%ref_sss_local = .true.
    atm%use_virt_salt = (trim(which_ALE)=='linfs')   ! M6a-3: zstar/zlevel -> real freshwater flux
    atm%l_snow        = .true.
    atm%surf_relax_S  = 1.929e-06_WP

    !===========================================================================
    ! M3f-3 native CORE2 forcing read setup (optional, FESOM3_FORCING_DIR). The 8 NCAR
    ! fields are read + bilinear-interpolated + g2r-rotated like M2.10a (proven on these
    ! exact files); the per-step rdate advances like FESOM2 sbc_do (clock 0 1 1948, dt=1800).
    call get_environment_variable('FESOM3_FORCING_DIR', forcing_dir)
    native_forcing = (len_trim(forcing_dir) > 0)
    if (.not. prescribe_atm .and. .not. native_forcing) then
        write(*,'(a)') 'fesom_lifecycle_native: FULLY NATIVE mode (no FESOM3_ATMFLUX_FILE) &
                       &requires FESOM3_FORCING_DIR'
        error stop 1
    end if
    if (native_forcing) then
        fyear = 1948
        frc%nfld = 8;  frc%nnod = mesh%nod2D
        frc%iyear = 1948; frc%imm = 1; frc%idd = 1; frc%freq = 1; frc%tmid = 1
        frc%ic_cyclic = .true.; frc%rotated_grid = .true.
        frc%i_xwind = 1; frc%i_ywind = 2
        frc%f(1)%file_base = trim(forcing_dir)//'/u_10.';        frc%f(1)%varname = 'U_10_MOD'
        frc%f(2)%file_base = trim(forcing_dir)//'/v_10.';        frc%f(2)%varname = 'V_10_MOD'
        frc%f(3)%file_base = trim(forcing_dir)//'/q_10.';        frc%f(3)%varname = 'Q_10_MOD'
        frc%f(4)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(4)%varname = 'SWDN_MOD'
        frc%f(5)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(5)%varname = 'LWDN_MOD'
        frc%f(6)%file_base = trim(forcing_dir)//'/t_10.';        frc%f(6)%varname = 'T_10_MOD'
        frc%f(7)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(7)%varname = 'RAIN'
        frc%f(8)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(8)%varname = 'SNOW'
        rdate_cold = real(forcing_julday(fyear,1,1,'noleap'),WP)
        call forcing_alloc(frc)
        do fld = 1, frc%nfld
            call forcing_read_grid(frc, fld, fyear)
            call forcing_build_bilin(frc, fld, mesh, partit)
        end do
        do fld = 1, frc%nfld
            call forcing_getcoeffld(frc, fld, fyear, rdate_cold, mesh, partit)
        end do
        call forcing_rotate_wind(frc, mesh, partit)
        allocate(nuw(mesh%nod2D), nvw(mesh%nod2D), nta(mesh%nod2D), nsh(mesh%nod2D), &
                 nswr(mesh%nod2D), nlw(mesh%nod2D), npr(mesh%nod2D), nps(mesh%nod2D), &
                 ncd(mesh%nod2D), nch(mesh%nod2D), nce(mesh%nod2D), &
                 nsx(mesh%nod2D), nsy(mesh%nod2D), nix(mesh%nod2D), niy(mesh%nod2D))
        ! M3f-3b native runoff + SSS climatology (read ONCE — both constant for a Jan run).
        ! runoff: read_other_NetCDF('Foxx_o_roff', rec 1, check_dummy=.false. -> missing=0) then
        ! /1000 (kg/s/m^2 -> m/s); Ssurf: read_other_NetCDF('SALT', rec month=1, check_dummy=.true.
        ! -> 30-neighbour fill). Both interpolate to vertices (do_onvert=.true.).
        call get_environment_variable('FESOM3_RUNOFF_FILE', runoff_file)
        if (len_trim(runoff_file) == 0) &
            runoff_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc'
        call get_environment_variable('FESOM3_SSS_FILE', sss_file)
        if (len_trim(sss_file) == 0) &
            sss_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc'
        allocate(nro(mesh%nod2D), nss(mesh%nod2D))
        ! 1-rank: omit partit (like ocean2ice/ice_timestep) -> num = mesh%nod2D. M3f-4 passes
        ! partit (with halo dims) so the interp also fills the eDim halo nodes.
        call read_other_NetCDF(trim(runoff_file), 'Foxx_o_roff', 1, nro, .false., .true., mesh)
        nro = nro / 1000.0_WP
        call read_other_NetCDF(trim(sss_file), 'SALT', 1, nss, .true., .true., mesh)
        write(*,'(a,a)') 'fesom_lifecycle_native: NATIVE forcing self-check ON, dir=', trim(forcing_dir)
        write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle_native: native runoff[', &
            minval(nro), maxval(nro), '] Ssurf[', minval(nss), maxval(nss), ']'
    end if

    !===========================================================================
    ! ocean-step surface BC arrays. heat_flux/water_flux/virtual_salt/relax_salt/stress_surf
    ! are now produced NATIVELY each step (oce_fluxes / oce_fluxes_mom write into atm%* and
    ! the local stress_surf); Ki/real_salt_flux stay 0 (M4). M6a-3: is_nonlinfs=1 for
    ! non-linfs (zstar) so bc_surface uses the real-salt / advective-heat freshwater terms.
    is_nonlinfs = merge(1.0_WP, 0.0_WP, trim(which_ALE)/='linfs')
    allocate(Ki(nl-1, mesh%nod2D), real_salt_flux(mesh%nod2D), stress_surf(2, mesh%elem2D))
    Ki = 0.0_WP; real_salt_flux = 0.0_WP; stress_surf = 0.0_WP

    !===========================================================================
    ! per-substep node dump (mod_dump).
    call init_dump_identity(partit%mype, mesh%nod2D, mesh%elem2D)

    !===========================================================================
    ! open the per-step prescribed atmospheric forcing dump (16 nod2D arrays/step), OR
    ! announce FULLY NATIVE mode (the whole atmosphere computed natively each step).
    if (prescribe_atm) then
        open(newunit=atm_unit, file=trim(atm_file), status='old', form='unformatted', &
             access='stream', action='read')
        write(*,'(a,a)') 'fesom_lifecycle_native: prescribing per-step atm forcing from ', trim(atm_file)
    else
        write(*,'(a)') 'fesom_lifecycle_native: FULLY NATIVE atmosphere (NCAR read + bulk + &
                       &stresses + runoff + Ssurf computed each step)'
    end if

    ! optional native-vs-oracle flux self-check.
    call get_environment_variable('FESOM3_FLUX_FILE', flux_file)
    selfcheck = (len_trim(flux_file) > 0)
    if (selfcheck) then
        open(newunit=flux_unit, file=trim(flux_file), status='old', form='unformatted', &
             access='stream', action='read')
        allocate(o_hf(mesh%nod2D), o_wf(mesh%nod2D), o_vs(mesh%nod2D), o_rs(mesh%nod2D), &
                 o_ss(2, mesh%elem2D))
        write(*,'(a)') 'fesom_lifecycle_native: per-step native-vs-oracle flux self-check ON'
    end if

    !===========================================================================
    ! runloop: ocean2ice -> ice_timestep -> oce_fluxes_mom -> oce_fluxes -> step_oce.
    ! This is the FESOM2 runloop order (ocean2ice 689 -> update_atm_forcing 700 ->
    ! ice_timestep 719 -> oce_fluxes_mom 727 -> oce_fluxes 728 -> oce_timestep_ale 750),
    ! with update_atm_forcing replaced by the prescribed read.
    do n = 1, nsteps
        if (prescribe_atm) call read_atm_record(atm_unit, atm, ice)

        call ocean2ice(ice, dyn, tracers, mesh)
        ! Native atmosphere (after ocean2ice -> srfoce is live; before EVP -> ice%uice/vice
        ! are the previous step's, as FESOM2 update_atm_forcing uses for the wind-on-ice
        ! stress). Two modes:
        !   prescribe_atm  : recompute natively + SELF-CHECK vs the prescribed atm (the step
        !                    still USES the prescribed values — the proven M3f-2 path).
        !   fully native   : compute natively + WRITE into atm%* / ice%stress_atmice (the step
        !                    USES the native values — M3f-3 complete).
        if (native_forcing) then
            if (prescribe_atm) then
                call forcing_selfcheck(n)
            else
                call apply_native_forcing(n)
            end if
        end if

        call ice_timestep(ice, mesh, atm)
        call oce_fluxes_mom(ice, atm, stress_surf, mesh)
        call oce_fluxes(ice, tracers, atm, mesh)
        ! M5c: shortwave penetration — the LAST step of FESOM2 oce_fluxes (ice_oce_coupling.F90:873,
        ! after the whole air-sea budget). Fills dyn%work%sw_3d + adds the visible band back to
        ! atm%heat_flux; matches the oracle flux dump (recorded after oce_fluxes) so the self-check
        ! and the tracer TDMA both see the modified heat_flux. albw = ice%thermo%albw (work_core 0.1).
        if (use_sw_pene) call cal_shortwave_rad(.true., ice%thermo%albw, atm%shortwave, chl, &
                              ice%data(1)%values(1:mesh%nod2D), atm%heat_flux, dyn%work%sw_3d, mesh)

        if (selfcheck) call flux_selfcheck(flux_unit, n, atm, stress_surf, mesh, &
                                           o_hf, o_wf, o_vs, o_rs, o_ss)

        ! M5c: stress_node_surf=atm%stress_node_surf (oce_fluxes_mom output) feeds KPP ustar;
        ! step_oce reads it only when mix_scheme_nmb==1 (PP path ignores it).
        ! M6a-3: pass the NATIVE atm%real_salt_flux (obudget output, the ice-brine salt flux,
        ! nonzero for zstar use_virt_salt=.false.) — NOT the local zeroed real_salt_flux. linfs
        ! ignores it (is_nonlinfs=0), so this is byte-neutral for the M3/M4/M5 linfs gates.
        call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                      atm%heat_flux, atm%water_flux, atm%virtual_salt, atm%relax_salt, &
                      atm%real_salt_flux, is_nonlinfs, stress_surf, stress_node_surf=atm%stress_node_surf)
        write(*,'(a,i0,a,es12.4,a,es12.4,a,es12.4)') 'fesom_lifecycle_native: step ', n, &
            '  max|eta_n|=', maxval(abs(dyn%eta_n)), '  max|uv|=', maxval(abs(dyn%uv)), &
            '  max|a_ice|=', maxval(abs(ice%data(1)%values(1:mesh%nod2D)))
    end do
    if (prescribe_atm) close(atm_unit)
    if (selfcheck) close(flux_unit)

    call dump_finalize()
    write(*,'(a,i0,a)') 'fesom_lifecycle_native: done (', nsteps, ' steps, NATIVE fluxes).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    ! per-substep node dump init with identity global ids (1-rank).
    subroutine init_dump_identity(mype, nn, ne)
        integer, intent(in) :: mype, nn, ne
        integer, allocatable :: idlist(:)
        integer :: i
        allocate(idlist(max(nn, ne)))
        do i = 1, size(idlist); idlist(i) = i; end do
        call dump_init(mype, nn, idlist(1:nn), ne, idlist(1:ne))
    end subroutine init_dump_identity

    ! read one per-step atm forcing record (the 16 nod2D arrays, the fesom_atmflux_dump order).
    subroutine read_atm_record(unit, a, ic)
        integer,         intent(in)    :: unit
        type(t_atmflux), intent(inout) :: a
        type(t_ice),     intent(inout) :: ic
        integer(int32) :: s, nn, ne
        read(unit) s, nn, ne
        read(unit) a%shortwave
        read(unit) a%longwave
        read(unit) a%Tair
        read(unit) a%shum
        read(unit) a%prec_rain
        read(unit) a%prec_snow
        read(unit) a%runoff
        read(unit) a%u_wind
        read(unit) a%v_wind
        read(unit) a%Ch_atm_oce_arr
        read(unit) a%Ce_atm_oce_arr
        read(unit) a%stress_atmoce_x
        read(unit) a%stress_atmoce_y
        read(unit) ic%stress_atmice_x
        read(unit) ic%stress_atmice_y
        read(unit) a%Ssurf
    end subroutine read_atm_record

    ! compare the NATIVE fluxes against the oracle flux dump (fesom_flux_dump order) and
    ! print max|delta| per field — localizes a coupling/flux mismatch before the full gate.
    subroutine flux_selfcheck(unit, n, a, ss, m, hf, wf, vs, rs, oss)
        integer,         intent(in)    :: unit, n
        type(t_atmflux), intent(in)    :: a
        real(kind=WP),   intent(in)    :: ss(:,:)
        type(t_mesh),    intent(in)    :: m
        real(kind=WP),   intent(inout) :: hf(:), wf(:), vs(:), rs(:), oss(:,:)
        integer(int32) :: s, nn, ne
        read(unit) s, nn, ne
        read(unit) hf
        read(unit) wf
        read(unit) vs
        read(unit) rs
        read(unit) oss
        write(*,'(a,i0,a,5es11.3)') 'fesom_lifecycle_native: [selfcheck step ', n, &
            '] max|d(hf,wf,vs,rs,ss)| = ', &
            maxval(abs(a%heat_flux(1:m%nod2D)    - hf)), &
            maxval(abs(a%water_flux(1:m%nod2D)   - wf)), &
            maxval(abs(a%virtual_salt(1:m%nod2D) - vs)), &
            maxval(abs(a%relax_salt(1:m%nod2D)   - rs)), &
            maxval(abs(ss - oss))
    end subroutine flux_selfcheck

    ! M3f-3 compute the native CORE2 atmosphere into the n* work arrays: the 8 NCAR fields
    ! (timeinterp at the per-step rdate) + the NCAR bulk (Ch/Ce, using the live srfoce) +
    ! stress_atmoce + the wind-on-ice stress (using the previous-step uice/vice). runoff (nro)
    ! + Ssurf (nss) were read once at setup (constant for a January run).
    subroutine compute_native_forcing(n)
        integer, intent(in) :: n
        integer :: nn
        real(kind=WP) :: tnew, rcur
        nn = mesh%nod2D
        tnew = real(n,WP)*dt                          ! daynew=1 for n<48 (within forcing day 1)
        rcur = real(forcing_julday(fyear,1,1,'noleap'),WP) + real(1-1,WP) &
             + tnew/86400._WP - dt/86400._WP/2._WP
        call forcing_timeinterp(frc, rcur, partit)
        nuw = frc%atmdata(1,1:nn);  nvw = frc%atmdata(2,1:nn);  nsh = frc%atmdata(3,1:nn)
        nswr = frc%atmdata(4,1:nn);  nlw = frc%atmdata(5,1:nn)
        nta = frc%atmdata(6,1:nn) - 273.15_WP
        npr = frc%atmdata(7,1:nn) / 1000._WP
        nps = frc%atmdata(8,1:nn) / 1000._WP
        ncd = 0.0_WP; nch = 0.0_WP; nce = 0.0_WP
        call forcing_bulk_ncar(10.0_WP, 10.0_WP, 10.0_WP, nta, nsh, nuw, nvw, &
                               ice%srfoce_temp, ice%srfoce_u, ice%srfoce_v, ncd, nch, nce, mesh)
        call forcing_wind_stress(0.0_WP, nuw, nvw, ice%srfoce_u, ice%srfoce_v, ncd, nsx, nsy, mesh)
        call forcing_ice_stress(0.0012_WP, nuw, nvw, ice%uice, ice%vice, nix, niy, mesh)
    end subroutine compute_native_forcing

    ! prescribe + self-check mode: compute natively, print max|delta| vs the prescribed atm.
    ! 0 => the native CORE2 forcing read + bulk + stress + runoff/SSS are byte-exact.
    subroutine forcing_selfcheck(n)
        integer, intent(in) :: n
        integer :: nn
        nn = mesh%nod2D
        call compute_native_forcing(n)
        write(*,'(a,i0,a,8es10.2)') 'fesom_lifecycle_native: [forcing step ', n, &
            '] max|d NCAR(sw,lw,Ta,sh,rn,sn,uw,vw)| = ', &
            maxval(abs(nswr-atm%shortwave)),  maxval(abs(nlw-atm%longwave)),  &
            maxval(abs(nta-atm%Tair)),       maxval(abs(nsh-atm%shum)),      &
            maxval(abs(npr-atm%prec_rain)),  maxval(abs(nps-atm%prec_snow)), &
            maxval(abs(nuw-atm%u_wind)),     maxval(abs(nvw-atm%v_wind))
        write(*,'(a,i0,a,6es10.2)') 'fesom_lifecycle_native: [forcing step ', n, &
            '] max|d bulk(Ch,Ce,sox,soy,six,siy)| = ', &
            maxval(abs(nch-atm%Ch_atm_oce_arr)), maxval(abs(nce-atm%Ce_atm_oce_arr)), &
            maxval(abs(nsx-atm%stress_atmoce_x)), maxval(abs(nsy-atm%stress_atmoce_y)), &
            maxval(abs(nix-ice%stress_atmice_x(1:nn))), maxval(abs(niy-ice%stress_atmice_y(1:nn)))
        ! M3f-3b: native runoff + SSS climatology vs the prescribed atm (both should be 0).
        write(*,'(a,i0,a,2es10.2)') 'fesom_lifecycle_native: [forcing step ', n, &
            '] max|d clim(runoff,Ssurf)| = ', &
            maxval(abs(nro-atm%runoff)), maxval(abs(nss-atm%Ssurf))
    end subroutine forcing_selfcheck

    ! fully-native mode: compute natively, WRITE into atm%* / ice%stress_atmice for the step.
    subroutine apply_native_forcing(n)
        integer, intent(in) :: n
        integer :: nn
        nn = mesh%nod2D
        call compute_native_forcing(n)
        atm%shortwave       = nswr
        atm%longwave        = nlw
        atm%Tair            = nta
        atm%shum            = nsh
        atm%prec_rain       = npr
        atm%prec_snow       = nps
        atm%u_wind          = nuw
        atm%v_wind          = nvw
        atm%Ch_atm_oce_arr  = nch
        atm%Ce_atm_oce_arr  = nce
        atm%stress_atmoce_x = nsx
        atm%stress_atmoce_y = nsy
        ice%stress_atmice_x(1:nn) = nix
        ice%stress_atmice_y(1:nn) = niy
        atm%runoff          = nro
        atm%Ssurf           = nss
    end subroutine apply_native_forcing

    subroutine alloc_atm(a, nn)
        type(t_atmflux), intent(inout) :: a
        integer,         intent(in)    :: nn
        allocate(a%shortwave(nn), a%longwave(nn), a%Tair(nn), a%shum(nn))
        allocate(a%prec_rain(nn), a%prec_snow(nn), a%runoff(nn))
        allocate(a%u_wind(nn), a%v_wind(nn), a%Ch_atm_oce_arr(nn), a%Ce_atm_oce_arr(nn))
        allocate(a%evaporation(nn), a%ice_sublimation(nn), a%flice(nn), a%real_salt_flux(nn))
        allocate(a%fw_ice(nn), a%fw_snw(nn))
        allocate(a%hf_Qlat(nn), a%hf_Qsen(nn), a%hf_Qradtot(nn))
        allocate(a%hf_Qswr(nn), a%hf_Qlwr(nn), a%hf_Qlwrout(nn))
        allocate(a%heat_flux(nn), a%water_flux(nn), a%heat_flux_in(nn))
        allocate(a%virtual_salt(nn), a%relax_salt(nn))
        allocate(a%stress_node_surf(2, nn))
        allocate(a%stress_atmoce_x(nn), a%stress_atmoce_y(nn), a%Ssurf(nn))
        a%shortwave = 0.0_WP; a%longwave = 0.0_WP; a%Tair = 0.0_WP; a%shum = 0.0_WP
        a%prec_rain = 0.0_WP; a%prec_snow = 0.0_WP; a%runoff = 0.0_WP
        a%u_wind = 0.0_WP; a%v_wind = 0.0_WP; a%Ch_atm_oce_arr = 0.0_WP; a%Ce_atm_oce_arr = 0.0_WP
        a%evaporation = 0.0_WP; a%ice_sublimation = 0.0_WP; a%flice = 0.0_WP; a%real_salt_flux = 0.0_WP
        a%fw_ice = 0.0_WP; a%fw_snw = 0.0_WP
        a%hf_Qlat = 0.0_WP; a%hf_Qsen = 0.0_WP; a%hf_Qradtot = 0.0_WP
        a%hf_Qswr = 0.0_WP; a%hf_Qlwr = 0.0_WP; a%hf_Qlwrout = 0.0_WP
        a%heat_flux = 0.0_WP; a%water_flux = 0.0_WP; a%heat_flux_in = 0.0_WP
        a%virtual_salt = 0.0_WP; a%relax_salt = 0.0_WP
        a%stress_node_surf = 0.0_WP
        a%stress_atmoce_x = 0.0_WP; a%stress_atmoce_y = 0.0_WP; a%Ssurf = 0.0_WP
    end subroutine alloc_atm

end program fesom_lifecycle_native
