program fesom_lifecycle_native_mr
    ! M3f-4: the MULTI-RANK native-flux forced lifecycle driver — the multi-rank analog of
    ! the 1-rank fesom_lifecycle_native (and the native-flux analog of fesom_lifecycle_mr).
    ! Runs on dist_<NP> (NP>=1) the REAL free-running coupled sea-ice + ocean lifecycle with
    ! the WHOLE air-sea forcing computed NATIVELY each step (NO prescribed atmosphere, NO
    ! prescribed fluxes):
    !
    !   native atmosphere  8 NCAR fields (mod_forcing_read) + NCAR bulk Ch/Ce + stress_atmoce
    !                      + wind-on-ice stress (mod_forcing_bulk) + runoff/Ssurf monthly
    !                      climatology (mod_forcing_other), ALL over owned+halo (nNodL)
    !   ocean2ice          ocean surface -> srfoce_u/v/temp/salt/ssh                  (M3b)
    !   ice_timestep       EVP -> FCT -> cut_off -> thermodynamics                     (M3f-1)
    !   oce_fluxes_mom     ice+atm momentum stress -> stress_surf                      (M3e)
    !   oce_fluxes         heat/freshwater/salt budget -> the tracer/SSH BCs           (M3e)
    !   step_oce           the whole ocean step on the NATIVE fluxes (through partit)  (M2.12)
    !
    ! Each kernel runs the M2.12 optional-partit path: at npes==1 the proven 1-rank code
    ! VERBATIM, at npes>1 owned/halo loop bounds + the FESOM2 halo exchanges. The whole atm
    ! forcing + ice + budget chain is partition-independent per-node (computed over owned+halo
    ! from halo-valid inputs, no extra exchange); ocean2ice exchanges u_w/v_w, the EVP exchanges
    ! uice/vice, the FCT exchanges its solve intermediates, and oce_fluxes' integrate_nod /
    ! mesh%ocean_area do the cross-rank allreduce_sum. step_oce is the byte-proven M2.12c MR step.
    !
    ! The gate (tools/run_lifecycle_fullynative_gate_multirank.sh) compares the per-substep
    ! gid-keyed NODE dumps (mod_dump) PER RANK vs the REAL FESOM2 forced lifecycle on the SAME
    ! dist_<NP> partition (the L8 same-partition rule). NOTE: the oracle flux dump
    ! (fesom_flux_dump) is skipped at npes/=1, so there is no per-step flux self-check at
    ! multi-rank — the node-substep gate is the validation (the 1-rank fully-native gate's
    ! self-check already proved the native fluxes byte-exact).
    !
    !   FESOM3_MESH_DIR     mesh dir (must contain dist_<NP>/)  default: CORE2
    !   FESOM3_IC_FILE      IC netcdf                           default: pool phc3.0_winter.nc
    !   FESOM3_FORCING_DIR  native NCAR forcing dir             (REQUIRED — fully native)
    !   FESOM3_RUNOFF_FILE  runoff climatology                  default: pool CORE2_runoff.nc
    !   FESOM3_SSS_FILE     SSS restoring climatology           default: pool PHC2_salx.nc
    !   FESOM3_WHICHEVP     0=EVP / 1=mEVP                       (default 0)
    !   FESOM_DUMP_FILE     per-rank node dump prefix (mod_dump; <prefix>.<mype5>)
    !   FESOM_DUMP_MAXSTEPS dump step cap
    !   FESOM3_NSTEPS       number of steps                     (default: 3)
    !   FESOM3_SHARED_SYNC  1 = one owner per shared element: update_vel syncs uv (default);
    !                       0 = off (unsynced copies, A/B only). Read in src/infra/mod_shared_elem.F90.
    use mpi
    use mod_timer           ! permanent per-component wall-clock timing (TMR_* ids, timer_start/stop/report)
    use, intrinsic :: ieee_arithmetic   ! M8c Step-3 fix: control flush-to-zero (denormal) underflow mode
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use mod_precision,      only: WP, MP, MPI_WP
    use mod_constants,      only: density_0
    use mod_halo,           only: shared_copy_max_elem   ! FESOM3_STEP_DIAG shared-element spreads
    use mod_param_phys,     only: N2smth_h, alpha, theta, scale_area, shear_splines, N2_splines, rvo_upwind
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    ! M4f: GM bolus (FESOM3_FER_GM) + Redi isopycnal diffusion (FESOM3_REDI) at MULTI-RANK —
    ! the work_core GM+Redi config (the M4 routines already thread the optional partit + exchanges).
    use mod_param_phys,     only: Fer_GM, Redi, K_GM_max, K_GM_min, K_GM_bvref, &
                                  K_GM_rampmax, K_GM_rampmin, K_GM_resscalorder, K_GM_cm, &
                                  K_GM_cmin, K_GM_Ktaper, scaling_Ferreira, scaling_Rossby, &
                                  scaling_resolution, scaling_FESOM14, scaling_GMzexp, &
                                  scaling_GINsea, GMzexp_zref, GMzexp_smin
    use mod_param_phys,     only: Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                                  scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    ! M5d: KPP vertical mixing + shortwave penetration + (gated-off) ghats nonlocal flux at MR.
    use mod_param_phys,     only: Ricr, concv, visc_sh_limit, diff_sh_limit, &
                                  use_kpp_nonlclflx, ref_sss, ref_sss_local
    use mod_param_phys,     only: tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, &
                                  tke_kappaM_min, tke_kappaM_max, tke_min, tke_surf_min, &
                                  tke_mxl_choice, tke_only, tke_use_ubound_dirichlet, &
                                  tke_use_lbound_dirichlet, tke_dolangmuir
    use mod_config,         only: use_sw_pene, which_ALE, &
                                  cfg_dt => dt, step_per_day, run_length, run_length_unit, &
                                  runid, RestartInPath, RestartOutPath, include_fleapyear
    use mod_clock,          only: clock, clock_init, clock_nsteps, clock_finish, r_restart, &
                                  yearnew, yearold, daynew, timenew, month, &
                                  yearstart, ndpyr, day_in_month, fleapyear, num_day_in_month
    use oce_mixing_kpp,     only: oce_mixing_kpp_init
    use oce_mixing_tke,     only: tke_init
    use oce_wsplit,         only: wsplit_assert_params
    use oce_shortwave_pene, only: cal_shortwave_rad
    use mod_mesh,           only: t_mesh
    use mod_partit,         only: t_partit
    use mod_partitioning,   only: par_init, par_ex, set_partition
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
    use mod_forcing_read    ! native CORE2 atmospheric forcing read
    use mod_forcing_bulk,   only: forcing_bulk_ncar, forcing_wind_stress, forcing_ice_stress
    use mod_forcing_other,  only: read_other_NetCDF   ! runoff + SSS climatology read
    use oce_initial_state,  only: t_ic3d_config, do_ic3d
    use oce_muscl_adv,      only: muscl_adv_init
    use oce_ssh_rhs,        only: init_stiff_mat_ale
    use mod_step_oce,       only: step_oce
    use mod_io_meshdiag,    only: meshdiag_write
    use mod_io_means,       only: t_io_means, t_io_config, t_io_entry, t_means_clock, MEANS_MAXF, &
                                  means_init, means_read_namelist, means_define_node2d, &
                                  means_define_node3d, means_define_vector3d, means_define_elem2d, &
                                  means_define_elem3d, means_define_vector3d_elem, means_accumulate, &
                                  means_has, means_output, means_finalize
    use mod_dump,           only: dump_init, dump_finalize
    ! Restart Stage 5: checkpoint registry + register/write/read hooks (the lifecycle wiring of mod_io_restart).
    use mod_io_restart,     only: t_restart, restart_init, restart_register_state, restart_write, &
                                  restart_read, restart_finalize, restart_resolve_latest
    implicit none

    ! Model timestep, from step_per_day (env FESOM3_STEP_PER_DAY; default 48 => dt = 1800 s = 30 min,
    ! the CORE2 namelist timestep). Set once in the env block below, BEFORE init_stiff_mat_ale and
    ! ice_setup bake it in. Leaving the env unset reproduces the old parameter exactly, so every
    ! byte-gate is unchanged.
    integer                  :: spd, spd_env
    real(kind=WP)            :: dt

    character(len=512) :: mesh_dir, ic_file, env, whichevp_str
    type(t_partit)        :: partit
    type(t_mesh),  target :: mesh        ! target: restart registers live pointers into these (Stage 5)
    type(t_dyn),   target :: dyn
    type(t_tracer),target :: tracers
    type(t_ice),   target :: ice
    type(t_atmflux)    :: atm
    type(t_ic3d_config):: ic
    integer :: n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len, nsw, whichevp
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    ! M9 Stage 2: field output (mod_io_means), env-gated (FESOM3_OUTPUT) + configured via namelist.io
    ! (&nml_general knobs + &nml_list rows; Task 2.6). FESOM3_OUTPUT_EVERY is the no-namelist fallback.
    type(t_io_means)    :: oio
    type(t_io_config)   :: iocfg
    type(t_io_entry)    :: iolist(MEANS_MAXF)
    type(t_means_clock) :: oclk
    logical          :: do_output, have_nml
    integer          :: out_every, n_io, ii
    character(len=4096) :: out_env, out_dir_o, nml_io_path
    ! M8a: model clock (cold start from FESOM3_START_CLOCK) + run-length driver
    integer            :: clk_d0, clk_y0, clk_unit, idx, ierr
    real(kind=WP)      :: clk_t0, dstat(9)   ! dstat: M8c Step-3 global per-step physical diagnostics
    real(kind=WP)      :: dsum(2)            ! [sum ssh_rhs_old, sum|ssh_rhs_old|] over owned nodes
    logical            :: ftz_supported   ! M8c Step-3: flush-to-zero (denormal) underflow control (match FESOM2)
    logical            :: step_diag        ! perf: per-step global MPI_MAX stability diagnostic (env FESOM3_STEP_DIAG; default OFF — FESOM2 has no per-step collective)
    integer            :: mon_every        ! perf monitoring: emit a per-component timing report every mon_every steps (env FESOM3_TIMING_EVERY; 0 = only the final report)
    character(len=512) :: start_clock, restart_in
    ! Restart Stage 5 (Tasks 5.1/5.3): checkpoint registry + global cadence + restart-mode detection.
    type(t_restart)     :: rst
    logical             :: do_restart      ! FESOM3_RESTART set => register + write the checkpoint state
    logical             :: restart_mode    ! restart launch: do_restart AND a complete checkpoint resolves
    logical             :: wrote_final     ! did the last in-loop step already checkpoint? (skip the final)
    logical             :: rl_ok
    character(len=4096) :: restart_dir     ! FESOM3_RESTART value (no trailing slash); checkpoint folder root
    character(len=:), allocatable :: rl_folder
    integer             :: restart_length  ! cadence length  (&nml_restart restart_length; FESOM3_RESTART_LENGTH)
    character(len=8)    :: restart_length_unit  ! cadence unit (y|m|d|h|s|off; FESOM3_RESTART_UNIT)
    integer             :: cp_year, cp_day      ! normalized checkpoint clock (year-rollover, like clock_finish)
    real(kind=real64)   :: cp_tsec
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    logical :: use_fer_gm, use_redi   ! M4f: GM bolus / Redi isopycnal diffusion toggles
    logical :: use_kpp, do_swpene, do_nonlcl   ! M5d: KPP / shortwave pene / ghats nonlocal flux
    logical :: use_tke      ! M7d: FESOM3_MIX_TKE -> cvmix_TKE producer at multi-rank (LOCAL nNodL)
    integer :: momadv_env   ! FESOM3_MOMADV_OPT: 1 = vector invariant, 2 = scalar (default)
    logical       :: use_wsplit                            ! FESOM3_WSPLIT (presence)
    real(kind=WP) :: wsplit_mincfl_env, wsplit_maxcfl_env  ! FESOM3_WSPLIT_MINCFL/_MAXCFL (values)
    real(kind=WP), allocatable :: chl(:)       ! M5d const 0.1 / M8c Sweeney monthly climatology
    logical :: use_chl_sweeney                 ! M8c: FESOM3_CHL_SWEENEY -> read Sweeney chl monthly
    ! native atmospheric forcing read (the whole atmosphere, over owned+halo).
    character(len=512)  :: forcing_dir, runoff_file, sss_file, chl_file
    ! M8c: forcing dataset select. CORE2 = the legacy regression substitute (noleap, NCAR fields);
    ! JRA55 = the production target (gregorian + leap, JRA55-do fields, 3-hourly, tmid=0).
    character(len=16)   :: forc_set, forc_calendar
    type(t_atm_forcing) :: frc
    integer             :: fld, fyear
    real(kind=WP)       :: rdate_cold
    real(kind=WP), allocatable :: nuw(:), nvw(:), nta(:), nsh(:), nswr(:), nlw(:), npr(:), nps(:)
    real(kind=WP), allocatable :: ncd(:), nch(:), nce(:), nsx(:), nsy(:), nix(:), niy(:)
    real(kind=WP), allocatable :: nro(:), nss(:)
    integer, allocatable :: idlist(:)
    ! M8c bisection: optional per-field atmosphere self-check vs the F2 atmflux dump (np=1 only).
    character(len=512) :: atmchk_file
    integer            :: atmchk_unit
    logical            :: do_atm_chk
    real(kind=WP), allocatable :: tscr(:)

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
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
    ! Vertical cubic-spline reconstruction of the shear (oce_vert_spline) and of N^2
    ! (oce_pressure_bv), FESOM's RI_SPLINES analogue. Presence-based like the rest.
    ! Set BOTH for like-by-like Ri and neutral-slope ratios; one alone re-creates the
    ! numerator/denominator asymmetry the pair exists to remove.
    call get_environment_variable('FESOM3_SHEAR_SPLINES', env, length=env_len, status=ios)
    shear_splines = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_N2_SPLINES', env, length=env_len, status=ios)
    N2_splines = (ios == 0 .and. env_len > 0)
    ! Momentum advection operator. VALUE-based, NOT presence-based like the switches
    ! above: momadv_opt is 1 (vector invariant, oce_dyn_vinv) or 2 (scalar/flux form),
    ! so a bare presence test could not distinguish them. Unset => 2, the default.
    call get_environment_variable('FESOM3_MOMADV_OPT', env, length=env_len, status=ios)
    momadv_env = 2
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) momadv_env
    if (momadv_env /= 1 .and. momadv_env /= 2) &
        error stop 'FESOM3_MOMADV_OPT must be 1 or 2'
    ! Upwind blending of the face relative vorticity (oce_dyn_vinv; only active with
    ! momadv_opt==1). VALUE-based: a fraction in [0,1], not a presence switch.
    call get_environment_variable('FESOM3_RVO_UPWIND', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) rvo_upwind
    if (rvo_upwind < 0.0_WP .or. rvo_upwind > 1.0_WP) &
        error stop 'FESOM3_RVO_UPWIND must be in [0,1]'
    ! Smooth Courant-dependent w = w_e + w_i split (oce_wsplit, compute_Wvel_split): the
    ! part of w above the explicit CFL_z cap is advected by the implicit upstream solves
    ! (momentum TDMA, tracer vertical-diffusion TDMA, FCT adv_tra_vert_impl). Presence-
    ! based switch like FESOM3_FER_GM; the two parameters are VALUE-based like
    ! FESOM3_MOMADV_OPT. Defaults: the t_dyn cap 1.0 and an onset of 0.5*maxcfl, so setting
    ! only FESOM3_WSPLIT_MAXCFL keeps the type's onset/cap ratio. Validated once here
    ! (wsplit_assert_params: maxcfl > 0, 0 <= mincfl <= maxcfl). Off unless set
    ! (production precedent).
    ! An unparsable value is an error, not a silent fall-back to the default: a production
    ! run must never carry a cap other than the one the job file asked for.
    call get_environment_variable('FESOM3_WSPLIT', env, length=env_len, status=ios)
    use_wsplit = (ios == 0 .and. env_len > 0)
    wsplit_maxcfl_env = 1.0_WP
    call get_environment_variable('FESOM3_WSPLIT_MAXCFL', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) then
        read(env, *, iostat=ios) wsplit_maxcfl_env
        if (ios /= 0) error stop 'FESOM3_WSPLIT_MAXCFL: not a real'
    end if
    wsplit_mincfl_env = 0.5_WP*wsplit_maxcfl_env
    call get_environment_variable('FESOM3_WSPLIT_MINCFL', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) then
        read(env, *, iostat=ios) wsplit_mincfl_env
        if (ios /= 0) error stop 'FESOM3_WSPLIT_MINCFL: not a real'
    end if
    call wsplit_assert_params(wsplit_mincfl_env, wsplit_maxcfl_env)   ! the t_dyn obligation
    if (partit%mype == 0 .and. (shear_splines .or. N2_splines)) write(*,'(a,l1,a,l1,a)') &
        'fesom_lifecycle_native_mr: vertical splines ENABLED (shear=', shear_splines, &
        ', N2=', N2_splines, ')'
    ! M8c: FESOM3_CHL_SWEENEY -> read the Sweeney monthly chlorophyll climatology (production) instead
    ! of the const-0.1 fallback. Only meaningful with shortwave penetration (use_sw_pene) on.
    call get_environment_variable('FESOM3_CHL_SWEENEY', env, length=env_len, status=ios)
    use_chl_sweeney = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_KPP_NONLCL', env, length=env_len, status=ios)
    do_nonlcl = (ios == 0 .and. env_len > 0)
    ! perf: the per-step GLOBAL MPI_Allreduce stability diagnostic (max|eta|/|uv|/a_ice/m_ice/T/S) is a
    ! sync point FESOM2 does NOT have (it logs only every logfile_outfreq). OFF by default so the production
    ! driver matches FESOM2's per-step work; the M8e free-running stability run sets FESOM3_STEP_DIAG=1.
    call get_environment_variable('FESOM3_STEP_DIAG', env, length=env_len, status=ios)
    step_diag = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_TIMING_EVERY', env, length=env_len, status=ios)
    mon_every = 0
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) mon_every
    ! M6a-4: ALE vertical coordinate (default linfs). 'zstar' -> full free surface + real
    ! freshwater flux (use_virt_salt=.false., is_nonlinfs=1; mirror of the 1-rank M6a-3 wiring).
    call get_environment_variable('FESOM3_WHICH_ALE', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) which_ALE = trim(env)
    ! Timestep: FESOM3_STEP_PER_DAY steps per day (default 48 = 1800 s = 30 min, the CORE2
    ! namelist timestep). A finer mesh needs a larger value (e.g. 96 => 900 s). dt must be
    ! resolved HERE: init_stiff_mat_ale and ice_setup below capture it at setup time.
    spd = 48
    call get_environment_variable('FESOM3_STEP_PER_DAY', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) then
        read(env, *, iostat=ios) spd_env
        if (ios == 0 .and. spd_env > 0) spd = spd_env
    end if
    dt = 86400.0_WP / real(spd, WP)

    !===========================================================================
    ! model_init: MR mesh remap + geometry (set_partition -> read_mesh dispatches to
    ! read_mesh_local at npes>1; npes==1 reads the global mesh).
    call par_init(partit)
    ! M8c Step-3 FIX (root cause of the day-107 byte-divergence): FESOM2 runs with flush-to-zero ON
    ! (abrupt underflow) but FESOM3's process had it OFF, so FESOM3 RETAINED a denormal m_snow (~1e-309)
    ! where FESOM2 flushed it to exactly 0.0 -> the if(hsn>0) ice-albedo branch flipped at day ~107 ->
    ! a global SSH divergence. Match FESOM2 by flushing denormals to zero (also the physically-correct
    ! behaviour — a snow thickness of 1e-309 m IS zero). Confirmed: with this on, the 2-year JRA55
    ! headline is byte-exact (max|delta|=0). Re-asserted each step below to survive any library MXCSR reset.
    ftz_supported = ieee_support_underflow_control(1.0_WP)
    if (ftz_supported) then
        call ieee_set_underflow_mode(gradual=.false.)
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: flush-to-zero (FTZ) ON to match FESOM2 (denormals -> 0)'
    else if (partit%mype == 0) then
        write(*,'(a)') 'fesom_lifecycle_native_mr: WARNING ieee underflow control unsupported — FTZ not forced'
    end if
    call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
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
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_lifecycle_native_mr: nod2D=', mesh%nod2D, &
            ' elem2D=', mesh%elem2D, ' nl=', nl, ' CW swaps=', nsw
    if (partit%mype == 0) &
        write(*,'(a,i0,a,f0.1,a)') 'fesom_lifecycle_native_mr: step_per_day=', spd, &
            ' dt=', dt, ' s'

    !===========================================================================
    ! ALE depth/thickness state (linfs full cells; local sizes) — as fesom_lifecycle_mr.
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
    mesh%hbar = 0.0_MP; mesh%hbar_old = 0.0_MP; mesh%dhe = 0.0_MP
    mesh%hnode_new = mesh%hnode

    !===========================================================================
    ! dynamics state — COLD START (local sizes).
    allocate(dyn%uv(2, nl-1, nElemF), dyn%uv_rhs(2, nl-1, nElemF))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, nElemF))
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
    dyn%momadv_opt    = momadv_env
    if (partit%mype == 0 .and. dyn%momadv_opt == 1) write(*,'(a)') &
        'fesom_lifecycle_native_mr: momentum advection = VECTOR INVARIANT (momadv_opt=1)'
    if (partit%mype == 0 .and. rvo_upwind > 0.0_WP) write(*,'(a,f6.3)') &
        'fesom_lifecycle_native_mr: face-vorticity upwind blend rvo_upwind = ', rvo_upwind
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    dyn%use_wsplit    = use_wsplit            ! .false. unless FESOM3_WSPLIT (production precedent)
    dyn%wsplit_mincfl = wsplit_mincfl_env
    dyn%wsplit_maxcfl = wsplit_maxcfl_env
    if (partit%mype == 0 .and. dyn%use_wsplit) write(*,'(a,f8.4,a,f8.4,a)') &
        'fesom_lifecycle_native_mr: implicit vertical advection split w = w_e + w_i ENABLED ' // &
        '(use_wsplit; wsplit_mincfl=', dyn%wsplit_mincfl, ', wsplit_maxcfl=', dyn%wsplit_maxcfl, ')'

    !===========================================================================
    ! 2-tracer ocean state (data(1)=T, data(2)=S; local sizes) + do_ic3d phc3.0 IC.
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, nNodL), tracers%data(2)%values(nl-1, nNodL))
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
    call do_ic3d(tracers, ic, mesh, partit)
    if (partit%mype == 0) &
        write(*,'(a,2es12.4,a,2es12.4)') 'fesom_lifecycle_native_mr: IC(rank0 owned) T=', &
            minval(tracers%data(1)%values(:,1:nNodO)), maxval(tracers%data(1)%values(:,1:nNodO)), &
            '  S=', minval(tracers%data(2)%values(:,1:nNodO)), maxval(tracers%data(2)%values(:,1:nNodO))

    !===========================================================================
    ! tracer advection machinery (cold start = values; local sizes; adv_flux_hor owned edges).
    do tr_num = 1, 2
        allocate(tracers%data(tr_num)%valuesAB(nl-1, nNodL))
        allocate(tracers%data(tr_num)%valuesold(2, nl-1, nNodL))
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
    mix_scheme_nmb = 2           ! reduced-M2 mixing = PP (KPP is M5c/M5d)
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.
    use_windmix   = .false.

    !===========================================================================
    ! M4f GM/Redi (FESOM3_FER_GM / FESOM3_REDI): enable Gent-McWilliams bolus advection
    ! (+ Redi isopycnal diffusion) at MULTI-RANK. Same work_core config + inits as the 1-rank
    ! fesom_lifecycle_native (M4e), but the GM/Redi work arrays are LOCAL-sized (nNodL/nElemF)
    ! like the rest of the MR state. step_oce threads the optional partit to every M4 routine
    ! (producers + init_Redi_GM/fer_solve_Gamma/fer_gamma2vel + the Redi diff terms), each of
    ! which already does owned-loop bounds + the FESOM2 halo exchanges (sigma_xy/neutral_slope/
    ! slope_tapered/fer_c/fer_K/Ki/fer_gamma exchange_nod, fer_uv exchange_elem, fer_w in
    ! vert_vel_ale). sw_alpha/sw_beta/tr_z are computed over owned+halo (no exchange needed).
    if (use_fer_gm .or. use_redi) then
        Fer_GM = .true.        ! Redi rides the GM coupling (Ki = max(fer_scal*Redi_Kmax, K_GM_min))
        K_GM_max = 1000.0_WP; K_GM_min = 2.0_WP; K_GM_bvref = 1
        K_GM_rampmax = -1.0_WP; K_GM_rampmin = -1.0_WP; K_GM_resscalorder = 2.0_WP
        K_GM_cm = 3.0_WP; K_GM_cmin = 0.1_WP; K_GM_Ktaper = .false.
        scaling_Ferreira = .false.; scaling_Rossby = .false.; scaling_resolution = .true.
        scaling_FESOM14 = .false.; scaling_GMzexp = .true.; scaling_GINsea = .false.
        GMzexp_zref = 500.0_WP; GMzexp_smin = 0.6_WP
        allocate(dyn%fer_uv(2, nl-1, nElemF), dyn%fer_w(nl, nNodL))
        allocate(dyn%work%sw_alpha(nl-1, nNodL), dyn%work%sw_beta(nl-1, nNodL))
        allocate(dyn%work%sigma_xy(2, nl-1, nNodL))
        allocate(dyn%work%fer_K(nl, nNodL), dyn%work%fer_c(nNodL), dyn%work%fer_scal(nNodL))
        allocate(dyn%work%fer_gamma(2, nl, nNodL))
        dyn%fer_uv = 0.0_WP; dyn%fer_w = 0.0_WP
        dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP; dyn%work%sigma_xy = 0.0_WP
        dyn%work%fer_K = 500.0_WP; dyn%work%fer_c = 1.0_WP; dyn%work%fer_scal = 0.0_WP
        dyn%work%fer_gamma = 0.0_WP
        if (use_redi) then
            Redi = .true.; Redi_Ktaper = .true.; Redi_Kmax = 0.0_WP; Redi_Kmin = 100.0_WP
            scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP; scaling_LDD97 = .false.
            K_hor = 0.0_WP
            allocate(dyn%work%Ki(nl-1, nNodL), dyn%work%fer_tapfac(nl-1, nNodL))
            allocate(dyn%work%neutral_slope(3, nl-1, nNodL), dyn%work%slope_tapered(3, nl-1, nNodL))
            allocate(tracers%work%tr_z(nl, nNodL))
            dyn%work%Ki = 0.0_WP; dyn%work%fer_tapfac = 1.0_WP
            dyn%work%neutral_slope = 0.0_WP; dyn%work%slope_tapered = 0.0_WP
            tracers%work%tr_z = 0.0_WP
            if (partit%mype == 0) write(*,'(a)') &
                'fesom_lifecycle_native_mr: Fer_GM + Redi ENABLED (work_core GM+Redi config)'
        else
            Redi = .false.
            if (partit%mype == 0) write(*,'(a)') &
                'fesom_lifecycle_native_mr: Fer_GM ENABLED (work_core GM config; Redi off)'
        end if
    end if

    !===========================================================================
    ! M5d KPP vertical mixing (FESOM3_MIX_KPP) at MULTI-RANK — the M5c block with LOCAL-sized
    ! (nNodL) dyn%work arrays. The KPP module is optional-`partit` from the start (oce_mixing_
    ! kpp_driver does the owned-loop + exchange_nod(blmc/diffK/ghats/viscA) + smooth_blmc), so
    ! this is pure WIRING (the M4f lesson). sw_alpha/sw_beta were allocated by the GM/Redi block.
    if (use_kpp) then
        mix_scheme_nmb = 1
        Ricr = 0.3_WP; concv = 1.6_WP
        visc_sh_limit = 5.0e-3_WP; diff_sh_limit = 5.0e-3_WP
        ref_sss = 34.0_WP; ref_sss_local = .true.
        if (.not. allocated(dyn%work%sw_alpha)) then
            allocate(dyn%work%sw_alpha(nl-1, nNodL), dyn%work%sw_beta(nl-1, nNodL))
            dyn%work%sw_alpha = 0.0_WP; dyn%work%sw_beta = 0.0_WP
        end if
        allocate(dyn%work%Kv_double(nl, nNodL, tracers%num_tracers))
        allocate(dyn%work%viscA_kpp(nl, nNodL), dyn%work%blmc(nl, nNodL, 3))
        allocate(dyn%work%ghats(nl-1, nNodL), dyn%work%dkm1(nNodL, 3))
        allocate(dyn%work%dbsfc(nl, nNodL), dyn%work%dVsq(nl, nNodL))
        allocate(dyn%work%sw_3d(nl, nNodL))
        allocate(dyn%work%hbl(nNodL), dyn%work%bfsfc(nNodL))
        allocate(dyn%work%stable(nNodL), dyn%work%caseA(nNodL))
        allocate(dyn%work%ustar(nNodL), dyn%work%Bo(nNodL), dyn%work%kbl(nNodL))
        dyn%work%Kv_double = 0.0_WP; dyn%work%viscA_kpp = 0.0_WP; dyn%work%blmc = 0.0_WP
        dyn%work%ghats = 0.0_WP; dyn%work%dkm1 = 0.0_WP; dyn%work%dbsfc = 0.0_WP
        dyn%work%dVsq = 0.0_WP; dyn%work%sw_3d = 0.0_WP
        dyn%work%hbl = 0.0_WP; dyn%work%bfsfc = 0.0_WP; dyn%work%stable = 0.0_WP
        dyn%work%caseA = 0.0_WP; dyn%work%ustar = 0.0_WP; dyn%work%Bo = 0.0_WP; dyn%work%kbl = 0
        call oce_mixing_kpp_init(Ricr, concv)   ! wmt/wst lookup tables + Vtc/cg (once)
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: KPP vertical mixing ENABLED (work_core KPP)'
    end if
    !===========================================================================
    ! M7d TKE vertical mixing (FESOM3_MIX_TKE) at MULTI-RANK — the M7c block with LOCAL-sized
    ! (nNodL) dyn%work arrays. oce_mixing_tke.calc_cvmix_tke is optional-`partit` from the start
    ! (owned-node loop over nNodO + exchange_nod(tke_Kv)/exchange_nod(tke_Av) BEFORE the element
    ! average; tke is NEVER exchanged — the recurrence is partition-local), so this is pure WIRING
    ! (the M4f/M5d lesson). Pair FESOM3_SW_PENE=1 (work_*_tke use_sw_pene=.true.).
    if (use_tke) then
        mix_scheme_nmb = 5
        allocate(dyn%work%tke(nl, nNodL), dyn%work%tke_Av(nl, nNodL), dyn%work%tke_Kv(nl, nNodL))
        dyn%work%tke = 0.0_WP; dyn%work%tke_Av = 0.0_WP; dyn%work%tke_Kv = 0.0_WP
        call tke_init(tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, tke_kappaM_min, &
                      tke_kappaM_max, tke_min, tke_surf_min, tke_mxl_choice, tke_only, &
                      tke_use_ubound_dirichlet, tke_use_lbound_dirichlet, tke_dolangmuir)
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: TKE vertical mixing ENABLED (cvmix_TKE)'
    end if
    !===========================================================================
    ! M5d shortwave penetration (FESOM3_SW_PENE); cal_shortwave_rad (partit) fills dyn%work%sw_3d over
    ! owned+halo after oce_fluxes. sw_3d allocated by the KPP block. chl is seeded to the const-0.1
    ! fallback here; the Sweeney monthly read (M8c, needs clock_init's `month`) overwrites it below.
    if (do_swpene) then
        use_sw_pene = .true.
        allocate(chl(nNodL)); chl = 0.1_WP
        if (.not. allocated(dyn%work%sw_3d)) then
            allocate(dyn%work%sw_3d(nl, nNodL)); dyn%work%sw_3d = 0.0_WP
        end if
        if (partit%mype == 0) then
            if (use_chl_sweeney) then
                write(*,'(a)') 'fesom_lifecycle_native_mr: shortwave penetration ENABLED, chl=Sweeney monthly'
            else
                write(*,'(a)') 'fesom_lifecycle_native_mr: shortwave penetration ENABLED, chl=const 0.1'
            end if
        end if
    end if
    !===========================================================================
    ! M5d KPP nonlocal counter-gradient flux (FESOM3_KPP_NONLCL) — DEAD in production; the gate
    ! variant turns it on BOTH sides to byte-verify the ghats term at multi-rank.
    if (do_nonlcl) then
        use_kpp_nonlclflx = .true.
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: KPP nonlocal counter-gradient flux ENABLED (ghats)'
    end if

    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt, partit)

    !===========================================================================
    ! ice setup (allocate local-sized + ice_mass_matrix_fill + cold-start ice_initial_state)
    ! + the CORE2 &ice_dyn / &ice_therm config overrides (same as fesom_lifecycle_native).
    call ice_setup(ice, tracers, mesh, dt, partit)
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
    ! atmflux: allocate (local nNodL) + scalar config (the CORE2 namelist scalars).
    call alloc_atm(atm, nNodL)
    atm%Ch_atm_ice    = 0.00175_WP
    atm%Ce_atm_ice    = 0.00175_WP
    atm%ref_sss       = 34.0_WP
    atm%ref_sss_local = .true.
    atm%use_virt_salt = (trim(which_ALE)=='linfs')   ! M6a-4: zstar -> real freshwater flux
    atm%l_snow        = .true.
    atm%surf_relax_S  = 1.929e-06_WP

    !===========================================================================
    ! M8a/M8b: model clock + run-length driver. Set the CORE2 timestep into mod_config so the
    ! ported `clock` advances mod_config%dt per step (== the local dt parameter) and forcing_sbc_do
    ! reads the same dt. FESOM3_START_CLOCK ("t d y", default "0 1 1948") positions a COLD start:
    ! rank 0 writes a 2-identical-line .clock file, then every rank reads it through the faithful
    ! clock_init (the M8b/c boundary gates set FESOM3_START_CLOCK to "0 31 1948" / "0 365 1948").
    ! nsteps comes from clock_nsteps(run_length/unit) with FESOM3_NSTEPS keeping the short-gate
    ! override. clock_init MUST precede the forcing cold-start build below: M8b seeds rdate_cold +
    ! the SSS month from the clock's yearnew/daynew/timenew/month (prereq #6).
    ! M8c: resolve the forcing dataset FIRST — include_fleapyear feeds check_fleapyr inside
    ! clock_init below, and forc_calendar drives the cold-start/per-step rdate. Default CORE2 keeps
    ! the legacy regression gates (noleap) byte-unchanged; FESOM3_FORCING=JRA55 = the production
    ! target (gregorian + leap year cycle, the work_zstar_tke config).
    call get_environment_variable('FESOM3_FORCING', forc_set, length=env_len, status=ios)
    if (ios /= 0 .or. env_len == 0) forc_set = 'CORE2'
    if (trim(forc_set) == 'JRA55') then
        forc_calendar = 'gregorian'; include_fleapyear = .true.
    else
        forc_calendar = 'noleap';    include_fleapyear = .false.
    end if
    step_per_day    = spd
    cfg_dt          = dt              ! 86400/spd (default 48 = 1800 s); clock + forcing_sbc_do read mod_config%dt
    run_length      = 2
    run_length_unit = 'y'
    call get_environment_variable('FESOM3_RUN_LENGTH', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) run_length
    call get_environment_variable('FESOM3_RUN_UNIT', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) run_length_unit = trim(env)
    ! RestartInPath for the .clock file: FESOM3_RESTART_IN, else dirname(FESOM_DUMP_FILE), else ./
    call get_environment_variable('FESOM3_RESTART_IN', restart_in, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) then
        if (restart_in(env_len:env_len) /= '/') restart_in = trim(restart_in)//'/'
    else
        call get_environment_variable('FESOM_DUMP_FILE', restart_in, length=env_len, status=ios)
        idx = 0
        if (ios == 0 .and. env_len > 0) idx = index(trim(restart_in), '/', back=.true.)
        if (idx > 0) then; restart_in = restart_in(1:idx); else; restart_in = './'; end if
    end if
    RestartInPath = trim(restart_in)
    ! ---- Restart Stage 5 (Task 5.1): restart-mode detection, BEFORE the cold .clock overwrite ----
    ! A restart LAUNCH iff FESOM3_RESTART=<dir> is set AND a COMPLETE checkpoint resolves under
    ! RestartInPath (restart.latest -> a folder carrying checkpoint.json). FESOM3_RESTART also enables
    ! the checkpoint WRITE path (do_restart) and sets RestartOutPath (where clock_finish + the folders go).
    call get_environment_variable('FESOM3_RESTART', restart_dir, length=env_len, status=ios)
    do_restart = (ios == 0 .and. env_len > 0)
    if (do_restart) then
        ! strip a trailing '/': restart_write joins <dir>/<folder>; RestartOutPath keeps a trailing slash.
        if (len_trim(restart_dir) > 1) then
            if (restart_dir(len_trim(restart_dir):len_trim(restart_dir)) == '/') &
                restart_dir = restart_dir(1:len_trim(restart_dir)-1)
        end if
        RestartOutPath = trim(restart_dir)//'/'
    end if
    restart_mode = .false.
    if (do_restart) then
        call restart_resolve_latest(trim(RestartInPath), rl_folder, rl_ok)
        restart_mode = rl_ok
    end if
    ! restart cadence knobs (default 1 'y'; &nml_restart-equivalent via FESOM3_* env per the plan's
    ! "env fallbacks are fine" allowance). restart_keep / compressor / n_writers enter via restart_init.
    restart_length = 1; restart_length_unit = 'y'
    call get_environment_variable('FESOM3_RESTART_LENGTH', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) restart_length
    call get_environment_variable('FESOM3_RESTART_UNIT', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) restart_length_unit = trim(env)
    ! cold-start clock values from FESOM3_START_CLOCK (default 0 1 <start year>); both lines equal
    ! => cold start. JRA55 forcing starts 1958 (no 1949), CORE2 NCAR at 1948.
    clk_t0 = 0.0_WP; clk_d0 = 1; clk_y0 = merge(1958, 1948, trim(forc_set) == 'JRA55')
    call get_environment_variable('FESOM3_START_CLOCK', start_clock, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(start_clock, *, iostat=ios) clk_t0, clk_d0, clk_y0
    ! In restart mode SKIP the unconditional cold .clock overwrite so clock_init reads the PREVIOUS
    ! segment's clock_finish-written .clock (two DIFFERING lines => r_restart=.true.); the overwrite
    ! would force two EQUAL lines => r_restart always false => the READ hook below would be dead
    ! (verified mod_clock.F90:127-132).
    if (.not. restart_mode .and. partit%mype == 0) then
        open(newunit=clk_unit, file=trim(RestartInPath)//trim(runid)//'.clock', &
             status='replace', action='write')
        write(clk_unit,*) clk_t0, clk_d0, clk_y0
        write(clk_unit,*) clk_t0, clk_d0, clk_y0
        close(clk_unit)
    end if
    call MPI_Barrier(partit%MPI_COMM_FESOM, ierr)
    if (restart_mode .and. partit%mype == 0) write(*,'(a)') &
        'fesom_lifecycle_native_mr: RESTART MODE — checkpoint found under '//trim(RestartInPath)// &
        ' (cold .clock overwrite skipped; clock_init reads the chained .clock)'
    call clock_init(partit)
    nsteps = clock_nsteps(partit)
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps   ! short-gate override

    !===========================================================================
    ! native CORE2 forcing read setup (REQUIRED — fully native). frc%nnod = nNodL so the
    ! 8 NCAR fields read + bilinear + g2r-rotate fill the OWNED+HALO nodes (each node
    ! interpolates from the full global raw grid — partition-independent, no exchange).
    call get_environment_variable('FESOM3_FORCING_DIR', forcing_dir)
    if (len_trim(forcing_dir) == 0) then
        if (partit%mype == 0) write(*,'(a)') &
            'fesom_lifecycle_native_mr: FULLY NATIVE requires FESOM3_FORCING_DIR'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if
    ! M8b (prereq #6): cold-start forcing from the clock (yearnew/daynew/timenew/month set by
    ! clock_init above) so a non-day-1 START_CLOCK initialises forcing at the SAME point as F2.
    fyear = yearnew
    frc%nfld = 8;  frc%nnod = nNodL
    frc%imm = 1; frc%idd = 1; frc%freq = 1
    frc%ic_cyclic = .true.; frc%rotated_grid = .true.
    frc%i_xwind = 1; frc%i_ywind = 2
    ! field->slot order matches FESOM2 update_atm_forcing's atmdata indexing + conversions
    ! (gen_forcing_couple.F90:681-694): 1,2=wind 3=humi 4=sw 5=lw 6=Tair(-273.15) 7=rain 8=snow(/1000).
    if (trim(forc_set) == 'JRA55') then
        ! JRA55-do-v1.4.0: 3-hourly, gregorian (days since 1900-01-01), start-of-interval => tmid=0
        ! (the mid-point shift fires). var-name == file-name.
        frc%iyear = 1900; frc%tmid = 0
        frc%f(1)%file_base = trim(forcing_dir)//'/uas.';  frc%f(1)%varname = 'uas'
        frc%f(2)%file_base = trim(forcing_dir)//'/vas.';  frc%f(2)%varname = 'vas'
        frc%f(3)%file_base = trim(forcing_dir)//'/huss.'; frc%f(3)%varname = 'huss'
        frc%f(4)%file_base = trim(forcing_dir)//'/rsds.'; frc%f(4)%varname = 'rsds'
        frc%f(5)%file_base = trim(forcing_dir)//'/rlds.'; frc%f(5)%varname = 'rlds'
        frc%f(6)%file_base = trim(forcing_dir)//'/tas.';  frc%f(6)%varname = 'tas'
        frc%f(7)%file_base = trim(forcing_dir)//'/prra.'; frc%f(7)%varname = 'prra'
        frc%f(8)%file_base = trim(forcing_dir)//'/prsn.'; frc%f(8)%varname = 'prsn'
    else
        ! CORE2 NCAR: 6-hourly winds/q/t (1460) + daily rad (365) + monthly precip (12), noleap, tmid=1.
        frc%iyear = 1948; frc%tmid = 1
        frc%f(1)%file_base = trim(forcing_dir)//'/u_10.';        frc%f(1)%varname = 'U_10_MOD'
        frc%f(2)%file_base = trim(forcing_dir)//'/v_10.';        frc%f(2)%varname = 'V_10_MOD'
        frc%f(3)%file_base = trim(forcing_dir)//'/q_10.';        frc%f(3)%varname = 'Q_10_MOD'
        frc%f(4)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(4)%varname = 'SWDN_MOD'
        frc%f(5)%file_base = trim(forcing_dir)//'/ncar_rad.';    frc%f(5)%varname = 'LWDN_MOD'
        frc%f(6)%file_base = trim(forcing_dir)//'/t_10.';        frc%f(6)%varname = 'T_10_MOD'
        frc%f(7)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(7)%varname = 'RAIN'
        frc%f(8)%file_base = trim(forcing_dir)//'/ncar_precip.'; frc%f(8)%varname = 'SNOW'
    end if
    rdate_cold = real(forcing_julday(yearnew,1,1,forc_calendar),WP) &  ! cold-start rdate, NO half-step
               + real(daynew-1,WP) + timenew/86400._WP                 ! (FESOM2 nc_sbc_ini:643-644)
    call forcing_alloc(frc)
    do fld = 1, frc%nfld
        call forcing_read_grid(frc, fld, fyear)
        call forcing_build_bilin(frc, fld, mesh, partit)
    end do
    do fld = 1, frc%nfld
        call forcing_getcoeffld(frc, fld, fyear, rdate_cold, mesh, partit)
    end do
    call forcing_rotate_wind(frc, mesh, partit)
    allocate(nuw(nNodL), nvw(nNodL), nta(nNodL), nsh(nNodL), &
             nswr(nNodL), nlw(nNodL), npr(nNodL), nps(nNodL), &
             ncd(nNodL), nch(nNodL), nce(nNodL), &
             nsx(nNodL), nsy(nNodL), nix(nNodL), niy(nNodL))
    ! native runoff (CORE2 single-slice climatology, read ONCE) + SSS at the cold-start month
    ! (M8b prereq #6: i=month, not hard-coded 1; M8c adds the in-loop monthly read-ahead). PASS
    ! partit so read_other_NetCDF interpolates over owned+halo (myDim+eDim; read_dist_partition
    ! sets myDim — the L42 trap is the 1-rank case where myDim=0, here npes>1 sets it).
    call get_environment_variable('FESOM3_RUNOFF_FILE', runoff_file)
    if (len_trim(runoff_file) == 0) &
        runoff_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/CORE2_runoff.nc'
    call get_environment_variable('FESOM3_SSS_FILE', sss_file)
    if (len_trim(sss_file) == 0) &
        sss_file = '/pool/data/AWICM/FESOM2/FORCING/JRA55-do-v1.4.0/PHC2_salx.nc'
    allocate(nro(nNodL), nss(nNodL))
    if (partit%npes == 1) then
        ! 1-rank: omit partit (read_other_NetCDF then counts mesh%nod2D — the L42 trap fix).
        call read_other_NetCDF(trim(runoff_file), 'Foxx_o_roff', 1, nro, .false., .true., mesh)
        call read_other_NetCDF(trim(sss_file), 'SALT', month, nss, .true., .true., mesh)
    else
        call read_other_NetCDF(trim(runoff_file), 'Foxx_o_roff', 1, nro, .false., .true., mesh, partit)
        call read_other_NetCDF(trim(sss_file), 'SALT', month, nss, .true., .true., mesh, partit)
    end if
    nro = nro / 1000.0_WP
    ! M8c: production Sweeney monthly chlorophyll, read at the cold-start month (overwrites the
    ! const-0.1 seed). Guard on use_sw_pene so chl is allocated; FESOM3_CHL_FILE overrides the path.
    if (use_sw_pene .and. use_chl_sweeney) then
        call get_environment_variable('FESOM3_CHL_FILE', chl_file)
        if (len_trim(chl_file) == 0) &
            chl_file = '/pool/data/AWICM/FESOM2/FORCING/Sweeney/Sweeney_2005.nc'
        if (partit%npes == 1) then
            call read_other_NetCDF(trim(chl_file), 'chl', month, chl, .true., .true., mesh)
        else
            call read_other_NetCDF(trim(chl_file), 'chl', month, chl, .true., .true., mesh, partit)
        end if
    end if
    if (partit%mype == 0) &
        write(*,'(a,2es12.4,a,2es12.4,a)') 'fesom_lifecycle_native_mr: native runoff[', &
            minval(nro(1:nNodO)), maxval(nro(1:nNodO)), '] Ssurf[', &
            minval(nss(1:nNodO)), maxval(nss(1:nNodO)), '] (rank0 owned)'

    !===========================================================================
    ! ocean-step surface BC arrays. heat_flux/water_flux/virtual_salt/relax_salt/stress_surf
    ! are produced NATIVELY each step into atm%* / stress_surf; Ki/real_salt_flux stay 0 (M4).
    ! M6a-4: is_nonlinfs=1 for zstar (bc_surface uses real_salt_flux + advective-heat). Local sizes.
    is_nonlinfs = merge(1.0_WP, 0.0_WP, trim(which_ALE)/='linfs')
    allocate(Ki(nl-1, nNodL), real_salt_flux(nNodL), stress_surf(2, nElemF))
    Ki = 0.0_WP; real_salt_flux = 0.0_WP; stress_surf = 0.0_WP

    !===========================================================================
    ! M8c bisection: FESOM3_ATMFLUX_CHECK = an F2 atmflux dump (fesom_atmflux_dump order, np=1) =>
    ! each step, read the F2 record and print per-field max|delta| between the NATIVE atmosphere and
    ! the oracle's, localizing a forcing-read mismatch BEFORE it propagates into ssh_rhs. np=1 only
    ! (the dump + the local mesh are both global-sized there).
    do_atm_chk = .false.
    call get_environment_variable('FESOM3_ATMFLUX_CHECK', atmchk_file, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0 .and. partit%npes == 1) then
        open(newunit=atmchk_unit, file=trim(atmchk_file), status='old', form='unformatted', &
             access='stream', action='read', iostat=ios)   ! fesom_atmflux_dump is STREAM
        if (ios == 0) then; do_atm_chk = .true.; allocate(tscr(nNodL)); end if
    end if

    !===========================================================================
    ! per-rank dump (mod_dump: gid-keyed probes, the rank owning a probe writes it).
    allocate(idlist(max(nNodO, nElemO)))
    do n = 1, size(idlist); idlist(n) = n; end do
    if (partit%npes == 1) then
        call dump_init(partit%mype, nNodO, idlist(1:nNodO), nElemO, idlist(1:nElemO))
    else
        call dump_init(partit%mype, nNodO, partit%myList_nod2D, nElemO, partit%myList_elem2D)
    end if

    ! M9 Stage 1: emit fesom.mesh.diag.zarr after setup (mesh + zbar_e_bot ready), env-gated.
    block
        character(len=4096) :: md_env, md_out
        call get_environment_variable('FESOM3_MESHDIAG', md_env)
        if (len_trim(md_env) > 0) then
            call get_environment_variable('FESOM3_MESHDIAG_OUT', md_out)
            if (len_trim(md_out) == 0) md_out = 'fesom.mesh.diag.zarr'
            call meshdiag_write(trim(md_out), mesh, partit)
            if (partit%mype == 0) write(*,'(a)') 'M9: wrote '//trim(md_out)
        end if
    end block

    ! M9 Stage 2: register field output, env-gated by FESOM3_OUTPUT=<dir>. Configuration (the variable
    ! list + global writer knobs) comes from namelist.io in the rundir (Task 2.6) when present; otherwise
    ! the fixed default set at a FESOM3_OUTPUT_EVERY step cadence. FESOM3_* env overrides namelist knobs.
    do_output = .false.; out_every = 1; have_nml = .false.; n_io = 0
    call get_environment_variable('FESOM3_OUTPUT', out_env)
    if (len_trim(out_env) > 0) then
        do_output = .true.
        out_dir_o = out_env
        call get_environment_variable('FESOM3_OUTPUT_EVERY', env)
        if (len_trim(env) > 0) then; read(env,*,iostat=ios) out_every; if (ios /= 0) out_every = 1; end if
        if (out_every < 1) out_every = 1
        ! namelist.io (rundir): &nml_general global knobs + &nml_list output variable rows.
        call get_environment_variable('FESOM3_NAMELIST_IO', nml_io_path)
        if (len_trim(nml_io_path) == 0) nml_io_path = 'namelist.io'
        call means_read_namelist(trim(nml_io_path), iocfg, iolist, n_io, have_nml)
        if (have_nml) then
            call means_init(oio, trim(out_dir_o), mesh, partit, calendar=trim(forc_calendar), &
                            chunk_horiz=iocfg%chunk_horiz, n_writers=iocfg%n_writers, &
                            chunk_time=iocfg%chunk_time, chunk_vert=iocfg%chunk_vert, &
                            compressor=trim(iocfg%compressor), filesplit_freq=iocfg%filesplit_freq, &
                            vec_frame=trim(iocfg%vec_frame))
            do ii = 1, n_io
                call register_output_var(oio, iolist(ii))
            end do
        else
            call means_init(oio, trim(out_dir_o), mesh, partit, calendar=trim(forc_calendar))
            call register_default_outputs(oio, out_every)
        end if
        if (partit%mype == 0) then
            if (have_nml) then
                write(*,'(a,i0,a)') 'M9: field output ON -> '//trim(out_dir_o)//' (namelist.io: ', &
                    n_io, ' streams)'
            else
                write(*,'(a,i0,a)') 'M9: field output ON -> '//trim(out_dir_o)// &
                    ' (default set, every ', out_every, ' steps)'
            end if
        end if
    end if

    !===========================================================================
    ! Restart Stage 5 (Task 5.3 REGISTER + Task 5.1 READ). When FESOM3_RESTART=<dir> is set, register the
    ! FULL prognostic state (live pointers — oce + ice incl. EVP sigma + MP hbar/hnode) so the in-loop /
    ! end-of-run hooks below can checkpoint it. When this is a RESTART launch (r_restart, set by clock_init
    ! from the chained .clock), restart_read then restores the newest checkpoint INTO those live arrays,
    ! OVERWRITING the cold state allocated above — placed after the whole cold init (dyn/tracers/ice +
    ! forcing) and before the loop; nothing between clock_init and here consumes the cold ocean/ice state
    ! (the forcing build reads only the clock + atmosphere). restart_read resolves the folder via
    ! restart.latest under RestartInPath; clock_year/day/time_sec drive the time-vs-clock safety warning
    ! (time_sec is the sec-of-day == timenew, matching the folder tag + checkpoint.json convention).
    if (do_restart) then
        call restart_init(rst, mesh, partit)
        call restart_register_state(rst, dyn, tracers, ice, mesh, mix_scheme=mix_scheme_nmb)
        if (partit%mype == 0) write(*,'(a,i0,a)') &
            'restart: registered ', rst%nf, ' prognostic fields -> '//trim(RestartOutPath)
        if (r_restart) then
            call restart_read(rst, trim(RestartInPath), mesh, partit, &
                              clock_year=yearnew, clock_day=daynew, clock_time_sec=real(timenew, real64))
            if (partit%mype == 0) write(*,'(a,f9.1,a,i0,a,i0)') &
                'restart: state restored (resumed run) at clock time=', timenew, ' day=', daynew, &
                ' year=', yearnew
        end if
    end if

    !===========================================================================
    ! runloop: ocean2ice -> [native atm] -> ice_timestep -> oce_fluxes_mom -> oce_fluxes ->
    ! step_oce (the FESOM2 runloop order; update_atm_forcing replaced by apply_native_forcing).
    call timer_init()
    wrote_final = .false.
    do n = 1, nsteps
        call clock                                   ! M8a: advance the model clock (top of step)
        if (ftz_supported) call ieee_set_underflow_mode(gradual=.false.)   ! M8c: keep FTZ on each step (match FESOM2; survive library MXCSR resets)
        call timer_start(TMR_OCEAN2ICE)
        call ocean2ice(ice, dyn, tracers, mesh, partit)
        call timer_stop(TMR_OCEAN2ICE)
        ! native atmosphere (after ocean2ice -> srfoce live; before EVP -> ice%uice/vice are
        ! the previous step's, as FESOM2 update_atm_forcing uses for the wind-on-ice stress).
        call timer_start(TMR_FORCING)
        call apply_native_forcing(n)
        call timer_stop(TMR_FORCING)
        if (do_atm_chk) call atm_selfcheck(n)
        call timer_start(TMR_ICE)
        call ice_timestep(ice, mesh, atm, partit)
        call timer_stop(TMR_ICE)
        call timer_start(TMR_FLUXES)
        call oce_fluxes_mom(ice, atm, stress_surf, mesh, partit)
        call oce_fluxes(ice, tracers, atm, mesh, partit)
        ! M5d: shortwave penetration after oce_fluxes (over owned+halo via partit) — fills
        ! dyn%work%sw_3d + adds the visible band back to atm%heat_flux. albw = ice%thermo%albw.
        if (use_sw_pene) call cal_shortwave_rad(.true., ice%thermo%albw, atm%shortwave, chl, &
                              ice%data(1)%values(1:nNodL), atm%heat_flux, dyn%work%sw_3d, mesh, partit)
        call timer_stop(TMR_FLUXES)
        ! M5d: stress_node_surf=atm%stress_node_surf (oce_fluxes_mom) feeds KPP ustar (PP ignores it).
        call timer_start(TMR_STEP_OCE)
        ! Restart Stage 5 (Task 5.2): first-resumed-step AB guard. lfirst = (n==1) .and. .not. r_restart
        ! so the first step of a RESUMED run blends AB2 (ff=ab2) over the RESTORED uv_rhsAB instead of
        ! forward-Euler (ff=1.0) discarding it (oce_dyn_velrhs.F90:62-64,135-136). A cold run keeps the
        ! Euler first step (r_restart=.false.).
        call step_oce(n, dt, (n == 1) .and. .not. r_restart, dyn, tracers, mesh, Ki, &
                      atm%heat_flux, atm%water_flux, atm%virtual_salt, atm%relax_salt, &
                      atm%real_salt_flux, is_nonlinfs, stress_surf, partit, &   ! M6a-4: native rsf (zstar)
                      stress_node_surf=atm%stress_node_surf)
        call timer_stop(TMR_STEP_OCE)
        ! M9 Stage 2: field output. Accumulate the live (post-step_oce) state EVERY step (FESOM2
        ! update_means) — snapshots overwrite (count=1), means sum — then means_output evaluates each
        ! field's freq/unit event to decide which write now. ABOVE the step_diag cycle so production
        ! steps output too. Only registered fields accumulate (means_has).
        if (do_output) then
            if (means_has(oio,'ssh'))    call means_accumulate(oio,'ssh',    dyn%eta_n(1:nNodO))
            if (means_has(oio,'sst'))    call means_accumulate(oio,'sst',    tracers%data(1)%values(1,1:nNodO))
            if (means_has(oio,'sss'))    call means_accumulate(oio,'sss',    tracers%data(2)%values(1,1:nNodO))
            if (means_has(oio,'a_ice'))  call means_accumulate(oio,'a_ice',  ice%data(1)%values(1:nNodO))
            if (means_has(oio,'m_ice'))  call means_accumulate(oio,'m_ice',  ice%data(2)%values(1:nNodO))
            if (means_has(oio,'m_snow')) call means_accumulate(oio,'m_snow', ice%data(3)%values(1:nNodO))
            if (means_has(oio,'temp'))   call means_accumulate(oio,'temp',   tracers%data(1)%values(1:nl-1, 1:nNodO))
            if (means_has(oio,'salt'))   call means_accumulate(oio,'salt',   tracers%data(2)%values(1:nl-1, 1:nNodO))
            if (means_has(oio,'w'))      call means_accumulate(oio,'w',      dyn%w(1:nl, 1:nNodO))
            if (means_has(oio,'unod'))   call means_accumulate(oio,'unod',   dyn%uvnode(1, 1:nl-1, 1:nNodO))
            if (means_has(oio,'vnod'))   call means_accumulate(oio,'vnod',   dyn%uvnode(2, 1:nl-1, 1:nNodO))
            ! Task 2.7 ELEMENT fields (owned elements nElemO): u/v (dyn%uv), Av (full levels), GM bolus.
            if (means_has(oio,'u'))      call means_accumulate(oio,'u',      dyn%uv(1, 1:nl-1, 1:nElemO))
            if (means_has(oio,'v'))      call means_accumulate(oio,'v',      dyn%uv(2, 1:nl-1, 1:nElemO))
            if (means_has(oio,'Av'))     call means_accumulate(oio,'Av',     dyn%work%Av(1:nl, 1:nElemO))
            if (means_has(oio,'bolus_u')) call means_accumulate(oio,'bolus_u', dyn%fer_uv(1, 1:nl-1, 1:nElemO))
            if (means_has(oio,'bolus_v')) call means_accumulate(oio,'bolus_v', dyn%fer_uv(2, 1:nl-1, 1:nElemO))
            oclk%year = yearnew; oclk%yearstart = yearstart; oclk%daynew = daynew; oclk%ndpyr = ndpyr
            oclk%month = month;  oclk%day_in_month = day_in_month
            oclk%ndim_month = num_day_in_month(fleapyear, month); oclk%timenew = real(timenew, real64)
            call means_output(oio, n, oclk)
        end if
        ! Restart Stage 5 (Task 5.3): periodic checkpoint of the live (post-step_oce) state — the same
        ! end-of-step instant M9 output sees. restart_due is the GLOBAL is_due cadence (restart_length/
        ! unit; the last step is always due). clock_finish writes RestartOutPath//runid//.clock after
        ! EVERY write (mirror io_restart.F90:573-578) so a mid-run checkpoint is resumable next launch —
        ! else the .clock stays at the cold-start time and r_restart stays false. The folder tag +
        ! checkpoint.json use clock_finish's own year-rollover normalization so the stamps match the
        ! .clock that the next segment's clock_init reads back (time_sec = sec-of-day = timenew).
        if (do_restart) then
            if (restart_due(n)) then
                cp_year = yearnew; cp_day = daynew; cp_tsec = real(timenew, real64)
                if (daynew == ndpyr .and. timenew == 86400._WP) then
                    cp_tsec = 0.0_real64; cp_day = 1; cp_year = yearold + 1
                end if
                call restart_write(rst, trim(restart_dir), cp_year, cp_day, cp_tsec, globalstep=n)
                if (partit%mype == 0) call clock_finish()
                if (n == nsteps) wrote_final = .true.
                if (partit%mype == 0) write(*,'(a,i0,a,i4.4,a,i3.3,a,i5.5)') &
                    'restart: checkpoint written at step ', n, ' -> fesom.', cp_year, '.', cp_day, &
                    '.', int(cp_tsec)
            end if
        end if
        ! perf monitoring: periodic (cumulative) per-component timing report every mon_every steps.
        ! Collective (all ranks call it — placed before the step_diag cycle); 0 => only the final report.
        ! NOT a single `.and.`: Fortran does not guarantee short-circuit evaluation, and
        ! mon_every defaults to 0, so `mod(n, mon_every)` would be an integer divide by zero
        ! on every step of a default run. Same hazard class as the leap-year forcing bug.
        if (mon_every > 0) then
            if (mod(n, mon_every) == 0) &
                call timer_report(partit%MPI_COMM_FESOM, partit%mype, partit%npes, n, 'monitor')
        end if
        ! M8c Step 3: GLOBAL per-step diagnostics. ALL ranks compute local extrema and reduce with ONE
        ! MPI_MAX; only rank 0 prints. Stability (max|eta|,max|uv|), cryosphere (max a_ice in [0,1], max
        ! m_ice), and warm/salty drift ceilings (global Tmax, Smax). Diagnostic-only: reads the prognostic
        ! state after step_oce, never writes it. GATED (perf): the per-step MPI_Allreduce is a global sync
        ! FESOM2 lacks — skip unless FESOM3_STEP_DIAG is set (the stability run sets it; benchmarks/production do not).
        if (.not. step_diag) cycle
        dstat(1) = maxval(abs(dyn%eta_n(1:nNodO)))
        dstat(2) = maxval(abs(dyn%uv(:,:,1:nElemO)))
        dstat(3) = maxval(ice%data(1)%values(1:nNodO))            ! a_ice  (area fraction, >=0)
        dstat(4) = maxval(ice%data(2)%values(1:nNodO))            ! m_ice  (effective thickness, >=0)
        dstat(5) = maxval(tracers%data(1)%values(:,1:nNodO))      ! T max  (tropical SST ceiling)
        dstat(6) = maxval(tracers%data(2)%values(:,1:nNodO))      ! S max  (evaporative-basin salinity)
        ! Shared elements (docs/plans/completed/20261012-shared-element-sync.md): the copy spreads
        ! max|copy - designated owner's copy| of uv (0 with FESOM3_SHARED_SYNC=1; the drift
        ! without it), uv_rhsAB and sigma (0 by construction, guarded by test_vinv V13 and
        ! test_ice_dyn_mixed D6), and the divergence-sum ratio |sum ssh_rhs_old| /
        ! sum|ssh_rhs_old| over owned nodes: a linfs leak shows as a ratio above ~1e-12
        ! (zstar adds the freshwater term to ssh_rhs_old, so there it is only reported).
        ! Appended AFTER Smax: tools/run_cutout_short.sh reads fields 5, 7 and 9.
        dstat(7) = shared_copy_max_elem(dyn%uv, partit)
        dstat(8) = shared_copy_max_elem(dyn%uv_rhsAB(1,:,:,:), partit)
        dstat(9) = 0.0_WP
        if (allocated(ice%work%sigma11)) dstat(9) = max(shared_copy_max_elem(ice%work%sigma11, partit), &
                                                        shared_copy_max_elem(ice%work%sigma12, partit), &
                                                        shared_copy_max_elem(ice%work%sigma22, partit))
        call MPI_Allreduce(MPI_IN_PLACE, dstat, 9, MPI_WP, MPI_MAX, partit%MPI_COMM_FESOM, ierr)
        dsum = [sum(dyn%ssh_rhs_old(1:nNodO)), sum(abs(dyn%ssh_rhs_old(1:nNodO)))]
        call MPI_Allreduce(MPI_IN_PLACE, dsum, 2, MPI_WP, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
        if (partit%mype == 0) &
            write(*,'(a,i0,10(a,es12.4))') 'fesom_lifecycle_native_mr: step ', n, &
                '  max|eta|=', dstat(1), '  max|uv|=', dstat(2), &
                '  a_ice=', dstat(3), '  m_ice=', dstat(4), &
                '  Tmax=', dstat(5), '  Smax=', dstat(6), &
                '  divsum=', abs(dsum(1))/max(dsum(2), tiny(1.0_WP)), &
                '  dUV=', dstat(7), '  dAB=', dstat(8), '  dsig=', dstat(9)
    end do
    ! Restart Stage 5 (Task 5.3): end-of-run final checkpoint — a safety net for the case the last
    ! in-loop step was NOT cadence-due. restart_due already treats the last step as due, so this is
    ! normally skipped (wrote_final); kept so a run can never end without a resumable checkpoint.
    ! clock_finish after the write keeps RestartOutPath//.clock at the final instant for the next launch.
    if (do_restart .and. .not. wrote_final) then
        cp_year = yearnew; cp_day = daynew; cp_tsec = real(timenew, real64)
        if (daynew == ndpyr .and. timenew == 86400._WP) then
            cp_tsec = 0.0_real64; cp_day = 1; cp_year = yearold + 1
        end if
        call restart_write(rst, trim(restart_dir), cp_year, cp_day, cp_tsec, globalstep=nsteps)
        if (partit%mype == 0) call clock_finish()
        if (partit%mype == 0) write(*,'(a)') 'restart: end-of-run final checkpoint written'
    end if
    if (do_restart) call restart_finalize(rst)
    if (do_output) then
        call means_finalize(oio)
        if (partit%mype == 0) write(*,'(a)') 'M9: field output finalized.'
    end if
    ! permanent end-of-run per-component diagnostics (mean/min/max ms/step + % of loop, ocean broken down)
    call timer_report(partit%MPI_COMM_FESOM, partit%mype, partit%npes, nsteps, 'FINAL')

    call dump_finalize()
    if (partit%mype == 0) &
        write(*,'(a,i0,a)') 'fesom_lifecycle_native_mr: done (', nsteps, ' steps, NATIVE fluxes).'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    ! Restart Stage 5 (Task 5.3): the GLOBAL restart cadence predicate. FESOM2 io_restart.F90:937 is_due
    ! + gen_events.F90 annual/monthly/daily/hourly/step_event ported FRESH (the M9 mod_io_means event_due
    ! is private + coupled to t_means_clock, not reusable). Reads the host's restart_length/unit + the
    ! live mod_clock state. The LAST step is ALWAYS due so end-of-run always leaves a resumable
    ! checkpoint; 'off' disables only the PERIODIC writes. Identical 'y'/'m'/'d'/'h'/'s' semantics to the
    ! M9 output event so a checkpoint and an annual output land on the same instant.
    logical function restart_due(istep)
        integer, intent(in) :: istep
        integer :: fr
        restart_due = .false.
        fr = max(1, restart_length)
        select case (restart_length_unit(1:1))
        case ('y')                                            ! annual_event
            restart_due = (mod(yearnew - yearstart + 1, fr) == 0 .and. &
                           daynew == ndpyr .and. timenew == 86400._WP)
        case ('m')                                            ! monthly_event
            restart_due = (mod(month, fr) == 0 .and. &
                           day_in_month == num_day_in_month(fleapyear, month) .and. timenew == 86400._WP)
        case ('d')                                            ! daily_event
            restart_due = (mod(daynew, fr) == 0 .and. timenew == 86400._WP)
        case ('h')                                            ! hourly_event
            restart_due = (mod(timenew, 3600._WP*real(fr, WP)) == 0._WP)
        case ('s')                                            ! step_event
            restart_due = (mod(istep, fr) == 0)
        case ('o')                                            ! 'off' — no periodic writes
            restart_due = .false.
        case default
            restart_due = .false.
        end select
        if (istep == nsteps) restart_due = .true.             ! end-of-run always checkpoints
    end function restart_due

    ! M9 Task 2.6: register one output stream from a parsed namelist.io row (id + freq/unit/precision/
    ! mean|snap), dispatching the FESOM3 field's fixed name/long_name/units/kind. Unknown ids warn+skip.
    ! 'unod' registers the unod/vnod vector PAIR; a lone 'vnod' row is then a no-op.
    subroutine register_output_var(io, e)
        type(t_io_means), intent(inout) :: io
        type(t_io_entry), intent(in)    :: e
        character(len=8) :: prec
        logical          :: ismean
        integer          :: fr
        character(len=1) :: un
        prec = '<f4'; if (e%precision == 8) prec = '<f8'
        ismean = (trim(e%op) == 'mean')
        fr = e%freq; un = e%unit
        select case (trim(e%id))
        case ('ssh')
            call means_define_node2d(io,'ssh','sea surface elevation','m', &
                 std='sea_surface_height_above_geoid', precision=prec, mean=ismean, freq=fr, unit=un)
        case ('sst')
            call means_define_node2d(io,'sst','sea surface temperature','C', &
                 std='sea_surface_temperature', precision=prec, mean=ismean, freq=fr, unit=un)
        case ('sss')
            call means_define_node2d(io,'sss','sea surface salinity','psu', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('a_ice')
            call means_define_node2d(io,'a_ice','ice concentration','', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('m_ice')
            call means_define_node2d(io,'m_ice','effective ice thickness','m', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('m_snow')
            call means_define_node2d(io,'m_snow','effective snow thickness','m', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('temp')
            call means_define_node3d(io,'temp','sea water potential temperature','C', &
                 on_full_levels=.false., std='sea_water_potential_temperature', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('salt')
            call means_define_node3d(io,'salt','sea water salinity','psu', &
                 on_full_levels=.false., std='sea_water_salinity', &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('w')
            call means_define_node3d(io,'w','vertical velocity','m/s', &
                 on_full_levels=.true., precision=prec, mean=ismean, freq=fr, unit=un)
        case ('unod')
            call means_define_vector3d(io,'unod','vnod','zonal velocity at nodes', &
                 'meridional velocity at nodes','m/s', on_full_levels=.false., &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('vnod')
            continue   ! registered together with 'unod' (the vector pair)
        ! Task 2.7 ELEMENT outputs (user-requested): element velocity u/v (dyn%uv), vertical viscosity
        ! Av (dyn%work%Av, full levels), GM bolus u/v (dyn%fer_uv, Fer_GM only).
        case ('u')
            call means_define_vector3d_elem(io,'u','v','zonal velocity at elements', &
                 'meridional velocity at elements','m/s', on_full_levels=.false., &
                 precision=prec, mean=ismean, freq=fr, unit=un)
        case ('v')
            continue   ! registered together with 'u' (the element vector pair)
        case ('Av')
            call means_define_elem3d(io,'Av','vertical viscosity (momentum)','m2/s', &
                 on_full_levels=.true., precision=prec, mean=ismean, freq=fr, unit=un)
        case ('bolus_u')
            if (use_fer_gm) then
                call means_define_vector3d_elem(io,'bolus_u','bolus_v','GM bolus velocity x', &
                     'GM bolus velocity y','m/s', on_full_levels=.false., &
                     precision=prec, mean=ismean, freq=fr, unit=un)
            else if (partit%mype == 0) then
                write(*,'(a)') 'M9 WARNING: bolus_u/bolus_v requested but Fer_GM off (skipped)'
            end if
        case ('bolus_v')
            continue   ! registered together with 'bolus_u'
        case default
            if (partit%mype == 0) &
                write(*,'(a)') 'M9 WARNING: unknown output var "'//trim(e%id)//'" (skipped)'
        end select
    end subroutine register_output_var

    ! the no-namelist fallback default set (FESOM3_OUTPUT without a namelist.io): the snapshot fields
    ! the env path always produced, on a per-step cadence (freq=every step, unit='s').
    subroutine register_default_outputs(io, every)
        type(t_io_means), intent(inout) :: io
        integer,          intent(in)    :: every
        call means_define_node2d(io,'ssh','sea surface elevation','m', &
             std='sea_surface_height_above_geoid', freq=every, unit='s')
        call means_define_node2d(io,'sst','sea surface temperature','C', &
             std='sea_surface_temperature', freq=every, unit='s')
        call means_define_node2d(io,'sss','sea surface salinity','psu', freq=every, unit='s')
        call means_define_node2d(io,'a_ice','ice concentration','', freq=every, unit='s')
        call means_define_node2d(io,'m_ice','effective ice thickness','m', freq=every, unit='s')
        call means_define_node2d(io,'m_snow','effective snow thickness','m', freq=every, unit='s')
        call means_define_node3d(io,'temp','sea water potential temperature','C', &
             on_full_levels=.false., std='sea_water_potential_temperature', freq=every, unit='s')
        call means_define_node3d(io,'salt','sea water salinity','psu', &
             on_full_levels=.false., std='sea_water_salinity', freq=every, unit='s')
        call means_define_node3d(io,'w','vertical velocity','m/s', on_full_levels=.true., &
             freq=every, unit='s')
        call means_define_vector3d(io,'unod','vnod','zonal velocity at nodes', &
             'meridional velocity at nodes','m/s', on_full_levels=.false., freq=every, unit='s')
        ! Task 2.7 element fields: u/v (dyn%uv) + Av (full levels); bolus only when Fer_GM is on.
        call means_define_vector3d_elem(io,'u','v','zonal velocity at elements', &
             'meridional velocity at elements','m/s', on_full_levels=.false., freq=every, unit='s')
        call means_define_elem3d(io,'Av','vertical viscosity (momentum)','m2/s', &
             on_full_levels=.true., freq=every, unit='s')
        if (use_fer_gm) &
            call means_define_vector3d_elem(io,'bolus_u','bolus_v','GM bolus velocity x', &
                 'GM bolus velocity y','m/s', on_full_levels=.false., freq=every, unit='s')
    end subroutine register_default_outputs

    ! compute the native CORE2 atmosphere over owned+halo (nNodL): the 8 NCAR fields
    ! (timeinterp at the per-step rdate) + the NCAR bulk (Ch/Ce, live srfoce) + stress_atmoce
    ! + the wind-on-ice stress (previous-step uice/vice). runoff (nro) + Ssurf (nss) were read
    ! once at setup. Every routine threads partit -> the OWNED+HALO nodes are filled.
    subroutine compute_native_forcing(n)
        integer, intent(in) :: n                       ! step index -> per-step forcing rdate (rcur)
        real(kind=WP) :: rcur
        ! M8b: refresh interp coefficients on record/day crossings (forcing_sbc_do reads the model
        ! clock yearnew/daynew/timenew), THEN evaluate atmdata at the per-step rdate. rcur MUST equal
        ! forcing_sbc_do's per-field rdate (same formula, same clock; forc_calendar matches the file
        ! calendar — CORE2 'noleap'=>365*yyyy, JRA55 'gregorian') — FESOM2 sbc_do uses one rdate for
        ! both the crossing test and the interp.
        call timer_start(TMR_FRC_SBC)
        call forcing_sbc_do(frc, mesh, partit)
        call timer_stop(TMR_FRC_SBC)
        rcur = real(forcing_julday(yearnew,1,1,forc_calendar),WP) + real(daynew-1,WP) &
             + timenew/86400._WP - dt/86400._WP/2._WP
        call timer_start(TMR_FRC_INTERP)
        call forcing_timeinterp(frc, rcur, partit)
        nuw = frc%atmdata(1,1:nNodL);  nvw = frc%atmdata(2,1:nNodL);  nsh = frc%atmdata(3,1:nNodL)
        nswr = frc%atmdata(4,1:nNodL);  nlw = frc%atmdata(5,1:nNodL)
        nta = frc%atmdata(6,1:nNodL) - 273.15_WP
        npr = frc%atmdata(7,1:nNodL) / 1000._WP
        nps = frc%atmdata(8,1:nNodL) / 1000._WP
        ncd = 0.0_WP; nch = 0.0_WP; nce = 0.0_WP
        call timer_stop(TMR_FRC_INTERP)
        call timer_start(TMR_FRC_BULK)
        call forcing_bulk_ncar(10.0_WP, 10.0_WP, 10.0_WP, nta, nsh, nuw, nvw, &
                               ice%srfoce_temp, ice%srfoce_u, ice%srfoce_v, ncd, nch, nce, mesh, partit)
        call timer_stop(TMR_FRC_BULK)
        call timer_start(TMR_FRC_STRESS)
        call forcing_wind_stress(0.0_WP, nuw, nvw, ice%srfoce_u, ice%srfoce_v, ncd, nsx, nsy, mesh, partit)
        call forcing_ice_stress(0.0012_WP, nuw, nvw, ice%uice, ice%vice, nix, niy, mesh, partit)
        call timer_stop(TMR_FRC_STRESS)
    end subroutine compute_native_forcing

    ! M8c monthly climatology read-ahead (FESOM2 sbc_do:1590-1618). At the last instant of a month
    ! (timenew==86400 on the month's final day) and at step 1, re-read the SSS restoring and (when the
    ! production Sweeney path is on) the chl for the NEXT month. update_monthly_flag + the i=month /
    ! (mstep>1 -> +1) / wrap-at-12 index logic is byte-faithful to the oracle; the read_other_NetCDF
    ! call is the same one the cold start used, so the slice bytes match. Called after the atmosphere
    ! (compute_native_forcing) so the per-step order mirrors sbc_do (atmosphere -> SSS -> chl).
    subroutine roll_monthly_clim(n)
        use mod_clock, only: month, day_in_month, num_day_in_month, fleapyear, timenew
        integer, intent(in) :: n
        logical :: update_monthly
        integer :: i
        update_monthly = ( (day_in_month == num_day_in_month(fleapyear, month) .and. &
                            timenew == 86400._WP) .or. n == 1 )
        if (.not. update_monthly) return
        ! SSS restoring (sss_data_source='CORE2'; always active in this native lifecycle)
        i = month; if (n > 1) i = i + 1; if (i > 12) i = 1
        if (partit%npes == 1) then
            call read_other_NetCDF(trim(sss_file), 'SALT', i, nss, .true., .true., mesh)
        else
            call read_other_NetCDF(trim(sss_file), 'SALT', i, nss, .true., .true., mesh, partit)
        end if
        ! Sweeney chl (production; same month index)
        if (use_sw_pene .and. use_chl_sweeney) then
            i = month; if (n > 1) i = i + 1; if (i > 12) i = 1
            if (partit%npes == 1) then
                call read_other_NetCDF(trim(chl_file), 'chl', i, chl, .true., .true., mesh)
            else
                call read_other_NetCDF(trim(chl_file), 'chl', i, chl, .true., .true., mesh, partit)
            end if
        end if
        if (partit%mype == 0) write(*,'(a,i0,a,i0)') &
            ' roll_monthly_clim: month update -> SSS/chl slice ', i, ' at step ', n
    end subroutine roll_monthly_clim

    ! compute natively + WRITE into atm%* / ice%stress_atmice for the step.
    subroutine apply_native_forcing(n)
        integer, intent(in) :: n
        call compute_native_forcing(n)
        call roll_monthly_clim(n)          ! M8c: SSS + Sweeney chl monthly read-ahead
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
        ice%stress_atmice_x(1:nNodL) = nix
        ice%stress_atmice_y(1:nNodL) = niy
        atm%runoff          = nro
        atm%Ssurf           = nss
    end subroutine apply_native_forcing

    ! M8c bisection: read one F2 atmflux record (fesom_atmflux_dump order) into tscr field-by-field
    ! and print max|delta| of each NATIVE atmosphere/bulk field vs the oracle's, over owned nodes.
    subroutine atm_selfcheck(n)
        integer, intent(in) :: n
        integer :: is, inn, ine, no
        real(kind=WP) :: d_uw, d_vw, d_ta, d_sh, d_sw, d_lw, d_pr, d_ps
        real(kind=WP) :: d_ch, d_ce, d_sx, d_sy, d_ix, d_iy
        read(atmchk_unit) is, inn, ine            ! inn = oracle myDim_nod2D (record array size)
        no = min(inn, nNodO)                      ! compare over owned nodes (np=1 => global order)
        read(atmchk_unit) tscr(1:inn); d_sw = maxval(abs(nswr(1:no) - tscr(1:no)))   ! shortwave
        read(atmchk_unit) tscr(1:inn); d_lw = maxval(abs(nlw (1:no) - tscr(1:no)))   ! longwave
        read(atmchk_unit) tscr(1:inn); d_ta = maxval(abs(nta (1:no) - tscr(1:no)))   ! Tair [degC]
        read(atmchk_unit) tscr(1:inn); d_sh = maxval(abs(nsh (1:no) - tscr(1:no)))   ! shum
        read(atmchk_unit) tscr(1:inn); d_pr = maxval(abs(npr (1:no) - tscr(1:no)))   ! prec_rain [m/s]
        read(atmchk_unit) tscr(1:inn); d_ps = maxval(abs(nps (1:no) - tscr(1:no)))   ! prec_snow [m/s]
        read(atmchk_unit) tscr(1:inn)                                                ! runoff (read once)
        read(atmchk_unit) tscr(1:inn); d_uw = maxval(abs(nuw (1:no) - tscr(1:no)))   ! u_wind
        read(atmchk_unit) tscr(1:inn); d_vw = maxval(abs(nvw (1:no) - tscr(1:no)))   ! v_wind
        read(atmchk_unit) tscr(1:inn); d_ch = maxval(abs(nch (1:no) - tscr(1:no)))   ! Ch_atm_oce
        read(atmchk_unit) tscr(1:inn); d_ce = maxval(abs(nce (1:no) - tscr(1:no)))   ! Ce_atm_oce
        read(atmchk_unit) tscr(1:inn); d_sx = maxval(abs(nsx (1:no) - tscr(1:no)))   ! stress_atmoce_x
        read(atmchk_unit) tscr(1:inn); d_sy = maxval(abs(nsy (1:no) - tscr(1:no)))   ! stress_atmoce_y
        read(atmchk_unit) tscr(1:inn); d_ix = maxval(abs(nix (1:no) - tscr(1:no)))   ! stress_atmice_x
        read(atmchk_unit) tscr(1:inn); d_iy = maxval(abs(niy (1:no) - tscr(1:no)))   ! stress_atmice_y
        read(atmchk_unit) tscr(1:inn)                                                ! Ssurf
        if (partit%mype == 0) then
            write(*,'(a,i0,a,8es10.2)') 'ATMCHK step ', n, &
                ' |d| uw,vw,ta,sh,sw,lw,pr,ps=', d_uw, d_vw, d_ta, d_sh, d_sw, d_lw, d_pr, d_ps
            write(*,'(a,6es10.2)') '            |d| ch,ce,sx,sy,ix,iy=', &
                d_ch, d_ce, d_sx, d_sy, d_ix, d_iy
        end if
    end subroutine atm_selfcheck

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

end program fesom_lifecycle_native_mr
