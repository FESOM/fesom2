program fesom_conserve
    ! CONSERVATION + NO-LEAKAGE GATE for the bottom-at-vertices change.
    !
    ! This is the numerical regression net that replaces the retired FESOM2 byte-gates.
    ! It is a copy of fesom_lifecycle_mr (deliberate: drivers in this tree are duplicated
    ! rather than sharing an init helper) with the mod_dump probes swapped for per-step
    ! invariants that the bottom change would break if it were wrong:
    !
    !   1. CONSERVATION. Total heat and salt content
    !          sum over owned n, nz in [ulevels_nod2D(n), nlevels_nod2D(n)-1] of
    !              tr(nz,n) * hnode(nz,n) * areasvol(nz,n)
    !      must not drift. The run is UNFORCED (all surface fluxes 0) and advection is
    !      conservative, so any drift beyond round-off means the flux areas and the cell
    !      volumes disagree -- which is exactly the risk when area() becomes
    !      depth-independent and nlevels() flips from element-defined to vertex-defined.
    !
    !   2. NO VELOCITY LEAKAGE (spec T10). UV(:,nz,e) must be exactly 0 for
    !      nz >= nlevels(e): velocity DOF in partly-land prisms must never be assembled.
    !
    !   3. THICKNESS CONSISTENCY. helem(nz,e) == sum(elem_w*hnode(nz,elnodes)) over the FULL
    !      element range nz in [ulevels(e), nlevels(e)-1], and zbar_e_bot(e) ==
    !      zbar(nlevels(e)). Catches the zstar helem bottom-layer bound: today the node
    !      hnode commit stops at nlevels_nod2D_min(n)-2 <= nlevels(e)-2 so the element's
    !      deepest layer is never stretched, but after the inversion nlevels(e) <=
    !      nlevels_nod2D(n) and that guarantee is gone.
    !
    !   4. FINITENESS of UV, w, hnode and the tracers, so a NaN is caught at its step
    !      rather than after it has spread.
    !
    ! Runs on dist_<NP> (NP>=1). Report the per-step contents on rank 0; error stop on any
    ! violated invariant. Drift is judged by tools/run_conserve_pi.sh.
    !
    !   FESOM3_MESH_DIR    mesh dir (must contain dist_<NP>/)  default: CORE2
    !   FESOM3_IC_FILE     IC netcdf                           default: pool phc3.0_winter.nc
    !   FESOM3_NSTEPS      number of steps                     default: 3
    !   FESOM3_WHICH_ALE   linfs | zlevel | zstar              default: linfs
    !   FESOM3_CONSERVE_TOL  max |relative drift| before error stop; unset/0 = report only
    !   FESOM3_WSPLIT      presence -> dyn%use_wsplit (smooth w = w_e + w_i split)
    !   FESOM3_WSPLIT_MAXCFL  explicit CFL_z cap (maxcfl)      default: 1.0
    !   FESOM3_WSPLIT_MINCFL  onset of the split (mincfl)      default: 0.5*maxcfl
    !   FESOM3_WSPLIT_EXPECT_SPLIT=N  error stop unless >= N owned faces split (hook comment)
    !   FESOM3_CONSERVE_ICE  presence -> sea ice ON (see below)
    !   FESOM3_WHICHEVP    0=EVP / 1=mEVP (with FESOM3_CONSERVE_ICE)    default: 0
    !   FESOM3_DIVSUM_TOL  shared-element gate (below); unset/0 = report only
    !
    ! SHARED ELEMENTS (docs/plans/completed/20261012-shared-element-sync.md). A border element is owned by
    ! every rank that owns one of its vertices; update_vel syncs uv to the designated owner's
    ! copy (FESOM3_SHARED_SYNC, src/infra/mod_shared_elem.F90). Every step the driver takes the
    ! divergence-sum ratio |sum ssh_rhs_old| / sum|ssh_rhs_old| over the owned nodes (a leak of
    ! the linear free surface shows as a ratio above ~1e-12) and the copy spreads of uv,
    ! uv_rhsAB and sigma11/12/22 (max|copy - designated copy|). No per-step line (run_tri_
    ! regression greps "step 20 "): the run maxima go into ONE end-of-run line. With
    ! FESOM3_DIVSUM_TOL > 0 the run stops with an error if the max ratio exceeds it, or if the
    ! sync is on and any spread is nonzero. With zstar ssh_rhs_old also carries the freshwater
    ! term, so gate the ratio with linfs.
    !
    ! SEA ICE (FESOM3_CONSERVE_ICE, hybrid-mesh Task 17). The ice dynamics (EVP/mEVP) and the
    ! ice TG/FCT advection run every step, coupled to the ocean both ways as in the production
    ! runloop: ocean2ice (surface velocity, hbar) -> EVPdynamics_solve -> ice_TG_rhs ->
    ! ice_fct_solve -> oce_fluxes_mom (ice-ocean drag -> stress_surf, stress_node_surf). The ice
    ! is SEEDED analytically (a cold start puts ice only where SST < 0, i.e. none on a tropical
    ! cut-out): a_ice = 0.45*(1 + sin(k*lon)*cos(k*lat)) with a 5 deg wavelength (geographic
    ! coordinates; ice-free troughs, so the FCT limiter acts), m_ice = a_ice*(1.5 + 0.5*cos(k*lon)),
    ! m_snow = 0.1*a_ice*(1 + sin(k*lat)), driven by a uniform wind-on-ice stress (0.1, -0.05)
    ! N/m^2 (model frame); the atmosphere-ocean stress stays 0. cut_off and the thermodynamics are
    ! NOT run: they exchange heat, salt and water with the ocean (a budget, not a conserved
    ! quantity) and cut_off removes ice by design. The ocean heat/salt/volume gate is unchanged
    ! (the drag moves momentum only), and the ice masses sum(fct_lumpedmass*f) of a_ice, m_ice
    ! and m_snow over owned nodes (the FCT lumped mass, Task 16) are gated with the same
    ! FESOM3_CONSERVE_TOL.
    !
    ! NOTE ON linfs. The linear free surface is NOT tracer-conserving by construction:
    ! hnode is frozen, so the surface vertical advective flux -w*T*area at nzmin is a real
    ! source/sink with no thickness change to balance it. Measured on pi over 20 steps at
    ! the IC: heat drifts -2.2e-04, salt -1.2e-06 (non-monotone -- a free-surface
    ! adjustment transient). Do NOT set FESOM3_CONSERVE_TOL for linfs; run it for the T10 /
    ! helem / finiteness invariants instead. zstar IS conserving (pi, 20 steps: heat
    ! -4.9e-15, salt -1.7e-14 at np=1; 0.0 and -7.8e-15 at np=2) and is the mode to gate on.
    use mpi
    use, intrinsic :: iso_fortran_env, only: int32
    use mod_precision,      only: WP, MP, MPI_WP
    use mod_constants,      only: density_0, pi, rad
    use mod_param_phys,     only: N2smth_h, alpha, theta, shear_splines, N2_splines, rvo_upwind
    use mod_param_phys,     only: mix_coeff_PP, A_ver, K_ver, Kv0_const, mix_scheme_nmb
    use mod_param_phys,     only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_param_phys,     only: Fer_GM, Redi, K_GM_max, K_GM_min, K_GM_bvref, &
                                  K_GM_rampmax, K_GM_rampmin, K_GM_resscalorder, K_GM_cm, &
                                  K_GM_cmin, K_GM_Ktaper, scaling_Ferreira, scaling_Rossby, &
                                  scaling_resolution, scaling_FESOM14, scaling_GMzexp, &
                                  scaling_GINsea, GMzexp_zref, GMzexp_smin
    use mod_param_phys,     only: Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                                  scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    use mod_param_phys,     only: Ricr, concv, visc_sh_limit, diff_sh_limit
    use mod_param_phys,     only: tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, &
                                  tke_kappaM_min, tke_kappaM_max, tke_min, tke_surf_min, &
                                  tke_mxl_choice, tke_only, tke_use_ubound_dirichlet, &
                                  tke_use_lbound_dirichlet, tke_dolangmuir
    use mod_config,         only: use_sw_pene
    use oce_mixing_kpp,     only: oce_mixing_kpp_init
    use oce_mixing_tke,     only: tke_init
    use mod_mesh,           only: t_mesh, MAX_NV
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
    use mod_config,         only: which_ALE
    use mod_halo,           only: allreduce_sum, allreduce_max, shared_copy_max_elem
    use oce_wsplit,         only: wsplit_assert_params
    use mod_ice,            only: t_ice
    use mod_ice_setup,      only: ice_setup
    use mod_ice_dyn,        only: ocean2ice, EVPdynamics_solve
    use mod_ice_fct,        only: ice_TG_rhs, ice_fct_solve
    use mod_ice_thermo,     only: t_atmflux
    use mod_ice_oce_coupling, only: oce_fluxes_mom
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
    integer :: n, nz, nl, nzmin, nzmax, e, tr_num, nsteps, ios, env_len, i_ice
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: Ki(:,:), heat_flux(:), water_flux(:), virtual_salt(:)
    real(kind=WP), allocatable :: relax_salt(:), real_salt_flux(:), stress_surf(:,:)
    real(kind=WP) :: is_nonlinfs
    integer       :: momadv_env   ! FESOM3_MOMADV_OPT (1 = vector invariant, 2 = scalar)
    logical       :: use_wsplit, wsplit_expect_split          ! FESOM3_WSPLIT, _EXPECT_SPLIT
    integer       :: wsplit_expect_min                        ! FESOM3_WSPLIT_EXPECT_SPLIT value
    real(kind=WP) :: wsplit_mincfl_env, wsplit_maxcfl_env     ! FESOM3_WSPLIT_MINCFL/_MAXCFL
    real(kind=WP) :: content0(2), content(2), conserve_tol, drift(2)
    real(kind=WP) :: volume0, volume, vdrift
    logical :: use_fer_gm, use_redi, use_kpp, use_tke
    real(kind=WP), allocatable :: stress_node_surf(:,:)
    logical         :: use_ice                 ! FESOM3_CONSERVE_ICE
    integer         :: whichevp                ! FESOM3_WHICHEVP
    type(t_ice)     :: ice
    type(t_atmflux) :: atm
    real(kind=WP)   :: icemass0(3), icemass(3), idrift(3), klon, kwave
    real(kind=WP)   :: divsum_tol, shmax(4)    ! FESOM3_DIVSUM_TOL; run maxima: ratio, dUV, dAB, dsig
    logical         :: shared_ok
    character(len=24) :: verdict

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/pool/data/AWICM/FESOM2/MESHES_FESOM2.1/core2'
    call get_environment_variable('FESOM3_IC_FILE', ic_file)
    if (len_trim(ic_file) == 0) &
        ic_file = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_winter.nc'
    nsteps = 3
    call get_environment_variable('FESOM3_NSTEPS', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) nsteps
    ! ALE vertical coordinate. 'zstar' routes step_oce through update_stiff_mat_ale +
    ! the vert_vel_ale / update_thickness_ale stretch -- the path invariant 3 exists for.
    call get_environment_variable('FESOM3_WHICH_ALE', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) which_ALE = trim(env)
    conserve_tol = 0.0_WP
    call get_environment_variable('FESOM3_CONSERVE_TOL', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) conserve_tol
    divsum_tol = 0.0_WP
    call get_environment_variable('FESOM3_DIVSUM_TOL', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) then
        read(env, *, iostat=ios) divsum_tol
        if (ios /= 0) error stop 'FESOM3_DIVSUM_TOL: not a real'
    end if
    ! GM bolus / Redi isopycnal diffusion. Redi matters here specifically: it is the only
    ! caller of diff_ver_part_redi_expl, which is the only place tr_xynodes is computed --
    ! and tr_xynodes is one of the two node-averaging denominators the depth-independent
    ! area changed. Without FESOM3_REDI that code path never executes.
    call get_environment_variable('FESOM3_FER_GM', env, length=env_len, status=ios)
    use_fer_gm = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_REDI', env, length=env_len, status=ios)
    use_redi = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_MIX_KPP', env, length=env_len, status=ios)
    use_kpp = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_MIX_TKE', env, length=env_len, status=ios)
    use_tke = (ios == 0 .and. env_len > 0)
    ! Spline shear reconstruction (oce_vert_spline, FESOM's RI_SPLINES analogue). Off by
    ! default the shared producer reproduces the historical two-point form bit for bit,
    ! so gating BOTH settings is what proves the spline path is conservative too.
    call get_environment_variable('FESOM3_SHEAR_SPLINES', env, length=env_len, status=ios)
    shear_splines = (ios == 0 .and. env_len > 0)
    call get_environment_variable('FESOM3_N2_SPLINES', env, length=env_len, status=ios)
    N2_splines = (ios == 0 .and. env_len > 0)
    ! Vector-invariant momentum advection (oce_dyn_vinv). VALUE-based, not presence-based:
    ! momadv_opt is 1 or 2, so a bare presence test could not select between them.
    call get_environment_variable('FESOM3_MOMADV_OPT', env, length=env_len, status=ios)
    momadv_env = 2
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) momadv_env
    ! Upwind face-vorticity blend (oce_dyn_vinv). VALUE-based like FESOM3_MOMADV_OPT.
    call get_environment_variable('FESOM3_RVO_UPWIND', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) rvo_upwind
    if (rvo_upwind < 0.0_WP .or. rvo_upwind > 1.0_WP) &
        error stop 'FESOM3_RVO_UPWIND must be in [0,1]'
    ! Smooth Courant-dependent w = w_e + w_i split (oce_wsplit, compute_Wvel_split): the
    ! part of w above the explicit CFL_z cap is advected by the implicit upstream solves
    ! (momentum TDMA, tracer vertical-diffusion TDMA, FCT adv_tra_vert_impl) instead of
    ! the explicit schemes. Presence-based switch; the two parameters are VALUE-based.
    ! Defaults: the t_dyn cap 1.0 and an onset of 0.5*maxcfl, so a config that sets only
    ! FESOM3_WSPLIT_MAXCFL keeps the type's onset/cap ratio (the gate configs set a cap
    ! far below 1.0 because the 20-step cold start never reaches CFL_z ~ 1). Validated
    ! once here (wsplit_assert_params: maxcfl > 0, 0 <= mincfl <= maxcfl). An unparsable
    ! value is an error, not a silent fall-back to the default (ctests
    ! fesom_conserve_wsplit_bad{params,value}).
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
    ! Gate guard (VALUE-based: the minimum number of split faces): a wsplit config in
    ! which no -- or only a handful of -- owned faces split at the last step, or none
    ! reaches the implicit-dominated range f >= 0.5, proves nothing about the split
    ! consumers, so it must FAIL, not pass.
    call get_environment_variable('FESOM3_CONSERVE_ICE', env, length=env_len, status=ios)
    use_ice = (ios == 0 .and. env_len > 0)
    whichevp = 0
    call get_environment_variable('FESOM3_WHICHEVP', env, length=env_len, status=ios)
    if (ios == 0 .and. env_len > 0) read(env, *, iostat=ios) whichevp
    call get_environment_variable('FESOM3_WSPLIT_EXPECT_SPLIT', env, length=env_len, status=ios)
    wsplit_expect_split = (ios == 0 .and. env_len > 0)
    wsplit_expect_min   = 1
    if (wsplit_expect_split) then
        read(env, *, iostat=ios) wsplit_expect_min
        if (ios /= 0 .or. wsplit_expect_min < 1) &
            error stop 'FESOM3_WSPLIT_EXPECT_SPLIT: not a positive integer'
    end if

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
        write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_conserve: nod2D=', mesh%nod2D, &
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
    dyn%momadv_opt    = momadv_env
    dyn%opt_visc      = 7
    dyn%visc_gamma0   = 0.003_WP
    dyn%visc_gamma1   = 0.1_WP
    dyn%visc_gamma2   = 0.285_WP
    dyn%visc_gamma0_h = 0.0_WP
    dyn%visc_gamma1_h = 0.0_WP
    ! CORE2 production pins use_wsplit=.false. (M1.4/M2.9b precedent); FESOM3_WSPLIT turns
    ! the split on here so the gate can run the w_e/w_i consumers end to end.
    dyn%use_wsplit    = use_wsplit
    dyn%wsplit_mincfl = wsplit_mincfl_env
    dyn%wsplit_maxcfl = wsplit_maxcfl_env
    if (partit%mype == 0 .and. dyn%use_wsplit) write(*,'(a,f8.4,a,f8.4,a)') &
        'fesom_conserve: implicit vertical advection split w = w_e + w_i ENABLED ' // &
        '(use_wsplit; wsplit_mincfl=', dyn%wsplit_mincfl, ', wsplit_maxcfl=', dyn%wsplit_maxcfl, ')'

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
        write(*,'(a,2es12.4,a,2es12.4)') 'fesom_conserve: IC(rank0 owned) T=', &
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
    is_nonlinfs = merge(1.0_WP, 0.0_WP, trim(which_ALE) /= 'linfs')
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

    !===========================================================================
    ! GM / Redi (work_core config), local sizes. Transcribed from fesom_lifecycle's
    ! 1-rank block. Redi rides the GM coupling, so enabling Redi enables Fer_GM too.
    if (use_fer_gm .or. use_redi) then
        Fer_GM = .true.
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
            scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP
            scaling_LDD97 = .false.
            K_hor = 0.0_WP
            allocate(dyn%work%Ki(nl-1, nNodL), dyn%work%fer_tapfac(nl-1, nNodL))
            allocate(dyn%work%neutral_slope(3, nl-1, nNodL))
            allocate(dyn%work%slope_tapered(3, nl-1, nNodL))
            allocate(tracers%work%tr_z(nl, nNodL))
            dyn%work%Ki = 0.0_WP; dyn%work%fer_tapfac = 1.0_WP
            dyn%work%neutral_slope = 0.0_WP; dyn%work%slope_tapered = 0.0_WP
            tracers%work%tr_z = 0.0_WP
        else
            Redi = .false.
        end if
    end if

    !===========================================================================
    ! KPP (FESOM3_MIX_KPP) / TKE (FESOM3_MIX_TKE), local sizes. Ported from the 1-rank
    ! fesom_lifecycle. These matter here because the vertical-diffusion TDMA and the KPP
    ! non-local / shortwave terms are the remaining consumers of area/areasvol, and they
    ! are only reachable with the corresponding scheme switched on.
    if (use_kpp) then
        mix_scheme_nmb = 1
        use_sw_pene    = .false.
        Ricr = 0.3_WP; concv = 1.6_WP
        visc_sh_limit = 5.0e-3_WP; diff_sh_limit = 5.0e-3_WP
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
        dyn%work%caseA = 0.0_WP; dyn%work%ustar = 0.0_WP; dyn%work%Bo = 0.0_WP
        dyn%work%kbl = 0
        if (.not. allocated(stress_node_surf)) then
            allocate(stress_node_surf(2, nNodL)); stress_node_surf = 0.0_WP
        end if
        call oce_mixing_kpp_init(Ricr, concv)
    end if
    if (use_tke) then
        mix_scheme_nmb = 5
        allocate(dyn%work%tke(nl, nNodL), dyn%work%tke_Av(nl, nNodL), &
                 dyn%work%tke_Kv(nl, nNodL))
        dyn%work%tke = 0.0_WP; dyn%work%tke_Av = 0.0_WP; dyn%work%tke_Kv = 0.0_WP
        if (.not. allocated(stress_node_surf)) then
            allocate(stress_node_surf(2, nNodL)); stress_node_surf = 0.0_WP
        end if
        call tke_init(tke_c_k, tke_c_eps, tke_cd, tke_alpha, tke_mxl_min, tke_kappaM_min, &
                      tke_kappaM_max, tke_min, tke_surf_min, tke_mxl_choice, tke_only, &
                      tke_use_ubound_dirichlet, tke_use_lbound_dirichlet, tke_dolangmuir)
    end if
    if (partit%mype == 0) write(*,'(a,l1,a,l1,a,l1,a,l1)') &
        'fesom_conserve: Fer_GM=', Fer_GM, ' Redi=', Redi, ' KPP=', use_kpp, ' TKE=', use_tke

    !===========================================================================
    ! SSH stiffness (built ONCE; dt = CORE2 namelist timestep).
    call init_stiff_mat_ale(mesh, dt, partit)

    !===========================================================================
    ! sea ice (FESOM3_CONSERVE_ICE): setup on the SSH CSR pattern (mass matrix + lumped mass),
    ! the production &ice_dyn config (fesom_lifecycle_native_mr), the analytic seed and the
    ! wind-on-ice stress described in the header. Pointwise formulas over owned+halo nodes, so
    ! the halo needs no exchange.
    if (use_ice) then
        call ice_setup(ice, tracers, mesh, dt, partit)
        ice%whichEVP      = whichevp
        ice%cd_oce_ice    = 0.0055_WP
        ice%delta_min     = 1.0e-11_WP
        ice%ice_diff      = 0.0_WP
        ice%ice_gamma_fct = 0.5_WP
        kwave = 2.0_WP*pi/(5.0_WP*rad)   ! 5 deg wavelength
        do n = 1, nNodL
            if (mesh%ulevels_nod2D(n) > 1) cycle
            klon = kwave*mesh%geo_coord_nod2D(1, n)
            ice%data(1)%values(n) = 0.45_WP*(1.0_WP + sin(klon)*cos(kwave*mesh%geo_coord_nod2D(2, n)))
            ice%data(2)%values(n) = ice%data(1)%values(n)*(1.5_WP + 0.5_WP*cos(klon))
            ice%data(3)%values(n) = 0.1_WP*ice%data(1)%values(n) &
                                    *(1.0_WP + sin(kwave*mesh%geo_coord_nod2D(2, n)))
        end do
        ice%stress_atmice_x = 0.10_WP
        ice%stress_atmice_y = -0.05_WP
        allocate(atm%stress_node_surf(2, nNodL), atm%stress_atmoce_x(nNodL), atm%stress_atmoce_y(nNodL))
        atm%stress_node_surf = 0.0_WP; atm%stress_atmoce_x = 0.0_WP; atm%stress_atmoce_y = 0.0_WP
        if (partit%mype == 0) write(*,'(a,i0,a)') 'fesom_conserve: sea ice ON (whichEVP=', whichevp, &
            '; seeded, dynamics + FCT advection + ice-ocean drag, no thermodynamics)'
    end if

    !===========================================================================
    ! runloop: N steps, state persists (multi-step AB2 evolution). lfirst=(n==1).
    ! Invariants are checked BEFORE the first step (the IC baseline) and after each step.
    call check_invariants(0)
    shmax = 0.0_WP
    content0 = content
    volume0  = volume
    if (use_ice) icemass0 = icemass
    do n = 1, nsteps
        if (use_ice) then
            call ocean2ice(ice, dyn, tracers, mesh, partit)
            call EVPdynamics_solve(ice, mesh, partit)
            call ice_TG_rhs(ice, mesh, partit)
            call ice_fct_solve(ice, mesh, partit)
            call oce_fluxes_mom(ice, atm, stress_surf, mesh, partit)
            if (allocated(stress_node_surf)) stress_node_surf = atm%stress_node_surf
        end if
        if (allocated(stress_node_surf)) then
            call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                          heat_flux, water_flux, virtual_salt, relax_salt, &
                          real_salt_flux, is_nonlinfs, stress_surf, partit, stress_node_surf)
        else
            call step_oce(n, dt, (n == 1), dyn, tracers, mesh, Ki, &
                          heat_flux, water_flux, virtual_salt, relax_salt, &
                          real_salt_flux, is_nonlinfs, stress_surf, partit)
        end if
        call check_invariants(n)
        call shared_element_diag()
    end do

    ! Last-step vertical Courant number + w split statistics (owned faces, all ranks).
    ! Printed for EVERY config: max CFL_z is what FESOM3_WSPLIT_MAXCFL for the gate
    ! configs is chosen from, the split counts their non-vacuity evidence.
    call wsplit_statistics()

    drift(1) = reldrift(content(1), content0(1))
    drift(2) = reldrift(content(2), content0(2))
    vdrift   = reldrift(volume, volume0)
    idrift   = 0.0_WP
    if (use_ice) then
        do i_ice = 1, 3
            idrift(i_ice) = reldrift(icemass(i_ice), icemass0(i_ice))
        end do
    end if
    if (partit%mype == 0) then
        write(*,'(a)') 'fesom_conserve: relative drift over the run'
        write(*,'(a,es24.16)') '  heat: ', drift(1)
        write(*,'(a,es24.16)') '  salt: ', drift(2)
        write(*,'(a,es24.16)') '  vol : ', vdrift
        if (use_ice) write(*,'(a,3es24.16)') '  ice a/m/m_snow: ', idrift
        write(*,'(a,i0,a,a,a)') 'fesom_conserve: done (', nsteps, ' steps, which_ALE=', &
            trim(which_ALE), ').'
    end if
    ! shared elements: one line with the run maxima (and the verdict when gated)
    shared_ok = shmax(1) <= divsum_tol
    if (partit%sync_shared_elem) shared_ok = shared_ok .and. all(shmax(2:4) == 0.0_WP)
    if (partit%mype == 0) then
        if (divsum_tol > 0.0_WP) then
            if (shared_ok) then
                verdict = ' -- OK within tol='
            else
                verdict = ' -- VIOLATED, tol='
            end if
            write(*,'(a,a,4(a,es10.3),a,es9.2)') 'fesom_conserve: shared elements (sync ', &
                trim(merge('ON ', 'OFF', partit%sync_shared_elem)), '): max divsum=', shmax(1), &
                '  dUV=', shmax(2), '  dAB=', shmax(3), '  dsig=', shmax(4), trim(verdict), divsum_tol
        else
            write(*,'(a,a,4(a,es10.3))') 'fesom_conserve: shared elements (sync ', &
                trim(merge('ON ', 'OFF', partit%sync_shared_elem)), '): max divsum=', shmax(1), &
                '  dUV=', shmax(2), '  dAB=', shmax(3), '  dsig=', shmax(4)
        end if
    end if
    if (conserve_tol > 0.0_WP) then
        ! written as .not.(<= tol) so that a NaN drift fails the gate
        if (.not. (abs(drift(1)) <= conserve_tol .and. abs(drift(2)) <= conserve_tol .and. &
                   abs(vdrift)   <= conserve_tol .and. all(abs(idrift) <= conserve_tol))) then
            if (partit%mype == 0) &
                write(*,'(a,es12.4)') 'fesom_conserve: CONSERVATION VIOLATED, tol=', conserve_tol
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            error stop 1
        end if
        if (partit%mype == 0) &
            write(*,'(a,es12.4)') 'fesom_conserve: CONSERVATION OK within tol=', conserve_tol
    end if
    if (divsum_tol > 0.0_WP .and. .not. shared_ok) then
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        error stop 1
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)

contains

    subroutine shared_element_diag()
        ! the step's divergence-sum ratio and shared-copy spreads -> run maxima shmax
        ! (header SHARED ELEMENTS); spreads are local, so one MPI_MAX
        real(kind=WP) :: s(2), d(4)
        integer :: ierr
        s = [sum(dyn%ssh_rhs_old(1:nNodO)), sum(abs(dyn%ssh_rhs_old(1:nNodO)))]
        if (partit%npes > 1) call allreduce_sum(s, partit)
        d(1) = abs(s(1)) / max(s(2), tiny(1.0_WP))
        d(2) = shared_copy_max_elem(dyn%uv, partit)
        d(3) = shared_copy_max_elem(dyn%uv_rhsAB(1,:,:,:), partit)
        d(4) = 0.0_WP
        if (use_ice) d(4) = max(shared_copy_max_elem(ice%work%sigma11, partit), &
                                shared_copy_max_elem(ice%work%sigma12, partit), &
                                shared_copy_max_elem(ice%work%sigma22, partit))
        if (partit%npes > 1) call MPI_Allreduce(MPI_IN_PLACE, d(2:4), 3, MPI_WP, MPI_MAX, &
                                                partit%MPI_COMM_FESOM, ierr)
        shmax = max(shmax, d)
    end subroutine shared_element_diag

    real(kind=WP) function reldrift(a, b)
        real(kind=WP), intent(in) :: a, b
        reldrift = 0.0_WP
        if (abs(b) > 0.0_WP) reldrift = (a - b) / abs(b)
    end function reldrift

    logical function finite(x)
        real(kind=WP), intent(in) :: x
        finite = (abs(x) <= huge(x))     ! false for both NaN and +/-Inf
    end function finite

    subroutine wsplit_statistics()
        ! Vertical Courant number and w = w_e + w_i split statistics over OWNED faces
        ! (levels ulevels_nod2D..nlevels_nod2D) as vert_vel_ale left them on the LAST step
        ! (cfl_z is computed with or without use_wsplit; w_i == 0 everywhere without it).
        ! f = |w_i|/|w| is the implicit share of the face velocity. One line: max CFL_z (what
        ! a cap is chosen from), the split faces w_i /= 0 (by the split function exactly the
        ! faces with CFL_z > wsplit_mincfl, the onset), the faces on the CAPPED branch
        ! CFL_z >= Cu_cut = 2*maxcfl - mincfl (where the explicit CFL_z saturates exactly at
        ! wsplit_maxcfl -- at CFL_z = maxcfl itself a face is still on the bend), the faces
        ! with f >= 0.5 (exactly CFL_z >= 2*maxcfl, the implicit-dominated range) and max f.
        ! MPI-summed/maxed, rank 0 prints. FESOM3_WSPLIT_EXPECT_SPLIT=N: error stop unless
        ! at least N faces split and at least one has f >= 0.5.
        integer       :: j, k
        integer       :: n_split, n_capped, n_half   ! faces with w_i /= 0; CFL_z >= Cu_cut; f >= 0.5
        real(kind=WP) :: cfl, w, cflmax, fmax, cu_cut

        cu_cut = 2.0_WP*dyn%wsplit_maxcfl - dyn%wsplit_mincfl
        n_split = 0; n_capped = 0; n_half = 0; cflmax = 0.0_WP; fmax = 0.0_WP
        do j = 1, nNodO
            if (mesh%nlevels_nod2D(j) <= 0) cycle
            do k = mesh%ulevels_nod2D(j), mesh%nlevels_nod2D(j)
                cfl    = dyn%cfl_z(k,j)
                cflmax = max(cflmax, cfl)
                if (cfl >= cu_cut) n_capped = n_capped + 1
                w = dyn%w(k,j)
                if (dyn%w_i(k,j) /= 0.0_WP) then
                    n_split = n_split + 1
                    if (abs(dyn%w_i(k,j)) >= 0.5_WP*abs(w)) n_half = n_half + 1
                end if
                if (w /= 0.0_WP) fmax = max(fmax, abs(dyn%w_i(k,j))/abs(w))
            end do
        end do
        if (partit%npes > 1) then
            call allreduce_sum(n_split, partit)
            call allreduce_sum(n_capped, partit)
            call allreduce_sum(n_half, partit)
            call allreduce_max(cflmax, partit)
            call allreduce_max(fmax, partit)
        end if
        if (partit%mype == 0) then
            write(*,'(a,es10.3,a,f8.4,a,i0,a,f8.4,a,i0,a,f8.4,a,i0,a,f8.5)') &
                'fesom_conserve: last-step w split (owned faces): max CFL_z=', cflmax, &
                ' ; w_i/=0 (CFL_z > mincfl=', dyn%wsplit_mincfl, '): ', n_split, &
                ' ; capped (CFL_z >= Cu_cut=', cu_cut, '): ', n_capped, &
                ' ; f>=0.5 (CFL_z >= 2*maxcfl=', 2.0_WP*dyn%wsplit_maxcfl, '): ', n_half, &
                ' ; max f=', fmax
        end if
        if (wsplit_expect_split .and. (n_split < wsplit_expect_min .or. n_half == 0)) then
            if (partit%mype == 0) write(*,'(a,i0,a,i0,a,i0,a)') &
                'fesom_conserve: FESOM3_WSPLIT_EXPECT_SPLIT=', wsplit_expect_min, &
                ' but only ', n_split, ' faces split, ', n_half, ' with f >= 0.5 (vacuous config)'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            error stop 1
        end if
    end subroutine wsplit_statistics

    subroutine check_invariants(step)
        ! Contents (1=heat, 2=salt) into `content`, plus the three hard invariants.
        integer, intent(in) :: step
        integer       :: i, j, k, kmin, kmax, elnodes(MAX_NV)
        real(kind=WP) :: acc, hmean, tol, uicemax

        !______________________________________________________________________
        ! 1. total tracer content over OWNED scalar cells (allreduced at npes>1)
        do i = 1, 2
            acc = 0.0_WP
            do j = 1, nNodO
                if (mesh%nlevels_nod2D(j) <= 0) cycle
                do k = mesh%ulevels_nod2D(j), mesh%nlevels_nod2D(j)-1
                    acc = acc + tracers%data(i)%values(k,j) &
                              * real(mesh%hnode(k,j), WP) * real(mesh%areasvol(j), WP)
                end do
            end do
            if (partit%npes > 1) call allreduce_sum(acc, partit)
            content(i) = acc
        end do

        !______________________________________________________________________
        ! 1b. total VOLUME, sum over owned cells of hnode*areasvol. Tracer content can be
        ! conserved while volume is not (an error in h and an opposite one in T would
        ! cancel in the product), so this is a separate check, not a corollary.
        acc = 0.0_WP
        do j = 1, nNodO
            if (mesh%nlevels_nod2D(j) <= 0) cycle
            do k = mesh%ulevels_nod2D(j), mesh%nlevels_nod2D(j)-1
                acc = acc + real(mesh%hnode(k,j), WP) * real(mesh%areasvol(j), WP)
            end do
        end do
        if (partit%npes > 1) call allreduce_sum(acc, partit)
        volume = acc

        !______________________________________________________________________
        ! 2. T10: no velocity below the element's last FULL prism
        do j = 1, nElemO
            if (mesh%nlevels(j) <= 0) cycle
            do k = mesh%nlevels(j), nl-1
                if (dyn%uv(1,k,j) /= 0.0_WP .or. dyn%uv(2,k,j) /= 0.0_WP) then
                    write(*,'(a,i0,a,i0,a,i0,a,i0,a,2es12.4)') &
                        'fesom_conserve: T10 VIOLATED at step ', step, ' rank ', partit%mype, &
                        ' elem ', j, ' nz ', k, ' UV=', dyn%uv(1,k,j), dyn%uv(2,k,j)
                    error stop 1
                end if
            end do
        end do

        !______________________________________________________________________
        ! 3. helem == kite-weighted (elem_w) mean of the vertex hnode over the FULL element
        !    range, and the element bottom depth is the derived one. Relative tolerance: the
        !    weighted sum is not exact even when all vertices agree.
        tol = 1.0e-13_WP
        do j = 1, nElemO
            if (mesh%nlevels(j) <= 0) cycle
            if (abs(real(mesh%zbar_e_bot(j), WP) - real(mesh%zbar(mesh%nlevels(j)), WP)) > 0.0_WP) then
                write(*,'(a,i0,a,i0,a,i0)') 'fesom_conserve: zbar_e_bot VIOLATED at step ', &
                    step, ' rank ', partit%mype, ' elem ', j
                error stop 1
            end if
            elnodes = mesh%elem2D_nodes(:, j)
            kmin = mesh%ulevels(j); kmax = mesh%nlevels(j)-1
            do k = kmin, kmax
                hmean = sum(mesh%elem_w(:, j)*real(mesh%hnode(k, elnodes), WP))
                if (abs(real(mesh%helem(k,j), WP) - hmean) > tol*max(abs(hmean), 1.0_WP)) then
                    write(*,'(a,i0,a,i0,a,i0,a,i0,a,2es24.16)') &
                        'fesom_conserve: helem VIOLATED at step ', step, ' rank ', partit%mype, &
                        ' elem ', j, ' nz ', k, ' helem/mean=', real(mesh%helem(k,j), WP), hmean
                    error stop 1
                end if
            end do
        end do

        !______________________________________________________________________
        ! 1c. sea-ice masses sum(fct_lumpedmass*f) over owned nodes (a_ice, m_ice, m_snow) --
        !     the quantity the ice FCT conserves (Task 16); finite ice velocity and tracers
        if (use_ice) then
            do i = 1, 3
                acc = 0.0_WP
                do j = 1, nNodO
                    acc = acc + ice%work%fct_lumpedmass(j)*ice%data(i)%values(j)
                end do
                if (partit%npes > 1) call allreduce_sum(acc, partit)
                icemass(i) = acc
            end do
            do j = 1, nNodO
                if (.not. finite(ice%uice(j)) .or. .not. finite(ice%vice(j))) then
                    write(*,'(a,i0,a,i0,a,i0)') 'fesom_conserve: non-finite ice velocity at step ', &
                        step, ' rank ', partit%mype, ' node ', j
                    error stop 1
                end if
                if (.not. finite(ice%data(1)%values(j)) .or. .not. finite(ice%data(2)%values(j)) .or. &
                    .not. finite(ice%data(3)%values(j))) then
                    write(*,'(a,i0,a,i0,a,i0)') 'fesom_conserve: non-finite a/m_ice/m_snow at step ', &
                        step, ' rank ', partit%mype, ' node ', j
                    error stop 1
                end if
            end do
            uicemax = maxval(sqrt(ice%uice(1:nNodO)**2 + ice%vice(1:nNodO)**2))
            if (partit%npes > 1) call allreduce_max(uicemax, partit)
        end if

        !______________________________________________________________________
        ! 4. finiteness of the live state (catch a NaN at its own step)
        do j = 1, nElemO
            do k = 1, nl-1
                if (.not. finite(dyn%uv(1,k,j)) .or. .not. finite(dyn%uv(2,k,j))) then
                    write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_conserve: non-finite UV at step ', &
                        step, ' rank ', partit%mype, ' elem ', j, ' nz ', k
                    error stop 1
                end if
            end do
        end do
        do j = 1, nNodO
            do k = 1, nl
                if (.not. finite(dyn%w(k,j))) then
                    write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_conserve: non-finite w at step ', &
                        step, ' rank ', partit%mype, ' node ', j, ' nz ', k
                    error stop 1
                end if
            end do
            do k = 1, nl-1
                if (.not. finite(real(mesh%hnode(k,j), WP)) .or. &
                    .not. finite(tracers%data(1)%values(k,j)) .or. &
                    .not. finite(tracers%data(2)%values(k,j))) then
                    write(*,'(a,i0,a,i0,a,i0,a,i0)') 'fesom_conserve: non-finite hnode/tracer at step ', &
                        step, ' rank ', partit%mype, ' node ', j, ' nz ', k
                    error stop 1
                end if
            end do
        end do

        !______________________________________________________________________
        if (partit%mype == 0) &
            write(*,'(a,i0,a,es24.16,a,es24.16)') 'fesom_conserve: step ', step, &
                '  heat=', content(1), '  salt=', content(2)
        if (use_ice .and. partit%mype == 0) &
            write(*,'(a,i0,a,3es24.16,a,es10.3)') 'fesom_conserve: ice step ', step, &
                '  a/m/m_snow=', icemass, '  max|uice|=', uicemax
    end subroutine check_invariants

end program fesom_conserve
