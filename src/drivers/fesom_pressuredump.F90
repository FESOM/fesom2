program fesom_pressuredump
    ! M2.1 pressure/EOS + M2.2 PGF + M2.3 vel_rhs operator byte-gate driver. Loads pi
    ! 1-rank with the SAME rotation as the FESOM2 pi run (alpha/beta/gamma=50/15/-90,
    ! cyclic 360), computes geometry (incl. M2.3 coriolis), builds the ALE depths
    ! Z_3d_n/zbar_3d_n + node thickness hnode (linfs full cells; M1.2/M1.3 proved
    ! them), PRESCRIBES analytic T/S + UV/eta_n/UV_rhsAB fields (identical formulas to
    ! the FESOM2 oracle src/fesom_pressure_dump.F90, from the byte-identical coords),
    ! then runs oce_pressure_bv (raw + smoothed N^2) -> oce_pgf -> oce_dyn_velrhs and
    ! dumps density_m_rho0 / hpressure / bvfreq / pgf_x/pgf_y / coriolis / uv_rhs* in
    ! the oracle's FADVHDMP format. tools/pressure_diff.py compares for max|delta|=0.
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): the FESOM2-oracle
    ! byte gate is retired (docs/HANDOFF.md "Gate migration"). REMAPPED: the driver's own
    ! node->element means (KPP viscAE, unstable-element count) run over the 4 slots with the
    ! kite weights elem_w like the kernels; it stays a diagnostic dump (roundoff vs FESOM2).
    !
    !   FESOM3_MESH_DIR      mesh dir   (default: pi)
    !   FESOM3_PRESSURE_OUT  out path   (default: pressure_f3.bin)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_constants,    only: density_0
    use mod_param_phys,   only: N2smth_h
    use mod_param_phys,   only: Fer_GM, Redi, Redi_Ktaper, scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_dyn,          only: t_dyn
    use oce_wsplit,       only: wsplit_assert_params
    use oce_pressure_bv,  only: pressure_bv, sw_alpha_beta, compute_sigma_xy, compute_neutral_slope
    use oce_fer_gm,       only: init_Redi_GM, fer_solve_Gamma, fer_gamma2vel
    use mod_param_phys,   only: K_GM_max, K_GM_min, K_GM_bvref, K_GM_rampmax, K_GM_rampmin, &
                                K_GM_resscalorder, K_GM_cm, K_GM_cmin, K_GM_Ktaper, &
                                scaling_Ferreira, scaling_Rossby, scaling_resolution, &
                                scaling_FESOM14, scaling_GMzexp, scaling_GINsea, GMzexp_zref, GMzexp_smin
    use oce_pgf,          only: pressure_force_4_linfs_fullcell, pressure_force_4_zxxxx_shchepetkin
    use oce_dyn_velrhs,   only: compute_vel_rhs
    use oce_dyn_visc,     only: viscosity_filter
    use oce_dyn_ivertvisc, only: impl_vert_visc_ale
    use oce_ssh_rhs,      only: init_stiff_mat_ale, compute_ssh_rhs_ale
    use oce_ssh_solve,    only: solve_ssh_ale
    use oce_ale,          only: update_vel, compute_hbar_ale, update_eta_n, vert_vel_ale
    use oce_ale_mixing_pp, only: oce_mixing_pp
    use oce_vert_spline,   only: compute_shear2
    use oce_mixing_kpp,   only: ri_iwmix, bldepth, oce_mixing_kpp_init, blmix_kpp, enhance
    use mod_constants,    only: density_0_r, vcpw, g
    use mod_param_phys,   only: Ricr, concv
    use mod_config,       only: use_sw_pene
    use oce_mo_conv,      only: mo_convect
    use mod_tracer,       only: t_tracer
    use oce_tracer_grad,  only: tracer_gradient_elements
    use oce_ale_tracer,   only: diff_tracers_ale
    use mod_param_phys,   only: alpha, theta
    use mod_param_phys,   only: mix_coeff_PP, A_ver, K_ver, Kv0_const
    use mod_param_phys,   only: visc_sh_limit, diff_sh_limit
    use mod_param_phys,   only: use_instabmix, instabmix_kv, use_momix, use_windmix
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r1, wr_r2, wr_r3
    implicit none

    ! Pinned shared constant (NOT pi's namelist dt=86400/36=2400 s): the vel_rhs gate
    ! only needs BOTH codes to use the same dt in the dt*(...)/elem_area scaling.
    real(kind=WP), parameter :: dt_velrhs = 1800.0_WP
    ! M2.6 SSH stiffness dt: the FESOM2 oracle built ssh_stiff at ocean_setup (line 140)
    ! with the pi NAMELIST dt = 86400/step_per_day, step_per_day=36 -> 2400 s, BEFORE the
    ! shim overrides dt to dt_velrhs. So the stiffness matrix uses this value; compute it
    ! the SAME way FESOM2 does (gen_model_setup.F90:92) — NOT a 2400.0 literal — so the
    ! -no-prec-div bits agree (the only dt-dependent part of M2.6).
    ! step_per_day is env-configurable (FESOM3_STEP_PER_DAY) so the SAME driver gates
    ! the SSH stiffness on pi (36 -> dt=2400) AND CORE2 (48 -> dt=1800); the oracle built
    ! ssh_stiff at ocean_setup with its namelist step_per_day, so they MUST match.
    real(kind=WP) :: dt_ssh
    integer       :: spd_ssh

    character(len=512) :: mesh_dir, out_path, env_spd
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    integer :: nsw, n, nz, nl, u, nzmin, nzmax, e
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    real(kind=WP), allocatable :: temp(:,:), salt(:,:), density_ref(:,:)
    real(kind=WP), allocatable :: density(:,:), hpressure(:,:), bvfreq(:,:)
    real(kind=WP), allocatable :: bvfreq_raw(:,:), pgf_x(:,:), pgf_y(:,:)
    ! --- M6a-1: Shchepetkin full-free-surface PGF (which_ALE='zlevel', which_pgf='shchepetkin') ---
    real(kind=WP), allocatable :: pgf_x_shchep(:,:), pgf_y_shchep(:,:)
    character(len=256) :: env_shchep
    logical            :: do_shchep
    ! --- M2.3 vel_rhs (Coriolis + AB2 + PGF + SSH gradient) + M2.4 momentum advection ---
    type(t_dyn) :: dyn
    real(kind=WP), allocatable :: uv_in(:,:,:), uv_rhsAB_prev(:,:,:)
    real(kind=WP), allocatable :: uv_rhs_eul(:,:,:), uv_rhs_ab2(:,:,:), uv_rhsAB_cor(:,:,:)
    real(kind=WP), allocatable :: uvnode_rhs_dump(:,:,:)   ! M2.4 momadv intermediate
    ! --- M2.4 biharmonic viscosity (opt_visc=7, visc_filt_bidiff) ---
    real(kind=WP), allocatable :: uv_rhs_visc(:,:,:), u_c_dump(:,:), v_c_dump(:,:)
    integer       :: ed, el(2), ndu, ng0, ng1, ng2
    real(kind=WP) :: du, dumax
    ! --- M2.5 implicit vertical viscosity (TDMA / Thomas solve, impl_vert_visc_ale) ---
    real(kind=WP), allocatable :: Av(:,:), stress_surf(:,:), uv_rhs_ivv(:,:,:)
    integer       :: nwp, nwm
    ! --- M2.6 SSH (stiffness matrix + ssh_rhs + CG solve) ---
    real(kind=WP), allocatable :: ssh_diag(:), ssh_Aeta(:)
    integer       :: n_ssh_iter, row, ni, ne2
    real(kind=WP) :: sum_ssh_rhs, resid_inf
    ! --- M2.7 ALE (linfs) velocity / SSH / thickness-W update ---
    ! Saved copies of the prescribed inputs the M2.7 kernels overwrite in place
    ! (eta_n by update_eta_n; w_e/w_i by compute_Wvel_split) so the M2.3/M2.4/M2.5
    ! INPUT dump records still echo the prescription, not the M2.7 output.
    real(kind=WP), allocatable :: hbar_in(:), eta_n_in(:), w_e_in(:,:), w_i_in(:,:)
    real(kind=WP), allocatable :: uv_upd(:,:,:)
    integer       :: n_split_faces
    ! --- M2.8 PP vertical mixing (oce_mixing_pp -> Kv nodes / Av elements) ---
    real(kind=WP) :: shear, dz_inv, factor, fmin, fmax
    integer       :: nf, nf_big
    ! --- M2.8b mo_convect convective adjustment (Kv/Av floored where bvfreq<0) ---
    ! pp_Kv_dump/pp_Av_dump save the post-PP (pre-mo_convect) coefficients so the M2.8
    ! pp_Kv/pp_Av records still echo the PP output (mo_convect overwrites Kv/Av in place).
    real(kind=WP), allocatable :: pp_Kv_dump(:,:), pp_Av_dump(:,:)
    integer       :: n_unstab_n, n_unstab_e
    ! --- M2.9a tracer-solve assembly (diff_tracers_ale: horiz diff + ALE reconstruct +
    ! implicit vertical-diffusion TDMA, the FIRST consumer of the M2.8 dyn%work%Kv) ---
    type(t_tracer) :: tracers
    real(kind=WP), allocatable :: Ki(:,:), del_ttf_in(:,:), tr_xy(:,:,:)
    ! --- M4a GM/Redi producers (sw_alpha_beta -> sigma_xy -> neutral_slope) ---
    real(kind=WP), allocatable :: sw_alpha(:,:), sw_beta(:,:), sigma_xy(:,:,:)
    real(kind=WP), allocatable :: neutral_slope(:,:,:), slope_tapered(:,:,:), fer_tapfac(:,:)
    ! --- M4b GM diffusivity / streamfunction / bolus velocity (fer_uv/fer_w ride dyn) ---
    real(kind=WP), allocatable :: fer_K(:,:), fer_c(:), fer_scal(:), fer_gamma(:,:,:)
    ! --- M5a-2 KPP ri_iwmix interior mixing (node viscA + per-tracer diffK) ---
    real(kind=WP), allocatable :: viscA_ri(:,:), diffK_ri(:,:,:)
    ! --- M5a-3 KPP boundary-layer depth (prestep dVsq/ustar/Bo + dbsfc + sw_3d + bldepth) ---
    real(kind=WP), allocatable :: dbsfc_kpp(:,:), dVsq_kpp(:,:), sw3d_kpp(:,:)
    real(kind=WP), allocatable :: ustar_kpp(:), bo_kpp(:), hbl_kpp(:), bfsfc_kpp(:)
    real(kind=WP), allocatable :: stable_kpp(:), caseA_kpp(:), hf_kpp(:), wf_kpp(:)
    real(kind=WP), allocatable :: stress_node_surf(:,:)
    integer,       allocatable :: kbl_kpp(:)
    ! --- M5a-4 KPP blmix_kpp + enhance + combine + node->elem viscAE average ---
    real(kind=WP), allocatable :: blmc_kpp(:,:,:), ghats_kpp(:,:), dkm1_kpp(:,:)
    real(kind=WP), allocatable :: viscA_fin(:,:), diffK_fin(:,:,:), viscAE_kpp(:,:)
    real(kind=WP)              :: minmix
    integer                    :: elem, elnodes(MAX_NV)
    real(kind=WP)              :: usurf, vsurf, u_loc, v_loc, sw0
    real(kind=WP), allocatable :: heat_flux(:), water_flux(:), virtual_salt(:), relax_salt(:)
    real(kind=WP), allocatable :: real_salt_flux(:)
    real(kind=WP), allocatable :: t_solved(:,:), s_solved(:,:), del_ttf_t(:,:), del_ttf_s(:,:)
    real(kind=WP) :: is_nonlinfs, kvmin, kvmax
    integer       :: tr_num

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    spd_ssh = 36
    call get_environment_variable('FESOM3_STEP_PER_DAY', env_spd)
    if (len_trim(env_spd) > 0) read(env_spd, *) spd_ssh
    dt_ssh = 86400.0_WP / real(spd_ssh, WP)
    call get_environment_variable('FESOM3_PRESSURE_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'pressure_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_pressuredump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_pressuredump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    ! hnode (per-node layer thickness, linfs full cells): hnode(nz,n)=zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels_nod2D(n), nlevels_nod2D(n)-1]. Proven by the M1.3 gate.
    allocate(mesh%hnode(nl-1, mesh%nod2D)); mesh%hnode = 0.0_MP
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do

    ! helem (per-element layer thickness, linfs full cells): helem(nz,e)=zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels(e), nlevels(e)-1]; the element analog of hnode (M1 advection gate).
    ! zbar_e_bot (partial-cell bottom depth) = zbar(nlevels(e)) for full cells (FESOM2
    ! init_bottom_elem_thickness, use_partial_cell=.false.). Both consumed by M2.5
    ! impl_vert_visc_ale to rebuild the per-column zbar_n/Z_n.
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    allocate(mesh%zbar_e_bot(mesh%elem2D)); mesh%zbar_e_bot = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
        mesh%zbar_e_bot(e) = mesh%zbar(mesh%nlevels(e))
    end do

    ! zbar_3d_n / Z_3d_n: per-node ALE interface/mid depths, built exactly as FESOM2
    ! init_ale (oce_ale.F90:531-566) at the initial state (no cavity, full cells).
    ! Proven field-by-field by the M1.2 gate.
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

    ! Prescribe analytic temperature + salinity (at nodes). MUST match the FESOM2
    ! oracle src/fesom_pressure_dump.F90 exactly:
    !   T(nz,n) = 12.0 + 8.0*cos(lat)*cos(lon) - 0.20*nz + dT_unstable(nz,n)
    !   S(nz,n) = 34.5 + 0.5*sin(2*lon)*cos(lat) + 0.03*nz   (always > 0 for sqrt(s))
    ! Strong horizontal (cos lat/lon) + vertical (nz) structure exercises the EOS,
    ! the top-down hpressure integration, the N^2 difference, and the horizontal
    ! N^2 smoother (which only changes a horizontally-varying field).
    !
    ! M2.8b: dT_unstable = 8.0*max(0,cos(lat)*cos(lon))*min(nz,8)/8 — a warm subsurface
    ! lens in the warm hemisphere whose downward warming (d/dnz > 0 over the top 8 levels)
    ! OVERCOMES the stabilising -0.20*nz (T) and +0.03*nz (S) there, so N^2 goes NEGATIVE
    ! (statically unstable) in a controlled region. This makes the M2.8b mo_convect
    ! convective adjustment (Kv/Av floored to instabmix_kv where bvfreq<0) NON-VACUOUS.
    ! It CASCADES (T -> density -> hpressure -> pgf -> ... -> every downstream field), but
    ! all M2.1-M2.8 records re-verify max|Δ|=0 (both sides use this identical T). Bounded
    ! at +8 C (T <= ~28 C) and S unchanged (>0) so the EOS stays well-posed.
    allocate(temp(nl-1, mesh%nod2D), salt(nl-1, mesh%nod2D), density_ref(nl-1, mesh%nod2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            temp(nz, n) = 12.0_WP + 8.0_WP*cos(lat)*cos(lon) - 0.20_WP*real(nz, WP) &
                        + 8.0_WP*max(0.0_WP, cos(lat)*cos(lon))*min(real(nz,WP),8.0_WP)/8.0_WP
            salt(nz, n) = 34.5_WP + 0.5_WP*sin(2.0_WP*lon)*cos(lat) + 0.03_WP*real(nz, WP)
        end do
    end do
    ! density_ref = density_0 everywhere (use_density_ref=.false. on pi; the
    ! density anomaly subtracts this ARRAY, not the scalar).
    density_ref = density_0

    allocate(density(nl-1, mesh%nod2D), hpressure(nl, mesh%nod2D), bvfreq(nl, mesh%nod2D))
    allocate(bvfreq_raw(nl, mesh%nod2D))
    allocate(dbsfc_kpp(nl, mesh%nod2D))   ! M5a-3: KPP buoyancy-diff-wrt-surface (filled by pressure_bv)
    dbsfc_kpp = 0.0_WP

    ! Pass 1: horizontal N^2 smoothing OFF -> raw bvfreq (+ density_m_rho0, hpressure).
    ! The caller pre-zeros the outputs (pressure_bv leaves below-bottom entries as-is).
    N2smth_h = .false.
    density = 0.0_WP; hpressure = 0.0_WP; bvfreq = 0.0_WP
    call pressure_bv(temp, salt, density_ref, mesh, density, hpressure, bvfreq)
    bvfreq_raw = bvfreq

    ! Pass 2: horizontal N^2 smoothing ON -> smoothed bvfreq. dbsfc (M5a-3) is density-only
    ! (smoothing-independent), filled here for the KPP bldepth gate below.
    N2smth_h = .true.
    bvfreq = 0.0_WP
    call pressure_bv(temp, salt, density_ref, mesh, density, hpressure, bvfreq, dbsfc=dbsfc_kpp)

    ! ============== M4a GM/Redi producers (feed-forward diagnostics) ==============
    ! sw_alpha_beta -> compute_sigma_xy -> compute_neutral_slope, from T/S + the smoothed
    ! bvfreq. Force the work_core GM/Redi taper config so compute_neutral_slope takes the
    ! PRODUCTION branch (slope_tapered = ns*sqrt(c1*c2), fer_tapfac = c1*c2); the oracle shim
    ! forces the SAME flags + allocates fer_tapfac. Restored after so the downstream reduced-M2
    ! step (28-field gate) is unperturbed. ODM95_Scr=0.2e-2 is the work_core override (not 1e-2).
    Fer_GM = .true.; Redi = .true.; Redi_Ktaper = .true.
    scaling_ODM95 = .true.; ODM95_Scr = 0.2e-2_WP; ODM95_Sd = 1.0e-3_WP; scaling_LDD97 = .false.
    allocate(sw_alpha(nl-1, mesh%nod2D), sw_beta(nl-1, mesh%nod2D))
    allocate(sigma_xy(2, nl-1, mesh%nod2D))
    allocate(neutral_slope(3, nl-1, mesh%nod2D), slope_tapered(3, nl-1, mesh%nod2D))
    allocate(fer_tapfac(nl-1, mesh%nod2D))
    sw_alpha = 0.0_WP; sw_beta = 0.0_WP; sigma_xy = 0.0_WP
    neutral_slope = 0.0_WP; slope_tapered = 0.0_WP; fer_tapfac = 0.0_WP
    call sw_alpha_beta(temp, salt, mesh, sw_alpha, sw_beta)
    call compute_sigma_xy(temp, salt, sw_alpha, sw_beta, mesh, sigma_xy)
    call compute_neutral_slope(sigma_xy, bvfreq, mesh, neutral_slope, slope_tapered, fer_tapfac)
    Fer_GM = .false.; Redi = .false.; Redi_Ktaper = .false.   ! restore reduced-M2 for downstream

    ! M2.2: hydrostatic PGF from the M2.1 hpressure (the smoothing pass leaves
    ! hpressure unchanged - it only touches bvfreq, L13). Caller pre-zeros pgf_x/pgf_y
    ! (oce_pgf writes only levels ule..nle; FESOM2 leaves below-bottom at alloc-zero).
    allocate(pgf_x(nl-1, mesh%elem2D), pgf_y(nl-1, mesh%elem2D))
    pgf_x = 0.0_WP; pgf_y = 0.0_WP
    call pressure_force_4_linfs_fullcell(hpressure, mesh, pgf_x, pgf_y)

    ! M6a-1: the full-free-surface (which_ALE='zlevel') Shchepetkin PGF, run on the
    ! SAME prescribed rest state (eta=0 -> Z_3d_n=initial, so non-vacuous for the
    ! ALGORITHM regardless of eta). Env-gated (FESOM3_PGF_SHCHEP=1) and dumped under
    ! separate field names so the proven linfs pgf_x/pgf_y gate is byte-unaffected.
    call get_environment_variable('FESOM3_PGF_SHCHEP', env_shchep)
    do_shchep = (len_trim(env_shchep) > 0)
    allocate(pgf_x_shchep(nl-1, mesh%elem2D), pgf_y_shchep(nl-1, mesh%elem2D))
    pgf_x_shchep = 0.0_WP; pgf_y_shchep = 0.0_WP
    if (do_shchep) then
        call pressure_force_4_zxxxx_shchepetkin(density, mesh, pgf_x_shchep, pgf_y_shchep)
    end if

    ! ============== M2.3 vel_rhs assembly + M2.4 momentum advection ==============
    ! Coriolis + AB2 + PGF + SSH-gradient + momentum advection (momadv_opt=2). Build a
    ! t_dyn, prescribe analytic UV (elements) / eta_n (nodes) / previous-step UV_rhsAB +
    ! w_e (nodes, the explicit vertical velocity momadv reads; all identical to the
    ! FESOM2 oracle, from the byte-identical coords), copy in the M2.2 pgf_x/pgf_y, and
    ! run oce_dyn_velrhs TWICE: lfirst=.true. (first Euler step, ff=1.0) then
    ! lfirst=.false. (AB2 step, ff=ab2=1.6). The oracle drives the REAL FESOM2
    ! compute_vel_rhs the same way (now with momadv_opt=2 -> momentum_adv_scalar).
    allocate(dyn%uv(2, nl-1, mesh%elem2D), dyn%uv_rhs(2, nl-1, mesh%elem2D))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, mesh%elem2D))     ! (AB_order-1, 2, nl-1, elem2D)
    allocate(dyn%eta_n(mesh%nod2D))
    allocate(dyn%w_e(nl, mesh%nod2D))                   ! explicit vertical velocity (nodes)
    allocate(dyn%w_i(nl, mesh%nod2D))                   ! M2.5 implicit vertical velocity (nodes)
    allocate(dyn%work%uvnode_rhs(2, nl-1, mesh%nod2D))  ! momadv nodal scratch
    allocate(dyn%work%vorticity(nl-1, mesh%nod2D))   ! momadv_opt==1 (vector-invariant) scratch
    dyn%work%vorticity = 0.0_WP
    allocate(dyn%work%pgf_x(nl-1, mesh%elem2D), dyn%work%pgf_y(nl-1, mesh%elem2D))
    allocate(dyn%work%u_c(nl-1, mesh%elem2D), dyn%work%v_c(nl-1, mesh%elem2D))  ! visc scratch
    allocate(Av(nl, mesh%elem2D), stress_surf(2, mesh%elem2D))  ! M2.5 prescribed inputs
    ! M2.8 PP mixing: nodal velocity (prescribed input), smoothed N^2 (copied from the
    ! M2.1 local bvfreq), and the Kv/Av output fields.
    allocate(dyn%uvnode(2, nl-1, mesh%nod2D))                   ! nodal velocity (PP shear input)
    allocate(dyn%work%bvfreq(nl, mesh%nod2D))                   ! smoothed N^2 (PP reads it)
    allocate(dyn%work%shear2(nl, mesh%nod2D))                   ! shared |du/dz|^2 (PP reads it)
    allocate(dyn%work%Kv(nl, mesh%nod2D), dyn%work%Av(nl, mesh%elem2D))
    dyn%AB_order          = 2
    dyn%momadv_opt        = 2     ! M2.4: enable momentum_adv_scalar (pi production value)
    ! M2.4 biharmonic viscosity (opt_visc=7) — pin the pi/reduced-M2 namelist values
    ! (visc_gamma0=0.003 OVERRIDES the t_dyn default 0.03; gamma_h=0 -> pure biharmonic).
    ! MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    dyn%opt_visc          = 7
    dyn%visc_gamma0       = 0.003_WP
    dyn%visc_gamma1       = 0.1_WP
    dyn%visc_gamma2       = 0.285_WP
    dyn%visc_gamma0_h     = 0.0_WP
    dyn%visc_gamma1_h     = 0.0_WP
    dyn%w_e               = 0.0_WP
    dyn%w_i               = 0.0_WP
    dyn%work%uvnode_rhs   = 0.0_WP
    dyn%work%u_c          = 0.0_WP
    dyn%work%v_c          = 0.0_WP
    dyn%work%pgf_x = pgf_x        ! reuse the proven M2.2 PGF (exact copy)
    dyn%work%pgf_y = pgf_y

    ! analytic SSH + explicit vertical velocity w_e at nodes (w_e: sign varies in
    ! space AND depth via cos(0.3*nz) -> non-trivial w*du/dz; MUST match the oracle).
    ! M2.5 also prescribes the IMPLICIT vertical velocity w_i (nodes): a DISTINCT analytic
    ! formula (sin(lon)*cos(2*lat)*cos(0.25*nz)) so its sign varies in space AND depth ->
    ! exercises both the min(0,wu/wd) and max(0,wu/wd) upwind branches of the vertical-
    ! advection update in impl_vert_visc_ale. MUST match the FESOM2 oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        dyn%eta_n(n) = 0.5_WP*cos(lat)*sin(lon) + 0.3_WP*sin(2.0_WP*lat)
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz,WP))
            dyn%w_i(nz,n) = 2.0e-4_WP*sin(lon)*cos(2.0_WP*lat)*cos(0.25_WP*real(nz,WP))
        end do
    end do
    ! analytic element velocity + previous-step AB array (from element node-1 coords).
    ! Horizontal amplitude 2.00/1.50 m/s (was 0.10/0.08 at M2.3/M2.4): a strong-current
    ! stress test so the across-edge velocity jump |du| spans ALL THREE biharmonic-
    ! viscosity flow-aware branches of max(gamma0, max(gamma1*|du|, gamma2*|du|^2)). With
    ! the pi gammas the winners cross at |du|=gamma0/gamma1=0.03 (gamma0->gamma1) and
    ! |du|=gamma1/gamma2=0.35 (gamma1->gamma2), so max|du|~0.71 (~2x the upper crossover)
    ! exercises gamma0, gamma1 AND gamma2 with a clear margin (L11/L17). The depth terms
    ! (-0.005*nz / +0.004*nz) cancel in du=UV(el1)-UV(el2). M2.3/M2.4 re-gate with the new
    ! values (identical formula on both sides). MUST match the FESOM2 oracle.
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1,e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1,e))
        do nz = 1, nl-1
            dyn%uv(1,nz,e) =  2.00_WP*cos(lat)*sin(lon)        - 0.005_WP*real(nz,WP)
            dyn%uv(2,nz,e) = -1.50_WP*sin(lat)*cos(2.0_WP*lon) + 0.004_WP*real(nz,WP)
            dyn%uv_rhsAB(1,1,nz,e) =  1.0e6_WP*sin(lon)*cos(lat)        + 1.0e4_WP*real(nz,WP)
            dyn%uv_rhsAB(1,2,nz,e) = -1.0e6_WP*cos(lon)*sin(2.0_WP*lat) - 1.0e4_WP*real(nz,WP)
        end do
        ! M2.5 prescribed inputs (elements): Av strictly POSITIVE vertical viscosity
        ! ([1e-3,1.4e-2] m^2/s, realistic) over all nz=1..nl (the TDMA reads Av(nzmax)),
        ! stress_surf sign-varying wind stress (~0.1 N/m^2 -> exercises the surface BC for
        ! both signs). MUST match the FESOM2 oracle.
        do nz = 1, nl
            Av(nz,e) = 5.0e-3_WP + 4.0e-3_WP*cos(lat)*cos(lon) + 1.0e-4_WP*real(nz,WP)
        end do
        stress_surf(1,e) =  0.10_WP*cos(lat)*sin(lon)
        stress_surf(2,e) = -0.08_WP*sin(lat)*cos(2.0_WP*lon)
    end do

    allocate(uv_in(2,nl-1,mesh%elem2D), uv_rhsAB_prev(2,nl-1,mesh%elem2D))
    allocate(uv_rhs_eul(2,nl-1,mesh%elem2D), uv_rhs_ab2(2,nl-1,mesh%elem2D))
    allocate(uv_rhsAB_cor(2,nl-1,mesh%elem2D), uvnode_rhs_dump(2,nl-1,mesh%nod2D))
    allocate(uv_rhs_visc(2,nl-1,mesh%elem2D))
    allocate(u_c_dump(nl-1,mesh%elem2D), v_c_dump(nl-1,mesh%elem2D))
    uv_in         = dyn%uv
    uv_rhsAB_prev = dyn%uv_rhsAB(1,:,:,:)

    dyn%uv_rhs = 0.0_WP
    call compute_vel_rhs(dyn, mesh, dt_velrhs, lfirst=.true.)    ! Euler start (ff=1.0)
    uv_rhs_eul      = dyn%uv_rhs
    uv_rhsAB_cor    = dyn%uv_rhsAB(1,:,:,:)         ! Coriolis + momentum advection (M2.4)
    uvnode_rhs_dump = dyn%work%uvnode_rhs          ! momadv nodal scratch (post-normalize)
    call compute_vel_rhs(dyn, mesh, dt_velrhs, lfirst=.false.)   ! AB2 step (ff=ab2)
    uv_rhs_ab2   = dyn%uv_rhs

    ! M2.4 biharmonic viscosity (opt_visc=7): a SEPARATE operator run AFTER
    ! compute_vel_rhs (FESOM2 oce_ale.F90:3822). It ADDS the biharmonic increment into
    ! dyn%uv_rhs (reads dyn%uv only; the incoming uv_rhs=uv_rhs_ab2 is just accumulated
    ! into). The oracle drives the REAL FESOM2 visc_filt_bidiff the same way.
    call viscosity_filter(7, dyn, mesh, dt_velrhs)
    uv_rhs_visc = dyn%uv_rhs                        ! post-viscosity UV_rhs (gate target)
    u_c_dump    = dyn%work%u_c                       ! first-stage Laplacian (intermediate)
    v_c_dump    = dyn%work%v_c

    ! gate-strength diagnostic (NOT dumped): over interior edges x levels, which term
    ! WINS the viscosity coefficient max(gamma0, max(gamma1*|du|, gamma2*|du|^2))? A
    ! non-zero share for each of gamma0/gamma1/gamma2 confirms all three flow-aware
    ! branches of visc_filt_bidiff are genuinely selected (L11/L17), not just computed.
    dumax = 0.0_WP; ndu = 0; ng0 = 0; ng1 = 0; ng2 = 0
    do ed = 1, mesh%edge2D
        if (ed > mesh%edge2D_in) cycle
        el    = mesh%edge_tri(:, ed)
        nzmin = maxval(mesh%ulevels(el))
        nzmax = minval(mesh%nlevels(el))
        do nz = nzmin, nzmax-1
            du = sqrt((dyn%uv(1,nz,el(1))-dyn%uv(1,nz,el(2)))**2 &
                    + (dyn%uv(2,nz,el(1))-dyn%uv(2,nz,el(2)))**2)
            dumax = max(dumax, du); ndu = ndu + 1
            if (0.003_WP >= max(0.1_WP*du, 0.285_WP*du*du)) then
                ng0 = ng0 + 1
            else if (0.1_WP*du >= 0.285_WP*du*du) then
                ng1 = ng1 + 1
            else
                ng2 = ng2 + 1
            end if
        end do
    end do
    write(*,'(a,es10.3,a,f5.1,a,f5.1,a,f5.1,a)') &
        'fesom_pressuredump: visc strength: max|du|=', dumax, &
        ' m/s ; gamma0/gamma1/gamma2 selected on ', &
        100.0_WP*real(ng0,WP)/real(max(ndu,1),WP), '/', &
        100.0_WP*real(ng1,WP)/real(max(ndu,1),WP), '/', &
        100.0_WP*real(ng2,WP)/real(max(ndu,1),WP), '% of edge-levels'

    ! ============== M2.5 implicit vertical viscosity (TDMA / Thomas solve) ==============
    ! A SEPARATE operator run AFTER viscosity_filter in the timestep (FESOM2
    ! oce_ale.F90:3874-3876). It CONSUMES the post-viscosity uv_rhs_visc as the explicit
    ! rhs, solves the per-column tridiagonal (implicit vertical viscosity Av + vertical
    ! advection w_i + wind-stress top BC + bottom drag), and OVERWRITES dyn%uv_rhs with the
    ! solution. The oracle drives the REAL FESOM2 impl_vert_visc_ale on the same inputs.
    allocate(uv_rhs_ivv(2,nl-1,mesh%elem2D))
    call impl_vert_visc_ale(dyn, mesh, dt_velrhs, Av, stress_surf)
    uv_rhs_ivv = dyn%uv_rhs                           ! post-TDMA UV_rhs (gate target)

    ! gate-strength diagnostic (NOT dumped): the TDMA must non-trivially change the rhs,
    ! and w_i must take BOTH signs so both upwind branches of the advection update fire.
    nwp = count(dyn%w_i >  0.0_WP)
    nwm = count(dyn%w_i <  0.0_WP)
    write(*,'(a,es10.3,a,es10.3,a,i0,a,i0)') &
        'fesom_pressuredump: ivertvisc: max|uv_rhs_ivv|=', maxval(abs(uv_rhs_ivv)), &
        ' ; max|d(uv_rhs)|=', maxval(abs(uv_rhs_ivv - uv_rhs_visc)), &
        ' ; w_i>0/<0 = ', nwp, '/', nwm

    ! ============== M2.6 SSH: stiffness matrix + ssh_rhs + CG solve ===================
    ! The free-surface implicit solve, run AFTER impl_vert_visc_ale and BEFORE the
    ! velocity/SSH update in the timestep (FESOM2 oce_ale.F90:3920-3930). For which_ale=
    ! 'linfs' the stiffness matrix is built ONCE (init_stiff_mat_ale; update_stiff_mat_ale
    ! is skipped, oce_ale.F90:3921), then compute_ssh_rhs_ale assembles the rhs from the
    ! prescribed UV + post-TDMA UV_rhs (dyn%uv_rhs = uv_rhs_ivv) and solve_ssh_ale CG-solves
    ! for d_eta. The FESOM2 oracle drives the REAL routines on the same inputs.
    alpha = 1.0_WP; theta = 1.0_WP      ! pi/reduced-M2: full implicitness (FESOM2 default)
    allocate(dyn%ssh_rhs(mesh%nod2D), dyn%ssh_rhs_old(mesh%nod2D), dyn%d_eta(mesh%nod2D))
    dyn%ssh_rhs_old = 0.0_WP            ! (1-alpha)*ssh_rhs_old = 0 since alpha=1
    dyn%d_eta       = 0.0_WP            ! CG initial guess x0 = 0

    call init_stiff_mat_ale(mesh, dt_ssh)              ! build ssh_stiff (sparsity + values)
    call compute_ssh_rhs_ale(dyn, mesh)               ! assemble dyn%ssh_rhs
    call solve_ssh_ale(dyn, mesh, n_iter=n_ssh_iter)  ! CG -> dyn%d_eta

    ! gate-strength diagnostics (NOT dumped):
    !  (1) the stiffness diagonal + a full matvec A*eta_n (DUMPED, gated) localise the
    !      matrix assembly; computed here only for the prints below,
    !  (2) sum(ssh_rhs) telescopes to ~1e-13 (a proper edge-divergence; plan M2.6),
    !  (3) the CG converged in n_ssh_iter iters and ||A*d_eta - ssh_rhs||_inf shows d_eta
    !      actually solves the system (non-vacuity, the L11 weak-gate guard).
    allocate(ssh_diag(mesh%nod2D), ssh_Aeta(mesh%nod2D))
    do row = 1, mesh%nod2D
        ni  = mesh%ssh_stiff%rowptr_loc(row)
        ne2 = mesh%ssh_stiff%rowptr_loc(row+1) - 1
        ssh_diag(row) = mesh%ssh_stiff%values(ni)
        ssh_Aeta(row) = sum(mesh%ssh_stiff%values(ni:ne2) * dyn%eta_n(mesh%ssh_stiff%colind_loc(ni:ne2)))
    end do
    sum_ssh_rhs = sum(dyn%ssh_rhs(1:mesh%nod2D))
    resid_inf   = 0.0_WP
    do row = 1, mesh%nod2D
        ni  = mesh%ssh_stiff%rowptr_loc(row)
        ne2 = mesh%ssh_stiff%rowptr_loc(row+1) - 1
        resid_inf = max(resid_inf, abs( &
            sum(mesh%ssh_stiff%values(ni:ne2) * dyn%d_eta(mesh%ssh_stiff%colind_loc(ni:ne2))) - dyn%ssh_rhs(row)))
    end do
    write(*,'(a,i0,a,i0,a,es10.3,a,es10.3,a,es10.3)') &
        'fesom_pressuredump: ssh: nza=', mesh%ssh_stiff%nza, ' ; CG iters=', n_ssh_iter, &
        ' ; max|d_eta|=', maxval(abs(dyn%d_eta)), ' ; sum(ssh_rhs)=', sum_ssh_rhs, &
        ' ; ||A d_eta - rhs||inf=', resid_inf

    ! ============== M2.7 ALE (linfs) velocity / SSH / thickness-W update ==============
    ! The post-CG tail of the timestep (FESOM2 oce_ale.F90:3946-4081), run AFTER the SSH
    ! solve. Build the ALE thickness state the update reads/writes: mesh%hbar (prescribed,
    ! the previous-step elevation), hbar_old/dhe (outputs), hnode_new=hnode (linfs never
    ! evolves it), and dyn%w/cfl_z. Force the pi wsplit config (use_wsplit=.true.,
    ! wsplit_maxcfl=1.0; wsplit_mincfl at its type default 0.5) so compute_Wvel_split runs
    ! here (the w_i /= 0 count below is its evidence). The split is FESOM3's smooth
    ! Shchepetkin function (oce_wsplit), NOT FESOM2's hard switch, so w_e/w_i are not
    ! byte-comparable and are NOT dumped; FESOM2's w_split_e/w_split_i records are listed
    ! in pressure_diff.py's SKIP_FIELDS (the FESOM2 byte-gates are retired since
    ! bottom-at-vertices anyway).
    ! alpha=theta=1.0 already set above. The chain is, in order:
    !   update_vel (UV += UV_rhs + SSH-grad) -> compute_hbar_ale (hbar/dhe/ssh_rhs_old)
    !   -> update_eta_n (eta_n=hbar) -> vert_vel_ale (w + cfl_z + Wvel split). dyn%uv_rhs
    ! is the post-TDMA uv_rhs_ivv; dyn%uv is the prescribed UV (update_vel overwrites it).
    allocate(mesh%hbar(mesh%nod2D), mesh%hbar_old(mesh%nod2D), mesh%dhe(mesh%elem2D))
    allocate(mesh%hnode_new(nl-1, mesh%nod2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%cfl_z(nl, mesh%nod2D))
    mesh%hbar_old  = 0.0_MP
    mesh%dhe       = 0.0_MP
    mesh%hnode_new = mesh%hnode            ! linfs: hnode_new == hnode (never evolves)
    dyn%w          = 0.0_WP
    dyn%cfl_z      = 0.0_WP
    dyn%use_wsplit    = .true.             ! pi production value (namelist.dyn)
    dyn%wsplit_maxcfl = 1.0_WP             ! the explicit CFL_z cap; wsplit_mincfl = 0.5 (type default)
    call wsplit_assert_params(dyn%wsplit_mincfl, dyn%wsplit_maxcfl)   ! the t_dyn obligation

    ! analytic previous-step elevation hbar (sign-varying, ~0.6 m). MUST match the oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        mesh%hbar(n) = 0.4_WP*sin(lon)*cos(lat) - 0.2_WP*cos(2.0_WP*lat)
    end do

    ! save the prescribed inputs the M2.7 kernels overwrite in place
    allocate(hbar_in(mesh%nod2D), eta_n_in(mesh%nod2D))
    allocate(w_e_in(nl,mesh%nod2D), w_i_in(nl,mesh%nod2D), uv_upd(2,nl-1,mesh%elem2D))
    hbar_in  = mesh%hbar
    eta_n_in = dyn%eta_n
    w_e_in   = dyn%w_e
    w_i_in   = dyn%w_i

    call update_vel(dyn, mesh, dt_velrhs)           ! UV += UV_rhs + [-g*theta*dt*grad(d_eta)]
    uv_upd = dyn%uv
    call compute_hbar_ale(dyn, mesh, dt_velrhs)     ! ssh_rhs_old, hbar_old, hbar, dhe
    call update_eta_n(dyn, mesh)                    ! eta_n = alpha*hbar + (1-alpha)*hbar_old

    ! ============== M4b GM diffusivity + streamfunction + bolus velocity ==============
    ! init_Redi_GM (fer_K/fer_c/fer_scal) -> fer_solve_Gamma (fer_gamma TDMA) ->
    ! fer_gamma2vel (fer_uv), then vert_vel_ale computes fer_w (vertical divergence of fer_uv)
    ! with Fer_GM on. Force the work_core GM config (Fer_GM=T, Redi=F isolating gate;
    ! scaling_resolution+GMzexp, K_GM_max=1000, K_GM_cm=3, K_GM_cmin=0.1, ramp off, Ktaper off).
    ! Inits MUST match FESOM2 oce_setup_step.F90:962-967 (fer_K=500/fer_c=1/fer_scal=0/
    ! fer_gamma=0) + :670/:688 (fer_uv/fer_w=0) so the below-bottom / unwritten regions
    ! byte-match. init_Redi_GM/fer_solve_Gamma read the M4a sigma_xy + the smoothed bvfreq.
    ! Restored after vert_vel_ale so the downstream reduced-M2 M2.8/M2.9 step is unperturbed.
    Fer_GM = .true.; Redi = .false.
    K_GM_max = 1000.0_WP; K_GM_min = 2.0_WP; K_GM_bvref = 1
    K_GM_rampmax = -1.0_WP; K_GM_rampmin = -1.0_WP; K_GM_resscalorder = 2.0_WP
    K_GM_cm = 3.0_WP; K_GM_cmin = 0.1_WP; K_GM_Ktaper = .false.
    scaling_Ferreira = .false.; scaling_Rossby = .false.; scaling_resolution = .true.
    scaling_FESOM14 = .false.; scaling_GMzexp = .true.; scaling_GINsea = .false.
    GMzexp_zref = 500.0_WP; GMzexp_smin = 0.6_WP
    allocate(fer_K(nl, mesh%nod2D), fer_c(mesh%nod2D), fer_scal(mesh%nod2D))
    allocate(fer_gamma(2, nl, mesh%nod2D))
    allocate(dyn%fer_uv(2, nl-1, mesh%elem2D), dyn%fer_w(nl, mesh%nod2D))
    fer_K = 500.0_WP; fer_c = 1.0_WP; fer_scal = 0.0_WP; fer_gamma = 0.0_WP
    dyn%fer_uv = 0.0_WP; dyn%fer_w = 0.0_WP
    call init_Redi_GM(mesh, bvfreq, fer_K, fer_c, fer_scal)
    call fer_solve_Gamma(mesh, sigma_xy, bvfreq, fer_c, fer_K, fer_gamma)
    call fer_gamma2vel(mesh, fer_gamma, dyn%fer_uv)

    call vert_vel_ale(dyn, mesh, dt_velrhs)         ! w + cfl_z + w_e/w_i split + (Fer_GM) fer_w
    Fer_GM = .false.; Redi = .false.   ! restore reduced-M2 for the downstream M2.8/M2.9 step

    ! gate-strength diagnostic (NOT dumped): the smooth Wvel split fires (w_i /= 0) only
    ! where CFL_z > wsplit_mincfl. A non-zero count confirms the split function is genuinely
    ! exercised (not just the trivial w_e=w, w_i=0 range), the L11/L17 guard.
    n_split_faces = 0
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)
            if (dyn%w_i(nz,n) /= 0.0_WP) n_split_faces = n_split_faces + 1
        end do
    end do
    write(*,'(a,es10.3,a,es10.3,a,es10.3,a,i0)') &
        'fesom_pressuredump: ale: max|uv_upd|=', maxval(abs(uv_upd)), &
        ' ; max|hbar|=', maxval(abs(mesh%hbar)), ' ; max|w|=', maxval(abs(dyn%w)), &
        ' ; w_i/=0 (split faces) on ', n_split_faces

    ! ============== M2.8 PP (Pacanowski-Philander) vertical mixing ====================
    ! The Richardson-number mixing coefficients (FESOM2 oce_ale.F90:3728 -> oce_mixing_PP),
    ! the FIRST operator of the dynamics part of the step (it runs BEFORE compute_vel_rhs;
    ! here it is appended as a STANDALONE operator on prescribed inputs — wiring it into the
    ! step order, and sourcing impl_vert_visc_ale's Av from dyn%work%Av, is M2.9). PP reads
    ! the SMOOTHED N^2 (the M2.1 bvfreq, copied into dyn%work%bvfreq) and the nodal velocity
    ! shear (dyn%uvnode, prescribed here with a STRONG vertical shear so the Ri factor spans
    ! [0, ~0.7] and f^2/f^3 are exercised — compute_vel_nodes, the real uvnode source, is an
    ! area-weighted elem->node average gated with the step at M2.9). Force the pi/reduced-M2
    ! mixing knobs (mix_coeff_PP=0.01, A_ver=1e-4 [pi NAMELIST, NOT the 1e-3 module default],
    ! K_ver=1e-5, Kv0_const=.true.). MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    mix_coeff_PP = 0.01_WP
    A_ver        = 1.0e-4_WP
    K_ver        = 1.0e-5_WP
    Kv0_const    = .true.

    ! nodal velocity with strong vertical shear (sin/cos(k*nz)) modulated horizontally so
    ! the Ri factor varies from ~0 (weak shear) to ~0.7 (strong shear). MUST match the oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            dyn%uvnode(1,nz,n) =  1.20_WP*cos(lat)*sin(lon)*sin(0.5_WP*real(nz,WP)) &
                                + 0.20_WP*sin(2.0_WP*lon)*cos(lat)
            dyn%uvnode(2,nz,n) = -0.90_WP*sin(lat)*cos(2.0_WP*lon)*cos(0.4_WP*real(nz,WP)) &
                                + 0.15_WP*cos(lon)
        end do
    end do

    dyn%work%bvfreq = bvfreq                 ! the M2.1 SMOOTHED N^2 (PP reads dyn%work%bvfreq)
    dyn%work%Kv = 0.0_WP                      ! caller pre-zeros (PP writes only nzmin+1..nzmax-1)
    dyn%work%Av = 0.0_WP
    dyn%work%shear2 = 0.0_WP
    call compute_shear2(dyn, mesh)
    call oce_mixing_pp(dyn, mesh)

    ! gate-strength diagnostic (NOT dumped): re-derive the Ri factor (the pass-1 Kv) to
    ! confirm it spans a meaningful range, so the f^2 (Av) and f^3 (Kv) maps are exercised
    ! across their dynamic range, not just at the near-zero background (the L11/L17 guard).
    fmin = 1.0_WP; fmax = 0.0_WP; nf = 0; nf_big = 0
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n)+1, mesh%nlevels_nod2D(n)-1
            dz_inv = 1.0_WP/(mesh%Z_3d_n(nz-1,n) - mesh%Z_3d_n(nz,n))
            shear  = ((dyn%uvnode(1,nz-1,n)-dyn%uvnode(1,nz,n))**2 &
                    + (dyn%uvnode(2,nz-1,n)-dyn%uvnode(2,nz,n))**2)*dz_inv*dz_inv
            factor = shear/(shear + 5.0_WP*max(bvfreq(nz,n),0.0_WP) + 1.0e-14_WP)
            fmin = min(fmin,factor); fmax = max(fmax,factor); nf = nf + 1
            if (factor > 0.1_WP) nf_big = nf_big + 1
        end do
    end do
    write(*,'(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,f5.1,a)') &
        'fesom_pressuredump: PP mix: factor=[', fmin, ',', fmax, &
        '] ; max|Kv|=', maxval(dyn%work%Kv), ' ; max|Av|=', maxval(dyn%work%Av), &
        ' ; factor>0.1 on ', 100.0_WP*real(nf_big,WP)/real(max(nf,1),WP), '% of node-levels'

    ! ============== M5a-2: KPP ri_iwmix interior mixing ==============================
    ! KPP's interior coefficients (Richardson-number shear instability + constant internal-
    ! wave background, oce_ale_mixing_kpp.F90:1008 ri_iwmix), run on the SAME prescribed
    ! strong-shear uvnode + smoothed bvfreq as PP (isolated producer gate). Force the work_core
    ! KPP shear-mixing knobs (visc_sh_limit/diff_sh_limit=5e-3; A_ver/K_ver/Kv0_const already
    ! set above). viscA_ri = node viscosity ; diffK_ri(:,:,1)=T, (:,:,2)=S diffusivity. MUST
    ! match the FESOM2 oracle src/fesom_pressure_dump.F90.
    visc_sh_limit = 5.0e-3_WP
    diff_sh_limit = 5.0e-3_WP
    allocate(viscA_ri(nl, mesh%nod2D), diffK_ri(nl, mesh%nod2D, 2))
    viscA_ri = 0.0_WP; diffK_ri = 0.0_WP
    call ri_iwmix(viscA_ri, diffK_ri, dyn, mesh)

    ! ============== M5a-3: KPP boundary-layer depth (bldepth) ========================
    ! The KPP driver pre-step (dVsq/ustar/Bo, oce_ale_mixing_kpp.F90:344-411) + bldepth
    ! (:746, the C-port HIGHEST-RISK routine: bulk-Ri accumulation + sw interpolation +
    ! ekman/monob limit). Prescribe the surface fluxes (heat/water/wind-stress on nodes) +
    ! a physical decaying shortwave profile sw_3d ANALYTICALLY (byte-identical formula both
    ! sides, M3b style) so Bo spans BOTH signs -> bldepth runs the stable AND unstable
    ! branches. dbsfc is the pressure_bv output (dbsfc_kpp, above). build the wmt/wst lookup
    ! tables first (oce_mixing_kpp_init, M5a-1) so bldepth's wscale reads them. MUST match the
    ! FESOM2 oracle src/fesom_pressure_dump.F90.
    use_sw_pene = .true.
    call oce_mixing_kpp_init(Ricr, concv)   ! Ricr=0.3, concv=1.6 (work_core; mod_param_phys defaults)
    allocate(dVsq_kpp(nl, mesh%nod2D), sw3d_kpp(nl, mesh%nod2D))
    allocate(ustar_kpp(mesh%nod2D), bo_kpp(mesh%nod2D), hbl_kpp(mesh%nod2D))
    allocate(bfsfc_kpp(mesh%nod2D), stable_kpp(mesh%nod2D), caseA_kpp(mesh%nod2D))
    allocate(kbl_kpp(mesh%nod2D), hf_kpp(mesh%nod2D), wf_kpp(mesh%nod2D))
    allocate(stress_node_surf(2, mesh%nod2D))
    dVsq_kpp = 0.0_WP; sw3d_kpp = 0.0_WP

    ! prescribe surface fluxes (heat_flux/water_flux +up; wind stress on nodes) + the
    ! in-water shortwave sw_3d (positive, exp-decaying with depth; surface ~2.4e-5 K m/s =
    ! (1-albw)*0.54*~200 W/m2 / vcpw). heat_flux sign-varying so Bo (= -g*(alpha*hf/vcpw +
    ! beta*S*wf)) spans both signs. MUST match the FESOM2 oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        hf_kpp(n)             =  200.0_WP*cos(lat)*sin(lon)
        wf_kpp(n)             =  2.0e-6_WP*sin(2.0_WP*lon)*cos(lat)
        stress_node_surf(1,n) =  0.15_WP*cos(lat)*sin(lon)
        stress_node_surf(2,n) = -0.10_WP*sin(lat)*cos(2.0_WP*lon)
        sw0 = 2.4e-5_WP*(0.5_WP + 0.5_WP*cos(lon)*cos(lat))             ! >= 0
        do nz = 1, nl
            sw3d_kpp(nz,n) = sw0*exp(real(mesh%zbar_3d_n(nz,n),WP)/20.0_WP)
        end do
    end do

    ! prestep: dVsq (velocity shear re surface @ Z) + ustar (friction velocity) + Bo
    ! (surface turbulent buoyancy forcing). dyn%uvnode is the strong-shear nodal velocity
    ! prescribed for PP/ri_iwmix above; sw_alpha/sw_beta are the M4a producers; salt the
    ! prescribed surface salinity. MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    do n = 1, mesh%nod2D
        nzmin = mesh%ulevels_nod2D(n)
        nzmax = mesh%nlevels_nod2D(n)
        dVsq_kpp(nzmin,n) = 0.0_WP
        usurf = dyn%uvnode(1,nzmin,n)
        vsurf = dyn%uvnode(2,nzmin,n)
        do nz = nzmin+1, nzmax-1
            u_loc = 0.5_WP*( dyn%uvnode(1,nz-1,n) + dyn%uvnode(1,nz,n) )
            v_loc = 0.5_WP*( dyn%uvnode(2,nz-1,n) + dyn%uvnode(2,nz,n) )
            dVsq_kpp(nz,n) = ( usurf - u_loc )**2 + ( vsurf - v_loc )**2
        end do
        dVsq_kpp(nzmax,n) = dVsq_kpp(nzmax-1,n)
        ustar_kpp(n) = sqrt( sqrt( stress_node_surf(1,n)**2 + stress_node_surf(2,n)**2 )*density_0_r )
        bo_kpp(n)    = -g*( sw_alpha(nzmin,n)*hf_kpp(n)/vcpw &
                          + sw_beta (nzmin,n)*wf_kpp(n)*salt(nzmin,n) )
    end do

    call bldepth(dVsq_kpp, dbsfc_kpp, ustar_kpp, bo_kpp, sw3d_kpp, sw_alpha, bvfreq, &
                 hbl_kpp, kbl_kpp, bfsfc_kpp, stable_kpp, caseA_kpp, mesh)

    ! gate-strength diagnostic (NOT dumped): hbl range + how many nodes hit each forcing
    ! branch (stable=1 vs unstable=0) so both wscale/Vtsq paths are genuinely exercised.
    nwp = count(stable_kpp >  0.5_WP)
    nwm = count(stable_kpp <= 0.5_WP)
    write(*,'(a,es10.3,a,es10.3,a,i0,a,i0)') &
        'fesom_pressuredump: bldepth: hbl=[', minval(hbl_kpp), ',', maxval(hbl_kpp), &
        '] m ; stable/unstable nodes = ', nwp, '/', nwm

    ! ============== M5a-4: KPP blmix + enhance + combine + viscAE average ============
    ! The boundary-layer mixing coefficients (oce_ale_mixing_kpp.F90 blmix_kpp:1228 +
    ! enhance:1419) + the driver tail (combine + node->elem viscAE average, :451-492).
    ! blmix reads the ri_iwmix interior viscA_ri/diffK_ri (M5a-2) + the bldepth prestep/
    ! outputs (M5a-3) — all byte-proven — and fills blmc(3)/dkm1(3)/ghats; enhance blends
    ! the kbl-1 interface. The combine takes max(interior,blmc) WITHIN the boundary layer
    ! (into the SEPARATE copies viscA_fin/diffK_fin so the M5a-2 ri_* dumps still echo the
    ! interior) and zeroes ghats outside it; viscAE = node->elem average with the minmix=3e-3
    ! surface floor. MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    allocate(blmc_kpp(nl, mesh%nod2D, 3), ghats_kpp(nl-1, mesh%nod2D), dkm1_kpp(mesh%nod2D, 3))
    allocate(viscA_fin(nl, mesh%nod2D), diffK_fin(nl, mesh%nod2D, 2))
    allocate(viscAE_kpp(nl, mesh%elem2D))
    blmc_kpp = 0.0_WP; ghats_kpp = 0.0_WP; dkm1_kpp = 0.0_WP   ! mirror oce_mixing_kpp_init zero
    viscAE_kpp = 0.0_WP

    call blmix_kpp(viscA_ri, diffK_ri, hbl_kpp, ustar_kpp, bfsfc_kpp, stable_kpp, &
                   caseA_kpp, kbl_kpp, blmc_kpp, ghats_kpp, dkm1_kpp, mesh)
    call enhance(viscA_ri, diffK_ri, hbl_kpp, caseA_kpp, kbl_kpp, &
                 blmc_kpp, ghats_kpp, dkm1_kpp, mesh)

    ! combine: within the boundary layer take max(interior, blmc); outside zero ghats.
    ! Operate on copies so ri_viscA/ri_diffKt/ri_diffKs (M5a-2) keep echoing the interior.
    viscA_fin = viscA_ri
    diffK_fin = diffK_ri
    do n = 1, mesh%nod2D
        nzmin = mesh%ulevels_nod2D(n)
        nzmax = mesh%nlevels_nod2D(n)
        do nz = nzmin+1, nzmax-1
            if (nz < kbl_kpp(n)) then
                viscA_fin(nz,n)   = max(viscA_fin(nz,n),   blmc_kpp(nz,n,1))
                diffK_fin(nz,n,1) = max(diffK_fin(nz,n,1), blmc_kpp(nz,n,2))
                diffK_fin(nz,n,2) = max(diffK_fin(nz,n,2), blmc_kpp(nz,n,3))
            else
                ghats_kpp(nz,n) = 0.0_WP
            end if
        end do
    end do

    ! node->elem viscAE average (minmix=3e-3 surface floor). 1-rank: no pre-average exchange.
    minmix = 3.0e-3_WP
    do elem = 1, mesh%elem2D
        elnodes = mesh%elem2D_nodes(:,elem)
        nzmin   = mesh%ulevels(elem)
        nzmax   = mesh%nlevels(elem)
        do nz = nzmin, nzmax-1
            viscAE_kpp(nz,elem) = sum(mesh%elem_w(:,elem)*viscA_fin(nz,elnodes))
        end do
        viscAE_kpp(nzmax,elem) = viscAE_kpp(nzmax-1,elem)
        if (viscAE_kpp(nzmin,elem) < minmix) viscAE_kpp(nzmin,elem) = minmix
    end do

    ! ============== M2.8b mo_convect convective adjustment ===========================
    ! Static-instability adjustment, run AFTER PP in the step (FESOM2 oce_ale.F90:3729 ->
    ! mo_convect): where N^2<0 floor Kv (nodes) / Av (elements) to instabmix_kv (=0.1).
    ! Save the post-PP Kv/Av first (mo_convect overwrites them in place — the L20 pattern;
    ! the pp_Kv/pp_Av records keep echoing the PP output). Force the gate knobs:
    ! use_instabmix=.true./instabmix_kv=0.1 (pi defaults), use_momix=.false. (CRITICAL —
    ! the TB04 path needs forcing/ice not present pre-forcing; deferred to M2.10),
    ! use_windmix=.false. (pi default; the branch is guarded off). The unstable T band
    ! prescribed above makes bvfreq<0 in the warm hemisphere's upper levels -> non-vacuous.
    ! MUST match the FESOM2 oracle src/fesom_pressure_dump.F90.
    allocate(pp_Kv_dump(nl, mesh%nod2D), pp_Av_dump(nl, mesh%elem2D))
    pp_Kv_dump = dyn%work%Kv                 ! save the PP output (pre-adjustment)
    pp_Av_dump = dyn%work%Av
    use_instabmix = .true.
    instabmix_kv  = 0.1_WP
    use_momix     = .false.                  ! TB04 Monin-Obukhov DEFERRED to M2.10 (forcing/ice)
    use_windmix   = .false.
    call mo_convect(dyn, mesh)

    ! gate-strength diagnostic (NOT dumped): how many node/elem levels were statically
    ! unstable (bvfreq<0) so the convective floor genuinely fired (the L11 weak-gate guard).
    n_unstab_n = 0; n_unstab_e = 0
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n)+1, mesh%nlevels_nod2D(n)-1
            if (bvfreq(nz,n) < 0.0_WP) n_unstab_n = n_unstab_n + 1
        end do
    end do
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e)+1, mesh%nlevels(e)-1
            if (any(bvfreq(nz, mesh%elem2D_nodes(:,e)) < 0.0_WP)) n_unstab_e = n_unstab_e + 1
        end do
    end do
    write(*,'(a,i0,a,i0,a,es10.3,a,es10.3)') &
        'fesom_pressuredump: mo_convect: bvfreq<0 on ', n_unstab_n, ' node-levels / ', &
        n_unstab_e, ' elem-levels ; max|d Kv|=', maxval(abs(dyn%work%Kv - pp_Kv_dump)), &
        ' ; max|d Av|=', maxval(abs(dyn%work%Av - pp_Av_dump))

    ! ============== M2.9a tracer-solve assembly (diff_tracers_ale) ====================
    ! The per-tracer diffusion + ALE reconstruct (FESOM2 oce_ale_tracer.F90:335-491), run
    ! per tracer inside solve_tracers_ale AFTER advection. Here it is ISOLATED on prescribed
    ! inputs (the M2.8 pattern): del_ttf enters with a prescribed advection tendency, the
    ! horizontal diffusivity Ki and the surface fluxes are prescribed (Ki needs mesh_resolution
    ! [M4]; heat_flux/water_flux/virtual_salt/relax_salt come from forcing [M2.10]). The
    ! implicit vertical-diffusion TDMA CONSUMES dyn%work%Kv — the LIVE post-mo_convect PP
    ! coefficient (the FIRST consumer of a PP output; no longer prescribed). MUST match the
    ! FESOM2 oracle src/fesom_pressure_dump.F90.
    !
    ! Gate config: Redi=.false. (isredi=0), tra_adv_lim='FCT' (do_wimpl=.false.),
    ! i_vert_diff=.true. (implicit path), PP (no KPP nonlocal), linfs (is_nonlinfs=0).
    is_nonlinfs = 0.0_WP
    allocate(Ki(nl-1, mesh%nod2D), del_ttf_in(nl-1, mesh%nod2D), tr_xy(2, nl-1, mesh%elem2D))
    allocate(heat_flux(mesh%nod2D), water_flux(mesh%nod2D), virtual_salt(mesh%nod2D), &
             relax_salt(mesh%nod2D), real_salt_flux(mesh%nod2D))
    allocate(t_solved(nl-1, mesh%nod2D), s_solved(nl-1, mesh%nod2D))
    allocate(del_ttf_t(nl-1, mesh%nod2D), del_ttf_s(nl-1, mesh%nod2D))

    ! prescribed horizontal diffusivity Ki (>0 so diff_part_hor_redi is non-vacuous; pi
    ! production K_hor=0 makes it a no-op, so exercise it synthetically — the L21 pattern).
    ! Surface fluxes: heat_flux (T BC), virtual_salt + relax_salt (S BC; the active surface
    ! salinity restoring). water_flux/real_salt_flux drop on linfs (x is_nonlinfs=0) but are
    ! referenced, so prescribe them too. MUST match the FESOM2 oracle.
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            Ki(nz,n)         = 300.0_WP + 150.0_WP*cos(lat)*cos(lon) + 10.0_WP*real(nz,WP)
            del_ttf_in(nz,n) = 0.50_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz,WP))
        end do
        heat_flux(n)      =  50.0_WP*cos(lat)*sin(lon)
        water_flux(n)     =  1.0e-6_WP*sin(2.0_WP*lon)*cos(lat)
        virtual_salt(n)   =  3.0e-5_WP*sin(lon)*cos(lat)
        relax_salt(n)     =  2.0e-5_WP*cos(2.0_WP*lon)*cos(lat)
        real_salt_flux(n) =  0.0_WP
    end do

    ! build the tracer state: T (ID=1) / S (ID=2), values = the prescribed (M2.1) T/S (a COPY,
    ! so the 'temp'/'salt' dump records keep echoing the prescription). FCT + implicit vertical
    ! diffusion (pi namelist.tra). use_wsplit stays .true. (M2.7) but FCT -> do_wimpl=.false.
    tracers%num_tracers = 2
    allocate(tracers%data(2))
    allocate(tracers%data(1)%values(nl-1, mesh%nod2D), tracers%data(2)%values(nl-1, mesh%nod2D))
    allocate(tracers%work%del_ttf(nl-1, mesh%nod2D))
    tracers%data(1)%values = temp;  tracers%data(1)%ID = 1
    tracers%data(2)%values = salt;  tracers%data(2)%ID = 2
    tracers%data(1)%tra_adv_lim = 'FCT';  tracers%data(1)%i_vert_diff = .true.
    tracers%data(2)%tra_adv_lim = 'FCT';  tracers%data(2)%i_vert_diff = .true.

    ! BASELINE FIX (M6a-1): mirror the oracle (fesom_pressure_dump.F90:653) — restore
    ! use_sw_pene=.false. BEFORE the tracer solve. The KPP bldepth above needed it .true.
    ! (line 647), but the M5c diff_ver_part_impl_ale sw_3d term (use_sw_pene .and. id==1)
    ! reads dyn%work%sw_3d, which this shim never allocates -> it must be OFF here (the
    ! oracle's tracer solve also runs use_sw_pene=.false.). Without this, the Release build
    ! read an unallocated allocatable (benign-by-luck at m5 until the heap layout shifted;
    ! debug -check all flags it deterministically). See the M5a-3 META-LESSON (uninit memory).
    use_sw_pene = .false.

    do tr_num = 1, 2
        tracers%work%del_ttf = del_ttf_in                                  ! prescribed advection tendency
        call tracer_gradient_elements(tracers%data(tr_num)%values, tr_xy, mesh)  ! M1.1 tr_xy (horiz diff)
        call diff_tracers_ale(tr_num, dt_velrhs, dyn, tracers, mesh, tr_xy, Ki, &
                              heat_flux, water_flux, virtual_salt, relax_salt, &
                              real_salt_flux, is_nonlinfs)
        if (tr_num == 1) then
            del_ttf_t = tracers%work%del_ttf;  t_solved = tracers%data(1)%values
        else
            del_ttf_s = tracers%work%del_ttf;  s_solved = tracers%data(2)%values
        end if
    end do

    ! gate-strength diagnostic (NOT dumped): the consumed Kv range + how much the solve moved
    ! T/S, confirming the TDMA + bc_surface + horizontal diffusion are non-vacuous.
    kvmin = minval(dyn%work%Kv); kvmax = maxval(dyn%work%Kv)
    write(*,'(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') &
        'fesom_pressuredump: tracer solve: Kv consumed=[', kvmin, ',', kvmax, &
        '] ; max|del_ttf_T|=', maxval(abs(del_ttf_t)), &
        ' ; max|dT|=', maxval(abs(t_solved-temp)), ' ; max|dS|=', maxval(abs(s_solved-salt))

    ! --- dump (same FADVHDMP format / names as the FESOM2 oracle) ---
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r2(u, 'temp',           real(temp,        MP))
    call wr_r2(u, 'salt',           real(salt,        MP))
    call wr_r2(u, 'density_ref',    real(density_ref, MP))
    call wr_r2(u, 'zbar_3d_n',      mesh%zbar_3d_n(1:nl,   1:mesh%nod2D))
    call wr_r2(u, 'Z_3d_n',         mesh%Z_3d_n(1:nl-1,    1:mesh%nod2D))
    call wr_r2(u, 'hnode',          mesh%hnode(1:nl-1,     1:mesh%nod2D))
    call wr_r2(u, 'density_m_rho0', real(density,            MP))
    call wr_r2(u, 'hpressure',      real(hpressure(1:nl-1,:), MP))
    call wr_r2(u, 'bvfreq_raw',     real(bvfreq_raw,         MP))
    call wr_r2(u, 'bvfreq',         real(bvfreq,             MP))
    call wr_r2(u, 'sw_alpha',       real(sw_alpha,           MP))
    call wr_r2(u, 'sw_beta',        real(sw_beta,            MP))
    call wr_r3(u, 'sigma_xy',       real(sigma_xy,           MP))
    call wr_r3(u, 'neutral_slope',  real(neutral_slope,      MP))
    call wr_r3(u, 'slope_tapered',  real(slope_tapered,      MP))
    call wr_r2(u, 'fer_tapfac',     real(fer_tapfac,         MP))
    call wr_r2(u, 'pgf_x',          real(pgf_x(1:nl-1, :),   MP))
    call wr_r2(u, 'pgf_y',          real(pgf_y(1:nl-1, :),   MP))
    call wr_r2(u, 'pgf_x_shchep',   real(pgf_x_shchep(1:nl-1, :), MP))   ! M6a-1
    call wr_r2(u, 'pgf_y_shchep',   real(pgf_y_shchep(1:nl-1, :), MP))   ! M6a-1
    ! M2.3 vel_rhs + M2.4 momadv: gated coriolis + prescribed inputs (incl. w_e) +
    ! the momadv nodal intermediate (uvnode_rhs) + the full-assembly outputs.
    call wr_r1(u, 'coriolis',       real(mesh%coriolis(1:mesh%elem2D), MP))
    call wr_r1(u, 'eta_n',          real(eta_n_in,         MP))   ! prescribed (pre-M2.7)
    call wr_r3(u, 'uv_in',          real(uv_in,            MP))
    call wr_r3(u, 'uv_rhsAB_prev',  real(uv_rhsAB_prev,    MP))
    call wr_r2(u, 'w_e',            real(w_e_in(1:nl, :),  MP))   ! prescribed (pre-split)
    call wr_r3(u, 'uvnode_rhs',     real(uvnode_rhs_dump,  MP))
    call wr_r3(u, 'uv_rhsAB_cor',   real(uv_rhsAB_cor,     MP))
    call wr_r3(u, 'uv_rhs_eul',     real(uv_rhs_eul,       MP))
    call wr_r3(u, 'uv_rhs_ab2',     real(uv_rhs_ab2,       MP))
    ! M2.4 biharmonic viscosity (opt_visc=7): the first-stage Laplacian intermediate
    ! u_c/v_c + the post-viscosity UV_rhs (gate target).
    call wr_r2(u, 'visc_u_c',       real(u_c_dump,         MP))
    call wr_r2(u, 'visc_v_c',       real(v_c_dump,         MP))
    call wr_r3(u, 'uv_rhs_visc',    real(uv_rhs_visc,      MP))
    ! M2.5 implicit vertical viscosity (TDMA): prescribed inputs (Av/stress_surf/w_i) +
    ! the post-solve UV_rhs (gate target).
    call wr_r2(u, 'Av',             real(Av(1:nl, :),      MP))
    call wr_r2(u, 'stress_surf',    real(stress_surf,      MP))
    call wr_r2(u, 'w_i',            real(w_i_in(1:nl, :),  MP))   ! prescribed (pre-split)
    call wr_r3(u, 'uv_rhs_ivv',     real(uv_rhs_ivv,       MP))
    ! M2.6 SSH: the stiffness diagonal + matvec A*eta_n (localise the matrix assembly),
    ! the assembled ssh_rhs, and the CG solution d_eta (the gate target).
    call wr_r1(u, 'ssh_stiff_diag', real(ssh_diag,    MP))
    call wr_r1(u, 'ssh_Aeta',       real(ssh_Aeta,    MP))
    call wr_r1(u, 'ssh_rhs',        real(dyn%ssh_rhs, MP))
    call wr_r1(u, 'd_eta',          real(dyn%d_eta,   MP))
    ! M2.7 ALE velocity/SSH/thickness-W update: the prescribed hbar input + the updated
    ! UV (update_vel) + the divergence ssh_rhs_old + new hbar/dhe (compute_hbar_ale) + the
    ! blended eta_n + the vertical velocity w (vert_vel_ale) + hnode_new (=hnode, linfs) +
    ! cfl_z (compute_CFLz). The explicit/implicit split w_e/w_i (compute_Wvel_split) is
    ! the smooth FESOM3 function and is NOT dumped: nothing reads it, and FESOM2's
    ! w_split_e/w_split_i records are skipped by pressure_diff.py (never byte-compared).
    call wr_r1(u, 'hbar_in',        real(hbar_in,          MP))
    call wr_r3(u, 'uv_upd',         real(uv_upd,           MP))
    call wr_r1(u, 'ssh_rhs_old',    real(dyn%ssh_rhs_old,  MP))
    call wr_r1(u, 'hbar',           real(mesh%hbar,        MP))
    call wr_r1(u, 'dhe',            real(mesh%dhe,         MP))
    call wr_r1(u, 'eta_n_upd',      real(dyn%eta_n,        MP))
    call wr_r2(u, 'w',              real(dyn%w(1:nl, :),   MP))
    call wr_r2(u, 'hnode_new',      real(mesh%hnode_new(1:nl-1, :), MP))
    call wr_r2(u, 'cfl_z',          real(dyn%cfl_z(1:nl, :), MP))
    ! M4b GM: the diffusivity fer_K (nl) + gravity-wave fer_c (cm^2) + scaling fer_scal +
    ! the streamfunction fer_gamma (fer_solve_Gamma TDMA) + bolus velocity fer_uv
    ! (fer_gamma2vel) + bolus vertical velocity fer_w (vert_vel_ale, Fer_GM branch).
    call wr_r2(u, 'fer_K',          real(fer_K(1:nl, :),    MP))
    call wr_r1(u, 'fer_c',          real(fer_c,             MP))
    call wr_r1(u, 'fer_scal',       real(fer_scal,          MP))
    call wr_r3(u, 'fer_gamma',      real(fer_gamma(1:2, 1:nl, :), MP))
    call wr_r3(u, 'fer_uv',         real(dyn%fer_uv,        MP))
    call wr_r2(u, 'fer_w',          real(dyn%fer_w(1:nl, :), MP))
    ! M2.8 PP vertical mixing: the prescribed nodal-velocity shear input + the gate
    ! targets pp_Kv (vertical diffusivity, nodes) and pp_Av (vertical viscosity, elements).
    call wr_r3(u, 'uvnode',         real(dyn%uvnode,            MP))
    call wr_r2(u, 'pp_Kv',          real(pp_Kv_dump(1:nl, :),   MP))   ! PP output (pre-adjustment)
    call wr_r2(u, 'pp_Av',          real(pp_Av_dump(1:nl, :),   MP))
    ! M2.8b mo_convect: the post-adjustment Kv/Av (floored to instabmix_kv where bvfreq<0).
    call wr_r2(u, 'moc_Kv',         real(dyn%work%Kv(1:nl, :),  MP))
    call wr_r2(u, 'moc_Av',         real(dyn%work%Av(1:nl, :),  MP))
    ! M5a-2 KPP ri_iwmix: node viscosity + per-tracer (T/S) interior diffusivity.
    call wr_r2(u, 'ri_viscA',       real(viscA_ri(1:nl, :),     MP))
    call wr_r2(u, 'ri_diffKt',      real(diffK_ri(1:nl, :, 1),  MP))
    call wr_r2(u, 'ri_diffKs',      real(diffK_ri(1:nl, :, 2),  MP))
    ! M5a-3 KPP bldepth: prescribed inputs (wind stress + sw_3d), the prestep
    ! (dVsq/ustar/Bo) + dbsfc (pressure_bv), and the boundary-layer-depth outputs
    ! (hbl/kbl/bfsfc/stable/caseA). kbl dumped as real (exactly representable).
    call wr_r2(u, 'kpp_stress',     real(stress_node_surf,      MP))
    call wr_r2(u, 'kpp_sw3d',       real(sw3d_kpp(1:nl, :),     MP))
    call wr_r2(u, 'kpp_dVsq',       real(dVsq_kpp(1:nl, :),     MP))
    call wr_r2(u, 'kpp_dbsfc',      real(dbsfc_kpp(1:nl, :),    MP))
    call wr_r1(u, 'kpp_ustar',      real(ustar_kpp,             MP))
    call wr_r1(u, 'kpp_Bo',         real(bo_kpp,                MP))
    call wr_r1(u, 'kpp_hbl',        real(hbl_kpp,               MP))
    call wr_r1(u, 'kpp_kbl',        real(kbl_kpp,               MP))
    call wr_r1(u, 'kpp_bfsfc',      real(bfsfc_kpp,             MP))
    call wr_r1(u, 'kpp_stable',     real(stable_kpp,            MP))
    call wr_r1(u, 'kpp_caseA',      real(caseA_kpp,             MP))
    ! M5a-4 KPP blmix+enhance (BL mixing coeffs blmc(1:3)=mom/T/S + kbl-1 diffs dkm1 +
    ! nonlocal ghats), then the combine (final node viscA/diffKt/diffKs) + node->elem viscAE.
    call wr_r2(u, 'kpp_blmc1',      real(blmc_kpp(1:nl, :, 1),  MP))
    call wr_r2(u, 'kpp_blmc2',      real(blmc_kpp(1:nl, :, 2),  MP))
    call wr_r2(u, 'kpp_blmc3',      real(blmc_kpp(1:nl, :, 3),  MP))
    call wr_r1(u, 'kpp_dkm1m',      real(dkm1_kpp(:, 1),        MP))
    call wr_r1(u, 'kpp_dkm1t',      real(dkm1_kpp(:, 2),        MP))
    call wr_r1(u, 'kpp_dkm1s',      real(dkm1_kpp(:, 3),        MP))
    call wr_r2(u, 'kpp_ghats',      real(ghats_kpp(1:nl-1, :),  MP))
    call wr_r2(u, 'kpp_viscA',      real(viscA_fin(1:nl, :),    MP))
    call wr_r2(u, 'kpp_diffKt',     real(diffK_fin(1:nl, :, 1), MP))
    call wr_r2(u, 'kpp_diffKs',     real(diffK_fin(1:nl, :, 2), MP))
    call wr_r2(u, 'kpp_viscAE',     real(viscAE_kpp(1:nl, :),   MP))
    ! M2.9a tracer solve: prescribed inputs (Ki, surface fluxes, the advection tendency
    ! del_ttf_in) + the gate targets per tracer — del_ttf after diffusion (del_ttf_T/S) and
    ! the solved temperature/salinity (tsol_T/tsol_S). The TDMA consumed the live moc_Kv.
    call wr_r2(u, 'tsol_Ki',           real(Ki(1:nl-1, :),     MP))
    call wr_r1(u, 'tsol_heat_flux',    real(heat_flux,         MP))
    call wr_r1(u, 'tsol_virtual_salt', real(virtual_salt,      MP))
    call wr_r1(u, 'tsol_relax_salt',   real(relax_salt,        MP))
    call wr_r2(u, 'tsol_del_ttf_in',   real(del_ttf_in,        MP))
    call wr_r2(u, 'tsol_del_ttf_T',    real(del_ttf_t,         MP))
    call wr_r2(u, 'tsol_del_ttf_S',    real(del_ttf_s,         MP))
    call wr_r2(u, 'tsol_T',            real(t_solved,          MP))
    call wr_r2(u, 'tsol_S',            real(s_solved,          MP))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_pressuredump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_pressuredump
