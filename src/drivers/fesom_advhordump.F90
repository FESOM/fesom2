program fesom_advhordump
    ! M1.1 horizontal-advection operator byte-gate driver. Loads pi 1-rank with the
    ! SAME rotation as the FESOM2 pi run (alpha/beta/gamma=50/15/-90, cyclic 360),
    ! computes geometry, builds helem (linfs full-cell: zbar(nz)-zbar(nz+1)),
    ! PRESCRIBES an analytic velocity + tracer field (identical formula to the
    ! FESOM2 oracle src/fesom_advhor_dump.F90, computed from the byte-identical
    ! coordinates), runs the full horizontal-advection pipeline, and dumps every
    ! field in the FESOM2 oracle's format. tools/advhor_diff.py then compares for
    ! max|delta|=0.
    !
    ! Pipeline (per FESOM2 do_oce_adv_tra horizontal path):
    !   muscl_adv_init -> edge_up_dn_tri
    !   tracer_gradient_elements -> tr_xy (the kernels take the up/dn reconstruction on
    !   the fly; FESOM2's edge_up_dn_grad record is not written)
    !   adv_tra_hor_upw1/_muscl -> adv_flux_hor (edge) ; scatter -> del_ttf_advhoriz
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 14): the FESOM2-oracle
    ! byte gate is retired (docs/HANDOFF.md "Gate migration"). The driver has no triangle-
    ! specific code of its own (it runs the converted kernels) and stays a diagnostic dump.
    !
    !   FESOM3_MESH_DIR    mesh dir   (default: pi)
    !   FESOM3_ADVHOR_OUT  out path   (default: advhor_f3.bin)
    use mpi
    use mod_precision,    only: WP, MP
    use mod_mesh,         only: t_mesh
    use mod_partit,       only: t_partit
    use mod_partitioning, only: par_init, par_ex
    use mod_mesh_read,    only: read_mesh
    use mod_mesh_areas,   only: compute_geometry
    use mod_tracer,       only: t_tracer_work, t_tracer
    use mod_dyn,          only: t_dyn
    use oce_tracer_grad,  only: tracer_gradient_elements
    use oce_muscl_adv,    only: muscl_adv_init
    use oce_adv_tra_hor,  only: adv_tra_hor_upw1, adv_tra_hor_muscl
    use oce_adv_tra_ver,  only: adv_tra_ver_upw1, adv_tra_ver_qr4c
    use oce_adv_tra_flux, only: oce_tra_adv_flux2dtracer
    use oce_adv_tra_fct,  only: oce_tra_adv_fct
    use oce_ale_tracer,   only: advect_tracer       ! M1.4 assembled step body
    use mod_advhor_dump,  only: advhor_dump_open, advhor_dump_close, wr_r1, wr_r2, wr_r3, wr_i1, wr_i2
    implicit none

    ! Controlled-input parameters (must equal the FESOM2 oracle's)
    real(kind=WP), parameter :: dt      = 1800.0_WP   ! time step [s]
    real(kind=WP), parameter :: num_ord = 0.75_WP     ! MUSCL/QR4C 4th-order fraction (exercises both terms)
    real(kind=WP), parameter :: opth    = 0.0_WP      ! MUSCL 4th-order fraction in the pi FCT config (hor.Ord=0 -> 3rd order)
    real(kind=WP), parameter :: optv    = 1.0_WP      ! QR4C 4th-order fraction in the pi FCT config (vert.Ord=1 -> 4th order)

    character(len=512) :: mesh_dir, out_path
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    type(t_tracer_work) :: twork
    integer :: nsw, e, n, nz, nl, u, nzmin, nzmax
    real(kind=WP), allocatable :: vel(:,:,:), ttf(:,:), tr_xy(:,:,:)
    real(kind=WP), allocatable :: aflux_u(:,:)
    real(kind=WP), allocatable :: dttf_u(:,:), dttf_v(:,:), flux_v(:,:)
    real(kind=WP), allocatable :: wvel(:,:), aflux_vu(:,:), aflux_vq(:,:)
    real(kind=WP), allocatable :: dttf_vu(:,:), dttf_vq(:,:), zero_h(:,:)
    ! --- M1.3 FCT path ---
    real(kind=WP), allocatable :: ttfAB(:,:), aflux_mfct(:,:)
    real(kind=WP), allocatable :: fct_LO(:,:), fct_ttf_max(:,:), fct_ttf_min(:,:)
    real(kind=WP), allocatable :: fct_plus(:,:), fct_minus(:,:)
    real(kind=WP), allocatable :: adf_h(:,:), adf_v(:,:)
    real(kind=WP), allocatable :: aflux_h_ho(:,:), aflux_v_ho(:,:)
    real(kind=WP), allocatable :: dttf_h_fct(:,:), dttf_v_fct(:,:)
    integer :: e2, enodes(2), el(2), nl1, nl2, nu1, nu2, nl12, nu12
    real(kind=WP) :: lon, lat
    real(kind=MP) :: zbar_srf, zbar_bot
    ! --- M1.4 assembled step (real init_tracers_AB + do_oce_adv_tra via advect_tracer) ---
    type(t_tracer) :: tr
    type(t_dyn)    :: dyn
    real(kind=WP), allocatable :: valuesAB_s(:,:)
    real(kind=WP), allocatable :: dttf_h_s(:,:),   dttf_v_s(:,:),   dttf_s(:,:)     ! FCT config
    real(kind=WP), allocatable :: dttf_h_sn(:,:),  dttf_v_sn(:,:),  dttf_sn(:,:)    ! non-FCT config

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call get_environment_variable('FESOM3_ADVHOR_OUT', out_path)
    if (len_trim(out_path) == 0) out_path = 'advhor_f3.bin'

    call par_init(partit)
    if (partit%npes /= 1) then
        if (partit%mype == 0) write(*,'(a)') 'fesom_advhordump: requires 1 rank'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype); error stop 1
    end if

    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    nl = mesh%nl
    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'fesom_advhordump: nod2D=', mesh%nod2D, &
        ' elem2D=', mesh%elem2D, ' edge2D=', mesh%edge2D, ' nl=', nl, ' CW swaps=', nsw

    ! helem (linfs, full cells): layer thickness at element = zbar(nz)-zbar(nz+1)
    ! for nz in [ulevels(elem), nlevels(elem)-1]. (FESOM2 init_thickness_ale linfs,
    ! bottom_elem_thickness = zbar(nle-1)-zbar(nle) for full cells; verified by gate.)
    allocate(mesh%helem(nl-1, mesh%elem2D)); mesh%helem = 0.0_MP
    do e = 1, mesh%elem2D
        do nz = mesh%ulevels(e), mesh%nlevels(e)-1
            mesh%helem(nz, e) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do

    ! hnode / hnode_new (per-node layer thickness): FESOM2 init_thickness_ale linfs
    ! (oce_ale.F90:1014-1033) sets hnode(nz,n)=zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n) for
    ! nz in [ulevels_nod2D(n), nlevels_nod2D(n)-2] and hnode(nlevels_nod2D(n)-1,n)=
    ! bottom_node_thickness; at the initial full-cell state both reduce to
    ! zbar(nz)-zbar(nz+1) (the node analog of helem). hnode_new=hnode (oce_ale.F90:1217).
    ! Used by the FCT low-order vertical update, the b2 limiter, and the use_lo scatter;
    ! dumped + gated below.
    allocate(mesh%hnode(nl-1, mesh%nod2D), mesh%hnode_new(nl-1, mesh%nod2D))
    mesh%hnode = 0.0_MP
    do n = 1, mesh%nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n)-1
            mesh%hnode(nz, n) = mesh%zbar(nz) - mesh%zbar(nz+1)
        end do
    end do
    mesh%hnode_new = mesh%hnode

    ! zbar_3d_n / Z_3d_n: per-node ALE interface/mid depths needed by QR4C, built
    ! exactly as FESOM2 init_ale (oce_ale.F90:531-566) at the initial state. pi has
    ! no cavity (ulevels_nod2D=1) and full cells, so zbar_n_srf=zbar(nzmin) and
    ! zbar_n_bot=zbar(nzmax); the surface/bottom Z_3d_n use the /2 boundary formula.
    ! Verified field-by-field by the gate (zbar_3d_n, Z_3d_n).
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

    ! Prescribe analytic velocity (at elements) and tracer (at nodes). MUST match
    ! FESOM2 src/fesom_advhor_dump.F90 exactly:
    !   ttf(nz,n)   = sin(2*lon_n)*cos(lat_n) + 0.02*nz
    !   vel(1,nz,e) = 0.20*cos(lat_e)*sin(lon_e) + 0.001*nz
    !   vel(2,nz,e) = 0.15*cos(2*lat_e)
    ! lon/lat are the rotated coords; lat_e/lon_e use the element's 1st node.
    allocate(ttf(nl-1, mesh%nod2D), vel(2, nl-1, mesh%elem2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            ttf(nz, n) = sin(2.0_WP*lon)*cos(lat) + 0.02_WP*real(nz, WP)
        end do
    end do
    do e = 1, mesh%elem2D
        lon = mesh%coord_nod2D(1, mesh%elem2D_nodes(1, e))
        lat = mesh%coord_nod2D(2, mesh%elem2D_nodes(1, e))
        do nz = 1, nl-1
            vel(1, nz, e) = 0.20_WP*cos(lat)*sin(lon) + 0.001_WP*real(nz, WP)
            vel(2, nz, e) = 0.15_WP*cos(2.0_WP*lat)
        end do
    end do

    ! ttfAB: the Adams-Bashforth (high-order) tracer for the FCT path. The FCT limiter
    ! clips the antidiffusive flux HO(ttfAB)-LO(ttf); a deliberately SHARP, large-
    ! amplitude oscillatory ttfAB (vs the smooth ttf) makes the high-order flux exceed
    ! the cluster bounds set by ttf, so the limiter actively clips (fct_plus/minus<1)
    ! and the b1/b2/b3 limiting logic is genuinely exercised, not just multiplied by 1.
    ! MUST match the FESOM2 oracle src/fesom_advhor_dump.F90 exactly.
    allocate(ttfAB(nl-1, mesh%nod2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl-1
            ttfAB(nz, n) = 20.0_WP*sin(8.0_WP*lon)*cos(6.0_WP*lat) &
                         + 10.0_WP*sin(5.0_WP*lon)*cos(1.5_WP*real(nz, WP))
        end do
    end do

    ! Prescribe analytic vertical velocity Wvel at nodes (nl interfaces). MUST match
    ! FESOM2 src/fesom_advhor_dump.F90 exactly. Sign varies in space (sin 2*lon) and
    ! depth (cos 0.3*nz) so both upwind branches (W>0 / W<0) and the QR4C 4th-order
    ! interior are exercised.
    allocate(wvel(nl, mesh%nod2D))
    do n = 1, mesh%nod2D
        lon = mesh%coord_nod2D(1, n); lat = mesh%coord_nod2D(2, n)
        do nz = 1, nl
            wvel(nz, n) = 1.0e-4_WP*sin(2.0_WP*lon)*cos(lat)*cos(0.3_WP*real(nz, WP))
        end do
    end do

    ! MUSCL setup: edge_up_dn_tri
    call muscl_adv_init(twork, mesh)

    ! Elemental tracer gradient
    allocate(tr_xy(2, nl-1, mesh%elem2D))
    call tracer_gradient_elements(ttf, tr_xy, mesh)

    allocate(aflux_u(nl-1, mesh%edge2D))
    allocate(dttf_u(nl-1, mesh%nod2D))
    allocate(dttf_v(nl-1, mesh%nod2D), flux_v(nl, mesh%nod2D))
    flux_v = 0.0_WP

    ! --- UPW1 horizontal flux + scatter -> del_ttf_advhoriz ---
    call adv_tra_hor_upw1(vel, ttf, mesh, aflux_u, o_init_zero=.true.)
    dttf_u = 0.0_WP; dttf_v = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_u, dttf_v, aflux_u, flux_v, mesh)

    ! FESOM2's clamped MUSCL (adv_flux_hor_muscl / del_ttf_advhoriz_muscl) is not ported:
    ! FESOM3's single MUSCL kernel is FESOM2's MFCT (adv_flux_hor_mfct below).

    ! --- VERTICAL advection: upw1 + qr4c, scatter -> del_ttf_advvert ---
    ! Scatter with a zero horizontal flux so only the vertical branch contributes;
    ! dttf_v here is the throwaway dttf_h argument (stays 0), the gated output is
    ! dttf_vu/dttf_vq.  oce_tra_adv_flux2dtracer:
    !   del_ttf_advvert(nz,n) += (flux_v(nz,n)-flux_v(nz+1,n))*dt/areasvol(nz,n)
    allocate(aflux_vu(nl, mesh%nod2D), aflux_vq(nl, mesh%nod2D))
    allocate(dttf_vu(nl-1, mesh%nod2D), dttf_vq(nl-1, mesh%nod2D))
    allocate(zero_h(nl-1, mesh%edge2D)); zero_h = 0.0_WP

    call adv_tra_ver_upw1(wvel, ttf, mesh, aflux_vu, o_init_zero=.true.)
    dttf_v = 0.0_WP; dttf_vu = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_v, dttf_vu, zero_h, aflux_vu, mesh)

    call adv_tra_ver_qr4c(wvel, ttf, mesh, num_ord, aflux_vq, o_init_zero=.true.)
    dttf_v = 0.0_WP; dttf_vq = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_v, dttf_vq, zero_h, aflux_vq, mesh)

    ! =====================================================================
    ! M1.3 FCT (Zalesak) path. pi namelist.tra tracers 1/2:
    !   tra_adv_hor='MFCT'(opth=0.0; FESOM3: 'MUSCL'), tra_adv_ver='QR4C'(optv=1.0), tra_adv_lim='FCT'.
    ! Mirror do_oce_adv_tra's FCT branch (oce_adv_tra_driver.F90:111-347):
    !   LO = UPW1(ttf) horiz + vert -> fct_LO ; HO antidiffusive = {MFCT(ttfAB)-LO,
    !   QR4C(ttfAB)-LO} ; oce_tra_adv_fct clips ; oce_tra_adv_flux2dtracer(use_lo).
    ! tr_xy is the gradient of TTF (values), NOT ttfAB — as FESOM2's init_tracers_AB
    ! (oce_tracer_mod.F90:127, the valuesAB variant is commented out); it is reused for
    ! the MUSCL(ttfAB) reconstruction.
    ! =====================================================================
    allocate(aflux_mfct(nl-1, mesh%edge2D))
    allocate(fct_LO(nl-1, mesh%nod2D), fct_ttf_max(nl-1, mesh%nod2D), fct_ttf_min(nl-1, mesh%nod2D))
    allocate(fct_plus(nl-1, mesh%nod2D), fct_minus(nl-1, mesh%nod2D))
    allocate(adf_h(nl-1, mesh%edge2D), adf_v(nl, mesh%nod2D))
    allocate(aflux_h_ho(nl-1, mesh%edge2D), aflux_v_ho(nl, mesh%nod2D))
    allocate(dttf_h_fct(nl-1, mesh%nod2D), dttf_v_fct(nl-1, mesh%nod2D))

    ! standalone MUSCL high-order horizontal flux from ttf (num_ord=0.75), dumped under
    ! FESOM2's name adv_flux_hor_mfct: gates the kernel independently of FCT.
    call adv_tra_hor_muscl(vel, ttf, mesh, num_ord, aflux_mfct, tr_xy, twork%edge_up_dn_tri, o_init_zero=.true.)

    ! (1) low-order upwind horizontal flux from ttf -> adf_h
    call adv_tra_hor_upw1(vel, ttf, mesh, adf_h, o_init_zero=.true.)
    ! (2) fct_LO = scatter of the LO horizontal flux over edges
    fct_LO = 0.0_WP
    do e2 = 1, mesh%edge2D
        enodes = mesh%edges(:, e2)
        el     = mesh%edge_tri(:, e2)
        nl1 = mesh%nlevels(el(1))-1; nu1 = mesh%ulevels(el(1)); nl2 = 0; nu2 = 0
        if (el(2) /= el(1)) then
            nl2 = mesh%nlevels(el(2))-1; nu2 = mesh%ulevels(el(2))
        end if
        nl12 = max(nl1, nl2); nu12 = nu1
        if (nu2 > 0) nu12 = min(nu1, nu2)
        do nz = nu12, nl12
            fct_LO(nz, enodes(1)) = fct_LO(nz, enodes(1)) + adf_h(nz, e2)
            fct_LO(nz, enodes(2)) = fct_LO(nz, enodes(2)) - adf_h(nz, e2)
        end do
    end do
    ! (3) low-order upwind vertical flux from ttf -> adf_v
    call adv_tra_ver_upw1(wvel, ttf, mesh, adf_v, o_init_zero=.true.)
    ! (4) finish the low-order solution (combine horiz accumulation + vert flux)
    do n = 1, mesh%nod2D
        nu1 = mesh%ulevels_nod2D(n); nl1 = mesh%nlevels_nod2D(n)
        do nz = nu1, nl1-1
            fct_LO(nz,n) = (ttf(nz,n)*mesh%hnode(nz,n) &
                          + (fct_LO(nz,n) + (adf_v(nz,n)-adf_v(nz+1,n)))*dt/mesh%areasvol(n)) &
                          / mesh%hnode_new(nz,n)
        end do
    end do
    ! (exchange_nod(fct_LO) — no-op at 1 rank)
    ! (5) high-order antidiffusive horizontal flux: MUSCL(ttfAB)-LO (opth, o_init_zero=.false.)
    call adv_tra_hor_muscl(vel, ttfAB, mesh, opth, adf_h, tr_xy, twork%edge_up_dn_tri, o_init_zero=.false.)
    ! (6) high-order antidiffusive vertical flux: QR4C(ttfAB)-LO (optv, o_init_zero=.false.)
    call adv_tra_ver_qr4c(wvel, ttfAB, mesh, optv, adf_v, o_init_zero=.false.)
    aflux_h_ho = adf_h; aflux_v_ho = adf_v          ! capture pre-clip antidiffusive fluxes
    ! (7) FCT limiter: clips adf_h/adf_v in place; fills fct_ttf_max/min + fct_plus/minus
    call oce_tra_adv_fct(dt, ttf, fct_LO, adf_h, adf_v, fct_ttf_min, fct_ttf_max, fct_plus, fct_minus, mesh)
    ! (8) scatter clipped fluxes with the low-order reconstruction -> del_ttf_*_fct
    dttf_h_fct = 0.0_WP; dttf_v_fct = 0.0_WP
    call oce_tra_adv_flux2dtracer(dt, dttf_h_fct, dttf_v_fct, adf_h, adf_v, mesh, &
                                  use_lo=.true., ttf=ttf, lo=fct_LO)

    ! =====================================================================
    ! M1.4 assembled step. Builds a real t_tracer + t_dyn and runs the REAL
    ! advect_tracer (init_tracers_AB -> do_oce_adv_tra -> del_ttf accumulation) —
    ! the same body FESOM2 adv_tracers_ale runs per tracer. Prescribes
    !   values    (smooth) = ttf    (the M1.1 field)
    !   valuesold (sharp)  = ttfAB  (the M1.3 field)
    ! so init_tracers_AB's AB(2) interpolation valuesAB = -(0.5+eps)*valuesold +
    ! (1.5+eps)*values (eps=0.1) is exercised + byte-gated. Two configs: FCT
    ! (MUSCL/QR4C/FCT, opth=0/optv=1, the pi production config) and non-FCT
    ! (MUSCL/QR4C/NON, ph=pv=0.75) which exercises do_zero_flux + the order knobs.
    ! =====================================================================
    tr%num_tracers = 1
    allocate(tr%data(1))
    allocate(tr%data(1)%values   (nl-1, mesh%nod2D))
    allocate(tr%data(1)%valuesAB (nl-1, mesh%nod2D))
    allocate(tr%data(1)%valuesold(2, nl-1, mesh%nod2D))
    allocate(tr%work%fct_LO          (nl-1, mesh%nod2D))
    allocate(tr%work%adv_flux_hor    (nl-1, mesh%edge2D))
    allocate(tr%work%adv_flux_ver    (nl,   mesh%nod2D))
    allocate(tr%work%fct_ttf_max     (nl-1, mesh%nod2D))
    allocate(tr%work%fct_ttf_min     (nl-1, mesh%nod2D))
    allocate(tr%work%fct_plus        (nl-1, mesh%nod2D))
    allocate(tr%work%fct_minus       (nl-1, mesh%nod2D))
    allocate(tr%work%del_ttf         (nl-1, mesh%nod2D))
    allocate(tr%work%del_ttf_advhoriz(nl-1, mesh%nod2D))
    allocate(tr%work%del_ttf_advvert (nl-1, mesh%nod2D))
    call muscl_adv_init(tr%work, mesh)     ! edge_up_dn_tri
    allocate(dyn%uv(2, nl-1, mesh%elem2D))
    allocate(dyn%w(nl, mesh%nod2D), dyn%w_e(nl, mesh%nod2D), dyn%w_i(nl, mesh%nod2D))
    dyn%uv = vel
    dyn%w = wvel; dyn%w_e = wvel; dyn%w_i = wvel     ! use_wsplit=.false. -> w==w_e

    allocate(valuesAB_s(nl-1, mesh%nod2D))
    allocate(dttf_h_s(nl-1, mesh%nod2D),  dttf_v_s(nl-1, mesh%nod2D),  dttf_s(nl-1, mesh%nod2D))
    allocate(dttf_h_sn(nl-1, mesh%nod2D), dttf_v_sn(nl-1, mesh%nod2D), dttf_sn(nl-1, mesh%nod2D))

    ! --- FCT config (pi tracer 1: MUSCL/QR4C/FCT, opth=0/optv=1, AB2) ---
    tr%data(1)%values(:,:)      = ttf
    tr%data(1)%valuesold(1,:,:) = ttfAB
    tr%data(1)%AB_order    = 2
    tr%data(1)%tra_adv_hor = 'MUSCL'
    tr%data(1)%tra_adv_ver = 'QR4C'
    tr%data(1)%tra_adv_lim = 'FCT'
    tr%data(1)%tra_adv_ph  = opth
    tr%data(1)%tra_adv_pv  = optv
    call advect_tracer(dt, 1, dyn, tr, mesh)
    valuesAB_s = real(tr%data(1)%valuesAB,      WP)
    dttf_h_s   = real(tr%work%del_ttf_advhoriz, WP)
    dttf_v_s   = real(tr%work%del_ttf_advvert,  WP)
    dttf_s     = real(tr%work%del_ttf,          WP)

    ! --- non-FCT config (MUSCL/QR4C/NON, ph=pv=0.75): do_zero_flux + order knobs ---
    tr%data(1)%values(:,:)      = ttf
    tr%data(1)%valuesold(1,:,:) = ttfAB
    tr%data(1)%tra_adv_hor = 'MUSCL'
    tr%data(1)%tra_adv_ver = 'QR4C'
    tr%data(1)%tra_adv_lim = 'NON'
    tr%data(1)%tra_adv_ph  = 0.75_WP
    tr%data(1)%tra_adv_pv  = 0.75_WP
    call advect_tracer(dt, 1, dyn, tr, mesh)
    dttf_h_sn = real(tr%work%del_ttf_advhoriz, WP)
    dttf_v_sn = real(tr%work%del_ttf_advvert,  WP)
    dttf_sn   = real(tr%work%del_ttf,          WP)

    ! --- dump (same order/names as the FESOM2 oracle) ---
    call advhor_dump_open(u, trim(out_path), mesh%nod2D, mesh%elem2D, mesh%edge2D, nl)
    call wr_r2(u, 'ttf',                    real(ttf, MP))
    call wr_r3(u, 'vel',                    real(vel, MP))
    call wr_r2(u, 'helem',                  mesh%helem(1:nl-1, 1:mesh%elem2D))
    call wr_i2(u, 'edge_up_dn_tri',         twork%edge_up_dn_tri(1:2, 1:mesh%edge2D))
    call wr_r3(u, 'tr_xy',                  real(tr_xy, MP))
    call wr_r2(u, 'adv_flux_hor_upw1',      real(aflux_u, MP))
    call wr_r2(u, 'del_ttf_advhoriz_upw1',  real(dttf_u, MP))
    call wr_r2(u, 'wvel',                   real(wvel, MP))
    call wr_r2(u, 'zbar_3d_n',              mesh%zbar_3d_n(1:nl,   1:mesh%nod2D))
    call wr_r2(u, 'Z_3d_n',                 mesh%Z_3d_n(1:nl-1,    1:mesh%nod2D))
    call wr_r1(u, 'area',                   mesh%area(1:mesh%nod2D))
    call wr_r2(u, 'adv_flux_ver_upw1',      real(aflux_vu, MP))
    call wr_r2(u, 'del_ttf_advvert_upw1',   real(dttf_vu, MP))
    call wr_r2(u, 'adv_flux_ver_qr4c',      real(aflux_vq, MP))
    call wr_r2(u, 'del_ttf_advvert_qr4c',   real(dttf_vq, MP))
    ! --- M1.3 FCT fields ---
    call wr_r2(u, 'adv_flux_hor_mfct',      real(aflux_mfct, MP))
    call wr_r2(u, 'ttfAB',                  real(ttfAB, MP))
    call wr_r2(u, 'hnode',                  mesh%hnode(1:nl-1, 1:mesh%nod2D))
    call wr_r2(u, 'hnode_new',              mesh%hnode_new(1:nl-1, 1:mesh%nod2D))
    call wr_r2(u, 'fct_LO',                 real(fct_LO, MP))
    call wr_r2(u, 'adv_flux_hor_fct_ho',    real(aflux_h_ho, MP))
    call wr_r2(u, 'adv_flux_ver_fct_ho',    real(aflux_v_ho, MP))
    call wr_r2(u, 'fct_ttf_max',            real(fct_ttf_max, MP))
    call wr_r2(u, 'fct_ttf_min',            real(fct_ttf_min, MP))
    call wr_r2(u, 'fct_plus',               real(fct_plus, MP))
    call wr_r2(u, 'fct_minus',              real(fct_minus, MP))
    call wr_r2(u, 'adv_flux_hor_fct',       real(adf_h, MP))
    call wr_r2(u, 'adv_flux_ver_fct',       real(adf_v, MP))
    call wr_r2(u, 'del_ttf_advhoriz_fct',   real(dttf_h_fct, MP))
    call wr_r2(u, 'del_ttf_advvert_fct',    real(dttf_v_fct, MP))
    ! --- M1.4 assembled-step fields (real advect_tracer) ---
    call wr_r2(u, 'valuesAB',                 real(valuesAB_s, MP))
    call wr_r2(u, 'del_ttf_advhoriz_step',    real(dttf_h_s,  MP))
    call wr_r2(u, 'del_ttf_advvert_step',     real(dttf_v_s,  MP))
    call wr_r2(u, 'del_ttf_step',             real(dttf_s,    MP))
    call wr_r2(u, 'del_ttf_advhoriz_stepnon', real(dttf_h_sn, MP))
    call wr_r2(u, 'del_ttf_advvert_stepnon',  real(dttf_v_sn, MP))
    call wr_r2(u, 'del_ttf_step_non',         real(dttf_sn,   MP))
    call advhor_dump_close(u)
    write(*,'(a)') 'fesom_advhordump: wrote '//trim(out_path)

    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program fesom_advhordump
