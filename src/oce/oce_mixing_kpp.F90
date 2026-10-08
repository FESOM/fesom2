module oce_mixing_kpp
    ! FESOM K-Profile Parameterization (Large et al. 1994) vertical mixing — the
    ! production CORE2 scheme (work_core mix_scheme='KPP', mix_scheme_nmb==1; NOT cvmix).
    ! Transcribed VERBATIM from FESOM2 v2.7.3 oce_ale_mixing_kpp.F90 (MODULE o_mixing_KPP_mod).
    !
    ! M5 milestone (KPP). Ported incrementally, byte-gated max|Δ|=0 vs the oracle, mirroring
    ! M2/M3/M4. Optional-`partit` multi-rank pattern from the start (M4f lesson).
    !
    ! M5a-1 SCOPE (this commit): oce_mixing_kpp_init (the rank-independent constants Vtc/cg/
    ! deltaz/deltau + the wmt/wst turbulent-velocity-scale lookup tables, eqn 13 & B1) + wscale
    ! (the lookup into those tables). These are PURE math — no mesh/partit/forcing state — so
    ! they gate standalone against the oracle's FESOM_KPP_DUMP_DIR init/wscale dumps.
    !
    ! STRUCTURAL NOTE (byte-neutral): the oracle's oce_mixing_kpp_init ALSO allocates+zeroes the
    ! per-node state arrays (ghats/hbl/blmc/...). FESOM3 keeps explicit dataflow: those arrays
    ! live in dyn%work (allocated in the KPP setup, M5b); this init computes ONLY the tables/
    ! constants. The arithmetic that fills wmt/wst/Vtc/cg is byte-identical to the oracle.
    !
    ! The lookup tables + constants are computed ONCE (at setup) and are RANK-INDEPENDENT, so
    ! they are module-saved (mirrors the oracle's module variables).
    use mod_precision, only: WP
    implicit none
    private

    public :: oce_mixing_kpp_driver   ! the KPP driver (FESOM2 oce_mixing_KPP; renamed to
                                      ! avoid the case-insensitive clash with the module name)
    public :: oce_mixing_kpp_init
    public :: wscale
    public :: ri_iwmix
    public :: bldepth
    public :: blmix_kpp
    public :: enhance
    ! Exposed (read-only init results) for the M5a-1 gate driver fesom_kppdump.
    public :: Vtc, cg, deltaz, deltau, wmt, wst, nni, nnj

    ! --- KPP fixed parameters (oce_ale_mixing_kpp.F90:50-72) -------------------
    real(kind=WP), parameter :: epsln       = 1.0e-40_WP   ! a small value
    real(kind=WP), parameter :: epsilon_kpp = 0.1_WP
    real(kind=WP), parameter :: vonk        = 0.4_WP       ! von Karman constant
    real(kind=WP), parameter :: conc1       = 5.0_WP

    real(kind=WP), parameter :: zmin = -4.e-7_WP   ! m3/s3 limit for lookup table of wm and ws
    real(kind=WP), parameter :: zmax =  0.0_WP     ! m3/s3 limit for lookup table of wm and ws
    real(kind=WP), parameter :: umin =  0.0_WP     ! m/s limit for lookup table of wm and ws
    real(kind=WP), parameter :: umax =  0.04_WP    ! m/s limit for lookup table of wm and ws

    integer, parameter :: nni = 890   ! number of values for zehat in the lookup table
    integer, parameter :: nnj = 480   ! number of values for ustar in the lookup table

    ! --- module-saved constants/tables (computed once in oce_mixing_kpp_init) ---
    real(kind=WP), save :: cg      ! non-dimensional coefficient for counter-gradient term
    real(kind=WP), save :: Vtc     ! non-dim coeff for turbulent velocity shear (eqn 23)
    real(kind=WP), save :: deltaz  ! delta zehat in table
    real(kind=WP), save :: deltau  ! delta ustar in table
    real(kind=WP), save, dimension(0:nni+1,0:nnj+1) :: wmt ! lookup table for wm (momentum)
    real(kind=WP), save, dimension(0:nni+1,0:nnj+1) :: wst ! lookup table for ws (scalars)

contains

    !#######################################################################
    ! The KPP vertical-mixing driver (= oce_mixing_KPP, oce_ale_mixing_kpp.F90:306).
    ! Assembles the byte-proven sub-kernels (M5a) into the production scheme:
    !   prestep   dVsq (surface-referenced velocity shear @ Z, eqn 21) + ustar/Bo
    !             (surface friction velocity + turbulent buoyancy forcing, eqn 2 &
    !             A2/A3); dbsfc is filled upstream by pressure_bv (optional output, M5a-3).
    !   ri_iwmix  interior viscA/diffK (shear instability + background + static, M5a-2)
    !   bldepth   OBL depth hbl/kbl + bfsfc/stable/caseA (M5a-3)
    !   blmix_kpp BL mixing coeffs blmc(3) + dkm1(3) + nonlocal ghats (M5a-4)
    !   enhance   kbl-1 interface blend (M5a-4)
    !   combine   within the BL take max(interior, blmc); outside zero ghats (in place)
    !   average   node viscA -> element viscAE (Av, kite-weighted elem_w mean), minmix=3e-3
    !             surface floor
    ! Outputs: dynamics%work%Av (element momentum viscosity; impl_vert_visc_ale reads it,
    ! UNCHANGED) + dynamics%work%Kv_double (node T/S diffusivity; step_oce copies (:,:,1)
    ! -> Kv for the tracer TDMA, UNCHANGED). Kv_double(:,:,2) (S channel) is a dead output
    ! in the reduced config (only feeds the gated dump / blmc(:,:,3)).
    !
    ! All scratch (dVsq/ustar/Bo/hbl/kbl/bfsfc/stable/caseA/blmc/ghats/dkm1/viscA_kpp)
    ! lives in dynamics%work (allocated in the KPP setup). Optional-partit from the start:
    ! owned-node loops + the FESOM2 end-of-driver halo exchanges (M5d-ready; at 1-rank the
    ! exchanges are skipped and owned == all). The prestep is the verbatim oce_ale_mixing_kpp
    ! .F90:347-413 dVsq/ustar/Bo and the combine + node->elem average is :452-492 — both
    ! already byte-proven inline in the M5a-3/a-4 CORE2+pi pressure gate.
    subroutine oce_mixing_kpp_driver(dynamics, tracers, stress_node_surf, heat_flux, water_flux, mesh, partit)
        use mod_mesh,        only: t_mesh, MAX_NV
        use mod_dyn,         only: t_dyn
        use mod_tracer,      only: t_tracer
        use mod_partit,      only: t_partit
        use mod_part_bounds, only: owned_bounds, is_multirank
        use mod_halo,        only: exchange_nod
        use mod_constants,   only: g, vcpw, density_0_r
        use oce_pressure_bv, only: smooth_nod
        type(t_dyn),    intent(inout), target :: dynamics
        type(t_tracer), intent(in),    target :: tracers
        real(kind=WP),  intent(in)            :: stress_node_surf(:,:)  ! (2, nnod) surface stress on nodes
        real(kind=WP),  intent(in)            :: heat_flux(:)           ! (nnod)   +up
        real(kind=WP),  intent(in)            :: water_flux(:)          ! (nnod)   +up
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in),    optional :: partit

        real(kind=WP), parameter :: minmix = 3.0e-3_WP  ! surface viscAE floor (avoids huge surf vel)
        integer :: node, nz, nzmin, nzmax, elem, elnodes(MAX_NV)
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: usurf, vsurf, u_loc, v_loc
        real(kind=WP), pointer :: viscAE(:,:), viscA(:,:), Kv_double(:,:,:)
        real(kind=WP), pointer :: sw_alpha(:,:), sw_beta(:,:), dbsfc(:,:), dVsq(:,:)
        real(kind=WP), pointer :: ustar(:), Bo(:), salt(:,:)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        viscAE    => dynamics%work%Av
        viscA     => dynamics%work%viscA_kpp
        Kv_double => dynamics%work%Kv_double
        sw_alpha  => dynamics%work%sw_alpha
        sw_beta   => dynamics%work%sw_beta
        dbsfc     => dynamics%work%dbsfc
        dVsq      => dynamics%work%dVsq
        ustar     => dynamics%work%ustar
        Bo        => dynamics%work%Bo
        salt      => tracers%data(2)%values

        !_______________________________________________________________________
        ! reset the node momentum viscosity over owned+halo (:342-344)
        do node = 1, nNodL
            viscA(:, node) = 0.0_WP
        end do

        !_______________________________________________________________________
        ! prestep: dVsq = squared velocity shear referenced to the surface @ Z (:367-383);
        ! dbsfc surface reference zeroed (its interior is filled by pressure_bv, M5a-3).
        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            dVsq (nzmin,node) = 0.0_WP
            dbsfc(nzmin,node) = 0.0_WP
            usurf = dynamics%uvnode(1,nzmin,node)
            vsurf = dynamics%uvnode(2,nzmin,node)
            do nz = nzmin+1, nzmax-1
                u_loc = 0.5_WP*( dynamics%uvnode(1,nz-1,node) + dynamics%uvnode(1,nz,node) )
                v_loc = 0.5_WP*( dynamics%uvnode(2,nz-1,node) + dynamics%uvnode(2,nz,node) )
                dVsq(nz,node) = ( usurf - u_loc )**2 + ( vsurf - v_loc )**2
            end do
            dVsq(nzmax,node) = dVsq(nzmax-1,node)
        end do

        !_______________________________________________________________________
        ! friction velocity ustar (eqn 2) + surface turbulent buoyancy forcing Bo
        ! (eqns A2c/A2d/A3b/A3d). heat_flux/water_flux positive up. :403-413.
        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            ustar(node) = sqrt( sqrt( stress_node_surf(1,node)**2 + stress_node_surf(2,node)**2 )*density_0_r )
            Bo(node)    = -g*( sw_alpha(nzmin,node)*heat_flux(node)/vcpw &
                             + sw_beta (nzmin,node)*water_flux(node)*salt(nzmin,node) )
        end do

        !_______________________________________________________________________
        ! interior mixing coefficients (Ri shear instability + constant background +
        ! static instability); writes viscA + Kv_double(:,:,1:2). :421.
        call ri_iwmix(viscA, Kv_double, dynamics, mesh, partit)

        !_______________________________________________________________________
        ! boundary-layer depth (hbl/kbl + bfsfc/stable/caseA), boundary-layer
        ! diffusivities (blmc/dkm1/ghats), and the kbl-1 interface enhancement. :430-437.
        call bldepth(dVsq, dbsfc, ustar, Bo, dynamics%work%sw_3d, sw_alpha, dynamics%work%bvfreq, &
                     dynamics%work%hbl, dynamics%work%kbl, dynamics%work%bfsfc, &
                     dynamics%work%stable, dynamics%work%caseA, mesh, partit)
        call blmix_kpp(viscA, Kv_double, dynamics%work%hbl, ustar, dynamics%work%bfsfc, &
                       dynamics%work%stable, dynamics%work%caseA, dynamics%work%kbl, &
                       dynamics%work%blmc, dynamics%work%ghats, dynamics%work%dkm1, mesh, partit)
        call enhance(viscA, Kv_double, dynamics%work%hbl, dynamics%work%caseA, dynamics%work%kbl, &
                     dynamics%work%blmc, dynamics%work%ghats, dynamics%work%dkm1, mesh, partit)

        !_______________________________________________________________________
        ! smooth_blmc=.true. (oracle default, oce_ale_mixing_kpp.F90:439-449): 3 area-weighted
        ! horizontal smoothing sweeps per BL-mixing channel (the lumped-P1 smooth_nod3D = the
        ! M2.1 bvfreq smoother). The pre-smooth exchange makes the blmix-written OWNED blmc
        ! halo-valid before the patch accumulation reads neighbours (multi-rank); at 1-rank
        ! both the exchange and smooth_nod's per-sweep exchange are no-ops but the smoothing
        ! itself runs. NOT exercised by the M5a-4 isolated gate (that called blmix/enhance
        ! directly) -> this is the live-driver step that the lifecycle needs.
        if (is_multirank(partit)) then
            call exchange_nod(dynamics%work%blmc(:,:,1), partit)
            call exchange_nod(dynamics%work%blmc(:,:,2), partit)
            call exchange_nod(dynamics%work%blmc(:,:,3), partit)
        end if
        call smooth_nod(dynamics%work%blmc(:,:,1), 3, mesh, partit)
        call smooth_nod(dynamics%work%blmc(:,:,2), 3, mesh, partit)
        call smooth_nod(dynamics%work%blmc(:,:,3), 3, mesh, partit)

        !_______________________________________________________________________
        ! combine: within the boundary layer take max(interior, blmc); outside zero the
        ! nonlocal ghats. In place on viscA / Kv_double (:452-466).
        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            do nz = nzmin+1, nzmax-1
                if (nz < dynamics%work%kbl(node)) then
                    viscA(nz,node)       = max(viscA(nz,node),       dynamics%work%blmc(nz,node,1))
                    Kv_double(nz,node,1) = max(Kv_double(nz,node,1), dynamics%work%blmc(nz,node,2))
                    Kv_double(nz,node,2) = max(Kv_double(nz,node,2), dynamics%work%blmc(nz,node,3))
                else
                    dynamics%work%ghats(nz,node) = 0.0_WP
                end if
            end do
        end do

        !_______________________________________________________________________
        ! halo exchange before the node->element average (multi-rank only; the owned
        ! loops above leave the halo stale). FESOM2 :470-475. 1-rank: skipped.
        if (is_multirank(partit)) then
            call exchange_nod(Kv_double(:,:,1), partit)
            call exchange_nod(Kv_double(:,:,2), partit)
            call exchange_nod(dynamics%work%ghats, partit)
            call exchange_nod(viscA, partit)
        end if

        !_______________________________________________________________________
        ! node viscA -> element viscAE (Av), with the minmix surface floor. :478-492.
        ! Element mean = kite-weighted sum over the 4 vertex slots (FESOM2: plain /3).
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:,elem)
            nzmin   = mesh%ulevels(elem)
            nzmax   = mesh%nlevels(elem)
            do nz = nzmin, nzmax-1
                viscAE(nz,elem) = sum(mesh%elem_w(:,elem)*viscA(nz,elnodes))
            end do
            viscAE(nzmax,elem) = viscAE(nzmax-1,elem)
            if (viscAE(nzmin,elem) < minmix) viscAE(nzmin,elem) = minmix
        end do
    end subroutine oce_mixing_kpp_driver

    !#######################################################################
    ! Initialization for the KPP vertical mixing scheme.
    ! output: Vtc, cg, wmt, wst (= oce_mixing_kpp_init, oce_ale_mixing_kpp.F90:111).
    ! Ricr/concv are the work_core namelist.oce DOUBLES (Ricr=0.3, concv=1.6).
    subroutine oce_mixing_kpp_init(Ricr, concv)
        real(kind=WP), intent(in) :: Ricr, concv

        real(kind=WP), parameter :: cstar = 10.0_WP    ! proportionality coeff for nonlocal transport
        real(kind=WP), parameter :: conam =  1.257_WP
        real(kind=WP), parameter :: concm =  8.380_WP
        real(kind=WP), parameter :: conc2 = 16.0_WP
        real(kind=WP), parameter :: zetam = -0.2_WP
        real(kind=WP), parameter :: conas = -28.86_WP
        real(kind=WP), parameter :: concs = 98.96_WP
        real(kind=WP), parameter :: conc3 = 16.0_WP
        real(kind=WP), parameter :: zetas = -1.0_WP

        real(kind=WP) :: zehat   ! = zeta * ustar**3
        real(kind=WP) :: zeta    ! = stability parameter d/L
        real(kind=WP) :: usta
        integer :: i, j

        ! Vtc used in eqn. 23
        Vtc = concv * sqrt(0.2_WP/concs/epsilon_kpp) / vonk**2 / Ricr

        ! cg = cs in eqn. 20
        cg = cstar * vonk * (concs * vonk * epsilon_kpp)**(1._WP/3._WP)

        ! Construct the wm and ws lookup tables (eqn. 13 & B1)
        deltaz = (zmax-zmin)/real(nni+1,WP)
        deltau = (umax-umin)/real(nnj+1,WP)

        do i=0,nni+1
            zehat = deltaz*(i) + zmin
            do j=0,nnj+1
                usta = deltau*(j) + umin
                zeta = zehat/(usta**3+epsln)

                if(zehat >= 0._WP) then
                    wmt(i,j) = vonk*usta/(1.+conc1*zeta)
                    wst(i,j) = wmt(i,j)
                else
                    if(zeta > zetam) then
                        wmt(i,j) = vonk* usta * (1._WP-conc2*zeta)**(1._WP/4._WP)
                    else
                        wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1._WP/3._WP)
                    endif
                    if(zeta > zetas) then
                        wst(i,j) = vonk* usta * (1._WP-conc3*zeta)**(1._WP/2._WP)
                    else
                        wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1._WP/3._WP)
                    endif
                endif
            enddo
        enddo
    end subroutine oce_mixing_kpp_init

    !#######################################################################
    ! Compute turbulent velocity scales wm, ws via the 2D lookup table for
    ! unstable conditions (zehat <= zmax); direct formula otherwise.
    ! (= wscale, oce_ale_mixing_kpp.F90:946.)
    subroutine wscale(zehat, us, wm, ws)
        real(kind=WP), intent(in)  :: zehat, us
        real(kind=WP), intent(out) :: wm, ws
        real(kind=WP) :: zdiff, udiff, zfrac, ufrac, fzfrac
        real(kind=WP) :: wam, wbm, was, wbs, u3
        integer :: iz, izp1, ju, jup1

        ! use lookup table for zehat < zmax only; otherwise use stable formulae
        ! zehat = vonk * sigma * abs(depth) * bfsfc
        IF (zehat <= zmax) THEN
            zdiff = zehat-zmin
            iz    = INT( zdiff/deltaz)
            iz    = MIN( iz , nni )
            iz    = MAX( iz , 0  )
            izp1  = iz + 1

            udiff = us-umin
            ju    = INT( MIN(udiff/deltau,real(nnj,WP)))
            ju    = MAX( ju , 0  )
            jup1  = ju+1

            zfrac = zdiff/deltaz - real(iz,WP)
            ufrac = udiff/deltau - real(ju,WP)

            fzfrac= 1._WP-zfrac
            wam   = (fzfrac)  * wmt(iz,jup1) + zfrac * wmt(izp1,jup1)
            wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac * wmt(izp1,ju  )
            wm    = (1._WP-ufrac)* wbm          + ufrac * wam

            was   = (fzfrac)  * wst(iz,jup1) + zfrac * wst(izp1,jup1)
            wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac * wst(izp1,ju  )
            ws    = (1._WP-ufrac)* wbs          + ufrac * was
        ELSE
            u3    = us*us*us
            wm    = vonk * us * u3 / ( u3 + conc1*zehat + epsln )
            ws    = wm
        ENDIF
    end subroutine wscale

    !#######################################################################
    ! Interior viscosity/diffusivity from shear instability (local Richardson
    ! number), constant internal-wave background, and static instability.
    ! (= ri_iwmix, oce_ale_mixing_kpp.F90:1008.) Outputs:
    !   viscA(:,node)      node viscosity (momentum, m^2/s)
    !   diffK(:,node,1)    T diffusivity ; diffK(:,node,2) = S diffusivity (= copy of T here)
    ! Loops OWNED nodes only (the oracle loops myDim_nod2D; the full KPP driver halo-exchanges
    ! viscA/diffK at the end). smooth_Ri_ver/smooth_Ri_hor = .false. -> the smoothing passes are
    ! dead under work_core and omitted (re-add if ever enabled).
    !
    ! L29: writes go through local POINTERS (vA/dK), mirroring oce_mixing_pp — the "compiler
    ! can't disprove aliasing -> scalar divide" idiom that keeps the Ri-number divide byte-exact
    ! vs the oracle (whose diffK is an explicit-shape dummy). AMAX1/AMIN1 kept VERBATIM (REAL
    ! intrinsics; the C port's 1e-9 drift came from mapping them to fmax — Fortran-to-Fortran is
    ! exact because both compile the same AMAX1 with the same ifort).
    subroutine ri_iwmix(viscA, diffK, dyn, mesh, partit)
        use mod_mesh,          only: t_mesh
        use mod_dyn,           only: t_dyn
        use mod_partit,        only: t_partit
        use mod_part_bounds,   only: owned_bounds
        use mod_param_phys,    only: A_ver, K_ver, Kv0_const, visc_sh_limit, diff_sh_limit
        use mod_constants,     only: rad
        use oce_ale_mixing_pp, only: Kv0_background_qiang
        real(kind=WP), dimension(:,:),   intent(inout), target :: viscA  ! (nl, nnod_alloc)
        real(kind=WP), dimension(:,:,:), intent(inout), target :: diffK  ! (nl, nnod_alloc, ntr)
        type(t_dyn),    intent(in), target   :: dyn
        type(t_mesh),   intent(in), target   :: mesh
        type(t_partit), intent(in), optional :: partit

        real(kind=WP), parameter :: Riinfty = 0.8_WP   ! local Ri limit for shear instability
        integer :: node, nz, nzmin, nzmax, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: shear, Rigg, ratio, frit, Kv0_b
        real(kind=WP), pointer :: vA(:,:), bvf(:,:)
        real(kind=WP), pointer :: dK(:,:,:), sh2(:,:)

        vA  => viscA
        dK  => diffK
        bvf => dyn%work%bvfreq
        sh2 => dyn%work%shear2
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !_______________________________________________________________________
        ! Richardson number, stored temporarily in dK(:,:,1) to save memory.
        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            do nz = nzmin+1, nzmax-1
                ! shared shear (oce_vert_spline::compute_shear2) -- this is the direct
                ! analogue of the Rig that ROMS reconstructs under RI_SPLINES.
                shear  = sh2(nz,node)
                dK(nz,node,1) = MAX( bvf(nz,node), 0.0_WP ) / (shear + epsln)  ! avoid NaN at start
            end do
            ! surface/bottom are not used by the model (diffK @ zbar)
            dK(nzmin,node,1) = dK(nzmin+1,node,1)
            dK(nzmax,node,1) = dK(nzmax-1,node,1)
        end do

        !_______________________________________________________________________
        ! viscA + diffK from the Ri factor (eqn 28b&c shear-instability shape).
        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            do nz = nzmin+1, nzmax-1
                Rigg  = AMAX1( dK(nz,node,1), 0.0_WP )
                ratio = AMIN1( Rigg/Riinfty, 1.0_WP )
                frit  = (1.0_WP - ratio*ratio)
                frit  = frit*frit*frit
                vA(nz,node) = visc_sh_limit * frit + A_ver
                if (Kv0_const) then
                    dK(nz,node,1) = diff_sh_limit * frit + K_ver
                else
                    call Kv0_background_qiang(Kv0_b, real(mesh%geo_coord_nod2D(2,node),WP)/rad, &
                                              abs(real(mesh%zbar_3d_n(nz,node),WP)))
                    dK(nz,node,1) = diff_sh_limit * frit + Kv0_b
                end if
                dK(nz,node,2) = dK(nz,node,1)
            end do
            vA(nzmin,node)   = vA(nzmin+1,node)
            dK(nzmin,node,1) = dK(nzmin+1,node,1)
            dK(nzmin,node,2) = dK(nzmin+1,node,2)
            vA(nzmax,node)   = vA(nzmax-1,node)
            dK(nzmax,node,1) = dK(nzmax-1,node,1)
            dK(nzmax,node,2) = dK(nzmax-1,node,2)
        end do
    end subroutine ri_iwmix

    !#######################################################################
    ! Diagnose the ocean-boundary-layer depth hbl + the first level below it
    ! kbl, plus the surface buoyancy forcing bfsfc, the stable flag, and caseA
    ! (= bldepth, oce_ale_mixing_kpp.F90:746). The bulk Richardson number
    !   Rib(z) = z*dbsfc(z) / ( dVsq(z) + Vtsq(z) + eps )
    ! (eqn 21) accumulates downward until it crosses Ricr; hbl is the linear-
    ! interpolated crossing, then limited by the Ekman/Monin-Obukhov depths
    ! (eqn 24) in stable forcing. Vtsq is the unresolved-shear term (eqn 23,
    ! velocity scale ws from wscale). With use_sw_pene the surface buoyancy
    ! forcing bfsfc absorbs the in-water shortwave (sw_3d) penetrating to z.
    !
    ! Inputs (prescribed/produced upstream): dVsq/dbsfc (driver prestep +
    ! pressure_bv), ustar/Bo (prestep), sw_3d (cal_shortwave_rad), sw_alpha
    ! (M4a, surface only), bvfreq (M2.3). Outputs: hbl/kbl/bfsfc/stable/caseA.
    ! Loops OWNED nodes (the oracle's myDim_nod2D; smooth_hbl=.false. so no
    ! exchange — the KPP driver halo-exchanges hbl/blmc later). Optional-partit
    ! from the start (M5d-ready). The C port marks this the HIGHEST-RISK routine
    ! (Rib accumulation + sw interp); every divide here is scalar (the inner
    ! nz-loop EXITs at the crossing + calls wscale, so no L29 SIMD trap).
    ! SIGN/AMIN1 kept VERBATIM (REAL intrinsics; Fortran-to-Fortran exact).
    subroutine bldepth(dVsq, dbsfc, ustar, Bo, sw_3d, sw_alpha, bvfreq, &
                       hbl, kbl, bfsfc, stable, caseA, mesh, partit)
        use mod_mesh,        only: t_mesh
        use mod_partit,      only: t_partit
        use mod_part_bounds, only: owned_bounds
        use mod_param_phys,  only: Ricr
        use mod_config,      only: use_sw_pene
        use mod_constants,   only: g
        real(kind=WP), dimension(:,:), intent(in)    :: dVsq      ! (nl,   nnod)
        real(kind=WP), dimension(:,:), intent(in)    :: dbsfc     ! (nl,   nnod)
        real(kind=WP), dimension(:),   intent(in)    :: ustar     ! (nnod)
        real(kind=WP), dimension(:),   intent(in)    :: Bo        ! (nnod)
        real(kind=WP), dimension(:,:), intent(in)    :: sw_3d     ! (nl,   nnod)
        real(kind=WP), dimension(:,:), intent(in)    :: sw_alpha  ! (nl-1, nnod)
        real(kind=WP), dimension(:,:), intent(in)    :: bvfreq    ! (nl,   nnod)
        real(kind=WP), dimension(:),   intent(inout) :: hbl       ! (nnod)
        integer,       dimension(:),   intent(inout) :: kbl       ! (nnod)
        real(kind=WP), dimension(:),   intent(inout) :: bfsfc     ! (nnod)
        real(kind=WP), dimension(:),   intent(inout) :: stable    ! (nnod)
        real(kind=WP), dimension(:),   intent(inout) :: caseA     ! (nnod)
        type(t_mesh),   intent(in), target   :: mesh
        type(t_partit), intent(in), optional :: partit

        real(kind=WP), parameter :: cekman = 0.7_WP   ! constant for Ekman depth
        real(kind=WP), parameter :: cmonob = 1.0_WP   ! constant for Monin-Obukhov depth
        integer :: node, nz, nzmin, nzmax, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: Ritop, bvsq, Vtsq, hekman, hmonob, hlimit
        real(kind=WP) :: Rib_km1, Rib_k, coeff_sw, zk, zkm1, sigma, zehat, wm, ws, dzup

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! initialize hbl and kbl to bottomed-out values
        do node = 1, nNodO
            kbl(node) = mesh%nlevels_nod2D(node)
            hbl(node) = ABS( mesh%zbar_3d_n( mesh%nlevels_nod2D(node), node ) )
        end do

        do node = 1, nNodO
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)

            if (use_sw_pene) coeff_sw = g * sw_alpha(nzmin,node)   ! @ the surface @ Z (m/s2/K)

            Rib_km1     = 0.0_WP
            bfsfc(node) = Bo(node)

            do nz = nzmin+1, nzmax
                zk   = ABS( mesh%zbar_3d_n(nz,  node) )
                zkm1 = ABS( mesh%zbar_3d_n(nz-1,node) )

                ! bfsfc = Bo + sw contribution (sw_3d: K m/s, positive downward)
                if (use_sw_pene) &
                    bfsfc(node) = Bo(node) + coeff_sw * ( sw_3d(nzmin,node) - sw_3d(nz,node) )

                stable(node) = 0.5_WP + SIGN( 0.5_WP, bfsfc(node) )
                sigma        = stable(node) + ( 1.0_WP - stable(node) ) * epsilon_kpp

                ! velocity scales at sigma, for z=-zbar(nz) (eqn 23)
                zehat = vonk * sigma * zk * bfsfc(node)
                call wscale(zehat, ustar(node), wm, ws)

                bvsq = bvfreq(nz,node)                  ! N^2 @ zbar
                Vtsq = zk * ws * SQRT(ABS(bvsq)) * Vtc

                ! bulk Richardson number at the new level (eqn 21)
                Ritop = zk    *   dbsfc( nz, node )
                Rib_k = Ritop / ( dVsq ( nz, node ) + Vtsq + epsln )
                dzup  = zk    -   zkm1

                if (Rib_k > Ricr) then
                    ! linearly interpolate to find hbl where Rib = Ricr
                    hbl(node) = zkm1 + dzup*(Ricr-Rib_km1)/(Rib_k-Rib_km1+epsln)
                    kbl(node) = nz
                    exit
                else
                    Rib_km1 = Rib_k
                end if

                ! stability + buoyancy forcing for the boundary layer
                if (use_sw_pene) then
                    ! linear interpolation of sw_3d to depth hbl
                    bfsfc(node) = Bo(node) + &
                                  coeff_sw * &
                                  ( sw_3d(nzmin,node) - &
                                                  ( sw_3d(nz-1,node) + &
                                                                     ( sw_3d(nz,node) - sw_3d(nz-1,node) ) * ( hbl(node) - zkm1 ) / dzup &
                                                  ) &
                                  )
                    stable(node) = 0.5_WP + SIGN( 0.5_WP, bfsfc(node) )
                    bfsfc (node) = bfsfc(node) + stable(node) * epsln   ! ensures bfsfc never = 0
                end if
            end do

            ! check hbl limits for hekman or hmonob (eqn 24); no ekman/monob under cavity
            if (bfsfc(node) > 0.0_WP .and. nzmin==1) then
                hekman = cekman * ustar(node) / MAX( ABS( mesh%coriolis_node(node) ), epsln)
                hmonob = cmonob * ustar(node) * ustar(node) * ustar(node)     &
                       /vonk / (bfsfc(node) + epsln)
                hlimit = stable(node) * AMIN1( hekman, hmonob )
                hbl(node) = AMIN1( hbl(node), hlimit )
                hbl(node) = MAX( hbl(node), ABS(mesh%zbar_3d_n(2,node)) )
            end if
        end do

        ! smooth_hbl = .false. -> no exchange/smoothing pass (dead under work_core)

        do node = 1, nNodO
            nzmax = mesh%nlevels_nod2D(node)
            nzmin = mesh%ulevels_nod2D(node)
            ! find new kbl
            kbl(node) = nzmax
            do nz = nzmin+1, nzmax
                if (ABS(mesh%zbar_3d_n(nz,node)) > hbl(node)) then
                    kbl(node) = nz
                    exit
                end if
            end do

            ! stability + buoyancy forcing for the final hbl
            if (use_sw_pene) then
                coeff_sw = g * sw_alpha(nzmin,node)   ! @ the surface @ Z (m/s2/K)
                ! linear interpolation of sw_3d to depth hbl
                bfsfc(node) = Bo(node) + &
                              coeff_sw * &
                              ( sw_3d(nzmin,node) - &
                                              ( sw_3d(kbl(node)-1, node) + &
                                                                         ( sw_3d(kbl(node), node) - sw_3d(kbl(node)-1, node) ) &
                                                                         * ( hbl(node) + mesh%zbar_3d_n( kbl(node)-1,node) ) &
                                                                         / ( mesh%zbar_3d_n( kbl(node)-1,node) - mesh%zbar_3d_n(kbl(node),node) ) ) )
                stable(node) = 0.5_WP + SIGN(0.5_WP, bfsfc(node))
                bfsfc(node)  = bfsfc(node) + stable(node) * epsln
            end if

            ! determine caseA (=1 if hbl is above the mid point of level kbl, else 0)
            dzup        = mesh%zbar_3d_n(kbl(node)-1,node) - mesh%zbar_3d_n(kbl(node),node)
            caseA(node) = 0.5_WP + SIGN( 0.5_WP, ABS( mesh%zbar_3d_n(kbl(node),node) ) - 0.5_WP * dzup - hbl(node) )
        end do
    end subroutine bldepth

    !#######################################################################
    ! Boundary-layer mixing coefficients blmc(:,:,1:3) (momentum/T/S) + the
    ! kbl-1 diffusivities dkm1(:,1:3) + the nonlocal-transport counter-gradient
    ! flux ghats (= blmix_kpp, oce_ale_mixing_kpp.F90:1228). Within the OBL the
    ! shape function G(sigma) (eqn 11) modulates the surface velocity scale
    ! (wscale) to give the diffusivity profile (eqn 10); ghats is the nonlocal
    ! term (eqn 20, nonzero only for scalars in unstable forcing).
    !
    ! Channel wiring: diff_col(:,1)=viscA (momentum) -> blmc(:,:,1)/dkm1(:,1);
    ! diff_col(:,2)=diffK(:,:,1) (T) -> blmc(:,:,2)/dkm1(:,2); diff_col(:,3)=
    ! diffK(:,:,2) (S) -> blmc(:,:,3)/dkm1(:,3). Inputs viscA/diffK are the
    ! ri_iwmix interior coeffs (M5a-2); hbl/bfsfc/stable/caseA/kbl/ustar are the
    ! bldepth prestep+outputs (M5a-3) — all byte-proven, so blmix reads them
    ! directly. Loops OWNED nodes (the oracle's myDim_nod2D); blmc is pre-zeroed
    ! over OWNED+HALO (matches the oracle); ghats/dkm1 rely on the caller's zero
    ! (the oracle's oce_mixing_kpp_init zero). Optional-partit from the start.
    !
    ! L29: every divide here (/dthick, /(hbl+epsln), /(wm+epsln), /(ws+epsln))
    ! is SCALAR per-node (the inner nz-loop EXITs at kbl + calls wscale) -> no
    ! SIMD-divide trap; explicit-shape dummies + AMIN1/MIN/INT/ABS kept VERBATIM.
    subroutine blmix_kpp(viscA, diffK, hbl, ustar, bfsfc, stable, caseA, kbl, &
                         blmc, ghats, dkm1, mesh, partit)
        use mod_mesh,        only: t_mesh
        use mod_partit,      only: t_partit
        use mod_part_bounds, only: owned_bounds
        real(kind=WP), dimension(:,:),   intent(in)    :: viscA   ! (nl,   nnod)
        real(kind=WP), dimension(:,:,:), intent(in)    :: diffK   ! (nl,   nnod, 2)
        real(kind=WP), dimension(:),     intent(in)    :: hbl     ! (nnod)
        real(kind=WP), dimension(:),     intent(in)    :: ustar   ! (nnod)
        real(kind=WP), dimension(:),     intent(in)    :: bfsfc   ! (nnod)
        real(kind=WP), dimension(:),     intent(in)    :: stable  ! (nnod)
        real(kind=WP), dimension(:),     intent(in)    :: caseA   ! (nnod)
        integer,       dimension(:),     intent(in)    :: kbl     ! (nnod)
        real(kind=WP), dimension(:,:,:), intent(inout) :: blmc    ! (nl,   nnod, 3)
        real(kind=WP), dimension(:,:),   intent(inout) :: ghats   ! (nl-1, nnod)
        real(kind=WP), dimension(:,:),   intent(inout) :: dkm1    ! (nnod, 3)
        type(t_mesh),   intent(in), target   :: mesh
        type(t_partit), intent(in), optional :: partit

        integer :: node, nz, kn, knm1, knp1, nl1, nu1
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: delhat, R, dvdzup, dvdzdn
        real(kind=WP) :: viscp, difsp, diftp, visch, difsh, difth, f1
        real(kind=WP) :: sig, a1, a2, a3, Gm, Gs, Gt
        real(kind=WP) :: sigma, zehat, wm, ws
        real(kind=WP) :: gat1m, gat1t, gat1s, dat1m, dat1s, dat1t
        real(kind=WP) :: dthick(mesh%nl), diff_col(mesh%nl,3)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do node = 1, nNodL
            blmc(:, node, :) = 0.0_WP
        end do

        do node = 1, nNodO
            nl1 = mesh%nlevels_nod2D(node)
            nu1 = mesh%ulevels_nod2D(node)

            if (nl1     < 3) cycle   ! a temporary solution
            if (nl1-nu1 < 2) cycle

            dthick(nu1+1:nl1-1) = 0.5_WP*( mesh%hnode(nu1:nl1-2,node) + mesh%hnode(nu1+1:nl1-1,node) )
            dthick(nu1)         = mesh%hnode(nu1,    node)*0.5_WP
            dthick(nl1)         = mesh%hnode(nl1-1,  node)*0.5_WP

            diff_col(nu1:nl1-1,1)   = viscA(nu1:nl1-1,node)
            diff_col(nu1:nl1-1,2:3) = diffK(nu1:nl1-1,node,:)
            diff_col(nl1,:)         = diff_col(nl1-1,:)

            ! velocity scales at hbl (recall epsilon_kpp=0.1)
            sigma = stable(node) * 1.0_WP + (1.0_WP-stable(node)) * epsilon_kpp
            zehat = vonk * sigma * hbl(node) * bfsfc(node)
            call wscale(zehat, ustar(node), wm, ws)

            kn   = INT(caseA(node)+epsln) *(kbl(node) -1) + &
                   (1-INT(caseA(node)+epsln)) * kbl(node)
            kn   = MIN(kn,nl1-1)
            knm1 = MAX(kn-1,nu1)
            knp1 = MIN(kn+1,nl1)

            ! interior viscosities + derivatives at hbl (eqn 18)
            delhat = ABS(mesh%Z_3d_n(kn,node)) - hbl(node)
            R      = 1.0_WP - delhat / dthick(kn)

            dvdzup = (diff_col(knm1,1) - diff_col(kn,1))/dthick(kn)
            dvdzdn = (diff_col(kn,1) - diff_col(knp1,1))/dthick(knp1)
            viscp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+ &
                                R  * (dvdzdn + abs(dvdzdn)) )

            dvdzup = (diff_col(knm1,3) - diff_col(kn,3))/dthick(kn)
            dvdzdn = (diff_col(kn,3) - diff_col(knp1,3))/dthick(knp1)
            difsp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+ &
                                R  * (dvdzdn + ABS(dvdzdn)) )

            dvdzup = (diff_col(knm1,2) - diff_col(kn,2))/dthick(kn)
            dvdzdn = (diff_col(kn,2) - diff_col(knp1,2))/dthick(knp1)
            diftp  = 0.5_WP * ( (1.0_WP - R) * (dvdzup + ABS(dvdzup))+ &
                                R  * (dvdzdn + ABS(dvdzdn)) )

            visch  = diff_col(kn,1) + viscp * delhat
            difsh  = diff_col(kn,3) + difsp * delhat
            difth  = diff_col(kn,2) + diftp * delhat

            f1 = stable(node) * conc1 * bfsfc(node) / (ustar(node)**4+epsln)

            gat1m = visch / (hbl(node) + epsln) / (wm + epsln)
            dat1m = -viscp / (wm+epsln) + f1 * visch
            dat1m = min(dat1m, 0.0_WP)

            gat1s = difsh  / (hbl(node) + epsln) / (ws + epsln)
            dat1s = -difsp / (ws+epsln) + f1 * difsh
            dat1s = min(dat1s, 0.0_WP)

            gat1t = difth /  (hbl(node) + epsln) / (ws + epsln)
            dat1t = -diftp / (ws+epsln) + f1 * difth
            dat1t = min(dat1t, 0.0_WP)

            do nz = nu1+1, nl1-1
                if (nz >= kbl(node)) exit

                ! turbulent velocity scales on the interfaces
                sig   = ABS(mesh%Z_3d_n(nz,node)) / (hbl(node)+epsln)
                sigma = stable(node) * sig &
                      + (1.0_WP - stable(node)) * AMIN1(sig, epsilon_kpp)
                zehat = vonk * sigma * hbl(node) * bfsfc(node)
                call wscale(zehat, ustar(node), wm, ws)

                ! dimensionless shape functions at the interfaces (eqn 11)
                a1 = sig    - 2.0_WP
                a2 = 3.0_WP - 2.0_WP * sig
                a3 = sig    - 1.0_WP

                Gm = a1 + a2 * gat1m + a3 * dat1m
                Gs = a1 + a2 * gat1s + a3 * dat1s
                Gt = a1 + a2 * gat1t + a3 * dat1t

                ! boundary layer diffusivities at the interfaces (eqn 10)
                blmc(nz,node,1) = hbl(node) * wm * sig * (1.0_WP + sig * Gm)
                blmc(nz,node,2) = hbl(node) * ws * sig * (1.0_WP + sig * Gt)
                blmc(nz,node,3) = hbl(node) * ws * sig * (1.0_WP + sig * Gs)

                ! nonlocal transport term = ghats * <ws>o (eqn 20)
                ghats(nz,node) = (1.0_WP - stable(node)) * cg &
                               / (ws * hbl(node) + epsln)
            end do

            ! diffusivities at kbl-1 grid level
            sig   = ABS(mesh%zbar_3d_n(kbl(node)-1,node)) / (hbl(node) + epsln)
            sigma = stable(node) * sig &
                  + (1.0_WP - stable(node)) * MIN(sig, epsilon_kpp)
            zehat = vonk * sigma * hbl(node) * bfsfc(node)
            call wscale(zehat, ustar(node), wm, ws)

            a1 = sig    - 2.0_WP
            a2 = 3.0_WP - 2.0_WP * sig
            a3 = sig    - 1.0_WP

            Gm = a1 + a2 * gat1m + a3 * dat1m
            Gs = a1 + a2 * gat1s + a3 * dat1s
            Gt = a1 + a2 * gat1t + a3 * dat1t

            dkm1(node,1) = hbl(node) * wm * sig * (1.0_WP + sig * Gm)
            dkm1(node,2) = hbl(node) * ws * sig * (1.0_WP + sig * Gt)
            dkm1(node,3) = hbl(node) * ws * sig * (1.0_WP + sig * Gs)
        end do
    end subroutine blmix_kpp

    !#######################################################################
    ! Enhance the boundary-layer diffusivity at the kbl-0.5 interface
    ! (= enhance, oce_ale_mixing_kpp.F90:1419). Blends the interior diffusivity
    ! at kbl-1 with the BL value via the fractional depth delta, and scales the
    ! nonlocal ghats by (1-caseA). Modifies blmc(kbl-1,:,1:3) + ghats(kbl-1,:).
    ! Loops OWNED nodes. One scalar divide (delta) -> no L29 SIMD trap.
    subroutine enhance(viscA, diffK, hbl, caseA, kbl, blmc, ghats, dkm1, mesh, partit)
        use mod_mesh,        only: t_mesh
        use mod_partit,      only: t_partit
        use mod_part_bounds, only: owned_bounds
        real(kind=WP), dimension(:,:),   intent(in)    :: viscA   ! (nl,   nnod)
        real(kind=WP), dimension(:,:,:), intent(in)    :: diffK   ! (nl,   nnod, 2)
        real(kind=WP), dimension(:),     intent(in)    :: hbl     ! (nnod)
        real(kind=WP), dimension(:),     intent(in)    :: caseA   ! (nnod)
        integer,       dimension(:),     intent(in)    :: kbl     ! (nnod)
        real(kind=WP), dimension(:,:,:), intent(inout) :: blmc    ! (nl,   nnod, 3)
        real(kind=WP), dimension(:,:),   intent(inout) :: ghats   ! (nl-1, nnod)
        real(kind=WP), dimension(:,:),   intent(in)    :: dkm1    ! (nnod, 3)
        type(t_mesh),   intent(in), target   :: mesh
        type(t_partit), intent(in), optional :: partit

        integer :: node, k, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: delta, dkmp5, dstar

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do node = 1, nNodO
            k     = kbl(node) - 1
            delta = (hbl(node) + mesh%zbar_3d_n(k,node)) &
                  / (mesh%zbar_3d_n(k,node) - mesh%zbar_3d_n(k+1,node))

            ! momentum
            dkmp5 = caseA(node) * viscA(k,node) &
                  + ( 1.0_WP - caseA(node) ) * blmc( k, node, 1 )
            dstar = ( 1.0_WP - delta )**2 * dkm1( node, 1 ) + delta**2 * dkmp5
            blmc( k, node, 1 ) = (1.0_WP - delta) * viscA(k,node) &
                  + delta * dstar

            ! temperature
            dkmp5 = caseA(node) * diffK(k,node,1) &
                  + ( 1.0_WP - caseA(node) ) * blmc( k, node, 2 )
            dstar = ( 1.0_WP - delta )**2 * dkm1( node, 2 ) + delta**2 * dkmp5
            blmc( k, node, 2 ) = ( 1.0_WP - delta ) * diffK( k, node, 1) &
                  + delta * dstar

            ! salinity
            dkmp5 = caseA(node) * diffK(k,node,2) &
                  + ( 1.0_WP - caseA(node) ) * blmc( k, node, 3 )
            dstar = ( 1.0_WP - delta )**2 * dkm1( node, 3 ) + delta**2 * dkmp5
            blmc( k, node, 3 ) = ( 1.0_WP - delta ) * diffK( k, node, 2 ) &
                  + delta * dstar

            ghats(k,node) = (1.0_WP-caseA(node)) * ghats(k,node)
        end do
    end subroutine enhance

end module oce_mixing_kpp
