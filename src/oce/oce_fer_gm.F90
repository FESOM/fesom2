module oce_fer_gm
    ! M4b — Gent-McWilliams eddy bolus parameterization (Ferrari et al. 2010).
    ! Faithful transcription of FESOM2 v2.7.3 oce_fer_gm.F90 (init_Redi_GM :215,
    ! fer_solve_Gamma :40, fer_gamma2vel :168), the explicit-dataflow + optional-partit
    ! pattern (absent partit => the proven 1-rank path verbatim; present+npes>1 => owned
    ! loops + the FESOM2 exchanges).
    !
    ! SCOPE: ONLY the work_core GM path is ported (Fer_GM=.true., Redi=.false.,
    ! scaling_resolution + scaling_GMzexp, K_GM_Ktaper=.false.). The dead branches
    ! (scaling_Ferreira/Rossby/GINsea/FESOM14, K_GM_Ktaper, the Redi Ki path) are NOT
    ! transcribed — they are .false. in work_core/namelist.oce, so unexercised and
    ! ungateable here. init_Redi_GM asserts those flags are off (error stop) so a future
    ! config that flips one fails loudly rather than silently computing a wrong K. The
    ! Redi Ki path lands in M4d.
    !
    ! Byte-exactness notes (L29 SIMD-divide watch-list, plan 2026-06-23-m4-gm-redi.md):
    !  - fer_solve_Gamma mirrors the oracle's `tr => fer_gamma(:,:,n)` POINTER so the
    !    Thomas-sweep divides keep the same (scalar) codegen as the oracle (the L29 trap:
    !    a derived-type/explicit array the compiler can prove distinct may vectorize a
    !    divide the oracle compiles scalar). The TDMA is a serial recurrence anyway.
    !  - fer_gamma2vel `zinv=1/helem` and init_Redi_GM `bvfreq/bvref` (dead) /
    !    GMzexp exp-divide are scalar per node/level; gated, add !DIR$ NOVECTOR only if a
    !    field drifts ~1 ULP.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_constants,   only: density_0, g, pi
    use mod_param_phys,  only: Fer_GM, Redi, Redi_Kmax, Redi_Kmin, Redi_Ktaper, K_hor, &
                               K_GM_max, K_GM_min, K_GM_cm, K_GM_cmin, K_GM_resscalorder, &
                               K_GM_rampmax, K_GM_rampmin, K_GM_Ktaper, &
                               scaling_Ferreira, scaling_Rossby, scaling_resolution, &
                               scaling_FESOM14, scaling_GMzexp, scaling_GINsea, &
                               GMzexp_zref, GMzexp_smin
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod, exchange_elem
    implicit none
    private
    public :: init_Redi_GM, fer_solve_Gamma, fer_gamma2vel

contains

    !===========================================================================
    subroutine init_Redi_GM(mesh, bvfreq, fer_K, fer_c, fer_scal, partit, Ki, fer_tapfac)
        ! GM diffusivity K = F1(x,y)*F2(z) (FESOM2 oce_fer_gm.F90:215). F1: baroclinic
        ! gravity-wave speed cm (-> fer_c=cm^2) + resolution/area scaling -> fer_scal,
        ! fer_K(nzmin)=fer_scal*K_GM_max (clamped to K_GM_min). F2: Ferreira-style vertical
        ! downscaling (work_core: scaling_GMzexp, exp(-|zbar|/zref)) -> fer_K(nz).
        ! fer_K MUST enter init'd to 500.0 (oce_setup_step.F90:964) so the below-bottom
        ! region (levels never written) matches the oracle byte-for-byte.
        type(t_mesh),   intent(in)    :: mesh
        real(kind=WP),  intent(in)    :: bvfreq(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(inout) :: fer_K(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(inout) :: fer_c(mesh%nod2D)
        real(kind=WP),  intent(inout) :: fer_scal(mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        ! M4d Redi: Ki (output) + fer_tapfac (input, from compute_neutral_slope). Present iff
        ! Redi=.true. (the GM-only callers omit them; the Redi blocks are if(Redi)-guarded).
        real(kind=WP),  intent(inout), optional :: Ki(mesh%nl-1, mesh%nod2D)
        real(kind=WP),  intent(in),    optional :: fer_tapfac(mesh%nl-1, mesh%nod2D)
        integer       :: n, k, nz, nzmax, nzmin, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: reso, cm, scaling
        real(kind=WP) :: zscaling(mesh%nl)
        real(kind=WP) :: refscalresol = 100000._WP   ! 100 km

        ! M4b/M4d port the work_core GM + Redi paths; assert the unsupported scalings are off.
        if (scaling_Ferreira .or. scaling_Rossby .or. scaling_GINsea .or. &
            scaling_FESOM14 .or. K_GM_Ktaper) then
            error stop 'oce_fer_gm::init_Redi_GM: only the work_core path &
                &(scaling_resolution + scaling_GMzexp, no GM Ktaper) is ported. &
                &Ferreira/Rossby/GINsea/FESOM14/Ktaper not implemented.'
        end if

        ! Redi: synchronise Redi_Kmax with K_GM_max when Redi_Kmax<=0 (FESOM2 oce_fer_gm.F90:240)
        if (Redi .and. Redi_Kmax <= 0) Redi_Kmax = K_GM_max

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !******************************* F1(x,y) ******************************
        do n = 1, nNodO
            ! min levels below / max ulevels above elements around node n (explicit loop,
            ! NOT minval/maxval — the oracle's deliberate Intel-vectorisation workaround :254)
            nzmax = mesh%nl
            nzmin = 1
            do k = 1, mesh%nod_in_elem2D_num(n)
                nzmax = min(nzmax, mesh%nlevels(mesh%nod_in_elem2D(k, n)))
                nzmin = max(nzmin, mesh%ulevels(mesh%nod_in_elem2D(k, n)))
            end do
            reso = mesh%mesh_resolution(n)
            if (Fer_GM) then
                ! baroclinic gravity-wave speed cm (abs() guards the -0 cray case, :266)
                cm = 0._WP
                do nz = nzmin, nzmax-1
                    cm = cm + mesh%hnode_new(nz,n)*(sqrt(abs(max(bvfreq(nz,n),   0._WP))) &
                                                  + sqrt(abs(max(bvfreq(nz+1,n), 0._WP))))/2._WP
                end do
                cm = max(K_GM_cmin, cm/pi/K_GM_cm)   ! mth baroclinic speed, floored
                scaling = 1._WP

                ! scale K_GM with resolution (referenced to 100 km); work_core resscalorder=2
                if     (scaling_resolution .and. K_GM_resscalorder==1) then
                    scaling = scaling*(reso/refscalresol)**K_GM_resscalorder
                elseif (scaling_resolution .and. K_GM_resscalorder==2) then
                    scaling = scaling*(mesh%area(n)/(refscalresol**2)*2)**(1/K_GM_resscalorder)
                end if

                ! resolution ramp (work_core K_GM_rampmax=-1 => always skipped; reso>0).
                ! denom floored at 1.e-12 — the oracle's deliberate Intel divide-by-0 guard :308
                if (reso/1000.0_WP < K_GM_rampmax) then
                    scaling = scaling*max((reso/1000.0_WP-K_GM_rampmin) &
                                          /max(K_GM_rampmax-K_GM_rampmin, 1.e-12_WP), 0._WP)
                end if

                fer_scal(n)    = min(scaling, 1.0_WP)
                fer_K(nzmin,n) = fer_scal(n)*K_GM_max     ! surface template
                fer_K(nzmin,n) = max(fer_K(nzmin,n), K_GM_min)
                fer_c(n)       = cm*cm
            end if

            ! Redi diffusivity surface template (FESOM2 oce_fer_gm.F90:337-346). work_core
            ! resscalorder=2 + K_hor=0 -> Ki(nzmin)=0**(..) =0 here; OVERWRITTEN by the
            ! Redi.and.Fer_GM coupling below. Note the oracle's literal K_hor** (exponentiation).
            if (Redi) then
                if     (K_GM_resscalorder==1) then
                    Ki(nzmin,n) = K_hor*(reso/refscalresol)**K_GM_resscalorder
                elseif (K_GM_resscalorder==2) then
                    Ki(nzmin,n) = K_hor**(mesh%area(n)/(refscalresol**2)*2)**(1/K_GM_resscalorder)
                end if
            end if
        end do

        ! Like FESOM 1.4, make Redi follow GM: Ki(surface) = fer_scal*Redi_Kmax (clamped).
        ! (FESOM2 oce_fer_gm.F90:351-359; the dominant Ki source — K_hor=0 zeroes the resol term.)
        ! nzmin = ulevels_nod2D_max(n) (= the F1 per-node template index; 1 with no cavity).
        if (Redi .and. Fer_GM) then
            do n = 1, nNodO
                nzmin = mesh%ulevels_nod2D_max(n)
                Ki(nzmin, n) = max(fer_scal(n)*Redi_Kmax, K_GM_min)
            end do
        end if

        !*************** F2(z) (Ferreira et al. 2005) — work_core: GMzexp ***************
        do n = 1, nNodO
            if (Redi .or. Fer_GM) then
                nzmax = mesh%nlevels_nod2D(n)
                nzmin = mesh%ulevels_nod2D(n)
                zscaling = 1.0_WP
                ! work_core: scaling_GMzexp -> exp(-|zbar|/zref) depth downscaling
                if (scaling_GMzexp) then
                    do nz = nzmin, nzmax
                        zscaling(nz) = GMzexp_smin + (1.0_WP-GMzexp_smin) &
                                       *exp(-abs(mesh%zbar_3d_n(nz,n)/GMzexp_zref))
                        zscaling(nz) = max(min(zscaling(nz), 1.0_WP), GMzexp_smin)
                    end do
                end if
            end if

            ! vertical Ferreira scaling of the GM diffusivity (K_GM_Ktaper=.false. => no taper)
            if (Fer_GM) then
                do nz = nzmin+1, nzmax
                    fer_K(nz,n) = fer_K(nzmin,n)*zscaling(nz)
                end do
                fer_K(nzmin,n) = fer_K(nzmin,n)*zscaling(nzmin)
            end if

            ! vertical Ferreira scaling of the Redi diffusivity + Redi_Ktaper sqrt-split of the
            ! slope taper between Ki and the neutral slope (FESOM2 oce_fer_gm.F90:479-498).
            if (Redi) then
                do nz = nzmin+1, nzmax-1
                    Ki(nz,n) = Ki(nzmin,n)*0.5_WP*(zscaling(nz)+zscaling(nz+1))
                end do
                Ki(nzmin,n) = Ki(nzmin,n)*0.5_WP*(zscaling(nzmin)+zscaling(nzmin+1))
                if (Redi_Ktaper) then
                    do nz = nzmin, nzmax-1
                        Ki(nz,n) = Ki(nz,n)*sqrt(fer_tapfac(nz,n)) &
                                 + Redi_Kmin*abs(sqrt(fer_tapfac(nz,n))-1)
                    end do
                end if
            end if
        end do

        if (present(partit)) then
            if (is_multirank(partit)) then
                if (Fer_GM) call exchange_nod(fer_c, partit)
                if (Fer_GM) call exchange_nod(fer_K, partit)
                if (Redi)   call exchange_nod(Ki, partit)
            end if
        end if
    end subroutine init_Redi_GM

    !===========================================================================
    subroutine fer_solve_Gamma(mesh, sigma_xy, bvfreq, fer_c, fer_K, fer_gamma, partit)
        ! Per-node Thomas (TDMA) solve of (c*Gamma_zz - N^2*Gamma) = (g/rho0)*K*grad(sigma)
        ! -> fer_gamma(1:2,nz,n), the GM streamfunction (FESOM2 oce_fer_gm.F90:40). zbar_n /
        ! Z_n are rebuilt bottom-up from hnode_new exactly as the oracle (the FP accumulation
        ! order matters); zbar_n_bot(n) == mesh%zbar_3d_n(nlevels_nod2D(n),n) (FESOM2
        ! oce_ale.F90:550). `tr` is a POINTER into fer_gamma (mirror the oracle so the sweep
        ! divides keep its scalar codegen, L29).
        type(t_mesh),   intent(in)            :: mesh
        real(kind=WP),  intent(in)            :: sigma_xy(2, mesh%nl-1, mesh%nod2D)
        real(kind=WP),  intent(in)            :: bvfreq(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(in)            :: fer_c(mesh%nod2D)
        real(kind=WP),  intent(in)            :: fer_K(mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(inout), target :: fer_gamma(2, mesh%nl, mesh%nod2D)
        type(t_partit), intent(in), optional  :: partit
        integer       :: nz, n, nzmax, nzmin, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: zinv1, zinv2, zinv, m, r
        real(kind=WP) :: a(mesh%nl), b(mesh%nl), c(mesh%nl)
        real(kind=WP) :: cp(mesh%nl), tp(2,mesh%nl)
        real(kind=WP) :: zbar_n(mesh%nl), z_n(mesh%nl-1)
        real(kind=WP), dimension(:,:), pointer :: tr

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodO
            tr => fer_gamma(:,:,n)
            nzmax = mesh%nlevels_nod2D(n)
            nzmin = mesh%ulevels_nod2D(n)

            ! Z_n / zbar_n for the current ALE step, bottom-up from the partial-cell bottom
            zbar_n = 0.0_WP
            z_n    = 0.0_WP
            zbar_n(nzmax) = mesh%zbar_3d_n(nzmax, n)                       ! = zbar_n_bot(n)
            z_n(nzmax-1)  = zbar_n(nzmax) + mesh%hnode_new(nzmax-1,n)/2.0_WP
            do nz = nzmax-1, nzmin+1, -1
                zbar_n(nz) = zbar_n(nz+1) + mesh%hnode_new(nz,n)
                z_n(nz-1)  = zbar_n(nz)   + mesh%hnode_new(nz-1,n)/2.0_WP
            end do
            zbar_n(nzmin) = zbar_n(nzmin+1) + mesh%hnode_new(nzmin,n)

            ! min levels below / max ulevels above elements around node n
            nzmax = mesh%nlevels_nod2D_min(n)
            nzmin = mesh%ulevels_nod2D_max(n)

            ! FESOM2 58ccb44: a node touching BOTH deep and shallow cavity elements can
            ! invert this range (nzmin >= nzmax). The sweeps below would then run zero
            ! times and the back substitution would publish an UNASSIGNED tp element into
            ! fer_gamma. Zero the column and skip. Inert without a cavity (ulevels_nod2D_max
            ! == 1 everywhere, so nzmin < nzmax), but bottom-at-vertices lowers
            ! nlevels_nod2D_min relative to nlevels_nod2D, which moves TOWARDS this case.
            if (nzmin >= nzmax) then
                tr(:,:) = 0.0_WP
                cycle
            end if

            ! tridiagonal coefficients
            c(nzmin) = 0.0_WP
            a(nzmin) = 0.0_WP
            b(nzmin) = 1.0_WP
            zinv2 = 1.0_WP/(zbar_n(nzmin)-zbar_n(nzmin+1))
            do nz = nzmin+1, nzmax-1
                zinv1 = zinv2
                zinv2 = 1.0_WP/(zbar_n(nz)-zbar_n(nz+1))
                zinv  = 1.0_WP/(z_n(nz-1)-z_n(nz))
                a(nz) = fer_c(n)*zinv1*zinv
                c(nz) = fer_c(n)*zinv2*zinv
                b(nz) = -a(nz)-c(nz)-max(bvfreq(nz,n), 1.e-8_WP)
            end do
            nz = nzmax
            c(nz) = 0.0_WP
            a(nz) = 0.0_WP
            b(nz) = 1.0_WP

            ! rhs
            tr(:, nzmin) = 0.0_WP
            tr(:, nzmax) = 0.0_WP
            do nz = nzmin+1, nzmax-1
                r = g/density_0
                tr(1, nz) = r*0.5_WP*sum(sigma_xy(1,nz-1:nz,n))*fer_K(nz, n)
                tr(2, nz) = r*0.5_WP*sum(sigma_xy(2,nz-1:nz,n))*fer_K(nz, n)
            end do

            ! Thomas sweep
            cp(nzmin)   = c(nzmin)/b(nzmin)
            tp(:,nzmin) = tr(:,nzmin)/b(nzmin)
            do nz = nzmin+1, nzmax
                m = b(nz)-cp(nz-1)*a(nz)
                cp(nz) = c(nz)/m
                tp(:,nz) = (tr(:,nz)-tp(:,nz-1)*a(nz))/m
            end do
            tr(:,nzmax) = tp(:,nzmax)
            do nz = nzmax-1, nzmin, -1
                tr(:,nz) = tp(:,nz)-cp(nz)*tr(:,nz+1)
            end do
        end do

        if (present(partit)) then
            if (is_multirank(partit)) call exchange_nod(fer_gamma, partit)
        end if
    end subroutine fer_solve_Gamma

    !===========================================================================
    subroutine fer_gamma2vel(mesh, fer_gamma, fer_uv, partit)
        ! GM bolus horizontal velocity fer_uv(1:2,nz,el) = (1/h)*sum_k elem_w(k)*(Gamma_nz-Gamma_nz+1)
        ! over the 4 vertex slots (FESOM2 oce_fer_gm.F90:168 had the plain 1/3 mean). fer_uv MUST
        ! enter zeroed below-bottom (only nzmin..nzmax-1 written). zinv=1/helem is the
        ! SIMD-divide watch item.
        type(t_mesh),   intent(in)    :: mesh
        real(kind=WP),  intent(in)    :: fer_gamma(2, mesh%nl, mesh%nod2D)
        real(kind=WP),  intent(inout) :: fer_uv(2, mesh%nl-1, mesh%elem2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: nz, nzmax, nzmin, el, elnod(MAX_NV), nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: zinv

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do el = 1, nElemO
            elnod = mesh%elem2D_nodes(:,el)
            nzmax = mesh%nlevels(el)
            nzmin = mesh%ulevels(el)
            do nz = nzmin, nzmax-1
                zinv = 1.0_WP/mesh%helem(nz,el)
                fer_uv(1,nz,el) = sum(mesh%elem_w(:,el)*(fer_gamma(1,nz,elnod)-fer_gamma(1,nz+1,elnod)))*zinv
                fer_uv(2,nz,el) = sum(mesh%elem_w(:,el)*(fer_gamma(2,nz,elnod)-fer_gamma(2,nz+1,elnod)))*zinv
            end do
        end do

        if (present(partit)) then
            if (is_multirank(partit)) call exchange_elem(fer_uv, partit)
        end if
    end subroutine fer_gamma2vel

end module oce_fer_gm
