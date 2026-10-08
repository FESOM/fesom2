module oce_mixing_tke
    ! FESOM TKE (Gaspar et al. 1990) prognostic turbulent-kinetic-energy vertical mixing,
    ! with the Blanke & Delecluse (1993) mixing length. In FESOM2 v2.7.3 this scheme is
    ! reachable ONLY through CVMix (oce_ale.F90:3749, #if defined(__cvmix)); M7 transcribes
    ! the CVMix TKE column math DIRECTLY into FESOM3 — library-free, no tke_type put/get, no
    ! diagnostics-output bloat. The prognostic column solver (integrate_tke) + the Thomas
    ! tridiagonal solve (solve_tridiag) are 1:1 from the oracle:
    !   port2/fesom2/src/cvmix_driver/cvmix_tke.F90        (integrate_tke, :415-987)
    !   port2/fesom2/src/cvmix_driver/cvmix_utils_addon.F90 (solve_tridiag, :116-149)
    !   port2/fesom2/src/cvmix_driver/gen_modules_cvmix_tke.F90 (init_cvmix_tke + calc_cvmix_tke)
    !
    ! Target config (work_*_tke): mix_scheme='cvmix_TKE' (mix_scheme_nmb==5) -> tke_only=.true.
    ! (pure standalone TKE; IDEMIX uncoupled, tidal mixing off). tke_mxl_choice==2 is the only
    ! ported mixing-length branch (Blanke-Delecluse). The namelist &param_tke DOUBLES are baked
    ! (L38): tke_cd=3.75 (namelist OVERRIDES the code default 1.0; used at the Neumann surface BC
    ! despite the "Dirichlet" comment), tke_c_k=0.1, tke_c_eps=0.7, tke_alpha=30, tke_mxl_min=1e-8,
    ! tke_kappaM_max=100, tke_surf_min=1e-4, tke_min=1e-6.
    !
    ! Byte-exactness (the M5a-2 Fortran<->Fortran lesson): under the -r8 anchor build the oracle
    ! literals (6.6, sqrt(2D0), **(3./2.), 0.d0, 1d-12, MAX/MIN) compile as doubles, so they are
    ! transcribed VERBATIM — do NOT "simplify" them. The C port's 6.6/pow/fmax traps do NOT recur
    ! here PROVIDED nothing is changed.
    !
    ! M7a-1 scope (this commit): tke_init (param store) + integrate_tke (the column core) +
    ! solve_tridiag. The driver (calc_cvmix_tke) is added M7a-2.
    use mod_precision, only: WP
    implicit none
    private

    public :: tke_init
    public :: integrate_tke
    public :: calc_cvmix_tke

    ! Module-saved TKE constants (set once by tke_init; mirror cvmix_tke.F90
    ! tke_constants_saved, minus the put/get string-dispatch plumbing).
    real(kind=WP), save :: tkep_c_k
    real(kind=WP), save :: tkep_c_eps
    real(kind=WP), save :: tkep_cd
    real(kind=WP), save :: tkep_alpha_tke
    real(kind=WP), save :: tkep_KappaM_min
    real(kind=WP), save :: tkep_KappaM_max
    real(kind=WP), save :: tkep_mxl_min
    real(kind=WP), save :: tkep_tke_min
    real(kind=WP), save :: tkep_tke_surf_min
    integer,       save :: tkep_mxl_choice
    logical,       save :: tkep_only_tke
    logical,       save :: tkep_l_lc
    logical,       save :: tkep_use_ubound_dirichlet
    logical,       save :: tkep_use_lbound_dirichlet

contains

!=================================================================================
! Store the TKE parameters (the &param_tke DOUBLES + the gate logicals) into the
! module-saved constants read by integrate_tke. Mirrors init_cvmix_tke -> init_tke
! (gen_modules_cvmix_tke.F90:357-371) without the CVMix put/get range-checks.
subroutine tke_init(c_k, c_eps, cd, alpha_tke, mxl_min, KappaM_min, KappaM_max, &
                    tke_min, tke_surf_min, tke_mxl_choice, only_tke, &
                    use_ubound_dirichlet, use_lbound_dirichlet, l_lc)
    real(kind=WP), intent(in) :: c_k, c_eps, cd, alpha_tke, mxl_min
    real(kind=WP), intent(in) :: KappaM_min, KappaM_max, tke_min, tke_surf_min
    integer,       intent(in) :: tke_mxl_choice
    logical,       intent(in) :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet, l_lc

    tkep_c_k                  = c_k
    tkep_c_eps                = c_eps
    tkep_cd                   = cd
    tkep_alpha_tke            = alpha_tke
    tkep_mxl_min              = mxl_min
    tkep_KappaM_min           = KappaM_min
    tkep_KappaM_max           = KappaM_max
    tkep_tke_min              = tke_min
    tkep_tke_surf_min         = tke_surf_min
    tkep_mxl_choice           = tke_mxl_choice
    tkep_only_tke             = only_tke
    tkep_use_ubound_dirichlet = use_ubound_dirichlet
    tkep_use_lbound_dirichlet = use_lbound_dirichlet
    tkep_l_lc                 = l_lc
end subroutine tke_init

!=================================================================================
! Prognostic TKE column solver. 1:1 transcription of cvmix_tke.F90 integrate_tke
! (:415-987): P1 mixing length (Blanke-Delecluse) -> P2 diffusivities -> P3 forcing
! -> P4 implicit tridiag diffusion+dissipation -> P5 floor -> P6 diagnostics.
! Per column, k=1..nlev+1 (nlev = nln-nun+1). dzw(nlev)=hnode(nun:nln); dzt(nlev+1)
! =dz_trr(nun:nln+1). iw_diss/E_iw/alpha_c/tke_plc are passed as length-(nlev+1) ZERO
! arrays (read unconditionally / in dead .not.only_tke / l_lc branches) — the C-port
! NULL-pointer trap. The debug args (i/j/tstep_count/cvmix_int_*), max_nlev,
! bottom_fric, and the vestigial old_KappaM/old_KappaH (never read) are dropped.
subroutine integrate_tke(tke_old, tke_new, KappaM_out, KappaH_out, &
                         dzw, dzt, nlev, Ssqr, Nsqr,               &
                         tke_Tbpr, tke_Tspr, tke_Tdif, tke_Tdis, tke_Twin, &
                         tke_Tiwf, tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr,   &
                         tke_plc, forc_tke_surf, E_iw, dtime,      &
                         iw_diss, forc_rho_surf, rho_ref, grav, alpha_c)

    integer, intent(in)                          :: nlev
    ! OLD values
    real(kind=WP), dimension(nlev+1), intent(in) :: tke_old, dzt
    real(kind=WP), dimension(nlev+1), intent(in) :: Ssqr, Nsqr
    real(kind=WP), dimension(nlev),   intent(in) :: dzw
    ! Langmuir + IDEMIX inputs (zero arrays here)
    real(kind=WP), dimension(nlev+1), intent(in) :: tke_plc
    real(kind=WP), dimension(nlev+1), intent(in) :: E_iw, alpha_c, iw_diss
    real(kind=WP), intent(in)                    :: forc_rho_surf, rho_ref, dtime, grav
    real(kind=WP), intent(in)                    :: forc_tke_surf
    ! NEW values
    real(kind=WP), dimension(nlev+1), intent(out) :: tke_new, KappaM_out, KappaH_out
    ! diagnostics
    real(kind=WP), dimension(nlev+1), intent(out) :: tke_Tbpr, tke_Tspr, tke_Tdif, &
                                                     tke_Tdis, tke_Twin, tke_Tiwf, &
                                                     tke_Tbck, tke_Ttot, tke_Lmix, tke_Pr

    ! local variables
    real(kind=WP), dimension(nlev+1) :: tke_unrest, tke_upd, mxl, sqrttke, prandtl, &
                                        Rinum, K_diss_v, P_diss_v, forc
    real(kind=WP) :: tke_surf, tke_bott
    real(kind=WP) :: alpha_tke, c_eps, cd, KappaM_max, mxl_min, c_k, tke_surf_min, tke_min
    integer       :: tke_mxl_choice
    logical       :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet, l_lc
    real(kind=WP) :: zzw, depth, diff_surf_forc, diff_bott_forc
    real(kind=WP), dimension(nlev+1) :: a_dif, b_dif, c_dif, a_tri, b_tri, c_tri, d_tri, ke
    integer       :: k, kk, kp1
    real(kind=WP) :: kappaM_min

    ! initialize diagnostics
    tke_Tbpr = 0.0
    tke_Tspr = 0.0
    tke_Tdif = 0.0
    tke_Tdis = 0.0
    tke_Twin = 0.0
    tke_Tiwf = 0.0
    tke_Tbck = 0.0
    tke_Ttot = 0.0

    tke_new = 0.0
    tke_upd = 0.0
    tke_surf= 0.0

    a_dif = 0.0
    b_dif = 0.0
    c_dif = 0.0
    a_tri = 0.0
    b_tri = 0.0
    c_tri = 0.0

    !---------------------------------------------------------------------------------
    ! set tke_constants locally
    !---------------------------------------------------------------------------------
    alpha_tke  = tkep_alpha_tke
    c_eps      = tkep_c_eps
    cd         = tkep_cd
    KappaM_max = tkep_KappaM_max
    mxl_min    = tkep_mxl_min
    c_k        = tkep_c_k
    tke_min    = tkep_tke_min
    tke_surf_min   = tkep_tke_surf_min
    tke_mxl_choice = tkep_mxl_choice
    only_tke = tkep_only_tke
    l_lc     = tkep_l_lc
    use_ubound_dirichlet = tkep_use_ubound_dirichlet
    use_lbound_dirichlet = tkep_use_lbound_dirichlet

    kappaM_min = 0.0

    !---------------------------------------------------------------------------------
    ! Part 1: calculate mixing length scale
    !---------------------------------------------------------------------------------
    sqrttke = sqrt(max(0d0,tke_old))

    ! turbulent mixing length
    mxl = sqrt(2D0)*sqrttke/sqrt(max(1d-12,Nsqr))

    ! constrain mixing length scale as in MITgcm
    if (tke_mxl_choice==2) then
      !FIXME: What should we do at the surface and bottom?
      mxl(1) = 0.d0
      mxl(nlev+1) = 0.d0
      do k=2,nlev
        mxl(k) = min(mxl(k), mxl(k-1)+dzw(k-1))
      enddo
      mxl(nlev) = min(mxl(nlev), mxl_min+dzw(nlev))
      do k=nlev-1,2,-1
        mxl(k) = min(mxl(k), mxl(k+1)+dzw(k))
      enddo
      mxl= max(mxl,mxl_min)
    ! bounded by the distance to surface/bottom
    elseif (tke_mxl_choice==3) then
      depth = sum(dzw(1:nlev))
      do k=2,nlev+1
       zzw = sum(dzw(1:k-1))
       mxl(k) = min(zzw,mxl(k),depth-zzw)
      enddo
      mxl(1) = mxl(2)
      mxl= max(mxl,mxl_min)
    else
      write(*,*) 'Wrong choice of tke_mxl_choice. Aborting...'
      stop
    endif

    !---------------------------------------------------------------------------------
    ! Part 2: calculate diffusivities
    !---------------------------------------------------------------------------------
    ! see. Blanke and Delecluse 1993, eq. 2.25
    KappaM_out = min(KappaM_max,c_k*mxl*sqrttke)
    Rinum = Nsqr/max(Ssqr,1d-12)

    if (.not.only_tke) then  !IDEMIX is on
      Rinum = min(Rinum,KappaM_out*Nsqr/max(1d-12,alpha_c*E_iw**2))
    end if

    ! Richardson number dependent expression is used for the P_rt (Osborn 1980,
    ! Crawford 1982) --> see. Blanke and Delecluse 1993, eq. 2.33 & 2.34
    prandtl=max(1d0,min(10d0,6.6*Rinum))

    ! see. Blanke and Delecluse 1993, eq. 2.26
    KappaH_out=KappaM_out/prandtl

    !---------------------------------------------------------------------------------
    ! Part 3: tke forcing
    !---------------------------------------------------------------------------------
    ! initialize forcing
    forc = 0.0

    ! --- forcing by shear and buoycancy production
    K_diss_v   = Ssqr*KappaM_out
    P_diss_v   = Nsqr*KappaH_out
    P_diss_v(1) = -forc_rho_surf*grav/rho_ref
    forc = forc + K_diss_v - P_diss_v

    ! --- additional langmuir turbulence term
    if (l_lc) then
      forc = forc + tke_plc
    endif

    ! --- forcing by internal wave dissipation
    if (.not.only_tke) then
      forc = forc + iw_diss
    endif

    !---------------------------------------------------------------------------------
    ! Part 4: vertical diffusion and dissipation is solved implicitely
    !---------------------------------------------------------------------------------
    ke = 0.d0
    do k = 1, nlev
      kp1 = min(k+1,nlev)
      kk  = max(k,2)
      ke(k) = alpha_tke*0.5*(KappaM_out(kp1)+KappaM_out(kk))
    enddo

    !--- c is lower diagonal of matrix
    do k=1,nlev
      c_dif(k) = ke(k)/( dzt(k)*dzw(k) )
    enddo
    c_dif(nlev+1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

    !--- b is main diagonal of matrix
    do k=2,nlev
      b_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) ) + ke(k)/( dzt(k)*dzw(k) )
    enddo

    !--- a is upper diagonal of matrix
    do k=2,nlev+1
      a_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) )
    enddo
    a_dif(1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

    ! copy tke_old
    tke_upd(1:nlev+1) = tke_old(1:nlev+1)

    ! upper boundary condition
    if (use_ubound_dirichlet) then
      sqrttke(1)      = 0.d0 ! to suppres dissipation for k=1
      forc(1)         = 0.d0 ! to suppres forcing for k=1
      tke_surf        = max(tke_surf_min, cd*forc_tke_surf)
      tke_upd(1)      = tke_surf
      ! add diffusive part that depends on tke_surf to forcing
      diff_surf_forc  = a_dif(2)*tke_surf
      forc(2)         = forc(2)+diff_surf_forc
      a_dif(2)        = 0.d0 ! and set matrix element to zero
      b_dif(1)        = 0.d0 ! 0 line in matrix for k=1
      c_dif(1)        = 0.d0 ! 0 line in matrix for k=1
    else
      ! add wind forcing
      forc(1) = forc(1) + (cd*forc_tke_surf**(3./2.))/(dzt(1))
      b_dif(1)        = ke(1)/( dzt(1)*dzw(1) )
      diff_surf_forc  = 0.0
    endif

    ! lower boundary condition
    if (use_lbound_dirichlet) then
      sqrttke(nlev+1) = 0.d0 ! to suppres dissipation for k=nlev+1
      forc(nlev+1)    = 0.d0 ! to suppres forcing for k=nlev+1
      tke_bott        = tke_min
      tke_upd(nlev+1) = tke_bott
      ! add diffusive part that depends on tke_bott to forcing
      diff_bott_forc  = c_dif(nlev)*tke_bott
      forc(nlev)      = forc(nlev)+diff_bott_forc
      c_dif(nlev)     = 0.d0 ! and set matrix element to zero
      b_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
      a_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
    else
      b_dif(nlev+1)   = ke(nlev)/( dzt(nlev+1)*dzw(nlev) )
      diff_bott_forc  = 0.0
    endif

    !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
    a_tri = -dtime*a_dif
    b_tri = 1+dtime*b_dif
    b_tri(2:nlev) = b_tri(2:nlev) + dtime*c_eps*sqrttke(2:nlev)/mxl(2:nlev)
    c_tri = -dtime*c_dif

    !--- d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
    d_tri(1:nlev+1)  = tke_upd(1:nlev+1) + dtime*forc(1:nlev+1)

    ! solve the tri-diag matrix
    call solve_tridiag(a_tri, b_tri, c_tri, d_tri, tke_new, nlev+1)

    ! --- diagnose implicite tendencies (only for diagnostics)
    ! vertical diffusion of TKE
    do k=2,nlev
      tke_Tdif(k) = a_dif(k)*tke_new(k-1) - b_dif(k)*tke_new(k) + c_dif(k)*tke_new(k+1)
    enddo
    tke_Tdif(1) = - b_dif(1)*tke_new(1) + c_dif(1)*tke_new(2)
    tke_Tdif(nlev+1) = a_dif(nlev+1)*tke_new(nlev) - b_dif(nlev+1)*tke_new(nlev+1)
    tke_Tdif(2) = tke_Tdif(2) + diff_surf_forc
    tke_Tdif(nlev) = tke_Tdif(nlev) + diff_bott_forc

    if (use_ubound_dirichlet) then
      tke_Tdif(1) = - ke(1)/dzw(1)/dzt(1) &
                      * (tke_surf-tke_new(2))
    endif
    if (use_lbound_dirichlet) then
      k = nlev+1
      tke_Tdif(k) = ke(k-1)/dzw(k-1)/dzt(k) &
                      * (tke_new(k-1)-tke_bott)
    endif

    ! dissipation of TKE
    tke_Tdis = 0.d0
    tke_Tdis(2:nlev) = -c_eps/mxl(2:nlev)*sqrttke(2:nlev)*tke_new(2:nlev)

    !---------------------------------------------------------------------------------
    ! Part 5: reset tke to bounding values
    !---------------------------------------------------------------------------------
    ! copy of unrestored tke to diagnose energy input by restoring
    tke_unrest = tke_new

    ! restrict values of TKE to tke_min, if IDEMIX is not used
    if (only_tke) then
      tke_new(1:nlev+1) = MAX(tke_new(1:nlev+1), tke_min)
    end if

    !---------------------------------------------------------------------------------
    ! Part 6: Assign diagnostic variables
    !---------------------------------------------------------------------------------
    tke_Tbpr(1:nlev+1) = -P_diss_v(1:nlev+1)
    tke_Tspr(1:nlev+1) = K_diss_v(1:nlev+1)
    !tke_Tdif is set above
    tke_Tbck = (tke_new-tke_unrest)/dtime
    if (use_ubound_dirichlet) then
      tke_Twin(1) = (tke_new(1)-tke_old(1))/dtime - tke_Tdif(1)
      tke_Tbck(1) = 0.0
    else
      tke_Twin(1) = (cd*forc_tke_surf**(3./2.))/(dzt(1))
    endif
    if (use_lbound_dirichlet) then
      tke_Twin(nlev+1) = (tke_new(nlev+1)-tke_old(nlev+1))/dtime - tke_Tdif(nlev+1)
      tke_Tbck(nlev+1) = 0.0
    else
      tke_Twin(nlev+1) = 0.0
    endif

    tke_Tiwf(1:nlev+1) = iw_diss(1:nlev+1)
    tke_Ttot = (tke_new-tke_old)/dtime
    tke_Lmix(nlev+1:) = 0.0
    tke_Lmix(1:nlev+1) = mxl(1:nlev+1)
    tke_Pr(nlev+1:) = 0.0
    tke_Pr(1:nlev+1) = prandtl(1:nlev+1)

end subroutine integrate_tke

!=================================================================================
! TKE Av/Kv PRODUCER for one ALE timestep. 1:1 transcription of calc_cvmix_tke
! (gen_modules_cvmix_tke.F90:378-660), library-free + optional-partit. Per OWNED node:
! assemble vshear2/bvfreq2/dz_trr/normstress + tke_old, call integrate_tke (the prognostic
! column solve), zero the tke_Av/tke_Kv endpoints; then exchange_nod(tke_Kv);Kv=tke_Kv and
! exchange_nod(tke_Av) BEFORE the node->elem kite-weighted (elem_w) average into Av. tke is the prognostic
! field carried step->step in dynamics%work%tke (tke_old copy -> integrate_tke -> tke_new in
! the same slab). IDEMIX/Langmuir are off (tke_only=.true.) -> iw_diss/E_iw/alpha_c/tke_plc
! are passed as length-(nlev+1) ZERO arrays; forc_rho_surf=0. The 10 diagnostics are computed
! into discarded local scratch (NOT load-bearing for tke/Av/Kv; no diagnostics-output bloat).
subroutine calc_cvmix_tke(dynamics, stress_node_surf, dt, mesh, partit)
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_dyn,         only: t_dyn
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    use mod_constants,   only: density_0, g
    type(t_dyn),    intent(inout), target   :: dynamics
    real(kind=WP),  intent(in)              :: stress_node_surf(:,:)   ! (2, nnod) node surface stress
    real(kind=WP),  intent(in)              :: dt                      ! time step [s]
    type(t_mesh),   intent(in),    target   :: mesh
    type(t_partit), intent(in),    optional :: partit

    integer :: node, nz, nln, nun, nlev, elem, elnodes(MAX_NV)
    integer :: nNodO, nNodL, nEdgeO, nElemO
    real(kind=WP) :: forc_surf
    real(kind=WP), dimension(:,:),   pointer :: bvfreq, tke, tke_Av, tke_Kv, Kv, Av, shear2
    real(kind=WP), dimension(mesh%nl) :: vshear2, bvfreq2, dz_trr, tke_old_col, dzw_col, zero_col
    ! discarded diagnostic scratch (integrate_tke requires them; not load-bearing here)
    real(kind=WP), dimension(mesh%nl) :: dg_tbpr, dg_tspr, dg_tdif, dg_tdis, dg_twin, &
                                         dg_tiwf, dg_tbck, dg_ttot, dg_lmix, dg_pr

    shear2 => dynamics%work%shear2
    bvfreq => dynamics%work%bvfreq
    tke    => dynamics%work%tke
    tke_Av => dynamics%work%tke_Av
    tke_Kv => dynamics%work%tke_Kv
    Kv     => dynamics%work%Kv
    Av     => dynamics%work%Av
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

    zero_col = 0.0_WP

    do node = 1, nNodO
        nln  = mesh%nlevels_nod2D(node)-1
        nun  = mesh%ulevels_nod2D(node)
        nlev = nln-nun+1

        ! TKE surface momentum forcing — norm of the nodal surface wind stress
        forc_surf = sqrt( stress_node_surf(1,node)**2 + stress_node_surf(2,node)**2)/density_0

        ! 3D vertical velocity shear (zero the full nl column, then fill nun+1..nln)
        ! shared shear (oce_vert_spline::compute_shear2). NOTE this differs from the
        ! previous TKE-local form in the LAST BITS even at shear_splines=.false.: TKE
        ! divided by (Z(nz-1)-Z(nz))**2 once, the shared producer multiplies by dz_inv
        ! twice (the PP / ri_iwmix operand order). Same value, different rounding.
        vshear2=0.0_WP
        vshear2(nun+1:nln) = shear2(nun+1:nln, node)

        ! square of Brunt-Vaisala frequency (bvfreq already holds N^2)
        bvfreq2        = 0.0_WP
        bvfreq2(nun+1:nln) = bvfreq(nun+1:nln,node)

        ! dz_trr distance between tracer points; surface/bottom = half the layer thickness
        dz_trr            = 0.0_WP
        dz_trr(nun+1:nln) = abs(mesh%Z_3d_n(nun:nln-1,node)-mesh%Z_3d_n(nun+1:nln,node))
        dz_trr(nun)       = mesh%hnode(nun,node)/2.0_WP
        dz_trr(nln+1)     = mesh%hnode(nln,node)/2.0_WP

        ! prognostic recurrence: tke_old = the prior-step tke slab; dzw=hnode(nun:nln)
        tke_old_col      = tke(:,node)
        dzw_col(1:nlev)  = mesh%hnode(nun:nln,node)

        call integrate_tke( &
             tke_old      = tke_old_col(nun:nln+1),  &
             tke_new      = tke(   nun:nln+1,node),  & ! out --> turbulent kinetic energy
             KappaM_out   = tke_Av(nun:nln+1,node),  & ! out
             KappaH_out   = tke_Kv(nun:nln+1,node),  & ! out
             dzw          = dzw_col(1:nlev),         &
             dzt          = dz_trr(nun:nln+1),       &
             nlev         = nlev,                    &
             Ssqr         = vshear2(nun:nln+1),      &
             Nsqr         = bvfreq2(nun:nln+1),      &
             tke_Tbpr     = dg_tbpr(1:nlev+1),       &
             tke_Tspr     = dg_tspr(1:nlev+1),       &
             tke_Tdif     = dg_tdif(1:nlev+1),       &
             tke_Tdis     = dg_tdis(1:nlev+1),       &
             tke_Twin     = dg_twin(1:nlev+1),       &
             tke_Tiwf     = dg_tiwf(1:nlev+1),       &
             tke_Tbck     = dg_tbck(1:nlev+1),       &
             tke_Ttot     = dg_ttot(1:nlev+1),       &
             tke_Lmix     = dg_lmix(1:nlev+1),       &
             tke_Pr       = dg_pr(1:nlev+1),         &
             tke_plc      = zero_col(1:nlev+1),      &
             forc_tke_surf= forc_surf,               &
             E_iw         = zero_col(1:nlev+1),      &
             dtime        = dt,                      &
             iw_diss      = zero_col(1:nlev+1),      &
             forc_rho_surf= 0.0_WP,                  &
             rho_ref      = density_0,               &
             grav         = g,                       &
             alpha_c      = zero_col(1:nlev+1))

        tke_Av(nln+1,node)=0.0_WP
        tke_Kv(nln+1,node)=0.0_WP
        tke_Av(nun  ,node)=0.0_WP
        tke_Kv(nun  ,node)=0.0_WP
    end do !--> do node = 1,nNodO

    !___________________________________________________________________________
    ! write out diffusivity (nodes) — Kv consumer UNCHANGED
    if (is_multirank(partit)) call exchange_nod(tke_Kv, partit)
    Kv = tke_Kv

    !___________________________________________________________________________
    ! write out viscosity (elements) — exchange tke_Av BEFORE the node->elem average
    if (is_multirank(partit)) call exchange_nod(tke_Av, partit)
    Av = 0.0_WP
    do elem=1, nElemO
        elnodes=mesh%elem2D_nodes(:,elem)
        do nz=mesh%ulevels(elem)+1,mesh%nlevels(elem)-1
            Av(nz,elem) = sum(mesh%elem_w(:,elem)*tke_Av(nz,elnodes))    ! (elementwise, kite-weighted)
        end do
    end do
end subroutine calc_cvmix_tke

!=================================================================================
! Thomas tridiagonal solve. VERBATIM from cvmix_utils_addon.F90:116-149 — row 1 is a
! DIRECT divide (cp(1)=c(1)/b(1), dp(1)=d(1)/b(1)); the elimination rows use the
! reciprocal-multiply form (fxa=1D0/m; cp(i)=c(i)*fxa) — NOT c(i)/m (the C-port
! last-bit lesson). Match the operand/operation order line-for-line.
subroutine solve_tridiag(a,b,c,d,x,n)
    implicit none
    integer, intent(in)                       :: n
    real(kind=WP), dimension(n), intent(in)   :: a,b,c,d
    real(kind=WP), dimension(n), intent(out)  :: x
    real(kind=WP), dimension(n) :: cp,dp
    real(kind=WP) :: m,fxa
    integer i

    ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do i = 2,n
      m = b(i)-cp(i-1)*a(i)
      fxa = 1D0/m
      cp(i) = c(i)*fxa
      dp(i) = (d(i)-dp(i-1)*a(i))*fxa
    enddo
    ! initialize x
    x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
    end do
end subroutine solve_tridiag

end module oce_mixing_tke
