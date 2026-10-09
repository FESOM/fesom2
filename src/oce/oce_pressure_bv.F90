module oce_pressure_bv
    ! Hydrostatic pressure + Brunt-Vaisala frequency in a single pass, using the
    ! SPLIT form of the Jackett-McDougall equation of state. Transcribed from
    ! FESOM2 v2.7.3 oce_ale_pressure_bv.F90:
    !   pressure_bv          (194-501)  density_m_rho0 + hpressure + bvfreq
    !   densityJM_components (2605-2669) split-form secant bulk modulus + rhopot
    !   smooth_nod3D         (gen_support.F90:99-198) horizontal N^2 smoother
    !
    ! M2.1 SCOPE — the first dynamics kernel; gate target = density_m_rho0 (PGF
    ! density anomaly), hpressure (top-down hydrostatic integration) and bvfreq
    ! (N^2, squared). OMITTED (pure diagnostics that do NOT feed the gated fields):
    ! the dMOC density (ldiag_dMOC), the KPP buoyancy diff dbsfc1/dbsfc, db_max,
    ! and the mixed-layer depths MLD1/2/3. Dropping them leaves density_m_rho0 /
    ! hpressure / bvfreq bit-for-bit unchanged (density_m_rho0 is computed BEFORE
    ! rho_surf/dbsfc1/db_max; the MLD logic only READS bvfreq/rhopot/Z_3d_n).
    !
    ! EOS NOTES (bit-identity critical):
    !  - NEVER linearize alpha/beta — the split form factorizes the secant bulk
    !    modulus (bulk_0 + pz*(bulk_pz + pz*bulk_pz2)) and the bits depend on that
    !    factorization. Transcribe densityJM_components verbatim.
    !  - density anomaly subtracts the density_ref(nz,node) ARRAY (here = density_0
    !    everywhere, since use_density_ref=.false. on pi; init_ref_density is the
    !    deferred use_density_ref=.true. path). N^2 divides by the SCALAR density_0.
    !  - real(state_equation): 1 (full EOS, case 1). Kept a variable so the
    !    0.1*Z*real(state_equation) terms transcribe verbatim (state_equation=0
    !    linear EOS would zero the pressure term — not gated here).
    !
    ! which_ALE='linfs' (pi) routes the hydrostatic pressure integration through
    ! THIS routine (the zstar/zlevel cases compute it in pressure_force_4_zxxxx,
    ! M2.x). use_cavity / partial cells (nzmin>1, the cavity branches below) are
    ! v1-dropped and UNGATED on pi (nzmin==1, full cells) — transcribed faithfully
    ! for the M2.11 cavity re-gate.
    !
    ! 1-rank only: smooth_nod's per-cycle exchange_nod(bvfreq) is a no-op here and
    ! is dropped (lifted with the other M1/M2 halo exchanges at M2.12). The three
    ! outputs are intent(inout) and NOT zeroed inside (FESOM2 pressure_bv does not
    ! initialize the below-bottom entries either) — the CALLER must zero them so the
    ! dump's below-bottom region is a deterministic 0 on both sides of the gate.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_constants,   only: density_0, g, pi
    use mod_config,      only: which_ALE
    use oce_vert_spline, only: spline_eval
    use mod_param_phys,  only: state_equation, N2smth_h, N2smth_hidx, N2_splines, &
                               scaling_ODM95, ODM95_Scr, ODM95_Sd, scaling_LDD97, &
                               LDD97_c, LDD97_rmin, LDD97_rmax, Fer_GM, Redi, Redi_Ktaper
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    implicit none
    private
    public :: pressure_bv, densityJM_components, insitu2pot
    public :: sw_alpha_beta, compute_sigma_xy, compute_neutral_slope
    public :: smooth_nod   ! area-weighted patch smoother (= FESOM2 smooth_nod3D); KPP smooth_blmc

contains

    !===========================================================================
    subroutine pressure_bv(temp, salt, density_ref, mesh, density_m_rho0, hpressure, bvfreq, partit, dbsfc, bvfreq_raw)
        ! temp/salt/density_ref: (nl-1, nod2D) inputs. density_m_rho0: (nl-1, nod2D)
        ! out. hpressure/bvfreq: (nl, nod2D) out (only 1..nzmax used). The caller
        ! pre-zeros the three outputs (see header).
        ! M2.12c: optional partit -> owned+halo node loop (FESOM2 :237-244 do node=1,
        ! myDim_nod2D+eDim_nod2D — the EOS is per-node so computing the halo here saves
        ! an exchange; the downstream pgf/smoothing then read it locally) + smooth_nod's
        ! per-sweep exchange_nod(bvfreq).
        ! M5a-3 (KPP): optional dbsfc (nl, nod2D) buoyancy-difference-wrt-surface output —
        ! the bulk-Ri numerator bldepth consumes. Only the KPP gate/driver passes it; the
        ! lifecycle/step callers omit it (a separate NOVECTOR loop, byte-neutral when absent).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: temp(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: salt(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: density_ref(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: density_m_rho0(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: hpressure(mesh%nl, mesh%nod2D)
        real(kind=WP), intent(inout) :: bvfreq(mesh%nl, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        real(kind=WP), intent(inout), optional :: dbsfc(mesh%nl, mesh%nod2D)
        ! M?: optional copy of N^2 BEFORE the horizontal smoother, for consumers whose
        ! answer depends on the true local sign (mo_convect). Absent -> nothing written.
        real(kind=WP), intent(inout), optional :: bvfreq_raw(mesh%nl, mesh%nod2D)

        integer       :: node, nz, nzmax, nzmin
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: zmean, dz_inv, a, rho_up, rho_dn, t, s, smin
        real(kind=WP) :: bulk_up, bulk_dn, rho_surf
        real(kind=WP) :: rhopot(mesh%nl), bulk_0(mesh%nl), bulk_pz(mesh%nl)
        real(kind=WP) :: bulk_pz2(mesh%nl), rho(mesh%nl), dbsfc1(mesh%nl)
        ! N2_splines scratch (one column)
        integer       :: n_sp, k_sp
        real(kind=WP) :: alpha_sp, beta_sp
        real(kind=WP) :: xc_sp(mesh%nl), xi_sp(mesh%nl), tc_sp(mesh%nl), sc_sp(mesh%nl)
        real(kind=WP) :: ti_sp(mesh%nl), si_sp(mesh%nl), dt_sp(mesh%nl), ds_sp(mesh%nl)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !_______________________________________________________________________
        ! Screen salinity (diagnostic; s<0 would break sqrt(s) in the EOS).
        smin = 0.0_WP
        do node=1, nNodL
            do nz = mesh%ulevels_nod2D(node), mesh%nlevels_nod2D(node)-1
                smin = min(smin, salt(nz,node))
            end do
        end do
        if (smin < 0.0_WP) write(*,*) ' --> oce_pressure_bv: s<0 happens!', smin

        !_______________________________________________________________________
        do node=1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)

            rho      = 0.0_WP
            bulk_0   = 0.0_WP
            bulk_pz  = 0.0_WP
            bulk_pz2 = 0.0_WP
            rhopot   = 0.0_WP

            !___________________________________________________________________
            ! apply equation of state
            do nz=nzmin, nzmax-1
                t=temp(nz, node)
                s=salt(nz, node)
                call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
            end do

            !___________________________________________________________________
            ! compute density for PGF
            do nz=nzmin, nzmax-1
                rho(nz) = bulk_0(nz) + mesh%Z_3d_n(nz,node)*(bulk_pz(nz) + mesh%Z_3d_n(nz,node)*bulk_pz2(nz))
                rho(nz) = rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*mesh%Z_3d_n(nz,node)*real(state_equation,WP))-density_ref(nz,node)
                density_m_rho0(nz,node) = rho(nz)
            end do

            !___________________________________________________________________
            ! M5a-3 (KPP): buoyancy difference wrt the surface (dbsfc, m/s2) — the bulk-Ri
            ! numerator bldepth consumes. Bring the surface water adiabatically to depth nz
            ! (rho_surf = surface bulk-modulus components evaluated at Z_3d_n(nz)), then
            ! dbsfc1 = -g*(rho_surf - rho_full(nz))/rho_full(nz). Faithful to FESOM2
            ! oce_ale_pressure_bv.F90:326-339 (filled there under mixing_kpp). !DIR$ NOVECTOR
            ! keeps the divide scalar (L29): the oracle's dbsfc1 sits in the db_max-reduction
            ! density loop (scalar); operands are the byte-matched density_m_rho0+density_ref,
            ! so scalar-vs-scalar is exact. dbsfc1(nzmin)=0 (rho_surf==rho_full there).
            if (present(dbsfc)) then
                !DIR$ NOVECTOR
                do nz=nzmin, nzmax-1
                    rho_surf = bulk_0(nzmin) + mesh%Z_3d_n(nz,node)*(bulk_pz(nzmin) + mesh%Z_3d_n(nz,node)*bulk_pz2(nzmin))
                    rho_surf = rho_surf*rhopot(nzmin)/(rho_surf+0.1_WP*mesh%Z_3d_n(nz,node)*real(state_equation,WP))
                    dbsfc1(nz) = -g*(rho_surf - (density_m_rho0(nz,node)+density_ref(nz,node))) / (density_m_rho0(nz,node)+density_ref(nz,node))
                end do
                dbsfc1(nzmax) = dbsfc1(nzmax-1)
                dbsfc(nzmin:nzmax, node) = dbsfc1(nzmin:nzmax)
            end if

            !___________________________________________________________________
            ! fill density levels occupied by the cavity (nzmin>1; v1-dropped, ungated
            ! on pi). Mass that corresponds to T/S at the cavity-ocean interface.
            if (nzmin>1) then
                t=temp(nzmin, node)
                s=salt(nzmin, node)
                do nz=1, nzmin-1
                    call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
                    rho(nz)= bulk_0(nz)   + mesh%Z_3d_n(nz,node)*(bulk_pz(nz)   + mesh%Z_3d_n(nz,node)*bulk_pz2(nz))
                    rho(nz)=rho(nz)*rhopot(nz)/(rho(nz)+0.1_WP*mesh%Z_3d_n(nz,node)*real(state_equation,WP))-density_ref(nz,node)
                    density_m_rho0(nz,node) = rho(nz)
                end do
            end if

            !___________________________________________________________________
            ! calculate pressure (linfs or cavity case)
            if (trim(which_ALE)=='linfs') then
                if (nzmin>1) then ! cavity case (v1-dropped, ungated on pi)
                    hpressure(nzmin, node)=0.5_WP*(mesh%zbar_3d_n(1,node)-mesh%zbar_3d_n(2,node))*rho(1)*g
                    do nz=2,nzmin
                        a=0.5_WP*g*(rho(nz-1)*(mesh%zbar_3d_n(nz-1,node)-mesh%zbar_3d_n(nz,node))+rho(nz)*(mesh%zbar_3d_n(nz,node)-mesh%zbar_3d_n(nz+1,node)))
                        hpressure(nzmin, node)=hpressure(nzmin, node)+a
                    end do
                else
                    hpressure(nzmin, node)=-mesh%Z_3d_n(nzmin,node)*rho(nzmin)*g
                end if

                ! pressure below surface boundary: integrate g*rho*dz over half the
                ! previous + half the actual layer thickness to mid-depth of the layer
                do nz=nzmin+1,nzmax-1
                    a=0.5_WP*g*(rho(nz-1)*mesh%hnode(nz-1,node)+rho(nz)*mesh%hnode(nz,node))
                    hpressure(nz, node)=hpressure(nz-1, node)+a
                end do
            end if

            !___________________________________________________________________
            ! squared Brunt-Vaisala frequency N^2 (defined on full levels except the
            ! first and the last). N^2>0 stable, N^2<0 unstable stratification.
            do nz=nzmin+1,nzmax-1
                zmean   = 0.5_WP*sum(mesh%Z_3d_n(nz-1:nz, node))
                bulk_up = bulk_0(nz-1) + zmean*(bulk_pz(nz-1) + zmean*bulk_pz2(nz-1))
                bulk_dn = bulk_0(nz)   + zmean*(bulk_pz(nz)   + zmean*bulk_pz2(nz))
                rho_up  = bulk_up*rhopot(nz-1) / (bulk_up + 0.1_WP*zmean*real(state_equation,WP))
                rho_dn  = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1_WP*zmean*real(state_equation,WP))
                dz_inv  = 1.0_WP/(mesh%Z_3d_n(nz-1,node)-mesh%Z_3d_n(nz,node))
                bvfreq(nz,node)  = -g*dz_inv*(rho_up-rho_dn)/density_0
            end do

            !___________________________________________________________________
            ! N2_splines: overwrite the interior levels with the spline form
            !     N2 = g*(beta*dS/dx - alpha*dT/dx),   x = depth (increasing downward)
            ! evaluated AT the interface zbar (the two-point form above is centred on the
            ! midpoint of the two layer centres, which on this stretched grid sits up to
            ! ~17% of the local dz away from the interface). There is no single rho(z)
            ! column to spline -- each interface re-references both parcels to its own
            ! zmean -- so the derivative is taken through the chain rule with alpha/beta at
            ! the interface state. Sign check: T decreasing with depth -> dT/dx<0 -> N2>0.
            ! Same natural cubic spline as the shear (oce_vert_spline), for like-by-like
            ! ratios downstream. Columns with a single layer have no interior level.
            if (N2_splines) then
                n_sp = nzmax - nzmin
                if (n_sp >= 2) then
                    do k_sp = 1, n_sp
                        xc_sp(k_sp) = -mesh%Z_3d_n(nzmin+k_sp-1, node)
                        tc_sp(k_sp) =  temp(nzmin+k_sp-1, node)
                        sc_sp(k_sp) =  salt(nzmin+k_sp-1, node)
                    end do
                    do k_sp = 1, n_sp-1
                        xi_sp(k_sp) = -mesh%zbar_3d_n(nzmin+k_sp, node)
                    end do
                    call spline_eval(n_sp, xc_sp(1:n_sp), tc_sp(1:n_sp), xi_sp(1:n_sp-1), ti_sp(1:n_sp-1), dt_sp(1:n_sp-1))
                    call spline_eval(n_sp, xc_sp(1:n_sp), sc_sp(1:n_sp), xi_sp(1:n_sp-1), si_sp(1:n_sp-1), ds_sp(1:n_sp-1))
                    do k_sp = 1, n_sp-1
                        call sw_alpha_beta_point(ti_sp(k_sp), si_sp(k_sp), xi_sp(k_sp), alpha_sp, beta_sp)
                        bvfreq(nzmin+k_sp, node) = g*(beta_sp*ds_sp(k_sp) - alpha_sp*dt_sp(k_sp))
                    end do
                end if
            end if

            bvfreq(nzmin,node)=bvfreq(nzmin+1,node)
            bvfreq(nzmax,node)=bvfreq(nzmax-1,node)
        end do

        !_______________________________________________________________________
        ! hand out the UNSMOOTHED N^2 first -- smooth_nod overwrites bvfreq in place, and
        ! its 1/3-own + 2/3-neighbour stencil can flip the SIGN of N^2 near convection,
        ! which is exactly what a static-instability test must not see.
        ! BOUNDED BY nNodL, NOT a whole-array assignment: the dummy is declared with the
        ! GLOBAL mesh%nod2D while the MR drivers allocate (nl, nNodL) with nNodL =
        ! myDim+eDim. The over-declared dummy is only a view, safe as long as every access
        ! stays inside nNodL like the rest of this routine -- `bvfreq_raw = bvfreq` did not,
        ! and wrote nl*nod2D elements into an nl*nNodL array (heap corruption at 512 ranks,
        ! invisible at np=1 where nNodL == nod2D).
        if (present(bvfreq_raw)) then
            do node = 1, nNodL
                bvfreq_raw(:, node) = bvfreq(:, node)
            end do
        end if

        ! apply horizontal smoothing of the N^2 buoyancy frequency (GM/Redi + the
        ! mixing schemes consume this smoothed field)
        if (N2smth_h) call smooth_nod(bvfreq, N2smth_hidx, mesh, partit)

    end subroutine pressure_bv

    !===========================================================================
    subroutine densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
        ! Split form of the Jackett-McDougall (1992, CSIRO) in-situ density EOS.
        ! Verbatim from FESOM2 oce_ale_pressure_bv.F90:2605-2669 (N. Rakowski split).
        ! rho_in-situ(pz) = bulk*rhopot/(bulk+0.1*pz) with bulk = bulk_0 + pz*(bulk_pz
        ! + pz*bulk_pz2). DO NOT pre-multiply/refactor — the bit pattern depends on
        ! the exact Horner factorization below.
        real(kind=WP), intent(in)  :: t, s
        real(kind=WP), intent(out) :: bulk_0, bulk_pz, bulk_pz2, rhopot
        real(kind=WP) :: s_sqrt

        real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
        real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
        real(kind=WP), parameter   :: at4   = -1.361629e-5
        real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
        real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4
        real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729
        real(kind=WP), parameter   :: asst2 = -1.909078e-2
        real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
        real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7
        real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4
        real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
        real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10
        real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
        real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9
        real(kind=WP), parameter   :: ap2st2= -3.508914e-9

        real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
        real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
        real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
        real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
        real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7
        real(kind=WP), parameter   :: bst4 = 5.38750e-9
        real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
        real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4

        s_sqrt = sqrt(s)

        bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
                + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
                     + s_sqrt*(ass  + t*(asst + t*asst2)))

        bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                        + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

        bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)                           &
                      + s *(ap2s + t*(ap2st + t*ap2st2))

        rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))  &
                     + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                        + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                             + s* bss2)
    end subroutine densityJM_components

    !===========================================================================
    subroutine smooth_nod(arr, N_smooth, mesh, partit)
        ! Mass-matrix horizontal smoother, transcribed from FESOM2 g_support
        ! smooth_nod3D (gen_support.F90:99-198): applies the lumped P1 mass matrix
        ! N_smooth times. Hybrid meshes: node value = the kite-share-weighted mean
        ! (elem_area*nod_in_elem2D_w) of the kite-weighted element means (elem_w) of the
        ! wet elements around the node -- FESOM2's sum(A*(a1+a2+a3))/(3*sum(A)) on
        ! triangles (weights 1/3), equal to roundoff. M2.12c: optional partit ->
        ! owned-node loops + the per-sweep exchange_nod(arr) (the patch accumulation
        ! reads arr at an owned node's owned-element nodes, which can be halo nodes, so
        ! arr must be halo-valid each sweep).
        ! Per-level patch areas vary with the bathymetry (a deep node's shallower
        ! neighbour elements drop out level by level). elem_area is geom-proven; the
        ! nod_in_elem2D / elem2D_nodes order is area-gate-proven (L9).
        real(kind=WP), intent(inout)       :: arr(:,:)
        integer,       intent(in)          :: N_smooth
        type(t_mesh),  intent(in), target  :: mesh
        type(t_partit), intent(in), optional :: partit
        integer :: n, q, el, nz, j, nlev, uln, nln, ule, nle
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: wa
        real(kind=WP), allocatable :: vol(:,:), work_array(:,:)

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        nlev=ubound(arr,1)
        allocate(vol(mesh%nl, nNodL), work_array(nlev, nNodL))

        ! First sweep: precompute the (inverse) patch areas, then smooth.
        do n=1, nNodO
            uln = mesh%ulevels_nod2d(n)
            nln = min(nlev, mesh%nlevels_nod2d(n))
            vol(       1:nln,n) = 0._WP
            work_array(1:nln,n) = 0._WP
            do j=1, mesh%nod_in_elem2D_num(n)
                el  = mesh%nod_in_elem2D(j,n)
                wa  = mesh%elem_area(el)*mesh%nod_in_elem2D_w(j,n)
                ule = max( uln, mesh%ulevels(el) )
                nle = min( nln, min(nlev, mesh%nlevels(el)) )
                do nz=ule, nle
                    vol(nz,n) = vol(nz,n) + wa
                    work_array(nz,n) = work_array(nz,n) + wa * elem_mean(nz, el)
                end do
            end do
            do nz=uln,nln
                vol(nz,n) = 1._WP / vol(nz,n)  ! inverse
            end do
        end do
        do n=1, nNodO
            uln = mesh%ulevels_nod2d(n)
            nln = min(nlev, mesh%nlevels_nod2d(n))
            do nz=uln,nln
                arr(nz, n) = work_array(nz, n) * vol(nz,n)
            end do
        end do
        if (is_multirank(partit)) call exchange_nod(arr, partit)

        ! Remaining sweeps reuse the precomputed inverse patch areas.
        do q=1,N_smooth-1
            do n=1, nNodO
                uln = mesh%ulevels_nod2d(n)
                nln = min(nlev, mesh%nlevels_nod2d(n))
                work_array(1:nln,n) = 0._WP
                do j=1, mesh%nod_in_elem2D_num(n)
                    el  = mesh%nod_in_elem2D(j,n)
                    wa  = mesh%elem_area(el)*mesh%nod_in_elem2D_w(j,n)
                    ule = max( uln, mesh%ulevels(el) )
                    nle = min( nln, min(nlev, mesh%nlevels(el)) )
                    do nz=ule, nle
                        work_array(nz,n) = work_array(nz,n) + wa * elem_mean(nz, el)
                    end do
                end do
            end do
            do n=1, nNodO
                uln = mesh%ulevels_nod2d(n)
                nln = min(nlev, mesh%nlevels_nod2d(n))
                do nz=uln,nln
                    arr(nz, n) = work_array(nz, n) * vol(nz,n)
                end do
            end do
            if (is_multirank(partit)) call exchange_nod(arr, partit)
        end do

        deallocate(vol, work_array)

    contains
        pure real(kind=WP) function elem_mean(nz, el)
            ! kite-weighted element mean of arr at level nz (4 slots; a triangle's slot 4
            ! has weight 0)
            integer, intent(in) :: nz, el
            elem_mean = mesh%elem_w(1,el)*arr(nz, mesh%elem2D_nodes(1,el)) &
                      + mesh%elem_w(2,el)*arr(nz, mesh%elem2D_nodes(2,el)) &
                      + mesh%elem_w(3,el)*arr(nz, mesh%elem2D_nodes(3,el)) &
                      + mesh%elem_w(4,el)*arr(nz, mesh%elem2D_nodes(4,el))
        end function elem_mean
    end subroutine smooth_nod

    !===========================================================================
    ! Convert in-situ temperature -> potential temperature, IN PLACE (temp), using
    ! salinity (read-only) at reference pressure pr=0. Transcribed from FESOM2
    ! oce_ale_pressure_bv.F90 insitu2pot (3074-3118): per-node loop, the in-situ T
    ! at each level is replaced by ptheta(S, T, |Z(nz)|, 0). FESOM2 uses the 1-D
    ! mid-depth Z(nz) for the pressure proxy pp=abs(Z(nz)) (NOT Z_3d_n — the comment
    ! at :3108 keeps Z for partial-cell stability at init). The IC driver feeds
    ! tracers%data(1)%values (in-situ T) / data(2)%values (S) here (M2.11b do_ic3d,
    ! t_insitu=.true.). M2.12-MVP: optional partit -> the FESOM2 owned+halo node loop
    ! (n=1..myDim+eDim); per-node independent so the halo potential T is computed locally
    ! (= the owner's value, no exchange needed). Absent partit -> n=1..nod2D (1-rank).
    subroutine insitu2pot(temp, salt, mesh, partit)
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(inout) :: temp(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: salt(mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, nzmin, nzmax, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: pp, pr, tt, ss
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        pr = 0.0_WP
        do n = 1, nNodL
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax-1
                tt = temp(nz, n)
                ss = salt(nz, n)
                pp = abs(mesh%Z(nz))
                temp(nz, n) = ptheta(ss, tt, pp, pr)
            end do
        end do
    end subroutine insitu2pot

    !===========================================================================
    ! Local potential temperature at reference pressure pr, via Bryden-1973
    ! adiabatic lapse rate + 4th-order Runge-Kutta. Verbatim from FESOM2
    ! oce_ale_pressure_bv.F90 ptheta (2674-2714). Args are mutated locally (t,p);
    ! the caller passes scalar temporaries so the mutation is harmless. checkvalue:
    ! theta = 36.89073 C for s=40, t=40, p=10000, pr=0.
    function ptheta(s, t, p, pr) result(theta)
        real(kind=WP) :: theta
        real(kind=WP) :: s, t, p, pr
        real(kind=WP) :: h, xk, q
        h  = pr - p
        xk = h*atg(s, t, p)
        t  = t + 0.5_WP*xk
        q  = xk
        p  = p + 0.5_WP*h
        xk = h*atg(s, t, p)
        t  = t + 0.29289322_WP*(xk-q)
        q  = 0.58578644_WP*xk + 0.121320344_WP*q
        xk = h*atg(s, t, p)
        t  = t + 1.707106781_WP*(xk-q)
        q  = 3.414213562_WP*xk - 4.121320344_WP*q
        p  = p + 0.5_WP*h
        xk = h*atg(s, t, p)
        theta = t + (xk-2.0_WP*q)/6.0_WP
    end function ptheta

    !===========================================================================
    ! Adiabatic temperature gradient deg C / decibar (Bryden 1973). Verbatim from
    ! FESOM2 oce_ale_pressure_bv.F90 atg (2719-2746). checkvalue: atg=3.255976e-4
    ! C/dbar for s=40, t=40, p=10000.
    function atg(s, t, p) result(g_atg)
        real(kind=WP) :: g_atg
        real(kind=WP) :: s, t, p, ds
        ds = s - 35.0_WP
        g_atg = (((-2.1687e-16_WP*t+1.8676e-14_WP)*t-4.6206e-13_WP)*p   &
              +((2.7759e-12_WP*t-1.1351e-10_WP)*ds+((-5.4481e-14_WP*t        &
              +8.733e-12_WP)*t-6.7795e-10_WP)*t+1.8741e-8_WP))*p             &
              +(-4.2393e-8_WP*t+1.8932e-6_WP)*ds                          &
              +((6.6228e-10_WP*t-6.836e-8_WP)*t+8.5258e-6_WP)*t+3.5803e-5_WP
    end function atg

    !===========================================================================
    subroutine sw_alpha_beta(temp, salt, mesh, sw_alpha, sw_beta, partit)
        ! Thermal expansion (sw_alpha) + saline contraction (sw_beta) coefficients,
        ! McDougall 1987 (FESOM2 oce_ale_pressure_bv.F90:2751). All-scalar per node/level
        ! -> NO SIMD-divide trap. t1 = T*1.00024 (ITS-90->IPTS-68), p1 = |Z_3d_n| [db],
        ! s35 = S-35. Per-node EOS -> loop owned+halo (nNodL), no exchange (FESOM3 idiom,
        ! same as pressure_bv). M4 GM/Redi producer; outputs feed compute_sigma_xy.
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: temp(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: salt(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: sw_alpha(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: sw_beta(mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, nzmin, nzmax, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: t1, t1_2, t1_3, t1_4, p1, p1_2, p1_3, s1, s35, s35_2, a_over_b
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do n = 1, nNodL
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n)
            do nz = nzmin, nzmax-1
                call sw_alpha_beta_point(temp(nz,n), salt(nz,n), abs(mesh%Z_3d_n(nz,n)), &
                                         sw_alpha(nz,n), sw_beta(nz,n))
            end do
        end do
    end subroutine sw_alpha_beta

    !===========================================================================
    pure subroutine sw_alpha_beta_point(t, s, p, alpha, beta)
        ! One (T, S, p) -> (alpha, beta). McDougall 1987, the polynomial sw_alpha_beta
        ! used inline before it was factored out here; the N2_splines path calls it at the
        ! INTERFACE state (spline-interpolated T,S, p = interface depth) so that the vertical
        ! density gradient is linearised exactly the way compute_sigma_xy linearises the
        ! horizontal one -- the neutral slope then divides like by like.
        ! t in ITS-90 degC (the 1.00024 IPTS-68 factor is applied here), s in psu, p in db
        ! (= depth in m, as the callers pass |Z|).
        real(kind=WP), intent(in)  :: t, s, p
        real(kind=WP), intent(out) :: alpha, beta
        real(kind=WP) :: t1, t1_2, t1_3, t1_4, p1, p1_2, p1_3, s35, s35_2, a_over_b
        t1 = t*1.00024_WP
        p1 = p
        t1_2 = t1*t1; t1_3 = t1_2*t1; t1_4 = t1_3*t1
        p1_2 = p1*p1; p1_3 = p1_2*p1
        s35  = s-35.0_WP; s35_2 = s35*s35
        beta = 0.785567e-3_WP - 0.301985e-5_WP*t1 &
             + 0.555579e-7_WP*t1_2 - 0.415613e-9_WP*t1_3 &
             + s35*(-0.356603e-6_WP + 0.788212e-8_WP*t1 &
             + 0.408195e-10_WP*p1 - 0.602281e-15_WP*p1_2) &
             + s35_2*(0.515032e-8_WP) &
             + p1*(-0.121555e-7_WP + 0.192867e-9_WP*t1 - 0.213127e-11_WP*t1_2) &
             + p1_2*(0.176621e-12_WP - 0.175379e-14_WP*t1) &
             + p1_3*(0.121551e-17_WP)
        a_over_b = 0.665157e-1_WP + 0.170907e-1_WP*t1 &
             - 0.203814e-3_WP*t1_2 + 0.298357e-5_WP*t1_3 &
             - 0.255019e-7_WP*t1_4 &
             + s35*(0.378110e-2_WP - 0.846960e-4_WP*t1 &
             - 0.164759e-6_WP*p1 - 0.251520e-11_WP*p1_2) &
             + s35_2*(-0.678662e-5_WP) &
             + p1*(0.380374e-4_WP - 0.933746e-6_WP*t1 + 0.791325e-8_WP*t1_2) &
             + p1_2*t1_2*(0.512857e-12_WP) &
             - p1_3*(0.302285e-13_WP)
        alpha = a_over_b*beta
    end subroutine sw_alpha_beta_point

    !===========================================================================
    subroutine compute_sigma_xy(temp, salt, sw_alpha, sw_beta, mesh, sigma_xy, partit)
        ! Density gradient sigma_xy (FESOM2 oce_ale_pressure_bv.F90:2851). Per owned node:
        ! accumulate element T/S gradients (gradient_sca, 4 slots) x the node's kite share
        ! elem_area*nod_in_elem2D_w over nod_in_elem2D (vol = the sum of the wet shares; on
        ! triangles FESOM2's elem_area weighting, both scaled by 1/3), then
        ! sigma_xy(1/2) = (-sw_alpha*Tgrad + sw_beta*Sgrad)/vol*density_0. The /vol is
        ! an array divide -> SIMD-divide trap (L29): gate then NOVECTOR if codegen vectorizes.
        ! The owned-node element neighbourhood is complete (M2.12a invariant) -> loop owned +
        ! exchange (oracle uses the 3x MPI_BARRIER aux dance; the rank-3 exchange_nod is the
        ! byte-equivalent owner->halo broadcast).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: temp(mesh%nl-1, mesh%nod2D), salt(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: sw_alpha(mesh%nl-1, mesh%nod2D), sw_beta(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: sigma_xy(2, mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, k, el, nln, uln, nle, ule, nNodO, nNodL, nEdgeO, nElemO, en(MAX_NV)
        real(kind=WP) :: wa
        real(kind=WP) :: tx(mesh%nl-1), ty(mesh%nl-1), sx(mesh%nl-1), sy(mesh%nl-1), vol(mesh%nl-1)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do n = 1, nNodO
            nln = mesh%nlevels_nod2D(n)-1
            uln = mesh%ulevels_nod2D(n)
            vol(uln:nln) = 0.0_WP
            tx(uln:nln) = 0.0_WP; ty(uln:nln) = 0.0_WP
            sx(uln:nln) = 0.0_WP; sy(uln:nln) = 0.0_WP
            do k = 1, mesh%nod_in_elem2D_num(n)
                el  = mesh%nod_in_elem2D(k, n)
                nle = mesh%nlevels(el)-1
                ule = mesh%ulevels(el)
                wa  = mesh%elem_area(el)*mesh%nod_in_elem2D_w(k, n)   ! kite share of el at n
                en  = mesh%elem2D_nodes(:, el)
                do nz = ule, nle
                    vol(nz) = vol(nz) + wa
                    tx(nz) = tx(nz) + (mesh%gradient_sca(1,el)*temp(nz,en(1)) + mesh%gradient_sca(2,el)*temp(nz,en(2)) &
                                     + mesh%gradient_sca(3,el)*temp(nz,en(3)) + mesh%gradient_sca(4,el)*temp(nz,en(4)))*wa
                    ty(nz) = ty(nz) + (mesh%gradient_sca(5,el)*temp(nz,en(1)) + mesh%gradient_sca(6,el)*temp(nz,en(2)) &
                                     + mesh%gradient_sca(7,el)*temp(nz,en(3)) + mesh%gradient_sca(8,el)*temp(nz,en(4)))*wa
                    sx(nz) = sx(nz) + (mesh%gradient_sca(1,el)*salt(nz,en(1)) + mesh%gradient_sca(2,el)*salt(nz,en(2)) &
                                     + mesh%gradient_sca(3,el)*salt(nz,en(3)) + mesh%gradient_sca(4,el)*salt(nz,en(4)))*wa
                    sy(nz) = sy(nz) + (mesh%gradient_sca(5,el)*salt(nz,en(1)) + mesh%gradient_sca(6,el)*salt(nz,en(2)) &
                                     + mesh%gradient_sca(7,el)*salt(nz,en(3)) + mesh%gradient_sca(8,el)*salt(nz,en(4)))*wa
                end do
            end do
            sigma_xy(1,uln:nln,n) = (-sw_alpha(uln:nln,n)*tx(uln:nln)+sw_beta(uln:nln,n)*sx(uln:nln))/vol(uln:nln)*density_0
            sigma_xy(2,uln:nln,n) = (-sw_alpha(uln:nln,n)*ty(uln:nln)+sw_beta(uln:nln,n)*sy(uln:nln))/vol(uln:nln)*density_0
        end do
        if (present(partit)) then
            if (is_multirank(partit)) call exchange_nod(sigma_xy, partit)
        end if
    end subroutine compute_sigma_xy

    !===========================================================================
    subroutine compute_neutral_slope(sigma_xy, bvfreq, mesh, neutral_slope, slope_tapered, fer_tapfac, partit)
        ! Neutral slope + ODM95/LDD97 tapering (FESOM2 oce_ale_pressure_bv.F90:2949).
        ! ro_z_inv = 2*g/density_0/max(bvfreq(nz)+bvfreq(nz+1), eps^2) -> variable-divisor
        ! SIMD trap. ODM95 c1 = 0.5*(1+tanh((Scr-|S|)/Sd)), 0 if N2<=0; LDD97 c2 (off at
        ! scaling_LDD97=.false. -> c2=1). work_core (Fer_GM.and.Redi.and.Redi_Ktaper) branch:
        ! fer_tapfac = c1*c2, slope_tapered = ns*sqrt(c1*c2); else slope_tapered = ns*c1*c2.
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: sigma_xy(2, mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: bvfreq(mesh%nl, mesh%nod2D)
        real(kind=WP), intent(inout) :: neutral_slope(3, mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: slope_tapered(3, mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(inout) :: fer_tapfac(mesh%nl-1, mesh%nod2D)
        type(t_partit), intent(in), optional :: partit
        integer       :: n, nz, nl1, ul1, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: ro_z_inv, eps, f_min, dep_scale, rssby
        real(kind=WP) :: c1(mesh%nl-1), c2(mesh%nl-1)
        eps   = 5.0e-6_WP
        f_min = 1.0e-6_WP
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do n = 1, nNodO
            slope_tapered(:, :, n) = 0.0_WP
            nl1 = mesh%nlevels_nod2D(n)-1
            ul1 = mesh%ulevels_nod2D(n)
            do nz = ul1, nl1
                ! N2 = -g*drho/dz; the minus is hidden in the buoyancy definition
                ro_z_inv = 2.0_WP*g/density_0/max(bvfreq(nz,n)+bvfreq(nz+1,n), eps**2)
                neutral_slope(1,nz,n) = sigma_xy(1,nz,n)*ro_z_inv
                neutral_slope(2,nz,n) = sigma_xy(2,nz,n)*ro_z_inv
                neutral_slope(3,nz,n) = sqrt(neutral_slope(1,nz,n)**2+neutral_slope(2,nz,n)**2)
            end do
            ! ODM95 (Danabasoglu & McWilliams 1995) tanh slope tapering
            c1 = 1.0_WP
            if (scaling_ODM95) then
                do nz = ul1, nl1
                    c1(nz) = 0.5_WP*(1.0_WP + tanh((ODM95_Scr - neutral_slope(3,nz,n))/ODM95_Sd))
                    if ((bvfreq(nz,n) <= 0.0_WP) .or. (bvfreq(nz+1,n) <= 0.0_WP)) c1(nz) = 0.0_WP
                end do
            end if
            ! LDD97 (Large et al. 1997) surface Rossby-radius taper (off in work_core)
            c2 = 1.0_WP
            if (scaling_LDD97) then
                rssby = LDD97_c/max(abs(mesh%coriolis_node(n)), f_min)
                rssby = min(LDD97_rmax, max(LDD97_rmin, rssby))
                do nz = ul1, nl1
                    dep_scale = rssby*neutral_slope(3,nz,n)
                    if (abs(mesh%Z_3d_n(nz,n)) < dep_scale) then
                        c2(nz) = 0.5_WP*(1.0_WP + sin(pi*abs(mesh%Z_3d_n(nz,n))/dep_scale - pi/2.0_WP))
                    end if
                end do
            end if
            ! taper slope with c1*c2; Redi_Ktaper splits sqrt between slope and Ki
            if (Fer_GM .and. Redi .and. Redi_Ktaper) then
                do nz = ul1, nl1
                    fer_tapfac(nz, n) = c1(nz) * c2(nz)
                    slope_tapered(:, nz, n) = neutral_slope(:, nz, n) * sqrt(c1(nz) * c2(nz))
                end do
            else
                do nz = ul1, nl1
                    slope_tapered(:, nz, n) = neutral_slope(:, nz, n) * c1(nz) * c2(nz)
                end do
            end if
        end do
        if (present(partit)) then
            if (is_multirank(partit)) then
                call exchange_nod(neutral_slope, partit)
                call exchange_nod(slope_tapered, partit)
            end if
        end if
    end subroutine compute_neutral_slope

end module oce_pressure_bv
