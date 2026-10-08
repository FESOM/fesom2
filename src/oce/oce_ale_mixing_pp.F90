module oce_ale_mixing_pp
    ! Pacanowski & Philander (1981) Richardson-number vertical mixing. Transcribed
    ! verbatim from FESOM2 v2.7.3 oce_ale_mixing_pp.F90:
    !   oce_mixing_pp        (  1- 93)  Av (elements) + Kv (nodes) from Ri = N^2/(du/dz)^2
    !   Kv0_background_qiang (100-134)  lat/depth background diffusivity (Kv0_const=.false.)
    !   Kv0_background       (141-177)  alt background (NOT called by oce_mixing_pp; kept
    !                                   faithful — same routine ships in the oracle file)
    !
    ! M2.8 SCOPE — the vertical mixing coefficients. Gate target = Kv (nodes) and Av
    ! (elements), the production fields the timestep feeds to impl_vert_visc_ale (Av) and
    ! to the tracer vertical diffusion (Kv). Three SEQUENTIAL passes, the ordering is
    ! load-bearing:
    !   1. node loop: Kv := factor = shear/(shear + 5*max(N^2,0) + 1e-14), the
    !      inverse-Richardson mixing factor 1/(1+5*Ri) rewritten (Ri = N^2/shear),
    !      shear = |d(uvnode)/dz|^2 (vertical velocity shear at the node).
    !   2. elem loop: Av = mix_coeff_PP * mean_over_3_nodes(factor^2) + A_ver.
    !   3. node loop: Kv = mix_coeff_PP * factor^3 + K_ver           (Kv0_const=.true., pi),
    !      or          Kv = mix_coeff_PP * factor^3 + Kv0_background_qiang(lat,depth)
    !                                                                (Kv0_const=.false.).
    ! Pass 2 reads the pass-1 factor BEFORE pass 3 overwrites Kv with the final
    ! diffusivity (Av uses factor^2, Kv uses factor^3 — the same scratch array Kv holds
    ! both the factor and the result, so the elem loop must run between them).
    !
    ! reads:  dyn%uvnode (2,nl-1,nod2D), dyn%work%bvfreq (nl,nod2D, the SMOOTHED N^2 the
    !         step computed in pressure_bv), mesh Z_3d_n / zbar_3d_n / geo_coord_nod2D +
    !         level arrays. params mix_coeff_PP / A_ver / K_ver / Kv0_const (mod_param_phys).
    ! writes: dyn%work%Kv (nl,nod2D), dyn%work%Av (nl,elem2D). The CALLER pre-zeros both
    ! (oce_mixing_pp writes only nz in nzmin+1..nzmax-1 — surface nzmin, bottom nzmax and
    ! below-bottom are left untouched), so the dump's unwritten region is a deterministic 0.
    !
    ! 1-rank: the node loops run 1..nod2D (FESOM2 myDim+eDim, eDim=0); the elem loop runs
    ! 1..elem2D (FESOM2 myDim, no halo elements on 1-rank). oce_mixing_pp has NO halo
    ! exchange (purely local per-node/per-element), so nothing is dropped here. Byte-
    ! identical to FESOM2 by L9 transitivity: every operand is pinned (bvfreq M2.1; uvnode/
    ! Z_3d_n/geometry geom-gated). The element mean of Kv**2 is the kite-weighted one over
    ! the 4 vertex slots (sum(elem_w*Kv(nz,elnodes)**2); FESOM2: plain /3), so triangles
    ! now agree with FESOM2 to roundoff only.
    !
    ! Kv0_const=.false. (Kv0_background_qiang) and the cavity nzmin>1 path are transcribed
    ! faithfully but UNGATED on pi (Kv0_const=.true., no cavity) — re-gated at M2.11.
    use mod_precision,  only: WP
    use mod_mesh,       only: t_mesh, MAX_NV
    use mod_dyn,        only: t_dyn
    use mod_constants,  only: rad
    use mod_param_phys, only: mix_coeff_PP, A_ver, K_ver, Kv0_const
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: oce_mixing_pp, Kv0_background_qiang, Kv0_background

contains

    !===========================================================================
    subroutine oce_mixing_pp(dyn, mesh, partit)
        ! M2.12c: optional partit -> the FESOM2 bounds (oce_ale_mixing_pp.F90): the two
        ! node passes run 1..myDim_nod2D+eDim_nod2D (the elem pass reads the factor at the
        ! element's vertices, which can be halo, so the factor must be valid on the
        ! halo — uvnode is halo-exchanged by compute_vel_nodes, bvfreq computed on the
        ! halo by pressure_bv, so no exchange is needed here); the elem pass runs
        ! 1..myDim_elem2D (Av is consumed at owned elements by impl_vert_visc_ale).
        type(t_mesh), intent(in)            :: mesh
        type(t_dyn),  intent(inout), target :: dyn
        type(t_partit), intent(in), optional :: partit
        integer       :: node, nz, nzmax, nzmin, elem, elnodes(MAX_NV)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: shear, Kv0_b
        real(kind=WP), dimension(:,:),   pointer :: Kv, Av, bvfreq, shear2

        shear2 => dyn%work%shear2
        Kv     => dyn%work%Kv
        Av     => dyn%work%Av
        bvfreq => dyn%work%bvfreq
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !_______________________________________________________________________
        ! Richardson mixing factor (stored temporarily in Kv).
        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            ! ALE: changing zlevel at every node (Z_3d_n is per-node).
            do nz = nzmin+1, nzmax-1
                ! shear now comes from the shared producer (oce_vert_spline::compute_shear2),
                ! which reproduces the previous inline two-point form bit for bit when
                ! shear_splines=.false. and splines it when .true.
                shear  = shear2(nz,node)
                Kv(nz,node) = shear/(shear + 5._WP*max(bvfreq(nz,node),0.0_WP) + 1.0e-14_WP)  ! avoid NaN at start
            end do
        end do

        !_______________________________________________________________________
        ! viscosity (elements): Av = mix_coeff_PP*mean(factor^2) + A_ver (kite-weighted mean).
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:,elem)
            nzmin = mesh%ulevels(elem)
            nzmax = mesh%nlevels(elem)
            do nz = nzmin+1, nzmax-1
                Av(nz,elem) = mix_coeff_PP*sum(mesh%elem_w(:,elem)*Kv(nz,elnodes)**2) + A_ver
            end do
        end do

        !_______________________________________________________________________
        ! diffusivity (nodes): Kv = mix_coeff_PP*factor^3 + Kv0.
        do node = 1, nNodL
            nzmin = mesh%ulevels_nod2D(node)
            nzmax = mesh%nlevels_nod2D(node)
            !___________________________________________________________________
            ! constant background diffusivity (namelist K_ver) — pi default.
            if (Kv0_const) then
                do nz = nzmin+1, nzmax-1
                    Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3 + K_ver
                end do
            !___________________________________________________________________
            ! latitude/depth-dependent background (Q. Wang FESOM1.4) — ungated on pi.
            else
                do nz = nzmin+1, nzmax-1
                    call Kv0_background_qiang(Kv0_b, &
                         real(mesh%geo_coord_nod2D(2,node),WP)/rad, &
                         abs(real(mesh%zbar_3d_n(nz,node),WP)))
                    Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3 + Kv0_b
                end do
            end if
        end do
    end subroutine oce_mixing_pp

    !===========================================================================
    ! Non-constant background diffusion coefficient (also used by KPP ri_iwmix),
    ! after Q. Wang FESOM1.4. Verbatim from FESOM2 oce_ale_mixing_pp.F90:100-134.
    subroutine Kv0_background_qiang(Kv0_b, lat, dep)
        real(kind=WP), intent(out) :: Kv0_b
        real(kind=WP), intent(in)  :: lat, dep      ! latitude [deg], depth [m, >=0]
        real(kind=WP)              :: aux, ratio

        ! latitude/depth-dependent background diffusivity (Q. Wang FESOM1.4).
        aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.0e-5_WP

        ! latitudinal equatorial scaling.
        if (abs(lat) < 5.0_WP) then
            ratio = 1.0_WP
        else
            ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
        end if

        ! latitudinal arctic scaling.
        if (lat > 70.0_WP) then
            if (dep <= 50.0_WP)     then
                ratio = 4.0_WP + 6.0_WP * (50.0_WP - dep) / 50.0_WP
            else
                ratio = 4.0_WP
            endif
        end if
        Kv0_b = aux*ratio
    end subroutine Kv0_background_qiang

    !===========================================================================
    ! Alternative non-constant background (FESOM1.4, Q. Wang). Verbatim from
    ! FESOM2 oce_ale_mixing_pp.F90:141-177. Not called by oce_mixing_pp; kept so the
    ! ported file mirrors the oracle file 1:1.
    subroutine Kv0_background(Kv0_b, lat, dep)
        real(kind=WP), intent(out) :: Kv0_b
        real(kind=WP), intent(in)  :: lat, dep
        real(kind=WP)              :: aux, ratio

        ! latitudinal equatorial scaling.
        if (abs(lat) < 5.0_WP) then
            ratio = 1.0_WP
        else
            ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
        end if

        ! latitudinal <70 deg vs arctic scaling.
        if (lat < 70.0_WP) then
            aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1e-5_WP
        else
            aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.e-6_WP
            ratio = 3.0_WP
            if (dep < 80.0_WP)     then
                ratio = 1.0_WP
            elseif(dep < 100.0_WP) then
                ratio = 2.0_WP
            endif
        end if
        Kv0_b = aux*ratio
    end subroutine Kv0_background

end module oce_ale_mixing_pp
