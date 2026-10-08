module mod_forcing_bulk
    ! Surface-forcing BULK transfer coefficients + wind stress (M2.10b). Transcribed
    ! from FESOM2 v2.7.3:
    !   ncar_ocean_fluxes_mode  (gen_bulk_formulae.F90:126-341)  the LIVE NCAR bulk
    !       (Large&Yeager 2004 + Large 2009 drag, n_itts=5, 3 measurement heights,
    !        Monin-Obukhov stability iteration) -> Cd/Ch/Ce transfer coefficients.
    !   wind stress             (gen_forcing_couple.F90:749-756) stress_atmoce =
    !       Cd*(rho_air*|du|)*du,  du = u_wind-(1-Swind)*u_w.
    !   node->element stress    (ice_oce_coupling.F90:128-149)   with a_ice=0 ->
    !       stress_node_surf=stress_atmoce; stress_surf(elem)=sum(elem_w*node) (kite-
    !       weighted mean over the 4 vertex slots; FESOM2: sum(node,elnodes)/3).
    !
    ! BIT-IDENTITY NOTES (every literal transcribed VERBATIM — same build flags):
    ! * inc_ratio=1.0e-4 (NO _WP), inv_rhoair=1./1.3, tmelt=273.15, rhoair=1.3 are
    !   written WITHOUT a _WP suffix in FESOM2 (MOD_ICE type defaults / the gen_bulk
    !   local) — keep them un-suffixed so the default-real literal rounds identically.
    ! * (ustar*ustar) NOT ustar**2 (line 233); atan(1.0_WP) stays a runtime call;
    !   cd=cd_n10/(1+cd_n10_rt*xx)**2; qs order ((0.98*q1)*inv_rhoair)*exp(q2/ts);
    !   bstar uses (q10+1.0_WP/0.608_WP); test uses (cd+1.0e-8_WP).
    ! * The drag is Large-2009 eq.11 (the u10**6 term + 33 m/s -> 2.34e-3 cap), NOT
    !   the L-Y2004 6a form (commented out). cd_n10_rt (neutral), not cd_rt, in 9a/10.
    ! * The stress + node->elem are gated against the trivial REAL formulas; Cd/Ch/Ce
    !   are gated against the REAL ncar_ocean_fluxes_mode. The node->element mean is the
    !   kite-weighted one (hybrid meshes), so it equals FESOM2 to roundoff only.
    !
    ! SCOPE: transfer coefficients + wind stress only. heat/water flux (obudget,
    ! ice_thermo_oce.F90) is M3. SST + surface ocean velocity are PRESCRIBED inputs
    ! (the M2.5/M2.8 pattern; they come from the ocean/ice state, not ported here).
    use mod_precision, only: WP
    use mod_mesh, only: t_mesh, MAX_NV
    use mod_partit, only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: forcing_bulk_ncar, forcing_wind_stress, forcing_ice_stress, forcing_stress_surf

    ! MULTI-RANK (M3f-4): the bulk transfer coefficients + wind stresses are a per-node
    ! computation (each node reads only its own atm forcing + ocean/ice surface state), so
    ! they are partition-independent — computed over OWNED+HALO (nNodL, FESOM2's
    ! myDim+eDim loop in gen_forcing_couple.F90:703/738) with NO exchange. The optional
    ! partit selects nNodL; absent/npes==1 -> mesh%nod2D (the proven 1-rank path verbatim).

    ! FESOM2 ice%thermo defaults (MOD_ICE.F90:53,67) — VERBATIM un-suffixed literals.
    real(kind=WP), parameter :: bulk_inv_rhoair = 1./1.3
    real(kind=WP), parameter :: bulk_tmelt      = 273.15
    real(kind=WP), parameter :: bulk_rhoair     = 1.3

contains

    ! ---- ncar_ocean_fluxes_mode (gen_bulk_formulae.F90:126-341) ------------------
    subroutine forcing_bulk_ncar(z_wind, z_tair, z_shum, tair, shum, u_wind, v_wind, &
                                 sst, u_w, v_w, cd_oce, ch_oce, ce_oce, mesh, partit)
        real(kind=WP), intent(in)  :: z_wind, z_tair, z_shum   ! ncar_bulk_z_* (10.0 pi)
        real(kind=WP), intent(in)  :: tair(:), shum(:), u_wind(:), v_wind(:)  ! atm (degC/kg-kg/m-s)
        real(kind=WP), intent(in)  :: sst(:), u_w(:), v_w(:)   ! ocean surface T/u/v (prescribed)
        real(kind=WP), intent(out) :: cd_oce(:), ch_oce(:), ce_oce(:)
        type(t_mesh),  intent(in)  :: mesh
        type(t_partit), intent(in), optional :: partit
        integer, parameter :: n_itts = 5
        integer :: i, j, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: cd_n10, ce_n10, ch_n10, cd_n10_rt, hl1
        real(kind=WP) :: cd, ce, ch, cd_rt
        real(kind=WP) :: x2, x, stab
        real(kind=WP) :: zeta_u, zeta_t, zeta_q
        real(kind=WP) :: psi_m_u, psi_h_u, psi_m_t, psi_h_t, psi_m_q, psi_h_q
        real(kind=WP) :: ts, qs, tv, xx, dux, dvy
        real(kind=WP) :: t, t10, q, q10, u, u10
        real(kind=WP) :: tstar, qstar, ustar, bstar
        real(kind=WP), parameter :: grav = 9.80_WP, vonkarm = 0.40_WP
        real(kind=WP), parameter :: q1 = 640380._WP, q2 = -5107.4_WP
        real(kind=WP), parameter :: u10min = 0.3_WP
        real(kind=WP) :: test, cd_prev, inc_ratio = 1.0e-4

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do i = 1, nNodL
            if (mesh%ulevels_nod2d(i) > 1) cycle
            t  = tair(i) + bulk_tmelt
            ts = sst(i)  + bulk_tmelt
            q  = shum(i)
            qs = 0.98_WP*q1*bulk_inv_rhoair*exp(q2/ts)                  ! L-Y eqn. 5
            tv = t*(1.0_WP + 0.608_WP*q)
            dux = u_wind(i) - u_w(i)
            dvy = v_wind(i) - v_w(i)
            u   = max(sqrt(dux**2 + dvy**2), u10min)
            u10 = u
            t10 = t
            q10 = q
            hl1    = (2.7_WP/u10 + 0.142_WP + 0.0764_WP*u10 - 3.14807e-10_WP*(u10**6)) / 1.0e3_WP
            cd_n10 = (0.5_WP - sign(0.5_WP, u10-33.0_WP)) * hl1 &
                   + (0.5_WP + sign(0.5_WP, u10-33.0_WP)) * 2.34e-3_WP
            cd_n10_rt = sqrt(cd_n10)
            ce_n10 = 34.6_WP*cd_n10_rt*1.0e-3_WP
            stab   = 0.5_WP + sign(0.5_WP, t-ts)
            ch_n10 = (18.0_WP*stab + 32.7_WP*(1.0_WP-stab))*cd_n10_rt*1.e-3_WP
            cd = cd_n10
            ch = ch_n10
            ce = ce_n10
            cd_prev = cd
            do j = 1, n_itts
                cd_rt = sqrt(cd)
                ustar = cd_rt*u                                        ! 7a
                tstar = (ch/cd_rt)*(t10-ts)                            ! 7b
                qstar = (ce/cd_rt)*(q10-qs)                            ! 7c
                bstar = grav*(tstar/tv + qstar/(q10 + 1.0_WP/0.608_WP))
                ! (2a) z_wind
                zeta_u = vonkarm*bstar*z_wind/(ustar*ustar)           ! 8a
                zeta_u = sign(min(abs(zeta_u), 10.0_WP), zeta_u)
                x2 = sqrt(abs(1._WP - 16._WP*zeta_u))                 ! 8b
                x2 = max(x2, 1.0_WP)
                x  = sqrt(x2)
                if (zeta_u > 0._WP) then
                    psi_m_u = -5._WP*zeta_u                           ! 8c
                    psi_h_u = -5._WP*zeta_u
                else
                    psi_m_u = log((1._WP+2._WP*x+x2)*(1.0_WP+x2)/8._WP) - 2._WP*(atan(x)-atan(1.0_WP)) ! 8d
                    psi_h_u = 2._WP*log((1._WP+x2)/2._WP)             ! 8e
                end if
                ! (2b) z_tair
                zeta_t = vonkarm*bstar*z_tair/(ustar*ustar)
                zeta_t = sign(min(abs(zeta_t), 10.0_WP), zeta_t)
                x2 = sqrt(abs(1._WP - 16._WP*zeta_t))
                x2 = max(x2, 1.0_WP)
                x  = sqrt(x2)
                if (zeta_t > 0._WP) then
                    psi_m_t = -5._WP*zeta_t
                    psi_h_t = -5._WP*zeta_t
                else
                    psi_m_t = log((1._WP+2._WP*x+x2)*(1.0_WP+x2)/8._WP) - 2._WP*(atan(x)-atan(1.0_WP))
                    psi_h_t = 2._WP*log((1._WP+x2)/2._WP)
                end if
                ! (2c) z_shum
                zeta_q = vonkarm*bstar*z_shum/(ustar*ustar)
                zeta_q = sign(min(abs(zeta_q), 10.0_WP), zeta_q)
                x2 = sqrt(abs(1._WP - 16._WP*zeta_q))
                x2 = max(x2, 1.0_WP)
                x  = sqrt(x2)
                if (zeta_q > 0._WP) then
                    psi_m_q = -5._WP*zeta_q
                    psi_h_q = -5._WP*zeta_q
                else
                    psi_m_q = log((1._WP+2._WP*x+x2)*(1.0_WP+x2)/8._WP) - 2._WP*(atan(x)-atan(1.0_WP))
                    psi_h_q = 2._WP*log((1._WP+x2)/2._WP)
                end if
                ! (3a) shift wind to 10m neutral; (3b) shift T,q to wind height
                u10 = u/(1.0_WP + cd_n10_rt*(log(z_wind/10._WP)-psi_m_u)/vonkarm) ! 9a
                u10 = max(u10, u10min)
                t10 = t - tstar/vonkarm*(log(z_tair/z_wind)+psi_h_u-psi_h_t)      ! 9b
                q10 = q - qstar/vonkarm*(log(z_shum/z_wind)+psi_h_u-psi_h_q)      ! 9b
                ! (3c) update tv
                tv = t10*(1.0_WP + 0.608_WP*q10)
                ! (4a) update neutral coeffs
                hl1    = (2.7_WP/u10 + 0.142_WP + 0.0764_WP*u10 - 3.14807e-10_WP*(u10**6)) / 1.0e3_WP
                cd_n10 = (0.5_WP - sign(0.5_WP, u10-33.0_WP)) * hl1 &
                       + (0.5_WP + sign(0.5_WP, u10-33.0_WP)) * 2.34e-3_WP
                cd_n10_rt = sqrt(cd_n10)
                ce_n10 = 34.6_WP*cd_n10_rt*1.e-3_WP
                stab   = 0.5_WP + sign(0.5_WP, zeta_u)
                ch_n10 = (18.0_WP*stab + 32.7_WP*(1.0_WP-stab))*cd_n10_rt*1.e-3_WP
                ! (4b) shift to measurement height + stability
                xx = (log(z_wind/10._WP)-psi_m_u)/vonkarm
                cd = cd_n10/(1.0_WP + cd_n10_rt*xx)**2                 ! 10a
                xx = (log(z_wind/10._WP)-psi_h_u)/vonkarm
                ch = ch_n10/(1.0_WP + ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)   ! 10b
                ce = ce_n10/(1.0_WP + ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)   ! 10c
                ! (5) convergence
                test = abs(cd - cd_prev) / (cd + 1.0e-8_WP)
                cd_prev = cd
                if (test < inc_ratio) exit
            end do
            cd_oce(i) = cd
            ch_oce(i) = ch
            ce_oce(i) = ce
        end do
    end subroutine forcing_bulk_ncar

    ! ---- wind stress on nodes (gen_forcing_couple.F90:749-756) -------------------
    subroutine forcing_wind_stress(swind, u_wind, v_wind, u_w, v_w, cd_oce, &
                                   stress_x, stress_y, mesh, partit)
        real(kind=WP), intent(in)  :: swind
        real(kind=WP), intent(in)  :: u_wind(:), v_wind(:), u_w(:), v_w(:), cd_oce(:)
        real(kind=WP), intent(out) :: stress_x(:), stress_y(:)
        type(t_mesh),  intent(in)  :: mesh
        type(t_partit), intent(in), optional :: partit
        integer :: i, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dux, dvy, aux
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do i = 1, nNodL
            if (mesh%ulevels_nod2d(i) > 1) then
                stress_x(i) = 0.0_WP; stress_y(i) = 0.0_WP; cycle
            end if
            dux = u_wind(i) - (1.0_WP-swind)*u_w(i)
            dvy = v_wind(i) - (1.0_WP-swind)*v_w(i)
            aux = sqrt(dux**2 + dvy**2)*bulk_rhoair
            stress_x(i) = cd_oce(i)*aux*dux
            stress_y(i) = cd_oce(i)*aux*dvy
        end do
    end subroutine forcing_wind_stress

    ! ---- wind-on-ICE stress on nodes (gen_forcing_couple.F90:759-763) ------------
    ! stress_atmice = Cd_atm_ice*(rhoair*|u_wind-u_ice|)*(u_wind-u_ice). Cd_atm_ice is the
    ! CONSTANT namelist drag (0.0012; AOMIP_drag_coeff=.false. => no cal_wind_drag_coeff),
    ! NOT the bulk Cd_atm_oce_arr. Cavity nodes -> 0 (matches the combined FESOM2 loop's
    ! ulevels guard at :740, which zeros BOTH stresses).
    subroutine forcing_ice_stress(cd_atm_ice, u_wind, v_wind, u_ice, v_ice, &
                                  stress_x, stress_y, mesh, partit)
        real(kind=WP), intent(in)  :: cd_atm_ice
        real(kind=WP), intent(in)  :: u_wind(:), v_wind(:), u_ice(:), v_ice(:)
        real(kind=WP), intent(out) :: stress_x(:), stress_y(:)
        type(t_mesh),  intent(in)  :: mesh
        type(t_partit), intent(in), optional :: partit
        integer :: i, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: dux, dvy, aux
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        do i = 1, nNodL
            if (mesh%ulevels_nod2d(i) > 1) then
                stress_x(i) = 0.0_WP; stress_y(i) = 0.0_WP; cycle
            end if
            dux = u_wind(i) - u_ice(i)
            dvy = v_wind(i) - v_ice(i)
            aux = sqrt(dux**2 + dvy**2)*bulk_rhoair
            stress_x(i) = cd_atm_ice*aux*dux
            stress_y(i) = cd_atm_ice*aux*dvy
        end do
    end subroutine forcing_ice_stress

    ! ---- node->element surface stress, a_ice=0 (ice_oce_coupling.F90:128-149) -----
    ! 1-rank (fesom_forcingdump): elem2D == the owned elements elem_w is allocated for.
    subroutine forcing_stress_surf(stress_x, stress_y, stress_surf, mesh)
        real(kind=WP), intent(in)  :: stress_x(:), stress_y(:)   ! = stress_node_surf (a_ice=0)
        real(kind=WP), intent(out) :: stress_surf(:,:)           ! (2, elem2D)
        type(t_mesh),  intent(in)  :: mesh
        integer :: elem, elnodes(MAX_NV)
        do elem = 1, mesh%elem2D
            if (mesh%ulevels(elem) > 1) cycle
            elnodes = mesh%elem2D_nodes(:, elem)
            stress_surf(1, elem) = sum(mesh%elem_w(:, elem)*stress_x(elnodes))
            stress_surf(2, elem) = sum(mesh%elem_w(:, elem)*stress_y(elnodes))
        end do
    end subroutine forcing_stress_surf

end module mod_forcing_bulk
