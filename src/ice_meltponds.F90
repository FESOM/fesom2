!===============================================================================
! FESOM-2.7 Melt Pond Parameterization Module
!
! Based on Icepack `meltpond_lvl` (Stefan/'hlid' branch) — energy-balance lid
! growth/melt; pond geometry follows the CESM aspect-ratio formulation. The
! pond-ice lid (`ipnd`) is bounded by physical constraints from the Stefan
! approximation; the pond depth is bounded by ice freeboard. This replaces an
! earlier home-grown variant in which `ipnd` accumulated `meltt*dt` without
! a thermal-conduction limit and could grow to multi-metre values.
!
! References:
!   Hunke et al., 2015: CICE the Los Alamos sea ice model documentation.
!   Icepack columnphysics/icepack_meltpond_lvl.F90 (hlid branch).
!===============================================================================

module ice_meltponds
    use o_param, only: WP
    implicit none

    private
    public :: meltpond_area, meltpond_albedo, init_meltponds

    ! Melt pond parameters (similar to Icepack ponds_nml)
    real(kind=WP), parameter :: hi_min     = 0.10_WP  ! Min ice thickness for ponds [m]
    real(kind=WP), parameter :: hs1        = 0.03_WP  ! Snow depth threshold for pond reset [m]
    real(kind=WP), parameter :: rfracmin   = 0.15_WP  ! Min retained melt fraction
    real(kind=WP), parameter :: rfracmax   = 0.50_WP  ! Max retained melt fraction
    real(kind=WP), parameter :: pndaspect  = 0.80_WP  ! Pond depth/area aspect ratio

    ! Pond albedo (CICE/CESM range)
    real(kind=WP), parameter :: albpnd     = 0.30_WP  ! Open melt pond albedo
    real(kind=WP), parameter :: albpnd_frz = 0.50_WP  ! Frozen pond albedo

    ! Physical constants
    real(kind=WP), parameter :: puny     = 1.e-11_WP
    real(kind=WP), parameter :: c0       = 0.0_WP
    real(kind=WP), parameter :: c1       = 1.0_WP
    real(kind=WP), parameter :: p5       = 0.5_WP
    real(kind=WP), parameter :: c2       = 2.0_WP
    real(kind=WP), parameter :: Tffresh  = 273.15_WP   ! K
    ! Pond / sea-ice / sea-water densities and physics constants
    real(kind=WP), parameter :: rhoice_p   = 917.0_WP    ! kg m^-3 (pond context)
    real(kind=WP), parameter :: rhosno_p   = 330.0_WP    ! kg m^-3
    real(kind=WP), parameter :: rhowat_p   = 1025.0_WP   ! kg m^-3
    real(kind=WP), parameter :: rhofresh   = 1000.0_WP   ! kg m^-3
    real(kind=WP), parameter :: kice       = 2.03_WP     ! W m^-1 K^-1, lid conductivity
    real(kind=WP), parameter :: Lfresh     = 3.34e5_WP   ! J kg^-1 latent heat fusion (fresh)

contains

    !===========================================================================
    subroutine init_meltponds()
        implicit none
    end subroutine init_meltponds

    !===========================================================================
    ! Update pond area/depth/lid for one node, one time step.
    !
    !   meltt, melts:  top-ice / snow melt rates (>=0, m s^-1 of fresh water)
    !   frain      :  rain rate (m s^-1 fresh water; converted to mass below)
    !   Tsfc_K     :  ice surface temperature (K); used for lid Stefan growth
    !   fsurfn     :  atm->ice surface heat flux (W m^-2, positive downward)
    !   dt         :  time step (s)
    !   apnd,hpnd,ipnd : pond area frac, pond depth (m), lid thickness (m)
    !   fpond      :  diagnostic pond freshwater flux (m s^-1)
    !===========================================================================
    subroutine meltpond_area(aice, hice, hsno, &
                             meltt, melts, frain, &
                             Tsfc_K, fsurfn,      &
                             dt, apnd, hpnd, ipnd, fpond)
        implicit none

        real(kind=WP), intent(in)    :: aice, hice, hsno
        real(kind=WP), intent(in)    :: meltt, melts, frain
        real(kind=WP), intent(in)    :: Tsfc_K, fsurfn
        real(kind=WP), intent(in)    :: dt
        real(kind=WP), intent(inout) :: apnd, hpnd, ipnd
        real(kind=WP), intent(out)   :: fpond

        ! local
        real(kind=WP) :: volpn, volpn_init, dvn
        real(kind=WP) :: hi, hs
        real(kind=WP) :: armor, rfrac
        real(kind=WP) :: apondn, hpondn
        real(kind=WP) :: hlid, dhlid, bdt, Ts
        real(kind=WP) :: alid, fbcap

        fpond = c0

        ! No ponds on thin ice / heavy snow / no ice
        if (aice < puny .or. hice < hi_min .or. hsno > hs1) then
            apnd = c0; hpnd = c0; ipnd = c0
            return
        endif

        ! Per-ice quantities (state passed as grid-cell-mean: hice = aice*hi)
        hi = hice / max(aice, puny)
        hs = hsno / max(aice, puny)

        ! Retained melt fraction (more retained on thicker ice)
        armor = min(c1, hi / 0.50_WP)
        rfrac = rfracmin + (rfracmax - rfracmin) * armor

        ! Initial pond volume per unit grid area  (m of fresh water)
        volpn_init = apnd * hpnd * aice
        volpn      = volpn_init
        hlid       = ipnd

        ! Melt-water input (ice melt -> fresh, snow melt -> fresh, rain mass)
        dvn = rfrac/rhofresh * ( meltt*rhoice_p  &
            +                    melts*rhosno_p  &
            +                    frain*rhofresh ) * dt * aice

        !-----------------------------------------------------------------------
        ! Stefan lid growth / melt  (Icepack `meltpond_lvl` hlid branch)
        ! Pond is at melt temperature (~0 °C); lid grows if surface is colder.
        !-----------------------------------------------------------------------
        Ts = Tsfc_K - Tffresh
        if (dvn <= c0) then
            ! Freeze branch: cold surface, no fresh melt water this step
            if (Ts < c0) then
                bdt   = -c2 * Ts * kice * dt / (rhoice_p * Lfresh)
                dhlid = p5 * sqrt(max(bdt, c0))            ! open water freezing
                if (hlid > dhlid) dhlid = p5 * bdt / hlid  ! conduction through existing lid
                dhlid = min(dhlid, hpnd * rhofresh/rhoice_p) ! can't freeze more than pond holds
                hlid  = hlid + dhlid
            else
                dhlid = c0
            endif
        else
            ! Melt branch: warm surface, energy thins the lid
            dhlid = max(fsurfn*dt / (rhoice_p * Lfresh), c0)  ! >= 0 magnitude
            dhlid = -min(dhlid, hlid)                          ! <= 0, capped at existing lid
            hlid  = max(hlid + dhlid, c0)
        endif
        ! Lid grow/melt redistributes pond water vs. ice mass
        alid = apnd * aice                          ! lid sits on pond area
        dvn  = dvn - dhlid * alid * rhoice_p/rhofresh

        volpn = volpn + dvn

        !-----------------------------------------------------------------------
        ! Pond geometry — CESM aspect-ratio formulation
        !-----------------------------------------------------------------------
        if (volpn > puny) then
            apondn = min(sqrt(volpn / (pndaspect * aice)), aice * 0.9_WP)
            apondn = max(apondn, puny)
            hpondn = volpn / max(apondn, puny)
            ! Limit pond depth to maintain positive freeboard:
            !   draft <= hi  =>  hpondn <= ((rhow-rhoi)*hi - rhos*hs) / rhofresh
            fbcap  = ((rhowat_p - rhoice_p) * hi - rhosno_p * hs) / rhofresh
            hpondn = min(hpondn, max(fbcap, c0))
        else
            apondn = c0
            hpondn = c0
            volpn  = c0
            hlid   = c0
        endif

        ! Lid cannot exceed pond depth (lid sits on water column)
        hlid = min(hlid, hpondn)

        ! Write back
        apnd  = apondn
        hpnd  = hpondn
        ipnd  = hlid
        fpond = (volpn - volpn_init) / dt
    end subroutine meltpond_area

    !===========================================================================
    ! Albedo with pond contribution (binary lid threshold at 1 cm)
    !===========================================================================
    subroutine meltpond_albedo(aice, hsno, apnd, ipnd, t_sfc, &
                               albi_dry, albsn_dry, albedo_eff)
        implicit none
        real(kind=WP), intent(in)  :: aice, hsno, apnd, ipnd, t_sfc
        real(kind=WP), intent(in)  :: albi_dry, albsn_dry
        real(kind=WP), intent(out) :: albedo_eff

        real(kind=WP) :: albedo_ice, albedo_pond
        real(kind=WP) :: frac_ice, frac_pond, frac_open
        real(kind=WP), parameter :: albw = 0.066_WP

        if (aice < puny) then
            albedo_eff = albw
            return
        endif

        if (hsno > 0.001_WP) then
            albedo_ice = albsn_dry
        else
            albedo_ice = albi_dry
        endif

        if (ipnd > 0.01_WP) then
            albedo_pond = albpnd_frz
        else
            albedo_pond = albpnd
        endif

        frac_pond = min(apnd, aice)
        frac_ice  = aice - frac_pond
        frac_open = c1 - aice

        albedo_eff = frac_ice  * albedo_ice  + &
                     frac_pond * albedo_pond + &
                     frac_open * albw
        albedo_eff = max(albw, min(0.9_WP, albedo_eff))
    end subroutine meltpond_albedo

end module ice_meltponds
