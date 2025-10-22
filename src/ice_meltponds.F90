!===============================================================================
! FESOM-2.6 Melt Pond Parameterization Module
! 
! Based on CESM melt pond scheme from ICEPACK
! Adapted for use with ice_thermo_cpl.F90
!
! Author: Extracted from ICEPACK for FESOM-2.6 coupled simulations
! Date: 2024
!===============================================================================

module ice_meltponds
    use o_param, only: WP
    implicit none
    
    private
    public :: meltpond_area, meltpond_albedo, init_meltponds
    
    ! Melt pond parameters (similar to ICEPACK ponds_nml)
    real(kind=WP), parameter :: hp1 = 0.01_WP        ! Critical pond depth [m]
    real(kind=WP), parameter :: hs0 = 0.0_WP         ! Snow depth controlling drainage [m] 
    real(kind=WP), parameter :: hs1 = 0.03_WP        ! Snow depth threshold [m]
    real(kind=WP), parameter :: rfracmin = 0.15_WP   ! Minimum retained melt fraction
    real(kind=WP), parameter :: rfracmax = 1.0_WP    ! Maximum retained melt fraction
    real(kind=WP), parameter :: pndaspect = 0.8_WP   ! Pond aspect ratio
    real(kind=WP), parameter :: dpscale = 1.e-3_WP   ! Drainage timescale parameter
    
    ! Pond albedo parameters
    real(kind=WP), parameter :: albpnd = 0.2_WP      ! Pond albedo (open water-like)
    real(kind=WP), parameter :: albpnd_frz = 0.36_WP ! Frozen pond albedo
    
    ! Physical constants
    real(kind=WP), parameter :: puny = 1.e-11_WP     ! Small number
    real(kind=WP), parameter :: c0 = 0.0_WP
    real(kind=WP), parameter :: c1 = 1.0_WP
    
    contains
    
    !===========================================================================
    ! Initialize melt pond variables
    !===========================================================================
    subroutine init_meltponds()
        implicit none
        ! Currently no initialization needed - parameters are set above
        ! This routine is available for future extensions
    end subroutine init_meltponds
    
    !===========================================================================
    ! Calculate melt pond area fraction
    ! 
    ! Based on CESM melt pond scheme from ICEPACK
    !===========================================================================
    subroutine meltpond_area(aice, hice, hsno, meltt, melts, dt, &
                            apnd, hpnd, ipnd, fpond)
        implicit none
        
        ! Input variables
        real(kind=WP), intent(in) :: aice     ! Ice concentration
        real(kind=WP), intent(in) :: hice     ! Ice thickness [m]
        real(kind=WP), intent(in) :: hsno     ! Snow thickness [m]  
        real(kind=WP), intent(in) :: meltt    ! Top melt rate [m/s]
        real(kind=WP), intent(in) :: melts    ! Snow melt rate [m/s]
        real(kind=WP), intent(in) :: dt       ! Time step [s]
        
        ! Input/Output variables
        real(kind=WP), intent(inout) :: apnd  ! Pond area fraction
        real(kind=WP), intent(inout) :: hpnd  ! Pond depth [m]
        real(kind=WP), intent(inout) :: ipnd  ! Pond ice thickness [m]
        
        ! Output variables
        real(kind=WP), intent(out) :: fpond   ! Pond fresh water flux [m/s]
        
        ! Local variables
        real(kind=WP) :: volpnd               ! Pond volume per unit area [m]
        real(kind=WP) :: volpnd_init          ! Initial pond volume
        real(kind=WP) :: hi_min               ! Minimum ice thickness for ponds
        real(kind=WP) :: rfrac                ! Retained fraction of melt water
        real(kind=WP) :: dvolpnd              ! Change in pond volume
        real(kind=WP) :: armor                ! Armor parameter
        real(kind=WP) :: pp                   ! Pond parameter
        real(kind=WP) :: apondn               ! New pond area
        real(kind=WP) :: hpondn               ! New pond depth
        
        ! Initialize
        fpond = c0
        hi_min = 0.1_WP  ! Minimum ice thickness for pond formation
        
        ! No ponds on thin ice or high snow cover
        if (aice < puny .or. hice < hi_min .or. hsno > hs1) then
            apnd = c0
            hpnd = c0 
            ipnd = c0
            return
        endif
        
        ! Calculate retained fraction of melt water
        ! More melt water retained on thicker ice
        armor = min(c1, hice / 0.5_WP)  ! Armor based on ice thickness
        rfrac = rfracmin + (rfracmax - rfracmin) * armor
        
        ! Current pond volume per unit area
        volpnd_init = apnd * hpnd
        
        ! Add melt water (both surface and snow melt)
        dvolpnd = rfrac * (meltt + melts) * dt
        volpnd = volpnd_init + dvolpnd
        
        ! Pond drainage - simple exponential decay when snow is thin
        if (hsno < hs0 .and. volpnd > puny) then
            volpnd = volpnd * exp(-dt * dpscale)
        endif
        
        ! Calculate new pond geometry
        if (volpnd > puny) then
            ! CESM pond area calculation
            pp = 0.5_WP * pndaspect * aice
            apondn = (-c1 + sqrt(c1 + 4.0_WP * pp * volpnd / hp1)) / (2.0_WP * pp)
            apondn = max(puny, min(aice * 0.9_WP, apondn))
            
            if (apondn > puny) then
                hpondn = volpnd / apondn
                hpondn = max(puny, min(c1, hpondn))  ! Limit pond depth
            else
                apondn = c0
                hpondn = c0
            endif
        else
            apondn = c0
            hpondn = c0
            volpnd = c0
        endif
        
        ! Update pond variables
        apnd = apondn
        hpnd = hpondn
        
        ! Calculate fresh water flux (pond volume change)
        fpond = (volpnd - volpnd_init) / dt
        
        ! Simple pond ice formation/melting (frozen lid)
        if (meltt < c0 .and. apnd > puny) then
            ! Freezing conditions - form/thicken pond ice
            ipnd = ipnd - meltt * dt  ! meltt is negative for freezing
            ipnd = max(c0, ipnd)
        else
            ! Melting conditions - melt pond ice
            ipnd = max(c0, ipnd + meltt * dt)
        endif
        
    end subroutine meltpond_area
    
    !===========================================================================
    ! Calculate effective albedo including melt pond effects
    !===========================================================================
    subroutine meltpond_albedo(aice, hsno, apnd, ipnd, t_sfc, &
                              albi_dry, albsn_dry, albedo_eff)
        implicit none
        
        ! Input variables
        real(kind=WP), intent(in) :: aice       ! Ice concentration
        real(kind=WP), intent(in) :: hsno       ! Snow thickness [m]
        real(kind=WP), intent(in) :: apnd       ! Pond area fraction
        real(kind=WP), intent(in) :: ipnd       ! Pond ice thickness [m]
        real(kind=WP), intent(in) :: t_sfc      ! Surface temperature [K]
        real(kind=WP), intent(in) :: albi_dry   ! Dry ice albedo
        real(kind=WP), intent(in) :: albsn_dry  ! Snow albedo
        
        ! Output
        real(kind=WP), intent(out) :: albedo_eff ! Effective albedo
        
        ! Local variables
        real(kind=WP) :: albedo_ice              ! Ice surface albedo
        real(kind=WP) :: albedo_pond             ! Pond albedo
        real(kind=WP) :: frac_ice                ! Fraction of ice not covered by ponds
        real(kind=WP) :: frac_pond               ! Effective pond fraction
        real(kind=WP) :: frac_open               ! Open water fraction
        real(kind=WP) :: albw                    ! Open water albedo
        
        ! Initialize
        albw = 0.066_WP  ! Open water albedo
        
        if (aice < puny) then
            albedo_eff = albw
            return
        endif
        
        ! Determine ice surface albedo (ice or snow)
        if (hsno > 0.001_WP) then
            albedo_ice = albsn_dry  ! Snow-covered ice
        else
            albedo_ice = albi_dry   ! Bare ice
        endif
        
        ! Determine pond albedo based on pond ice thickness
        if (ipnd > 0.01_WP) then
            ! Thick pond ice - use frozen pond albedo
            albedo_pond = albpnd_frz
        else
            ! Thin or no pond ice - use open pond albedo
            albedo_pond = albpnd
        endif
        
        ! Calculate area fractions
        frac_pond = min(apnd, aice)  ! Ponds only exist on ice
        frac_ice = aice - frac_pond   ! Ice not covered by ponds
        frac_open = c1 - aice         ! Open water
        
        ! Weight albedos by area fractions
        albedo_eff = frac_ice * albedo_ice + &
                     frac_pond * albedo_pond + &
                     frac_open * albw
        
        ! Ensure reasonable bounds
        albedo_eff = max(albw, min(0.9_WP, albedo_eff))
        
    end subroutine meltpond_albedo
    
end module ice_meltponds
