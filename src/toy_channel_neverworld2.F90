MODULE Toy_Neverworld2
    use mod_mesh
    use o_ARRAYS
    use o_PARAM
    use MOD_PARSUP
    use MOD_PARTIT
    use MOD_TRACER
    use MOD_DYN
    use g_config
    use g_comm_auto
    use g_support
    use g_clock, only: daynew, timenew, ndpyr
    implicit none
    SAVE 
    private
    public            :: initial_state_neverworld2, relax_2_tsurf, relax_2_ssurf, oce_mixing_TOY,  &
                        do_wind, wind_opt, tau_inv, do_Trelax, do_Tpert, trelax_opt, &
                        gamma_restore, do_north_cold_patch, north_pole_target, north_cold_patch_power, &
                        do_south_cold_patch, south_pole_target, south_cold_patch_power, &
                        equator_target, &
                        thermal_alpha, &
                        do_cabbeling, cabbeling_Cb, &
                        do_thermobar, thermobaric_Th, &
                        do_Sstrat, Sstrat_amp, &
                        do_haline, haline_beta, &
                        do_Srelax, s_north, s_south, s_equator, &
                        do_Tseasonal, north_seasonal_amp, south_seasonal_amp, &
                        oce_neverworld2
    logical           :: do_wind  = .True.    ! apply surface windstress
    integer           :: wind_opt = 2         ! 1: interpolate tau from profile data, 2: read already to elem interp profile data

    logical           :: do_Tpert = .False.    ! apply temp. perturbation to trigger instabilities

    ! Two alternative SST restoring implementations, picked by trelax_opt (same pattern as wind_opt):
    !   1: original -- directly nudges the temperature tracer: T += dt*tau_inv*(Ttarget-T)
    !   2: Genevieve's -- applies restoring as a surface heat flux instead, so it is visible
    !      in the heat budget: heat_flux = gamma_restore*(Tsurface-Ttarget) [W/m^2/K].
    !      FESOM heat_flux convention: >0 cools the ocean (heat leaving), <0 warms it.
    logical           :: do_Trelax= .True.   ! apply surface temp relaxation
    integer           :: Trelax_opt = 1
    real(kind=WP)     :: tau_inv  =1.0/50.0/24.0/3600.0   ! used by trelax_opt=1, i.e. tau=1/tau_inv=50 days
    !
    ! gamma_restore = 40.0 W/m^2/K is Genevieve's original default, carried over as-is.
    ! It is NOT independently derived/tuned for this mesh/config -- treat it as a
    ! placeholder and revisit before relying on it scientifically.
    !
    ! To pick a value from a target restoring timescale instead, note that a heat
    ! flux of gamma_restore*(Tsurface-Ttarget) applied over a layer of thickness h
    ! is equivalent to a Newtonian relaxation (trelax_opt=1 style) with timescale
    !     tau_eff = h*vcpw / gamma_restore                 (vcpw = 4.2e6 J/m^3/K)
    ! i.e.
    !     gamma_restore = h*vcpw / tau_eff
    ! For comparison against trelax_opt=1's tau=50 days (the tau_inv default above):
    !   h=50m, tau_eff=50 days -> gamma_restore = 50*4.2e6/(50*86400) = 48.6 W/m^2/K
    ! i.e. Genevieve's 40 W/m^2/K already corresponds to tau_eff=h*vcpw/40=61 days at
    ! h=50m -- in the same ballpark as the original 50-day default, not wildly off.
    real(kind=WP)     :: gamma_restore = 40.0_WP           ! used by trelax_opt=2

    ! Optional asymmetric cold patch near the northern and/or southern edge, only used with
    ! trelax_opt=2: blends the symmetric cosine restoring target towards north_pole_target /
    ! south_pole_target as ynorm=lat/Ly -> +-1. Each edge is independently switchable/tunable.
    ! Off by default.
    logical           :: do_north_cold_patch    = .False.
    real(kind=WP)     :: north_pole_target      = 0.0_WP
    real(kind=WP)     :: north_cold_patch_power = 2.0_WP
    
    logical           :: do_south_cold_patch    = .False.
    real(kind=WP)     :: south_pole_target      = 0.0_WP
    real(kind=WP)     :: south_cold_patch_power = 2.0_WP

    ! Equatorial value of trelax_opt=2's symmetric cosine restoring target,
    ! t_base = 5 + (equator_target-5)*cos(0.5*pi*lat/Ly) -- the poleward endpoint
    ! (5degC, before any cold-patch blending) stays hardcoded, only the equatorial
    ! value is exposed. Was hardcoded to 33degC (5+28); default unchanged, now
    ! namelist-tunable -- e.g. to bring it closer to DINO's Theta*_eq=27degC
    ! (Kamm et al. 2025), which is closer to realistic tropical SST than 33degC.
    real(kind=WP)     :: equator_target = 33.0_WP

    ! Linear thermal expansion coefficient for the neverworld2 branch of density_linear()
    ! (rho = rho0 - rho0*thermal_alpha*(T-10)) and the matching constant sw_alpha in
    ! sw_alpha_beta() -- same normalized convention as dbgyre's 0.0002052. Was hardcoded
    ! to 0.0002 before; default unchanged, now namelist-tunable. Also the base slope used
    ! by the do_cabbeling quadratic term below (a0 = thermal_alpha*rho0), and the base
    ! sw_alpha that the do_cabbeling per-node alpha is offset from in sw_alpha_beta().
    real(kind=WP)     :: thermal_alpha = 0.0002_WP   ! [1/K], normalized like sw_alpha

    ! Optional thermal cabbeling term for the neverworld2 branch of density_linear()
    ! (oce_ale_pressure_bv.F90), only used when state_equation==0 -- state_equation itself
    ! is NOT a free case selector there (it's used as a literal 0/1 pressure-correction
    ! multiplier throughout that file), so this stays a toy-level flag instead of a new
    ! state_equation value. Adds curvature (d^2rho/dT^2 /= 0) that the plain linear EOS
    ! lacks by construction, following the quadratic-in-Theta cabbeling term of the DINO
    ! configuration's simplified nonlinear EOS (Roquet et al. 2015, Eq.6; parameters from
    ! Caneill et al. 2022): rho = rho0 - (a0 + 0.5*Cb*Ta)*Ta, with Ta = T-10degC. The
    ! existing linear slope (a0 = thermal_alpha*rho0) is kept unchanged; only Cb is added
    ! on top, so do_cabbeling=.False. reproduces the old formula exactly. DINO's haline term (b0)
    ! is deliberately NOT added -- neverworld2 keeps salinity dynamically inert (beta=0),
    ! this isolates thermal cabbeling only. sw_alpha_beta() must mirror this formula's local
    ! d(rho)/dT when do_cabbeling=.True. (see there) -- same consistency requirement noted
    ! for the plain linear case above it. cabbeling_Cb is rho0-NORMALIZED [K^-2], same
    ! convention as thermal_alpha [K^-1] -- NOT DINO/Roquet's raw absolute Cb [kg/m^3/K^2]
    ! (their Table 1 value 9.9e-3 kg/m^3/K^2 converts to 9.9e-3/density_0 = 9.6e-6 K^-2).
    logical           :: do_cabbeling = .False.
    real(kind=WP)     :: cabbeling_Cb = 1.0e-5_WP   ! [K^-2], TEOS-10/GSW-range realistic value

    ! Optional thermobaric term for the neverworld2 branch of density_linear(), independent
    ! of/combinable with do_cabbeling. Following DINO's EOS (Roquet et al. 2015, Eq.6):
    ! rho = rho0 - (a0 + 0.5*Cb*Ta + Th*depth)*Ta, depth positive-downward [m] (NOT FESOM's
    ! own negative-below-surface Z -- DINO's Th is only physically sensible, i.e. increasing
    ! effective thermal expansion with depth, when multiplied by a positive-downward depth).
    ! Mechanistically distinct from cabbeling: depth-dependent rather than uniform, so it
    ! only strengthens sinking/already-dense water rather than also lightening the warm
    ! equatorial surface the way do_cabbeling's quadratic-in-T term does. thermobaric_Th is
    ! rho0-normalized [K^-1 m^-1], same convention as thermal_alpha/cabbeling_Cb -- default
    ! converts DINO/Roquet's raw Table 1 value (Th=2.47e-5 kg/m^3/K/m) via /density_0.
    logical           :: do_thermobar   = .False.
    real(kind=WP)     :: thermobaric_Th = 2.4e-8_WP   ! [K^-1 m^-1], DINO's 2.47e-5/1030

    ! Optional initial salinity stratification, off by default -- when .False. (the old
    ! behavior), initial_state_neverworld2 sets salinity flat/uniform (35 psu everywhere,
    ! horizontally and with depth), exactly as before. When .True., it instead uses a
    ! single-exponential-decay profile anchored at 35 psu at depth (same functional style
    ! as the temperature stratification below, but offset so it doesn't decay toward the
    ! unphysical zero-salinity limit): S(z) = 35 + Sstrat_amp*exp(-|Z|/1000). Dynamically
    ! inert unless do_haline is also switched on (see below).
    logical           :: do_Sstrat   = .False.
    real(kind=WP)     :: Sstrat_amp  = 0.5_WP    ! [psu], surface-minus-deep salinity anomaly

    ! Optional linear haline (salt) term in the neverworld2 branch of density_linear(),
    ! independent of/combinable with do_cabbeling. Off by default (beta=0), matching the
    ! original neverworld2 design (salinity dynamically inert regardless of do_Sstrat).
    ! When .True., adds density_0*haline_beta*(S-35) on top of whichever temperature term
    ! is active (plain linear or cabbeling), same normalized-beta convention already used
    ! for the dbgyre toy branch just above (rho = rho0 - rho0*alpha*(T-10) +
    ! rho0*beta*(S-35)) -- default value reuses dbgyre's existing beta=0.00079 as a
    ! same-codebase, already-sane order of magnitude rather than converting DINO's raw
    ! b0 (different units/convention). sw_alpha_beta() must mirror this with a matching
    ! sw_beta when do_haline=.True. (see there). Note: turning this on without also
    ! giving salinity real forcing/restoring (do_Sstrat is IC-only, no restoring exists
    ! yet) means whatever salinity structure exists will still just get mixed away over
    ! a long run, since nothing sustains the gradient against diffusion/advection.
    logical           :: do_haline    = .False.
    real(kind=WP)     :: haline_beta  = 0.00079_WP   ! [1/psu], normalized like sw_alpha (dbgyre's value)

    ! Optional surface salinity restoring, off by default. Provides the missing salt-
    ! restoring *target* that do_haline needs to have any real, sustained effect (do_Sstrat
    ! is IC-only and gets mixed away over a long run otherwise). Target follows DINO's Eq.
    ! B2 (Kamm et al. 2025, Table 1/Appendix B): constant in time (unlike Theta*, no
    ! seasonal cycle), same cos(0.5*pi*lat/Ly) shape as trelax_opt=2's t_base above, no
    ! cold-patch-style polar blend (DINO's S* doesn't have one either). Applied as a flux
    ! via relax_2_ssurf() -> relax_salt, reimplementing FESOM2's generic SSS-relaxation
    ! mechanism (ice_oce_coupling.F90) directly in this module -- that mechanism only ever
    ! runs inside the use_ice branch of the main timestep loop, which this toy config never
    ! takes (use_ice=.false.), the same reason trelax_opt=2 above has to set heat_flux_in
    ! itself rather than relying on the generic (ice-coupled) flux path.
    logical           :: do_Srelax  = .False.
    real(kind=WP)     :: s_north    = 35.00_WP   ! [psu], DINO S*_n
    real(kind=WP)     :: s_south    = 35.10_WP   ! [psu], DINO S*_s
    real(kind=WP)     :: s_equator  = 37.25_WP   ! [psu], DINO S*_eq

    ! Optional seasonal cycle in the trelax_opt=2 north/south cold-patch pole targets,
    ! off by default -- follows DINO's Eq. B3/B4 (Kamm et al. 2025, Appendix B): a plain
    ! annual cosine, evaluated from FESOM's own clock (daynew/timenew/ndpyr, g_clock)
    ! rescaled onto DINO's 1-360 day convention so their exact cos(pi*(d-201)/180) phasing
    ! applies unchanged regardless of FESOM's own (365-day, include_fleapyear=.false. for
    ! this config) calendar length. Composes with (adds on top of) north_pole_target/
    ! south_pole_target rather than replacing them with DINO's own absolute means (5degC
    ! north / -0.5degC south) -- so every existing sweep case's pole target stays its own
    ! tuned mean, now with an optional seasonal wobble on top, rather than introducing a
    ! second, DINO-specific mean-value convention alongside the existing one. Hemispheres
    ! are deliberately out of phase (north: +amp*cos(...), south: -amp*cos(...)), matching
    ! DINO's own sign convention (real NH-summer/SH-winter asymmetry, not two independent
    ! wobbles). Recomputed every relax_2_tsurf() call (cheap -- pure per-node function of
    ! latitude + day-of-year, no communication) rather than once at init time, since the
    ! whole point is that the target itself now varies through the run.
    logical           :: do_Tseasonal       = .False.
    real(kind=WP)     :: north_seasonal_amp = 3.0_WP   ! [degC], DINO Eq. B3 amplitude
    real(kind=WP)     :: south_seasonal_amp = 0.5_WP   ! [degC], DINO Eq. B4 amplitude

    ! Meridional half-extent of the domain [rad], computed once in initial_state_neverworld2
    ! (MPI reduction over mesh latitude) and reused here as a module-level constant -- needed
    ! by relax_2_tsurf()'s do_Tseasonal recompute, which cannot afford to redo that reduction
    ! every timestep.
    real(kind=WP), save :: Ly

    ! All toy_neverworld2 parameters above are namelist-controlled: read from the
    ! &oce_neverworld2 group in namelist.oce by gen_model_setup.F90 (only when
    ! which_toy=='neverworld2'), tolerant of the group being missing (falls back to
    ! the compiled-in defaults declared above).
    NAMELIST /oce_neverworld2/ do_wind, wind_opt, do_Tpert, do_Trelax, Trelax_opt, tau_inv, &
                                gamma_restore, &
                                do_north_cold_patch, north_pole_target, north_cold_patch_power, &
                                do_south_cold_patch, south_pole_target, south_cold_patch_power, &
                                equator_target, &
                                thermal_alpha, &
                                do_cabbeling, cabbeling_Cb, &
                                do_thermobar, thermobaric_Th, &
                                do_Sstrat, Sstrat_amp, &
                                do_haline, haline_beta, &
                                do_Srelax, s_north, s_south, s_equator, &
                                do_Tseasonal, north_seasonal_amp, south_seasonal_amp
    contains
    !
    !
    !_______________________________________________________________________________
    subroutine initial_state_neverworld2(dynamics, tracers, partit, mesh)
        
        implicit none
        type(t_mesh)  , intent(inout), target :: mesh
        type(t_partit), intent(inout), target :: partit
        type(t_tracer), intent(inout), target :: tracers
        type(t_dyn)   , intent(inout), target :: dynamics 

        integer                               :: elem, node, ii, nz, nzmin, nzmax, elnodes(3), taul, idx
        real(kind=WP)                         :: lon, lat, dlat_tau, loc_Ly
        ! Ly itself is now module-level (see declaration above) -- reused by relax_2_tsurf()'s
        ! do_Tseasonal recompute, so it must survive past this subroutine's return.
        real(kind=WP)                         :: ynorm, t_base, cold_weight   ! used by trelax_opt=2
        real(kind=WP), allocatable            :: lat_tau(:), val_tau(:)
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

        ! default values
        stress_surf   = 0.0_WP !make sure how it works (y-component is 0)
        heat_flux     = 0.0_WP
        water_flux    = 0.0_WP        
        relax2clim    = 0.0_WP
        ! Initial salinity: flat/uniform by default (old behavior, unchanged), or an
        ! optional single-exponential-decay stratification -- same functional style as
        ! the temperature profile below, but anchored at 35 psu at depth rather than
        ! decaying toward the unphysical zero-salinity limit -- see do_Sstrat.
        if (do_Sstrat) then
            do node = 1, myDim_nod2D+eDim_nod2D
                nzmin = ulevels_nod2D(node)
                nzmax = nlevels_nod2D(node)-1
                do nz = nzmin, nzmax
                    tracers%data(2)%values(nz,node) = 35.0_WP + Sstrat_amp*exp(-abs(Z(nz))/1000.0_WP)
                end do
            end do
        else
            tracers%data(2)%values(:,:) = 35.0_WP
        end if
        Tsurf         = tracers%data(1)%values(1,:)
        Ssurf         = tracers%data(2)%values(1,:)
        
        !___________________________________________________________________________
        ! determine latitudinal domain size from mesh
        loc_Ly=omp_min_max_sum1(coord_nod2D(2,:), 1, myDim_nod2D, 'max', partit) 
        call MPI_AllREDUCE(loc_Ly, Ly, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
        
        !___________________________________________________________________________
        ! Do windprofile input (momentum flux):
        ! switch on wind 
        if (do_wind) then 
            if (wind_opt==1) then 
                ! Option (1):
                ! wind forcing (momentum flux) from windstress.out file with interpolation  
                open(20, file=trim(meshpath)//'windstress.out', status='old')
                
                ! number of windprofile  points 
                read(20,*) taul
                
                ! read profile 
                allocate(lat_tau(taul), val_tau(taul))
                do ii = 1,taul
                    read(20,*) lat_tau(ii), val_tau(ii)
                end do
                close(20)
                
                ! lat resolution of wind stress profile 
                dlat_tau = (lat_tau(2)-lat_tau(1))
                
                ! linear interpolate windstress profile data to our elem locations
                do elem=1, myDim_elem2D
                    elnodes = elem2d_nodes(:,elem)
                    lat     = sum(coord_nod2D(2,elnodes))/3.0_WP
                    lat     = lat*180/pi    ! We return to degrees
                    idx     = int((lat+lat_tau(1))/dlat_tau)+1
                    stress_surf(1,elem) = val_tau(idx) + (val_tau(idx+1)-val_tau(idx))* (lat-lat_tau(idx))/dlat_tau
                end do
                deallocate(val_tau, lat_tau)
                
            ! Option (2) read already to elements interpoalted windprofile file 
            elseif (wind_opt == 2) then 
                allocate(val_tau(elem2d))
                open(20, file=trim(meshpath)//'windstress@elem.out', status='old')
                read(20, *) val_tau
                stress_surf(1,:)=val_tau(myList_elem2D)
                deallocate(val_tau)
                
            else
                write(*,*) " -ERROR-> This wind_opt is not supported !"
                call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
            end if 
        
        ! switch off wind 
        else 
            stress_surf = 0.0_WP
        
        end if
        
        !___________________________________________________________________________
        !  Initial temperature stratification
        do node = 1, myDim_nod2D+eDim_nod2D
            nzmin = ulevels_nod2D(node)
            nzmax = nlevels_nod2D(node)-1
            do nz = nzmin, nzmax
                tracers%data(1)%values(nz,node) = 28.0 * exp(-abs(Z(nz))/1000.0)
            end do
        end do
        ! In order to be consistent with Neverworld2, we need to take a linear equation 
        ! of state with temperature the only tracer. In the equation of state, alpha should be 
        ! 0.0002, or 0.2 when multiplied with rho_0   pho=0.2(T_ref-T), T_ref can be any.
        ! --> see oce_ale_pressure_bv.F90 --> subroutine density_linear(...)
        
        !___________________________________________________________________________
        ! Surface temperature relaxation target:
        ! The Neverworld2 is adiabatic, but I suspect some surface forcing will be
        ! needed for long runs. This is an example of some surface temperature
        ! distribution. Relaxation of surface temperature to this distribution will
        ! be introducing some temperature forcing.
        ! Some experiments are needed to adjust Tsurf
        if (do_Trelax) then
            select case(trelax_opt)
            case(1)
                ! original: symmetric cosine target, ~18-28C, lat in degrees
                do node = 1, myDim_nod2D+eDim_nod2D
                    lat          = coord_nod2D(2,node)/rad
                    Tsurf(node) = 18.0+10.0*cos(pi*lat/Ly)
                end do
            case(2)
                ! Genevieve: symmetric cosine target, ~5-28C, consistent radian units
                ! throughout (lat and Ly both in radians), plus optional asymmetric
                ! cold patch(es) that pull the target down to north_pole_target /
                ! south_pole_target near the northern (ynorm=lat/Ly -> 1) and/or
                ! southern (ynorm -> -1) edge, independently switchable.
                do node = 1, myDim_nod2D+eDim_nod2D
                    lat   = coord_nod2D(2,node)
                    ynorm = lat/Ly
                    t_base = 5.0_WP + (equator_target-5.0_WP)*cos(0.5_WP*pi*lat/Ly)
                    if (do_north_cold_patch .and. ynorm>0.0_WP) then
                        cold_weight = ynorm**north_cold_patch_power
                        Tsurf(node) = (1.0_WP-cold_weight)*t_base + cold_weight*north_pole_target
                    else if (do_south_cold_patch .and. ynorm<0.0_WP) then
                        cold_weight = abs(ynorm)**south_cold_patch_power
                        Tsurf(node) = (1.0_WP-cold_weight)*t_base + cold_weight*south_pole_target
                    else
                        Tsurf(node) = t_base
                    end if
                end do
            case default
                write(*,*) " -ERROR-> This trelax_opt is not supported !"
                call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
            end select
        end if
        ! In order to introduce forcing we need to allow surface relaxation to
        ! climatology.

        !___________________________________________________________________________
        ! Surface salinity relaxation target (DINO Eq. B2): constant-in-time, same
        ! cos(0.5*pi*lat/Ly) shape as t_base above, no cold-patch blend -- see do_Srelax
        ! declaration above for the full rationale. Actual restoring flux is applied
        ! every timestep by relax_2_ssurf(), not here.
        if (do_Srelax) then
            do node = 1, myDim_nod2D+eDim_nod2D
                lat = coord_nod2D(2,node)
                if (lat>=0.0_WP) then
                    Ssurf(node) = s_north + (s_equator-s_north)*cos(0.5_WP*pi*lat/Ly)
                else
                    Ssurf(node) = s_south + (s_equator-s_south)*cos(0.5_WP*pi*lat/Ly)
                end if
            end do
        end if

        !___________________________________________________________________________
        ! Temperature perturbation to trigger instabilities
        if (do_Tpert) then
            do node = 1, myDim_nod2D+eDim_nod2D
                ! inject surface perturbations (in the upper two layers) in the 
                ! middle of the ocean @[30°E, -50°S] in a
                ! distance radius of 5deg around tha point
                lat=coord_nod2D(2,node)+50.0*rad
                lon=coord_nod2D(1,node)-30.0*rad
                lat=sqrt(lat*lat+lon*lon)
                if (lat<=5*rad) then
                    tracers%data(1)%values(1:2,node)=tracers%data(1)%values(1:2,node)+0.5*cos(pi*lat/2.0/5.0/rad)
                end if
            end do 
        end if
        
    end subroutine initial_state_neverworld2



    !
    !
    !_______________________________________________________________________________
    subroutine relax_2_tsurf(tdata, partit, mesh)
        implicit none
        type(t_mesh),        intent(in),    target  :: mesh
        type(t_partit),      intent(inout), target  :: partit
        type(t_tracer_data), intent(inout), target  :: tdata
        integer                                     :: node
        real(kind=WP)                               :: t_surface, t_target   ! used by trelax_opt=2
        real(kind=WP)                               :: lat, ynorm, t_base, cold_weight  ! used by do_Tseasonal
        real(kind=WP)                               :: frac_year, d_dino, seasonal_phase
        real(kind=WP)                               :: north_pole_target_now, south_pole_target_now
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
        if (do_Trelax) then
            select case(trelax_opt)
            case(1)
                ! original: nudge the tracer directly towards Tsurf
                do node=1, myDim_nod2D
                    tdata%values(1,node)  = tdata%values(1,node)+dt*tau_inv*(Tsurf(node)-tdata%values(1,node))
                end do
            case(2)
                ! Genevieve: apply restoring as a surface heat flux instead of nudging
                ! the tracer directly, so it shows up in the heat budget. Overwrites
                ! heat_flux every call since, for this adiabatic toy config, restoring
                ! is its only source.
                ! Recompute Tsurf every call (not just once at init) so relax_2_tsurf is the
                ! single source of truth for the trelax_opt=2 target under BOTH do_Tseasonal
                ! states -- avoids keeping two separate copies of the same cold-patch formula
                ! (this one, and initial_state_neverworld2's) that could drift out of sync.
                ! initial_state_neverworld2's own version now only supplies a harmless initial
                ! value, immediately overwritten here on the first call.
                if (do_Tseasonal) then
                    ! DINO Eq. B3/B4 (Kamm et al. 2025): a plain annual cosine on top of
                    ! the existing (namelist-tuned) pole-target means, hemispheres out of
                    ! phase (real NH-summer/SH-winter asymmetry). Rescale FESOM's own
                    ! day-of-year (daynew/timenew, out of ndpyr) onto DINO's 1-360 day
                    ! convention so cos(pi*(d-201)/180) applies with its original phasing
                    ! regardless of FESOM's calendar length.
                    frac_year = (real(daynew-1,WP) + timenew/86400.0_WP) / real(ndpyr,WP)
                    d_dino    = 1.0_WP + frac_year*360.0_WP
                    seasonal_phase = cos(pi*(d_dino-201.0_WP)/180.0_WP)
                    north_pole_target_now = north_pole_target +  north_seasonal_amp*seasonal_phase
                    south_pole_target_now = south_pole_target -  south_seasonal_amp*seasonal_phase
                else
                    north_pole_target_now = north_pole_target
                    south_pole_target_now = south_pole_target
                end if
                do node = 1, myDim_nod2D+eDim_nod2D
                    lat   = coord_nod2D(2,node)
                    ynorm = lat/Ly
                    t_base = 5.0_WP + (equator_target-5.0_WP)*cos(0.5_WP*pi*lat/Ly)
                    if (do_north_cold_patch .and. ynorm>0.0_WP) then
                        cold_weight = ynorm**north_cold_patch_power
                        Tsurf(node) = (1.0_WP-cold_weight)*t_base + cold_weight*north_pole_target_now
                    else if (do_south_cold_patch .and. ynorm<0.0_WP) then
                        cold_weight = abs(ynorm)**south_cold_patch_power
                        Tsurf(node) = (1.0_WP-cold_weight)*t_base + cold_weight*south_pole_target_now
                    else
                        Tsurf(node) = t_base
                    end if
                end do
                heat_flux = 0.0_WP
                do node=1, myDim_nod2D
                    t_surface = tdata%values(1,node)
                    t_target  = Tsurf(node)
                    heat_flux(node) = gamma_restore*(t_surface-t_target)
                end do
                ! heat_flux_in is what diag_densMOC actually reads for its heat-driven
                ! buoyancy-flux term (dens_flux_e, std_heat_flux) -- it's normally set in
                ! oce_fluxes/ice_oce_coupling.F90, a path this toy config never runs
                ! (use_ice=.false.), so without this line heat_flux_in would stay 0 forever
                ! and the DMOC diagnostic would silently miss this restoring flux entirely.
                heat_flux_in = heat_flux
                ! The loop above only fills owned nodes (1:myDim_nod2D); halo nodes are
                ! still whatever heat_flux=0.0_WP reset them to. diag_densMOC reads these
                ! per-element (elnodes), and elements straddling a partition boundary have
                ! halo-node vertices -- without this exchange those would silently read 0
                ! instead of the neighbouring rank's value, biasing dens_flux_e/std_heat_flux
                ! low near partition boundaries.
                call exchange_nod(heat_flux, heat_flux_in, partit)
            end select
        end if

    end subroutine relax_2_tsurf



    !
    !
    !_______________________________________________________________________________
    subroutine relax_2_ssurf(tdata, partit, mesh)
        implicit none
        type(t_mesh),        intent(in),    target  :: mesh
        type(t_partit),      intent(inout), target  :: partit
        type(t_tracer_data), intent(inout), target  :: tdata
        integer                                     :: node
        real(kind=WP)                               :: net
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
        if (do_Srelax) then
            ! Apply salinity restoring as a flux into relax_salt, mirroring FESOM2's
            ! generic SSS-relaxation mechanism (ice_oce_coupling.F90:504) including its
            ! domain-mean subtraction (enforces zero net domain-integrated salt input at
            ! every timestep) -- per Kamm et al. 2025 Appendix B, this is exactly what
            ! makes DINO's own two-phase E-P-climatology spin-up trick unnecessary here.
            do node = 1, myDim_nod2D+eDim_nod2D
                relax_salt(node) = surf_relax_S*(Ssurf(node)-tdata%values(1,node))
            end do
            call integrate_nod(relax_salt, net, partit, mesh)
            net = net/ocean_area
            do node = 1, myDim_nod2D+eDim_nod2D
                relax_salt(node) = relax_salt(node)-net
            end do
        end if

    end subroutine relax_2_ssurf


    !
    !
    !_______________________________________________________________________________
    ! simple toy mixing only increase vertical mixing where buoyancy is negativ !!!
    subroutine oce_mixing_TOY(partit, mesh)
        use MOD_MESH
        use MOD_PARTIT
        use MOD_PARSUP
        use o_PARAM
        use o_ARRAYS
        use g_config
        implicit none

        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(inout), target :: partit
        integer                               :: node, nz, nzmax, nzmin, elem, elnodes(3)    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
        !
        !___________________________________________________________________________
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(node, nz, nzmin, nzmax)
        do node=1, myDim_nod2D+eDim_nod2D
            nzmax = nlevels_nod2d(node)
            nzmin = ulevels_nod2d(node)
            do nz=nzmin+1, nzmax-1
                ! set vertical background diffusivity
                Kv(nz,node) = K_ver
                
                ! force on convection if unstable 
                if (bvfreq(nz, node) < 0._WP)  Kv(nz,node)=max(Kv(nz,node),instabmix_kv)
            end do
        end do
    !$OMP END PARALLEL DO
        
        !
        !___________________________________________________________________________
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem, elnodes, nz, nzmin, nzmax)
        do elem=1, myDim_elem2D
            elnodes=elem2D_nodes(:,elem)
            nzmax = nlevels(elem)
            nzmin = ulevels(elem)
            do nz=nzmin+1,nzmax-1
                 ! set vertical background viscosity
                Av(nz,elem) = A_ver
                
                ! force on convection if unstable 
                if (any(bvfreq(nz, elnodes) < 0._WP)) Av(nz,elem)=max(Av(nz,elem),instabmix_kv)
            end do
        end do
    !$OMP END PARALLEL DO    
    end subroutine oce_mixing_TOY

END MODULE Toy_Neverworld2
