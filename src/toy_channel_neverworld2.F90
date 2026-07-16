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
    implicit none
    SAVE 
    private
    public            :: initial_state_neverworld2, relax_2_tsurf, oce_mixing_TOY,  &
                        do_wind, wind_opt, tau_inv, do_Trelax, do_Tpert, trelax_opt, &
                        gamma_restore, do_north_cold_patch, north_pole_target, north_cold_patch_power, &
                        neverworld2_forcing
    logical           :: do_wind  = .True.    ! apply surface windstress
    integer           :: wind_opt = 2         ! 1: interpolate tau from profile data, 2: read already to elem interp profile data

    logical           :: do_Trelax= .False.   ! apply surface temp relaxation
    logical           :: do_Tpert = .True.    ! apply temp. perturbation to trigger instabilities

    ! Two alternative SST restoring implementations, picked by trelax_opt (same pattern as wind_opt):
    !   1: original -- directly nudges the temperature tracer: T += dt*tau_inv*(Ttarget-T)
    !   2: Genevieve's -- applies restoring as a surface heat flux instead, so it is visible
    !      in the heat budget: heat_flux = gamma_restore*(Tsurface-Ttarget) [W/m^2/K].
    !      FESOM heat_flux convention: >0 cools the ocean (heat leaving), <0 warms it.
    integer           :: trelax_opt = 1
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

    ! Optional asymmetric cold patch near the northern edge, only used with trelax_opt=2:
    ! blends the symmetric cosine restoring target towards north_pole_target as
    ! ynorm=lat/Ly -> 1. Off by default.
    logical           :: do_north_cold_patch    = .False.
    real(kind=WP)     :: north_pole_target      = 0.0_WP
    real(kind=WP)     :: north_cold_patch_power = 2.0_WP

    ! trelax_opt and gamma_restore are namelist-controlled: read from the &neverworld2_forcing
    ! group in namelist.oce by gen_model_setup.F90 (only when which_toy=='neverworld2').
    ! The other toy parameters above stay compiled-in defaults, same as before.
    NAMELIST /neverworld2_forcing/ trelax_opt, gamma_restore
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
        real(kind=WP)                         :: lon, lat, dlat_tau, loc_Ly, Ly
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
        tracers%data(2)%values(:,:) = 35.0_WP
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
                ! cold patch that pulls the target down to north_pole_target near the
                ! northern edge (ynorm=lat/Ly -> 1).
                do node = 1, myDim_nod2D+eDim_nod2D
                    lat   = coord_nod2D(2,node)
                    ynorm = lat/Ly
                    t_base = 5.0_WP + 28.0_WP*cos(0.5_WP*pi*lat/Ly)
                    if (do_north_cold_patch .and. ynorm>0.0_WP) then
                        cold_weight = ynorm**north_cold_patch_power
                        Tsurf(node) = (1.0_WP-cold_weight)*t_base + cold_weight*north_pole_target
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
