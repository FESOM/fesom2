module ocean2ice_interface
    interface
        subroutine ocean2ice(ice, dynamics, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_DYN
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_dyn)   , intent(in)   , target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module oce_fluxes_interface
    interface
        subroutine oce_fluxes(ice, dynamics, tracers, partit, mesh, total_nsteps, nstep)
        USE MOD_ICE
        USE MOD_DYN
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_dyn)   , intent(in)   , target :: dynamics
        type(t_tracer), intent(inout), target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        integer       , intent(in)            :: total_nsteps
        integer       , intent(in)            :: nstep
        end subroutine
        
        subroutine oce_fluxes_mom(ice, dynamics, partit, mesh)
        USE MOD_ICE
        USE MOD_DYN
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_dyn)   , intent(in)   , target :: dynamics
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine fw_surf_anomaly(hSv, tracers, partit, mesh)
        use o_PARAM
        use o_ARRAYS
        USE MOD_TRACER
        USE MOD_PARTIT
        use MOD_MESH
        use MOD_PARSUP
        USE g_CONFIG
        use g_comm_auto
        use g_support
        type(t_tracer), intent(in),    target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        real(kind=WP) , intent(in)            :: hSv
        end subroutine

    end interface
end module

!
!
!_______________________________________________________________________________
! transmits the relevant fields from the ice to the ocean model
subroutine oce_fluxes_mom(ice, dynamics, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use o_ARRAYS
    USE g_CONFIG
    use g_comm_auto
    use cavity_interfaces    
#if defined (__icepack)
    use icedrv_main,   only: icepack_to_fesom
#endif
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                  :: n, elem, elnodes(3),n1
    real(kind=WP)            :: aux
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice, a_ice, u_w, v_w
    real(kind=WP), dimension(:), pointer  :: stress_iceoce_x, stress_iceoce_y  
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice(:)
    v_ice           => ice%vice(:)
    a_ice           => ice%data(1)%values(:)
    u_w             => ice%srfoce_u(:)
    v_w             => ice%srfoce_v(:)
    stress_iceoce_x => ice%stress_iceoce_x(:)
    stress_iceoce_y => ice%stress_iceoce_y(:)
  
    ! ==================
    ! momentum flux:
    ! ==================
    !___________________________________________________________________________

#if defined (__icepack)
     call icepack_to_fesom(nx_in=(myDim_nod2D+eDim_nod2D), &
                           aice_out=a_ice)
#endif
    !___________________________________________________________________________
    ! compute total surface stress (iceoce+atmoce) on nodes 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, elem, elnodes, n1, aux)
!$OMP DO
    do n=1,myDim_nod2D+eDim_nod2D   
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(n)>1) cycle
        
        !_______________________________________________________________________
        if(a_ice(n)>0.001_WP) then
            aux=sqrt((u_ice(n)-u_w(n))**2+(v_ice(n)-v_w(n))**2)*density_0*ice%cd_oce_ice
            stress_iceoce_x(n) = aux * (u_ice(n)-u_w(n))
            stress_iceoce_y(n) = aux * (v_ice(n)-v_w(n))
        else
            stress_iceoce_x(n)=0.0_WP
            stress_iceoce_y(n)=0.0_WP
        end if
        
        stress_node_surf(1,n) = stress_iceoce_x(n)*a_ice(n) + stress_atmoce_x(n)*(1.0_WP-a_ice(n))
        stress_node_surf(2,n) = stress_iceoce_y(n)*a_ice(n) + stress_atmoce_y(n)*(1.0_WP-a_ice(n))
    end do
!$OMP END DO
    !___________________________________________________________________________
    ! compute total surface stress (iceoce+atmoce) on elements
!$OMP DO
    DO elem=1,myDim_elem2D
        !_______________________________________________________________________
        ! if cavity element skip it 
        if (ulevels(elem)>1) cycle
        
        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        stress_surf(1,elem)=sum(stress_iceoce_x(elnodes)*a_ice(elnodes) + &
                                stress_atmoce_x(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
        stress_surf(2,elem)=sum(stress_iceoce_y(elnodes)*a_ice(elnodes) + &
                                stress_atmoce_y(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
    END DO
!$OMP END DO
!$OMP END PARALLEL
    !___________________________________________________________________________
    if (use_cavity) call cavity_momentum_fluxes(dynamics, partit, mesh)
  
end subroutine oce_fluxes_mom
!
!
!_______________________________________________________________________________
! transmits the relevant fields from the ocean to the ice model
subroutine ocean2ice(ice, dynamics, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    integer :: n, elem, k
    real(kind=WP) :: uw, vw, vol
    real(kind=WP), dimension(:,:)  , pointer :: temp, salt
    real(kind=WP), dimension(:,:,:), pointer :: UV
    real(kind=WP), dimension(:)    , pointer :: S_oc_array, T_oc_array, u_w, v_w, elevation
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    temp       => tracers%data(1)%values(:,:)
    salt       => tracers%data(2)%values(:,:)
    UV         => dynamics%uv(:,:,:)
    u_w        => ice%srfoce_u(:)
    v_w        => ice%srfoce_v(:)
    T_oc_array => ice%srfoce_temp(:)
    S_oc_array => ice%srfoce_salt(:)
    elevation  => ice%srfoce_ssh(:)
    
    !___________________________________________________________________________
    ! the arrays in the ice model are renamed
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, elem, k, uw, vw, vol)
     if (ice%ice_update) then
!$OMP DO
        do n=1, myDim_nod2d+eDim_nod2d  
            if (ulevels_nod2D(n)>1) cycle 
            T_oc_array(n) = temp(1,n)
            S_oc_array(n) = salt(1,n)
            elevation(n)  = hbar(n)
        end do
!$OMP END DO
    else
!$OMP DO
        do n=1, myDim_nod2d+eDim_nod2d
            if (ulevels_nod2D(n)>1) cycle 
             T_oc_array(n) = (T_oc_array(n)*real(ice%ice_steps_since_upd,WP)+temp(1,n))/real(ice%ice_steps_since_upd+1,WP)
             S_oc_array(n) = (S_oc_array(n)*real(ice%ice_steps_since_upd,WP)+salt(1,n))/real(ice%ice_steps_since_upd+1,WP)
             elevation(n)  = (elevation(n) *real(ice%ice_steps_since_upd,WP)+  hbar(n))/real(ice%ice_steps_since_upd+1,WP)
        end do
!$OMP END DO
    end if

!$OMP DO
    do n=1, myDim_nod2d+eDim_nod2d
       u_w(n) = 0.0_WP
       v_w(n) = 0.0_WP
    end do
!$OMP END DO

!$OMP DO
    do n=1, myDim_nod2d  
        if (ulevels_nod2D(n)>1) cycle 
        uw  = 0.0_WP
        vw  = 0.0_WP
        vol = 0.0_WP
        do k=1, nod_in_elem2D_num(n)
            elem=nod_in_elem2D(k,n)
            if (ulevels(elem)>1) cycle
            vol = vol + elem_area(elem)
            uw  = uw+ UV(1,1,elem)*elem_area(elem)
            vw  = vw+ UV(2,1,elem)*elem_area(elem)
        end do
        uw = uw/vol
        vw = vw/vol
        
        if (ice%ice_update) then
            u_w(n)=uw
            v_w(n)=vw
        else
            u_w(n)=(u_w(n)*real(ice%ice_steps_since_upd,WP)+uw)/real(ice%ice_steps_since_upd+1,WP)
            v_w(n)=(v_w(n)*real(ice%ice_steps_since_upd,WP)+vw)/real(ice%ice_steps_since_upd+1,WP)
        endif
    end do
!$OMP END DO
!$OMP END PARALLEL
    call exchange_nod(u_w, v_w, partit)
end subroutine ocean2ice
!
!
!_______________________________________________________________________________
subroutine oce_fluxes(ice, dynamics, tracers, partit, mesh, total_nsteps, nstep)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_CONFIG
    use o_ARRAYS
    use g_comm_auto
    use g_forcing_param, only: use_virt_salt
    use g_forcing_arrays
    use g_support
    use cavity_interfaces
#if defined (__icepack)
    use icedrv_main,   only: icepack_to_fesom,    &
                            init_flux_atm_ocn
#endif
    use cavity_interfaces
    use g_clock
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    integer       , intent(in)            :: total_nsteps
    integer       , intent(in)            :: nstep
    !___________________________________________________________________________
    integer                    :: n, elem, elnodes(3),n1
    real(kind=WP)              :: rsss, net, hSv, fw0
    real(kind=WP), allocatable :: flux(:)
    !___________________________________________________________________________
    real(kind=WP), dimension(:,:), pointer :: temp, salt
    real(kind=WP), dimension(:)  , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)  , pointer :: a_ice_old
    real(kind=WP), dimension(:)  , pointer :: thdgr, thdgrsn
    real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
    real(kind=WP)                , pointer :: rhoice, rhosno, inv_rhowat
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    temp          => tracers%data(1)%values(:,:)
    salt          => tracers%data(2)%values(:,:)
    a_ice         => ice%data(1)%values(:)
    m_ice         => ice%data(2)%values(:)
    m_snow        => ice%data(3)%values(:)
    a_ice_old     => ice%data(1)%values_old(:)
    thdgr         => ice%thermo%thdgr(:)
    thdgrsn       => ice%thermo%thdgrsn(:)
    fresh_wa_flux => ice%flx_fw(:)
    net_heat_flux => ice%flx_h(:)
    rhoice        => ice%thermo%rhoice
    rhosno        => ice%thermo%rhosno
    inv_rhowat    => ice%thermo%inv_rhowat
    
    !___________________________________________________________________________
    allocate(flux(myDim_nod2D+eDim_nod2D))

!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2d  
       flux(n) = 0.0_WP
    end do
!$OMP END PARALLEL DO

    ! ==================
    ! heat and freshwater
    ! ==================   
    !___________________________________________________________________________
    ! from here on: 
    !    (-)  (+)
    !     |    ^
    ! ~~~~|~~~~|~~~~
    !     V    |
    !     
#if defined (__icepack)

    call icepack_to_fesom (nx_in=(myDim_nod2D+eDim_nod2D), &
                           aice_out=a_ice,                 &
                           vice_out=m_ice,                 &
                           vsno_out=m_snow,                &
                           fhocn_tot_out=net_heat_flux,    &
                           fresh_tot_out=fresh_wa_flux,    &
                           fsalt_out=real_salt_flux,       &
                           dhi_dt_out=thdgrsn,             &
                           dhs_dt_out=thdgr,               &
                           evap_ocn_out=evaporation        )

    heat_flux(:)   = - net_heat_flux(:)
    water_flux(:)  = - (fresh_wa_flux(:)/1000.0_WP) - runoff(:)

    ! Evaporation
    evaporation(:) = - evaporation(:) / 1000.0_WP
    ice_sublimation(:) = 0.0_WP

    call init_flux_atm_ocn()

#else
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2d  
       heat_flux(n)   = -net_heat_flux(n)
       water_flux(n)  = -fresh_wa_flux(n)
    end do
!$OMP END PARALLEL DO
#endif
    
    !___freshwater hosing routine_______________________________________________
    !
    !1992-2020 hosing experiment increasing 1.1*10^-3 Sv/yr
    if (yearnew<1992) then
        fw0 = 0.0
        hSv = 0.0 
    else
        fw0 = 0.0011*(yearnew-1992)
        hSv = fw0 + (0.0011/total_nsteps)*nstep
    endif
    !
    !constant flux
    !hSv=0.1
    !
    call fw_surf_anomaly(hSv, tracers, partit, mesh)

    !___________________________________________________________________________
    ! add heat and fresh water flux from cavity 
    if (use_cavity) then
        call cavity_heat_water_fluxes_3eq(ice, dynamics, tracers, partit, mesh)   
        call exchange_nod(heat_flux, water_flux, partit)
!$OMP BARRIER
    end if 
    
    !___________________________________________________________________________
    ! save total heat flux (heat_flux_in) since heat_flux will be alternated by 
    ! sw_pene
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2d  
       heat_flux_in(n)=heat_flux(n) ! sw_pene will change the heat_flux
    end do
!$OMP END PARALLEL DO
    
    !___________________________________________________________________________
    ! on freshwater inflow/outflow or virtual salinity:
    ! 1. In zlevel & zstar the freshwater flux is applied in the update of the 
    ! ssh matrix when solving the continuity equation of vertically 
    ! integrated flow. The salt concentration in the first layer will 
    ! be then adjusted according to the change in volume.
    ! In this case rsss is forced to be zero by setting ref_sss=0. and ref_sss_local=.false.
    ! in routines above.
    ! 2. In cases where the volume of the upper layer is fixed (i.e. linfs)  the freshwater flux 
    ! 'rsss*water_flux(n)' is applied as a virtual salt boundary condition via the vertical 
    ! diffusion operator.
    ! --> rsss*water_flux(n) : virtual salt flux 
    !___________________________________________________________________________
    ! balance virtual salt flux
    if (use_virt_salt) then ! will remain zero otherwise
        rsss=ref_sss
!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
            if (ref_sss_local) rsss = salt(ulevels_nod2d(n),n)
            virtual_salt(n)=rsss*water_flux(n) 
        end do
!$OMP END PARALLEL DO        

        call integrate_nod(virtual_salt, net, partit, mesh)
        
        ! we try not to change the virtual_salt flux values under the cavity from
        ! the balancing --> but we balance the contribution under the cavity over 
        ! the rest of the ocean
        ! ocean_area in case of cavity contains only the open ocean area
        net = net/ocean_area
!$OMP PARALLEL DO 
        do n=1, myDim_nod2D+eDim_nod2D
            if (ulevels_nod2d(n) > 1) cycle ! --> is cavity node 
            virtual_salt(n)=virtual_salt(n)-net
        end do
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    ! do virtual salt flux under the cavity also when zstar is switched on 
    ! since under the cavity nothing is allowed to move  --> linfs --> and linfs
    ! requires virtual salt
    elseif ( (.not. use_virt_salt) .and. (use_cavity) ) then ! will remain zero otherwise
        ! Introducing here in zstar a virtual salt flux in the cavity , will mess 
        ! up our global salinity conservation, we also cant do our usual global 
        ! virtual salt balancing since this will mess our local virtual_salt 
        ! flux, since we have no counter salt fluxes anywhere else in the ocean. 
        ! So what we try is to introduce an artificial very small counter virtual 
        ! saltflux at every open ocean vertice, to counter balance the virtual 
        ! salt flux in the cavity and to conserve the global salt budget   
        rsss=ref_sss
        
        ! compute virtual salt flux within the cavity
!$OMP PARALLEL DO         
        do n=1, myDim_nod2D+eDim_nod2D
            virtual_salt(n)=0.0_WP
            if (ulevels_nod2d(n) == 1) cycle ! --> is open ocean node 
            if (ref_sss_local) rsss = salt(ulevels_nod2d(n),n)
            virtual_salt(n)=rsss*water_flux(n) 
        end do
!$OMP END PARALLEL DO        
        
        ! integrate salt flux in the cavity(outside cavity virtual_salt is 0.0)
        call integrate_nod(virtual_salt, net, partit, mesh)
        
        ! counter balance the integrated cavity salt flux only in the open ocean 
        ! --> ensure global salt conservation !!!
        net = net/ocean_area
!$OMP PARALLEL DO         
        do n=1, myDim_nod2D+eDim_nod2D
            if (ulevels_nod2d(n) > 1) cycle ! --> is cavity node 
            virtual_salt(n)=virtual_salt(n)-net
        end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    ! balance SSS restoring to climatology
    if (use_cavity) then
!$OMP PARALLEL DO    
        do n=1, myDim_nod2D+eDim_nod2D
            relax_salt(n) = 0.0_WP
            if (ulevels_nod2d(n) > 1) cycle ! --> is cavity node --> only do salt relaxation in open ocean
            relax_salt(n)=surf_relax_S*(Ssurf(n)-salt(ulevels_nod2d(n),n))
        end do
!$OMP END PARALLEL DO
    else
!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
            relax_salt(n)=surf_relax_S*(Ssurf(n)-salt(ulevels_nod2d(n),n))
        end do
!$OMP END PARALLEL DO
    end if 
    
    ! --> if use_cavity=.true. relax_salt anyway zero where is cavity see above
    call integrate_nod(relax_salt, net, partit, mesh)
    net = net/ocean_area
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
        !--> only balance salt_relaxation in open ocean under the cavity it remains zero
        if (ulevels_nod2d(n) > 1) cycle ! --> is cavity node
        relax_salt(n)=relax_salt(n)-net
    end do
!$OMP END PARALLEL DO
    
    !___________________________________________________________________________
    ! enforce the total freshwater/salt flux be zero
    ! 1. water flux ! if (.not. use_virt_salt) can be used!
    ! we conserve only the fluxes from the database plus evaporation.
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
        flux(n) = evaporation(n)                      &
                  -ice_sublimation(n)                 & ! the ice2atmos subplimation does not contribute to the freshwater flux into the ocean
                  +prec_rain(n)                       &
                  +prec_snow(n)*(1.0_WP-a_ice_old(n)) &
#if defined (__oasis) || defined (__ifsinterface)                  
                  +residualifwflx(n)                  & ! balance residual ice flux only in coupled case
#endif
                  +runoff(n)
#if defined (__oasis) || defined (__ifsinterface)
! in the coupled mode the computation of freshwater flux takes into account the ratio between freshwater & salt water
        flux(n) = flux(n)*ice%thermo%rhofwt/ice%thermo%rhowat
#endif                  
    end do
!$OMP END PARALLEL DO
    ! --> In case of zlevel and zstar and levitating sea ice, sea ice is just sitting 
    ! on top of the ocean without displacement of water, there the thermodynamic 
    ! growth rates of sea ice have to be taken into account to preserve the fresh water 
    ! flux. In the case of floating sea ice, water is displaced by 
    ! sea ice and flux conservation from ocean sea ice transformation follows from 
    ! the conservation of volume
    ! --> In case of linfs ocean sea ice transformation is balanced by the virtual 
    ! salinity flux
    !!PS   if ( .not. use_floatice .and. .not. use_virt_salt) then
    if (.not. use_virt_salt) then
!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
           flux(n) = flux(n)-thdgr(n)*rhoice*inv_rhowat-thdgrsn(n)*rhosno*inv_rhowat
        end do
!$OMP END PARALLEL DO
    end if     
    
    !___________________________________________________________________________
    if (use_cavity) then
        ! with zstar we do not balance the freshwater flux under the cavity since its
        ! not contributing to the ocean volume increase/decrease since under the
        ! cavity ist linfs 
        if (.not. use_virt_salt) then !zstar, zlevel
            where (ulevels_nod2d > 1) flux = 0.0_WP
            
        ! with linfs we balance the freshwater flux globally, thereby ot changeing
        ! the freshwater flux under the cavity
        else ! linfs 
            where (ulevels_nod2d > 1) flux = -water_flux
        end if 
    end if 
    
    !___________________________________________________________________________
    ! compute total global net freshwater flux into the ocean 
    call integrate_nod(flux, net, partit, mesh)
    
    !___________________________________________________________________________
    ! here the + sign must be used because we switched up the sign of the 
    ! water_flux with water_flux = -fresh_wa_flux, but evap, prec_... and runoff still
    ! have there original sign
    ! if use_cavity=.false. --> ocean_area == ocean_areawithcav
    !! water_flux=water_flux+net/ocean_area
    if (use_cavity) then
        ! due to rigid lid approximation under the cavity we to not add freshwater
        ! under the cavity for the freshwater balancing we do this only for the open
        ! ocean
        net = net/ocean_area
!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
            if (ulevels_nod2d(n) > 1) cycle
            water_flux(n)=water_flux(n)+net
        end do
!$OMP END PARALLEL DO
    else
        net = net/ocean_area
!$OMP PARALLEL DO
        do n=1, myDim_nod2D+eDim_nod2D
           water_flux(n)=water_flux(n)+net
        end do
!$OMP END PARALLEL DO
    end if 
    
    !___________________________________________________________________________
    ! use the balanced water_flux and relax_salt flux (same as in the tracer 
    ! boundary condition) to compute the dens_flux for MOC diagnostic
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D    
        if (ulevels_nod2d(n) == 1) then ! --> is open ocean node 
            dens_flux(n)=sw_alpha(1,n) * heat_flux_in(n) / vcpw + sw_beta(1, n) * (relax_salt(n) + water_flux(n) * salt(1,n))
        else
            dens_flux(n)=0.0_WP
        end if
    end do
!$OMP END PARALLEL DO

    !___________________________________________________________________________
    if (use_sw_pene) call cal_shortwave_rad(ice, partit, mesh)
    
    !___________________________________________________________________________
    deallocate(flux)
    
end subroutine oce_fluxes
!
!
!_______________________________________________________________________________
subroutine fw_surf_anomaly(hSv, tracers, partit, mesh) 

    use o_PARAM
    use o_ARRAYS
    USE MOD_TRACER
    USE MOD_PARTIT
    use MOD_MESH
    use MOD_PARSUP
    USE g_CONFIG
    use g_comm_auto
    use g_support
    implicit none
    type(t_partit), intent(inout), target  :: partit
    type(t_mesh),   intent(in),    target  :: mesh
    type(t_tracer), intent(in),    target  :: tracers
    real(kind=WP),  intent(in)             :: hSv
    real(kind=WP)                          :: x, y, net
    integer                                :: n, ed(2)
    real(kind=WP), dimension(:,:), pointer :: temp
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
      temp          => tracers%data(1)%values(:,:)

      hosing_flux=0._WP
      hosing_heat_flux=0._WP

     !all nodes south of 60S 
     !do n=1, myDim_nod2D+eDim_nod2D
     !   y = geo_coord_nod2D(2,n)/rad
     !   if (y<-60._WP) hosing_flux(n)=1._WP
     !end do

     !only coastal nodes south of 60S
     do n=1, myDim_edge2D
        ed=mesh%edges(:, n)
        if (myList_edge2D(n) <= mesh%edge2D_in) cycle
        y = sum(geo_coord_nod2D(2,ed))/2._WP/rad
        if (y<-60._WP) hosing_flux(ed)=1._WP
     end do

     !call smooth_nod (hosing_flux, 8, partit, mesh)
     !do n=1,myDim_nod2d+eDim_nod2D
     !   if (hosing_flux(n)>0.0_WP) hosing_flux(n)=1.0_WP
     !end do

     call integrate_nod(hosing_flux, net, partit, mesh)

     if (abs(net)>1.e-6) then
        hosing_flux=hosing_flux/net*hSv*1.e6 ! hSv*1.e6 in m/s
     end if
     water_flux=water_flux-hosing_flux

     hosing_heat_flux=(temp(1,:)+1.8)*hosing_flux*vcpw !(tf - temp(1,:)) with tf=-1.8 
     !write(*,*) mype, 'heat_flux=', heat_flux, 'temp=', temp(1,:)
     heat_flux=heat_flux+hosing_heat_flux
     !write(*,*) mype, 'hosing_heat_flux=', hosing_heat_flux, 'heat_flux=', heat_flux
!write(*,*) mype, 'hSv=', hSv, 'net=', net, 'hf=', minval(hosing_flux, 1), maxval(hosing_flux, 1)

end subroutine fw_surf_anomaly
!
!
!_______________________________________________________________________________
subroutine fw_depth_anomaly(tts, ttt, hSv, partit, mesh)
    use o_PARAM
    use o_ARRAYS
    USE MOD_PARTIT
    use MOD_MESH
    use MOD_PARSUP
    use g_comm_auto
    use g_support
    use o_tracers
    use g_config,         only: dt
    implicit none

    integer                               :: n, ed(2), row, k, nz, nzmin, nzmax, elem1, elem2
    real(kind=WP)                         :: x, y, net, spar(100)
    real(kind=WP),  intent(in)            :: hSv
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP),  intent (inout)        :: tts(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP),  intent (inout)        :: ttt(mesh%nl-1,partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), save, allocatable      :: mask(:), mask3D(:,:)
    logical, save                         :: lfirst=.true.
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    hosing_flux3D=0._WP
    if (lfirst) then
       allocate(mask(myDim_nod2D+eDim_nod2D), mask3D(mesh%nl-1, myDim_nod2D+eDim_nod2D))
       mask =0._WP
       mask3D=0._WP
       spar=0._WP
       do n=1, myDim_edge2D
          ed=mesh%edges(:, n)
          if (myList_edge2D(n) <= mesh%edge2D_in) cycle
             y = sum(geo_coord_nod2D(2,ed))/2._WP/rad
          if (y<-60._WP) mask(ed)=1._WP
       end do
       
       call smooth_nod (mask, 5, partit, mesh)

       where (mask>1.e-12)
             mask=1.0_WP
       end where

       do n=1, myDim_edge2D
          ed = mesh%edges(:, n)
          if (any(mask(ed)<1.e-12)) cycle
          elem1 = mesh%edge_tri(1, n)
          elem2 = mesh%edge_tri(2, n)
          if (elem2<=0) then
             nzmin=1
             nzmax=mesh%nlevels(elem1)-1
          else
             nzmin = min(mesh%nlevels(elem1), mesh%nlevels(elem2))-1 !!remove -1 depending on what is desired
             nzmax = max(mesh%nlevels(elem1), mesh%nlevels(elem2))-1
          end if
          if (nzmax<nzmin) nzmax=nzmin
          nzmin=max(16, nzmin) !circa 200m
          nzmax=min(19, nzmax) !circa 400m
          mask3D(nzmin:nzmax, ed)=1._WP
       end do
       
!      mask=0.0_WP
       do row=1,myDim_nod2d+eDim_nod2D
          if (sum(mask3D(:, row))<1.e-12) mask(row)=0._WP
       end do
       
       call integrate_nod(mask, net, partit, mesh)

       if (abs(net)>1.e-6) then
          mask=mask/net*1.e6 ! hSv will be applied later as it changes in time
       end if
       
       do row=1,myDim_nod2d+eDim_nod2D 
          spar(:mesh%nl)=0._WP
          do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
             spar(k)=areasvol(k,row)*hnode(k,row)*mask3D(k, row) ! *any profile function(z) or 1.!                           
          end do
          !spar(:15)=0._WP
          !spar(16:19)=1._WP
          !spar(20:)=0._WP
          spar(:nlevels_nod2D(row)-1)=spar(:nlevels_nod2D(row)-1)/max(sum(spar(:nlevels_nod2D(row)-1)), 1.e-12)
          do k=ulevels_nod2D(row), nlevels_nod2D(row)-1
             spar(k)=spar(k)/areasvol(k,row)/hnode(k,row)*mask(row)*area(1,row)*dt
          end do
          !if ((mask(row)>1.e-12) .AND. (nlevels_nod2D(row)-1<=15)) then
          !    write(*,*) 'ALARM: (mask(row)>1.e-12) .AND. (nlevels_nod2D(row)-1<15)', nlevels_nod2D(row)
          !end if
          mask3D(:nlevels_nod2D(row)-1,row)=spar(:nlevels_nod2D(row)-1)
       end do
       lfirst=.false.
    !call integrate_nod(mask3D, net, partit, mesh)
    !write(*,*) 'integrated mask3D=', net
    !call integrate_nod(mask, net, partit, mesh)
    !write(*,*) 'integrated mask2D=', net
    end if
    tts=tts-mask3D*hSv*tts     
    ttt=ttt-mask3D*hSv*(ttt+1.8)
    call integrate_nod(tts, net, partit, mesh)
    write(*,*) 'integrated salinity=', net
    call integrate_nod(ttt, net, partit, mesh)
    write(*,*) 'integrated temperature=', net

    hosing_flux3D=hosing_flux3D+mask3D*hSv

!    do row=1,myDim_nod2d+eDim_nod2D
!       do k=1, mesh%nl-1
!          if (mask3(k, row)<1.e-12) CYCLE
!          ttf(k,row)=ttf(k,row)-mask3(k,row)*hSv*ttf(1, row)
!       end do
!    end do

!___ original routine___
!    hosing_flux=0._WP
!
!    do n=1, myDim_edge2D
!       ed=mesh%edges(:, n)
!       if (myList_edge2D(n) <= mesh%edge2D_in) cycle
!          y = sum(geo_coord_nod2D(2,ed))/2._WP/rad
!       if (y<-60._WP) hosing_flux(ed)=1._WP
!    end do
!
!    call smooth_nod (hosing_flux, 5, partit, mesh)
!    do row=1,myDim_nod2d+eDim_nod2D
!      if (hosing_flux(row)>0.0_WP) hosing_flux(row)=1.0_WP
!    end do
!
!    call integrate_nod(hosing_flux, net, partit, mesh)
!
!    if (abs(net)>1.e-6) then
!       hosing_flux=hosing_flux/net*hSv*1.e6 ! hSv*1.e6 in m/s
!    end if
!
!    do row=1, myDim_nod2d+eDim_nod2D
!        hosing_flux(row)=hosing_flux(row)*ttf(1, row)
!    end do
!
!    do row=1,myDim_nod2d+eDim_nod2D   ! myDim is sufficient
!     hosing_flux(row)=hosing_flux(row)*area(1,row)
!     if (ulevels_nod2D(row)>1) cycle
!        nzmin = ulevels_nod2D(row)
!        nzmax = min(21, nlevels_nod2D(row)-1)
!        do k=nzmin, nzmax
!          spar(k)=areasvol(k,row)*hnode(k,row) ! *any profile function(z)
!        end do
!        spar(nzmin:nzmax)=spar(nzmin:nzmax)/sum(spar(nzmin:nzmax))
!        do k=nzmin+1,nzmax
!           ttf(k,row)=ttf(k,row)-hosing_flux(row)*ttf(1, row)*spar(k)/areasvol(k,row)/hnode(k,row)*dt
!        end do
!    end do

end subroutine fw_depth_anomaly
!
!
!_______________________________________________________________________________
