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
        subroutine oce_fluxes(ice, dynamics, tracers, partit, mesh)
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
        ! total surface stress (iceoce+atmoce) on elements 
        elnodes=elem2D_nodes(:,elem)
        
        !!PS stress_surf(1,elem)=sum(stress_iceoce_x(elnodes)*a_ice(elnodes) + &
        !!PS                         stress_atmoce_x(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
        !!PS stress_surf(2,elem)=sum(stress_iceoce_y(elnodes)*a_ice(elnodes) + &
        !!PS                         stress_atmoce_y(elnodes)*(1.0_WP-a_ice(elnodes)))/3.0_WP
        stress_surf(1,elem)=sum(stress_node_surf(1,elnodes))/3.0_WP
        stress_surf(2,elem)=sum(stress_node_surf(2,elnodes))/3.0_WP

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
subroutine oce_fluxes(ice, dynamics, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_DYN
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_CONFIG
    use o_ARRAYS
    use g_comm_auto
    use g_forcing_param, only: use_virt_salt, use_landice_water, use_age_tracer, use_age_mask, age_start_year !---fwf-code, age-code
    use g_forcing_arrays
    use g_support
    use cavity_interfaces
#if defined (__icepack)
    use icedrv_main,   only: icepack_to_fesom,    &
                            init_flux_atm_ocn
#endif
    use iceberg_params
    use iceberg_ocean_coupling
    use cavity_interfaces
    !---fwf-code
    use g_clock
    !---fwf-code-end

    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                    :: n, elem, elnodes(3),n1
    real(kind=WP)              :: rsss, net
    real(kind=WP), allocatable :: flux(:)
    !___________________________________________________________________________
    real(kind=WP), dimension(:,:), pointer :: temp, salt
    real(kind=WP), dimension(:)  , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)  , pointer :: a_ice_old
    real(kind=WP), dimension(:)  , pointer :: thdgr, thdgrsn
    real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
    real(kind=WP)                , pointer :: rhoice, rhosno, inv_rhowat

    !---wiso-code
    integer                    :: nt
    real(kind=WP), dimension(:,:), pointer :: wiso_oce1, wiso_oce2, wiso_oce3
    real(kind=WP), dimension(3) :: zfrac_freezing
    real(kind=WP), parameter   :: zwisomin = 1.e-6_WP
    real(kind=WP), allocatable :: snmelt(:), icemelt(:)
    real(kind=WP), allocatable :: wiso_prec_o16(:)
    real(kind=WP), allocatable :: wiso_rain(:,:), wiso_snow(:,:), wiso_melt(:,:)
    real(kind=WP), allocatable :: wiso_delta_rain(:,:),wiso_delta_snow(:,:),wiso_delta_ocean(:,:),wiso_delta_seaice(:,:)
    !---wiso-code-end

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

    !---wiso-code
    if (lwiso) then
      wiso_oce1 => tracers%data(index_wiso_tracers(1))%values(:,:)
      wiso_oce2 => tracers%data(index_wiso_tracers(2))%values(:,:)
      wiso_oce3 => tracers%data(index_wiso_tracers(3))%values(:,:)
      allocate(snmelt(myDim_nod2D+eDim_nod2D))
      allocate(icemelt(myDim_nod2D+eDim_nod2D))
      allocate(wiso_prec_o16(myDim_nod2D+eDim_nod2D))
      allocate(wiso_rain(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_snow(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_melt(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_delta_rain(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_delta_snow(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_delta_ocean(myDim_nod2D+eDim_nod2D,3))
      allocate(wiso_delta_seaice(myDim_nod2D+eDim_nod2D,3))
    end if
    !---wiso-code-end

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

    call icepack_to_fesom (nx_in         = (myDim_nod2D+eDim_nod2D), &
                           aice_out      = a_ice,                    &
                           vice_out      = m_ice,                    &
                           vsno_out      = m_snow,                   &
                           fhocn_tot_out = net_heat_flux,            &
                           fresh_tot_out = fresh_wa_flux,            &
                           fsalt_out     = real_salt_flux,           &
                           dhs_dt_out    = thdgrsn,                  &
                           dhi_dt_out    = thdgr,                    &
                           evap_ocn_out  = evaporation,              &
                           evap_out      = ice_sublimation           )

!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2d  
        ! Heat flux 
        heat_flux(n)       = - net_heat_flux(n)

        ! Freshwater flux (convert units from icepack to fesom)
        water_flux(n)      = - (fresh_wa_flux(n) * inv_rhowat) - runoff(n)

        ! Evaporation (convert units from icepack to fesom)
        evaporation(n)     = - evaporation(n) * (1.0_WP - a_ice(n)) * inv_rhowat

        ! Ice-Sublimation is added to to the freshwater in icepack --> see 
        ! icepack_therm_vertical.90 --> subroutine thermo_vertical(...): Line: 453
        ! freshn = freshn + evapn - (rhoi*dhi + rhos*dhs) / dt , evapn==sublimation
        ice_sublimation(n) = - ice_sublimation(n) * inv_rhowat
    end do
!$OMP END PARALLEL DO

    call init_flux_atm_ocn()

#else
!$OMP PARALLEL DO
    do n=1, myDim_nod2d+eDim_nod2d  
       heat_flux(n)   = -net_heat_flux(n)
       water_flux(n)  = -fresh_wa_flux(n)
    end do
!$OMP END PARALLEL DO
#endif

    if (use_icebergs) then
        call icb2fesom(mesh, partit, ice)
    end if

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
    ! ICEPACK: adds rain, snow and evap is based on the newly formed ice 
    !          concentration (a_ice). In our standard ice model rain, snow and evap is
    !          added based on the ice concentration of the previous time step (a_ice_old)  
    !          So for the proper balancing of snow the proper aice has to be choosen 
    !          -icepack_therm_itd.F90 --> subroutine icepack_step_therm2(...)
    !           fresh  = fresh + frain*aice
    !          -icedrv_step.F90: subroutine ocn_mixed_layer_icepack(...
    !           fresh_tot = fresh + (-evap_ocn + frain + fsnow)*(c1-aice)
    !          At the end all rain is added to the ocean, only snow needs to be
    !          scaled with (1-aice )
    !          -Ice-Sublimation is not added to evap in icepack, therefor we dont need
    !           to compensate for it the ice2atmos subplimation does not contribute 
    !           to the freshwater flux into the ocean
                   
!$OMP PARALLEL DO
    do n=1, myDim_nod2D+eDim_nod2D
        flux(n) = evaporation(n)                     &
                  -ice_sublimation(n)                 & ! the ice2atmos subplimation does not contribute to the freshwater flux into the ocean
                  +prec_rain(n)                       &                  
#if defined (__icepack)
                  +prec_snow(n)*(1.0_WP-a_ice(n))     &
#else
                  +prec_snow(n)*(1.0_WP-a_ice_old(n)) &
#endif

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

    !---wiso-code
    if (lwiso) then
      ! calculate snow melt (> 0.) and sea ice melt/growth (melt: > 0.; growth < 0.) as fraction of freshwater flux into the ocean
      snmelt = 0._WP
      icemelt = 0._WP
      where (abs(flux) > zwisomin) snmelt = -(thdgrsn*rhosno*inv_rhowat)/abs(flux)
      where (abs(flux) > zwisomin) icemelt = -(thdgr*rhoice*inv_rhowat)/abs(flux)
    end if
    !---wiso-code-end

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
    if (use_icebergs) then
        if (lbalance_fw .and. (.not. turn_off_fw)) then
            flux = flux + (ibfwb + ibfwe + ibfwl + ibfwbv) !* steps_per_ib_step
        end if
        
        call integrate_nod(ibfwb + ibfwe + ibfwl + ibfwbv, net, partit, mesh)
        if (mype==0) write(*,*) " * total iceberg fw flux: ", net
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
!---fwf-code-begin
    if(use_landice_water) then
!$OMP PARALLEL DO
      do n=1, myDim_nod2D+eDim_nod2D
         water_flux(n)=water_flux(n)-runoff_landice(n)*landice_season(month)
      end do
!$OMP END PARALLEL DO
    end if
    
    !___________________________________________________________________________
    if(lwiso .and. use_landice_water) then
!$OMP PARALLEL DO
      do n=1, myDim_nod2D+eDim_nod2D
         wiso_flux_oce(n,1)=wiso_flux_oce(n,1)+runoff_landice(n)*1000.0*wiso_smow(1)*(1-30.0/1000.0)*landice_season(month)
         wiso_flux_oce(n,2)=wiso_flux_oce(n,2)+runoff_landice(n)*1000.0*wiso_smow(2)*(1-240.0/1000.0)*landice_season(month)
         wiso_flux_oce(n,3)=wiso_flux_oce(n,3)+runoff_landice(n)*1000.0*landice_season(month)
      end do
!$OMP END PARALLEL DO
    end if
!---fwf-code-end

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

    !---wiso-code

    if (lwiso) then

    ! *** Important ***: The following wiso tracer order is assumed: nt=1: H218O, nt=2: HDO, nt=3: H216O


    ! (i) calculate isotope fluxes (received from coupled atmosphere model) into open water and onto sea ice

    ! atmospheric total H216O flux over open water and sea ice
    wiso_prec_o16 = (www3+iii3)*1000._WP
    ! integrate total H2O fluxes over all nodes for following flux corrections
    call integrate_nod(wiso_prec_o16, net, partit, mesh)

    nt=3
    wiso_rain(:,nt)=www3*1000._WP                              ! atmospheric H216O flux over open water
    wiso_flux_oce(:,nt)=wiso_rain(:,nt)-net/ocean_area         ! correction to enforce total H216O flux be zero

    do nt=1,2
      if (nt .EQ. 1) wiso_rain(:,nt)=www1*1000._WP/((20._WP/18._WP)*100._WP)           ! atmospheric H218O flux over open water
      if (nt .EQ. 2) wiso_rain(:,nt)=www2*1000._WP/((19._WP/18._WP)*2._WP*1000._WP)    ! atmospheric HDO flux over open water
      wiso_delta_rain(:,nt)=wiso_smow(nt)
      where (abs(wiso_rain(:,3)).gt.zwisomin) wiso_delta_rain(:,nt) = wiso_rain(:,nt)/wiso_rain(:,3)
      wiso_flux_oce(:,nt) = wiso_delta_rain(:,nt) * wiso_flux_oce(:,3) ! flux into ocean: assume same delta as atmospheric flux
    end do

    nt=3
    wiso_snow(:,nt)=iii3*1000._WP                              ! atmospheric H216O flux over sea ice
    wiso_flux_ice(:,nt)=wiso_snow(:,nt)-net/ocean_area         ! correction to enforce total H216O flux be zero
    where (a_ice(:).le.0.001_WP) wiso_flux_ice(:,nt)=0.0_WP    ! limit corrected H216O flux to sea ice areas, only

    do nt=1,2
      if (nt .EQ. 1) wiso_snow(:,nt)=iii1*1000._WP/((20._WP/18._WP)*100._WP)           ! atmospheric H218O flux over sea ice
      if (nt .EQ. 2) wiso_snow(:,nt)=iii2*1000._WP/((19._WP/18._WP)*2._WP*1000._WP)    ! atmospheric HDO flux over sea ice
      wiso_delta_snow(:,nt)=wiso_smow(nt)
      where (abs(wiso_snow(:,3)).gt.zwisomin) wiso_delta_snow(:,nt) = wiso_snow(:,nt)/wiso_snow(:,3)
      wiso_flux_ice(:,nt) = wiso_delta_snow(:,nt) * wiso_flux_ice(:,3)                 ! flux onto sea ice: assume same delta as atmospheric flux
    end do


    ! (ii) balance isotope fluxes for growing/melting of sea ice and melting of snow on sea ice

    ! set delta values of various H216O water masses to SMOW(O16) (=1.)
    nt=3
    wiso_delta_ocean(:,nt)=wiso_smow(nt)
    wiso_delta_seaice(:,nt)=wiso_smow(nt)
    wiso_delta_snow(:,nt)=wiso_smow(nt)

    ! calculate delta values of ocean surface water and sea ice for H218O and HDO
    nt=1
    do n=1, myDim_nod2D+eDim_nod2D
       ! calculate delta of open water (top ocean level)
       wiso_delta_ocean(n,nt)=wiso_smow(nt)
       if (wiso_oce3(1,n).gt.zwisomin) wiso_delta_ocean(n,nt) = wiso_oce1(1,n)/wiso_oce3(1,n)
       ! calculate delta of sea ice
       wiso_delta_seaice(n,nt)=wiso_smow(nt)
       if (tr_arr_ice(n,3).gt.zwisomin) wiso_delta_seaice(n,nt) = tr_arr_ice(n,nt)/tr_arr_ice(n,3)
    end do

    nt=2
    do n=1, myDim_nod2D+eDim_nod2D
       ! calculate delta of open water (top ocean level)
       wiso_delta_ocean(n,nt)=wiso_smow(nt)
       if (wiso_oce3(1,n).gt.zwisomin) wiso_delta_ocean(n,nt) = wiso_oce2(1,n)/wiso_oce3(1,n)
       ! calculate delta of sea ice
       wiso_delta_seaice(n,nt)=wiso_smow(nt)
       if (tr_arr_ice(n,3).gt.zwisomin) wiso_delta_seaice(n,nt) = tr_arr_ice(n,nt)/tr_arr_ice(n,3)
    end do

    ! for melting of snow on seaice (snmelt > 0.): assume no fractionation during melting process
    wiso_melt(:,:) = 0.0_WP

    nt=3
    wiso_melt(:,nt) = snmelt(:) * abs(wiso_flux_oce(:,nt))                      ! H216O melt amount = snow melt fraction x H216O freshwater flux over ocean
    wiso_melt(:,nt) = max(min(wiso_melt(:,nt),wiso_flux_ice(:,nt)),0._WP)       ! limit snow melt amount to range (0...wiso_flux_ice)
    where (a_ice(:).le.0.001_WP) wiso_melt(:,nt)=0.0_WP                         ! limit isotope changes by snow melt to sea ice areas, only

    do nt=1,2
       ! H218O and HDO meltwater has the same isotope ratio as snow on sea ice; no fractionation during melting
       wiso_melt(:,nt) = wiso_delta_snow(:,nt) * wiso_melt(:,3)
    end do

    wiso_flux_oce(:,:)= wiso_flux_oce(:,:) + wiso_melt(:,:)
    wiso_flux_ice(:,:)= wiso_flux_ice(:,:) - wiso_melt(:,:)

    ! for melting of seaice (icemelt > 0.): assume no fractionation during melting process
    ! for growing of seaice (icemelt < 0.): assume fractionation during freezing process
    ! (use equilibrium fractionation factors by Lehmann & Siegenthaler, JofGlaciology, 1991)
    zfrac_freezing = (/1.00291_WP, 1.0212_WP, 1.0_WP/)
    wiso_melt(:,:) = 0.0_WP

    nt=3
    wiso_melt(:,nt) = icemelt(:) * abs(wiso_flux_oce(:,nt))                     ! H216O melt (or growth) amount = sea ice melt fraction x H216O freshwater flux over ocean
    where (a_ice(:).le.0.001_WP) wiso_melt(:,nt)=0.0_WP                         ! limit isotope changes by melting/growing of sea ice to sea ice areas, only

    do nt=1,2
       ! H218O and HDO meltwater has the same isotope ratio as sea ice; no fractionation during melting
       where (wiso_melt(:,3) > 0.0_WP) wiso_melt(:,nt) = wiso_delta_seaice(:,nt) * wiso_melt(:,3)
       ! newly formed H218O and HDO sea ice has isotope ratio of ocean water; fractionation during growing considered
       where (wiso_melt(:,3) < 0.0_WP) wiso_melt(:,nt) = wiso_delta_ocean(:,nt)  * wiso_melt(:,3) * zfrac_freezing(nt)
    end do

    wiso_flux_oce(:,:)= wiso_flux_oce(:,:) + wiso_melt(:,:)
    wiso_flux_ice(:,:)= wiso_flux_ice(:,:) - wiso_melt(:,:)

    ! here: update sea ice isotope tracer concentration, only
    ! sea ice isotope tracer concentration are limited to sea ice areas, only
    ! (ocean water isotope tracers are updated in routine *oce_ale_tracer*)
    do n=1, myDim_nod2D+eDim_nod2D
       tr_arr_ice(n,:) = tr_arr_ice(n,:) + dt*wiso_flux_ice(n,:)
    end do

    do n=1, myDim_nod2D+eDim_nod2D
       if (tr_arr_ice(n,3) .gt. 1500._WP) then                           ! check if H216O tracer concentration reaches (arbitrary) limit of 1500.
          tr_arr_ice(n,1) = 1500._WP*tr_arr_ice(n,1)/tr_arr_ice(n,3)     ! reduce H2O18 based on the original ratio between H2O18 and H2O16 (i.e. the delta values are not changed)
          tr_arr_ice(n,2) = 1500._WP*tr_arr_ice(n,2)/tr_arr_ice(n,3)     ! reduce HDO16 based on the original ratio between HDO16 and H2O16 (i.e. the delta values are not changed)
          tr_arr_ice(n,3) = 1500._WP                                     ! reduce H216O to (arbitrary) upper limit
       endif
    end do

    do n=1, myDim_nod2D+eDim_nod2D
       if (tr_arr_ice(n,3) .le. 1._WP) then              ! check if H216O tracer concentration becomes too small or even negative
          tr_arr_ice(n,1) = wiso_smow(1)                 ! set delta H2O18 to SMOW value
          tr_arr_ice(n,2) = wiso_smow(2)                 ! set delta HDO16 to SMOW value
          tr_arr_ice(n,3) = 1._WP                        ! set H216O to lower limit
       endif
    end do

    end if  ! lwiso end

    !---wiso-code-end

    !---age-code-begin
    if (use_age_tracer) then
       tracers%data(index_age_tracer)%values(:,:) = tracers%data(index_age_tracer)%values(:,:) + dt/(86400.0*(365+fleapyear))

       if (use_age_mask) then
          tracers%data(index_age_tracer)%values(1,:) = tracers%data(index_age_tracer)%values(1,:) * (1-age_tracer_loc_index(:))
       else
          tracers%data(index_age_tracer)%values(1,:) = 0.0
       end if

       where (tracers%data(index_age_tracer)%values(:,:) .gt. yearnew-age_start_year+1)
          tracers%data(index_age_tracer)%values(:,:) = yearnew-age_start_year+1
       end where
    end if
    !---age-code-end

    !___________________________________________________________________________
    if (use_sw_pene) call cal_shortwave_rad(ice, partit, mesh)
    
    !___________________________________________________________________________
    deallocate(flux)
    
end subroutine oce_fluxes
!
!
!_______________________________________________________________________________
