
!
!_______________________________________________________________________________
subroutine oce_fluxes_mom(mesh)
    ! transmits the relevant fields from the ice to the ocean model
    !
    use o_PARAM
    use o_ARRAYS
    use MOD_MESH
    use i_ARRAYS
    use g_PARSUP
    use i_PARAM
    USE g_CONFIG
    use g_comm_auto

#if defined (__icepack)
    use icedrv_main,   only: icepack_to_fesom
#endif

    implicit none
    
    integer                  :: n, elem, elnodes(3),n1
    real(kind=WP)            :: aux, aux1
    type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"

    ! ==================
    ! momentum flux:
    ! ==================
    !___________________________________________________________________________

#if defined (__icepack)
     call icepack_to_fesom(nx_in=(myDim_nod2D+eDim_nod2D), &
                           aice_out=a_ice)
#endif

    do n=1,myDim_nod2D+eDim_nod2D   
        !_______________________________________________________________________
        ! if cavity node skip it 
        if (ulevels_nod2d(n)>1) cycle
        
        !_______________________________________________________________________
        if(a_ice(n)>0.001_WP) then
            aux=sqrt((u_ice(n)-u_w(n))**2+(v_ice(n)-v_w(n))**2)*density_0*Cd_oce_ice
            stress_iceoce_x(n) = aux * (u_ice(n)-u_w(n))
            stress_iceoce_y(n) = aux * (v_ice(n)-v_w(n))
        else
            stress_iceoce_x(n)=0.0_WP
            stress_iceoce_y(n)=0.0_WP
        end if
        
        ! total surface stress (iceoce+atmoce) on nodes 
        stress_node_surf(1,n) = stress_iceoce_x(n)*a_ice(n) + stress_atmoce_x(n)*(1.0_WP-a_ice(n))
        stress_node_surf(2,n) = stress_iceoce_y(n)*a_ice(n) + stress_atmoce_y(n)*(1.0_WP-a_ice(n))
    end do
    
    !___________________________________________________________________________
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
        !!PS stress_surf(1,elem)=sum(stress_node_surf(1,elnodes))/3.0_WP
        !!PS stress_surf(2,elem)=sum(stress_node_surf(2,elnodes))/3.0_WP
    END DO
    
    !___________________________________________________________________________
    if (use_cavity) call cavity_momentum_fluxes(mesh)
  
end subroutine oce_fluxes_mom
!
!
!_______________________________________________________________________________
subroutine ocean2ice(mesh)
  
    ! transmits the relevant fields from the ocean to the ice model

    use o_PARAM
    use o_ARRAYS
    use i_ARRAYS
    use MOD_MESH
    use g_PARSUP
    USE g_CONFIG
    use g_comm_auto
    implicit none

    type(t_mesh), intent(in) , target :: mesh
    integer :: n, elem, k
    real(kind=WP) :: uw, vw, vol

#include  "associate_mesh.h"

    ! the arrays in the ice model are renamed
        
    if (ice_update) then
        do n=1, myDim_nod2d+eDim_nod2d  
            if (ulevels_nod2D(n)>1) cycle 
            T_oc_array(n) = tr_arr(1,n,1)
            S_oc_array(n) = tr_arr(1,n,2)
            elevation(n)  = hbar(n)
        end do
    else
        do n=1, myDim_nod2d+eDim_nod2d
            if (ulevels_nod2D(n)>1) cycle 
            T_oc_array(n) = (T_oc_array(n)*real(ice_steps_since_upd,WP)+tr_arr(1,n,1))/real(ice_steps_since_upd+1,WP)
            S_oc_array(n) = (S_oc_array(n)*real(ice_steps_since_upd,WP)+tr_arr(1,n,2))/real(ice_steps_since_upd+1,WP)
            elevation(n)  = (elevation(n) *real(ice_steps_since_upd,WP)+      hbar(n))/real(ice_steps_since_upd+1,WP)
        !NR !PS      elevation(n)=(elevation(n)*real(ice_steps_since_upd)+eta_n(n))/real(ice_steps_since_upd+1,WP)
        !NR     elevation(n)=(elevation(n)*real(ice_steps_since_upd)+hbar(n))/real(ice_steps_since_upd+1,WP) !PS
        end do
!!PS         elevation(:)= (elevation(:)*real(ice_steps_since_upd)+hbar(:))/real(ice_steps_since_upd+1,WP)
    end if
    
    u_w = 0.0_WP
    v_w = 0.0_WP
    do n=1, myDim_nod2d  
        if (ulevels_nod2D(n)>1) cycle 
        uw  = 0.0_WP
        vw  = 0.0_WP
        vol = 0.0_WP
        do k=1, nod_in_elem2D_num(n)
            elem=nod_in_elem2D(k,n)
            if (ulevels(elem)>1) cycle
            !uw = uw+ UV(1,1,elem)*elem_area(elem)
            !vw = vw+ UV(2,1,elem)*elem_area(elem)
            vol = vol + elem_area(elem)
            uw  = uw+ UV(1,1,elem)*elem_area(elem)
            vw  = vw+ UV(2,1,elem)*elem_area(elem)
        end do
        !!PS uw = uw/area(1,n)/3.0_WP	  
        !!PS vw = vw/area(1,n)/3.0_WP
        uw = uw/vol
        vw = vw/vol
        
        if (ice_update) then
            u_w(n)=uw
            v_w(n)=vw
        else
            u_w(n)=(u_w(n)*real(ice_steps_since_upd,WP)+uw)/real(ice_steps_since_upd+1,WP)
            v_w(n)=(v_w(n)*real(ice_steps_since_upd,WP)+vw)/real(ice_steps_since_upd+1,WP)
        endif
    end do
    call exchange_nod(u_w, v_w)
end subroutine ocean2ice
!
!
!_______________________________________________________________________________
subroutine oce_fluxes(mesh)

  use MOD_MESH
  USE g_CONFIG
  use o_ARRAYS
  use i_ARRAYS
  use g_comm_auto
  use g_forcing_param, only: use_virt_salt,use_landice_water !---fwf-code, add use_landice_water
  use g_forcing_arrays
  use g_PARSUP
  use g_support
  use i_therm_param
  !---wiso-code
  use o_PARAM, only: num_wiso_tracers
  !---wiso-code-end
  !---fwf-code
  use g_clock 
  !---fwf-code-end
  
#if defined (__icepack)
  use icedrv_main,   only: icepack_to_fesom,    &
                           init_flux_atm_ocn
#endif

  implicit none
  type(t_mesh), intent(in)   , target :: mesh
  integer                    :: n, elem, elnodes(3),n1
  real(kind=WP)              :: rsss, net
  real(kind=WP), allocatable :: flux(:)

  !---wiso-code
  integer                    :: nt
  real(kind=WP), dimension(3) :: zfrac_freezing
  real(kind=WP)              :: zwisomin
  real(kind=WP), allocatable :: snmelt(:), icemelt(:)
  real(kind=WP), allocatable :: wiso_prec_o16(:)
  real(kind=WP), allocatable :: wiso_rain(:,:), wiso_snow(:,:), wiso_melt(:,:)
  real(kind=WP), allocatable :: wiso_delta_rain(:,:),wiso_delta_snow(:,:),wiso_delta_ocean(:,:),wiso_delta_seaice(:,:)
  !---wiso-code-end

#include  "associate_mesh.h"
    
    allocate(flux(myDim_nod2D+eDim_nod2D))
    flux = 0.0_WP

  !---wiso-code
  if (lwiso) then
    allocate(snmelt(myDim_nod2D+eDim_nod2D))
    allocate(icemelt(myDim_nod2D+eDim_nod2D))
    allocate(wiso_prec_o16(myDim_nod2D+eDim_nod2D))
    allocate(wiso_rain(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_snow(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_melt(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_delta_rain(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_delta_snow(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_delta_ocean(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    allocate(wiso_delta_seaice(myDim_nod2D+eDim_nod2D,num_wiso_tracers))
    zwisomin= 1.e-6_WP
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
    heat_flux   = -net_heat_flux 
    water_flux  = -fresh_wa_flux
#endif 
   
    if (use_cavity) call cavity_heat_water_fluxes_3eq(mesh)
    !!PS if (use_cavity) call cavity_heat_water_fluxes_2eq(mesh)
    
!!PS     where(ulevels_nod2D>1) heat_flux=0.0_WP
!!PS     where(ulevels_nod2D>1) water_flux=0.0_WP

    heat_flux_in=heat_flux ! sw_pene will change the heat_flux
    
    !___________________________________________________________________________
    call exchange_nod(heat_flux, water_flux) 

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
        do n=1, myDim_nod2D+eDim_nod2D
            !!PS if (ref_sss_local) rsss = tr_arr(1,n,2)
            if (ref_sss_local) rsss = tr_arr(ulevels_nod2d(n),n,2)
            virtual_salt(n)=rsss*water_flux(n) 
        end do
        
        if (use_cavity) then
            flux = virtual_salt
            where (ulevels_nod2d > 1) flux = 0.0_WP
            call integrate_nod(flux, net, mesh)
        else   
            call integrate_nod(virtual_salt, net, mesh)
        end if    
        virtual_salt=virtual_salt-net/ocean_area
    end if

    where (ulevels_nod2d == 1)
          dens_flux=sw_alpha(1,:) * heat_flux_in / vcpw + sw_beta(1, :) * (relax_salt + water_flux * tr_arr(1,:,2))
    elsewhere
          dens_flux=0.0_WP
    end where
    !___________________________________________________________________________
    ! balance SSS restoring to climatology
    if (use_cavity) then 
        do n=1, myDim_nod2D+eDim_nod2D
            relax_salt(n) = 0.0_WP
            if (ulevels_nod2d(n)>1) cycle
            !!PS relax_salt(n)=surf_relax_S*(Ssurf(n)-tr_arr(1,n,2))
            relax_salt(n)=surf_relax_S*(Ssurf(n)-tr_arr(ulevels_nod2d(n),n,2))
        end do
    else
        do n=1, myDim_nod2D+eDim_nod2D
            !!PS relax_salt(n)=surf_relax_S*(Ssurf(n)-tr_arr(1,n,2))
            relax_salt(n)=surf_relax_S*(Ssurf(n)-tr_arr(ulevels_nod2d(n),n,2))
        end do
    end if 
    
    ! --> if use_cavity=.true. relax_salt anyway zero where is cavity see above
    call integrate_nod(relax_salt, net, mesh)
    relax_salt=relax_salt-net/ocean_area
    
    !___________________________________________________________________________
    ! enforce the total freshwater/salt flux be zero
    ! 1. water flux ! if (.not. use_virt_salt) can be used!
    ! we conserve only the fluxes from the database plus evaporation.
    flux = evaporation-ice_sublimation     & ! the ice2atmos subplimation does not contribute to the freshwater flux into the ocean
            +prec_rain                       &
            +prec_snow*(1.0_WP-a_ice_old)    &
            +runoff    
    
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
        flux = flux-thdgr*rhoice*inv_rhowat-thdgrsn*rhosno*inv_rhowat
    end if     
    
    ! Also balance freshwater flux that come from ocean-cavity boundary
    if (use_cavity) then
        if (.not. use_virt_salt) then
            ! only for full-free surface approach otherwise total ocean volume will drift
            where (ulevels_nod2d > 1) flux = -water_flux
        else
            where (ulevels_nod2d > 1) flux =  0.0_WP
        end if 
    end if 
            
    call integrate_nod(water_flux, net, mesh)
    ! here the + sign must be used because we switched up the sign of the 
    ! water_flux with water_flux = -fresh_wa_flux, but evap, prec_... and runoff still
    ! have there original sign
    ! if use_cavity=.false. --> ocean_area == ocean_areawithcav
    water_flux=water_flux-net/ocean_areawithcav
    
!---wiso-code

    if (lwiso) then
    
    ! *** Important ***: The following wiso tracer order is assumed: nt=1: H218O, nt=2: HDO, nt=3: H216O
    

    ! (i) calculate isotope fluxes (received from coupled atmosphere model) into open water and onto sea ice

    ! atmospheric total H216O flux over open water and sea ice
    wiso_prec_o16 = (www3+iii3)*1000._WP
    ! integrate total H2O fluxes over all nodes for following flux corrections
    call integrate_nod(wiso_prec_o16, net, mesh)

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
    do n=1, myDim_nod2D+eDim_nod2D
      do nt=1,2
        ! calculate delta of open water (top ocean level)
        wiso_delta_ocean(n,nt)=wiso_smow(nt)
        if (tr_arr(1,n,num_tracers+3).gt.zwisomin) wiso_delta_ocean(n,nt) = tr_arr(1,n,num_tracers+nt)/tr_arr(1,n,num_tracers+3)
        !calculate delta of sea ice
        wiso_delta_seaice(n,nt)=wiso_smow(nt)
        if (tr_arr_ice(n,3).gt.zwisomin) wiso_delta_seaice(n,nt) = tr_arr_ice(n,nt)/tr_arr_ice(n,3)
      end do
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

!---fwf-code-begin
   if(use_landice_water) then
     do n=1, myDim_nod2D+eDim_nod2D
        water_flux(n)=water_flux(n)-runoff_landice(n)*landice_season(month)
     end do
   end if

   if(lwiso .and. use_landice_water) then
     do n=1, myDim_nod2D+eDim_nod2D
     wiso_flux_oce(n,1)=wiso_flux_oce(n,1)+runoff_landice(n)*1000.0*wiso_smow(1)*(1-30.0/1000.0)*landice_season(month)
     wiso_flux_oce(n,2)=wiso_flux_oce(n,2)+runoff_landice(n)*1000.0*wiso_smow(2)*(1-240.0/1000.0)*landice_season(month)
     wiso_flux_oce(n,3)=wiso_flux_oce(n,3)+runoff_landice(n)*1000.0*landice_season(month)
     end do
   end if
!---fwf-code-end

    !___________________________________________________________________________
    if (use_sw_pene) call cal_shortwave_rad(mesh)
    
    !___________________________________________________________________________
    deallocate(flux)
    
end subroutine oce_fluxes
!
!
!_______________________________________________________________________________
