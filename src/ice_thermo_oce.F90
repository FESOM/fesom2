module ice_thermodynamics_interfaces
    interface
        subroutine thermodynamics(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine thermodynamics
        
        subroutine cut_off(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine cut_off
    end interface  
end module ice_thermodynamics_interfaces

module ice_therm_interface
    interface
        subroutine therm_ice(ithermp, h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
        ug,ustar,T_oc,S_oc,H_ML,t,ice_dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
        rsf, dhgrowth, dhsngrowth, iflice, hflatow, hfsenow, hflwrdout, hfswrow, &
        hflwrow, hfradow, lid_clo,geolon, geolat, subli)
        USE MOD_ICE
        type(t_ice_thermo), intent(in), target :: ithermp
        real(kind=WP)   h, hsn, A, fsh, flo, Ta, qa, rain, snow, runo, rsss, &
                        ug, ustar, T_oc, S_oc, H_ML, t, ice_dt, ch, ce, ch_i, ce_i, evap_in, fw, ehf, &
                        dhgrowth, dhsngrowth, ahf, prec, subli, subli_i, rsf, &
                        rhow, show, rhice, shice, sh, thick, thact, lat, &
                        rh, rA, qhst, sn, hsntmp, o2ihf, evap, iflice, hflatow, &
                        hfsenow, hflwrdout, hfswrow, hflwrow, hfradow, lid_clo, geolon, geolat
        end subroutine therm_ice
    end interface
end module ice_therm_interface

module ice_budget_interfaces
    interface
        subroutine budget(ithermp, hice, hsn, t, ta, qa, fsh, flo, ug, S_oc, ch_i, ce_i, fh, subli)
        USE MOD_ICE
        type(t_ice_thermo), intent(in), target :: ithermp
        real(kind=WP)  hice, hsn, t, ta, qa, fsh, flo, ug, S_oc, ch_i, ce_i, fh, subli
        end subroutine budget
        
        subroutine obudget(ithermp, qa, fsh, flo, t, ug, ta, ch, ce, geolon, & 
                           geolat, fh, evap, hflatow, hfsenow, hflwrdout, hfswrow, &
                           hflwrow, hfradow) 
        USE MOD_ICE
        type(t_ice_thermo), intent(in), target :: ithermp
        real(kind=WP)   qa, t, ta, fsh, flo, ug, ch, ce, geolon, geolat, fh, evap, &
                        hfsenow, hflatow, hflwrdout, hfswrow, hflwrow, hfradow
        end subroutine obudget
        
        subroutine flooding(ithermp, h, hsn)
        USE MOD_ICE
        type(t_ice_thermo), intent(in), target :: ithermp
        real(kind=WP)   h, hsn
        end subroutine flooding
    end interface
end module ice_budget_interfaces
!
!
!_______________________________________________________________________________
subroutine cut_off(ice, partit, mesh)
    use o_param
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_ICE
    use g_config, only: use_cavity
    implicit none
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_ice),    intent(inout), target :: ice
    integer                               :: n
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
#if defined (__oifs) || defined (__ifsinterface)
    real(kind=WP), dimension(:), pointer  :: ice_temp
#endif /* (__oifs) */
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice    => ice%data(1)%values(:)
    m_ice    => ice%data(2)%values(:)
    m_snow   => ice%data(3)%values(:)
#if defined (__oifs) || defined (__ifsinterface)
    ice_temp => ice%data(4)%values(:)
#endif /* (__oifs) */

    !___________________________________________________________________________
    ! upper cutoff: a_ice
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n)
DO n=1, myDim_nod2D+eDim_nod2D
   if (a_ice(n) > 1.0_WP)   a_ice(n)=1.0_WP
    ! lower cutoff: a_ice
   if (a_ice(n) < .1e-8_WP) then
       a_ice(n)=0.0_WP
       m_ice(n)   =0.0_WP
       m_snow(n)  =0.0_WP
#if defined (__oifs) || defined (__ifsinterface)
        
        ice_temp(n)=273.15_WP
#endif /* (__oifs) */
   end if
    !___________________________________________________________________________
    ! lower cutoff: m_ice
   if (m_ice(n) < .1e-8_WP) then
        m_ice(n)=0.0_WP 
        m_snow(n)  =0.0_WP
        a_ice(n)   =0.0_WP
#if defined (__oifs) || defined (__ifsinterface)
        ice_temp(n)=273.15_WP
#endif /* (__oifs) */
   end if
     
    !___________________________________________________________________________
#if defined (__oifs) || defined (__ifsinterface)
    if (ice_temp(n) > 273.15_WP) ice_temp(n)=273.15_WP
#endif /* (__oifs) */

#if defined (__oifs) || defined (__ifsinterface)
    if (ice_temp(n) < 173.15_WP .and. a_ice(n) >= 0.1e-8_WP) ice_temp(n)=271.35_WP
#endif /* (__oifs) */
END DO
!$OMP END PARALLEL DO
end subroutine cut_off

#if !defined (__oasis) && !defined (__ifsinterface)
!_______________________________________________________________________________
! Sea-ice thermodynamics routines
!
! Coded by N. Yakovlev and S. Danilov.
! Adjusted for upgraded model physics (NCEP forcing data; parameterization
! of ice-ocean heat flux considering friction velocity, etc) 
! by Ralph Timmermann.
! Adjusted for general forcing data and NlFs option, cleaned up, bug fixing,
! by Qiang Wang, 13.01.2009
!_______________________________________________________________________________
subroutine thermodynamics(ice, partit, mesh)
  !
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodydnamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !------------------------------------------------------------------------
  
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_config
    use g_forcing_param
    use g_forcing_arrays
    use g_comm_auto
    use g_sbf, only: l_snow
    use ice_therm_interface
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    !_____________________________________________________________________________
    integer        :: i, j, elem
    real(kind=WP)  :: h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,rsf,evap_in
    real(kind=WP)  :: ug,ustar,T_oc,S_oc,h_ml,t,ch,ce,ch_i,ce_i,fw,ehf,evap
    real(kind=WP)  :: ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout, subli, hfswrow, hflwrow, hfradow
    real(kind=WP)  :: lid_clo, o2ihf
    real(kind=WP)  :: lat
    real(kind=WP)  :: geolon, geolat
    
    !_____________________________________________________________________________
    ! pointer on necessary derived types
    integer                      , pointer :: myDim_nod2D, eDim_nod2D
    integer      , dimension(:)  , pointer :: ulevels_nod2D
    real(kind=WP), dimension(:,:), pointer :: geo_coord_nod2D
    real(kind=WP), dimension(:)  , pointer :: u_ice, v_ice
    real(kind=WP), dimension(:)  , pointer :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:)  , pointer :: a_ice_old, m_ice_old, m_snow_old
    real(kind=WP), dimension(:)  , pointer :: thdgr, thdgrsn, thdgr_old, t_skin, ustar_aux
    real(kind=WP), dimension(:)  , pointer :: S_oc_array, T_oc_array, u_w, v_w
    real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
    real(kind=WP), external  :: TFrez  ! Sea water freeze temperature
    myDim_nod2d   => partit%myDim_nod2D
    eDim_nod2D    => partit%eDim_nod2D
    ulevels_nod2D  (1    :myDim_nod2D+eDim_nod2D) => mesh%ulevels_nod2D(:)
    geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D) => mesh%geo_coord_nod2D(:,:)    
    u_ice         => ice%uice(:)
    v_ice         => ice%vice(:)
    a_ice         => ice%data(1)%values(:)
    m_ice         => ice%data(2)%values(:)
    m_snow        => ice%data(3)%values(:)
    a_ice_old     => ice%data(1)%values_old(:)
    m_ice_old     => ice%data(2)%values_old(:)
    m_snow_old    => ice%data(3)%values_old(:)
    thdgr         => ice%thermo%thdgr
    thdgrsn       => ice%thermo%thdgrsn
    thdgr_old     => ice%thermo%thdgr_old
    t_skin        => ice%thermo%t_skin
    ustar_aux     => ice%thermo%ustar
    u_w           => ice%srfoce_u(:)
    v_w           => ice%srfoce_v(:)
    T_oc_array    => ice%srfoce_temp(:)
    S_oc_array    => ice%srfoce_salt(:)
    net_heat_flux => ice%flx_h(:)
    fresh_wa_flux => ice%flx_fw(:)
  
    !___________________________________________________________________________
    rsss=ref_sss
    
    ! u_ice and v_ice are at nodes
    ! u_w, v_w are at nodes (interpolated from elements)
    ! u_wind and v_wind are always at nodes
    !___________________________________________________________________________
    ! Friction velocity 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, elem, h, hsn, A, fsh, flo, Ta, qa, rain, snow, runo, rsss, rsf, evap_in, ug, ustar, T_oc, S_oc, &
!$OMP                                  h_ml, t, ch, ce, ch_i, ce_i, fw, ehf, evap, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout,    &
!$OMP                                  subli, lid_clo, lat, geolon, geolat, o2ihf)
!$OMP DO
    do i=1, myDim_nod2D
        ustar=0.0_WP
        if(ulevels_nod2d(i)>1) cycle 
        ustar=((u_ice(i)-u_w(i))**2 + (v_ice(i)-v_w(i))**2)
        ustar_aux(i)=sqrt(ustar*ice%cd_oce_ice)
    end do
!$OMP END DO
!$OMP MASTER
    call exchange_nod(ustar_aux, partit)
!$OMP END MASTER
!$OMP BARRIER

    !___________________________________________________________________________
!$OMP DO
    do i=1, myDim_nod2d+eDim_nod2D
        !_______________________________________________________________________
        ! if there is a cavity no sea ice thermodynamics is apllied
        if(ulevels_nod2d(i)>1) cycle 
        
        !_______________________________________________________________________
        ! prepare inputs for ice thermodynamics step
        h       = m_ice(i)
        hsn     = m_snow(i)
        A       = a_ice(i)
        fsh     = shortwave(i)
        flo     = longwave(i)
        Ta      = Tair(i)
        qa      = shum(i)

           
        if (.not. l_snow) then
            if (Ta>=0.0_WP) then
            rain=prec_rain(i)
            snow=0.0_WP
            else
            rain=0.0_WP
            snow=prec_rain(i)
            endif
            evap_in=evaporation(i) !evap_in: positive up
        else
            rain = prec_rain(i)
            snow = prec_snow(i)
            evap_in=0.0_WP
        end if
        runo    = runoff(i)
        ug      = sqrt(u_wind(i)**2+v_wind(i)**2)
        ustar   = ustar_aux(i)
        T_oc    = T_oc_array(i)      
        S_oc    = S_oc_array(i)
        if(ref_sss_local) rsss = S_oc
        t       = t_skin(i)   
        ch	    = Ch_atm_oce_arr(i)
        ce	    = Ce_atm_oce_arr(i)
        ch_i    = Ch_atm_ice
        ce_i    = Ce_atm_ice
        h_ml    = ice%thermo%h_ml               ! 10.0 or 30. used previously
        fw      = 0.0_WP
        ehf     = 0.0_WP
        geolon = geo_coord_nod2D(1, i)
        geolat = geo_coord_nod2D(2, i)
        
        if (geolat>0) then !TODO 2 separate pars for each hemisphere
            lid_clo=ice%thermo%h0
        else
            lid_clo=ice%thermo%h0_s
        endif
        

        !_______________________________________________________________________
        ! do ice thermodynamics
        call therm_ice(ice%thermo,h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
            ug,ustar,T_oc,S_oc,h_ml,t,ice%ice_dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
            rsf, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout, hfswrow, & 
            hflwrow, hfradow, lid_clo, geolon, geolat, subli)
        
        !_______________________________________________________________________
        ! write ice thermodyn. results into arrays
        ! backup of old values
        m_ice_old(i)      = m_ice(i) !PS
        m_snow_old(i)     = m_snow(i) !PS
        a_ice_old(i)      = a_ice(i) !PS
        thdgr_old(i)      = thdgr(i) !PS
        
        ! new values
        m_ice(i)          = h
        m_snow(i)         = hsn
        a_ice(i)          = A
        
        t_skin(i)         = t
        fresh_wa_flux(i)  = fw      !positive down
        net_heat_flux(i)  = ehf     !positive down
        evaporation(i)    = evap    !negative up
        ice_sublimation(i)= subli 
        
        thdgr(i)          = ithdgr
        thdgrsn(i)        = ithdgrsn
        flice(i)          = iflice
        
        
        ! Add minus sign to make from ... 
        !           ^(-)                      ^(+)
        !    |      |                  |      |
        ! ~~~|~~~~~~|~~~    to      ~~~|~~~~~~|~~~
        !    |      |                  |      | 
        !    v(+)                      v(-)   
        hf_Qlat(i)        = - hflatow  ! latent heat flux
        hf_Qsen(i)        = - hfsenow  ! sensible heat flux 
        hf_Qradtot(i)     = - hfradow  ! total radiation heat flux
        hf_Qswr(i)        = - hfswrow  ! shortwave radiation heat flux incoming
        hf_Qlwr(i)        = - hflwrow  ! longwave radiation heatflux incoming 
        hf_Qlwrout(i)     = - hflwrdout! longwave radiation heat flux outgoing
        ! --> the minus sign for net_heat_flux and fresh_wa_flux is added in 
        !     ice_oce_coupling.F90 in subroutine oce_fluxes(...=)
        
        ! real salt flux due to salinity that is contained in the sea ice 4-5 psu
        real_salt_flux(i) = rsf !PS
        
        ! if snow file is not given snow computed from prec_rain --> but prec_snow 
        ! array needs to be filled --> so that the freshwater balancing adds up
        if (.not. l_snow) then
            prec_rain(i)     = rain
            prec_snow(i)     = snow
        end if
    end do
!$OMP END DO
!$OMP END PARALLEL 
end subroutine thermodynamics
!
!
!_______________________________________________________________________________
subroutine therm_ice(ithermp, h, hsn, A, fsh, flo, Ta, qa, rain, snow, runo, rsss, &
                    ug, ustar, T_oc, S_oc, H_ML, t, ice_dt, ch, ce, ch_i, ce_i,    &
                    evap_in, fw, ehf, evap, rsf, dhgrowth, dhsngrowth, iflice,     &
                    hflatow, hfsenow, hflwrdout, hfswrow, hflwrow, hfradow, lid_clo, geolon, geolat, subli)
    ! Ice Thermodynamic growth model     
    !
    ! Input parameters:
    !------------------
    ! h - ice mass [m]
    ! hsn - snow mass [m]
    ! A - ice compactness
    ! fsh - shortwave radiation
    ! flo - longwave radiation
    ! Ta - air temperature
    ! qa - specific humidity
    ! rain - precipitation rain
    ! snow - precipitation snow
    ! runo - runoff
    ! ug - wind speed
    ! ustar - friction velocity
    ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer)
    ! H_ML - mixed layer depth - should be specified.
    ! t - temperature of snow/ice top surface
    ! ice_dt - time step [s]
    ! ch - transfer coefficient for sensible heat (for open ocean)
    ! ce - transfer coefficient for evaporation   (for open ocean)
    ! ch_i - transfer coefficient for sensible heat (for ice)
    ! ce_i - transfer coefficient for evaporation   (for ice)  
    ! lid_clo - lid closing parameter
    ! geolon - geographical longitude
    ! geolat - geographical latitude
    ! Output parameters:
    !-------------------
    ! h - ice mass
    ! hsn - snow mass
    ! A - ice compactness
    ! t - temperature of snow/ice top surface
    ! fw - freshwater flux due to ice melting [m water/ice_dt]
    ! ehf - net heat flux at the ocean surface [W/m2]        !RTnew
    ! subli - sublimatione over ice
    ! o2ihf - ocean to ice heat flux [W/m2] 

    USE MOD_ICE
    use g_forcing_param,  only: use_virt_salt  
    use o_param
    use ice_budget_interfaces
    implicit none
    type(t_ice_thermo), intent(in), target :: ithermp
    integer k
    real(kind=WP)  h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,evap_in
    real(kind=WP)  ug,ustar,T_oc,S_oc,H_ML,t,ice_dt,ch,ce,ch_i,ce_i,fw,ehf
    real(kind=WP)  dhgrowth,dhsngrowth,ahf,prec,subli,subli_i,rsf
    real(kind=WP)  rhow,show,rhice,shice,sh,snthick,thick,thact,lat
    real(kind=WP)  rh,rA,qhst,sn,hsntmp,o2ihf,evap
    real(kind=WP)  iflice, hflatow, hfsenow, hflwrdout, hfswrow, hflwrow, hfradow
    real(kind=WP), external  :: TFrez  ! Sea water freeze temperature.
    real(kind=WP)  lid_clo, geolon, geolat
    !___________________________________________________________________________
    logical      , pointer :: snowdist, new_iclasses
    integer      , pointer :: iclasses, open_water_albedo
    real(kind=WP), pointer :: hmin, Sice, Armin, cc, cl, con, consn, rhosno, rhoice, inv_rhowat, inv_rhosno, c_melt, h_cutoff
    real(kind=WP), pointer, dimension (:) :: hpdf
    snowdist          => ithermp%snowdist
    new_iclasses      => ithermp%new_iclasses
    iclasses          => ithermp%iclasses
    open_water_albedo => ithermp%open_water_albedo
    hmin       => ithermp%hmin
    Armin      => ithermp%Armin
    Sice       => ithermp%Sice
    cc         => ithermp%cc
    cl         => ithermp%cl
    con        => ithermp%con
    consn      => ithermp%consn
    iclasses   => ithermp%iclasses
    rhosno     => ithermp%rhosno
    rhoice     => ithermp%rhoice
    inv_rhowat => ithermp%inv_rhowat
    inv_rhosno => ithermp%inv_rhosno
    c_melt            => ithermp%c_melt
    h_cutoff          => ithermp%h_cutoff
    hpdf              => ithermp%hpdf
    
    !___________________________________________________________________________
    ! Store ice thickness at start of growth routine
    dhgrowth=h  	  

    ! determine h(i,j)/a(i,j) = actual ice thickness.
    ! if snow layer is present, add hsn weighted with quotient
    ! of conductivities of ice and snow, according to 0-layer approach
    ! of Semtner (1976).   	    
    ! thickness at the ice covered part
    snthick=hsn*(con/consn)/max(A,Armin)  ! Effective snow thickness
    thick=h/max(A,Armin)            ! Effective ice thickness
    if (snowdist) thick=snthick+thick   ! Effective ice and ice thickness

    ! Growth rate for ice in open ocean
    rhow=0.0_WP
    evap=0.0_WP
    call obudget(ithermp, qa,fsh,flo,T_oc,ug,ta,ch,ce,geolon, geolat, rhow, evap, &
                 hflatow, hfsenow, hflwrdout, hfswrow, hflwrow, hfradow) 
    hflatow  = hflatow  *(1.0_WP-A)   ! latent heatflux 
    hfsenow  = hfsenow  *(1.0_WP-A)   ! sensible heatflux 
    hfradow  = hfradow  *(1.0_WP-A)   ! total radiation hfswrow+hflwrow+hflwrdout
    hfswrow  = hfswrow  *(1.0_WP-A)   ! incoming shortwave radiation
    hflwrow  = hflwrow  *(1.0_WP-A)   ! incoming longwave radiation
    hflwrdout= hflwrdout*(1.0_WP-A)   ! outgoing long wave radiation 
    
    ! add heat loss at open ocean due to melting snow fall
    !rhow=rhow+snow*1000.0/rhoice !qiang
    ! ice_dt and (1-A) will be multiplied afterwards

    ! growth rate of ice in ice covered part
    ! following Hibler 1984
    ! assuming ice thickness has an euqal, 7-level distribution from zero to two times h 
    rhice=0.0_WP                      
    subli=0.0_WP
    if (thick.gt.hmin) then
        do k=1,iclasses
            thact = real((2*k-1),WP)*thick/real(iclasses,WP) ! Thicknesses of actual ice class
            if(new_iclasses) thact=h_cutoff/2.*thact       ! h_cutoff is variable (originally hcutoff was 2*h => factor 2 singles out)
            if(.not. snowdist) thact=thact+snthick         ! if snowdist=.true. snow depth is the same on every ice class
            call budget(ithermp, thact, hsn,t,Ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i) 
            !Thick ice K-class growth rate
            if(new_iclasses) then
                rhice=rhice+shice*hpdf(k)
                subli=subli+subli_i*hpdf(k)
             else
                rhice=rhice+shice
                subli=subli+subli_i
             end if
        end do
        if(.not. new_iclasses) then
            rhice=rhice/real(iclasses,WP)      	! Add to average heat flux
            subli=subli/real(iclasses,WP)
        end if
    end if
    
    ! Convert growth rates [m ice/sec] into growth per time step DT.
    rhow=rhow*ice_dt
    rhice=rhice*ice_dt

    ! Multiply ice growth of open water and ice
    ! with the corresponding areal fractions of grid cell
    show =rhow*(1.0_WP-A)
    shice=rhice*A
    sh   =show+shice

    ! Store atmospheric heat flux, average over grid cell [W/m**2]
    ahf=-cl*sh/ice_dt   

    ! precipitation (into the ocean)
    prec=rain+runo+snow*(1.0_WP-A)  	        ! m water/s

    ! snow fall above ice
    hsn=hsn+snow*ice_dt*A*1000.0_WP*inv_rhosno	! Add snow fall to temporary snow thickness    !!!
    dhsngrowth=hsn   		        ! Store snow thickness after snow fall 

    evap=evap*(1.0_WP-A)    		! m water/s
    subli=subli*A

    ! If there is atmospheric melting, first melt any snow that is present.
    ! Atmospheric heat flux available for melting
    ! sh = MINUS atm. heat flux / specific latent heat of sea ice
    ! Note: (sh<0) for melting, (sh>0) for freezing
    hsntmp= -min(sh,0.0_WP)*rhoice*inv_rhosno

    ! hsntmp is the decrease in snow thickness due to atmospheric melting
    ! [m/DT]. Do not melt more snow than available
    hsntmp=min(hsntmp,hsn)
    hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt

    ! Negative atmospheric heat flux left after melting of snow
    ! Note: (sh<0) and (hsntmp>0) for melting conditions
    ! hsntmp=0 for non-snow-melting conditions
    rh=sh+hsntmp*rhosno/rhoice
    h=max(h,0.0_WP)

    ! Compute heat content qhst of mixed layer - sea ice system
    !
    ! Total heat content is the sum of
    !	h	ice thickness after calculation of dynamic effects
    !	178418rh	change in ice thickness due to atmospheric forcing
    ! and heat available in mixed layer, with
    !	T_oc	temperature of ocean surface layer
    !	Tfrez	freezing point of sea water
    !	H_ML	thickness of uppermost layer
    !
    !RT:
    ! There are three possibilities to do this.
    ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
    !     Any heat available is then instantaneously used to melt ice.
    !     (so-called ice-bath approach)
    !     This is what used to be used in the Lemke sea ice-mixed layer model.
    ! rh=rh-(T_oc-TFrez(S_oc))*H_ML*cc/cl
    ! qhst=h+rh 
    !
    ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
    !     as a function of temperature difference. For a first step 
    !     we can assume a constant exchange coefficient gamma_t:
    ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
    !        +(T_oc-Tfrez(S_oc))*H_ML/ice_dt*cc*(1.0-A) ! [W/m2]
    ! rh=rh-o2ihf*ice_dt/cl
    ! qhst=h+rh		                      	! [m]
    !
    ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
    !     as a function of temperature difference and the
    !     friction velocity:
    o2ihf= (T_oc-TFrez(S_oc))*0.006_WP*ustar*cc*A  &
        +(T_oc-Tfrez(S_oc))*H_ML/ice_dt*cc*(1.0_WP-A)  	! [W/m2]
    rh=rh-o2ihf*ice_dt/cl
    qhst=h+rh		              		! [m]

    ! Melt snow if there is any ML heat content left (qhst<0).
    ! This may be the case if advection moves ice (with snow) to regions
    ! with a warm mixed layer.
    sn=hsn+min(qhst,0.0_WP)*rhoice*inv_rhosno

    ! New temporary snow thickness must not be negative:
    sn=max(sn,0.0_WP)

    ! Update snow and ice depth
    hsn=sn
    h=max(qhst,0.0_WP)
    if (h.lt.1E-6_WP) h=0._WP        ! Avoid very small ice thicknesses

    ! heat and fresh water fluxes
    dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
    dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting

    ! (without snow fall). This is a negative value (MINUS snow melt)

    dhgrowth=dhgrowth/ice_dt       ! Conversion: 'per time step' -> 'per second'
    dhsngrowth=dhsngrowth/ice_dt   ! Conversion: 'per time step' -> 'per second'
    ! (radiation+turbulent) + freezing(-melting) sea-ice&snow 

    ehf = ahf + cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth)
        
    ! (prec+runoff)+evap - freezing(+melting) ice&snow
    if (.not. use_virt_salt) then
        fw= prec+evap - dhgrowth*rhoice*inv_rhowat - dhsngrowth*rhosno*inv_rhowat
        rsf= -dhgrowth*rhoice*inv_rhowat*Sice
    else
        fw= prec+evap - dhgrowth*rhoice*inv_rhowat*(rsss-Sice)/rsss - dhsngrowth*rhosno*inv_rhowat 
        rsf = 0.0_WP
    end if
    
    ! Changes in compactnesses (equation 16 of Hibler 1979)
    rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
    rA= rhow - o2ihf*ice_dt/cl !Qiang: it was -(T_oc-TFrez(S_oc))*H_ML*cc/cl, changed in June 2010
    !rA= rhow - (T_oc-TFrez(S_oc))*H_ML*cc/cl*(1.0-A)
    A=A + c_melt*min(rh,0.0_WP)*A/max(h,hmin) + max(rA,0.0_WP)*(1._WP-A)/lid_clo  !/h0   
    !meaning:           melting                         freezing
    
    A=min(A,h*1.e6_WP)     ! A -> 0 for h -> 0
    A=min(max(A,0.0_WP),1._WP) ! A >= 0, A <= 1

    ! Flooding (snow to ice conversion)
    iflice=h
    call flooding(ithermp, h, hsn)     
    iflice=(h-iflice)/ice_dt
    
    ! to maintain salt conservation for the current model version
    !(a way to avoid producing net salt from snow-type-ice) 
    if (.not. use_virt_salt) then
        rsf=rsf-iflice*rhoice*inv_rhowat*Sice
    else
        fw=fw+iflice*rhoice*inv_rhowat*Sice/rsss
    end if
    
    evap=evap+subli
    
end subroutine therm_ice
!
!
!_______________________________________________________________________________
subroutine budget (ithermp, hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh,subli)
    ! Thick ice growth rate [m ice/sec]
    !
    ! INPUT:
    ! hice - actual ice thickness [m]
    ! hsn - snow thickness, used for albedo parameterization [m]
    ! t - temperature of snow/ice surface [C]
    ! ta - air temperature [C]
    ! qa - specific humidity [Kg/Kg]
    ! fsh - shortwave radiation [W/m**2]
    ! flo - longwave radiation  [W/m**2]
    ! ug - wind speed [m/sec]
    ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
    ! ch_i - transfer coefficient for sensible heat (for ice)
    ! ce_i - transfer coefficient for evaporation   (for ice) 
    !
    ! OUTPUT: fh - growth rate
    !
    ! qiang: The formular for saturated humidity was modified according to Large/Yeager2004
    ! to allow direct comparison with the CORE results (Griffies et al. 2009). The new
    ! formular does not require sea level pressure.
    ! A similar change was also made for the obudget routine.
    ! It was found through experiments that the results are quite similar to that from the
    ! original code, and the simulated ice volume is only slightly larger after modification. 
    use MOD_ICE
    use o_param, only: WP
    implicit none
    type(t_ice_thermo), intent(in), target :: ithermp
    integer iter, imax      ! Number of iterations
    real(kind=WP)  hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh
    real(kind=WP)  hfsen,hfrad,hflat,hftot,subli         
    real(kind=WP)  alb             ! Albedo of sea ice
    real(kind=WP)  q1, q2	  ! coefficients for saturated specific humidity
    real(kind=WP)  A1,A2,A3,B,C, d1, d2, d3   
    real(kind=WP), external :: TFrez
    !___________________________________________________________________________
    real(kind=WP), pointer :: boltzmann, emiss_ice, tmelt, cl, clhi, con, cpair, &
                              inv_rhowat, inv_rhoair, rhoair, albim, albi, albsn, albsnm
    boltzmann  => ithermp%boltzmann
    emiss_ice  => ithermp%emiss_ice
    tmelt      => ithermp%tmelt
    cl         => ithermp%cl
    clhi       => ithermp%clhi
    con        => ithermp%con
    cpair      => ithermp%cpair
    inv_rhowat => ithermp%inv_rhowat
    inv_rhoair => ithermp%inv_rhoair
    rhoair     => ithermp%rhoair
    albim      => ithermp%albim
    albi       => ithermp%albi
    albsn      => ithermp%albsn
    albsnm     => ithermp%albsnm
    
    !___________________________________________________________________________
    q1   = 11637800.0_WP
    q2   = -5897.8_WP
    imax = 5
    
    !___________________________________________________________________________
    ! set albedo
    ! ice and snow, freezing and melting conditions are distinguished.
    if (t<0.0_WP) then ! --> freezing condition    
        if (hsn.gt.0.0_WP) then ! --> snow cover present  
            alb=albsn         
        else                    ! --> no snow cover       
            alb=albi      
        endif
    else               ! --> melting condition     
        if (hsn.gt.0.0_WP) then ! --> snow cover present  
            alb=albsnm	    
        else                    ! --> no snow cover       
            alb=albim
        endif
    endif

    !___________________________________________________________________________
    d1=rhoair*cpair*Ch_i
    d2=rhoair*Ce_i
    d3=d2*clhi

    ! total incoming atmospheric heat flux
    A1=(1.0_WP-alb)*fsh + flo + d1*ug*ta + d3*ug*qa   ! in LY2004 emiss is multiplied wiht flo
    ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER

    do iter=1,imax
        B=q1*inv_rhoair*exp(q2/(t+tmelt))       ! (saturated) specific humidity over ice
        A2=-d1*ug*t-d3*ug*B &
            -emiss_ice*boltzmann*((t+tmelt)**4) ! sensible and latent heat,and outward radiation
        A3=-d3*ug*B*q2/((t+tmelt)**2)           ! gradient coefficient for the latent heat part
        C=con/hice                              ! gradient coefficient for downward heat conductivity
        A3=A3+C+d1*ug &                         ! gradient coefficient for sensible heat and radiation 
            +4.0_WP*emiss_ice*boltzmann*((t+tmelt)**3)    
        C=C*(TFrez(S_oc)-t)                     ! downward conductivity term
        
        t=t+(A1+A2+C)/A3                        ! NEW ICE TEMPERATURE AS THE SUM OF ALL COMPONENTS
    end do
    t=min(0.0_WP,t)
    
    !___________________________________________________________________________
    ! heat fluxes [W/m**2]:
    hfrad= (1.0_WP-alb)*fsh &               ! absorbed short wave radiation
        +flo &                              ! long wave radiation coming in  ! in LY2004 emiss is multiplied
        -emiss_ice*boltzmann*((t+tmelt)**4) ! long wave radiation going out

    hfsen=d1*ug*(ta-t)                    ! sensible heat 
    subli=d2*ug*(qa-B)                    ! sublimation
    hflat=clhi*subli                      ! latent heat

    hftot=hfrad+hfsen+hflat               ! total heat

    fh= -hftot/cl                         ! growth rate [m ice/sec]
                                          ! +: ML gains energy, ice melts
                                          ! -: ML loses energy, ice grows
    subli=subli*inv_rhowat                ! negative upward

    return
end subroutine budget
!
!
!_______________________________________________________________________________
subroutine obudget (ithermp, qa,fsh,flo,t,ug,ta,ch,ce,geolon, geolat, fh, evap, & 
                    hflatow, hfsenow, hflwrdout, hfswrow, hflwrow, hfradow)  
    ! Ice growth rate for open ocean [m ice/sec]
    !
    ! INPUT:
    ! t - temperature of open water [C]
    ! fsh - shortwave radiation
    ! flo - longwave radiation
    ! ta - air temperature [C]
    ! qa  - specific humidity             
    ! ug - wind speed [m/sec]
    ! ch - transfer coefficient for sensible heat
    ! ce - transfer coefficient for evaporation
    ! geolon - geographical longitude
    ! geolat - geographical latitude
    !
    ! OUTPUT: fh - growth rate
    !         evap - evaporation
    use MOD_ICE  
    use MOD_MESH
    use o_param, only: WP
    use g_clock
    implicit none
    type(t_ice_thermo), intent(in), target :: ithermp
    real(kind=WP) qa,t,ta,fsh,flo,ug,ch,ce,fh,evap
    real(kind=WP) hfsenow, hfswrow, hflwrow, hfradow, hflatow, hftotow, hflwrdout,b
    real(kind=WP) q1, q2 		! coefficients for saturated specific humidity
    real(kind=WP) c1, c4, c5, coszen, geolon, geolat
    real(kind=WP), external  :: compute_solar_zenith_angle, albw_taylor, albw_briegleb
    logical :: standard_saturation_shum_formula = .true.
    integer :: ii
    !___________________________________________________________________________
    integer, pointer :: open_water_albedo
    real(kind=WP), pointer :: boltzmann, emiss_wat, inv_rhowat, inv_rhoair, rhoair, &
                              tmelt, cl, clhw, cpair, albw
    boltzmann  => ithermp%boltzmann
    emiss_wat  => ithermp%emiss_wat
    inv_rhowat => ithermp%inv_rhowat
    inv_rhoair => ithermp%inv_rhoair
    rhoair     => ithermp%rhoair
    tmelt      => ithermp%tmelt
    cl         => ithermp%cl
    clhw       => ithermp%clhw
    cpair      => ithermp%cpair
    albw       => ithermp%albw
    open_water_albedo => ithermp%open_water_albedo
    
    !___________________________________________________________________________
    c1 = 3.8e-3_WP    
    c4 = 17.27_WP
    c5 = 237.3_WP
    q1 = 640380._WP
    q2 = -5107.4_WP
    if(open_water_albedo > 0)then
        coszen=compute_solar_zenith_angle(daynew, timenew/3600., geolon, geolat)
        if(open_water_albedo > 1)then
           albw=albw_briegleb(coszen)
        else
           albw=albw_taylor(coszen)
        endif
     endif

    ! (saturated) surface specific humidity
    if(standard_saturation_shum_formula) then
        b=c1*exp(c4*t/(t+c5))                      ! a standard one
    else
        b=0.98_WP*q1*inv_rhoair*exp(q2/(t+tmelt)) 	! LY2004 NCAR version 
    end if
    
    ! radiation heat fluxe [W/m**2]:
    hfswrow  = (1.0_WP-albw)*fsh
    hflwrow  = flo             	                        ! long wave radiation coming in !put emiss/check
    hflwrdout= -emiss_wat*boltzmann*((t+tmelt)**4) ! long wave radiation going out !in LY2004 emiss=1
    hfradow  = hfswrow + hflwrow + hflwrdout

    ! sensible heat fluxe [W/m**2]:
    hfsenow  = rhoair*cpair*ch*ug*(ta-t)             ! sensible heat 
    
    ! latent heat fluxe [W/m**2]:
    evap     = rhoair*ce*ug*(qa-b)  ! evaporation kg/m2/s
    hflatow  = clhw*evap                             ! latent heat W/m2

    ! total heat fluxe [W/m**2]:
    hftotow  = hfradow+hfsenow+hflatow               ! total heat W/m2
    
    fh= -hftotow/cl                             	! growth rate [m ice/sec]
    !                                           	+: ML gains energy, ice melts
    !                                           	-: ML loses energy, ice grows
    evap=evap*inv_rhowat 	 			! evaporation rate [m water/s],negative up !!!

    return
end subroutine obudget
!
!
!_______________________________________________________________________________
subroutine flooding (ithermp, h, hsn)
    use MOD_ICE
    type(t_ice_thermo), intent(in), target :: ithermp
    real(kind=WP) h,hsn,hdraft,hflood
    !___________________________________________________________________________
    real(kind=WP), pointer :: inv_rhowat, inv_rhosno, rhoice, rhosno
    inv_rhowat => ithermp%inv_rhowat
    inv_rhosno => ithermp%inv_rhosno
    rhoice     => ithermp%rhoice
    rhosno     => ithermp%rhosno

    !___________________________________________________________________________
    hdraft=(rhosno*hsn+h*rhoice)*inv_rhowat ! Archimedes: displaced water
    hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding
    h=h+hflood                          ! Add converted snow to ice volume
    hsn=hsn-hflood*rhoice*inv_rhosno        ! Subtract snow from snow layer

    !RT   This is what all AWI sea ice models do, but
    !RT   I wonder whether it really is correct for the heat budget.
    !RT   I suggest we initially keep it to allow for a comparison with BRIOS results
    !RT   and rethink it at a later stage.

    return
end subroutine flooding
!
!
!_______________________________________________________________________________
function TFrez(S)
    ! Nonlinear correlation for the water freezing temperature.
    ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
    use o_param, only: WP
    implicit none
    real(kind=WP) :: S, TFrez

    TFrez= -0.0575_WP*S+1.7105e-3_WP *sqrt(S**3)-2.155e-4_WP *S*S

end function TFrez


function compute_solar_zenith_angle(day_of_year, hour_utc, longitude, latitude) result(cos_zenith)
    !-----------------------------------------------------------------------
    ! Purpose: Compute the cosine of the solar zenith angle given the day of the year,
    !          time in UTC, and geographic coordinates.
    !
    ! Inputs:
    !   - day_of_year: Integer, day of the year (1 to 365 or 366)
    !   - hour_utc: Real, time in UTC (decimal hours)
    !   - longitude: Real, longitude in radians
    !   - latitude: Real, latitude in radians
    !
    ! Output:
    !   - cos_zenith: Real, cosine of the solar zenith angle
    !
    ! Notes:
    !   - Based on standard solar position calculations
    !   - Uses an approximation for solar declination
    !
    ! Written for Apache 2.0 licensed projects.
    !-----------------------------------------------------------------------

    use o_param, only: WP  ! Ensure precision consistency

    implicit none

    ! Input variables
    integer, intent(in)        :: day_of_year
    real(kind=WP), intent(in)  :: hour_utc
    real(kind=WP), intent(in)  :: longitude, latitude

    ! Output variable
    real(kind=WP)              :: cos_zenith

    ! Constants
    real(kind=WP), parameter   :: PI = 3.141592653589793_WP
    real(kind=WP), parameter   :: DEG_TO_RAD = PI / 180.0_WP
    real(kind=WP), parameter   :: RAD_TO_DEG = 180.0_WP / PI
    real(kind=WP), parameter   :: DAYS_PER_YEAR = 365.25_WP

    ! Computed values
    real(kind=WP)              :: solar_declination, hour_angle
    real(kind=WP)              :: solar_fraction

    ! Compute solar position parameters
    solar_fraction = 2.0_WP * PI * day_of_year / DAYS_PER_YEAR

    ! Approximate solar declination using Fourier terms
    solar_declination = 0.006918_WP - 0.399912_WP * cos(solar_fraction) &
                      + 0.070257_WP * sin(solar_fraction) &
                      - 0.006758_WP * cos(2.0_WP * solar_fraction) &
                      + 0.000907_WP * sin(2.0_WP * solar_fraction) &
                      - 0.002697_WP * cos(3.0_WP * solar_fraction) &
                      + 0.001480_WP * sin(3.0_WP * solar_fraction)

    ! Compute hour angle (relative to noon, converted to radians)
    hour_angle = (hour_utc - 12.0_WP) * 15.0_WP * DEG_TO_RAD + longitude

    ! Compute cosine of the solar zenith angle
    cos_zenith = sin(latitude) * sin(solar_declination) &
               + cos(latitude) * cos(solar_declination) * cos(hour_angle)

    ! Ensure cosine is non-negative (no below-horizon values)
    if (cos_zenith < 0.0_WP) cos_zenith = 0.0_WP

    return
end function compute_solar_zenith_angle

  
  function albw_taylor(coszen)
    !     purpose:  zenith depending open water albedo
    !     method:   taylor et al. (1996)
    !     author:   frank kauker
    !     date:     26. april 2017
    use o_param, only: WP
    implicit none
  
    real(kind=WP), intent(in)  :: coszen
    real(kind=WP) :: albw_taylor
  
    albw_taylor = 0.037_WP/(1.1_WP * coszen**1.4_WP + 0.15_WP)
  
    return
  end function albw_taylor
  
  function albw_briegleb(coszen)
    !     purpose:  zenith depending open water albedo
    !     method:   briegleb et al. (1986)
    !     author:   frank kauker
    !     date:     26. april 2017
    
    use o_param, only: WP
    implicit none
  
    real(kind=WP), intent(in)  :: coszen
    real(kind=WP)              :: albw_briegleb
  
    albw_briegleb = 0.026_WP/(1.1_WP * coszen**1.7_WP + 0.065_WP) + 0.15_WP * (coszen-1._WP)**2 * (coszen-0.5_WP)
  
    return
  end function albw_briegleb
  !
!
!
!_______________________________________________________________________________
#endif /* #if !defined (__oasis) && !defined (__ifsinterface) */
