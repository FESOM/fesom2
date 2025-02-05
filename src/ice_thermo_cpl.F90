#if defined (__coupled) || defined (__ifsinterface)
subroutine thermodynamics(ice, partit, mesh)

  !===================================================================
  !
  ! This subroutine computes thermodynamic changes of ice and snow
  ! for coupled simulations of FESOM2.
  ! It replaces the original FESOM scheme for uncoupled simulations.
  ! Note that atmospheric fluxes already need to be available.
  !
  ! Reference: Dorn et al. (2009), Ocean Modelling 29, 103-114.
  !
  ! Author: Wolfgang Dorn (AWI), Aug-2012
  !         Wolfgang Dorn (AWI), Oct-2012 (h0min adapted)
  !
  !===================================================================

  use o_param
  USE MOD_ICE
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_MESH
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_comm_auto
  use g_rotate_grid
  implicit none
  type(t_ice)   , intent(inout), target :: ice
  type(t_partit), intent(inout), target :: partit
  type(t_mesh)  , intent(in)   , target :: mesh
  !_____________________________________________________________________________
  integer :: inod
  !---- prognostic variables (updated in `ice_growth')
  real(kind=WP)  :: A, h, hsn, alb, t
  !---- atmospheric heat fluxes (provided by ECHAM)
  real(kind=WP)  :: a2ohf, a2ihf, qres, qcon
  !---- evaporation and sublimation (provided by ECHAM)
  real(kind=WP)  :: evap, subli
  !---- add residual freshwater flux over ice to freshwater (setted in ice_growth)
  real(kind=WP)  :: resid
  !---- precipitation and runoff (provided by ECHAM)
  real(kind=WP)  :: rain, snow, runo
  !---- ocean variables (provided by FESOM)
  real(kind=WP)  :: T_oc, S_oc, ustar
  !---- local variables (set in this subroutine)
  real(kind=WP)  :: rsss
  !---- output variables (computed in `ice_growth')
  real(kind=WP)  :: ehf, fw, rsf, dhgrowth, dhsngrowth, dhflice

  !---- geographical coordinates
  real(kind=WP)  :: geolon, geolat
  !---- minimum and maximum of the lead closing parameter
  real(kind=WP)  :: h0min = 0.5, h0max = 1.5

  real(kind=WP), parameter :: Aimin = 0.001, himin = 0.005

  !_____________________________________________________________________________
  ! pointer on necessary derived types
  integer      ,                 pointer :: myDim_nod2D, eDim_nod2D
  integer      , dimension(:)  , pointer :: ulevels_nod2D
  real(kind=WP), dimension(:,:), pointer :: geo_coord_nod2D
  real(kind=WP), dimension(:)  , pointer :: u_ice, v_ice
  real(kind=WP), dimension(:)  , pointer :: a_ice, m_ice, m_snow
  real(kind=WP), dimension(:)  , pointer :: thdgr, thdgrsn
  real(kind=WP), dimension(:)  , pointer :: a_ice_old, m_ice_old, m_snow_old, thdgr_old 
  real(kind=WP), dimension(:)  , pointer :: S_oc_array, T_oc_array, u_w, v_w
  real(kind=WP), dimension(:)  , pointer :: fresh_wa_flux, net_heat_flux
#if defined (__oifs) || defined (__ifsinterface)
  real(kind=WP), dimension(:) , pointer  :: ice_temp, ice_alb, enthalpyoffuse, ice_heat_qres, ice_heat_qcon
#endif
#if defined (__coupled) || defined (__ifsinterface)
  real(kind=WP), dimension(:)  , pointer ::  oce_heat_flux, ice_heat_flux 
#endif 
  real(kind=WP)                , pointer :: rhoice, rhosno, rhowat, rhofwt, Sice, cl, cc, cpice, consn, con 
  myDim_nod2d=>partit%myDim_nod2D
  eDim_nod2D =>partit%eDim_nod2D
  ulevels_nod2D  (1    :myDim_nod2D+eDim_nod2D) => mesh%ulevels_nod2D
  geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D) => mesh%geo_coord_nod2D(:,:)
  u_ice         => ice%uice(:)
  v_ice         => ice%vice(:)
  a_ice         => ice%data(1)%values(:)
  m_ice         => ice%data(2)%values(:)
  m_snow        => ice%data(3)%values(:)  
  thdgr         => ice%thermo%thdgr(:)
  thdgrsn       => ice%thermo%thdgrsn(:)
  a_ice_old     => ice%data(1)%values_old(:)
  m_ice_old     => ice%data(2)%values_old(:)
  m_snow_old    => ice%data(3)%values_old(:)
  thdgr_old     => ice%thermo%thdgr_old
  T_oc_array    => ice%srfoce_temp(:)
  S_oc_array    => ice%srfoce_salt(:)
  u_w           => ice%srfoce_u(:)
  v_w           => ice%srfoce_v(:)
  fresh_wa_flux => ice%flx_fw(:)
  net_heat_flux => ice%flx_h(:)
#if defined (__oifs) || defined (__ifsinterface)
  ice_temp      => ice%data(4)%values(:)
  ice_alb       => ice%atmcoupl%ice_alb(:)
  enthalpyoffuse=> ice%atmcoupl%enthalpyoffuse(:)
  ice_heat_qres => ice%atmcoupl%flx_qres(:)
  ice_heat_qcon => ice%atmcoupl%flx_qcon(:)
#endif 
#if defined (__coupled) || defined (__ifsinterface)
  oce_heat_flux => ice%atmcoupl%oce_flx_h(:)
  ice_heat_flux => ice%atmcoupl%ice_flx_h(:)
#endif
  rhoice        => ice%thermo%rhoice  
  rhosno        => ice%thermo%rhosno
  rhowat        => ice%thermo%rhowat
  rhofwt        => ice%thermo%rhofwt
  Sice          => ice%thermo%Sice
  cl            => ice%thermo%cl
  cc            => ice%thermo%cc
  cpice         => ice%thermo%cpice
  consn         => ice%thermo%consn
  con           => ice%thermo%con
  rhoice        => ice%thermo%rhoice
  !_____________________________________________________________________________  
  rsss = ref_sss

  !---- loop over all surface node
  do inod=1,myDim_nod2d+eDim_nod2D

     if (ulevels_nod2D(inod) > 1) cycle

     A       = a_ice(inod)
     h       = m_ice(inod)
     hsn     = m_snow(inod)

#if defined (__oifs) || defined (__ifsinterface)
     a2ohf   = oce_heat_flux(inod) + shortwave(inod) + enthalpyoffuse(inod)
#else
     a2ohf   = oce_heat_flux(inod) + shortwave(inod)
#endif
     a2ihf   = ice_heat_flux(inod)
     evap    = evap_no_ifrac(inod)
     subli   = sublimation(inod)
     rain    = prec_rain(inod)
     snow    = prec_snow(inod)
     runo    = runoff(inod)     

     ustar   = sqrt(ice%cd_oce_ice)*sqrt((u_ice(inod)-u_w(inod))**2+(v_ice(inod)-v_w(inod))**2)
     T_oc    = T_oc_array(inod)      
     S_oc    = S_oc_array(inod)
     if (ref_sss_local) rsss = S_oc

     ehf     = 0._WP
     fw      = 0._WP
     if (.not. use_virt_salt) then
        rsf     = 0._WP
     end if

#if defined (__oifs) || defined (__ifsinterface)

     !---- For AWI-CM3 we calculate ice surface temp and albedo in fesom,
     ! then send those to OpenIFS where they are used to calucate the 
     ! energy fluxes ---!
     t                   = ice_temp(inod)
     qres     = 0.0_WP
     qcon     = 0.0_WP
     if(A>Aimin) then
        call ice_surftemp(ice%thermo, max(h/(max(A,Aimin)),0.05), hsn/(max(A,Aimin)), a2ihf, t)
        ice_temp(inod)  = t
     else
        ! Freezing temp of saltwater in K
        ice_temp(inod) = -0.0575_WP*S_oc_array(inod) + 1.7105e-3_WP*sqrt(S_oc_array(inod)**3) -2.155e-4_WP*(S_oc_array(inod)**2)+273.15_WP        
     endif
     call ice_albedo(ice%thermo, h, hsn, t, alb)
     ice_alb(inod)       = alb
     ice_heat_qres(inod) = qres
     ice_heat_qcon(inod) = qcon
#endif
     call ice_growth

     !__________________________________________________________________________
     ! save old ice variables
     m_ice_old(inod)      = m_ice(inod) 
     m_snow_old(inod)     = m_snow(inod)
     a_ice_old(inod)      = a_ice(inod) 
     thdgr_old(inod)      = thdgr(inod) 
     
     !__________________________________________________________________________
     ! save new ice variables
     a_ice(inod)          = A
     m_ice(inod)          = h
     m_snow(inod)         = hsn
     net_heat_flux(inod)  = ehf
     fresh_wa_flux(inod)  = fw
     if (.not. use_virt_salt) then
        real_salt_flux(inod)= rsf
     end if
     thdgr(inod)          = dhgrowth
     thdgrsn(inod)        = dhsngrowth
     flice(inod)          = dhflice
     
     !---- total evaporation (needed in oce_salt_balance.F90) = evap+subli
     evaporation(inod)    = evap + subli
     ice_sublimation(inod)= subli
#if defined (__oasis) || defined (__ifsinterface)
     residualifwflx(inod) = resid
#endif     
  enddo
  return

contains

  !===================================================================
  ! Thermodynamic ice growth model     
  !===================================================================

  subroutine ice_growth
      
    implicit none
 
    !---- thermodynamic production rates (pos.: growth; neg.: melting)
    real(kind=WP)  :: dsnow, dslat, dhice, dhiow, dcice, dciow

    !---- heat fluxes (positive upward, negative downward)
    real(kind=WP)  :: Qatmice, Qatmocn, Qocnice, Qocnatm, Qicecon
    real(kind=WP)  :: ahf, ohf

    !---- atmospheric freshwater fluxes (precipitation minus evaporation)
    real(kind=WP)  :: PmEice, PmEocn

    !---- local variables and dummys
    real(kind=WP)  :: hold, hsnold, htmp, hsntmp, heff, h0cur, hdraft, hflood

    !---- cut-off ice thickness (hcutoff) used to avoid very small ice
    !---- thicknesses as well as division by zero. NOTE: the standard
    !---- cut-off ice thickness hmin=0.05 is set in `i_therm_parms'
    !---- and is questionable in terms of conservation of energy.
    real(kind=WP), parameter :: hcutoff = 1.e-6

    !---- minimum ice concentration (Aimin) and ice thickness (himin)
    real(kind=WP), parameter :: Aimin = 0.001, himin = 0.005

    !---- an arbitrary big value, but note that bigval*hcutoff should
    !---- be greater than one (= maximum ice concentration)
    real(kind=WP), parameter :: bigval = 1.e10

    !---- heat transfer rate (gamma_t = h_ml/tau0, where h_ml is the
    !---- mixed layer depth and tau0 is a damping time constant for a
    !---- delayed adaptation of the mixed layer temperature. We assume
    !---- this rate to be 10 meters per day. NOTE: tau0 should be
    !---- significantly greater than the time step dt
    real(kind=WP), parameter :: gamma_t = 10./86400.

    !---- freezing temperature of freshwater [deg C]
    real(kind=WP), parameter :: Tfrez0 = 0.

    !---- freezing temperature of sea-water [deg C]
    real(kind=WP)  :: Tfrezs

    !---- compute freezing temperature of sea-water from salinity
    TFrezs = -0.0575_WP*S_oc + 1.7105e-3_WP*sqrt(S_oc**3) - 2.155e-4_WP*(S_oc**2)

    !---- effective thermodynamic thickness of the snow/ice layer
    heff = h + hsn*con/consn
    heff = heff/max(A,Aimin)

    !---- conductive heat flux through the snow/ice layer for melting
    !---- conditions at the top of the layer (Tice = Tfrez0)
    Qicecon = (Tfrezs-Tfrez0)*con/max(heff,himin)

    !---- atmospheric heat fluxes (provided by the atmosphere model)

#if defined (__oifs) || defined (__ifsinterface)
    Qatmice = -qres-qcon
#else
    Qatmice = -a2ihf
#endif
    Qatmocn = -a2ohf

    !---- oceanic heat fluxes
    !---- NOTE: for freezing conditions: Qocnatm < Qatmocn (due to
    !---- latent heat release), otherwise the fluxes must be balanced:
    !---- Qocnatm = Qatmocn (no ice melt over open water)
    !!Qocnice = (T_oc-Tfrezs)*0.006*ustar*cc
    !!Qocnatm = (T_oc-Tfrezs)*h_ml/dt*cc
    Qocnice = (T_oc-Tfrezs)*gamma_t*cc
    Qocnatm = min(Qocnice,Qatmocn)

    !---- total atmospheric and oceanic heat fluxes
    !---- average over grid cell [W/m**2]
    ahf = A*Qatmice + (1._WP-A)*Qatmocn
    ohf = A*Qocnice + (1._WP-A)*Qocnatm



    !---- convert heat fluxes [W/m**2] into growth per time step dt [m]
    Qicecon = Qicecon*dt/cl
    Qatmice = Qatmice*dt/cl
    Qatmocn = Qatmocn*dt/cl
    Qocnice = Qocnice*dt/cl
    Qocnatm = Qocnatm*dt/cl

    !---- atmospheric freshwater fluxes (provided by the atmosphere model)
    !---- NOTE: evaporation and sublimation represent potential fluxes and
    !---- must be area-weighted (like the heat fluxes); in contrast,
    !---- precipitation (snow and rain) and runoff are effective fluxes
!already weighted in IFS coupling
#if !defined (__ifsinterface)
    subli  = A*subli
    evap   = (1._WP-A)*evap
#endif
    PmEice = A*snow + subli
    PmEocn = evap + rain + (1._WP-A)*snow + runo

    !---- convert freshwater fluxes [m/s] into growth per time step dt [m]
    PmEice = PmEice*dt
    PmEocn = PmEocn*dt

    !---- add snowfall minus sublimation to temporary snow thickness
    hsn = hsn + PmEice*rhofwt/rhosno

    !---- residual freshwater flux over ice
    if (hsn.lt.0) then
       PmEice = hsn*rhosno/rhofwt
       hsn = 0._WP
    else
       PmEice = 0._WP
    endif

    !---- subtract sublimation from ice thickness (PmEice <= 0)
    h = h + PmEice*rhofwt/rhoice

    !---- residual freshwater flux over ice
    if (h.lt.0) then
       PmEice = h*rhoice/rhofwt
       h = 0._WP
    else
       PmEice = 0._WP
    endif
    resid  = PmEice/dt
    
    !---- add residual freshwater flux over ice to freshwater flux over ocean
    PmEocn = PmEocn + PmEice
    PmEice = 0._WP

    !---- store snow and ice thickness after snowfall and sublimation
    hsnold = hsn
    hold = h

    !---- snow melt rate over sea ice (dsnow <= 0)
    !---- if there is atmospheric melting over sea ice, first melt any
    !---- snow that is present, but do not melt more snow than available
#if defined (__oifs) || defined (__ifsinterface)
    !---- new condition added - surface temperature must be
    !----                       larger than 273K to melt snow
    if (t.gt.273_WP) then
        dsnow = A*min(Qatmice-Qicecon,0._WP)
        dsnow = max(dsnow*rhoice/rhosno,-hsn)
    else
        dsnow = 0.0_WP
    endif
#else
    dsnow = A*min(Qatmice-Qicecon,0._WP)
    dsnow = max(dsnow*rhoice/rhosno,-hsn)
#endif 

    !---- update snow thickness after atmospheric snow melt
    hsn = hsn + dsnow

    !---- ice growth/melt rate over sea ice (dhice)
    dhice = A*(Qatmice-Qocnice)

    !---- subtract atmospheric heat already used for snow melt
    dhice = dhice - dsnow*rhosno/rhoice

    !---- ice growth rate over open water (dhiow >= 0)
    dhiow = (1._WP-A)*max(Qatmocn-Qocnatm,0._WP)

    !---- temporary new ice thickness [m]
    htmp = h + dhice + dhiow

    if (htmp.lt.0._WP) then
       !---- all ice melts; now try to melt snow if still present,
       !---- but do not melt more snow than available
       hsntmp = max(htmp*rhoice/rhosno,-hsn)

       !---- update snow thickness after snow melt
       hsn = hsn + hsntmp

       !---- new ice thickness
       h = 0._WP
    else
       h = htmp
    endif

    !---- avoid very small ice thicknesses
    if (h.lt.hcutoff) h = 0._WP

    !---- ice thickness before any thermodynamic change
    !---- (for h=0 use cut-off ice thickness to avoid division by zero)
    htmp = max(hold,hcutoff)

    !---- ice concentration change by melting of ice (dhice <= 0)
    dcice = 0.5_WP*A*min(dhice,0._WP)/htmp

    !---- lateral snow melt if lateral ice melt exceeds snow melt
    !---- due to atmospheric forcing ( dcice*hsn/A - dsnow < 0 )
    if (A.le.0._WP) then
       dslat = -hsn
    else
       hsntmp = hsnold/max(A,Aimin)
       dslat = min(dcice*hsntmp-dsnow,0._WP)
       dslat = max(dslat,-hsn)
    endif

    !---- update snow thickness after lateral snow melt
    hsn = hsn + dslat

    !---- subtract heat required to melt this additional amount of
    !---- snow from total oceanic heat flux (dslat <= 0)
    ohf = ohf + dslat*rhosno/rhoice*cl/dt

    !---- lead closing parameter
    h0cur = max(h0min,min(h0max,hold))

    !---- alternative lead closing parameter when h0max is negative.
    !---- h0min is then interpreted as 'Phi_F' according to the ice
    !---- model by Mellor and Kantha (1989)
    if (h0max.le.0._WP) then
       htmp = hold/max(A,Aimin)
       h0cur = max(htmp,himin)/h0min
    endif

    !---- ice concentration change by freezing in open leads (dhiow >= 0)
    !---- NOTE: dhiow already represents an areal fraction
    dciow = max(dhiow,0._WP)/h0cur

    !---- new ice concentration
    A = A + dcice + dciow

    !---- set a=0 for h=0
    A = min(A,h*bigval)

    !---- restrict ice concentration to values between zero and one
    A = max(0._WP,min(1._WP,A))

    !---- change in snow and ice thickness due to thermodynamic effects [m/s]
    dhsngrowth = (hsn-hsnold)/dt
    dhgrowth = (h-hold)/dt

    !---- convert growth per time step dt [m] into freshwater fluxes [m/s]
    PmEocn = PmEocn/dt

    !---- total freshwater mass flux into the ocean [kg/m**2/s]
    if (.not. use_virt_salt) then
       fw  = PmEocn*rhofwt - dhgrowth*rhoice - dhsngrowth*rhosno
       rsf = -dhgrowth*rhoice*Sice/rhowat
    else
       fw = PmEocn*rhofwt - dhgrowth*rhoice*(rsss-Sice)/rsss - dhsngrowth*rhosno 
    end if

    !---- total energie flux into the ocean [W/m**2] (positive downward)
    !---- NOTE: ehf = -ohf (in case of no cut-off)
    ehf = -ahf + cl*(dhgrowth + dhsngrowth*rhosno/rhoice) 

    !---- store ice thickness before flooding (snow to ice conversion)
    htmp = h

    !---- Archimedes: displaced water
    hdraft = (h*rhoice + hsn*rhosno)/rhowat

    !---- increase in mean ice thickness due to flooding
    hflood = hdraft - min(h,hdraft)

    !---- add converted snow to ice thickness
    h = h + hflood

    !---- subtract converted snow from snow thickness
    hsn = hsn - hflood*rhoice/rhosno

    !---- rate of flooding snow to ice
    dhflice = (h-htmp)/dt

    !---- to maintain salt conservation for the current model version
    !---- (a way to avoid producing net salt from snow-type-ice) 
    if (.not. use_virt_salt) then
       rsf = rsf - dhflice*rhoice*Sice/rhowat
    else
       fw = fw + dhflice*rhoice*Sice/rsss
    end if

    !---- convert freshwater mass flux [kg/m**2/s] into sea-water volume flux [m/s]
    fw   = fw/rhowat
    ! keep in mind that for computation of FW all imposed fluxes were accounted with the ratio rhofwt/rhowat:
    !evap = evap *rhofwt/rhowat
    !rain = rain *rhofwt/rhowat
    !snow = snow *rhofwt/rhowat
    !runo = runo *rhofwt/rhowat
    !subli= subli*rhofwt/rhowat
    !resid= resid*rhofwt/rhowat
    return
  end subroutine ice_growth

 subroutine ice_surftemp(ithermp, h,hsn,a2ihf,t)
  ! INPUT:
  ! a2ihf - Total atmo heat flux to ice
  ! A  - Ice fraction
  ! h  - Ice thickness
  ! hsn   - Snow thickness
  ! 
  ! INPUT/OUTPUT:
  ! t     - Ice surface temperature

  implicit none
  type(t_ice_thermo), intent(in), target :: ithermp
  !---- atmospheric heat net flux into to ice (provided by OpenIFS)
  real(kind=WP)  a2ihf
  !---- ocean variables (provided by FESOM)
  real(kind=WP)  h
  real(kind=WP)  hsn
  real(kind=WP)  t
  !---- local variables
  real(kind=WP)  snicecond
  real(kind=WP)  zsniced
  real(kind=WP)  zicefl
  real(kind=WP)  hcapice
  real(kind=WP)  zcpdt
  real(kind=WP)  zcpdte
  real(kind=WP)  zcprosn
  !---- local parameters
  real(kind=WP), parameter :: dice  = 0.10_WP                       ! Thickness for top ice "layer"
  !---- freezing temperature of sea-water [K]
  real(kind=WP)  :: TFrezs
  
  real(kind=WP), pointer :: con, consn, cpsno, rhoice, rhosno
  con    => ice%thermo%con
  consn  => ice%thermo%consn
  cpsno  => ice%thermo%cpsno
  rhoice => ice%thermo%rhoice
  rhosno => ice%thermo%rhosno

  !---- compute freezing temperature of sea-water from salinity
  TFrezs = -0.0575_WP*S_oc + 1.7105e-3_WP*sqrt(S_oc**3) - 2.155e-4_WP*(S_oc**2)+273.15

  snicecond = con/consn                 ! equivalence fraction thickness of ice/snow
  zsniced=h+snicecond*hsn               ! Ice + Snow-Ice-equivalent thickness [m]
  zicefl=con*TFrezs/zsniced             ! Conductive heat flux through sea ice [W/m²]
  hcapice=rhoice*cpice*dice             ! heat capacity of upper 0.05 cm sea ice layer [J/(m²K)]
  zcpdt=hcapice/dt                      ! Energy required to change temperature of top ice "layer" [J/(sm²K)]
  zcprosn=rhosno*cpsno/dt               ! Specific Energy required to change temperature of 1m snow on ice [J/(sm³K)]
  zcpdte=zcpdt !+zcprosn*hsn            ! Combined Energy required to change temperature of snow + 0.05m of upper ice
  t=(zcpdte*t+a2ihf+zicefl)/(zcpdte+con/zsniced) ! New sea ice surf temp [K]
  if (t>273.15_WP) then
     qres=(con/zsniced+zcpdte)*(t-273.15_WP)
     t=273.15_WP
  endif
  qcon=con*(t-TFrezs)/max(zsniced, himin)
! t=min(273.15_WP,t)
 end subroutine ice_surftemp

 subroutine ice_albedo(ithermp, h, hsn, t, alb)
  ! INPUT:
  ! h      - ice thickness [m]
  ! hsn    - snow thickness [m]
  ! t      - temperature of snow/ice surface [C]
  ! 
  ! OUTPUT:
  ! alb    - selected broadband albedo
  implicit none
  type(t_ice_thermo), intent(in), target :: ithermp
  real(kind=WP) :: h
  real(kind=WP) :: hsn    
  real(kind=WP) :: t    
  real(kind=WP) :: alb
  real(kind=WP) :: geolat
  real(kind=WP), pointer :: albsn, albi, albsnm, albim
  albsn  => ice%thermo%albsn
  albi   => ice%thermo%albi
  albsnm => ice%thermo%albsnm
  albim  => ice%thermo%albim
  
  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished
  if (h>0.0_WP) then
     if (t<273.15_WP) then         ! freezing condition    
        if (hsn.gt.0.001_WP) then !   snow cover present  
           alb=albsn       
        else                    !   no snow cover       
           alb=albi           
        endif
     else                               ! melting condition     
        if (hsn.gt.0.001_WP) then !   snow cover present  
           alb=albsnm          
        else                    !   no snow cover       
           alb=albim
        endif
     endif
   else
      alb=0.066_WP            !  ocean albedo
   endif
 end subroutine ice_albedo

end subroutine thermodynamics
#endif
