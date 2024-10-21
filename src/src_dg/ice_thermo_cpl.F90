#if defined (__oasis)
subroutine thermodynamics(mesh)

  !===================================================================
  !
  ! This subroutine computes thermodynamic changes of ice and snow
  ! for coupled simulations of FESOM and ECHAM.
  ! It replaces the original FESOM scheme for uncoupled simulations.
  ! Note that atmospheric fluxes already need to be available.
  !
  ! Reference: Dorn et al. (2009), Ocean Modelling 29, 103-114.
  !
  ! Author: Wolfgang Dorn (AWI), Aug-2012
  !         Wolfgang Dorn (AWI), Oct-2012 (h0min adapted)
  !
  !===================================================================
  !---- variables from oce_modules.F90
#if 0
  use o_param,          only: ref_sss, ref_sss_local
#ifdef use_fullfreesurf
  use o_array,          only: real_salt_flux
#endif
  use g_parsup,         only: myDim_nod2D, eDim_nod2D
#ifdef use_cavity
  use o_mesh,           only: coord_nod2D, cavity_flag_nod2d
#else
  use o_mesh,           only: coord_nod2D
#endif

  !---- variables from ice_modules.F90
  use i_dyn_parms,      only: Cd_oce_ice
  use i_therm_parms,    only: rhowat, rhoice, rhosno, cc, cl, con, consn, Sice
  use i_array,          only: a_ice, m_ice, m_snow, u_ice, v_ice, u_w, v_w  &
       , fresh_wa_flux, net_heat_flux, oce_heat_flux, ice_heat_flux, S_oc_array, T_oc_array

  !---- variables from gen_modules_config.F90
  use g_config,         only: dt

  !---- variables from gen_modules_forcing.F90
  use g_forcing_arrays, only: shortwave, evap_no_ifrac, sublimation  &
       , prec_rain, prec_snow, runoff, evaporation, thdgr, thdgrsn, flice

  !---- variables from gen_modules_rotate_grid.F90
  use g_rotate_grid,    only: r2g
#endif

  use o_param
  use mod_mesh
  use i_therm_param
  use i_param
  use i_arrays
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_parsup
  use g_comm_auto
  use g_rotate_grid
  implicit none

  integer :: inod
  !---- prognostic variables (updated in `ice_growth')
  real(kind=WP)  :: A, h, hsn, alb, t
  !---- atmospheric heat fluxes (provided by ECHAM)
  real(kind=WP)  :: a2ohf, a2ihf
  !---- evaporation and sublimation (provided by ECHAM)
  real(kind=WP)  :: evap, subli
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
  real(kind=WP)  :: h0min = 0.50, h0max = 1.5
  type(t_mesh), intent(in)   , target :: mesh  

#include  "associate_mesh.h"


  rsss = ref_sss

  !---- total evaporation (needed in oce_salt_balance.F90)
  evaporation = evap_no_ifrac*(1.-a_ice) + sublimation*a_ice

  !---- loop over all surface node
  do inod=1,myDim_nod2d+eDim_nod2d

#ifdef use_cavity
     if (cavity_flag_nod2d(inod).eq.1) cycle
#endif

     A       = a_ice(inod)
     h       = m_ice(inod)
     hsn     = m_snow(inod)

     a2ohf   = oce_heat_flux(inod) + shortwave(inod)
     a2ihf   = ice_heat_flux(inod)
     evap    = evap_no_ifrac(inod)
     subli   = sublimation(inod)
     rain    = prec_rain(inod)
     snow    = prec_snow(inod)
     runo    = runoff(inod)     

     ustar   = sqrt(Cd_oce_ice)*sqrt((u_ice(inod)-u_w(inod))**2+(v_ice(inod)-v_w(inod))**2)
     T_oc    = T_oc_array(inod)      
     S_oc    = S_oc_array(inod)
     if (ref_sss_local) rsss = S_oc

     ehf     = 0._WP
     fw      = 0._WP
#ifdef use_fullfreesurf
     rsf     = 0._WP
#endif

     !---- different lead closing parameter for NH and SH
     call r2g(geolon, geolat, coord_nod2d(1,inod), coord_nod2d(2,inod))
!     if (geolat.lt.0.) then
!        h0min = 1.0
!        h0max = 1.5
!     else
!        h0min = 0.75
!        h0max = 1.0
!     endif

     call ice_growth
#if defined (__oifs)
     call ice_albedo(hsn,t,alb)
     ice_alb(inod)	 = alb
     ice_temp(inod)      = t+tmelt
#endif


     a_ice(inod)         = A
     m_ice(inod)         = h
     m_snow(inod)        = hsn
     net_heat_flux(inod) = ehf
     fresh_wa_flux(inod) = fw
#ifdef use_fullfreesurf
     real_salt_flux(inod)= rsf
#endif
     thdgr(inod)         = dhgrowth
     thdgrsn(inod)       = dhsngrowth
     flice(inod)         = dhflice

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

    !---- density of freshwater [kg/m**3].
    real(kind=WP), parameter :: rhofwt = 1000.

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
    Qatmice = -a2ihf
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
    PmEice = A*snow + A*subli
    PmEocn = rain + runo + (1._WP-A)*snow + (1._WP-A)*evap

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

    !---- add residual freshwater flux over ice to freshwater flux over ocean
    PmEocn = PmEocn + PmEice
    PmEice = 0._WP

    !---- store snow and ice thickness after snowfall and sublimation
    hsnold = hsn
    hold = h

    !---- snow melt rate over sea ice (dsnow <= 0)
    !---- if there is atmospheric melting over sea ice, first melt any
    !---- snow that is present, but do not melt more snow than available
    dsnow = A*min(Qatmice-Qicecon,0._WP)
    dsnow = max(dsnow*rhoice/rhosno,-hsn)

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
#ifdef use_fullfreesurf
    fw = PmEocn*rhofwt - dhgrowth*rhoice - dhsngrowth*rhosno 
    rsf = -dhgrowth*rhoice*Sice
#else
    fw = PmEocn*rhofwt - dhgrowth*rhoice*(rsss-Sice)/rsss - dhsngrowth*rhosno 
#endif

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
#ifdef use_fullfreesurf
    rsf = rsf - dhflice*rhoice*Sice
#else
    fw = fw + dhflice*rhoice*Sice/rsss
#endif

    !---- convert freshwater mass flux [kg/m**2/s] into sea-water volume flux [m/s]
    fw = fw/rhowat
#ifdef use_fullfreesurf
    rsf = rsf/rhowat
#endif

    return
  end subroutine ice_growth


 subroutine ice_albedo (hsn,t,alb)
  ! INPUT:
  ! hsn - snow thickness, used for albedo parameterization [m]
  use i_therm_param
  implicit none

  real(kind=WP)  hsn    
  real(kind=WP)  t    
  real(kind=WP)  alb             ! Albedo of sea ice

  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  if (hsn.gt.0.0) then	!   snow cover present  
     alb=albsn         	
  else              		!   no snow cover       
     alb=albi       	
  endif
 end subroutine ice_albedo

end subroutine thermodynamics
#endif
