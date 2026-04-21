!===============================================================================
! Specific declarations related to carbon isotope simulations
!-------------------------------------------------------------------------------
module REcoM_ciso
  implicit none
  save


! Options for carbon isotope simulations (see namelist.recom)
  logical :: ciso_init             = .false.  ! Initial fractionation of bulk organic matter
  logical :: ciso_14               = .false.  ! Include radiocarbon (-> 31 or 38 tracers)
  logical :: ciso_organic_14       = .false.  ! Include organic radiocarbon (-> 38 tracers)
  real(kind=8) :: delta_co2_13        = -6.61
  real(kind=8) :: big_delta_co2_14(3) = (/0., 0., 0./) 
  real(kind=8) :: lambda_14 = 3.8561e-12      ! Decay constant of carbon-14
! for revised atbox 14CO2 implementation
  logical      :: atbox_spinup    = .true.
  real(kind=8) :: cosmic_14_init  = 2.0       ! Initial 14C production flux (atoms / s / cm**2)


  namelist / paciso / ciso_init, ciso_14, ciso_organic_14, &
                      lambda_14, delta_co2_13, big_delta_co2_14, &
                      atbox_spinup, cosmic_14_init

! Extensions of other modules or subroutines
! Module REcoM_constants: ciso tracer indices
  integer :: idic_13, iphyc_13, idetc_13, ihetc_13, idoc_13, idiac_13, iphycal_13, idetcal_13, &
             idic_14, iphyc_14, idetc_14, ihetc_14, idoc_14, idiac_14, iphycal_14, idetcal_14

! Module REcoM_declarations:
  real(kind=8)  :: quota_13, quota_14, quota_dia_13, quota_dia_14,                     &  ! quotas
                   recipQuota_13, recipQuota_14, recipQuota_dia_13, recipQuota_dia_14, &  ! reciprocal quotas
                   recipQZoo_13, recipQZoo_14
  real(kind=8)  :: HetRespFlux_13, HetRespFlux_14           ! zooplankton respiration fluxes
  real(kind=8)  :: calcification_13, calcification_14, &    ! calcification
                   calc_loss_agg_13, calc_loss_agg_14, &
                   calc_loss_gra_13, calc_loss_gra_14, &
                   calc_diss_13, calc_diss_14

! Module REcoM_GloVar:
  real(kind=8),dimension(12)              :: AtmCO2_13                       ! [uatm] Atmospheric 13CO2 partial pressure. One value for the whole planet for each month
  real(kind=8),dimension(3,12)            :: AtmCO2_14                       ! [uatm] Atmospheric 14CO2 partial pressure. Three latitude zones for each month
  real(kind=8),allocatable,dimension(:)   :: GloPCO2surf_13, GloPCO2surf_14  ! [uatm] Surface ocean 13|14CO2 partial pressure
  real(kind=8),allocatable,dimension(:)   :: GloCO2flux_13, GloCO2flux_14    ! [mmol/m2/day] Positive downwards
  real(kind=8),allocatable,dimension(:)   :: GloCO2flux_seaicemask_13, GloCO2flux_seaicemask_14
  real(kind=8),allocatable,dimension(:)   :: RiverineDOCOrig_13, RiverineDOCOrig_14, RiverDOC2D_13, RiverDOC2D_14 

! Module REcoM_LocVar:
  real(kind=8) :: pCO2surf_13(1), pCO2surf_14(1), & ! [uatm] Partial pressure of 13|14CO2 in surface layer at current 2D node
                  co2flux_13(1), co2flux_14(1), &      ! mocsy output: air-to-sea flux of 13|14CO2 [mol/(m^2 * s)]
                  co2flux_seaicemask_13(1), co2flux_seaicemask_14(1) ! air-to-sea flux of CO2 [mmol/m2/s]
  real(kind=8) :: LocAtmCO2_13(1), LocAtmCO2_14(1)  ! [uatm]
  real(kind=8) :: LocRiverDOC_13, LocRiverDOC_14    ! CHECK

! Subroutines REcoM_main & REcoM_extra:
  real(kind=8) :: lat_val                           ! nodal latitude (of atmospheric input)

! Subroutine REcoM_extra:
  real(kind=8) :: delta_co2_14                      ! atmospheric Delta14CO2

! Subroutine REcoM_forcing:
  real(kind=8) :: co2sat,                         & ! dissolved CO2 at saturation (CO2*air) [mol / m**3]
                  kwco2                             ! piston velocity of CO2

! Subroutine REcoM_sms:
  real(kind=8) :: DIC_13, DIC_14,                 & ! [mmol/m3] Dissolved Inorganic 13|14Carbon
                  PhyC_13, PhyC_14,               & ! [mmol/m3] Intracellular conc of 13|14Carbon in small phytoplankton
                  DetC_13, DetC_14,               & ! [mmol/m3] Conc of 13|14C in Detritus
                  HetC_13, HetC_14,               & ! [mmol/m3] Conc of 13|14C in heterotrophs
                  EOC_13, EOC_14,                 & ! [mmol/m3] Extracellular Organic 13|14C conc
                  DiaC_13, DiaC_14,               & ! [mmol/m3] Intracellular conc of 13|14Carbon in diatoms
                  PhyCalc_13, PhyCalc_14,         & ! [mmol/m3] Conc of 13|14C in calcite of phytoplankton
                  DetCalc_13, DetCalc_14            ! [mmol/m3] Conc of 13|14C in calcite of detritus
  real(kind=8),allocatable,dimension(:) :: Cphot_z, Cphot_dia_z ! Vertical profiles of photosynthesis rates, fesom1: 46 -> 47 in fesom2

! Subroutine REcoM_init:
  real(kind=8),allocatable,dimension(:,:) :: delta_dic_13_init, &     ! auxiliary initial
                                             delta_dic_14_init, &     ! d|Delta13|14C
                                             big_delta_dic_14_init    ! fields

! Atmospheric box model (global variables):
  real(kind=8),allocatable,dimension(:) ::   x_co2atm_13, x_co2atm_14, & ! atmospheric CO2 mixing ratio (mole fraction)
                                             cosmic_14                   ! cosmogenic 14 production (mol / s)
  real(kind=8) :: production_rate_to_flux_14, &                          ! conversion factor
                  r_atm_spinup_13, r_atm_spinup_14                       ! 13|14CO2 / 12CO2 spinup ratios

! Specific factors related the carbon-isotopic composition
! Isotopic ratios
  real(kind=8) :: r_atm_13, r_atm_14,             & ! atmospheric CO2
                  r_co2s_13, r_co2s_14,           & ! dissolved CO2
                  r_dic_13, r_dic_14,             & ! DIC in seawater
                  r_phyc_13, r_phyc_14,           & ! nanophytoplankton
                  r_diac_13, r_diac_14,           & ! diatoms
                  r_iorg_13 = 0.975,              & ! initial ratios of organic matter
                  r_iorg_14 = 0.950
                  
! Fractionation factors
  real(kind=8) :: alpha_k_13 = 0.99912,           & ! gas transfer (kinetic fractionation,
                  alpha_k_14 = 0.99824,           & ! mean values for 5-21C by Zhang et al., 1995)
                  alpha_aq_13, alpha_aq_14,       & ! dissolution of CO2 in sewater (equilibrium fractionation)
                  alpha_dic_13, alpha_dic_14,     & ! hydrolysis / dissociation of CO2 <-> DIC (equilibrium fract.)
                  alpha_p_13, alpha_p_14,         & ! photosynthesis of nanophytoplankton
                  alpha_p_dia_13, alpha_p_dia_14, & ! photosynthesis of diatoms
                  alpha_calc_13 = 1.000,          & ! calcification (Romanek et al., 1992: 1.001, 1.002)
                  alpha_calc_14 = 1.000,          &
                  alpha_dcal_13 = 1.000,          & ! dissolution of calcite (Romanek et al., 1992: 0.999, 0.998)
                  alpha_dcal_14 = 1.000
! Radioactive decay constant of carbon-14
! t1/2 = 5700 years (Bé et al., 2013; recommended by Orr et al., 2017, for OMIP-BGC)
! if 1 year := 365.25 days:  lambda_14 = 3.8534e-12 / second
! if 1 year := 365.00 days:  lambda_14 = 3.8561e-12 / second
! if 1 year := 360    days:  lambda_14 = 3.9096e-12 / second

! Tracer IDs to be considered in decay calculations (oce_ale_tracer.F90)
  integer, dimension(8) :: c14_tracer_id = (/1402, 1405, 1408, 1410, 1412, 1414, 1420, 1421/)

  contains


    subroutine recom_ciso_airsea(tempc, co3, dic)
!   ----------------------------------------------------------------------------------
!     Subroutine to calculate carbon-isotopic fractionation during air-sea exchange
!   ----------------------------------------------------------------------------------
!
!     Input variables:
!     tempc              lokal temperature in C
!     co3                carbonate ion concentration
!     dic                total carbon concentration
!
!     Output variables, defined in module REcoM_ciso:
!     alpha_k_13,14      kinetic fract. factors for gas transfer
!     alpha_aq_13,14     equilib. fract. factors for dissolution
!     alpha_dic_13,14    equilib. fract. factors for DIC <-> CO2
!
!     Internal variables:
!     epsilon_aq_13,14   equilib. fractionation for dissolution
!     epsilon_dic_13,14  equilib. fractionation for DIC <-> CO2
!     fco3               total carbon fraction
!
!     mbutzin, 2016 - 2019.


!     Declarations
      implicit none

      real(kind=8), intent(in) ::  tempc, co3, dic
      real(kind=8) ::  epsilon_aq_13, epsilon_dic_13, fco3

!     Calculation of carbon-isotopic fractionation factors, where
!
!     alpha_xy   = Rx / Ry               = fractionation factor
!     epsilon_xy = (alpha_xy - 1) * 1000 = fractionation (in per mill)
!     epsilon_14 = 2 * epsilon_13 => alpha_14 = 2 * alpha_13 - 1.

!     We use parametrisations and numerical values determined for carbon-13
!     by Zhang et al. (1995).

!     Kinetic fractionation during gas transfer, mean values between 5 and 21C
!     (values are defined in module REcoM_ciso)
!     epsilon_k_13 = -0.86 => alpha_k_13 =  0.99914, alpha_k_14 =  0.99828

!     Equilibrium fractionation during gas dissolution
!
      epsilon_aq_13 = 0.0049 * tempc - 1.31
      alpha_aq_13   = 1. + 0.001 * epsilon_aq_13

!     Equilibrium fractionation between DIC and CO2
!
!     The equilibrium fractionation between DIC and CO2 cannot be simply
!     calculated from the fractionation factors for HCO3, CO3 and CO2star.
!     Here, we employ an empirical function involving fCO3 = [CO3] / DIC
!     assuming that fCO3 is the same for all carbon isotopes
      fco3 = co3 / dic
      epsilon_dic_13 = (0.014 * fco3 - 0.107 ) * tempc + 10.53
      alpha_dic_13 = 1. + 0.001 * epsilon_dic_13

!     Fractionation of radiocarbon
      if (ciso_organic_14) then
        alpha_aq_14  = 2. * alpha_aq_13  - 1.
        alpha_dic_14 = 2. * alpha_dic_13 - 1.
      else
!       no fractionation in the inorganic approximation
        alpha_aq_14  = 1.
        alpha_dic_14 = 1.
      end if

    return
    end subroutine recom_ciso_airsea
!   ----------------------------------------------------------------------------------

    subroutine recom_ciso_photo (co2st)
!   ----------------------------------------------------------------------------------
!        Subroutine calculating carbon-isotopic fractionation during photosynthesis
!   ----------------------------------------------------------------------------------
!     Input:
!     dissolved CO2 (co2st) in mol / m**3
!
!     Output:
!     isotopic fractionation factors for phytoplankton and diatoms due to
!     photosynthesis (alpha_p_13|14, declared at the head of the module)
!
!     Note that we are interested in effective values (implictly including the 
!     fractionation of dissolved CO2) which are actually derived in field studies
!     or lab experiments. Young et al. 2013, eq. (5) with values from paragraph [35]
!
!     Here, we follow Young et al. 2013, eq. (5) with values from paragraph [35]
!     eps_p = eps_pm * (1. - rho / co2aq) = 17.6 * (1 - 2.02 / co2aq)
!     where co2aq is in umol / L
!
!     mbutzin, 2017 - 2021.

      implicit none
      real(kind=8), intent(in):: co2st
      real(kind=8)            :: co2aq

!     Convert dissolved CO2 from mol / m**3 to umol / L and prevent from division by zero
      co2aq = max(1.d-8, co2st * 1000.)

!     Fractionation wrt carbon-13
      alpha_p_13     = max(1., 1. + 0.001 * (17.6 * (1 - 2.02 / co2aq)))
      alpha_p_dia_13 = alpha_p_13 

!     Fractionation wrt carbon-14
      alpha_p_14     = 2. * alpha_p_13 - 1.
      alpha_p_dia_14 = 2. * alpha_p_dia_13 - 1.

    return
    end subroutine recom_ciso_photo
!   ----------------------------------------------------------------------------------  
 
 
    function lat_zone(lat_n)
!   ----------------------------------------------------------------------------------
!   Assign latitude zones from nodal latitude values 
!   ----------------------------------------------------------------------------------
    
      implicit none
      integer                  :: lat_zone
    
!     Input: Latitude value corresponding to node n
      real(kind=8), intent(in) :: lat_n

!     Binning of latitudes to three zones
      if (lat_n > 30.)  then       ! Northern Hemisphere polewards of 30°N
        lat_zone = 1
      else if (lat_n <- 30.) then  ! Southern Hemisphere polewards of 30°S
        lat_zone = 3
      else                         ! (Sub-) Tropical zone
        lat_zone = 2
      end if
    
      return
    end function lat_zone


    function wind_10(windstr_x, windstr_y)
!   ----------------------------------------------------------------------------------
!    computes wind speed at 10 m height "wind10" from wind stress fields tau_x, tau_y
!    as long as wind10 is not properly passed from ECHAM in coupled simulations.
!    We follow Peixoto & Oort (1992, Eq. (10.28), (10,29)) and Charnock (1955); 
!    also see MPI report 349 (2003), Eq. (5.7).
!   ----------------------------------------------------------------------------------
      implicit none
     
      real(kind=8) :: wind_10

!     Input
      real(kind=8), intent(in) :: windstr_x, windstr_y

!     Internal variables and parameters
!     Zonal and meridional velocities at 10 m height
      real(kind=8) :: u_10, v_10
!     Zonal and meridional friction velocities
      real(kind=8) :: u_fric, v_fric
!     Zonal and meridional roughness lengths
      real(kind=8) :: l_rough_x, l_rough_y
!     Inverse von-Karman constant (0.4), Charnock constant (0.018) divided by g, inverse density of air (1.3), log(10)
      real(kind=8), parameter :: inv_karm = 2.5, charn_g = 0.00173, inv_dens_air = 0.76923, log_10 = 2.30258
     
!     Calculate friction velocities (Peixoto & Oort, 1992, Eq. (10.28))
      u_fric = sqrt(abs(windstr_x) * inv_dens_air)
      v_fric = sqrt(abs(windstr_y) * inv_dens_air)

!     Calculate roughness lengths (MPI report 349, 2003, Eq. (5.7), quoting Charnock, 1955)
      l_rough_x = max((charn_g * u_fric**2), 1.5e-5)
      l_rough_y = max((charn_g * v_fric**2), 1.5e-5)

!     Calculate wind speed at 10 m (Peixoto & Oort, 1992, Eq. (10.29))
      u_10 = inv_karm * u_fric * (log_10 - log(l_rough_x))
      v_10 = inv_karm * v_fric * (log_10 - log(l_rough_y))
     
      wind_10 = sqrt(u_10**2 + v_10**2)
     
      return
    end function wind_10
!   ----------------------------------------------------------------------------------
   
   
end module REcoM_ciso
