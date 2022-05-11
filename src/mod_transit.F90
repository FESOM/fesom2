!==========================================================
MODULE mod_transit
! Parameters, variables and functions for transient tracer simulations.
! By mbutzin, 2019-2021

  implicit none
  save

! Normalized and fractionation-corrected atmospheric 14CO2 / 12CO2 ratios
  real(kind=8) :: r14c_a  = 1.0, & ! Value passed in air-sea flux calculation
                  r14c_nh = 1.0, & ! Northern Hemisphere
                  r14c_tz = 1.0, & ! Tropics
                  r14c_sh = 1.0    ! Southern Hemisphere
! Normalized and fractionation-corrected atmospheric 39Ar/40Ar ratio
  real(kind=8) :: r39ar_a  = 1.0   ! Global average and/or value in air-sea flux calculation
! Atmospheric CO2 concentration (mole fraction in dry air)
! CMIP6 & OMIP-BGC: xCO2_a = 284.32 ppm for 1700-1850 CE
! PMIP4:            xCO2_a = 190.00 ppm for 21 kcal BP
  real(kind=8) :: xCO2_a = 284.32e-6
! Atmospheric concentrations of CFC-12 and SF6 (ppt in dry air)
  real(kind=8) :: xf12_a  = 0.00, &  ! CFC-12, value passed in air-sea flux calculations
                  xf12_nh = 0.00, &  ! CFC-12, Northern Hemisphere
                  xf12_sh = 0.00, &  ! CFC-12, Southern Hemisphere
                  xsf6_a  = 0.00, &  ! SF6, value passed in air-sea flux calculations
                  xsf6_nh = 0.00, &  ! SF6, Northern Hemisphere
                  xsf6_sh = 0.00     ! SF6, Southern Hemisphere
! Atmospheric concentration trends of atmospheric CFC-12 and SF6 (ppt / year)
  real(kind=8) :: f12t_a  = 0.00, &  ! CFC-12, value passed in air-sea flux calculations
                  f12t_nh = 0.00, &  ! CFC-12, Northern Hemisphere
                  f12t_sh = 0.00, &  ! CFC-12, Southern Hemisphere
                  sf6t_a  = 0.00, &  ! SF6, value passed in air-sea flux calculations
                  sf6t_nh = 0.00, &  ! SF6, Northern Hemisphere
                  sf6t_sh = 0.00     ! SF6, Southern Hemisphere
! Atmospheric Argon concentration (mole fraction in dry air)
  real(kind=8) :: xarg_a  = 9.34e-3      ! value passed in air-sea flux calculation
! Global-mean concentrations of DIC and Argon in the mixed layer (mol / m**3)
  real(kind=8) :: dic_0 = 2.00, &        ! GLODAPv2, 0-50 m: TCO2 ~ 2050 umol / kg
                  arg_0 = 0.01           ! Hamme et al. 2019, doi:10.1146/annurev-marine-121916-063604
! Radioactive decay constants (1 / s; default values assume that 1 year = 365.00 days)
  real(kind=8) :: decay14 = 3.8561e-12 , & ! 14C; t1/2 = 5700 a following OMIP-BGC
                  decay39 = 8.1708e-11     ! 39Ar; t1/2 = 269 a
! Latitude of atmospheric boundary conditions and latitudinal interpolation weight
  real(kind=8) :: y_abc, yy_nh
! Tracer indices of transient tracers
  integer ::      id_r14c, id_r39ar, id_f12, id_sf6
  
! Namelist to modify default parameter settings
  namelist / transit_param / r14c_nh, r14c_tz, r14c_sh, &  ! atmospheric F14C
                             r39ar_a, &                    ! atmospheric 39Ar/Ar ratio
                             xarg_a, &                     ! atmospheric mole fraction of Argon
                             xco2_a, &                     ! atmospheric mole fraction of CO2
                             xf12_nh, xf12_sh, &           ! atmospheric mole fractions of CFC-12
                             xsf6_nh, xsf6_sh, &           ! atmospheric mole fractions of SF6
                             f12t_nh, f12t_sh, &           ! atmospheric trends of CFC-12
                             sf6t_nh, sf6t_sh, &           ! atmospheric trends of SF6
                             dic_0, arg_0, &               ! mixed layer values of DIC and Argon
                             decay14, decay39              ! decay constants of 14C and 39Ar


  contains


    function iso_flux(which_gas, temp_c, sal, u_10, v_10, f_ice, p_atm, x_gas, r_air, r_sea, c_surf)
!     Calculate isotopic air-sea exchange fluxes in 1 / (m**2 * s) assuming local solubility equilibrium
!     for the abundant isotopologue. Positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: iso_flux
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! trace gas name

      real(kind=8), intent(in) :: temp_c, sal, &   ! SST (deg C) and SSS ("PSU" or permil)
                                  u_10, v_10,  &   ! wind speed at 10 m height (m / s)
                                  f_ice, &         ! sea-ice fractional coverage
                                  p_atm, &         ! total atmospheric pressure (Pa)
                                  x_gas, &         ! atmospheric mole fraction of the abundant isotope
                                  r_air, r_sea, &  ! isotopic ratios in atmosphere and ocean
                                  c_surf           ! surface water concentration of the abundant isotope (mol / m**3)

      iso_flux = transfer_vel(which_gas, temp_c, u_10, v_10) * &
                 solub(which_gas, temp_c, sal) * p_atm / 1.01325e5 * x_gas * &
                 (r_air - r_sea) * (1. - f_ice) / c_surf
      return
    end function iso_flux


    function gas_flux(which_gas, temp_c, sal, u_10, v_10, f_ice, p_atm, x_gas, c_surf)
!     Computes air-sea exchange gas fluxes in mol / (m**2 * s) , positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: gas_flux
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! trace gas name
      real(kind=8), intent(in) :: temp_c, sal, & ! SST (deg C) and SSS ("PSU" or permil)
                                  u_10, v_10,  & ! wind speed at 10 m height (m / s)
                                  f_ice, &       ! sea-ice fractional coverage
                                  p_atm, &       ! total atmospheric pressure (Pa)
                                  x_gas, &       ! atmospheric mole fraction 
                                  c_surf         ! marine surface water concentration (mol / m**3)
!     Internal variables
      real(kind=8) :: c_sat                      ! marine saturation concentration (mol / m**3)
      c_sat = solub(which_gas, temp_c, sal) * p_atm / 1.01325e5 * x_gas
      gas_flux = transfer_vel(which_gas, temp_c, u_10, v_10) * (c_sat - c_surf) * (1. - f_ice)

      return
    end function gas_flux


    function solub(which_gas, temp_c, sal)
!     Computes the solubility of trace gases in seawater.
!     This parametrization includes the effect of water vapor.
      implicit none
      real(kind=8) :: solub                  ! solubility ((p)mol / (m**3 * atm))
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! tracer name
      real(kind=8), intent(in) :: temp_c, &      ! temperature (deg C)
                                  sal            ! salinity ("PSU" or permil)
      real(kind=8) :: a1, a2, a3, a4, &          ! polynomial coefficients of the
                      b1, b2, b3, b4, c1, &      ! solubility function
                      temp_k100, &               ! water temperature in K / 100
                      con2con                    ! concentration units conversion factor
      integer ::      pow                        ! power in solubility function

      temp_k100 = (temp_c + 273.15) * 0.01

      select case (which_gas)
      case ("co2")
!       CO2 in mol / (L * atm) (Weiss & Price 1985, doi:10.1016/0304-4203(80)90024-9, Table VI) 
        a1 = -160.7333;  a2 = 215.4152;   a3 = 89.8920;   a4 = -1.47759;  pow = 2
        b1 =  0.029941;  b2 = -0.027455;  b3 = 0.0053407; c1 = 0.
        con2con = 1000.  ! convert to mol / (m**3 * atm)
      case ("f12") 
!       CFC-12 in mol / (L * atm) (Warner & Weiss 1985, doi:10.1016/0198-0149(85)90099-8, Table 5)
        a1 = -218.0971;  a2 = 298.9702;   a3 = 113.8049;   a4 = -1.39165; pow = 2
        b1 = -0.143566;  b2 = 0.091015;   b3 = -0.0153924; c1 = 0.
        con2con = 1000. ! convert to mol / (m**3 * atm)
      case ("sf6") 
!       SF6 in mol / (L * atm) (Bullister et al. 2002, doi:10.1016/S0967-0637(01)00051-6, Table 3)
        a1 = -80.0343;   a2 = 117.232;    a3 = 29.5817;    a4 = 0.;       pow = 2
        b1 =  0.0335183; b2 = -0.0373942; b3 = 0.00774862; c1 = 0.
        con2con = 1000. ! convert to mol / (m**3 * atm)
      case("arg")
!       Argon in mol / kg (Jenkins et al. 2019, doi:10.1016/j.marchem.2019.03.007, Table 4)
        a1 = -227.4607; a2 = 305.4347;   a3 = 180.5278;   a4 = -27.99450; pow = 1
        b1 = -0.066942; b2 = 0.037201;   b3 = -0.0056364; c1 = -5.30e-6
        con2con = 1024.5  ! convert to mol / m**3 assuming homogeneous density of surface water
      end select

      solub = exp(       a1 + a2 / temp_k100 + a3 * log(temp_k100) + a4 * temp_k100 **pow + & 
                  sal * (b1 + b2 * temp_k100 + b3 * temp_k100**2   + c1 * sal))
      solub = solub * con2con

      return
    end function solub


    function sc_660(which_gas, temp_c)
!     Schmidt numbers of trace gases in sea water with S = 35 
!     normalized to 20 degC (Sc(CO2) ~660; Wanninkhof 2014, tab. 1)).
      implicit none

!     Result
      real(kind=8) :: sc_660                       ! Schmidt number
!     Input parameters
      character(len=3), intent(in) :: which_gas    ! tracer name
      real(kind=8), intent(in) :: temp_c           ! temperature (deg C)
!     Internal parameters and/or variables
      real(kind=8) :: as, bs, cs, ds, es           ! polynomial coefficients

      select case (which_gas)
      case ("co2") ! CO2
        as = 2116.8; bs = -136.25; cs = 4.7353; ds = -0.092307; es = 0.0007555
      case ("f12") ! CFC-12
        as = 3828.1; bs = -249.86; cs = 8.7603; ds = -0.171600; es = 0.0014080
      case ("sf6") ! SF6
        as = 3177.5; bs = -200.57; cs = 6.8865; ds = -0.133350; es = 0.0010877
      case ("arg") ! Ar-39
        as = 2078.1; bs = -146.74; cs = 5.6403; ds = -0.118380; es = 0.0010148
      end select
      
      sc_660 = (as + bs *temp_c + cs * temp_c**2 + ds * temp_c**3 + es * temp_c**4) / 660.
      
      return
    end function sc_660


    function transfer_vel(which_gas, temp_c, u_10, v_10)
!     Compute gas transfer velocities of / for tracers
!     Result
      real(kind=8) :: transfer_vel                 ! transfer velocity (m / s)
!     Input parameters
      character(len=3), intent(in) :: which_gas    ! tracer name
      real(kind=8), intent(in) :: temp_c, &        ! temperature (deg C)
                                  u_10, v_10       ! wind speed at 10 m height (m / s)

!     Wanninkhof (2014), eq. (4) with a = 0.251 (cm / h) / (m / s)**2 -> 6.9722e-7 s / m 
!     to obtain the gas transfer velocity in m / s
      transfer_vel = 6.9722e-7 * sc_660(which_gas, temp_c)**(-0.5) * (u_10**2 + v_10**2) 

      return
    end function transfer_vel

END MODULE mod_transit
!==========================================================