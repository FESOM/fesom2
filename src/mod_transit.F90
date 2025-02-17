!==========================================================
MODULE mod_transit
! Parameters, variables and functions for transient tracer simulations.
! By mbutzin, 2019-2021

  implicit none
  save

! Atmospheric pressure, local (dummy) variable, global-mean SLP, and local wind speed at 10 m
  real(kind=8) :: press_a, mean_slp = 1.01325e5, wind_2

! Atmospheric trace gas (dummy) values used in air-sea flux calculations
! Isotopic ratios are normalized and fractionation-corrected,
! volume mixing ratios are mole fractions in dry air.
  real(kind=8) :: r14c_a  = 1.0, &       ! 14CO2 / 12CO2 (may vary with latitude in transient runs)
                  r39ar_a = 1.0, &       ! 39Ar / 40 Ar (homogeneous)
                  xarg_a  = 9.34e-3, &   ! Argon (homogeneous)
                  xCO2_a  = 284.32e-6, & ! CO2 (CMIP6 & OMIP-BGC: 284.32e-6 for 1700-1850, PMIP4: 190.00e-6 for 21 ka BP)
                  xf11_a  = 0.0, &       ! CFC-11 (latitude dependent)
                  xf12_a  = 0.0, &       ! CFC-12 (latitude dependent)
                  xsf6_a  = 0.0          ! SF6 (latitude dependent)

! Transient values of atmospheric trace gases (1d-arrays of variable length to be specified in namelist.config -> length_transit)
  real(kind=8), allocatable, dimension(:) :: r14c_nh, r14c_tz, r14c_sh, & ! 14CO2 / 12CO2, latitude-dependent (e.g., bomb 14C)
                                             r14c_ti, &                   ! 14CO2 / 12CO2, homogenous (e.g., IntCal)
                                             xCO2_ti, &                   ! CO2
                                             xf11_nh, xf11_sh, &          ! CFC-11, latitude-dependent
                                             xf12_nh, xf12_sh, &          ! CFC-12, latitude-dependent
                                             xsf6_nh, xsf6_sh             ! SF6, latitude-dependent
  integer, allocatable, dimension(:)      :: year_ce                      ! current year in anthropenic runs (control output)
  integer                                 :: length_transit = 1, &        ! length (years) of transient tracer input
                                             ti_start_transit = 1         ! index of the first tracer input year in ifile_transit
  logical                                 :: l_r14c = .false., &          ! switch on / off specific tracerss
                                             l_r39ar = .false. , &
                                             l_f11 = .false., &
                                             l_f12 = .false., &
                                             l_sf6 = .false., &
                                             anthro_transit = .false., &  
                                             paleo_transit = .false.      ! specify tracer input scenario
  character(300)                          :: ifile_transit ='Table_CO2_isoC_CFCs_SF6.txt'! tracer input file; not neccessary for steady state simulations


! Parameters which can be changed via namelist.oce (-> transit_param)
! Global-mean concentrations of DIC and Argon in the mixed layer (mol / m**3)
  real(kind=8) :: dic_0 = 2.00, &        ! GLODAPv2, 0-50 m: TCO2 ~ 2050 umol / kg
                  arg_0 = 0.01           ! Hamme et al. 2019, doi:10.1146/annurev-marine-121916-063604
                  
! Radioactive decay constants (1 / s; default values assume that 1 year = 365.00 days)
  real(kind=8) :: decay14 = 3.8561e-12 , & ! 14C; t1/2 = 5700 a following OMIP-BGC
                  decay39 = 8.1708e-11     ! 39Ar; t1/2 = 269 a

! Further internal parameters
! Latitude of atmospheric boundary conditions and latitudinal interpolation weight
  real(kind=8) :: y_abc, yy_nh
! Tracer indices of transient tracers
  integer ::      id_r14c, id_r39ar, id_f11, id_f12, id_sf6, index_transit_r14c, index_transit_r39ar, index_transit_f11, index_transit_f12, index_transit_sf6
! Time index (=year) in transient simulations
  integer ::      ti_transit

! Namelist to modify default parameter settings
  namelist / transit_param / &
                             l_r14c, &                     ! switch on R14C
                             l_r39ar, &                    ! switch on R39Ar
                             l_f11, &                      ! switch on CFC-11
                             l_f12, &                      ! switch on CFC-12
                             l_sf6, &                      ! switch on SF6
                             anthro_transit, &             ! anthropogenic transient tracers
                             paleo_transit, &              ! paleo transient tracers
                             length_transit, &             ! 166 if (anthro_transit=.true.)
                             ti_start_transit, &           ! 1 for D14C, 80 for CFC-12
                             ifile_transit, &              ! forcing file
                             r14c_a, &                     ! atm. 14C/C ratio, global mean
                             r39ar_a, &                    ! atmospheric 39Ar/Ar ratio
                             xarg_a, &                     ! atmospheric mole fraction of Argon
                             xco2_a, &                     ! atmospheric mole fraction of CO2
                             dic_0, arg_0, &               ! mixed layer values of DIC and Argon
                             decay14, decay39              ! decay constants of 14C and 39Ar

  contains

    function iso_flux(which_gas, temp_c, sal, wind_2, f_ice, p_atm, x_gas, r_air, r_sea, c_surf)
!     Calculate isotopic air-sea exchange fluxes in 1 / (m**2 * s) assuming local solubility equilibrium
!     for the abundant isotopologue. Positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: iso_flux
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! trace gas name

      real(kind=8), intent(in) :: temp_c, sal, &   ! SST (deg C) and SSS ("PSU" or permil)
                                  wind_2, &        ! wind speed at 10 m heigth squared
                                  f_ice, &         ! sea-ice fractional coverage
                                  p_atm, &         ! total atmospheric pressure (Pa)
                                  x_gas, &         ! atmospheric mole fraction of the abundant isotope
                                  r_air, r_sea, &  ! isotopic ratios in atmosphere and ocean
                                  c_surf           ! surface water concentration of the abundant isotope (mol / m**3)

      iso_flux = transfer_vel(which_gas, temp_c, wind_2) * &
                 solub(which_gas, temp_c, sal) * p_atm / 1.01325e5 * x_gas * &
                 (r_air - r_sea) * (1. - f_ice) / c_surf
      return
    end function iso_flux


    function gas_flux(which_gas, temp_c, sal, wind_2, f_ice, p_atm, x_gas, c_surf)
!     Computes air-sea exchange gas fluxes in mol / (m**2 * s) , positive values mean oceanic uptake.
      implicit none

      real(kind=8) :: gas_flux
!     Input parameters
      character(len=3), intent(in) :: which_gas  ! trace gas name
      real(kind=8), intent(in) :: temp_c, sal, & ! SST (deg C) and SSS ("PSU" or permil)
                                  wind_2,     &  ! wind speed at 10 m heigth squared
                                  f_ice, &       ! sea-ice fractional coverage
                                  p_atm, &       ! total atmospheric pressure (Pa)
                                  x_gas, &       ! atmospheric mole fraction
                                  c_surf         ! marine surface water concentration (mol / m**3)
!     Internal variables
      real(kind=8) :: c_sat                      ! marine saturation concentration (mol / m**3)
      c_sat = solub(which_gas, temp_c, sal) * p_atm / 1.01325e5 * x_gas
      gas_flux = transfer_vel(which_gas, temp_c, wind_2) * (c_sat - c_surf) * (1. - f_ice)

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
      case ("f11")
!       CFC-11 in mol / (L * atm) (Warner & Weiss 1985, doi:10.1016/0198-0149(85)90099-8, Table 5)
        a1 = -229.9261;  a2 = 319.6552;   a3 = 119.4471;   a4 = -1.39165; pow = 2
        b1 = -0.142382;  b2 = 0.091459;   b3 = -0.0157274; c1 = 0.
        con2con = 1000. ! convert to mol / (m**3 * atm)
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
      case ("f11") ! CFC-11
        as = 3579.2; bs = -222.63; cs = 7.5749; ds = -0.145950; es = 0.0011870
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


    function transfer_vel(which_gas, temp_c, wind_2)
!     Compute gas transfer velocities of / for tracers
      implicit none
!     Result
      real(kind=8) :: transfer_vel                 ! transfer velocity (m / s)
!     Input parameters
      character(len=3), intent(in) :: which_gas    ! tracer name
      real(kind=8), intent(in) :: temp_c, &        ! temperature (deg C)
                                  wind_2           ! wind speed squared at 10 m height (m / s)

!     Wanninkhof (2014), eq. (4) with a = 0.251 (cm / h) / (m / s)**2 -> 6.9722e-7 s / m
!     to obtain the gas transfer velocity in m / s
      transfer_vel = 6.9722e-7 * sc_660(which_gas, temp_c)**(-0.5) * wind_2

      return
    end function transfer_vel


    function speed_2(windstr_x, windstr_y)
!     Computes the square of wind speed at 10 m height from wind stress fields
!     in coupled simulations as long as it is not provided by the AGCM / OASIS.
!     We follow Peixoto & Oort (1992, Eq. (10.28), (10,29)) and Charnock (1955);
!     also see MPI report 349 (2003), Eq. (5.7).
      implicit none
      real(kind=8) :: speed_2

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

      speed_2 = u_10**2 + v_10**2

      return
    end function speed_2


    subroutine read_transit_input(ifileunit)
!   Read atmospheric input of isoCO2 and / or other tracers
      implicit none

!     Internal variables
      integer, intent(in) :: ifileunit
      integer :: jj
      real(kind=8), allocatable, dimension(:) :: d14c_nh, d14c_tz, d14c_sh, d14c_ti, d13c_dummy

      if (anthro_transit) then
!       Anthropogenic input for 1850 - 2015 CE
        allocate(d14c_nh(length_transit))
        allocate(d14c_tz(length_transit))
        allocate(d14c_sh(length_transit))
        allocate(r14c_nh(length_transit))
        allocate(r14c_tz(length_transit))
        allocate(r14c_sh(length_transit))
        allocate(xCO2_ti(length_transit))
        allocate(xf11_nh(length_transit))
        allocate(xf11_sh(length_transit))
        allocate(xf12_nh(length_transit))
        allocate(xf12_sh(length_transit))
        allocate(xsf6_nh(length_transit))
        allocate(xsf6_sh(length_transit))
        allocate(year_ce(length_transit))
        allocate(d13c_dummy(length_transit))

!       Skip header lines
        do jj = 1,15
          read (ifileunit, fmt=*)
        end do
!       Read input values
        do jj = 1, length_transit
          read (ifileunit, fmt=*) year_ce(jj), &
                           xCO2_ti(jj), &
                           d14c_nh(jj), d14c_tz(jj), d14c_sh(jj), &
                           d13c_dummy(jj), &
                           xf11_nh(jj), xf11_sh(jj), &
                           xf12_nh(jj), xf12_sh(jj), &
                           xsf6_nh(jj), xsf6_sh(jj)
        end do

!       Convert Delta14C to F14C
        r14c_nh = 1. + 0.001 * d14c_nh
        r14c_tz = 1. + 0.001 * d14c_tz
        r14c_sh = 1. + 0.001 * d14c_sh
!       Convert volume mixing ratios
        xCO2_ti = xCO2_ti * 1.e-6
        xf11_nh = xf11_nh * 1.e-12
        xf11_sh = xf11_sh * 1.e-12
        xf12_nh = xf12_nh * 1.e-12
        xf12_sh = xf12_sh * 1.e-12
        xsf6_nh = xsf6_nh * 1.e-12
        xsf6_sh = xsf6_sh * 1.e-12
      elseif (paleo_transit) then
!       UNDER CONSTRUCTION
        allocate(d14c_ti(length_transit))
        allocate(r14c_ti(length_transit))
        allocate(xCO2_ti(length_transit))
!
      else
!       Read constant parameter values from namelist.oce.
!       This is done in subroutine gen_model_setup.
        allocate(d14c_ti(1))
        allocate(r14c_ti(1))
        allocate(xCO2_ti(1))
      end if

      return
      end subroutine read_transit_input

END MODULE mod_transit
!==========================================================
