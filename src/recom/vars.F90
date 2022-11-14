!> \file vars.f90
!! \BRIEF 
!> Module with vars subroutine - compute carbonate system vars from DIC,Alk,T,S,P,nuts
MODULE mvars
CONTAINS
!>    Computes standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
SUBROUTINE vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                      &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon, verbose    )

  !   Purpose:
  !     Computes other standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     as 1D arrays
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations (all 1-D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = also be used to convert absolute to practical salinity, when optS='Sabs' 
  !             = dummy array (unused when optP='db' and optS='Sprc')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !             = conservative temperature [degrees C] (with optT='Tcsv')
  !     sal     = practical [psu] or absolute [g/kg] salinity
  !     alk     = total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =               [eq/kg]  with optCON = 'mol/kg'
  !     dic     = dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = silica    [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     phos    = phosphate [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input & output concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC, ALK, sil, & phos given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC, ALK, sil, & phos given in mol/m^3 (std MODEL units)
  !     -----------
  !     optT: choose in-situ, potential or conservative temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tcsv' means input is Conservative Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (thus p = depth)
  !     ---------
  !     optB: choose total boron formulation - Uppström (1974) vs. Lee et al. (2010)
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
  !     ---------
  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     ----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2
  !     ---------
  !       PRESSURE corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ----------
  !     optS: choose practical [psu] or absolute [g/kg] salinity as input
  !     ----------
  !       -> 'Sprc' means input is practical salinity according to EOS-80 convention
  !       -> 'Sabs' means input is absolute salinity according to TEOS-10 convention (practical sal. will be computed)
  !     ---------
  !     lon:  longitude in degrees East
  !     ----------
  !        Optional, it may be used along with latitude when optS is "Sabs".
  !        Then, they are parameters for conversion from Absolute to Practical Salinity.
  !
  !        When seawater is not of standard composition, Practical Salinity alone is not sufficient 
  !        to compute Absolute Salinity and vice-versa. One needs to know the chemical composition, 
  !        mainly silicate and nitrate concentration. When those concentrations are unknown and 'lon' and 'lat' 
  !        are given, absolute salinity conversion is based on WOA silicate concentration at given location. 
  !
  !        Alternative when optS is 'Sabs' :
  !        -------------------------------
  !        When silicate and phosphate concentrations are known, nitrate concentration is inferred from phosphate
  !        (using Redfield ratio), then practical salinity is computed from absolute salinity, 
  !        total alcalinity (alk), DIC (dic), silicate (sil) and phosphate (phos).
  !        In that case, do not pass optional parameter 'lon'.
  !
  !        When neither chemical composition nor location are known, an arbitrary geographic point is chosen:
  !        mid equatorial Atlantic. Note that this implies an error on computed practical salinity up to 0.02 psu.
  !        In that case, do pass parameter 'lon' and set each of its elements to 1.e20.
  !     ---------
  !     verbose: turn on or off the print statements  (.true. or .false ; default is .true.)
  !     ----------
  !
  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state
  !     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  !     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  !     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  !     tempis  = in-situ temperature [degrees C]

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif
  USE msingledouble
  USE mconstants
  USE mp80
  USE mrho
  USE meos
  USE msw_temp
  USE mvarsolver

  IMPLICIT NONE

! Input variables
  !> N: number of records
!f2py intent(hide) n
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: temp
  !> salinity <b>[psu] or [g/kg]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: phos
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: lat

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature
  !>  \b 'Tcsv" for conservative temperature (in two last cases, in-situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) [default] or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) [default] or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(*), OPTIONAL, INTENT(in) :: optGAS
  !> choose \b 'Sprc' for practical sal. (EOS-80, default) [default] or \b 'Sabs' for absolute salinity (TEOS-10)
!f2py character(4) optional, intent(in) :: optS = 'Sprc'
  CHARACTER(4), OPTIONAL, INTENT(in) :: optS
  !> longitude <b>[degrees east]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lon = -25.
  REAL(kind=rx), OPTIONAL, INTENT(in), DIMENSION(N) :: lon
  !> to print warnings when input out of bounds, use .true.; for no warnings, use .false.
!f2py logical optional, intent(in) :: verbose = .true.
  LOGICAL, OPTIONAL, INTENT(in) :: verbose

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: tempis

! Local variables
  ! practical salinity (psu)
  REAL(kind=rx), DIMENSION(N) :: salprac

  
  ! Call the subroutine that actually computes
  call vars_sprac (ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                         &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon, salprac, verbose   )
                
END SUBROUTINE vars


!>    This is the subroutine that does the actual computations
!!    This subroutine is only called internaly (i.e. by other Mocsy subroutines)
!!    Its output parameter is "Practical Salinity", when Absolute Salinity is passed in,
!!    is used by those internal calling routines.
!!
SUBROUTINE vars_sprac (ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                             &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon, salprac, verbose    )

  USE msingledouble
  USE mconstants
  USE mp80
  USE mrho
  USE meos
  USE gsw_mod_toolbox, only: gsw_t_from_ct, gsw_ct_from_t, gsw_rho
  USE msw_temp
  USE mvarsolver

  IMPLICIT NONE

! Input variables
  !>     number of records
!f2py intent(hide) n
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: temp
  !> salinity <b>[psu] or [g/kg]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: phos
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: lat

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature
  !>  \b 'Tcsv" for conservative temperature (in two last cases, in-situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) [default] or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='u74'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) [default] or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) [default] or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(*), OPTIONAL, INTENT(in) :: optGAS
  !> choose \b 'Sprc' for practical sal. (EOS-80, default) or \b 'Sabs' for absolute salinity (TEOS-10)
!f2py character*4 optional, intent(in) :: optS='Sprc'
  CHARACTER(4), OPTIONAL, INTENT(in) :: optS
  !> longitude <b>[degrees east]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lon = -25.
  REAL(kind=rx), OPTIONAL, INTENT(in), DIMENSION(N) :: lon
!f2py logical optional, intent(in) :: verbose=.true.
  LOGICAL, OPTIONAL, INTENT(in) :: verbose

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: tempis
  !> practical salinity \b <b>[psu]</b>
  REAL(kind=rx), OPTIONAL, INTENT(out), DIMENSION(N) :: salprac

! Local variables
  REAL(kind=rx) :: ssal, salk, sdic, ssil, sphos
  REAL(kind=r8) :: tempot, tempis68, tempot68, tempis90, tempcsv
  REAL(kind=r8) :: drho

  ! local 1-long array version of scalar variables
  REAL(kind=r8), DIMENSION(1) :: aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(kind=r8), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(kind=r8), DIMENSION(1) :: aSt, aFt, aBt

  REAL(kind=rx), DIMENSION(1) :: sabs1, spra1, p1, lon1, lat1
  REAL(kind=rx), DIMENSION(1) :: tc1, ta1, sit1, nt1
  REAL(kind=rx), DIMENSION(1) :: sal1
  
  REAL(kind=r8) :: Patmd
  REAL(kind=r8) :: Ptot
  REAL(kind=r8) :: Phydro_atm

  INTEGER :: i, icount

  REAL(kind=r8) :: prb
  REAL(kind=r8) :: s
  REAL(kind=r8) :: tc, ta
  REAL(kind=r8) :: sit, pt
  
  REAL(kind=r8), DIMENSION(2) :: dicdel, pco2del
  REAL(kind=r8) :: dx, Rf
  REAL(kind=r8) :: dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC

  INTEGER :: kcomp
  INTEGER :: j, minusplus

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS
  LOGICAL      :: verbosity

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
!   Default is Uppstrom (1974) for total boron
    opB = 'u74'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF
  IF (PRESENT(optS)) THEN
    opS = optS
  ELSE
    opS = 'Sprc'
  ENDIF
  IF (PRESENT(verbose)) THEN
    verbosity = verbose
  ELSE
    verbosity = .true.
  ENDIF

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!    - depth [m] => convert to pressure [db]
!    - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     !print *,"optP =", optP, "end"
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
!        write(*,*), 'i,depth,lat,p',i,depth(i),lat(i),p(i)
        p(i) = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p(i) = depth(i)
     ELSE
        !print *,"optP =", optP, "end"
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = DBLE(temp(i))
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002_rx) / 0.99975_rx
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p(i), 0.d0 )
!jh        tempis68 = sw_temp(sal(i), SGLE(tempot68), p(i), SGLE(0.d0) )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis90 = 0.99975*tempis68 + 0.0002_r8
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis90 = DBLE(temp(i))
        tempis68  = (tempis90 - 0.0002_r8) / 0.99975_r8
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
!       dtempot   = 0.99975*dtempot68 + 0.0002
     ELSEIF (trim(optT) == 'Tcsv' .OR. trim(optT) == 'tcsv') THEN
!       Convert given conservative temperature to in-situ temperature
        ! First convert salinity to absolute sal., if necessary
        IF (trim(opS) == 'Sprc')  THEN
            ! conversion will use default geographic location
            spra1(1) = sal(i)
            p1(1) = p(i)
            CALL sp2sa_geo (spra1, 1, sabs1, p1)
        ELSE
            sabs1(1) = sal(i)
        END IF
        ! Then convert temperature
        tempis90 = gsw_t_from_ct (DBLE(sabs1(1)), DBLE(temp(i)), DBLE(p(i)))
        tempis68  = (tempis90 - 0.0002_r8) / 0.99975_r8
     ELSE
        PRINT *,"optT must be either 'Tpot, 'Tinsitu' or 'Tcsv'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF
     tempis(i) = SGLE(tempis90)

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (verbosity .EQV. .true.) THEN
            IF (       sal(i) < 0.0_rx  &
                 .OR.  alk(i) < 0.0_rx  &
                 .OR.  dic(i) < 0.0_rx  &
                 .OR.  sil(i) < 0.0_rx  &
                 .OR. phos(i) < 0.0_rx  &
                 .OR.  sal(i) > 1e+3_rx &
                 .OR.  alk(i) > 1e+3_rx &
                 .OR.  dic(i) > 1e+3_rx &
                 .OR.  sil(i) > 1e+3_rx &
                 .OR. phos(i) > 1e+3_rx) THEN
               PRINT *, 'i, icount, tempot, sal,    alk,    dic,    sil,    phos =', &
                         i, icount, tempot, sal(i), alk(i), dic(i), sil(i), phos(i)
            ENDIF
        ENDIF
!       Zero out any negative salinity, phosphate, silica, dic, and alk
        IF (sal(i) < 0.0_rx) THEN
           ssal = 0.0_rx
        ELSE
           ssal = sal(i)
        ENDIF
        IF (phos(i) < 0.0_rx) THEN
           sphos = 0.0_rx
        ELSE
           sphos = phos(i)
        ENDIF
        IF (sil(i) < 0.0_rx) THEN
           ssil = 0.0_rx
        ELSE
           ssil = sil(i)
        ENDIF
        IF (dic(i) < 0.0_rx) THEN
          sdic = 0.0_rx
        ELSE
          sdic = dic(i)
        ENDIF
        IF (alk(i) < 0.0_rx) THEN
          salk = 0.0_rx
        ELSE
          salk = alk(i)
        ENDIF

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))
!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p(i) / SGLE(10.0d0))
        Phydro_atm = prb / 1.01325d0  ! convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
!       Total pressure [atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           Ptot = Patmd + Phydro_atm   ! total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF

!       Compute in-situ density [kg/m^3]
        ! If Absolute salinity is given
        IF (trim(opS) == 'Sabs')  THEN
           ! If in-situ or potential temperature is given
           IF (trim(optT) /= 'Scsv') THEN
              ! First compute conservative temperature
              tempcsv = gsw_ct_from_t (DBLE(ssal), tempis90, DBLE(p(i)))
           ELSE
              tempcsv = DBLE(temp(i))
           ENDIF
           rhoSW(i) = gsw_rho(DBLE(ssal), tempcsv, prb)
        ELSE
           rhoSW(i) = rho(ssal, tempis68, prb)
!jh           rhoSW(i) = rho(ssal, SGLE(tempis68), SGLE(prb))
        ENDIF

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (trim(optCON) == 'mol/kg') THEN
!          No conversion:  drho = 1.
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = 1.0
        ELSEIF (trim(optCON) == 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = DBLE(rhoSW(i))
        ELSE
           PRINT *,"optCON must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

!       Initialise ta, tc, pt and sit
        ta  = DBLE(salk)  / drho
        tc  = DBLE(sdic)  / drho
        pt  = DBLE(sphos) / drho
        sit = DBLE(ssil)  / drho

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Convert from Absolute to Practical salinity if needed
        IF (trim(opS) == 'Sabs')  THEN
           sabs1(1) = ssal   ! absolute sal.
           ! If longitude is passed in
           IF (PRESENT(lon)) THEN
               p1(1) = p(i)
               IF (lon(i) .NE. 1.e20_rx) THEN
                  ! longitude and latitude are defined
                  lon1(1) = lon(i)
                  lat1(1) = lat(i)
                  CALL sa2sp_geo (sabs1, 1, spra1, p1, lon1, lat1)
               ELSE
                  ! will use default geographic location
                  CALL sa2sp_geo (sabs1, 1, spra1, p1)
               ENDIF
           ELSE
               tc1(1) = tc
               ta1(1) = ta
               sit1(1) = sit
               ! Nitrate total : from Phosphate using Redfield ratio (16)
               nt1(1) = 16.0 * pt
               CALL sa2sp_chem(sabs1, ta1, tc1, nt1, sit1, 1, spra1)
           ENDIF
           s = DBLE(spra1(1))
           IF (PRESENT(salprac)) salprac(i) = SGLE(s)
        ENDIF
          
!       Get all equilibrium constants and total concentrations of SO4, F, B
        sal1(1) = SGLE(s)
        CALL constants (aK0, aK1, aK2, aKb, aKw, aKs, aKf,            &
                    aKspc, aKspa, aK1p, aK2p, aK3p, aKsi,             &
                    aSt, aFt, aBt,                                    &
                    temp(i), sal1, Patm(i),                           &
                    depth(i), lat(i), 1,                              &
                    optT, optP, opB, opK1K2, opKf, opGAS              )

!       Solve for pH and all other variables
!       ------------------------------------
        CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC,   &
                    tempis90, s, ta, tc, pt, sit,                                &
                    aBt(1), aSt(1), aFt(1),                                      &
                    aK0(1), aK1(1), aK2(1), aKb(1), aKw(1), aKs(1), aKf(1),      &
                    aKspc(1), aKspa(1), aK1p(1), aK2p(1), aK3p(1), aKsi(1),      & 
                    Patmd, prb, drho, opGAS                   )

!       Convert all output variables from double to single precision, if necessary
        pH(i)     = SGLE(dph)
        co2(i)    = SGLE(dco2)
        hco3(i)   = SGLE(dhco3)
        co3(i)    = SGLE(dco3)
        fCO2(i)   = SGLE(dfCO2)
        pCO2(i)   = SGLE(dpCO2)
        OmegaA(i) = SGLE(dOmegaA)
        OmegaC(i) = SGLE(dOmegaC)

!       Compute Revelle factor numerically (derivative using centered-difference scheme)
        DO j=1,2
           minusplus = (-1)**j
           dx = 0.1 * 1e-6         ! Numerical tests found for DIC that optimal dx = 0.1 umol/kg (0.1e-6 mol/kg)
           dicdel(j) = tc + DBLE(minusplus)*dx/2.0d0
           CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC, &
                tempis90, s, ta, dicdel(j), pt, sit,                         &
                aBt(1), aSt(1), aFt(1),                                      &
                aK0(1), aK1(1), aK2(1), aKb(1), aKw(1), aKs(1), aKf(1),      &
                aKspc(1), aKspa(1), aK1p(1), aK2p(1), aK3p(1), aKsi(1),      & 
                Patmd, prb, drho, opGAS                        )
            pco2del(j) = dpco2
        END DO
       !Classic finite centered difference formula for derivative (2nd order accurate)
        Rf = (pco2del(2) - pco2del(1)) / (dicdel(2) - dicdel(1))       ! dpCO2/dDIC
       !Rf = (pco2del(2) - pco2del(1)) / (dx)                          ! dpCO2/dDIC (same as just above)
        Rf = Rf * tc / dpco2                                           ! R = (dpCO2/dDIC) * (DIC/pCO2)

        BetaD(i) = SGLE(Rf)

     ELSE

        ph(i)     = 1.e20_rx
        pco2(i)   = 1.e20_rx
        fco2(i)   = 1.e20_rx
        co2(i)    = 1.e20_rx
        hco3(i)   = 1.e20_rx
        co3(i)    = 1.e20_rx
        OmegaA(i) = 1.e20_rx
        OmegaC(i) = 1.e20_rx
        BetaD(i)  = 1.e20_rx
        rhoSW(i)  = 1.e20_rx
        p(i)      = 1.e20_rx
        tempis(i) = 1.e20_rx

     ENDIF

  END DO

  RETURN
END SUBROUTINE vars_sprac


!>    Same routine as vars() above
!!    except that it perturbs slightly one dissociation constant 
SUBROUTINE vars_pertK(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC,       &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,        &
                var_index, abs_delta,                                       &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon, verbose )

  !   Purpose:
  !     
  !     Same routine as vars() above except that it perturbs slightly one dissociation constant 
  !     before proceeding to computation of carbonate system variables (call to varsolver())
  !
  !     This routine is intended only to be called internally by derivnum.f90
  !     which computes numerical derivatives.
  

  !     INPUT variables:
  !     ================
  !     same as routine vars() above,     plus :
  !     var_index  =   numerical id of dissociation constant (and total boron) to perturb
  !                    an index (1...8) from the list : K0, K1, K2, Kb, Kw, Kspa, Kspc, Bt
  !     abs_delta  =   perturbation value
  !
  !     INPUT options:
  !     ==============
  !     same as routine vars() above

  !     OUTPUT variables:
  !     =================
  !     same as routine vars() above except :  BetaD, rhoSW, p, tempis

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE mconstants
  USE mp80
  USE mrho
  USE meos
  USE gsw_mod_toolbox, only: gsw_t_from_ct, gsw_ct_from_t, gsw_rho
  USE msw_temp
  USE mvarsolver

  IMPLICIT NONE

! Input variables
  !>     number of records
!f2py optional , depend(sal) :: n=len(sal)
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: temp
  !> salinity <b>[psu]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: phos
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: lat
  !> Numeric id of dissociation constant
  INTEGER, INTENT(in)  ::  var_index
  !> Perturbation value
  REAL(kind=rx), INTENT(in)  :: abs_delta
  
  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS
  !> choose \b 'Sprc' for practical sal. (EOS-80, default) or \b 'Sabs' for absolute salinity (TEOS-10)
!f2py character*4 optional, intent(in) :: optS='Sprc'
  CHARACTER(4), OPTIONAL, INTENT(in) :: optS
  !> longitude <b>[degrees east]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lon = -25.
  REAL(kind=rx), OPTIONAL, INTENT(in), DIMENSION(N) :: lon
!f2py logical optional, intent(in) :: verbose
  LOGICAL, OPTIONAL, INTENT(in) :: verbose

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: OmegaC

! Local variables
  !> in-situ density of seawater in <b>[kg/m3]</b>
  REAL(kind=rx)  :: rhoSW
  !> pressure <b>[decibars]</b>
  REAL(kind=rx) :: p

  REAL(kind=rx) :: ssal, salk, sdic, ssil, sphos

  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=r8) :: tempot, tempis68, tempot68, tempis90, tempcsv
  REAL(kind=r8) :: drho

  ! local 1-long array version of scalar variables
  REAL(kind=r8), DIMENSION(1) :: aK0, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(kind=r8), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(kind=r8), DIMENSION(1) :: aSt, aFt, aBt

  REAL(kind=rx), DIMENSION(1) :: sabs1, spra1, p1, lon1, lat1
  REAL(kind=rx), DIMENSION(1) :: tc1, ta1, sit1, nt1
  REAL(kind=rx), DIMENSION(1) :: sal1 
  
  REAL(kind=r8) :: Patmd
  REAL(kind=r8) :: Ptot
  REAL(kind=r8) :: Phydro_atm

  INTEGER :: i, icount

  REAL(kind=r8) :: prb
  REAL(kind=r8) :: s
  REAL(kind=r8) :: tc, ta
  REAL(kind=r8) :: sit, pt
  
  REAL(kind=r8) :: dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS
  LOGICAL      :: verbosity

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
!   print *,"optB present:"
!   print *,"optB = ", optB 
    opB = optB
  ELSE
!   Default is Lee et al (2010) for total boron
!   print *,"optB NOT present:"
    opB = 'l10'
!   print *,"opB = ", opB 
  ENDIF
  IF (PRESENT(optKf)) THEN
!   print *,"optKf = ", optKf
    opKf = optKf
  ELSE
!   print *,"optKf NOT present:"
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
!   print *,"opKf = ", opKf
  ENDIF
  IF (PRESENT(optK1K2)) THEN
!   print *,"optK1K2 = ", optK1K2
    opK1K2 = optK1K2
  ELSE
!   print *,"optK1K2 NOT present:"
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
!   print *,"opK1K2 = ", opK1K2
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF
  IF (PRESENT(optS)) THEN
    opS = optS
  ELSE
    opS = 'Sprc'
  ENDIF
  IF (PRESENT(verbose)) THEN
    verbosity = verbose
  ELSE
    verbosity = .true.
  ENDIF

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!    - depth [m] => convert to pressure [db]
!    - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     !print *,"optP =", optP, "end"
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = depth(i)
     ELSE
        !print *,"optP =", optP, "end"
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = DBLE(temp(i))
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p, 0.d0 )
!jh        tempis68 = sw_temp(sal(i), SGLE(tempot68), p, SGLE(0.d0) )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis90 = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis90 = DBLE(temp(i))
        tempis68  = (tempis90 - 0.0002_r8) / 0.99975_r8
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
!       dtempot   = 0.99975*dtempot68 + 0.0002
     ELSEIF (trim(optT) == 'Tcsv' .OR. trim(optT) == 'tcsv') THEN
!       Convert given conservative temperature to in-situ temperature
        ! First convert salinity to absolute sal., if necessary
        IF (trim(opS) == 'Sprc')  THEN
            ! conversion will use default geographic location
            spra1(1) = sal(i)
            p1(1) = p
            CALL sp2sa_geo (spra1, 1, sabs1, p1)
        ELSE
            sabs1(1) = sal(i)
        END IF
        ! Then convert temperature
        tempis90 = gsw_t_from_ct (DBLE(sabs1(1)), DBLE(temp(i)), DBLE(p))
        tempis68  = (tempis90 - 0.0002_r8) / 0.99975_r8
     ELSE
        PRINT *,"optT must be either 'Tpot, 'Tinsitu' or 'Tcsv'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (verbosity .EQV. .true.) THEN
            IF (       sal(i) < 0.   &
                 .OR.  alk(i) < 0.   &
                 .OR.  dic(i) < 0.   &
                 .OR.  sil(i) < 0.   &
                 .OR. phos(i) < 0.   &
                 .OR.  sal(i) > 1e+3 &
                 .OR.  alk(i) > 1e+3 &
                 .OR.  dic(i) > 1e+3 &
                 .OR.  sil(i) > 1e+3 &
                 .OR. phos(i) > 1e+3) THEN
               PRINT *, 'i, icount, tempot, sal,    alk,    dic,    sil,    phos =', &
                         i, icount, tempot, sal(i), alk(i), dic(i), sil(i), phos(i)
            ENDIF
        ENDIF
!       Zero out any negative salinity, phosphate, silica, dic, and alk
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF
        IF (phos(i) < 0.0) THEN
           sphos = 0.0
        ELSE
           sphos = phos(i)
        ENDIF
        IF (sil(i) < 0.0) THEN
           ssil = 0.0
        ELSE
           ssil = sil(i)
        ENDIF
        IF (dic(i) < 0.0) THEN
          sdic = 0.0
        ELSE
          sdic = dic(i)
        ENDIF
        IF (alk(i) < 0.0) THEN
          salk = 0.0
        ELSE
          salk = alk(i)
        ENDIF

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))
!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p / SGLE(10.0d0))
        Phydro_atm = prb / 1.01325d0  ! convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
!       Total pressure [atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           Ptot = Patmd               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           Ptot = Patmd + Phydro_atm   ! total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF

!       Compute in-situ density [kg/m^3]
        ! If Absolute salinity is given
        IF (trim(opS) == 'Sabs')  THEN
           ! If in-situ or potential temperature is given
           IF (trim(optT) /= 'Scsv') THEN
              ! First compute conservative temperature
              tempcsv = gsw_ct_from_t (DBLE(ssal), tempis90, DBLE(p))
           ELSE
              tempcsv = DBLE(temp(i))
           ENDIF
           rhoSW = gsw_rho(DBLE(ssal), tempcsv, prb)
        ELSE
           rhoSW = rho(ssal, tempis68, prb)
!jh           rhoSW = rho(ssal, SGLE(tempis68), SGLE(prb))
        ENDIF

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (trim(optCON) == 'mol/kg') THEN
!          No conversion:  drho = 1.
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = 1.0
        ELSEIF (trim(optCON) == 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = DBLE(rhoSW)
        ELSE
           PRINT *,"optCON must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

!       Initialise ta, tc, pt and sit
        ta  = DBLE(salk)  / drho
        tc  = DBLE(sdic)  / drho
        pt  = DBLE(sphos) / drho
        sit = DBLE(ssil)  / drho

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Convert from Absolute to Practical salinity if needed
        IF (trim(opS) == 'Sabs')  THEN
           sabs1(1) = ssal   ! absolute sal.
           ! If longitude is passed in
           IF (PRESENT(lon)) THEN
               p1(1) = p
               IF (lon(i) .NE. 1e20_rx) THEN
                  ! longitude and latitude are defined
                  lon1(1) = lon(i)
                  lat1(1) = lat(i)
                  CALL sa2sp_geo (sabs1, 1, spra1, p1, lon1, lat1)
               ELSE
                  ! will use default geographic location
                  CALL sa2sp_geo (sabs1, 1, spra1, p1)
               ENDIF
           ELSE
               tc1(1) = tc
               ta1(1) = ta
               sit1(1) = sit
               ! Nitrate total : from Phosphate using Redfield ratio (16)
               nt1(1) = 16.0 * pt
               CALL sa2sp_chem(sabs1, ta1, tc1, nt1, sit1, 1, spra1)
           ENDIF
           s = DBLE(spra1(1))
        ENDIF
          
!       Get all equilibrium constants and total concentrations of SO4, F, B
        sal1(1) = SGLE(s)
        CALL constants (aK0, aK1, aK2, aKb, aKw, aKs, aKf,            &
                    aKspc, aKspa, aK1p, aK2p, aK3p, aKsi,             &
                    aSt, aFt, aBt,                                    &
                    temp(i), sal1, Patm(i),                           &
                    depth(i), lat(i), 1,                              &
                    optT, optP, opB, opK1K2, opKf, opGAS              )

!       Apply perturbation on K
        SELECT CASE (var_index)
            CASE (1)
                aK0(1) = aK0(1) + abs_delta
            CASE (2)
                aK1(1) = aK1(1) + abs_delta
            CASE (3)
                aK2(1) = aK2(1) + abs_delta
            CASE (4)
                aKb(1) = aKb(1) + abs_delta
            CASE (5)
                aKw(1) = aKw(1) + abs_delta
            CASE (6)
                aKspa(1) = aKspa(1) + abs_delta
            CASE (7)
                aKspc(1) = aKspc(1) + abs_delta
            CASE (8)
                aBt(1) = aBt(1) + abs_delta
        END SELECT

!       Solve for pH and all other variables
!       ------------------------------------
        CALL varsolver(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC,   &
                    tempis90, s, ta, tc, pt, sit,                                &
                    aBt(1), aSt(1), aFt(1),                                      &
                    aK0(1), aK1(1), aK2(1), aKb(1), aKw(1), aKs(1), aKf(1),      &
                    aKspc(1), aKspa(1), aK1p(1), aK2p(1), aK3p(1), aKsi(1),      & 
                    Patmd, prb, drho, opGAS                   )

!       Convert all output variables from double to single precision, if necessary
        pH(i)     = SGLE(dph)
        co2(i)    = SGLE(dco2)
        hco3(i)   = SGLE(dhco3)
        co3(i)    = SGLE(dco3)
        fCO2(i)   = SGLE(dfCO2)
        pCO2(i)   = SGLE(dpCO2)
        OmegaA(i) = SGLE(dOmegaA)
        OmegaC(i) = SGLE(dOmegaC)

     ELSE

        ph(i)     = 1.e20_rx
        pco2(i)   = 1.e20_rx
        fco2(i)   = 1.e20_rx
        co2(i)    = 1.e20_rx
        hco3(i)   = 1.e20_rx
        co3(i)    = 1.e20_rx
        OmegaA(i) = 1.e20_rx
        OmegaC(i) = 1.e20_rx

     ENDIF

  END DO

  RETURN
END SUBROUTINE vars_pertK
END MODULE mvars
