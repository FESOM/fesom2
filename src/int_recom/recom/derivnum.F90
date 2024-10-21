!> \file derivnum.f90
!! \BRIEF 
!> Module with derivnum subroutine - compute numerical derivatives of carbonate system vars 
!> with respect to DIC, Alk, total phosphorus, total silicon, T, S
MODULE mderivnum
CONTAINS


!>    Computes numerical derivatives of standard carbonate system variables  
!>    (H+, pCO2, fCO2, CO2*, HCO3- and CO3--, OmegaA, OmegaC) with respect to one given input variable
!!    Input variables are :
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations
SUBROUTINE derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,                   &
                     dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                     temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, derivar,   &
                     optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon     )
  !   Purpose:
  !     Computes numerical derivatives of standard carbonate system variables 
  !     (H+, pCO2, fCO2, CO2*, HCO3- and CO3--, OmegaA, OmegaC) with respect to one given input variable
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !             = conservative temperature [degrees C] (with optT='Tcsv')
  !     sal     = salinity in [psu]
  !     alk     = total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =               [eq/kg]  with optCON = 'mol/kg'
  !     dic     = dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = silica    [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     phos    = phosphate [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     derivar = 3-lowercase-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', 'sal'
  !               and dissociation constants 'k0', 'k1', 'k2', 'kw', 'kb', 'ka', 'kc'
  !               and total dissolved inorganic boron 'bt'
  !               The last two of the k's (constants) are the solubility products for Aragonite and Calcite
  !
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input & output concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC, ALK, sil, & phos given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC, ALK, sil, & phos given in mol/m^3 (std MODEL units)
  !     -----------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tcsv' means input is Conservative Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
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
  !     ---------
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

  !     OUTPUT variables:
  !     =================
  !     dh_dx      = derivative of ion [H+] concentration on total scale in [mol/kg]
  !     dpco2_dx   = derivative of CO2 partial pressure (uatm)
  !     dfco2_dx   = derivative of CO2 fugacity (uatm)
  !     dco2_dx    = derivative of aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dhco3_dx   = derivative of bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dco3_dx    = derivative of carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dOmegaA_dx = derivative of Omega for aragonite, i.e., the aragonite saturation state
  !     dOmegaC_dx = derivative of Omega for calcite, i.e., the   calcite saturation state


  USE msingledouble
  USE mvars
  
  IMPLICIT NONE

! Input variables
  !>     number of records
!f2py intent(hide) n
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
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
  ! 3-character identifier of input variable with respect to which derivative is requested
  ! ('alk', 'dic', 'pho', 'sil', 'tem', 'sal'
  CHARACTER(3), INTENT(in) ::  derivar

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
!f2py character*3 optional, intent(in) :: optB='u74'
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
!f2py character*7 optional, intent(in) :: optS='Sprc'
  CHARACTER(4), OPTIONAL, INTENT(in) :: optS
  !> longitude <b>[degrees east]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lon = -25.
  REAL(kind=rx), OPTIONAL, INTENT(in), DIMENSION(N) :: lon

! Output variables:
  !> derivative of H on the <b>total scale</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dh_dx
  !> derivative of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dpco2_dx
  !> derivative of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dfco2_dx
  !> derivative of aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dco2_dx
  !> derivative of (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dhco3_dx
  !> derivative of (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dco3_dx
  !> derivative of Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dOmegaA_dx
  !> derivative of Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dOmegaC_dx

! Local variables
  !> pH on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(1,2) :: ph
  !> ion (H+) concentration on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(1,2) :: h
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(1,2) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(1,2) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), DIMENSION(1,2) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), DIMENSION(1,2) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), DIMENSION(1) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=rx), DIMENSION(1) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx), DIMENSION(1) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=rx), DIMENSION(1) :: tempis

  ! Arrays that are copies of input parameters
  REAL(kind=rx), DIMENSION(1) :: atemp
  REAL(kind=rx), DIMENSION(1) :: asal
  REAL(kind=rx), DIMENSION(1) :: aalk
  REAL(kind=rx), DIMENSION(1) :: adic
  REAL(kind=rx), DIMENSION(1) :: asil
  REAL(kind=rx), DIMENSION(1) :: aphos
  REAL(kind=rx), DIMENSION(1) :: aPatm
  REAL(kind=rx), DIMENSION(1) :: adepth
  REAL(kind=rx), DIMENSION(1) :: alat
 
  ! value of small delta to apply to input variable when computing numerical derivative
  ! it is actually the ratio of delta relative to input variable value
  REAL(kind=rx) :: rel_delta_x
  ! Value of input variable with respect to which we derive
  REAL(kind=rx) :: input_value
  REAL(kind=rx), DIMENSION(1) :: ainput1, ainput2
  REAL(kind=rx) :: abs_delta, dX
  
! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(3) :: opK1K2
  CHARACTER(2) :: opKf
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS

  INTEGER :: i
  LOGICAL :: deriv_K    !  if one derives with respect to a dissoc constant
  INTEGER :: var_index  ! index of variable with respect to which one derives

! Approximate values for K0, K1, K2, Kb, Kspa and Kspc
! They will be used to compute an absolute perturbation value on these constants
  REAL,PARAMETER :: K_values(7) = (/0.034, 1.2e-06, 8.3e-10, 2.1e-09, 6.1e-14, 6.7e-07, 4.3e-07/)

! Multiplicative conversion factor to convert from umol/kg to optCON units (mol/kg or mol/m3)
  REAL(kind=rx) :: xfac
  
! Reference value (typically global surface means, for calculating consistent absolute derivatives)
  REAL(kind=rx) :: ref_value
  
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

  SELECT CASE (derivar)
      CASE ('alk')
          var_index = 1
          deriv_K   = .FALSE.
      CASE ('dic')
          var_index = 2
          deriv_K   = .FALSE.
      CASE ('pho')
          var_index = 3
          deriv_K   = .FALSE.
      CASE ('sil')
          var_index = 4
          deriv_K   = .FALSE.
      CASE ('tem')
          var_index = 5
          deriv_K   = .FALSE.
      CASE ('sal')
          var_index = 6
          deriv_K   = .FALSE.
      CASE ('k0')
          var_index = 1
          deriv_K   = .TRUE.
      CASE ('k1')
          var_index = 2
          deriv_K   = .TRUE.
      CASE ('k2')
          var_index = 3
          deriv_K   = .TRUE.
      CASE ('kb')
          var_index = 4
          deriv_K   = .TRUE.
      CASE ('kw')
          var_index = 5
          deriv_K   = .TRUE.
      CASE ('ka')
          var_index = 6
          deriv_K   = .TRUE.
      CASE ('kc')
          var_index = 7
          deriv_K   = .TRUE.
      CASE ('bt')
          ! For Bt (propotional to S): call vars_pertK but treat slightly differently than other constants   
          var_index = 8
          deriv_K   = .FALSE.
      CASE DEFAULT
          PRINT *,"derivar must be 3-char variable: 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'"
          PRINT *," or a 2-char dissoc. constant name : 'k0', 'k1', 'k2', 'kb', 'kw', 'ka', 'kc', 'bt'"
          STOP
  END SELECT
  
  DO i = 1, N
  
    ! Copy input values
    atemp(1)  = temp(i)
    asal(1)   = sal(i)
    aalk(1)   = alk(i)
    adic(1)   = dic(i)
    asil(1)   = sil(i)
    aphos(1)  = phos(i)
    aPatm(1)  = Patm(i)
    adepth(1) = depth(i)
    alat(1)   = lat(i)

    IF (deriv_K) THEN
        ! Choose value of absolute perturbation
        abs_delta = K_values(var_index) * 1.d-3   ! 0.1 percent of Kx value
        dx = 2 * abs_delta
    
        ! Call routine vars_pertK() with request to increase one K value
        call vars_pertK(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),      &
                        OmegaA(:,1), OmegaC(:,1),                                          &
                        atemp, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                        var_index, -abs_delta,                                             &
                        optCON, optT, optP, opB, opK1K2, opKf, opGAS,  optS, lon           ) 
        ! Call routine vars_pertK() with request to decrease the same K value
        call vars_pertK(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),      &
                        OmegaA(:,2), OmegaC(:,2),                                          &
                        atemp, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                        var_index, abs_delta,                                              &
                        optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon            )

    ELSE    ! we do NOT derive w/ respect to one of the  Ks

        ! Simple conversion factors
       IF (optCON .EQ. 'mol/kg') THEN
          xfac = 1.d-6    !Factor to convert umol/kg to mol/kg
       ELSE
          xfac = 1.025d-3 !Factor to convert umol/kg to mol/m3
       ENDIF
       
        ! Select input value to vary slightly
        ! values for individual rel_delta_x determined by minimizing difference
        ! between numerical derivatives and mocsy-dnad automatic derivatives
        ! (mocsy_dnad), which gives results as good as analytic solutions
        SELECT CASE (var_index)
            CASE (1)
                input_value = alk(i)
                ref_value = 2300.0d0 * xfac ! global surface mean in optCON units
                rel_delta_x = 1.0d-6
            CASE (2)
                input_value = dic(i)
                ref_value = 2000.0d0 * xfac     ! global surface mean in optCON units
                rel_delta_x = 1.0d-6
            CASE (3)
                input_value = phos(i)
                ref_value = 0.5d0 * xfac   ! global surface mean (Orr et al., 2017) in optCON units
                rel_delta_x = 1.0d-3
            CASE (4)
                input_value = sil(i)
                ref_value = 7.5d0 * xfac   ! global surface mean (Orr et al., 2017) in optCON units
                rel_delta_x = 1.0d-3
            CASE (5)
                input_value = temp(i)
                ref_value = 18.0d0 ! global surface mean (C)
                rel_delta_x = 1.0d-4
            CASE (6)
                input_value = sal(i)
                ref_value = 35.0d0 ! global surface mean (C)
                rel_delta_x = 1.0d-4
            CASE (8)
                ! For Bt (propotional to salinity)    
                input_value = 0.0002414d0 * (sal(i)  / 1.80655d0) / 10.811d0
                ref_value   = 0.0002414d0 * (35.0d0 / 1.80655d0) / 10.811d0
                rel_delta_x = 1.0d-3
        END SELECT

        ! Determine two slightly different values of selected input value
        abs_delta  = ref_value * rel_delta_x
        ainput1(1) = input_value - abs_delta
        ainput2(1) = input_value + abs_delta
        ! Compute total absolue delta
        dx = ainput2(1) - ainput1(1)
    
        ! Call routine vars() twice with two slightly different values of selected input
        SELECT CASE (var_index)
            CASE (1)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, ainput1, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. ) 
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, ainput2, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon               ) 
            CASE (2)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, ainput1, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. ) 
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, ainput2, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon               )
            CASE (3)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, adic, asil, ainput1, aPatm, adepth, alat, 1,       &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. )
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, adic, asil, ainput2, aPatm, adepth, alat, 1,       &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon             ) 
            CASE (4)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, adic, ainput1, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. )
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                atemp, asal, aalk, adic, ainput2, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon             ) 
            CASE (5)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                ainput1, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,       &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. )
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                ainput2, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,       &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon             ) 
            CASE (6)
                call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                                OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                                atemp, ainput1, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon, verbose=.false. )
                call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                                OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                                atemp, ainput2, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon             ) 
            CASE (8)
                ! For Bt (must treat with vars_pertK, but not quite like other 'constants' to account for proportionality to S)
                ! Call routine vars_pertK() with request to increase one Bt value
                call vars_pertK(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),      &
                                OmegaA(:,1), OmegaC(:,1),                                          &
                                atemp, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                var_index, -abs_delta,                                             &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon            ) 
                ! Call routine vars_pertK() with request to decrease the same Bt value
                call vars_pertK(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),      &
                                OmegaA(:,2), OmegaC(:,2),                                          &
                                atemp, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                                var_index, abs_delta,                                              &
                                optCON, optT, optP, opB, opK1K2, opKf, opGAS, optS, lon            ) 
        END SELECT
    ENDIF
    
    ! H+ concentration from ph
    h(1,1) = 10**(-ph(1,1))
    h(1,2) = 10**(-ph(1,2))
    
    ! Compute derivatives (centered differences)
    dh_dx(i)      = (h(1,2)    - h(1,1))    / dx
    dpco2_dx(i)   = (pco2(1,2) - pco2(1,1)) / dx
    dfco2_dx(i)   = (fco2(1,2) - fco2(1,1)) / dx
    dco2_dx(i)    = (co2(1,2)  - co2(1,1))  / dx
    dhco3_dx(i)   = (hco3(1,2) - hco3(1,1)) / dx
    dco3_dx(i)    = (co3(1,2)  - co3(1,1))  / dx
    dOmegaA_dx(i) = (OmegaA(1,2) - OmegaA(1,1)) / dx
    dOmegaC_dx(i) = (OmegaC(1,2) - OmegaC(1,1)) / dx
  END DO

  RETURN
END SUBROUTINE derivnum
END MODULE mderivnum
