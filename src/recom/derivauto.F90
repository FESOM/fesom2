!> \file derivauto.f90
!! \BRIEF 
!> Module with derivauto subroutine - compute partial derivatives of carbonate system vars from DIC,Alk,T,S,P,nuts
MODULE mderivauto
CONTAINS
!>    Computes partial derivatives of standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
!!    WITH RESPECT TO
!!    Alk, DIC, phosphate, silicate, temperature and salinity.
SUBROUTINE derivauto(ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv,   &
                OmegaA_deriv, OmegaC_deriv,                                                &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                       &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS                           )
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
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
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
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
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

  !     OUTPUT variables:
  !     =================
  !
  !     For each chemical OUTPUT variable, there is an OPTIONAL one that contains 
  !     partial derivatives with respect to 6 input variables : 
  !         alk, dic, phosphate, silicate, temperature and salinity.
  !
  !     ph_deriv   = derivatives of pH on total scale
  !     pco2_deriv = derivatives of CO2 partial pressure (uatm)
  !     fco2_deriv = derivatives of CO2 fugacity (uatm)
  !     co2_deriv  = derivatives of aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     hco3_deriv = derivatives of bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     co3_deriv  = derivatives of carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     OmegaA_deriv = derivatives of Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC_deriv = derivatives of Omega for calcite, i.e., the   calcite saturation state
  !

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE mconstants
  USE mp80
  USE mrho
  USE msw_temp
  USE mvarsolver
  USE Dual_Num_Auto_Diff
  
  IMPLICIT NONE

! Input variables
  !>     number of records
!f2py intent(hide) n
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

! Output variables:
  !> derivatives of pH on the <b>total scale</b>
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: ph_deriv
  !> derivatives of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: pco2_deriv
  !> derivatives of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: fco2_deriv
  !> derivatives of aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: co2_deriv
  !> derivatives of bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: hco3_deriv
  !> derivatives of carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: co3_deriv
  !> derivatives of Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: OmegaA_deriv
  !> derivatives of Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(out), DIMENSION(6,N) :: OmegaC_deriv

! Local variables
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx) :: press

  REAL(kind=rx) :: ssal, salk, sdic, ssil, sphos

  TYPE(DUAL_NUM) :: tempot, tempis68, tempot68, tempis90
! REAL(kind=r8) :: dtempot, dtempot68
  TYPE(DUAL_NUM) :: drho_sw, drho

  TYPE(DUAL_NUM), DIMENSION(1) :: K0, K1, K2, Kb, Kw, Ks, Kf, Kspc
  TYPE(DUAL_NUM), DIMENSION(1) :: Kspa, K1p, K2p, K3p, Ksi
  TYPE(DUAL_NUM), DIMENSION(1) :: St, Ft, Bt

  TYPE(DUAL_NUM), DIMENSION(1) :: Patmd
  TYPE(DUAL_NUM) :: pd
  REAL(kind=r8) :: Ptot
  REAL(kind=r8) :: Phydro_atm

  INTEGER :: i, icount

  TYPE(DUAL_NUM), DIMENSION(1) :: t
  TYPE(DUAL_NUM), DIMENSION(1) :: s
  TYPE(DUAL_NUM) :: prb
  TYPE(DUAL_NUM) :: tc, ta
  TYPE(DUAL_NUM) :: sit, pt

  ! Carbonate system output variables
  TYPE(DUAL_NUM) :: dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC

  TYPE (DUAL_NUM),PARAMETER:: zero = DUAL_NUM(0d0 ,(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

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
        press = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        press = depth(i)
     ELSE
        !print *,"optP =", optP, "end"
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = DUAL_NUM(DBLE(temp(i)),(/0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0/))
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
!       Salinity (in double precision and with derivatives)
        s(1) = DUAL_NUM(DBLE(sal(i)),(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0/))
!       Pressure (in double precision and with derivatives)
        pd = DUAL_NUM(DBLE(press),(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
        tempis68 = sw_temp_DNAD(s(1), tempot68, pd, zero )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis90 = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis90 = DUAL_NUM(DBLE(temp(i)),(/0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0/))
        tempis68  = (tempis90 - 0.0002) / 0.99975
!       dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
!       dtempot   = 0.99975*dtempot68 + 0.0002
     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
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
        Patmd(1) = DUAL_NUM(DBLE(Patm(i)),(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
!       Hydrostatic pressure (prb is in bars)
        prb = DUAL_NUM(DBLE(press) / 10.0d0,(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
        Phydro_atm = prb%x_ad_ / 1.01325d0  ! convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
!       Total pressure [atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           Ptot = Patmd(1)%x_ad_               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           Ptot = Patmd(1)%x_ad_               ! total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           Ptot = Patmd(1)%x_ad_ + Phydro_atm   ! total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF

!       Salinity (in double precision and with derivative w/ respect to itelf)
        s(1) = DUAL_NUM(DBLE(ssal),(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0/))

!       Temperature (Celsius) and its derivatives relative to input variables
!       (variable 't' used here just for calculation of constants)
        t(1) = DUAL_NUM(DBLE(temp(i)),(/0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0/))

!       Get all equilibrium constants and total concentrations of SO4, F, B
!       and their derivatives with respect to temperature and salinity
        CALL constants_DNAD( K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,     &
                              K1p, K2p, K3p, Ksi, St, Ft, Bt,            &
                              t, s, Patmd,                               &
                              depth(i), lat(i), 1,                       &
                              optT, optP, opB, opK1K2, opKf, opGAS       )

!       Compute in-situ density [kg/m^3]
!       and derivatives with respect to temperature and salinity
        drho_sw = rho_DNAD(s(1), tempis68, prb)

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (trim(optCON) == 'mol/kg') THEN
!          No conversion:  drho = 1.
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = DUAL_NUM(1.0,(/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
        ELSEIF (trim(optCON) == 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = drho_sw
        ELSE
           PRINT *,"optCON must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

!       Initialise ta, tc, pt and sit and their derivatives 
        ta  = DUAL_NUM(DBLE(salk),  (/1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/)) / drho
        tc  = DUAL_NUM(DBLE(sdic),  (/0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0/)) / drho
        pt  = DUAL_NUM(DBLE(sphos), (/0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0/)) / drho
        sit = DUAL_NUM(DBLE(ssil),  (/0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0/)) / drho

!       Solve for pH and all other variables
!       ------------------------------------

!       Compute chemical variables and their derivatives
        CALL varsolver_DNAD(dph, dpco2, dfco2, dco2, dhco3, dco3, dOmegaA, dOmegaC, &
                    tempis90, s(1), ta, tc, pt, sit,                                &
                    Bt(1), St(1), Ft(1),                                            &
                    K0(1), K1(1), K2(1), Kb(1), Kw(1), Ks(1), Kf(1),                &
                    Kspc(1), Kspa(1), K1p(1), K2p(1), K3p(1), Ksi(1),               & 
                    Patmd(1), prb, drho, opGAS                                      )

!       Save derivatives
        ph_deriv(:,i)     = SGLE(dph%xp_ad_)
        pCO2_deriv(:,i)   = SGLE(dpCO2%xp_ad_)
        fCO2_deriv(:,i)   = SGLE(dfCO2%xp_ad_)
        co2_deriv(:,i)    = SGLE(dco2%xp_ad_)
        hco3_deriv(:,i)   = SGLE(dhco3%xp_ad_)
        co3_deriv(:,i)    = SGLE(dco3%xp_ad_)
        OmegaA_deriv(:,i) = SGLE(dOmegaA%xp_ad_)
        OmegaC_deriv(:,i) = SGLE(dOmegaC%xp_ad_)

     ELSE

        ph_deriv(:,i)     = 1.e20_rx
        pCO2_deriv(:,i)   = 1.e20_rx
        fCO2_deriv(:,i)   = 1.e20_rx
        co2_deriv(:,i)    = 1.e20_rx
        hco3_deriv(:,i)   = 1.e20_rx
        co3_deriv(:,i)    = 1.e20_rx
        OmegaA_deriv(:,i) = 1.e20_rx
        OmegaC_deriv(:,i) = 1.e20_rx

     ENDIF

  END DO

  RETURN
END SUBROUTINE derivauto
END MODULE mderivauto
