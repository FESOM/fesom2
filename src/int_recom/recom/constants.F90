!> \file constants.f90
!! \BRIEF 
!> Module with contants subroutine - computes carbonate system constants
!! from S,T,P 
MODULE mconstants
CONTAINS
!> Compute thermodynamic constants
!! FROM temperature, salinity, and pressure (1D arrays)
SUBROUTINE constants(K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
                     K1p, K2p, K3p, Ksi,                      &
                     St, Ft, Bt,                              &
                     temp, sal, Patm,                         &
                     depth, lat, N,                           &
                     optT, optP, optB, optK1K2, optKf, optGAS, optS, lon )

  !   Purpose:
  !     Compute thermodynamic constants
  !     FROM: temperature, salinity, and pressure (1D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     ---------
  !     optT: choose in-situ, potential or conservative temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp is computed)
  !       -> 'Tcsv' means input is Conservative Temperature (in situ Temp is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)
  !     ---------
  !     optB:
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
  !       -> 'w14' means use Waters (2014) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     -----------
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
  !        Optional, it may be used along with latitude when optS is "Sabs" as conversion parameters 
  !           from Absolute to Practical Salinity.
  !
  !        When seawater is not of standard composition, Practical Salinity alone is not sufficient 
  !        to compute Absolute Salinity and vice-versa. One needs to know the chemical composition, 
  !        mainly silicate and nitrate concentration. When parameters 'lon' and 'lat' are given, 
  !        absolute salinity conversion is based on WOA silicate concentration at given location. 
  !
  !        When 'lon' and 'lat' are unknown, an arbitrary geographic point is chosen:
  !        which is mid equatorial Atlantic. Note that this implies an error on computed practical salinity up to 0.02 psu.
  !        In that case, do not pass parameter 'lon' or set its elements to 1.e20.
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     St, Ft, Bt

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE mp80
  USE meos
  USE gsw_mod_toolbox, only: gsw_t_from_ct
  USE msw_temp
  USE msw_ptmp
  IMPLICIT NONE

! Input variables
  !>     number of records
!f2py intent(hide) n
  INTEGER, INTENT(in) :: N
  !> in <b>situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: temp
  !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: lat
  !> salinity <b>[psu]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal

  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm

  !> for temp input, choose \b 'Tinsitu' for in situ Temp or 
  !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
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
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) or \b 'w14' (Waters et al., 2014)
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
  REAL(kind=rx), OPTIONAL, INTENT(in),    DIMENSION(N) :: lon

! Ouput variables
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Bt

! Local variables
  REAL(kind=rx) :: ssal
  REAL(kind=rx) :: p
  REAL(kind=rx) :: tempot, tempis68, tempot68
  REAL(kind=rx) :: tempis
  REAL(kind=r8) :: is, invtk, dlogtk, is2, s2, sqrtis
  REAL(kind=r8) :: Ks_0p, Kf_0p
  REAL(kind=r8) :: total2free, free2SWS, total2SWS, SWS2total
  REAL(kind=r8) :: total2free_0p, free2SWS_0p, total2SWS_0p
! REAL(kind=r8) :: free2SWS, free2SWS_0p

  REAL(kind=r8) :: dtempot, dtempot68
  REAL(kind=r8) :: R

  REAL(kind=r8) :: pK1o, ma1, mb1, mc1, pK1
  REAL(kind=r8) :: pK2o, ma2, mb2, mc2, pK2

  REAL(kind=r8), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
  REAL(kind=r8), DIMENSION(12) :: deltav, deltak, lnkpok0
  REAL(kind=r8) :: tmp, nK0we74

  INTEGER :: i, icount, ipc

  REAL(kind=r8) :: t, tk, tk0, prb
  REAL(kind=r8) :: s, sqrts, s15, scl

  REAL(kind=r8) :: Phydro_atm, Patmd, Ptot, Rgas_atm, vbarCO2

! local 1-long array version of scalar variables
  REAL(kind=rx), DIMENSION(1) :: lon1, lat1
  REAL(kind=rx), DIMENSION(1) :: p1, spra1, sabs1

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS

  ! CONSTANTS
  ! =========
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! with corrected coefficients for Kb, Kw, Ksi, etc.

  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi

  DATA a0 /-25.5_r8, -15.82_r8, -29.48_r8, -20.02_r8, &
          -18.03_r8,  -9.78_r8, -48.76_r8, -45.96_r8, &
          -14.51_r8, -23.12_r8, -26.57_r8, -29.48_r8/
  DATA a1 /0.1271_r8, -0.0219_r8, 0.1622_r8, 0.1119_r8, &
           0.0466_r8, -0.0090_r8, 0.5304_r8, 0.5304_r8, &
           0.1211_r8, 0.1758_r8, 0.2020_r8, 0.1622_r8/
  DATA a2 /     0.0_r8,       0.0_r8, -2.608e-3_r8, -1.409e-3_r8, &
           0.316e-3_r8, -0.942e-3_r8,  0.0_r8,       0.0_r8, &
          -0.321e-3_r8, -2.647e-3_r8, -3.042e-3_r8, -2.6080e-3_r8/
  DATA b0 /-3.08e-3_r8, 1.13e-3_r8,  -2.84e-3_r8,   -5.13e-3_r8, &
           -4.53e-3_r8, -3.91e-3_r8, -11.76e-3_r8, -11.76e-3_r8, &
           -2.67e-3_r8, -5.15e-3_r8,  -4.08e-3_r8,  -2.84e-3_r8/
  DATA b1 /0.0877e-3_r8, -0.1475e-3_r8, 0.0_r8,       0.0794e-3_r8, &
           0.09e-3_r8,    0.054e-3_r8,  0.3692e-3_r8, 0.3692e-3_r8, &
           0.0427e-3_r8,  0.09e-3_r8,   0.0714e-3_r8, 0.0_r8/
  DATA b2 /12*0.0_r8/

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
    opB = 'u74'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
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

  R = 83.14472_r8

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
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = p80(depth(i), lat(i))
     ELSEIF (trim(optP) == 'db' ) THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = depth(i)
     ELSE
        PRINT *,"optP must be 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002_rx) / 0.99975_rx
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p, 0.0d0) !! JH changed: SGLE(0.0D0) ) to 0.0d0)
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis = 0.99975_rx*tempis68 + 0.0002_rx
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis    = temp(i)
        tempis68  = (temp(i) - 0.0002_rx) / 0.99975_rx
        dtempot68 = sw_ptmp(DBLE(sal(i)), DBLE(tempis68), DBLE(p), 0.0d0)
        dtempot   = 0.99975_rx*dtempot68 + 0.0002_rx
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
        tempis = SGLE(gsw_t_from_ct (DBLE(sabs1(1)), DBLE(temp(i)), DBLE(p)))
        tempis68  = (tempis - 0.0002_rx) / 0.99975_rx
     ELSE
        PRINT *,"optT must be either 'Tpot, 'Tinsitu' or 'Tcsv'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    Compute constants:
     IF (temp(i) >= -5.0_rx .AND. temp(i) < 1.0e+2_rx) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (      sal(i) < 0.  .OR.  sal(i) > 1e+3) THEN
           PRINT *, 'i, icount, temp, sal =', i, icount, temp(i), sal(i)
        ENDIF
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF

!       Absolute temperature (Kelvin) and related values
        t = DBLE(tempis)
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Atmospheric pressure
        Patmd = DBLE(Patm(i))

!       Hydrostatic pressure (prb is in bars)
        prb = DBLE(p) / 10.0d0

!       Convert from Absolute to Practical salinity if needed
        IF (trim(opS) == 'Sabs')  THEN
           IF (PRESENT(lon)) THEN
               ! longitude is passed in
               sabs1(1) = ssal
               p1(1) = p
               lon1(1) = lon(i)
               lat1(1) = lat(i)
           ELSE
               lon1(1) = 1.e20_rx
               lat1(1) = 1.e20_rx
           ENDIF
           IF (lon1(1) .NE. 1.e20_rx .AND. lat1(1) .NE. 1.e20_rx) THEN
              ! longitude and latitude are defined
              CALL sa2sp_geo (sabs1, 1, spra1, p1, lon1, lat1)
           ELSE
              ! use default geographic location
              CALL sa2sp_geo (sabs1, 1, spra1, p1)
           ENDIF
           s = DBLE(spra1(1))
        ELSE
           s = DBLE(ssal)
        ENDIF

        ! Salinity and simply related values
        s2=s*s
        sqrts=SQRT(s)
        s15=s**1.5d0
        scl=s/1.80655d0

!       Ionic strength:
        is = 19.924d0*s/(1000.0d0 - 1.005d0*s)
        is2 = is*is
        sqrtis = SQRT(is)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
        St(i) = 0.14d0 * scl/96.062d0

!       Fluoride:  Riley (1965)
        Ft(i) = 0.000067d0 * scl/18.9984d0

!       Boron:
        IF (trim(opB) == 'l10') THEN
!          New formulation from Lee et al (2010)
           Bt(i) = 0.0002414d0 * scl/10.811d0
        ELSEIF (trim(opB) == 'u74') THEN
!          Classic formulation from Uppström (1974)
           Bt(i) = 0.000232d0  * scl/10.811d0
        ELSE
           PRINT *,"optB must be 'l10' or 'u74'"
           STOP
        ENDIF

!       K0 (K Henry)
!       CO2(g) <-> CO2(aq.)
!       K0  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           tk0 = tk                   !in situ temperature (K) for K0 calculation
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           tk0 = dtempot + 273.15d0   !potential temperature (K) for K0 calculation as needed for potential fCO2 & pCO2
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           tk0 = tk                     !in situ temperature (K) for K0 calculation
           Phydro_atm = prb / 1.01325d0 !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
           Ptot = Patmd + Phydro_atm    !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF
        tmp = 9345.17d0/tk0 - 60.2409d0 + 23.3585d0 * LOG(tk0/100.0d0)
        nK0we74 = tmp + s*(0.023517d0 - 0.00023656d0*tk0 + 0.0047036e-4_r8*tk0*tk0)
        K0(i) = EXP(nK0we74)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
        IF (trim(opK1K2) == 'l') THEN
!         Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
          K1(i) = 10.0d0**(-1.0d0*(3633.86d0*invtk - 61.2172d0 + 9.6777d0*dlogtk  &
                  - 0.011555d0*s + 0.0001152d0*s2))
          K2(i) = 10.0d0**(-1*(471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk      &
                  - 0.01781d0*s + 0.0001122d0*s2))
        ELSEIF (trim(opK1K2) == 'm10') THEN
!         Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.4038d0*sqrts + 0.03206d0*s - (5.242e-5_rx)*s2
          mb1 = -530.659d0*sqrts - 5.8210d0*s
          mc1 = -2.0664d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.3728d0*sqrts + 0.1218d0*s - (3.688e-4_r8)*s2
          mb2 = -788.289d0*sqrts - 19.189d0*s
          mc2 = -3.374d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSEIF (trim(opK1K2) == 'w14') THEN
!         Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.409160d0*sqrts + 0.031646d0*s - (5.1895e-5)*s2
          mb1 = -531.3642d0*sqrts - 5.713d0*s
          mc1 = -2.0669166d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.225890d0*sqrts + 0.12450870d0*s - (3.7243e-4_r8)*s2
          mb2 = -779.3444d0*sqrts - 19.91739d0*s
          mc2 = -3.3534679d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSE
           PRINT *, "optK1K2 must be either 'l', 'm10', or 'w14'"
           STOP
        ENDIF

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
        Kb(i) = EXP((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +  &
                1.728d0*s15 - 0.0996d0*s2)*invtk +              &
                (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +   &
                (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *      &
                dlogtk + 0.053105d0*sqrts*tk)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equation's 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K1p(i) = EXP(-4576.752d0*invtk + 115.540d0 - 18.453d0*dlogtk +  &
                 (-106.736d0*invtk + 0.69171d0) * sqrts +             &
                 (-0.65643d0*invtk - 0.01844d0) * s)

!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equation's 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K2p(i) = EXP(-8814.715d0*invtk + 172.1033d0 - 27.927d0*dlogtk +  &
                 (-160.340d0*invtk + 1.3566d0)*sqrts +                 &
                 (0.37335d0*invtk - 0.05778d0)*s)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equation's 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K3p(i) = EXP(-3070.75d0*invtk - 18.126d0 +            &
                 (17.27039d0*invtk + 2.81197d0) *             &
                 sqrts + (-44.99486d0*invtk - 0.09984d0) * s)

!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equation's 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Ksi(i) = EXP(-8904.2d0*invtk  + 117.400d0 - 19.334d0*dlogtk +  &
                 (-458.79d0*invtk + 3.5913d0) * sqrtis +             &
                 (188.74d0*invtk - 1.5998d0) * is +                  &
                 (-12.1652d0*invtk + 0.07871d0) * is2 +              &
                 LOG(1.0 - 0.001005d0*s))

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equation's 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Kw(i) = EXP(-13847.26d0*invtk + 148.9802d0 - 23.6521d0*dlogtk +  &
               (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *          &
               sqrts - 0.01615d0 * s)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks_0p = EXP(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk          &
                + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk) * sqrtis  &
                + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is      &
                - 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is2              &
                + LOG(1.0d0 - 0.001005d0*s))

!       Kf = [H][F]/[HF]
!       (total scale)
        IF (trim(opKf) == 'dg') THEN
!          Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
           Kf_0p = EXP(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +  &
                   LOG(1.0d0 - 0.001005d0*s) +                     &
                   LOG(1.0d0 + St(i)/Ks_0p))
        ELSEIF (trim(opKf) == 'pf') THEN
!          Perez and Fraga (1987) - Already on Total scale (no need for last line above)
!          Formulation as given in Dickson et al. (2007)
           Kf_0p = EXP(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)
        ELSE
           PRINT *, "optKf must be either 'dg' or 'pf'"
           STOP
        ENDIF

!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
        Kspc(i) = 10d0**(-171.9065d0 - 0.077993d0*tk + 2839.319d0/tk    &
                 + 71.595d0*LOG10(tk)                             &
                 + (-0.77712d0 + 0.0028426d0*tk + 178.34d0/tk)*sqrts  &
                 -0.07711d0*s + 0.0041249d0*s15 )


!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
        Kspa(i) = 10.d0**(-171.945d0 - 0.077993d0*tk + 2903.293d0/tk &
             +71.595d0*LOG10(tk) &
             +(-0.068393d0 + 0.0017276d0*tk + 88.135d0/tk)*sqrts &
             -0.10018d0*s + 0.0059415d0*s15 )

!       Pressure effect on K0 based on Weiss (1974, equation 5)
        Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
        vbarCO2 = 32.3_r8           ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
        K0(i) = K0(i) * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk0) )   ! Weiss (1974, equation 5)

!       Pressure effect on all other K's (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
        DO ipc = 1, 12
           deltav(ipc)  =  a0(ipc) + a1(ipc) *t + a2(ipc) *t*t
           deltak(ipc)   = (b0(ipc)  + b1(ipc) *t + b2(ipc) *t*t)
           lnkpok0(ipc)  = (-(deltav(ipc)) &
                +(0.5d0*deltak(ipc) * prb) &
                )                         * prb/(R*tk)
        END DO

!       Pressure correction on Ks (Free scale)
        Ks(i) = Ks_0p*EXP(lnkpok0(5))
!       Conversion factor total -> free scale
        total2free     = 1.d0/(1.d0 + St(i)/Ks(i))   ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
        total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
        Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
        Kf(i) = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
        Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total hydrogen (pH) scales
        free2SWS  = 1.d0 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
        total2SWS = total2free * free2SWS                          ! KSWS = Ktotal*total2SWS
        SWS2total = 1.d0 / total2SWS
!       Conversion at pressure zero
        free2SWS_0p  = 1.d0 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  ! using Kf on free scale
        total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
        IF (trim(optK1K2) == 'l') THEN
          K1(i)  = K1(i)*total2SWS_0p
          K2(i)  = K2(i)*total2SWS_0p
          !This conversion is unnecessary for the K1,K2 from Millero (2010),
          !since we use here the formulation already on the seawater scale
        ENDIF
        Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw

!       Other contants (keep on another scale):
!          - K0         (independent of pH scale, already pressure corrected)
!          - Ks         (already on Free scale;   already pressure corrected)
!          - Kf         (already on Total scale;  already pressure corrected)
!          - Kspc, Kspa (independent of pH scale; pressure-corrected below)

!       Perform actual pressure correction (on seawater scale)
        K1(i)   = K1(i)*EXP(lnkpok0(1))
        K2(i)   = K2(i)*EXP(lnkpok0(2))
        Kb(i)   = Kb(i)*EXP(lnkpok0(3))
        Kw(i)   = Kw(i)*EXP(lnkpok0(4))
        Kspc(i) = Kspc(i)*EXP(lnkpok0(7))
        Kspa(i) = Kspa(i)*EXP(lnkpok0(8))
        K1p(i)  = K1p(i)*EXP(lnkpok0(9))
        K2p(i)  = K2p(i)*EXP(lnkpok0(10))
        K3p(i)  = K3p(i)*EXP(lnkpok0(11))
        Ksi(i)  = Ksi(i)*EXP(lnkpok0(12))

!       Convert back to original total scale:
        K1(i)  = K1(i) *SWS2total
        K2(i)  = K2(i) *SWS2total
        K1p(i) = K1p(i)*SWS2total
        K2p(i) = K2p(i)*SWS2total
        K3p(i) = K3p(i)*SWS2total
        Kb(i)  = Kb(i) *SWS2total
        Ksi(i) = Ksi(i)*SWS2total
        Kw(i)  = Kw(i) *SWS2total

     ELSE

        K0(i)   = 1.e20_r8
        K1(i)   = 1.e20_r8
        K2(i)   = 1.e20_r8
        Kb(i)   = 1.e20_r8
        Kw(i)   = 1.e20_r8
        Ks(i)   = 1.e20_r8
        Kf(i)   = 1.e20_r8
        Kspc(i) = 1.e20_r8
        Kspa(i) = 1.e20_r8
        K1p(i)  = 1.e20_r8
        K2p(i)  = 1.e20_r8
        K3p(i)  = 1.e20_r8
        Ksi(i)  = 1.e20_r8
        Bt(i)   = 1.e20_r8
        Ft(i)   = 1.e20_r8
        St(i)   = 1.e20_r8

     ENDIF

  END DO

  RETURN
END SUBROUTINE constants

!> Compute thermodynamic constants
!! and their derivatives with respect to temperature and salinity
!! FROM temperature, salinity, and pressure (1D arrays)
SUBROUTINE constants_DNAD(K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
                     K1p, K2p, K3p, Ksi,                      &
                     St, Ft, Bt,                              &
                     temp, sal, Patm,                         &
                     depth, lat, N,                           &
                     optT, optP, optB, optK1K2, optKf, optGAS, optS, lon )

  !   Purpose:
  !     Compute thermodynamic constants
  !     FROM: temperature, salinity, and pressure (1D arrays)
  !
  !     It is similar to subroutine 'constants' above except that it also computes
  !     partial derivatives of all output thermodynamic constants
  !     with respect to temperature and salinity.
  !
  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     ---------
  !     optT: choose in-situ, potential or conservative temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp is computed)
  !       -> 'Tcsv' means input is Conservative Temperature (in situ Temp is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)
  !     ---------
  !     optB:
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
  !     -----------
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
  !        Optional, it may be used along with latitude when optS is "Sabs" as conversion parameters 
  !           from Absolute to Practical Salinity.
  !
  !        When seawater is not of standard composition, Practical Salinity alone is not sufficient 
  !        to compute Absolute Salinity and vice-versa. One needs to know the chemical composition, 
  !        mainly silicate and nitrate concentration. When parameters 'lon' and 'lat' are given, 
  !        absolute salinity conversion is based on WOA silicate concentration at given location. 
  !
  !        When 'lon' and 'lat' are unknown, an arbitrary geographic point is chosen:
  !        which is mid equatorial Atlantic. Note that this implies an error on computed practical salinity up to 0.02 psu.
  !        In that case, do not pass parameter 'lon' or set its elements to 1.e20.
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     St, Ft, Bt

  USE msingledouble
  USE mp80
  USE meos
  USE gsw_mod_toolbox, only: gsw_t_from_ct
  USE msw_temp
  USE msw_ptmp
  USE Dual_Num_Auto_Diff
  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> in <b>situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
  TYPE(DUAL_NUM), INTENT(in),    DIMENSION(N) :: temp
  !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in),    DIMENSION(N) :: lat
  !> salinity <b>[psu]</b>
  TYPE(DUAL_NUM), INTENT(in), DIMENSION(N) :: sal
!f2py optional , depend(sal) :: n=len(sal)

  !> atmospheric pressure <b>[atm]</b>
  TYPE(DUAL_NUM), INTENT(in), DIMENSION(N) :: Patm

  !> for temp input, choose \b 'Tinsitu' for in situ Temp or 
  !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
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
  REAL(kind=rx), OPTIONAL, INTENT(in),    DIMENSION(N) :: lon

! Ouput variables
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K2
  !> equilibrium constant for dissociation of boric acid 
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: St
  !> total fluoride  (Riley, 1965)
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  TYPE(DUAL_NUM), INTENT(out), DIMENSION(N) :: Bt

! Local variables
  TYPE(DUAL_NUM) ::ssal
  TYPE(DUAL_NUM) ::p
  TYPE(DUAL_NUM) ::tempot, tempis68, tempot68
  TYPE(DUAL_NUM) ::tempis
  TYPE(DUAL_NUM) ::is, invtk, dlogtk, is2, s2, sqrtis
  TYPE(DUAL_NUM) ::Ks_0p, Kf_0p
  TYPE(DUAL_NUM) ::total2free, free2SWS, total2SWS, SWS2total
  TYPE(DUAL_NUM) ::total2free_0p, free2SWS_0p, total2SWS_0p
! TYPE(DUAL_NUM) ::free2SWS, free2SWS_0p

  TYPE(DUAL_NUM) ::dtempot, dtempot68
  TYPE(DUAL_NUM) ::R

  TYPE(DUAL_NUM) ::pK1o, ma1, mb1, mc1, pK1
  TYPE(DUAL_NUM) ::pK2o, ma2, mb2, mc2, pK2

  REAL(kind=r8), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
  TYPE(DUAL_NUM), DIMENSION(12) :: deltav, deltak, lnkpok0
  TYPE(DUAL_NUM) ::tmp, nK0we74

  INTEGER :: i, icount, ipc

  TYPE(DUAL_NUM) ::t, tk, tk0, prb
  TYPE(DUAL_NUM) ::s, sqrts, s15, scl

  TYPE(DUAL_NUM) ::Phydro_atm, Patmd, Ptot, Rgas_atm, vbarCO2

  TYPE (DUAL_NUM),PARAMETER:: zero=DUAL_NUM(0d0,0.D0), ten=DUAL_NUM(10d0,0.D0)

! local 1-long array version of scalar variables
  REAL(kind=rx), DIMENSION(1) :: lon1, lat1
  REAL(kind=rx), DIMENSION(1) :: p1, spra1, sabs1

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS

  ! CONSTANTS
  ! =========
  ! Constants in formulation for Pressure effect on K's (Millero, 95)
  ! with corrected coefficients for Kb, Kw, Ksi, etc.

  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi

  DATA a0 /-25.5_r8, -15.82_r8, -29.48_r8, -20.02_r8, &
          -18.03_r8,  -9.78_r8, -48.76_r8, -45.96_r8, &
          -14.51_r8, -23.12_r8, -26.57_r8, -29.48_r8/
  DATA a1 /0.1271_r8, -0.0219_r8, 0.1622_r8, 0.1119_r8, &
           0.0466_r8, -0.0090_r8, 0.5304_r8, 0.5304_r8, &
           0.1211_r8, 0.1758_r8, 0.2020_r8, 0.1622_r8/
  DATA a2 /     0.0_r8,       0.0_r8, -2.608e-3_r8, -1.409e-3_r8, &
           0.316e-3_r8, -0.942e-3_r8,  0.0_r8,       0.0_r8, &
          -0.321e-3_r8, -2.647e-3_r8, -3.042e-3_r8, -2.6080e-3_r8/
  DATA b0 /-3.08e-3_r8, 1.13e-3_r8,  -2.84e-3_r8,   -5.13e-3_r8, &
           -4.53e-3_r8, -3.91e-3_r8, -11.76e-3_r8, -11.76e-3_r8, &
           -2.67e-3_r8, -5.15e-3_r8,  -4.08e-3_r8,  -2.84e-3_r8/
  DATA b1 /0.0877e-3_r8, -0.1475e-3_r8, 0.0_r8,       0.0794e-3_r8, &
           0.09e-3_r8,    0.054e-3_r8,  0.3692e-3_r8, 0.3692e-3_r8, &
           0.0427e-3_r8,  0.09e-3_r8,   0.0714e-3_r8, 0.0_r8/
  DATA b2 /12*0.0_r8/

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
    opB = optB
  ELSE
    opB = 'l10'
  ENDIF
  IF (PRESENT(optKf)) THEN
    opKf = optKf
  ELSE
    opKf = 'pf'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
    opK1K2 = optK1K2
  ELSE
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

  R = 83.14472_r8

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
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (trim(optP) == 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = DUAL_NUM(p80(depth(i), lat(i)), (/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
     ELSEIF (trim(optP) == 'db' ) THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = DUAL_NUM(depth(i), (/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/))
     ELSE
        PRINT *,"optP must be 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (trim(optT) == 'Tpot' .OR. trim(optT) == 'tpot') THEN
        tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot - 0.0002) / 0.99975
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp_DNAD(sal(i), tempot68, p, zero )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis = 0.99975*tempis68 + 0.0002
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (trim(optT) == 'Tinsitu' .OR. trim(optT) == 'tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis    = temp(i)
        tempis68  = (temp(i) - 0.0002) / 0.99975
        dtempot68 = sw_ptmp_DNAD(sal(i), tempis68, p, zero )
        dtempot   = 0.99975*dtempot68 + 0.0002
     ELSEIF (trim(optT) == 'Tcsv' .OR. trim(optT) == 'tcsv') THEN
!       Convert given conservative temperature to in-situ temperature
        ! First convert salinity to absolute sal., if necessary
        IF (trim(opS) == 'Sprc')  THEN
            ! conversion will use default geographic location
            spra1(1) = SGLE(sal(i)%x_ad_)
            p1(1) = SGLE(p%x_ad_)
            CALL sp2sa_geo (spra1, 1, sabs1, p1)
        ELSE
            sabs1(1) = SGLE(sal(i)%x_ad_)
        END IF
        ! Then convert temperature
        tempis%x_ad_ = gsw_t_from_ct (DBLE(sabs1(1)), temp(i)%x_ad_, p%x_ad_)
!       Sorry but no computation of derivatives because gsw_t_from_ct does not support it.
!          instead, assume that any derivative d(tempis)/dx is as d(temp)/dx
        tempis%xp_ad_(:) = temp(i)%xp_ad_(:)
        tempis68  = (tempis - 0.0002_r8) / 0.99975_r8
     ELSE
        PRINT *,"optT must be either 'Tpot, 'Tinsitu' or 'Tcsv'"
        PRINT *,"you specified optT =", trim(optT) 
        STOP
     ENDIF

!    Compute constants:
     IF (temp(i) >= -5. .AND. temp(i) < 1.0e+2) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (      sal(i) < 0.  .OR.  sal(i) > 1e+3) THEN
           PRINT *, 'i, icount, temp, sal =', i, icount, temp(i), sal(i)
        ENDIF
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF

!       Absolute temperature (Kelvin) and related values
        t = tempis
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Atmospheric pressure
        Patmd = Patm(i)

!       Hydrostatic pressure (prb is in bars)
        prb = p / 10.0d0

!       Convert from Absolute to Practical salinity if needed
        IF (trim(opS) == 'Sabs')  THEN
           IF (PRESENT(lon)) THEN
               ! longitude is passed in
               sabs1(1) = SGLE(ssal%x_ad_)
               p1(1) = SGLE(p%x_ad_)
               lon1(1) = lon(i)
               lat1(1) = lat(i)
           ELSE
               lon1(1) = 1.e20_rx
               lat1(1) = 1.e20_rx
           ENDIF
           IF (lon1(1) .NE. 1.e20_rx .AND. lat1(1) .NE. 1.e20_rx) THEN
              ! longitude and latitude are defined
              CALL sa2sp_geo (sabs1, 1, spra1, p1, lon1, lat1)
           ELSE
              ! use default geographic location
              CALL sa2sp_geo (sabs1, 1, spra1, p1)
           ENDIF
           s%x_ad_ = DBLE(spra1(1))
!          Sorry but no computation of derivatives because sa2sp_geo does not support it.
!            instead, assume that any derivative d(s)/dx is as d(sal)/dx
           s%xp_ad_(:) = ssal%xp_ad_(:)
        ELSE
           s = ssal
        ENDIF

!       Salinity and simply related values
        s = ssal
        s2=s*s
        sqrts=SQRT(s)
        s15=s**1.5d0
        scl=s/1.80655d0

!       Ionic strength:
        is = 19.924d0*s/(1000.0d0 - 1.005d0*s)
        is2 = is*is
        sqrtis = SQRT(is)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
        St(i) = 0.14d0 * scl/96.062d0

!       Fluoride:  Riley (1965)
        Ft(i) = 0.000067d0 * scl/18.9984d0

!       Boron:
        IF (trim(opB) == 'l10') THEN
!          New formulation from Lee et al (2010)
           Bt(i) = 0.0002414d0 * scl/10.811d0
        ELSEIF (trim(opB) == 'u74') THEN
!          Classic formulation from Uppström (1974)
           Bt(i) = 0.000232d0  * scl/10.811d0
        ELSE
           PRINT *,"optB must be 'l10' or 'u74'"
           STOP
        ENDIF

!       K0 (K Henry)
!       CO2(g) <-> CO2(aq.)
!       K0  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
        IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
           tk0 = tk                   !in situ temperature (K) for K0 calculation
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Ppot'    .OR. trim(opGAS) == 'ppot') THEN
           tk0 = dtempot + 273.15d0   !potential temperature (K) for K0 calculation as needed for potential fCO2 & pCO2
           Ptot = Patmd               !total pressure (in atm) = atmospheric pressure ONLY
        ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
           tk0 = tk                     !in situ temperature (K) for K0 calculation
           Phydro_atm = prb / 1.01325d0 !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
           Ptot = Patmd + Phydro_atm    !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
        ELSE
           PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
           STOP
        ENDIF
        tmp = 9345.17d0/tk0 - 60.2409d0 + 23.3585d0 * LOG(tk0/100.0d0)
        nK0we74 = tmp + s*(0.023517d0 - 0.00023656d0*tk0 + 0.0047036e-4_r8*tk0*tk0)
        K0(i) = EXP(nK0we74)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
        IF (trim(opK1K2) == 'l') THEN
!         Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total pH scale)
          K1(i) = ten**(-1.0d0*(3633.86d0*invtk - 61.2172d0 + 9.6777d0*dlogtk  &
                  - 0.011555d0*s + 0.0001152d0*s2))
          K2(i) = ten**(-1*(471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk      &
                  - 0.01781d0*s + 0.0001122d0*s2))
        ELSEIF (trim(opK1K2) == 'm10') THEN
!         Millero (2010, Mar. Fresh Wat. Res.) (total pH scale)
!         pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
!         ma1 = 13.4051d0*sqrts + 0.03185d0*s - (5.218e-5)*s2
!         mb1 = -531.095d0*sqrts - 5.7789d0*s
!         mc1 = -2.0663d0*sqrts
!         pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
!         K1(i) = ten**(-pK1) 

!         pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
!         ma2 = 21.5724d0*sqrts + 0.1212d0*s - (3.714e-4)*s2
!         mb2 = -798.292d0*sqrts - 18.951d0*s
!         mc2 = -3.403d0*sqrts
!         pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
!         K2(i) = ten**(-pK2)

!         Millero (2010, Mar. Fresh Wat. Res.) (seawater pH scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.4038d0*sqrts + 0.03206d0*s - (5.242e-5)*s2
          mb1 = -530.659d0*sqrts - 5.8210d0*s
          mc1 = -2.0664d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = ten**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.3728d0*sqrts + 0.1218d0*s - (3.688e-4)*s2
          mb2 = -788.289d0*sqrts - 19.189d0*s
          mc2 = -3.374d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = ten**(-pK2)
        ELSE
           PRINT *, "optK1K2 must be either 'l' or 'm10'"
           STOP
        ENDIF

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
        Kb(i) = EXP((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +  &
                1.728d0*s15 - 0.0996d0*s2)*invtk +              &
                (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +   &
                (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *      &
                dlogtk + 0.053105d0*sqrts*tk)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equation's 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K1p(i) = EXP(-4576.752d0*invtk + 115.540d0 - 18.453d0*dlogtk +  &
                 (-106.736d0*invtk + 0.69171d0) * sqrts +             &
                 (-0.65643d0*invtk - 0.01844d0) * s)

!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equation's 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K2p(i) = EXP(-8814.715d0*invtk + 172.1033d0 - 27.927d0*dlogtk +  &
                 (-160.340d0*invtk + 1.3566d0)*sqrts +                 &
                 (0.37335d0*invtk - 0.05778d0)*s)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equation's 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K3p(i) = EXP(-3070.75d0*invtk - 18.126d0 +            &
                 (17.27039d0*invtk + 2.81197d0) *             &
                 sqrts + (-44.99486d0*invtk - 0.09984d0) * s)

!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equation's 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Ksi(i) = EXP(-8904.2d0*invtk  + 117.400d0 - 19.334d0*dlogtk +  &
                 (-458.79d0*invtk + 3.5913d0) * sqrtis +             &
                 (188.74d0*invtk - 1.5998d0) * is +                  &
                 (-12.1652d0*invtk + 0.07871d0) * is2 +              &
                 LOG(1.0 - 0.001005d0*s))

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equation's 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Kw(i) = EXP(-13847.26d0*invtk + 148.9802d0 - 23.6521d0*dlogtk +  &
               (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *          &
               sqrts - 0.01615d0 * s)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks_0p = EXP(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk          &
                + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk) * sqrtis  &
                + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is      &
                - 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is2              &
                + LOG(1.0d0 - 0.001005d0*s))

!       Kf = [H][F]/[HF]
!       (total scale)
        IF (trim(opKf) == 'dg') THEN
!          Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
           Kf_0p = EXP(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +  &
                   LOG(1.0d0 - 0.001005d0*s) +                     &
                   LOG(1.0d0 + St(i)/Ks_0p))
        ELSEIF (trim(opKf) == 'pf') THEN
!          Perez and Fraga (1987) - Already on Total scale (no need for last line above)
!          Formulation as given in Dickson et al. (2007)
           Kf_0p = EXP(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)
        ELSE
           PRINT *, "optKf must be either 'dg' or 'pf'"
           STOP
        ENDIF

!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
        Kspc(i) = ten**(-171.9065d0 - 0.077993d0*tk + 2839.319d0/tk    &
                 + 71.595d0*LOG10(tk)                             &
                 + (-0.77712d0 + 0.0028426d0*tk + 178.34d0/tk)*sqrts  &
                 -0.07711d0*s + 0.0041249d0*s15 )


!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
        Kspa(i) = ten**(-171.945d0 - 0.077993d0*tk + 2903.293d0/tk &
             +71.595d0*LOG10(tk) &
             +(-0.068393d0 + 0.0017276d0*tk + 88.135d0/tk)*sqrts &
             -0.10018d0*s + 0.0059415d0*s15 )

!       Pressure effect on K0 based on Weiss (1974, equation 5)
        Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
        vbarCO2 = 32.3_r8           ! partial molal volume (cm3 / mol) from Weiss (1974, Appendix, paragraph 3)
        K0(i) = K0(i) * exp( ((1-Ptot)*vbarCO2)/(Rgas_atm*tk0) )   ! Weiss (1974, equation 5)

!       Pressure effect on all other K's (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
        DO ipc = 1, 12
           deltav(ipc)  =  a0(ipc) + a1(ipc) *t + a2(ipc) *t*t
           deltak(ipc)   = (b0(ipc)  + b1(ipc) *t + b2(ipc) *t*t)
           lnkpok0(ipc)  = (-(deltav(ipc)) &
                +(0.5d0*deltak(ipc) * prb) &
                )                         * prb/(R*tk)
        END DO

!       Pressure correction on Ks (Free scale)
        Ks(i) = Ks_0p*EXP(lnkpok0(5))
!       Conversion factor total -> free scale
        total2free     = 1.d0/(1.d0 + St(i)/Ks(i))   ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
        total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
        Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
        Kf(i) = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
        Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total hydrogen (pH) scales
        free2SWS  = 1.d0 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
        total2SWS = total2free * free2SWS                          ! KSWS = Ktotal*total2SWS
        SWS2total = 1.d0 / total2SWS
!       Conversion at pressure zero
        free2SWS_0p  = 1.d0 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  ! using Kf on free scale
        total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
        IF (trim(optK1K2) == 'l') THEN
          K1(i)  = K1(i)*total2SWS_0p
          K2(i)  = K2(i)*total2SWS_0p
          !This conversion is unnecessary for the K1,K2 from Millero (2010),
          !since we use here the formulation already on the seawater scale
        ENDIF
        Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw

!       Other contants (keep on another scale):
!          - K0         (independent of pH scale, already pressure corrected)
!          - Ks         (already on Free scale;   already pressure corrected)
!          - Kf         (already on Total scale;  already pressure corrected)
!          - Kspc, Kspa (independent of pH scale; pressure-corrected below)

!       Perform actual pressure correction (on seawater scale)
        K1(i)   = K1(i)*EXP(lnkpok0(1))
        K2(i)   = K2(i)*EXP(lnkpok0(2))
        Kb(i)   = Kb(i)*EXP(lnkpok0(3))
        Kw(i)   = Kw(i)*EXP(lnkpok0(4))
        Kspc(i) = Kspc(i)*EXP(lnkpok0(7))
        Kspa(i) = Kspa(i)*EXP(lnkpok0(8))
        K1p(i)  = K1p(i)*EXP(lnkpok0(9))
        K2p(i)  = K2p(i)*EXP(lnkpok0(10))
        K3p(i)  = K3p(i)*EXP(lnkpok0(11))
        Ksi(i)  = Ksi(i)*EXP(lnkpok0(12))

!       Convert back to original total scale:
        K1(i)  = K1(i) *SWS2total
        K2(i)  = K2(i) *SWS2total
        K1p(i) = K1p(i)*SWS2total
        K2p(i) = K2p(i)*SWS2total
        K3p(i) = K3p(i)*SWS2total
        Kb(i)  = Kb(i) *SWS2total
        Ksi(i) = Ksi(i)*SWS2total
        Kw(i)  = Kw(i) *SWS2total

     ELSE

        K0(i)   = 1.e20_r8
        K1(i)   = 1.e20_r8
        K2(i)   = 1.e20_r8
        Kb(i)   = 1.e20_r8
        Kw(i)   = 1.e20_r8
        Ks(i)   = 1.e20_r8
        Kf(i)   = 1.e20_r8
        Kspc(i) = 1.e20_r8
        Kspa(i) = 1.e20_r8
        K1p(i)  = 1.e20_r8
        K2p(i)  = 1.e20_r8
        K3p(i)  = 1.e20_r8
        Ksi(i)  = 1.e20_r8
        Bt(i)   = 1.e20_r8
        Ft(i)   = 1.e20_r8
        St(i)   = 1.e20_r8

     ENDIF

  END DO

  RETURN
END SUBROUTINE constants_DNAD
END MODULE mconstants
