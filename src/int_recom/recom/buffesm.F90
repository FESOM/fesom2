!> \file buffesm2.f90
!! \BRIEF 
!> Module with buffesm subroutine - compute carbonate system vars from DIC,Alk,T,S,P,nuts
MODULE mbuffesm
CONTAINS
!>    Computes buffer factors of the seawater carbonate system as defined by Egleston et al. (2010)
!!    (corrected for sign error & modified to account for effects of total dissolved inorganic P and Si acid systems) 
!!    as 1D arrays FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations (all 1-D arrays)
SUBROUTINE buffesm(gammaDIC, betaDIC, omegaDIC, gammaALK, betaALK, omegaALK, Rf, &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,             &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon      )

  !   Purpose:
  !     Computes buffer factors for seawater carbonate system as defined by Egleston et al. (2010), corrected & modified,
  !     as 1D arrays
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations (also, all 1-D arrays)

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
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

  !     OUTPUT variables:
  !     =================
  !     gammaDIC = (d ln CO2 / dDIC)^-1
  !                in units of [mol/kg] or [mol/m^3] depending on optCON
  !     betaDIC =  (d ln H+ / dDIC)^-1
  !                in units of [mol/kg] or [mol/m^3] depending on optCON
  !     omegaDIC = (d ln Omega / dDIC)^-1
  !                in units of [mol/kg] or [mol/m^3] depending on optCON
  !     gammaALK = (d ln CO2 / dAlk)^-1
  !                in units of [mol/kg] or [mol/m^3] depending on optCON
  !     betaALK  = (d ln H+ / dAlk)^-1
  !                in units of [mol/kg] or [mol/m^3] depending on optCON
  !     omegaALK = (d ln Omega / dAlk)^-1
  !     Rf =   dpCO2/pCO2 / dDIC/DIC    (i.e., the Revelle factor, unitless)


  USE msingledouble
  USE mconstants
  USE mvars

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
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature
  !>  \b 'Tcsv" for conservative temperature (in two last cases, in-situ Temp is computed, needed for models)
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
  !> choose \b 'Sprc' for practical sal. (EOS-80, default) or \b 'Sabs' for absolute salinity (TEOS-10)
!f2py character*4 optional, intent(in) :: optS='Sprc'
  CHARACTER(4), OPTIONAL, INTENT(in) :: optS
  !> longitude <b>[degrees east]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lon = -25.
  REAL(kind=rx), OPTIONAL, INTENT(in),    DIMENSION(N) :: lon

! Output variables:
! -----------------
  !> (d ln CO2 / dDIC)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: gammaDIC
  !> (d ln H+ / dDIC)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: betaDIC
  !> (d ln Omega / dDIC)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: omegaDIC
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  !> (d ln CO2 / dAlk)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: gammaALK
  !> (d ln H+ / dAlk)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: betaALK
  !> (d ln Omega / dAlk)^-1
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: omegaALK
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: Rf

! Local variables
! -----------------
  ! 1) Output from constants.f90:
  REAL(kind=r8), DIMENSION(N) :: K0
  REAL(kind=r8), DIMENSION(N) :: K1
  REAL(kind=r8), DIMENSION(N) :: K2
  REAL(kind=r8), DIMENSION(N) :: Kb
  REAL(kind=r8), DIMENSION(N) :: Kw
  REAL(kind=r8), DIMENSION(N) :: Ks
  REAL(kind=r8), DIMENSION(N) :: Kf
  REAL(kind=r8), DIMENSION(N) :: Kspc
  REAL(kind=r8), DIMENSION(N) :: Kspa
  REAL(kind=r8), DIMENSION(N) :: K1p
  REAL(kind=r8), DIMENSION(N) :: K2p
  REAL(kind=r8), DIMENSION(N) :: K3p
  REAL(kind=r8), DIMENSION(N) :: Ksi
  REAL(kind=r8), DIMENSION(N) :: St
  REAL(kind=r8), DIMENSION(N) :: Ft
  REAL(kind=r8), DIMENSION(N) :: Bt

  ! 2) Output from vars.f90:
  REAL(kind=rx), DIMENSION(N) :: ph
  REAL(kind=rx), DIMENSION(N) :: pco2
  REAL(kind=rx), DIMENSION(N) :: fco2
  REAL(kind=rx), DIMENSION(N) :: co2
  REAL(kind=rx), DIMENSION(N) :: hco3
  REAL(kind=rx), DIMENSION(N) :: co3
  REAL(kind=rx), DIMENSION(N) :: OmegaA
  REAL(kind=rx), DIMENSION(N) :: OmegaC
  REAL(kind=rx), DIMENSION(N) :: BetaD
  REAL(kind=rx), DIMENSION(N) :: rhoSW
  REAL(kind=rx), DIMENSION(N) :: p
  REAL(kind=rx), DIMENSION(N) :: tempis
  ! practical salinity [psu] computed when absolute saliniry is given 
  REAL(kind=rx), DIMENSION(N) :: salprac

  ! 3) Other Local variables (needed to compute buffer factors)
  REAL(kind=r8) :: Alkc, Borate, h, oh
  REAL(kind=r8) :: h3po4, h2po4, hpo4, po4
  REAL(kind=r8) :: sioh4, sioh3
  REAL(kind=r8) :: Pt, Sit
  REAL(kind=r8) :: Pegle, Segle, Qegle
  REAL(kind=r8) :: SegleCBW, SeglePt, SegleSit
  
  INTEGER :: i
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
!   Default is Lueker et al. (2000) for K1 & K2
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

  print *, 'Calling vars_sprac'
!  Compute carbonate system variables from DIC, ALK, T, S, nutrients, etc
!  -------------------------------------------------------------------------------------
  CALL vars_sprac(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                      &
             optCON, optT, optP, opB, opK1K2, opKf, opGAS, opS, lon, salprac      )

  IF (trim(opS) == 'Sprc') THEN
     salprac = sal
  ENDIF

! Get equilibrium constants and total concentrations of SO4, F, B
! ------------------------------------------------------------------
  CALL constants(K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,           &
                 K1p, K2p, K3p, Ksi,                               &
                 St, Ft, Bt,                                       &
                 temp, salprac, Patm,                              &
                 depth, lat, N,                                    &
                 optT, optP, opB, opK1K2, opKf, opGAS, opS, lon    )

!  Compute buffer factors
!  ----------------------
   DO i=1,N
      IF (ph(i) .EQ. 1.e20_rx) THEN
         ! Missing input data, thus mask values
         gammaDIC(i)   = 1.e20_rx
         betaDIC(i)    = 1.e20_rx
         omegaDIC(i)   = 1.e20_rx
         gammaALK(i)   = 1.e20_rx
         betaALK(i)    = 1.e20_rx
         omegaALK(i)   = 1.e20_rx
         Rf(i)         = 1.e20_rx
      ELSE
         ! Compute buffer factors

         ! Total concentrations
         Pt = phos(i)
         Sit = sil(i)

         ! Carbonate alkalinity & concentrations of H+, OH-, Borate 
         h = 10**(-ph(i))          !Total scale
         Alkc = 2*co3(i) + hco3(i)
         Borate  = Bt(i) / (1 + h/Kb(i) ) 
         oh = Kw(i) / h

         ! Phosphate species
         ! [h2po4-] = K1p * [h3po4] / [H+]
         ! [hpo4--] = K1p * K2p * [h3po4] / [H+]²
         ! [po4---] = K1p * K2p * K3p * [h3po4] / [H+]³
         ! Pt       = [h3po4] * (1 + K1p/[H+] + K1p*K2p/[H+]² + K1p*K2p*K3p/[H+]³) 

         h3po4 = Pt / (1 + K1p(i)/h + K1p(i)*K2p(i)/h**2 + K1p(i)*K2p(i)*K3p(i)/h**3)
         h2po4 = K1p(i) * h3po4 / h
         hpo4  = K2p(i) * h2po4 / h
         po4   = K3p(i) * hpo4  / h

         ! Silicate
         ! [SiO(OH)3-] = Ksi * [Si(OH)4] / [H+]
         ! Sit = [Si(OH)4] * (1 + Ksi / [H+])

         sioh4 = Sit /(1 + Ksi(i)/h)
         sioh3 = Ksi(i) * Sit/(h+Ksi(i))
      
        !Segle: (1) corrected from sign error in Egleston et al., and (2) modified to account for P & Si acid systems
        !------------------------------------------------------------------------------------------------------------
         !Segle  =   hco3(i) + 4*co3(i) + (h*Borate/(Kb(i) + h)) + h + oh   !Formula from Egleston et al (2010) with corrected sign error
                                                                            !Egleston et al ignore P & Si acid systems!
        !REMEDY:
        !Separate Segle into 3 parts:
        ! SegleCBW: contribution of C, B, and water systems (only part accounted for by Egleston et al.)
        ! SeglePt:  contribution of total phosphate     (from J.-M. Epitalon)
        ! SegleSit: contribution of total silicon       (from J.-M. Epitalon)

!         Segle  = hco3(i) + 4*co3(i) + (h*Borate/(Kb(i) + h)) + h + oh     & 
!                    + (- h3po4 * (-h2po4 - 2*hpo4 - 3*po4)                 &
!                       + hpo4  * (2*h3po4 + h2po4 - po4)                   &
!                       + 2*po4 * (3*h3po4 + 2*h2po4 + hpo4)                &
!                      ) / Pt                                               &
!                    + h * sioh3/(Ksi(i) + h)                               &   
!                  )

         ! Part 1: Carbon, Boron, Water
         SegleCBW = hco3(i) + 4*co3(i) + (h*Borate/(Kb(i) + h)) + h + oh

         ! Part 2: Phosphorus
         IF (Pt .eq. 0.0d0) THEN
            SeglePt = 0.0
         ELSE
            SeglePt = (- h3po4 * (-h2po4 - 2*hpo4 - 3*po4)                 &
                       + hpo4  * (2*h3po4 + h2po4 - po4)                   &
                       + 2*po4 * (3*h3po4 + 2*h2po4 + hpo4)                &
                      ) / Pt
         ENDIF

         ! Part 3: Silicon
         IF (Sit .EQ. 0.0d0) THEN
            SegleSit = 0.0
         ELSE
            SegleSit = h * sioh3/(Ksi(i) + h)                           
         ENDIF

!        Add the 3 parts of Segle
         Segle  = SegleCBW + SeglePt + SegleSit

         !print *, 'hco3, co3, h, Borate, Kb, oh = ', hco3(i), co3(i), h, Borate, Kb(i), oh 
         !print *, 'h3po4, h2po4, hpo4, po4 =', h3po4, h2po4, hpo4, po4
         !print *, 'Pt = ', Pt
         !print *, 'sioh3, Ksi =', sioh3, Ksi(i)

         Pegle  = 2*co2(i) + hco3(i)                                  
         Qegle  = 2*Alkc - Segle                                 

         ! Compute 6 buffer factors: 
         !*** NOTE - units of buffer factors depend on optCON (mol/m3 OR mol/kg) 
         !         - to compare with Egleston, must convert to mmol/kg, e.g., multiply factor in mol/kg by 1000 

         gammaDIC(i) = dic(i) - (Alkc*Alkc)/Segle
         gammaALK(i) = (Alkc*Alkc - dic(i)*Segle) / Alkc
         betaDIC(i)  = (dic(i)*Segle - Alkc*Alkc) / Alkc
         betaALK(i)  = Alkc*Alkc/dic(i) - Segle
         omegaDIC(i) = dic(i) - (Alkc * Pegle/Qegle)                   
         omegaALK(i) = Alkc - dic(i) * Qegle/Pegle                     

         Rf = dic(i) / gammaDIC
         !print *, 'Segle, Pegle, Qegle = ', Segle, Pegle, Qegle
         !print *, 'i, gammaDIC(i) = ', i, gammaDIC(i)
      ENDIF
   END DO

  RETURN
END SUBROUTINE buffesm
END MODULE mbuffesm
