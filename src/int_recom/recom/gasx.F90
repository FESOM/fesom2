!> \file gasx.f90
!! \BRIEF 
!> Module with routines needed to compute gas exchange (flxco2, scco2, atmospheric xCO2 and pCO2)
MODULE gasx
CONTAINS
!>    Computes air-sea CO2 flux & surface-ocean carbonate system vars (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!!    from T, S, P, ALK, DIC, total inorganic silicon, total inorganic phosphorus, all as 1-D arrays
SUBROUTINE flxco2(co2flux, co2ex, dpco2,                                                    &
                  ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis, K0,  &
                  temp, sal, alk, dic, sil, phos, kw660, xco2, Patm, dz1, N, lon, lat,      &
                  optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS          )
  !   Purpose:
  !     Computes air-sea CO2 flux & surface ocean carbonate system vars (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
  !     from T, S, P, ALK, DIC, total inorganic silicon, total inorganic phosphorus, all as 1-D arrays
  !
  !     INPUT variables:
  !     ================
  !     kw660   = gas transfer velocity (piston velocity) for CO2 [m/s] 
  !               without T-dependant Schmidt number correction
  !               BUT accounting for sea ice fraction (See OCMIP2 design & HOWTO documents for details)
  !     xco2    = atmospheric mole fraction of CO2 [ppm]
  !     Patm    = atmospheric pressure at surface [atm]
  !     dz1     = depth of the top vertical layer of the model [m]
  !     temp    = surface potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = surface in situ temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = surface salinity in [psu] (practical sal.) or [g/kg] (absolute sal.)
  !     alk     = surface total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =                             [eq/kg]  with optCON = 'mol/kg'
  !     dic     = surface dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = surface total inorganic silicon   [mol/m^3] with optCON = 'mol/m3'
  !             =                                   [mol/kg]  with optCON = 'mol/kg'
  !     phos    = surface total inorganinc phosphorus [mol/m^3] with optCON = 'mol/m3'
  !             =                                     [mol/kg]  with optCON = 'mol/kg'
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC and ALK given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC and ALK given in mol/m^3 (std MODEL units)
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
  !                **** Orr et al. (GMDD, 2014) identify large discrepancies between packages w/ this option
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
  !     lat:  latitude in degrees North
  !     ----------
  !        lat and lon are optional; they may be used when optS is "Sabs" as conversion parameters 
  !           from Absolute to Practical Salinity.
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
  !        In that case, do not pass optional parameter 'lon' an 'lat'
  !
  !        When neither chemical composition nor location are known, an arbitrary geographic point is chosen:
  !        mid equatorial Atlantic. Note that this implies an error on computed practical salinity up to 0.02 psu.
  !        In that case, do pass parameter 'lon' and 'lat' and set each of their elements to 1.e20.
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     co2flux = air-to-sea flux of CO2 [mol/(m^2 * s)]
  !     co2ex = time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
  !     dpco2 = difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
  !     ph   = pH on total scale
  !     pco2 = oceanic partial pressure of CO2 (uatm)
  !     fco2 = oceanic fugacity of CO2 (uatm)
  !     co2  = aqueous CO2 concentration [mol/m^3]
  !     hco3 = bicarbonate (HCO3-) concentration [mol/m^3]
  !     co3  = carbonate (CO3--) concentration [mol/m^3]
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state
  !     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  !     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  !     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  !     tempis  = in-situ temperature [degrees C]
  !     K0 = CO2 solubility [(mol/kg) / atm] 


#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE meos
  USE gsw_mod_toolbox, only: gsw_t_from_ct
  USE mvars
  USE mp2fCO2

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
  !> gas transfer velocity (piston velocity) at a Schmidt number of 660 <b>[m/s]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: kw660
  !> atmospheric mole fraction of CO2 <b>[ppm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: xco2
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> thickness of the surface layer of the model <b>[m]</b> 
  REAL(kind=rx), INTENT(in) :: dz1

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
!f2py character*6 optional, intent(in) :: optCON='mol/m3'
  CHARACTER(6), OPTIONAL, INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
!f2py character*7 optional, intent(in) :: optT='Tinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
!f2py character*2 optional, intent(in) :: optP='m'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optP
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
  !> latitude <b>[degrees north]</b>
!f2py real(8) optional, intent(in), dimension(n) :: lat = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in),    DIMENSION(N) :: lat

! Output variables:
  !> air-to-sea CO2 flux <b>[mol/(m^2 * s)]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co2flux
  !> rate of change of surface DIC concentration <b>[mol/(m^3 * s)]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: co2ex
  !> difference of surface ocean pCO2 minus atmospheric pCO2 <b>[uatm]</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: dpCO2
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
  !> K0, solubility of CO2  \b <b>[(mol/kg) / atm]</b>                                               
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: K0

! Local variables
  REAL(kind=r8) :: tk, invtk, dtemp
  REAL(kind=r8) :: tmp, co2star, co2starair, kwco2 ! K0
  REAL(kind=rx), DIMENSION(N) :: pCO2atm, fCO2atm
  REAL(kind=rx), DIMENSION(N) :: depth0, lat0
 
  ! local 1-long array version of scalar variables
  REAL(kind=r8), DIMENSION(1) :: sa1, s1, p1, lon1, lat1
  REAL(kind=r8), DIMENSION(1) :: tc1, ta1, sit1, nt1
  
  ! practical salinity [psu] computed when absolute saliniry is given 
  REAL(kind=rx), DIMENSION(N) :: salprac

  INTEGER :: i
  INTEGER :: kcomp

! Optional arguments: if not present use defaults (defined below)
  CHARACTER(6) :: opCON
  CHARACTER(7) :: opT
  CHARACTER(2) :: opP
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaration
!        Those !f2py statements should be consistent w/ defaults below
  IF (PRESENT(optCON)) THEN
    opCON = optCON
  ELSE
!   Default (typical for models, not data)
    opCON = 'mol/m3'
  ENDIF
  IF (PRESENT(optT)) THEN
    opT = optT
  ELSE
!   Default (this option is irrelevant for surface values, 0 m)
    opT = 'Tinsitu'
  ENDIF
  IF (PRESENT(optP)) THEN
    opP = optP
  ELSE
!   Default (this option is irrelevant for surface values, 0 m)
    opP = 'm'
  ENDIF
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

  depth0 = fco2 * 0.0
  lat0   = depth0

! Compute derived variables from input (DIC, ALK, ...)
  IF (PRESENT(lat)) THEN
    CALL vars_sprac(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
              temp, sal, alk, dic, sil, phos, Patm, depth0, lat, N,                     &
              opCON, opT, opP, opB, opK1K2, opKf, opGAS, opS, lon, salprac              )
  ELSE
    CALL vars_sprac(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
              temp, sal, alk, dic, sil, phos, Patm, depth0, lat0, N,                    &
              opCON, opT, opP, opB, opK1K2, opKf, opGAS, opS, lon, salprac              )
  ENDIF
  IF (TRIM(opS) == 'Sprc') THEN
    salprac(:) = sal(:)
  ENDIF
  
! Compute pCO2atm [uatm] from xCO2 [ppm], atmospheric pressure [atm], & vapor pressure of seawater
! pCO2atm = (Patm - pH20(i)) * xCO2,   where pH20 is the vapor pressure of seawater [atm]
  CALL x2pCO2atm(xco2, temp, salprac, Patm, N, pco2atm)

! Compute fCO2atm [uatm] from pCO2atm [uatm] & fugacity coefficient [unitless]
! fCO2atm = pCO2atm * fugcoeff,   where fugcoeff= exp(Patm*(B + 2.0*xc2*Del)/(R*tk) )
  CALL p2fCO2(pco2atm, temp, Patm, depth0, N, fco2atm)

! Compute flux, absolute rate of change of surface DIC, & Delta pCO2
  DO i = 1, N
     IF (co3(i) .eq. 1.e20_rx) THEN
!       Masked values (land)
        co2flux(i) = 1.e20_rx
        co2ex(i)   = 1.e20_rx
        dpco2(i)   = 1.e20_rx
     ELSE
        dtemp = DBLE(temp(i))
        tk = dtemp + 273.15d0
        invtk = 1.0d0/tk

!       Transfer velocity for CO2 in m/s (see equation [4] in OCMIP2 design document & OCMIP2 Abiotic HOWTO)
        kwco2 = DBLE(kw660(i)) * (660/scco2(dtemp))**0.5

!       Surface K0 [(mol/kg) / atm] at T, S of surface water
        tmp = 9345.17d0*invtk - 60.2409d0 + 23.3585d0 * LOG(tk/100.0d0)
        K0 = EXP( tmp + DBLE(salprac(i))*(0.023517d0 - 0.00023656d0*tk + 0.0047036e-4_r8*tk*tk) )

!       "Atmospheric" [CO2*], air-sea CO2 flux, sfc DIC rate of change, & Delta pCO2
        co2starair = K0(i) * DBLE(fco2atm(i)) * 1.0e-6_r8 * DBLE(rhoSW(i)) !Equil. [CO2*] for atm CO2 at Patm & sfc-water T,S [mol/m3]
        co2star = DBLE(co2(i))                                          !Oceanic [CO2*] in [mol/m3] from vars.f90
        co2flux(i) = SGLE(kwco2 * (co2starair - co2star))               !Air-sea CO2 flux [mol/(m2 * s)]
!       the conversion from co2flux to impact on dic is done in recom_forcing/recom_sms 
        co2ex(i) = co2flux(i) ! / dz1                                     !Change in sfc DIC due to gas exchange [mol/[m3 * s)]
        dpco2(i) = pco2(i) - pco2atm(i)                                 !Delta pCO2 (oceanic - atmospheric pCO2) [uatm]
     ENDIF

  END DO

  RETURN
END SUBROUTINE flxco2

!>    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure
SUBROUTINE pCO2atm2xCO2(pCO2atm, temp, salt, Patm, N, xCO2)
  !    Purpose:
  !    Compute xCO2 from arrays of pCO2atm, in situ T, S, & atm pressure

  USE msingledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> atmospheric partial pressure of CO2 [uatm] 
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: pCO2atm
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> mole fraction of CO2 [ppm]
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: xCO2

! LOCAL variables:
  REAL(kind=r8) :: dpCO2atm, dPatm
  REAL(kind=r8), DIMENSION(N) :: pH20
  REAL(kind=r8) :: dxCO2

  INTEGER :: i

  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dpCO2atm  = DBLE(pCO2atm(i))
     dPatm     = DBLE(Patm(i))
     dxCO2     = dpCO2atm / (dPatm - pH20(i))
     xCO2(i) = SGLE(dxCO2)
  END DO

  RETURN
END SUBROUTINE pCO2atm2xCO2

!>    Compute piston velolicty kw660 (at 25 C) from wind speed
SUBROUTINE pistonvel(windspeed, Fice, N, kw660)
  !    Purpose:
  !    Compute piston velocity from wind speed, BUT without Schmidt number temperature correction (Sc differs each gas)

  USE msingledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> wind speed at 10-m height
  REAL(kind=8), INTENT(in), DIMENSION(N) :: windspeed
!jh  REAL(kind=r8), INTENT(in), DIMENSION(N) :: windspeed
  !> modeled sea-ice cover: fraction of grid cell, varying between 0.0 (no ice) and 1.0 (full cover)
  REAL(kind=8), INTENT(in), DIMENSION(N) :: Fice
!  REAL(kind=r8), INTENT(in), DIMENSION(N) :: Fice
!f2py optional , depend(windspeed) :: n=len(windspeed)

! OUTPUT variables:
  !> piston velocity at 25°C [m/s], uncorrected by the Schmidt number for different temperatures
  REAL(kind=8), INTENT(out), DIMENSION(N) :: kw660
!jh  REAL(kind=rx), INTENT(out), DIMENSION(N) :: kw660

! LOCAL variables:
  REAL(kind=r8) :: a, xfac

  INTEGER :: i

! Coefficient from Wanninkhof (2014, Limnol. Oceanograph. Methods, 12, 351-362)
  a = 0.251d0
! where, by convention, the units of k are in cm h–1 and U is in m s–1. Thus, the units of the coefficient 0.251 are (cm h–1) (ms–1)–2

! Conversion factor to convert from cm/h to m/s
  xfac = 0.01d0 / 3600d0
  
  DO i = 1,N
     kw660(i) = a * windspeed(i)**2 * (1.0d0 - Fice(i)) * xfac
  END DO

  RETURN
END SUBROUTINE pistonvel

!>    Compute Schmidt number for CFC11 in seawater from temperature
FUNCTION sccfc11(temp)

!  Compute Schmidt number of CFC11 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: temp
   REAL(kind=r8) :: sccfc11

   sccfc11 = 3579.2_r8 - 222.63_r8*temp + 7.5749_r8*temp**2 - 0.14595_r8*temp**3  + 0.0011874_r8*temp**4
   
   RETURN
END FUNCTION sccfc11

!>    Compute Schmidt number for CFC12 in seawater from temperature
FUNCTION sccfc12(Tc)

!  Compute Schmidt number of CFC12 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: Tc
   REAL(kind=r8) :: sccfc12

   sccfc12 = 3828.1_r8 - 249.86_r8*Tc + 8.7603_r8*Tc**2 - 0.1716_r8*Tc**3   + 0.001408_r8*Tc**4
      
   RETURN
END FUNCTION sccfc12

!>    Compute Schmidt number for SF6 in seawater from temperature
FUNCTION scsf6(Tc)

!  Compute Schmidt number of SF6 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: Tc
   REAL(kind=r8) :: scsf6

   scsf6 = 3177.5_r8 - 200.57_r8*Tc + 6.8865_r8*Tc**2 - 0.13335_r8*Tc**3 + 0.0010877_r8*Tc**4
      
   RETURN
END FUNCTION scsf6

!>    Compute Schmidt number for CO2 in seawater from temperature
FUNCTION scco2(Tc)

!  Compute Schmidt number of CO2 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: Tc
   REAL(kind=r8) :: scco2

   scco2 = 2116.8_r8 - 136.25_r8*Tc + 4.7353_r8*Tc**2 - 0.092307_r8*Tc**3 + 0.0007555_r8*Tc**4

   RETURN
END FUNCTION scco2

!>    Compute Schmidt number for O2 in seawater from temperature
FUNCTION sco2(Tc)

!  Compute Schmidt number of O2 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

   USE msingledouble
   IMPLICIT NONE

!  Input & output variables:
   REAL(kind=r8), INTENT(in) :: Tc
   REAL(kind=r8) :: sco2

   sco2 = 1920.4_r8 - 135.6_r8*Tc  + 5.2122_r8*Tc**2 - 0.10939_r8*Tc**3  + 0.00093777_r8*Tc**4

   RETURN
END FUNCTION sco2

!>    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure
SUBROUTINE x2pCO2atm(xCO2, temp, salt, Patm, N, pCO2atm)
  !    Purpose:
  !    Compute pCO2atm from arrays of xCO2, in situ T, S, & atm pressure

  USE msingledouble

  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> mole fraction of CO2 [ppm]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: xCO2
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !> atmospheric pressure [atm]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm] 
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pCO2atm

! LOCAL variables:
  REAL(kind=r8) :: dxCO2, dPatm
  REAL(kind=r8), DIMENSION(N) :: pH20
  REAL(kind=r8) :: dpCO2atm

  INTEGER :: i

! Compute vapor pressure of seawater [in atm]
  call vapress(temp, salt, N, pH20)

  DO i = 1,N
     dxCO2     = DBLE(xCO2(i))
     dPatm     = DBLE(Patm(i))
     dpCO2atm = (dPatm - pH20(i)) * dxCO2
     pCO2atm(i) = SGLE(dpCO2atm)
  END DO

  RETURN
END SUBROUTINE x2pCO2atm

!>    Compute solubilities of CFC-11, CFC-12, SF6, CO2, and N20 at 1 atm pressure, i.e., Phi0 (atm), in mol L-1 atm-1 
SUBROUTINE phizero(gasname, temp, salt, N, phi0)
  !    Purpose:
  !    Compute solubilities of CFC-11, CFC-12, SF6, CO2, and N20 at 1 atm pressure, i.e., Phi0 (atm), in mol L-1 atm-1 
  !    Phi0 = K0 * Cf * (Pa0 - pH20),
  !    where
  !    * K0 is the solubility of a given gas,
  !    * Cf is its fugacity coefficient,
  !    * Pa0 is 1 atm of total atmospheric pressure, and
  !    * pH20 is the water vapor pressure at saturation (also in atm), all described in Orr et al. (2016, GMDD).

  !    Usage: this routine must be called once for each gas, e.g.,
  !           call phizero('co2', temp, salt, 20, phi0_co2)
  !           call phizero('sf6', temp, salt, 20, phi0_sf6)
  
  ! James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, Université Paris Saclay, France
  ! 5 August 2016
  
  USE msingledouble
  IMPLICIT NONE

  INTEGER, PARAMETER :: ngas = 5
  
  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !f2py optional , depend(temp) :: n=len(temp)
  !> name of gas: 'cfc11', 'cfc12', 'sf6', 'co2', or 'n2o'
  CHARACTER(*), INTENT(in) :: gasname

! OUTPUT variables:
  !> solubility of gas in seawater, Phi0 = K0 * Cf * (Pa0 - pH20), with units of [mol L-1 atm-1] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) ::  phi0

! LOCAL variables:
  ! Coefficients to compute Phi0
  REAL(kind=r8), DIMENSION(7,ngas) :: Phi0coeffs  !Array of coeffs (a1, a2, a3, a4, b1, b2, b3) times 3 gases
  REAL(kind=r8), DIMENSION(ngas,7) :: a            !Transpose of above array
  ! Absolute temperature and salinity
  REAL(kind=r8) :: dsalt, tk, tk100, ln_phi0

  INTEGER :: i, ig

! Coefficients to compute solubilities (Phi0) of CFC-11, CFC-12, SF6, CO2, and N20
  ! ------------------------------------------------------------------------------------------
  ! Coefficients a1, a2, a3, a4, b1, b2, & b3 from Table 3 of Orr et al. (2016).
  ! These are used to compute Phi0 for each gas following Eq. (15) in same paper.
  ! ------------------------------------------------------------------------------------------
  !                     a1        a2        a3       a4        b1          b2         b3
  ! ------------------------------------------------------------------------------------------
  DATA Phi0coeffs / -229.9261, 319.6552, 119.4471, -1.39165, -0.142382,   0.091459, -0.0157274,  & !CFC-11
                    -218.0971, 298.9702, 113.8049, -1.39165, -0.143566,   0.091015, -0.0153924,  & !CFC-12
                     -80.0343,  117.232,  29.5817,  0.0,      0.0335183, -0.0373942, 0.00774862, & !SF6
                    -160.7333, 215.4152,  89.8920, -1.47759,  0.029941,  -0.027455,  0.0053407,  & !CO2
                    -165.8806, 222.8743,  92.0792, -1.48425, -0.056235,   0.031619, -0.0048472   / !N2O
  ! -------------------------------------------------------------------------------
  ! * The originial publications providing these coefficients are
  ! * - Warner and Weiss (1985) for CFC-11 and CFC-12,
  ! * - Bullister et al. (2002) for SF6, and
  ! * - Weiss and Price (1980) for CO2 and N2O.

  a = TRANSPOSE(Phi0coeffs)

  ! Set value of the index that refers to to each gasname
  SELECT CASE (trim(gasname))
      CASE ('cfc11')
          ig = 1
      CASE ('cfc12')
          ig = 2
      CASE ('sf6')
          ig = 3
      CASE ('co2')
          ig = 4
      CASE ('n2o')
          ig = 5
      CASE DEFAULT
          PRINT *,"ERROR in 'phizero' routine in gasx.f90:"
          PRINT *,"'gasname' input var must be one of the following: 'cfc11', 'cfc12', 'sf6', 'co2', or 'n2o'"
          STOP
  END SELECT

! Compute phi0  
  DO i = 1, N
      tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
      dsalt = DBLE(salt(i))             
      tk100 = tk/100.d0
      ! Phi0 for the chosen gas:
      ln_phi0 = a(ig,1) + a(ig,2)/tk100 + a(ig,3)*log(tk100) + a(ig,4)*tk100**2  &
              + dsalt * (a(ig,5) + a(ig,6)*tk100 + a(ig,7)*tk100**2 )
      phi0(i) = exp(ln_phi0)
  END DO

  RETURN
END SUBROUTINE phizero

!>    Compute K' of CFC-11, CFC-12, and SF6; units are in mol L-1 atm-1 
SUBROUTINE kprime(gasname, temp, salt, N, kp)
  !    Purpose:
  !    Compute K' of CFC-11, CFC-12, and SF6; units are in mol L-1 atm-1 
  !    To avoid confusion this can be used directly in the ocean solubility equation
  !    (Eq 21 from Orr et al., 2016) and in the atmospheric saturation equation (Eq 14),
  !    The atmospheric saturation would then need to separately account for atmospheric pressure 
  !    and humidity (vapress routine).

  !    Usage: this routine should be called once for each gas, e.g.,
  !           call kprime('cfc11', temp, salt, 20, kp_cfc11)
  !           call kprime('cfc12', temp, salt, 20, kp_cfc12)
  !           call kprime('sf6',   temp, salt, 20, kp_sf6)
  
  ! James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, Université Paris Saclay, France
  ! 8 August 2016
  
  USE msingledouble
  IMPLICIT NONE

  INTEGER, PARAMETER :: ngas = 3
  
  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !f2py optional , depend(temp) :: n=len(temp)

  !> name of gas: 'cfc11', 'cfc12', or 'sf6'
  CHARACTER(*), INTENT(in) :: gasname

! OUTPUT variables:
  !> solubility of gas in seawater, in units of [mol L-1 atm-1] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) ::  kp

! LOCAL variables:
  ! Coefficients to compute K' 
  REAL(kind=r8), DIMENSION(6,ngas) :: Kpcoeffs  !Array of coeffs (a1, a2, a3, b1, b2, b3) times 3 gases
  REAL(kind=r8), DIMENSION(ngas,6) :: a         !Transpose of above array
  ! Absolute temperature and salinity
  REAL(kind=r8) :: dsalt, tk, tk100, ln_kp

  INTEGER :: i, ig

  ! -------------------------------------------------------------------------------
  ! Input array of coefficients (top half of Table 3 of Orr et al. (2016, GMDD)
  ! These are used to compute K' for each gas following Eq. (24) in same paper.
  ! -------------------------------------------------------------------------------
  !                   a1        a2       a3       b1          b2         b3
  ! -------------------------------------------------------------------------------
  DATA Kpcoeffs / -134.1536, 203.2156, 56.2320, -0.144449,   0.092952, -0.0159977,   &  !CFC-11
                  -122.3246, 182.5306, 50.5898, -0.145633,   0.092509, -0.0156627,   &  !CFC-12
                   -96.5975, 139.883,  37.8193,  0.0310693, -0.0356385, 0.00743254   /  !SF6
  ! -------------------------------------------------------------------------------
  ! * The originial publications providing these coefficients are
  ! * - Warner and Weiss (1985) for CFC-11 and CFC-12, and
  ! * - Bullister et al. (2002) for SF6, and

  a = TRANSPOSE(Kpcoeffs)

  write(*,*)'gasname = ', gasname
! Set value of the index that refers to to each gasname
  SELECT CASE (trim(gasname))
      CASE ('cfc11')
          ig = 1
      CASE ('cfc12')
          ig = 2
      CASE ('sf6')
          ig = 3
      CASE DEFAULT
          PRINT *,"ERROR in 'kprime' routine in gasx.f90:"
          PRINT *,"'gasname' input var must be one of the following: 'cfc11', 'cfc12', or 'sf6'"
          STOP
  END SELECT

! Compute K' 
  DO i = 1, N
      tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
      dsalt = DBLE(salt(i))             
      tk100 = tk/100.d0
      ! Phi0 for the chosen gas:
      ln_kp = a(ig,1) + a(ig,2)/tk100 + a(ig,3)*log(tk100)                   &
            + dsalt * (a(ig,4) + a(ig,5)*tk100 + a(ig,6)*tk100**2 )
      kp(i) = exp(ln_kp)
  END DO

  RETURN
END SUBROUTINE kprime

!>    Compute K0 of CO2 and N20; units are in mol L-1 atm-1 
SUBROUTINE kzero(gasname, temp, salt, N, k0)
  !    Purpose:
  !    Compute K0 of CO2 and N20; units are in mol L-1 atm-1 
  !    Uses same form of the equation 24 (given for K') in Orr et al. (2016)

  !    Usage: this routine should be called once for each gas, e.g.,
  !           call kzero('co2', temp, salt, 20, k0_co2)
  !           call kzero('n2o', temp, salt, 20, kp_n2o)
  
  ! James Orr, LSCE/IPSL, CEA-CNRS-UVSQ, Université Paris Saclay, France
  ! 8 August 2016
  
  USE msingledouble
  IMPLICIT NONE

  INTEGER, PARAMETER :: ngas = 2
  
  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !f2py optional , depend(temp) :: n=len(temp)

  !> name of gas: 'co2' or 'n2o'
  CHARACTER(*), INTENT(in) :: gasname

! OUTPUT variables:
  !> solubility of gas in seawater, in units of [mol L-1 atm-1] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) ::  k0

! LOCAL variables:
  ! Coefficients to compute K' 
  REAL(kind=r8), DIMENSION(6,ngas) :: K0coeffs  !Array of coeffs (a1, a2, a3, b1, b2, b3) times 3 gases
  REAL(kind=r8), DIMENSION(ngas,6) :: a         !Transpose of above array
  ! Absolute temperature and salinity
  REAL(kind=r8) :: dsalt, tk, tk100, ln_k0

  INTEGER :: i, ig

  ! -------------------------------------------------------------------------------
  ! Input array of coefficients (bottom half of Table 3 of Orr et al. (2016, GMDD)
  ! These are used to compute K' for each gas following Eq. (24) in same paper.
  ! -------------------------------------------------------------------------------
  !                   a1        a2       a3       b1          b2         b3
  ! -------------------------------------------------------------------------------
  DATA K0coeffs /  -58.0931,  90.5069, 22.2940,  0.027766,  -0.025888,  0.0050578,   &  ! CO2
                   -62.7062,  97.3066, 24.1406, -0.058420,   0.033193,  -0.0051313   /  !N2O
  ! -------------------------------------------------------------------------------
  ! * The originial publications providing these coefficients are
  ! * - Weiss (1974, Table 1, column 1) for CO2
  ! * - Weiss and Price (1980, Table 2, column 1) for N2O

  a = TRANSPOSE(K0coeffs)

! Set value of the index that refers to to each gasname
  SELECT CASE (trim(gasname))
      CASE ('co2')
          ig = 1
      CASE ('n2o')
          ig = 2
      CASE DEFAULT
          PRINT *,"ERROR in 'kzero' routine in gasx.f90:"
          PRINT *,"'gasname' input var must be one of the following: 'co2' or 'n2o'"
          STOP
  END SELECT

! Compute K0
  DO i = 1, N
      tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
      dsalt = DBLE(salt(i))             
      tk100 = tk/100.d0
      ! k0 for the chosen gas:
      ln_k0 = a(ig,1) + a(ig,2)/tk100 + a(ig,3)*log(tk100)                   &
            + dsalt * (a(ig,4) + a(ig,5)*tk100 + a(ig,6)*tk100**2 )
      k0(i) = exp(ln_k0)
  END DO

  RETURN
END SUBROUTINE kzero

!>    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)
SUBROUTINE vapress(temp, salt, N, vpsw)
  !    Purpose:
  !    Compute vapor pressure of seawater (atm) following preocedure from Weiss & Price (1980)

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> vapor pressure of seawater [atm] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: vpsw

! LOCAL variables:
  REAL(kind=r8) :: tk, dsalt

  INTEGER :: i

  DO i = 1,N
     dsalt = DBLE(salt(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     vpsw(i) = exp(24.4543d0 - 67.4509d0*(100.0d0/tk) - 4.8489d0*log(tk/100) - 0.000544d0*dsalt)
  END DO

  RETURN
END SUBROUTINE vapress

!>    Compute O2 saturation concentration of surface seawater (mol/m3) at 1 atm (Garcia & Gordon, 1992)
SUBROUTINE o2sato(T, S, N, o2sat_molm3)
  !    Purpose:
  !    Compute O2 saturation concentration of surface seawater (mol/m3) at 1 atm (Garcia & Gordon, 1992)
  !
  !    ********************************************************************
  !    Computes the oxygen saturation concentration at 1 atm total pressure
  !    in mol/m^3 given sea surface temperature T (deg C) and salinity S (permil) 
  !
  !    From: Garcia & Gordon (1992) Oxygen solubility in seawater: better fitting equations,
  !          Limnol. Oceanogr., 37(6), 1307-1312.
  !          This routine uses:
  !          - equation (8) on page 1310
  !          - coefficients from Table 1, column 2 (Benson & Krause, [cm3/dm3], i.e, same as [ml/L])
  !
  !    *** NOTE: The "A3*ts^2" term in the equation (8) in the paper is a TYPO.    ***
  !    *** It shouldn't be there. It is not used in this routine.                  ***
  !
  !    'o2sat' is fit between T(freezing) <= T <= 40(deg C)  and  0 <= S <= 42 permil
  !
  !    CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
  !    o2sat_molm3 = 0.282015 mol/m^3
  !    ********************************************************************

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> surface temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: T
  !> surface salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: S
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> O2 saturation concentration of seawater [mol/m3] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: o2sat_molm3

! LOCAL variables:
  REAL(kind=r8) :: A0, A1, A2, A3, A4, A5
  REAL(kind=r8) :: B0, B1, B2, B3
  REAL(kind=r8) :: C0
  REAL(kind=r8) :: tmp
  REAL(kind=r8) :: o2sat_mlL 
  REAL(kind=r8) :: tt, tk, ts, ts2, ts3, ts4, ts5
  INTEGER :: i
  
  DATA A0/ 2.00907_r8   /, A1/ 3.22014_r8   /, A2/ 4.05010_r8 /,  &
       A3/ 4.94457_r8   /, A4/-2.56847E-1_r8/, A5/ 3.88767_r8 /
  DATA B0/-6.24523E-3_r8/, B1/-7.37614E-3_r8/, B2/-1.03410E-2_r8/, B3/-8.17083E-3_r8/
  DATA C0/-4.88682E-7_r8/
      
  DO i = 1, N
      tt  = 298.15_r8 - T(i)
      tk  = 273.15_r8 + T(i)
      ts  = LOG(tt/tk)

      ts2 = ts**2
      ts3 = ts**3
      ts4 = ts**4
      ts5 = ts**5

!     O2 saturation concentration (ml/L) 
      tmp  = A0 + A1*ts + A2*ts2 + A3*ts3 + A4*ts4 + A5*ts5  &
               + S(i)*(B0 + B1*ts + B2*ts2 + B3*ts3)         &
               + C0*(S(i)*S(i))
      o2sat_mlL = EXP(tmp)

!     Convert from ml/L to mol/m^3
      o2sat_molm3(i) = o2sat_mlL / 22391.6_r8*1000.0_r8
  END DO
   
  RETURN
END SUBROUTINE o2sato

!>    Compute time rate of change of O2 in the surface layer due to air-sea gas exchange [mol/(m^3 *s)].
SUBROUTINE o2flux(T, S, kw660, ppo, o2,  N, o2ex)

  !    **********************************************************************
  !    Purpose: Compute time rate of change of O2 in the surface layer due to air-sea gas exchange [mol/(m^3 *s)]
  !
  !    Input:
  !      T       model surface temperature (deg C)
  !      S       model surface salinity (permil)
  !      kw660   gas transfer velocity at a Schmidt number of 660, accounting
  !              for sea ice fraction (m/s)
  !      ppo     surface pressure divided by 1 atm.
  !      o2      surface ocean O2 concentration (mol/m^3)
  !      dz1     thickness of surface grid box (m) - taken out

  !    Output:
  !      o2ex    time rate of change of oxygen in the surface layer due
  !              to air-sea exchange (mol/m^3/s)
  !
  !    Two routines are called:
  !      sco2   - function to compute Schmidt number of oxygen
  !      o2sato - subroutine to compute oxygen saturation concentration at 1 atm (mol/m^3)
  !
  !    Numbers in brackets refer to equation numbers in OCMIP2 simulation design document
  !
  !    Original for OCMIP2: Ray Najjar, 29 January 1999
  !    Modified for OMIP:   James Orr, LSCE/IPSL France, 14 March 2015
  !    **********************************************************************

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> sea surface temperature [C]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: T
  !> sea surface salinity [psu]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: S
  !> gas transfer velocity at a Schmidt number of 660, accounting for sea ice fraction [m/s]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: kw660
  !> surface atmospheric pressure [atm]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: ppo
  !> modeled surface ocean dissolved O2 concentration [mol/m^3]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: o2
  !> thickness of surface grid box [m]
!  REAL(kind=r8), INTENT(in), DIMENSION(N) :: dz1
!f2py optional , depend(temp) :: n=len(temp)

! OUTPUT variables:
  !> rate of change of dissolved O2 in the surface layer due to air-sea O2 exchange [mol/(m^3*s)]
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: o2ex

! LOCAL variables:
  REAL(kind=r8) :: kwo2, o2sat
  REAL(kind=r8), DIMENSION(N) :: o2sat_1atm
  INTEGER :: i
  
! Dissolved O2 saturation concentraion [mol/m^3] (in equilibrium with atmosphere) at 1 atm pressure 
  CALL o2sato(T, S, N, o2sat_1atm)
!jh  CALL o2sato(SGLE(T), SGLE(S), N, o2sat_1atm)

  DO i = 1, N
!     Transfer velocity for O2 in m/s [4]
      kwo2 = (kw660(i) * (660._r8/sco2(T(i)))**0.5)
      
!     O2 saturation concentration at given atm pressure [3]
      o2sat = o2sat_1atm(i) * ppo(i)

!     Time rate of change of surface dissolved O2 due to gas exchange (mol/(m3 * s) [1]
      o2ex(i) = kwo2*(o2sat - o2(i)) !/ dz1(i) mol/m2/s as not divided by thickness of box
  END DO

  RETURN
END SUBROUTINE o2flux

END MODULE gasx

! REFERENCES
!
!  Bullister, J. L., Wisegarver, D. P., and Menzia, F. A.: The solubility of
!  sulfur hexafluoride in water and seawater, Deep-Sea Res. I, 49, 175 –187, 2002.
!
!  Orr, J. C., Najjar, R. G., Aumount, O., Bopp, L., Bullister, J. L.,
!  Danabasoglu, G., Doney, S. C., Dunne, J. P., Dutay, J.-C., Graven,
!  H., Griffies, S. M., John, J. G., Joos, F., Levin, I., Lindsay, K.,
!  Matear, R. J., McKinley, G. A., Mouchet, A., Oschlies, A., Romanou,
!  A., Schlitzer, R., Tagliabue, A., Tanhua, T., and Yool, A.:
!  Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model
!  Intercomparison Project (OMIP), Geosci. Model Dev. Discuss.,
!  doi:10.5194/gmd-2016-155, in review, 2016.
!
!  Warner, M. J. and Weiss, R. F.: Solubilities of chlorofluorocarbons 11 and 12
! in water and seawater, Deep-Sea Res. Part A., 32, 1485–1497, 1985 .
!
!  Weiss, R. F. and Price, B. A.: Nitrous oxide solubility in water and seawater,
!  Mar. Chem., 8, 347–359, 1980.

