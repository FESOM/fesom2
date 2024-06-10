!> \file varsolver.f90
!! \BRIEF 
!> Module with varsolver subroutine - solve for pH and other carbonate system variables
MODULE mvarsolver
CONTAINS
!>    Solve for pH and other carbonate system variables (with input from vars routine)
SUBROUTINE varsolver(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC,             &
                    temp, salt, ta, tc, pt, sit,                                 &
                    Bt, St, Ft,                                                  &
                    K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  & 
                    Patm, Phydro_bar, rhodum, optGAS                             )

  !   Purpose: Solve for pH and other carbonate system variables (with input from vars routine)

  !     INPUT variables:
  !     ================
  !     temp    = in situ temperature [degrees C]
  !     ta      = total alkalinity                     in [eq/m^3] or [eq/kg]   based on optCON in calling routine (vars)
  !     tc      = dissolved inorganic carbon           in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     pt      = total dissolved inorganic phosphorus in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     sit     = total dissolved inorganic silicon    in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     Bt      = total dissolved inorganic boron      computed in calling routine (vars)
  !     St      = total dissolved inorganic sulfur     computed in calling routine (vars)
  !     Ft      = total dissolved inorganic fluorine   computed in calling routine (vars)
  !     K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi 
  !     Patm    = atmospheric pressure [atm]
  !     Phydro_bar = hydrostatic pressure [bar]
  !     rhodum  = density factor as computed in calling routine  (vars)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2 (default optGAS = 'Pinsitu')
  !     ---------
  !       PRESSURE & T corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] determined by rhodum (depends on optCON in calling routine)
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state

  USE msingledouble
  USE mphsolvers
  USE msw_ptmp

  IMPLICIT NONE

! Input variables
  !> <b>in situ temperature</b> [degrees C]
  REAL(kind=r8), INTENT(in) :: temp
  !> <b>salinity</b> [on the practical salinity scale, dimensionless]
  REAL(kind=r8), INTENT(in) :: salt
  !> total alkalinity in <b>[eq/m^3]</b> OR in <b>[eq/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: ta
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: tc
  !> phosphate concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: pt
  !> total dissolved inorganic silicon concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  REAL(kind=r8), INTENT(in) :: sit
  !> total boron from either Uppstrom (1974) or Lee et al. (2010), depending on optB in calling routine
  REAL(kind=r8), INTENT(in) :: Bt
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), INTENT(in) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), INTENT(in) :: Ft
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), INTENT(in) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(in) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(in) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=r8), INTENT(in) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), INTENT(in) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! from Dickson and Riley (1979) or Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), INTENT(in) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), INTENT(in) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), INTENT(in) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), INTENT(in) :: Ksi
  !> total atmospheric pressure <b>[atm]</b>
  REAL(kind=r8), INTENT(in) :: Patm
  !> total hydrostatic pressure <b>[bar]</b>
  REAL(kind=r8), INTENT(in) :: Phydro_bar
  !> density factor as computed incalling routine  (vars)
  REAL(kind=r8), INTENT(in) :: rhodum
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(kind=r8), INTENT(out) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=r8), INTENT(out) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=r8), INTENT(out) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=r8), INTENT(out) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=r8), INTENT(out) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=r8), INTENT(out) :: OmegaC

! Local variables
  REAL(kind=r8) :: Phydro_atm, Ptot
  REAL(kind=r8) :: Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: tk, tk0
  real(kind=r8) :: temp68, tempot, tempot68
  REAL(kind=r8) :: Hinit, H
  REAL(kind=r8) :: HSO4, HF, HSI, HPO4
  REAL(kind=r8) :: ab, aw, ac
  REAL(kind=r8) :: cu, cb, cc
  REAL(kind=r8) :: Ca
! Array to pass optional arguments
  CHARACTER(7) :: opGAS

  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

! Compute pH from constants and total concentrations
! - use SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) modified to use mocsy's Ks instead of its own
! 1) Solve for H+ using above result as the initial H+ value
  H = solve_at_general(ta, tc, Bt,                                         & 
                       pt,     sit,                                        &
                       St, Ft,                                             &
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi )
  
! 2) Calculate pH from from H+ concentration (mol/kg)
  IF (H > 0.d0) THEN
     pH = -1.*LOG10(H)
  ELSE
     pH = 1.e20_r8
  ENDIF

! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0d0 + Ks/(H/(1.0d0 + St/Ks)))
  HF = 1.0d0/(1.0d0 + Kf/H)
  HSI = 1.0d0/(1.0d0 + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3) /                &
  (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0d0 + H/Kb)
  aw = Kw/H - H/(1.0d0 + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H

! When optCON = 'mol/m3' in calling routine (vars), then:
! convert output var concentrations from mol/kg to mol/m^3
! e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
  co2  = cu * rhodum
  hco3 = cb * rhodum
  co3  = cc * rhodum

! Determine CO2 fugacity [uatm]
! NOTE: equation just below requires CO2* in mol/kg
  fCO2 = cu * 1.e6_r8/K0

! Determine CO2 partial pressure from CO2 fugacity [uatm]
  tk = 273.15d0 + temp
  !Compute EITHER "potential pCO2" OR "in situ pCO2" (T and P used for calculations will differ)
  IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
     tk0 = tk                 !in situ temperature (K) for K0 calculation
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Ppot' .OR. trim(opGAS) == 'ppot') THEN
     !Use potential temperature and atmospheric pressure (water parcel adiabatically brought back to surface)
     !temp68 = (temp - 0.0002d0) / 0.99975d0          !temp = in situ T; temp68 is same converted to ITPS-68 scale
     !tempot68 = sw_ptmp(salt, temp68, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     !tempot   = 0.99975*tempot68 + 0.0002
     !tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     tempot = sw_ptmp(salt, temp, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
     !Use in situ temperature and total pressure 
     tk0 = tk                             !in situ temperature (K) for fugacity coefficient calculation
     Phydro_atm = Phydro_bar / 1.01325d0  !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
     Ptot = Patm + Phydro_atm            !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
  ELSE
     PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
     STOP
  ENDIF

! Now that we have T and P in the right form, continue with calculation of fugacity coefficient (and pCO2)
  Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
! To compute fugcoeff, we need 3 other terms (B, Del, xc2) in addition to 3 others above (tk, Ptot, Rgas_atm)
  B = -1636.75d0 + 12.0408d0*tk0 - 0.0327957d0*(tk0*tk0) + 0.0000316528d0*(tk0*tk0*tk0)
  Del = 57.7d0 - 0.118d0*tk0
! "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
! x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
! Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
  xCO2approx = fCO2 * 1.e-6_r8
  IF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
!    xCO2approx = 400.0e-6_r8      !a simple test (gives about same result as seacarb for pCO2insitu)
!    approximate surface xCO2 ~ surface fCO2 (i.e., in situ fCO2 d by exponential pressure correction)
     xCO2approx = xCO2approx * exp( ((1-Ptot)*32.3_r8)/(82.05736_r8*tk0) )   ! of K0 press. correction, see Weiss (1974, equation 5)
  ENDIF
  xc2 = (1.0d0 - xCO2approx)**2 
  fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk0) )
  pCO2 = fCO2 / fugcoeff

! Determine Omega Calcite et Aragonite
! OmegaA = ((0.01028d0*salt/35.0d0)*cc) / Kspa
! OmegaC = ((0.01028d0*salt/35.0d0)*cc) / Kspc
! - see comments from Munhoven on the best value "0.02128" which differs slightly from the best practices guide (0.02127)
  Ca = (0.02128d0/40.078d0) * salt/1.80655d0
  OmegaA = (Ca*cc) / Kspa
  OmegaC = (Ca*cc) / Kspc

  RETURN
END SUBROUTINE varsolver

!>    Solve for pH and other carbonate system variables (with input from vars routine)
!>    and compute partial derivatives (buffer factors) using dual numbers (DNAD) technique
SUBROUTINE varsolver_DNAD (ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC,       &
                    temp, salt, ta, tc, pt, sit,                                 &
                    Bt, St, Ft,                                                  &
                    K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi,  & 
                    Patm, Phydro_bar, rhodum, optGAS                             )

  !   Purpose: Solve for pH and other carbonate system variables (with input from vars routine)
  !   and compute partial derivatives (buffer factors) using dual numbers (DNAD) technique
  !
  ! It is similar to subroutine varsolver above except that it computes
  ! partial derivatives of all output (ph, pco2, fco2, co2, hco3, co3, OmegaA and OmegaC)
  ! with respect to four input variables : ta, tc pt and sit
  !     

  !     INPUT variables:
  !     ================
  !     temp    = in situ temperature [degrees C]
  !     ta      = total alkalinity                     in [eq/m^3] or [eq/kg]   based on optCON in calling routine (vars)
  !     tc      = dissolved inorganic carbon           in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     pt      = total dissolved inorganic phosphorus in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     sit     = total dissolved inorganic silicon    in [mol/m^3] or [mol/kg] based on optCON in calling routine (vars)
  !     Bt      = total dissolved inorganic boron      computed in calling routine (vars)
  !     St      = total dissolved inorganic sulfur     computed in calling routine (vars)
  !     Ft      = total dissolved inorganic fluorine   computed in calling routine (vars)
  !     K's     = K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi 
  !     Patm    = atmospheric pressure [atm]
  !     Phydro_bar = hydrostatic pressure [bar]
  !     rhodum  = density factor as computed in calling routine  (vars)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2 (default optGAS = 'Pinsitu')
  !     ---------
  !       PRESSURE & T corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     ph   = pH on total scale
  !     pco2 = CO2 partial pressure (uatm)
  !     fco2 = CO2 fugacity (uatm)
  !     co2  = aqueous CO2 concentration in [mol/kg] or [mol/m^3] determined by rhodum (depends on optCON in calling routine)
  !     hco3 = bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     co3  = carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] determined by rhodum
  !     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
  !     OmegaC = Omega for calcite, i.e., the   calcite saturation state

  USE msingledouble
  USE mphsolvers
  USE msw_ptmp
  USE Dual_Num_Auto_Diff

  IMPLICIT NONE

! Input variables
  !> <b>in situ temperature</b> [degrees C]
  TYPE(DUAL_NUM), INTENT(in) :: temp
  !> <b>salinity</b> [on the practical salinity scale, dimensionless]
  TYPE(DUAL_NUM), INTENT(in) :: salt
  !> total alkalinity in <b>[eq/m^3]</b> OR in <b>[eq/kg]</b>, depending on optCON in calling routine
  TYPE(DUAL_NUM), INTENT(in) :: ta
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  TYPE(DUAL_NUM), INTENT(in) :: tc
  !> phosphate concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  TYPE(DUAL_NUM), INTENT(in) :: pt
  !> total dissolved inorganic silicon concentration in <b>[mol/m^3]</b> OR in <b>[mol/kg]</b>, depending on optCON in calling routine
  TYPE(DUAL_NUM), INTENT(in) :: sit
  !> total boron from either Uppstrom (1974) or Lee et al. (2010), depending on optB in calling routine
  TYPE(DUAL_NUM), INTENT(in) :: Bt
  !> total sulfate (Morris & Riley, 1966)
  TYPE(DUAL_NUM), INTENT(in) :: St
  !> total fluoride  (Riley, 1965)
  TYPE(DUAL_NUM), INTENT(in) :: Ft
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  TYPE(DUAL_NUM), INTENT(in) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  TYPE(DUAL_NUM), INTENT(in) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  TYPE(DUAL_NUM), INTENT(in) :: K2
  !> equilibrium constant for dissociation of boric acid 
  TYPE(DUAL_NUM), INTENT(in) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(in) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  TYPE(DUAL_NUM), INTENT(in) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! from Dickson and Riley (1979) or Perez and Fraga (1987), depending on optKf
  TYPE(DUAL_NUM), INTENT(in) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  TYPE(DUAL_NUM), INTENT(in) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  TYPE(DUAL_NUM), INTENT(in) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(in) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(in) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(in) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  TYPE(DUAL_NUM), INTENT(in) :: Ksi
  !> total atmospheric pressure <b>[atm]</b>
  TYPE(DUAL_NUM), INTENT(in) :: Patm
  !> total hydrostatic pressure <b>[bar]</b>
  TYPE(DUAL_NUM), INTENT(in) :: Phydro_bar
  !> density factor as computed incalling routine  (vars)
  TYPE(DUAL_NUM), INTENT(in) :: rhodum
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> pH on the <b>total scale</b>
  TYPE(DUAL_NUM), INTENT(out) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  TYPE(DUAL_NUM), INTENT(out) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  TYPE(DUAL_NUM), INTENT(out) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  TYPE(DUAL_NUM), INTENT(out) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  TYPE(DUAL_NUM), INTENT(out) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  TYPE(DUAL_NUM), INTENT(out) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  TYPE(DUAL_NUM), INTENT(out) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  TYPE(DUAL_NUM), INTENT(out) :: OmegaC

! Local variables
  TYPE(DUAL_NUM) :: Phydro_atm, Ptot
  TYPE(DUAL_NUM) :: Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  TYPE(DUAL_NUM) :: tk, tk0
  real(kind=r8) :: temp68, tempot, tempot68
  TYPE(DUAL_NUM) :: Hinit, H
  TYPE(DUAL_NUM) :: HSO4, HF, HSI, HPO4
  TYPE(DUAL_NUM) :: ab, aw, ac
  TYPE(DUAL_NUM) :: cu, cb, cc
  TYPE(DUAL_NUM) :: Ca
! Array to pass optional arguments
  CHARACTER(7) :: opGAS

  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

! Compute pH from constants and total concentrations
! - use SolveSAPHE v1.0.1 routines from Munhoven (2013, GMD) 
!    * modified to use mocsy's Ks instead of its own
!    * modified to support DNAD (Dual Numbers Automatic Derivation)
! 1) Solve for H+ using above result as the initial H+ value
  H = solve_at_general_DNAD(ta, tc, Bt,                                    & 
                       pt,     sit,                                        &
                       St, Ft,                                             &
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi )
  
! 2) Calculate pH from from H+ concentration (mol/kg)
  IF (H > 0.d0) THEN
     pH = -1.*LOG10(H)
  ELSE
     pH = 1.e20_r8
  ENDIF

! Compute carbonate Alk (Ac) by difference: from total Alk and other Alk components
  HSO4 = St/(1.0d0 + Ks/(H/(1.0d0 + St/Ks)))
  HF = 1.0d0/(1.0d0 + Kf/H)
  HSI = 1.0d0/(1.0d0 + H/Ksi)
  HPO4 = (K1p*K2p*(H + 2.*K3p) - H**3) /                &
  (H**3 + K1p*H**2 + K1p*K2p*H + K1p*K2p*K3p)
  ab = Bt/(1.0d0 + H/Kb)
  aw = Kw/H - H/(1.0d0 + St/Ks)
  ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4

! Calculate CO2*, HCO3-, & CO32- (in mol/kg soln) from Ct, Ac, H+, K1, & K2
  cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / H)
  cb = K1 * cu / H
  cc = K2 * cb / H

! When optCON = 'mol/m3' in calling routine (vars), then:
! convert output var concentrations from mol/kg to mol/m^3
! e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
  co2  = cu * rhodum
  hco3 = cb * rhodum
  co3  = cc * rhodum

! Determine CO2 fugacity [uatm]
! NOTE: equation just below requires CO2* in mol/kg
  fCO2 = cu * 1.e6_r8/K0

! Determine CO2 partial pressure from CO2 fugacity [uatm]
  tk = 273.15d0 + temp
  !Compute EITHER "potential pCO2" OR "in situ pCO2" (T and P used for calculations will differ)
  IF     (trim(opGAS) == 'Pzero'   .OR. trim(opGAS) == 'pzero') THEN
     tk0 = tk                 !in situ temperature (K) for K0 calculation
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Ppot' .OR. trim(opGAS) == 'ppot') THEN
     !Use potential temperature and atmospheric pressure (water parcel adiabatically brought back to surface)
     !temp68 = (temp - 0.0002d0) / 0.99975d0          !temp = in situ T; temp68 is same converted to ITPS-68 scale
     !tempot68 = sw_ptmp(salt, temp68, Phydro_bar*10d0, 0.0d0) !potential temperature (C)
     !tempot   = 0.99975*tempot68 + 0.0002
     !tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     tempot = sw_ptmp(salt%x_ad_, temp%x_ad_, Phydro_bar%x_ad_*10d0, 0.0d0) !potential temperature (C)
     tk0 = tempot + 273.15d0  !potential temperature (K) for fugacity coeff. calc as needed for potential fCO2 & pCO2
     Ptot = Patm              !total pressure (in atm) = atmospheric pressure ONLY
  ELSEIF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
     !Use in situ temperature and total pressure 
     tk0 = tk                             !in situ temperature (K) for fugacity coefficient calculation
     Phydro_atm = Phydro_bar / 1.01325d0  !convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
     Ptot = Patm + Phydro_atm            !total pressure (in atm) = atmospheric pressure + hydrostatic pressure
  ELSE
     PRINT *, "optGAS must be 'Pzero', 'Ppot', or 'Pinsitu'"
     STOP
  ENDIF

! Now that we have T and P in the right form, continue with calculation of fugacity coefficient (and pCO2)
  Rgas_atm = 82.05736_r8      ! (cm3 * atm) / (mol * K)  CODATA (2006)
! To compute fugcoeff, we need 3 other terms (B, Del, xc2) in addition to 3 others above (tk, Ptot, Rgas_atm)
  B = -1636.75d0 + 12.0408d0*tk0 - 0.0327957d0*(tk0*tk0) + 0.0000316528d0*(tk0*tk0*tk0)
  Del = 57.7d0 - 0.118d0*tk0
! "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
! x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
! Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
  xCO2approx = fCO2 * 1.e-6_r8
  IF (trim(opGAS) == 'Pinsitu' .OR. trim(opGAS) == 'pinsitu') THEN
!    xCO2approx = 400.0e-6_r8      !a simple test (gives about same result as seacarb for pCO2insitu)
!    approximate surface xCO2 ~ surface fCO2 (i.e., in situ fCO2 d by exponential pressure correction)
     xCO2approx = xCO2approx * exp( ((1-Ptot)*32.3_r8)/(82.05736_r8*tk0) )   ! of K0 press. correction, see Weiss (1974, equation 5)
  ENDIF
  xc2 = (1.0d0 - xCO2approx)**2 
  fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk0) )
  pCO2 = fCO2 / fugcoeff

! Determine Omega Calcite et Aragonite
! OmegaA = ((0.01028d0*salt/35.0d0)*cc) / Kspa
! OmegaC = ((0.01028d0*salt/35.0d0)*cc) / Kspc
! - see comments from Munhoven on the best value "0.02128" which differs slightly from the best practices guide (0.02127)
  Ca = (0.02128d0/40.078d0) * salt/1.80655d0
  OmegaA = (Ca*cc) / Kspa
  OmegaC = (Ca*cc) / Kspc

  RETURN
END SUBROUTINE varsolver_DNAD
END MODULE mvarsolver
