!> \file errors.f90
!! \BRIEF 
!> Module with errors subroutine - Propagate standard error (or uncertainty) in carbonate system vars 
MODULE merrors
CONTAINS


!>    Compute standard errors (or uncertainties) of standard carbonate system variables  
!>    (H+, pCO2, fCO2, CO2*, HCO3- and CO3--, OmegaA, OmegaC) based on uncertainties
!!    in input variables, which are Alk, DIC, Total Silicate and Phosphate, T and S
!!    and in dissociations constants  K0, K1, K2, Kb, Kw, Kspa, Kspc
SUBROUTINE errors  (eH, epCO2, efCO2, eCO2, eHCO3, eCO3, eOmegaA, eOmegaC,    &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,          &
                etemp, esal, ealk, edic, esil, ephos,                         &
                optCON, optT, optP, optB, optK1K2, optKf, optGAS, optS, lon,  & 
                r, epK, ebt )
         
!     This subroutine does error propagation on the computation of carbonate system variables 
!     from errors (or uncertainties) on six input 
!      - Alk and DIC
!      - nutrients (silicate and phosphate concentrations)
!      - temperature and salinity
!     plus errors on dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
!    
!     It computes numerical derivatives then applies error propagation using the method of moments.
!     The method of moments is a very general technique for estimating the second moment of a variable z
!     (variance or standard deviation) based on a first-order approximation to z.
!    
!     INPUT variables:
!     ================
!     - temp, sal, alk, dic, sil, phos, Patm, depth, lat :  same as input of subroutine  vars() : vectors
!    
!     - N              :  number of data points (= size of input and output vectors)
!     - ealk, edic     :  standard error (or uncertainty) on alk and DIC
!     - esal, etemp    :  standard error (or uncertainty) on Salinity and Temperature
!     - ephos, esil    :  standard error (or uncertainty) on Phosphate and Silicate total concentrations
!     - epK            :  standard error (or uncertainty) on all seven dissociation constants (a vector) [pK units]
!     - ebt            :  standard error (or uncertainty) on total boron fractional error [default=0.01], i.e., 1%
!    
!       All parameters are vectors of length N except epK 
!         epK must be vector of seven values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
!         these errors are assumed to be equal for all data points.
!         In constrast, ealk, edic, esal, etemp, ephos and esilt are errors associated with each data point.
!
!       epK is optional :  if not given, following default error values will be taken :
!                   pK0     pK1    pK2,   pKb    pKw    pKspa   pKspc
!                   0.004   0.015  0.03   0.01   0.01   0.02    0.02 
!
!       ebt is optional :  if not given, the default error is 0.01 (i.e., 1%)
!                   Bt (rel. err)
!                   0.01
!    
!     INPUT options:
!     ==============
!    
!     - optCON, optT, optP, opB, optK1K2, optKf, optGAS, optS, lon :  same as for input of subroutine vars() 
!    
!    
!     OUTPUT variables:
!     =================
!
!     - eH       total error to [H+] concentration  (mol/kg) on the <b>total scale</b>
!     - epCO2    total error to "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (µatm)
!     - efCO2    total error to "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (µatm)
!     - eCO2     total error to CO2 concentration (mol/kg)
!     - eHCO3    total error to HCO3 concentration (mol/kg)
!     - eCO3     total error to CO3 concentration (mol/kg)
!     - eOmegaA  total error of Omega aragonite (aragonite saturation state)
!     - eOmegaC  total error of Omega calcite   (calcite saturation state)
!
!*****************************************************************************************************

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  USE mconstants
  USE mderivnum

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

  !> standard error (or uncertainty) on Temperature
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: etemp
  !> standard error (or uncertainty) on Salinity
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: esal
  !> standard error (or uncertainty) on alk
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: ealk
  !> standard error (or uncertainty) on DIC
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: edic
  !> standard error (or uncertainty) on Silicate total concentrations
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: esil
  !> standard error (or uncertainty) on Phosphate total concentrations
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: ephos
  !> standard error (or uncertainty) on all seven dissociation constants (a vector)
!f2py real(8) intent(in), optional, dimension(7) :: epK=(0.004,0.015,0.03,0.01,0.01,0.02,0.02)
  REAL(kind=rx), INTENT(in), OPTIONAL, DIMENSION(7) :: epK
  !> correlation coefficient (-1 < r < 1) for correlation between ALK and DIC (zero by default)
!f2py  real intent(in), optional :: r = 0.0
  REAL(kind=rx), OPTIONAL, INTENT(in) :: r
  !> standard error (or uncertainty) on total boron (Bt) - a single number, not a vector like other errors
!f2py  real intent(in), optional :: ebt=0.01
  REAL(kind=rx), OPTIONAL, INTENT(in) :: ebt 

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
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
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
  REAL(kind=rx), OPTIONAL, INTENT(in),    DIMENSION(N) :: lon

! Output variables:
  !> total error to [H+] concentration  (mol/kg) on the <b>total scale</b>
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eH
  !> total error to "standard" pCO2, CO2 partial pressure computed at in situ temperature and atmospheric pressure (µatm)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: epCO2
  !> total error to "standard" fCO2, CO2 fugacity computed at in situ temperature and atmospheric pressure (µatm)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: efCO2
  !> total error to CO2 concentration (mol/kg)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eCO2
  !> total error to HCO3 concentration (mol/kg)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eHCO3
  !> total error to CO3 concentration (mol/kg)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eCO3
  !> total error of Omega aragonite (aragonite saturation state)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eOmegaA
  !> total error of Omega calcite   (calcite saturation state)
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: eOmegaC

! Local variables

  ! Default value for errors on pK
  REAL(kind=rx), DIMENSION(7) :: epK_local = (/0.004_r8, 0.015_r8, 0.03_r8, 0.01_r8, 0.01_r8, 0.02_r8, 0.02_r8/)
  REAL(kind=rx), DIMENSION(7) :: epKstd = (/0.004_r8, 0.015_r8, 0.03_r8, 0.01_r8, 0.01_r8, 0.02_r8, 0.02_r8/)
  REAL(kind=rx), DIMENSION(7) :: epKzero
! REAL(kind=rx), DIMENSION(7) :: epKstd
! Extend epK_local by 1 to later include error for Bt (simplifies coding)
  REAL(kind=rx), DIMENSION(8) :: epK_local8
  CHARACTER*3, DIMENSION(8) ::  Kid = (/'k0 ','k1 ','k2 ','kb ','kw ','ka ','kc ', 'bt '/)
! CHARACTER*3, DIMENSION(7) :: Kid = (/'k0 ','k1 ','k2 ','kb ','kw ','ka ','kc '/)

  ! Default value for error on Total Boron (ebt)
  REAL(kind=rx) :: ebt_local = 0.01_rx
! REAL(kind=rx) :: ebt_local 

  ! Default value for correlation between ALK & DIC
  REAL(kind=rx) :: r_local = 0.0_r8

  ! derivative of H on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(N) :: dh_dx
  ! derivative of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(N) :: dpco2_dx
  ! derivative of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(N) :: dfco2_dx
  ! derivative of aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: dco2_dx
  ! derivative of (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: dhco3_dx
  ! derivative of (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: dco3_dx
  ! derivative of Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), DIMENSION(N) :: dOmegaA_dx
  ! derivative of Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), DIMENSION(N) :: dOmegaC_dx

  ! covariance tmp array for derivative of H on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(N) :: r_dh_dx
  ! covariance tmp array for derivative of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(N) :: r_dpco2_dx
  ! covariance tmp array for derivative of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(N) :: r_dfco2_dx
  ! covariance tmp array for derivative of aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: r_dco2_dx
  ! covariance tmp array for derivative of (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: r_dhco3_dx
  ! covariance tmp array for derivative of (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(N) :: r_dco3_dx
  ! covariance tmp array for derivative of Omega for aragonite, i.e., the aragonite saturation state 
  REAL(kind=rx), DIMENSION(N) :: r_dOmegaA_dx
  ! covariance tmp array for Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), DIMENSION(N) :: r_dOmegaC_dx

  REAL(kind=r8), DIMENSION(N) :: K0, K1, K2, Kb, Kw, Ks, Kf, Kspc
  REAL(kind=r8), DIMENSION(N) :: Kspa, K1p, K2p, K3p, Ksi
  REAL(kind=r8), DIMENSION(N) :: St, Ft, Bt
  ! Standard error on dissociation constant Kx
  REAL(kind=r8), DIMENSION(N) :: eK

! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(3) :: opK1K2
  CHARACTER(2) :: opKf
  CHARACTER(7) :: opGAS
  CHARACTER(4) :: opS

  INTEGER :: i

! *******************************************************************

! Set defaults for optional arguments (in Fortran90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion

  epKstd = (/0.004d0, 0.015d0, 0.03d0, 0.01d0, 0.01d0, 0.02d0, 0.02d0/)
  epK_local(:) = epKstd(:)
  epKzero = epKstd * 0.0

  !------------------------------------------------------------------------------------------------------
  ! IMPORTANT: in python (after f2py), the result of "PRESENT()" clause below is always .TRUE.
  !            Thus the corresponding ELSE blocks below are executed only when in pure FORTRAN;
  !            They are not reached in a python call of the same fortran routines. Thus
  !            for f2py, the same defaults are assigned with "!f2py" statements above
  IF (PRESENT(optB)) THEN
     opB = optB
  ELSE
!   Default is Uppstrom (1974) for total boron
     opB = 'u74'
  ENDIF
  IF (PRESENT(optK1K2)) THEN
     opK1K2 = optK1K2
  ELSE
!    Default is Lueker et al. 2000) for K1 & K2
     opK1K2 = 'l'
  ENDIF
  IF (PRESENT(optKf)) THEN
     opKf = optKf
  ELSE
!    Default is Perez & Fraga (1987) for Kf
     opKf = 'pf'
  ENDIF
  IF (PRESENT(optGAS)) THEN
     opGAS = optGAS
  ELSE
     opGAS = 'Pinsitu'
  ENDIF
  IF (.NOT. PRESENT(optS)) THEN
     opS = optS
  ELSE
     opS = 'Sprc'
  ENDIF
! IF (PRESENT(lon)) THEN
!    lon_local = lon
! ELSE
!    lon_local = 1.e20
! ENDIF
  IF (PRESENT(r)) THEN
     r_local = r
  ELSE
     r_local = 0.0_r8
  ENDIF
  IF (PRESENT(epK)) THEN
     epK_local(:) = epK(:)
     ! Kludge: when used in python, the init. of epK fails (does not work with an array)
     ! i.e;, all 7 epK members are equal to last specified element on !f2py statement)
     IF (epK(1)==epK(7) .and. epK(2)==epK(7) .and. epK(3)==epK(7) .AND. &
         epK(4)==epK(7) .and. epK(5)==epK(7) .and. epK(6)==epK(7)         ) THEN
         ! if so, initialize to the default
         IF (epK(1)==0.0) THEN
            epK_local = epKzero
         ELSE
            epK_local = epKstd
         ENDIF
     ENDIF
  ELSE
     epK_local = (/0.004_r8, 0.015_r8, 0.03_r8, 0.01_r8, 0.01_r8, 0.02_r8, 0.02_r8/)
  ENDIF
  IF (PRESENT(ebt)) THEN
     ebt_local = ebt
  ELSE
     ebt_local = 0.01_rx
  ENDIF
  
  ! initialise total square error
  eH(:) = 0
  epCO2(:) = 0
  efCO2(:) = 0
  eCO2(:) = 0
  eHCO3(:) = 0
  eCO3(:) = 0
  eOmegaA(:) = 0
  eOmegaC(:) = 0
    
  ! Contribution of Alk to squared standard error
  ! ---------------------------------------------
  
  ! Compute sensitivities (partial derivatives)
  CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
            dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
            temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'alk',     &
            optCON, optT, optP, opB, opK1K2, opKf, opGAS                    )

  ! Covariance (only when R is not zero)
  IF (r_local .NE. 0.0_r8) THEN
     r_dh_dx(:)      = dh_dx(:)      * ealk(:)
     r_dpco2_dx(:)   = dpco2_dx(:)   * ealk(:)
     r_dfco2_dx(:)   = dfco2_dx(:)   * ealk(:)
     r_dco2_dx(:)    = dco2_dx(:)    * ealk(:)
     r_dhco3_dx(:)   = dhco3_dx(:)   * ealk(:)
     r_dco3_dx(:)    = dco3_dx(:)    * ealk(:)
     r_dOmegaA_dx(:) = dOmegaA_dx(:) * ealk(:)
     r_dOmegaC_dx(:) = dOmegaC_dx(:) * ealk(:)
  ENDIF

  ! multiply derivatives by error
  dh_dx(:)      = dh_dx(:)      * ealk(:)
  dpco2_dx(:)   = dpco2_dx(:)   * ealk(:)
  dfco2_dx(:)   = dfco2_dx(:)   * ealk(:)
  dco2_dx(:)    = dco2_dx(:)    * ealk(:)
  dhco3_dx(:)   = dhco3_dx(:)   * ealk(:)
  dco3_dx(:)    = dco3_dx(:)    * ealk(:)
  dOmegaA_dx(:) = dOmegaA_dx(:) * ealk(:)
  dOmegaC_dx(:) = dOmegaC_dx(:) * ealk(:)
  ! Square 
  dh_dx(:)      = dh_dx(:) * dh_dx(:)
  dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
  dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
  dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
  dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
  dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
  dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
  dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
  ! sum up
  eH(:)    = eH(:)    + dh_dx(:)
  epCO2(:) = epCO2(:) + dpco2_dx(:)
  efCO2(:) = efCO2(:) + dfco2_dx(:)
  eCO2(:)  = eCO2(:)  + dco2_dx(:)
  eHCO3(:) = eHCO3(:) + dhco3_dx(:)
  eCO3(:)  = eCO3(:)  + dco3_dx(:)
  eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
  eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)

  ! Contribution of DIC to squared standard error
  ! ---------------------------------------------

  ! Compute sensitivities (partial derivatives)
  CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
            dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
            temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'dic',     &
            optCON, optT, optP, opB, opK1K2, opKf, opGAS                    )
  
  ! Covariance (only when R is not zero) = 2*r*(dvar/dAt)*(dvar/dCt)*sigma1*sigma2
  IF (r_local .NE. 0.0_r8) THEN
     eH(:)      = eH(:)      + 2.0_r8 * r_local * r_dh_dx(:)      * dh_dx(:)      * edic(:)
     epco2(:)   = epco2(:)   + 2.0_r8 * r_local * r_dpco2_dx(:)   * dpco2_dx(:)   * edic(:)
     efco2(:)   = efco2(:)   + 2.0_r8 * r_local * r_dfco2_dx(:)   * dfco2_dx(:)   * edic(:)
     eco2(:)    = eco2(:)    + 2.0_r8 * r_local * r_dco2_dx(:)    * dco2_dx(:)    * edic(:)
     ehco3(:)   = ehco3(:)   + 2.0_r8 * r_local * r_dhco3_dx(:)   * dhco3_dx(:)   * edic(:)
     eco3(:)    = eco3(:)    + 2.0_r8 * r_local * r_dco3_dx(:)    * dco3_dx(:)    * edic(:)
     eOmegaA(:) = eOmegaA(:) + 2.0_r8 * r_local * r_dOmegaA_dx(:) * dOmegaA_dx(:) * edic(:)
     eOmegaC(:) = eOmegaC(:) + 2.0_r8 * r_local * r_dOmegaC_dx(:) * dOmegaC_dx(:) * edic(:)
  ENDIF

  ! multiply derivatives by error
  dh_dx(:)      = dh_dx(:)      * edic(:)
  dpco2_dx(:)   = dpco2_dx(:)   * edic(:)
  dfco2_dx(:)   = dfco2_dx(:)   * edic(:)
  dco2_dx(:)    = dco2_dx(:)    * edic(:)
  dhco3_dx(:)   = dhco3_dx(:)   * edic(:)
  dco3_dx(:)    = dco3_dx(:)    * edic(:)
  dOmegaA_dx(:) = dOmegaA_dx(:) * edic(:)
  dOmegaC_dx(:) = dOmegaC_dx(:) * edic(:)
  ! Square 
  dh_dx(:)      = dh_dx(:) * dh_dx(:)
  dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
  dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
  dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
  dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
  dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
  dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
  dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
  ! sum up
  eH(:)    = eH(:)    + dh_dx(:)
  epCO2(:) = epCO2(:) + dpco2_dx(:)
  efCO2(:) = efCO2(:) + dfco2_dx(:)
  eCO2(:)  = eCO2(:)  + dco2_dx(:)
  eHCO3(:) = eHCO3(:) + dhco3_dx(:)
  eCO3(:)  = eCO3(:)  + dco3_dx(:)
  eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
  eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)

  ! Contribution of Silicon (total dissolved inorganic concentration) to squared standard error
  ! -------------------------------------------------------------------------------------------

  IF (any (esil .ne. 0.0)) THEN
      ! Compute sensitivities (partial derivatives)
      CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
                dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'sil',     &
                optCON, optT, optP, opB, opK1K2, opKf, opGAS                )

      ! Cancels derivative where Sil = 0 
      !   because computation of derivative w/ respect to sil fails in that case
      WHERE (sil(:) .EQ. 0)
          dh_dx(:)      = 0
          dpco2_dx(:)   = 0
          dfco2_dx(:)   = 0
          dco2_dx(:)    = 0
          dhco3_dx(:)   = 0
          dco3_dx(:)    = 0
          dOmegaA_dx(:) = 0
          dOmegaC_dx(:) = 0
      END WHERE
      
      ! multiply derivatives by error
      dh_dx(:)      = dh_dx(:)      * esil(:)
      dpco2_dx(:)   = dpco2_dx(:)   * esil(:)
      dfco2_dx(:)   = dfco2_dx(:)   * esil(:)
      dco2_dx(:)    = dco2_dx(:)    * esil(:)
      dhco3_dx(:)   = dhco3_dx(:)   * esil(:)
      dco3_dx(:)    = dco3_dx(:)    * esil(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * esil(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * esil(:)
      ! Square 
      dh_dx(:)      = dh_dx(:) * dh_dx(:)
      dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
      dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
      dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
      dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
      dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
      ! sum up
      eH(:)    = eH(:)    + dh_dx(:)
      epCO2(:) = epCO2(:) + dpco2_dx(:)
      efCO2(:) = efCO2(:) + dfco2_dx(:)
      eCO2(:)  = eCO2(:)  + dco2_dx(:)
      eHCO3(:) = eHCO3(:) + dhco3_dx(:)
      eCO3(:)  = eCO3(:)  + dco3_dx(:)
      eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
      eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)
  ENDIF

  ! Contribution of Phosphorus (total dissoloved inorganic concentration) to squared standard error
  ! -----------------------------------------------------------------------------------------------

  IF (any (ephos .ne. 0.0)) THEN

      ! Compute sensitivities (partial derivatives)
      CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
                dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'pho',     &
                optCON, optT, optP, opB, opK1K2, opKf, opGAS                )

      ! Cancels derivative where Phos = 0 
      !   because computation of derivative w/ respect to phos fails in that case
      WHERE (phos(:) .EQ. 0)
          dh_dx(:)      = 0
          dpco2_dx(:)   = 0
          dfco2_dx(:)   = 0
          dco2_dx(:)    = 0
          dhco3_dx(:)   = 0
          dco3_dx(:)    = 0
          dOmegaA_dx(:) = 0
          dOmegaC_dx(:) = 0
      END WHERE
      
      ! multiply derivatives by error
      dh_dx(:)      = dh_dx(:)      * ephos(:)
      dpco2_dx(:)   = dpco2_dx(:)   * ephos(:)
      dfco2_dx(:)   = dfco2_dx(:)   * ephos(:)
      dco2_dx(:)    = dco2_dx(:)    * ephos(:)
      dhco3_dx(:)   = dhco3_dx(:)   * ephos(:)
      dco3_dx(:)    = dco3_dx(:)    * ephos(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * ephos(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * ephos(:)
      ! Square 
      dh_dx(:)      = dh_dx(:) * dh_dx(:)
      dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
      dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
      dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
      dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
      dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
      ! sum up
      eH(:)    = eH(:)    + dh_dx(:)
      epCO2(:) = epCO2(:) + dpco2_dx(:)
      efCO2(:) = efCO2(:) + dfco2_dx(:)
      eCO2(:)  = eCO2(:)  + dco2_dx(:)
      eHCO3(:) = eHCO3(:) + dhco3_dx(:)
      eCO3(:)  = eCO3(:)  + dco3_dx(:)
      eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
      eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)
  ENDIF

  ! Contribution of temperature to squared standard error
  ! -----------------------------------------------------
  
  IF (any (etemp .ne. 0.0)) THEN
      ! Compute sensitivities (partial derivatives)
      CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
                dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'tem',     &
                optCON, optT, optP, opB, opK1K2, opKf, opGAS                )

      ! multiply derivatives by error
      dh_dx(:)      = dh_dx(:)      * etemp(:)
      dpco2_dx(:)   = dpco2_dx(:)   * etemp(:)
      dfco2_dx(:)   = dfco2_dx(:)   * etemp(:)
      dco2_dx(:)    = dco2_dx(:)    * etemp(:)
      dhco3_dx(:)   = dhco3_dx(:)   * etemp(:)
      dco3_dx(:)    = dco3_dx(:)    * etemp(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * etemp(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * etemp(:)
      ! Square 
      dh_dx(:)      = dh_dx(:) * dh_dx(:)
      dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
      dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
      dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
      dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
      dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
      ! sum up
      eH(:)    = eH(:)    + dh_dx(:)
      epCO2(:) = epCO2(:) + dpco2_dx(:)
      efCO2(:) = efCO2(:) + dfco2_dx(:)
      eCO2(:)  = eCO2(:)  + dco2_dx(:)
      eHCO3(:) = eHCO3(:) + dhco3_dx(:)
      eCO3(:)  = eCO3(:)  + dco3_dx(:)
      eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
      eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)
  ENDIF


  ! Contribution of salinity to squared standard error
  ! --------------------------------------------------
  
  IF (any (esal .ne. 0.0)) THEN
      ! Compute sensitivities (partial derivatives)
      CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
                dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, 'sal',     &
                optCON, optT, optP, opB, opK1K2, opKf, opGAS                )

      ! multiply derivatives by error
      dh_dx(:)      = dh_dx(:)      * esal(:)
      dpco2_dx(:)   = dpco2_dx(:)   * esal(:)
      dfco2_dx(:)   = dfco2_dx(:)   * esal(:)
      dco2_dx(:)    = dco2_dx(:)    * esal(:)
      dhco3_dx(:)   = dhco3_dx(:)   * esal(:)
      dco3_dx(:)    = dco3_dx(:)    * esal(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * esal(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * esal(:)
      ! Square 
      dh_dx(:)      = dh_dx(:) * dh_dx(:)
      dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
      dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
      dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
      dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
      dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
      dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
      dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
      ! sum up
      eH(:)    = eH(:)    + dh_dx(:)
      epCO2(:) = epCO2(:) + dpco2_dx(:)
      efCO2(:) = efCO2(:) + dfco2_dx(:)
      eCO2(:)  = eCO2(:)  + dco2_dx(:)
      eHCO3(:) = eHCO3(:) + dhco3_dx(:)
      eCO3(:)  = eCO3(:)  + dco3_dx(:)
      eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
      eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)
  ENDIF

  ! Contribution of the equil constants and Bt to squared standard error
  ! --------------------------------------------------------------------
  
  ! Make vector of (epK_local, ebt_local), i.e., adding 1 value at end of epK (for testing below)
  epK_local8(1:7) = epK_local
  epK_local8(8)   = ebt_local
  
  ! Preliminary calculations for dissociation constants and total boron
  IF (any (epK_local8 .NE. 0) ) THEN
  
      ! Get all equilibrium constants and total concentrations of SO4, F, B
      CALL constants (K0, K1, K2, Kb, Kw, Ks, Kf,                &
                  Kspc, Kspa, K1p, K2p, K3p, Ksi, St, Ft, Bt,    &
                  temp, sal, Patm, depth, lat, N,                &
                  optT, optP, opB, opK1K2, opKf, opGAS        )


      ! Contribution of all pKi to squared standard error
      DO i = 1,8
          ! if error on Ki is given
          IF (epK_local8(i) .ne. 0.0) THEN

              ! compute error on Ki from that on pKi
              SELECT CASE (i)
              CASE(1)
                  eK(:) = - epK_local(i) * K0(:) * log(10.0)
              CASE(2)
                  eK(:) = - epK_local(i) * K1(:) * log(10.0)
              CASE(3)
                  eK(:) = - epK_local(i) * K2(:) * log(10.0)
              CASE(4)
                  eK(:) = - epK_local(i) * Kb(:) * log(10.0)
              CASE(5)
                  eK(:) = - epK_local(i) * Kw(:) * log(10.0)
              CASE(6)
                  eK(:) = - epK_local(i) * Kspa(:) * log(10.0)
              CASE(7)
                  eK(:) = - epK_local(i) * Kspc(:) * log(10.0)
              CASE(8)
                  ! for convenience we tag on ebt here, even though Bt is not an equilibrium constant
                  ! For Bt, we start from % error in Bt (not from absolute delta pK, so no log(10) conversion)
                  eK(:) =  ebt_local    * Bt(:) 
              END SELECT
              
              ! Compute sensitivities (partial derivatives)
              CALL derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,              &
                        dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                        temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, Kid(i),    &
                        optCON, optT, optP, opB, opK1K2, opKf, opGAS                    )

              ! multiply derivatives by error
              dh_dx(:)      = dh_dx(:)      * eK(:)
              dpco2_dx(:)   = dpco2_dx(:)   * eK(:)
              dfco2_dx(:)   = dfco2_dx(:)   * eK(:)
              dco2_dx(:)    = dco2_dx(:)    * eK(:)
              dhco3_dx(:)   = dhco3_dx(:)   * eK(:)
              dco3_dx(:)    = dco3_dx(:)    * eK(:)
              dOmegaA_dx(:) = dOmegaA_dx(:) * eK(:)
              dOmegaC_dx(:) = dOmegaC_dx(:) * eK(:)
              ! Square 
              dh_dx(:)      = dh_dx(:) * dh_dx(:)
              dpco2_dx(:)   = dpco2_dx(:) * dpco2_dx(:) 
              dfco2_dx(:)   = dfco2_dx(:) * dfco2_dx(:)
              dco2_dx(:)    = dco2_dx(:)  * dco2_dx(:)
              dhco3_dx(:)   = dhco3_dx(:) * dhco3_dx(:)
              dco3_dx(:)    = dco3_dx(:)  * dco3_dx(:)
              dOmegaA_dx(:) = dOmegaA_dx(:) * dOmegaA_dx(:)
              dOmegaC_dx(:) = dOmegaC_dx(:) * dOmegaC_dx(:)
              ! sum up
              eH(:)    = eH(:)    + dh_dx(:)
              epCO2(:) = epCO2(:) + dpco2_dx(:)
              efCO2(:) = efCO2(:) + dfco2_dx(:)
              eCO2(:)  = eCO2(:)  + dco2_dx(:)
              eHCO3(:) = eHCO3(:) + dhco3_dx(:)
              eCO3(:)  = eCO3(:)  + dco3_dx(:)
              eOmegaA(:) = eOmegaA(:) + dOmegaA_dx(:)
              eOmegaC(:) = eOmegaC(:) + dOmegaC_dx(:)
          ENDIF
      END DO
  ENDIF

  ! Compute resulting total error (or uncertainty)
  eH       = sqrt(eH)   
  epCO2    = sqrt(epCO2)
  efCO2    = sqrt(efCO2)
  eCO2     = sqrt(eCO2) 
  eHCO3    = sqrt(eHCO3)
  eCO3     = sqrt(eCO3) 
  eOmegaA  = sqrt(eOmegaA)
  eOmegaC  = sqrt(eOmegaC)

  RETURN
END SUBROUTINE errors
END MODULE merrors
