!> \file phsolvers.f90
!! \BRIEF 
!> Module with routines needed to solve pH-total alkalinity equation (Munhoven, 2013, GMD)
MODULE mphsolvers
!   Module of fastest solvers from Munhoven (2013, Geosci. Model Dev., 6, 1367-1388)
!   ! Taken from SolveSAPHE (mod_phsolvers.F90) & adapted very slightly for use with mocsy
!   ! SolveSaphe is distributed under the GNU Lesser General Public License
!   ! mocsy is distributed under the MIT License
!
! Modifications J. C. Orr, LSCE/IPSL, CEA-CNRS-UVSQ, France, 11 Sep 2014:
! 1) kept only the 3 fastest solvers (atgen, atsec, atfast) and routines which they call
! 2) reduced vertical white space: deleted many blank lines & comment lines that served as divisions
! 3) converted name from .F90 to .f90, deleting a few optional preprocesse if statements
! 4) read in mocsy computed equilibrium constants (as arguments) instead of USE MOD_CHEMCONST
! 5) converted routine names from upper case to lower case
! 6) commented out arguments and equations for NH4 and H2S acid systems

USE msingledouble
USE Dual_Num_Auto_Diff
IMPLICIT NONE

! General parameters
REAL(KIND=wp), PARAMETER :: pp_rdel_ah_target = 1.E-8_wp
REAL(KIND=wp), PARAMETER :: pp_ln10 = 2.302585092994045684018_wp

! Maximum number of iterations for each method
INTEGER, PARAMETER :: jp_maxniter_atgen    = 50
INTEGER, PARAMETER :: jp_maxniter_atsec    = 50
INTEGER, PARAMETER :: jp_maxniter_atfast   = 50

! Bookkeeping variables for each method
! - SOLVE_AT_GENERAL and SOLVE_AT_GENERAL_DNAD
INTEGER :: niter_atgen    = jp_maxniter_atgen

! - SOLVE_AT_GENERAL_SEC
INTEGER :: niter_atsec    = jp_maxniter_atsec

! - SOLVE_AT_FAST (variant of SOLVE_AT_GENERAL w/o bracketing
INTEGER :: niter_atfast   = jp_maxniter_atfast

! Keep the following functions private to avoid conflicts with
! other modules that provide similar ones.
!PRIVATE AHINI_FOR_AT

CONTAINS
!===============================================================================
SUBROUTINE anw_infsup(p_dictot, p_bortot,                                     &
                      p_po4tot, p_siltot,                                     &
                      p_so4tot, p_flutot,                                     &
                      p_alknw_inf, p_alknw_sup)

! Subroutine returns the lower and upper bounds of "non-water-selfionization"
! contributions to total alkalinity (the infimum and the supremum), i.e
! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

USE msingledouble
IMPLICIT NONE

! Argument variables
REAL(KIND=wp), INTENT(IN)  :: p_dictot
REAL(KIND=wp), INTENT(IN)  :: p_bortot
REAL(KIND=wp), INTENT(IN)  :: p_po4tot
REAL(KIND=wp), INTENT(IN)  :: p_siltot
!REAL(KIND=wp), INTENT(IN)  :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)  :: p_h2stot
REAL(KIND=wp), INTENT(IN)  :: p_so4tot
REAL(KIND=wp), INTENT(IN)  :: p_flutot
REAL(KIND=wp), INTENT(OUT) :: p_alknw_inf
REAL(KIND=wp), INTENT(OUT) :: p_alknw_sup

p_alknw_inf =  -p_po4tot - p_so4tot - p_flutot
p_alknw_sup =   p_dictot + p_dictot + p_bortot &
              + p_po4tot + p_po4tot + p_siltot !&
!             + p_nh4tot + p_h2stot

RETURN
END SUBROUTINE anw_infsup

!===============================================================================

FUNCTION equation_at(p_alktot, p_h,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                      &
                     p_so4tot, p_flutot,                                      &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                     p_deriveqn)

  ! Purpose: Compute total alkalinity from ion concentrations and equilibrium constants

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: equation_at

! Argument variables
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_h
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_deriveqn

! Local variables 
!-----------------
REAL(KIND=wp) :: znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
REAL(KIND=wp) :: znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
REAL(KIND=wp) :: znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
REAL(KIND=wp) :: znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
REAL(KIND=wp) :: znumer_nh4, zdnumer_nh4, zdenom_nh4, zalk_nh4, zdalk_nh4
REAL(KIND=wp) :: znumer_h2s, zdnumer_h2s, zdenom_h2s, zalk_h2s, zdalk_h2s
REAL(KIND=wp) :: znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
REAL(KIND=wp) :: znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
REAL(KIND=wp) ::                                      zalk_wat, zdalk_wat
REAL(KIND=wp) :: aphscale

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

! H2CO3 - HCO3 - CO3 : n=2, m=0
znumer_dic = 2._wp*K1*K2 + p_h*       K1
zdenom_dic =       K1*K2 + p_h*(      K1 + p_h)
zalk_dic   = p_dictot * (znumer_dic/zdenom_dic)

! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       Kb
zdenom_bor =       Kb + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._wp*K1p*K2p*K3p + p_h*(2._wp*K1p*K2p + p_h* K1p)
zdenom_po4 =       K1p*K2p*K3p + p_h*(      K1p*K2p + p_h*(K1p + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._wp) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       Ksi
zdenom_sil =       Ksi + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! NH4 - NH3 : n=1, m=0
!znumer_nh4 =       api1_nh4
!zdenom_nh4 =       api1_nh4 + p_h
!zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)
! Note: api1_nh4 = Knh4

! H2S - HS : n=1, m=0
!znumer_h2s =       api1_h2s
!zdenom_h2s =       api1_h2s + p_h
!zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)
! Note: api1_h2s = Kh2s

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       Ks
zdenom_so4 =       Ks + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._wp)

! HF - F : n=1, m=1
znumer_flu =       Kf
zdenom_flu =       Kf + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._wp)

! H2O - OH
zalk_wat   = Kw/p_h - p_h/aphscale

equation_at =    zalk_dic + zalk_bor + zalk_po4 + zalk_sil &
               + zalk_so4 + zalk_flu                       &
               + zalk_wat - p_alktot

IF(PRESENT(p_deriveqn)) THEN
   ! H2CO3 - HCO3 - CO3 : n=2
   zdnumer_dic = K1*K1*K2 + p_h*(4._wp*K1*K2               &
                          + p_h*       K1    )
   zdalk_dic   = -p_dictot*(zdnumer_dic/zdenom_dic**2)

   ! B(OH)3 - B(OH)4 : n=1
   zdnumer_bor = Kb
   zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor**2)

   ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3
   zdnumer_po4 = K1p*K2p*K1p*K2p*K3p + p_h*(4._wp*K1p*K1p*K2p*K3p                &
                                     + p_h*(9._wp*K1p*K2p*K3p + K1p*K1p*K2p      &
                                     + p_h*(4._wp*K1p*K2p                        &
                                     + p_h*       K1p)))
   zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4**2)

   ! H4SiO4 - H3SiO4 : n=1
   zdnumer_sil = Ksi
   zdalk_sil   = -p_siltot * (zdnumer_sil/zdenom_sil**2)

!  ! NH4 - NH3 : n=1
!  zdnumer_nh4 = Knh4
!  zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4**2)

!  ! H2S - HS : n=1
!  zdnumer_h2s = api1_h2s
!  zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s**2)

   ! HSO4 - SO4 : n=1
   zdnumer_so4 = Ks
   zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4**2)

   ! HF - F : n=1
   zdnumer_flu = Kf
   zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu**2)

!  p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
!               + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu &
!               - Kw/p_h**2 - 1._wp/aphscale
   p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
                                        + zdalk_so4 + zdalk_flu &
                - Kw/p_h**2 - 1._wp/aphscale
ENDIF
RETURN
END FUNCTION equation_at

!===============================================================================

FUNCTION equation_at_DNAD(p_alktot, p_h,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                      &
                     p_so4tot, p_flutot,                                      &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                     p_deriveqn)
  ! Purpose: Compute total alkalinity from ion concentrations and equilibrium constants
  ! and compute its partial derivatives using dual numbers (DNAD) technique

  ! It is similar to subroutine equation_at above except that it computes
  ! partial derivatives of total alkalinity
  ! with respect to four input variables : ta, tc pt and sit

USE msingledouble
IMPLICIT NONE
TYPE(DUAL_NUM) ::equation_at_DNAD

! Argument variables
TYPE(DUAL_NUM), INTENT(IN)            :: p_alktot
TYPE(DUAL_NUM), INTENT(IN)            :: p_h
TYPE(DUAL_NUM), INTENT(IN)            :: p_dictot
TYPE(DUAL_NUM), INTENT(IN)            :: p_bortot
TYPE(DUAL_NUM), INTENT(IN)            :: p_po4tot
TYPE(DUAL_NUM), INTENT(IN)            :: p_siltot
!TYPE(DUAL_NUM), INTENT(IN)            :: p_nh4tot
!TYPE(DUAL_NUM), INTENT(IN)            :: p_h2stot
TYPE(DUAL_NUM), INTENT(IN)            :: p_so4tot
TYPE(DUAL_NUM), INTENT(IN)            :: p_flutot
TYPE(DUAL_NUM), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
TYPE(DUAL_NUM), INTENT(IN)            :: K1p, K2p, K3p, Ksi
TYPE(DUAL_NUM), INTENT(OUT), OPTIONAL :: p_deriveqn

! Local variables 
!-----------------
TYPE(DUAL_NUM) ::znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
TYPE(DUAL_NUM) ::znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
TYPE(DUAL_NUM) ::znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
TYPE(DUAL_NUM) ::znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
TYPE(DUAL_NUM) ::znumer_nh4, zdnumer_nh4, zdenom_nh4, zalk_nh4, zdalk_nh4
TYPE(DUAL_NUM) ::znumer_h2s, zdnumer_h2s, zdenom_h2s, zalk_h2s, zdalk_h2s
TYPE(DUAL_NUM) ::znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
TYPE(DUAL_NUM) ::znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
TYPE(DUAL_NUM) ::                                     zalk_wat, zdalk_wat
TYPE(DUAL_NUM) ::aphscale

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

! H2CO3 - HCO3 - CO3 : n=2, m=0
znumer_dic = 2._wp*K1*K2 + p_h*       K1
zdenom_dic =       K1*K2 + p_h*(      K1 + p_h)
zalk_dic   = p_dictot * (znumer_dic/zdenom_dic)

! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       Kb
zdenom_bor =       Kb + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._wp*K1p*K2p*K3p + p_h*(2._wp*K1p*K2p + p_h* K1p)
zdenom_po4 =       K1p*K2p*K3p + p_h*(      K1p*K2p + p_h*(K1p + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._wp) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       Ksi
zdenom_sil =       Ksi + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! NH4 - NH3 : n=1, m=0
!znumer_nh4 =       api1_nh4
!zdenom_nh4 =       api1_nh4 + p_h
!zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)
! Note: api1_nh4 = Knh4

! H2S - HS : n=1, m=0
!znumer_h2s =       api1_h2s
!zdenom_h2s =       api1_h2s + p_h
!zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)
! Note: api1_h2s = Kh2s

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       Ks
zdenom_so4 =       Ks + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._wp)

! HF - F : n=1, m=1
znumer_flu =       Kf
zdenom_flu =       Kf + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._wp)

! H2O - OH
zalk_wat   = Kw/p_h - p_h/aphscale

equation_at_DNAD =  zalk_dic + zalk_bor + zalk_po4 + zalk_sil &
                  + zalk_so4 + zalk_flu                       &
                  + zalk_wat - p_alktot

IF(PRESENT(p_deriveqn)) THEN
   ! H2CO3 - HCO3 - CO3 : n=2
   zdnumer_dic = K1*K1*K2 + p_h*(4._wp*K1*K2               &
                          + p_h*       K1    )
   zdalk_dic   = -p_dictot*(zdnumer_dic/zdenom_dic**2)

   ! B(OH)3 - B(OH)4 : n=1
   zdnumer_bor = Kb
   zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor**2)

   ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3
   zdnumer_po4 = K1p*K2p*K1p*K2p*K3p + p_h*(4._wp*K1p*K1p*K2p*K3p                &
                                     + p_h*(9._wp*K1p*K2p*K3p + K1p*K1p*K2p      &
                                     + p_h*(4._wp*K1p*K2p                        &
                                     + p_h*       K1p)))
   zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4**2)

   ! H4SiO4 - H3SiO4 : n=1
   zdnumer_sil = Ksi
   zdalk_sil   = -p_siltot * (zdnumer_sil/zdenom_sil**2)

!  ! NH4 - NH3 : n=1
!  zdnumer_nh4 = Knh4
!  zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4**2)

!  ! H2S - HS : n=1
!  zdnumer_h2s = api1_h2s
!  zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s**2)

   ! HSO4 - SO4 : n=1
   zdnumer_so4 = Ks
   zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4**2)

   ! HF - F : n=1
   zdnumer_flu = Kf
   zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu**2)

!  p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
!               + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu &
!               - Kw/p_h**2 - 1._wp/aphscale
   p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
                                        + zdalk_so4 + zdalk_flu &
                - Kw/p_h**2 - 1._wp/aphscale
ENDIF
RETURN
END FUNCTION equation_at_DNAD

!===============================================================================

SUBROUTINE ahini_for_at(p_alkcb, p_dictot, p_bortot, K1, K2, Kb, p_hini)

! Subroutine returns the root for the 2nd order approximation of the
! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
! around the local minimum, if it exists.

! Returns * 1E-03_wp if p_alkcb <= 0
!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
!                    and the 2nd order approximation does not have a solution

!USE MOD_CHEMCONST, ONLY : api1_dic, api2_dic, api1_bor

USE msingledouble
IMPLICIT NONE

! Argument variables 
!--------------------
REAL(KIND=wp), INTENT(IN)   ::  p_alkcb, p_dictot, p_bortot
REAL(KIND=wp), INTENT(IN)   ::  K1, K2, Kb
REAL(KIND=wp), INTENT(OUT)  ::  p_hini

! Local variables 
!-----------------
REAL(KIND=wp)  ::  zca, zba
REAL(KIND=wp)  ::  zd, zsqrtd, zhmin
REAL(KIND=wp)  ::  za2, za1, za0

IF (p_alkcb <= 0._wp) THEN
  p_hini = 1.e-3_wp
ELSEIF (p_alkcb >= (2._wp*p_dictot + p_bortot)) THEN
  p_hini = 1.e-10_wp
ELSE
  zca = p_dictot/p_alkcb
  zba = p_bortot/p_alkcb

  ! Coefficients of the cubic polynomial
  za2 = Kb*(1._wp - zba) + K1*(1._wp-zca)
  za1 = K1*Kb*(1._wp - zba - zca) + K1*K2*(1._wp - (zca+zca))
  za0 = K1*K2*Kb*(1._wp - zba - (zca+zca))
                                        ! Taylor expansion around the minimum
  zd = za2*za2 - 3._wp*za1              ! Discriminant of the quadratic equation
                                        ! for the minimum close to the root

  IF(zd > 0._wp) THEN                   ! If the discriminant is positive
    zsqrtd = SQRT(zd)
    IF(za2 < 0) THEN
      zhmin = (-za2 + zsqrtd)/3._wp
    ELSE
      zhmin = -za1/(za2 + zsqrtd)
    ENDIF
    p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
  ELSE
    p_hini = 1.e-7_wp
  ENDIF

ENDIF
RETURN
END SUBROUTINE ahini_for_at

!===============================================================================

FUNCTION solve_at_general(p_alktot, p_dictot, p_bortot,                       &
                          p_po4tot, p_siltot,                                 &
                          p_so4tot, p_flutot,                                 &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                          p_hini,   p_val)

! Purpose: Compute [H°] ion concentration from sea-water ion concentrations,
!          alkalinity, DIC, and equilibrium constants
! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_GENERAL

! Argument variables 
!--------------------
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables 
!-----------------
REAL(KIND=wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=wp)  ::  zalknw_inf, zalknw_sup
REAL(KIND=wp)  ::  zh_min, zh_max
REAL(KIND=wp)  ::  zdelta, zh_delta
REAL(KIND=wp)  ::  zeqn, zdeqndh, zeqn_absmin
REAL(KIND=wp)  ::  aphscale
LOGICAL        :: l_exitnow
REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL ahini_for_at(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF

 CALL anw_infsup(p_dictot, p_bortot,                                           &
                 p_po4tot, p_siltot,                                           &
                 p_so4tot, p_flutot,                                           &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
   zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atgen        = 0                 ! Reset counters of iterations
zeqn_absmin        = HUGE(1._wp)

DO
   IF(niter_atgen >= jp_maxniter_atgen) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zh_prev = zh
   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                      P_DERIVEQN = zdeqndh)

   ! Adapt bracketing interval
   IF(zeqn > 0._wp) THEN
      zh_min = zh_prev
   ELSEIF(zeqn < 0._wp) THEN
      zh_max = zh_prev
   ELSE
      ! zh is the root; unlikely but, one never knows
      EXIT
   ENDIF

   ! Now determine the next iterate zh
   niter_atgen = niter_atgen + 1

   IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
      ! if the function evaluation at the current point is
      ! not decreasing faster than with a bisection step (at least linearly)
      ! in absolute value take one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      !
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
      zh = SQRT(zh_max * zh_min)
      zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
   ELSE
      ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
      !           = -zdeqndh * LOG(10) * [H]
      ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
      !
      ! pH_new = pH_old + \deltapH
      !
      ! [H]_new = 10**(-pH_new)
      !         = 10**(-pH_old - \Delta pH)
      !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

      zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

      IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
         zh          = zh_prev*EXP(zh_lnfactor)
      ELSE
         zh_delta    = zh_lnfactor*zh_prev
         zh          = zh_prev + zh_delta
      ENDIF

      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)
         zh                = SQRT(zh_prev * zh_min)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)
         zh                = SQRT(zh_prev * zh_max)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF
   ENDIF

   zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <-- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT
ENDDO

solve_at_general = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)    
   ELSE
      p_val = HUGE(1._wp)
   ENDIF
ENDIF
RETURN
END FUNCTION solve_at_general

!===============================================================================

FUNCTION solve_at_general_DNAD (p_alktot, p_dictot, p_bortot,                 &
                          p_po4tot, p_siltot,                                 &
                          p_so4tot, p_flutot,                                 &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,     &
                          p_hini,   p_val)

! Purpose: Compute [H+] ion concentration from sea-water ion concentrations,
! alkalinity, DIC, and equilibrium constants
! It also computes partial derivatives of [H+] using dual numbers (DNAD) technique

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

! It is similar to subroutine solve_at_general above except that it computes
! partial derivatives of [H°] ion concentration
! with respect to four input variables : ta, tc pt and sit

! NOTE : all input and output parameters are of type DUAL_NUM (dual number)
!        where they are of type REAL in standard version of this subroutine.
!        But only four input parameters are required to be passed as such :
!         p_alktot, p_dictot, p_po4tot and p_siltot

USE msingledouble
IMPLICIT NONE
! value and partial derivatives of [H+] are returned in one objet of type DUAL_NUM
TYPE(DUAL_NUM) ::SOLVE_AT_GENERAL_DNAD

! Argument variables 
!--------------------
TYPE(DUAL_NUM), INTENT(IN)            :: p_alktot
TYPE(DUAL_NUM), INTENT(IN)            :: p_dictot
TYPE(DUAL_NUM), INTENT(IN)            :: p_bortot
TYPE(DUAL_NUM), INTENT(IN)            :: p_po4tot
TYPE(DUAL_NUM), INTENT(IN)            :: p_siltot
!TYPE(DUAL_NUM), INTENT(IN)            :: p_nh4tot
!TYPE(DUAL_NUM), INTENT(IN)            :: p_h2stot
TYPE(DUAL_NUM), INTENT(IN)            :: p_so4tot
TYPE(DUAL_NUM), INTENT(IN)            :: p_flutot
TYPE(DUAL_NUM), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
TYPE(DUAL_NUM), INTENT(IN)            :: K1p, K2p, K3p, Ksi
TYPE(DUAL_NUM), INTENT(IN), OPTIONAL  :: p_hini
TYPE(DUAL_NUM), INTENT(OUT), OPTIONAL :: p_val

! Local variables 
!-----------------
REAL(KIND=wp)  ::  zh_ini
TYPE(DUAL_NUM) ::  zh, zh_prev, zh_lnfactor
REAL(KIND=wp)  ::  zalknw_inf, zalknw_sup
TYPE(DUAL_NUM) ::  zh_min, zh_max
TYPE(DUAL_NUM) ::  zdelta, zh_delta
TYPE(DUAL_NUM) ::  zeqn, zdeqndh, zeqn_absmin
TYPE(DUAL_NUM) ::  aphscale
LOGICAL        :: l_exitnow
REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini%x_ad_
ELSE
   CALL ahini_for_at(p_alktot%x_ad_, p_dictot%x_ad_, p_bortot%x_ad_, K1%x_ad_, K2%x_ad_, Kb%x_ad_, zh_ini)
ENDIF

 CALL anw_infsup(p_dictot%x_ad_, p_bortot%x_ad_,                               &
                 p_po4tot%x_ad_, p_siltot%x_ad_,                               &
                 p_so4tot%x_ad_, p_flutot%x_ad_,                               &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
   zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atgen        = 0                 ! Reset counters of iterations
zeqn_absmin        = HUGE(1._wp)

DO
   IF(niter_atgen >= jp_maxniter_atgen) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zh_prev = zh
   zeqn = equation_at_DNAD(p_alktot, zh,       p_dictot, p_bortot,             &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                      P_DERIVEQN = zdeqndh)

   ! Adapt bracketing interval
   IF(zeqn > 0._wp) THEN
      zh_min = zh_prev
   ELSEIF(zeqn < 0._wp) THEN
      zh_max = zh_prev
   ELSE
      ! zh is the root; unlikely but, one never knows
      EXIT
   ENDIF

   ! Now determine the next iterate zh
   niter_atgen = niter_atgen + 1

   IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
      ! if the function evaluation at the current point is
      ! not decreasing faster than with a bisection step (at least linearly)
      ! in absolute value take one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      !
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
      zh = SQRT(zh_max * zh_min)
      zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
   ELSE
      ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
      !           = -zdeqndh * LOG(10) * [H]
      ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
      !
      ! pH_new = pH_old + \deltapH
      !
      ! [H]_new = 10**(-pH_new)
      !         = 10**(-pH_old - \Delta pH)
      !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

      zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

      IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
         zh          = zh_prev*EXP(zh_lnfactor)
      ELSE
         zh_delta    = zh_lnfactor*zh_prev
         zh          = zh_prev + zh_delta
      ENDIF

      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)
         zh                = SQRT(zh_prev * zh_min)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)
         zh                = SQRT(zh_prev * zh_max)
         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
      ENDIF
   ENDIF

   zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)

   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <-- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT
ENDDO

solve_at_general_DNAD = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at_DNAD(p_alktot, zh,       p_dictot, p_bortot,         &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)    
   ELSE
      p_val = HUGE(1._wp)
   ENDIF
ENDIF
RETURN
END FUNCTION solve_at_general_DNAD

!===============================================================================

FUNCTION solve_at_general_sec(p_alktot, p_dictot, p_bortot,                   &
                              p_po4tot, p_siltot,                             &
                              p_so4tot, p_flutot,                             &
                              K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi, &
                              p_hini,   p_val) 

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

!USE MOD_CHEMCONST, ONLY: api1_wat, aphscale
USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_GENERAL_SEC

! Argument variables 
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
REAL(KIND=wp)  ::  zh_ini, zh, zh_1, zh_2, zh_delta
REAL(KIND=wp)  ::  zalknw_inf, zalknw_sup
REAL(KIND=wp)  ::  zh_min, zh_max
REAL(KIND=wp)  ::  zeqn, zeqn_1, zeqn_2, zeqn_absmin
REAL(KIND=wp)  ::  zdelta
REAL(KIND=wp)  ::  aphscale
LOGICAL        ::  l_exitnow

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL ahini_for_at(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF

 CALL anw_infsup(p_dictot, p_bortot,                                      &
                 p_po4tot, p_siltot,                                      &
                 p_so4tot, p_flutot,                                      &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._wp*Kw/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._wp*Kw /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._wp
ENDIF

zdelta = (p_alktot-zalknw_sup)**2 + 4._wp*Kw/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._wp
ELSE
   zh_max = 2._wp*Kw /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

zh = MAX(MIN(zh_max, zh_ini), zh_min)
niter_atsec        = 0                 ! Reset counters of iterations

! Prepare the secant iterations: two initial (zh, zeqn) pairs are required
! We have the starting value, that needs to be completed by the evaluation
! of the equation value it produces.

! Complete the initial value with its equation evaluation
! (will take the role of the $n-2$ iterate at the first secant evaluation)

niter_atsec = 0                        ! zh_2 is the initial value;

zh_2   = zh
zeqn_2 = equation_at(p_alktot, zh_2,     p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                     &
                     p_so4tot, p_flutot,                                     &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

zeqn_absmin        = ABS(zeqn_2)

! Adapt bracketing interval and heuristically set zh_1
IF(zeqn_2 < 0._wp) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_max:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_max = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) below the current zh_max,
                                       ! but stay above SQRT(zh_min*zh_max), which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MAX(zh_max/1.25_wp, SQRT(zh_min*zh_max))
ELSEIF(zeqn_2 > 0._wp) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_min:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_min = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) above the current zh_min,
                                       ! but stay below SQRT(zh_min*zh_max) which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MIN(zh_min*1.25_wp, SQRT(zh_min*zh_max))
ELSE ! we have got the root; unlikely, but one never knows
   solve_at_general_sec = zh_2
   IF(PRESENT(p_val)) p_val = zeqn_2
   RETURN
ENDIF

! We now have the first pair completed (zh_2, zeqn_2).
! Define the second one (zh_1, zeqn_1), which is also the first iterate.
! zh_1 has already been set above
niter_atsec = 1                        ! Update counter of iterations

zeqn_1 = equation_at(p_alktot, zh_1,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot,                                       &
                     p_so4tot, p_flutot,                                       &
                     K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

! Adapt bracketing interval: we know that zh_1 <= zh <= zh_max (if zeqn_1 > 0)
! or zh_min <= zh <= zh_1 (if zeqn_1 < 0), so this can always be done
IF(zeqn_1 > 0._wp) THEN
   zh_min = zh_1
ELSEIF(zeqn_1 < 0._wp) THEN
   zh_max = zh_1
ELSE ! zh_1 is the root
   solve_at_general_sec = zh_1
   IF(PRESENT(p_val)) p_val = zeqn_1
ENDIF

IF(ABS(zeqn_1) > zeqn_absmin) THEN     ! Swap zh_2 and zh_1 if ABS(zeqn_2) < ABS(zeqn_1)
                                       ! so that zh_2 and zh_1 lead to decreasing equation
                                       ! values (in absolute value)
   zh     = zh_1
   zeqn   = zeqn_1
   zh_1   = zh_2
   zeqn_1 = zeqn_2
   zh_2   = zh
   zeqn_2 = zeqn
ELSE
   zeqn_absmin = ABS(zeqn_1)
ENDIF

! Pre-calculate the first secant iterate (this is the second iterate)
niter_atsec = 2

zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))
zh = zh_1 + zh_delta

! Make sure that zh_min < zh < zh_max (if not,
! bisect around zh_1 which is the best estimate)
IF (zh > zh_max) THEN                  ! this can only happen if zh_2 < zh_1
                                       ! and zeqn_2 > zeqn_1 > 0
   zh = SQRT(zh_1*zh_max)
ENDIF

IF (zh < zh_min) THEN                  ! this can only happen if zh_2 > zh_1
                                       ! and zeqn_2 < zeqn_1 < 0
   zh = SQRT(zh_1*zh_min)
ENDIF

DO
   IF(niter_atsec >= jp_maxniter_atsec) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      &
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)

   ! Adapt bracketing interval: since initially, zh_min <= zh <= zh_max
   ! we are sure that zh will improve either bracket, depending on the sign
   ! of zeqn
   IF(zeqn > 0._wp) THEN
     zh_min = zh
   ELSEIF(zeqn < 0._wp) THEN
     zh_max = zh
   ELSE
     ! zh is the root
     EXIT
   ENDIF

   ! start calculation of next iterate
   niter_atsec = niter_atsec + 1

   zh_2   = zh_1
   zeqn_2 = zeqn_1
   zh_1   = zh
   zeqn_1 = zeqn

   IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin) THEN
      ! if the function evaluation at the current point
      ! is not decreasing faster in absolute value than
      ! we may expect for a bisection step, then take
      ! one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
      zh                = SQRT(zh_max * zh_min)
      zh_delta           = zh - zh_1
   ELSE
      ! \Delta H = -zeqn_1*(h_2 - h_1)/(zeqn_2 - zeqn_1) 
      ! H_new = H_1 + \Delta H
      zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))
      zh       = zh_1 + zh_delta

      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)
         zh                = SQRT(zh_1 * zh_min)
         zh_delta          = zh - zh_1
      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)
         zh                 = SQRT(zh_1 * zh_max)
         zh_delta           = zh - zh_1
      ENDIF
   ENDIF

   zeqn_absmin = MIN(ABS(zeqn), zeqn_absmin)

   ! Stop iterations once |([H]-[H_1])/[H_1]| < rdel
   l_exitnow = (ABS(zh_delta) < pp_rdel_ah_target*zh_1)

   IF(l_exitnow) EXIT
ENDDO

SOLVE_AT_GENERAL_SEC = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
   ELSE
     p_val = HUGE(1._wp)
   ENDIF
ENDIF

RETURN
END FUNCTION SOLVE_AT_GENERAL_SEC

!===============================================================================

FUNCTION SOLVE_AT_FAST(p_alktot, p_dictot, p_bortot,                          &
                       p_po4tot, p_siltot,                                    &
                       p_so4tot, p_flutot,                                    & 
                       K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,        &
                       p_hini,   p_val)

! Fast version of SOLVE_AT_GENERAL, without any bounds checking.

USE msingledouble
IMPLICIT NONE
REAL(KIND=wp) :: SOLVE_AT_FAST

! Argument variables 
REAL(KIND=wp), INTENT(IN)            :: p_alktot
REAL(KIND=wp), INTENT(IN)            :: p_dictot
REAL(KIND=wp), INTENT(IN)            :: p_bortot
REAL(KIND=wp), INTENT(IN)            :: p_po4tot
REAL(KIND=wp), INTENT(IN)            :: p_siltot
!REAL(KIND=wp), INTENT(IN)            :: p_nh4tot
!REAL(KIND=wp), INTENT(IN)            :: p_h2stot
REAL(KIND=wp), INTENT(IN)            :: p_so4tot
REAL(KIND=wp), INTENT(IN)            :: p_flutot
REAL(KIND=wp), INTENT(IN)            :: K0, K1, K2, Kb, Kw, Ks, Kf
REAL(KIND=wp), INTENT(IN)            :: K1p, K2p, K3p, Ksi
REAL(KIND=wp), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=wp), INTENT(OUT), OPTIONAL :: p_val

! Local variables
REAL(KIND=wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=wp)  ::  zhdelta
REAL(KIND=wp)  ::  zeqn, zdeqndh
!REAL(KIND=wp)  ::  aphscale
LOGICAL        :: l_exitnow
REAL(KIND=wp), PARAMETER :: pz_exp_threshold = 1.0_wp

! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
!aphscale = 1._wp + p_so4tot/Ks

IF(PRESENT(p_hini)) THEN
   zh_ini = p_hini
ELSE
   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, K1, K2, Kb, zh_ini)
ENDIF
zh = zh_ini

niter_atfast    = 0                 ! Reset counters of iterations
DO
   niter_atfast = niter_atfast + 1
   IF(niter_atfast > jp_maxniter_atfast) THEN
      zh = -1._wp
      EXIT
   ENDIF

   zh_prev = zh

   zeqn = equation_at(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot,                                      &
                      p_so4tot, p_flutot,                                      & 
                      K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi,          &
                      P_DERIVEQN = zdeqndh)

   IF(zeqn == 0._wp) EXIT               ! zh is the root

   zh_lnfactor = -zeqn/(zdeqndh*zh_prev)
   IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
      zh          = zh_prev*EXP(zh_lnfactor)
   ELSE
      zhdelta     = zh_lnfactor*zh_prev
      zh          = zh_prev + zhdelta
   ENDIF

   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT
ENDDO

SOLVE_AT_FAST = zh

IF(PRESENT(p_val)) THEN
   IF(zh > 0._wp) THEN
      p_val = equation_at(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot,                                  &
                          p_so4tot, p_flutot,                                  &
                          K0, K1, K2, Kb, Kw, Ks, Kf, K1p, K2p, K3p, Ksi)
   ELSE
      p_val = HUGE(1._wp)
   ENDIF
ENDIF

RETURN
END FUNCTION solve_at_fast
!===============================================================================

END MODULE mphsolvers
