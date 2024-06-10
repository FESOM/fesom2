!> \file f2pCO2.f90
!! \BRIEF 
!>    Module with f2pCO2 subroutine - compute pCO2 from fCO2, in situ T, atm pressure, hydrostatic pressure
MODULE mf2pCO2
CONTAINS
!>    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure
SUBROUTINE f2pCO2(fCO2, temp, Patm, p, N, pCO2)
  !    Purpose:
  !    Compute pCO2 from arrays of fCO2, in situ temp, atm pressure, & hydrostatic pressure

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  IMPLICIT NONE

  !> number of records
!f2py intent(hide) n
  INTEGER, intent(in) :: N

! INPUT variables
  !> oceanic fugacity of CO2 [uatm]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: fCO2
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> atmospheric pressure [atm]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> hydrostatic pressure [db]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: p

! OUTPUT variables:
  !> oceanic partial pressure of CO2 [uatm] 
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pCO2

! LOCAL variables:
  REAL(kind=r8) :: dfCO2, dtemp, tk, dPatm, prb
  REAL(kind=r8) :: Ptot, Rgas_atm, B, Del, xCO2approx, xc2, fugcoeff
  REAL(kind=r8) :: dpCO2

  INTEGER :: i

! REAL(kind=r8) :: sw_ptmp
! EXTERNAL sw_ptmp

  DO i = 1,N
     dfCO2     = DBLE(fCO2(i))
     dtemp     = DBLE(temp(i))
     dPatm     = DBLE(Patm(i))
     tk = 273.15d0 + DBLE(temp(i))     !Absolute temperature (Kelvin)
     prb = DBLE(p(i)) / 10.0d0         !Pressure effect (prb is in bars)
     Ptot = dPatm + prb/1.01325d0      !Total pressure (atmospheric + hydrostatic) [atm]
     Rgas_atm = 82.05736_r8            !R in (cm3 * atm) / (mol * K)  from CODATA (2006)
!    To compute fugcoeff, we need 3 other terms (B, Del, xc2) as well as 3 others above (tk, Ptot, Rgas_atm)
     B = -1636.75d0 + 12.0408d0*tk - 0.0327957d0*(tk*tk) + 0.0000316528d0*(tk*tk*tk)
     Del = 57.7d0 - 0.118d0*tk
!    "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
!    x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
!    Let's assume that xCO2 = fCO2. Resulting fugcoeff is identical to 8th digit after the decimal.
     xCO2approx = dfCO2 * 1.e-6_r8
     xc2 = (1.0d0 - xCO2approx)**2 
     fugcoeff = exp( Ptot*(B + 2.0d0*xc2*Del)/(Rgas_atm*tk) )
     dpCO2 = dfCO2 / fugcoeff
     pCO2(i) = SGLE(dpCO2)
  END DO

  RETURN
END SUBROUTINE f2pCO2
END MODULE mf2pCO2
