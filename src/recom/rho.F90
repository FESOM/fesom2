!> \file rho.f90
!! \BRIEF 
!> Module with rho function - computes in situ density from S, T, P
MODULE mrho
CONTAINS
!> Function to compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
FUNCTION rho(salt, temp, pbar)

  ! Compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  USE msingledouble
  IMPLICIT NONE

  !> salinity [psu]
  REAL(kind=rx) :: salt
  !> in situ temperature (C)
  REAL(kind=rx) :: temp
  !> pressure (bar) [Note units: this is NOT in decibars]
  REAL(kind=rx) :: pbar

  REAL(kind=r8) :: s, t, p
! REAL(kind=r8) :: t68
  REAL(kind=r8) :: X
  REAL(kind=r8) :: rhow, rho0
  REAL(kind=r8) :: a, b, c
  REAL(kind=r8) :: Ksbmw, Ksbm0, Ksbm
  REAL(kind=r8) :: drho

  REAL(kind=rx) :: rho

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = in situ temperature [degree C (IPTS-68)]
  !     p  = pressure            [bar] !!!!  (not in [db]

  s = DBLE(salt)
  t = DBLE(temp)
  p = DBLE(pbar)

! Convert the temperature on today's "ITS 90" scale to older IPTS 68 scale
! (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
! According to Best-Practices guide, line above should be commented & 2 lines below should be uncommented
! Guide's answer of rho (s=35, t=25, p=0) = 1023.343 is for temperature given on ITPS-68 scale
! t68 = (T - 0.0002) / 0.99975
! X = t68
! Finally, don't do the ITS-90 to IPTS-68 conversion (T input var now already on IPTS-68 scale)
  X = T

! Density of pure water
  rhow = 999.842594d0 + 6.793952e-2_r8*X          &
       -9.095290e-3_r8*X*X + 1.001685e-4_r8*X**3  &
       -1.120083e-6_r8*X**4 + 6.536332e-9_r8*X**5

! Density of seawater at 1 atm, P=0
  A = 8.24493e-1_r8 - 4.0899e-3_r8*X                         &
       + 7.6438e-5_r8*X*X - 8.2467e-7_r8*X**3 + 5.3875e-9_r8*X**4
  B = -5.72466e-3_r8 + 1.0227e-4_r8*X - 1.6546e-6_r8*X*X
  C = 4.8314e-4_r8

  rho0 = rhow + A*S + B*S*SQRT(S) + C*S**2.0d0

! Secant bulk modulus of pure water
! The secant bulk modulus is the average change in pressure
! divided by the total change in volume per unit of initial volume.
  Ksbmw = 19652.21d0 + 148.4206d0*X - 2.327105d0*X*X &
       + 1.360477e-2_r8*X**3 - 5.155288e-5_r8*X**4

! Secant bulk modulus of seawater at 1 atm
  Ksbm0 = Ksbmw + S*( 54.6746d0 - 0.603459d0*X + 1.09987e-2_r8*X**2 &
       - 6.1670e-5_r8*X**3) &
       + S*SQRT(S)*( 7.944e-2_r8 + 1.6483e-2_r8*X - 5.3009e-4_r8*X**2)

! Secant bulk modulus of seawater at S,T,P
  Ksbm = Ksbm0 &
       + P*(3.239908d0 + 1.43713e-3_r8*X + 1.16092e-4_r8*X**2 - 5.77905e-7_r8*X**3) &
       + P*S*(2.2838e-3_r8 - 1.0981e-5_r8*X - 1.6078e-6_r8*X**2) &
       + P*S*SQRT(S)*1.91075e-4_r8 &
       + P*P*(8.50935e-5_r8 - 6.12293e-6_r8*X + 5.2787e-8_r8*X**2) &
       + P*P*S*(-9.9348e-7_r8 + 2.0816e-8_r8*X + 9.1697e-10_r8*X**2)

! Density of seawater at S,T,P
  drho = rho0/(1.0d0 - P/Ksbm)
  rho = SGLE(drho)

  RETURN
END FUNCTION rho


!> Function to compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
!! and its derivatives with respect to salinity and temperature.
FUNCTION rho_DNAD(salt, temp, pbar)

  ! Compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
  !
  ! It is similar to subroutine 'rho' above except that it computes
  ! partial derivatives of in situ density
  ! with respect to salinity and temperature.

  USE msingledouble
  USE Dual_Num_Auto_Diff
  IMPLICIT NONE

  !> salinity [psu]
  TYPE (DUAL_NUM) :: salt
  !> in situ temperature (C)
  TYPE (DUAL_NUM) :: temp
  !> pressure (bar) [Note units: this is NOT in decibars]
  TYPE (DUAL_NUM) :: pbar

  TYPE (DUAL_NUM) :: s, p
! REAL(kind=r8) :: t68
  TYPE (DUAL_NUM) :: X
  TYPE (DUAL_NUM) :: rhow, rho0
  TYPE (DUAL_NUM) :: a, b, c
  TYPE (DUAL_NUM) :: Ksbmw, Ksbm0, Ksbm
  
  TYPE (DUAL_NUM) :: rho_DNAD

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = in situ temperature [degree C (IPTS-68)]
  !     p  = pressure            [bar] !!!!  (not in [db]

  s = salt
  p = pbar

! Convert the temperature on today's "ITS 90" scale to older IPTS 68 scale
! (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
! According to Best-Practices guide, line above should be commented & 2 lines below should be uncommented
! Guide's answer of rho (s=35, t=25, p=0) = 1023.343 is for temperature given on ITPS-68 scale
! t68 = (T - 0.0002) / 0.99975
! X = t68
! Finally, don't do the ITS-90 to IPTS-68 conversion (T input var now already on IPTS-68 scale)
  X = temp

! Density of pure water
  rhow = 999.842594d0 + 6.793952e-2_r8*X          &
       -9.095290e-3_r8*X*X + 1.001685e-4_r8*X**3  &
       -1.120083e-6_r8*X**4 + 6.536332e-9_r8*X**5

! Density of seawater at 1 atm, P=0
  A = 8.24493e-1_r8 - 4.0899e-3_r8*X                         &
       + 7.6438e-5_r8*X*X - 8.2467e-7_r8*X**3 + 5.3875e-9_r8*X**4
  B = -5.72466e-3_r8 + 1.0227e-4_r8*X - 1.6546e-6_r8*X*X
  C = 4.8314e-4_r8

  rho0 = rhow + A*S + B*S*SQRT(S) + C*S**2.0d0

! Secant bulk modulus of pure water
! The secant bulk modulus is the average change in pressure
! divided by the total change in volume per unit of initial volume.
  Ksbmw = 19652.21d0 + 148.4206d0*X - 2.327105d0*X*X &
       + 1.360477e-2_r8*X**3 - 5.155288e-5_r8*X**4

! Secant bulk modulus of seawater at 1 atm
  Ksbm0 = Ksbmw + S*( 54.6746d0 - 0.603459d0*X + 1.09987e-2_r8*X**2 &
       - 6.1670e-5_r8*X**3) &
       + S*SQRT(S)*( 7.944e-2_r8 + 1.6483e-2_r8*X - 5.3009e-4_r8*X**2)

! Secant bulk modulus of seawater at S,T,P
  Ksbm = Ksbm0 &
       + P*(3.239908d0 + 1.43713e-3_r8*X + 1.16092e-4_r8*X**2 - 5.77905e-7_r8*X**3) &
       + P*S*(2.2838e-3_r8 - 1.0981e-5_r8*X - 1.6078e-6_r8*X**2) &
       + P*S*SQRT(S)*1.91075e-4_r8 &
       + P*P*(8.50935e-5_r8 - 6.12293e-6_r8*X + 5.2787e-8_r8*X**2) &
       + P*P*S*(-9.9348e-7_r8 + 2.0816e-8_r8*X + 9.1697e-10_r8*X**2)

! Density of seawater at S,T,P
  rho_DNAD = rho0/(1.0d0 - P/Ksbm)

  RETURN
END FUNCTION rho_DNAD
END MODULE mrho
