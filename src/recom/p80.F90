!> \file p80.f90
!! \BRIEF 
!> Module with p80 function - compute pressure from depth
MODULE mp80
CONTAINS
!>     Function to compute pressure from depth using Saunder's (1981) formula with eos80.
FUNCTION p80(dpth,xlat)

  !     Compute Pressure from depth using Saunder's (1981) formula with eos80.

  !     Reference:
  !     Saunders, Peter M. (1981) Practical conversion of pressure
  !     to depth, J. Phys. Ooceanogr., 11, 573-574, (1981)

  !     Coded by:
  !     R. Millard
  !     March 9, 1983
  !     check value: p80=7500.004 dbars at lat=30 deg., depth=7321.45 meters

  !     Modified (slight format changes + added ref. details):
  !     J. Orr, 16 April 2009

  USE msingledouble
  IMPLICIT NONE

! Input variables:
  !> depth [m]
  REAL(kind=rx), INTENT(in) :: dpth
  !> latitude [degrees]
  REAL(kind=rx), INTENT(in) :: xlat

! Output variable:
  !> pressure [db]
  REAL(kind=rx) :: p80

! Local variables:
  REAL(kind=rx) :: pi
  REAL(kind=rx) :: plat, d, c1

  pi=3.141592654

  plat = ABS(xlat*pi/180.)
  d  = SIN(plat)
  c1 = 5.92e-3+d**2 * 5.25e-3

  p80 = ((1-c1)-SQRT(((1-c1)**2)-(8.84e-6*dpth))) / 4.42e-6

  RETURN
END FUNCTION p80
END MODULE mp80
