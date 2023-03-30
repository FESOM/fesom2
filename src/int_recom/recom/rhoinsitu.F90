!> \file rhoinsitu.f90
!! \BRIEF 
!> Module with rhoinsitu subroutine - compute in situ density from S, Tis, P
MODULE mrhoinsitu
CONTAINS
!>     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db).
!!     This subroutine is needed because rho is a function (using scalars not arrays)
SUBROUTINE rhoinsitu(salt, tempis, pdbar, N, rhois)

  !     Purpose:
  !     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db)
  !     Needed because rho is a function

  USE msingledouble
  USE mrho
  IMPLICIT NONE

  INTEGER :: N

! INPUT variables
  ! salt   = salinity [psu]
  ! tempis = in situ temperature [C]
  ! pdbar  = pressure [db]

  !> salinity [psu]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: pdbar
!f2py optional , depend(salt) :: n=len(salt)

! OUTPUT variables:
  ! rhois  = in situ density

  !> in situ density [kg/m3]
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: rhois

! Local variables
  INTEGER :: i

! REAL(kind=rx) ::  rho
! EXTERNAL rho

  DO i = 1,N
     rhois(i) = rho(salt(i), tempis(i), pdbar(i)/10.)
  END DO

  RETURN
END SUBROUTINE rhoinsitu
END MODULE mrhoinsitu
