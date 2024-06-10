
!> \file singledouble.f90
!! \BRIEF 
!> Module that defines single and double precision - used by all other modules
MODULE msingledouble
! INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)
! INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12)

  INTEGER, PARAMETER :: rx = KIND(1.0d0)

  INTEGER, PARAMETER :: r8 = KIND(1.0d0)
  INTEGER, PARAMETER :: wp = KIND(1.0d0)
END MODULE msingledouble