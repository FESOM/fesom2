!> \file depth2press.f90
!! \BRIEF 
!> Module with depth2press subroutine - converts depth to pressure
!! with Saunders (1981) formula
MODULE mdepth2press
CONTAINS
!>     Compute pressure [db] from depth [m] & latitude [degrees north].
!!     This subroutine is needed because p80 is a function (using scalars not arrays)
SUBROUTINE depth2press(depth, lat, pdbar, N)

  !     Purpose:
  !     Compute pressure [db] from depth [m] & latitude [degrees north].
  !     Needed because p80 is a function 

  USE msingledouble
  USE mp80
  IMPLICIT NONE

  !> number of records
!f2py intent(hide) n
  INTEGER, intent(in) :: N

! INPUT variables
  !> depth [m]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: depth
  !> latitude [degrees]
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: lat

! OUTPUT variables:
  !> pressure [db]
  REAL(kind=rx), INTENT(out), DIMENSION(N) :: pdbar

  !     Local variables
  INTEGER :: i

! REAL(kind=rx) ::  p80
! EXTERNAL p80

  DO i = 1,N
     pdbar(i) = p80(depth(i), lat(i))
  END DO

  RETURN
END SUBROUTINE depth2press
END MODULE mdepth2press
