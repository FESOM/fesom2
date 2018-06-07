! Set of modules needed by the interface to IFS.
!
! -Original code by Kristian Mogensen, ECMWF.

MODULE par_kind
  IMPLICIT NONE
  INTEGER, PUBLIC, PARAMETER ::          &  !: Floating point section
       sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
       dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
       wp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
       ik = SELECTED_INT_KIND(6)            !: integer precision 
END MODULE par_kind

MODULE interinfo

   ! Parallel regridding information

   USE parinter

   IMPLICIT NONE

   SAVE

   ! IFS to NEMO

   TYPE(parinterinfo) :: gausstoT,gausstoUV

   ! NEMO to IFS

   TYPE(parinterinfo) :: Ttogauss, UVtogauss

   ! Read parinterinfo on task 0 only and broadcast.

   LOGICAL :: lparbcast = .FALSE.

END MODULE interinfo


