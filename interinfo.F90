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
