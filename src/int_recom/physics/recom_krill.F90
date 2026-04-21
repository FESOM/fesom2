!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
 subroutine krill_resp(n, partit, mesh)
    use REcoM_declarations
    use REcoM_LocVar
    use REcoM_GloVar
    use g_clock
    use o_PARAM
    use mod_mesh
    use MOD_PARTIT
    use MOD_PARSUP
    use g_comm_auto

    implicit none

    ! Input parameters
    integer                                          :: n
    type(t_partit), intent(inout),   target          :: partit
    type(t_mesh)  , intent(inout),   target          :: mesh

#include "../../associate_part_def.h"
#include "../../associate_mesh_def.h"
#include "../../associate_part_ass.h"
#include "../../associate_mesh_ass.h"
 
   ! Values from FESOM                                                                                                 

   if (geo_coord_nod2D(2,n)<0.0_WP) then  !SH
      if ((daynew .le. 105)) then
       res_zoo2_a = 0.d0
      else if((105 .le. daynew).and.(daynew .le. 150)) then
       res_zoo2_a = (-1./90.)*daynew +7./6.
      else if((150 .lt. daynew).and.(daynew .lt. 250)) then
       res_zoo2_a = -0.5
      else if((250 .le. daynew).and.(daynew .le. 295)) then
       res_zoo2_a = (1/90.)*daynew - 59./18.
      else if((daynew .gt. 295)) then
       res_zoo2_a = 0.d0
      end if
   else
      if ((daynew .le. 65)) then
       res_zoo2_a = -0.5
      else if((285 .le. daynew).and.(daynew .le. 330)) then
       res_zoo2_a = (-1./90.)*daynew +57./18.
      else if((330 .lt. daynew)) then
       res_zoo2_a = -0.5
      else if((65 .le. daynew).and.(daynew .le. 110)) then
       res_zoo2_a = (1/90.)*daynew - 22./18.
      else if((110 .lt. daynew).and.(daynew .lt. 285)) then
       res_zoo2_a = 0.d0
      end if
   endif
 end subroutine krill_resp

