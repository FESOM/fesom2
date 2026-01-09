!===============================================================================
! Subroutine for calculating flux-depth and thickness of control volumes
!===============================================================================
subroutine Depth_calculations(n,Nn,wF,zF,thick,recipthick, partit, mesh)
  use recom_config
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
  USE o_PARAM
  USE o_ARRAYS
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto

  use g_clock
  use g_rotate_grid

  implicit none
! Input
    type(t_partit), intent(inout),   target         :: partit
    type(t_mesh)  , intent(inout),   target         :: mesh
 
  Integer,                           intent(in)     :: Nn	     ! Total number of nodes
! Output
  real(kind=8),dimension(mesh%nl),   intent(out)    :: zF             ! [m] Depth of vertical fluxes

  real(kind=8),dimension(mesh%nl-1), intent(out)    :: thick          ! [m] Distance between two nodes = thickness
  real(kind=8),dimension(mesh%nl-1), intent(out)    :: recipthick     ! [1/m] reciprocal of thickness

  real(kind=8),dimension(mesh%nl,5), intent(out)    :: wF             ! [m/day] Velocities of fluxes at the border of the control volumes
  Integer                                           :: k, n           ! Index for depth      

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
! ======================================================================================
!! zbar(nl) allocate the array for storing the standard depths (depth of layers)
!! zbar is negative 
!! Z(nl-1)  mid-depths of layers

    ! zbar_n: depth of layers due to ale thinkness variations at every node n 
!    allocate(zbar_n(nl))
!    allocate(zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))
    
    ! Z_n: mid depth of layers due to ale thinkness variations at every node n 
!    allocate(Z_n(nl-1))
!    allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
! ============================================================================== modular

!! Background sinking speed

  wF(2:Nn,ivphy)   = VPhy  
  wF(2:Nn,ivdia)   = VDia
  wF(2:Nn,ivdet)   = VDet
  wF(2:Nn,ivdetsc) = VDet_zoo2
  wF(2:Nn,ivcoc)   = VCocco

  wF(1,:)          = 0.d0
  wF(Nn+1,:)       = 0.d0

!if (allow_var_sinking) then
!!    wF(2:Nn+1,ivdet) = Vdet_a * abs(zbar_n(2:Nn+1)) + VDet
!!YY: use Vcalc instead of Vdet_a, only needed for calculating calc_diss
!    wF(2:Nn+1,ivdet) = Vcalc * abs(zbar_3d_n(2:Nn+1, n)) + VDet
!end if

!----------------------------------------------------
! Vertical layers thickness

    thick   =0.0_WP
    recipthick=0.0_WP
    zF=0.0_WP

!do n=1,myDim_nod2D+eDim_nod2D
   do k=1, Nn !nlevels_nod2D(n)-1
      thick(k)=hnode(k,n)
      if (hnode(k,n) > 0._WP) then
         recipthick(k)=1./hnode(k,n)
      else
         recipthick(k)=0._WP
      end if
   end do
!end do

! layer depth (negative)

!do n=1,myDim_nod2D+eDim_nod2D
   do k=1, Nn+1 !nlevels_nod2D(n)
      zF(k)=zbar_3d_n(k,n)
   end do
!end do
  
end subroutine Depth_calculations

!===============================================================================
! Subroutine for calculating cos(AngleOfIncidence)
!===============================================================================
subroutine Cobeta(partit, mesh)
    use REcoM_GloVar
    use g_clock
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use o_PARAM
    use g_comm_auto
    Implicit none
  	
! Declarations
    real(kind=8)                          :: yearfrac              ! The fraction of the year that has passed [0 1]
    real(kind=8)                          :: yDay                  ! yearfrac in radians [0 2*pi]
    real(kind=8)                          :: declination   = 0.d0  ! Declination of the sun at present lat and time
    real(kind=8)                          :: CosAngleNoon  = 0.d0  ! Cos(Angle of Incidence) at Noon ?
    integer                               :: n

! Constants
    real(kind=8)		          :: nWater        = 1.33

    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

!! find day (****NOTE for year starting in winter*****)  
!! Paltridge, G. W. and C. M. R. Platt, Radiative Processes 
!! in Meteorology and Climatology, Developments in 
!! Atmospheric Sciences, vol. 5, Elsevier Scientific 
!! Publishing Company, Amsterdam, Oxford, 
!! New York, 1976, ISBN 0-444-41444-4.

  yearfrac    = mod(real(daynew),real(ndpyr))/real(ndpyr)
  yDay        = 2 * pi * yearfrac
  declination = 0.006918                   &
         	- 0.399912 * cos(    yDay) &
     	        + 0.070257 * sin(    yDay) &
          	- 0.006758 * cos(2 * yDay) &
        	+ 0.000907 * sin(2 * yDay) &
        	- 0.002697 * cos(3 * yDay) &
        	+ 0.001480 * sin(3 * yDay) 

  do n=1, myDim_nod2D!+eDim_nod2D 

        cosAngleNoon   =   sin(geo_coord_nod2D(2,n)) * sin(declination) &
     		         + cos(geo_coord_nod2D(2,n)) * cos(declination)
        cosAI(n)       = sqrt(1-(1-cosAngleNoon**2)/nWater**2)

  end do
end subroutine Cobeta
!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
 subroutine krill_resp(n, partit, mesh)
   use REcoM_declarations
   use REcoM_LocVar
   use REcoM_GloVar
   use g_clock
   use o_PARAM
!   use g_PARSUP
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
   use g_comm_auto
   implicit none
   integer                          :: n

    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"
 
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
