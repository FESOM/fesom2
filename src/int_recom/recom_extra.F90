!===============================================================================
! Subroutine for calculating flux-depth and thickness of control volumes
!===============================================================================

subroutine Depth_calculations(n,Nn,wF,zF,thick,recipthick,mesh)
  use recom_config
  USE mod_MESH
  USE o_PARAM
  USE o_ARRAYS
  USE g_PARSUP
  USE g_CONFIG
  use g_forcing_arrays
  use g_comm_auto
  use i_param
  use i_arrays
  use i_therm_param
  use g_clock
  use g_rotate_grid
  use g_comm
  implicit none
! Input
  type(t_mesh), intent(in) , target                :: mesh
  Integer,                     intent(in)          :: Nn	     ! Total number of nodes
! Output
  real(kind=8),dimension(mesh%nl),intent(out)      :: zF             ! [m] Depth of vertical fluxes

  real(kind=8),dimension(mesh%nl-1),intent(out)    :: thick          ! [m] Distance between two nodes = thickness
  real(kind=8),dimension(mesh%nl-1),intent(out)    :: recipthick     ! [1/m] reciprocal of thickness

  real(kind=8),dimension(mesh%nl,4),intent(out)    :: wF             ! [m/day] Velocities of fluxes at the border of the control volumes
  Integer                                          :: k, n           ! Index for depth      
#include "../associate_mesh.h"
! ======================================================================================
!! zbar (depth of layers) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(nl) allocate the array for storing the standard depths
!! Z(nl-1)  mid-depths of cells

    ! zbar_n: depth of layers due to ale thinkness variations at every node n 
!    allocate(zbar_n(nl))
!    allocate(zbar_3d_n(nl,myDim_nod2D+eDim_nod2D))
    
    ! Z_n: mid depth of layers due to ale thinkness variations at every node n 
!    allocate(Z_n(nl-1))
!    allocate(Z_3d_n(nl-1,myDim_nod2D+eDim_nod2D)) 
! ============================================================================== modular

  wF(2:Nn,ivphy) = VPhy  
  wF(2:Nn,ivdia) = VDia
  wF(2:Nn,ivdet) = VDet
  wF(2:Nn,ivdetsc) = VDet_zoo2

  wF(1,:)         = 0.d0
  wF(Nn+1,:)      = 0.d0

!----------------------------------------------------
! calculate thickness of vertical layers

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
subroutine Cobeta(mesh)
  use REcoM_GloVar
  use g_clock
  use o_PARAM
  use g_PARSUP
  use mod_MESH
  use g_comm_auto
  Implicit none
  	
! Declarations
  Real(kind=8)                     :: yearfrac              ! The fraction of the year that has passed [0 1]
  Real(kind=8)                     :: yDay                  ! yearfrac in radians [0 2*pi]
  Real(kind=8)                     :: declination   = 0.d0  ! Declination of the sun at present lat and time
  Real(kind=8)                     :: CosAngleNoon  = 0.d0  ! Cos(Angle of Incidence) at Noon ?
  integer                          :: n

! Constants
  Real(kind=8)		           :: nWater        = 1.33
  type(t_mesh), intent(in), target :: mesh  
#include  "../associate_mesh.h"
!
! find day (****NOTE for year starting in winter*****)  
! Paltridge, G. W. and C. M. R. Platt, Radiative Processes in Meteorology and Climatology, Developments in Atmospheric Sciences, vol. 5, Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York, 1976, ISBN 0-444-41444-4.
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
! Subroutine controlling and calculating atm. dep. of Fe
!================================================================================

subroutine Atm_input(mesh)
  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  use g_clock
  use g_read_other_NetCDF
  use g_PARSUP
  use o_PARAM, only : mstep, WP
  use mod_MESH
  use REcoM_ciso
  implicit none
  type(t_mesh), intent(in) , target :: mesh  
#include "netcdf.inc"

  real(kind=8), allocatable :: ncdata(:)
  integer	                :: status, ncid, varid
  character(300)            :: Ironname, CO2filename, DustNfilename
  character(15)             :: Ironvari, CO2vari, Nvari
  integer, dimension(2)     :: istart, icount
  integer                   :: CO2start, CO2count
  integer                   :: firstyearofcurrentCO2cycle, totnumyear
  character(4)              :: currentCO2year_char
  logical                   :: do_read=.false.
  integer                   :: i
#include "../associate_mesh.h" 
 
!-Checking if files need to be opened---------------------------------------------
  if (mstep == 1) then ! The year has changed
     if (mype==0) write(*,*), month
     i=month

     if(UseFeDust) then
  
!-Testing if files exist for the year in question---------------------------------
        if (UseDustClimAlbani) then
            Ironname = trim(REcoMDataPath)//'DustClimMonthlyAlbani.nc'
            Ironvari     = 'DustClim'
        else  
            if((yearnew .LT. 1979) .OR.(yearnew .GT. 2008)) then
            Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
            Ironvari     = 'DustClim'
            else
               if (UseDustClim) then
                  Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
                  Ironvari     = 'DustClim'
               else
                  Ironname = trim(REcoMDataPath)//'DustMonthlyMahowald.nc'
                  Ironvari     = 'Dust'//cyearnew
               end if
            end if 
        end if    
!-Reading Iron--------------------------------------------------------------------
        if (mype==0) write(*,*) 'Updating Iron restoring data for month ', i,' from', Ironname     
        call read_2ddata_on_grid_NetCDF(Ironname, Ironvari, i, GloFeDust, mesh)  
    endif ! FeDust

else
  call monthly_event(do_read)
  if(do_read) then ! file is opened and read every month
     i=month+1
     if (i > 12) i=1             
     if(UseFeDust) then
  
!-Testing if files exist for the year in question---------------------------------
        if (UseDustClimAlbani) then
            Ironname = trim(REcoMDataPath)//'DustClimMonthlyAlbani.nc'
            Ironvari     = 'DustClim'
        else  
            if((yearnew .LT. 1979) .OR.(yearnew .GT. 2008)) then
            Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
            Ironvari     = 'DustClim'
            else
               if (UseDustClim) then
                  Ironname = trim(REcoMDataPath)//'DustClimMonthlyMahowald.nc'
                  Ironvari     = 'DustClim'
               else
                  Ironname = trim(REcoMDataPath)//'DustMonthlyMahowald.nc'
                  Ironvari     = 'Dust'//cyearnew
               end if
            end if 
        end if    
!-Reading Iron--------------------------------------------------------------------
        if (mype==0) write(*,*) 'Updating Iron restoring data for month ', i,' from', Ironname     
        call read_2ddata_on_grid_NetCDF(Ironname, Ironvari, i, GloFeDust, mesh)  
    endif ! FeDust
  endif ! time
endif    
!-Reading CO2----------------------------------------------------------------------
  if (mstep == 1) then ! The year has changed
    if (constant_CO2) then
      AtmCO2(:) = CO2_for_spinup
      !if (mype==0) write(*,*),'in atm_input: Atm CO2=',AtmCO2
      if (ciso) then
        AtmCO2_13(:) = CO2_for_spinup_13
        AtmCO2_14(:) = CO2_for_spinup_14
        if (mype==0) write(*,*),'in atm_input: Atm CO2_13=',AtmCO2_13
        if (mype==0) write(*,*),'in atm_input: Atm CO2_14=',AtmCO2_14
      end if
    else  

!     CO2filename = trim(REcoMDataPath)//'MonthlyAtmCO2_2019.nc'
!     CO2filename = trim(REcoMDataPath)//'MonthlyAtmCO2_gcb2020.nc'
     CO2filename = trim(REcoMDataPath)//'MonthlyAtmCO2_gcb2021.nc'

     totnumyear                 = lastyearoffesomcycle-firstyearoffesomcycle+1
     firstyearofcurrentCO2cycle = lastyearoffesomcycle-numofCO2cycles*totnumyear+(currentCO2cycle-1)*totnumyear
    
     currentCO2year = firstyearofcurrentCO2cycle + (yearnew-firstyearoffesomcycle)+1
     if(mype==0) write(*,*),currentCO2year, firstyearofcurrentCO2cycle, yearnew, firstyearoffesomcycle
     write(currentCO2year_char,'(i4)') currentCO2year
     CO2vari     = 'AtmCO2_'//currentCO2year_char

     ! open file
     status=nf_open(CO2filename, nf_nowrite, ncid)
     if (status.ne.nf_noerr)then
        print*,'ERROR: CANNOT READ CO2 FILE CORRECTLY !!!!!'
        print*,'Error in opening netcdf file '//CO2filename
        stop
     call par_ex
     endif    
	
!	! data
     allocate(ncdata(12))
     status=nf_inq_varid(ncid, CO2vari, varid)
     CO2start = 1
     CO2count = 12
     status=nf_get_vara_double(ncid,varid,CO2start,CO2count,ncdata)
     AtmCO2(:)=ncdata(:)
     deallocate(ncdata)
     if (mype==0) write(*,*),'Current carbon year=',currentCO2year
     if (mype==0) write(*,*),'Atm CO2=', AtmCO2

    status=nf_close(ncid)

   endif  ! constant CO2 or changing
    
 ! Aeolian nitrogen deposition
    if (useAeolianN) then
      i=1 ! A single time entry
      DustNfilename = trim(REcoMDataPath)//'AeolianNitrogenDep.nc'
      if (yearnew .lt. 2010) then
         Nvari      = 'NDep'//cyearnew
      else
         Nvari      = 'NDep2009'
      endif

      if (mype==0) write(*,*) 'Updating Nitrogen deposition data for month ', i     
      call read_2ddata_on_grid_NetCDF(DustNfilename, Nvari, i, GloNDust, mesh) 
    else
      GloNDust(:) = 0.d0 ! no aeolian N input 
    end if
  else
    return
  end if
end subroutine Atm_input
!================================================================================
! Calculating second zooplankton respiration rates
!================================================================================
 subroutine krill_resp(n,mesh)
   use REcoM_declarations
   use REcoM_LocVar
   use REcoM_GloVar
   use g_clock
   use o_PARAM
   use g_PARSUP
   use mod_MESH
   use g_comm_auto
   implicit none
   integer                          :: n
   type(t_mesh), intent(in), target :: mesh  
#include  "../associate_mesh.h"
 
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



!  if ((Latr .lt. 0).and.(daynew .le. 105)) then
!       res_zoo2_a = 0.d0
!   else if((Latr .lt. 0).and.(105 .le. daynew).and.(daynew .le. 150)) then
!       res_zoo2_a = (-1./90.)*daynew +7./6.
!   else if((Latr .lt. 0).and.(150 .lt. daynew).and.(daynew .lt. 250)) then
!       res_zoo2_a = -0.5
!   else if((Latr .lt. 0).and.(250 .le. daynew).and.(daynew .le. 295)) then
!       res_zoo2_a = (1/90.)*daynew - 59./18.
!   else if((Latr .lt. 0).and.(daynew .gt. 295)) then
!       res_zoo2_a = 0.d0
!  end if

!!For Northern Hemisphere


!  if ((Latr .ge. 0).and.(daynew .le. 65)) then
!       res_zoo2_a = -0.5
!   else if((Latr .ge. 0).and.(285 .le. daynew).and.(daynew .le. 330)) then
!       res_zoo2_a = (-1./90.)*daynew +57./18.
!   else if((Latr .ge. 0).and.(330 .lt. daynew)) then
!       res_zoo2_a = -0.5
!   else if((Latr .ge. 0).and.(65 .le. daynew).and.(daynew .le. 110)) then
!       res_zoo2_a = (1/90.)*daynew - 22./18.
!   else if((Latr .ge. 0).and.(110 .lt. daynew).and.(daynew .lt. 285)) then
!       res_zoo2_a = 0.d0
!  end if

 
 end subroutine krill_resp

!================================================================================
! Subroutine controlling and calculating river deposition of Nutrients
!================================================================================

subroutine River_input(mesh)

  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  use REcoM_declarations
  use g_clock
  use g_read_other_NetCDF
  use g_PARSUP
  use o_PARAM, only : mstep, WP
  use mod_MESH
  use REcoM_ciso
  implicit none
  type(t_mesh), intent(in) , target :: mesh  
#include "netcdf.inc"

  real(kind=8), allocatable :: ncdata(:)
  integer	            :: status, ncid, varid
  character(300)            :: Riverfilename
  integer, dimension(2)     :: istart, icount
  integer                   :: CO2start, CO2count
  logical                   :: do_read=.false.
  integer                   :: i
  integer                   :: node_size, num_sec_in_month

#include "../associate_mesh.h" 

  node_size=myDim_nod2D+eDim_nod2D 

!-Switch river nutrients on-------------------------------------------------------

  if (useRivers) then

  ! River inputs are in mmol/m2/s

    ! add river nutrients as surface boundary condition (surface_bc function in oce_ale_tracers)
     is_riverinput = 1.0d0

     if (mstep == 1) then ! The year has changed
        if (mype==0) write(*,*), month
        i=month
!        num_sec_in_month = num_day_in_month(fleapyear,month)
!        num_sec_in_month=num_sec_in_month*86400


           Riverfilename = trim(REcoMDataPath)//'RiverineInput.nc'
           if (mype==0) write(*,*) 'Updating riverine restoring data for month ', i,' from ', Riverfilename     
           call read_2ddata_on_grid_NetCDF(Riverfilename, 'Alkalinity', i, RiverAlk2D, mesh)
!           write(*,*) mype, 'RiverAlk2D', maxval(RiverAlk2D(:)), minval(RiverAlk2D(:))        
!           molar convertion of [CaCo3] * 2  -> [total Alkalinity]   
           RiverAlk2D = RiverAlk2D * 2    

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DIC', i, RiverDIC2D, mesh) 
!           write(*,*) mype, 'RiverDIC2D', maxval(RiverDIC2D(:)), minval(RiverDIC2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DIN', i, RiverDIN2D, mesh) 
!           write(*,*) mype, 'RiverDIN2D', maxval(RiverDIN2D(:)), minval(RiverDIN2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DOC', i, RiverDOC2D, mesh) 
!           write(*,*) mype, 'RiverDOC2D', maxval(RiverDOC2D(:)), minval(RiverDOC2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DON', i, RiverDON2D, mesh) 
!           write(*,*) mype, 'RiverDON2D', maxval(RiverDON2D(:)), minval(RiverDON2D(:))     

!           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DSi', i, RiverDSi2D, mesh) 
!           write(*,*) mype, 'RiverDSi2D', maxval(RiverDSi2D(:)), minval(RiverDSi2D(:))  
            RiverDSi2D = RiverDIN2D * (16/15)   


     else

!-Checking if files need to be opened---------------------------------------------
        call monthly_event(do_read)
!        if (mype==0) write(*,*), do_read, month
        if(do_read) then ! file is opened and read every month
           i=month+1
           if (i > 12) i=1        
!           num_sec_in_month = num_day_in_month(fleapyear,month)
!           num_sec_in_month=num_sec_in_month*86400

           Riverfilename = trim(REcoMDataPath)//'RiverineInput.nc'
           if (mype==0) write(*,*) 'Updating riverine restoring data for month ', i,' from ', Riverfilename     
           call read_2ddata_on_grid_NetCDF(Riverfilename, 'Alkalinity', i, RiverAlk2D, mesh) 
!           write(*,*) mype, 'RiverAlk2D', maxval(RiverAlk2D(:)), minval(RiverAlk2D(:))   
!           molar convertion of [CaCo3] * 2  -> [total Alkalinity]   
           RiverAlk2D = RiverAlk2D * 2        

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DIC', i, RiverDIC2D, mesh) 
!           write(*,*) mype, 'RiverDIC2D', maxval(RiverDIC2D(:)), minval(RiverDIC2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DIN', i, RiverDIN2D, mesh) 
!           write(*,*) mype, 'RiverDIN2D', maxval(RiverDIN2D(:)), minval(RiverDIN2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DOC', i, RiverDOC2D, mesh) 
!           write(*,*) mype, 'RiverDOC2D', maxval(RiverDOC2D(:)), minval(RiverDOC2D(:))     

           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DON', i, RiverDON2D, mesh) 
!           write(*,*) mype, 'RiverDON2D', maxval(RiverDON2D(:)), minval(RiverDON2D(:))     

!           call read_2ddata_on_grid_NetCDF(Riverfilename, 'DSi', i, RiverDSi2D, mesh) 
!           write(*,*) mype, 'RiverDSi2D', maxval(RiverDSi2D(:)), minval(RiverDSi2D(:))  
            RiverDSi2D = RiverDIN2D * (16/15)
        end if
     end if

  else 

     is_riverinput = 0.0d0
     RiverAlk2D = 0.0d0
     RiverDIC2D = 0.0d0
     RiverDIN2D = 0.0d0
     RiverDOC2D = 0.0d0
     RiverDON2D = 0.0d0
     RiverDSi2D = 0.0d0

  end if 
end subroutine River_input

!================================================================================
! Subroutine controlling and calculating erosion input
!================================================================================

subroutine Erosion_input(mesh)

  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  use REcoM_declarations
  use g_clock
  use g_read_other_NetCDF
  use g_PARSUP
  use o_PARAM, only : mstep, WP
  use mod_MESH
  use REcoM_ciso
  implicit none
  type(t_mesh), intent(in) , target :: mesh  
#include "netcdf.inc"

  real(kind=8), allocatable :: ncdata(:)
  integer	            :: status, ncid, varid
  character(300)            :: Erosionfilename
  integer, dimension(2)     :: istart, icount
  integer                   :: CO2start, CO2count
  logical                   :: do_read=.false.
  integer                   :: i
  integer                   :: node_size, num_sec_in_month

#include "../associate_mesh.h" 

  node_size=myDim_nod2D+eDim_nod2D 
!-Switch Erosion on-------------------------------------------------------

  if (useErosion) then

  ! River inputs are in mmol/m2/s

     ! add erosion as surface boundary condition (surface_bc function in oce_ale_tracers)
     is_erosioninput = 1.0d0


     if (mstep == 1) then ! The year has changed
        if (mype==0) write(*,*), month
        i=month
!        num_sec_in_month = num_day_in_month(fleapyear,month)
!        num_sec_in_month=num_sec_in_month*86400

           Erosionfilename = trim(REcoMDataPath)//'ErosionInput.nc'
           if (mype==0) write(*,*) 'Updating erosion restoring data for month ', i,' from ', Erosionfilename     
           call read_2ddata_on_grid_NetCDF(Erosionfilename, 'POC', i, ErosionTOC2D, mesh) 
!           write(*,*) mype, 'ErosionTOC2D', maxval(ErosionTOC2D(:)), minval(ErosionTOC2D(:))        

           call read_2ddata_on_grid_NetCDF(Erosionfilename, 'PON', i, ErosionTON2D, mesh) 
!           write(*,*) mype, 'ErosionTON2D', maxval(ErosionTON2D(:)), minval(ErosionTON2D(:))

           ! No silicates in erosion, we convert from nitrogen with redfieldian ratio     
	   ErosionTSi2D=ErosionTON2D * 16/15
!           write(*,*) mype, 'ErosionTSi2D', maxval(ErosionTSi2D(:)), minval(ErosionTSi2D(:))        
     else

!-Checking if files need to be opened---------------------------------------------
        call monthly_event(do_read)
!        if (mype==0) write(*,*), do_read, month
        if(do_read) then ! file is opened and read every month
           i=month+1
           if (i > 12) i=1        
!           num_sec_in_month = num_day_in_month(fleapyear,month)
!           num_sec_in_month=num_sec_in_month*86400

           Erosionfilename = trim(REcoMDataPath)//'ErosionInput.nc'
           if (mype==0) write(*,*) 'Updating erosion restoring data for month ', i,' from ', Erosionfilename     
           call read_2ddata_on_grid_NetCDF(Erosionfilename, 'POC', i, ErosionTOC2D, mesh) 
!           write(*,*) mype, 'ErosionTOC2D', maxval(ErosionTOC2D(:)), minval(ErosionTOC2D(:))        

           call read_2ddata_on_grid_NetCDF(Erosionfilename, 'PON', i, ErosionTON2D, mesh) 
!           write(*,*) mype, 'ErosionTON2D', maxval(ErosionTON2D(:)), minval(ErosionTON2D(:))        

            ! No silicates in erosion, we convert from nitrogen with redfieldian ratio     
	    ErosionTSi2D=ErosionTON2D * 16/15 
!           write(*,*) mype, 'ErosionTSi2D', maxval(ErosionTSi2D(:)), minval(ErosionTSi2D(:))        
        end if
     end if
  else
     is_erosioninput = 0.0d0

     ErosionTOC2D = 0.0d0
     ErosionTON2D = 0.0d0
     ErosionTSi2D = 0.0d0
  end if 
end subroutine Erosion_input

