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

  real(kind=8),dimension(mesh%nl,5),intent(out)    :: wF             ! [m/day] Velocities of fluxes at the border of the control volumes
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

!! Background sinking speed

  wF(2:Nn,ivphy)   = VPhy  
  wF(2:Nn,ivdia)   = VDia
  wF(2:Nn,ivdet)   = VDet
  wF(2:Nn,ivdetsc) = VDet_zoo2
  wF(2:Nn,ivcoc)   = VCocco            ! NEW

  wF(1,:)          = 0.d0
  wF(Nn+1,:)       = 0.d0

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

    if (use_atbox) then  
!     Atmospheric box model CO2 values
      AtmCO2(:)                   = x_co2atm(1)
      if (ciso) then 
        AtmCO2_13(:)              = x_co2atm_13(1)
        if (ciso_14) AtmCO2_14(:,1) = x_co2atm_14(1)
      end if
    else 
!     Prescribed atmospheric CO2 values
    if (constant_CO2) then
      AtmCO2(:) = CO2_for_spinup
      !if (mype==0) write(*,*),'in atm_input: Atm CO2=',AtmCO2
        if (ciso) then
          AtmCO2_13          = CO2_for_spinup * (1. + 0.001 * delta_co2_13)
          if (ciso_14) then
!           Atmospheric 14C varies with latitude
            do i=1, myDim_nod2D
!             Latitude of atmospheric input data
              lat_val = geo_coord_nod2D(2,i) / rad
!             Binning to latitude zones
              if (ciso_organic_14) then
!               Convert Delta_14C to delta_14C
                delta_co2_14 = (big_delta_co2_14(lat_zone(lat_val)) + 2. * delta_co2_13 + 50.) / &
                               (0.95 - 0.002 * delta_co2_13)
              else
!               "Inorganic" 14C approximation: delta_14C := Delta_14C 
                delta_co2_14 = big_delta_co2_14(lat_zone(lat_val))
              end if
              AtmCO2_14(lat_zone(lat_val),:) = CO2_for_spinup * (1. + 0.001 * delta_co2_14)
            end do
          end if
        end if

    else  

     CO2filename = trim(REcoMDataPath)//'MonthlyAtmCO2_gcb2023.nc'

     totnumyear                 = lastyearoffesomcycle-firstyearoffesomcycle+1
     firstyearofcurrentCO2cycle = lastyearoffesomcycle-numofCO2cycles*totnumyear+(currentCO2cycle-1)*totnumyear
    
     currentCO2year = firstyearofcurrentCO2cycle + (yearnew-firstyearoffesomcycle)+1
     if(mype==0) write(*,*),currentCO2year, firstyearofcurrentCO2cycle, yearnew, firstyearoffesomcycle
     write(currentCO2year_char,'(i4)') currentCO2year
     CO2vari     = 'AtmCO2_'//currentCO2year_char

! fesom2.1_recom:
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

! fesom-2.1-recom-paleo:
!        call read_2ddata_on_grid_NetCDF(CO2filename, CO2vari, i, AtmCO2, mesh) 
      end if ! constant or file
    end if   ! atmospheric box model or prescribed CO2 values   

!   Control output of atmospheric CO2 values
    if (mype==0 .and. my_fesom_group == 0) then !OG
      print *,                "In atm_input: AtmCO2    = ", AtmCO2(1)
      if (ciso) then
        print *,              "              AtmCO2_13 = ", AtmCO2_13(1)
        if (ciso_14) print *, "              AtmCO2_14 = ", AtmCO2_14(:,1)
      end if
      if (use_atbox) print *, "              use_atbox = .true."
    end if
 
 ! Aeolian nitrogen deposition
    if (useAeolianN) then
      i=1 ! A single time entry
      DustNfilename = trim(REcoMDataPath)//'AeolianNitrogenDep.nc'
      if (yearnew .gt. 2009) then
         Nvari      = 'NDep2009'
      else if (yearnew .lt. 1850) then
         Nvari      = 'NDep1850'
      else
         Nvari      =  'NDep'//cyearnew
      endif

      if (mype==0) write(*,*) 'Updating Nitrogen deposition data for month ', i     
      call read_2ddata_on_grid_NetCDF(DustNfilename, Nvari, i, GloNDust, mesh) 
    else
      GloNDust(:) = 0.d0 ! no aeolian N input 
    end if
  else
    return
  end if ! mstep == 1 for CO2
end subroutine Atm_input

!================================================================================
! Subroutine controlling and calculating sed. input of DIC, Alk, O2 and
! nutrients 
!================================================================================
subroutine Sed_input(mesh)

  use REcoM_GloVar
  use REcoM_locVar
  use recom_config
  use g_clock
  use g_read_other_NetCDF
  use g_PARSUP
  use o_PARAM, only : mstep, WP
  use mod_MESH
  use o_mesh
  use REcoM_ciso
  use REcoM_declarations
  USE g_forcing_arrays,    only: runoff

  implicit none
  type(t_mesh), intent(in) , target :: mesh
#include "netcdf.inc"

  character(300)            :: sedfilename
  integer                   :: i
!  integer                   :: num_sec_in_month
  logical                   :: do_read=.false.
  real(kind=8)              :: ncdata(9)
  integer                   :: status, ncid, varid, n_lb
  integer, dimension(2)     :: istart, icount
  real(kind=8)              :: total_runoff
#include "../associate_mesh.h"

!-Checking if files need to be opened---------------------------------------------

  if(use_MEDUSA .and. (sedflx_num .ne. 0)) then
! MEDUSA input needs to be renamed via jobscript
!   sedfilename  = trim(ResultPath)//'medusa_flux2fesom.'//cyearold//'.nc'
   sedfilename  = trim(ResultPath)//'medusa_flux2fesom.nc'

!   num_sec_in_month = num_day_in_month(fleapyear,month)
!   num_sec_in_month=num_sec_in_month*86400

! medusa output contains one annual mean of several diffusitive fluxes
! and should be read once in the beginning of the run and then as annual_event

   if (mstep == 1) then ! The year has changed
      if (mype==0 .and. my_fesom_group == 0) write(*,*) 'Updating sedimentary input first time from', sedfilename !OG

!-Opening files--------------------------------------------------------------------

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_din', 1, GloSed(:,1), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DIN flux:', maxval(GloSed(:,1)), minval(GloSed(:,1))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic', 1, GloSed(:,2), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DIC flux:', maxval(GloSed(:,2)), minval(GloSed(:,2))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_alk', 1, GloSed(:,3), mesh)
!      if (mype==0) write(*,*) mype, 'sediment Alk flux:', maxval(GloSed(:,3)), minval(GloSed(:,3))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dsi', 1, GloSed(:,4), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DSi flux:', maxval(GloSed(:,4)), minval(GloSed(:,4))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_o2', 1, GloSed(:,5), mesh)
!      if (mype==0) write(*,*) mype, 'sediment O2 flux:', maxval(GloSed(:,5)), minval(GloSed(:,5))

      if(ciso) then
        call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic13', 1, GloSed(:,6), mesh)
!        if (mype==0) write(*,*) mype, 'sediment DIC13 flux:', maxval(GloSed(:,6)), minval(GloSed(:,6))
        if(ciso_14) then
          call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic14', 1, GloSed(:,7), mesh)
!        if (mype==0) write(*,*) mype, 'sediment DIC14 flux:', maxval(GloSed(:,7)), minval(GloSed(:,7))
        end if ! ciso_14
      end if ! ciso

! unit conversion
      GloSed(:,:)=GloSed(:,:)/86400

! read loopback fluxes from the same file
      if(add_loopback) then
        if (mype==0 .and. my_fesom_group == 0) write(*,*) 'adding loopback fluxes through runoff for the first time' !OG

        istart = (/1,1/)
        icount = (/1,1/)
        ncdata = 0.d0

!        total_runoff = sum(runoff*area(1,:))*86400
!       Does 'area' only contain values on one node? sum of area not equal total
!       ocean surface area!
        total_runoff = 8.76d5*86400
!        if (mype==0) write(*,*) mype, 'total runoff (m3/day):', total_runoff
!        if (mype==0) write(*,*) mype, 'runoff:', maxval(runoff),minval(runoff)
!        if (mype==0) write(*,*) mype, 'surface area',maxval(area(1,:)),minval(area(1,:))

        status=nf_open(sedfilename, nf_nowrite, ncid)
        if(status.ne.nf_noerr) call handle_err(status)

        status=nf_inq_varid(ncid, 'loopback_orgm_din', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(1))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_din (mmolN/day):', ncdata(1) !OG

        status=nf_inq_varid(ncid, 'loopback_orgm_dic', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(2))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_dic (mmolC/day):', ncdata(2) !OG

        status=nf_inq_varid(ncid, 'loopback_orgm_alk', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(3))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_alk (mmolAlk/day):', ncdata(3) !OG

        status=nf_inq_varid(ncid, 'loopback_opal', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(4))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_opal (mmolSi/day):', ncdata(4) !OG

        status=nf_inq_varid(ncid, 'loopback_caco3', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(5))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco3 (mmolC/day):', ncdata(5) !OG

      if(ciso) then
        status=nf_inq_varid(ncid, 'loopback_orgm_dic13', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(6))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_dic13:', ncdata(6)      !OG   

        status=nf_inq_varid(ncid, 'loopback_caco313', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(7))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco313:', ncdata(7)!OG

       if(ciso_14 .and. ciso_organic_14) then
        status=nf_inq_varid(ncid, 'loopback_orgm_dic14', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(8))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_dic14:', ncdata(8) !OG

        status=nf_inq_varid(ncid, 'loopback_caco314', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(9))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco314:', ncdata(9) !OG

       end if ! ciso_14 .and. ciso_organic_14
      end if ! ciso
        status=nf_close(ncid)

! calculating fluxes back to ocean surface through rivers (mmol/m2/s)
! converting from fluxes out of sediment to fluxes into the ocean 
        do n_lb = 1,9
           lb_flux(:,n_lb) = -runoff*ncdata(n_lb)/total_runoff*lb_tscale
        end do

!        if (mype==0) write(*,*) mype, 'sum of surface area (m2)',
!        sum(area(1,:))
!        if (mype==0) write(*,*) mype, 'total ocean area (m2)', ocean_area        
!        if (mype==0) write(*,*) mype, 'DSi concentration in rivers',ncdata(4)/total_runoff
!        if (mype==0) write(*,*) mype, 'DIC concentration in rivers',ncdata(2)/total_runoff
!        if (mype==0) write(*,*) mype, 'Alk concentration in rivers',ncdata(3)/total_runoff
!        if (mype==0) write(*,*) mype, 'DIN concentration in rivers',ncdata(1)/total_runoff
!        if (mype==0) write(*,*) mype, 'DIN surface input:',minval(lb_flux(:,1)),maxval(lb_flux(:,1))
!        if (mype==0) write(*,*) mype, 'DIC surface input:',minval(lb_flux(:,2)),maxval(lb_flux(:,2))
!        if (mype==0) write(*,*) mype, 'Alk surface input:',minval(lb_flux(:,3)),maxval(lb_flux(:,3))
!        if (mype==0) write(*,*) mype, 'DSi surface input:',minval(lb_flux(:,4)),maxval(lb_flux(:,4))
!        if (mype==0) write(*,*) mype, 'DIC(calcite) surface input:',minval(lb_flux(:,5)),maxval(lb_flux(:,5))

      end if ! add_loopback

   else

!-Checking if files need to be opened---------------------------------------------
     call monthly_event(do_read)
     if(do_read) then ! file is opened and read every year
      i=month
      if (i > 12) i=1
      if (mype==0 .and. my_fesom_group == 0) write(*,*) 'Updating sedimentary input for month', i, 'from', sedfilename !OG

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_din', 1, GloSed(:,1), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DIN flux:', maxval(GloSed(:,1)), minval(GloSed(:,1))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic', 1, GloSed(:,2), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DIC flux:', maxval(GloSed(:,2)), minval(GloSed(:,2))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_alk', 1, GloSed(:,3), mesh)
!      if (mype==0) write(*,*) mype, 'sediment Alk flux:', maxval(GloSed(:,3)), minval(GloSed(:,3))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dsi', 1, GloSed(:,4), mesh)
!      if (mype==0) write(*,*) mype, 'sediment DSi flux:', maxval(GloSed(:,4)), minval(GloSed(:,4))

      call read_2ddata_on_grid_NetCDF(sedfilename, 'df_o2', 1, GloSed(:,5), mesh)
!      if (mype==0) write(*,*) mype, 'sediment O2 flux:', maxval(GloSed(:,5)), minval(GloSed(:,5))

      if(ciso) then
        call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic13', 1, GloSed(:,6), mesh)
!        if (mype==0) write(*,*) mype, 'sediment DIC13 flux:', maxval(GloSed(:,6)), minval(GloSed(:,6))
        if(ciso_14) then
          call read_2ddata_on_grid_NetCDF(sedfilename, 'df_dic14', 1, GloSed(:,7), mesh)
!          if (mype==0) write(*,*) mype, 'sediment DIC14 flux:', maxval(GloSed(:,7)), minval(GloSed(:,7))
        end if ! ciso_14
      end if ! ciso

!to mmol/m2/s
      GloSed(:,:)=GloSed(:,:)/86400

! read loopback fluxes from the same file
      if(add_loopback) then
        if (mype==0 .and. my_fesom_group == 0) write(*,*) 'adding loopback fluxes into the ocean monthly' !OG

        istart = (/1,1/)
        icount = (/1,1/)
        ncdata = 0.d0

        total_runoff = 8.76d5*86400

        status=nf_open(sedfilename, nf_nowrite, ncid)
        if(status.ne.nf_noerr) call handle_err(status)

        status=nf_inq_varid(ncid, 'loopback_orgm_din', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(1))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_din (mmolN/day):', ncdata(1) !OG

        status=nf_inq_varid(ncid, 'loopback_orgm_dic', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(2))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_dic (mmolC/day):', ncdata(2) !OG

        status=nf_inq_varid(ncid, 'loopback_orgm_alk', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(3))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_orgm_alk (mmolAlk/day):', ncdata(3) !OG

        status=nf_inq_varid(ncid, 'loopback_opal', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(4))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_opal (mmolSi/day):', ncdata(4) !OG

        status=nf_inq_varid(ncid, 'loopback_caco3', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(5))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco3 (mmolC/day):', ncdata(5) !OG

      if(ciso) then
        status=nf_inq_varid(ncid, 'loopback_orgm_dic13', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(6))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_dic13:', ncdata(6)     !OG   

        status=nf_inq_varid(ncid, 'loopback_caco313', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(7))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco313:', ncdata(7) !OG

       if(ciso_14 .and. ciso_organic_14) then
        status=nf_inq_varid(ncid, 'loopback_orgm_dic14', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(8))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_dic14:', ncdata(8) !OG

        status=nf_inq_varid(ncid, 'loopback_caco314', varid)
        if(status.ne.nf_noerr) call handle_err(status)
        status=nf_get_vara_double(ncid,varid,istart,icount,ncdata(9))
        if(status.ne.nf_noerr) call handle_err(status)
        if (mype==0 .and. my_fesom_group == 0) write(*,*) mype, 'loopback_caco314:', ncdata(9) !OG

       end if ! ciso_14 .and. ciso_organic_14
      end if ! ciso
        status=nf_close(ncid)

! calculating fluxes back to ocean surface through rivers (mmol/m2/s)
! converting from fluxes out of sediment to fluxes into the ocean 
        do n_lb = 1,9
           lb_flux(:,n_lb) = -runoff*ncdata(n_lb)/total_runoff*lb_tscale
        end do

      end if ! add_loopback

    end if ! do_read
   end if ! mstep==1
  else
    if (mype==0 .and. my_fesom_group == 0) write(*,*) 'sedimentary input from MEDUSA not used!' !OG
  endif ! use_MEDUSA and sedflx_num not 0

end subroutine Sed_input

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
  USE g_forcing_arrays,    only: runoff
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

  if (useRivFe) then

  ! River runoff is converted in gen_surface_forcing.F90 into m/s
  ! here multiplied with Fe concentration * muemolFe/m3 -> muemolFe/m2/s

  ! add river nutrients as surface boundary condition (surface_bc function in
  ! oce_ale_tracers)
     RiverFe = runoff * RiverFeConc
  else
     RiverFe = 0.0d0
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


