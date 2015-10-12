! Routines for initializing and updating the ocean surface forcing fields
! for the case without the ice module.

subroutine init_atm_forcing_OnlyOcean
  ! initialize the atmospheric forcing data for the ocean-alone model
  ! assume forcing data on T62 NCEP/NCAR grid
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_PARAM
  use o_MESH
  use o_arrays
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parsup
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, i, k, n2
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:) 

  n2=myDim_nod2D+eDim_nod2D       

  ! predefinition/correction
  ! for the CORE case:
  if(wind_data_source=='CORE1' .or. wind_data_source=='CORE2') wind_ttp_ind=1

  if(mype==0) write(*,*) 'Forcing data which are constant in time are initialized'

end subroutine init_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------ 
!
subroutine update_atm_forcing_OnlyOcean(istep)
  ! update atmospheric forcing data for ocean-alone cases
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_PARAM
  use o_MESH
  use o_arrays
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_parsup
  use g_clock
  use g_config
  implicit none

  integer		:: i, istep,elnodes(3)
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy
  real(kind=8)          :: rhoair_local
  real              	:: t1, t2

  data rhoair_local /1.3/

  t1=MPI_Wtime()  
  ! first, read forcing data
  call read_new_atm_forcing_OnlyOcean
  ! second, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call cal_wind_drag_coeff
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
     call ncar_ocean_fluxes_mode
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if

  ! third, compute wind stress
  do i=1,myDim_elem2D
     elnodes=elem2D_nodes(:,i)
     dux=sum(u_wind(elnodes))/3.0_WP-UV(1,1,i)
     dvy=sum(v_wind(elnodes))/3.0_WP-UV(2,1,i)
     !dux=sum(u_wind(elnodes))/3.0_WP
     !dvy=sum(v_wind(elnodes))/3.0_WP
     aux=sqrt(dux**2+dvy**2)*rhoair_local*sum(Cd_atm_oce_arr(elnodes))/3.0_WP 
     stress_surf(1,i) = aux*dux
     stress_surf(2,i) = aux*dvy
    
    
  end do

  t2=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if

end subroutine update_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------
!
subroutine read_new_atm_forcing_OnlyOcean
  ! read the second record of atmospheric forcing data for the ocean alone cases 
  ! assume forcing data on T62 NCEP/NCAR grid
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_PARAM
  use o_MESH
  use o_arrays
   use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parsup
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, m, i, k, n2
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:)       

  n2=myDim_nod2D+eDim_nod2D 

  !==========================================================================
  ! wind u and v
  if(wind_data_source=='CORE2') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind) 
       
        ! 10-m wind m/s ----------------------------------------

        filevari='u_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        !if(rotated_grid) 
	call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

     end if

  elseif(wind_data_source=='CORE1') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        itime=forcing_rec(wind_ttp_ind)
        
        ! 10-m wind m/s ----------------------------------------

        filevari='u_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        !if(rotated_grid) 
	call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

     end if
  endif

  !==========================================================================
  ! climatology of ssT/ssS
  ! in this case we use WOA05 (NetCDF files)

  if(update_forcing_flag(sss_ttp_ind)==1) then

     itime=forcing_rec(sss_ttp_ind)
     
     filevari='t0112an1'
     file=trim(ClimateDataPath)//trim(filevari)//'.nc'
     vari='t0112an1'
     call read_surf_hydrography_NetCDF(file, vari, itime, Tsurf)

     filevari='s0112an1'
     file=trim(ClimateDataPath)//trim(filevari)//'.nc'
     vari='s0112an1'
     call read_surf_hydrography_NetCDF(file, vari, itime, Ssurf)
  end if

end subroutine read_new_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------
!
