! Routines for initializing and updating ocean surface forcing fields
!-------------------------------------------------------------------------

subroutine init_atm_forcing
  ! read in forcing data that are constant in time
  ! the time varying forcing fields will be read in by read_new_atm_forcing
  ! assume atmosphere forcing data on T62 NCEP/NCAR grid
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??

  use o_PARAM
  use o_MESH
  use o_arrays
  use i_therm_param
  use i_arrays
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
  integer                   		:: itime, i, k, n2, year_first_rec
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
  if(rad_data_source=='CORE1' .or. rad_data_source=='CORE2') rad_ttp_ind=2
  if(precip_data_source=='CORE1' .or. precip_data_source=='CORE2') precip_ttp_ind=3
  if(runoff_data_source=='CORE1' .or. runoff_data_source=='CORE2') runoff_ttp_ind=0
  if(sss_data_source=='CORE1' .or. sss_data_source=='CORE2') sss_ttp_ind=4


  !==========================================================================
  ! runoff    

  if(runoff_data_source=='CORE1' .or. runoff_data_source=='CORE2' ) then

     ! runoff in CORE is constant in time

     ! Warning: For a global mesh, conservative scheme is to be updated!!

     file=trim(ForcingDataPath)//trim(runoff_data_source)//'/runoff.nc'
     vari='Foxx_o_roff'
     check_dummy=.false.

     itime=1
     call read_other_NetCDF(file, vari, itime, runoff, check_dummy) 
     runoff=runoff/1000.0  ! Kg/s/m2 --> m/s
  end if
  

  !==========================================================================
  ! sss restoring

  if(surf_relax_S>0.) then
     if(sss_data_source=='AAOMIP' .OR. sss_data_source=='ECHAM5') then

        ! taking the annual mean of PHC2 SSS

        file=trim(ForcingDataPath)//'CORE2'//'/PHC2_salx.nc'
        vari='SALT'
        check_dummy=.true.

        Ssurf=0.0

        allocate(aux(n2))
        do itime=1,12
           call read_other_NetCDF(file, vari, itime, aux, check_dummy) 
           Ssurf=Ssurf+aux
        end do
        Ssurf=Ssurf/12.0
        deallocate(aux)
     endif
  end if

  if(mype==0) write(*,*) 'Parts of forcing data (only constant in time fields) are read'

end subroutine init_atm_forcing
!
!------------------------------------------------------------------------------------------
!
subroutine update_atm_forcing(istep)
  ! update atmospheric forcing data
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !---------------------------------------------------------
  
  use o_PARAM
  use o_MESH
  use o_arrays
  use i_arrays
  use i_param
  use i_therm_param
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_parsup
  use g_clock
  use g_config
  use g_forcing_interp
  use g_comm_auto
  use g_read_CORE_NetCDF
  implicit none

  integer		:: i, istep,itime,n2,n,nz,k,elem
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy,tx,ty,tvol
  real              	:: t1, t2
  character(15)                         :: vari, filevari
  character(4)                          :: fileyear
  integer, parameter                    :: nci=192, ncj=94 ! T62 grid
  real(kind=8), dimension(nci,ncj)      :: array_nc, array_nc2,array_nc3,x
  character(80)                         :: file

  t1=MPI_Wtime()  
  ! first, read forcing data
  call read_new_atm_forcing
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
  do i=1,myDim_nod2d+eDim_nod2d     
     dux=u_wind(i)-u_w(i) 
     dvy=v_wind(i)-v_w(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
     stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
     dux=u_wind(i)-u_ice(i) 
     dvy=v_wind(i)-v_ice(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
     stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
  end do

stress_atmoce_x=stress_atmoce_x!*clim_growth
stress_atmoce_y=stress_atmoce_y!*clim_growth

  ! heat and fresh water fluxes are treated in i_therm and ice2ocean

  t2=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update of forcing data took', t2-t1
  end if

end subroutine update_atm_forcing
!
!------------------------------------------------------------------------------------
!
subroutine read_new_atm_forcing
  ! update atmospheric forcing, SSS, Chl etc. 
  ! assume forcing data on T62 NCEP/NCAR grid
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  
  use o_PARAM
  use o_MESH
  use o_arrays
  use i_therm_param
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
  use i_ARRAYS
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
  ! wind u and v, Tair, and shum              
  if(wind_data_source=='NCEP') then

     if(update_forcing_flag(wind_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind)

        ! three temporal types (6 hourly, daily and monthly) are possible 

        if(wind_ttp_ind==1) then ! 6 hourly data

           ! 10-m wind m/s ----------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'uwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'vwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'air.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)   
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'shum.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2)    

        elseif(wind_ttp_ind==2) then ! daily data      

           ! 10-m wind --------------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'uwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'vwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)  
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'air.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)   
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'shum.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

        elseif(wind_ttp_ind==3) then ! monthly data

           ! 10-m wind m/s ----------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'uwnd10m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'vwnd10m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'air2m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)  
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'shum2m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

        end if

     end if

  elseif(wind_data_source=='CORE2') then

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
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

        ! 10-m temperature -------------------------------------

        filevari='t_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='T_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2) 
        Tair=Tair-tmelt  ! Kelvin --> degree celcium

        ! 10 m specific humdity  Kg/Kg -------------------------

        filevari='q_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='Q_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

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
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

        ! 10-m temperature -------------------------------------

        filevari='t_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='T_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2) 
        Tair=Tair-tmelt  ! Kelvin --> Degree Celcius

        ! 10 m specific humdity  Kg/Kg -------------------------

        filevari='q_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='Q_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2)  

     end if

  endif


  !==========================================================================
  ! radiation 

  if(rad_data_source=='NCEP') then

     if(update_forcing_flag(rad_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)

        ! two temporal types (6 hourly, daily) are possible 

        if(rad_ttp_ind==1) then ! 6 hourly data

           ! short wave W/m2 --------------------------------------

           vari='dswrf'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'dswrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

           ! long wave W/m2 ---------------------------------------

           vari='dlwrf'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'dlwrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2)  

        elseif(rad_ttp_ind==2) then ! daily data

           ! short wave W/m2 --------------------------------------

           vari='dswrf'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'dswrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

           ! long wave W/m2 ---------------------------------------

           vari='dlwrf'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'dlwrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2) 

        end if

     end if

  elseif(rad_data_source=='CORE2') then

     ! in CORE daily radiation fluxes are used 

     if(update_forcing_flag(rad_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)

        ! short wave W/m2 --------------------------------------

        filevari='ncar_rad.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='SWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

        ! long wave W/m2 ---------------------------------------

        vari='LWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2) 

     end if

  elseif(rad_data_source=='CORE1') then

     ! in CORE daily radiation fluxes are used 

     if(update_forcing_flag(rad_ttp_ind)==1) then

        itime=forcing_rec(rad_ttp_ind)

        ! short wave W/m2 --------------------------------------

        filevari='ncar_rad'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='SWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

        ! long wave W/m2 ---------------------------------------

        vari='LWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2)  

     end if

  end if


  !==========================================================================
  ! precipitation

  if(precip_data_source=='NCEP') then

     if(update_forcing_flag(precip_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)

        ! four temporal types (6 hourly, daily and monthly, monthly ltm) are possible 

        if(precip_ttp_ind==1) then ! 6 hourly data

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'prate.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)  
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==2) then ! daily data      

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'prate.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==3) then ! monthly data 

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'prate.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==4) then ! monthly ltm data 

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_mon_ltm/'//'prate.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        end if

     end if

  elseif(precip_data_source=='CORE2') then

     ! in CORE monthly precipitation is used; 
     ! And rain and snow are separated.

     if(update_forcing_flag(precip_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)

        ! rain mm/s --------------------------------------------

        filevari='ncar_precip.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='RAIN'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
        prec_rain=prec_rain/1000.  ! mm/s --> m/s

        ! snow mm/s --------------------------------------------

        vari='SNOW'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_snow,n2)  
        prec_snow=prec_snow/1000.  ! mm/s --> m/s

     end if

  elseif(precip_data_source=='CORE1') then

     ! in CORE monthly precipitation is used; 
     ! And rain and snow are separated.

     if(update_forcing_flag(precip_ttp_ind)==1) then

        itime=forcing_rec(precip_ttp_ind)

        ! rain mm/s --------------------------------------------

        filevari='ncar_precip'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='RAIN'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)  
        prec_rain=prec_rain/1000.  ! mm/s --> m/s

        ! snow mm/s --------------------------------------------

        vari='SNOW'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_snow,n2)  
        prec_snow=prec_snow/1000.  ! mm/s --> m/s

     end if

  end if


  !==========================================================================
  ! runoff  

  if(runoff_data_source=='Dai09') then

     if(update_forcing_flag(runoff_ttp_ind)==1) then
        if(runoff_ttp_ind==4) then
           !climatology monthly mean

           itime=forcing_rec(runoff_ttp_ind)
           file=trim(MeshPath)//'runoff_on_grid//runoff_clim.nc' 
           vari='runoff'

           call read_2ddata_on_grid_NetCDF(file,vari,itime,runoff)

           !kg/m2/s -> m/s
           runoff=runoff/1000.

        elseif(runoff_ttp_ind==3) then
           !monthly data

           write(*,*) 'Monthly runoff need to be updated. Forced to stop.'
           call par_ex
           stop
        end if
     end if

  elseif(runoff_data_source=='AAOMIP') then

     ! runoff is monthly ltm in AOMIP/AAOMIP

     if(update_forcing_flag(runoff_ttp_ind)==1) then

        allocate(aux(nod2D))

        itime=forcing_rec(runoff_ttp_ind)

        file=trim(ForcingDataPath)//'AAOMIP'//'/river_runoff.dat'
        if(system==1) then
           readtype=2
        else
           readtype=8
        end if

        open(51,file=trim(file),form='unformatted', access='direct',recl=readtype*nod2d)
        read(51,rec=itime) aux                
        runoff=aux(myList_nod2D)        
        close(51)

        deallocate(aux)
     end if

  end if


  !==========================================================================
  ! sss restoring

  if(surf_relax_S > 0.) then
     if(sss_data_source=='CORE1' .or. sss_data_source=='CORE2') then

        ! sss is monthly ltm in CORE cases

        if(update_forcing_flag(sss_ttp_ind)==1) then

           itime=forcing_rec(sss_ttp_ind)

           file=trim(ForcingDataPath)//trim(sss_data_source)//'/PHC2_salx.nc'
           vari='SALT'
           check_dummy=.true.
           call read_other_NetCDF(file, vari, itime, Ssurf, check_dummy)  

        end if

     end if
  end if

end subroutine read_new_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_T62_wind(xarray, yarray)
  ! rotate wind on T62 grid from geographical coord. to rotated coordinates.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  
  use o_param
  use g_config
  use g_rotate_grid
  implicit none

  integer, parameter 	:: ni=192, nj=94  ! NCEP and CORE are on the same grid.
  integer               :: i, j
  real(kind=8)      	:: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj) 

  ! NCEP/CORE latitude
  cy=(/-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
       -77.1394, -75.2351, -73.3307, -71.4262, -69.5217, -67.6171, &  
       -65.7125, -63.8079, -61.9033, -59.9986, -58.0939, -56.1893, &  
       -54.2846, -52.3799, -50.4752, -48.5705, -46.6658, -44.7611,&  
       -42.8564, -40.9517, -39.0470, -37.1422, -35.2375, -33.3328, &  
       -31.4281, -29.5234, -27.6186, -25.7139, -23.8092, -21.9044, &  
       -19.9997, -18.0950, -16.1902, -14.2855, -12.3808, -10.47604, &  
       -8.57131, -6.66657, -4.76184, -2.8571, -0.952368, 0.952368, &  
       2.8571, 4.76184, 6.66657, 8.57131, 10.47604, 12.3808, &  
       14.2855, 16.1902, 18.095, 19.9997, 21.9044, 23.8092, &  
       25.7139, 27.6186, 29.5234, 31.4281, 33.3328, 35.2375,&  
       37.1422, 39.047,  40.9517, 42.8564, 44.7611, 46.6658,&  
       48.5705, 50.4752, 52.3799, 54.2846, 56.1893, 58.0939,&  
       59.9986, 61.9033, 63.8079, 65.7125, 67.6171, 69.5217, &  
       71.4262, 73.3307, 75.2351, 77.1394, 79.0435, 80.9473, &  
       82.8508, 84.7532, 86.6531, 88.542 /)*rad

  ! NCEP/CORE longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+1.875*rad
  enddo

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_T62_wind
!
!-----------------------------------------------------------------------------------------
!
subroutine shortwave_radiation(day,nhi,nhf,n2)
  !
  ! calculates the incoming shortwave radiation following
  ! Parkinson and Washington, 1979
  !
  ! input data:
  ! ndoyr               day of year  (1..365 or 366)
  ! nhi                 initial hour (0...23)
  ! nhf                 final hour   (0...23)
  !                     (last hour to be completely considered)
  !
  ! For instance, for a 6h-timestep:     nhi    nhf
  !                                      0      5
  !                                      6      11
  !                                      12     17
  !                                      18     23
  !
  !            for a daily timestep:     0      23
  !
  !
  ! Derived from subroutine swr of BRIOS sea ice model by Ralph Timmermann, 23.9.2004
  ! Copied to the new FESOM version by Qiang Wang
  ! Reviewed by ??
  !=========================================================================
  USE o_MESH, only: coord_nod2d, geo_coord_nod2D
  USE o_PARAM
  USE g_forcing_arrays
  
  IMPLICIT NONE

  INTEGER                        :: day, ndoyrd, nhi, nhf
  INTEGER                        :: it, iday, itd, n2
  real(kind=WP), dimension(0:23) :: cosha
  real(kind=WP), dimension(n2)   :: sinlat, coslat, cosz
  real(kind=WP), dimension(n2)   :: cf, evap, fsh
  real(kind=WP)                  :: dec, sindec, cosdec


  ! ----------------------------------------------------------------------
  ! Cosine of hour angle for every hour
  ! [Parkinson and Washington, 1979, eq. (3)]
  ! ----------------------------------------------------------------------
  do it=0,23
     cosha(it)=cos(2*pi*(12-it)/24.)             ! correct for GM
  enddo

  ! ----------------------------------------------------------------------
  ! Declination of the sun for every day
  ! -----------------------------------------------------------------------
  ndoyrd=min(day,365)    !no fleap year
  dec=rad*23.44*cos((174.-(ndoyrd-1))*2*pi/365.) ! 365 days/year
  sindec=sin(dec)
  cosdec=cos(dec)

  !-----------------------------------------------------------------------
  ! Sine and cosine of geographical latitude
  !-----------------------------------------------------------------------

  coslat=cos(geo_coord_nod2d(2,:))
  sinlat=sin(geo_coord_nod2d(2,:))

  !-----------------------------------------------------------------------
  !       fsh             shortwave radiation
  !       cf              cloud factor
  !-----------------------------------------------------------------------
  fsh  = 0.
  evap = shum*Pair*0.01/0.622*1.e-3  ! Pair with unit Pa, *0.01->hPa  
  !ecmwf e_vapor(i,j)=611.e-5*exp(19.*tdew(i,j)/(tdew(i,j)+250.))
  cf = 1.-0.6*cloudiness**3


  !-----------------------------------------------------------------------
  !LOOP for every hour
  !-----------------------------------------------------------------------
  do it=nhi,nhf
     if (it.lt.0) then
        itd=it+24
     else
        itd=it
     endif
     !-----------------------------------------------------------------------
     ! Cosine of zenith distance
     ! Parkinson and Washington, 1979, eq. (2)
     !-----------------------------------------------------------------------
     cosz=sinlat*sindec+coslat*cosdec*cosha(itd)
     !-----------------------------------------------------------------------
     ! At night, when no sun is present and cosz < 0.,
     ! the incoming solar radiation is zero (and not negative).
     !-----------------------------------------------------------------------
     cosz=max(cosz,0.)
     !-----------------------------------------------------------------------
     ! Add up incoming shortwave radiation for each and every hour
     ! Parkinson and Washington, 1979, eq. (1)
     !
     ! Solar constant = 1353 [W/m**2]
     !         56.375 = 1353 / 24 for contributions from 24 hours
     !         1353/float(nhf-nhi+1) for contributions from nhf-nhi hours
     !-----------------------------------------------------------------------
     fsh=fsh+1353./float(nhf-nhi+1)*cosz*cosz*cf(:) &
          /((cosz+2.7)*evap+1.085*cosz+0.1)
     !-----------------------------------------------------------------------
  enddo

  shortwave=fsh

end subroutine shortwave_radiation


