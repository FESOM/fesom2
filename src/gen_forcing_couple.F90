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
  character(500)             		:: file
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
#if defined (__oasis)
  use cpl_config
  use cpl_driver
  use g_rotate_grid
#endif

  implicit none

  integer		:: i, istep,itime,n2,n,nz,k,elem
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy,tx,ty,tvol
  real              	:: t1, t2
#ifdef __oasis
  real(kind=8)        				  :: flux_global(2), flux_local(2), eff_vol(2)
  real(kind=8), dimension(:), allocatable , save  :: exchange
  real(kind=8), dimension(:), allocatable , save  :: mask !, weight   
  logical, save                                   :: firstcall=.true.
  logical                                         :: action  
  logical                                         :: do_rotate_oce_wind=.false.
  logical                                         :: do_rotate_ice_wind=.false.    
  INTEGER                                         :: my_global_rank, ierror  
  INTEGER 					  :: status(MPI_STATUS_SIZE)
#endif
  character(15)                         :: vari, filevari
  character(4)                          :: fileyear
  integer, parameter                    :: nci=192, ncj=94 ! T62 grid
  real(kind=8), dimension(nci,ncj)      :: array_nc, array_nc2,array_nc3,x
  character(500)                        :: file

  t1=MPI_Wtime()
#ifdef __oasis
     if (firstcall) then
     allocate(exchange(myDim_nod2D+eDim_nod2D), mask(myDim_nod2D+eDim_nod2D))!, weight(myDim_nod2D+eDim_nod2D))
     allocate(a2o_fcorr_stat(nrecv,6))
     a2o_fcorr_stat=0.
     exchange=0.
     mask=0.
     !weight=0.     
     firstcall=.false.     
     end if
     do i=1,nsend
         exchange =0.0
         if (i.eq.1) then
            do n=1,myDim_nod2D+eDim_nod2D
               exchange(n)=tr_arr(1, n, 1)+tmelt                          ! sea surface temperature deg C
            end do
            elseif (i.eq.2) then
            exchange(:) = m_ice(:)                                  ! ice thickness [m]
            elseif (i.eq.3) then
            exchange(:) = a_ice(:)                                  ! ice concentation [%]
            elseif (i.eq.4) then
            exchange(:) = m_snow(:)                                 ! snow thickness
            elseif (i.eq.5) then
            exchange(:) = ice_temp(:)                               ! ice temperature
            elseif (i.eq.6) then
            exchange(:) = ice_alb(:)                                ! ice albedo
            else	    
            print *, 'not installed yet or error in cpl_oasis3mct_send', mype
         endif
         call cpl_oasis3mct_send(i, exchange, action)
      enddo
#ifdef VERBOSE
      do i=1, nsend 
        if (mype==0) write(*,*) 'SEND: field ', i, ' max val:', maxval(exchange), ' . ACTION? ', action 
      enddo
#endif
      mask=1.
      do i=1,nrecv
         exchange =0.0
#ifdef VERBOSE
	 if (mype==0) write(*,*) 'Trying to RECV: flux ', i 	  
#endif
         call cpl_oasis3mct_recv (i,exchange,action)
	 !if (.not. action) cycle
	 !Do not apply a correction at first time step!
	 if (i==1 .and. action .and. istep/=1) call net_rec_from_atm(action)
         if (i.eq.1) then
     	     if (.not. action) cycle
             stress_atmoce_x(:) =  exchange(:)                    ! taux_oce
	     do_rotate_oce_wind=.true.
         elseif (i.eq.2) then
     	     if (.not. action) cycle	 
             stress_atmoce_y(:) =  exchange(:)                    ! tauy_oce
	     do_rotate_oce_wind=.true.
         elseif (i.eq.3) then
     	     if (.not. action) cycle	 
             stress_atmice_x(:) =  exchange(:)                    ! taux_ice
	     do_rotate_ice_wind=.true.
         elseif (i.eq.4) then
     	     if (.not. action) cycle	 
             stress_atmice_y(:) =  exchange(:)                    ! tauy_ice
	     do_rotate_ice_wind=.true.	     
         elseif (i.eq.5) then
             if (action) then 
	        prec_rain(:)    =  exchange(:)	                  ! tot_prec
		mask=1.
		call force_flux_consv(prec_rain, mask, i, 0,action)
	     end if
         elseif (i.eq.6) then 
	     if (action) then
	        prec_snow(:)    =  exchange(:)                    ! snowfall
		mask=1.
		call force_flux_consv(prec_snow, mask,i,1,action) ! Northern hemisphere
		call force_flux_consv(prec_snow, mask,i,2,action) ! Southern Hemisphere
             end if
         elseif (i.eq.7) then
             if (action) then
	     evap_no_ifrac(:)     =  exchange(:)        	  ! tot_evap
	     tmp_evap_no_ifrac(:) =  exchange(:) 		  ! to reset for flux 
	     							  ! correction
	     end if
 	     mask=1.-a_ice
	     evap_no_ifrac(:)     =  tmp_evap_no_ifrac(:)
	     call force_flux_consv(evap_no_ifrac,mask,i,0,action)
         elseif (i.eq.8) then
             if (action) then
	     sublimation(:)       =  exchange(:)        	  ! tot_subl
	     tmp_sublimation(:)   =  exchange(:) 		  ! to reset for flux 
	     							  ! correction
	     end if
	     mask=a_ice 
	     sublimation(:)       =  tmp_sublimation(:)
	     call force_flux_consv(sublimation,mask,i,1,action) ! Northern hemisphere
	     call force_flux_consv(sublimation,mask,i,2,action) ! Southern Hemisphere
         elseif (i.eq.9) then
             if (action) then
	     oce_heat_flux(:)     =  exchange(:)        	  ! heat_oce
	     tmp_oce_heat_flux(:) =  exchange(:) 		  ! to reset for flux 
	     							  ! correction
	     end if
	     mask=1.-a_ice
	     oce_heat_flux(:)     =  tmp_oce_heat_flux(:)
	     call force_flux_consv(oce_heat_flux, mask, i, 0,action)
         elseif (i.eq.10) then
             if (action) then
	     ice_heat_flux(:)     =  exchange(:)        	  ! heat_ice
	     tmp_ice_heat_flux(:) =  exchange(:) 		  ! to reset for flux 
	     							  ! correction
	     end if
	     mask=a_ice
	     ice_heat_flux(:)     =  tmp_ice_heat_flux(:)
	     call force_flux_consv(ice_heat_flux, mask, i, 1,action) ! Northern hemisphere
	     call force_flux_consv(ice_heat_flux, mask, i, 2,action) ! Southern Hemisphere	     
         elseif (i.eq.11) then
             if (action) then
	     shortwave(:)         =  exchange(:)		  ! heat_swr
	     tmp_shortwave(:)     =  exchange(:) 		  ! to reset for flux 
	     							  ! correction
	     end if
	     mask=1.-a_ice
	     shortwave(:)   =  tmp_shortwave(:)
	     call force_flux_consv(shortwave, mask, i, 0,action)
         elseif (i.eq.12) then
             if (action) then
	        runoff(:)                   =  exchange(:)        ! runoff + calving
    	        mask=1.
		call force_flux_consv(runoff, mask, i, 0,action)
		!runoff=0.
             end if
	  end if  	  
!#ifdef VERBOSE
	  if (mype==0) then
		write(*,*) 'RECV: flux ', i, ', max val: ', maxval(stress_atmoce_x)
	  end if
!#endif
      end do

      if ((do_rotate_oce_wind .AND. do_rotate_ice_wind) .AND. rotated_grid) then
         do n=1, myDim_nod2D+eDim_nod2D
	    call vector_g2r(stress_atmoce_x(n), stress_atmoce_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
	    call vector_g2r(stress_atmice_x(n), stress_atmice_y(n), coord_nod2D(1, n), coord_nod2D(2, n), 0)
	 end do
	 do_rotate_oce_wind=.false.
         do_rotate_ice_wind=.false.
      end if
#else	
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

  ! heat and fresh water fluxes are treated in i_therm and ice2ocean
#endif /* (__oasis) */

  t2=MPI_Wtime()

#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
#ifdef __oasis     
     write(*,*) 'oasis part + exchange mesh:', time_recv(1)+time_send(1)
     write(*,*) 'exchange mesh             :', time_recv(2)+time_send(2)
#endif /* (__oasis) */     
  end if
#endif
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
  character(500)            		:: file
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
        if(system_arch==1) then
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
#if defined (__oasis)
!
!=================================================================
!
! Modifies the flux 'field2d' in order to conserve 
! the net fluxes.
! We distinguish between NH and SH fluxes (integer 'h') 
! where possible.
!
!================================================================= 
! History :
!  07-11  (D.Sidorenko,	AWI Germany) first routine
!  10-12  (T.Rackow, 	AWI Germany) code reordering and cleanup  
!-----------------------------------------------------------------
!
SUBROUTINE force_flux_consv(field2d, mask, n, h, do_stats)

  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
				flux_correction_north, flux_correction_south,	&
				flux_correction_total
  use g_parsup,	         only : myDim_nod2D, eDim_nod2D, mype
  use o_mesh,		 only :	geo_coord_nod2D
  use cpl_driver,	 only : nrecv, cpl_recv
  use cpl_config, only : a2o_fcorr_stat   
  use o_PARAM,           only : mstep
  IMPLICIT NONE
  
  real(kind=8), INTENT (INOUT) 	:: field2d(myDim_nod2D+eDim_nod2D)
  real(kind=8), INTENT (IN)	:: mask(myDim_nod2D+eDim_nod2D)
  INTEGER,      INTENT (IN)	:: n
  INTEGER,      INTENT (IN)	:: h !hemisphere: 0=GL, 1=NH, 2=SH
  logical,      INTENT (IN)	:: do_stats
  
  real(kind=8)			:: rmask(myDim_nod2D+eDim_nod2D)   
  real(kind=8)			:: weight(myDim_nod2D+eDim_nod2D)    
  real(kind=8)			:: flux_global(2), flux_local(2)
  real(kind=8)			:: eff_vol(2)
return !!!OIFS
  if (mstep==1) then
	if (mype == 0) write(*,*) 'Do not apply a correction at first time step for ', trim(cpl_recv(n))
	return
  end if

  
  !just keep NH or SH part in mask
  rmask=mask
  SELECT CASE(h)
	   CASE( 1 ) ; where (geo_coord_nod2D(2, :)< 0) !just NH 
	                     rmask=0. 
		       end where
	   CASE( 2 ) ; where (geo_coord_nod2D(2, :)>=0) !just SH 
	                     rmask=0. 
		       end where
	   CASE DEFAULT 
	   !keep global mask
  END SELECT
 
  !residual (net) fluxes; computes also oce_net_fluxes_*
  call compute_residual(field2d, rmask, n)
  
#ifdef VERBOSE
  if (mype == 0) then
  !atm net fluxes, oce net fluxes before modification
  write(*,'(3A,3e15.7)') 'atm NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    atm_net_fluxes_north(n), atm_net_fluxes_south(n), 	&
			    atm_net_fluxes_north(n)+ atm_net_fluxes_south(n)
			    
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			    oce_net_fluxes_north(n), oce_net_fluxes_south(n), 	&
			    oce_net_fluxes_north(n)+ oce_net_fluxes_south(n)
  end if			    
#endif
  if (do_stats)	then
     if (h==0 .or. h==1) then
        a2o_fcorr_stat(n, 1) =a2o_fcorr_stat(n, 1)+atm_net_fluxes_north(n)
        a2o_fcorr_stat(n, 4) =a2o_fcorr_stat(n, 4)+oce_net_fluxes_north(n)
     end if
     if (h==0 .or. h==2) then
        a2o_fcorr_stat(n, 2) =a2o_fcorr_stat(n, 2)+atm_net_fluxes_south(n)
        a2o_fcorr_stat(n, 5) =a2o_fcorr_stat(n, 5)+oce_net_fluxes_south(n)
     end if     
  end if

  !integrate (masked) abs(field2d) to get positive weights
  call integrate_2D(flux_global, flux_local, eff_vol, abs(field2d), rmask)
  
  !get weight pattern with integral 1
  if (abs(sum(flux_global))>1.e-10) then
    weight=abs(field2d/sum(flux_global))
  else
    !should rarely happen
    weight=1.0_8 / sum(eff_vol)
    write(*,*) 'Warning: Constant redistribution for flux ', trim(cpl_recv(n))
  end if
  
  !weight is still global 2D field, just keep NH or SH part
  where (rmask<1.e-10)
    weight=0.
  end where
  
  !redistribute the residual according to the mask
  SELECT CASE(h)
	   CASE( 0 ) ; field2d=field2d+weight*flux_correction_total(n) ! GL	       
	   CASE( 1 ) ; field2d=field2d+weight*flux_correction_north(n) ! NH
	   CASE( 2 ) ; field2d=field2d+weight*flux_correction_south(n) ! SH
  END SELECT
  
  !check conservation
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, rmask)
#ifdef VERBOSE
  if (mype == 0) then
  write(*,'(3A,3e15.7)') 'oce NH SH GL / ', trim(cpl_recv(n)), ': ', 		&
  			   flux_global(1), flux_global(2), sum(flux_global)
  end if
#endif
  
  !last flux			   
  if (n==nrecv .AND. mype==0) write(*,*) 'Fluxes have been modified.'  
END SUBROUTINE force_flux_consv

!
! Compute the difference between the net fluxes seen by the atmosphere
! and ocean component (residual flux) for flux n.
!
SUBROUTINE compute_residual(field2d, mask, n)

  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
				flux_correction_north, flux_correction_south,	&
				flux_correction_total
  use g_parsup, 	only : 	myDim_nod2D, eDim_nod2D
  
  IMPLICIT NONE
  
  real(kind=8), INTENT(IN)   :: field2d(myDim_nod2D+eDim_nod2D)
  real(kind=8), INTENT(IN)   :: mask(myDim_nod2D+eDim_nod2D)  
  INTEGER,      INTENT(IN)   :: n
  
  real(kind=8)               :: flux_global(2), flux_local(2)
  real(kind=8)               :: eff_vol(2)
  write(*,*) 'c1'
  !compute net flux (for flux n) on ocean side
  call integrate_2D(flux_global, flux_local, eff_vol, field2d, mask)
  oce_net_fluxes_north(n)=flux_global(1)
  oce_net_fluxes_south(n)=flux_global(2)
  
  !compute the residual fluxes for NH, SH and Global
  flux_correction_north(n)= atm_net_fluxes_north(n) - oce_net_fluxes_north(n)
  flux_correction_south(n)= atm_net_fluxes_south(n) - oce_net_fluxes_south(n)
  flux_correction_total(n)= flux_correction_north(n) + flux_correction_south(n)
  write(*,*) 'c2'
END SUBROUTINE compute_residual

!
! -field_2d (input) is any (partitioned) 2D field
! -flux_local  (returned) is the net local flux (for current pc)
! -flux_global (returned) is the communicated and summarized flux_local  
!
SUBROUTINE integrate_2D(flux_global, flux_local, eff_vol, field2d, mask)
 

  use g_parsup !myDim_nod2D, eDim_nod2D, MPI stuff
  use o_MESH,	only :	lump2d_north, lump2d_south
  
  IMPLICIT NONE

  real(kind=8), INTENT(OUT)  :: flux_global(2), flux_local(2)
  real(kind=8), INTENT(OUT)  :: eff_vol(2)
  real(kind=8), INTENT(IN)   :: field2d(myDim_nod2D+eDim_nod2D)
  real(kind=8), INTENT(IN)   :: mask(myDim_nod2D   +eDim_nod2D) 
   
  real(kind=8)               :: eff_vol_local(2)

  flux_local(1)=sum(lump2d_north*field2d(1:myDim_nod2D)*mask(1:myDim_nod2D))
  flux_local(2)=sum(lump2d_south*field2d(1:myDim_nod2D)*mask(1:myDim_nod2D))
  call MPI_AllREDUCE(flux_local, flux_global, 2, &
  		     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		     
		     
  eff_vol_local(1)=sum(lump2d_north*mask(1:myDim_nod2D))
  eff_vol_local(2)=sum(lump2d_south*mask(1:myDim_nod2D))
  call MPI_AllREDUCE(eff_vol_local, eff_vol,  2, & 
  		     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		     
END SUBROUTINE integrate_2D

!
!---------------------------------------------------------------------------------------------------
!  Receive atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
!SUBROUTINE net_rec_from_atm(action)
!  
!
!  use g_forcing_arrays,	only : 	atm_net_fluxes_north, atm_net_fluxes_south, 	&
!  				oce_net_fluxes_north, oce_net_fluxes_south, 	&
!				flux_correction_north, flux_correction_south,	&
!				flux_correction_total
!  use cpl_driver,	only : 	nrecv
!  
!  IMPLICIT NONE
!
!  LOGICAL, INTENT(IN)   :: action
!
!  if (action) then
!  atm_net_fluxes_north=0.
!  atm_net_fluxes_south=0.
!
!  atm_net_fluxes_north=(/ 1.,2.,3.,4. , 0.5e7, -0.1e10, -0.3e8, 0.5e-13, -0.7e17, -0.4e-6, 0.5e5, 0.0/)
!  atm_net_fluxes_south=(/ 1.,2.,3.,4. , 0.4e7, 0.3e06, -0.3e8, -0.3,-0.7e17, -0.4e-6, 0.2e6, 0.0/)
!  end if
!END SUBROUTINE net_rec_from_atm
!
!
!---------------------------------------------------------------------------------------------------
!  Receieve atmospheric net fluxes (atm_net_fluxes_north and atm_net_fluxes_south)
!
SUBROUTINE net_rec_from_atm(action)
!
  use g_forcing_arrays
  use g_parsup
  use cpl_driver

  IMPLICIT NONE

  LOGICAL,      INTENT (IN)   		          :: action
  INTEGER                                         :: my_global_rank, ierror
  INTEGER                                         :: n  
  INTEGER 					  :: status(MPI_STATUS_SIZE,npes) 
  INTEGER                                         :: request(2)
  real(kind=8)                 			  :: aux(nrecv)
RETURN !!!OIFS
  if (action) then
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_global_rank, ierror)
     atm_net_fluxes_north=0.
     atm_net_fluxes_south=0.
     if (my_global_rank==target_root) then
	CALL MPI_IRecv(atm_net_fluxes_north(1), nrecv, MPI_DOUBLE_PRECISION, source_root, 111, MPI_COMM_WORLD, request(1), MPIerr)
        CALL MPI_IRecv(atm_net_fluxes_south(1), nrecv, MPI_DOUBLE_PRECISION, source_root, 112, MPI_COMM_WORLD, request(2), MPIerr)
        CALL MPI_Waitall(2, request, status, MPIerr)
     end if
  call MPI_Barrier(MPI_COMM_FESOM, MPIerr)     
  call MPI_AllREDUCE(atm_net_fluxes_north(1), aux, nrecv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  atm_net_fluxes_north=aux
  call MPI_AllREDUCE(atm_net_fluxes_south(1), aux, nrecv, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
  atm_net_fluxes_south=aux
  end if
END SUBROUTINE net_rec_from_atm
#endif
