! ==============================================================
subroutine setup_model(partit)
  USE MOD_PARTIT
  USE MOD_PARSUP
  use o_param
!   use i_therm_param
  use g_forcing_param
  use g_config
  use diagnostics, only: ldiag_solver,lcurt_stress_surf,lcurt_stress_surf, ldiag_Ri, ldiag_TurbFlux, ldiag_trflx, &
                         ldiag_dMOC, ldiag_DVD, diag_list
  use g_clock,     only: timenew, daynew, yearnew
  use g_ic3d
#ifdef __recom
  use recom_config
  use recom_ciso
#endif
  use mod_transit

  implicit none
  type(t_partit), intent(inout), target :: partit
  character(len=MAX_PATH)               :: nmlfile
  integer fileunit

  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=modelname)
  read (fileunit, NML=timestep)
  read (fileunit, NML=clockinit) 
  read (fileunit, NML=paths)
  read (fileunit, NML=restart_log)
  read (fileunit, NML=ale_def)
  read (fileunit, NML=geometry)
  read (fileunit, NML=calendar)
  read (fileunit, NML=run_config)
  read (fileunit,NML=icebergs)

!!$  read (fileunit, NML=machine)
  close (fileunit)
  
  
  ! ==========
  ! compute dt
  ! ========== 
  dt=86400._WP/real(step_per_day,WP)
  
  ! ==========
  ! degree2radian
  ! ==========
  cyclic_length=cyclic_length*rad
  alphaEuler=alphaEuler*rad 	
  betaEuler=betaEuler*rad
  gammaEuler=gammaEuler*rad

! =================================
 
  nmlfile ='namelist.oce'    ! name of ocean namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=oce_dyn)
  close (fileunit)

  nmlfile ='namelist.tra'    ! name of ocean namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=tracer_phys)
  close (fileunit)

  nmlfile ='namelist.forcing'    ! name of forcing namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=forcing_exchange_coeff)
  read (fileunit, NML=forcing_bulk)
  read (fileunit, NML=land_ice)
  read (fileunit, NML=age_tracer) !---age-code
  close (fileunit)

!   if(use_ice) then
!   nmlfile ='namelist.ice'    ! name of ice namelist file
!   open (newunit=fileunit, file=nmlfile)
! !   read (fileunit, NML=ice_dyn)
!   read (fileunit, NML=ice_therm)
!   close (fileunit)
!   endif
  
  nmlfile ='namelist.io'    ! name of forcing namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=diag_list)
  close (fileunit)

#if defined (__recom)
  nmlfile ='namelist.recom'    ! name of recom namelist file
  open (newunit=fileunit, file=nmlfile)
  read (fileunit, NML=pavariables)
  read (fileunit, NML=pasinking)
  read (fileunit, NML=painitialization_N)
  read (fileunit, NML=paArrhenius)
  read (fileunit, NML=palimiter_function)
  read (fileunit, NML=palight_calculations)
  read (fileunit, NML=paphotosynthesis)
  read (fileunit, NML=paassimilation)
  read (fileunit, NML=pairon_chem)
  read (fileunit, NML=pazooplankton)
  read (fileunit, NML=pasecondzooplankton)
  read (fileunit, NML=pathirdzooplankton)
  read (fileunit, NML=pagrazingdetritus)
  read (fileunit, NML=paaggregation)
  read (fileunit, NML=padin_rho_N)
  read (fileunit, NML=padic_rho_C1)
  read (fileunit, NML=paphytoplankton_N)
  read (fileunit, NML=paphytoplankton_C)
  read (fileunit, NML=paphytoplankton_ChlA)
  read (fileunit, NML=padetritus_N)
  read (fileunit, NML=padetritus_C)
  read (fileunit, NML=paheterotrophs)
  read (fileunit, NML=paseczooloss)
  read (fileunit, NML=pathirdzooloss)
  read (fileunit, NML=paco2lim)
  read (fileunit, NML=pairon)
  read (fileunit, NML=pacalc)
  read (fileunit, NML=pabenthos_decay_rate)
  read (fileunit, NML=paco2_flux_param)
  read (fileunit, NML=paalkalinity_restoring)
  read (fileunit, NML=paballasting)
  read (fileunit, NML=paciso)
  close (fileunit)
#endif

  if (use_transit) then
! Transient tracer input, input file names have to be specified in
! namelist.config, nml=run_config
    if(partit%mype==0) print *, "Transient tracers are ON. Tracer input file: ", ifile_transit
    open (20,file=ifile_transit)
    if (anthro_transit .or. paleo_transit) then
      call read_transit_input
    else
!     Spinup / equilibrium runs with constant tracer input,
!     read parameter values from namelist.oce
      read (20,nml=transit_param)
    end if
    close (20)
  end if

  if(partit%mype==0) write(*,*) 'Namelist files are read in'
  
  !_____________________________________________________________________________
  ! Check for namelist parameter consistency
  if(partit%mype==0) then
    
    ! check for valid step per day number
    if (mod(86400,step_per_day)==0) then
        write(*,*) 'time step size is set to ', dt, 'sec'
    else
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: The used step_per_day variable is not valid, model'
        write(*,*) '        simulation STOPS here. The variable step_per_day must be'
        write(*,*) '        an integer multiple of 86400 (mod(86400,step_per_day)==0).'
        write(*,*) '        That means valid steps per day are:'
        write(*,*) '        ... 32(45min), 36(40min), 40, 45, 48(30min), 50, 54, 60(24min), '
        write(*,*) '        64, 72(20min), 75, 80, 90, 96(15min), 100, 108, 120, 128,  '
        write(*,*) '        135, 144(10min), 150, 160(9min), 180(8min), 192, 200, 216, '
        write(*,*) '        225, 240(6min), 270, 288(5min), 300, 320, 360(4min), 384,  '
        write(*,*) '        400, 432, 450, 480(3min), 540, 576, 600, 640, 675, 720(2min),'
        write(*,*) '        800, 864, 900, 960, 1080, 1152, 1200, 1350, 1440(1min), ...'
        write(*,*)
        write(*,*) '        --> check your namelist.config !!!'
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
    endif
    

  endif
! if ((output_length_unit=='s').or.(int(real(step_per_day)/24.0)<=1)) use_means=.false.
end subroutine setup_model


#if defined(__usetp)
! kh 11.11.21 read num_fesom_groups for multi FESOM group loop parallelization
! =================================================================
subroutine read_namelist_run_config(partit)
  ! Reads run_config namelist and overwrite default parameters.
  !
  ! kh 11.11.21 Copied by Kai Himstedt (based on read_namelist)
  !--------------------------------------------------------------
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_config
  implicit none
  type(t_partit), intent(inout), target :: partit

  character(len=100)   :: nmlfile
  integer fileunit

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (newunit=fileunit, file=nmlfile)

  open (fileunit,file=nmlfile)
!  read (fileunit,NML=run_config)
  read (fileunit,NML=run_config_tp)
  close (fileunit)
end subroutine read_namelist_run_config
#endif
! =================================================================
subroutine get_run_steps(nsteps, partit)
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------  
  use g_clock
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none

  type(t_partit), intent(inout) :: partit
  integer,        intent(inout) :: nsteps 
  integer                       :: i, temp_year, temp_mon, temp_fleapyear

  ! clock should have been inialized before calling this routine

  if(run_length_unit=='s') then
     nsteps=run_length
  elseif(run_length_unit=='d') then
     nsteps=step_per_day*run_length
  elseif(run_length_unit=='m') then
     nsteps=0
     temp_mon=month-1
     temp_year=yearnew
     temp_fleapyear=fleapyear
     do i=1,run_length
        temp_mon=temp_mon+1
        if(temp_mon>12) then
           temp_year=temp_year+1
           temp_mon=1
           call check_fleapyr(temp_year, temp_fleapyear)
        end if
        nsteps=nsteps+step_per_day*num_day_in_month(temp_fleapyear,temp_mon)
     end do
  elseif(run_length_unit=='y') then
     nsteps=0
     do i=1,run_length
        temp_year=yearnew+i-1
        call check_fleapyr(temp_year, temp_fleapyear)
        nsteps=nsteps+step_per_day*(365+temp_fleapyear)
     end do
  else
     write(*,*) 'Run length unit ', run_length_unit, ' is not defined.'
     write(*,*) 'Please check and update the code.'
     call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
     stop
  end if

  if(partit%mype==0) write(*,*) nsteps, ' steps to run for ', runid, ' job submission'
end subroutine get_run_steps

    
! ==============================================================
