! ==============================================================
subroutine setup_model
  implicit none
  call read_namelist    ! should be before clock_init
end subroutine setup_model
! ==============================================================
subroutine read_namelist
  ! Reads namelist files and overwrites default parameters.
  !
  ! Coded by Lars Nerger
  ! Modified by Qiang Wang, SD
  !--------------------------------------------------------------
  use o_param
  use i_param
  use i_therm_param
  use g_forcing_param
  use g_parsup
  use g_config
  use diagnostics, only: ldiag_solver,lcurt_stress_surf,lcurt_stress_surf, ldiag_energy, &
                         ldiag_dMOC, ldiag_DVD, diag_list
  use g_clock, only: timenew, daynew, yearnew
  use g_ic3d 
  implicit none

  character(len=MAX_PATH)   :: nmlfile
  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (20,file=nmlfile)
  read (20,NML=modelname)
  read (20,NML=timestep)
  read (20,NML=clockinit) 
  read (20,NML=paths)
  read (20,NML=restart_log)
  read (20,NML=ale_def)
  read (20,NML=geometry)
  read (20,NML=calendar)
  read (20,NML=run_config)
!!$  read (20,NML=machine)
  close (20)
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
  open (20,file=nmlfile)
  read (20,NML=oce_dyn)
  read (20,NML=oce_tra)
  read (20,NML=oce_init3d)
  close (20)

  nmlfile ='namelist.forcing'    ! name of forcing namelist file
  open (20,file=nmlfile)
  read (20,NML=forcing_exchange_coeff)
  read (20,NML=forcing_bulk)
  read (20,NML=land_ice)
  close (20)

  if(use_ice) then
  nmlfile ='namelist.ice'    ! name of ice namelist file
  open (20,file=nmlfile)
  read (20,NML=ice_dyn)
  read (20,NML=ice_therm)
  close (20)
  endif
  
  nmlfile ='namelist.io'    ! name of forcing namelist file
  open (20,file=nmlfile)
  read (20,NML=diag_list)
  close (20)

  if(mype==0) write(*,*) 'Namelist files are read in'
  
  !_____________________________________________________________________________
  ! Check for namelist parameter consistency
  if(mype==0) then
    
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
        call par_ex(0)
    endif
    

  endif

! if ((output_length_unit=='s').or.(int(real(step_per_day)/24.0)<=1)) use_means=.false.
end subroutine read_namelist
! =================================================================
subroutine get_run_steps(nsteps)
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------
  
  use g_clock
  use g_parsup
  implicit none

  integer      :: i, temp_year, temp_mon, temp_fleapyear, nsteps

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
     call par_ex(1)
     stop
  end if

  if(mype==0) write(*,*) nsteps, ' steps to run for ', runid, ' job submission'
end subroutine get_run_steps

    
! ==============================================================
