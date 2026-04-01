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
  integer                               :: fileunit, istat

  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=modelname, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'modelname', nmlfile, partit)
  
  read (fileunit, NML=timestep, iostat=istat)  
  if (istat /= 0) call check_namelist_read(fileunit, 'timestep', nmlfile, partit)

  read (fileunit, NML=clockinit, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'clockinit', nmlfile, partit)

  read (fileunit, NML=paths, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paths', nmlfile, partit)

  ! Set defaults for restart paths if not specified (backward compatibility)
  if (len_trim(RestartInPath) == 0) RestartInPath = ResultPath
  if (len_trim(RestartOutPath) == 0) RestartOutPath = ResultPath
  
  ! Ensure paths have trailing slash for consistent concatenation
  call ensure_trailing_slash(RestartInPath)
  call ensure_trailing_slash(RestartOutPath)
  call ensure_trailing_slash(ResultPath)
  
  ! Report the configuration to user
  if(partit%mype==0) then
    write(*,*) 'Restart path configuration:'
    write(*,*) '  Input path:  ', trim(RestartInPath)
    write(*,*) '  Output path: ', trim(RestartOutPath)
    if (trim(RestartInPath) == trim(RestartOutPath)) then
      write(*,*) '  Note: Using same directory for input and output (legacy mode)'
    else
      write(*,*) '  Note: Separate input/output directories configured for reproducibility'
    end if
  end if

  read (fileunit, NML=restart_log, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'restart_log', nmlfile, partit)

  read (fileunit, NML=ale_def, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'ale_def', nmlfile, partit)

  read (fileunit, NML=geometry, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'geometry', nmlfile, partit)

  read (fileunit, NML=calendar, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'calendar', nmlfile, partit)

  read (fileunit, NML=run_config, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'run_config', nmlfile, partit)

  read (fileunit, NML=icebergs, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'icebergs', nmlfile, partit)

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
  open (newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=oce_dyn, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'oce_dyn', nmlfile, partit)
  close (fileunit)

  nmlfile ='namelist.tra'    ! name of ocean namelist file
  open (newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=tracer_phys, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'tracer_phys', nmlfile, partit)
  close (fileunit)

  nmlfile ='namelist.forcing'    ! name of forcing namelist file
  open (newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=forcing_exchange_coeff, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'forcing_exchange_coeff', nmlfile, partit)

  read (fileunit, NML=forcing_bulk, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'forcing_bulk', nmlfile, partit)

  read (fileunit, NML=land_ice, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'land_ice', nmlfile, partit)

  read (fileunit, NML=age_tracer, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'age_tracer', nmlfile, partit)
  close (fileunit)

!   if(use_ice) then
!   nmlfile ='namelist.ice'    ! name of ice namelist file
!   open (newunit=fileunit, file=nmlfile)
! !   read (fileunit, NML=ice_dyn)
!   read (fileunit, NML=ice_therm)
!   close (fileunit)
!   endif
  
  nmlfile ='namelist.io'    ! name of forcing namelist file
  open (newunit=fileunit, file=nmlfile, status='OLD', iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=diag_list, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'diag_list', nmlfile, partit)
  close (fileunit)

#if defined (__recom)
  nmlfile ='namelist.recom'    ! name of recom namelist file
  open (newunit=fileunit, file=nmlfile, iostat=istat)
  if (istat /= 0) then
    if(partit%mype==0) then
      write(*,*) 'ERROR: Could not open namelist file ', trim(nmlfile)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
  endif
  read (fileunit, NML=pavariables, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pavariables', nmlfile, partit)

  read (fileunit, NML=pasinking, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pasinking', nmlfile, partit)

  read (fileunit, NML=painitialization_N, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'painitialization_N', nmlfile, partit)

  read (fileunit, NML=paArrhenius, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paArrhenius', nmlfile, partit)

  read (fileunit, NML=palimiter_function, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'palimiter_function', nmlfile, partit)

  read (fileunit, NML=palight_calculations, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'palight_calculations', nmlfile, partit)

  read (fileunit, NML=paphotosynthesis, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paphotosynthesis', nmlfile, partit)

  read (fileunit, NML=paassimilation, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paassimilation', nmlfile, partit)

  read (fileunit, NML=pairon_chem, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pairon_chem', nmlfile, partit)

  read (fileunit, NML=pazooplankton, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pazooplankton', nmlfile, partit)

  read (fileunit, NML=pasecondzooplankton, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pasecondzooplankton', nmlfile, partit)

  read (fileunit, NML=pathirdzooplankton, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pathirdzooplankton', nmlfile, partit)

  read (fileunit, NML=pagrazingdetritus, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pagrazingdetritus', nmlfile, partit)

  read (fileunit, NML=paaggregation, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paaggregation', nmlfile, partit)

  read (fileunit, NML=padin_rho_N, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'padin_rho_N', nmlfile, partit)

  read (fileunit, NML=padic_rho_C1, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'padic_rho_C1', nmlfile, partit)

  read (fileunit, NML=paphytoplankton_N, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paphytoplankton_N', nmlfile, partit)

  read (fileunit, NML=paphytoplankton_C, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paphytoplankton_C', nmlfile, partit)

  read (fileunit, NML=paphytoplankton_ChlA, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paphytoplankton_ChlA', nmlfile, partit)

  read (fileunit, NML=padetritus_N, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'padetritus_N', nmlfile, partit)

  read (fileunit, NML=padetritus_C, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'padetritus_C', nmlfile, partit)

  read (fileunit, NML=paheterotrophs, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paheterotrophs', nmlfile, partit)

  read (fileunit, NML=paseczooloss, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paseczooloss', nmlfile, partit)

  read (fileunit, NML=pathirdzooloss, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pathirdzooloss', nmlfile, partit)

  read (fileunit, NML=paco2lim, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paco2lim', nmlfile, partit)

  read (fileunit, NML=pairon, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pairon', nmlfile, partit)

  read (fileunit, NML=pacalc, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pacalc', nmlfile, partit)

  read (fileunit, NML=pabenthos_decay_rate, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'pabenthos_decay_rate', nmlfile, partit)

  read (fileunit, NML=paco2_flux_param, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paco2_flux_param', nmlfile, partit)

  read (fileunit, NML=paalkalinity_restoring, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paalkalinity_restoring', nmlfile, partit)

  read (fileunit, NML=paballasting, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paballasting', nmlfile, partit)

  read (fileunit, NML=paciso, iostat=istat)
  if (istat /= 0) call check_namelist_read(fileunit, 'paciso', nmlfile, partit)

  close (fileunit)
#endif

  if (use_transit) then
    if(partit%mype==0) print *, "Transient tracers are ON."
    nmlfile = 'namelist.transit'  ! name of transient tracers namelist file
    open (newunit=fileunit, file=nmlfile)
    read (fileunit, nml=transit_param)
    close (fileunit)
    if (anthro_transit) then
!     transient values of historical CO2, bomb radiocarbon, CFC-12, and SF6
      if(partit%mype==0) print *, "Reading transient input values from file: ", ifile_transit
      open (newunit=fileunit, file=ifile_transit)
      call read_transit_input(fileunit)
      close (fileunit)
    elseif (paleo_transit) then
!     transient values of atmospheric CO2 and radiocarbon reconstructions
!     under construction / not yet realized
    else
!     Spinup / equilibrium runs with constant tracer input specified in namelist.transit
      if(partit%mype==0) print *, "Reading constant input values from file: ", nmlfile
    end if
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

subroutine check_namelist_read(fileunit, nml_name, nmlfile, partit)
    use MOD_PARTIT
    use MOD_PARSUP
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    integer,          intent(in) :: fileunit
    character(len=*), intent(in) :: nml_name
    character(len=*), intent(in) :: nmlfile
    type(t_partit),   intent(in) :: partit
    character(len=256) :: line

    backspace(fileunit)
    read(fileunit, fmt='(A)') line
    if(partit%mype==0) then
        write(error_unit,'(A)') 'ERROR: Could not read namelist '//trim(nml_name)//' from '//trim(nmlfile)
        write(error_unit,'(A)') 'Invalid line in namelist: '//trim(line)
    endif
    call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
end subroutine check_namelist_read

! ==============================================================
! Helper subroutine to ensure paths have trailing slash
subroutine ensure_trailing_slash(path)
    implicit none
    character(len=*), intent(inout) :: path
    integer :: path_len
    
    path_len = len_trim(path)
    if (path_len > 0 .and. path(path_len:path_len) /= '/') then
        path = trim(path) // '/'
    end if
end subroutine ensure_trailing_slash
