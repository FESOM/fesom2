! ==============================================================
subroutine setup_model
  implicit none
#ifdef __ifsinterface
  call read_namelist_ifs    ! should be before clock_init
#else
  call read_namelist    ! should be before clock_init
#endif
  call define_prog_tracer

end subroutine setup_model
! ==============================================================
#ifdef __ifsinterface
subroutine read_namelist_ifs
#else
subroutine read_namelist
#endif
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
  use g_clock, only: timenew, daynew, yearnew
  implicit none

  character(len=100)   :: nmlfile
#ifdef __ifsinterface
  integer, parameter :: iunit = 21
#else
  integer, parameter :: iunit = 20
#endif
  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (iunit,file=nmlfile)
  read (iunit,NML=modelname)
  read (iunit,NML=timestep)
  read (iunit,NML=clockinit) 
  read (iunit,NML=paths)
  read (iunit,NML=initialization)  
  read (iunit,NML=inout)
  read (iunit,NML=mesh_def)
  read (iunit,NML=geometry)
  read (iunit,NML=calendar)
  read (iunit,NML=run_config)
!!$  read (iunit,NML=machine)
  close (iunit)
  ! ==========
  ! compute dt
  ! ========== 
  dt=86400./float(step_per_day)
  if(mype==0) write(*,*) 'time step size is set to ', real(dt,4), 'sec'
  ! ==========
  ! degree2radian
  ! ==========
  cyclic_length=cyclic_length*rad
  alphaEuler=alphaEuler*rad 	
  betaEuler=betaEuler*rad
  gammaEuler=gammaEuler*rad

! =================================
 
  nmlfile ='namelist.oce'    ! name of ocean namelist file
  open (iunit,file=nmlfile)
  read (iunit,NML=oce_dyn)
  read (iunit,NML=oce_tra)
  close (iunit)

  nmlfile ='namelist.forcing'    ! name of forcing namelist file
  open (iunit,file=nmlfile)
  read (iunit,NML=forcing_exchange_coeff)
  read (iunit,NML=forcing_source)
  read (iunit,NML=forcing_bulk)
  read (iunit,NML=land_ice)
  close (iunit)

  if(use_ice) then
  nmlfile ='namelist.ice'    ! name of ice namelist file
  open (iunit,file=nmlfile)
  read (iunit,NML=ice_dyn)
  read (iunit,NML=ice_therm)
  close (iunit)
  endif

  if(mype==0) write(*,*) 'Namelist files are read in'

! if ((output_length_unit=='s').or.(int(real(step_per_day)/24.0)<=1)) use_means=.false.
#ifdef __ifsinterface
end subroutine read_namelist_ifs
#else
end subroutine read_namelist
#endif
! =================================================================
subroutine define_prog_tracer
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------
  
  use o_param
  use o_arrays, only : prog_tracer_name
  use g_parsup
  implicit none

  integer      :: j, num
  character(1) :: cageind
  character(4) :: tr_name

  ! allocate prog_tracer_name

  num_tracer=2  ! t and s

  allocate(prog_tracer_name(num_tracer))

  ! fill prog_tracer_name

  num=2  ! t and s
  prog_tracer_name(1)='temp'
  prog_tracer_name(2)='salt'

  if(mype==0) write(*,*) 'Number of prognostic ocean tracers: ',num
end subroutine define_prog_tracer
! ==============================================================
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
