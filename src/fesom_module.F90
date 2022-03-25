! synopsis: save any derived types we initialize
!           so they can be reused after fesom_init
module fesom_main_storage_module
  USE MOD_MESH
  USE MOD_ICE
  USE MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE o_ARRAYS
  USE o_PARAM
  use g_clock
  use g_config
  use g_comm_auto
  use g_forcing_arrays
  use io_RESTART
  use io_MEANDATA
  use io_mesh_info
  use diagnostics
  use mo_tidal
  use tracer_init_interface
  use ocean_setup_interface
  use ice_setup_interface
  use ocean2ice_interface
  use oce_fluxes_interface
  use update_atm_forcing_interface
  use before_oce_step_interface
  use oce_timestep_ale_interface
  use read_mesh_interface
  use fesom_version_info_module
  use command_line_options_module
  ! Define icepack module
#if defined (__icepack)
  use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif

#if defined (__oasis)
  use cpl_driver
#endif

  implicit none
    
  type :: fesom_main_storage_type

    integer           :: n, from_nstep, offset, row, i, provided
    integer           :: which_readr ! read which restart files (0=netcdf, 1=core dump,2=dtype)
    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESOM
    real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t0_ice, t1_ice, t0_frc, t1_frc
    real(kind=WP)     :: rtime_fullice,    rtime_write_restart, rtime_write_means, rtime_compute_diag, rtime_read_forcing
    real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing 
    real(kind=real32) :: rtime_setup_ice,  rtime_setup_other, rtime_setup_restart
    real(kind=real32) :: runtime_alltimesteps


    type(t_mesh)   mesh
    type(t_tracer) tracers
    type(t_dyn)    dynamics
    type(t_partit) partit
    type(t_ice)    ice


    character(LEN=256)               :: dump_dir, dump_filename
    logical                          :: L_EXISTS
    type(t_mesh)   mesh_copy
    type(t_tracer) tracers_copy
    type(t_dyn)    dynamics_copy
    type(t_ice)    ice_copy

    character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version_txt
    integer mpi_version_len
    logical fesom_did_mpi_init
    
  end type
  type(fesom_main_storage_type), save, target :: f

end module


! synopsis: main FESOM program split into 3 parts
!           this way FESOM can e.g. be used as a library with an external time loop driver
!           used with IFS-FESOM
module fesom_module
  implicit none
  public fesom_init, fesom_runloop, fesom_finalize
  private

contains
 
  subroutine fesom_init(fesom_total_nsteps)
      use fesom_main_storage_module
      integer, intent(out) :: fesom_total_nsteps
      ! EO parameters
      logical mpi_is_initialized

#if !defined  __ifsinterface
      if(command_argument_count() > 0) then
        call command_line_options%parse()
        stop
      end if
#endif
      
      mpi_is_initialized = .false.
      f%fesom_did_mpi_init = .false.

#ifndef __oifs
        !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
        !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
        call MPI_Initialized(mpi_is_initialized, f%i)
        if(.not. mpi_is_initialized) then
            ! do not initialize MPI here if it has been initialized already, e.g. via IFS when fesom is called as library (__ifsinterface is defined)
            call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, f%provided, f%i)
            f%fesom_did_mpi_init = .true.
        end if
#endif
    

#if defined (__oasis)
        call cpl_oasis3mct_init(f%partit,f%partit%MPI_COMM_FESOM)
#endif
        f%t1 = MPI_Wtime()

        call par_init(f%partit)

        f%mype          =>f%partit%mype
        f%MPIerr        =>f%partit%MPIerr
        f%MPI_COMM_FESOM=>f%partit%MPI_COMM_FESOM
        f%npes          =>f%partit%npes
        if(f%mype==0) then
            write(*,*)
            print *,"FESOM2 git SHA: "//fesom_git_sha()
            call MPI_Get_library_version(f%mpi_version_txt, f%mpi_version_len, f%MPIERR)
            print *,"MPI library version: "//trim(f%mpi_version_txt)
            print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
            print *, achar(27)//'[7;32m'//' --> FESOM BUILDS UP MODEL CONFIGURATION                    '//achar(27)//'[0m'
        end if
        !=====================
        ! Read configuration data,  
        ! load the mesh and fill in 
        ! auxiliary mesh arrays
        !=====================
        call setup_model(f%partit)  ! Read Namelists, always before clock_init
        call clock_init(f%partit)   ! read the clock file 
        call get_run_steps(fesom_total_nsteps, f%partit)
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call mesh_setup'//achar(27)//'[0m'
        call mesh_setup(f%partit, f%mesh)

        if (f%mype==0) write(*,*) 'FESOM mesh_setup... complete'
    
        !=====================
        ! Allocate field variables 
        ! and additional arrays needed for 
        ! fancy advection etc.  
        !=====================
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call check_mesh_consistency'//achar(27)//'[0m'
        call check_mesh_consistency(f%partit, f%mesh)
        if (f%mype==0) f%t2=MPI_Wtime()

        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
        call dynamics_init(f%dynamics, f%partit, f%mesh)
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
        call tracer_init(f%tracers, f%partit, f%mesh)                ! allocate array of ocean tracers (derived type "t_tracer")
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
        call arrays_init(f%tracers%num_tracers, f%partit, f%mesh)    ! allocate other arrays (to be refactured same as tracers in the future)
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean_setup'//achar(27)//'[0m'
        call ocean_setup(f%dynamics, f%tracers, f%partit, f%mesh)

        if (f%mype==0) then
           write(*,*) 'FESOM ocean_setup... complete'
           f%t3=MPI_Wtime()
        endif
        call forcing_setup(f%partit, f%mesh)

        if (f%mype==0) f%t4=MPI_Wtime()
        if (use_ice) then 
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_setup'//achar(27)//'[0m'
            call ice_setup(f%ice, f%tracers, f%partit, f%mesh)
            f%ice%ice_steps_since_upd = f%ice%ice_ave_steps-1
            f%ice%ice_update=.true.
            if (f%mype==0) write(*,*) 'EVP scheme option=', f%ice%whichEVP
        else 
            ! create a dummy ice derived type with only a_ice, m_ice, m_snow and 
            ! uvice since oce_timesteps still needs in moment
            ! ice as an input for mo_convect(ice, partit, mesh), call 
            ! compute_vel_rhs(ice, dynamics, partit, mesh),  
            ! call write_step_info(...) and call check_blowup(...)
            call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)
        endif
        
        if (f%mype==0) f%t5=MPI_Wtime()
        call compute_diagnostics(0, f%dynamics, f%tracers, f%partit, f%mesh) ! allocate arrays for diagnostic
#if defined (__oasis)
        call cpl_oasis3mct_define_unstr(f%partit, f%mesh)
        if(f%mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
#endif

#if defined (__icepack)
        !=====================
        ! Setup icepack
        !=====================
        if (f%mype==0) write(*,*) 'Icepack: reading namelists from namelist.icepack'
        call set_icepack(f%ice, f%partit)
        call alloc_icepack
        call init_icepack(f%ice, f%tracers%data(1), f%mesh)
        if (f%mype==0) write(*,*) 'Icepack: setup complete'
#endif
        call clock_newyear                        ! check if it is a new year
        if (f%mype==0) f%t6=MPI_Wtime()
        !___CREATE NEW RESTART FILE IF APPLICABLE___________________________________
        call restart(0, r_restart, f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        if (f%mype==0) f%t7=MPI_Wtime()
        ! store grid information into netcdf file
        if (.not. r_restart) call write_mesh_info(f%partit, f%mesh)

        !___IF RESTART WITH ZLEVEL OR ZSTAR IS DONE, ALSO THE ACTUAL LEVELS AND ____
        !___MIDDEPTH LEVELS NEEDS TO BE CALCULATET AT RESTART_______________________
        if (r_restart .and. .not. f%which_readr==2) then
            call restart_thickness_ale(f%partit, f%mesh)
        end if
        if (f%mype==0) then
           f%t8=MPI_Wtime()
    
           f%rtime_setup_mesh    = real( f%t2 - f%t1              ,real32)
           f%rtime_setup_ocean   = real( f%t3 - f%t2              ,real32)
           f%rtime_setup_forcing = real( f%t4 - f%t3              ,real32)
           f%rtime_setup_ice     = real( f%t5 - f%t4              ,real32)
           f%rtime_setup_restart = real( f%t7 - f%t6              ,real32)
           f%rtime_setup_other   = real((f%t8 - f%t7) + (f%t6 - f%t5) ,real32)

           write(*,*) '=========================================='
           write(*,*) 'MODEL SETUP took on mype=0 [seconds]      '
           write(*,*) 'runtime setup total      ',real(f%t8-f%t1,real32)      
           write(*,*) ' > runtime setup mesh    ',f%rtime_setup_mesh   
           write(*,*) ' > runtime setup ocean   ',f%rtime_setup_ocean  
           write(*,*) ' > runtime setup forcing ',f%rtime_setup_forcing
           write(*,*) ' > runtime setup ice     ',f%rtime_setup_ice    
           write(*,*) ' > runtime setup restart ',f%rtime_setup_restart
           write(*,*) ' > runtime setup other   ',f%rtime_setup_other 
            write(*,*) '============================================' 
        endif

    !    f%dump_dir='DUMP/'
    !    INQUIRE(file=trim(f%dump_dir), EXIST=f%L_EXISTS)
    !    if (.not. f%L_EXISTS) call system('mkdir '//trim(f%dump_dir))

    !    write (f%dump_filename, "(A7,I7.7)") "t_mesh.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%mesh
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%mesh_copy
    !    close (f%mype+300)
         
    !    write (f%dump_filename, "(A9,I7.7)") "t_tracer.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%tracers
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%dynamics_copy
    !    close (f%mype+300)

    !    write (f%dump_filename, "(A9,I7.7)") "t_dynamics.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%dynamics
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%tracers_copy
    !    close (f%mype+300)
    
    !call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype)
    !stop
    !         
    !    if (f%mype==10) write(,) f%mesh1%ssh_stiff%values-f%mesh%ssh_stiff%value    
  
    ! Initialize timers
    f%rtime_fullice       = 0._WP
    f%rtime_write_restart = 0._WP
    f%rtime_write_means   = 0._WP
    f%rtime_compute_diag  = 0._WP
    f%rtime_read_forcing  = 0._WP

    f%from_nstep = 1
  end subroutine


  subroutine fesom_runloop(current_nsteps)
    use fesom_main_storage_module
    integer, intent(in) :: current_nsteps 
    ! EO parameters
    integer n

    !=====================
    ! Time stepping
    !=====================

    if (f%mype==0) write(*,*) 'FESOM start iteration before the barrier...'
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    
    if (f%mype==0) then
       write(*,*) 'FESOM start iteration after the barrier...'
       f%t0 = MPI_Wtime()
    endif
    if(f%mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESOM STARTS TIME LOOP                                 '//achar(27)//'[0m'
    end if
    !___MODEL TIME STEPPING LOOP________________________________________________
    if (use_global_tides) then
       call foreph_ini(yearnew, month, f%partit)
    end if
    do n=f%from_nstep, f%from_nstep-1+current_nsteps        
        if (use_global_tides) then
           call foreph(f%partit, f%mesh)
        end if
        mstep = n
        if (mod(n,logfile_outfreq)==0 .and. f%mype==0) then
            write(*,*) 'FESOM ======================================================='
!             write(*,*) 'FESOM step:',n,' day:', n*dt/24./3600.,
            write(*,*) 'FESOM step:',n,' day:', daynew,' year:',yearnew 
            write(*,*)
        end if
#if defined (__oifs) || defined (__oasis)
            seconds_til_now=INT(dt)*(n-1)
#endif
        call clock      
        !___compute horizontal velocity on nodes (originaly on elements)________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_vel_nodes'//achar(27)//'[0m'
        call compute_vel_nodes(f%dynamics, f%partit, f%mesh)
        
        !___model sea-ice step__________________________________________________
        f%t1 = MPI_Wtime()
        if(use_ice) then
            !___compute fluxes from ocean to ice________________________________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean2ice(n)'//achar(27)//'[0m'
            call ocean2ice(f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
            
            !___compute update of atmospheric forcing____________________________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
            f%t0_frc = MPI_Wtime()
            call update_atm_forcing(n, f%ice, f%tracers, f%partit, f%mesh)
            f%t1_frc = MPI_Wtime()       
            !___compute ice step________________________________________________
            if (f%ice%ice_steps_since_upd>=f%ice%ice_ave_steps-1) then
                f%ice%ice_update=.true.
                f%ice%ice_steps_since_upd = 0
            else
                f%ice%ice_update=.false.
                f%ice%ice_steps_since_upd=f%ice%ice_steps_since_upd+1
            endif
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
            if (f%ice%ice_update) call ice_timestep(n, f%ice, f%partit, f%mesh)  
            !___compute fluxes to the ocean: heat, freshwater, momentum_________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom...'//achar(27)//'[0m'
            call oce_fluxes_mom(f%ice, f%dynamics, f%partit, f%mesh) ! momentum only
            call oce_fluxes(f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        end if
        call before_oce_step(f%dynamics, f%tracers, f%partit, f%mesh) ! prepare the things if required
        f%t2 = MPI_Wtime()
        
        !___model ocean step____________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_timestep_ale'//achar(27)//'[0m'
        call oce_timestep_ale(n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)

        f%t3 = MPI_Wtime()
        !___compute energy diagnostics..._______________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1)'//achar(27)//'[0m'
        call compute_diagnostics(1, f%dynamics, f%tracers, f%partit, f%mesh)

        f%t4 = MPI_Wtime()
        !___prepare output______________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
        call output (n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)

        f%t5 = MPI_Wtime()
        call restart(n, .false., f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        f%t6 = MPI_Wtime()
        
        f%rtime_fullice       = f%rtime_fullice       + f%t2 - f%t1
        f%rtime_compute_diag  = f%rtime_compute_diag  + f%t4 - f%t3
        f%rtime_write_means   = f%rtime_write_means   + f%t5 - f%t4   
        f%rtime_write_restart = f%rtime_write_restart + f%t6 - f%t5
        f%rtime_read_forcing  = f%rtime_read_forcing  + f%t1_frc - f%t0_frc
    end do

    f%from_nstep = f%from_nstep+current_nsteps
  end subroutine


  subroutine fesom_finalize()
    use fesom_main_storage_module
    ! EO parameters
    real(kind=real32) :: mean_rtime(15), max_rtime(15), min_rtime(15)

    call finalize_output()
    call finalize_restart()

    !___FINISH MODEL RUN________________________________________________________

    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    if (f%mype==0) then
       f%t1 = MPI_Wtime()
       f%runtime_alltimesteps = real(f%t1-f%t0,real32)
       write(*,*) 'FESOM Run is finished, updating clock'
    endif

    mean_rtime(1)  = rtime_oce         
    mean_rtime(2)  = rtime_oce_mixpres 
    mean_rtime(3)  = rtime_oce_dyn     
    mean_rtime(4)  = rtime_oce_dynssh  
    mean_rtime(5)  = rtime_oce_solvessh
    mean_rtime(6)  = rtime_oce_GMRedi  
    mean_rtime(7)  = rtime_oce_solvetra
    mean_rtime(8)  = rtime_ice         
    mean_rtime(9)  = rtime_tot  
    mean_rtime(10) = f%rtime_fullice - f%rtime_read_forcing 
    mean_rtime(11) = f%rtime_compute_diag
    mean_rtime(12) = f%rtime_write_means
    mean_rtime(13) = f%rtime_write_restart
    mean_rtime(14) = f%rtime_read_forcing   
    
    max_rtime(1:14) = mean_rtime(1:14)
    min_rtime(1:14) = mean_rtime(1:14)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
    mean_rtime(1:14) = mean_rtime(1:14) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)

    if(f%fesom_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype) ! finalize MPI before FESOM prints its stats block, otherwise there is sometimes output from other processes from an earlier time in the programm AFTER the starts block (with parastationMPI)
    if (f%mype==0) then
        41 format (a35,a10,2a15) !Format for table heading
        42 format (a30,3f15.4)   !Format for table content

        print 41, '___MODEL RUNTIME per task [seconds]','_____mean_','___________min_', '___________max_'
        print 42, '  runtime ocean:              ',    mean_rtime(1),     min_rtime(1),      max_rtime(1)
        print 42, '    > runtime oce. mix,pres. :',    mean_rtime(2),     min_rtime(2),      max_rtime(2)
        print 42, '    > runtime oce. dyn. u,v,w:',    mean_rtime(3),     min_rtime(3),      max_rtime(3)
        print 42, '    > runtime oce. dyn. ssh  :',    mean_rtime(4),     min_rtime(4),      max_rtime(4)
        print 42, '    > runtime oce. solve ssh :',    mean_rtime(5),     min_rtime(5),      max_rtime(5)
        print 42, '    > runtime oce. GM/Redi   :',    mean_rtime(6),     min_rtime(6),      max_rtime(6)
        print 42, '    > runtime oce. tracer    :',    mean_rtime(7),     min_rtime(7),      max_rtime(7)
        print 42, '  runtime ice  :              ',    mean_rtime(10),    min_rtime(10),     max_rtime(10)
        print 42, '    > runtime ice step :      ',    mean_rtime(8),     min_rtime(8),      max_rtime(8)
        print 42, '  runtime diag:               ',    mean_rtime(11),    min_rtime(11),     max_rtime(11)
        print 42, '  runtime output:             ',    mean_rtime(12),    min_rtime(12),     max_rtime(12)
        print 42, '  runtime restart:            ',    mean_rtime(13),    min_rtime(13),     max_rtime(13)
        print 42, '  runtime forcing:            ',    mean_rtime(14),    min_rtime(14),     max_rtime(14)
        print 42, '  runtime total (ice+oce):    ',    mean_rtime(9),     min_rtime(9),      max_rtime(9)

        43 format (a33,i15)        !Format Ncores
        44 format (a33,i15)        !Format OMP threads
        45 format (a33,f15.4,a4)   !Format runtime

        write(*,*)
        write(*,*) '======================================================'
        write(*,*) '================ BENCHMARK RUNTIME ==================='
        print 43, '    Number of cores :            ',f%npes
#if defined(_OPENMP)
        print 44, '    Max OpenMP threads :         ',OMP_GET_MAX_THREADS()
#endif
        print 45, '    Runtime for all timesteps :  ',f%runtime_alltimesteps,' sec'
        write(*,*) '======================================================'
        write(*,*)
    end if    
!   call clock_finish  
  end subroutine

end module
