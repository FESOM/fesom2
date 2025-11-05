! synopsis: Simplified toy ocean model module
!           Minimal version of fesom_module.F90 for ocean-only toy configurations
!           Removes: ice model, icebergs, OASIS coupling, RECOM, transient tracers,
!                    tidal forcing, land ice coupling, age tracers, OpenIFS coupling
module fesom_toy_storage_module
  ! Core data structures
  USE MOD_MESH
  USE MOD_ICE       ! Keep for dummy ice structure
  USE MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE o_ARRAYS
  USE o_PARAM

  ! Configuration and control
  use g_clock
  use g_config
  use g_comm_auto
  use g_forcing_arrays

  ! I/O
  use io_RESTART
  use io_MEANDATA
  use io_mesh_info

  ! Diagnostics
  use diagnostics

  ! Initialization and setup interfaces
  use tracer_init_interface
  use ocean_setup_interface
  use read_mesh_interface

  ! Ocean timestep and dynamics
  use before_oce_step_interface
  use oce_timestep_ale_interface

  ! Version and utilities
  use fesom_version_info_module
  use command_line_options_module
  use, intrinsic :: iso_fortran_env, only : real32

  ! Toy ocean specific
  use Toy_Channel_Soufflet, only: compute_zonal_mean

#if defined (FESOM_PROFILING)
  use fesom_profiler
#endif

  implicit none

  type :: fesom_toy_storage_type
    integer           :: n, from_nstep, offset, row, i, provided, id
    integer           :: which_readr ! read which restart files (0=netcdf, 1=core dump, 2=dtype)
    integer           :: total_nsteps
    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESOM, MPI_COMM_WORLD
    real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8
    real(kind=WP)     :: rtime_write_restart, rtime_write_means, rtime_compute_diag
    real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing
    real(kind=real32) :: rtime_setup_other, rtime_setup_restart
    real(kind=real32) :: runtime_alltimesteps

    type(t_mesh)   mesh
    type(t_tracer) tracers
    type(t_dyn)    dynamics
    type(t_partit) partit
    type(t_ice)    ice        ! Dummy ice structure for toy ocean

    character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version_txt
    integer mpi_version_len
    logical fesom_did_mpi_init

  end type fesom_toy_storage_type
  type(fesom_toy_storage_type), save, target :: f

end module fesom_toy_storage_module


! synopsis: Toy ocean FESOM module with init/runloop/finalize
!           Simplified for ocean-only configurations (e.g., Soufflet channel test case)
module fesom_toy_module
#if defined (FESOM_PROFILING)
  use fesom_profiler
#endif
  implicit none
  public toy_fesom_init, toy_fesom_runloop, toy_fesom_finalize
  private

contains

  subroutine toy_fesom_init(fesom_total_nsteps)
      use fesom_toy_storage_module
      integer, intent(out) :: fesom_total_nsteps
      ! EO parameters
      logical mpi_is_initialized
      integer              :: tr_num

#if !defined  __ifsinterface
      if(command_argument_count() > 0) then
        call command_line_options%parse()
        stop
      end if
#endif

      mpi_is_initialized = .false.
      f%fesom_did_mpi_init = .false.

      ! MPI initialization (unless already initialized by coupling framework)
      call MPI_Initialized(mpi_is_initialized, f%i)
      if(.not. mpi_is_initialized) then
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, f%provided, f%i)
          f%fesom_did_mpi_init = .true.
      end if

      f%t1 = MPI_Wtime()

      ! Initialize enhanced profiler
#if defined (FESOM_PROFILING)
      call fesom_profiler_init(.true.)
      call fesom_profiler_start("toy_fesom_init_total")
#endif

      ! Parallel environment setup
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("par_init")
#endif
      call par_init(f%partit)
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("par_init")
#endif

      f%mype          =>f%partit%mype
      f%MPIerr        =>f%partit%MPIerr
      f%MPI_COMM_FESOM=>f%partit%MPI_COMM_FESOM
      f%MPI_COMM_WORLD=>f%partit%MPI_COMM_WORLD
      f%npes          =>f%partit%npes

      if(f%mype==0) then
          call plot_fesomlogo()
          write(*,*)
          print *,"FESOM2 TOY OCEAN MODEL"
          print *,"FESOM2 git SHA: "//fesom_git_sha()
          call MPI_Get_library_version(f%mpi_version_txt, f%mpi_version_len, f%MPIERR)
          print *,"MPI library version: "//trim(f%mpi_version_txt)
          print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
          print *, achar(27)//'[7;32m'//' --> TOY OCEAN: BUILDING MODEL CONFIGURATION                '//achar(27)//'[0m'
      end if

      !=====================
      ! Read configuration data, load mesh
      !=====================
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call setup_model'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("setup_model")
#endif
      call setup_model(f%partit)  ! Read Namelists
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("setup_model")
#endif

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call clock_init'//achar(27)//'[0m'
      call clock_init(f%partit)   ! read the clock file

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call get_run_steps'//achar(27)//'[0m'
      call get_run_steps(fesom_total_nsteps, f%partit)
      f%total_nsteps=fesom_total_nsteps
#if defined (FESOM_PROFILING)
      call fesom_profiler_set_timesteps(fesom_total_nsteps)
      call fesom_profiler_set_timestep_size(86400.0d0 / real(step_per_day, kind=8))
#endif

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call mesh_setup'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("mesh_setup")
#endif
      call mesh_setup(f%partit, f%mesh)
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("mesh_setup")
#endif

      if (f%mype==0) write(*,*) 'TOY OCEAN mesh_setup... complete'

      if (f%mype==0) f%t2=MPI_Wtime()

      !=====================
      ! Allocate field variables
      !=====================
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call check_mesh_consistency'//achar(27)//'[0m'
      call check_mesh_consistency(f%partit, f%mesh)

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
      call dynamics_init(f%dynamics, f%partit, f%mesh)

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
      call tracer_init(f%tracers, f%partit, f%mesh)

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
      call arrays_init(f%tracers%num_tracers, f%partit, f%mesh)

      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean_setup'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
      call fesom_profiler_start("ocean_setup")
#endif
      call ocean_setup(f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
      call fesom_profiler_end("ocean_setup")
#endif

      if (f%mype==0) then
         write(*,*) 'TOY OCEAN ocean_setup... complete'
         f%t3=MPI_Wtime()
      endif

      ! For toy ocean, we may still need forcing arrays (even if zero for pure toy cases)
      call forcing_setup(f%partit, f%mesh)

      if (f%mype==0) f%t4=MPI_Wtime()

      ! Create dummy ice structure (no actual ice model)
      if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_init_toyocean_dummy'//achar(27)//'[0m'
      call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)

      if (f%mype==0) f%t5=MPI_Wtime()

      call compute_diagnostics(0, f%dynamics, f%tracers, f%ice, f%partit, f%mesh)

      call clock_newyear  ! check if it is a new year
      if (f%mype==0) f%t6=MPI_Wtime()

      !___READ INITIAL CONDITIONS IF THIS IS A RESTART RUN________________________
      if (r_restart) then
          call read_initial_conditions(f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
      end if
      if (f%mype==0) f%t7=MPI_Wtime()

      ! Recompute zonal profiles for toy ocean restart (Soufflet case)
      if (toy_ocean .and. r_restart) then
          SELECT CASE (TRIM(which_toy))
              CASE ("soufflet")
                  call compute_zonal_mean(f%dynamics, f%tracers, f%partit, f%mesh)
          END SELECT
      end if

      ! Store grid information into netcdf file
      if (.not. r_restart) call write_mesh_info(f%partit, f%mesh)

      !___IF RESTART WITH ZLEVEL OR ZSTAR, CALCULATE LEVELS______________________
      if (r_restart .and. .not. f%which_readr==2) then
          call restart_thickness_ale(f%partit, f%mesh)
      end if

      if (f%mype==0) then
         f%t8=MPI_Wtime()

         f%rtime_setup_mesh    = real( f%t2 - f%t1              ,real32)
         f%rtime_setup_ocean   = real( f%t3 - f%t2              ,real32)
         f%rtime_setup_forcing = real( f%t4 - f%t3              ,real32)
         f%rtime_setup_restart = real( f%t7 - f%t6              ,real32)
         f%rtime_setup_other   = real((f%t8 - f%t7) + (f%t6 - f%t5) + (f%t5 - f%t4) ,real32)

         write(*,*) '=========================================='
         write(*,*) 'TOY OCEAN MODEL SETUP [seconds]           '
         write(*,*) 'runtime setup total      ',real(f%t8-f%t1,real32)
         write(*,*) ' > runtime setup mesh    ',f%rtime_setup_mesh
         write(*,*) ' > runtime setup ocean   ',f%rtime_setup_ocean
         write(*,*) ' > runtime setup forcing ',f%rtime_setup_forcing
         write(*,*) ' > runtime setup restart ',f%rtime_setup_restart
         write(*,*) ' > runtime setup other   ',f%rtime_setup_other
         write(*,*) '============================================'
      endif

      ! Initialize timers
      f%rtime_write_restart = 0._WP
      f%rtime_write_means   = 0._WP
      f%rtime_compute_diag  = 0._WP

      f%from_nstep = 1

#if defined (FESOM_PROFILING)
      call fesom_profiler_end("toy_fesom_init_total")
#endif

  end subroutine toy_fesom_init


  subroutine toy_fesom_runloop(current_nsteps)
    use fesom_toy_storage_module
    integer, intent(in) :: current_nsteps
    ! EO parameters
    integer n, nstart, ntotal

    !=====================
    ! Time stepping
    !=====================
    if (f%mype==0) write(*,*) 'TOY OCEAN start iteration before the barrier...'
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
    if (f%mype==0) then
       write(*,*) 'TOY OCEAN start iteration after the barrier...'
       f%t0 = MPI_Wtime()
    endif
    if(f%mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> TOY OCEAN: TIME LOOP STARTS                            '//achar(27)//'[0m'
    end if

#if defined (FESOM_PROFILING)
    call fesom_profiler_start("toy_fesom_runloop_total")
#endif

    !___MODEL TIME STEPPING LOOP________________________________________________
    nstart=f%from_nstep
    ntotal=f%from_nstep-1+current_nsteps

    do n=nstart, ntotal
        mstep = n
        if (mod(n,logfile_outfreq)==0 .and. f%mype==0) then
            write(*,*) 'TOY OCEAN ==============================================='
            write(*,*) 'TOY OCEAN step:',n,' day:', daynew,' year:',yearnew
            write(*,*)
        end if

        call clock

        !___compute horizontal velocity on nodes (originally on elements)________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_vel_nodes'//achar(27)//'[0m'
        call compute_vel_nodes(f%dynamics, f%partit, f%mesh)

        !___prepare ocean step__________________________________________________
        call before_oce_step(f%dynamics, f%tracers, f%partit, f%mesh)

        !___model ocean step (CORE COMPUTATION)_________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_timestep_ale'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("oce_timestep_ale")
#endif
        call oce_timestep_ale(n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("oce_timestep_ale")
#endif

        f%t3 = MPI_Wtime()
        !___compute energy diagnostics__________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("compute_diagnostics")
#endif
        call compute_diagnostics(1, f%dynamics, f%tracers, f%ice, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("compute_diagnostics")
#endif

        f%t4 = MPI_Wtime()
        !___prepare output______________________________________________________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call output'//achar(27)//'[0m'
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("output")
#endif
        call output (n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("output")
#endif

        f%t5 = MPI_Wtime()
#if defined (FESOM_PROFILING)
        call fesom_profiler_start("restart")
#endif
        call write_initial_conditions(n, nstart, f%total_nsteps, f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
#if defined (FESOM_PROFILING)
        call fesom_profiler_end("restart")
#endif
        f%t6 = MPI_Wtime()

        ! Accumulate timing statistics
        f%rtime_compute_diag  = f%rtime_compute_diag  + f%t4 - f%t3
        f%rtime_write_means   = f%rtime_write_means   + f%t5 - f%t4
        f%rtime_write_restart = f%rtime_write_restart + f%t6 - f%t5

    end do

    f%from_nstep = f%from_nstep+current_nsteps

#if defined (FESOM_PROFILING)
    call fesom_profiler_end("toy_fesom_runloop_total")
#endif

  end subroutine toy_fesom_runloop


  subroutine toy_fesom_finalize()
    use fesom_toy_storage_module
    ! EO parameters
    real(kind=real32) :: mean_rtime(10), max_rtime(10), min_rtime(10)

#if defined (FESOM_PROFILING)
    call fesom_profiler_start("toy_fesom_finalize_total")
#endif

    call finalize_output()
    call finalize_restart()

    !___FINISH MODEL RUN________________________________________________________
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)

    if (f%mype==0) then
       f%t1 = MPI_Wtime()
       f%runtime_alltimesteps = real(f%t1-f%t0,real32)
       write(*,*) 'TOY OCEAN Run is finished'
    endif

    ! Gather runtime statistics
    mean_rtime(1)  = rtime_oce
    mean_rtime(2)  = rtime_oce_mixpres
    mean_rtime(3)  = rtime_oce_dyn
    mean_rtime(4)  = rtime_oce_dynssh
    mean_rtime(5)  = rtime_oce_solvessh
    mean_rtime(6)  = rtime_oce_GMRedi
    mean_rtime(7)  = rtime_oce_solvetra
    mean_rtime(8)  = f%rtime_compute_diag
    mean_rtime(9)  = f%rtime_write_means
    mean_rtime(10) = f%rtime_write_restart

    max_rtime(1:10) = mean_rtime(1:10)
    min_rtime(1:10) = mean_rtime(1:10)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 10, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
    mean_rtime(1:10) = mean_rtime(1:10) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  10, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  10, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)

#if defined (FESOM_PROFILING)
    call fesom_profiler_end("toy_fesom_finalize_total")
    call fesom_profiler_report(f%MPI_COMM_FESOM, f%mype)
#endif

    if(f%fesom_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype)

    if (f%mype==0) then
        41 format (a35,a10,2a15)
        42 format (a30,3f15.4)

        write(*,*)
        write(*,*) '======================================================'
        write(*,*) '============ TOY OCEAN MODEL RUNTIME ================='
        print 41, '___MODEL RUNTIME per task [seconds]','_____mean_','___________min_', '___________max_'
        print 42, '  runtime ocean:              ',    mean_rtime(1),     min_rtime(1),      max_rtime(1)
        print 42, '    > runtime oce. mix,pres. :',    mean_rtime(2),     min_rtime(2),      max_rtime(2)
        print 42, '    > runtime oce. dyn. u,v,w:',    mean_rtime(3),     min_rtime(3),      max_rtime(3)
        print 42, '    > runtime oce. dyn. ssh  :',    mean_rtime(4),     min_rtime(4),      max_rtime(4)
        print 42, '    > runtime oce. solve ssh :',    mean_rtime(5),     min_rtime(5),      max_rtime(5)
        print 42, '    > runtime oce. GM/Redi   :',    mean_rtime(6),     min_rtime(6),      max_rtime(6)
        print 42, '    > runtime oce. tracer    :',    mean_rtime(7),     min_rtime(7),      max_rtime(7)
        print 42, '  runtime diag:               ',    mean_rtime(8),     min_rtime(8),      max_rtime(8)
        print 42, '  runtime output:             ',    mean_rtime(9),     min_rtime(9),      max_rtime(9)
        print 42, '  runtime restart:            ',    mean_rtime(10),    min_rtime(10),     max_rtime(10)

        43 format (a33,i15)
        45 format (a33,f15.4,a4)

        write(*,*)
        write(*,*) '================ BENCHMARK RUNTIME ==================='
        print 43, '    Number of cores :            ',f%npes
        print 45, '    Runtime for all timesteps :  ',f%runtime_alltimesteps,' sec'
        write(*,*) '======================================================'
        write(*,*)
    end if

  end subroutine toy_fesom_finalize

end module fesom_toy_module
