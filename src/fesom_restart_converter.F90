!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!                   Raw Restart to NetCDF Converter
!
!=============================================================================!
!  This utility converts raw (binary) restart files to portable NetCDF format.
!  It allows continuing a simulation with a different number of MPI tasks.
!  
!  Usage:
!    1. Run this tool with the SAME number of MPI tasks as the raw restart
!    2. It will read raw restart files and write NetCDF restart files
!    3. Continue your simulation with any number of MPI tasks using NetCDF
!
!  Example:
!    # Original simulation used 288 tasks and created raw restart
!    mpirun -n 288 ./fesom_restart_converter.x
!    
!    # Now continue with 144 tasks using the generated NetCDF restart
!    mpirun -n 144 ./fesom.x
!=============================================================================!    

program fesom_restart_converter
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use MOD_DYN
  use MOD_TRACER
  use MOD_ICE
  use o_PARAM
  use o_ARRAYS
  use g_config
  use g_clock
  use io_RESTART
#if defined(__recom)
  use recom_glovar
  use recom_config
#endif
  use, intrinsic :: iso_fortran_env, only : real32

  implicit none
  
  integer           :: provided, MPIerr
  logical           :: mpi_is_initialized
  real(kind=WP)     :: t0, t1, t2, t3
  type(t_mesh)      :: mesh
  type(t_partit)    :: partit
  type(t_dyn)       :: dynamics
  type(t_tracer)    :: tracers
  type(t_ice)       :: ice
  integer           :: which_readr
  integer           :: expected_npes, actual_npes
  character(500)    :: raw_restart_infopath
  logical           :: restart_exists
  
  ! Initialize MPI
  mpi_is_initialized = .false.
  call MPI_Initialized(mpi_is_initialized, MPIerr)
  if(.not. mpi_is_initialized) then
    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, MPIerr)
  end if
  
  ! Set up partitioning
  call par_init(partit)
  
  if (partit%mype == 0) then
    write(*,*) '============================================================'
    write(*,*) 'FESOM2 Raw Restart to NetCDF Converter'
    write(*,*) '============================================================'
    write(*,*) ''
    t0 = MPI_Wtime()
  endif
  
  ! Setup model configuration (read namelists)
  if (partit%mype == 0) then
    write(*,*) 'Step 1: Reading configuration...'
    t1 = MPI_Wtime()
  endif
  call setup_model(partit)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    write(*,*) '  Configuration loaded in ', real(t2-t1, real32), ' seconds'
    write(*,*) '  RunID: ', trim(runid)
    write(*,*) '  Input path: ', trim(RestartInPath)
    write(*,*) '  Output path: ', trim(RestartOutPath)
    write(*,*) ''
  endif
  
  ! Validate MPI task count matches raw restart
  if (partit%mype == 0) then
    write(*,*) 'Step 2: Validating MPI task count...'
    write(raw_restart_infopath, '(3A,I0,A)') trim(RestartInPath), &
          trim(runid), '_raw_restart/np', partit%npes, '.info'
    inquire(file=trim(raw_restart_infopath), exist=restart_exists)
    
    if (.not. restart_exists) then
      write(*,*) ''
      write(*,*) 'ERROR: Raw restart file not found!'
      write(*,*) '  Looking for: ', trim(raw_restart_infopath)
      write(*,*) ''
      write(*,*) 'Make sure you run this converter with the SAME number of'
      write(*,*) 'MPI tasks as the original simulation that created the raw restart.'
      write(*,*) ''
      write(*,*) 'Current MPI tasks: ', partit%npes
      write(*,*) ''
      call MPI_ABORT(MPI_COMM_WORLD, 1, MPIerr)
    endif
    
    write(*,*) '  Raw restart found for np=', partit%npes
    write(*,*) '  MPI task count matches: OK'
    write(*,*) ''
  endif
  
  ! Setup mesh
  if (partit%mype == 0) then
    write(*,*) 'Step 3: Setting up mesh and domain decomposition...'
    t1 = MPI_Wtime()
  endif
  call mesh_setup(partit, mesh)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    write(*,*) '  Mesh setup complete in ', real(t2-t1, real32), ' seconds'
    write(*,*) '  Total nodes: ', mesh%nod2D
    write(*,*) '  Total elements: ', mesh%elem2D
    write(*,*) '  Vertical levels: ', mesh%nl
    write(*,*) ''
  endif
  
  ! Initialize dynamics and tracers
  if (partit%mype == 0) then
    write(*,*) 'Step 4: Initializing data structures...'
    t1 = MPI_Wtime()
  endif
  
  call dynamics_init(dynamics, partit, mesh)
  call tracer_init(tracers, partit, mesh)
  call arrays_init(tracers%num_tracers, partit, mesh)
  
  if (use_ice) then
    call ice_init(ice, partit, mesh)
  endif
  
#if defined(__recom)
  if (REcoM_restart .or. use_REcoM) then
    call recom_init(tracers, partit, mesh)
  endif
#endif
  
  ! Complete ocean setup to ensure all arrays are properly initialized
  call ocean_setup(dynamics, tracers, partit, mesh)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    write(*,*) '  Data structures initialized in ', real(t2-t1, real32), ' seconds'
    write(*,*) '  Number of tracers: ', tracers%num_tracers
    write(*,*) '  Ice module: ', use_ice
    write(*,*) ''
  endif
  
  ! Initialize clock to get correct year for output filename
  if (partit%mype == 0) then
    write(*,*) 'Step 5: Initializing clock...'
  endif
  call clock_init(partit)
  
  if (partit%mype == 0) then
    write(*,*) '  Year from clock: ', yearnew
    write(*,*) ''
  endif
  
  ! Read raw restart files
  if (partit%mype == 0) then
    write(*,*) 'Step 6: Reading raw restart files...'
    write(*,*) '  This may take several minutes for large simulations...'
    t1 = MPI_Wtime()
  endif
  
  call read_initial_conditions(which_readr, ice, dynamics, tracers, partit, mesh)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    if (which_readr == 1) then
      write(*,*) '  Raw restart files read successfully in ', real(t2-t1, real32), ' seconds'
    else
      write(*,*) ''
      write(*,*) 'WARNING: Expected to read raw restart but got format ', which_readr
      write(*,*) 'Continuing anyway, but output may not be what you expected.'
    endif
    write(*,*) ''
  endif
  
  ! Write NetCDF restart files
  if (partit%mype == 0) then
    write(*,*) 'Step 7: Writing NetCDF restart files...'
    write(*,*) '  This may take several minutes for large simulations...'
    t1 = MPI_Wtime()
  endif
  
  ! Force write of portable restart by calling write_initial_conditions with istep=1
  ! The function will write NetCDF output since we're explicitly triggering it
  call write_initial_conditions(1, 1, 1, which_readr, ice, dynamics, tracers, partit, mesh)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    write(*,*) '  NetCDF restart files written in ', real(t2-t1, real32), ' seconds'
    write(*,*) ''
  endif
  
  ! Finalize restart I/O (close files, join threads)
  if (partit%mype == 0) then
    write(*,*) 'Step 8: Finalizing and cleaning up...'
    t1 = MPI_Wtime()
  endif
  
  call finalize_restart()
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    t3 = MPI_Wtime()
    write(*,*) '  Cleanup complete in ', real(t2-t1, real32), ' seconds'
    write(*,*) ''
    write(*,*) '============================================================'
    write(*,*) 'CONVERSION COMPLETE!'
    write(*,*) '============================================================'
    write(*,*) ''
    write(*,*) 'Output files written to: ', trim(RestartOutPath)
    write(*,*) '  Ocean: ', trim(runid), '.', yearnew, '.oce.restart.nc/'
    if (use_ice) then
      write(*,*) '  Ice:   ', trim(runid), '.', yearnew, '.ice.restart.nc/'
    endif
#if defined(__recom)
    if (REcoM_restart .or. use_REcoM) then
      write(*,*) '  Bio:   ', trim(runid), '.', yearnew, '.bio.restart.nc/'
    endif
#endif
    write(*,*) ''
    write(*,*) 'You can now restart FESOM with a different number of MPI tasks'
    write(*,*) 'using these portable NetCDF restart files.'
    write(*,*) ''
    write(*,*) 'Total runtime: ', real(t3-t0, real32), ' seconds'
    write(*,*) '============================================================'
  endif
  
  ! Finalize MPI
  call MPI_FINALIZE(MPIerr)
  
end program fesom_restart_converter
