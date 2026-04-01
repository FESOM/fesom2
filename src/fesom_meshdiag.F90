!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!                     Mesh Diagnostics Utility
!
!=============================================================================!
!  This utility generates mesh diagnostics NetCDF file and exits.
!  It performs minimal initialization required for mesh diagnostics generation.
!=============================================================================!    

program fesom_meshdiag
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use MOD_DYN
  use MOD_TRACER
  use o_PARAM
  use o_ARRAYS
  use g_config
  use g_comm_auto
  use io_mesh_info
  use, intrinsic :: iso_fortran_env, only : real32

  implicit none
  
  ! Local variables
  integer           :: provided, MPIerr
  logical           :: mpi_is_initialized
  real(kind=WP)     :: t0, t1, t2
  type(t_mesh)      :: mesh
  type(t_partit)    :: partit
  type(t_dyn)       :: dynamics
  type(t_tracer)    :: tracers
  
  ! Initialize MPI (following FESOM2 approach)
  mpi_is_initialized = .false.
  call MPI_Initialized(mpi_is_initialized, MPIerr)
  if(.not. mpi_is_initialized) then
    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, MPIerr)
  end if
  
  ! Set up global pointers
  call par_init(partit)
  
  if (partit%mype == 0) then
    write(*,*) '=========================================='
    write(*,*) 'FESOM2 Mesh Diagnostics Utility'
    write(*,*) '=========================================='
    t0 = MPI_Wtime()
  endif
  
  ! Setup model (read namelists)
  if (partit%mype == 0) write(*,*) 'Reading namelists...'
  call setup_model(partit)
  
  ! Setup mesh
  if (partit%mype == 0) then
    write(*,*) 'Setting up mesh...'
    t1 = MPI_Wtime()
  endif
  call mesh_setup(partit, mesh)
  
  if (partit%mype == 0) then
    t2 = MPI_Wtime()
    write(*,*) 'FESOM mesh_setup... complete'
    write(*,*) 'Mesh setup time: ', real(t2-t1, real32), ' seconds'
  endif
  
  ! Check mesh consistency
  if (partit%mype == 0) write(*,*) 'Checking mesh consistency...'
  call check_mesh_consistency(partit, mesh)
  
  ! Initialize minimal dynamics and tracers for mesh diagnostics
  if (partit%mype == 0) write(*,*) 'Initializing dynamics and tracers...'
  call dynamics_init(dynamics, partit, mesh)
  call tracer_init(tracers, partit, mesh)
  call arrays_init(tracers%num_tracers, partit, mesh)
  
  ! Complete ocean setup to ensure all arrays are properly initialized
  if (partit%mype == 0) write(*,*) 'Completing ocean setup...'
  call ocean_setup(dynamics, tracers, partit, mesh)
  
  ! Generate mesh diagnostics NetCDF file
  if (partit%mype == 0) write(*,*) 'Writing mesh diagnostics...'
  call write_mesh_info(partit, mesh)
  
  if (partit%mype == 0) then
    write(*,*) 'Mesh diagnostics written successfully!'
    write(*,*) 'Output file: ', trim(ResultPath)//trim(runid)//'.mesh.diag.nc'
    write(*,*) 'Total runtime: ', real(MPI_Wtime()-t0, real32), ' seconds'
    write(*,*) '=========================================='
  endif
  
  ! Finalize MPI
  call MPI_FINALIZE(MPIerr)
  
end program fesom_meshdiag

