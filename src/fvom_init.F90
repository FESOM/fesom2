!> @author
!> Sergei Danilov
!> @brief
!> fvom_init.F90
!! A set of routines that does domain partitioning for a distributed setup.
!! Cell-vertex finite-volume discretization
!! routines used here are the stripped down versions of the main setup.
!! Do no mix them up!

!=============================================================================
!> @brief
!> Main driver routine for initialization
program MAIN

  use o_PARAM, only: MAX_PATH
  use MOD_MESH, only: t_mesh
  use MOD_PARTIT, only: t_partit
  use g_CONFIG, only: paths, geometry, run_config, machine, cyclic_length, rad, &
                      alphaEuler, betaEuler, gammaEuler, use_cavity
  use g_rotate_grid, only: set_mesh_transform_matrix
  use iso_fortran_env, only: error_unit
  use mod_mesh_utils, only: read_mesh_ini, read_mesh_cavity, test_tri_ini, &
                            find_edges_ini, find_elem_neighbors_ini, find_levels, &
                            find_levels_cavity, stiff_mat_ini, set_par_support_ini, &
                            communication_ini
  
  implicit none


  character(len=MAX_PATH)         :: nmlfile  !> name of configuration namelist file
  integer                     :: start_t, interm_t, finish_t, rate_t
  integer                     :: file_unit
  integer                     :: ierr
  character(512)              :: errmsg
  type(t_mesh),   target, save :: mesh
  type(t_partit), target, save :: partit

  call system_clock(start_t, rate_t)
  interm_t = start_t
  
  nmlfile ='namelist.config'

  open(newunit=file_unit, file=nmlfile, action='read', &
      status='old', iostat=ierr, iomsg=errmsg)
  if (ierr /= 0) then
    write (unit=error_unit, fmt='(3A)') &
      '### error: can not open file ', trim(nmlfile), &
      ', error: ' // trim(errmsg)
    error stop
  end if

  read (file_unit, NML=paths)         ! We need MeshPath
  read (file_unit, NML=geometry)      ! We need cyclic_length and cartesian
  read (file_unit, NML=run_config)    ! We need use_cavity=true/false
  read (file_unit, NML=machine)       ! We need partitioning hierarchy
  close (file_unit)

  cyclic_length=cyclic_length*rad 
  alphaEuler=alphaEuler*rad 	
  betaEuler=betaEuler*rad
  gammaEuler=gammaEuler*rad
  call set_mesh_transform_matrix  !(rotated grid) 
  call read_mesh_ini(mesh)
  if (use_cavity) call read_mesh_cavity(mesh)
  
  call system_clock(finish_t)
  print '("**** Reading initial mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call test_tri_ini(mesh)
  call find_edges_ini(mesh)
  call find_elem_neighbors_ini(mesh)
  call find_levels(mesh)
  if (use_cavity) call find_levels_cavity(mesh)
    
! NR Some arrays are not needed for partitioning, after setting up the grid
  deallocate(mesh%coord_nod2D)
  deallocate(mesh%edge_tri)
  deallocate(mesh%zbar,mesh%Z,mesh%nlevels,mesh%depth)
  call system_clock(finish_t)
  print '("**** Checking and initializing mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t

  call stiff_mat_ini(mesh)
  call set_par_support_ini(partit, mesh)
  call system_clock(finish_t)
  print '("**** Partitioning time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call communication_ini(partit, mesh)
  call system_clock(finish_t)
  print '("**** Storing partitioned mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  print '("**** Total time = ",f12.3," seconds. ****")', &
       real(finish_t-start_t)/real(rate_t)
end program MAIN
