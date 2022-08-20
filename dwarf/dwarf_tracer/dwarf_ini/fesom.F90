!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main
USE MOD_MESH
USE MOD_PARTIT
USE MOD_TRACER
USE MOD_DYN
USE MOD_ICE
USE MOD_PARSUP
USE g_comm_auto
USE par_support_interfaces
USE restart_derivedtype_module
USE fortran_utils
IMPLICIT NONE

character(LEN=500)               :: resultpath, npepath
character(LEN=256)               :: npes_string
logical                          :: dir_exist
logical                          :: L_EXISTS
type(t_mesh),       target, save :: mesh
type(t_tracer),     target, save :: tracers
type(t_partit),     target, save :: partit
type(t_dyn),        target, save :: dyn
type(t_ice),        target, save :: ice
integer                          :: i, n, nzmax, nzmin


call MPI_INIT(i)
call par_init(partit)

resultpath='../'
! check if resultpath exist
INQUIRE(directory=trim(resultpath), EXIST=dir_exist)
if (.not. dir_exist) then 
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(resultpath)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if 

npepath =trim(resultpath)//"/fesom_bin_restart/np"//int_to_txt(partit%npes)
INQUIRE(directory=trim(npepath), EXIST=dir_exist)
if (.not. dir_exist) then 
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(npepath)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if 

!_______________________________________________________________________________
! read derived type binary restart files
call read_all_bin_restarts(npepath, partit=partit, mesh=mesh, dynamics=dyn, tracers=tracers)
! even though the partitioning has been read some things regarding MPI shall be computed during the runtime
! these include: MPI_TYPE_COMMIT etc.
! used to be call set_par_support(partit, mesh)
call init_mpi_types(partit, mesh)
call init_gatherLists(partit)

do i=1, 10
    !___________________________________________________________________________
    ! ale tracer advection 
    tracers%work%del_ttf_advhoriz = 0.0_WP
    tracers%work%del_ttf_advvert  = 0.0_WP
!   if (mype==0) write(*,*) 'start advection part.......'
    call do_oce_adv_tra(1.e-3, dyn%uv, dyn%w, dyn%w_i, dyn%w_e, 1, dyn, tracers, partit, mesh)
!   if (mype==0) write(*,*) 'advection part completed...'
    if (partit%mype==0) write(*,*) minval(tracers%data(1)%values), maxval(tracers%data(1)%values), sum(tracers%data(1)%values)
    !_____________________________________________________
    !___________________________________________________________________________
    ! update array for total tracer flux del_ttf with the fluxes from horizontal
    ! and vertical advection
    tracers%work%del_ttf=tracers%work%del_ttf+tracers%work%del_ttf_advhoriz+tracers%work%del_ttf_advvert

   do n=1, partit%myDim_nod2D 
      nzmax=mesh%nlevels_nod2D(n)-1
      nzmin=mesh%ulevels_nod2D(n)
      tracers%data(1)%values(nzmin:nzmax,n)=tracers%data(1)%values(nzmin:nzmax,n)+tracers%work%del_ttf(nzmin:nzmax,n)/mesh%hnode_new(nzmin:nzmax,n) ! LINFS
   end do
   call exchange_nod(tracers%data(1)%values(:,:), partit)
   call exchange_nod(tracers%data(2)%values(:,:), partit)
end do
call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program main
