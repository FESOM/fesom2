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
USE MOD_ICE
USE MOD_PARSUP
USE g_comm_auto
USE par_support_interfaces
USE restart_derivedtype_module
use fortran_utils
use nvtx
IMPLICIT NONE

character(LEN=500)           :: resultpath, npepath, meta
character(LEN=256)           :: npes_string
logical                      :: dir_exist
type(t_mesh)  , target, save :: mesh
type(t_partit), target, save :: partit
type(t_ice)   , target, save :: ice
integer                      :: i, n, nzmax, nzmin
real(kind=WP) , allocatable  :: UV(:,:,:), wvel(:,:), wvel_i(:,:), wvel_e(:,:)
integer                      :: node_size, elem_size
integer, dimension(3)        :: time=(/0, 0, 0/)

!_______________________________________________________________________________
resultpath='../'

!_______________________________________________________________________________
call MPI_INIT(i)
call par_init(partit)

!_______________________________________________________________________________
! check if resultpath exist
#if defined(__PGI)
INQUIRE(file=trim(resultpath), EXIST=dir_exist)
#else
INQUIRE(directory=trim(resultpath), EXIST=dir_exist)
#endif
if (.not. dir_exist) then
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(resultpath)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if

npepath =trim(resultpath)//"/fesom_bin_restart/np"//int_to_txt(partit%npes)
#if defined(__PGI)
INQUIRE(file=trim(npepath), EXIST=dir_exist)
#else
INQUIRE(directory=trim(npepath), EXIST=dir_exist)
#endif
if (.not. dir_exist) then
    if (partit%mype==0) print *, achar(27)//'[1;31m'//' -ERROR-> could not find:'//trim(npepath)//achar(27)//'[0m'
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    stop
end if

!_______________________________________________________________________________
! read derived type binary restart files
call read_all_bin_restarts(npepath, ice=ice, partit=partit, mesh=mesh)

!_______________________________________________________________________________
! even though the partitioning has been read some things regarding MPI shall be computed during the runtime
! these include: MPI_TYPE_COMMIT etc.
! used to be call set_par_support(partit, mesh)
call init_mpi_types(partit, mesh)
call init_gatherLists(partit)

!_______________________________________________________________________________
node_size=partit%myDim_nod2D +partit%eDim_nod2D
elem_size=partit%myDim_elem2D+partit%eDim_elem2D

!_______________________________________________________________________________
do i=1, 10
    if (partit%mype==0) write(*,*) i
    !___________________________________________________________________________
    ! Dynamics
    select case (ice%whichEVP)
        case (0)
            if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'
            call nvtxStartRange('EVP')
            call EVPdynamics  (ice, partit, mesh)
            call nvtxEndRange()
        case (1)
            if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_m...'//achar(27)//'[0m'
            call EVPdynamics_m(ice, partit, mesh)
        case (2)
            if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_a...'//achar(27)//'[0m'
            call EVPdynamics_a(ice, partit, mesh)
        case default
            if (partit%mype==0) write(*,*) 'a non existing EVP scheme specified!'
            call par_ex(partit%MPI_COMM_FESOM, partit%mype)
            stop
    end select

    !___________________________________________________________________________
    ! Advection
    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_TG_rhs_div...'//achar(27)//'[0m'
    call nvtxStartRange('ice_TG')
    call ice_TG_rhs_div    (ice, partit, mesh)
    call nvtxEndRange()

    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m'
    call nvtxStartRange('ice_fct')
    call ice_fct_solve     (ice, partit, mesh)
    call nvtxEndRange()

    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
    call nvtxStartRange('ice_update')
    call ice_update_for_div(ice, partit, mesh)
    call nvtxEndRange()

end do

meta=trim(resultpath)//"/fesom_bin_restart/meta.time"

!_______________________________________________________________________________
! write derived type binary restart files
call write_all_bin_restarts(time, npepath, meta, partit=partit, mesh=mesh, ice=ice)

call par_ex(partit%MPI_COMM_FESOM, partit%mype)

end program main
