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
USE g_comm
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

!$ACC DATA COPY(mesh, mesh%coriolis_node, mesh%nn_num, mesh%nn_pos) &
!$ACC      COPY(mesh%ssh_stiff, mesh%ssh_stiff%rowptr) &
!$ACC      COPY(mesh%gradient_sca, mesh%metric_factor, mesh%elem_area, mesh%area, mesh%edge2D_in) &
!$ACC      COPY(mesh%elem2D_nodes, mesh%ulevels, mesh%ulevels_nod2d, mesh%edges, mesh%edge_tri) &
!$ACC      COPY(partit, partit%eDim_nod2D, partit%myDim_edge2D) &
!$ACC      COPY(partit%myDim_elem2D, partit%myDim_nod2D, partit%myList_edge2D) &
!$ACC      COPY(ice, ice%data, ice%work, ice%work%fct_massmatrix) &
!$ACC      COPY(ice%delta_min, ice%Tevp_inv, ice%cd_oce_ice) &
!$ACC      COPY(ice%work%fct_tmax, ice%work%fct_tmin) &
!$ACC      COPY(ice%work%fct_fluxes, ice%work%fct_plus, ice%work%fct_minus) &
!$ACC      COPY(ice%work%eps11, ice%work%eps12, ice%work%eps22) &
!$ACC      COPY(ice%work%sigma11, ice%work%sigma12, ice%work%sigma22) &
!$ACC      COPY(ice%work%ice_strength, ice%stress_atmice_x, ice%stress_atmice_y) &
!$ACC      COPY(ice%srfoce_ssh, ice%thermo%rhosno, ice%thermo%rhoice, ice%pstar, ice%c_pressure) &
!$ACC      COPY(ice%work%inv_areamass, ice%work%inv_mass, ice%uice_rhs, ice%vice_rhs) &
!$ACC      COPY(ice%uice, ice%vice, ice%srfoce_u, ice%srfoce_v, ice%uice_old, ice%vice_old) &
!$ACC      COPY(ice%data(1)%values, ice%data(2)%values, ice%data(3)%values) &
!$ACC      COPY(ice%data(1)%valuesl, ice%data(2)%valuesl, ice%data(3)%valuesl) &
!$ACC      COPY(ice%data(1)%dvalues, ice%data(2)%dvalues, ice%data(3)%dvalues) &
!$ACC      COPY(ice%data(1)%values_rhs, ice%data(2)%values_rhs, ice%data(3)%values_rhs) &
!$ACC      COPY(ice%data(1)%values_div_rhs, ice%data(2)%values_div_rhs, ice%data(3)%values_div_rhs)

!_______________________________________________________________________________
do i=1, 10
    if (partit%mype==0) write(*,*) i
    !___________________________________________________________________________
    ! Dynamics
    select case (ice%whichEVP)
        case (0)
            if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'
            call EVPdynamics  (ice, partit, mesh)
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
    call ice_TG_rhs_div    (ice, partit, mesh)

    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m'
    call ice_fct_solve     (ice, partit, mesh)

    if (partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
    call ice_update_for_div(ice, partit, mesh)

end do

!$ACC END DATA

meta=trim(resultpath)//"/fesom_bin_restart/meta.time"

!_______________________________________________________________________________
! write derived type binary restart files
call write_all_bin_restarts(time, npepath, meta, partit=partit, mesh=mesh, ice=ice)

call par_ex(partit%MPI_COMM_FESOM, partit%mype)

end program main
