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

character(LEN=500)               :: resultpath, npepath, meta
character(LEN=256)               :: npes_string
logical                          :: dir_exist
logical                          :: L_EXISTS
type(t_mesh),       target, save :: mesh
type(t_tracer),     target, save :: tracers
type(t_partit),     target, save :: partit
type(t_dyn),        target, save :: dyn
type(t_ice),        target, save :: ice
integer                          :: i, n, nz, nzmax, nzmin
integer, dimension(3)            :: time=(/0, 0, 0/)

call MPI_INIT(i)
call par_init(partit)

resultpath='../'
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
call read_all_bin_restarts(npepath, partit=partit, mesh=mesh, dynamics=dyn, tracers=tracers)
! even though the partitioning has been read some things regarding MPI shall be computed during the runtime
! these include: MPI_TYPE_COMMIT etc.
! used to be call set_par_support(partit, mesh)
call init_mpi_types(partit, mesh)
call init_gatherLists(partit)

    !$ACC DATA COPY(mesh, partit, tracers, dyn) &
    !$ACC      COPY(mesh%helem, mesh%elem_cos, mesh%edge_cross_dxdy, mesh%elem2d_nodes, mesh%nl) &
    !$ACC      COPY(mesh%nlevels_nod2D, mesh%ulevels_nod2D, mesh%nod_in_elem2D, mesh%nod_in_elem2D_num) &
    !$ACC      COPY(mesh%ulevels, mesh%edge_tri, mesh%edge_dxdy, mesh%edges, mesh%nlevels, mesh%hnode, mesh%hnode_new) &
    !$ACC      COPY(mesh%zbar_3d_n, mesh%z_3d_n, mesh%area, mesh%areasvol) &
    !$ACC      COPY(partit%mydim_elem2d, partit%myDim_nod2D, partit%edim_nod2d, partit%myDim_edge2D) &
    !$ACC      COPY(dyn%w, dyn%w_e, dyn%uv) &
    !$ACC      COPY(tracers%data(1), tracers%work) &
    !$ACC      COPY(tracers%data(1)%values, tracers%data(1)%valuesAB) &
    !$ACC      COPY(tracers%work%fct_ttf_min, tracers%work%fct_ttf_max, tracers%work%fct_plus, tracers%work%fct_minus) &
    !$ACC      COPY(tracers%work%adv_flux_hor, tracers%work%adv_flux_ver, tracers%work%fct_LO) &
    !$ACC      COPY(tracers%work%del_ttf_advvert, tracers%work%del_ttf_advhoriz, tracers%work%edge_up_dn_grad) &
    !$ACC      COPY(tracers%work%del_ttf)

do i=1, 10
    !___________________________________________________________________________
    ! ale tracer advection
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, partit%myDim_nod2D + partit%eDim_nod2D
        do nz = 1, mesh%nl - 1
        tracers%work%del_ttf_advhoriz(nz, n) = 0.0_WP
        tracers%work%del_ttf_advvert(nz, n)  = 0.0_WP
        end do
    end do
    !$ACC END PARALLEL LOOP

    call do_oce_adv_tra(1.e-3, dyn%uv, dyn%w, dyn%w_i, dyn%w_e, 1, dyn, tracers, partit, mesh)

    if ( partit%mype == 0 ) write(*,*) minval(tracers%data(1)%values), maxval(tracers%data(1)%values), sum(tracers%data(1)%values)

    !_____________________________________________________
    !___________________________________________________________________________
    ! update array for total tracer flux del_ttf with the fluxes from horizontal
    ! and vertical advection
    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n = 1, partit%myDim_nod2D + partit%eDim_nod2D
        nzmax=mesh%nlevels_nod2D(n)-1
        nzmin=mesh%ulevels_nod2D(n)
        !$ACC LOOP VECTOR
        do nz = nzmin, nzmax
            tracers%work%del_ttf(nz, n) = tracers%work%del_ttf(nz, n) + tracers%work%del_ttf_advhoriz(nz, n) + tracers%work%del_ttf_advvert(nz, n)
        end do
    end do
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
    do n=1, partit%myDim_nod2D
        nzmax=mesh%nlevels_nod2D(n)-1
        nzmin=mesh%ulevels_nod2D(n)
        !$ACC LOOP VECTOR
        do nz = nzmin, nzmax
            tracers%data(1)%values(nz ,n) = tracers%data(1)%values(nz ,n) + tracers%work%del_ttf(nz ,n) / mesh%hnode_new(nz ,n) ! LINFS
        end do
    end do
    !$ACC END PARALLEL LOOP

    call exchange_nod(tracers%data(1)%values(:,:), partit, luse_g2g = .true.)
end do

!$ACC END DATA

meta=trim(resultpath)//"/fesom_bin_restart/meta.time"
!_______________________________________________________________________________
! write derived type binary restart files
call write_all_bin_restarts(time, npepath, meta, partit=partit, mesh=mesh, dynamics=dyn, tracers=tracers)

call par_ex(partit%MPI_COMM_FESOM, partit%mype)
end program main
