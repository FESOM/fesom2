#if defined(__yac)
module cpl_yac_driver

  use mo_yac_finterface
  USE o_PARAM
  implicit none
  save

  character(len=*), PARAMETER   :: comp_name = "fesom2"
  integer :: comp_id, grid_id, points_id
  integer :: send_field_id(4), recv_field_id(12)
  real(kind=WP), dimension(:,:),   allocatable   :: a2o_fcorr_stat  !flux correction statistics for the output
  integer, parameter         :: nsend = 2
  integer, parameter         :: nrecv = 5
  character(len=32)          :: cpl_send(nsend), cpl_recv(nrecv)
  integer                    :: cpl_send_collection_size(nsend), cpl_recv_collection_size(nrecv)
  integer                    :: source_root, target_root   !this root/source in MPI_COMM_WORLD
  logical                    :: commRank       ! true for ranks doing YAC communication

  public nsend, nrecv
  public cpl_send, cpl_recv
  public source_root, target_root, commRank
  public a2o_fcorr_stat

contains

  subroutine cpl_yac_init( localCommunicator )
    implicit none

    integer, intent(out) :: localCommunicator

#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_init : coupler initialization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */

    CALL yac_finit()
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fread_config_yaml ("coupling.yaml")
    CALL yac_fdef_comp(comp_name, comp_id)

    CALL yac_fget_comp_comm(comp_id, localCommunicator)

  end subroutine cpl_yac_init

  subroutine compute_midpoint(geo_a, geo_b, mid)
    implicit none

    REAL(kind=WP), intent(in) :: geo_a(2), geo_b(2)
    REAL(kind=WP), intent(out) :: mid(2)
    REAL(kind=WP) :: cos_lat_a, cos_lat_b, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    x_mid = 0.5*(cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)))
    y_mid = 0.5*(cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)))
    z_mid = 0.5*(sin(geo_a(2)) + sin(geo_b(2)))

    mid(1) = atan2(y_mid , x_mid)
    mid(2) = PI/2 - acos(z_mid/sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_midpoint

  subroutine compute_center(geo_a, geo_b, geo_c, mid)
    implicit none

    REAL(kind=WP), intent(in) :: geo_a(2), geo_b(2), geo_c(2)
    REAL(kind=WP), intent(out) :: mid(2)
    REAL(kind=WP) :: cos_lat_a, cos_lat_b, cos_lat_c, x_mid, y_mid, z_mid

    cos_lat_a = cos(geo_a(2))
    cos_lat_b = cos(geo_b(2))
    cos_lat_c= cos(geo_c(2))
    x_mid = (cos_lat_a * cos(geo_a(1)) + cos_lat_b * cos(geo_b(1)) + cos_lat_c * cos(geo_c(1))) / 3.
    y_mid = (cos_lat_a * sin(geo_a(1)) + cos_lat_b * sin(geo_b(1)) + cos_lat_c * sin(geo_c(1)) ) / 3.
    z_mid = (sin(geo_a(2)) + sin(geo_b(2)) + sin(geo_c(2))) /3.

    mid(1) = atan2(y_mid , x_mid)
    mid(2) = PI/2 - acos(z_mid/sqrt(x_mid*x_mid + y_mid*y_mid + z_mid*z_mid))
  end subroutine compute_center

  subroutine cpl_yac_define_unstr(partit, mesh)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_rotate_grid
    use g_config, only: dt

    implicit none

    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP), allocatable :: x_vertices(:), y_vertices(:)
    real(kind=WP) :: mid(2)
    integer, allocatable :: nbr_vertices_per_cell(:), cell_to_vertex(:)
    integer :: ierr, i, j, k, nbr_vertices, nbr_boundary_nodes, nbr_connections, vtx_idx, c2v_idx
    integer :: curr_elem, curr_edge
    logical, allocatable :: node_is_boundary(:)
    character(len=4)           :: dt_str

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    ! find boundary nodes
    ALLOCATE(node_is_boundary(myDim_nod2D))
    node_is_boundary = .FALSE.
    DO i=1,myDim_edge2D
       IF (myList_edge2D(i) > edge2D_in) THEN
          IF (edges(1,i) <= myDim_nod2D) &
               node_is_boundary(edges(1,i)) = .TRUE.
          IF (edges(2,i) <= myDim_nod2D) &
               node_is_boundary(edges(2,i)) = .TRUE.
       END IF
    END DO
    nbr_boundary_nodes = COUNT(node_is_boundary)

    nbr_vertices = myDim_elem2D + myDim_edge2D + nbr_boundary_nodes
    nbr_connections = 2*SUM(nod_in_elem2D_num(1:myDim_nod2D)) + 2*nbr_boundary_nodes

    WRITE (0,*) "nbr_boundary_nodes: ", nbr_boundary_nodes
    WRITE (0,*) "nbr_vertices: " , nbr_vertices
    WRITE (0,*) "nbr_cells: " , myDim_nod2D
    WRITE (0,*) "nbr_connections: " , nbr_connections

    ALLOCATE(x_vertices(nbr_vertices))
    ALLOCATE(y_vertices(nbr_vertices))
    ALLOCATE(nbr_vertices_per_cell(myDim_nod2D))
    ALLOCATE(cell_to_vertex(nbr_connections))

    ! compute vertices
    ! 1 element centers
    DO i=1,myDim_elem2D
       CALL compute_center(geo_coord_nod2D(1:2,elem2D_nodes(1, i)), &
            geo_coord_nod2D(1:2,elem2D_nodes(2, i)), &
            geo_coord_nod2D(1:2,elem2D_nodes(3, i)), mid)
       x_vertices(i) = mid(1)
       y_vertices(i) = mid(2)
    END DO
    ! 2 edges midpoints
    DO i=1,myDim_edge2D
       CALL compute_midpoint(geo_coord_nod2D(1:2,edges(1,i)), geo_coord_nod2D(1:2,edges(2,i)), mid)
       x_vertices(myDim_elem2D + i) = mid(1)
       y_vertices(myDim_elem2D + i) = mid(2)
    END DO
    ! 3 boundary nodes -> are added on the fly when the cells are computed

    ! compute cells
    vtx_idx = myDim_elem2D + myDim_edge2D
    c2v_idx = 0
    DO i=1,myDim_nod2D
       nbr_vertices_per_cell(i) = 2*nod_in_elem2D_num(i)
       curr_elem = nod_in_elem2D(1,i)
       curr_edge = elem_edges(1, curr_elem)
       IF (ALL(edges(1:2,curr_edge) /= i)) curr_edge = elem_edges(2, curr_elem)
       IF (node_is_boundary(i)) THEN
          nbr_vertices_per_cell(i) = nbr_vertices_per_cell(i) + 2
          ! we're starting with the boundary node
          vtx_idx = vtx_idx+1
          x_vertices(vtx_idx) = geo_coord_nod2D(1,i)
          y_vertices(vtx_idx) = geo_coord_nod2D(2,i)
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = vtx_idx
          elem_loop: DO j=1,nod_in_elem2D_num(i)
             DO k=1,3
                IF (myList_edge2D(elem_edges(k,nod_in_elem2D(j,i))) > edge2D_in .AND. &
                     ANY(edges(1:2,elem_edges(k,nod_in_elem2D(j,i))) == i)) THEN
                   curr_elem = nod_in_elem2D(j,i)
                   curr_edge = elem_edges(k, curr_elem)
                   EXIT elem_loop
                END IF
             END DO
          END DO elem_loop
       END IF

       DO j=1,nod_in_elem2D_num(i)
          IF (curr_elem == 0) THEN
             WRITE (0,*) "Error: cells with two or more coast lines are currently not supported by the yac setup code"
             STOP
          ENDIF
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = myDim_elem2D + curr_edge
          ! add the center of curr_elem
          c2v_idx = c2v_idx + 1
          IF (curr_elem < 1) THEN
             WRITE (0,*) "elem 0 detected in cell ", i, " with ", nbr_vertices_per_cell(i), " vertices - j is ", j, " is_boundary ", node_is_boundary(i)
          END IF
          cell_to_vertex(c2v_idx) = curr_elem
          ! find next edge
          edge_loop: DO k=1,3
             IF (elem_edges(k, curr_elem) /= curr_edge .AND. &
                  ANY(edges(1:2,elem_edges(k, curr_elem)) == i)) THEN
                curr_edge = elem_edges(k, curr_elem)
                EXIT edge_loop
             END IF
          END DO edge_loop
          IF (edge_tri(1, curr_edge) /= curr_elem) THEN
             curr_elem = edge_tri(1, curr_edge)
          ELSE
             curr_elem = edge_tri(2, curr_edge)
          END IF
       END DO

       ! in case of a boundary volume, we add the others edges midpoint
       IF (node_is_boundary(i)) THEN
          ! add the midpoint of curr_edge
          c2v_idx = c2v_idx + 1
          cell_to_vertex(c2v_idx) = myDim_elem2D + curr_edge
       END IF
    END DO

    CALL yac_fdef_grid_nonuniform_r8( &
         "fesom_grid", & ! grid_name
         nbr_vertices, &   ! nbr_vertices
         myDim_nod2D, &  ! nbr_cells
         nbr_connections, & ! nbr_connections
         nbr_vertices_per_cell, & ! nbr_vertices_per_cell
         x_vertices, & ! x_vertices
         y_vertices, & ! y_vertices
         cell_to_vertex, &
         grid_id)

    CALL yac_fdef_points(grid_id, &
         myDim_nod2D, &
         YAC_LOCATION_CELL, &
         geo_coord_nod2D(1,1:myDim_nod2D), &
         geo_coord_nod2D(2,1:myDim_nod2D), &
         points_id)

    cpl_send( 1)='sst_feom' ! 1. sea surface temperature [°C]      ->
    cpl_send_collection_size(1) = 1
    cpl_send( 2)='ocean_sea_ice_bundle'
    cpl_send_collection_size(2) = 3

    WRITE (dt_str, '(I4)') int(dt)
    DO i=1,nsend
       CALL yac_fdef_field(cpl_send(i), comp_id, [points_id], 1, &
            cpl_send_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, send_field_id(i))
    END DO

    cpl_recv(1)  = 'taux'
    cpl_recv_collection_size(1) = 2
    cpl_recv(2)  = 'tauy'
    cpl_recv_collection_size(2) = 2
    cpl_recv(3)  = 'surface_fresh_water_flux'
    cpl_recv_collection_size(3) = 3
    cpl_recv(4) = 'total_heat_flux'
    cpl_recv_collection_size(4) = 4
    cpl_recv(5) = 'atmosphere_sea_ice_bundle'
    cpl_recv_collection_size(5) = 2

    DO i=1,nrecv
       CALL yac_fdef_field(cpl_recv(i), comp_id, [points_id], 1, cpl_recv_collection_size(i), &
            dt_str, YAC_TIME_UNIT_SECOND, recv_field_id(i))
    END DO

    CALL yac_fsearch(ierr)

  end subroutine cpl_yac_define_unstr


  subroutine cpl_yac_send(ind, data_array, action)
    implicit none

    integer, intent( IN )          :: ind       ! variable Id
    logical, intent( OUT )         :: action    !
    real(kind=WP),  intent(IN)     :: data_array(:, :)
    integer :: info, ierr

    call yac_fput(send_field_id(ind), SIZE(data_array, 1), SIZE(data_array, 2), &
         data_array, info, ierr)
    action = info .eq. YAC_ACTION_COUPLING
  end subroutine cpl_yac_send

  subroutine cpl_yac_recv(ind, data_array, action)
    implicit none

    integer, intent( IN )  :: ind       ! variable Id
    logical, intent( OUT ) :: action    ! 
    real(kind=WP), intent( INOUT )    :: data_array(:,:)
    integer :: info, ierr

    call yac_fget(recv_field_id(ind), SIZE(data_array, 1), SIZE(data_array, 2),&
         data_array, info, ierr)
    action = info .eq. YAC_ACTION_COUPLING
  end subroutine cpl_yac_recv

  subroutine cpl_yac_finalize ()
    implicit none
#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_finalize : coupler finalization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */
    CALL yac_ffinalize()
  end subroutine cpl_yac_finalize

end module cpl_yac_driver
#endif
