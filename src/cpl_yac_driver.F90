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
  integer, parameter         :: nsend = 4
  integer, parameter         :: nrecv = 12
  character(len=32)          :: cpl_send(4), cpl_recv(12)
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
    CALL yac_fdef_comp(comp_name, comp_id)
    CALL yac_fread_config_yaml ("coupling.yaml")

    CALL yac_fget_comp_comm(comp_id, localCommunicator)

  end subroutine cpl_yac_init

  subroutine cpl_yac_define_unstr(partit, mesh)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_rotate_grid
    use g_config, only: dt

    implicit none

    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    integer :: ierr, elnodes(3), elem, k, i, j, dimsize, nbr_connections
    real(kind=WP) :: ax(3), amin, x, y
    real(kind=WP), allocatable :: elem_x(:), elem_y(:)
    character(len=4)           :: dt_str
    integer, allocatable :: cell_to_vertex(:)

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    allocate(elem_x(myDim_elem2D))
    allocate(elem_y(myDim_elem2D))
    DO elem=1,myDim_elem2D
       elnodes=mesh%elem2D_nodes(:,elem)
       ax=mesh%coord_nod2D(1, elnodes)
       amin=minval(ax)
       DO k=1,3
          if(ax(k)-amin>=cyclic_length/2.0_WP) ax(k)=ax(k)-cyclic_length
          if(ax(k)-amin<-cyclic_length/2.0_WP) ax(k)=ax(k)+cyclic_length
       END DO
       x=sum(ax)/3.0_WP
       y=sum(mesh%coord_nod2D(2,elnodes))/3.0_WP
       CALL r2g(elem_x(elem), elem_y(elem), x, y)
    END DO

    nbr_connections = SUM(nod_in_elem2D_num(1:myDim_nod2D))

    ALLOCATE(cell_to_vertex(nbr_connections))
    k = 1
    DO i = 1, myDim_nod2D
       DO j = 1, nod_in_elem2D_num(i)
          cell_to_vertex(k) = nod_in_elem2D(j, i)
          k = k + 1
       END DO
    END DO

    CALL yac_fdef_grid_nonuniform_r8("fesom_grid", &
         partit%myDim_elem2D, &  ! integer
         partit%myDim_nod2D, &   ! integer
         nbr_connections, &
         nod_in_elem2D_num(1:partit%myDim_nod2D), &
         elem_x, &
         elem_y, &
         cell_to_vertex, &
         grid_id)

    CALL yac_fdef_points(grid_id, &
         partit%myDim_nod2D, &
         YAC_LOCATION_CELL, &
         geo_coord_nod2D(1, 1:partit%myDim_nod2D), &
         geo_coord_nod2D(2, 1:partit%myDim_nod2D), &
         points_id)

    cpl_send( 1)='sst_feom' ! 1. sea surface temperature [°C]      ->
    cpl_send( 2)='sit_feom' ! 2. sea ice thickness [m]             ->
    cpl_send( 3)='sie_feom' ! 3. sea ice extent [%-100]            ->
    cpl_send( 4)='snt_feom' ! 4. snow thickness [m]                ->

    WRITE (dt_str, '(I4)') int(dt)
    DO i=1,SIZE(cpl_send)
       CALL yac_fdef_field(cpl_send(i), comp_id, [points_id], 1, 1, &
            dt_str, YAC_TIME_UNIT_SECOND, send_field_id(i))
    END DO

    cpl_recv(1)  = 'taux_oce'
    cpl_recv(2)  = 'tauy_oce'
    cpl_recv(3)  = 'taux_ico'
    cpl_recv(4)  = 'tauy_ico'    
    cpl_recv(5)  = 'prec_oce'
    cpl_recv(6)  = 'snow_oce'    
    cpl_recv(7)  = 'evap_oce'
    cpl_recv(8)  = 'subl_oce'
    cpl_recv(9)  = 'heat_oce'
    cpl_recv(10) = 'heat_ico'
    cpl_recv(11) = 'heat_swo'    
    cpl_recv(12) = 'hydr_oce'

    DO i=1,SIZE(cpl_recv)
       CALL yac_fdef_field(cpl_recv(i), comp_id, [points_id], 1, 1, &
            dt_str, YAC_TIME_UNIT_SECOND, recv_field_id(i))
    END DO

    CALL yac_fsearch(ierr)

  end subroutine cpl_yac_define_unstr


  subroutine cpl_yac_send(ind, data_array, action, partit)
    USE MOD_PARTIT
    
    implicit none

    integer, intent( IN )          :: ind       ! variable Id
    logical, intent( OUT )         :: action    !
    type(t_partit), intent(in)     :: partit
    real(kind=WP),  intent(IN)     :: data_array(partit%myDim_nod2D+partit%eDim_nod2D)
    integer :: info, ierr
    
    call yac_fput(send_field_id(1), SIZE(data_array), 1, 1, RESHAPE(data_array, [1, 1, SIZE(data_array)]), info, ierr)
    action = info .eq. YAC_ACTION_COUPLING
  end subroutine cpl_yac_send

  subroutine cpl_yac_recv(ind, data_array, action, partit)
    USE MOD_PARTIT

    implicit none

    integer, intent( IN )  :: ind       ! variable Id
    logical, intent( OUT ) :: action    ! 
    real(kind=WP), intent( OUT ), TARGET    :: data_array(:)
    type(t_partit), intent(inout), target :: partit
    integer :: info, ierr
    type(yac_r8_ptr) :: recv_field(1)

    recv_field(1)%p => data_array

    call yac_fget(send_field_id(1), 1, &
         recv_field, info, ierr)
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
