#if defined (__oasis)
module cpl_driver
  !======================================================================
  !
  ! for coupling between the FESOM ocean ECHAM6 atmosphere using OASIS3-MCT
  !
  !=====================================================================
  ! History :
  !  09-09  (R. Redler, Germany)  Original code
  !  09-09  (K.Fieg, AWI Germany) Adjustment for FESOM
  !  07-12  (D.Barbi, AWI Germany) Switch to ECHAM6.1 and OASIS3-MCT
  !----------------------------------------------------------------------
  ! Modules used
  !
  use mod_oasis                    ! oasis module
  use g_config, only : dt, use_icebergs, lwiso
#if defined(__usetp)
  use g_config, only : num_fesom_groups ! kh 03.12.21 OG 08.09.23
#endif
  use o_param,  only : rad
  USE MOD_PARTIT
  implicit none
  save   
  !
  ! Exchange parameters for coupling FESOM with ECHAM6
  !

  !---wiso-code
  ! define nsend and nrecv as variables instead of fixed parameters
  ! (final number of fields depends now on lwiso switch and is set in subroutine cpl_oasis3mct_define_unstr)

#if defined (__oifs)
  integer                    :: nsend = 7
  integer                    :: nrecv = 13
#else
  integer                    :: nsend = 4
  integer                    :: nrecv = 12
#endif
  
  ! define send_id and recv_id with variable dimension as nsend and nrecv are now variables)
  integer, allocatable, dimension(:) :: send_id
  integer, allocatable, dimension(:) :: recv_id

  ! define cpl_send and cpl_recv with variable dimension as nsend and nrecv are now variables)
  character(len=32), allocatable, dimension(:) :: cpl_send
  character(len=32), allocatable, dimension(:) :: cpl_recv

  !---wiso-code-end

  character(len=16)          :: appl_name      ! application name for OASIS use
  character(len=16)          :: comp_name      ! name of this component
  character(len=11)          :: grid_name      ! name of the grid

  private

  integer                    :: source_root, target_root   !this root/source in MPI_COMM_WORLD
  integer, parameter         :: localRoot  = 0

  integer                    :: localRank      ! local MPI rank
  integer                    :: localSize      ! local MPI size
  integer                    :: localComm      ! local MPI size
  logical                    :: commRank       ! true for ranks doing OASIS communication
  integer                    :: comp_id        ! id returned by oasis_init_comp

  logical, save              :: oasis_was_initialized
  logical, save              :: oasis_was_terminated
  integer, save              :: write_grid

  integer, save              :: seconds_til_now=0
  integer                    :: ierror              ! return error code
  logical                    :: rootexchg   =.true. ! logical switch 


  integer                    :: o2a_call_count=0
  integer                    :: a2o_call_count=0

  REAL(kind=WP), POINTER                          :: exfld(:)          ! buffer for receiving global exchange fields
  real(kind=WP), allocatable, dimension(:,:)      :: cplsnd


  integer, dimension(1,3)    :: iextent    
  integer, dimension(1,3)    :: ioffset   

 !
  real(kind=WP), dimension(:,:),   allocatable   :: a2o_fcorr_stat  !flux correction statistics for the output

  !
  ! Routine accessibility
  !
  public cpl_oasis3mct_init
  public cpl_oasis3mct_define_unstr
  public cpl_oasis3mct_send
  public cpl_oasis3mct_recv
  public cpl_oasis3mct_finalize
  public fesom_flush
  public exchange_roots
  public seconds_til_now 
  public send_id, recv_id
  public nsend, nrecv, cpl_send, cpl_recv
  public source_root, target_root, commRank
  public a2o_fcorr_stat

contains

    subroutine node_contours(my_x_corners, my_y_corners, partit, mesh)
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE o_PARAM
        use g_comm_auto
        use o_ARRAYS
        use g_rotate_grid

        IMPLICIT NONE
        type(t_mesh),   intent(in), target :: mesh
        type(t_partit), intent(inout), target :: partit
        real(kind=WP), allocatable, intent(inout) :: my_x_corners(:,:)     ! longitude node corners
        real(kind=WP), allocatable, intent(inout) :: my_y_corners(:,:)     ! latitude node corners    
        integer                               :: bEdge_left, bEdge_right
        integer,              dimension(2)    :: belem_left, belem_right
        integer                               :: edge_left, edge_right
        integer                               :: n, ee, elem, nn, el(2), flag, nn1, nn2
        integer                               :: current_pos
        integer                               :: pos_increment=-1 ! counter clockwise is negative, otherwise +1!
        integer, allocatable, dimension(:)    :: nedges, nelems, nedges1, nelems1, nedges2, nelems2
        real(kind=WP)                         :: this_x_coord, this_y_coord

include "associate_part_def.h"
include "associate_mesh_def.h"
include "associate_part_ass.h"
include "associate_mesh_ass.h"

    if (.not. allocated(my_x_corners)) then
        ALLOCATE(my_x_corners(myDim_nod2D, 25)) !maxval(nod_in_elem2D_num, 1)*2+2))
    endif
    if (.not. allocated(my_y_corners)) then
        ALLOCATE(my_y_corners(myDim_nod2D, 25)) !maxval(nod_in_elem2D_num, 1)*2+2))
    endif
    do n=1, myDim_nod2D
        ! find the type/of node: internal or at boundary
        bEdge_left =0
        belem_left =0
        bEdge_right=0
        belem_right=0

        do ee=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(ee,n)
           if (elem2D_nodes(1,elem)==n) then
              edge_left=elem_edges(3,elem)
              edge_right=elem_edges(2,elem)
           elseif (elem2D_nodes(2,elem)==n) then
              edge_left=elem_edges(1,elem)
              edge_right=elem_edges(3,elem)
           else
              edge_left=elem_edges(2,elem)
              edge_right=elem_edges(1,elem)
           end if
           if (myList_edge2D(edge_left)>edge2D_in) then
              bEdge_left=bEdge_left+1
              belem_left(bEdge_left)=elem
           end if
           if (myList_edge2D(edge_right)>edge2D_in) then
              bEdge_right=bEdge_right+1
              belem_right(bEdge_right)=elem
           end if
        end do

    ! now we have three cases
       if (bEdge_left==0) then      ! inner contour
          elem=nod_in_elem2D(1, n)  ! we can start from any
          allocate(nedges(nod_in_elem2D_num(n)))
          nedges=0
          allocate(nelems(nod_in_elem2D_num(n)))
          nelems=0
          !!!!!!! inner_node_contour
include "node_contour_inner.h"
          if (pos_increment<0) then 
             current_pos=2*nod_in_elem2D_num(n)
          else
             current_pos =1
          end if
          do nn=1, nod_in_elem2D_num(n)
             call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
             current_pos=current_pos+pos_increment
             call elem_center(nelems(nn), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
             current_pos=current_pos+pos_increment
          end do
          current_pos=2*nod_in_elem2D_num(n)+1
          do nn=current_pos, size(my_x_corners, 2)
             my_x_corners(n, nn)=my_x_corners(n, current_pos-1)
             my_y_corners(n, nn)=my_y_corners(n, current_pos-1)
          end do
          deallocate(nedges, nelems)
       end if


       if (bEdge_left==1) then ! standard boundary node
          elem=belem_left(1)
          allocate(nedges(nod_in_elem2D_num(n)+1))
          nedges=0
          allocate(nelems(nod_in_elem2D_num(n)))
          nelems=0
          !!!!!!!boundary_node_contour
include "node_contour_boundary.h"
          if (pos_increment<0) then 
             current_pos=2*nod_in_elem2D_num(n)+2 !one more for the node n itself also we have 2 boundary edges
          else
             current_pos =1
          end if
          do nn=1, nod_in_elem2D_num(n)
             call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
             current_pos=current_pos+pos_increment
             call elem_center(nelems(nn), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
             current_pos=current_pos+pos_increment
          end do
          nn=nod_in_elem2D_num(n)+1
          call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
          current_pos=current_pos+pos_increment
          my_x_corners(n, current_pos)=coord_nod2D(1,n)
          my_y_corners(n, current_pos)=coord_nod2D(2,n)
          current_pos=2*nod_in_elem2D_num(n)+3
          do nn=current_pos, size(my_x_corners, 2)
             my_x_corners(n, nn)=my_x_corners(n, current_pos-1)
             my_y_corners(n, nn)=my_y_corners(n, current_pos-1)
          end do
          !!!!!!!
          deallocate(nedges, nelems)
       end if

       if (bEdge_left==2) then  ! strange boundary node
           elem=belem_left(1)
           allocate(nedges (nod_in_elem2D_num(n)+1))
           allocate(nedges1(nod_in_elem2D_num(n)+1))
           nedges =0
           nedges1=0
           allocate(nelems (nod_in_elem2D_num(n)))
           allocate(nelems1(nod_in_elem2D_num(n)))
           nelems=0
           nelems1=0
           if (pos_increment<0) then 
            current_pos=2*nod_in_elem2D_num(n)+4 !two more for the node n itself also we have 4 boundary edges
           else
            current_pos =1
           end if
           !!!!!!!boundary_node_contour
include "node_contour_boundary.h"
           where (nedges>0)
                 nedges1=nedges
           end where
           where (nelems>0)
                 nelems1=nelems
           end where
           nn1=nn
           do nn=1, nn1
              call edge_center(edges(1, nedges1(nn)), edges(2, nedges1(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
              current_pos=current_pos+pos_increment
              call elem_center(nelems1(nn), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
              current_pos=current_pos+pos_increment
           end do
           nn=nn1+1
           call edge_center(edges(1, nedges1(nn)), edges(2, nedges1(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
           current_pos=current_pos+pos_increment
           nn=nn1+2
           my_x_corners(n, current_pos)=coord_nod2D(1,n)
           my_y_corners(n, current_pos)=coord_nod2D(2,n)
           current_pos=current_pos+pos_increment
           !!!!!!!
           elem=belem_left(2)
           allocate(nedges2(nod_in_elem2D_num(n)+1))
           nedges =0
           nedges2=0
           allocate(nelems2(nod_in_elem2D_num(n)))
           nelems =0
           nelems2=0
           !!!!!!!boundary_node_contour
include "node_contour_boundary.h"
           where (nedges>0)
                nedges2=nedges
           end where
           where (nelems>0)
                 nelems2=nelems
           end where
           nn2=nn
           do nn=nn1+3, nn1+nn2+2
              call edge_center(edges(1, nedges2(nn)), edges(2, nedges2(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
              current_pos=current_pos+pos_increment
              call elem_center(nelems2(nn), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
              current_pos=current_pos+pos_increment
           end do
           nn=nn1+nn2+3
           call edge_center(edges(1, nedges2(nn)), edges(2, nedges2(nn)), my_x_corners(n, current_pos), my_y_corners(n, current_pos), mesh)
           current_pos=current_pos+pos_increment
           nn=nn1+nn2+4
           my_x_corners(n, nn)=coord_nod2D(1,n)
           my_y_corners(n, nn)=coord_nod2D(2,n)
           current_pos=2*nod_in_elem2D_num(n)+5
           do nn=current_pos, size(my_x_corners, 2)
              my_x_corners(n, nn)=my_x_corners(n, current_pos-1)
              my_y_corners(n, nn)=my_y_corners(n, current_pos-1)
           end do
           !!!!!!!
           deallocate(nedges, nelems, nedges1, nelems1, nedges2, nelems2)
       end if
    end do
    do n=1, myDim_nod2D
       do nn=1, size(my_x_corners, 2)
          this_x_coord=my_x_corners(n, nn)
          this_y_coord=my_y_corners(n, nn)
          call r2g(my_x_corners(n, nn), my_y_corners(n, nn), this_x_coord, this_y_coord)
       end do
    end do
    my_x_corners=my_x_corners/rad
    my_y_corners=my_y_corners/rad
    end subroutine node_contours

! kh 02.12.21
#if defined(__usetp)
  subroutine cpl_oasis3mct_init(partit, localCommunicator, num_fesom_groups)
#else
  subroutine cpl_oasis3mct_init(partit, localCommunicator)
#endif

    USE MOD_PARTIT
    implicit none   
    save

    !-------------------------------------------------------------------
    ! Initialize coupled mode communication for ocean
    ! exchange between AGCM, OGCM and COUPLER. (OASIS3-MCT software)
    !
    !--------------------------------------------------------------------
    ! Arguments
    !
    integer, intent(OUT)       :: localCommunicator
    type(t_partit), intent(inout), target :: partit
#if defined(__usetp)
! kh 02.12.21
    integer, intent(inout)     :: num_fesom_groups
#endif
    !
    ! Local declarations
    !
    !--------------------------------------------------------------------
    !

#ifdef VERBOSE
      print *, '================================================='
      print *, 'cpl_oasis3mct_init : coupler initialization for OASIS3-MCT'
      print *, '*************************************************'
      call fesom_flush
#endif /* VERBOSE */

    appl_name = 'ocean'
    comp_name = 'fesom'
    
    grid_name = 'feom'

    !------------------------------------------------------------------
    ! 1st Initialize the OASIS3-MCT coupling system for the application
    !------------------------------------------------------------------
! kh 02.12.21
#if defined(__usetp)
    CALL oasis_init_comp(comp_id, comp_name, ierror, num_program_groups = num_fesom_groups)
#else
    CALL oasis_init_comp(comp_id, comp_name, ierror )
#endif
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'Init_comp failed.')
    ENDIF

    ! Unit for output messages : one file for each process
    CALL MPI_Comm_Rank ( MPI_COMM_WORLD, commRank, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_rank failed.')
    ENDIF

! kh 02.12.21
#if defined(__usetp)
    CALL oasis_get_localcomm_all_groups( localCommunicator, ierror )
#else
    CALL oasis_get_localcomm( localCommunicator, ierror )
#endif
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'get_local_comm failed.')
    ENDIF

    ! Get MPI size and rank
    CALL MPI_Comm_Size ( localCommunicator, partit%npes, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_size failed.')
    ENDIF
    
    CALL MPI_Comm_Rank ( localCommunicator, partit%mype, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_rank failed.')
    ENDIF

  end subroutine cpl_oasis3mct_init

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_define_unstr(partit, mesh)
   
#ifdef __oifs
    use mod_oasis_auxiliary_routines, ONLY:	oasis_get_debug, oasis_set_debug
#else
    use mod_oasis_method, ONLY:	oasis_get_debug, oasis_set_debug
#endif
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_rotate_grid
    use mod_oasis, only: oasis_write_area, oasis_write_mask
    implicit none
    save
    type(t_mesh),   intent(in),    target :: mesh
    type(t_partit), intent(inout), target :: partit
    !-------------------------------------------------------------------
    ! Definition of grid and field information for ocean
    ! exchange between FESOM, ECHAM6 and OASIS3-MCT.
    !--------------------------------------------------------------------
    ! Arguments
    !
    ! Local declarations

    integer, parameter         :: nbr_grids   = 1
    integer, parameter         :: nbr_parts   = 1
    integer, parameter         :: nbr_masks   = 1
    integer, parameter         :: nbr_methods = 1

    integer                    :: grid_id(nbr_grids)    ! ids returned by oasis_def_grid
    integer                    :: part_id(nbr_parts)    ! ids returned by oasis_def_grid

    integer                    :: point_id(nbr_methods) ! ids returned by oasis_set_points

    integer                    :: mask_id(nbr_masks)    ! ids returned by oasis_set_mask

    integer                    :: shape(2,3)     ! shape of arrays passed to OASIS4
    integer                    :: total_shape(4) ! shape of arrays passed to OASIS3MCT
    integer                    :: nodim(2)
    integer                    :: data_type      ! data type of transients

    integer                    :: ig_paral(3)

    integer                    :: il_flag
    logical                    :: new_points

    integer                    :: i, j, k        ! local loop indicees
    integer                    :: l,m            ! local loop indicees

    character(len=32)          :: point_name     ! name of the grid points

    integer                    :: my_number_of_points
    integer                    :: number_of_all_points
    integer                    :: my_displacement

    integer,allocatable        :: counts_from_all_pes(:)
    integer,allocatable        :: displs_from_all_pes(:)
    integer,allocatable        :: unstr_mask(:,:)
    real(kind=WP)              :: this_x_coord          ! longitude coordinates
    real(kind=WP)              :: this_y_coord          ! latitude coordinates
    !
    ! Corner data structure for a OASIS3-MCT Reglonlatvrt grid
    !
    real(kind=WP), allocatable :: my_x_coords(:)     ! longitude coordinates
    real(kind=WP), allocatable :: my_y_coords(:)     ! latitude  coordinates

    real(kind=WP), allocatable :: all_x_coords(:, :)     ! longitude coordinates
    real(kind=WP), allocatable :: all_y_coords(:, :)     ! latitude  coordinates
    real(kind=WP), allocatable :: all_area(:,:)    

    real(kind=WP), allocatable :: my_x_corners(:,:)     ! local longitude node corners
    real(kind=WP), allocatable :: my_y_corners(:,:)     ! local latitude node corners
    real(kind=WP), allocatable :: all_x_corners(:,:,:)     ! global longitude node corners
    real(kind=WP), allocatable :: all_y_corners(:,:,:)     ! global latitude node corners

include "associate_part_def.h"
include "associate_mesh_def.h"
include "associate_part_ass.h"
include "associate_mesh_ass.h"


#ifdef VERBOSE
      print *, '=============================================================='
      print *, 'cpl_oasis3mct_define_unstr : coupler definition for OASIS3-MCT'
      print *, '**************************************************************'
      call fesom_flush
#endif /* VERBOSE */

    ! -----------------------------------------------------------------
    ! ... Some initialisation
    ! -----------------------------------------------------------------

!---wiso-code
    ALLOCATE(cpl_send(nsend))
    ALLOCATE(cpl_recv(nrecv))

    ALLOCATE(send_id(nsend))
    ALLOCATE(recv_id(nrecv))
!---wiso-code-end
    
    ALLOCATE(displs_from_all_pes(partit%npes))
    ALLOCATE(counts_from_all_pes(partit%npes))

    send_id = 0
    recv_id = 0

    ! -----------------------------------------------------------------
    ! ... Some MPI stuff relevant for optional exchange via root only
    ! -----------------------------------------------------------------

    localRank = mype 
    localSize = npes

#ifdef VERBOSE
    print *, 'This is FESOM local rank ', localRank, commRank
#endif /* VERBOSE */

! Note: We assume here that each process sends its local field.
! -----------------------------------------------------------------
! ... Define the partition
! -----------------------------------------------------------------

    my_number_of_points = myDim_nod2d
    number_of_all_points = nod2d
    if (mype .eq. 0) then 
      print *, 'FESOM Before ALLGATHERV'
    endif
    CALL MPI_ALLGATHER(my_number_of_points, 1, MPI_INTEGER, & 
                       counts_from_all_pes, 1, MPI_INTEGER, MPI_COMM_FESOM, ierror)
    if (mype .eq. 0) then
      print *, 'FESOM after ALLGATHERV'
    endif

    if (mype .eq. 0) then
      my_displacement = 0
    else
      my_displacement = SUM(counts_from_all_pes(1:mype))
    endif

    ig_paral(1) = 1                       ! Apple Partition
    ig_paral(2) = my_displacement         ! Global Offset
    ig_paral(3) = my_number_of_points     ! Local Extent

    if (mype .eq. 0) then
      print *, 'FESOM before def partition'
    endif
    CALL oasis_def_partition( part_id(1), ig_paral, ierror )
    if (mype .eq. 0) then
      print *, 'FESOM after def partition'
    endif
    if ( ierror /= 0 ) then
       print *, 'FESOM commRank def_partition failed'
       call oasis_abort(comp_id, 'cpl_oasis3mct_define_unstr', 'def_partition failed')
    endif
      
    ALLOCATE(my_x_coords(my_number_of_points))
    ALLOCATE(my_y_coords(my_number_of_points))
    ALLOCATE(my_x_corners(myDim_nod2D, 25))
    ALLOCATE(my_y_corners(myDim_nod2D, 25))

    do i = 1, my_number_of_points
      this_x_coord = coord_nod2D(1, i)
      this_y_coord = coord_nod2D(2, i)
      call r2g(my_x_coords(i), my_y_coords(i), this_x_coord, this_y_coord)
    end do   

    my_x_coords=my_x_coords/rad
    my_y_coords=my_y_coords/rad

    if (mype .eq. 0) then
      print *, 'FESOM before corner computation'
    endif
    call node_contours(my_x_corners, my_y_corners, partit, mesh)
    if (mype .eq. 0) then
      print *, 'FESOM after corner computation'
    endif

    if (mype .eq. localroot) then
      ALLOCATE(all_x_coords(number_of_all_points, 1))
      ALLOCATE(all_y_coords(number_of_all_points, 1))
      ALLOCATE(all_area(number_of_all_points, 1))
      ALLOCATE(all_x_corners(number_of_all_points, 1, 25))
      ALLOCATE(all_y_corners(number_of_all_points, 1, 25))
    else 
      ALLOCATE(all_x_coords(1, 1))
      ALLOCATE(all_y_coords(1, 1))
      ALLOCATE(all_area(1, 1))
      ALLOCATE(all_x_corners(1, 1, 1))
      ALLOCATE(all_y_corners(1, 1, 1))
    endif
   


    displs_from_all_pes(1) = 0
    do i = 2, npes
      displs_from_all_pes(i) = SUM(counts_from_all_pes(1:(i-1)))
    enddo  

    if (mype .eq. 0) then 
      print *, 'FESOM before 1st GatherV', displs_from_all_pes(npes), counts_from_all_pes(npes), number_of_all_points
    endif
    CALL MPI_GATHERV(my_x_coords, my_number_of_points, MPI_DOUBLE_PRECISION, all_x_coords,  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)

    if (mype .eq. 0) then 
      print *, 'FESOM before 2nd GatherV'
    endif
    CALL MPI_GATHERV(my_y_coords, my_number_of_points, MPI_DOUBLE_PRECISION, all_y_coords,  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)

    if (mype .eq. 0) then 
      print *, 'FESOM before 3rd GatherV'
    endif
    CALL MPI_GATHERV(area(1,:), my_number_of_points, MPI_DOUBLE_PRECISION, all_area,  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)

    do j = 1, 25
      CALL MPI_GATHERV(my_x_corners(:,j), myDim_nod2D, MPI_DOUBLE_PRECISION, all_x_corners(:,:,j),  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)
      CALL MPI_GATHERV(my_y_corners(:,j), myDim_nod2D, MPI_DOUBLE_PRECISION, all_y_corners(:,:,j),  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)
    end do

    CALL MPI_Barrier(MPI_COMM_FESOM, ierror)
    if (mype .eq. 0) then 
      print *, 'FESOM after Barrier'
    endif

! kh 30.11.21
#if defined(__usetp)
    if(my_fesom_group == 0) then
#endif

    if (mype .eq. localroot) then
      print *, 'FESOM before grid writing to oasis grid files'
       CALL oasis_start_grids_writing(il_flag)
       IF (il_flag .NE. 0) THEN

          print *, 'FESOM before write grid centers'
          CALL oasis_write_grid (grid_name, number_of_all_points, 1, all_x_coords(:,:), all_y_coords(:,:))

          print *, 'FESOM before write corner'
          CALL oasis_write_corner (grid_name, number_of_all_points, 1, 25, all_x_corners(:,:,:), all_y_corners(:,:,:))

          ALLOCATE(unstr_mask(number_of_all_points, 1))
          unstr_mask=0
          print *, 'FESOM before write mask'
          CALL oasis_write_mask(grid_name, number_of_all_points, 1, unstr_mask)
          DEALLOCATE(unstr_mask)

          print *, 'FESOM before write area'
          CALL oasis_write_area(grid_name, number_of_all_points, 1, all_area)

       end if
      print *, 'FESOM before terminate_grids_writing'
      call oasis_terminate_grids_writing()
      print *, 'FESOM after terminate_grids_writing'
    endif !localroot
     
#if defined(__usetp)
    end if !(my_fesom_group == 0) then     
#endif


    DEALLOCATE(all_x_coords, all_y_coords, my_x_coords, my_y_coords, displs_from_all_pes, counts_from_all_pes)
!------------------------------------------------------------------
! 3rd Declare the transient variables
!------------------------------------------------------------------
!
!
! ... Define symbolic names for the transient fields send by the ocean
!     These must be identical to the names specified in the SMIOC file.
!
#if defined (__oifs)
    cpl_send( 1)='sst_feom' ! 1. sea surface temperature [K]       ->
    cpl_send( 2)='sie_feom' ! 2. sea ice extent [%-100]            ->
    cpl_send( 3)='snt_feom' ! 3. snow thickness [m]                ->
    cpl_send( 4)='ist_feom' ! 4. sea ice surface temperature [K]   ->
    cpl_send( 5)='sia_feom' ! 5. sea ice albedo [%-100]            ->
    cpl_send( 6)='u_feom'   ! 6. eastward  surface velocity [m/s]  ->
    cpl_send( 7)='v_feom'   ! 7. northward surface velocity [m/s]  ->
#else
    cpl_send( 1)='sst_feom' ! 1. sea surface temperature [°C]      ->
    cpl_send( 2)='sit_feom' ! 2. sea ice thickness [m]             ->
    cpl_send( 3)='sie_feom' ! 3. sea ice extent [%-100]            ->
    cpl_send( 4)='snt_feom' ! 4. snow thickness [m]                ->
!---wiso-code
! add isotope coupling fields
    IF (lwiso) THEN
      cpl_send( 5)='o18w_oce' !                 -> h2o18 of ocean water
      cpl_send( 6)='hdow_oce' !                 -> hdo16 of ocean water
      cpl_send( 7)='o16w_oce' !                 -> h2o16 of ocean water
      cpl_send( 8)='o18i_oce' !                 -> h2o18 of sea ice
      cpl_send( 9)='hdoi_oce' !                 -> hdo16 of sea ice
      cpl_send(10)='o16i_oce' !                 -> h2o16 of sea ice
    END IF
!---wiso-code-end
#endif


    
!
! ...  Define symbolic names for transient fields received by the ocean.
!      These must be identical to the names specified in the SMIOC file.
!
#if defined (__oifs)
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
    cpl_recv(13) = 'enth_oce'
#else
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
! --- icebergs ---
    IF (lwiso) THEN
      cpl_recv(13) = 'w1_oce'
      cpl_recv(14) = 'w2_oce'
      cpl_recv(15) = 'w3_oce'
      cpl_recv(16) = 'i1_oce'
      cpl_recv(17) = 'i2_oce'
      cpl_recv(18) = 'i3_oce'
      IF (use_icebergs) THEN
        cpl_recv(19) = 'u10w_oce'
        cpl_recv(20) = 'v10w_oce'
      END IF
    ELSE IF (use_icebergs) THEN
      cpl_recv(13) = 'u10w_oce'
      cpl_recv(14) = 'v10w_oce'
    END IF
! --- icebergs ---
#endif

    if (mype .eq. 0) then 
       print *, 'FESOM after declaring the transient variables'
    endif

    data_type = oasis_DOUBLE
!     nodim(1): the overall number of dimensions of array that is going
!               to be transferred with var_id, i.e. same than grid_nodims
!               except for bundle variables for which it is one more.

    nodim(1) = 1 ! check, 3 for OASIS4
    nodim(2) = 1
!
! ... Announce send variables, all on T points. 
!
    total_shape(1) = number_of_all_points
    total_shape(2) = 1

    do i = 1, nsend
       call oasis_def_var(send_id(i), cpl_send(i), part_id(1), &
            nodim, oasis_Out, total_shape, oasis_Double, ierror)
       if (ierror /= oasis_Ok) then
          print *, 'Failed to define transient ', i, trim(cpl_send(i))
          call oasis_abort(comp_id, 'cpl_oasis3mct_define_unstr', 'def_var failed')
       endif
    enddo

    if (mype .eq. 0) then 
       print *, 'FESOM after announcing send variables'
    endif

!
!
! ... Announce recv variables, all on T points
!
    do i = 1, nrecv
       call oasis_def_var(recv_id(i), cpl_recv(i), part_id(1), &
            nodim, oasis_In, total_shape, oasis_Double, ierror)
       if (ierror /= oasis_Ok) then
          print *, 'Failed to define transient ', i, trim(cpl_recv(i))
          call oasis_abort(comp_id, 'cpl_oasis3mct_define_unstr', 'def_var failed')
       endif
    enddo

    if (mype .eq. 0) then 
       print *, 'FESOM after announcing revieved variables'
    endif

!------------------------------------------------------------------
! 4th End of definition phase
!------------------------------------------------------------------

   call oasis_enddef(ierror)
   if (commRank) print *, 'fesom oasis_enddef: COMPLETED'
#ifndef __oifs
   if (commRank) print *, 'FESOM: calling exchange_roots'
   call exchange_roots(source_root, target_root, 1, partit%MPI_COMM_FESOM, MPI_COMM_WORLD)
   if (commRank) print *, 'FESOM source/target roots: ', source_root, target_root
#endif

   ! WAS VOM FOLGENDEN BRAUCHE ICH NOCH ??? 

   allocate(cplsnd(nsend, myDim_nod2D+eDim_nod2D))
   allocate(exfld(myDim_nod2D))
   cplsnd=0.
   o2a_call_count=0

   CALL MPI_BARRIER(MPI_COMM_FESOM, ierror)
   if (mype .eq. 0) then 
      print *, 'Survived cpl_oasis3mct_define'
   endif   

  end subroutine cpl_oasis3mct_define_unstr

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_send(ind, data_array, action, partit)
    use o_param
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none
    save
    !---------------------------------------------------------------------
    !
    ! At each coupling time-step this routine 
    ! checks whether fields have to be send.
    ! If fields have to be send it
    !
    !  - collects to data from the FESOM unstructured grid
    !  - interpolates onto a regular grid
    !  - and sends fields
    !
    !----------------------------------------------------------------------
    ! Arguments
    !
    integer, intent( IN )          :: ind       ! variable Id
    logical, intent( OUT )         :: action    !
    type(t_partit), intent(in)     :: partit
    real(kind=WP),  intent(IN)     :: data_array(partit%myDim_nod2D+partit%eDim_nod2D)
    !
    ! Local declarations
    !
    integer                :: info
    !
    real (kind=WP)          :: t1, t2, t3
    !
    !--------------------------------------------------------------------
    !
    t1=MPI_Wtime()
    action = .false.

    if (ind==1) o2a_call_count=o2a_call_count+1
    cplsnd(ind, :)=cplsnd(ind, :)+data_array
    ! call do_oce_2_atm(cplsnd(ind, :)/real(o2a_call_count), atm_fld, 1)

    exfld = cplsnd(ind, 1:partit%myDim_nod2D)/real(o2a_call_count)

    t2=MPI_Wtime()
#ifdef VERBOSE
    if (partit%mype==0) then
        print *, 'FESOM oasis_send: ', cpl_send(ind)   
    endif     
#endif
    call oasis_put(send_id(ind), seconds_til_now, exfld, info)
    action=(info==4 .OR. info==8)
    if (action) then
       if (ind==nsend) then
          cplsnd=0.
          o2a_call_count=0
       end if
    end if
    t3=MPI_Wtime()
    
  end subroutine cpl_oasis3mct_send

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_recv(ind, data_array, action, partit)
    use o_param
    use g_comm_auto
    USE MOD_PARTIT
    USE MOD_PARSUP
    implicit none
    save
    !---------------------------------------------------------------------
    !
    ! At each coupling time-step,this routine receives fields
    !
    !----------------------------------------------------------------------
    ! Arguments
    !
    integer, intent( IN )  :: ind       ! variable Id
    logical, intent( OUT ) :: action    ! 
    real(kind=WP), intent( OUT )    :: data_array(:)
    type(t_partit), intent(inout), target :: partit
    !
    ! Local declarations
    !
    integer                :: info
    integer                :: j
    integer, save          :: ncount = 0
    real (kind=WP)         :: t1, t2, t3        
    !
    !--------------------------------------------------------------------
    !    
    t1=MPI_Wtime()        
    action=.false.
    !
    ! receive data from OASIS3-MCT on local root
    !
#ifdef VERBOSE
    if (partit%mype==0) then
        print *, 'oasis_get: ', cpl_recv(ind)
    endif    
#endif

#if defined(__usetp)
! kh 06.12.21 the coupling is in principle as it was before, i.e. the fesom processes - in group 0 - receive their data from echam
    if(my_fesom_group == 0) then
#endif

    call oasis_get(recv_id(ind), seconds_til_now, exfld,info)

#if defined(__usetp)
    else

! kh 06.12.21 defensive: assignment statement "action=(info==3 ..." below is "don't care" in this case, because the actual value for action
! is received via MPI_Bcast anyway
        info = 0

    end if
#endif

    t2=MPI_Wtime()
 !
 ! FESOM's interpolation routine interpolates structured
 ! VarStrLoc coming from OASIS3MCT to local unstructured data_array
 ! and delivered back to FESOM.
   action=(info==3 .OR. info==10 .OR. info==11 .OR. info==12 .OR. info==13)

#if defined(__usetp)
! kh 03.12.21
   if(num_fesom_groups > 1) then
      call MPI_Bcast(action, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, MPIerr)
   end if
#endif 

   if (action) then
#if defined(__usetp)
! kh 03.12.21
      if(my_fesom_group == 0) then
#endif
      data_array(1:partit%myDim_nod2d) = exfld
#if defined(__usetp)
      end if

! kh 03.12.21
      if(num_fesom_groups > 1) then
          call MPI_Bcast(data_array, myDim_nod2d, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, MPIerr)
      end if
#endif

      call exchange_nod(data_array, partit)
   end if   
   t3=MPI_Wtime()
  end subroutine cpl_oasis3mct_recv
  
!
!  Exchange the roots between components
!  



SUBROUTINE exchange_roots(source_root, target_root, il_side, &
                                local_comm, global_comm) 
!source root (returned) is the 0 process of the source component (i.e. comm_echam) in comm_psmile
!target root (returned) is the 0 process of the target component (i.e. comm_fesom) in comm_psmile
!is side(input) is 1 on the source side (comm_echam) and 0 on the target side (comm_fesom)
!local_comm  (i.e. comm_echam here)
!global_comm (i.e. comm_psmile here)

        IMPLICIT NONE

        INTEGER, INTENT(IN)        :: il_side
        INTEGER, INTENT(IN)        :: local_comm, global_comm
        INTEGER, INTENT(OUT)       :: source_root, target_root

        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: local_rank, my_global_rank, ierror


        source_root = 500000
        target_root = 500000

        CALL MPI_COMM_RANK(local_comm, local_rank, ierror)
        CALL MPI_COMM_RANK(global_comm, my_global_rank, ierror)
        write(*,*) 'local_rank=', local_rank, ' ; global_rank=', my_global_rank
        IF (local_rank == 0) THEN
          IF (il_side == 0) THEN
            source_root = my_global_rank
            IF (source_root .ne. 0) THEN
              CALL MPI_Send(source_root, 1, MPI_INTEGER, 0, 1, global_comm, ierror)
            END IF
          ELSE IF (il_side == 1) THEN
            target_root = my_global_rank
            IF (target_root .ne. 0) THEN
              CALL MPI_Send(target_root, 1, MPI_INTEGER, 0, 2, global_comm, ierror)
            END IF
          END IF
        END IF
        IF (my_global_rank == 0) THEN
          IF (my_global_rank .ne. source_root) THEN
            CALL MPI_Recv(source_root, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, global_comm, status, ierror)
          END IF
          IF (my_global_rank .ne. target_root) THEN
            CALL MPI_Recv(target_root, 1, MPI_INTEGER, MPI_ANY_SOURCE, 2, global_comm, status, ierror)
          END IF
        END IF
        CALL MPI_BCast(source_root, 1, MPI_INTEGER, 0, global_comm, ierror)
        CALL MPI_BCast(target_root, 1, MPI_INTEGER, 0, global_comm, ierror)
END SUBROUTINE exchange_roots

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_finalize

    implicit none
    save

    !---------------------------------------------------------------------
    !
    ! Finalizes the coupling. If MPI_init has not been called
    ! explicitly before cpl_oasis3mct_finalize it will also close
    ! MPI communication.
    !
    !----------------------------------------------------------------------

 
!!!!#ifdef VERBOSE
      print *, '================================================='
      print *, 'cpl_oasis3mct_finalize : ending coupling process via oasis3mct'
      print *, '*************************************************'
      call fesom_flush
!!!!#endif /* VERBOSE */

      call oasis_terminate( ierror )

  end subroutine cpl_oasis3mct_finalize

  subroutine fesom_flush ()
  end subroutine fesom_flush

end module cpl_driver
#else
module cpl_driver
end module cpl_driver
#endif
