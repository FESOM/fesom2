#if defined (__oasis)
module cpl_driver
  !======================================================================
  !
  ! for coupling between the FESOM and an AOGCM using OASIS3-MCT
  !
  !=====================================================================
  ! History :
  !  09-09  (R. Redler, Germany)  Original code
  !  09-09  (K.Fieg, AWI Germany) Adjustment for FESOM
  !  07-12  (D.Barbi, AWI Germany) Switch to ECHAM6.1 and OASIS3-MCT
  !  01-19  (J.Streffing, AWI Germany) Added OpenIFS coupling
  !  03-23  (J.Streffing, AWI Germany) Added corner point computation
  !                                    for 1st order conserv remapping  
  !----------------------------------------------------------------------
  ! Modules used
  !
  use mod_oasis                    ! oasis module
  use g_config, only : dt, use_icebergs, lwiso
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

  subroutine cpl_oasis3mct_init(partit, localCommunicator )
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
    CALL oasis_init_comp(comp_id, comp_name, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'Init_comp failed.')
    ENDIF

    ! Unit for output messages : one file for each process
    CALL MPI_Comm_Rank ( MPI_COMM_WORLD, commRank, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_rank failed.')
    ENDIF

    CALL oasis_get_localcomm( localCommunicator, ierror )
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

    integer                    :: i, j, k            ! local loop indicees
    integer                    :: l,m,n, done        ! local loop indicees

    character(len=32)          :: point_name     ! name of the grid points

    integer                    :: my_number_of_points
    integer                    :: number_of_all_points
    integer                    :: counts_from_all_pes(partit%npes)
    integer                    :: displs_from_all_pes(partit%npes)
    integer                    :: my_displacement
    integer                    :: my_max_elem(partit%npes)
    integer                    :: my_max_edge(partit%npes)
    integer                    :: all_max_elem, all_max_edge, n_neg, n_pos
    integer                    :: el(2), enodes(2), edge

    integer,allocatable        :: unstr_mask(:,:), coastal_edge_list(:,:)
    real(kind=WP)              :: max_x          ! max longitude on corners of control volume
    real(kind=WP)              :: min_x          ! min longitude on corners of control volume
    real(kind=WP)              :: temp                  ! temp storage for corner sorting
    real(kind=WP)              :: this_x_coord          ! longitude coordinates
    real(kind=WP)              :: this_y_coord          ! latitude coordinates
    real(kind=WP)              :: this_x_corners        ! longitude node corners
    real(kind=WP)              :: this_y_corners        ! latitude node corners
    !
    ! Corner data structure for a OASIS3-MCT Reglonlatvrt grid
    !
    real(kind=WP), allocatable :: pos_x(:)       ! longitude to the right of dateline
    real(kind=WP), allocatable :: pos_y(:)       ! latitude to the right of dateline
    real(kind=WP), allocatable :: neg_x(:)       ! longitude to the left of dateline
    real(kind=WP), allocatable :: neg_y(:)       ! latitude to the left of dateline
    real(kind=WP), allocatable :: temp_x_coord(:)    ! longitude coordinates
    real(kind=WP), allocatable :: temp_y_coord(:)    ! longitude coordinates
    real(kind=WP), allocatable :: my_x_coords(:)     ! longitude coordinates
    real(kind=WP), allocatable :: my_y_coords(:)     ! latitude  coordinates
    real(kind=WP), allocatable :: angle(:,:)         ! array for holding corner angle for sorting
    real(kind=WP), allocatable :: my_x_corners(:,:)     ! longitude node corners
    real(kind=WP), allocatable :: my_y_corners(:,:)     ! latitude node corners
    real(kind=WP), allocatable :: coord_e_edge_center(:,:,:)   ! edge center coords
    real(kind=WP), allocatable :: all_x_coords(:, :)     ! longitude coordinates
    real(kind=WP), allocatable :: all_y_coords(:, :)     ! latitude  coordinates
    real(kind=WP), allocatable :: all_x_corners(:,:,:)    ! longitude node corners
    real(kind=WP), allocatable :: all_y_corners(:,:,:)    ! latitude node corners
    real(kind=WP), allocatable :: all_area(:,:)    
    logical, allocatable       :: coastal_nodes(:)    


#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


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

    CALL MPI_BARRIER(MPI_COMM_FESOM, ierror)

    my_max_elem=0
    my_max_elem = maxval(nod_in_elem2D_num(1:myDim_nod2D))
    all_max_elem = 0
    call MPI_Allreduce(my_max_elem, all_max_elem, &
          1, MPI_INTEGER, MPI_MAX, &
          MPI_COMM_FESOM, MPIerr)

    my_max_edge=0
    my_max_edge=maxval(nn_num)
    all_max_edge=0
    call MPI_AllREDUCE( my_max_edge, all_max_edge, &
          1, MPI_INTEGER,MPI_MAX, &
          MPI_COMM_FESOM, MPIerr)

    CALL MPI_BARRIER(MPI_COMM_FESOM, ierror)

    if (mype .eq. 0) then
      print *, 'Max elements per node:', all_max_elem, 'Max edges per node:', all_max_edge
      print *, 'FESOM before def partition'
    endif

    ig_paral(1) = 1                       ! Apple Partition
    ig_paral(2) = my_displacement         ! Global Offset
    ig_paral(3) = my_number_of_points     ! Local Extent

    ! For MPI_GATHERV we need the location of the local segment in the global vector
    displs_from_all_pes(1) = 0
    do i = 2, npes
      displs_from_all_pes(i) = SUM(counts_from_all_pes(1:(i-1)))
    enddo  

    CALL oasis_def_partition( part_id(1), ig_paral, ierror )
    if (mype .eq. 0) then
      print *, 'FESOM after def partition'
    endif
    if ( ierror /= 0 ) then
       print *, 'FESOM commRank def_partition failed'
       call oasis_abort(comp_id, 'cpl_oasis3mct_define_unstr', 'def_partition failed')
    endif

    ALLOCATE(coastal_nodes(number_of_all_points))
    ALLOCATE(angle(my_number_of_points,all_max_elem+all_max_edge))
    ALLOCATE(my_x_corners(my_number_of_points,all_max_elem+all_max_edge))
    ALLOCATE(my_y_corners(my_number_of_points,all_max_elem+all_max_edge))
    ALLOCATE(coord_e_edge_center(2,my_number_of_points, all_max_edge))
   
    ! We need to know for every node if any of it's edges are coastal, because 
    ! in case they are the center point will be a corner of the nodal area
    coastal_nodes=.False.
    allocate (coastal_edge_list(my_number_of_points*2,my_number_of_points*2))
    do edge=1, myDim_edge2D
      ! local indice of nodes that span up edge
      enodes=edges(:,edge)      
      ! local index of element that contribute to edge
      el=edge_tri(:,edge)
      if(el(2)>0) then
        ! Inner edge
        continue
      else   
        ! Boundary/coastal edge
        coastal_nodes(enodes(1))=.True.
        coastal_nodes(enodes(2))=.True.
        coastal_edge_list(enodes(1),enodes(2))=edge
        coastal_edge_list(enodes(2),enodes(1))=edge
      end if  
    end do


    ! For every node, loop over neighbours, calculate edge center as mean of node center and neighbour node center.
    coord_e_edge_center=0
    do i = 1, my_number_of_points
      ! if we are on coastal node, include node center n=1 as corner
      if (coastal_nodes(i)==.True.) then 
        do n = 1, nn_num(i)
          call edge_center(i, nn_pos(n,i), this_x_coord, this_y_coord, mesh)
          call r2g(coord_e_edge_center(1,i,n), coord_e_edge_center(2,i,n), this_x_coord, this_y_coord)
        end do
      ! else we skip n=1 and use only the edge centers n=2:nn_num(i)
      else
        do n = 2, nn_num(i)
          call edge_center(i, nn_pos(n,i), this_x_coord, this_y_coord, mesh)
          call r2g(coord_e_edge_center(1,i,n-1), coord_e_edge_center(2,i,n-1), this_x_coord, this_y_coord)
        end do
      end if
    end do

    ALLOCATE(my_x_coords(my_number_of_points))
    ALLOCATE(my_y_coords(my_number_of_points))

    ! Obtain center coordinates as node center on open ocean and as mean of corners at coastline
    do i = 1, my_number_of_points
      ! Center coord as mean of corner coordiantes along coastline
      if (coastal_nodes(i)==.True.) then
        ! So we define temp_corner coordiantes 
        allocate(temp_x_coord(nod_in_elem2D_num(i)+nn_num(i)))
        allocate(temp_y_coord(nod_in_elem2D_num(i)+nn_num(i)))
        temp_x_coord=0
        temp_y_coord=0
        do j = 1, nod_in_elem2D_num(i)
          temp_x_coord(j) = x_corners(i,j)*rad
          temp_y_coord(j) = y_corners(i,j)*rad
        end do
        ! Loop over edges
        do j = 1, nn_num(i)
          ! We skip coastal edge center points for the new center point calculation
          ! such that 1 element islands have the node center at the right angle
          ! We only do so if n elements is > 2, to avoid having only 3 corners
          if ((j>1) .and. (nod_in_elem2D_num(i) > 2)) then
            edge = coastal_edge_list(i,nn_pos(j,i))
            ! if edge is coastal, we leave it out of the mean equation, replaced by the node center
            if (edge>0) then
              this_x_coord = coord_nod2D(1, i)
              this_y_coord = coord_nod2D(2, i)
              ! unrotate grid
              call r2g(my_x_coords(i), my_y_coords(i), this_x_coord, this_y_coord)
              temp_x_coord(j+nod_in_elem2D_num(i))=my_x_coords(i)
              temp_y_coord(j+nod_in_elem2D_num(i))=my_y_coords(i)
            ! case for only two elements, we need the real edge centers to ensure center coord
            ! is inside polygon
            else
              temp_x_coord(j+nod_in_elem2D_num(i)) = coord_e_edge_center(1,i,j)
              temp_y_coord(j+nod_in_elem2D_num(i)) = coord_e_edge_center(2,i,j)
            end if
          ! Open ocean case, we just use the corner coords
          else
            temp_x_coord(j+nod_in_elem2D_num(i)) = coord_e_edge_center(1,i,j)
            temp_y_coord(j+nod_in_elem2D_num(i)) = coord_e_edge_center(2,i,j)
          end if
        end do
        min_x = minval(temp_x_coord)
        max_x = maxval(temp_x_coord)
        ! if we are at dateline (fesom cell larger than pi)
        if (max_x-min_x > pi) then

          ! set up separate data structures for the two hemispheres
          n_pos=count(temp_x_coord>=0)
          n_neg=count(temp_x_coord<0)
          allocate(pos_x(n_pos))
          allocate(pos_y(n_pos))
          allocate(neg_x(n_neg))
          allocate(neg_y(n_neg))
          pos_x = 0
          pos_y = 0
          neg_x = 0
          neg_x = 0
          n=1
          do j = 1, size(temp_x_coord)
            ! build separate corner vectors for the hemispheres
            if (temp_x_coord(j) >= 0) then
              pos_x(n) = temp_x_coord(j)
              pos_y(n) = temp_y_coord(j)
              n=n+1
            end if
          end do
          n=1
          do j = 1, size(temp_x_coord)
            if (temp_x_coord(j) < 0) then
              neg_x(n) = temp_x_coord(j)
              neg_y(n) = temp_y_coord(j)
              n=n+1
            end if
          end do
          ! if sum on right side of dateline are further from the dateline we shift the negative sum over to the right
          if (-sum(pos_x)+pi*n_pos >= sum(neg_x)+pi*n_neg) then
            this_x_coord = (sum(pos_x) + sum(neg_x) + 2*pi*n_neg) / (n_pos + n_neg)
            this_y_coord = (sum(pos_y) + sum(neg_y)) / (n_pos + n_neg)
          ! else we shift the positive sum over to the left side
          else
            this_x_coord = (sum(pos_x) - 2*pi*n_pos + sum(neg_x)) / (n_pos + n_neg)
            this_y_coord = (sum(pos_y) + sum(neg_y)) / (n_pos + n_neg)
          end if
          deallocate(pos_x,pos_y,neg_x,neg_y)
        ! max_x-min_x > pi -> we are not at dateline, just a normal mean is enough
        else
          this_x_coord = sum(temp_x_coord)/(size(temp_x_coord))
          this_y_coord = sum(temp_y_coord)/(size(temp_y_coord))
        end if
        my_x_coords(i)=this_x_coord
        my_y_coords(i)=this_y_coord
        deallocate(temp_x_coord, temp_y_coord)
      ! coastal_nodes(i)==.True. -> Node center on open ocean, we can use node center
      else
        this_x_coord = coord_nod2D(1, i)
        this_y_coord = coord_nod2D(2, i)
        ! unrotate grid
        call r2g(my_x_coords(i), my_y_coords(i), this_x_coord, this_y_coord)
      end if
    end do

    ! Add the different corner types to single array in preparation for angle calculation
    do i = 1, my_number_of_points
      ! First for element center based corners
      do j = 1, nod_in_elem2D_num(i)
        my_x_corners(i,j) = x_corners(i,j)*rad ! atan2 takes radian and elem corners come in grad
        my_y_corners(i,j) = y_corners(i,j)*rad
      end do
      ! Then we repeat for edge center coordinate
      ! The the coast j=1 is the node center
      if (coastal_nodes(i)==.True.) then
        do j = 1, nn_num(i)
          my_x_corners(i,j+nod_in_elem2D_num(i)) = coord_e_edge_center(1,i,j)
          my_y_corners(i,j+nod_in_elem2D_num(i)) = coord_e_edge_center(2,i,j)
        end do
      ! On open ocean we dont use the node center as corner, and thus have one less corner
      else
        do j = 1, nn_num(i)-1
          my_x_corners(i,j+nod_in_elem2D_num(i)) = coord_e_edge_center(1,i,j)
          my_y_corners(i,j+nod_in_elem2D_num(i)) = coord_e_edge_center(2,i,j)
        end do
      end if
    end do

    ! calculate angle between corners and center
    do i = 1, my_number_of_points
      if (coastal_nodes(i)==.True.) then
        n=0
      else
        n=1
      end if
      do j = 1, nod_in_elem2D_num(i)+nn_num(i)-n
        ! If they have different sign we are near the dateline and need to bring the corner onto
        ! the same hemisphere as the center (only for angle calc, the coord for oasis remains as before)
        ! Default: same sign -> normal atan2
        if (my_x_coords(i) <=0 .and. my_x_corners(i,j) <=0 .or. my_x_coords(i) >0 .and. my_x_corners(i,j) >0) then 
          angle(i,j) = atan2(my_x_corners(i,j) - my_x_coords(i), my_y_corners(i,j) - my_y_coords(i))
        else
          ! at dateline center is on the right side
          if (my_x_coords(i) >=pi/2) then 
            angle(i,j) = atan2(my_x_corners(i,j) + 2*pi - my_x_coords(i), my_y_corners(i,j) - my_y_coords(i))
          ! at dateline center is on the left side
          else if (my_x_coords(i) <=-pi/2) then
            angle(i,j) = atan2(my_x_corners(i,j) - 2*pi - my_x_coords(i), my_y_corners(i,j) - my_y_coords(i))
          ! at prime meridan -> also default
          else
            angle(i,j) = atan2(my_x_corners(i,j) - my_x_coords(i), my_y_corners(i,j) - my_y_coords(i))
          end if
        end if
      end do
    end do   

    ! Oasis requires corners sorted counterclockwise, so we sort by angle
    do i = 1, my_number_of_points
      if (coastal_nodes(i)==.True.) then
        n=0
      else
        n=1
      end if
      do l = 1, nod_in_elem2D_num(i)+nn_num(i)-1-n
        do m = l+1, nod_in_elem2D_num(i)+nn_num(i)-n
          if (angle(i,l) < angle(i,m)) then
            ! Swap angle
            temp = angle(i,m)
            angle(i,m) = angle(i,l)
            angle(i,l) = temp
            ! Swap lon
            temp = my_x_corners(i,m)
            my_x_corners(i,m) = my_x_corners(i,l)
            my_x_corners(i,l) = temp
            ! Swap lat
            temp = my_y_corners(i,m)
            my_y_corners(i,m) = my_y_corners(i,l)
            my_y_corners(i,l) = temp
          end if
        end do
      end do
    end do

    ! We can have a variable number of corner points.
    ! Luckly oasis can deal with that by just repeating the last one.
    ! Note, we are only allowed to repeat one coordinate and 
    ! the last one is not an element center, but an edge center
    do i = 1, my_number_of_points
      do j = 1, all_max_elem+all_max_edge
        if (coastal_nodes(i)==.True.) then
          if (j < nod_in_elem2D_num(i)+nn_num(i)) then
            my_y_corners(i,j)=my_y_corners(i,j)
            my_x_corners(i,j)=my_x_corners(i,j)
          else
            my_y_corners(i,j)=my_y_corners(i,nod_in_elem2D_num(i)+nn_num(i))
            my_x_corners(i,j)=my_x_corners(i,nod_in_elem2D_num(i)+nn_num(i))
          end if
        else
          if (j < nod_in_elem2D_num(i)+nn_num(i)-1) then
            my_y_corners(i,j)=my_y_corners(i,j)
            my_x_corners(i,j)=my_x_corners(i,j)
          else
            my_y_corners(i,j)=my_y_corners(i,nod_in_elem2D_num(i)+nn_num(i)-1)
            my_x_corners(i,j)=my_x_corners(i,nod_in_elem2D_num(i)+nn_num(i)-1)
          end if
        end if
      end do
    end do

    ! Oasis takes grad angles
    my_x_coords=my_x_coords/rad
    my_y_coords=my_y_coords/rad
    my_x_corners=my_x_corners/rad
    my_y_corners=my_y_corners/rad

    
    if (mype .eq. localroot) then
      ALLOCATE(all_x_coords(number_of_all_points, 1))
      ALLOCATE(all_y_coords(number_of_all_points, 1))
      ALLOCATE(all_x_corners(number_of_all_points, 1, all_max_elem+all_max_edge))
      ALLOCATE(all_y_corners(number_of_all_points, 1, all_max_elem+all_max_edge))
      ALLOCATE(all_area(number_of_all_points, 1))
    else 
      ALLOCATE(all_x_coords(1, 1))
      ALLOCATE(all_y_coords(1, 1))
      ALLOCATE(all_x_corners(1, 1, all_max_elem+all_max_edge))
      ALLOCATE(all_y_corners(1, 1, all_max_elem+all_max_edge))
      ALLOCATE(all_area(1, 1))
    endif


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
      print *, 'FESOM before 3rd GatherV', displs_from_all_pes(npes), counts_from_all_pes(npes), number_of_all_points
    endif

    do j = 1, all_max_elem+all_max_edge
      CALL MPI_GATHERV(my_x_corners(:,j), my_number_of_points, MPI_DOUBLE_PRECISION, all_x_corners(:,:,j),  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)
      CALL MPI_GATHERV(my_y_corners(:,j), my_number_of_points, MPI_DOUBLE_PRECISION, all_y_corners(:,:,j),  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)
    end do

    if (mype .eq. 0) then 
      print *, 'FESOM before 4th GatherV'
    endif
    CALL MPI_GATHERV(area(1,:), my_number_of_points, MPI_DOUBLE_PRECISION, all_area,  &
                    counts_from_all_pes, displs_from_all_pes, MPI_DOUBLE_PRECISION, localroot, MPI_COMM_FESOM, ierror)

    if (mype .eq. 0) then 
      print *, 'FESOM after 4th GatherV'
    endif


    CALL MPI_Barrier(MPI_COMM_FESOM, ierror)
    if (mype .eq. 0) then 
      print *, 'FESOM after Barrier'
    endif

    if (mype .eq. localroot) then
      print *, 'FESOM before start_grids_writing'
       CALL oasis_start_grids_writing(il_flag)
       IF (il_flag .NE. 0) THEN

          print *, 'FESOM before write grid'
          CALL oasis_write_grid (grid_name, number_of_all_points, 1, all_x_coords(:,:), all_y_coords(:,:))

          print *, 'FESOM before write corner'
          CALL oasis_write_corner (grid_name, number_of_all_points, 1, all_max_elem+all_max_edge, all_x_corners(:,:,:), all_y_corners(:,:,:))

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

    DEALLOCATE(all_x_coords, all_y_coords, my_x_coords, my_y_coords) 
    DEALLOCATE(all_x_corners, all_y_corners, my_x_corners, my_y_corners, angle) 
    DEALLOCATE(coastal_nodes, coord_e_edge_center) 
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
   if (mype .eq. 0) then 
      print *, 'After enddef'
   endif   

   ! WAS VOM FOLGENDEN BRAUCHE ICH NOCH ??? 

   allocate(cplsnd(nsend, myDim_nod2D+eDim_nod2D))
   allocate(exfld(myDim_nod2D))
   cplsnd=0.
   o2a_call_count=0
   if (mype .eq. 0) then 
      print *, 'Before last barrier'
   endif   

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

    call oasis_get(recv_id(ind), seconds_til_now, exfld,info)
    t2=MPI_Wtime()
 !
 ! FESOM's interpolation routine interpolates structured
 ! VarStrLoc coming from OASIS3MCT to local unstructured data_array
 ! and delivered back to FESOM.
   action=(info==3 .OR. info==10 .OR. info==11 .OR. info==12 .OR. info==13)
   if (action) then
      data_array(1:partit%myDim_nod2d) = exfld
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
