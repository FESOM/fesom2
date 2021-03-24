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
  use g_config, only : dt
  use o_param,  only : rad
  use g_PARSUP
  implicit none
  save   
  !
  ! Exchange parameters for coupling FESOM with ECHAM6
  !

#if defined (__oifs)
  integer, parameter         :: nsend = 5
  integer, parameter         :: nrecv = 13
#else
  integer, parameter         :: nsend = 4
  integer, parameter         :: nrecv = 12
#endif
  
  integer, dimension(nsend)  :: send_id
  integer, dimension(nrecv)  :: recv_id

  character(len=32)          :: cpl_send(nsend)
  character(len=32)          :: cpl_recv(nrecv)

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

  real(kind=WP)              :: time_send(2), time_recv(2)

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

  subroutine cpl_oasis3mct_init( localCommunicator )
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
    CALL MPI_Comm_Size ( localCommunicator, npes, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_size failed.')
    ENDIF
    
    CALL MPI_Comm_Rank ( localCommunicator, mype, ierror )
    IF (ierror /= 0) THEN
        CALL oasis_abort(comp_id, 'cpl_oasis3mct_init', 'comm_rank failed.')
    ENDIF

  end subroutine cpl_oasis3mct_init

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_define_unstr(mesh)
   
#ifdef __oifs
    use mod_oasis_auxiliary_routines, ONLY:	oasis_get_debug, oasis_set_debug
#else
    use mod_oasis_method, ONLY:	oasis_get_debug, oasis_set_debug
#endif
    use mod_mesh
    use g_rotate_grid
    use mod_oasis, only: oasis_write_area, oasis_write_mask
    implicit none
    save
    type(t_mesh), intent(in), target :: mesh
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
    integer                    :: counts_from_all_pes(npes)
    integer                    :: displs_from_all_pes(npes)
    integer                    :: my_displacement

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

#include "associate_mesh.h"

#ifdef VERBOSE
      print *, '=============================================================='
      print *, 'cpl_oasis3mct_define_unstr : coupler definition for OASIS3-MCT'
      print *, '**************************************************************'
      call fesom_flush
#endif /* VERBOSE */

    ! -----------------------------------------------------------------
    ! ... Some initialisation
    ! -----------------------------------------------------------------

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

    do i = 1, my_number_of_points
      this_x_coord = coord_nod2D(1, i)
      this_y_coord = coord_nod2D(2, i)
      call r2g(my_x_coords(i), my_y_coords(i), this_x_coord, this_y_coord)
    end do   

    my_x_coords=my_x_coords/rad
    my_y_coords=my_y_coords/rad

    if (mype .eq. localroot) then
      ALLOCATE(all_x_coords(number_of_all_points, 1))
      ALLOCATE(all_y_coords(number_of_all_points, 1))
      ALLOCATE(all_area(number_of_all_points, 1))
    else 
      ALLOCATE(all_x_coords(1, 1))
      ALLOCATE(all_y_coords(1, 1))
      ALLOCATE(all_area(1, 1))
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

    if (mype .eq. 0) then 
      print *, 'FESOM after 3rd GatherV'
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
#else
    cpl_send( 1)='sst_feom' ! 1. sea surface temperature [°C]      ->
    cpl_send( 2)='sit_feom' ! 2. sea ice thickness [m]             ->
    cpl_send( 3)='sie_feom' ! 3. sea ice extent [%-100]            ->
    cpl_send( 4)='snt_feom' ! 4. snow thickness [m]                ->
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
   call exchange_roots(source_root, target_root, 1, MPI_COMM_FESOM, MPI_COMM_WORLD)
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

  subroutine cpl_oasis3mct_send(ind, data_array, action)
    use o_param
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
    real(kind=WP), intent(IN)       :: data_array(myDim_nod2D+eDim_nod2D)
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

    exfld = cplsnd(ind, 1:myDim_nod2D)/real(o2a_call_count)

    t2=MPI_Wtime()
#ifdef VERBOSE
    if (mype==0) then
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
    
    if (ind==1) then
       time_send(1)=t3-t1       
       time_send(2)=t2-t1
    else
       time_send(1)=time_send(1)+t3-t1
       time_send(2)=time_send(2)+t2-t1
    endif	         
  end subroutine cpl_oasis3mct_send

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine cpl_oasis3mct_recv(ind, data_array, action)
    use o_param
    use g_comm_auto
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
    if (mype==0) then
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
      data_array(1:myDim_nod2d) = exfld
      call exchange_nod(data_array)
   end if   
   t3=MPI_Wtime()
   if (ind==1) then
      time_recv(1)=t3-t1
      time_recv(2)=t3-t2
   else      
      time_recv(1)=time_recv(1)+t3-t1
      time_recv(2)=time_recv(2)+t3-t2
   endif
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
		
        INTEGER, INTENT(IN) :: il_side
        INTEGER, INTENT(IN) :: local_comm, global_comm
        INTEGER, INTENT(OUT) :: source_root, target_root

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
