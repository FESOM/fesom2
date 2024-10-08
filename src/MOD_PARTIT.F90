!==========================================================
! Variables to organize parallel work
module MOD_PARTIT
USE O_PARAM
USE, intrinsic :: ISO_FORTRAN_ENV, only : int32
USE MOD_WRITE_BINARY_ARRAYS
USE MOD_READ_BINARY_ARRAYS
#if defined(_OPENMP)
    USE OMP_LIB
#endif
IMPLICIT NONE
SAVE
include 'mpif.h'
integer, parameter   :: MAX_LAENDERECK=16
integer, parameter   :: MAX_NEIGHBOR_PARTITIONS=32


type com_struct
     integer                                       :: rPEnum ! the number of PE I receive info from
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: rPE    ! their list
     integer, dimension(MAX_NEIGHBOR_PARTITIONS+1) :: rptr   ! allocatables to the list of nodes
     integer, dimension(:), allocatable            :: rlist  ! the list of nodes
     integer                                       :: sPEnum ! send part
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sPE
     integer, dimension(MAX_NEIGHBOR_PARTITIONS)   :: sptr
     integer, dimension(:), allocatable            :: slist
     integer, dimension(:), allocatable            :: req    ! request for MPI_Wait
     integer                                       :: nreq   ! number of requests for MPI_Wait
                                                             ! (to combine halo exchange of several fields)
     contains
     procedure READ_T_COM_STRUCT
     procedure WRITE_T_COM_STRUCT
end type com_struct

TYPE T_PARTIT
  
  !---------------------------------------------------
  !LA 2023-01-31 add asynchronous icebergs
  ! kh 10.02.21 communicator for async iceberg computations based on OpenMP
  integer              :: MPI_COMM_FESOM_IB
  !---------------------------------------------------

  type(com_struct) :: com_nod2D
  type(com_struct) :: com_elem2D
  type(com_struct) :: com_elem2D_full

  !---------------------------------------------------
  !LA 2023-01-31 add asynchronous icebergs
  ! kh 11.02.21
  integer            :: MPIERR_IB
  !---------------------------------------------------
  integer              :: npes
  integer              :: mype
  integer              :: maxPEnum=100
  integer, allocatable, dimension(:)  :: part

  integer                             :: myDim_nod2D
  integer                             :: eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D

  integer                             :: myDim_elem2D, myDim_elem2D_shrinked
  integer                             :: eDim_elem2D
  integer                             :: eXDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D, myInd_elem2D_shrinked

  integer                             :: myDim_edge2D
  integer                             :: eDim_edge2D
  integer, allocatable, dimension(:)  :: myList_edge2D
  integer :: pe_status = 0 ! if /=0 then something is wrong

  integer              :: MPI_COMM_FESOM ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)
  integer              :: MPI_COMM_WORLD ! FESOM communicator (for ocean only runs if often a copy of MPI_COMM_WORLD)

  ! MPI Datatypes for interface exchange
  ! Element fields (2D; 2D integer; 3D with nl-1 or nl levels, 1 - 4 values)
  !                 small halo and / or full halo
  !!! s(r)_mpitype_* are constructed during the runtime ans shall not be dumped!!!
  integer, allocatable :: s_mpitype_elem2D(:,:),       r_mpitype_elem2D(:,:)
  integer, allocatable :: s_mpitype_elem2D_full_i(:),  r_mpitype_elem2D_full_i(:)
  integer, allocatable :: s_mpitype_elem2D_full(:,:),  r_mpitype_elem2D_full(:,:)
  integer, allocatable :: s_mpitype_elem3D(:,:,:),     r_mpitype_elem3D(:,:,:)
  integer, allocatable :: s_mpitype_elem3D_full(:,:,:),r_mpitype_elem3D_full(:,:,:)

  ! Nodal fields (2D; 2D integer; 3D with nl-1 or nl levels, one, two, or three values)
  integer, allocatable       :: s_mpitype_nod2D(:),     r_mpitype_nod2D(:)
  integer, allocatable       :: s_mpitype_nod2D_i(:),   r_mpitype_nod2D_i(:)
  integer, allocatable       :: s_mpitype_nod3D(:,:,:), r_mpitype_nod3D(:,:,:)

  integer            :: MPIERR
  
  !!! remPtr_* are constructed during the runtime and shall not be dumped!!!
  integer, allocatable ::  remPtr_nod2D(:),  remList_nod2D(:)
  integer, allocatable ::  remPtr_elem2D(:), remList_elem2D(:)

  logical :: elem_full_flag
#if defined(_OPENMP)
  !!! plock is constructed during the runtime and shall not be dumped!!!
    integer(omp_lock_kind), allocatable :: plock(:)
#endif

  contains
#if defined(__PGI)
  procedure,private  READ_T_PARTIT
  procedure,private  WRITE_T_PARTIT
#else
  procedure READ_T_PARTIT
  procedure WRITE_T_PARTIT
#endif
  generic :: read(unformatted)  =>  READ_T_PARTIT
  generic :: write(unformatted) =>  WRITE_T_PARTIT
END TYPE T_PARTIT
contains

! Unformatted writing for COM_STRUCT TYPE
subroutine WRITE_T_COM_STRUCT(tstruct, unit)
    IMPLICIT NONE
    class(COM_STRUCT),    intent(in)     :: tstruct
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg

    write(unit, iostat=iostat, iomsg=iomsg) tstruct%rPEnum
    call write1d_int_static(tstruct%rPE,   unit, iostat, iomsg)
    call write1d_int_static(tstruct%rptr,  unit, iostat, iomsg)
    call write_bin_array(tstruct%rlist, unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) tstruct%sPEnum
    call write1d_int_static(tstruct%sPE,   unit, iostat, iomsg)
    call write1d_int_static(tstruct%sptr,  unit, iostat, iomsg)
    call write_bin_array(tstruct%slist, unit, iostat, iomsg)
    ! req is constructed during the runtime
    ! call write_bin_array(tstruct%req,   unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) tstruct%nreq
end subroutine WRITE_T_COM_STRUCT

subroutine READ_T_COM_STRUCT(tstruct, unit)
    IMPLICIT NONE
    class(COM_STRUCT),    intent(inout)  :: tstruct
    integer,              intent(in)     :: unit
    integer                              :: iostat
    character(len=1024)                  :: iomsg

    read(unit, iostat=iostat, iomsg=iomsg) tstruct%rPEnum
    call read1d_int_static(tstruct%rPE,   unit, iostat, iomsg)
    call read1d_int_static(tstruct%rptr,  unit, iostat, iomsg)
    call read_bin_array(tstruct%rlist, unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) tstruct%sPEnum
    call read1d_int_static(tstruct%sPE,   unit, iostat, iomsg)
    call read1d_int_static(tstruct%sptr,  unit, iostat, iomsg)
    call read_bin_array(tstruct%slist, unit, iostat, iomsg)
!   req is constructed during the runtime
!   call read_bin_array(tstruct%req,   unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) tstruct%nreq
end subroutine READ_T_COM_STRUCT

! Unformatted writing for T_PARTIT
subroutine WRITE_T_PARTIT(partit, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_PARTIT),      intent(in)     :: partit
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call partit%com_nod2D%WRITE_T_COM_STRUCT(unit)
    call partit%com_elem2D%WRITE_T_COM_STRUCT(unit)
    call partit%com_elem2D_full%WRITE_T_COM_STRUCT(unit)

    write(unit, iostat=iostat, iomsg=iomsg) partit%npes
    write(unit, iostat=iostat, iomsg=iomsg) partit%mype
    write(unit, iostat=iostat, iomsg=iomsg) partit%maxPEnum
!PS     write(unit, iostat=iostat, iomsg=iomsg) partit%flag_debug
    call write_bin_array(partit%part,           unit, iostat, iomsg)

    write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_nod2D
    write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_nod2D
    call write_bin_array(partit%myList_nod2D,   unit, iostat, iomsg)

    write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_elem2D
    write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_elem2D
    write(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_elem2D
    call write_bin_array(partit%myList_elem2D,  unit, iostat, iomsg)

    write(unit, iostat=iostat, iomsg=iomsg) partit%myDim_edge2D
    write(unit, iostat=iostat, iomsg=iomsg) partit%eDim_edge2D
    call write_bin_array(partit%myList_edge2D,  unit, iostat, iomsg)
    write(unit, iostat=iostat, iomsg=iomsg) partit%pe_status
end subroutine WRITE_T_PARTIT
! Unformatted reading for T_PARTIT
subroutine READ_T_PARTIT(partit, unit, iostat, iomsg)
    IMPLICIT NONE
    class(T_PARTIT),      intent(inout)  :: partit
    integer,              intent(in)     :: unit
    integer,              intent(out)    :: iostat
    character(*),         intent(inout)  :: iomsg

    call partit%com_nod2D%READ_T_COM_STRUCT(unit)
    call partit%com_elem2D%READ_T_COM_STRUCT(unit)
    call partit%com_elem2D_full%READ_T_COM_STRUCT(unit)

    read(unit, iostat=iostat, iomsg=iomsg) partit%npes
    read(unit, iostat=iostat, iomsg=iomsg) partit%mype
    read(unit, iostat=iostat, iomsg=iomsg) partit%maxPEnum
!PS     read(unit, iostat=iostat, iomsg=iomsg) partit%flag_debug
    call read_bin_array(partit%part,           unit, iostat, iomsg)

    read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_nod2D
    read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_nod2D
    call read_bin_array(partit%myList_nod2D,   unit, iostat, iomsg)

    read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_elem2D
    read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_elem2D
    read(unit, iostat=iostat, iomsg=iomsg) partit%eXDim_elem2D
    call read_bin_array(partit%myList_elem2D,  unit, iostat, iomsg)

    read(unit, iostat=iostat, iomsg=iomsg) partit%myDim_edge2D
    read(unit, iostat=iostat, iomsg=iomsg) partit%eDim_edge2D
    call read_bin_array(partit%myList_edge2D,  unit, iostat, iomsg)
    read(unit, iostat=iostat, iomsg=iomsg) partit%pe_status
end subroutine READ_T_PARTIT

end module MOD_PARTIT
