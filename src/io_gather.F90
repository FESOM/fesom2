module io_gather_module
  implicit none
  public init_io_gather, gather_nod2D, gather_real4_nod2D, gather_elem2D, gather_real4_elem2D
  private

  logical, save :: nod2D_lists_initialized = .false.
  integer, save :: rank0Dim_nod2D
  integer, save, allocatable, dimension(:) :: rank0List_nod2D

  logical, save :: elem2D_lists_initialized = .false.
  integer, save :: rank0Dim_elem2D
  integer, save, allocatable, dimension(:) :: rank0List_elem2D
  
contains


  subroutine init_io_gather()
    integer err

    if(.not. nod2D_lists_initialized) call init_nod2D_lists()
    if(.not. elem2D_lists_initialized) call init_elem2D_lists()
  end subroutine


  subroutine init_nod2D_lists()
    use g_PARSUP
    implicit none
    ! EO args

    ! todo: initialize with the other comm arrays, probably in "init_gatherLists" subroutine 
    if(mype /= 0) then
      if(.not. allocated(remPtr_nod2D)) allocate(remPtr_nod2D(npes))
    end if
    call MPI_Bcast(remPtr_nod2D, size(remPtr_nod2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    if(mype /= 0) then
      if(.not. allocated(remList_nod2D)) allocate(remList_nod2D(remPtr_nod2D(npes)))
    end if
    call MPI_Bcast(remList_nod2D, size(remList_nod2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)

    if(mype == 0) rank0Dim_nod2D = myDim_nod2D
    call mpi_bcast(rank0Dim_nod2D, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)

    if(mype == 0) then
      allocate(rank0List_nod2D(0))
      call mpi_bcast(myList_nod2D, myDim_nod2D, MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    else
      allocate(rank0List_nod2D(rank0Dim_nod2D))
      call mpi_bcast(rank0List_nod2D, size(rank0List_nod2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    end if
    
    nod2D_lists_initialized = .true.
  end subroutine


  subroutine init_elem2D_lists()
    use g_PARSUP
    implicit none
    ! EO args

    ! todo: initialize with the other comm arrays, probably in "init_gatherLists" subroutine 
    if(mype /= 0) then
      if(.not. allocated(remPtr_elem2D)) allocate(remPtr_elem2D(npes))
    end if
    call MPI_Bcast(remPtr_elem2D, size(remPtr_elem2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    if(mype /= 0) then
      if(.not. allocated(remList_elem2D)) allocate(remList_elem2D(remPtr_elem2D(npes)))
    end if
    call MPI_Bcast(remList_elem2D, size(remList_elem2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)

    if(mype == 0) rank0Dim_elem2D = myDim_elem2D
    call mpi_bcast(rank0Dim_elem2D, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)

    if(mype == 0) then
      allocate(rank0List_elem2D(0))
      call mpi_bcast(myList_elem2D, myDim_elem2D, MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    else
      allocate(rank0List_elem2D(rank0Dim_elem2D))
      call mpi_bcast(rank0List_elem2D, size(rank0List_elem2D), MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
    end if

    elem2D_lists_initialized = .true.
  end subroutine
  

  ! thread-safe procedure
  subroutine gather_nod2D(arr2D, arr2D_global, root_rank, tag, io_comm)
    use g_PARSUP
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in)  :: arr2D(:)
    real(real64), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    integer, intent(in) :: tag
    integer io_comm
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_node_count = -1
    real(real64), allocatable :: sendbuf(:)
    real(real64), allocatable :: recvbuf(:) ! todo: alloc only for root_rank
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_DOUBLE_PRECISION

    if(.not. nod2D_lists_initialized) stop "io_gather_module has not been initialized"

    include "io_gather_nod.inc"  
  end subroutine


  ! thread-safe procedure
  subroutine gather_real4_nod2D(arr2D, arr2D_global, root_rank, tag, io_comm)
    use g_PARSUP
    use, intrinsic :: iso_fortran_env, only: real32
    implicit none
    real(real32), intent(in)  :: arr2D(:)
    real(real32), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    integer, intent(in) :: tag
    integer io_comm
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_node_count = -1
    real(real32), allocatable :: sendbuf(:)
    real(real32), allocatable :: recvbuf(:) ! todo: alloc only for root_rank
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_REAL

    if(.not. nod2D_lists_initialized) stop "io_gather_module has not been initialized"

    include "io_gather_nod.inc"  
  end subroutine


  ! thread-safe procedure
  subroutine gather_elem2D(arr2D, arr2D_global, root_rank, tag, io_comm)
    use g_PARSUP
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in)  :: arr2D(:)
    real(real64), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    integer, intent(in) :: tag
    integer io_comm
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_elem_count = -1
    real(real64), allocatable              :: sendbuf(:)
    real(real64), allocatable              :: recvbuf(:)
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_DOUBLE_PRECISION

    if(.not. elem2D_lists_initialized) stop "io_gather_module has not been initialized"

    include "io_gather_elem.inc"
  end subroutine


  ! thread-safe procedure
  subroutine gather_real4_elem2D(arr2D, arr2D_global, root_rank, tag, io_comm)
    use g_PARSUP
    use, intrinsic :: iso_fortran_env, only: real32
    implicit none
    real(real32), intent(in)  :: arr2D(:)
    real(real32), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    integer, intent(in) :: tag
    integer io_comm
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_elem_count = -1
    real(real32), allocatable              :: sendbuf(:)
    real(real32), allocatable              :: recvbuf(:)
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_REAL

    if(.not. elem2D_lists_initialized) stop "io_gather_module has not been initialized"

    include "io_gather_elem.inc"
  end subroutine

end module

