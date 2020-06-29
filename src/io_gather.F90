module io_gather_module
  implicit none
  public gather_nod2D, gather_real4_nod2D, gather_elem2D, gather_real4_elem2D
  private

  logical, save :: nod2D_lists_initialized = .false.
  integer, save :: rank0Dim_nod2D
  integer, save, allocatable, dimension(:) :: rank0List_nod2D

  logical, save :: elem2D_lists_initialized = .false.
  integer, save :: rank0Dim_elem2D
  integer, save, allocatable, dimension(:) :: rank0List_elem2D
  
contains


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
  

  subroutine gather_nod2D(arr2D, arr2D_global, root_rank, mesh)
    use g_PARSUP
    use mod_mesh
    implicit none
    type(t_mesh), intent(in) :: mesh
    real(real64), intent(in)  :: arr2D(:)
    real(real64), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_node_count = -1
    real(real64), allocatable :: sendbuf(:)
    real(real64), allocatable :: recvbuf(:) ! todo: alloc only for root_rank
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_DOUBLE_PRECISION

    if(.not. nod2D_lists_initialized) call init_nod2D_lists()

    include "io_gather_nod.inc"  
  end subroutine


subroutine gather_real4_nod2D(arr2D, arr2D_global)

! Make nodal information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, n2D

 if (npes> 1) then

CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then
      allocate(recvbuf(ubound(arr2D_global,1)))
      do  n = 1, npes-1
         n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
         start = remPtr_nod2D(n)
         call MPI_IRECV(recvbuf(start), n2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                       = recvbuf(1 : remPtr_nod2D(npes)-1)
      deallocate(recvbuf)
   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_nod2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF

endif
end subroutine gather_real4_nod2D


  subroutine gather_elem2D(arr2D, arr2D_global, root_rank, mesh)
    use g_PARSUP
    use mod_mesh
    implicit none
    type(t_mesh), intent(in) :: mesh
    real(real64), intent(in)  :: arr2D(:)
    real(real64), intent(out) :: arr2D_global(:)
    integer, intent(in) :: root_rank ! rank of receiving process
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_elem_count = -1
    real(real64), allocatable              :: sendbuf(:)
    real(real64), allocatable              :: recvbuf(:)
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: mpi_precision = MPI_DOUBLE_PRECISION

    if(.not. elem2D_lists_initialized) call init_elem2D_lists()

    include "io_gather_elem.inc"
  end subroutine


subroutine gather_real4_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real32)  ::  arr2D(:)
real(real32)  ::  arr2D_global(:)
real(real32), allocatable :: recvbuf(:)
integer        :: req(npes-1)
integer        :: start, e2D


 if (npes> 1) then
CALL MPI_BARRIER(MPI_COMM_FESOM,MPIerr)

! Consider MPI-datatypes to recv directly into arr2D_global!

IF ( mype == 0 ) THEN
   
   if (npes>1) then

      allocate(recvbuf(remPtr_elem2D(npes)))

      do  n = 1, npes-1
         e2D   = remPtr_elem2D(n+1) - remPtr_elem2D(n)
         start = remPtr_elem2D(n)
         call MPI_IRECV(recvbuf(start), e2D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
      enddo
      
      arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D)
   
      call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)
   
      arr2D_global(remList_elem2D(1 : remPtr_elem2D(npes)-1)) &
                       = recvbuf(1 : remPtr_elem2D(npes)-1)

      deallocate(recvbuf)

   else

      arr2D_global(:) = arr2D(:)
     
   endif

ELSE
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF
end if

end subroutine gather_real4_elem2D

end module

