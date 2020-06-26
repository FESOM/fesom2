module io_gather_module
  implicit none
  public gather_nod2D, gather_real4_nod2D, gather_elem2D, gather_real4_elem2D
  private


contains


  subroutine gather_nod2D(arr2D, arr2D_global, root_rank, mesh)
    use g_PARSUP
    use mod_mesh
    implicit none
    type(t_mesh), intent(in) :: mesh
    real(real64), intent(in)  :: arr2D(myDim_nod2D+eDim_nod2D)
    real(real64), intent(out) :: arr2D_global(mesh%nod2D)
    integer, intent(in) :: root_rank ! rank of receiving process
    ! EO args
    integer  ::  remote_rank = -1
    integer :: remote_node_count = -1
    real(real64)              :: sendbuf(myDim_nod2D)
    real(real64)              :: recvbuf(mesh%nod2D) ! todo: alloc only for root_rank
    integer                   :: req(npes-1)
    integer :: request_index
    integer :: rank0Dim_nod2D
    integer, allocatable, dimension(:) :: rank0List_nod2D

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

    if( mype == root_rank ) then
      request_index = 1
      rank0Dim_nod2D = mesh%nod2D - remPtr_nod2D(npes) +1
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_node_count = rank0Dim_nod2D
        else
          remote_node_count = remPtr_nod2D(remote_rank+1) - remPtr_nod2D(remote_rank)
        endif

        if(remote_rank == 0) then
          call mpi_irecv(recvbuf(remPtr_nod2D(npes)), remote_node_count, MPI_DOUBLE_PRECISION, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
        else
          call mpi_irecv(recvbuf(remPtr_nod2D(remote_rank)), remote_node_count, MPI_DOUBLE_PRECISION, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
        endif
        request_index = request_index + 1
      end do
    
      call mpi_waitall(size(req), req, MPI_STATUSES_IGNORE, MPIerr)    
    
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) then
          arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D) ! local data
        else if(remote_rank == 0) then
          arr2D_global(rank0List_nod2D(1:rank0Dim_nod2D)) = recvbuf(remPtr_nod2D(npes):remPtr_nod2D(npes)+rank0Dim_nod2D-1) ! rank 0 data
        else
          arr2D_global(remList_nod2D(remPtr_nod2D(remote_rank):remPtr_nod2D(remote_rank+1)-1)) = recvbuf(remPtr_nod2D(remote_rank):remPtr_nod2D(remote_rank+1)-1) ! data of any rank!=0
        end if
      end do
    else
      sendbuf(1:myDim_nod2D) = arr2D(1:myDim_nod2D)
      call mpi_send(sendbuf, myDim_nod2D, MPI_DOUBLE_PRECISION, root_rank, 2, MPI_COMM_FESOM, MPIerr)
    end if
  
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


subroutine gather_elem2D(arr2D, arr2D_global)

! Make element information available to master PE 

use g_PARSUP
USE o_MESH

IMPLICIT NONE

integer      :: n

real(real64) ::  arr2D(:)
real(real64) ::  arr2D_global(:)
real(real64), allocatable :: recvbuf(:)
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
         call MPI_IRECV(recvbuf(start), e2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
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
   
   call MPI_SEND( arr2D, myDim_elem2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
   
ENDIF
end if

end subroutine gather_elem2D


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

