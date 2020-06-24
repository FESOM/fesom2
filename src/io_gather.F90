module io_gather_module
  implicit none
  public gather_nod2D, gather_real4_nod2D, gather_elem2D, gather_real4_elem2D
  private


contains


  subroutine gather_nod2D(arr2D, arr2D_global)
    use g_PARSUP
    use o_MESH
    implicit none
    integer      :: n
    real(real64) ::  arr2D(:)
    real(real64) ::  arr2D_global(:)
    real(real64), allocatable :: recvbuf(:)
    integer        :: req(npes-1)
    integer        :: start, n2D
  
    if ( mype == 0 ) then
      if (npes>1) then
    
        allocate(recvbuf(size(arr2D_global)))
        do  n = 1, npes-1
          n2D   = remPtr_nod2D(n+1) - remPtr_nod2D(n)
          start = remPtr_nod2D(n)
          call mpi_irecv(recvbuf(start), n2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
        end do

        arr2D_global(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)

        call mpi_waitall(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

        arr2D_global(remList_nod2D(1 : remPtr_nod2D(npes)-1)) = recvbuf(1 : remPtr_nod2D(npes)-1)
        deallocate(recvbuf)
    
      else
        arr2D_global(:) = arr2D(:)
      end if
    else
  
      call mpi_send( arr2D, myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )
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

