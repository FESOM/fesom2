    if( partit%mype == root_rank ) then
      allocate(recvbuf(partit%remPtr_elem2D(partit%npes)+rank0Dim_elem2D))
      request_index = 1
      do remote_rank = 0, partit%npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_elem_count = rank0Dim_elem2D
        else
          remote_elem_count = partit%remPtr_elem2D(remote_rank+1) - partit%remPtr_elem2D(remote_rank)
        endif

        if(remote_rank == 0) then
          call mpi_irecv(recvbuf(partit%remPtr_elem2D(partit%npes)), remote_elem_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), partit%MPIerr)
        else
          call mpi_irecv(recvbuf(partit%remPtr_elem2D(remote_rank)), remote_elem_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), partit%MPIerr)
        endif
        request_index = request_index + 1
      end do
    
      call mpi_waitall(size(req), req, MPI_STATUSES_IGNORE, partit%MPIerr)    
    
      do remote_rank = 0, partit%npes-1
        if(remote_rank == root_rank) then
          arr2D_global(partit%myList_elem2D(1:partit%myDim_elem2D)) = arr2D(1:partit%myDim_elem2D) ! local data
        else if(remote_rank == 0) then
          arr2D_global(rank0List_elem2D(1:rank0Dim_elem2D)) = recvbuf(partit%remPtr_elem2D(partit%npes):partit%remPtr_elem2D(partit%npes)+rank0Dim_elem2D-1) ! rank 0 data
        else
          arr2D_global(partit%remList_elem2D(partit%remPtr_elem2D(remote_rank):partit%remPtr_elem2D(remote_rank+1)-1)) = recvbuf(partit%remPtr_elem2D(remote_rank):partit%remPtr_elem2D(remote_rank+1)-1) ! data of any rank!=0
        end if
      end do
    else
      allocate(sendbuf(partit%myDim_elem2D))
      sendbuf(1:partit%myDim_elem2D) = arr2D(1:partit%myDim_elem2D)
      call mpi_send(sendbuf, partit%myDim_elem2D, mpi_precision, root_rank, tag, io_comm, partit%MPIerr)
    end if
