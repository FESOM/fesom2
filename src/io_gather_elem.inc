    if( mype == root_rank ) then
      allocate(recvbuf(remPtr_elem2D(npes)+rank0Dim_elem2D))
      request_index = 1
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_elem_count = rank0Dim_elem2D
        else
          remote_elem_count = remPtr_elem2D(remote_rank+1) - remPtr_elem2D(remote_rank)
        endif

        if(remote_rank == 0) then
          call mpi_irecv(recvbuf(remPtr_elem2D(npes)), remote_elem_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
        else
          call mpi_irecv(recvbuf(remPtr_elem2D(remote_rank)), remote_elem_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
        endif
        request_index = request_index + 1
      end do
    
      call mpi_waitall(size(req), req, MPI_STATUSES_IGNORE, MPIerr)    
    
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) then
          arr2D_global(myList_elem2D(1:myDim_elem2D)) = arr2D(1:myDim_elem2D) ! local data
        else if(remote_rank == 0) then
          arr2D_global(rank0List_elem2D(1:rank0Dim_elem2D)) = recvbuf(remPtr_elem2D(npes):remPtr_elem2D(npes)+rank0Dim_elem2D-1) ! rank 0 data
        else
          arr2D_global(remList_elem2D(remPtr_elem2D(remote_rank):remPtr_elem2D(remote_rank+1)-1)) = recvbuf(remPtr_elem2D(remote_rank):remPtr_elem2D(remote_rank+1)-1) ! data of any rank!=0
        end if
      end do
    else
      allocate(sendbuf(myDim_elem2D))
      sendbuf(1:myDim_elem2D) = arr2D(1:myDim_elem2D)
      call mpi_send(sendbuf, myDim_elem2D, mpi_precision, root_rank, tag, io_comm, MPIerr)
    end if
