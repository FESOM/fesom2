    if( mype == root_rank ) then
      allocate(recvbuf(remPtr_nod2D(npes)+rank0Dim_nod2D))
      request_index = 1
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_node_count = rank0Dim_nod2D
        else
          remote_node_count = remPtr_nod2D(remote_rank+1) - remPtr_nod2D(remote_rank)
        endif

        if(remote_rank == 0) then
          call mpi_irecv(recvbuf(remPtr_nod2D(npes)), remote_node_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
        else
          call mpi_irecv(recvbuf(remPtr_nod2D(remote_rank)), remote_node_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
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
      allocate(sendbuf(myDim_nod2D))
      sendbuf(1:myDim_nod2D) = arr2D(1:myDim_nod2D)
      call mpi_send(sendbuf, myDim_nod2D, mpi_precision, root_rank, tag, io_comm, MPIerr)
    end if
