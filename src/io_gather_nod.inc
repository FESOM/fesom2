    if( partit%mype == root_rank ) then
      allocate(recvbuf(partit%remPtr_nod2D(partit%npes)+rank0Dim_nod2D))
      request_index = 1
      do remote_rank = 0, partit%npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_node_count = rank0Dim_nod2D
        else
          remote_node_count = partit%remPtr_nod2D(remote_rank+1) - partit%remPtr_nod2D(remote_rank)
        endif

        if(remote_rank == 0) then
          call mpi_irecv(recvbuf(partit%remPtr_nod2D(partit%npes)), remote_node_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), partit%MPIerr)
        else
          call mpi_irecv(recvbuf(partit%remPtr_nod2D(remote_rank)), remote_node_count, mpi_precision, remote_rank, tag, io_comm, req(request_index), partit%MPIerr)
        endif
        request_index = request_index + 1
      end do
    
      call mpi_waitall(size(req), req, MPI_STATUSES_IGNORE, partit%MPIerr)    
    
      do remote_rank = 0, partit%npes-1
        if(remote_rank == root_rank) then
          arr2D_global(partit%myList_nod2D(1:partit%myDim_nod2D)) = arr2D(1:partit%myDim_nod2D) ! local data
        else if(remote_rank == 0) then
          arr2D_global(rank0List_nod2D(1:rank0Dim_nod2D)) = recvbuf(partit%remPtr_nod2D(partit%npes):partit%remPtr_nod2D(partit%npes)+rank0Dim_nod2D-1) ! rank 0 data
        else
          arr2D_global(partit%remList_nod2D(partit%remPtr_nod2D(remote_rank):partit%remPtr_nod2D(remote_rank+1)-1)) = recvbuf(partit%remPtr_nod2D(remote_rank):partit%remPtr_nod2D(remote_rank+1)-1) ! data of any rank!=0
        end if
      end do
    else
      allocate(sendbuf(partit%myDim_nod2D))
      sendbuf(1:partit%myDim_nod2D) = arr2D(1:partit%myDim_nod2D)
      call mpi_send(sendbuf, partit%myDim_nod2D, mpi_precision, root_rank, tag, io_comm, partit%MPIerr)
    end if
