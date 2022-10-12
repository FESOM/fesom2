    
    
    if( mype == root_rank ) then
!      allocate(recvbuf(remPtr_elem2D(npes)+rank0Dim_elem2D)) ! -1?
      allocate(recvbuf3D(nlvl,remPtr_elem2D(npes)+rank0Dim_elem2D-1))

!      allocate(recvbufs(remote_elem_count, npes-1))
      request_index = 1
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle
        if(remote_rank == 0) then
          remote_elem_count = rank0Dim_elem2D
        else
          remote_elem_count = remPtr_elem2D(remote_rank+1) - remPtr_elem2D(remote_rank)
        endif

        if(remote_rank == 0) then
          count = remote_elem_count*nlvl

           !rb_ptr2 => recvbuf3D(:,remPtr_elem2D(npes)+remote_elem_count-1)
           call c_f_pointer(c_loc(recvbuf3D(1,remPtr_elem2D(npes))), rb_ptr, [count])
           call mpi_irecv(rb_ptr, count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
!           call mpi_irecv(recvbuf3D(1,remPtr_elem2D(npes)), count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
!          call mpi_irecv(recvbufs(1,request_index), count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
        else
          count = remote_elem_count*nlvl
 
          call c_f_pointer(c_loc(recvbuf3D(1,remPtr_elem2D(remote_rank))), rb_ptr, [count])
          call mpi_irecv(rb_ptr, count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)

!         call mpi_irecv(recvbuf3D(1,remPtr_elem2D(remote_rank)), count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
!          call mpi_irecv(recvbufs(1,request_index), count, mpi_precision, remote_rank, tag, io_comm, req(request_index), MPIerr)
        endif
        request_index = request_index + 1
      end do
    
      call mpi_waitall(size(req), req, MPI_STATUSES_IGNORE, MPIerr)    
    
      request_index = 1
      do remote_rank = 0, npes-1
         do lvl = 1, nlvl      
           if(remote_rank == root_rank) then
             arr3D_global(lvl,myList_elem2D(1:myDim_elem2D)) = arr3D(lvl,1:myDim_elem2D) ! local data
           else if(remote_rank == 0) then
            arr3D_global(lvl,rank0List_elem2D(1:rank0Dim_elem2D)) = recvbuf3D(lvl, remPtr_elem2D(npes):remPtr_elem2D(npes)+rank0Dim_elem2D-1) ! rank 0 data

!             arr3D_global(:,rank0List_elem2D(1:rank0Dim_elem2D)) = recvbufs(remPtr_elem2D(npes):remPtr_elem2D(npes)+rank0Dim_elem2D-1, request_index) ! rank 0 data
           else
             arr3D_global(lvl,remList_elem2D(remPtr_elem2D(remote_rank):remPtr_elem2D(remote_rank+1)-1)) = recvbuf3D(lvl,remPtr_elem2D(remote_rank):remPtr_elem2D(remote_rank+1)-1) ! data of any rank!=0

           end if
         end do
         if(remote_rank /= root_rank) request_index = request_index + 1
      end do
    else
      if(product(shape(arr3D)) /= myDim_elem2D*nlvl) stop "error in arr size"
      call c_f_pointer(c_loc(arr3D), sb_ptr, [product(shape(arr3D))])
      call mpi_send(sb_ptr, myDim_elem2D*nlvl, mpi_precision, root_rank, tag, io_comm, MPIerr)

    end if