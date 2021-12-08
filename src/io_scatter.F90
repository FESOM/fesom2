module io_scatter_module
  implicit none
  public scatter_nod2D, scatter_elem2D
  private
  
contains
  

  ! thread-safe procedure
  subroutine scatter_nod2D(arr2D_global, arr2D_local, root_rank, comm)
    use g_PARSUP
    use o_mesh
    use, intrinsic :: iso_fortran_env, only: real64
    real(real64), intent(in) :: arr2D_global(:)
    real(real64), intent(out)  :: arr2D_local(:)
    integer, intent(in) :: root_rank ! rank of sending process
    integer, intent(in) :: comm
    ! EO args
    integer :: tag = 0
    integer :: mpi_precision = MPI_DOUBLE_PRECISION
    integer status(MPI_STATUS_SIZE)
    integer :: remote_rank
    integer, allocatable :: remote_list_nod2d(:)
    real(real64), allocatable :: sendbuf(:)
    integer node_size

    call assert(size(arr2D_local) == size(mylist_nod2d), __LINE__) ! == mydim_nod2d+edim_nod2d, i.e. partition nodes + halo nodes

    if(mype == root_rank) then
      arr2D_local = arr2D_global(mylist_nod2d)
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle
        
        ! receive remote partition 2D size
        call mpi_recv(node_size, 1, mpi_integer, remote_rank, tag+0, comm, status, mpierr)

        ! receive remote mylist_nod2d
        allocate(remote_list_nod2d(node_size))
        call mpi_recv(remote_list_nod2d(1), node_size, mpi_integer, remote_rank, tag+1, comm, status, mpierr)

        allocate(sendbuf(node_size))
        sendbuf = arr2D_global(remote_list_nod2d)
        deallocate(remote_list_nod2d)

        call mpi_send(sendbuf(1), node_size, mpi_precision, remote_rank, tag+2, comm, mpierr)
        deallocate(sendbuf)
      end do
    
    else  
      node_size = size(mylist_nod2d)
      call mpi_send(node_size, 1, mpi_integer, root_rank, tag+0, comm, mpierr)
      call mpi_send(mylist_nod2d(1), node_size, mpi_integer, root_rank, tag+1, comm, mpierr)
      
      call mpi_recv(arr2D_local(1), node_size, mpi_precision, root_rank, tag+2, comm, status, mpierr) ! aleph blocks here
    end if

    ! without a barrier, we get wrong results in arr2D_local
    ! todo: not sure why this happens (probably because the 3D levels have the same send/recv signature), get rid of the barrier if possible
    call mpi_barrier(comm, mpierr)
  end subroutine


  ! thread-safe procedure
  subroutine scatter_elem2D(arr2D_global, arr2D_local, root_rank, comm)
    use g_PARSUP
    use o_mesh
    use, intrinsic :: iso_fortran_env, only: real64
    real(real64), intent(in) :: arr2D_global(:)
    real(real64), intent(out)  :: arr2D_local(:)
    integer, intent(in) :: root_rank ! rank of sending process
    integer, intent(in) :: comm
    ! EO args
    integer :: tag = 0
    integer :: mpi_precision = MPI_DOUBLE_PRECISION
    integer status(MPI_STATUS_SIZE)
    integer :: remote_rank
    integer, allocatable :: remote_list_elem2d(:)
    real(real64), allocatable :: sendbuf(:)
    integer elem_size

    elem_size = size(arr2D_local)
    call assert(elem_size == mydim_elem2d+edim_elem2d, __LINE__) ! mylist_elem2d is larger and can not be used for comparison here

    if(mype == root_rank) then
      arr2D_local = arr2D_global(myList_elem2D(1:elem_size))
      do remote_rank = 0, npes-1
        if(remote_rank == root_rank) cycle

        ! receive remote partition 2D size
        call mpi_recv(elem_size, 1, mpi_integer, remote_rank, tag+0, comm, status, mpierr)

        ! receive remote mylist_elem2d
        allocate(remote_list_elem2d(elem_size))
        call mpi_recv(remote_list_elem2d(1), elem_size, mpi_integer, remote_rank, tag+1, comm, status, mpierr)

        allocate(sendbuf(elem_size))
        sendbuf = arr2D_global(remote_list_elem2d)
        deallocate(remote_list_elem2d)

        call mpi_send(sendbuf(1), elem_size, mpi_precision, remote_rank, tag+2, comm, mpierr)
        deallocate(sendbuf)
      end do
    
    else  
      call mpi_send(elem_size, 1, mpi_integer, root_rank, tag+0, comm, mpierr)
      call mpi_send(mylist_elem2d(1), elem_size, mpi_integer, root_rank, tag+1, comm, mpierr)
      
      call mpi_recv(arr2D_local(1), elem_size, mpi_precision, root_rank, tag+2, comm, status, mpierr)
    end if
    
    ! without a barrier, we get wrong results in arr2D_local
    ! todo: not sure why this happens (probably because the 3D levels have the same send/recv signature), get rid of the barrier if possible
    call mpi_barrier(comm, mpierr)
  end subroutine


  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO parameters
    if(.not. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine

end module

