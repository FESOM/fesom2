module io_scatter_module
  implicit none
  public scatter_nod2D
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
    integer :: n, sender_rank
    integer, allocatable :: remote_list_nod2d(:)
    real(real64), allocatable :: sendbuf(:)
    integer node_size

    call assert(size(arr2D_local) == size(mylist_nod2d), __LINE__) ! == mydim_nod2d+edim_nod2d, i.e. partition nodes + halo nodes

    if(mype == root_rank) then
      arr2D_local = arr2D_global(mylist_nod2d)
      do  n = 1, npes-1
        ! receive remote partition 2D size
        call mpi_recv(node_size, 1, mpi_integer, MPI_ANY_SOURCE, tag+0, comm, status, mpierr)
        sender_rank = status(mpi_source)

        ! receive remote mylist_nod2d
        allocate(remote_list_nod2d(node_size))
        call mpi_recv(remote_list_nod2d(1), node_size, mpi_integer, sender_rank, tag+1, comm, status, mpierr)

        allocate(sendbuf(node_size))
        sendbuf = arr2D_global(remote_list_nod2d)
        deallocate(remote_list_nod2d)

        call mpi_send(sendbuf(1), node_size, mpi_double_precision, sender_rank, tag+2, comm, mpierr)
        deallocate(sendbuf)
      end do
    
    else  
      node_size = size(mylist_nod2d)
      call mpi_send(node_size, 1, mpi_integer, root_rank, tag+0, comm, mpierr)
      call mpi_send(mylist_nod2d(1), node_size, mpi_integer, root_rank, tag+1, comm, mpierr)
      
      call mpi_recv(arr2D_local(1), node_size, mpi_precision, root_rank, tag+2, comm, status, mpierr)
    end if
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

