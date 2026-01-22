module io_netcdf_workaround_module
  private
  public next_io_rank
  
  integer, parameter :: SEQUENTIAL_IO_RANK = 0
contains


  integer function next_io_rank(communicator, async_netcdf_allowed) result(result)
    use g_PARSUP
    use mpi_topology_module
    integer, intent(in) :: communicator
    logical, intent(out) :: async_netcdf_allowed
    ! EO args
    integer rank_use_count
    integer rank
    
    result = next_io_rank_helper(communicator, rank_use_count)
    if(rank_use_count > 1) then
      if(mype == SEQUENTIAL_IO_RANK) print *,"rejecting additional async NetCDF for process:",result, "use count:", rank_use_count, "falling back to sequential I/O on process ",SEQUENTIAL_IO_RANK
      result = SEQUENTIAL_IO_RANK
      async_netcdf_allowed = .false.
    else
      async_netcdf_allowed = .true.
    end if    
  end function


  integer recursive function next_io_rank_helper(communicator, rank_use_count) result(result)
    use mpi_topology_module
    integer, intent(in) :: communicator
    integer, intent(out) :: rank_use_count
    ! EO args
    integer rank
    
    rank = mpi_topology%next_host_head_rank(communicator, rank_use_count)

    ! we want to skip rank SEQUENTIAL_IO_RANK as we already do serial netcdf I/O calls from it
    if(rank == SEQUENTIAL_IO_RANK) then
      ! todo: recursion will never end if SEQUENTIAL_IO_RANK is our only rank
      rank = next_io_rank_helper(communicator, rank_use_count)
    end if
    
    result = rank
  end function

end module
