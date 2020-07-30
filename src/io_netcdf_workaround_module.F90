module io_netcdf_workaround_module
  private
  public next_io_rank
contains


  integer function next_io_rank(communicator) result(result)
    use mpi_topology_module
    integer, intent(in) :: communicator
    ! EO args
    integer rank_use_count
    integer rank
    
    result = next_io_rank_helper(communicator, rank_use_count)    
  end function


  integer recursive function next_io_rank_helper(communicator, rank_use_count) result(result)
    use mpi_topology_module
    integer, intent(in) :: communicator
    integer, intent(out) :: rank_use_count
    ! EO args
    integer rank
    
    rank = mpi_topology%next_host_head_rank(communicator, rank_use_count)

    ! we want to skip rank 0 as we already do serial netcdf I/O calls from it
    if(rank == 0) then
      rank = next_io_rank_helper(communicator, rank_use_count)    
    else if(rank_use_count > 1) then
      ! todo: automatically switch off multithreading for all I/O tasks where rank_use_count > 1 and process them on rank 0
      print *, "error in line ",__LINE__, __FILE__, " We can not re-use FESOM2 process",rank," All available FESOM2 processes are being used for a NetCDF I/O task already. Using multiple simultaneous I/O tasks per process may trigger an error in the NetCDF library."
      stop 1
    end if
    
    result = rank
  end function

end module
