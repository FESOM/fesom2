module hostname_sys_module
  contains

  subroutine hostname_sys(hostname)
    character(len=:), allocatable, intent(out) :: hostname
    integer*4 status, hostnm
    
    allocate(character(32) :: hostname) ! platform dependent length in limits.h or call `getconf HOST_NAME_MAX`
    status = hostnm(hostname)
  end subroutine
end module


module mpi_topology_module
! synopsis:
! collectively call mpi_topology%next_host_head_rank to get the first mpi rank of the next compute node (host) within the given communicator
! collectively call mpi_topology%reset_host_head_rank if there is need to start over with the first compute node

  use hostname_sys_module
  implicit none  
  public mpi_topology
  private

  type :: mpi_topology_type
  contains
    procedure, nopass :: next_host_head_rank, reset_host_head_rank, am_i_host_head_rank, set_hostname_strategy
  end type
  type(mpi_topology_type) mpi_topology

  integer, save :: MAXRANK = 0
  integer, save :: STEP = 0
  integer, save :: count = 0
  integer, save :: COMM = -1
  procedure(hostname_interface), pointer, save :: hostname_strategy => hostname_sys
  abstract interface
    subroutine hostname_interface(hostname)
    character(len=:), allocatable, intent(out) :: hostname
    end subroutine
  end interface

contains

  subroutine set_hostname_strategy(strategy)
    procedure(hostname_interface) strategy
    hostname_strategy => strategy
  end subroutine


  ! must be called collectively
  logical function am_i_host_head_rank(communicator) result(result)
    integer, intent(in) :: communicator
    integer rank, ierror
    
    if(communicator .ne. COMM) COMM = learn_topology(communicator)
    
    call MPI_COMM_RANK(communicator, rank, ierror)
    result = mod(rank, STEP)==0
  end function


  subroutine reset_host_head_rank()
    count = 0
  end


  integer recursive function next_host_head_rank(communicator) result(result)
    integer, intent(in) :: communicator
    
    if(communicator .ne. COMM) COMM = learn_topology(communicator)
    
    result = count*STEP
    if(result > MAXRANK) then
      call reset_host_head_rank()
      result = next_host_head_rank(communicator)
    else
      count = count + 1
    end if
  end function


  integer function learn_topology(communicator) result(result)
    use mpi ! should prefer mpi_f08, but it is not available on some older mpi installations
    integer rank, rank_count
    integer ierror
    character(len=:), allocatable :: hostname
    character(len=:), allocatable :: names(:)
    integer i
    integer ranks_per_host
    integer, intent(in) :: communicator

    call reset_host_head_rank()
    result = communicator
  
    call MPI_COMM_RANK(communicator, rank, ierror)
    call MPI_COMM_SIZE(communicator, rank_count, ierror)
    MAXRANK = rank_count-1
    
    call hostname_strategy(hostname)
    if(rank==0) then
      allocate(character(len(hostname)) :: names(rank_count))
    else
      allocate(character(0) :: names(0))
    end if
    call MPI_GATHER(hostname, len(hostname), MPI_CHAR, names, len(hostname), MPI_CHAR, 0, communicator, ierror)
    if(rank==0) then
        ranks_per_host = 1
        do i=1+1, size(names)
          if(hostname == names(i)) then
            ranks_per_host = ranks_per_host+1
          else
            exit
          end if
        end do    
    end if
    deallocate(names)

    call MPI_BCAST(ranks_per_host, 1, MPI_INT, 0, communicator, ierror)
    STEP = ranks_per_host
  end function
end module
