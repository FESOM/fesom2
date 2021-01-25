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
! after all hosts have been used up, we will start over at the first host and increment the rank by 1
! after all ranks have been used up, we will start over at the first host with rank 0
! optional second argument will return the number of times this rank has been returned

  use hostname_sys_module
  implicit none  
  public mpi_topology
  private

  type :: mpi_topology_type
  contains
    procedure, nopass :: next_host_head_rank, set_hostname_strategy, reset_state
  end type
  type(mpi_topology_type) mpi_topology

  logical, save :: IS_STATE_INITIALIZED = .false.

  integer, save :: MAXRANK
  integer, save :: STEP
  integer, save :: count
  integer, save :: lap
  integer, save :: host_use_count
  integer, save :: COMM
  procedure(hostname_interface), pointer, save :: hostname_strategy
  abstract interface
    subroutine hostname_interface(hostname)
    character(len=:), allocatable, intent(out) :: hostname
    end subroutine
  end interface

contains

  subroutine reset_state()
    MAXRANK = 0
    STEP = 0
    count = 0
    lap = 1
    host_use_count = lap
    COMM = -1
    hostname_strategy => hostname_sys
    
    IS_STATE_INITIALIZED = .true.
  end subroutine


  subroutine set_hostname_strategy(strategy)
    procedure(hostname_interface) strategy
    hostname_strategy => strategy
  end subroutine


  integer recursive function next_host_head_rank(communicator, rank_use_count) result(result)
    integer, intent(in) :: communicator
    integer, optional, intent(out) :: rank_use_count
    
    if(.not. IS_STATE_INITIALIZED) call reset_state()
    if(communicator .ne. COMM) COMM = learn_topology(communicator)
    
    result = count*STEP + lap-1
    if(result > MAXRANK) then ! start a new lap
      count = 0
      host_use_count = host_use_count + 1
      lap = lap + 1
      if(lap > STEP) lap = 1 ! start over with the first rank on a host
      result = next_host_head_rank(communicator)
    else
      count = count + 1
    end if
    
    if(present(rank_use_count)) then
      rank_use_count = (host_use_count-1)/STEP +1
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

    count = 0
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
