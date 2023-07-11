module forcing_provider_async_module
  use forcing_lookahead_reader_module
  use async_threads_module
  implicit none
  public forcing_provider
  private
  

  type forcing_provider_type
    private
    contains
    procedure, public, nopass :: get_forcingdata
  end type
  type(forcing_provider_type), save :: forcing_provider ! handle to singleton instance

  
  type, private :: forcing_async_reader_type
    character(:), allocatable :: varname
    integer fileyear
    type(forcing_lookahead_reader_type) reader_a, reader_b
    type(forcing_lookahead_reader_type), pointer :: reader_current , reader_next
    integer :: netcdf_timestep_size = -1
    type(thread_type) thread
    character(:), allocatable :: thread_filepath
    integer thread_timeindex
    contains
    procedure initialize_async_reader
  end type
  type(forcing_async_reader_type), allocatable, save, target :: all_readers(:) ! we can not put this inside of the forcing_provider_type as we must have it as a target to assign the current/next pointers (:sic:)
  
  
  contains


  subroutine get_forcingdata(varcount, varindex, async_netcdf_allowed, filepath, fileyear, varname, time_index, forcingdata)
    integer, intent(in) :: varcount
    integer, intent(in) :: varindex ! todo: remove this arg and just use a hashmap for varname
    logical async_netcdf_allowed
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_index
    real(4), intent(out) :: forcingdata(:,:)
    ! EO args
    type(forcing_lookahead_reader_type) new_reader_a, new_reader_b
    
    ! init our all_readers array if not already done
    ! dynamically increasing this array did not work because somehow the move_alloc messes with our pointer assignments (maybe it would work with a manual deep-copy instead?)
    ! as a workaround, we allocate the array once to a fixed size
    if(.not. allocated(all_readers)) then
      allocate(all_readers(varcount))
    else
      call assert(size(all_readers) == varcount, __LINE__)
    end if
    
    if( all_readers(varindex)%netcdf_timestep_size == -1 ) then ! reader has never been initialized
      all_readers(varindex)%netcdf_timestep_size = 0
      call all_readers(varindex)%initialize_async_reader(filepath, fileyear, varname)
      
      ! attach thread for this forcing field to the c++ library
      all_readers(varindex)%thread_filepath = filepath
      all_readers(varindex)%thread_timeindex = -1
      call all_readers(varindex)%thread%initialize(thread_callback, varindex)
      if(.not. async_netcdf_allowed) call all_readers(varindex)%thread%disable_async()
      
    else if(fileyear /= all_readers(varindex)%fileyear) then
      ! stop the thread, close our reader and create a new one
      if(all_readers(varindex)%thread_timeindex == time_index) then
        call all_readers(varindex)%thread%join()
        all_readers(varindex)%thread_timeindex = -1
      end if
      
      call all_readers(varindex)%reader_a%finalize()
      call all_readers(varindex)%reader_b%finalize()
      all_readers(varindex)%reader_a = new_reader_a
      all_readers(varindex)%reader_b = new_reader_b

      call all_readers(varindex)%initialize_async_reader(filepath, fileyear, varname)
    end if

    call assert(size(all_readers)>=varindex, __LINE__)
    
    ! join thread
    if(all_readers(varindex)%thread_timeindex == time_index) then
      call all_readers(varindex)%thread%join()
      all_readers(varindex)%thread_timeindex = -1
    end if

    call all_readers(varindex)%reader_current%yield_data(time_index, forcingdata)

    ! todo: kick off the thread before we fill the forcingdata array
    if(all_readers(varindex)%thread_timeindex == -1) then ! we did not already kick off a thread
      if(.not. all_readers(varindex)%reader_next%timeindex_in_cache_bounds(time_index+PREFETCH_SIZE)) then
        if(all_readers(varindex)%netcdf_timestep_size >= time_index+PREFETCH_SIZE) then
          ! prefetch the next timestep asynchronously
          call assert(all_readers(varindex)%thread_timeindex == -1, __LINE__)
          all_readers(varindex)%thread_timeindex = time_index+PREFETCH_SIZE
          call all_readers(varindex)%thread%run()
        end if
      end if
    end if
    
  end subroutine


  subroutine thread_callback(index)
    integer, intent(in) :: index
    ! EO args
    call all_readers(index)%reader_next%timeindex_hint(all_readers(index)%thread_timeindex)
  end subroutine


  subroutine initialize_async_reader(this, filepath, fileyear, varname)
    class(forcing_async_reader_type), target, intent(inout) :: this
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    character(len=*), intent(in) :: varname
    ! EO args
 
    this%varname = varname
    this%fileyear = fileyear

    call this%reader_a%initialize(filepath, fileyear, varname)
    call this%reader_b%initialize(filepath, fileyear, varname)
          
    this%netcdf_timestep_size = this%reader_a%netcdf_timestep_size()

    this%reader_current => this%reader_a
    this%reader_next => this%reader_b
    
    ! todo: this%thread%init should be called here, but we can not call it multiple times for the same forcing field
  end subroutine

  
  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine
end module
