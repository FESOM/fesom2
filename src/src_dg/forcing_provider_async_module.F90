module forcing_provider_async_module
  use iso_c_binding
  use forcing_lookahead_reader_module
  implicit none
  public forcing_provider
  private

    
  type cpp_thread
    private
    character(:), allocatable :: filepath
    character(:), allocatable :: varname
    integer timeindex
    contains
    procedure, public :: cpp_thread_begin
    procedure, public :: cpp_thread_end
    procedure, private :: init
  end type
  

  type forcing_provider_type
    private
    contains
    procedure, public :: get_forcingdata
  end type
  type(forcing_provider_type), save :: forcing_provider ! handle to singleton instance

  
  type, private :: forcing_async_reader_type
    character(:), allocatable :: varname
    integer fileyear
    type(forcing_lookahead_reader_type) reader_a, reader_b
    type(forcing_lookahead_reader_type), pointer :: reader_current , reader_next
    integer :: netcdf_timestep_size = -1
    type(cpp_thread) thread
    contains
    procedure initialize_async_reader
  end type
  type(forcing_async_reader_type), allocatable, save, target :: all_readers(:) ! we can not put this inside of the forcing_provider_type as we must have it as a target to assign the current/next pointers (:sic:)

  character(len=*), parameter :: FILENAMESUFFIX = ".nc"
  
  
  contains


  subroutine get_forcingdata(this, varcount, varindex, filepath, fileyear, varname, time_index, forcingdata)
    class(forcing_provider_type), intent(inout) :: this
    integer, intent(in) :: varcount
    integer, intent(in) :: varindex ! todo: remove this arg and just use a hashmap for varname
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_index
    real(4), intent(out) :: forcingdata(:,:)
    ! EO args
    type(forcing_async_reader_type), allocatable :: tmparr(:)
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
      call all_readers(varindex)%initialize_async_reader(varindex, filepath, fileyear, varname)
      
      ! attach thread for this forcing field to the c++ library
      call all_readers(varindex)%thread%init(varindex, filepath, varname)
    else if(fileyear /= all_readers(varindex)%fileyear) then
      ! stop the thread, close our reader and create a new one
      if(all_readers(varindex)%thread%timeindex == time_index) then
        call all_readers(varindex)%thread%cpp_thread_end(varindex)
      end if
      
      call all_readers(varindex)%reader_a%finalize()
      call all_readers(varindex)%reader_b%finalize()
      all_readers(varindex)%reader_a = new_reader_a
      all_readers(varindex)%reader_b = new_reader_b

      call all_readers(varindex)%initialize_async_reader(varindex, filepath, fileyear, varname)
    end if

call assert(allocated(all_readers), __LINE__)
call assert(size(all_readers)>=varindex, __LINE__)
if(.not. associated(all_readers(varindex)%reader_current, all_readers(varindex)%reader_a)) stop __LINE__
if(.not. associated(all_readers(varindex)%reader_next, all_readers(varindex)%reader_b)) stop __LINE__

    ! join thread
    if(all_readers(varindex)%thread%timeindex == time_index) then
      call all_readers(varindex)%thread%cpp_thread_end(varindex)
    end if

    call all_readers(varindex)%reader_current%yield_data(time_index, forcingdata)

    ! todo: kick off the thread before we fill the forcingdata array
    if(all_readers(varindex)%thread%timeindex == -1) then ! we did not already kick off a thread
      if(.not. all_readers(varindex)%reader_next%timeindex_in_cache_bounds(time_index+PREFETCH_SIZE)) then
        if(all_readers(varindex)%netcdf_timestep_size >= time_index+PREFETCH_SIZE) then
          ! prefetch the next timestep asynchronously
          call all_readers(varindex)%thread%cpp_thread_begin(varindex, time_index+PREFETCH_SIZE)
        end if
      end if
    end if
    
  end subroutine
  
  
  subroutine init(this, varindex, filepath, varname)
    class(cpp_thread), target :: this
    integer, intent(in) :: varindex
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    ! EO args
    
    this%filepath = filepath
    this%varname = varname
    this%timeindex = -1

   call init_ccall(varindex)
  end subroutine


  subroutine cpp_thread_begin(this, varindex, timeindex)
    class(cpp_thread) :: this
    integer, intent(in) :: varindex
    integer, intent(in) :: timeindex
    ! EO args
    call assert(this%timeindex == -1, __LINE__)
    this%timeindex = timeindex
    call begin_ccall(varindex)    
  end subroutine


   subroutine cpp_thread_end(this, varindex)
    class(cpp_thread), intent(inout) :: this
    integer, intent(in) :: varindex
    ! EO args  
    call end_ccall(varindex)
    this%timeindex = -1
  end subroutine


  subroutine fortran_call(index) bind (C, name="fortran_call")
    integer(c_int), intent(in), value :: index
    ! EO args
    call all_readers(index)%reader_next%timeindex_hint(all_readers(index)%thread%timeindex)
  end subroutine


  subroutine initialize_async_reader(this, varindex, filepath, fileyear, varname)
    class(forcing_async_reader_type), target, intent(inout) :: this
    integer, intent(in) :: varindex
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
