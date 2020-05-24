#define HGMODETHREADS

module forcing_provider_async_module
  use iso_c_binding
  use forcing_provider_netcdf_module
  implicit none
  public forcing_provider
  private

    
  type cpp_thread
    private
    character(:), allocatable :: filepath
    character(:), allocatable :: varname
    integer timeindex_first, timeindex_last
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

  
  type, private :: forcing_reader_type
    character(:), allocatable :: basepath
    character(:), allocatable :: varname
    integer fileyear
    type(netcdf_reader_handle) filehandle
    integer stored_timeindex_aa, stored_timeindex_bb
    real(4), allocatable :: stored_values_aa(:,:)
    real(4), allocatable :: stored_values_bb(:,:)
    integer, pointer :: stored_timeindex_current, stored_timeindex_next
    real(4), pointer :: stored_values_current(:,:)
    real(4), pointer :: stored_values_next(:,:)
    integer netcdf_timestep_size
    type(cpp_thread) thread
    contains
    final destructor
  end type
  type(forcing_reader_type), allocatable, save, target :: all_readers(:) ! we can not put this inside of the forcing_provider_type as we must have it as a target to assign the current/next pointers (:sic:)

  character(len=*), parameter :: FILENAMESUFFIX = ".nc"
  
  
  contains


  subroutine destructor(this)
    type(forcing_reader_type), intent(inout) :: this
    ! EO args
print *,"destructor forcing_reader_type ",this%varname, __LINE__
  end subroutine


  subroutine get_forcingdata(this, varindex, filepath, fileyear, varname, time_index, forcingdata)
    class(forcing_provider_type), intent(inout) :: this
    integer, intent(in) :: varindex ! todo: remove this arg and just use a hashmap for varname
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_index
    real(4), intent(out) :: forcingdata(:,:)
    ! EO args
    type(forcing_reader_type), allocatable :: tmparr(:)
    character(:), allocatable :: basepath
    real(4), allocatable :: values(:,:,:)
    integer, allocatable, dimension(:) :: dim_sizes
    
    ! init our all_readers array if not already done
!     if(.not. allocated(all_readers)) then
!       allocate(all_readers(varindex))
!     end if
!         
!     if(size(all_readers) < varindex) then
!       allocate( tmparr(varindex) )
!       tmparr(1:size(all_readers)) = all_readers
!       deallocate(all_readers)
!       call move_alloc(tmparr, all_readers)      
!     end if
! somehow the move_alloc messes with our pointer assignments, alloc only once
if(.not. allocated(all_readers)) then
  allocate(all_readers(10)) ! todo: pass max size as argument?
end if
if(size(all_readers) < varindex) stop __LINE__

    
    ! todo: this is probally not a save test, as the default value must not be 0
    if( len(all_readers(varindex)%varname) == 0 ) then ! reader has never been initialized
      basepath=basepath_from_path(filepath, fileyear)
      
      all_readers(varindex)%basepath = basepath
      all_readers(varindex)%varname = varname
      all_readers(varindex)%fileyear = fileyear
      all_readers(varindex)%stored_timeindex_aa = -1
      all_readers(varindex)%stored_timeindex_bb = -1
      call all_readers(varindex)%filehandle%initialize(filepath, varname) ! finalize() to close the file
      
      dim_sizes = all_readers(varindex)%filehandle%varshape      
      
      all_readers(varindex)%netcdf_timestep_size = dim_sizes(3)

      if( allocated(all_readers(varindex)%stored_values_aa) ) then
        if(all(dim_sizes(1:2) /= shape(all_readers(varindex)%stored_values_aa))) then
          deallocate(all_readers(varindex)%stored_values_aa)
        end if
      end if
      if(.not. allocated(all_readers(varindex)%stored_values_aa)) then
        allocate(all_readers(varindex)%stored_values_aa(dim_sizes(1),dim_sizes(2)))
      end if

      if( allocated(all_readers(varindex)%stored_values_bb) ) then
        if(all(dim_sizes(1:2) /= shape(all_readers(varindex)%stored_values_bb))) then
          deallocate(all_readers(varindex)%stored_values_bb)
        end if
      end if
      if(.not. allocated(all_readers(varindex)%stored_values_bb)) then
        allocate(all_readers(varindex)%stored_values_bb(dim_sizes(1),dim_sizes(2)))
      end if
      
      all_readers(varindex)%stored_values_current => all_readers(varindex)%stored_values_aa
      all_readers(varindex)%stored_timeindex_current => all_readers(varindex)%stored_timeindex_aa
      all_readers(varindex)%stored_values_next => all_readers(varindex)%stored_values_bb
      all_readers(varindex)%stored_timeindex_next => all_readers(varindex)%stored_timeindex_bb
      
      call all_readers(varindex)%thread%init(varindex, filepath, varname)      
    end if


call assert(allocated(all_readers), __LINE__)
call assert(size(all_readers)>=varindex, __LINE__)
if(.not. associated(all_readers(varindex)%stored_values_current)) stop __LINE__
if(.not. associated(all_readers(varindex)%stored_values_next)) stop __LINE__
if(.not. associated(all_readers(varindex)%stored_timeindex_current)) stop __LINE__
if(.not. associated(all_readers(varindex)%stored_timeindex_next)) stop __LINE__

    ! join thread
    if(all_readers(varindex)%thread%timeindex_first == time_index) then
#ifdef HGMODETHREADS
print *, "HGMODETHREADS: YES ",__LINE__
      call all_readers(varindex)%thread%cpp_thread_end(varindex)
#else
print *,"HGMODETHREADS: NO ",__LINE__
call fortran_call(varindex)
#endif
    end if
            
    if(all_readers(varindex)%stored_timeindex_current < time_index) then
      ! directly load the requested timestep
      call all_readers(varindex)%filehandle%read_netcdf_timesteps(time_index, time_index, values)
      all_readers(varindex)%stored_values_current = values(:,:,1)

      ! check if the outgoing array has the same shape as our data
      call assert(associated(all_readers(varindex)%stored_values_current), __LINE__)
      call assert( all( shape(forcingdata)==shape(all_readers(varindex)%stored_values_current(:,:)) ), __LINE__ )
      forcingdata = all_readers(varindex)%stored_values_current(:,:)
      all_readers(varindex)%stored_timeindex_current = time_index
           
    else if(all_readers(varindex)%stored_timeindex_current == time_index) then
      ! check if the outgoing array has the same shape as our data
      call assert(associated(all_readers(varindex)%stored_values_current), __LINE__)
      call assert( all( shape(forcingdata)==shape(all_readers(varindex)%stored_values_current(:,:)) ), __LINE__ )
      forcingdata = all_readers(varindex)%stored_values_current(:,:)
    else
      stop __LINE__ ! forcingdata not set
    end if


    if(all_readers(varindex)%stored_timeindex_next /= time_index+1) then
      if(all_readers(varindex)%netcdf_timestep_size > time_index+1) then
        ! prefetch the next timestep asynchronously
        call all_readers(varindex)%thread%cpp_thread_begin(varindex, time_index+1, time_index+1)
        all_readers(varindex)%stored_timeindex_next = time_index+1
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
    this%timeindex_first = -1
    this%timeindex_last = -1
    this%varname = varname

#ifdef HGMODETHREADS
   call init_ccall(varindex)
#endif
  end subroutine


  subroutine cpp_thread_begin(this, varindex, timeindex_first, timeindex_last)
    class(cpp_thread) :: this
    integer, intent(in) :: varindex
    integer, intent(in) :: timeindex_first, timeindex_last
    ! EO args
    this%timeindex_first = timeindex_first
    this%timeindex_last = timeindex_last
#ifdef HGMODETHREADS
   call begin_ccall(varindex)    
#endif
  end subroutine


   subroutine cpp_thread_end(this, varindex)
    class(cpp_thread), intent(in) :: this
    integer, intent(in) :: varindex
    ! EO args  
#ifdef HGMODETHREADS
   call end_ccall(varindex)
#endif
  end subroutine


  subroutine fortran_call(index) bind (C, name="fortran_call")
    integer(c_int), intent(in), value :: index
    ! EO args

call assert(all_readers(index)%thread%timeindex_first == all_readers(index)%thread%timeindex_last, __LINE__) ! we currently read only a single timestep
    if(all_readers(index)%thread%timeindex_first == all_readers(index)%stored_timeindex_aa) then
      call all_readers(index)%filehandle%read_netcdf_timestep_2d(all_readers(index)%thread%timeindex_first, all_readers(index)%stored_values_aa)
    else if(all_readers(index)%thread%timeindex_first == all_readers(index)%stored_timeindex_bb) then
      call all_readers(index)%filehandle%read_netcdf_timestep_2d(all_readers(index)%thread%timeindex_first, all_readers(index)%stored_values_bb)
    else
      stop __LINE__
    end if
  end subroutine
  
  
  function basepath_from_path(filepath, fileyear) result(r)
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    ! EO args
    integer number_of_yeardigits, suffix_size
    character(:), allocatable :: r

    number_of_yeardigits = int(log10(real(fileyear)))+1
    suffix_size = len(FILENAMESUFFIX)
    
    r = filepath(1:len(filepath)-number_of_yeardigits-suffix_size)
  end function

  
  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
      print *, "error in line ",line
      stop 1
    end if
  end subroutine
end module
