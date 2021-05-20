module forcing_lookahead_reader_module
  use forcing_provider_netcdf_module
  implicit none
  public forcing_lookahead_reader_type
  private
    
  type forcing_lookahead_reader_type
    private
    integer fileyear_
    integer :: first_stored_timeindex = -1, last_stored_timeindex = -1
    real(4), allocatable :: stored_values(:,:,:)
    integer :: netcdf_timestep_size_ = -1
    type(netcdf_reader_handle) filehandle
    contains
    procedure, public :: initialize, finalize, yield_data, is_initialized, fileyear, netcdf_timestep_size, timeindex_hint, timeindex_in_cache_bounds
    procedure, private :: read_data
  end type

  character(len=*), parameter :: FILENAMESUFFIX = ".nc"
  integer, public, parameter :: PREFETCH_SIZE = 10
  
  contains
  

  pure function fileyear(this) result(x)
    class(forcing_lookahead_reader_type), intent(in) :: this
    integer x
    ! EO args
    x = this%fileyear_
  end function
  

  pure function is_initialized(this) result(x)
    class(forcing_lookahead_reader_type), intent(in) :: this
    logical x
    ! EO args
    x = (this%netcdf_timestep_size_ /= -1)
  end function
  
  
  function netcdf_timestep_size(this) result(t)
    class(forcing_lookahead_reader_type), intent(in) :: this
    integer t
    t = this%netcdf_timestep_size_
  end function
  
  
  subroutine initialize(this, filepath, fileyear, varname)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    integer fileyear
    character(len=*), intent(in) :: varname
    ! EO args
          
    this%fileyear_ = fileyear
    call this%filehandle%initialize(filepath, varname) ! finalize() to close the file
    this%netcdf_timestep_size_ = this%filehandle%timestep_size() 
  end subroutine


  subroutine finalize(this)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    ! EO args
    call this%filehandle%finalize()
  end subroutine
  
  
  pure function timeindex_in_cache_bounds(this, time_index) result(x)
    class(forcing_lookahead_reader_type), intent(in) :: this
    integer, intent(in) :: time_index
    logical x
    ! EO args
    x = all([ this%first_stored_timeindex <= time_index, time_index <= this%last_stored_timeindex ])
  end function

  
  subroutine timeindex_hint(this, time_index)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    integer, intent(in) :: time_index
    ! EO args
      
    if(time_index <= this%netcdf_timestep_size_) then
      call this%read_data(time_index)
    end if
  end subroutine


  subroutine yield_data(this, time_index, values)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    integer, intent(in) :: time_index
    real(4), intent(out) :: values(:,:)
    ! EO args
    integer reader_time_index
      
    if(time_index <= this%netcdf_timestep_size_) then
      call this%read_data(time_index)
    end if
  
    ! check if the outgoing array has the same shape as our data
    reader_time_index = time_index-this%first_stored_timeindex+1
    call assert(allocated(this%stored_values), __LINE__)
    call assert( all( shape(values)==shape(this%stored_values(:,:,reader_time_index)) ), __LINE__ )
    values = this%stored_values(:,:,reader_time_index)    
  end subroutine


  subroutine read_data(this, time_index)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    integer, intent(in) :: time_index
    ! EO args
    
!    assert(time_index >= reader%first_stored_timeindex) ! we do not go back in time
    if(this%last_stored_timeindex < time_index) then
      
      this%first_stored_timeindex = time_index
      this%last_stored_timeindex = time_index+PREFETCH_SIZE
      if(this%last_stored_timeindex > this%netcdf_timestep_size_) then
        this%last_stored_timeindex = this%netcdf_timestep_size_
      end if
      call this%filehandle%read_netcdf_timesteps(this%first_stored_timeindex, this%last_stored_timeindex, this%stored_values)
      call assert(allocated(this%stored_values), __LINE__)    
    end if

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
