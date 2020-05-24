module forcing_lookahead_reader_module
  use forcing_provider_netcdf_module
  implicit none
  public forcing_lookahead_reader_type
  private
    
  type forcing_lookahead_reader_type
    private
    character(:), allocatable :: basepath
    integer fileyear_
    integer first_stored_timeindex, last_stored_timeindex
    real(4), allocatable :: stored_values(:,:,:)
    integer netcdf_timestep_size
    type(netcdf_reader_handle) filehandle
    contains
    procedure, public :: initialize_lookahead, yield_data_lookahead, read_data_lookahead, is_initialized, fileyear
  end type

  character(len=*), parameter :: FILENAMESUFFIX = ".nc"
  integer, parameter :: PREFETCH_SIZE = 10
  
  contains
  

  function fileyear(this) result(x)
    class(forcing_lookahead_reader_type), intent(in) :: this
    integer x
    ! EO args
    x = this%fileyear_
  end function
  

  function is_initialized(this) result(x)
    class(forcing_lookahead_reader_type), intent(in) :: this
    logical x
    ! EO args
    x = (len(this%basepath) /= 0)
  end function
  
  
  subroutine initialize_lookahead(this, filepath, fileyear, varname)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    integer fileyear
    character(len=*), intent(in) :: varname
    ! EO args
          
    this%basepath = basepath_from_path(filepath, fileyear)
    this%fileyear_ = fileyear
    this%first_stored_timeindex = -1
    this%last_stored_timeindex = -1
    call this%filehandle%initialize(filepath, varname) ! finalize() to close the file
    this%netcdf_timestep_size = this%filehandle%timestep_size() 
  end subroutine
  
  
  subroutine yield_data_lookahead(this, time_index, values)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    integer, intent(in) :: time_index
    real(4), intent(out) :: values(:,:)
    ! EO args
    integer reader_time_index
      
    if(time_index <= this%netcdf_timestep_size) then
      call this%read_data_lookahead(time_index)
    end if
  
    ! check if the outgoing array has the same shape as our data
    reader_time_index = time_index-this%first_stored_timeindex+1
    call assert(allocated(this%stored_values), __LINE__)
    call assert( all( shape(values)==shape(this%stored_values(:,:,reader_time_index)) ), __LINE__ )
    values = this%stored_values(:,:,reader_time_index)    
  end subroutine


  subroutine read_data_lookahead(this, time_index)
    class(forcing_lookahead_reader_type), intent(inout) :: this
    integer, intent(in) :: time_index
    ! EO args
    real(4), allocatable :: values(:,:,:)
    
!    assert(time_index >= reader%first_stored_timeindex) ! we do not go back in time
    if(this%last_stored_timeindex < time_index) then
      
      this%first_stored_timeindex = time_index
      this%last_stored_timeindex = time_index+PREFETCH_SIZE
      if(this%last_stored_timeindex > this%netcdf_timestep_size) then
        this%last_stored_timeindex = this%netcdf_timestep_size
      end if
      call this%filehandle%read_netcdf_timesteps(this%first_stored_timeindex, this%last_stored_timeindex, values)

      call assert(allocated(values), __LINE__)
      if( allocated(this%stored_values) ) then
        if(all(shape(values) /= shape(this%stored_values))) then
          deallocate(this%stored_values)
        end if
      end if
      if(.not. allocated(this%stored_values)) then
        allocate(this%stored_values(size(values, DIM=1),size(values, DIM=2),size(values, DIM=3)))
      end if
      this%stored_values = values
      call assert(allocated(this%stored_values), __LINE__)    
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
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine
end module
