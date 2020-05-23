module forcing_provider_module
  use forcing_provider_netcdf_module
  implicit none
  public forcing_provider
  private
  
  type forcing_provider_type
    private
    type(forcing_reader_type), allocatable :: all_readers(:)
    contains
    procedure, public :: get_forcingdata
  end type
  type(forcing_provider_type), save :: forcing_provider ! handle to singleton instance
  
  type, private :: forcing_reader_type
    character(:), allocatable :: basepath
    integer fileyear
    integer first_stored_timeindex, last_stored_timeindex
    real(4), allocatable :: stored_values(:,:,:)
    integer netcdf_timestep_size
    type(netcdf_reader_handle) filehandle
  end type

  character(len=*), parameter :: FILENAMESUFFIX = ".nc"
  integer, parameter :: PREFETCH_SIZE = 10
  
  
  contains

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
    integer reader_time_index
    real(4), allocatable :: values(:,:,:)
    
    ! init our all_readers array if not already done
    if(.not. allocated(this%all_readers)) then
      allocate(this%all_readers(varindex))
    end if
        
    if(size(this%all_readers) < varindex) then
      allocate( tmparr(varindex) )
      tmparr(1:size(this%all_readers)) = this%all_readers
      deallocate(this%all_readers)
      call move_alloc(tmparr, this%all_readers)      
    end if
    
    if( len(this%all_readers(varindex)%basepath) == 0 ) then ! reader has never been initialized ! todo: change this as it is probably compiler dependent
      basepath = basepath_from_path(filepath, fileyear)
      
      this%all_readers(varindex)%basepath = basepath
      this%all_readers(varindex)%fileyear = fileyear
      this%all_readers(varindex)%first_stored_timeindex = -1
      this%all_readers(varindex)%last_stored_timeindex = -1
      call this%all_readers(varindex)%filehandle%initialize(filepath, varname) ! finalize() to close the file
      this%all_readers(varindex)%netcdf_timestep_size = this%all_readers(varindex)%filehandle%timestep_size()
    else if(fileyear /= this%all_readers(varindex)%fileyear) then
      print *,"can not change years", __LINE__, __FILE__
      stop 1
    end if
    
!    assert(time_index >= reader%first_stored_timeindex) ! we do not go back in time
    if(this%all_readers(varindex)%last_stored_timeindex < time_index) then
      
      this%all_readers(varindex)%first_stored_timeindex = time_index
      this%all_readers(varindex)%last_stored_timeindex = time_index+PREFETCH_SIZE
      if(this%all_readers(varindex)%last_stored_timeindex > this%all_readers(varindex)%netcdf_timestep_size) then
        this%all_readers(varindex)%last_stored_timeindex = this%all_readers(varindex)%netcdf_timestep_size
      end if
      call this%all_readers(varindex)%filehandle%read_netcdf_timesteps(this%all_readers(varindex)%first_stored_timeindex, this%all_readers(varindex)%last_stored_timeindex, values)

      call assert(allocated(values), __LINE__)
      if( allocated(this%all_readers(varindex)%stored_values) ) then
        if(all(shape(values) /= shape(this%all_readers(varindex)%stored_values))) then
          deallocate(this%all_readers(varindex)%stored_values)
        end if
      end if
      if(.not. allocated(this%all_readers(varindex)%stored_values)) then
        allocate(this%all_readers(varindex)%stored_values(size(values, DIM=1),size(values, DIM=2),size(values, DIM=3)))
      end if
      this%all_readers(varindex)%stored_values = values
      call assert(allocated(this%all_readers(varindex)%stored_values), __LINE__)    
    end if
  
    ! check if the outgoing array has the same shape as our data
    reader_time_index = time_index-this%all_readers(varindex)%first_stored_timeindex+1
    call assert(allocated(this%all_readers(varindex)%stored_values), __LINE__)
    call assert( all( shape(forcingdata)==shape(this%all_readers(varindex)%stored_values(:,:,reader_time_index)) ), __LINE__ )
    forcingdata = this%all_readers(varindex)%stored_values(:,:,reader_time_index)    
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
