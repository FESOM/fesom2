module forcing_provider_module
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
    character(:), allocatable :: varname
    integer fileyear
    integer first_stored_timeindex, last_stored_timeindex
    real(4), allocatable :: stored_values(:,:,:)
    integer netcdf_timestep_size
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
    real(8), intent(out) :: forcingdata(:,:)
    ! EO args
    type(forcing_reader_type), allocatable :: tmparr(:)
    character(:), allocatable :: basepath
    integer i
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

    if( len(this%all_readers(varindex)%varname) == 0 ) then ! reader has never been initialized
      basepath=basepath_from_path(filepath, fileyear)
      
      this%all_readers(varindex)%basepath = basepath
      this%all_readers(varindex)%varname = varname
      this%all_readers(varindex)%fileyear = fileyear
      this%all_readers(varindex)%first_stored_timeindex = 0
      this%all_readers(varindex)%last_stored_timeindex = 0
      this%all_readers(varindex)%netcdf_timestep_size = read_netcdf_timestep_size(filepath, varname)
    end if
    
!    assert(time_index >= reader%first_stored_timeindex) ! we do not go back in time
    if(this%all_readers(varindex)%last_stored_timeindex < time_index) then
      
      this%all_readers(varindex)%first_stored_timeindex = time_index
      this%all_readers(varindex)%last_stored_timeindex = time_index+PREFETCH_SIZE
      if(this%all_readers(varindex)%last_stored_timeindex > this%all_readers(varindex)%netcdf_timestep_size) then
        this%all_readers(varindex)%last_stored_timeindex = this%all_readers(varindex)%netcdf_timestep_size
      end if
      call read_netcdf_timesteps(filepath, varname, this%all_readers(varindex)%first_stored_timeindex, this%all_readers(varindex)%last_stored_timeindex, values)

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
  
  
  subroutine read_netcdf_timesteps(filepath, varname, timeindex_first, timeindex_last, values)
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    integer, intent(in) :: timeindex_first, timeindex_last
    real(4), allocatable, intent(out) :: values(:,:,:)
    ! EO args
    integer, parameter :: timedim_index = 3
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer status
    integer fileid
    integer varid
    integer dim_size
    integer, allocatable, dimension(:) :: dimids
    integer, allocatable, dimension(:) :: dim_sizes
    integer i
    integer, allocatable, dimension(:) :: starts, sizes
   
    ! assume netcdf variable like: float q(time, lat, lon)
    status = nf_open(filepath,NF_NOWRITE,fileid)
    status = nf_inq_varid(fileid, varname, varid)
    status = nf_inq_varndims(fileid, varid, dim_size)
    allocate(dimids(dim_size))
    status = nf_inq_vardimid(fileid, varid, dimids)
    
    allocate(dim_sizes(dim_size))
    do i=1, dim_size
      status = nf_inq_dimlen(fileid, dimids(i), dim_sizes(i))
    end do
    
    ! todo: check if variable datatype is single precision (f77 real)
    
    allocate(starts(dim_size))
    allocate(sizes(dim_size))
    ! assert timedim_index == 3 && dim_size == 3
    call assert(timeindex_first <= timeindex_last, __LINE__)
    call assert(timeindex_last <= dim_sizes(timedim_index), __LINE__)
    
    ! todo: make this work if we have more than 3 dimensions and also if timedim_index != 3
    ! todo: check if values is already allocated
    starts(1) = 1
    starts(2) = 1
    starts(3) = timeindex_first
    sizes(1) = dim_sizes(1)
    sizes(2) = dim_sizes(2)
    sizes(3) = timeindex_last-timeindex_first+1
    allocate(values(sizes(1),sizes(2),sizes(3)))
    
    status = nf_get_vara_real(fileid, varid, starts, sizes, values)    
    status = nf_close(fileid)
    ! manual deallocation not required, the runtime environment will deallocate the arrays when out of scope
  end subroutine
  
  
  function read_netcdf_timestep_size(filepath, varname) result(r)
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    integer r
    ! EO args
    integer, parameter :: timedim_index = 3
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer status
    integer fileid
    integer varid
    integer dim_size
    integer, allocatable, dimension(:) :: dimids
   
    ! assume netcdf variable like: float q(time, lat, lon)
    status = nf_open(filepath,NF_NOWRITE,fileid)
    status = nf_inq_varid(fileid, varname, varid)
    status = nf_inq_varndims(fileid, varid, dim_size)
    allocate(dimids(dim_size))
    status = nf_inq_vardimid(fileid, varid, dimids)
    
    status = nf_inq_dimlen(fileid, dimids(timedim_index), r)
           
    status = nf_close(fileid)
  end function
  
  
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
