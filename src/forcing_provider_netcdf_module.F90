module forcing_provider_netcdf_module
  implicit none
  public netcdf_reader_handle
  private

  integer, parameter :: TIMEDIM_INDEX = 3

  type netcdf_reader_handle
    private
    character(:), allocatable :: filepath
    character(:), allocatable :: varname
    integer fileid
    integer varid
    integer, allocatable :: varshape(:)
    contains
    procedure, public :: initialize, finalize, read_netcdf_timesteps, read_netcdf_timestep_2d, timestep_size
    procedure open_netcdf_variable
  end type
    

  contains
  
  
  subroutine initialize(this, filepath, varname)
    class(netcdf_reader_handle), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    ! EO args
    this%filepath = filepath
    this%varname = varname
    
    ! assert varshape is not allocated, i.e. init has not been called
    call assert(.not. allocated(this%varshape), __LINE__)
    
    call this%open_netcdf_variable()
  end subroutine
  
  
  function timestep_size(this) result(t)
    class(netcdf_reader_handle), intent(in) :: this
    integer t
    t = this%varshape(TIMEDIM_INDEX)
  end function


  subroutine read_netcdf_timesteps(this, timeindex_first, timeindex_last, values)
    class(netcdf_reader_handle), intent(in) :: this
    integer, intent(in) :: timeindex_first
    integer, intent(in) :: timeindex_last
    real(4), allocatable, intent(inout) :: values(:,:,:) ! must be inout or the allocation is screwed
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer, allocatable, dimension(:) :: starts, sizes

    call assert(allocated(this%varshape), __LINE__)
   
    ! assume netcdf variable like: float q(time, lat, lon)    
    ! todo: check if variable datatype is single precision (f77 real)
    
    allocate(starts(size(this%varshape)))
    allocate(sizes(size(this%varshape)))
    ! assert TIMEDIM_INDEX == 3 && var_dim_size == 3
    call assert(0 < timeindex_first, __LINE__)
    call assert(timeindex_first <= timeindex_last, __LINE__)
    call assert(timeindex_last <= this%varshape(timedim_index), __LINE__)
    
    ! todo: make this work if we have more than 3 dimensions and also if TIMEDIM_INDEX != 3
    starts(1) = 1
    starts(2) = 1
    starts(3) = timeindex_first
    sizes(1) = this%varshape(1)
    sizes(2) = this%varshape(2)
    sizes(3) = timeindex_last-timeindex_first+1
    if(allocated(values)) then
      if(all(sizes /= shape(values))) deallocate(values)
    end if
    if(.not. allocated(values)) allocate(values(sizes(1),sizes(2),sizes(3)))
    call assert(allocated(values), __LINE__)
    
    call assert_nc(nf_get_vara_real(this%fileid, this%varid, starts, sizes, values), __LINE__)  
  end subroutine


  subroutine read_netcdf_timestep_2d(this, timeindex, values)
    class(netcdf_reader_handle), intent(in) :: this
    integer, intent(in) :: timeindex
    real(4), allocatable, intent(inout) :: values(:,:) ! must be inout or the allocation is screwed
    ! EO args
    integer, parameter :: TIMEDIM_INDEX = 3
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer, allocatable, dimension(:) :: starts, sizes

    call assert(allocated(this%varshape), __LINE__)
   
    ! assume netcdf variable like: float q(time, lat, lon)    
    ! todo: check if variable datatype is single precision (f77 real)
    
    allocate(starts(size(this%varshape)))
    allocate(sizes(size(this%varshape)))
    ! assert TIMEDIM_INDEX == 3 && var_dim_size == 3
    call assert(0 < timeindex, __LINE__)
    call assert(timeindex <= this%varshape(timedim_index), __LINE__)
    
    ! todo: make this work if we have more than 3 dimensions and also if TIMEDIM_INDEX != 3
    starts(1) = 1
    starts(2) = 1
    starts(3) = timeindex
    sizes(1) = this%varshape(1)
    sizes(2) = this%varshape(2)
    sizes(3) = 1 !timeindex_last-timeindex_first+1
    if(allocated(values)) then
      if(all(sizes /= shape(values))) deallocate(values)
    end if
    if(.not. allocated(values)) allocate(values(sizes(1),sizes(2)))
    call assert(allocated(values), __LINE__)
    
    call assert_nc(nf_get_vara_real(this%fileid, this%varid, starts, sizes, values), __LINE__)  
  end subroutine
    

  subroutine open_netcdf_variable(this) !, fileid, varid, varshape)
    class(netcdf_reader_handle), intent(inout) :: this
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    integer var_dim_size
    integer, allocatable, dimension(:) :: dimids
    integer i
   
    ! assume netcdf variable like: float q(time, lat, lon)
    call assert_nc( nf_open(this%filepath, NF_NOWRITE, this%fileid) , __LINE__)
    call assert_nc( nf_inq_varid(this%fileid, this%varname, this%varid) , __LINE__)
    call assert_nc( nf_inq_varndims(this%fileid, this%varid, var_dim_size) , __LINE__)
    allocate(dimids(var_dim_size))
    call assert_nc( nf_inq_vardimid(this%fileid, this%varid, dimids) , __LINE__)

    allocate(this%varshape(var_dim_size))
    do i=1, var_dim_size
      call assert_nc( nf_inq_dimlen(this%fileid, dimids(i), this%varshape(i)) , __LINE__)
    end do
  end subroutine  


  ! do not implicitly close the file (e.g. upon deallocation via destructor), as we might have a copy of this object with access to the same fileid
  subroutine finalize(this)
    class(netcdf_reader_handle), intent(inout) :: this
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    if(allocated(this%varshape)) then
      call assert_nc( nf_close(this%fileid) , __LINE__)
    end if
  end subroutine


  subroutine assert_nc(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    if(status /= NF_NOERR) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf_strerror(status))
      stop 1
    endif   
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
