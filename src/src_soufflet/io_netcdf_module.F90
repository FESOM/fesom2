module io_netcdf_module
  implicit none
  public netcdf_variable_handle
  private

  type netcdf_variable_handle
    private
    character(:), allocatable :: filepath
    character(:), allocatable :: varname
    integer fileid
    integer varid
    integer timedim_index
    integer, allocatable :: varshape(:)
    contains
    procedure, public :: initialize, finalize, number_of_timesteps, number_of_dims, dimsize_at
    generic, public :: read_values => read_values_r4,read_values_r8
    procedure, private :: open_netcdf_variable, read_values_r4, read_values_r8
  end type


  contains


  subroutine initialize(this, filepath, varname)
    class(netcdf_variable_handle), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    ! EO args
    include "netcdf.inc"
    integer mode
    
    this%filepath = filepath
    this%varname = varname
    
    ! assert varshape is not allocated, i.e. initialize has not been called
    call assert(.not. allocated(this%varshape), __LINE__)
    call this%open_netcdf_variable(NF_NOWRITE)

    ! assume the last dimension for this variable is the time dimension (i.e. first in ncdump)
    call assert(size(this%varshape) > 0, __LINE__)
    this%timedim_index = size(this%varshape)
  end subroutine


  subroutine open_netcdf_variable(this, mode)
    class(netcdf_variable_handle), intent(inout) :: this
    integer, intent(in) :: mode
    ! EO args
    include "netcdf.inc"
    integer var_dim_size
    integer, allocatable, dimension(:) :: dimids
    integer i
   
    call assert_nc( nf_open(this%filepath, mode, this%fileid) , __LINE__)
    call assert_nc( nf_inq_varid(this%fileid, this%varname, this%varid) , __LINE__)
    call assert_nc( nf_inq_varndims(this%fileid, this%varid, var_dim_size) , __LINE__)
    allocate(dimids(var_dim_size))
    call assert_nc( nf_inq_vardimid(this%fileid, this%varid, dimids) , __LINE__)

    allocate(this%varshape(var_dim_size))
    do i=1, var_dim_size
      call assert_nc( nf_inq_dimlen(this%fileid, dimids(i), this%varshape(i)) , __LINE__)
    end do
  end subroutine  


  subroutine finalize(this)
    ! do not implicitly close the file (e.g. upon deallocation via destructor), as we might have a copy of this object with access to the same fileid
    class(netcdf_variable_handle), intent(inout) :: this
    ! EO args
    include "netcdf.inc"
    if(allocated(this%varshape)) then
      call assert_nc( nf_close(this%fileid) , __LINE__)
    end if
  end subroutine

  
  function number_of_timesteps(this) result(t)
    class(netcdf_variable_handle), intent(in) :: this
    integer t
    ! EO args
    t = this%varshape(this%timedim_index)
  end function

  
  function number_of_dims(this) result(d)
    class(netcdf_variable_handle), intent(in) :: this
    integer d
    ! EO args
    d = size(this%varshape)
  end function

  
  function dimsize_at(this,index) result(s)
    class(netcdf_variable_handle), intent(in) :: this
    integer, intent(in) :: index
    integer s
    ! EO args
    call assert(index <= size(this%varshape), __LINE__)
    s = this%varshape(index)
  end function


  subroutine read_values_r8(this, timeindex, values)
    use io_netcdf_nf_interface
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_variable_handle), intent(in) :: this
    integer, intent(in) :: timeindex
    real(8), intent(inout), target :: values(:) ! must be inout or the allocation might be screwed
    ! EO args
    real(8), pointer :: values_ptr(:)
    integer, allocatable, dimension(:) :: starts, sizes

    call read_values_preflight(this, timeindex, starts, sizes)

    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_x(this%fileid, this%varid, starts, sizes, values_ptr), __LINE__)
  end subroutine
  
  
  subroutine read_values_r4(this, timeindex, values)
    use io_netcdf_nf_interface
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_variable_handle), intent(in) :: this
    integer, intent(in) :: timeindex
    real(4), intent(inout), target :: values(:) ! must be inout or the allocation might be screwed
    ! EO args
    real(4), pointer :: values_ptr(:)
    integer, allocatable, dimension(:) :: starts, sizes
    
    call read_values_preflight(this, timeindex, starts, sizes)

    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_x(this%fileid, this%varid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  subroutine read_values_preflight(this, timeindex, starts, sizes)
    class(netcdf_variable_handle), intent(in) :: this
    integer, intent(in) :: timeindex
    ! EO args

    integer, allocatable, dimension(:) :: starts, sizes

    call assert(allocated(this%varshape), __LINE__)
       
    allocate(starts(size(this%varshape)))
    allocate(sizes(size(this%varshape)))

    call assert(0 < timeindex, __LINE__)
    call assert(timeindex <= this%number_of_timesteps(), __LINE__)
    
    starts = 1
    sizes = this%varshape
    starts(this%timedim_index) = timeindex
    sizes(this%timedim_index) = 1 !timeindex_last-timeindex_first+1
  end subroutine


  subroutine assert_nc(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    include "netcdf.inc"
    if(status /= nf_noerr) then
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
