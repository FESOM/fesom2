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
    procedure, public :: initialize, finalize, number_of_timesteps
    procedure open_netcdf_variable
  end type


  contains


  subroutine initialize(this, filepath, varname)
    use netcdf
    class(netcdf_variable_handle), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    character(len=*), intent(in) :: varname
    ! EO args
    integer mode
    
    this%filepath = filepath
    this%varname = varname
    
    ! assert varshape is not allocated, i.e. initialize has not been called
    call assert(.not. allocated(this%varshape), __LINE__)
    call this%open_netcdf_variable(NF90_NOWRITE)

    ! assume the last dimension for this variable is the time dimension (i.e. first in ncdump)
    call assert(size(this%varshape) > 0, __LINE__)
    this%timedim_index = size(this%varshape)
  end subroutine


  subroutine open_netcdf_variable(this, mode)
    use netcdf
    class(netcdf_variable_handle), intent(inout) :: this
    integer, intent(in) :: mode
    ! EO args
    integer var_dim_size
    integer, allocatable, dimension(:) :: dimids
    integer i
   
    call assert_nc( nf90_open(this%filepath, mode, this%fileid) , __LINE__)
    call assert_nc( nf90_inq_varid(this%fileid, this%varname, this%varid) , __LINE__)
    call assert_nc( nf90_inquire_variable(this%fileid, this%varid, ndims=var_dim_size) , __LINE__)
    allocate(dimids(var_dim_size))
    call assert_nc( nf90_inquire_variable(this%fileid, this%varid, dimids=dimids) , __LINE__)

    allocate(this%varshape(var_dim_size))
    do i=1, var_dim_size
      call assert_nc( nf90_inquire_dimension(this%fileid, dimids(i), len=this%varshape(i)) , __LINE__)
    end do
  end subroutine  


  subroutine finalize(this)
    ! do not implicitly close the file (e.g. upon deallocation via destructor), as we might have a copy of this object with access to the same fileid
    use netcdf
    class(netcdf_variable_handle), intent(inout) :: this
    ! EO args
    if(allocated(this%varshape)) then
      call assert_nc( nf90_close(this%fileid) , __LINE__)
    end if
  end subroutine

  
  function number_of_timesteps(this) result(t)
    class(netcdf_variable_handle), intent(in) :: this
    integer t
    ! EO args
    t = this%varshape(this%timedim_index)
  end function


  subroutine assert_nc(status, line)
    use netcdf
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    if(status /= nf90_noerr) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf90_strerror(status))
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
