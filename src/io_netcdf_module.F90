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
    integer, allocatable :: varshape(:)
    contains
    procedure, public :: initialize
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
