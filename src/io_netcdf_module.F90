module io_netcdf_module
  implicit none
  private

  type netcdf_variable_handle
    private
    character(:), allocatable :: filepath
    character(:), allocatable :: varname
    integer fileid
    integer varid
    integer, allocatable :: varshape(:)
    contains
  end type

  contains


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
