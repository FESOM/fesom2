module io_fesom_file_module
  implicit none
  public fesom_file_type
  private


  type fesom_file_type
    private
    type(dim_type), allocatable :: dims(:)
  contains
    procedure, public :: initialize, add_dim
  end type
  
  
  type dim_type
    character(:), allocatable :: name
    integer len
  end type


contains


  subroutine initialize(this)
    class(fesom_file_type), intent(inout) :: this
  end subroutine


  function add_dim(this, name, len) result(dimindex)
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: len
    integer dimindex
    ! EO parameters
    type(dim_type), allocatable :: tmparr(:)
    
    if( .not. allocated(this%dims)) then
      allocate(this%dims(1))
    else
      allocate( tmparr(size(this%dims)+1) )
      tmparr(1:size(this%dims)) = this%dims
      deallocate(this%dims)
      call move_alloc(tmparr, this%dims)
    end if
    
    dimindex = size(this%dims)
    this%dims(dimindex) = dim_type(name=name, len=len)
  end function
  



  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO parameters
    if(.not. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine


  subroutine assert_nc(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO parameters
    include "netcdf.inc"
    if(status /= nf_noerr) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf_strerror(status))
      stop 1
    endif   
  end subroutine

end module
