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
  

end module
