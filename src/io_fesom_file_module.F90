module io_fesom_file_module
  implicit none
  public fesom_file
  private


  type fesom_file
  contains
    procedure, public :: initialize
  end type


contains


  subroutine initialize(this)
    class(fesom_file), intent(inout) :: this
  end subroutine  
    

end module
