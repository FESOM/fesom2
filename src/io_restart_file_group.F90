module restart_file_group_module
  use io_fesom_file_module
  implicit none
  private

  
  type, extends(fesom_file_type) :: restart_file_type
    private
    integer iter_varindex
    character(:), allocatable :: varname
    character(:), allocatable :: path
  end type


  type restart_file_group
    private
    type(restart_file_type) files(20); integer :: nfiles = 0 ! todo: allow dynamically allocated size without messing with shallow copied pointers
  contains
  end type

contains



end module
