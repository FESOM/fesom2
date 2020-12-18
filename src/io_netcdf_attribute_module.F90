module io_netcdf_attribute_module
  implicit none
  public att_type, att_type_text, att_type_int
  private
  
  
  type, abstract :: att_type
    character(:), allocatable :: name
  contains
    procedure(define_in_var), deferred :: define_in_var
  end type
  
  
  interface
    subroutine define_in_var(this, fileid, varid)
      import att_type
      class(att_type), intent(inout) :: this
      integer, intent(in) :: fileid
      integer, intent(in) :: varid
    end subroutine
  end interface


  type, extends(att_type) :: att_type_text
    character(:), allocatable :: text
  contains
    procedure :: define_in_var => define_in_var_text
  end type


  type, extends(att_type) :: att_type_int
    integer :: val
  contains
    procedure :: define_in_var => define_in_var_int
  end type


contains


  subroutine define_in_var_text(this, fileid, varid)
    class(att_type_text), intent(inout) :: this
    integer, intent(in) :: fileid
    integer, intent(in) :: varid
    ! EO parameters
    include "netcdf.inc"

    call assert_nc( nf_put_att_text(fileid, varid, this%name, len(this%text), this%text) , __LINE__)
  end subroutine


  subroutine define_in_var_int(this, fileid, varid)
    class(att_type_int), intent(inout) :: this
    integer, intent(in) :: fileid
    integer, intent(in) :: varid
    ! EO parameters
    include "netcdf.inc"

    call assert_nc( nf_put_att_int(fileid, varid, this%name, nf_int, 1, this%val) , __LINE__)
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
