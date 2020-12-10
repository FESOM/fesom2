module io_fesom_file_module
  implicit none
  public fesom_file_type
  private


  type fesom_file_type
    private
    type(dim_type), allocatable :: dims(:)

    character(:), allocatable :: filepath
    integer mode
    integer ncid
  contains
    procedure, public :: initialize, add_dim, open_readmode, close_file
  end type
  
  
  type dim_type
    character(:), allocatable :: name
    integer len
    
    integer ncid
  end type


contains


  subroutine initialize(this)
    class(fesom_file_type), intent(inout) :: this
    
    this%filepath = ""
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
    this%dims(dimindex) = dim_type(name=name, len=len, ncid=-1)
  end function
  

  subroutine open_readmode(this, filepath)
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    ! EO parameters
    include "netcdf.inc"
    integer i
    integer actual_len

    this%mode = nf_nowrite
    this%filepath = filepath
    
    call assert_nc( nf_open(this%filepath, this%mode, this%ncid) , __LINE__)
    
    ! attach our dims to their counterparts in the file
    do i=1, size(this%dims)
      call assert_nc( nf_inq_dimid(this%ncid, this%dims(i)%name, this%dims(i)%ncid) , __LINE__)
      call assert_nc( nf_inq_dimlen(this%ncid, this%dims(i)%ncid, actual_len) , __LINE__)
      call assert(this%dims(i)%len == actual_len, __LINE__)
    end do
  end subroutine


  subroutine close_file(this)
    ! do not implicitly close the file (e.g. upon deallocation via destructor), as we might have a copy of this object with access to the same ncid
    class(fesom_file_type), intent(inout) :: this
    ! EO parameters
    include "netcdf.inc"
    call assert_nc( nf_close(this%ncid) , __LINE__)
  end subroutine


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
