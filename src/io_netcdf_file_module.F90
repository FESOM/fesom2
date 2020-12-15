module io_netcdf_file_module
  implicit none
  public fesom_file_type
  private


  type fesom_file_type
    private
    type(dim_type), allocatable :: dims(:)
    type(var_type), allocatable :: vars(:)

    character(:), allocatable :: filepath
    integer ncid
  contains
    procedure, public :: initialize, add_dim, add_dim_unlimited, add_var, add_var_att, open_read, close_file
    generic, public :: read_var => read_var_r4, read_var_r8
    procedure, private :: read_var_r4, read_var_r8
  end type
  
  
  type dim_type
    character(:), allocatable :: name
    integer len
    
    integer ncid
  end type


  type var_type ! todo: use variable type from io_netcdf_module here
    character(:), allocatable :: name
    integer, allocatable :: dim_indices(:)
    type(att_type), allocatable :: atts(:)
    
    integer ncid
  end type
  
  
  type att_type
    character(:), allocatable :: name
    character(:), allocatable :: text
    ! todo: make this work for other data types like int
  end type


contains


  subroutine initialize(this)
    class(fesom_file_type), intent(inout) :: this
    
    this%filepath = ""
  end subroutine

  
  function add_dim_unlimited(this, name) result(dimindex)
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer dimindex
    ! EO parameters
    include "netcdf.inc"

    dimindex = this%add_dim(name, nf_unlimited)
  end function


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
  

  ! the sizes of the dims define the global shape of the var
  function add_var(this, name, dim_indices) result(varindex)
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer varindex
    ! EO parameters
    type(var_type), allocatable :: tmparr(:)
    
    if( .not. allocated(this%vars)) then
      allocate(this%vars(1))
    else
      allocate( tmparr(size(this%vars)+1) )
      tmparr(1:size(this%vars)) = this%vars
      deallocate(this%vars)
      call move_alloc(tmparr, this%vars)
    end if
    
    varindex = size(this%vars)
    this%vars(varindex) = var_type(name, dim_indices, ncid=-1)
  end function


  subroutine add_var_att(this, varindex, att_name, att_text)
    class(fesom_file_type), intent(inout) :: this
    integer, intent(in) :: varindex
    character(len=*), intent(in) :: att_name
    character(len=*), intent(in) :: att_text
    ! EO parameters
    type(att_type), allocatable :: tmparr(:)
    
    if( .not. allocated(this%vars(varindex)%atts)) then
      allocate(this%vars(varindex)%atts(1))
    else
      allocate( tmparr(size(this%vars(varindex)%atts)+1) )
      tmparr(1:size(this%vars(varindex)%atts)) = this%vars(varindex)%atts
      deallocate(this%vars(varindex)%atts)
      call move_alloc(tmparr, this%vars(varindex)%atts)
    end if
    
    this%vars(varindex)%atts( size(this%vars(varindex)%atts) ) = att_type(name=att_name, text=att_text)
  end subroutine


  subroutine open_read(this, filepath)
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    ! EO parameters
    include "netcdf.inc"
    integer i, ii
    integer actual_len
    integer actual_dimcount
    integer, allocatable :: actual_dimids(:)
    integer exp_dimid, act_dimid
    integer mode

    mode = nf_nowrite
    this%filepath = filepath
    
    call assert_nc( nf_open(this%filepath, mode, this%ncid) , __LINE__)
    
    ! attach our dims and vars to their counterparts in the file
    do i=1, size(this%dims)
      call assert_nc( nf_inq_dimid(this%ncid, this%dims(i)%name, this%dims(i)%ncid) , __LINE__)
      call assert_nc( nf_inq_dimlen(this%ncid, this%dims(i)%ncid, actual_len) , __LINE__)
      if(this%dims(i)%len .ne. nf_unlimited) call assert(this%dims(i)%len == actual_len, __LINE__)
    end do
    do i=1, size(this%vars)
      call assert_nc( nf_inq_varid(this%ncid, this%vars(i)%name, this%vars(i)%ncid) , __LINE__)
      ! see if this var has the expected dims
      call assert_nc( nf_inq_varndims(this%ncid, this%vars(i)%ncid, actual_dimcount) , __LINE__)
      call assert(size(this%vars(i)%dim_indices) == actual_dimcount, __LINE__)
      allocate(actual_dimids(actual_dimcount))
      call assert_nc( nf_inq_vardimid(this%ncid, this%vars(i)%ncid, actual_dimids) , __LINE__)
      do ii=1, actual_dimcount
        exp_dimid = this%dims( this%vars(i)%dim_indices(ii) )%ncid
        call assert(exp_dimid == actual_dimids(ii), __LINE__)
      end do
    end do
  end subroutine


  ! values array is not required to have the same shape as the variable but must fit the product of all items of the sizes array
  ! this way we can retrieve e.g. data from a 3D variable to a 2D array with one size set to 1 (e.g. to get a single timestep)
  ! starts and sizes must have the same rank as the variable has dimensions
  subroutine read_var_r8(this, varindex, starts, sizes, values)
    use io_netcdf_nf_interface
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(8), intent(inout), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    real(8), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%dims), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_x(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  ! see read_var_r8 for usage comment
  subroutine read_var_r4(this, varindex, starts, sizes, values)
    use io_netcdf_nf_interface
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(4), intent(inout), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    real(4), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%dims), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_x(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
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
