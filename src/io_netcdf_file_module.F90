module io_netcdf_file_module
  use io_netcdf_attribute_module
  implicit none
  public netcdf_file_type
  private


  type netcdf_file_type
    private
    type(dim_type), allocatable :: dims(:)
    type(var_type), allocatable :: vars(:)
    type(att_type_wrapper), allocatable :: gatts(:)

    character(:), allocatable :: filepath
    integer ncid
  contains
    procedure, public :: initialize, add_dim, add_dim_unlimited, add_var_double, add_var_real, add_var_int, open_read, flush_file, close_file, open_write_create, open_write_append
    procedure, public :: is_attached, read_var_shape
    procedure, public :: ndims
    generic, public :: read_var => read_var_r4, read_var_r8, read_var_integer
    generic, public :: write_var => write_var_r4, write_var_r8, write_var_integer
    generic, public :: read_var1 => read_var1_r4, read_var1_r8, read_var1_integer
    generic, public :: add_var_att => add_var_att_text, add_var_att_int
    generic, public :: add_global_att => add_global_att_text, add_global_att_int
    procedure, private :: read_var_r4, read_var_r8, read_var_integer, attach_dims_vars_to_file, add_var_x, write_var_r4, write_var_r8, write_var_integer, add_var_att_text, add_var_att_int
    procedure, private :: read_var1_r4, read_var1_r8, read_var1_integer
    procedure, private :: add_global_att_text, add_global_att_int
  end type
  
  
  type dim_type
    character(:), allocatable :: name
    integer len
    
    integer ncid
  end type


  type var_type ! todo: use variable type from io_netcdf_module here
    character(:), allocatable :: name
    integer, allocatable :: dim_indices(:)
    integer datatype
    type(att_type_wrapper), allocatable :: atts(:)
    
    integer ncid
  end type


  type att_type_wrapper ! work around Fortran not being able to have polymorphic types in the same array
    class(att_type), allocatable :: it
  end type
 
contains


  subroutine initialize(this)
    class(netcdf_file_type), intent(inout) :: this
    
    this%filepath = ""
    allocate(this%dims(0))
    allocate(this%vars(0))
    allocate(this%gatts(0))
  end subroutine

  
  function add_dim_unlimited(this, name) result(dimindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer dimindex
    ! EO parameters
    include "netcdf.inc"

    dimindex = this%add_dim(name, nf_unlimited)
  end function


  function add_dim(this, name, len) result(dimindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: len
    integer dimindex
    ! EO parameters
    type(dim_type), allocatable :: tmparr(:)
    
    ! assume the dims array is allocated
    allocate( tmparr(size(this%dims)+1) )
    tmparr(1:size(this%dims)) = this%dims
    deallocate(this%dims)
    call move_alloc(tmparr, this%dims)
    
    dimindex = size(this%dims)
    this%dims(dimindex) = dim_type(name=name, len=len, ncid=-1)
  end function


  ! return number of specified dimensions (which might be less dimensions than an attached file has)
  function ndims(this)
    class(netcdf_file_type), intent(inout) :: this
    integer ndims
    ! EO parameters
    
    ndims = size(this%dims)
  end function
  

  ! the sizes of the dims define the global shape of the var
  function add_var_double(this, name, dim_indices) result(varindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer varindex
    ! EO parameters
    include "netcdf.inc"
    
    varindex = this%add_var_x(name, dim_indices, nf_double)
  end function


  ! the sizes of the dims define the global shape of the var
  function add_var_real(this, name, dim_indices) result(varindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer varindex
    ! EO parameters
    include "netcdf.inc"
    
    varindex = this%add_var_x(name, dim_indices, nf_real)
  end function


  ! the sizes of the dims define the global shape of the var
  function add_var_int(this, name, dim_indices) result(varindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer varindex
    ! EO parameters
    include "netcdf.inc"
    
    varindex = this%add_var_x(name, dim_indices, nf_int)
  end function


  function add_var_x(this, name, dim_indices, netcdf_datatype) result(varindex)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer netcdf_datatype
    integer varindex
    ! EO parameters
    include "netcdf.inc"
    type(var_type), allocatable :: tmparr(:)
    type(att_type_wrapper), allocatable :: empty_atts(:)
    
    allocate(empty_atts(0)) ! if we use a fixed size array with size 0 there is a segfault at runtime when compiled with cray ftn

    ! assume the vars array is allocated
    allocate( tmparr(size(this%vars)+1) )
    tmparr(1:size(this%vars)) = this%vars
    deallocate(this%vars)
    call move_alloc(tmparr, this%vars)
    
    varindex = size(this%vars)
    this%vars(varindex) = var_type(name=name, dim_indices=dim_indices, datatype=netcdf_datatype, atts=empty_atts, ncid=-1)
  end function


  subroutine add_global_att_text(this, att_name, att_text)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    character(len=*), intent(in) :: att_text
    ! EO parameters
    type(att_type_wrapper), allocatable :: tmparr(:)
    
    allocate( tmparr(size(this%gatts)+1) )
    tmparr(1:size(this%gatts)) = this%gatts
    deallocate(this%gatts)
    call move_alloc(tmparr, this%gatts)
    
    this%gatts( size(this%gatts) )%it = att_type_text(name=att_name, text=att_text)
  end subroutine


  subroutine add_global_att_int(this, att_name, att_val)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer, intent(in) :: att_val
    ! EO parameters
    type(att_type_wrapper), allocatable :: tmparr(:)
    
    allocate( tmparr(size(this%gatts)+1) )
    tmparr(1:size(this%gatts)) = this%gatts
    deallocate(this%gatts)
    call move_alloc(tmparr, this%gatts)
    
    this%gatts( size(this%gatts) )%it = att_type_int(name=att_name, val=att_val)
  end subroutine


  subroutine add_var_att_text(this, varindex, att_name, att_text)
    class(netcdf_file_type), intent(inout) :: this
    integer, intent(in) :: varindex
    character(len=*), intent(in) :: att_name
    character(len=*), intent(in) :: att_text
    ! EO parameters
    type(att_type_wrapper), allocatable :: tmparr(:)

    allocate( tmparr(size(this%vars(varindex)%atts)+1) )
    tmparr(1:size(this%vars(varindex)%atts)) = this%vars(varindex)%atts
    deallocate(this%vars(varindex)%atts)
    call move_alloc(tmparr, this%vars(varindex)%atts)
    
    this%vars(varindex)%atts( size(this%vars(varindex)%atts) )%it = att_type_text(name=att_name, text=att_text)
  end subroutine


  subroutine add_var_att_int(this, varindex, att_name, att_val)
    class(netcdf_file_type), intent(inout) :: this
    integer, intent(in) :: varindex
    character(len=*), intent(in) :: att_name
    integer, intent(in) :: att_val
    ! EO parameters
    type(att_type_wrapper), allocatable :: tmparr(:)
    
    allocate( tmparr(size(this%vars(varindex)%atts)+1) )
    tmparr(1:size(this%vars(varindex)%atts)) = this%vars(varindex)%atts
    deallocate(this%vars(varindex)%atts)
    call move_alloc(tmparr, this%vars(varindex)%atts)
    
    this%vars(varindex)%atts( size(this%vars(varindex)%atts) )%it = att_type_int(name=att_name, val=att_val)
  end subroutine
  
  
  function is_attached(this) result(x)
    class(netcdf_file_type), intent(in) :: this
    logical x
    ! EO parameters

    x = (this%filepath .ne. "")
  end function


  subroutine open_read(this, filepath)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    ! EO parameters
    include "netcdf.inc"
    integer mode

    mode = nf_nowrite
    this%filepath = filepath
    
    call assert_nc( nf_open(this%filepath, mode, this%ncid) , __LINE__)
    
    ! attach our dims and vars to their counterparts in the file
    call this%attach_dims_vars_to_file()
  end subroutine


  ! return an array with the dimension sizes for all dimensions of the given variable
  ! this currently only makes sense for variables with unlimited dimensions,
  ! as all other dimensions must be known when adding the variable to the file specification, e.g before reading the file
  subroutine read_var_shape(this, varindex, varshape)
    class(netcdf_file_type), target, intent(in) :: this
    integer, intent(in) :: varindex
    integer, allocatable, intent(out) :: varshape(:)
    ! EO parameters
    include "netcdf.inc"
    type(var_type), pointer :: var
    integer var_ndims
    integer i
    
    var => this%vars(varindex)
    var_ndims = size(var%dim_indices)
    
    if(allocated(varshape)) deallocate(varshape)
    allocate(varshape(var_ndims))
         
    do i=1, var_ndims
      if(this%dims( var%dim_indices(i) )%len == nf_unlimited) then
        ! actually read from the file
        call assert_nc( nf_inq_dimlen(this%ncid, this%dims( var%dim_indices(i) )%ncid, varshape(i)) , __LINE__)
      else
        ! use the dim size which has been set without the file and is thus known anyway to the user
        varshape(i) = this%dims( var%dim_indices(i) )%len
      end if
    end do
  end subroutine


  ! values array is not required to have the same shape as the variable but must fit the product of all items of the sizes array
  ! this way we can retrieve e.g. data from a 3D variable to a 2D array with one size set to 1 (e.g. to get a single timestep)
  ! starts and sizes must have the same rank as the variable has dimensions
  subroutine read_var_r8(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(8), intent(inout), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    real(8), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_double(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  ! see read_var_r8 for usage comment
  subroutine read_var_r4(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(4), intent(inout), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    real(4), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_real(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  ! see read_var_r8 for usage comment
  subroutine read_var_integer(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    integer, intent(inout), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    integer, pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_get_vara_int(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  ! retrieve a single value specified via the indices array
  subroutine read_var1_r8(this, varindex, indices, value)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: indices
    real(8), intent(out) :: value
    ! EO parameters
    include "netcdf.inc"

    call assert(size(indices) == size(this%vars(varindex)%dim_indices), __LINE__)

    call assert_nc(nf_get_var1_double(this%ncid, this%vars(varindex)%ncid, indices, value), __LINE__)
  end subroutine


  ! see read_var1_r8 for usage comment
  subroutine read_var1_r4(this, varindex, indices, value)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: indices
    real(4), intent(out) :: value
    ! EO parameters
    include "netcdf.inc"

    call assert(size(indices) == size(this%vars(varindex)%dim_indices), __LINE__)

    call assert_nc(nf_get_var1_real(this%ncid, this%vars(varindex)%ncid, indices, value), __LINE__)
  end subroutine


  ! see read_var1_r8 for usage comment
  subroutine read_var1_integer(this, varindex, indices, value)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: indices
    integer, intent(out) :: value
    ! EO parameters
    include "netcdf.inc"

    call assert(size(indices) == size(this%vars(varindex)%dim_indices), __LINE__)

    call assert_nc(nf_get_var1_int(this%ncid, this%vars(varindex)%ncid, indices, value), __LINE__)
  end subroutine


  subroutine open_write_create(this, filepath)
    class(netcdf_file_type), target, intent(inout) :: this
    character(len=*), intent(in) :: filepath
    ! EO parameters
    include "netcdf.inc"
    integer cmode
    integer i, ii
    integer var_ndims
    integer, allocatable :: var_dimids(:)
    character(:), pointer :: att_name
    character(:), pointer :: att_text

    this%filepath = filepath

    cmode = ior(nf_noclobber, ior(nf_netcdf4, nf_classic_model))
    call assert_nc( nf_create(filepath, cmode, this%ncid) , __LINE__)
        
    ! create our dims in the file
    do i=1, size(this%dims)
      call assert_nc( nf_def_dim(this%ncid, this%dims(i)%name, this%dims(i)%len, this%dims(i)%ncid) , __LINE__)
    end do
    
    do i=1, size(this%gatts)
      call this%gatts(i)%it%define_in_var(this%ncid, nf_global)
    end do

    ! create our vars in the file
    do i=1, size(this%vars)
      var_ndims = size(this%vars(i)%dim_indices)
      if(allocated(var_dimids)) deallocate(var_dimids)
      allocate(var_dimids(var_ndims))
      do ii=1, var_ndims
        var_dimids(ii) = this%dims( this%vars(i)%dim_indices(ii) )%ncid
      end do
      call assert_nc( nf_def_var(this%ncid, this%vars(i)%name, this%vars(i)%datatype, var_ndims, var_dimids, this%vars(i)%ncid) , __LINE__)
      
      do ii=1, size(this%vars(i)%atts)
        call this%vars(i)%atts(ii)%it%define_in_var(this%ncid, this%vars(i)%ncid)
      end do
    end do

    call assert_nc( nf_enddef(this%ncid), __LINE__ )
  end subroutine


  ! open an existing file and prepare to write data to it
  subroutine open_write_append(this, filepath)
    class(netcdf_file_type), intent(inout) :: this
    character(len=*), intent(in) :: filepath
    ! EO parameters
    include "netcdf.inc"
    integer cmode

    this%filepath = filepath

    cmode = nf_write
    call assert_nc( nf_open(filepath, cmode, this%ncid) , __LINE__)

    ! make sure that all our dims and vars exist in this file and get hold of them
    call this%attach_dims_vars_to_file()
  end subroutine


  subroutine write_var_r8(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(8), intent(in), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    real(8), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_put_vara_double(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  subroutine write_var_r4(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    real(4), intent(in), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    real(4), pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_put_vara_real(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  subroutine write_var_integer(this, varindex, starts, sizes, values)
    use, intrinsic :: ISO_C_BINDING
    class(netcdf_file_type), intent(in) :: this
    integer, intent(in) :: varindex
    integer, dimension(:) :: starts, sizes
    integer, intent(in), target :: values(..) ! must be inout or the allocation might be screwed
    ! EO parameters
    include "netcdf.inc"
    integer, pointer :: values_ptr(:)

    call assert(size(sizes) == size(starts), __LINE__)
    call assert(size(starts) == size(this%vars(varindex)%dim_indices), __LINE__)
    call assert(product(sizes) == product(shape(values)), __LINE__)

    call c_f_pointer(c_loc(values), values_ptr, [product(shape(values))])
    call assert_nc(nf_put_vara_int(this%ncid, this%vars(varindex)%ncid, starts, sizes, values_ptr), __LINE__)
  end subroutine


  subroutine flush_file(this)
    class(netcdf_file_type), intent(inout) :: this
    ! EO parameters
    include "netcdf.inc"

    call assert_nc( nf_sync(this%ncid), __LINE__ ) ! flush the file to disk
  end subroutine


  subroutine close_file(this)
    ! do not implicitly close the file (e.g. upon deallocation via destructor), as we might have a copy of this object with access to the same ncid
    class(netcdf_file_type), intent(inout) :: this
    ! EO parameters
    include "netcdf.inc"
    call assert_nc( nf_close(this%ncid) , __LINE__)
    
    this%filepath = ""
  end subroutine


  ! connect our dims and vars to their counterparts in the NetCDF file, bail out if they do not match
  ! ignore any additional dims and vars the file might contain
  subroutine attach_dims_vars_to_file(this)
    class(netcdf_file_type), intent(inout) :: this
    ! EO parameters
    include "netcdf.inc"
    integer i, ii
    integer actual_len
    integer actual_dimcount
    integer, allocatable :: actual_dimids(:)
    integer exp_dimid, act_dimid
    integer actual_datatype

    do i=1, size(this%dims)
      call assert_nc( nf_inq_dimid(this%ncid, this%dims(i)%name, this%dims(i)%ncid) , __LINE__)
      call assert_nc( nf_inq_dimlen(this%ncid, this%dims(i)%ncid, actual_len) , __LINE__)
      if(this%dims(i)%len .ne. nf_unlimited) call assert(this%dims(i)%len == actual_len, __LINE__)
    end do
    do i=1, size(this%vars)
      call assert_nc( nf_inq_varid(this%ncid, this%vars(i)%name, this%vars(i)%ncid) , __LINE__)
      ! see if this var has the expected datatype
      call assert_nc( nf_inq_vartype(this%ncid, this%vars(i)%ncid, actual_datatype) , __LINE__)
      call assert(this%vars(i)%datatype == actual_datatype, __LINE__)
      ! see if this var has the expected dims
      call assert_nc( nf_inq_varndims(this%ncid, this%vars(i)%ncid, actual_dimcount) , __LINE__)
      call assert(size(this%vars(i)%dim_indices) == actual_dimcount, __LINE__)
      if(allocated(actual_dimids)) deallocate(actual_dimids)
      allocate(actual_dimids(actual_dimcount))
      call assert_nc( nf_inq_vardimid(this%ncid, this%vars(i)%ncid, actual_dimids) , __LINE__)
      do ii=1, actual_dimcount
        exp_dimid = this%dims( this%vars(i)%dim_indices(ii) )%ncid
        call assert(exp_dimid == actual_dimids(ii), __LINE__)
      end do
    end do
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
