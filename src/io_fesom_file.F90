 ! synopsis: generic implementation to asynchronously read/write FESOM mesh variable(s) with distributed cell or element data in 2D or 3D to/from a NetCDF file
module io_fesom_file_module
  use io_netcdf_file_module
  implicit none
  public fesom_file_type
  private
  
  
  type var_info
    integer var_index
    real(kind=8), pointer :: local_data_ptr3(:,:) => null()
    real(kind=8), allocatable :: global_level_data(:)
    integer :: global_level_data_size = 0
    logical is_elem_based
  end type
  
  
  type dim_info
    integer idx
    integer len ! better query the len from the netcdf_file_type ?
  end type

  
  type, extends(netcdf_file_type) :: fesom_file_type ! todo maybe: do not inherit but use composition to have different implementations for the iorank and non-io ranks
    private
    integer time_dimidx
    integer time_varidx
    type(var_info) var_infos(20); integer :: nvar_infos = 0 ! todo: allow dynamically allocated size without messing with shallow copied pointers
    type(dim_info), allocatable :: used_mesh_dims(:) ! the dims we add for our variables, we need to identify them when adding our mesh related variables
    integer :: rec_cnt = -1
    integer :: iorank = 0
    integer :: fesom_file_index
  contains
    procedure, public :: read_and_scatter_variables, gather_and_write_variables, init, specify_node_var, is_iorank, rec_count, time_varindex, time_dimindex
    procedure, public :: close_file ! inherited procedures we overwrite
    generic, public :: specify_elem_var => specify_elem_var_2d, specify_elem_var_3d
    procedure, private :: specify_elem_var_2d, specify_elem_var_3d
  end type
  
  
  integer, save :: m_nod2d
  integer, save :: m_elem2d
  integer, save :: m_nl
  

  type fesom_file_type_ptr
    class(fesom_file_type), pointer :: ptr
  end type
  type(fesom_file_type_ptr), allocatable, save :: all_fesom_files(:)


contains


  function is_iorank(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(in) :: this
    logical x
    x = (mype == this%iorank)
  end function


  ! return the number of timesteps of the file if a file is attached or return the default value of -1
  function rec_count(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(inout) :: this
    integer x
    ! EO parameters
    integer, allocatable :: time_shape(:)
    
    if(this%rec_cnt == -1 .and. this%is_attached()) then
      ! update from file if rec_cnt has never been used before
      call this%read_var_shape(this%time_varidx, time_shape)
      this%rec_cnt = time_shape(1)
    end if
    
    x = this%rec_cnt
  end function


  function time_varindex(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_varidx
  end function


  function time_dimindex(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_dimidx
  end function
  
  
  subroutine init(f, mesh_nod2d, mesh_elem2d, mesh_nl) ! todo: would like to call it initialize but Fortran is rather cluncky with overwriting base type procedures
    class(fesom_file_type), target, intent(inout) :: f
    integer mesh_nod2d
    integer mesh_elem2d
    integer mesh_nl
    ! EO parameters
    type(fesom_file_type_ptr), allocatable :: tmparr(:)

    ! get hold of our mesh data for later use (assume the mesh instance will not change)
    m_nod2d = mesh_nod2d
    m_elem2d = mesh_elem2d
    m_nl = mesh_nl
    call f%netcdf_file_type%initialize()

    allocate(f%used_mesh_dims(0))

    f%time_dimidx = f%add_dim_unlimited('time')

    f%time_varidx = f%add_var_double('time', [f%time_dimidx])

    ! add this instance to global array
    ! the array is being used to identify the instance in an async call
    if( .not. allocated(all_fesom_files)) then
      allocate(all_fesom_files(1))
    else
      allocate( tmparr(size(all_fesom_files)+1) )
      tmparr(1:size(all_fesom_files)) = all_fesom_files
      deallocate(all_fesom_files)
      call move_alloc(tmparr, all_fesom_files)
    end if
    all_fesom_files(size(all_fesom_files))%ptr => f
    f%fesom_file_index = size(all_fesom_files)
  end subroutine
  
  
  subroutine read_and_scatter_variables(f)
    use g_PARSUP
    use io_scatter_module
    class(fesom_file_type), target :: f
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    integer last_rec_idx
    type(var_info), pointer :: var
    real(kind=8), allocatable :: laux(:)
  
    last_rec_idx = f%rec_count()
    
    do i=1, f%nvar_infos
      var => f%var_infos(i)
    
      nlvl = size(var%local_data_ptr3,dim=1)
      is_2d = (nlvl == 1)
      allocate(laux( size(var%local_data_ptr3,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D

      if(mype == f%iorank) then
        ! todo: choose how many levels we read at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
        if(mype == f%iorank) then
          if(is_2d) then
            call f%read_var(var%var_index, [1,last_rec_idx], [size(var%global_level_data),1], var%global_level_data)
          else
            ! z,nod,time
            call f%read_var(var%var_index, [lvl,1,last_rec_idx], [1,size(var%global_level_data),1], var%global_level_data)
          end if
        end if

        if(var%is_elem_based) then
          call scatter_elem2D(var%global_level_data, laux, f%iorank, MPI_comm_fesom)
        else
          call scatter_nod2D(var%global_level_data, laux, f%iorank, MPI_comm_fesom)
        end if
        ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
        var%local_data_ptr3(lvl,:) = laux ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)
      end do
      deallocate(laux)
    end do
  end subroutine


  subroutine gather_and_write_variables(f)
    use g_PARSUP
    use io_gather_module
    class(fesom_file_type), target :: f
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    real(kind=8), allocatable :: laux(:)
    type(var_info), pointer :: var

    if(f%is_iorank()) f%rec_cnt = f%rec_count()+1
    
    do i=1, f%nvar_infos
      var => f%var_infos(i)

      nlvl = size(var%local_data_ptr3,dim=1)
      is_2d = (nlvl == 1)
      allocate(laux( size(var%local_data_ptr3,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D

      if(mype == f%iorank) then
        ! todo: choose how many levels we write at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
        ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
        laux = var%local_data_ptr3(lvl,:) ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)

        if(var%is_elem_based) then
          call gather_elem2D(laux, var%global_level_data, f%iorank, 42, MPI_comm_fesom)
        else
          call gather_nod2D (laux, var%global_level_data, f%iorank, 42, MPI_comm_fesom)
        end if

        if(mype == f%iorank) then
          if(is_2d) then
            call f%write_var(var%var_index, [1,f%rec_cnt], [size(var%global_level_data),1], var%global_level_data)
          else
            ! z,nod,time
            call f%write_var(var%var_index, [lvl,1,f%rec_cnt], [1,size(var%global_level_data),1], var%global_level_data)
          end if
        end if
      end do
      deallocate(laux)
    end do
  end subroutine


  subroutine specify_node_var(f, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    use g_PARSUP
    class(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(..) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: local_data_ptr3(:,:)
    type(dim_info) level_diminfo, depth_diminfo

    level_diminfo = obtain_diminfo(f, m_nod2d)
   
    if(size(shape(local_data)) == 1) then ! 1D data
      call c_f_pointer(c_loc(local_data), local_data_ptr3, [1,size(local_data)])
      call specify_variable(f, name, [level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data_ptr3, .false., longname, units)
    
    else if(size(shape(local_data)) == 2) then ! 2D data
      depth_diminfo = obtain_diminfo(f, size(local_data, dim=1))
      call c_f_pointer(c_loc(local_data), local_data_ptr3, [size(local_data, dim=1),size(local_data, dim=2)])
      call specify_variable(f, name, [depth_diminfo%idx, level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data_ptr3, .false., longname, units)
    end if        
  end subroutine


  subroutine specify_elem_var_2d(f, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    use g_PARSUP
    class(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: local_data_ptr3(:,:)
    type(dim_info) level_diminfo

    level_diminfo = obtain_diminfo(f, m_elem2d)

    local_data_ptr3(1:1,1:size(local_data)) => local_data
    call specify_variable(f, name, [level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data_ptr3, .true., longname, units)    
  end subroutine


  subroutine specify_elem_var_3d(f, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    use g_PARSUP
    class(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, depth_diminfo

    level_diminfo = obtain_diminfo(f, m_elem2d)
    depth_diminfo = obtain_diminfo(f, size(local_data, dim=1))
    
    call specify_variable(f, name, [depth_diminfo%idx, level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data, .true., longname, units)
  end subroutine
  
  
  function obtain_diminfo(f, len) result(info)
    type(fesom_file_type), intent(inout) :: f
    type(dim_info) info
    integer len
    ! EO parameters
    integer i
    type(dim_info), allocatable :: tmparr(:)
    
    do i=1, size(f%used_mesh_dims)
      if(f%used_mesh_dims(i)%len == len) then
        info = f%used_mesh_dims(i)
        return
      end if
    end do
    
    ! the dim has not been added yet, see if it is one of our allowed mesh related dims
    if(len == m_nod2d) then
      info = dim_info( idx=f%add_dim('node', len), len=len)
    else if(len == m_elem2d) then
      info = dim_info( idx=f%add_dim('elem', len), len=len)
    else if(len == m_nl-1) then
      info = dim_info( idx=f%add_dim('nz_1', len), len=len)
    else if(len == m_nl) then
      info = dim_info( idx=f%add_dim('nz', len), len=len)
    else
      print *, "error in line ",__LINE__, __FILE__," can not find dimension with size",len
      stop 1
    end if
    
    ! append the new dim to our list of used dims, i.e. the dims we use for the mesh based variables created via #specify_variable
    ! assume the used_mesh_dims array is allocated
    allocate( tmparr(size(f%used_mesh_dims)+1) )
    tmparr(1:size(f%used_mesh_dims)) = f%used_mesh_dims
    deallocate(f%used_mesh_dims)
    call move_alloc(tmparr, f%used_mesh_dims)
    f%used_mesh_dims( size(f%used_mesh_dims) ) = info
  end function


  subroutine specify_variable(f, name, dim_indices, global_level_data_size, local_data, is_elem_based, longname, units)
    use g_PARSUP
    type(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer global_level_data_size
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision?
    logical, intent(in) :: is_elem_based
    character(len=*), intent(in) :: units, longname
    ! EO parameters
    integer var_index

    var_index = f%add_var_double(name, dim_indices)
    call f%add_var_att(var_index, "units", units)
    call f%add_var_att(var_index, "long_name", longname)
    
    call assert(f%nvar_infos < size(f%var_infos), __LINE__)
    f%nvar_infos = f%nvar_infos+1
    f%var_infos(f%nvar_infos)%var_index = var_index
    f%var_infos(f%nvar_infos)%local_data_ptr3 => local_data
    f%var_infos(f%nvar_infos)%global_level_data_size = global_level_data_size
    f%var_infos(f%nvar_infos)%is_elem_based = is_elem_based
  end subroutine
  
  
  subroutine close_file(this)
    class(fesom_file_type), intent(inout) :: this
   
    this%rec_cnt = -1 ! reset state (should probably be done in all the open_ procedures, not here)
    call this%netcdf_file_type%close_file()
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
