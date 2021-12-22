 ! synopsis: generic implementation to asynchronously read/write FESOM mesh variable(s) with distributed cell or element data in 2D or 3D to/from a NetCDF file
module io_fesom_file_module
  use io_netcdf_file_module
  use async_threads_module
  implicit none
  public fesom_file_type
  private
  
  
  type var_info
    integer var_index
    real(kind=8), pointer :: external_local_data_ptr(:,:) => null()
    real(kind=8), allocatable, dimension(:,:) :: local_data_copy
    real(kind=8), allocatable :: global_level_data(:)
    integer :: global_level_data_size = 0
    logical is_elem_based
    character(:), allocatable :: varname ! todo: maybe use a getter in netcdf_file_type to get the name
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
    type(thread_type) thread
    logical :: thread_running = .false.
    integer :: comm
    logical gather_and_write
  contains
    procedure, public :: async_read_and_scatter_variables, async_gather_and_write_variables, join, init, is_iorank, rec_count, time_varindex, time_dimindex
    procedure, public :: read_variables_raw, write_variables_raw
    procedure, public :: close_file ! inherited procedures we overwrite
    generic, public :: specify_node_var => specify_node_var_2d, specify_node_var_3d
    generic, public :: specify_elem_var => specify_elem_var_2d, specify_elem_var_3d
    procedure, private :: specify_node_var_2d, specify_node_var_3d
    procedure, private :: specify_elem_var_2d, specify_elem_var_3d
    procedure, private :: read_and_scatter_variables, gather_and_write_variables
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
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_varidx
  end function


  function time_dimindex(this) result(x)
    class(fesom_file_type), intent(in) :: this
    integer x
    x = this%time_dimidx
  end function
  
  
  subroutine init(this, mesh_nod2d, mesh_elem2d, mesh_nl) ! todo: would like to call it initialize but Fortran is rather cluncky with overwriting base type procedures
    use g_PARSUP
    use io_netcdf_workaround_module
    use io_gather_module
    class(fesom_file_type), target, intent(inout) :: this
    integer mesh_nod2d
    integer mesh_elem2d
    integer mesh_nl
    ! EO parameters
    type(fesom_file_type_ptr), allocatable :: tmparr(:)
    logical async_netcdf_allowed
    integer err
    integer provided_mpi_thread_support_level

    call init_io_gather()

    ! get hold of our mesh data for later use (assume the mesh instance will not change)
    m_nod2d = mesh_nod2d
    m_elem2d = mesh_elem2d
    m_nl = mesh_nl
    call this%netcdf_file_type%initialize()

    allocate(this%used_mesh_dims(0))

    this%time_dimidx = this%add_dim_unlimited('time')

    this%time_varidx = this%add_var_double('time', [this%time_dimidx])

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
    all_fesom_files(size(all_fesom_files))%ptr => this
    this%fesom_file_index = size(all_fesom_files)

    ! set up async output
    
    this%iorank = next_io_rank(MPI_COMM_FESOM, async_netcdf_allowed)

    call MPI_Comm_dup(MPI_COMM_FESOM, this%comm, err)

    call this%thread%initialize(async_worker, this%fesom_file_index)
    if(.not. async_netcdf_allowed) call this%thread%disable_async()
  
    ! check if we have multi thread support available in the MPI library
    ! tough MPI_THREAD_FUNNELED should be enough here, at least on cray-mpich 7.5.3 async mpi calls fail if we do not have support level 'MPI_THREAD_MULTIPLE'
    ! on cray-mpich we only get level 'MPI_THREAD_MULTIPLE' if 'MPICH_MAX_THREAD_SAFETY=multiple' is set in the environment
    call MPI_Query_thread(provided_mpi_thread_support_level, err)
    if(provided_mpi_thread_support_level < MPI_THREAD_MULTIPLE) call this%thread%disable_async()    
  end subroutine
  
  
  subroutine read_and_scatter_variables(this)
    use io_scatter_module
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    integer last_rec_idx
    type(var_info), pointer :: var
    real(kind=8), allocatable :: laux(:)
    integer mpierr

    last_rec_idx = this%rec_count()
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
    
      nlvl = size(var%external_local_data_ptr,dim=1)
      is_2d = (nlvl == 1)
      allocate(laux( size(var%external_local_data_ptr,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D

      if(this%is_iorank()) then
        ! todo: choose how many levels we read at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
#ifdef ENABLE_ALEPH_CRAYMPICH_WORKAROUNDS
        ! aleph cray-mpich workaround
        call MPI_Barrier(this%comm, mpierr)
#endif
        if(this%is_iorank()) then
          if(is_2d) then
            call this%read_var(var%var_index, [1,last_rec_idx], [size(var%global_level_data),1], var%global_level_data)
          else
            ! z,nod,time
#ifdef UNTRANSPOSE_RESTART
            call this%read_var(var%var_index, [lvl,1,last_rec_idx], [1,size(var%global_level_data),1], var%global_level_data) ! untransposed
#else
            call this%read_var(var%var_index, [1,lvl,last_rec_idx], [size(var%global_level_data),1,1], var%global_level_data)
#endif
          end if
        end if

        if(var%is_elem_based) then
          call scatter_elem2D(var%global_level_data, laux, this%iorank, this%comm)
        else
          call scatter_nod2D(var%global_level_data, laux, this%iorank, this%comm)
        end if
        ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
       var%external_local_data_ptr(lvl,:) = laux ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)
      end do
      deallocate(laux)
    end do
  end subroutine


  subroutine gather_and_write_variables(this)
    use io_gather_module
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i,lvl, nlvl
    logical is_2d
    real(kind=8), allocatable :: laux(:)
    type(var_info), pointer :: var
    integer mpierr

    if(this%is_iorank()) this%rec_cnt = this%rec_count()+1
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)

      nlvl = size(var%local_data_copy,dim=1)
      is_2d = (nlvl == 1)
      allocate(laux( size(var%local_data_copy,dim=2) )) ! i.e. myDim_elem2D+eDim_elem2D or myDim_nod2D+eDim_nod2D

      if(this%is_iorank()) then
        ! todo: choose how many levels we write at once
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( var%global_level_data_size ))
      else
        if(.not. allocated(var%global_level_data)) allocate(var%global_level_data( 0 ))
      end if

      do lvl=1, nlvl
#ifdef ENABLE_ALEPH_CRAYMPICH_WORKAROUNDS
        ! aleph cray-mpich workaround
        call MPI_Barrier(this%comm, mpierr)
#endif
        ! the data from our pointer is not contiguous (if it is 3D data), so we can not pass the pointer directly to MPI
        laux = var%local_data_copy(lvl,:) ! todo: remove this buffer and pass the data directly to MPI (change order of data layout to be levelwise or do not gather levelwise but by columns)

        if(var%is_elem_based) then
          call gather_elem2D(laux, var%global_level_data, this%iorank, 42, this%comm)
        else
          call gather_nod2D (laux, var%global_level_data, this%iorank, 42, this%comm)
        end if

        if(this%is_iorank()) then
          if(is_2d) then
            call this%write_var(var%var_index, [1,this%rec_cnt], [size(var%global_level_data),1], var%global_level_data)
          else
            ! z,nod,time
#ifdef UNTRANSPOSE_RESTART
           call this%write_var(var%var_index, [lvl,1,this%rec_cnt], [1,size(var%global_level_data),1], var%global_level_data) ! untransposed
#else
           call this%write_var(var%var_index, [1,lvl,this%rec_cnt], [size(var%global_level_data),1,1], var%global_level_data)
#endif
          end if
        end if
      end do
      deallocate(laux)
    end do
    
    if(this%is_iorank()) call this%flush_file() ! flush the file to disk after each write
  end subroutine


  subroutine read_variables_raw(this, fileunit)
    class(fesom_file_type), target :: this
    integer, intent(in) :: fileunit
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    integer status
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
      read(fileunit) var%external_local_data_ptr ! directly use external_local_data_ptr, use the local_data_copy only when called asynchronously
    end do
  end subroutine


  subroutine write_variables_raw(this, fileunit)
    class(fesom_file_type), target :: this
    integer, intent(in) :: fileunit
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    
    do i=1, this%nvar_infos
      var => this%var_infos(i)
      write(fileunit) var%external_local_data_ptr ! directly use external_local_data_ptr, use the local_data_copy only when called asynchronously
    end do
  end subroutine


  subroutine join(this)
    class(fesom_file_type) this
    ! EO parameters
    
    if(this%thread_running) call this%thread%join()
    this%thread_running = .false.    
  end subroutine


  subroutine async_read_and_scatter_variables(this)
    class(fesom_file_type), target :: this

    call assert(.not. this%thread_running, __LINE__)

    this%gather_and_write = .false.
    call this%thread%run()
    this%thread_running = .true.
  end subroutine


  subroutine async_gather_and_write_variables(this)
    class(fesom_file_type), target :: this
    ! EO parameters
    integer i
    type(var_info), pointer :: var
    
    call assert(.not. this%thread_running, __LINE__)

    ! copy data so we can write the current values asynchronously
    do i=1, this%nvar_infos
      var => this%var_infos(i)
      if(.not. allocated(var%local_data_copy)) allocate( var%local_data_copy(size(var%external_local_data_ptr,dim=1), size(var%external_local_data_ptr,dim=2)) )
      var%local_data_copy = var%external_local_data_ptr
    end do
    
    this%gather_and_write = .true.
    call this%thread%run()
    this%thread_running = .true.
  end subroutine


  subroutine async_worker(fesom_file_index)
    integer, intent(in) :: fesom_file_index
    ! EO parameters
    type(fesom_file_type), pointer :: f

    f => all_fesom_files(fesom_file_index)%ptr

    if(f%gather_and_write) then
      call f%gather_and_write_variables()
    else
      call f%read_and_scatter_variables()
    end if
  end subroutine


  ! use separate procedures to specify node based or element based variables
  ! if we would otherwise specify the vars only via the sizes of their dimensions,
  ! we have to assign the corresponding dimindx somewhere else, which would be error prone
  subroutine specify_node_var_2d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: external_local_data_ptr(:,:)
    type(dim_info) level_diminfo

    level_diminfo = obtain_diminfo(this, m_nod2d)
   
    external_local_data_ptr(1:1,1:size(local_data)) => local_data
    call specify_variable(this, name, [level_diminfo%idx, this%time_dimidx], level_diminfo%len, external_local_data_ptr, .false., longname, units)    
  end subroutine


  subroutine specify_node_var_3d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, depth_diminfo

    level_diminfo = obtain_diminfo(this, m_nod2d)    
    depth_diminfo = obtain_diminfo(this, size(local_data, dim=1))

#ifdef UNTRANSPOSE_RESTART
   call specify_variable(this, name, [depth_diminfo%idx, level_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .false., longname, units) ! untransposed
#else
     call specify_variable(this, name, [level_diminfo%idx, depth_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .false., longname, units)
#endif
  end subroutine


  subroutine specify_elem_var_2d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    ! EO parameters
    real(8), pointer :: external_local_data_ptr(:,:)
    type(dim_info) level_diminfo

    level_diminfo = obtain_diminfo(this, m_elem2d)

    external_local_data_ptr(1:1,1:size(local_data)) => local_data
    call specify_variable(this, name, [level_diminfo%idx, this%time_dimidx], level_diminfo%len, external_local_data_ptr, .true., longname, units)    
  end subroutine


  subroutine specify_elem_var_3d(this, name, longname, units, local_data)
    use, intrinsic :: ISO_C_BINDING
    class(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    ! EO parameters
    type(dim_info) level_diminfo, depth_diminfo

    level_diminfo = obtain_diminfo(this, m_elem2d)
    depth_diminfo = obtain_diminfo(this, size(local_data, dim=1))
    
#ifdef UNTRANSPOSE_RESTART
   call specify_variable(this, name, [depth_diminfo%idx, level_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .true., longname, units) ! untransposed
#else
   call specify_variable(this, name, [level_diminfo%idx, depth_diminfo%idx, this%time_dimidx], level_diminfo%len, local_data, .true., longname, units)
#endif
  end subroutine
  
  
  function obtain_diminfo(this, len) result(info)
    type(fesom_file_type), intent(inout) :: this
    type(dim_info) info
    integer len
    ! EO parameters
    integer i
    type(dim_info), allocatable :: tmparr(:)
    
    do i=1, size(this%used_mesh_dims)
      if(this%used_mesh_dims(i)%len == len) then
        info = this%used_mesh_dims(i)
        return
      end if
    end do
    
    ! the dim has not been added yet, see if it is one of our allowed mesh related dims
    if(len == m_nod2d) then
      info = dim_info( idx=this%add_dim('node', len), len=len)
    else if(len == m_elem2d) then
      info = dim_info( idx=this%add_dim('elem', len), len=len)
    else if(len == m_nl-1) then
      info = dim_info( idx=this%add_dim('nz_1', len), len=len)
    else if(len == m_nl) then
      info = dim_info( idx=this%add_dim('nz', len), len=len)
    else
      print *, "error in line ",__LINE__, __FILE__," can not find dimension with size",len
      stop 1
    end if
    
    ! append the new dim to our list of used dims, i.e. the dims we use for the mesh based variables created via #specify_variable
    ! assume the used_mesh_dims array is allocated
    allocate( tmparr(size(this%used_mesh_dims)+1) )
    tmparr(1:size(this%used_mesh_dims)) = this%used_mesh_dims
    deallocate(this%used_mesh_dims)
    call move_alloc(tmparr, this%used_mesh_dims)
    this%used_mesh_dims( size(this%used_mesh_dims) ) = info
  end function


  subroutine specify_variable(this, name, dim_indices, global_level_data_size, local_data, is_elem_based, longname, units)
    type(fesom_file_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim_indices(:)
    integer global_level_data_size
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision?
    logical, intent(in) :: is_elem_based
    character(len=*), intent(in) :: units, longname
    ! EO parameters
    integer var_index

    var_index = this%add_var_double(name, dim_indices)
    call this%add_var_att(var_index, "units", units)
    call this%add_var_att(var_index, "long_name", longname)
    
    call assert(this%nvar_infos < size(this%var_infos), __LINE__)
    this%nvar_infos = this%nvar_infos+1
    this%var_infos(this%nvar_infos)%var_index = var_index
    this%var_infos(this%nvar_infos)%external_local_data_ptr => local_data
    this%var_infos(this%nvar_infos)%global_level_data_size = global_level_data_size
    this%var_infos(this%nvar_infos)%is_elem_based = is_elem_based
    this%var_infos(this%nvar_infos)%varname = name
  end subroutine
  
  
  subroutine close_file(this)
    class(fesom_file_type), intent(inout) :: this

    if(this%thread_running) call this%thread%join()
    this%thread_running = .false.    
   
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
