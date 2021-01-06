 ! synopsis: generic implementation to asynchronously read/write FESOM mesh variable(s) with distributed cell or element data in 2D or 3D to/from a NetCDF file
module io_fesom_file_module
  use mod_mesh
  use io_netcdf_file_module
  implicit none
  public fesom_file_type
  private
  
  
  type var_info
    integer var_index
    real(kind=8), pointer :: local_data_ptr3(:,:) => null()
    real(kind=8), allocatable :: global_level_data(:)
    integer :: global_level_data_size = 0
  end type
  
  
  type dim_info
    integer idx
    integer len ! better query the len from the netcdf_file_type ?
  end type

  
  type, extends(netcdf_file_type) :: fesom_file_type
    private
    integer time_dimidx
    integer time_varidx
    type(var_info) var_infos(20); integer :: nvar_infos = 0 ! todo: allow dynamically allocated size without messing with shallow copied pointers
    type(dim_info), allocatable :: dim_infos(:)
    integer :: rec_cnt = 0
    integer :: iorank = 0
  contains
    procedure, public :: gather_and_write, init, specify_node_var, is_iorank, rec_count, time_varindex, time_dimindex
  end type
  
  
  type(t_mesh), save :: mesh


contains


  function is_iorank(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(in) :: this
    logical x
    x = (mype == this%iorank)
  end function


  function rec_count(this) result(x)
    use g_PARSUP
    class(fesom_file_type), intent(in) :: this
    integer x
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


  subroutine init(f, mesh_) ! todo: would like to call it initialize but Fortran is rather cluncky with overwriting base type procedures
    use mod_mesh
    use o_arrays
    class(fesom_file_type), intent(inout) :: f
    type(t_mesh), intent(in) :: mesh_
    ! EO parameters

    mesh = mesh_ ! get hold of our mesh for later use (assume the mesh instance will never change)
    
    f%rec_cnt = 0

    call f%netcdf_file_type%initialize()

    ! add the dimensions we intend to use to the file spec and also store here so we can use them when creating the variables
    ! todo: store in a separate "dim pool" without calling f%add_dim and add only if a variable requires it
    allocate(f%dim_infos(4))
    f%dim_infos(1) = dim_info( idx=f%add_dim('node', mesh%nod2d), len=mesh%nod2d)
    f%dim_infos(2) = dim_info( idx=f%add_dim('elem', mesh%elem2d), len=mesh%elem2d)
    f%dim_infos(3) = dim_info( idx=f%add_dim('nz_1', mesh%nl-1), len=mesh%nl-1)
    f%dim_infos(4) = dim_info( idx=f%add_dim('nz', mesh%nl), len=mesh%nl)

    f%time_dimidx = f%add_dim_unlimited('time')

    f%time_varidx = f%add_var_double('time', [f%time_dimidx])
  end subroutine


  subroutine gather_and_write(f)
    use g_PARSUP
    use io_gather_module
    class(fesom_file_type), intent(inout) :: f
    ! EO parameters
    integer i,lvl, nlvl, nodes_per_lvl
    logical is_2d
  
    f%rec_cnt = f%rec_cnt+1
    
    do i=1, f%nvar_infos
    
call assert(associated(f%var_infos(i)%local_data_ptr3), __LINE__)

      nlvl = size(f%var_infos(i)%local_data_ptr3,dim=1)
      nodes_per_lvl = f%var_infos(i)%global_level_data_size
      is_2d = (nlvl == 1)

      if(mype == f%iorank) then
        ! todo: choose how many levels we write at once
        if(.not. allocated(f%var_infos(i)%global_level_data)) allocate(f%var_infos(i)%global_level_data(nodes_per_lvl))
      end if

      lvl=1 ! todo: loop lvls
      call gather_nod2D(f%var_infos(i)%local_data_ptr3(lvl,:), f%var_infos(i)%global_level_data, f%iorank, 42, MPI_comm_fesom) 
      if(mype == f%iorank) then
        if(is_2d) then
          call f%write_var(f%var_infos(i)%var_index, [1,f%rec_cnt], [size(f%var_infos(i)%global_level_data),1], f%var_infos(i)%global_level_data)
        else
          ! z,nod,time
          call f%write_var(f%var_infos(i)%var_index, [lvl,1,f%rec_cnt], [1,size(f%var_infos(i)%global_level_data),1], f%var_infos(i)%global_level_data)
        end if
     end if
    end do
  end subroutine


  subroutine specify_node_var(f, name, local_data, longname, units)
    use, intrinsic :: ISO_C_BINDING
    use g_PARSUP
    class(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
    real(kind=8), target, intent(inout) :: local_data(..) ! todo: be able to set precision
    character(len=*), intent(in) :: units, longname
    ! EO parameters
    real(8), pointer :: local_data_ptr3(:,:)
    type(dim_info) level_diminfo, depth_diminfo

    level_diminfo = find_diminfo(f, mesh%nod2d)
   
    if(size(shape(local_data)) == 1) then ! 1D data
      call c_f_pointer(c_loc(local_data), local_data_ptr3, [1,size(local_data)])
      call specify_variable(f, name, [level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data_ptr3, longname, units)
    
    else if(size(shape(local_data)) == 2) then ! 2D data
      depth_diminfo = find_diminfo(f, size(local_data, dim=1))
      call c_f_pointer(c_loc(local_data), local_data_ptr3, [size(local_data, dim=1),size(local_data, dim=2)])
      call specify_variable(f, name, [depth_diminfo%idx, level_diminfo%idx, f%time_dimidx], level_diminfo%len, local_data_ptr3, longname, units)
    end if        
  end subroutine
  
  
  function find_diminfo(f, len) result(info)
    type(fesom_file_type), intent(inout) :: f
    type(dim_info) info
    integer len
    ! EO parameters
    integer i
    
    do i=1, size(f%dim_infos)
      if(f%dim_infos(i)%len == len) then
        info = f%dim_infos(i)
        return
      end if
    end do
    
    print *, "error in line ",__LINE__, __FILE__," can not find dimension with size",len
    stop 1
  end function


  subroutine specify_variable(f, name, dim_indices, global_level_data_size, local_data, longname, units)
    use g_PARSUP
    type(fesom_file_type), intent(inout) :: f
    character(len=*), intent(in) :: name
!    integer, intent(in) :: global_shape(:)
    integer, intent(in) :: dim_indices(:)
    integer global_level_data_size
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision?
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
