! helper module to treat split restart files similar as the previous single-file ones
module restart_file_group_module
  use io_fesom_file_module
  implicit none
  public restart_file_group
  private

  
  type, extends(fesom_file_type) :: restart_file_type
    integer iter_varindex
    character(:), allocatable :: varname
    character(:), allocatable :: path
    logical must_exist_on_read
  end type


  type restart_file_group
    private
    type(restart_file_type), public :: files(20)
    integer, public :: nfiles = 0 ! todo: allow dynamically allocated size without messing with shallow copied pointers
  contains
    generic, public :: def_node_var => def_node_var_2d, def_node_var_3d
    generic, public :: def_elem_var => def_elem_var_2d, def_elem_var_3d
    procedure, private :: def_node_var_2d, def_node_var_3d
    procedure, private :: def_elem_var_2d, def_elem_var_3d
    ! def_*_optional procedures create a restart variable which does not have to exist when reading the restart file
    generic, public :: def_node_var_optional => def_node_var_2d_optional, def_node_var_3d_optional
    generic, public :: def_elem_var_optional => def_elem_var_2d_optional, def_elem_var_3d_optional
    procedure, private :: def_node_var_2d_optional, def_node_var_3d_optional
    procedure, private :: def_elem_var_2d_optional, def_elem_var_3d_optional
  end type
  
contains


  subroutine def_node_var_2d(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .true., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_node_var(name, longname, units, local_data)
  end subroutine


  subroutine def_node_var_3d(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .true., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_node_var(name, longname, units, local_data)
  end subroutine


  subroutine def_elem_var_2d(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .true., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_elem_var(name, longname, units, local_data)
  end subroutine


  subroutine def_elem_var_3d(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .true., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_elem_var(name, longname, units, local_data)
  end subroutine


  subroutine add_file(g, name, must_exist_on_read, mesh_nod2d, mesh_elem2d, mesh_nl)
    class(restart_file_group), target, intent(inout) :: g
    character(len=*), intent(in) :: name
    logical must_exist_on_read
    integer mesh_nod2d, mesh_elem2d, mesh_nl
    ! EO parameters
    type(restart_file_type), pointer :: f

    call assert(g%nfiles < size(g%files), __LINE__)
    g%nfiles = g%nfiles+1
    f => g%files(g%nfiles)
    
    f%path = ""
    f%varname = name
    f%must_exist_on_read = must_exist_on_read
    call f%fesom_file_type%init(mesh_nod2d, mesh_elem2d, mesh_nl)
    ! this is specific for a restart file
    f%iter_varindex = f%add_var_int('iter', [f%time_dimindex()])    
  end subroutine


  subroutine def_node_var_2d_optional(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .false., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_node_var(name, longname, units, local_data)
  end subroutine


  subroutine def_node_var_3d_optional(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .false., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_node_var(name, longname, units, local_data)
  end subroutine


  subroutine def_elem_var_2d_optional(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .false., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_elem_var(name, longname, units, local_data)
  end subroutine


  subroutine def_elem_var_3d_optional(this, name, longname, units, local_data, mesh)
    use mod_mesh
    class(restart_file_group), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units, longname
    real(kind=8), target, intent(inout) :: local_data(:,:) ! todo: be able to set precision
    type(t_mesh), intent(in) :: mesh
    ! EO parameters

    call add_file(this, name, .false., mesh%nod2d, mesh%elem2d, mesh%nl)
    call this%files(this%nfiles)%specify_elem_var(name, longname, units, local_data)
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

end module
