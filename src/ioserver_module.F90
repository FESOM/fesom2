module ioserver_module
  use MOD_MESH

  implicit none  
  public ioserver, comm_fesom_with_ioserver, ioserver_rank, ioserver_mesh, comm_fesom_npes
  private
  
  logical, save :: is_ioserver_ = .false.
  integer, save :: comm_fesom_with_ioserver
  integer, save :: comm_fesom_npes
  integer, save :: ioserver_rank
  type(t_mesh), save :: ioserver_mesh

  type :: ioserver_type
  contains
    procedure, nopass :: enable
    procedure, nopass :: is_ioserver
  end type
  type(ioserver_type) ioserver

contains

  subroutine enable()
    is_ioserver_ = .true.
  end subroutine
  
  pure function is_ioserver() result(x)
    logical x
    x = is_ioserver_
  end function


end module
