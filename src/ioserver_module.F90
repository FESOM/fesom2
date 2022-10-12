module ioserver_module

  implicit none  
  public ioserver
  private
  
  logical, save :: is_ioserver_ = .false.

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
