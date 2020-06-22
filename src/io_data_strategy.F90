module io_data_strategy_module
  implicit none
  public data_strategy_type, data_strategy_nf_float_type, data_strategy_nf_double_type
  private

  ! there seems to be no way to get this from the netcdf library except a call to nf_inq_type, which requires a netcdf file handle
  integer, parameter :: nf_float_precision = 4
  integer, parameter :: nf_double_precision = 8


  type, abstract :: data_strategy_type
  contains
    procedure(initialize), deferred, nopass :: initialize
    procedure(netcdf_type), deferred, nopass :: netcdf_type
  end type
  interface
    subroutine initialize()
    end subroutine
    
    function netcdf_type()
      integer netcdf_type
    end function
  end interface


  type, extends(data_strategy_type) :: data_strategy_nf_float_type
  contains
    procedure, nopass :: initialize => initialize_float
    procedure, nopass :: netcdf_type => netcdf_type_float
  end type


  type, extends(data_strategy_type) :: data_strategy_nf_double_type
  contains
    procedure, nopass :: initialize => initialize_double
    procedure, nopass :: netcdf_type => netcdf_type_double
  end type


contains


  subroutine initialize_float()
  end subroutine


  subroutine initialize_double()
  end subroutine


  function netcdf_type_float()
    include "netcdf.inc"
    integer netcdf_type_float
    netcdf_type_float = nf_float
  end function


  function netcdf_type_double()
    include "netcdf.inc"
    integer netcdf_type_double
    netcdf_type_double = nf_double
  end function


end module
