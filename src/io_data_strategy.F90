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
    procedure(add_values), deferred :: add_values
  end type
  interface
    subroutine initialize()
    end subroutine
    
    function netcdf_type()
      integer netcdf_type
    end function

    subroutine add_values(this, values)
      use o_PARAM, only : WP
      import data_strategy_type
      class(data_strategy_type) this
      real(kind=WP), dimension(:,:) :: values
    end subroutine
  end interface


  type, extends(data_strategy_type) :: data_strategy_nf_float_type
    private
    real(nf_float_precision), allocatable, dimension(:,:) :: local_values
  contains
    procedure, nopass :: initialize => initialize_float
    procedure, nopass :: netcdf_type => netcdf_type_float
    procedure :: add_values => add_values_float
  end type


  type, extends(data_strategy_type) :: data_strategy_nf_double_type
    private
    real(nf_double_precision), allocatable, dimension(:,:) :: local_values
  contains
    procedure, nopass :: initialize => initialize_double
    procedure, nopass :: netcdf_type => netcdf_type_double
    procedure :: add_values => add_values_double
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


  subroutine add_values_float(this, values)
    use o_PARAM, only : WP
    class(data_strategy_nf_float_type) this
    real(kind=WP), dimension(:,:) :: values
    ! EO args
    this%local_values = this%local_values + real(values) ! todo: check if casting the array to real here creates a copy of the array
  end subroutine


  subroutine add_values_double(this, values)
    use o_PARAM, only : WP
    class(data_strategy_nf_double_type) this
    real(kind=WP), dimension(:,:) :: values
    ! EO args
    this%local_values = this%local_values + values
  end subroutine

end module
