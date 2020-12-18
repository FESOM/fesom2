module io_netcdf_attribute_module
  implicit none
  public att_type
  private
  
  
  type att_type
    character(:), allocatable :: name
    character(:), allocatable :: text
    ! todo: make this work for other data types like int
  end type


contains


end module
