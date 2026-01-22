module g_oce_2_reg
  use o_DATA_TYPES
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS  
  implicit none
  save

  integer                  	::   reg_nx
  integer                  	::   reg_ny
  real(kind=8)             	::   LonMin= -180.  !-75.
  real(kind=8)             	::   LonMax=  180.  !-30.
  real(kind=8)             	::   LatMin=  -80.  !30.
  real(kind=8)             	::   LatMax=   90.   !65.
  real(kind=8)		       	::   RegDx=.25, RegDy=.25
  real(kind=8), allocatable     ::   reg_lon(:), reg_lat(:)
  type(sparse_matrix)           ::   oce_2_reg, oce_2_reg_el
  character(MAX_PATH)           ::   o2r_filename='oce_2_reg.bin_25'
  logical                  	::   do_mesh=.false.
namelist /regular_mesh/ LonMin, LonMax, LatMin, LatMax, RegDx, RegDy, o2r_filename, do_mesh

end module g_oce_2_reg
!
!----------------------------------------------------------------------------
!
