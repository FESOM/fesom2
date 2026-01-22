module o_param
  implicit none
  save
  !    
  ! *** Fixed parameters ***
  real(kind=8), parameter  	:: pi=3.141592653589793, rad=pi/180.0  
  real(kind=8), parameter  	:: omega=2.0*pi/(24.0*60.0*60.0)
  real(kind=8), parameter  	:: g=9.81                       ![m/s^2]
  real(kind=8), parameter  	:: r_earth=6.3675e6             ![m]
  real(kind=8), parameter  	:: rho0=1030.                   ![kg/m^3]
  real(kind=8), parameter 	:: rho0r=1.0/rho0 
  real(kind=8), parameter  	:: vcpw=4.2e6                   ![J/m^3/K]volum. heat cap. of water
  real(kind=8), parameter	:: small=1.0e-8                 !small value
end module o_param
!
!----------------------------------------------------------------------
!
module o_array
  implicit none
  save
  real(kind=8), allocatable, dimension(:)           :: density_ref, density_insitu  
  real(kind=8), allocatable, dimension(:)           :: hpressure
  real(kind=8), allocatable, target, dimension(:,:) :: tracer_rhs, tracer, dtracer
end module o_array
!
!----------------------------------------------------------------------
!
module o_DATA_TYPES
  implicit none
  save
  !
  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=8), pointer, dimension(:)      :: values
     integer(KIND=4), pointer,   dimension(:) :: colind
     integer(KIND=4), pointer,   dimension(:) :: rowptr
  end type sparse_matrix
  !
  type addresstype
     integer                                :: nmb
     integer(KIND=4), dimension(:), pointer :: addresses
  end type addresstype
  !
end module o_DATA_TYPES
!
!----------------------------------------------------------------------
!
module o_mesh
  !
  use o_DATA_TYPES
  implicit none
  save
  !
  integer, allocatable, dimension(:)           :: col_pos
  !
  integer                                      :: nod2D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod2D  
  integer, allocatable, dimension(:)           :: index_nod2D  
  real(kind=8), allocatable, dimension(:)      :: cos_elem2D 
  real(kind=8), allocatable, dimension(:)      :: depths
  !
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D     
  type(addresstype), allocatable, dimension(:) :: nod_in_opbnd_tri
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D
  type(addresstype), allocatable, dimension(:) :: nghbr_elem2D
  !
  integer, allocatable, dimension(:)           :: num_layers_below_nod2D
  !
  integer                                    ::   edge2D     ! the number of 2D edges
  integer                                    ::   edge2D_in  
  integer, allocatable, dimension(:,:)       ::   edges
  integer, allocatable, dimension(:,:)       ::   edge_tri
  integer, allocatable, dimension(:,:)       ::   elem_edges
  real(kind=8), allocatable, dimension(:,:)  ::   edge_dxdy, edge_cross_dxdy
  real(kind=8),allocatable,dimension(:,:)    ::   gradient_vec
  real(kind=8),allocatable,dimension(:,:)    ::   gradient_sca
  real(kind=8), allocatable, dimension(:)    ::   elem_cos, metric_factor
  integer                                    ::   nl_1, max_num_layers
  integer(KIND=4), allocatable, dimension(:) ::   nlvls, elvls

end module o_mesh
!
!----------------------------------------------------------------------------
!
module o_elements
  implicit none
  save
  integer                                      :: elem2D
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nodes 
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nghbrs 
  !
  real(kind=8)                                 :: Vol2D
  real(kind=8)                                 :: sProd_2Di
  real(kind=8), allocatable, dimension(:,:)    :: sProd_2Dij
  real(kind=8), allocatable, dimension(:,:)    :: derivative_stdbafu_x_2D 
  real(kind=8), allocatable, dimension(:,:)    :: bafux_2d, bafuy_2d
  real(kind=8), allocatable, dimension(:)      :: voltriangle
  real(kind=8), allocatable, dimension(:,:)    :: area
  integer,allocatable,dimension(:,:)           :: elem_neighbors
  !
  real(kind=8)                                 :: ocean_area
  real(kind=8),allocatable, dimension(:)       :: cluster_area_2D
 
end module o_elements
!
