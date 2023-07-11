!==========================================================
MODULE MOD_MESH
USE O_PARAM
USE,     intrinsic    :: ISO_FORTRAN_ENV
IMPLICIT NONE
SAVE
integer, parameter    :: MAX_ADJACENT=32 ! Max allowed number of adjacent nodes

TYPE SPARSE_MATRIX 
     integer :: nza
     integer :: dim
     real(kind=WP), allocatable,    dimension(:) :: values
     integer(int32), allocatable,   dimension(:) :: colind
     integer(int32), allocatable,   dimension(:) :: rowptr
     integer(int32), allocatable,   dimension(:) :: colind_loc
     integer(int32), allocatable,   dimension(:) :: rowptr_loc
END TYPE SPARSE_MATRIX

TYPE T_MESH
integer                                     :: nod2D           ! the number of 2D nodes
real(kind=WP)                               :: ocean_area, ocean_areawithcav
real(kind=WP), allocatable, dimension(:,:)  :: coord_nod2D, geo_coord_nod2D
integer                                     :: edge2D       ! the number of 2D edges
integer                                     :: edge2D_in    ! the number of internal 2D edges
integer                                     :: elem2D       ! the number of 2D elements
integer, allocatable, dimension(:,:)        :: elem2D_nodes ! elem2D_nodes(:,n) lists; 3 nodes of element n   
integer, allocatable, dimension(:,:)        :: edges        ! edge(:,n) lists 2 nodes; edge n
integer, allocatable, dimension(:,:)        :: edge_tri     ! edge_tri(:,n) lists 2 
                                                            ! elements containing edge n: the first one is to left 
                                                            ! of the line directed to the second node
integer,       allocatable, dimension(:,:)  :: elem_edges   ! elem_edges(:,n) are edges of element n.  
real(kind=WP), allocatable, dimension(:)    :: elem_area
real(kind=WP), allocatable, dimension(:,:)  :: edge_dxdy, edge_cross_dxdy
real(kind=WP), allocatable, dimension(:)    :: elem_cos, metric_factor
integer,       allocatable, dimension(:,:)  :: elem_neighbors
integer,       allocatable, dimension(:,:)  :: nod_in_elem2D
real(kind=WP), allocatable, dimension(:,:)  :: x_corners, y_corners ! cornes for the scalar points
integer,       allocatable, dimension(:)    :: nod_in_elem2D_num
real(kind=WP), allocatable, dimension(:)    :: depth                ! depth(n) is the depths at node n 
real(kind=WP), allocatable, dimension(:,:)  :: gradient_vec 
                                                           ! coefficients of linear reconstruction
                                                           ! of velocities on elements
real(kind=WP), allocatable, dimension(:,:)  :: gradient_sca ! Coefficients to compute gradient of scalars
                                                           ! on elements
INTEGER,       ALLOCATABLE, DIMENSION(:)    :: bc_index_nod2D(:)
                                                           ! vertical structure     
!
!
!___vertical mesh info__________________________________________________________                                                           
! total number of layers
integer                                     :: nl

! initial layer, mid-depth layer and element depth
real(kind=WP), allocatable, dimension(:)    :: zbar, Z,elem_depth

! upper boudnary index of all vertical vertice/element loops, default==1 but when 
! cavity is used becomes index of cavity-ocean boundary at vertices and elements
integer,       allocatable, dimension(:)    :: ulevels, ulevels_nod2D, ulevels_nod2D_max

! number of levels at elem and vertices considering bottom topography
integer,       allocatable, dimension(:)    :: nlevels, nlevels_nod2D, nlevels_nod2D_min

!
!
!___horizontal mesh info________________________________________________________
real(kind=WP), allocatable, dimension(:,:)  :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), allocatable, dimension(:)    :: mesh_resolution

!
!
!___cavity mesh info____________________________________________________________
! level index of cavity-ocean boundary at vertices and elements
! --> see: ulevels, ulevels_nod2D (fvom_main) 

! vertice/element yes=1/no=0 flag if cavity exists
integer,       allocatable, dimension(:)    :: cavity_flag_n, cavity_flag_e

! depth of cavity-ocean interface
real(kind=WP), allocatable, dimension(:)    :: cavity_depth


real(kind=WP), allocatable, dimension(:,:)  :: cavity_nrst_cavlpnt_xyz

!
!
!___Elevation stiffness matrix__________________________________________________
type(sparse_matrix)                         :: ssh_stiff

!#if defined (__oasis)
real(kind=WP), allocatable, dimension(:)    :: lump2d_south, lump2d_north  
integer,       allocatable, dimension(:)    :: ind_south, ind_north    
!#endif  

  character(:), allocatable :: representative_checksum
END TYPE T_MESH
end module MOD_MESH
!==========================================================

