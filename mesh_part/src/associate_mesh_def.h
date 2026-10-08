integer         , pointer :: nod2D
integer         , pointer :: elem2D
integer         , pointer :: edge2D
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
real(kind=WP)   , pointer :: ocean_areawithcav
integer         , pointer :: nl
integer         , pointer :: nn_size
real(kind=WP), dimension(:,:), pointer :: coord_nod2D, geo_coord_nod2D
integer, dimension(:,:)      , pointer :: elem2D_nodes
integer, dimension(:,:)      , pointer :: edges
integer, dimension(:,:)      , pointer :: edge_tri
integer, dimension(:,:)      , pointer :: elem_edges
real(kind=WP), dimension(:)  , pointer :: elem_area
real(kind=WP), dimension(:,:), pointer :: edge_dxdy, edge_cross_dxdy
real(kind=WP), dimension(:)  , pointer :: elem_cos, metric_factor
integer,       dimension(:,:), pointer :: elem_neighbors
integer,       dimension(:,:), pointer :: nod_in_elem2D
real(kind=WP), dimension(:,:), pointer :: x_corners, y_corners
integer,       dimension(:)  , pointer :: nod_in_elem2D_num
real(kind=WP), dimension(:)  , pointer :: depth
real(kind=WP), dimension(:,:), pointer :: gradient_vec
real(kind=WP), dimension(:,:), pointer :: gradient_sca
integer,       dimension(:)  , pointer :: bc_index_nod2D
real(kind=WP), dimension(:)  , pointer :: zbar, Z, elem_depth
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D, nlevels_nod2D_min
real(kind=WP), dimension(:,:), pointer :: area, area_inv, areasvol, areasvol_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n_ib
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
!real(kind=WP), dimension(:)  , pointer :: zbar_n
!real(kind=WP), dimension(:)  , pointer :: Z_n
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf
