integer         , pointer :: nod2D,  myDim_nod2D, eDim_nod2D
integer         , pointer :: elem2D, myDim_elem2D, eDim_elem2D, eXDim_elem2D
integer         , pointer :: edge2D, myDim_edge2D, eDim_edge2D
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
!!$integer,       dimension(:)  , pointer :: cavity_lev_nod2D, cavity_lev_elem2D
integer,       dimension(:)  , pointer :: cavity_flag_n, cavity_flag_e
real(kind=WP), dimension(:)  , pointer :: cavity_depth
integer,       dimension(:)  , pointer :: ulevels, ulevels_nod2D, ulevels_nod2D_max
integer,       dimension(:)  , pointer :: nn_num
integer,       dimension(:,:), pointer :: nn_pos

real(kind=WP), dimension(:,:), pointer :: hnode
real(kind=WP), dimension(:,:), pointer :: hnode_new
real(kind=WP), dimension(:,:), pointer :: zbar_3d_n
real(kind=WP), dimension(:,:), pointer :: Z_3d_n
real(kind=WP), dimension(:,:), pointer :: helem
real(kind=WP), dimension(:)  , pointer :: bottom_elem_thickness
real(kind=WP), dimension(:)  , pointer :: bottom_node_thickness
real(kind=WP), dimension(:)  , pointer :: dhe
real(kind=WP), dimension(:)  , pointer :: hbar
real(kind=WP), dimension(:)  , pointer :: hbar_old
real(kind=WP), dimension(:)  , pointer :: zbar_n_bot
real(kind=WP), dimension(:)  , pointer :: zbar_e_bot
real(kind=WP), dimension(:)  , pointer :: zbar_n_srf
real(kind=WP), dimension(:)  , pointer :: zbar_e_srf

nod2D              => mesh%nod2D              
elem2D             => mesh%elem2D             
edge2D             => mesh%edge2D             
edge2D_in          => mesh%edge2D_in          
ocean_area         => mesh%ocean_area
ocean_areawithcav  => mesh%ocean_areawithcav         
nl                 => mesh%nl  
nn_size            => mesh%nn_size


myDim_nod2D  => p_partit%myDim_nod2D
eDim_nod2D   => p_partit%eDim_nod2D
myDim_elem2D => p_partit%myDim_elem2D
eDim_elem2D  => p_partit%eDim_elem2D
eXDim_elem2D => p_partit%eXDim_elem2D
myDim_edge2D => p_partit%myDim_edge2D
eDim_edge2D  => p_partit%eDim_edge2D


coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D        
geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)              => mesh%geo_coord_nod2D    
elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes
edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges              
edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri           
elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges         
elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area 
edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy          
edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy    
elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos           
metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor      
elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors     
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)    
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax)) 
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))       
nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num  
depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth              
gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec       
gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca       
bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D     
zbar(1:mesh%nl)                                            => mesh%zbar               
Z(1:mesh%nl-1)                                             => mesh%Z
elem_depth         => mesh%elem_depth      ! never used, not even allocated
nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels            
nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D
nlevels_nod2D_min(1:myDim_nod2D+eDim_nod2D)                => mesh%nlevels_nod2D_min
area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                   => mesh%area   
areasvol(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%areasvol
area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)               => mesh%area_inv  
areasvol_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)           => mesh%areasvol_inv
mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution    
ssh_stiff                                                  => mesh%ssh_stiff
lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north
lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south
cavity_flag_n(1:myDim_nod2D+eDim_nod2D)                    => mesh%cavity_flag_n
cavity_flag_e(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%cavity_flag_e  
!!$cavity_lev_nod2D(1:myDim_nod2D+eDim_nod2D)                 => mesh%cavity_lev_nod2D  
!!$cavity_lev_elem2D(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%cavity_lev_elem2D  
cavity_depth(1:myDim_nod2D+eDim_nod2D)                     => mesh%cavity_depth  
ulevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%ulevels            
ulevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%ulevels_nod2D
ulevels_nod2D_max(1:myDim_nod2D+eDim_nod2D)                => mesh%ulevels_nod2D_max
nn_num(1:myDim_nod2D)                                      => mesh%nn_num
nn_pos(1:mesh%nn_size, 1:myDim_nod2D)                      => mesh%nn_pos
hnode(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)               => mesh%hnode
hnode_new(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)           => mesh%hnode_new
zbar_3d_n(1:mesh%nl, 1:myDim_nod2D+eDim_nod2D)             => mesh%zbar_3d_n
Z_3d_n(1:mesh%nl-1, 1:myDim_nod2D+eDim_nod2D)              => mesh%Z_3d_n
helem(1:mesh%nl-1, 1:myDim_elem2D)                         => mesh%helem
bottom_elem_thickness(1:myDim_elem2D)                      => mesh%bottom_elem_thickness
bottom_node_thickness(1:myDim_nod2D+eDim_nod2D)            => mesh%bottom_node_thickness
dhe(1:myDim_elem2D)                                        => mesh%dhe
hbar(1:myDim_nod2D+eDim_nod2D)                             => mesh%hbar
hbar_old(1:myDim_nod2D+eDim_nod2D)                         => mesh%hbar_old
zbar_n_bot(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_bot
zbar_e_bot(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_bot
zbar_n_srf(1:myDim_nod2D+eDim_nod2D)                       => mesh%zbar_n_srf
zbar_e_srf(1:myDim_elem2D+eDim_elem2D)                     => mesh%zbar_e_srf


