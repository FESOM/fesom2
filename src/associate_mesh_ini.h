integer         , pointer :: nod2D 
integer         , pointer :: elem2D   
integer         , pointer :: edge2D   
integer         , pointer :: edge2D_in
real(kind=WP)   , pointer :: ocean_area
integer         , pointer :: nl
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
integer,       dimension(:)  , pointer :: nlevels, nlevels_nod2D
real(kind=WP), dimension(:,:), pointer :: area, area_inv
real(kind=WP), dimension(:)  , pointer :: mesh_resolution
integer,       dimension(:)  , pointer :: cavity_flag, ulevels_nod2D, ulevels
real(kind=WP), dimension(:)  , pointer :: cavity_depth
type(sparse_matrix)          , pointer :: ssh_stiff

nod2D              => mesh%nod2D              
elem2D             => mesh%elem2D             
edge2D             => mesh%edge2D             
edge2D_in          => mesh%edge2D_in          
ocean_area         => mesh%ocean_area         
nl                 => mesh%nl  

coord_nod2D     => mesh%coord_nod2D        
geo_coord_nod2D => mesh%geo_coord_nod2D    
elem2D_nodes    => mesh%elem2D_nodes       
edges              => mesh%edges              
edge_tri           => mesh%edge_tri           
elem_edges         => mesh%elem_edges         
elem_area          => mesh%elem_area          
edge_dxdy          => mesh%edge_dxdy          
edge_cross_dxdy    => mesh%edge_cross_dxdy    
elem_cos           => mesh%elem_cos           
metric_factor      => mesh%metric_factor      
elem_neighbors     => mesh%elem_neighbors     
nod_in_elem2D      => mesh%nod_in_elem2D      
x_corners          => mesh%x_corners          
y_corners          => mesh%y_corners          
nod_in_elem2D_num  => mesh%nod_in_elem2D_num  
depth              => mesh%depth              
gradient_vec       => mesh%gradient_vec       
gradient_sca       => mesh%gradient_sca       
bc_index_nod2D     => mesh%bc_index_nod2D     
zbar               => mesh%zbar               
Z                  => mesh%Z
elem_depth         => mesh%elem_depth      
nlevels            => mesh%nlevels            
nlevels_nod2D      => mesh%nlevels_nod2D
area               => mesh%area     
area_inv           => mesh%area_inv     
mesh_resolution    => mesh%mesh_resolution    
ssh_stiff          => mesh%ssh_stiff        
!!$cavity_flag_n      => mesh%cavity_flag_n  
!!$cavity_flag_e      => mesh%cavity_flag_e
ulevels_nod2D      => mesh%ulevels_nod2D  
ulevels            => mesh%ulevels  
cavity_depth       => mesh%cavity_depth  
