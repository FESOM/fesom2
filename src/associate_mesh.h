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
real(kind=WP), dimension(:)  , pointer :: lump2d_north, lump2d_south
type(sparse_matrix)          , pointer :: ssh_stiff

nod2D              => mesh%nod2D              
elem2D             => mesh%elem2D             
edge2D             => mesh%edge2D             
edge2D_in          => mesh%edge2D_in          
ocean_area         => mesh%ocean_area         
nl                 => mesh%nl  

!!$coord_nod2D     => mesh%coord_nod2D        
!!$geo_coord_nod2D => mesh%geo_coord_nod2D    
!!$elem2D_nodes    => mesh%elem2D_nodes       
!!$edges              => mesh%edges              
!!$edge_tri           => mesh%edge_tri           
!!$elem_edges         => mesh%elem_edges         
!!$elem_area          => mesh%elem_area          
!!$edge_dxdy          => mesh%edge_dxdy          
!!$edge_cross_dxdy    => mesh%edge_cross_dxdy    
!!$elem_cos           => mesh%elem_cos           
!!$metric_factor      => mesh%metric_factor      
!!$elem_neighbors     => mesh%elem_neighbors     
!!$nod_in_elem2D      => mesh%nod_in_elem2D      
!!$x_corners          => mesh%x_corners          
!!$y_corners          => mesh%y_corners          
!!$nod_in_elem2D_num  => mesh%nod_in_elem2D_num  
!!$depth              => mesh%depth              
!!$gradient_vec       => mesh%gradient_vec       
!!$gradient_sca       => mesh%gradient_sca       
!!$bc_index_nod2D     => mesh%bc_index_nod2D     
!!$zbar               => mesh%zbar               
!!$Z                  => mesh%Z
!!$elem_depth         => mesh%elem_depth      
!!$nlevels            => mesh%nlevels            
!!$nlevels_nod2D      => mesh%nlevels_nod2D
!!$area               => mesh%area     
!!$area_inv           => mesh%area_inv     
!!$mesh_resolution    => mesh%mesh_resolution    
!!$ssh_stiff          => mesh%ssh_stiff          

if (allocated(mesh%coord_nod2D )) & 
     coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D)                  => mesh%coord_nod2D        
if (allocated(mesh%geo_coord_nod2D)) &
                   geo_coord_nod2D(1:2,1:myDim_nod2D+eDim_nod2D) => mesh%geo_coord_nod2D    
if (allocated(mesh%elem2D_nodes )) & 
     elem2D_nodes(1:3, 1:myDim_elem2D+eDim_elem2D+eXDim_elem2D) => mesh%elem2D_nodes
if (allocated(mesh%edges )) & 
     edges(1:2,1:myDim_edge2D+eDim_edge2D)                      => mesh%edges              
if (allocated(mesh%edge_tri )) & 
     edge_tri(1:2,1:myDim_edge2D+eDim_edge2D)                   => mesh%edge_tri           
if (allocated(mesh%elem_edges )) & 
     elem_edges(1:3,1:myDim_elem2D)                             => mesh%elem_edges         
if (allocated(mesh%elem_area )) & 
     elem_area(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)         => mesh%elem_area          
if (allocated(mesh%edge_dxdy )) & 
     edge_dxdy(1:2,1:myDim_edge2D+eDim_edge2D)                  => mesh%edge_dxdy          
if (allocated(mesh%edge_cross_dxdy )) & 
     edge_cross_dxdy(1:4,1:myDim_edge2D+eDim_edge2D)            => mesh%edge_cross_dxdy    
if (allocated(mesh%elem_cos )) & 
     elem_cos(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)          => mesh%elem_cos           
if (allocated(mesh%metric_factor )) & 
     metric_factor(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)     => mesh%metric_factor      
if (allocated(mesh%elem_neighbors )) & 
     elem_neighbors(1:3,1:myDim_elem2D)                         => mesh%elem_neighbors     
nod_in_elem2D      => mesh%nod_in_elem2D   ! (maxval(rmax),myDim_nod2D+eDim_nod2D)    
x_corners          => mesh%x_corners   ! (myDim_nod2D, maxval(rmax)) 
y_corners          => mesh%y_corners   ! (myDim_nod2D, maxval(rmax))       
if (allocated(mesh%nod_in_elem2D_num )) & 
     nod_in_elem2D_num(1:myDim_nod2D+eDim_nod2D)                => mesh%nod_in_elem2D_num  
if (allocated(mesh%depth )) & 
     depth(1:myDim_nod2D+eDim_nod2D)                            => mesh%depth              
if (allocated(mesh%gradient_vec )) & 
     gradient_vec(1:6,1:myDim_elem2D)                           => mesh%gradient_vec       
if (allocated(mesh%gradient_sca )) & 
     gradient_sca(1:6,1:myDim_elem2D)                           => mesh%gradient_sca       
if (allocated(mesh%bc_index_nod2D )) & 
     bc_index_nod2D(1:myDim_nod2D+eDim_nod2D)                   => mesh%bc_index_nod2D     
if (allocated(mesh%zbar )) & 
     zbar(1:mesh%nl)                                            => mesh%zbar               
if (allocated(mesh%Z )) & 
     Z(1:mesh%nl-1)                                             => mesh%Z
elem_depth         => mesh%elem_depth      ! never used, not even allocated
if (allocated(mesh%nlevels )) & 
     nlevels(1:myDim_elem2D+eDim_elem2D+eXDim_elem2D)           => mesh%nlevels            
if (allocated(mesh%nlevels_nod2D )) & 
     nlevels_nod2D(1:myDim_nod2D+eDim_nod2D)                    => mesh%nlevels_nod2D
if (allocated(mesh%area )) & 
     area(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                     => mesh%area     
if (allocated(mesh%area_inv )) & 
     area_inv(1:mesh%nl,1:myDim_nod2d+eDim_nod2D)                 => mesh%area_inv     
if (allocated(mesh%mesh_resolution )) & 
     mesh_resolution(1:myDim_nod2d+eDim_nod2D)                  => mesh%mesh_resolution    
ssh_stiff                                                  => mesh%ssh_stiff
if (allocated(mesh%lump2d_north )) & 
     lump2d_north(1:myDim_nod2d)                                => mesh%lump2d_north
if (allocated(mesh%lump2d_south )) & 
     lump2d_south(1:myDim_nod2d)                                => mesh%lump2d_south
