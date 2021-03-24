# Patrick Scholz, 14.12.2017
import time
import numpy as np
import pandas as pa
import math
from numpy.linalg import inv
from set_inputarray import *
import copy as cp

#+___INITIALISE/LOAD FESOM2.0 MESH CLASS IN MAIN PROGRAMM______________________+
#| IMPORTANT!!!:                                                               |                                         
#| only when mesh is initialised with this routine, the main programm is able  |
#| to recognize if mesh object already exist or not and if it exist to use it  |
#| and not to load it again !!!                                                |                    |
#+_____________________________________________________________________________+
def fesom_init_mesh(inputarray):
#     global mesh
    mesh = fesom_mesh(inputarray)
    return(mesh)
    
    
    
#+_____________________________________________________________________________+
#|                                                                             |
#|                        *** FESOM2.0 MESH CLASS ***                          |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_mesh:
    
    #____meshinfo, path, id...____________________
    metadata, path, id         = '' ,'', ''
    which_obj                = 'mesh'
    
    #____mesh euler angles_______________________
    alpha, beta, gamma        = 0, 0, 0
    focus                     = 0
    focus_old                 = 0
    rotate_grid               = True
    
    #____mesh vertical info______________________
    nlev, zlev, zmid        = 0, [], []
    
    #____mesh nodes & elem size__________________
    n2dn, n2dna                = 0, 0
    n2de, n2dna                = 0, 0
    
    #____mesh nodes info_________________________
    # original coordinates
    nodes_2d_x , nodes_2d_y   = [], []
    nodes_2d_z , nodes_2d_iz  = [], []
    nodes_2d_i                = []
    
    # coordinates rotated to geo
    nodes_2d_xg, nodes_2d_yg  = [], []
    nodes_2d_zg, nodes_2d_izg = [], []
    nodes_2d_ig               = []
    
    # cavity info 
    use_cavity                = False
    nodes_2d_cg, nodes_2d_icg = [], []
    nodes_2d_c,  nodes_2d_ic  = [], []
    elem0_2d_ic               = []
    
    # rotated coordiantes
    nodes_2d_xr, nodes_2d_yr  = [], []
        
    pbndn_2d_i                = []
    
    #____mesh elem info__________________________
    elem0_2d_i, elem_2d_i     = [], []
    elem0_2d_iz,pbndtri_2d_i  = [], []
    
    abnormtri_2d_i            = []
    
    #____fesom lsmask polygon____________________
    polygon_xy                = []
    
    #____triangle area in km^2___________________
    elem0_2d_area             = []
    elem_2d_resol             = []
    
    #____triangle area in km^2___________________
    nodes_2d_resol            = []
    nodes_2d_area             = []
    
    #+___INIT FESOM2.0 MESH OBJECT_____________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def __init__(self,inputarray):
        
        #_______________________________________________________________________
        # define euler angles
        self.id          = inputarray['mesh_id'   ]
        self.path        = inputarray['mesh_dir'  ]
        self.alpha       = inputarray['mesh_alpha']
        self.beta        = inputarray['mesh_beta' ]
        self.gamma       = inputarray['mesh_gamma']
        self.focus       = inputarray['mesh_focus']
        self.focus_old   = self.focus    
        self.rotate_grid = inputarray['mesh_rotate'	 ]
        
        print('')
        print('___LOAD FESOM MESH {}_________________________________________'.format(self.id))
        
        #_______________________________________________________________________
        # load mesh from file
        self.fesom_load_mesh()
        
        if (inputarray['use_cavity']==True):
            self.use_cavity = True
            self.fesom_load_cavity()
        
        #_______________________________________________________________________
        # rotate mesh from rot to geo coordinates
        if (inputarray['mesh_rotate'     ]==True):
            self.fesom_grid_rot_r2g('r2g')
        else:
            self.nodes_2d_xg  = self.nodes_2d_x
            self.nodes_2d_yg  = self.nodes_2d_y
            self.nodes_2d_zg  = self.nodes_2d_z
            self.nodes_2d_izg = self.nodes_2d_iz
            self.nodes_2d_ig  = self.nodes_2d_i
            if (inputarray['use_cavity']==True):
                self.nodes_2d_cg  = self.nodes_2d_c
                self.nodes_2d_icg = self.nodes_2d_ic
        #_______________________________________________________________________
        # change grid focus from -180...180 to e.g. 0...360
        if (inputarray['mesh_focus']!=0):
            self.fesom_grid_rot_r2g('focus')
        
        #_______________________________________________________________________
        # remove+augment periodic boundary
        self.fesom_remove_pbnd()
        
        # calculate fesom land mask interactivly
        #self.fesom_calc_landmask()
        
        #self.fesom_grid_rot_g2r()
        
        # calc triangle area in km^2 
        #self.fesom_calc_triarea()
        # calc triangle resol in km interpolated to nodes
        #mesh = fesom_calc_triresol(mesh)
        
        self.metadata = inputarray
        
        print('_______________________________________________________________')
    
    #+___LOAD FESOM MESH FROM FILE_____________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_load_mesh(self):
        
        print(' --> read grid files')
        
        #____load 2d node matrix________________________________________________
        print('     > nod2d.out',end='')
        file_content    = pa.read_csv(self.path+'nod2d.out', delim_whitespace=True, skiprows=1, \
                                      names=['node_number','x','y','flag'] )
        self.nodes_2d_x = file_content.x.values.astype('float32')
        self.nodes_2d_y = file_content.y.values.astype('float32')
        self.nodes_2d_i = file_content.flag.values.astype('uint16')   
        self.n2dn       = len(self.nodes_2d_x)
        print('  : #2dn={:d}'.format(self.n2dn))
        
        #____load 2d element matrix_____________________________________________
        print('     > elem2d.out',end='')
        file_content    = pa.read_csv(self.path+'elem2d.out', delim_whitespace=True, skiprows=1, \
                                       names=['1st_node_in_elem','2nd_node_in_elem','3rd_node_in_elem'])
        self.elem0_2d_i = file_content.values.astype('uint32') - 1
        self.n2de       = np.shape(self.elem0_2d_i)[0]
        print(' : #2de={:d}'.format(self.n2de))
        
        #____load 3d nodes alligned under 2d nodes______________________________
        print('     > aux3d.out') 
        with open(self.path+'aux3d.out') as f:
            self.nlev = int(next(f))
            self.zlev = np.array([next(f).rstrip() for x in range(self.nlev)]).astype(float)
        self.zmid     = (self.zlev[:-1]+self.zlev[1:])/2.
        
        #____load number of levels at each node_________________________________
        print('     > nlvls.out') 
        file_content    = pa.read_csv(self.path+'nlvls.out', delim_whitespace=True, skiprows=0, \
                                       names=['numb_of_lev'])
        self.nodes_2d_iz= file_content.values.astype('uint16') - 1
        self.nodes_2d_iz= self.nodes_2d_iz.squeeze()
        self.nodes_2d_z = np.float32(self.zlev[self.nodes_2d_iz])
        
        #____load number of levels at each elem_________________________________
        print('     > elvls.out') 
        file_content    = pa.read_csv(self.path+'elvls.out', delim_whitespace=True, skiprows=0, \
                                       names=['numb_of_lev'])
        self.elem0_2d_iz= file_content.values.astype('uint16') - 1
        self.elem0_2d_iz= self.elem0_2d_iz.squeeze()
        
    #+___LOAD FESOM CAVITY INFORMATION FROM FILE_______________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_load_cavity(self):
        
        print(' --> read cavity files')
        
        #____load number of levels at each node_________________________________
        print('     > cavity_nlvls.out') 
        file_content    = pa.read_csv(self.path+'cavity_nlvls.out', delim_whitespace=True, skiprows=0, \
                                       names=['numb_of_lev'])
        self.nodes_2d_ic= file_content.values.astype('uint16') - 1
        self.nodes_2d_ic= self.nodes_2d_ic.squeeze()
        self.nodes_2d_c = np.float32(self.zlev[self.nodes_2d_ic])
        
        #____load number of levels at each elem_________________________________
        print('     > cavity_elvls.out') 
        file_content    = pa.read_csv(self.path+'cavity_elvls.out', delim_whitespace=True, skiprows=0, \
                                       names=['numb_of_lev'])
        self.elem0_2d_ic= file_content.values.astype('uint16') - 1
        self.elem0_2d_ic= self.elem0_2d_ic.squeeze()
        
        
    #+___ROTATE GRID ROT-->GEO_________________________________________________+
    #| input : r2g (default) change coordinate from rotated-->geo              |
    #|         g2r           change coordinate from geo-->rotated              |
    #|         focus         change longitudinal focus                         |
    #+_________________________________________________________________________+
    def fesom_grid_rot_r2g(self,str_mode='r2g'):
        
        if (str_mode == 'focus'):
            print(' --> change mesh focus');
            if self.rotate_grid == True:
                alpha = self.alpha-self.focus
                beta  = self.beta
                gamma = self.gamma
            else:    
                alpha = -self.focus
                beta  = 0.
                gamma = 0.
        else:
            if (str_mode == 'r2g'):
                print(' --> rotate mesh rot2geo')
            else:
                print(' --> rotate mesh geo2rot')
                
            alpha = self.alpha
            beta  = self.beta
            gamma = self.gamma
        #_______________________________________________________________________
        # build euler rotation matrix
        rotate_matrix = np.zeros((3,3))
        rotate_matrix[0,0] =( math.cos(math.radians(gamma)) * 
                            math.cos(math.radians(alpha)) -
                            math.sin(math.radians(gamma)) *
                            math.cos(math.radians(beta )) *
                            math.sin(math.radians(alpha)))
        
        rotate_matrix[0,1] =( math.cos(math.radians(gamma)) *
                            math.sin(math.radians(alpha)) +
                            math.sin(math.radians(gamma)) *
                            math.cos(math.radians(beta )) *
                            math.cos(math.radians(alpha)))
        
        rotate_matrix[0,2] =( math.sin(math.radians(gamma)) *
                            math.sin(math.radians(beta )))
        
        rotate_matrix[1,0] =(-math.sin(math.radians(gamma)) *
                            math.cos(math.radians(alpha)) -
                            math.cos(math.radians(gamma)) *
                            math.cos(math.radians(beta )) *
                            math.sin(math.radians(alpha)))
        
        rotate_matrix[1,1] =(-math.sin(math.radians(gamma)) *
                            math.sin(math.radians(alpha)) +
                            math.cos(math.radians(gamma)) *
                            math.cos(math.radians(beta )) *
                            math.cos(math.radians(alpha)))
        
        rotate_matrix[1,2] =( math.cos(math.radians(gamma)) *
                            math.sin(math.radians(beta )))
        
        rotate_matrix[2,0] =( math.sin(math.radians(beta )) *
                            math.sin(math.radians(alpha)))
        
        rotate_matrix[2,1] =(-math.sin(math.radians(beta )) *
                            math.cos(math.radians(alpha)))
        
        rotate_matrix[2,2] =( math.cos(math.radians(beta )))
        
        #_______________________________________________________________________
        # make inverse of rotation matrix
        if (str_mode == 'r2g') or (str_mode == 'focus'):
            from numpy.linalg import inv
            rotate_matrix=inv(rotate_matrix);
            
        #____3D_________________________________________________________________
        # calculate Cartesian coordinates
        #if (str_mode == 'focus') or (str_mode == 'g2r'):
        if (str_mode == 'g2r'):
            xr=np.cos(np.radians(self.nodes_2d_yg)) * np.cos(np.radians(self.nodes_2d_xg));
            yr=np.cos(np.radians(self.nodes_2d_yg)) * np.sin(np.radians(self.nodes_2d_xg));
            zr=np.sin(np.radians(self.nodes_2d_yg));
        else:
            xr=np.cos(np.radians(self.nodes_2d_y)) * np.cos(np.radians(self.nodes_2d_x));
            yr=np.cos(np.radians(self.nodes_2d_y)) * np.sin(np.radians(self.nodes_2d_x));
            zr=np.sin(np.radians(self.nodes_2d_y));
        
        #_______________________________________________________________________
        # rotate to geographical cartesian coordinates:
        xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr;
        yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr;
        zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr;
        
        ##______________________________________________________________________
        #self.nodes_2d_yg = np.degrees(np.arcsin(zg));
        #self.nodes_2d_xg = np.degrees(np.arctan2(yg,xg));
        
        #_______________________________________________________________________
        if (str_mode == 'focus') or (str_mode == 'r2g'):
            self.nodes_2d_yg=np.degrees(np.arcsin(zg));
            self.nodes_2d_xg=np.degrees(np.arctan2(yg,xg));
            if (str_mode == 'focus'):
                self.nodes_2d_xg = self.nodes_2d_xg+self.focus
                
            self.nodes_2d_zg = self.nodes_2d_z
            self.nodes_2d_izg = self.nodes_2d_iz
            self.nodes_2d_ig = self.nodes_2d_i
            if (self.use_cavity == True):
                self.nodes_2d_cg = self.nodes_2d_c
                self.nodes_2d_icg = self.nodes_2d_ic
                
            if (str_mode == 'focus') and self.focus!=self.focus_old:
                self.fesom_remove_pbnd()
                self.focus_old = self.focus
                
        elif (str_mode == 'g2r'):
            self.nodes_2d_yr=np.degrees(np.arcsin(zg));
            self.nodes_2d_xr=np.degrees(np.arctan2(yg,xg));
            #if (str_mode == 'focus'):
                #self.nodes_2d_xr = self.nodes_2d_xr+self.focus
        #if (str_mode == 'focus'):
            #self.fesom_remove_pbnd()
            #self.fesom_calc_landmask()
    
    #+___ROTATE GRID GEO-->ROT_________________________________________________+
    #| input : g2r           change coordinate from geo-->rotated              |
    #+_________________________________________________________________________+
    def fesom_grid_rot_g2r(self):
        self.fesom_grid_rot_r2g('g2r')
    
    
    #+___SEARCH AND REMOVE PERIODIC BOUNDARIES_________________________________+
    #| identify the periodic boundary (nodind2d ==5000) and deletes it. Augment| 
    #| the left and right part                                                 |
    #+_________________________________________________________________________+
    def fesom_remove_pbnd(self):
        
        print(' --> remove cyclic boundary')
        #_______________________________________________________________________
        # find out 1st which element contribute to periodic boundary and 2nd
        # which nodes are involed in periodic boundary
        dx = self.nodes_2d_xg[self.elem0_2d_i].max(axis=1)-self.nodes_2d_xg[self.elem0_2d_i].min(axis=1)
        ind_elem_cyclbnd = np.where(dx >= 180)[0]
        ind_node_cyclbnd = np.array(self.elem0_2d_i[ind_elem_cyclbnd,:]).flatten()
        ind_node_cyclbnd = np.unique(ind_node_cyclbnd)
        del dx
        
        #_______________________________________________________________________
        # find out node indices that contribute to the left side of the periodic 
        # boundary (pbndn_l_2d_i) and to the right side (pbndn_r_2d_i)
        xmin, xmax = self.nodes_2d_xg.min(), self.nodes_2d_xg.max()
        aux_l_i    = np.where(self.nodes_2d_xg[ind_node_cyclbnd]<(xmin+(xmax-xmin)/2) )[0]
        aux_r_i    = np.where(self.nodes_2d_xg[ind_node_cyclbnd]>(xmin+(xmax-xmin)/2) )[0]
        pbndn_l_2d_i, pbndn_r_2d_i = ind_node_cyclbnd[aux_l_i], ind_node_cyclbnd[aux_r_i]
        npbnd_r_2d  ,npbnd_l_2d    = pbndn_r_2d_i.size, pbndn_l_2d_i.size
        self.pbndn_2d_i = np.concatenate((pbndn_r_2d_i,pbndn_l_2d_i))
        del aux_l_i, aux_r_i, ind_node_cyclbnd
        
        #_______________________________________________________________________
        # calculate augmentation positions for new left and right periodic boundaries
        aux_pos    = np.zeros(self.n2dn,dtype='uint32')
        aux_nnodesr= np.linspace(self.n2dn,self.n2dn+npbnd_r_2d-1,npbnd_r_2d,dtype='uint32')
        aux_nnodesl= np.linspace(self.n2dn+npbnd_r_2d,self.n2dn+npbnd_r_2d+npbnd_l_2d-1,npbnd_l_2d,dtype='uint32')
        aux_pos[pbndn_r_2d_i],aux_pos[pbndn_l_2d_i]=aux_nnodesr,aux_nnodesl
        del aux_nnodesl, aux_nnodesr, pbndn_l_2d_i, pbndn_r_2d_i
        
        #_______________________________________________________________________
        # Augment the nodes on the right and left side 
        xmin, xmax= np.floor(xmin),np.ceil(xmax)
        self.nodes_2d_xg = np.concatenate((self.nodes_2d_xg,
                                            np.zeros(npbnd_r_2d)+xmin,
                                            np.zeros(npbnd_l_2d)+xmax    ))
        self.nodes_2d_yg  = np.concatenate((self.nodes_2d_yg,self.nodes_2d_yg[self.pbndn_2d_i]))
        self.nodes_2d_zg  = np.concatenate((self.nodes_2d_z,self.nodes_2d_z[self.pbndn_2d_i  ]))
        self.nodes_2d_ig  = np.concatenate((self.nodes_2d_i,self.nodes_2d_i[self.pbndn_2d_i  ]))
        self.nodes_2d_izg = np.concatenate((self.nodes_2d_iz,self.nodes_2d_iz[self.pbndn_2d_i]))
        if (self.use_cavity==True):
            self.nodes_2d_cg  = np.concatenate((self.nodes_2d_c,self.nodes_2d_c[self.pbndn_2d_i  ]))
            self.nodes_2d_icg = np.concatenate((self.nodes_2d_ic,self.nodes_2d_ic[self.pbndn_2d_i]))
            
        self.n2dna=self.n2dn+self.pbndn_2d_i.size # new (augmented) n2d
        
        #____loop over all periodic bnd. segments___________________________________
        copy_elem = np.zeros(self.n2de,dtype='bool')
        self.elem_2d_i = np.copy(self.elem0_2d_i);
        
        #___________________________________________________________________________
        # (ii.a) 2d elements:
        # List all triangles that touch the cyclic boundary segments
        #_______________________________________________________________________
        elem_pbnd_l = np.copy(self.elem0_2d_i[ind_elem_cyclbnd,:])
        elem_pbnd_r = np.copy(elem_pbnd_l)
        for ei in range(0,ind_elem_cyclbnd.size):
            # node indices of periodic boundary triangle
            tri  = np.array(self.elem0_2d_i[ind_elem_cyclbnd[ei],:]).squeeze()
            
            # which triangle points belong to left periodic bnde or right periodic
            # boundary
            idx_l = np.where(self.nodes_2d_xg[tri].squeeze()<xmin+(xmax-xmin)/2)[0]
            idx_r = np.where(self.nodes_2d_xg[tri].squeeze()>xmin+(xmax-xmin)/2)[0]
            
            # change indices to left and right augmented boundary points
            elem_pbnd_l[ei,idx_r]=aux_pos[tri[idx_r]]
            elem_pbnd_r[ei,idx_l]=aux_pos[tri[idx_l]]
        del idx_l, idx_r, tri, aux_pos
        
        #_______________________________________________________________________
        # change existing periodic boundary triangles in elem_2d_i to augmented 
        # left boundary triangles
        #self.elem_2d_i = np.copy(self.elem0_2d_i)
        self.elem_2d_i[ind_elem_cyclbnd,:] = elem_pbnd_l
        
        # add additional augmented right periodic boundary triangles
        self.elem_2d_i = np.concatenate((self.elem_2d_i,elem_pbnd_r))
        self.n2dea     = self.elem_2d_i.shape[0]
        
        # indices in elem0_2d_i of periodic boudaries
        self.pbndtri_2d_i = ind_elem_cyclbnd
        del elem_pbnd_l, elem_pbnd_r, ind_elem_cyclbnd
        
        #_______________________________________________________________________
        
        
        """
        #_______________________________________________________________________
        # build land boundary edge matrix
        edge    = np.concatenate((self.elem_2d_i[:,[0,1]], \
                                self.elem_2d_i[:,[0,2]], \
                                self.elem_2d_i[:,[1,2]]))
        edge    = np.sort(edge,axis=1) 
        
        # python  sortrows algorythm --> matlab equivalent
        edge    = edge.tolist()
        edge.sort()
        edge    = np.array(edge)
        idx     = np.diff(edge,axis=0)==0
        idx     = np.all(idx,axis=1)
        idx     = np.logical_or(np.concatenate((idx,np.array([False]))),\
                                np.concatenate((np.array([False]),idx)))
        bnde    = edge[idx==False,:]
        
        #_______________________________________________________________________
        # kill periodic boundary edges--> periodic boundary edges have identical 
        # lon/lat points
        edge_x  = self.nodes_2d_xg[bnde]
        edge_y  = self.nodes_2d_yg[bnde]
        aux     = np.logical_or(edge_x==np.min(self.nodes_2d_xg),\
                                edge_x==np.max(self.nodes_2d_xg))
        aux     = np.logical_or(aux,edge_y>89)
        idx     = np.all(aux,axis=1)
        del aux
        del edge_x 
        del edge_y 
            
        bnde    = bnde[idx==False,:]
        del idx
        self.nodes_2d_i = np.zeros((self.n2dna),dtype='uint8')
        self.nodes_2d_i[bnde.reshape(bnde.size)]=1
        del bnde; del edge 
        """
    
    #+___CALCULATE MESH RESOLUTION_____________________________________________+
    #| calculate mean length of the three triangle sides                       |
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_calc_triresol(self):
        
        #_______________________________________________________________________
        print(' --> calc. triangle resolution in km')
        
        #_______________________________________________________________________
        nodes_2d_x     = np.concatenate((self.nodes_2d_xg[:self.n2dn] , \
                                      self.nodes_2d_xg[self.pbndn_2d_i]))
        nodes_2d_y     = np.concatenate((self.nodes_2d_yg[:self.n2dn] , \
                                      self.nodes_2d_yg[self.pbndn_2d_i]))
        elem_2d_x      = nodes_2d_x[self.elem_2d_i]
        elem_2d_y      = nodes_2d_y[self.elem_2d_i]
        elem_2d_xy     = np.array([elem_2d_x,elem_2d_y])
        del nodes_2d_x
        del nodes_2d_y
        del elem_2d_x
        
        #_______________________________________________________________________
        # calc jacobi matrix for all triangles 
        # | dx_12 dy_12 |
        # | dx_23 dy_23 |
        # | dx_31 dy_31 |_i , i=1....n2dea
        jacobian     = elem_2d_xy[:,:,1]-elem_2d_xy[:,:,0]
        jacobian     = np.array([jacobian,elem_2d_xy[:,:,2]-elem_2d_xy[:,:,1],\
                                          elem_2d_xy[:,:,0]-elem_2d_xy[:,:,2] ])
        
        # account for triangles with periodic bounaries
        for ii in range(3):
            idx = np.where(jacobian[ii,0,:]>180); 
            jacobian[ii,0,idx] = jacobian[ii,0,idx]-360;
            idx = np.where(jacobian[ii,0,:]<-180); 
            jacobian[ii,0,idx] = jacobian[ii,0,idx]+360;
            del idx
        
        # calc from geocoord to cartesian coord
        R_earth       = 12735/2;
        jacobian      = jacobian*R_earth*2*np.pi/360
        cos_theta     = np.cos(np.radians(elem_2d_y)).mean(axis=1)
        del elem_2d_y
        for ii in range(3):    
            jacobian[ii,0,:] = jacobian[ii,0,:]*cos_theta;
        del cos_theta
        
        #_______________________________________________________________________
        # calc vector length dr = sqrt(dx^2+dy^2)
        jacobian     = np.power(jacobian,2);
        jacobian     = np.sqrt(jacobian.sum(axis=1));
        jacobian     = jacobian.transpose()
        
        #_______________________________________________________________________
        # mean resolutiuon per element
        self.elem_2d_resol = jacobian.mean(axis=1)
        #self.nodes_2d_resol= self.fesom_interp_e2n(self.elem_2d_resol)
    
    
    #+___CALCULATE MESH RESOLUTION_____________________________________________+
    #| calculate mean length of the three triangle sides                       |
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_calc_noderesol(self):
        
        #_______________________________________________________________________
        print(' --> calc. node resolution in km')
        if len(self.elem_2d_resol)==0: self.fesom_calc_triresol()
        self.nodes_2d_resol= self.fesom_interp_e2n(self.elem_2d_resol)
        
    
    #___CALCULATE TRIANGLE AREA("volume")_______________________________________
    # calculate TRIANGLE AREA IN M^2
    #
    #___________________________________________________________________________
    def fesom_calc_triarea(self):
        
        #_______________________________________________________________________
        print(' --> calc. triangle area m^2',end='')
        t1 = time.time()
        #_______________________________________________________________________
        nodes_2d_x     = np.concatenate((self.nodes_2d_xg[:self.n2dn] , \
                                      self.nodes_2d_xg[self.pbndn_2d_i]))
        nodes_2d_y     = np.concatenate((self.nodes_2d_yg[:self.n2dn] , \
                                      self.nodes_2d_yg[self.pbndn_2d_i]))
        elem_2d_x      = nodes_2d_x[self.elem0_2d_i]
        elem_2d_y      = nodes_2d_y[self.elem0_2d_i]
        elem_2d_xy     = np.array([elem_2d_x,elem_2d_y])
        del nodes_2d_x, nodes_2d_y, elem_2d_x
        
        #_______________________________________________________________________
        # calc jacobi matrix for all triangles 
        # | dx_12 dy_12 |
        # | dx_13 dy_13 |_i , i=1....n2dea
        jacobian     = elem_2d_xy[:,:,1]-elem_2d_xy[:,:,0]
        jacobian     = np.array([jacobian,elem_2d_xy[:,:,2]-elem_2d_xy[:,:,0] ])
        
        # account for triangles with periodic bounaries
        for ii in range(2):
            idx = np.where(jacobian[ii,0,:]>180); 
            jacobian[ii,0,idx] = jacobian[ii,0,idx]-360;
            idx = np.where(jacobian[ii,0,:]<-180); 
            jacobian[ii,0,idx] = jacobian[ii,0,idx]+360;
            del idx
        
        # calc from geocoord to cartesian coord
        R_earth        =12735/2;
        jacobian     = jacobian*R_earth*2*np.pi/360
        cos_theta     = np.cos(np.radians(elem_2d_y)).mean(axis=1)
        del elem_2d_y
        for ii in range(2):    
            jacobian[ii,0,:] = jacobian[ii,0,:]*cos_theta;
        del cos_theta
        
        #_______________________________________________________________________
        # calc triangle area through vector product
        self.elem0_2d_area = 0.5*np.abs(jacobian[0,0,:]*jacobian[1,1,:]-\
                                       jacobian[0,1,:]*jacobian[1,0,:])
        #self.elem0_2d_area = np.concatenate((self.elem_2d_area,\
                                            #self.elem_2d_area[self.pbndtri_2d_i]))
        del jacobian
        
        #_______________________________________________________________________
        t2 = time.time()
        print(' >> time:{:.3f} s'.format(t2-t1))
        
    
    #___CALCULATE TRIANGLE AREA("volume")_______________________________________
    # calculate TRIANGLE AREA IN M^2
    #
    #___________________________________________________________________________
    def fesom_calc_nodearea(self):
        #_______________________________________________________________________
        print(' --> calc. triangle area on nodes m^2',end='')
        t1 = time.time()
        #_______________________________________________________________________
        if len(self.elem0_2d_area)==0: self.fesom_calc_triarea()
        #_______________________________________________________________________
        # calc node cluster area 
        self.nodes_2d_area=np.zeros(self.n2dna)
        aux_area = np.concatenate((self.elem0_2d_area, self.elem0_2d_area[self.pbndtri_2d_i]))
        aux_area = aux_area.reshape(aux_area.size,1)
        aux_area = np.concatenate((aux_area,aux_area,aux_area),axis=1).flatten()
        count    = 0
        # single loop over self.elem0_2d_i.flat is ~4 times faster than douple loop 
        # over for i in range(3): ,for j in range(self.n2de):
        for idx in self.elem0_2d_i.flat:
            self.nodes_2d_area[idx]=self.nodes_2d_area[idx]+aux_area[count]
            count=count+1 # count triangle index for aux_area[count] --> aux_area =[n2de*3,]
        del aux_area, count
        self.nodes_2d_area=self.nodes_2d_area/3.
        self.nodes_2d_area[self.n2dn:self.n2dna] = self.nodes_2d_area[self.pbndn_2d_i]
        #_______________________________________________________________________
        t2 = time.time()
        print(' >> time:{:.3f} s'.format(t2-t1))
        
        
    #+___CALCULATE INTERPOLATION MATRIX FROM ELEMENTS TO NODE POINTS RESOLUTION+
    #| interpolate from triangle values to node values                         |
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_interp_e2n(self,data_e):
        if len(self.elem0_2d_area)==0: self.fesom_calc_triarea()
        #if len(self.nodes_2d_area)==0: self.fesom_calc_nodearea()
        print(' --> calc. compute e2n interpolation',end='')
        t1 = time.time()
        #_______________________________________________________________________
        # case of one dimensional elemental data
        if data_e.ndim==1:
            if data_e.size==self.n2de:
                data_e=data_e*self.elem0_2d_area
                # elemental area of non bottom element
                data_e_area = self.elem0_2d_area*np.invert(np.isnan(data_e))
                data_e[np.isnan(data_e)]=0
                
            elif data_e.size==self.n2dea:
                data_e=data_e*np.concatenate((self.elem0_2d_area,self.elem0_2d_area[self.pbndtri_2d_i]))
                # elemental area of non bottom element
                data_e_area = np.concatenate((self.elem0_2d_area,self.elem0_2d_area[self.pbndtri_2d_i]))*np.invert(np.isnan(data_e))
                data_e[np.isnan(data_e)]=0
            
            data_n     =np.zeros((self.n2dna,))
            data_n_area=np.zeros((self.n2dna,))
            for ii in range(self.n2de):     
                data_n[     self.elem0_2d_i[ii,:]] = data_n[     self.elem0_2d_i[ii,:]] + np.array([1,1,1])*data_e[ii]
                data_n_area[self.elem0_2d_i[ii,:]] = data_n_area[self.elem0_2d_i[ii,:]] + np.array([1,1,1])*data_e_area[ii] 
            data_n[        0:self.n2dn ] = data_n[0:self.n2dn]/data_n_area[0:self.n2dn]
            data_n[self.n2dn:self.n2dna] = data_n[self.pbndn_2d_i]
        
        #_______________________________________________________________________
        # case of two dimensional elemental data
        elif data_e.ndim==2:
            nd1 = data_e.shape[0]
            nd2 = data_e.shape[1]
            
            if data_e.shape[0]==self.n2de:
                data_e      = data_e*np.matlib.repmat(self.elem0_2d_area,nd2,1).transpose()
                data_e_area = np.matlib.repmat(self.elem0_2d_area,nd2,1).transpose()*np.invert(data_e==0)
                
            elif data_e.shape[0]==self.n2dea:
                area_di     = np.concatenate((self.elem0_2d_area,self.elem0_2d_area[self.pbndtri_2d_i]))
                data_e      = data_e*np.matlib.repmat(area_di,nd2,1).transpose()
                data_e_area = np.matlib.repmat(area_di,nd2,1).transpose()*np.invert(data_e==0)
                del area_di
            
            data_n      = np.zeros((self.n2dna,nd2))
            data_n_area = np.zeros((self.n2dna,nd2))
            for ii in range(self.n2de):     
                data_n[     self.elem0_2d_i[ii,:], :] = data_n[     self.elem0_2d_i[ii,:], :] + np.matlib.repmat(     data_e[ii,:], 3, 1)
                data_n_area[self.elem0_2d_i[ii,:], :] = data_n_area[self.elem0_2d_i[ii,:], :] + np.matlib.repmat(data_e_area[ii,:], 3, 1)  
            #print(np.matlib.repmat(self.nodes_2d_area[0:self.n2dn]/3.,nd2,1).transpose())    
            data_n[        0:self.n2dn ,:] = data_n[0:self.n2dn,:]/data_n_area[0:self.n2dn,:]
            data_n[self.n2dn:self.n2dna,:] = data_n[self.pbndn_2d_i,:]
            
        t2 = time.time()
        print(' >> time:{:.3f} s'.format(t2-t1))
        
        return data_n
    
    
    #+___CALCULATE LAND MASK CONOURLINE________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def fesom_calc_landmask(self):

        #np.set_printoptions(edgeitems=10)
        print(' --> calc landmask contourline')
        #_______________________________________________________________________
        # build land boundary edge matrix
        
        t1 = time.time()
        edge    = np.concatenate((self.elem_2d_i[:,[0,1]], \
                                                self.elem_2d_i[:,[0,2]], \
                                                self.elem_2d_i[:,[1,2]]),axis=0)
        edge    = np.sort(edge,axis=1) 
        
        # python  sortrows algorythm --> matlab equivalent
        # twice as fast as list sorting
        #sortidx = np.lexsort((edge[:,0],edge[:,1]))
        #edge    = edge[sortidx,:].squeeze()
        #edge    = np.array(edge)
        
        ## python  sortrows algorythm --> matlab equivalent
        edge    = edge.tolist()
        edge.sort()
        edge    = np.array(edge)
        
        idx     = np.diff(edge,axis=0)==0
        idx     = np.all(idx,axis=1)
        idx     = np.logical_or(np.concatenate((idx,np.array([False]))),\
                                np.concatenate((np.array([False]),idx)))
        
        # all edges that belong to boundary
        bnde    = edge[idx==False,:]
        del edge, idx
        
        #_______________________________________________________________________
        # kill periodic boundary edges--> periodic boundary edges have identical 
        # lon/lat points
        edge_x  = self.nodes_2d_xg[bnde]
        edge_y  = self.nodes_2d_yg[bnde]
        
        aux     = np.logical_or(edge_x==self.nodes_2d_xg.min(),\
                                edge_x==self.nodes_2d_xg.max())
        aux     = np.logical_or(aux,edge_y>89)
        idx     = np.all(aux,axis=1)
        del aux, edge_x, edge_y 
        
        # periodic boundary edges
        pbnde   = bnde[idx==True ,:]
        
        # inner boundary edges
        bnde    = bnde[idx==False,:]
        del idx
        
        #_______________________________________________________________________
        #link correctly remaining points of periodic boundary
        lonly_p   = []
        aux_pbnde = np.sort(pbnde.reshape(pbnde.size))
        for ii in range(0,aux_pbnde.size):
            if np.sum(aux_pbnde==aux_pbnde[ii])==1:
                lonly_p.append(aux_pbnde[ii])
                
        lonly_p = np.array(lonly_p)
        sortvec = np.argsort(self.nodes_2d_yg[np.int32(lonly_p)])
        sortvec = np.array(sortvec[::-1]) #descent order
        lonly_p = lonly_p[sortvec]
        
        #_______________________________________________________________________
        addp_edge = [];
        xmin = self.nodes_2d_xg[np.int32(lonly_p)].min()
        xmax = self.nodes_2d_xg[np.int32(lonly_p)].max()
        
        #____START WHILE LOOP___________________________________________________
        while lonly_p.size>0:
            
            #___________________________________________________________________
            P1_i  = lonly_p[1];
            P1_x  = self.nodes_2d_xg[P1_i];
            P1_y  = self.nodes_2d_yg[P1_i];
            aux_i = np.squeeze(np.array(np.where(np.logical_and(\
                                self.nodes_2d_xg[lonly_p]==P1_x,lonly_p!=P1_i))))
            
            if np.size(aux_i)!=0:
                # find land mnask points on either the left or right
                # side of pbnd
                aux2_i = np.argmin(P1_y-self.nodes_2d_yg[lonly_p[aux_i]])
                
            else:
                # antarctica land mask points 
                if P1_x == xmin:
                    aux_i = np.array(np.where(self.nodes_2d_xg[lonly_p]==xmax)[0])
                    
                else:
                    aux_i = np.array(np.where(self.nodes_2d_xg[lonly_p]==xmin)[0])
                aux2_i = np.array(np.argmin(np.abs(P1_y-self.nodes_2d_yg[lonly_p[aux_i]])))
                
            P2_i      = lonly_p[aux_i[aux2_i]]
            addp_edge.append([P1_i,P2_i])
            #print(lonly_p) 
            lonly_p   = np.array(lonly_p[np.where(lonly_p!=P1_i)])
            lonly_p   = np.array(lonly_p[np.where(lonly_p!=P2_i)])
        
        #____END WHILE LOOP_____________________________________________________
        addp_edge = np.array(addp_edge)
        addp_edge = np.sort(addp_edge,axis=1)
        bnde      = np.concatenate((bnde,addp_edge),axis=0)
        del aux_i; del aux2_i; del P1_i; del P1_x; del P1_y; del P2_i; 
        del lonly_p; del addp_edge
        
        #_______________________________________________________________________
        # elimiate edges that go from -180 -->180 or 180--> -180
        edge_x  = np.sort(self.nodes_2d_xg[bnde],axis=1) ;
        idx     = (edge_x[:,1]-edge_x[:,0])>180;
        pbnde   = bnde[idx==True ,:];
        bnde    = bnde[idx==False,:];
        del edge_x; del idx
            
        #_______________________________________________________________________
        aux_x     = self.nodes_2d_xg[pbnde]
        aux_y     = self.nodes_2d_yg[pbnde]
        
        sortvec   = np.argsort(aux_y[:,0])
        aux_y     = aux_y[sortvec,:]
        aux_x     = aux_x[sortvec,:]
        aux_pbnde = pbnde[sortvec,:]
        
        sortvec = np.argsort(aux_x,axis=1)
        for ii in range(0,aux_y.shape[0]):
            aux_y[ii,:]     = aux_y[ii,sortvec[ii,:]]
            aux_pbnde[ii,:] = aux_pbnde[ii,sortvec[ii,:]]
            
        del sortvec
        
        #_______________________________________________________________________
        count_add_p  = np.array([self.n2dna])
        add_points_x = []
        add_points_y = [] 
        add_edge = []
        for ii in range(0,np.int(np.ceil(np.float16(aux_pbnde.shape[0])/2))) :
            if ii==0: # Antarctica
                y_m = np.sum(aux_y[0,:])/2
                
                dum_x = np.linspace(xmin,xmax,361);
                dum_y = dum_x*0-90
                count_add_p  = count_add_p[-1]+np.arange(0,dum_x.size+2,1)
                
                add_points_x.append([xmin])
                add_points_y.append([y_m])
                for jj in range(0,dum_x.size):
                    add_points_x.append([dum_x[jj]])
                    add_points_y.append([dum_y[jj]])
                add_points_x.append([xmax])
                add_points_y.append([y_m])
                
                add_edge.append([aux_pbnde[0,0], count_add_p[0]])
                for jj in range(0,count_add_p.size-1):
                    add_edge.append([count_add_p[jj], count_add_p[jj+1]])
                add_edge.append([count_add_p[-1], aux_pbnde[0,1]])
                
            else:
                count = ii+(ii-2);
                y_m = sum(aux_y[count:count+1,:])/2
                # add left additional edges and points
                count_add_p  = count_add_p[-1]+[1,2]
                #add_points_x.append([xmin, xmin])
                #add_points_y.append([y_m[1] y_m[2]])
                add_points_x.append([xmin])
                add_points_x.append([xmin])
                add_points_y.append([y_m[1]])
                add_points_y.append([y_m[2]])
                add_edge.append([aux_pbnde[count,0], count_add_p[0]])
                add_edge.append([count_add_p[0], count_add_p[1]])
                add_edge.append([count_add_p[1], aux_pbnde[count+1,0]])
                
                        
                # add right additional edges and points
                count_add_p  = count_add_p[-1]+[1,2]
                #add_points_x.append([xmax, xmax])
                #add_points_y.append([y_m[1] y_m[2]])
                add_points_x.append([xmax])
                add_points_x.append([xmax])
                add_points_y.append([y_m[1]])
                add_points_y.append([y_m[2]])
                add_edge.append([aux_pbnde[count,1], count_add_p[0]])
                add_edge.append([count_add_p[0], count_add_p[1]])
                add_edge.append([count_add_p[1], aux_pbnde[count+1,1]])
        #t5 = time.time()
        #_______________________________________________________________________
        # add new edge to all edge array
        add_edge = np.array(add_edge)
        add_edge = np.sort(add_edge,axis=1)
        bnde     = np.concatenate((bnde,add_edge),axis=0)
        nbnde    = bnde.shape[0];
                
        add_points_x = np.squeeze(np.array(add_points_x))
        add_points_y = np.squeeze(np.array(add_points_y))
        add_nodes_2d_xg = np.concatenate((self.nodes_2d_xg,add_points_x),axis=0)
        add_nodes_2d_yg = np.concatenate((self.nodes_2d_yg,add_points_y),axis=0)
            
        del y_m; del add_edge; del add_points_x; del add_points_y; del count_add_p
        del aux_x; del aux_y
        
        #_______________________________________________________________________
        run_cont        = np.zeros((1,nbnde))*np.nan
        run_cont[0,:2]  = bnde[0,:] # initialise the first landmask edge
        run_bnde        = bnde[1:,:] # remaining edges that still need to be distributed
        count_init      = 1;
        init_ind        = run_cont[0,0];
        ind_lc_s        = 0;
        
        polygon_xy = []
        for ii in range(0,nbnde):
            #___________________________________________________________________
            # search for next edge that contains that contains the last node index from 
            # run_cont
            kk_rc = np.column_stack(np.where( run_bnde==np.int(run_cont[0,count_init]) ))
            kk_r  = kk_rc[:,0]
            kk_c  = kk_rc[:,1]
            count_init  = count_init+1
            
            #___________________________________________________________________
            if kk_c[0] == 0 :
                run_cont[0,count_init] = run_bnde[kk_r[0],1]
            else:
                run_cont[0,count_init] = run_bnde[kk_r[0],0]
            
            #___________________________________________________________________
            # if a land sea mask polygon is closed
            if  np.any(run_bnde[kk_r[0],:] == init_ind):
                
                count_init  = count_init+1
                
                aux_lx = add_nodes_2d_xg[np.int64(run_cont[0,0:count_init])];
                aux_ly = add_nodes_2d_yg[np.int64(run_cont[0,0:count_init])];
                aux_xy = np.zeros((count_init,2))
                aux_xy[:,0] = aux_lx
                aux_xy[:,1] = aux_ly
                polygon_xy.append(aux_xy)
                del aux_lx; del aux_ly; del aux_xy
                
                ind_lc_s = ind_lc_s+count_init+1;
                
                count_init = count_init+1
                
                aux_ind  = np.arange(0,run_bnde.shape[0],1)
                run_bnde = run_bnde[aux_ind!=kk_r[0],:]
                if np.size(run_bnde)==0:
                    break
                
                #_______________________________________________________________
                run_cont        = np.zeros((1,nbnde))*np.nan
                run_cont[0,:2]  = run_bnde[0,:]
                run_bnde        = run_bnde[1:,:]
                count_init=1;
                init_ind        = run_cont[0,0]
            else:
                aux_ind =np.arange(0,run_bnde.shape[0],1)
                run_bnde=run_bnde[aux_ind!=kk_r[0],:]
            
        #_______________________________________________________________________
        self.polygon_xy = polygon_xy
        #t6 = time.time()
        #print('t6-t5:',t6-t5)
    
#+___ROTATE VECTOR GRID ROT-->GEO______________________________________________+
#|                                                                             |
#+_____________________________________________________________________________+
def fesom_vector_rot(mesh,u,v,do_output=True):
    
    if do_output==True: print(' --> do vector rotation rot2geo')
    #___________________________________________________________________________
    alpha = mesh.alpha
    beta  = mesh.beta
    gamma = mesh.gamma
    
    #___________________________________________________________________________
    if len(mesh.nodes_2d_yr)==0:
        mesh.fesom_grid_rot_g2r()
    
    #___________________________________________________________________________
    # build euler rotation matrix
    rotate_matrix = np.zeros((3,3))
    rotate_matrix[0,0] =( math.cos(math.radians(gamma)) * 
                          math.cos(math.radians(alpha)) -
                          math.sin(math.radians(gamma)) *
                          math.cos(math.radians(beta )) *
                          math.sin(math.radians(alpha)))
    
    rotate_matrix[0,1] =( math.cos(math.radians(gamma)) *
                          math.sin(math.radians(alpha)) +
                          math.sin(math.radians(gamma)) *
                          math.cos(math.radians(beta )) *
                          math.cos(math.radians(alpha)))
    
    rotate_matrix[0,2] =( math.sin(math.radians(gamma)) *
                          math.sin(math.radians(beta )))
    
    rotate_matrix[1,0] =(-math.sin(math.radians(gamma)) *
                          math.cos(math.radians(alpha)) -
                          math.cos(math.radians(gamma)) *
                          math.cos(math.radians(beta )) *
                          math.sin(math.radians(alpha)))
    
    rotate_matrix[1,1] =(-math.sin(math.radians(gamma)) *
                          math.sin(math.radians(alpha)) +
                          math.cos(math.radians(gamma)) *
                          math.cos(math.radians(beta )) *
                          math.cos(math.radians(alpha)))
    
    rotate_matrix[1,2] =( math.cos(math.radians(gamma)) *
                          math.sin(math.radians(beta )))
    
    rotate_matrix[2,0] =( math.sin(math.radians(beta )) *
                          math.sin(math.radians(alpha)))
    
    rotate_matrix[2,1] =(-math.sin(math.radians(beta )) *
                          math.cos(math.radians(alpha)))
    
    rotate_matrix[2,2] =( math.cos(math.radians(beta )))
    
    #___________________________________________________________________________
    # make inverse of rotation matrix
    #from numpy.linalg import inv
    rotate_matrix=inv(rotate_matrix);
    
    #___________________________________________________________________________
    u = np.squeeze(np.array(u))
    v = np.squeeze(np.array(v))
    #___________________________________________________________________________
    #if u.shape[0]==mesh.n2dna:
    if any(x in u.shape for x in [mesh.n2dna]):
        which_size = mesh.n2dna
        rlat = np.radians(mesh.nodes_2d_yr)
        rlon = np.radians(mesh.nodes_2d_xr)
        lat  = np.radians(mesh.nodes_2d_yg)
        lon  = np.radians(mesh.nodes_2d_xg)
        
    #elif u.shape[0]==mesh.n2dn:
    elif any(x in u.shape for x in [mesh.n2dn]):
        which_size = mesh.n2dn
        rlat = np.radians(mesh.nodes_2d_yr[0:mesh.n2dn])
        rlon = np.radians(mesh.nodes_2d_xr[0:mesh.n2dn])
        lat  = np.radians(mesh.nodes_2d_yg[0:mesh.n2dn])
        lon  = np.radians(mesh.nodes_2d_xg[0:mesh.n2dn])
        
    #elif u.shape[0]==mesh.n2dea:
    elif any(x in u.shape for x in [mesh.n2dea]):
        which_size = mesh.n2dea
        rlat = mesh.nodes_2d_yr[mesh.elem_2d_i]
        rlon = mesh.nodes_2d_xr[mesh.elem_2d_i]
        aux  =np.where( mesh.nodes_2d_x[mesh.elem0_2d_i].max(axis=1)-mesh.nodes_2d_x[mesh.elem0_2d_i].min(axis=1)>=180)[0]
        for ii in aux:
            for jj in range(0,3):
                if rlon[ii,jj]<rlon[ii,:].max()-180: rlon[ii,jj] = rlon[ii,jj]+360
        rlat = np.radians(rlat.sum(axis=1)/3.0)
        rlon = np.radians(rlon.sum(axis=1)/3.0)
        lat  = np.radians(mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0)
        lon  = np.radians(mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3.0)
        
    #elif u.shape[0]==mesh.n2de:
    elif any(x in u.shape for x in [mesh.n2de]):
        which_size = mesh.n2de
        rlat = mesh.nodes_2d_yr[mesh.elem0_2d_i]
        rlon = mesh.nodes_2d_xr[mesh.elem0_2d_i]
        aux  =np.where( mesh.nodes_2d_x[mesh.elem0_2d_i].max(axis=1)-mesh.nodes_2d_x[mesh.elem0_2d_i].min(axis=1)>=180)[0]
        for ii in aux:
            for jj in range(0,3):
                if rlon[ii,jj]<rlon[ii,:].max()-180: rlon[ii,jj] = rlon[ii,jj]+360
        rlat = np.radians(rlat.sum(axis=1)/3.0)
        rlon = np.radians(rlon.sum(axis=1)/3.0)
        lat  = np.radians(mesh.nodes_2d_yg[mesh.elem0_2d_i].sum(axis=1)/3.0)
        lon  = np.radians(mesh.nodes_2d_xg[mesh.elem0_2d_i].sum(axis=1)/3.0)
        
    ndims = u.shape
    if len(ndims)==1:
        #___________________________________________________________________________
        txg = -v*np.sin(rlat)*np.cos(rlon) - u*np.sin(rlon)
        tyg = -v*np.sin(rlat)*np.sin(rlon) + u*np.cos(rlon)
        tzg =  v*np.cos(rlat)
        
        txr = rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg;
        tyr = rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg;
        tzr = rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg;
        
        #___________________________________________________________________________
        v = txr*(-np.sin(lat))*np.cos(lon) - tyr*np.sin(lat)*np.sin(lon) + tzr*np.cos(lat)
        u = txr*(-np.sin(lon)) + tyr*np.cos(lon)
        
        #___________________________________________________________________________
        u = np.transpose(np.array(u,ndmin=2)).squeeze();
        v = np.transpose(np.array(v,ndmin=2)).squeeze();
    elif len(ndims)    == 2:
        
        which_index = list(u.shape).index(which_size)
        txg,tyg,tzg = np.zeros(u.shape), np.zeros(u.shape), np.zeros(u.shape)
            
        if which_index==0:
            #___________________________________________________________________________
            for ii in range(0,u.shape[1]):
                txg[:,ii] = -v[:,ii]*np.sin(rlat)*np.cos(rlon) - u[:,ii]*np.sin(rlon)
                tyg[:,ii] = -v[:,ii]*np.sin(rlat)*np.sin(rlon) + u[:,ii]*np.cos(rlon)
                tzg[:,ii] =  v[:,ii]*np.cos(rlat)
            
            txr = rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg;
            tyr = rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg;
            tzr = rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg;
            
            #___________________________________________________________________________
            for ii in range(0,u.shape[1]):
                v[:,ii] = txr[:,ii]*(-np.sin(lat))*np.cos(lon) - \
                          tyr[:,ii]*np.sin(lat)*np.sin(lon) + \
                          tzr[:,ii]*np.cos(lat)
                u[:,ii] = txr[:,ii]*(-np.sin(lon)) + tyr[:,ii]*np.cos(lon)
                
        elif which_index==1:
            #___________________________________________________________________________
            for ii in range(0,u.shape[0]):
                txg[ii,:] = -v[ii,:]*np.sin(rlat)*np.cos(rlon) - u[ii,:]*np.sin(rlon)
                tyg[ii,:] = -v[ii,:]*np.sin(rlat)*np.sin(rlon) + u[ii,:]*np.cos(rlon)
                tzg[ii,:] =  v[ii,:]*np.cos(rlat)
            
            txr = rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg;
            tyr = rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg;
            tzr = rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg;
            
            #___________________________________________________________________________
            for ii in range(0,u.shape[0]):
                v[ii,:] = txr[ii,:]*(-np.sin(lat))*np.cos(lon) - \
                         tyr[ii,:]*np.sin(lat)*np.sin(lon) + \
                         tzr[ii,:]*np.cos(lat)
                u[ii,:] = txr[ii,:]*(-np.sin(lon)) + tyr[ii,:]*np.cos(lon)
            
    return(u,v)
    
    
#+___EQUIVALENT OF MATLAB ISMEMBER FUNCTION____________________________________+
#|                                                                             |
#+_____________________________________________________________________________+
def ismember_rows(a, b):
    return np.flatnonzero(np.in1d(b[:,0], a[:,0]) & np.in1d(b[:,1], a[:,1]))
    
    
#+___TRANSFORM GEOCOORD to CARTESIAN COORD_____________________________________+
#|                                                                             |
#+_____________________________________________________________________________+
def geo2cart(glon,glat):
    Rearth=6371.0;
    x_cart=Rearth*np.cos(np.radians(glat))*np.cos(np.radians(glon))
    y_cart=Rearth*np.cos(np.radians(glat))*np.sin(np.radians(glon))
    z_cart=Rearth*np.sin(np.radians(glat))
    return np.array([x_cart,y_cart,z_cart])
