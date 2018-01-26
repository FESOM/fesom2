# Patrick Scholz, 14.12.2017
import time
import numpy as np
import pandas as pa
import math
from numpy.linalg import inv
from set_inputarray import *


#___LOAD FESOM MESH FROM FILE ALL IN ONE________________________________________
# input : inputarray
#_______________________________________________________________________________
def fesom_init_mesh(inputarray):
	global mesh
	#___________________________________________________________________________
	# initialse mesh object
	mesh = fesom_mesh
	
	#___________________________________________________________________________
	# define euler angles
	mesh.id    = inputarray['mesh_id'   ]
	mesh.path  = inputarray['mesh_dir'  ]
	mesh.alpha = inputarray['mesh_alpha']
	mesh.beta  = inputarray['mesh_beta' ]
	mesh.gamma = inputarray['mesh_gamma']
	mesh.focus = inputarray['mesh_focus']
	
	#___________________________________________________________________________
	# load mesh from file
	mesh = fesom_load_mesh(mesh)
	
	#___________________________________________________________________________
	# rotate mesh from rot to geo coordinates
	mesh = fesom_grid_rot_r2g(mesh,'r2g')
	
	#___________________________________________________________________________
	# change grid focus from -180...180 to e.g. 0...360
	if (inputarray['mesh_focus']!=0):
		mesh = fesom_grid_rot_r2g('focus')
	
	#___________________________________________________________________________
	# remove+augment periodic boundary
	mesh = fesom_remove_pbnd(mesh)
	
	# calculate fesom land mask interactivly
	mesh = fesom_calc_landmask(mesh)
	
	mesh = fesom_grid_rot_g2r(mesh)
	
	# calc triangle area in km^2 
	mesh = fesom_calc_triarea(mesh)
	
	# calc triangle resol in km interpolated to nodes
	#mesh = fesom_calc_triresol(mesh)
	
	mesh.metadata = inputarray
	
	return mesh
	
	
class fesom_mesh:
	"""existing instances are: path, n2d, e2d, nlev, zlevs, x2, y2, elem, n32, no_cyclic_elem, alpha, beta, gamma"""
	def __init__(self):
		self.metadata	    = ''
		self.path 	 	    = ''
		self.id 		    = ''
		self.path 	 	    = ''
		#____mesh euler angles_______________________
		self.alpha		    = 0
		self.beta		    = 0
		self.gamma 		    = 0
		self.focus 		    = 0
		#____mesh nodes & elem size__________________
		self.n2dn		    = 0
		self.n2dna		    = 0
		self.n2de 		    = 0
		self.n2dea		    = 0
		#____mesh vertical info______________________
		self.nlev		    = 0
		self.zlev		    = []
		self.zmid		    = []
		#____mesh nodes info_________________________
		self.nodes_2d_x     = []
		self.nodes_2d_y     = []
		self.nodes_2d_z     = []
		self.nodes_2d_i     = []
		self.nodes_2d_xg    = []
		self.nodes_2d_yg    = []
		self.nodes_2d_xr    = []
		self.nodes_2d_yr    = []
		self.pbndn_2d_i     = []
		self.pbndn_2d_iz    = []
		#____mesh elem info__________________________
		self.elem0_2d_i     = []
		self.elem_2d_i      = []
		self.elem0_2d_iz    = []
		self.pbndtri_2d_i   = []
		self.abnormtri_2d_i = []
		#____fesom lsmask polygon____________________
		self.polygon_xy     = []
		#____triangle area in km^2___________________
		self.elem_2d_area   = []
		self.elem_2d_resol  = []
		#____triangle area in km^2___________________
		self.nodes_2d_resol = []
		self.nodes_2d_area  = []
		#____mesh parameter for vertical interpolation
		
	def __str__(self):
		return "mesh path=%s" % self.path
	
	
#___LOAD FESOM MESH FROM FILE___________________________________________________
# input : inputarray
#_______________________________________________________________________________
def fesom_load_mesh(mesh):
	
	print('')
	print('___LOAD FESOM MESH_________________________________________')
	print(' --> read grid files')
	
	#____load 2d node matrix____________________________________________________
	print('     > nod2d.out')
	fid_n2d          = open(mesh.path+'nod2d.out', 'r')
	mesh.n2dn        = np.uint32(fid_n2d.readline().strip())
	dum              = np.matrix(pa.read_table(fid_n2d, header=-1,delim_whitespace=True))
	fid_n2d.close()
	del fid_n2d
	mesh.nodes_2d_x  = np.float32(np.array(dum[:,1]).reshape(mesh.n2dn))
	mesh.nodes_2d_y  = np.float32(np.array(dum[:,2]).reshape(mesh.n2dn))
	mesh.nodes_2d_i  = np.uint16( np.array(dum[:,3]).reshape(mesh.n2dn))
	del dum
	
	
	#____load 2d element matrix_________________________________________________
	print('     > elem2d.out')
	fid_e2d          = open(mesh.path+'elem2d.out', 'r')
	mesh.n2de      	 = np.int64(fid_e2d.readline().strip())
	## pandas module fastest option but not everywherre available
	mesh.elem0_2d_i  = np.matrix(pa.read_table(fid_e2d, header=-1,delim_whitespace=True),dtype='uint32')-1
	fid_e2d.close()
	
	
	#____load 3d nodes alligned under 2d nodes__________________________________
	print('     > aux3d.out') 
	fid_aux3d        = open(mesh.path+'aux3d.out', 'r')
	mesh.nlev        = np.uint16(fid_aux3d.readline().strip())
	dum			     = np.array(pa.read_table(fid_aux3d, 
												  header=-1,
												  delim_whitespace=True),
									dtype='int16')
	mesh.zlev        = dum[0:mesh.nlev].squeeze()
	
	mesh.zmid		 = (mesh.zlev[:-1]+mesh.zlev[1:])/2.
	
	#mesh.nodes_2d_z  = np.squeeze(dum[mesh.nlev:mesh.nlev+mesh.n2dn])
	fid_aux3d.close()
	del fid_aux3d
	del dum
	
	#____load number of levels at each node_____________________________________
	print('     > nlvls.out') 
	fid				 = open(mesh.path+'nlvls.out', 'r')
	mesh.nodes_2d_iz = np.array(pa.read_table(fid, header=-1,delim_whitespace=True))
	mesh.nodes_2d_iz = mesh.nodes_2d_iz.squeeze()
	mesh.nodes_2d_z  = np.float32(mesh.zlev[mesh.nodes_2d_iz])
	fid.close()
	
	#____load number of levels at each elem_____________________________________
	print('     > elvls.out') 
	fid				 = open(mesh.path+'elvls.out', 'r')
	mesh.elem0_2d_iz  = np.array(pa.read_table(fid, header=-1,delim_whitespace=True))
	mesh.elem0_2d_iz  = mesh.elem0_2d_iz.squeeze()
	fid.close()
	
	#___________________________________________________________________________
	return mesh
	
	
#___ROTATE GRID ROT-->GEO_______________________________________________________
# input : r2g (default) change coordinate from rotated-->geo
#         g2r           change coordinate from geo-->rotated
#         focus         change longitudinal focus 
#_______________________________________________________________________________
def fesom_grid_rot_r2g(mesh,str_mode):
	
	if (str_mode == 'focus'):
		print(' --> change mesh focus');
		alpha = -mesh.focus
		beta  = 0.
		gamma = 0.
	else:
		if (str_mode == 'r2g'):
			print(' --> rotate mesh rot2geo')
		else:
			print(' --> rotate mesh geo2rot')
			
		alpha = mesh.alpha
		beta  = mesh.beta
		gamma = mesh.gamma
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
	if (str_mode == 'r2g') or (str_mode == 'focus'):
		from numpy.linalg import inv
		rotate_matrix=inv(rotate_matrix);
		
	#____3D_____________________________________________________________________
	# calculate Cartesian coordinates
	if (str_mode == 'focus') or (str_mode == 'g2r'):
		xr=np.cos(np.radians(mesh.nodes_2d_yg)) * np.cos(np.radians(mesh.nodes_2d_xg));
		yr=np.cos(np.radians(mesh.nodes_2d_yg)) * np.sin(np.radians(mesh.nodes_2d_xg));
		zr=np.sin(np.radians(mesh.nodes_2d_yg));
	else:
		xr=np.cos(np.radians(mesh.nodes_2d_y)) * np.cos(np.radians(mesh.nodes_2d_x));
		yr=np.cos(np.radians(mesh.nodes_2d_y)) * np.sin(np.radians(mesh.nodes_2d_x));
		zr=np.sin(np.radians(mesh.nodes_2d_y));
	
	#___________________________________________________________________________
	# rotate to geographical cartesian coordinates:
	xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr;
	yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr;
	zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr;
	
	##___________________________________________________________________________
	#mesh.nodes_2d_yg = np.degrees(np.arcsin(zg));
	#mesh.nodes_2d_xg = np.degrees(np.arctan2(yg,xg));
	
	#___________________________________________________________________________
	if (str_mode == 'focus') or (str_mode == 'r2g'):
		mesh.nodes_2d_yg=np.degrees(np.arcsin(zg));
		mesh.nodes_2d_xg=np.degrees(np.arctan2(yg,xg));
		if (str_mode == 'focus'):
			mesh.nodes_2d_xg = mesh.nodes_2d_xg+mesh.focus
	elif (str_mode == 'g2r'):
		mesh.nodes_2d_yr=np.degrees(np.arcsin(zg));
		mesh.nodes_2d_xr=np.degrees(np.arctan2(yg,xg));
		if (str_mode == 'focus'):
			mesh.nodes_2d_xr = mesh.nodes_2d_xr+mesh.focus
		
	#___________________________________________________________________________
	return mesh
	
	
#___ROTATE GRID GEO-->ROT_______________________________________________________
# input : g2r           change coordinate from geo-->rotated
#_______________________________________________________________________________
def fesom_grid_rot_g2r(mesh):
	mesh = fesom_grid_rot_r2g(mesh,'g2r')
	return mesh
	
	
#___SEARCH AND REMOVE PERIODIC BOUNDARIES_______________________________________
# identify the periodic boundary (nodind2d ==5000) and deletes it. Augment 
# the left and right part
#_______________________________________________________________________________
def fesom_remove_pbnd(mesh):
	
	print(' --> remove cyclic boundary')
	#___________________________________________________________________________
	# edge list
	edge_1_i = mesh.elem0_2d_i[:,[0,1]]
	edge_2_i = mesh.elem0_2d_i[:,[1,2]]
	edge_3_i = mesh.elem0_2d_i[:,[2,0]]
	
	#___________________________________________________________________________
	# length of edges
	ledge_1_i= np.abs(mesh.nodes_2d_xg[edge_1_i][:,0]-mesh.nodes_2d_xg[edge_1_i][:,1])
	ledge_2_i= np.abs(mesh.nodes_2d_xg[edge_2_i][:,0]-mesh.nodes_2d_xg[edge_2_i][:,1])
	ledge_3_i= np.abs(mesh.nodes_2d_xg[edge_3_i][:,0]-mesh.nodes_2d_xg[edge_3_i][:,1])
	
	#___________________________________________________________________________
	# indice of period boundary edges [n x 2]
	ind_edge_1 = np.where( 
					np.squeeze( 
						np.array( [ledge_1_i >= 180] ) ) )[0]
	ind_edge_2 = np.where( 
					np.squeeze( 
						np.array( [ledge_2_i >= 180] ) ) )[0]
	ind_edge_3 = np.where( 
					np.squeeze( 
						np.array( [ledge_3_i >= 180] ) ) )[0]
	del ledge_1_i
	del ledge_2_i
	del ledge_3_i
	
	#___________________________________________________________________________
	# all unique indice of period boundary edges [n*2 x 0]
	ind_node_e1  = np.array(edge_1_i[ind_edge_1].reshape([np.size(ind_edge_1)*2]),dtype='uint32')
	ind_node_e2  = np.array(edge_2_i[ind_edge_2].reshape([np.size(ind_edge_2)*2]),dtype='uint32')
	ind_node_e3  = np.array(edge_3_i[ind_edge_3].reshape([np.size(ind_edge_3)*2]),dtype='uint32')
	ind_node_all = np.unique(np.concatenate((ind_node_e1,ind_node_e2,ind_node_e3),axis=1))
	del ind_node_e1
	del ind_node_e2
	del ind_node_e3
	del ind_edge_1
	del ind_edge_2
	del ind_edge_3
	del edge_1_i
	del edge_2_i
	del edge_3_i
	
	#___________________________________________________________________________
	# node indices of the 2d left and right periodic boundary
	xmin = np.min(mesh.nodes_2d_xg)
	xmax = np.max(mesh.nodes_2d_xg)
	
	aux_l_i = np.where(
				np.squeeze(
					np.array(
						[ (mesh.nodes_2d_xg[ind_node_all] < (xmin+(xmax-xmin)/2) ) #& 
						])))[0]
	aux_r_i = np.where(
				np.squeeze(
					np.array(
						[ (mesh.nodes_2d_xg[ind_node_all] > (xmin+(xmax-xmin)/2) ) #& 
						  #(mesh['nodes_2d_yg'][ind_node_all] < 89.75) &
						  #((180-np.abs(mesh['nodes_2d_xg'][ind_node_all]))<55 )
						])))[0]
	pbndn_l_2d_i 	= ind_node_all[aux_l_i]
	npbnd_l_2d   	= pbndn_l_2d_i.size
	pbndn_r_2d_i 	= ind_node_all[aux_r_i]
	mesh.pbndn_2d_i = np.concatenate((pbndn_r_2d_i,pbndn_l_2d_i))
	npbnd_r_2d   	= pbndn_r_2d_i.size
	
	del aux_l_i
	del aux_r_i
	del ind_node_all
	del xmin
	del xmax
	
	#___________________________________________________________________________
	#remove right & left periodic boundary
	aux_pos    = np.zeros(mesh.n2dn,dtype='uint32')
	aux_nnodesr= np.linspace(mesh.n2dn,
						 	 mesh.n2dn+npbnd_r_2d-1,
						 	 npbnd_r_2d,dtype='uint32')
	
	aux_pos[pbndn_r_2d_i]=aux_nnodesr
	aux_nnodesl= np.linspace(mesh.n2dn+npbnd_r_2d,
							 mesh.n2dn+npbnd_r_2d+npbnd_l_2d-1,
							 npbnd_l_2d,dtype='uint32')
	aux_pos[pbndn_l_2d_i]=aux_nnodesl
	del aux_nnodesl
	del aux_nnodesr
	
	#___________________________________________________________________________
	# Augment the mesh on the right periodic boundary
	# Straightforward part
	aux_xcyc  = np.transpose(np.matrix(mesh.nodes_2d_xg[mesh.pbndn_2d_i]))
	aux_ycyc  = np.transpose(np.matrix(mesh.nodes_2d_yg[mesh.pbndn_2d_i]))
	
	aux_xycyc = np.concatenate((aux_xcyc,
							 	aux_ycyc),axis=1)
	
	wall_pbnd_l = np.floor(np.min(aux_xcyc))
	wall_pbnd_r = np.ceil(np.max(aux_xcyc))
	
	mesh.nodes_2d_xg = np.concatenate((mesh.nodes_2d_xg,
										  np.zeros(npbnd_r_2d)+wall_pbnd_l,
										  np.zeros(npbnd_l_2d)+wall_pbnd_r    ))
	mesh.nodes_2d_yg = np.concatenate((mesh.nodes_2d_yg,
										  mesh.nodes_2d_yg[mesh.pbndn_2d_i]   ))
	mesh.nodes_2d_z = np.concatenate((mesh.nodes_2d_z,
										  mesh.nodes_2d_z[mesh.pbndn_2d_i]    ))
	mesh.nodes_2d_i = np.concatenate((mesh.nodes_2d_i,
										  mesh.nodes_2d_i[mesh.pbndn_2d_i]    ))
	mesh.nodes_2d_iz = np.concatenate((mesh.nodes_2d_iz,
										  mesh.nodes_2d_iz[mesh.pbndn_2d_i]   ))
	
	#____loop over all periodic bnd. segments___________________________________
	copy_elem = np.zeros(mesh.n2de,dtype='bool')
	mesh.elem_2d_i = mesh.elem0_2d_i;
	
	#___________________________________________________________________________
	# (ii.a) 2d elements:
	# List all triangles that touch the cyclic boundary segments
	#___________________________________________________________________________
	nn = np.array([1,2,3,1,2],dtype='uint8')-1
	for nn_i in range(0, 5):
		#_______________________________________________________________________
		aux_xc  = mesh.nodes_2d_xg[mesh.elem_2d_i[:,nn[nn_i]]]
		aux_yc  = mesh.nodes_2d_yg[mesh.elem_2d_i[:,nn[nn_i]]]
		aux_xyc = np.concatenate((aux_xc,aux_yc),axis=1)
		
		#_______________________________________________________________________
		# list all triangles that touch the periodic boundary segments
		n     = ismember_rows(aux_xycyc,aux_xyc)
		tlist = mesh.elem_2d_i[n,:]
		maxt  = np.max(mesh.nodes_2d_xg[tlist],axis=1)
		mint  = np.min(mesh.nodes_2d_xg[tlist],axis=1)
		mind  = np.squeeze(np.array([maxt-mint>=180]))
		n     = n[np.where(mind)[0]]
		del tlist
		del maxt
		del mint 
		del mind
		
		#_______________________________________________________________________
		# dublicate triangles that span over the periodic boundary
		n_copy = n[np.where( (n<=mesh.n2de-1) & 
					  		 (np.squeeze((aux_xc[n]>wall_pbnd_l+(wall_pbnd_r-wall_pbnd_l)/2))) )]
		n_copy = n_copy[np.where(copy_elem[n_copy]==False)]
		
		# duplicate triangles that have connection to the right pbndn
		if nn_i==0:
			mesh.pbndtri_2d_i = n_copy
		else:
			mesh.pbndtri_2d_i = np.concatenate((mesh.pbndtri_2d_i,n_copy))
		mesh.elem_2d_i=np.concatenate((mesh.elem_2d_i,mesh.elem_2d_i[n_copy,:]))
		
		# new aux_xc and aux_yc because of duplicated triangles
		aux_xc  = mesh.nodes_2d_xg[mesh.elem_2d_i[:,nn[nn_i]]]
		aux_yc  = mesh.nodes_2d_yg[mesh.elem_2d_i[:,nn[nn_i]]]
		aux_xyc = np.concatenate((aux_xc,aux_yc),axis=1)
				
		# mark already copied triangles
		copy_elem[n_copy]=True
		del n_copy
		
		#_______________________________________________________________________
		# rearange element list for the right side periodic triangles
		n_nodechange = n[np.where( (n<=mesh.n2de-1) & 
								   (np.squeeze((aux_xc[n]>wall_pbnd_l+(wall_pbnd_r-wall_pbnd_l)/2) )) & 
								   (np.squeeze(np.array(mesh.elem_2d_i[n,nn[nn_i]]<=mesh.n2dn-1)))
								)]
		mesh.elem_2d_i[n_nodechange,nn[nn_i]]=np.squeeze(aux_pos[mesh.elem_2d_i[n_nodechange,nn[nn_i]]])
		del n 
		del n_nodechange
		
		#_______________________________________________________________________
		# rearange element list for the left side periodic triangles
		n     = ismember_rows(aux_xycyc,aux_xyc)
		tlist = mesh.elem_2d_i[n,:]
		maxt  = np.max(mesh.nodes_2d_xg[tlist],axis=1)
		mint  = np.min(mesh.nodes_2d_xg[tlist],axis=1)
		mind  = np.squeeze(np.array([maxt-mint>=180]))
		n     = n[np.where(mind)[0]]
		
		n_nodechange = n[np.where( (n>mesh.n2de-1) & 
									np.squeeze((aux_xc[n]<wall_pbnd_l+(wall_pbnd_r-wall_pbnd_l)/2) ) & 
									np.squeeze(np.array(mesh.elem_2d_i[n,nn[nn_i]]<=mesh.n2dn-1))
								)]
		mesh.elem_2d_i[n_nodechange,nn[nn_i]]=np.squeeze(aux_pos[mesh.elem_2d_i[n_nodechange,nn[nn_i]]])
		
		del tlist
		del maxt
		del mint 
		del mind
		del n 
		del n_nodechange
	
	#___________________________________________________________________________
	# check if there are still any abnormal triangles an kick them out
	mesh.abnormtri_2d_i=[]
	#maxt = np.max(mesh.nodes_2d_xg[mesh.elem_2d_i],axis=1)
	#mint = np.min(mesh.nodes_2d_xg[mesh.elem_2d_i],axis=1)
	#mind = np.squeeze(np.array([maxt-mint<=180]))
	#if any(mind==False):
		#mesh.abnormtri_2d_i = np.where(mind==False)
	#mesh.elem_2d_i= mesh.elem_2d_i[np.where(mind==True)[0],:]
	#del maxt
	#del mint
	#del mind	
	mesh.n2dea    = mesh.elem_2d_i.shape[0]
	
	#___________________________________________________________________________
	mesh.n2dna=mesh.n2dn+mesh.pbndn_2d_i.size # new (augmented) n2d
	del aux_pos
	
	#___________________________________________________________________________
	# build land boundary edge matrix
	edge    = np.concatenate((mesh.elem_2d_i[:,[0,1]], \
							  mesh.elem_2d_i[:,[0,2]], \
							  mesh.elem_2d_i[:,[1,2]]))
	edge    = np.sort(edge,axis=1) 
	
	# python  sortrows algorythm --> matlab equivalent
	edge    = edge.tolist()
	edge.sort()
	edge    = np.array(edge)
	idx     = np.diff(edge,axis=0)==0
	idx     = np.all(idx,axis=1)
	idx     = np.logical_or(np.concatenate((idx,np.array([False]))),np.concatenate((np.array([False]),idx)))
	bnde    = edge[idx==False,:]
	
	#___________________________________________________________________________
	# kill periodic boundary edges--> periodic boundary edges have identical lon/lat points
	edge_x  = mesh.nodes_2d_xg[bnde]
	edge_y  = mesh.nodes_2d_yg[bnde]
	aux     = np.logical_or(edge_x==np.min(mesh.nodes_2d_xg),edge_x==np.max(mesh.nodes_2d_xg))
	aux     = np.logical_or(aux,edge_y>89)
	idx     = np.all(aux,axis=1)
	del aux
	del edge_x 
	del edge_y 
		
	bnde    = bnde[idx==False,:]
	del idx
	mesh.nodes_2d_i = np.zeros((mesh.n2dna),dtype='uint8')
	mesh.nodes_2d_i[bnde.reshape(bnde.size)]=1
	del bnde; del edge 
	
	#___________________________________________________________________________
	return mesh
	
	
#___CALCULATE MESH RESOLUTION___________________________________________________
# calculate mean length of the three triangle sides 
#
#_______________________________________________________________________________
def fesom_calc_triresol(mesh):
	
	#___________________________________________________________________________
	print(' --> calc. triangle resolution in km')
	
	#___________________________________________________________________________
	nodes_2d_x 	= np.concatenate((mesh.nodes_2d_xg[:mesh.n2dn] , mesh.nodes_2d_xg[mesh.pbndn_2d_i]))
	nodes_2d_y 	= np.concatenate((mesh.nodes_2d_yg[:mesh.n2dn] , mesh.nodes_2d_yg[mesh.pbndn_2d_i]))
	elem_2d_x  	= nodes_2d_x[mesh.elem_2d_i]
	elem_2d_y  	=  nodes_2d_y[mesh.elem_2d_i]
	elem_2d_xy 	= np.array([elem_2d_x,elem_2d_y])
	del nodes_2d_x
	del nodes_2d_y
	del elem_2d_x
	
	#___________________________________________________________________________
	# calc jacobi matrix for all triangles 
	# | dx_12 dy_12 |
	# | dx_23 dy_23 |
	# | dx_31 dy_31 |_i , i=1....n2dea
	jacobian 	= elem_2d_xy[:,:,1]-elem_2d_xy[:,:,0]
	jacobian 	= np.array([jacobian,elem_2d_xy[:,:,2]-elem_2d_xy[:,:,1],elem_2d_xy[:,:,0]-elem_2d_xy[:,:,2] ])
	
	# account for triangles with periodic bounaries
	for ii in range(3):
		idx = np.where(jacobian[ii,0,:]>180); 
		jacobian[ii,0,idx] = jacobian[ii,0,idx]-360;
		idx = np.where(jacobian[ii,0,:]<-180); 
		jacobian[ii,0,idx] = jacobian[ii,0,idx]+360;
		del idx
	
	# calc from geocoord to cartesian coord
	R_earth		= 12735/2;
	jacobian 	= jacobian*R_earth*2*np.pi/360
	cos_theta 	= np.cos(np.radians(elem_2d_y)).mean(axis=1)
	del elem_2d_y
	for ii in range(3):	
		jacobian[ii,0,:] = jacobian[ii,0,:]*cos_theta;
	del cos_theta
	
	#___________________________________________________________________________
	# calc vector length dr = sqrt(dx^2+dy^2)
	jacobian 	= np.power(jacobian,2);
	jacobian 	= np.sqrt(jacobian.sum(axis=1));
	jacobian 	= jacobian.transpose()
	
	#___________________________________________________________________________
	# mean resolutiuon per element
	mesh.elem_2d_resol = jacobian.mean(axis=1)
	mesh.nodes_2d_resol= fesom_interp_e2n(mesh,mesh.elem_2d_resol)
	
	##___________________________________________________________________________
	#jacobian 	= jacobian.reshape((mesh.n2dea*3))
	#interp_e2n_mat,interp_e2n_corector = fesom_e2n_interp(mesh)
	##___________________________________________________________________________
	## mesh resolution interpolated to node points
	#mesh.nodes_2d_resol = jacobian[interp_e2n_mat] * interp_e2n_corector
	#mesh.nodes_2d_resol = np.sum(mesh.nodes_2d_resol,axis=1) / np.sum(interp_e2n_corector,axis=1)
	
	return mesh
	
	
#___CALCULATE TRIANGLE AREA("volume")___________________________________________
# calculate TRIANGLE AREA IN M^2
#
#_______________________________________________________________________________
def fesom_calc_triarea(mesh):
	
	#___________________________________________________________________________
	print(' --> calc. triangle area m^2')
	
	#___________________________________________________________________________
	nodes_2d_x 	= np.concatenate((mesh.nodes_2d_xg[:mesh.n2dn] , mesh.nodes_2d_xg[mesh.pbndn_2d_i]))
	nodes_2d_y 	= np.concatenate((mesh.nodes_2d_yg[:mesh.n2dn] , mesh.nodes_2d_yg[mesh.pbndn_2d_i]))
	elem_2d_x 	= nodes_2d_x[mesh.elem0_2d_i]
	elem_2d_y 	= nodes_2d_y[mesh.elem0_2d_i]
	elem_2d_xy 	= np.array([elem_2d_x,elem_2d_y])
	del nodes_2d_x
	del nodes_2d_y
	del elem_2d_x
	
	#___________________________________________________________________________
	# calc jacobi matrix for all triangles 
	# | dx_12 dy_12 |
	# | dx_13 dy_13 |_i , i=1....n2dea
	jacobian 	= elem_2d_xy[:,:,1]-elem_2d_xy[:,:,0]
	jacobian 	= np.array([jacobian,elem_2d_xy[:,:,2]-elem_2d_xy[:,:,0] ])
	
	# account for triangles with periodic bounaries
	for ii in range(2):
		idx = np.where(jacobian[ii,0,:]>180); 
		jacobian[ii,0,idx] = jacobian[ii,0,idx]-360;
		idx = np.where(jacobian[ii,0,:]<-180); 
		jacobian[ii,0,idx] = jacobian[ii,0,idx]+360;
		del idx
	
	# calc from geocoord to cartesian coord
	R_earth		=12735/2;
	jacobian 	= jacobian*R_earth*2*np.pi/360
	cos_theta 	= np.cos(np.radians(elem_2d_y)).mean(axis=1)
	del elem_2d_y
	for ii in range(2):	
		jacobian[ii,0,:] = jacobian[ii,0,:]*cos_theta;
	del cos_theta
	
	#___________________________________________________________________________
	# calc triangle area through vector product
	mesh.elem_2d_area = 0.5*np.abs(jacobian[0,0,:]*jacobian[1,1,:]-jacobian[0,1,:]*jacobian[1,0,:])
	mesh.elem_2d_area = np.concatenate((mesh.elem_2d_area,mesh.elem_2d_area[mesh.pbndtri_2d_i]))
	
	#___________________________________________________________________________
	# calc node cluster area 
	mesh.nodes_2d_area=np.zeros(mesh.n2dna)
	for i in range(3):
		for j in range(mesh.n2de):
			n=mesh.elem0_2d_i[j,i]
			mesh.nodes_2d_area[n]=mesh.nodes_2d_area[n]+mesh.elem_2d_area[j]
	mesh.nodes_2d_area=mesh.nodes_2d_area/3.
	mesh.nodes_2d_area[mesh.n2dn:mesh.n2dna] = mesh.nodes_2d_area[mesh.pbndn_2d_i]
	
	#___________________________________________________________________________
	#if len(mesh.abnormtri_2d_i)!=0:
			#mesh.elem_2d_area = np.delete(mesh.elem_2d_area,mesh.abnormtri_2d_i)
	return mesh
	
	
#___CALCULATE INTERPOLATION MATRIX FROM ELEMENTS TO NODE POINTS RESOLUTION______
# interpolate from triangle values to node values
#	
#_______________________________________________________________________________
def fesom_interp_e2n(mesh,data_e):
	if data_e.size==mesh.n2de:
		data_e=data_e*mesh.elem_2d_area[0:mesh.n2de]
	elif data_e.size==mesh.n2dea:
		data_e=data_e*mesh.elem_2d_area 
	
	data_n=np.zeros((mesh.n2dna,))
	for ii in range(mesh.n2de):     
		data_n[mesh.elem0_2d_i[ii,:]]=data_n[mesh.elem0_2d_i[ii,:]]+np.array([1,1,1])*data_e[ii] 
	data_n[0:mesh.n2dn]=data_n[0:mesh.n2dn]/mesh.nodes_2d_area[0:mesh.n2dn]/3. 
	data_n[mesh.n2dn:mesh.n2dna] = data_n[mesh.pbndn_2d_i]
	
	return data_n
	
##___CALCULATE INTERPOLATION MATRIX FROM ELEMENTS TO NODE POINTS RESOLUTION______
## interpolate from triangle values to node values
##
##_______________________________________________________________________________
#def fesom_e2n_interp(mesh):
	
	#print(' --> calc. interpolation matrix elem2d-->nodes')
	#outstep        		= 10.0;
	##___________________________________________________________________________
	#aux					= np.array(mesh.elem_2d_i)
	#aux					= np.squeeze(aux.reshape((aux.shape[0]*aux.shape[1],1)))
	
	##___________________________________________________________________________
	#sort_index 			= np.argsort(aux)
	#sort_aux			= aux[sort_index]
	#kk_s				= 0;
	#del aux
	
	###__________________________________________________________________________
	## calc max number of neighbouring elemts of a node points
	#daux 				= sort_aux[1:mesh.n2dea*3]-sort_aux[0:mesh.n2dea*3-1]
	#idaux 				= np.linspace(1,mesh.n2dea*3-1,mesh.n2dea*3-1,dtype='uint32')
	#jdaux 				= idaux[daux==1]
	#kk_max 				= np.max(jdaux[1:]-jdaux[:-1])
	#del daux
	#del idaux
	#del jdaux
	
	##___________________________________________________________________________
	## allocate interoplation matrices
	#interp_e2n_mat      = np.zeros((mesh.n2dna,10),'uint32')
	#interp_e2n_corector = np.zeros((mesh.n2dna,10),'uint8')
	
	##___________________________________________________________________________
	## for loop over surface nodes
	##print( '     |',end=''),
	#print( '     |'),
	#for ii in range(0, mesh.n2dna):
		##_______________________________________________________________________
		#if np.mod(ii,np.floor(  mesh.n2dna*(outstep/100) ))==0:
			##print('{:2.0f}%|'.format(ii/np.floor(mesh['n2dna']*(outstep/100))*outstep),end='')
			#print('{:2.0f}%|'.format(ii/np.floor(mesh.n2dna*(outstep/100))*outstep)),
		##_______________________________________________________________________
		#if (kk_s+kk_max < 3*mesh.n2dea):
			#sect_index = sort_index[kk_s:kk_s+kk_max]
			#sect_aux   = sort_aux[kk_s:kk_s+kk_max]
		#else:
			#sect_index = sort_index[kk_s:]
			#sect_aux   = sort_aux[kk_s:]
		##_______________________________________________________________________
		#elemofnode  = sect_index[sect_aux==ii]
		
		##_______________________________________________________________________
		#interp_e2n_mat[ii,0:elemofnode.shape[0]]=elemofnode
		
		##_______________________________________________________________________
		#kk_s = kk_s + elemofnode.shape[0]
		#del elemofnode
		#del sect_aux
		#del sect_index
		
	#print('')
	##___________________________________________________________________________
	#interp_e2n_corector[np.where(interp_e2n_mat!=0)]=1
	
	##___________________________________________________________________________
	#return(interp_e2n_mat,interp_e2n_corector) 
	
	
#___CALCULATE LAND MASK CONOURLINE______________________________________________
#
#
#_______________________________________________________________________________
def fesom_calc_landmask(mesh):
	
	np.set_printoptions(edgeitems=10)
	print(' --> calc landmask contourline')
	#___________________________________________________________________________
	# build land boundary edge matrix
	edge    = np.concatenate((mesh.elem_2d_i[:,[0,1]], \
	  										 mesh.elem_2d_i[:,[0,2]], \
											 mesh.elem_2d_i[:,[1,2]]),axis=0)
	edge    = np.sort(edge,axis=1) 
	
	# python  sortrows algorythm --> matlab equivalent
	edge    = edge.tolist()
	edge.sort()
	edge    = np.array(edge)
	
	idx     = np.diff(edge,axis=0)==0
	idx     = np.all(idx,axis=1)
	idx     = np.logical_or(np.concatenate((idx,np.array([False]))),np.concatenate((np.array([False]),idx)))
	
	bnde    = edge[idx==False,:]
	del edge 
	del idx
	
	#___________________________________________________________________________
	# kill periodic boundary edges--> periodic boundary edges have identical lon/lat points
	edge_x  = mesh.nodes_2d_xg[bnde]
	edge_y  = mesh.nodes_2d_yg[bnde]
	
	aux     = np.logical_or(edge_x==np.min(mesh.nodes_2d_xg),edge_x==np.max(mesh.nodes_2d_xg))
	aux     = np.logical_or(aux,edge_y>89)
	idx     = np.all(aux,axis=1)
	del aux
	del edge_x 
	del edge_y 
		
	pbnde   = bnde[idx==True ,:]
	bnde    = bnde[idx==False,:]
	del idx
	
	#___________________________________________________________________________
	#link correctly remaining points of periodic boundary
	lonly_p   = []
	aux_pbnde = np.sort(pbnde.reshape(pbnde.size))
	for ii in range(0,aux_pbnde.size):
		if np.sum(aux_pbnde==aux_pbnde[ii])==1:
			lonly_p.append(aux_pbnde[ii])
			
	lonly_p = np.array(lonly_p)
	sortvec = np.argsort(mesh.nodes_2d_yg[lonly_p])
	sortvec = np.array(sortvec[::-1]) #descent order
	lonly_p = lonly_p[sortvec]
	
	#___________________________________________________________________________
	addp_edge = [];
	xmin = np.min(mesh.nodes_2d_xg[lonly_p])
	xmax = np.max(mesh.nodes_2d_xg[lonly_p])
	
	#____START WHILE LOOP_______________________________________________________
	while lonly_p.size>0:
		
		#_______________________________________________________________________
		P1_i  = lonly_p[1];
		P1_x  = mesh.nodes_2d_xg[P1_i];
		P1_y  = mesh.nodes_2d_yg[P1_i];
		aux_i = np.squeeze(np.array(np.where(np.logical_and(mesh.nodes_2d_xg[lonly_p]==P1_x,lonly_p!=P1_i))))
		
		if np.size(aux_i)!=0:
			# find land mnask points on either the left or right
			# side of pbnd
			aux2_i = np.argmin(P1_y-mesh.nodes_2d_yg[lonly_p[aux_i]])
		else:
			# antarctica land mask points 
			if P1_x == xmin:
				aux_i = np.array(np.where(mesh.nodes_2d_xg[lonly_p]==xmax))
			else:
				aux_i = np.array(np.where(mesh.nodes_2d_xg[lonly_p]==xmin))
			
			aux2_i = np.array(np.argmin(np.abs(P1_y-mesh.nodes_2d_yg[lonly_p[aux_i]])))
			
		P2_i      = lonly_p[aux_i[aux2_i]]
		addp_edge.append([P1_i,P2_i])
		lonly_p   = np.array(lonly_p[np.where(lonly_p!=P1_i)])
		lonly_p   = np.array(lonly_p[np.where(lonly_p!=P2_i)])
		
	#____END WHILE LOOP_________________________________________________________
	addp_edge = np.array(addp_edge)
	addp_edge = np.sort(addp_edge,axis=1)
	bnde      = np.concatenate((bnde,addp_edge),axis=0)
	del aux_i; del aux2_i; del P1_i; del P1_x; del P1_y; del P2_i; del lonly_p; del addp_edge
	
	#___________________________________________________________________________
	# elimiate edges that go from -180 -->180 or 180--> -180
	edge_x  = np.sort(mesh.nodes_2d_xg[bnde],axis=1) ;
	idx     = (edge_x[:,1]-edge_x[:,0])>180;
	pbnde   = bnde[idx==True ,:];
	bnde    = bnde[idx==False,:];
	del edge_x; del idx
		
	#___________________________________________________________________________
	aux_x     = mesh.nodes_2d_xg[pbnde]
	aux_y     = mesh.nodes_2d_yg[pbnde]
	
	sortvec   = np.argsort(aux_y[:,0])
	aux_y     = aux_y[sortvec,:]
	aux_x     = aux_x[sortvec,:]
	aux_pbnde = pbnde[sortvec,:]
	
	sortvec = np.argsort(aux_x,axis=1)
	for ii in range(0,aux_y.shape[0]):
		aux_y[ii,:]     = aux_y[ii,sortvec[ii,:]]
		aux_pbnde[ii,:] = aux_pbnde[ii,sortvec[ii,:]]
		
	del sortvec
	
	#___________________________________________________________________________
	count_add_p  = np.array([mesh.n2dna])
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
			
	#___________________________________________________________________________
	# add new edge to all edge array
	add_edge = np.array(add_edge)
	add_edge = np.sort(add_edge,axis=1)
	bnde     = np.concatenate((bnde,add_edge),axis=0)
	nbnde    = bnde.shape[0];
			
	add_points_x = np.squeeze(np.array(add_points_x))
	add_points_y = np.squeeze(np.array(add_points_y))
	add_nodes_2d_xg = np.concatenate((mesh.nodes_2d_xg,add_points_x),axis=0)
	add_nodes_2d_yg = np.concatenate((mesh.nodes_2d_yg,add_points_y),axis=0)
		
	del y_m; del add_edge; del add_points_x; del add_points_y; del count_add_p
	del aux_x; del aux_y
	
	#___________________________________________________________________________
	run_cont        = np.zeros((1,nbnde))*np.nan
	run_cont[0,:2]  = bnde[0,:];
	run_bnde        = bnde[1:,:];
	count_init      = 1;
	init_ind        = run_cont[0,0];
	
	#land_cont_x     = np.zeros((nbnde*2,))*np.nan
	#land_cont_y     = np.zeros((nbnde*2,))*np.nan
	ind_lc_s        = 0;
	
	polygon_xy = []
	for ii in range(0,nbnde):
		#_______________________________________________________________________
		kk_rc = np.column_stack(np.where( run_bnde==np.int(run_cont[0,count_init]) ))
		kk_r  = kk_rc[:,0]
		kk_c  = kk_rc[:,1]
		count_init  = count_init+1
		
		#_______________________________________________________________________
		if kk_c[0] == 0 :
			run_cont[0,count_init] = run_bnde[kk_r[0],1]
		else:
			run_cont[0,count_init] = run_bnde[kk_r[0],0]
		
		#_______________________________________________________________________
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
			
			#___________________________________________________________________
			run_cont        = np.zeros((1,nbnde))*np.nan
			run_cont[0,:2]  = run_bnde[0,:]
			run_bnde        = run_bnde[1:,:]
			count_init=1;
			init_ind        = run_cont[0,0]
		else:
			aux_ind =np.arange(0,run_bnde.shape[0],1)
			run_bnde=run_bnde[aux_ind!=kk_r[0],:]
		
	#___________________________________________________________________________
	mesh.polygon_xy = polygon_xy
	
	#___________________________________________________________________________
	return mesh
	
	
#___ROTATE VECTOR GRID ROT-->GEO________________________________________________
#
#
#_______________________________________________________________________________
def fesom_vector_rot(mesh,u,v):
	
	print(' --> do vector rotation rot2geo')
	#___________________________________________________________________________
	alpha = mesh.alpha
	beta  = mesh.beta
	gamma = mesh.gamma
	
	#___________________________________________________________________________
	#mesh = fesom_grid_rot_g2r()
	
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
	if u.shape[0]==mesh.n2dna:
		rlat = np.radians(mesh.nodes_2d_yr)
		rlon = np.radians(mesh.nodes_2d_xr)
		lat  = np.radians(mesh.nodes_2d_yg)
		lon  = np.radians(mesh.nodes_2d_xg)
	elif u.shape[0]==mesh.n2dn:
		rlat = np.radians(mesh.nodes_2d_yr[0:mesh.n2dn])
		rlon = np.radians(mesh.nodes_2d_xr[0:mesh.n2dn])
		lat  = np.radians(mesh.nodes_2d_yg[0:mesh.n2dn])
		lon  = np.radians(mesh.nodes_2d_xg[0:mesh.n2dn])
	elif u.shape[0]==mesh.n2dea:
		rlat = np.radians(np.sum(mesh.nodes_2d_yr[mesh.elem_2d_i],axis=1)/3)
		rlon = np.radians(np.sum(mesh.nodes_2d_xr[mesh.elem_2d_i],axis=1)/3)
		lat  = np.radians(np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3)
		lon  = np.radians(np.sum(mesh.nodes_2d_xg[mesh.elem_2d_i],axis=1)/3)
	elif u.shape[0]==mesh.n2de:
		rlat = np.radians(np.sum(mesh.nodes_2d_yr[mesh.elem0_2d_i],axis=1)/3)
		rlon = np.radians(np.sum(mesh.nodes_2d_xr[mesh.elem0_2d_i],axis=1)/3)
		lat  = np.radians(np.sum(mesh.nodes_2d_yg[mesh.elem0_2d_i],axis=1)/3)
		lon  = np.radians(np.sum(mesh.nodes_2d_xg[mesh.elem0_2d_i],axis=1)/3)
		
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
	elif len(ndims)	== 2:
		
		#___________________________________________________________________________
		txg,tyg,tzg = np.zeros(u.shape), np.zeros(u.shape), np.zeros(u.shape)
		for ii in range(0,u.shape[1]):
			txg[:,ii] = -v[:,ii]*np.sin(rlat)*np.cos(rlon) - u[:,ii]*np.sin(rlon)
			tyg[:,ii] = -v[:,ii]*np.sin(rlat)*np.sin(rlon) + u[:,ii]*np.cos(rlon)
			tzg[:,ii] =  v[:,ii]*np.cos(rlat)
		
		txr = rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg;
		tyr = rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg;
		tzr = rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg;
		
		#___________________________________________________________________________
		for ii in range(0,u.shape[1]):
			v[:,ii] = txr[:,ii]*(-np.sin(lat))*np.cos(lon) - tyr[:,ii]*np.sin(lat)*np.sin(lon) + tzr[:,ii]*np.cos(lat)
			u[:,ii] = txr[:,ii]*(-np.sin(lon)) + tyr[:,ii]*np.cos(lon)
		
	
	return(u,v)
	
	
#___EQUIVALENT OF MATLAB ISMEMBER FUNCTION______________________________________
#............
#............
#_______________________________________________________________________________
def ismember_rows(a, b):
	return np.flatnonzero(np.in1d(b[:,0], a[:,0]) & np.in1d(b[:,1], a[:,1]))
	
	
#___TRANSFORM GEOCOORD to CARTESIAN COORD_______________________________________
#............
#............
#_______________________________________________________________________________
def geo2cart(glon,glat):
	Rearth=6371.0;
	x_cart=Rearth*np.cos(np.radians(glat))*np.cos(np.radians(glon))
	y_cart=Rearth*np.cos(np.radians(glat))*np.sin(np.radians(glon))
	z_cart=Rearth*np.sin(np.radians(glat))
	return np.array([x_cart,y_cart,z_cart])