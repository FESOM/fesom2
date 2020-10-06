# Patrick Scholz, 14.12.2017
def set_inputarray():
	
	global inputarray
	
	import numpy as np
	
	inputarray = dict()
	#____FESOM RUN INFO_________________________________________________________
	inputarray['data_id'		 ] = 'fesom'
	# inputarray['data_dir1'		 ] = '/data_big/result_fvom_test/'
	#inputarray['data_dir1'		 ] = '/media/pscholz/data_ext_bckp/data_big/data_big/result_fvom_test/'
	inputarray['data_dir1'		 ] = '/scratch/users/pscholz/AWI_DATA/result_fvom_test/withPC-1/'    
	#inputarray['data_dir1'		 ] = '/scratch/users/pscholz/AWI_DATA/result_fvom_test/'
	#inputarray['data_dir1'		 ] = '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'    
	
	#____FESOM MESH INFO________________________________________________________
# 	inputarray['mesh_id'],inputarray['mesh_dir'] = 'COREv2','/work/ollie/pscholz/mesh_fesom2.0/mesh_CORE2_dsidorenko_meanval/'
	inputarray['mesh_id'],inputarray['mesh_dir'] = 'COREv2','/work/ollie/pscholz/mesh_fesom2.0/core2_meanz/'
	
	# Euler-Angle:for poles at (phi,theta)=(-40,75)
	inputarray['mesh_alpha'	 	 ] = np.float( 50.0)
	inputarray['mesh_beta'		 ] = np.float( 15.0)
	inputarray['mesh_gamma'	 	 ] = np.float(-90.0)
	inputarray['mesh_rotate'	 ] = True
	inputarray['mesh_remove_cyc' ] = True
	
	# Change mesh focus from atlantic to pacific values between [-180...180]
	inputarray['mesh_focus'	 	 ] = 0
	
	inputarray['which_mask'	 	 ] = 'fesom' # 'fesom', 'bluemarble', 'etopo'
	#inputarray['which_mask'	 	 ] = 'bluemarble' # 'fesom', 'bluemarble', 'etopo'
	#inputarray['which_mask'	 	 ] = 'etopo' # 'fesom', 'bluemarble', 'etopo'
	
	inputarray['save_fig'        ] = False
	inputarray['save_figpath'    ] = ''
	
	#inputarray['which_plot'     ]  = 'pcolor' #contourf/pcolor
	inputarray['which_plot'     ]  = 'contourf' #contourf/pcolor
	
	inputarray['which_box']        = [-180,180,-90,90]
	inputarray['proj'     ] 	   = 'cyl' # 'ortho', 'cyl', 'npstere'  
	inputarray['proj_lon' ] 	   = 0 #only for ortho
	inputarray['proj_lat' ] 	   = 75 #only for ortho
	
	inputarray['use_cavity' ] 	   = False #only for ortho
	#___________________________________________________________________________
	return inputarray

