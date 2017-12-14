# Patrick Scholz, 14.12.2017
#+____IMPORT RELATED LIBRARYS__________________________________________________+
#get_ipython().magic('matplotlib notebook')
#get_ipython().magic('load_ext autoreload')
#get_ipython().magic('autoreload 2')
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray  import *
from sub_fesom_mesh  import * 
from sub_fesom_data  import * 
from sub_fesom_plot  import *
from colormap_c2c    import *

#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** SET INPUT PARAMETER ***                         |
#|                                                                             |
#+_____________________________________________________________________________+
inputarray=set_inputarray()
#inputarray['save_fig'] = False

# set plot box for cyl projection (default: [-180,180,-90,90])
inputarray['which_box'] = [-180,180,50,90]

# set projection variable --> the lon/lat projection plot ranges are set via 
# inputarray['which_box'] = [lonmin,lonmax,latmin,latmax]
inputarray['proj'     ] = 'npstere' # 'ortho', 'cyl', 'npstere' , 'spstere'
inputarray['proj_lon' ] = 0 #only for ortho
inputarray['proj_lat' ] = 75 #only for ortho

#+_____________________________________________________________________________+
# setup variable name, runid and data path
data 			= fesom_init_data(inputarray) # init fesom2.0 data object
data.var 		= 'a_ice'
#data.crange     = []
#data.crange     = [0.0,5000.0,1000.0] # [cmin, cmax, cref]
#data.cmap       ='wbgyr'

#+_____________________________________________________________________________+
# select year, records, and depth range of which should be average
data.year		= [1948,1948]
# selct month to average over
#data.month		= [1,2,3,4,5,6,7,8,9,10,11,12]
#data.month		= [1,2,12]
data.month		= [9]
# select interpolate depth layers to average over
data.depth		= [0,10,20,30,40,50,75,100]

#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** LOAD FVSOM MESH ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
try:
	mesh
except NameError:
	mesh = fesom_init_mesh(inputarray)
else:
	print " --> ___FOUND FESOM MESH --> will use it!___________________________"
	
#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** LOAD FVSOM DATA ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
#_______________________________________________________________________________
# plot topography
if data.var=='depth':
	data.value 	= -mesh.nodes_2d_z
	data.sname 	= 'depth'
	data.lname 	= 'Depth'
	data.unit 	= 'm'
	data.levels = np.arange(0,np.max(data.value),200)# set resolution levels
	#data.cmap 	= 'blue2red'
	#data.cmap 	= 'cmocean.cm.balance'
	data.cmap 	= 'wbgyr'
#_______________________________________________________________________________
# plot triangle resolution interpolated to node
elif data.var=='triresol':
	mesh = fesom_calc_triresol(mesh)
	data.value 	= mesh.nodes_2d_resol
	data.sname 	= 'triresol'
	data.lname 	= 'Resolution'
	data.unit 	= 'km'
	data.levels = np.arange(0,np.max(data.value),10)# set resolution levels
	#data.cmap 	= 'blue2red'
	#data.cmap 	= 'cmocean.cm.balance'
	data.cmap 	= 'odv'
#_______________________________________________________________________________
# plot triangle area interpolated to node
elif data.var=='triarea':
	mesh = fesom_calc_triarea(mesh)
	data.value 	= mesh.nodes_2d_area
	data.sname 	= 'triarea'
	data.lname 	= 'Area'
	data.unit 	= 'km^2'
	data.levels = np.arange(0,np.max(data.value),10)# set resolution levels
	#data.cmap 	= 'blue2red'
	data.cmap 	= 'cmocean.cm.balance'
	#data.cmap 	= 'odv'
#_______________________________________________________________________________
# load all other 2d and 3d variables
else:
	fesom_load_data_horiz(mesh,data)
	
#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** PLOT FVSOM DATA ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
# plot 2d and 2dvec data
if len(data.value2)==0:
	test=fesom_plot2d_data(mesh,data)
else:
	fesom_plot2dvec_data(mesh,data)
	
###fesom_plot3d_earth(mesh,data)
###fesom_plot2d_geomesh(mesh)
###fesom_plot2d_rotmesh(mesh)
