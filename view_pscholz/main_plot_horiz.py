# Patrick Scholz, 23.01.2018
#+____IMPORT RELATED LIBRARYS__________________________________________________+
#get_ipython().magic('matplotlib notebook')
#get_ipython().magic('load_ext autoreload')
#get_ipython().magic('autoreload 2')
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import copy as cp

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray  import *
from sub_fesom_mesh  import * 
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
inputarray['save_fig'] = False
#inputarray['save_fig'] = True
inputarray['save_figpath'] = '/scratch/users/pscholz/AWI_PAPER/PAPER_FESOM2.0_evaluation/figures/withoutPC-2/zlevel-linfs/'

# set plot box for cyl projection (default: [-180,180,-90,90])
#inputarray['which_box'] = [0,360,-90,90] ; inputarray['mesh_focus']=180
inputarray['which_box'] = [-180,180,-90,90]
#inputarray['which_box'] = [-90,35,20,85]
#inputarray['which_box'] = [-75,0,30,55]
#inputarray['which_box'] = [-180,180,35,90]
#inputarray['which_box'] = [-180,180,-90,-50]

# set projection variable --> the lon/lat projection plot ranges are set via 
# inputarray['which_box'] = [lonmin,lonmax,latmin,latmax]
inputarray['proj'     ] = 'cyl' # 'ortho', 'cyl', 'npstere' , 'spstere'
inputarray['proj_lon' ] = -45 #only for ortho
inputarray['proj_lat' ] = 45 #only for ortho

#+_____________________________________________________________________________+
# setup variable name, runid and data path
data 			= fesom_data(inputarray) # init fesom2.0 data object
data.descript,data.path = 'linfs' , '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
#data.descript,data.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
#data.descript,data.path = 'test'  ,'/scratch/users/pscholz/AWI_DATA/result_fvom_test/withPC-1/'
data.var 		= 'triresol'

#data.crange     = []
#data.crange     = [0.0,5.0,2.5] # [cmin, cmax, cref]
#data.cmap       = 'rygbw'
#data.cnumb      = 25

#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
data.year		= [1960,2000]
#data.year		= [1990,2000]
#data.year		= [2000,2000]

# select month to average over
data.month		= [1,2,3,4,5,6,7,8,9,10,11,12]
#data.month		= [1,2,12]
#data.month		= [6,7,8]
#data.month		= [3]

# select linear interpolated depth layers to average over
#data.depth		= [0,10,20,30,40,50,75,100]
data.depth		= np.arange(   0, 200+1,10)
#data.depth		= np.arange( 200, 500+1,20)
#data.depth		= np.arange( 500,1000+1,50)
#data.depth		= np.arange(1000,1500+1,50)
#data.depth		= [50]

#+_____________________________________________________________________________+
# make anomaly
#do_anomaly      = True
do_anomaly      = False
if do_anomaly==True:
	data2 			= cp.copy(data) # init fesom2.0 data object
	data2.descript,data2.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
	#data2.descript,data2.path = 'zstar' ,'/scratch/users/pscholz/ollie_work/results_ALE/zstar/withoutPC-2/'
	
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
# plot topography
if data.var=='depth':
	data.value 	= -mesh.nodes_2d_z
	#data.value 	= mesh.nodes_2d_iz
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
	if len(mesh.nodes_2d_resol)==0: mesh.fesom_calc_triresol()
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
	if len(mesh.nodes_2d_area)==0: mesh.fesom_calc_triarea()
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
	
	if do_anomaly==True:
		fesom_load_data_horiz(mesh,data2)
		anom = fesom_data_anom(data,data2)
	
	
#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** PLOT FVSOM DATA ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
# plot 2d and 2dvec data
if len(data.value2)==0:
	
	#+_________________________________________________________________________+
	# plot anomaly
	if do_anomaly==False:
		#_______________________________________________________________________
		if data.value.size == mesh.n2dea: data.value = mesh.fesom_interp_e2n(np.array(data.value))
		fig,ax,map,cbar=fesom_plot2d_data(mesh,data)
	else:
		#_______________________________________________________________________
		# interpolate elemental values to nodal vlaues
		if data.value.size  == mesh.n2dea: data.value  = mesh.fesom_interp_e2n(np.array(data.value))
		if data2.value.size == mesh.n2dea: data2.value = mesh.fesom_interp_e2n(np.array(data2.value))
		if anom.value.size  == mesh.n2dea: anom.value  = mesh.fesom_interp_e2n(np.array(anom.value))
		
		#_______________________________________________________________________
		# do common crange for line and lin2
		cmax = np.max([np.nanmax(data.value),np.nanmax(data2.value)])
		cmin = np.min([np.nanmin(data.value),np.nanmin(data2.value)])
		cref = cmin + (cmax-cmin)/2
		cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		data.crange = [cmin,cmax,cref]	
		data2.crange= [cmin,cmax,cref]
		
		#_______________________________________________________________________
		fig,ax,map,cbar=fesom_plot2d_data(mesh,data)
		fig,ax,map,cbar=fesom_plot2d_data(mesh,data2)
		#anom.crange=[-0.25,0.25,0.0]
		fig,ax,map,cbar=fesom_plot2d_data(mesh,anom)
		
else:
	
	fesom_plot2dvec_data(mesh,data)
	
###fesom_plot3d_earth(mesh,data)
###fesom_plot2d_geomesh(mesh)
###fesom_plot2d_rotmesh(mesh)
