# Patrick Scholz, 23.01.2018
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
from sub_climatology  import *
from colormap_c2c    import *

#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** SET INPUT PARAMETER ***                         |
#|                                                                             |
#+_____________________________________________________________________________+
inputarray=set_inputarray()
inputarray['save_fig'] = False
#inputarray['save_fig'] = True
inputarray['save_figpath'] = '/scratch/users/pscholz/AWI_PAPER/PAPER_FESOM2.0_evaluation/figures/withoutPC-2/linfs-clim/'

# set plot box for cyl projection (default: [-180,180,-90,90])
#inputarray['which_box'] = [-180,180,-90,90]
#inputarray['which_box'] = [-90,35,20,85]
inputarray['which_box'] = [-180,180,15,90]
#inputarray['which_box'] = [-180,180,-90,-50]

# set projection variable --> the lon/lat projection plot ranges are set via 
# inputarray['which_box'] = [lonmin,lonmax,latmin,latmax]
inputarray['proj'     ] = 'npstere' # 'ortho', 'cyl', 'npstere' , 'spstere'
inputarray['proj_lon' ] = -45 #only for ortho
inputarray['proj_lat' ] = 45 #only for ortho

#+_____________________________________________________________________________+
# setup variable name, runid and data path
#data			= fesom_data()
data 			= fesom_data(inputarray) # init fesom2.0 data object
data.descript,data.path = 'linfs', '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
#data.descript,data.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
#data.path,data.descript = '/scratch/users/pscholz/AWI_DATA/result_fvom_test/withPC-1/' , '..-^-..-^-..-><(((*>'
data.var 		= 'temp'

#data.crange     = []
data.crange     = [-7.0,7.0,0.0] # [cmin, cmax, cref]
#data.crange     = [-2.0,2.0,0.0] # [cmin, cmax, cref]
#data.crange     = [-1.0,1.0,0.0] # [cmin, cmax, cref]
data.cmap       = 'blue2red'
#data.cnumb      = 25

#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
data.year		= [1998,2007]
#data.year		= [1948,1948]

# select month to average over
data.month		= [1,2,3,4,5,6,7,8,9,10,11,12]
#data.month		= [1,2,12]
#data.month		= [9]

# select linear interpolated depth layers to average over
#data.depth		= []
#data.depth		= np.arange(   0, 200+1,10)
#data.depth		= np.arange( 200, 500+1,20)
#data.depth		= np.arange( 500,1000+1,50)
#data.depth		= np.arange(1000,1500+1,50)
data.depth		= [300]
	
	
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
fesom_load_data_horiz(mesh,data)
#fig,ax,map,cbar=fesom_plot2d_data(mesh,data)

#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** LOAD CLIMATOLOGY ***                            |
#|                                                                             |
#+_____________________________________________________________________________+
path  = '/scratch/users/pscholz/AWI_FESOM2.0_git/view/woa2005/'
fname = 'woa2005TS.nc'
clim  = clim_data(path,fname,data.var)
clim.descript = data.descript+'-'+clim.descript
clim.crange,clim.cmap, clim.cnumb      = data.crange,data.cmap, data.cnumb
clim.str_time, clim.str_dep            = data.str_time, data.str_dep
clim.sname, clim.lname, clim.unit      = data.sname, data.lname, data.unit
clim.proj, clim.proj_lon, clim.proj_lat= data.proj, data.proj_lon, data.proj_lat

#+_____________________________________________________________________________+
#| interpolate and average clim data verticaly                                          
#+_____________________________________________________________________________+
data_clim=clim_vinterp(clim,data.depth)

#+_____________________________________________________________________________+
#| interpolate fesom data horizontalyto climatology grid                       |                                       
#+_____________________________________________________________________________+
from sub_regriding_adapt import *
data_fesom      = np.zeros(data_clim.shape)
mlon,mlat       = np.meshgrid(clim.lon, clim.lat)
distances, inds = create_indexes_and_distances(mesh, mlon, mlat,k=10, n_jobs=2)
data_fesom      = fesom2regular(data.value, mesh, mlon, mlat, distances=distances, inds=inds, radius_of_influence=100000)
data_fesom      = data_fesom.data
data_fesom[np.isnan(data_clim)]=np.nan

    
#+_____________________________________________________________________________+
#| calc anomaly between fesom and climatology data                             |                                       
#+_____________________________________________________________________________+
clim.anom = data_fesom-data_clim
#clim.anom = data_fesom
#clim.anom = data_clim
#clim.anom = data_clim-clim.value[1,:,:]

#+_____________________________________________________________________________+
#|                                                                             |
#|                         *** PLOT FVSOM DATA ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
fig,ax,map,cbar=clim_plot_anom(clim)
