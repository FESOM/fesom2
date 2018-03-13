# Patrick Scholz, 23.01.2018
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import cmocean
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import copy as cp

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray 		import *
from sub_fesom_mesh 		import * 
from sub_fesom_data 		import * 
from sub_fesom_plot 		import *
from sub_fesom_selectbox 	import *
from sub_regriding_adapt 	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
inputarray  = set_inputarray()
#inputarray['save_fig'] = True
inputarray['save_fig'] = False
inputarray['save_figpath'] = '/scratch/users/pscholz/AWI_PAPER/PAPER_FESOM2.0_evaluation/figures/withoutPC-2/zlevel-linfs/boxes/'

#+_____________________________________________________________________________+
# predefine section lines or leave it empty in this case interactive selection 
box 		      = fesom_box()
#box.box_define    =  list([[]]) # --> more boxes ,  list([[],[],...,...]) 
#box.box_define[0] = [ [-53.156, -49.688, -47.0, -49.719, -51.281, -53.719, -55.812, -53.375],
					  #[60.219, 59.656, 57.594, 53.531, 54.156, 55.406, 56.062, 60.438],
					   #'Test Polygon' ]
#box.box_define[0] = [ [-55.00, -47.00],
					 #[ 55.00,  60.00],
					   #'LabSea Test' ]
#box.box_define[1] = [ [-40.00, -30.00],
					 #[ 30.00,  40.00 ],
					   #'central Atl' ]
#+_____________________________________________________________________________+
box.descript,box.path = 'linfs', '/scratch/users/pscholz/AWI_DATA/fesom2.0_results/linfs/withoutPC-2/'
#box.descript,box.path = 'linfs', '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
box.var 	    = 'temp'
box.which_mean = 'None'

#box.crange     = []
#box.crange     = [2.5,6.0,4.5] # [cmin, cmax, cref]
#box.crange     = [-1.0,1.0,0.0] # [cmin, cmax, cref]
#box.cmap       = 'blue2red'
#box.cnumb      = 25
#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
#box.year	= [1990,2000]
box.year	= [1948,2009]

# select month to average over
box.month	= [1,2,3,4,5,6,7,8,9,10,11,12]

# select linear interpolated depth layers to average over, empty mean use all layers
# don't touch here, should be empty !
box.depth	= []

#+_____________________________________________________________________________+
# make anomaly
#do_anomaly      = True
do_anomaly      = False
if do_anomaly==True:
	box2  = cp.deepcopy(box)
	box2.descript,box2.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
	

#+_____________________________________________________________________________+
#|                         *** LOAD FVSOM MESH ***                             |
#+_____________________________________________________________________________+
#inputarray['which_plot']  = 'pcolor'
inputarray['which_plot']  = 'contourf'
try:
	mesh
except NameError:
	mesh = fesom_init_mesh(inputarray)
else:
	print " --> ___FOUND FESOM MESH --> will use it!___________________________"
	
#+_____________________________________________________________________________+
#|                    *** LINE SELECTION PROCESS ***                           |
#+_____________________________________________________________________________+
# --> interactive figure line selection 
#		[left mouse]  ... preselect area, move window
#		[right mouse] ... zoom out completly
#		button [+] ... zoom in
#		button [-] ... zoom out
#		button [1] ... first time push [1] - start line selection, line points can 
#					   now be selected with left mouse button, second time push of [1]
#					   finishes selection of a single line. If this is repeated 
#					   several lines can be selected
#		button [2] ... select polygon
#		button [d] ... delete point or rectangle that was just selected
#		button [q] ... finshes entire selection process, goes on with the calculation 
#					   of the cross-section
#	   
# --> BE AWARE !!!:  IN ORDER TO PROCEED IN MOMENT AFTER THE INTERACTIVE SELECTION
#					 ALL FIGURES THAT ARE STILL OPENED NEED TO BE CLOSED. OTHERWISE
#					 THE PROGRAMM CAN NOT PROCEED!
if len(box.box_define)==0:
	#___________________________________________________________________________
	fig, ax = plt.figure(figsize=(13, 13)), plt.gca()
	map 	= Basemap(projection = 'cyl',resolution = 'c',
				llcrnrlon = -180, urcrnrlon = 180, llcrnrlat = -90, urcrnrlat = 90)
	mx,my 	= map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
	#tri     = Triangulation(mx, my,mesh.elem_2d_i)
	#hp1		=plt.tripcolor(tri,-mesh.nodes_2d_z,cmap=cmocean.cm.deep)
	map.drawmapboundary(fill_color='0.9',linewidth=1.0)
	map.bluemarble()
	#map.etopo()
	fesom_plot_lmask(map,mesh,ax,'none','r')
	ax.grid(color='k', linestyle='-', linewidth=0.5)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.show(block=False)
	fig.canvas.draw()
	
	#___________________________________________________________________________
	# interactively
	#box._cid_pressb = fig.canvas.mpl_connect('button_press_event', box._anybutton_)
	#box._cid_pressk = fig.canvas.mpl_connect('key_press_event', box._anykey_)
	box._connect_(fig,ax,map)


#+_____________________________________________________________________________+
# analyse selected box, which points are inside box
print(' --> calculate points inside box')
box.select_pointsinbox(mesh)
if do_anomaly==True:
	box2.box_idx = box.box_idx
box.plot_index_position(mesh)

#+_____________________________________________________________________________+
#|                            *** LOAD 3D DATA ***                             |
#+_____________________________________________________________________________+
# load fesom2.0 data over time and depth --> average over points in box 
fesom_load_data_overtime(mesh,box)
if do_anomaly==True:
	fesom_load_data_overtime(mesh,box2)
	anom = cp.deepcopy(box)
	anom.data_anom(box,box2)


#+_____________________________________________________________________________+
#|               *** PLOT INDEXBOX DATA OVER DEPTH AND TIME ***                |
#+_____________________________________________________________________________+
if do_anomaly==False:
	fig2,ax,cax = box.plot_index_t_x_z()
if do_anomaly==True:
	#___________________________________________________________________________
	# do common crange for line and lin2
	cmax = np.max([np.nanmax(box.value),np.nanmax(box2.value)])
	cmin = np.min([np.nanmin(box.value),np.nanmin(box2.value)])
	cref = cmin + (cmax-cmin)/2
	cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
	box.crange=[cmin,cmax,cref]	
	box2.crange=[cmin,cmax,cref]
	#___________________________________________________________________________
	fig2,ax,cax = box.plot_index_t_x_z()
	fig2,ax,cax = box2.plot_index_t_x_z()
	fig2,ax,cax = anom.plot_index_t_x_z()
