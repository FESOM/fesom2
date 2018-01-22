import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import cmocean
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray 		import *
from sub_fesom_mesh 		import * 
from sub_fesom_data 		import * 
from sub_fesom_plot 		import *
from sub_fesom_selectbox 	import *
from sub_regriding_adapt 	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
#+_____________________________________________________________________________+
inputarray  = set_inputarray()
data 		= fesom_data(inputarray) # init fesom2.0 data object

#data.descript,data.path = 'linfs', '/scratch/users/pscholz/AWI_DATA/fesom2.0_results/linfs/withoutPC-2/'
data.descript,data.path = 'linfs', '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
data.var 	    = 'temp'
data.which_mean = 'None'

#data.crange     = []
#data.crange     = [2.5,6.0,4.5] # [cmin, cmax, cref]
#data.crange     = [-1.0,1.0,0.0] # [cmin, cmax, cref]
#data.cmap       = 'blue2red'
#data.cnumb      = 25
#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
#data.year	= [1990,2000]
data.year	= [1948,2009]

# select month to average over
data.month	= [1,2,3,4,5,6,7,8,9,10,11,12]

# select linear interpolated depth layers to average over, empty mean use all layers
# don't touch here, should be empty !
data.depth	= []

#+_____________________________________________________________________________+
# make anomaly
do_anomaly      = True
#do_anomaly      = False
if do_anomaly==True:
	data2 = fesom_data(inputarray) # init fesom2.0 data object
	data2.descript,data2.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
	data2.var, data2.year, data2.month, data2.depth = data.var, data.year, data.month, data.depth
	data2.crange, data2.cmap, data2.cnumb = data.crange, data.cmap, data.cnumb
	data2.which_mean = data.which_mean
	

#+_____________________________________________________________________________+
#+_____________________________________________________________________________+
# predefine section lines or leave it empty in this case interactive selection 
box 		     = fesom_init_box()
box.box_refxy    =  list([[]]) 
box.box_refxy[0] = [ [-55.00, -47.00, -47.00, -55.00],
					 [ 55.00,  55.00,  60.00,  60.00],
					   'LabSea Test' ]

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

#+_____________________________________________________________________________+
# --> interactive figure line selection 
#		[left mouse]  ... preselect area, move window
#		[right mouse] ... zoom out completly
#		button [+] ... zoom in
#		button [-] ... zoom out
#		button [b] ... first time push [l] - start line selection, line points can 
#					   now be selected with left mouse button, second time push of [L]
#					   finishes selection of a single line. If this is repeated 
#					   several lines can be selected
#		button [d] ... delete point or rectangle that was just selected
#		button [q] ... finshes entire selection process, goes on with the calculation 
#					   of the cross-section
#	   
# you need to close window to proceed
if len(box.box_refxy)==0:
	# interactively
	box.cid_pressb = fig.canvas.mpl_connect('button_press_event', box.anybutton)
	box.cid_pressk = fig.canvas.mpl_connect('key_press_event', box.anykey)
	box.connect(fig,ax,map)
	
	# redraw map
	fig, ax = plt.figure(figsize=(13, 13)), plt.gca()
	map 	= Basemap(projection = 'cyl',resolution = 'c',
				llcrnrlon = -180, urcrnrlon = 180, llcrnrlat = -90, urcrnrlat = 90)
	mx,my 	= map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
	map.drawmapboundary(fill_color='0.9',linewidth=1.0)
	map.bluemarble()
	fesom_plot_lmask(map,mesh,ax,'none','r')
	ax.grid(color='k', linestyle='-', linewidth=0.5)
	xmax,xmin,ymax,ymin = -180.0,180.0,-90.0,90.0
	patch=[]
	for ii in range(0,len(box.box_refxy)):
		patch.append(Polygon(zip(box.box_refxy[ii][0],box.box_refxy[ii][1]),closed=True,clip_on=True) )
		ax.plot(box.box_refxy[ii][0]    ,box.box_refxy[ii][1] ,linestyle='None'   ,color='w',linewidth=2.0,marker='o',mfc='w',mec='k',axes=ax)
		xmax = np.max([xmax,np.max(box.box_refxy[ii][0])])
		xmin = np.min([xmin,np.min(box.box_refxy[ii][0])])
		ymax = np.max([ymax,np.max(box.box_refxy[ii][1])])
		ymin = np.min([ymin,np.min(box.box_refxy[ii][1])])
	ax.add_collection(PatchCollection(patch, alpha=1.0,facecolor='none',edgecolor='w',zorder=1,linewidth=2.0,hatch='/'))	
	ax.set_xlim(np.max([-180,xmin-20.0]),np.min([180,xmax+20.0]))
	ax.set_ylim(np.max([ -90,ymin-20.0]),np.min([ 90,ymax+20.0]))
	ax.grid(color='k', linestyle='-', linewidth=0.5)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.show(block=False)
	fig.canvas.draw()
else:
	# predefined
	xmax,xmin,ymax,ymin = -180.0,180.0,-90.0,90.0
	patch=[]
	for ii in range(0,len(box.box_refxy)):
		patch.append(Polygon(zip(box.box_refxy[ii][0],box.box_refxy[ii][1]),closed=True,clip_on=True) )
		ax.plot(box.box_refxy[ii][0]    ,box.box_refxy[ii][1] ,linestyle='None'   ,color='w',linewidth=2.0,marker='o',mfc='w',mec='k',axes=ax)
		xmax = np.max([xmax,np.max(box.box_refxy[ii][0])])
		xmin = np.min([xmin,np.min(box.box_refxy[ii][0])])
		ymax = np.max([ymax,np.max(box.box_refxy[ii][1])])
		ymin = np.min([ymin,np.min(box.box_refxy[ii][1])])
	ax.add_collection(PatchCollection(patch, alpha=1.0,facecolor='none',edgecolor='w',zorder=1,linewidth=2.0,hatch='/'))
	ax.set_xlim(np.max([-180,xmin-20.0]),np.min([180,xmax+20.0]))
	ax.set_ylim(np.max([ -90,ymin-20.0]),np.min([ 90,ymax+20.0]))
	ax.grid(color='k', linestyle='-', linewidth=0.5)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.show(block=False)
	fig.canvas.draw()


#+_____________________________________________________________________________+
# analyse selected box, which points are inside box
print(' --> calculate points inside box')
box.select_pointsinbox(mesh)


#+_____________________________________________________________________________+
# load fesom2.0 data over time and depth --> average over points in box 
data=fesom_load_data_overtime(mesh,data,box.box_idx[0])
if do_anomaly==True:
	data2=fesom_load_data_overtime(mesh,data2,box.box_idx[0])
	anom = fesom_data_anom(data,data2)

#+_____________________________________________________________________________+
#|               *** PLOT INDEXBOX DATA OVER DEPTH AND TIME ***                |
#+_____________________________________________________________________________+
if do_anomaly==False:
	fig,ax,cax = plot_index_t_x_z(mesh,data,box)
if do_anomaly==True:	
	#fig,ax,cax = plot_index_t_x_z(mesh,data,box)
	#fig,ax,cax = plot_index_t_x_z(mesh,data2,box)
	fig,ax,cax = plot_index_t_x_z(mesh,anom,box)
