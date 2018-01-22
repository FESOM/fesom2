import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import cmocean
import matplotlib.gridspec as gridspec
import copy as cp

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray 		import *
from sub_fesom_mesh 		import * 
from sub_fesom_data 		import * 
from sub_fesom_plot 		import *
from sub_fesom_selectline 	import *
from sub_regriding_adapt 	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
inputarray  = set_inputarray()
inputarray['save_fig'] = False
#inputarray['save_fig'] = True
inputarray['save_figpath'] = '/scratch/users/pscholz/AWI_PAPER/PAPER_FESOM2.0_evaluation/figures/withoutPC-2/zlevel-linfs/lines/'

#+_____________________________________________________________________________+
data 		= fesom_data(inputarray) # init fesom2.0 data object
data.descript,data.path = 'linfs' , '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
#data.descript,data.path = 'linfs' , '/scratch/users/pscholz/AWI_DATA/result_fvom_test/withPC-1/'
data.var 	= 'salt'  # 
#data.var 	= 'vec_uv'  # --> volume flux
#data.var 	= 'vec_tuv' # --> heat flux
#data.var 	= 'vec_suv' # --> liquid freshwater flux
#data.var 	= '.....'   # --> all other 3d fesom2.0 variable (temp, salt, u,v,w,...)
#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
data.year	= [1960,2009]
#data.year	= [2000,2009]
#data.year	= [1948,1948]

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
	data2 			= fesom_data(inputarray) # init fesom2.0 data object
	data2.descript,data2.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'
	data2.var, data2.year, data2.month, data2.depth = data.var, data.year, data.month, data.depth
	data2.crange, data2.cmap, data2.cnumb = data.crange, data.cmap, data.cnumb
	#data2.var, data2.year, data2.month, data2.depth = data.var, data.year, [6,7,8], data.depth


#+_____________________________________________________________________________+
#+_____________________________________________________________________________+
# predefine section lines or leave it empty in this case interactive selection 
line 		= fesom_line()

#line.line_define    =  list([[],[]]) 
line.line_define    =  list([[]])

#line.line_define[0] = [ [-75.0, -5.0],[42.0,42.0],'E-W' ]
#line.line_define[0] = [ [-75.0, -0.0],[40.0,50.0],'E-W' ]
#line.line_define[0] = [ [ -40.0,-50.0],[50.0,40.0],'N-S' ]
line.line_define[0] = [[-72.5, -64.625], [41.781, 32.469], 'GS']

#+_____________________________________________________________________________+
#|                         *** LOAD FVSOM MESH ***                             |
#+_____________________________________________________________________________+
inputarray['which_plot']  = 'pcolor'
#inputarray['which_plot']  = 'contourf'
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
#hp1		= plt.tripcolor(tri,-mesh.nodes_2d_z,cmap=cmocean.cm.deep)
#hp1		=plt.tripcolor(tri,data_anom.value,cmap=cmocean.cm.balance)
#plt.clim(-0.50,0.50)

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
#		button [L] ... first time push [l] - start line selection, line points can 
#					   now be selected with left mouse button, second time push of [L]
#					   finishes selection of a single line. If this is repeated 
#					   several lines can be selected
#		button [d] ... delete point that was just selected
#		button [q] ... finshes entire selection process, goes on with the calculation 
#					   of the cross-section
#	   
# you need to close window to proceed
if len(line.line_define)==0:
	# interactively
	line._cid_pressb = fig.canvas.mpl_connect('button_press_event', line._anybutton_)
	line._cid_pressk = fig.canvas.mpl_connect('key_press_event', line._anykey_)
	line._connect_(fig,ax,map)
	
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
	for ii in range(0,len(line.line_define)):
		ax.plot(line.line_define[ii][0],line.line_define[ii][1] ,color='w',linewidth=2,marker='o',mfc='w',mec='k',axes=ax)
		ax.plot(line.line_define[ii][0][0],line.line_define[ii][1][0] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.0,0.8,0.0],mec='k',axes=ax)
		ax.plot(line.line_define[ii][0][-1],line.line_define[ii][1][-1] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.8,0.0,0.0],mec='k',axes=ax)
		xmax = np.max([xmax,np.max(line.line_define[ii][0])])
		xmin = np.min([xmin,np.min(line.line_define[ii][0])])
		ymax = np.max([ymax,np.max(line.line_define[ii][1])])
		ymin = np.min([ymin,np.min(line.line_define[ii][1])])
	ax.set_xlim(xmin-10.0,xmax+10.0)
	ax.set_ylim(ymin-10.0,ymax+10.0)
	plt.show(block=False)
	fig.canvas.draw()
else:
	# predefined
	xmax,xmin,ymax,ymin = -180.0,180.0,-90.0,90.0
	for ii in range(0,len(line.line_define)):
		ax.plot(line.line_define[ii][0],line.line_define[ii][1] ,color='w',linewidth=2,marker='o',mfc='w',mec='k',axes=ax)
		ax.plot(line.line_define[ii][0][0],line.line_define[ii][1][0] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.0,0.8,0.0],mec='k',axes=ax)
		ax.plot(line.line_define[ii][0][-1],line.line_define[ii][1][-1] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.8,0.0,0.0],mec='k',axes=ax)
		xmax = np.max([xmax,np.max(line.line_define[ii][0])])
		xmin = np.min([xmin,np.min(line.line_define[ii][0])])
		ymax = np.max([ymax,np.max(line.line_define[ii][1])])
		ymin = np.min([ymin,np.min(line.line_define[ii][1])])
	ax.set_xlim(xmin-10.0,xmax+10.0)
	ax.set_ylim(ymin-10.0,ymax+10.0)
	plt.show(block=False)
	fig.canvas.draw()

#+_____________________________________________________________________________+
# analyse selected lines, calculate interpolation points
print(' --> calculate interpolation points')
# which   ... res/npoints
# res     ... interpolation resolution in deg
# npoints ... use npoints per line segment for interpolation
line.analyse_lines(which='res',res=0.1)

for ii in range(0,len(line.line_define)):
	ax.plot(line.line_interp_pm[ii][0],line.line_interp_pm[ii][1] ,
			color='r',linewidth=4,marker='+',linestyle='None',axes=ax,markersize=6,mfc='r')
	step=10
	ax.quiver(line.line_interp_pm[ii][0][0::step],line.line_interp_pm[ii][1][0::step],
			  line.line_interp_pm_nvec[ii][0][0::step],line.line_interp_pm_nvec[ii][1][0::step],color='r')
plt.show(block=False)
fig.canvas.draw()


#+_____________________________________________________________________________+
#|                  *** LOAD AND INTERPOLATE 3D DATA ***                       |
#+_____________________________________________________________________________+
# load 3d fesom2.0 data on original depth layers --> data.value
print(' --> load 3d fesom2.0 data')
data=fesom_load_data_horiz(mesh,data)
if do_anomaly==True:
	data2=fesom_load_data_horiz(mesh,data2)
	#if data.var.find('vec')==-1: anom = fesom_data_anom(data,data2)


#+_____________________________________________________________________________+
# interpolate points on fesom2.0 data, use routine of Nikolay fesom2regular(...)
print(' --> do interpolation on data')
if do_anomaly==False:
	line.interp_lines(mesh,data)
	if data.var.find('vec')!=-1: line.calc_flux()
else:
	if data.var.find('vec')!=-1:
		
		# calculate data on line and line2 data
		line.interp_lines(mesh,data)
		line2 = cp.deepcopy(line)
		line2.interp_lines(mesh,data2)
		
		# calculate flux through line and line2
		line.calc_flux()
		line2.calc_flux()
		
		# calculate anomaly line2-line ... 
		#anom  = cp.deepcopy(line)
		anom  = cp.deepcopy(line)
		anom.data_anom(line,line2)
	else:
		line.interp_lines(mesh,data)
		
		line2 = cp.deepcopy(line)
		line2.interp_lines(mesh,data2)
		
		anom  = cp.deepcopy(line)
		anom.data_anom(line,line2)


#+_____________________________________________________________________________+
#|                  *** PLOT LINE OVER DEPTH ***                               |
#+_____________________________________________________________________________+
print(' --> plot line data')
if  do_anomaly==True and data.var.find('vec')!=-1:
	line.plot_lines_dist_x_z()
	line2.plot_lines_dist_x_z()
	anom.plot_lines_dist_x_z()
else:
	
	if do_anomaly==True:
		# do common crange for line and lin2
		cmax = np.max([np.nanmax(line.value[0]),np.nanmax(line2.value[0])])
		cmin = np.min([np.nanmin(line.value[0]),np.nanmin(line2.value[0])])
		cref = cmin + (cmax-cmin)/2
		cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		line.crange=[cmin,cmax,cref]	
		line2.crange=[cmin,cmax,cref]
		
		line.plot_lines_dist_x_z()
		line2.plot_lines_dist_x_z()
		anom.plot_lines_dist_x_z()
	else:
		line.plot_lines_dist_x_z()
