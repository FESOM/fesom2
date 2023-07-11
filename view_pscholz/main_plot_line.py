# Patrick Scholz, 23.01.2018
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
from sub_fesom_selectline	import *
from sub_regriding_adapt	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
inputarray  = set_inputarray()
inputarray['save_fig'] = False
#inputarray['save_fig'] = True
inputarray['save_figpath'] = '/scratch/users/pscholz/AWI_PAPER/PAPER_FESOM2.0_evaluation/figures/withoutPC-2/zlevel-linfs/lines/'

#+_____________________________________________________________________________+
# predefine section lines or leave it empty in this case interactive selection 
line				= fesom_line()
line.line_define    =  list([[]]) # --> more line ,  list([[],[],...,...]) 
#line.line_define[0]    = [[-81.75, -1.3398], [37.969, 32.156], 'default']
#line.line_define[0] = [[-76.5, -64.625], [47.781, 32.469], 'GS']
#line.line_define[0] = [[170.0, 130.0], [-30.0, 30.0], 'Test']
#line.line_define[0] = [[-72.5, -64.625,-58.5], [43.781, 32.469, 28.5], 'GS']
line.line_define[0] = [[-72.5, -64.625], [43.781, 32.469], 'GS']
#+_____________________________________________________________________________+
line.descript,line.path = 'linfs' , '/media/pscholz/data_ext_2/DATA_FESOM2.0/linfs/withoutPC-2/'
#line.descript,line.path = 'linfs' , '/scratch/users/pscholz/AWI_DATA/fesom2.0_results/linfs/withoutPC-2/'
line.var 	= 'vec_uv'  # 
#line.var 	= 'vec_uv'  # --> volume flux
#line.var 	= 'vec_tuv' # --> heat flux
#line.var 	= 'vec_suv' # --> liquid freshwater flux (not tested yet)
#line.var 	= '.....'   # --> all other 3d fesom2.0 variable (temp, salt, u,v,w,...)
#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
#line.year	= [1960,2009]
line.year	= [2000,2009]
#line.year	= [1948,1948]

# select month to average over
line.month	= [1,2,3,4,5,6,7,8,9,10,11,12]

# select linear interpolated depth layers to average over, empty mean use all layers
# don't touch here, should be empty !
line.depth	= []

#+_____________________________________________________________________________+
# make anomaly
#do_anomaly      = True
do_anomaly      = False
if do_anomaly==True:
	line2 			= cp.copy(line) # init fesom2.0 data object
	line2.descript,line2.path = 'zlevel','/media/pscholz/data_ext_2/DATA_FESOM2.0/zlevel/withoutPC-2/'


#+_____________________________________________________________________________+
#|                         *** LOAD FVSOM MESH ***                             |
#+_____________________________________________________________________________+
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
	#___________________________________________________________________________
	fig, ax = plt.figure(figsize=(13, 13)), plt.gca()
	map 	= Basemap(projection = 'cyl',resolution = 'c',
				llcrnrlon = -180, urcrnrlon = 180, llcrnrlat = -90, urcrnrlat = 90)
	mx,my 	= map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
	
	#___________________________________________________________________________
	#tri     = Triangulation(mx, my,mesh.elem_2d_i)
	#hp1		= plt.tripcolor(tri,-mesh.nodes_2d_z,cmap=cmocean.cm.deep)
	#hp1		=plt.tripcolor(tri,data_anom.value,cmap=cmocean.cm.balance)
	#plt.clim(-0.50,0.50)
	
	#___________________________________________________________________________
	map.drawmapboundary(fill_color='0.9',linewidth=1.0)
	map.bluemarble(scale=0.25)
	#map.etopo()
	fesom_plot_lmask(map,mesh,ax,'none','r')
	ax.grid(color='k', linestyle='-', linewidth=0.5)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.show(block=False)
	fig.canvas.draw()
	
	#___________________________________________________________________________
	# interactively
	line._cid_pressb = fig.canvas.mpl_connect('button_press_event', line._anybutton_)
	line._cid_pressk = fig.canvas.mpl_connect('key_press_event', line._anykey_)
	line._connect_(fig,ax,map)


#+_____________________________________________________________________________+
# analyse selected lines, calculate interpolation points
print(' --> calculate interpolation points')
# which   ... res/npoints
# res     ... interpolation resolution in deg
# npoints ... use npoints per line segment for interpolation
line.analyse_lines(which='res',res=0.1)
#line.analyse_lines(which='res',res=1.0)
if do_anomaly==True:
	line2.line_define = line.line_define
	line2.analyse_lines(which='res',res=0.1)

line.plot_lines_position(mesh)

#+_____________________________________________________________________________+
#|                            *** LOAD 3D DATA ***                             |
#+_____________________________________________________________________________+
# load 3d fesom2.0 data on original depth layers --> line.value
print(' --> load 3d fesom2.0 data')
fesom_load_data_horiz(mesh,line)
if do_anomaly==True:
	fesom_load_data_horiz(mesh,line2)

#+_____________________________________________________________________________+
#|                *** INTERPOLATE 3D DATA ON LINESEGMENT ***                   |
#+_____________________________________________________________________________+
# interpolate points on fesom2.0 data, use routine of Nikolay fesom2regular(...)
print(' --> do interpolation on data')
if do_anomaly==False:
	line.interp_lines(mesh)
	if line.var.find('vec')!=-1: line.calc_flux()
else:
	if line.var.find('vec')!=-1:
		# calculate data on line and line2 data
		line.interp_lines(mesh)
		line2.interp_lines(mesh)
		
		# calculate flux through line and line2
		line.calc_flux()
		line2.calc_flux()
		
		# calculate anomaly line2-line ... 
		anom = cp.copy(line)
		anom.data_anom(line,line2)
	else:
		line.interp_lines(mesh)
		line2.interp_lines(mesh)
		
		anom = cp.copy(line)
		anom.data_anom(line,line2)


#+_____________________________________________________________________________+
#|                        *** PLOT LINE OVER DEPTH ***                         |
#+_____________________________________________________________________________+
print(' --> plot line data')
if  do_anomaly==True and line.var.find('vec')!=-1:
	#___________________________________________________________________________
	# do common crange for line and lin2
	cmax = np.max([np.nanmax(line.value[0]),np.nanmax(line2.value[0])])
	cmin = np.min([np.nanmin(line.value[0]),np.nanmin(line2.value[0])])
	cref = cmin + (cmax-cmin)/2
	cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
	line.crange=[cmin,cmax,cref]	
	line2.crange=[cmin,cmax,cref]
	#___________________________________________________________________________
	line.plot_lines_dist_x_z()
	line2.plot_lines_dist_x_z()
	anom.plot_lines_dist_x_z()
else:
	if do_anomaly==True:
		#_______________________________________________________________________
		# do common crange for line and lin2
		cmax = np.max([np.nanmax(line.value[0]),np.nanmax(line2.value[0])])
		cmin = np.min([np.nanmin(line.value[0]),np.nanmin(line2.value[0])])
		cref = cmin + (cmax-cmin)/2
		cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		line.crange=[cmin,cmax,cref]	
		line2.crange=[cmin,cmax,cref]
		#_______________________________________________________________________
		line.plot_lines_dist_x_z()
		line2.plot_lines_dist_x_z()
		anom.plot_lines_dist_x_z()
	else:
		line.plot_lines_dist_x_z()
