import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import cmocean
import matplotlib.gridspec as gridspec

#+____IMPORT FESOM RELATET ROUTINES____________________________________________+
from set_inputarray 		import *
from sub_fesom_mesh 		import * 
from sub_fesom_data 		import * 
from sub_fesom_plot 		import *
from sub_fesom_selectline 	import *
from sub_regriding_adapt 	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
#+_____________________________________________________________________________+
inputarray  = set_inputarray()
data 		= fesom_init_data(inputarray) # init fesom2.0 data object
data.var 	= 'temp'
#+_____________________________________________________________________________+
# select year to average over [start_yr, end_yr]
data.year	= [1948,1948]

# select month to average over
data.month	= [1,2,12]

# select linear interpolated depth layers to average over, empty mean use all layers
data.depth	= []

#+_____________________________________________________________________________+
#+_____________________________________________________________________________+
# predefine section lines or leave it empty in this case interactive selection 
line 		= fesom_init_line()

#line.line_refxy    =  list([[],[]]) 
##line.line_refxy    =  list([[]]) 
#line.line_refxy[0] = [ [-75.188, -72.875, -65.875, -57.75, -50.156],
					   #[39.531, 36.594, 34.094, 34.312, 32.094],
					   #'test1' ]
#line.line_refxy[1] = [ [-25.8,-29.3,-32.2,-35.9,-38.2,-40.3,-41.1,-42.5,-45.9,-50.6,-53.4,-55.3,-57.6,-60.2],
					   #[67.3,65.4,65.1,63.7,62.9,61.9,59.6,58.6,58.4,60.8,62.5,62.7,61.91,60.8],
					   #'test2' ] 

#+_____________________________________________________________________________+
#|                         *** LOAD FVSOM MESH ***                             |
#+_____________________________________________________________________________+
inputarray['which_plot']  = 'pcolor'
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
#fesom_plot_lmask(map,mesh,ax,'0.6')
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
if len(line.line_refxy)==0:
	# interactively
	line.cid_pressb = fig.canvas.mpl_connect('button_press_event', line.anybutton)
	line.cid_pressk = fig.canvas.mpl_connect('key_press_event', line.anykey)
	line.connect(fig,ax,map)
else:
	# predefined
	for ii in range(0,len(line.line_refxy)):
		ax.plot(line.line_refxy[ii][0],line.line_refxy[ii][1] ,color='w',linewidth=2,marker='o',mfc='w',mec='k',axes=ax)
		ax.plot(line.line_refxy[ii][0][0],line.line_refxy[ii][1][0] ,color='k',linewidth=0.5,marker='o',mfc=[0.0,0.8,0.0],mec='k',axes=ax)
		ax.plot(line.line_refxy[ii][0][-1],line.line_refxy[ii][1][-1] ,color='k',linewidth=0.5,marker='o',mfc=[0.8,0.0,0.0],mec='k',axes=ax)
	plt.show(block=False)
	fig.canvas.draw()

#+_____________________________________________________________________________+
# analyse selected lines, calculate interpolation points
print(' --> calculate interpolation points')
line.analyse_lines(npoints=100)
#print(input.line_ipxy)


#+_____________________________________________________________________________+
#|                  *** LOAD AND INTERPOLATE 3D DATA ***                       |
#+_____________________________________________________________________________+
# load 3d fesom2.0 data on original depth layers --> data.value
print(' --> load 3d fesom2.0 data')
data=fesom_load_data_horiz(mesh,data)

#+_____________________________________________________________________________+
# interpolate points on fesom2.0 data, use routine of Nikolay fesom2regular(...)
print(' --> do interpolation on data')
line.interp_lines(mesh,data)

#+_____________________________________________________________________________+
#|                  *** PLOT LINE OVER DEPTH ***                               |
#+_____________________________________________________________________________+
for ii in range(0,len(line.line_refxy)):
	
	fsize=14
	fig = plt.figure(figsize=(20, 13))
	ax1  = plt.gca()
	# duplicate x-axes
	ax2  = ax1.twiny()
	
	#--------------------------------------------------------------------------#
	cnumb= 25
	cmin = np.nanmin(line.line_data[ii])
	cmax = np.nanmax(line.line_data[ii])
	cref = cmin + (cmax-cmin)/2
	#cref =0.0
	cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
	cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,'grads')
	do_drawedges=True
	if clevel.size>30: do_drawedges=False
	
	# overwrite discrete colormap
	#cmap0 = cmocean.cm.balance
	#do_drawedges=False
	
	#--------------------------------------------------------------------------#
	# make pcolor or contour plot 
	xxx = np.array(line.line_ipdist[ii])
	xxx1= np.array(line.line_ipxy[ii][0])
	yyy = mesh.zlev[0:-1] + (mesh.zlev[1::]-mesh.zlev[0:-1])/2.0
	yyy = -yyy
	yyylim = np.sum(~np.isnan(line.line_data[ii]),axis=1).max()+1
	if yyylim==yyy.shape: yyylim=yyylim-1
	yy,xx = np.meshgrid(yyy,xxx)
	
	if inputarray['which_plot']=='pcolor':
		hp=plt.pcolormesh(xx,yy,line.line_data[ii],
					shading='flat',#flat
					antialiased=False,
					edgecolor='None',
					cmap=cmap0,
					vmin=np.nanmin(line.line_data[ii]), vmax=np.nanmax(line.line_data[ii]))
	else: 
		hp=plt.contourf(xx,yy,line.line_data[ii],levels=clevel,
					antialiased=False,
					edgecolor='None',
					cmap=cmap0,
					vmin=np.nanmin(line.line_data[ii]), vmax=np.nanmax(line.line_data[ii]))
	hp.cmap.set_under([0.4,0.4,0.4])
	
	#--------------------------------------------------------------------------#
	# plot position of section points
	plt.plot(line.line_refdist[ii][0] ,0.0,        color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
	plt.plot(line.line_refdist[ii][0] ,yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
	plt.plot(line.line_refdist[ii][-1],0.0,        color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
	plt.plot(line.line_refdist[ii][-1],yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
	if len(line.line_refdist[ii])>2:
		for jj in range(1,len(line.line_refdist[ii])-1):
			plt.plot(np.ones((2,))*line.line_refdist[ii][jj],np.array([0.0,yyy[yyylim]]),color='k',linewidth=0.5,linestyle='--',marker=None)
			plt.plot(line.line_refdist[ii][jj],0.0,        color='k',linewidth=0.5,marker="^",markersize=8,mfc='k', clip_on=False)
			plt.plot(line.line_refdist[ii][jj],yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=8,mfc='k', clip_on=False)
	
	#--------------------------------------------------------------------------#
	# set main axes
	ax1.set_xlim(xxx.min(),xxx.max())
	ax1.set_ylim(0,yyy[yyylim])
	ax1.invert_yaxis()
	ax1.set_axis_bgcolor([0.25,0.25,0.25])
	ax1.tick_params(axis='both',which='major',direction='out',length=8)
	ax1.minorticks_on()
	ax1.tick_params(axis='both',which='minor',direction='out',length=4)
	ax1.set_xlabel('Distance from start point [km]',fontdict=dict(fontsize=12))
	ax1.set_ylabel('Depth [km]',fontdict=dict(fontsize=12))
	
	#--------------------------------------------------------------------------#
	# set upper secondary x-axes
	ax2.set_xlim(xxx.min(),xxx.max())
	ax2.set_ylim(0,yyy[yyylim])
	ax2.invert_yaxis()
	
	#lonticks = np.arange(-180.0,180.0,5.0)
	#lonticks = lonticks[lonticks>=xxx1.min()]
	#lonticks = lonticks[lonticks<=xxx1.max()]
	#lonticklabels=[]
	#for jj in lonticks: lonticklabels.append(str(jj))
	#xticks   = np.interp(lonticks,line.line_ipxy[ii][0],line.line_ipdist[ii])
	#ax2.set_xticks(xticks)
	#ax2.set_xticklabels(lonticklabels)
	#ax2.set_xlabel('Longitude [deg]',fontdict=dict(fontsize=12))
	#ax2.tick_params(axis='both',which='major',direction='out',length=8)
	
	lonticklabels=[]
	for jj in range(0,len(line.line_refxy[ii][0])): 
		strlon='E' 
		if line.line_refxy[ii][0][jj]<0:strlon='W'
		strlat='N' ; 
		if line.line_refxy[ii][1][jj]<0:strlon='S'
		lonticklabels.append('{:2.2f}$^{{\\circ}}${:s}\n{:2.2f}$^{{\\circ}}${:s}'.format(np.abs(line.line_refxy[ii][0][jj]),strlon,np.abs(line.line_refxy[ii][1][jj]),strlat))
	ax2.set_xticks(line.line_refdist[ii])
	ax2.set_xticklabels(lonticklabels,fontdict=dict(fontsize=10))
	ax2.set_xlabel('Longitude/Latitude [deg]',fontdict=dict(fontsize=12))
	ax2.tick_params(axis='both',which='major',direction='out',length=8)
	
	#--------------------------------------------------------------------------#
	# draw colorbar
	divider = make_axes_locatable(ax1)
	cax     = divider.append_axes("right", size="2.5%", pad=0.5)
	plt.clim(clevel[0],clevel[-1])
	
	cbar = plt.colorbar(hp,ax=[ax1,ax2],cax=cax,ticks=clevel,drawedges=do_drawedges)
	cbar.set_label(data.lname+' '+data.unit+'\n'+data.str_time, size=fsize+2)
	
	cl = plt.getp(cbar.ax, 'ymajorticklabels')
	plt.setp(cl, fontsize=fsize)
	
	# kickout some colormap labels if there are to many
	ncbar_l=len(cbar.ax.get_yticklabels()[:])
	idx_cref = np.where(clevel==cref)[0]
	idx_cref = np.asscalar(idx_cref)
	nmax_cbar_l = 10
	nstep = ncbar_l/nmax_cbar_l
	plt.setp(cbar.ax.get_yticklabels()[:], visible=False)
	plt.setp(cbar.ax.get_yticklabels()[idx_cref::nstep], visible=True)
	plt.setp(cbar.ax.get_yticklabels()[idx_cref::-nstep], visible=True)
	
	#--------------------------------------------------------------------------#
	# bug fix workaround to proper draw secondary axes
	fig.canvas.draw()
	ax2.set_position(ax1.get_position())
	
	#--------------------------------------------------------------------------#
	plt.show(block=False)
	print('finish')
	