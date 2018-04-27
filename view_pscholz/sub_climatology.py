# Patrick Scholz, 23.01.2018
import numpy as np
import time
import sys
from netCDF4 import Dataset
from sub_regriding_adapt import *
from mpl_toolkits.basemap import shiftgrid
import seawater as sw
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from colormap_c2c import colormap_c2c

#+_____________________________________________________________________________+
#|                                                                             |
#|                        *** CLIMATOLOGY DATA CLASS ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
class clim_data(object):
	#____data run variables______________________
	var                         = ''
	path, fname, descript       = '', '', ''
	
	#____data time variables_____________________
	#year, month, record         = [], [], []
	lon, lat, depth             = [], [], []
	str_time, str_dep           = '', ''
	
	#____data projection variables_______________
	proj, proj_lon, proj_lat    = '', 0.0, 90.0
	
	#____data description info___________________
	sname, lname, unit	        = '', '', ''
	
	#____data plot varaibles_____________________
	cmap, crange, cnumb         = 'grads', [], []
	which_plot                  = ''
	
	#____data varaible___________________________
	value, anom                 = [], []
	
	
	#___INIT CLIM OBJECT_______________________________________________________#
	#
	#__________________________________________________________________________#
	def __init__(self,path,fname,var='temp'):
		self.path = path 
		self.fname= fname
		self.var  = var
		
		#______________________________________________________________________#
		# load WOA 2005 data
		if self.fname=='woa2005TS.nc':
			self.descript = 'WOA2005'
			#__________________________________________________________________#
			ncid  = Dataset(path+fname, 'r') 
			
			#__________________________________________________________________#
			self.lon   = np.copy(ncid.variables['lon'][:])
			self.lat   = np.copy(ncid.variables['lat'][:])
			self.depth = np.copy(ncid.variables['depth'][:])
			valueT= np.copy(ncid.variables['t00an1'][:,:,:,:]).squeeze()
			valueS= np.copy(ncid.variables['s00an1'][:,:,:,:]).squeeze()
			ncid.close()
			
			#__________________________________________________________________#
			# convert missing value to nan
			valueT[valueT>=100]=np.nan
			valueS[valueS>=100]=np.nan
			
			#__________________________________________________________________#
			# shift coordinates from 0...360 --> -180...180
			valueT, dum      = shiftgrid(180.0,valueT, self.lon,start=False,cyclic=359)
			valueS, self.lon = shiftgrid(180.0,valueS, self.lon,start=False,cyclic=359)
			
			#__________________________________________________________________#
			# calculate potential temperature
			if var=='temp':
				self.value = valueT
				
			elif var=='ptemp':
				depth3d = np.zeros(valueS.shape)
				for di in range(0,self.depth.size):
					depth3d[di,:,:] = self.depth[di]
				self.value = sw.ptmp(valueS, valueT, depth3d)
				
			elif var=='salt':
				self.value = valueS
			
			else:
				print(' --> this variable is not supported for this climatology')
	
#___DO VERTICAL INTERPOLATE AVERAGE OVER CERTAIN LAYERS_________________________
#
#_______________________________________________________________________________
def clim_vinterp(data_in,levels):
	data_out = np.zeros((data_in.lat.size,data_in.lon.size))
	aux_div  = np.zeros((data_in.lat.size,data_in.lon.size))
	for di in range(0,len(levels)):
		# find upper and lower layer indices
		#idx_dwn = np.array(np.where( data.depth[di]<=abs(mesh.zlev))).squeeze()
		idx_dwn = np.array(np.where( levels[di]<=abs(data_in.depth))).squeeze()
		idx_dwn = idx_dwn[0]
		idx_up  = idx_dwn-1
		if idx_up<0: idx_up=0
		
		# linear vertical interpolant
		deltaz   = abs(data_in.depth[idx_dwn])-abs(data_in.depth[idx_up])
		#deltaz_i = abs(mesh.zlev[idx_dwn])-data.depth[di]
		deltaz_i = abs(data_in.depth[idx_dwn])-levels[di]
		
		# interpoalte verticaly and sum up
		if deltaz_i==0:
			auxval = data_in.value[idx_dwn,:,:]
		else:
			auxval = data_in.value[idx_dwn,:,:]-(data_in.value[idx_dwn,:,:]-data_in.value[idx_up,:,:])*deltaz_i/deltaz
		aux_div[~np.isnan(auxval)]=aux_div[~np.isnan(auxval)]+1.0
		auxval[np.isnan(auxval)]=0.0
		data_out = data_out + auxval
	
	# calculate mean over depth levels
	idx = aux_div!=0.0
	data_out[idx==True] = data_out[idx==True]/aux_div[idx==True]
	data_out[idx==False] = np.nan
	return(data_out)


#___DO VERTICAL INTERPOLATE AVERAGE OVER CERTAIN LAYERS_________________________
#
#_______________________________________________________________________________
def clim_plot_anom(clim,figsize=[]):
	from set_inputarray import inputarray
	
	if len(figsize)==0 : figsize=[14,10]
	#___________________________________________________________________________
	fig = plt.figure(figsize=figsize)
	#fig.patch.set_alpha(0.0)
	ax  = plt.gca()
	resolution = 'c'
	fsize = 12
	#+_________________________________________________________________________+
	#| SET PROJECTION PARAMETERS                                               |
	#+_________________________________________________________________________+
	if (clim.proj=='ortho' or clim.proj=='stere'):
		map = Basemap(projection = clim.proj, resolution = resolution,
					lat_0      = clim.proj_lat, lon_0 = clim.proj_lon,
					area_thresh=1000)
		ylabels= [0,0,0,0]
		xlabels= [1,1,1,1]
		xticks = np.arange(-180,180,20)
		yticks = np.arange(-90,90,15)
	elif clim.proj=='npstere':
		setbndlat = np.max([0.0,inputarray['which_box'][2]])
		map = Basemap(projection = clim.proj, resolution = resolution,
					lon_0      = 0, boundinglat = setbndlat,
					round      = True)
		ylabels=[0,0,0,0]
		xlabels=[1,1,1,1]
		xticks = np.arange(0.,360.,20.)
		yticks = np.arange(-90.,90.,15.)
	elif clim.proj=='spstere':
		setbndlat = np.min([0.0,inputarray['which_box'][3]])
		map = Basemap(projection = clim.proj, resolution = resolution,
					boundinglat= setbndlat, lon_0      = 180,
					round      = True)
		ylabels=[0,0,0,0]
		xlabels=[1,1,1,1]
		xticks = np.arange(0.,360.,20.)
		yticks = np.arange(-90.,90.,15.)
	else:	
		map = Basemap(projection = clim.proj,resolution = resolution,
						llcrnrlon  = inputarray['which_box'][0], urcrnrlon  = inputarray['which_box'][1],
						llcrnrlat  = inputarray['which_box'][2], urcrnrlat  = inputarray['which_box'][3],
						area_thresh=10000
						)
		ylabels=[1,0,0,0]
		xlabels=[0,0,0,1]
		ticknr   = 8
		tickstep = np.array([1.0,2.0,2.5,5.0,10.0,15.0,20.0,30.0])
		idx      = (inputarray['which_box'][1]-inputarray['which_box'][0])/ticknr
		idx1      = np.array(np.where(tickstep>=idx))
		if idx1.size==0 : 
			xticks   = np.arange(0.,360.,tickstep[-1])
		else:
			xticks   = np.arange(0.,360.,tickstep[idx1[0,0]])
		
		del idx
		idx      = (inputarray['which_box'][3]-inputarray['which_box'][2])/ticknr
		idx1      = np.array(np.where(tickstep>=idx))
		if idx1.size==0 : 
			yticks   = np.arange(-90.,90.,tickstep[-1])
		else:	
			yticks   = np.arange(-90.,90.,tickstep[idx1[0,0]])
		
	#___________________________________________________________________________
	# go from geo coord. to projection coord.
	mlon,mlat = np.meshgrid(clim.lon,clim.lat)
	mx,my     = map(mlon, mlat)

	#___________________________________________________________________________
	# calculate index of which data points are within boi
	#idx_box = ~np.isnan(clim.anom)
	idx_box = mlon<inputarray['which_box'][0]
	idx_box = np.logical_or(idx_box,mlon>inputarray['which_box'][1])
	idx_box = np.logical_or(idx_box,mlat<inputarray['which_box'][2])
	idx_box = np.logical_or(idx_box,mlat>inputarray['which_box'][3])
	idx_box = idx_box==False # true index for triangles that are within box 

	#+_________________________________________________________________________+
	#| set minimum, maximum and reference values for the creation of the       |
	#| adjustable colormap                                                     |
	#+_________________________________________________________________________+
	cnumb = 20; # minimum number of colors
	if len(clim.cnumb)!=0: cnumb=clim.cnumb[0]
	cmax  = np.nanmax(clim.anom[idx_box])
	cmin  = np.nanmin(clim.anom[idx_box])
	cref  = 0
	#clim.cmap='blue2red'

	#+_____________________________________________________________________________+
	# if predefined color range
	if len(clim.crange)!=0:
		if len(clim.crange)==2:
			cmin = np.float(clim.crange[0])
			cmax = np.float(clim.crange[1])
			cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		elif len(clim.crange)==3:
			cmin = np.float(clim.crange[0])
			cmax = np.float(clim.crange[1])
			cref = np.float(clim.crange[2])
		else:
			print(' this colorrange definition is not supported !!!')
			print('data.crange=[cmin,cmax] or data.crange=[cmin,cmax,cref]')

	#+_____________________________________________________________________________+
	# make custom colormap
	print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
	cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,clim.cmap)
	print('clevel = ',clevel)
	do_drawedges=True
	if clevel.size>30: do_drawedges=False
	
	data_plot = np.copy(clim.anom)
	data_plot[data_plot<clevel[ 0]]  = clevel[ 0]+np.finfo(np.float32).eps
	data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
		
	#+_____________________________________________________________________________+
	#| plot data on triangluar grid                                                |
	#+_____________________________________________________________________________+
	#hp1=map.pcolormesh(mlon,mlat,data_plot,
	#hp1=map.pcolormesh(mlon,mlat,data_plot,latlon=True,
	hp1=map.contourf(mlon,mlat,data_plot,clevel,#				shading='gouraud',#flat
					latlon=True,
					antialiased=False,
					edgecolor='None',
					cmap=cmap0,
					vmin=clevel[0], vmax=clevel[-1])
	#hp1.cmap.set_under([1.,1.,1.])

	#______________________________________________________________________________+
	# arange zonal & meriodional gridlines and labels
	map.drawparallels(yticks,#np.arange(-90.,90.,30.),
						labels=ylabels,
						linewidth=0.25,
						dashes=[1,1e-10],
						fontsize=fsize)
	map.drawmeridians(xticks,#np.arange(0.,360.,30.),
						linewidth=0.25,
						labels=xlabels,
						dashes=[1,1e-10],
						fontsize=fsize)
	map.drawmapboundary(fill_color='0.9',linewidth=1.0)
	map.drawcoastlines(color='k',linewidth=0.5)
	map.fillcontinents(color='0.6')

	#___________________________________________________________________________
	# arange colorbar position and labels
	if (clim.proj=='ortho'   or clim.proj=='stere' or \
		clim.proj=='npstere' or clim.proj=='spstere'):
		divider = make_axes_locatable(ax)
		cax     = divider.append_axes("right", size="2.5%", pad=0.70)
	else:
		divider = make_axes_locatable(ax)
		cax     = divider.append_axes("right", size="2.5%", pad=0.1)
	plt.clim(clevel[0],clevel[-1])
	cbar = plt.colorbar(hp1,cax=cax,ticks=clevel,drawedges=do_drawedges,
						extend='neither')

	cbar.set_label(clim.lname+' '+clim.unit+'\n'+clim.str_time+clim.str_dep, size=fsize+2)
	cl = plt.getp(cbar.ax, 'ymajorticklabels')
	plt.setp(cl, fontsize=fsize)

	# kickout some colormap labels if there are to many
	ncbar_l=len(cbar.ax.get_yticklabels()[:])
	idx_cref = np.where(clevel==cref)[0]
	idx_cref = np.asscalar(idx_cref)

	nmax_cbar_l = 10
	nstep = ncbar_l/nmax_cbar_l
	nstep = np.int(np.floor(nstep))
	if nstep==0: nstep=1
   
	plt.setp(cbar.ax.get_yticklabels()[:], visible=False)
	#plt.setp(cbar.ax.get_yticklabels()[::nstep], visible=True)
	plt.setp(cbar.ax.get_yticklabels()[idx_cref::nstep], visible=True)
	plt.setp(cbar.ax.get_yticklabels()[idx_cref::-nstep], visible=True)
	
	#______________________________________________________________________________+
	ax.tick_params(axis='both', direction='out')
	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
	ax.get_yaxis().tick_left()
	plt.sca(ax)
	plt.title(clim.descript+'\n',fontdict= dict(fontsize=24),verticalalignment='bottom')
	
	#+_________________________________________________________________________+
	#| SAVE FIGURE                                                             |
	#+_________________________________________________________________________+
	if inputarray['save_fig']==True:
		print(' --> save figure: png')
		#print(fig.dpi)
		str_times= clim.str_time.replace(' ','').replace(':','')
		str_deps = clim.str_dep.replace(' ','').replace(',','').replace(':','')
		sfname = 'plot_'+clim.descript+'_'+clim.proj+'_'+clim.sname+'_'+str_times+'_'+str_deps+'.png'
		sdname = inputarray['save_figpath'] 
		#if os.path.isdir(sdname)==False: os.mkdir(sdname)
		plt.savefig(sdname+sfname, \
		            format='png', dpi=600, \
					bbox_inches='tight', pad_inches=0,\
					transparent=True,frameon=True)
	
	#___________________________________________________________________________
	plt.show(block=False)
	
	#___________________________________________________________________________
	return(fig,ax,map,cbar)
	#return map