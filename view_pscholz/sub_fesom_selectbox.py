# Patrick Scholz, 23.01.2018
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sub_fesom_mesh 		import * 
from sub_fesom_plot 		import *
from colormap_c2c			import *


#+_____________________________________________________________________________+
class fesom_box:
	
	which_obj					= 'box'
	
	#____data run variables______________________
	var                         = ''
	runid, path, descript       = 'fesom', '', ''
	
	#____data time variables_____________________
	year,  month, record, depth = [], [], [], []
	str_time, str_dep           = '', ''
	
	#____data projection variables_______________
	proj, proj_lon, proj_lat    = '', 0.0, 90.0
	
	#____data description info___________________
	sname, lname, unit	        = '', '', ''
	
	#____data plot variables_____________________
	cmap, crange, cnumb         = 'grads', [], []
	which_plot                  = ''
	
	#____data variable___________________________
	value, value2, value3, time = [], [], [], [] 
	valueflx					= []; # store calculated fluxes (volume, heat, salt)
	which_mean                  = 'monthly'
	anom                        = False
	
	#____box variable____________________________
	box_define                  = []
	box_idx                     = []
	
	#___INIT BOX OBJECT________________________________________________________#
	#
	#__________________________________________________________________________#
	def __init__(self):
		#____selection variable______________________
		self._x, self._y, self._xc, self._yc  = [], [], [], []
		self._b, self._k, 					  = [] ,[]
		self._cid_pressb, self._cid_pressk	  = [], [],
		self._zoomfac, self._press			  = 20.0, 'None'
		self._fig, self._ax, self._map		  = [], [], []
		self._xlim, self._ylim				  = [], [],
		self._ptsx, self._ptsy, self._drawline= [], [], []
		self._text							  = []
	
	#+_________________________________________________________________________+
	def _connect_(self,fig,ax,map):
		print(' --> start interactive line selection')
		#self._figure = plt.gcf()
		#self._ax     = plt.gca()
		self._figure = fig
		self._ax     = ax
		self._map    = map
		#self.canvas = FigureCanvas(self._figure)
		self._xlim   = self._ax.get_xlim()
		self._ylim   = self._ax.get_ylim()
		self._text   = plt.figtext(0.5,0.05,'Interactive Selection',
							fontsize=24,
							fontweight='bold',
							horizontalalignment='center',
							bbox=dict(facecolor='w',edgecolor='k'),
							)
		# This line stops proceeding of the code...
		plt.show(block=True)
		
	#+_________________________________________________________________________+
	# change zoom factor to zoom in 
	def _drawzoom_in_(self):	
		print(' zoom_in')
		self._zoomfac = self._zoomfac/2.0
		
	#+_________________________________________________________________________+
	# change zoom factor to zoom out 
	def _drawzoom_out_(self):	
		print(' zoom_out')
		self._zoomfac = self._zoomfac*2.0
		
	#+_________________________________________________________________________+
	# zoom full out
	def _drawzoom_outfull_(self):	
		print(' zoom_fullout')
		self._ax.set_xlim(self._xlim[0],self._xlim[1])
		self._ax.set_ylim(self._ylim[0],self._ylim[1])
		#plt.xscale('linear')
		#plt.yscale('linear')
		self._figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# aplly actual zoom factor
	def _zoom_(self, k):
		self._k = k
		if self._k == '+':
			self._drawzoom_in_()
		elif self._k == '-':
			self._drawzoom_out_()
		self._ax.set_xlim(np.max([self._xc-self._zoomfac,-180.0]), np.min([self._xc+self._zoomfac,180.0]))
		self._ax.set_ylim(np.max([self._yc-self._zoomfac, -90.0]), np.min([self._yc+self._zoomfac, 90.0]))
		#plt.xscale('linear')
		#plt.yscale('linear')
		self._figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# update center position with left mouse click 
	def _update_center_(self, xc, yc):
		self._xc, self._yc, = xc, yc
		
	#+_________________________________________________________________________+
	# center new leftclick point
	def _move_center_(self):
		self._ax.set_xlim(np.max([self._xc-self._zoomfac,-180.0]), np.min([self._xc+self._zoomfac,180.0]))
		self._ax.set_ylim(np.max([self._yc-self._zoomfac, -90.0]), np.min([self._yc+self._zoomfac, 90.0]))
		#plt.xscale('linear')
		#plt.yscale('linear')
		self._figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# BoxBuilder
	def _Boxbuilder_(self,x,y,mode='add'):
		if np.size(x)==0 or np.size(y)==0: 
			self._ptsx = list(self._drawbox.get_xdata())
			self._ptsy = list(self._drawbox.get_ydata())
			plt.xscale('linear')
			plt.yscale('linear')
			return
		if mode == 'add':
			self._ptsx.append(np.float16(x))
			self._ptsy.append(np.float16(y))
			if len(self._ptsx)>=4:
				self._ptsx = [np.min(self._ptsx),np.max(self._ptsx),np.max(self._ptsx),np.min(self._ptsx),np.min(self._ptsx)]
				self._ptsy = [np.min(self._ptsy),np.min(self._ptsy),np.max(self._ptsy),np.max(self._ptsy),np.min(self._ptsy)]
				self._drawbox.set_linestyle('-')
		elif mode == 'remove':
			self._ptsx=[]
			self._ptsy=[]
			self._drawbox.set_linestyle('None')
			
		mboxx,mboxy = self._map(self._ptsx,self._ptsy)
		self._drawbox.set_data(mboxx, mboxy)
		self._drawbox.set_axes(self._ax)
		self._drawbox.figure.canvas.draw()
		plt.xscale('linear')
		plt.yscale('linear')
		self._figure.canvas.draw()
	
	#+_________________________________________________________________________+
	# what should happen if a mouse button is clicked when cursor is over axis 
	def _anybutton_(self,event):
		#global self
		if event.inaxes:
			x, y, b, k = event.xdata, event.ydata, event.button, []
			#print('button=%d, xdata=%f, ydata=%f' % ( event.button, event.xdata, event.ydata))
			#+___left mouse button_____________________________________________+
			if event.button==1 : 
				xc, yc = x, y
				if self._press=='None':
					self._update_center_(xc, yc)
					self._move_center_()
				elif self._press=='Box':
					self._Boxbuilder_(x,y)
			#+___middle mouse button___________________________________________+
			if event.button==2 : return
			if (event.xdata is None): return
			
			#+___right mouse button____________________________________________+
			if event.button==3 : 
				self._drawzoom_outfull_()
		
	#+_________________________________________________________________________+
	# what should happen if a keboard button is pressed when cursor is over axis 
	def _anykey_(self,event):
		#global self
		if event.inaxes:
			x, y , b, k = event.xdata, event.ydata, [], event.key
			
			#___ZOOM IN [+]_____________________________________________________
			if k=='+':
				self._zoom_(k)
				
			#___ZOOM OUT [-]____________________________________________________
			elif k=='-':
				self._zoom_(k)
				
			#___FINISH INTERACTION [q]__________________________________________
			elif k=='q':
				self._drawzoom_outfull_()
				
				# dis_connect_s button and key events
				self._figure.canvas.mpl_dis_connect_(self._cid_pressb)
				self._figure.canvas.mpl_dis_connect_(self.cid_pressk)
				
				# This line allows proceeding of the code...
				plt.close(self._figure)
				
			#___BOX MODE [l]___________________________________________________
			elif k=='b':
				#---------------------------------------------------------------
				# activate line mode, start line 
				if self._press == 'None':
					print('draw box: [ON]')
					self._text.set_text('draw box: [ON]')
					plt.draw(),self._figure.canvas.draw()
					
					self._press = 'Box'
					self._drawbox, = plt.plot([], [],color='w',linestyle='None',linewidth=2,marker='o',axes=self._ax)
					self._Boxbuilder_([],[])
					#plt.xscale('linear'),
					#plt.yscale('linear')
					self._figure.canvas.draw()
				#---------------------------------------------------------------
				# switch off line mode, end line
				elif self._press =='Box': 
					print('draw line: [OFF]')
					self._text.set_text('draw line: [OFF]')
					plt.draw(),self._figure.canvas.draw()
					self._press = 'None'
					self._drawbox.set_color([0.0,0.8,0.0])
					#self.line.set_color('w')
					self._drawbox.set_linewidth(2.0)
					plt.xscale('linear')
					plt.yscale('linear')
					self._figure.canvas.draw()
					#ask for name 
					#self._text.set_text('enter name in command line')
					#plt.draw(),self._figure.canvas.draw()
					#name = raw_input("Enter Name of line: ")
					name = 'default'
					self.box_define.append([])
					self.box_define[len(self.box_define)-1]=[self._ptsx[0:-1],self._ptsy[0:-1],name]
					
			#___DELETE DRAW POINTS [d]__________________________________________
			elif k=='d':
				if self._press=='Box':
					self._Boxbuilder_(x,y,mode='remove')
	
	
	#+_________________________________________________________________________+
	# interpolate line points on fesom data
	def select_pointsinbox(self,mesh):
		self.box_idx = []
		for ii in range(0,len(self.box_define)):
			self.box_idx.append([])
			self.box_idx[ii-1] = mesh.nodes_2d_xg >= np.min(self.box_define[ii][0])
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_xg <= np.max(self.box_define[ii][0]))
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_yg >= np.min(self.box_define[ii][1]))
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_yg <= np.max(self.box_define[ii][1]))
		self.box_idx[ii-1] = np.array(self.box_idx[ii-1])
		
	#+_________________________________________________________________________+
	#|																		   |
	#+_________________________________________________________________________+
	# calculate fluxes through section
	def data_anom(self, box, box2):
		self.descript     = box2.descript+'-'+box.descript
		if box.str_time!=box2.str_time:
			self.str_time = box2.str_time+'-'+box.str_time
		else:
			self.str_time = box.str_time
		self.str_dep = box.str_dep
		#_______________________________________________________________________
		#self.value   = []
		#self.value2  = []
		#self.value3  = []
		#self.valueflx= []
		#if len(box2.value)!=0:
			#for ii in range(0,len(box2.value)):
				#self.value.append([])
				#self.value[ii] = box2.value[ii]-box.value[ii]
				#if len(box2.value2)!=0:
					#self.value2.append([])
					#self.value2[ii]  = box2.value2[ii]-box.value2[ii]
				#if len(box2.value3)!=0:
					#self.value3.append([])
					#self.value3[ii]  = box2.value3[ii]-box.value3[ii]
				#if len(box2.valueflx)!=0:
					#self.valueflx.append([])
					#self.valueflx[ii]= box2.valueflx[ii]-box.valueflx[ii]
		self.value = box2.value-box.value
		self.anom  = True
		self.crange= []
		self.cmap  = 'blue2red'
		self.depth = box.depth
		self.sname = box.sname 
		self.lname = box.lname
		self.unit  = box.unit
	
	
	#___PLOT FESOM2.0 DATA IN INDEX BOX OVER DEPTH AND TIME DATA________________
	# 
	#___________________________________________________________________________
	def plot_index_t_x_z(self):
		from set_inputarray import inputarray
		
		fsize=14
		fig = plt.figure(figsize=(20, 13))
		ax1  = plt.gca()
		# duplicate x-axes
		
		#+______________________________________________________________________
		cnumb= 25
		cmin = np.nanmin(self.value)
		cmax = np.nanmax(self.value)
		cref = cmin + (cmax-cmin)/2
		cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		#cref =0.0
		
		#+______________________________________________________________________
		# if anomaly data	
		if self.anom==True: 
			cmax,cmin,cref = np.nanmax(self.value), np.nanmin(self.value), 0.0
			self.cmap='blue2red'
		#+______________________________________________________________________
		# if predefined color range	
		if len(self.crange)!=0:
			if len(self.crange)==2:
				cmin = np.float(self.crange[0])
				cmax = np.float(self.crange[1])
				cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
			elif len(self.crange)==3:
				cmin = np.float(self.crange[0])
				cmax = np.float(self.crange[1])
				cref = np.float(self.crange[2])
			else:
				print(' this colorrange definition is not supported !!!')
				print('data.crange=[cmin,cmax] or data.crange=[cmin,cmax,cref]')
		
		print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
		cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,self.cmap)
		print('clevel = ',clevel)
		
		do_drawedges=True
		if clevel.size>30: do_drawedges=False
		
		# overwrite discrete colormap
		#cmap0 = cmocean.cm.balance
		#do_drawedges=False
		
		#+______________________________________________________________________
		# make pcolor or contour plot 
		depth = self.depth[:-1] + (self.depth[1:]-self.depth[:-1])/2.0
		depth = -depth
		depthlim = np.sum(~np.isnan(self.value[0,:])).max()
		if depthlim==depth.shape: depthlim=depthlim-1
		yy,xx = np.meshgrid(depth,self.time)
		
		#+______________________________________________________________________
		data_plot = np.copy(self.value)
		data_plot[data_plot<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
		data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
		
		#+______________________________________________________________________
		if inputarray['which_plot']=='pcolor':
			hp=plt.pcolormesh(xx[:,0:depthlim],yy[:,0:depthlim],data_plot[:,0:depthlim],
						shading='flat',#flat
						antialiased=False,
						edgecolor='None',
						cmap=cmap0),
						#vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
		else: 
			hp=plt.contourf(xx[:,0:depthlim],yy[:,0:depthlim],data_plot[:,0:depthlim],levels=clevel,
						antialiased=False,
						cmap=cmap0,
						vmin=clevel[0], vmax=clevel[-1])
			#plt.contour(xx,yy,data.value,levels=clevel,
						#antialiased=True,
						#colors='k',
						#linewidths=0.5,
						#linestyles='solid',
						#vmin=np.nanmin(data.value), vmax=np.nanmax(data.value))
		#hp.cmap.set_under([0.4,0.4,0.4])
		
		#+______________________________________________________________________
		# set main axes
		ax1.set_xlim(self.time.min(),self.time.max())
		ax1.set_ylim(0,depth[depthlim-1])
		ax1.invert_yaxis()
		ax1.set_axis_bgcolor([0.25,0.25,0.25])
		ax1.tick_params(axis='both',which='major',direction='out',length=8)
		ax1.minorticks_on()
		ax1.tick_params(axis='both',which='minor',direction='out',length=4)
		ax1.set_xlabel('Time [years]',fontdict=dict(fontsize=12))
		ax1.set_ylabel('Depth [km]',fontdict=dict(fontsize=12))
		plt.title(self.descript+' - '+self.box_define[0][2]+'\n',fontdict= dict(fontsize=24),verticalalignment='bottom')
		
		#+______________________________________________________________________
		# draw colorbar
		divider = make_axes_locatable(ax1)
		cax     = divider.append_axes("right", size="2.5%", pad=0.5)
		plt.clim(clevel[0],clevel[-1])
		
		cbar = plt.colorbar(hp,ax=ax1,cax=cax,ticks=clevel,drawedges=do_drawedges)
		cbar.set_label(self.lname+' '+self.unit+'\n'+self.str_time, size=fsize+2)
		
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
		
		#___________________________________________________________________
		# save figure
		if inputarray['save_fig']==True:
			print(' --> save figure: png')
			str_times= self.str_time.replace(' ','').replace(':','')
			str_deps = self.str_dep.replace(' ','').replace(',','').replace(':','')
			sfname = 'boxplot_'+self.box_define[0][2]+'_'+self.descript+'_'+self.sname+'_'+str_times+'_'+str_deps+'.png'
			sdname = inputarray['save_figpath']
			if os.path.isdir(sdname)==False: os.mkdir(sdname)
			plt.savefig(sdname+sfname, \
						format='png', dpi=600, \
						bbox_inches='tight', pad_inches=0,\
						transparent=True,frameon=True)
		#+______________________________________________________________________
		plt.show(block=False)
		
		return(fig,ax1,cax)
		print('finish')
		
		
	#___PLOT FESOM2.0 DATA IN INDEX BOX POSITION________________________________
	# 
	#___________________________________________________________________________
	def plot_index_position(self,mesh):
		from set_inputarray import inputarray
		# draw position of box 
		for ii in range(0,len(self.box_define)):
			#___________________________________________________________________
			xmin,xmax = np.min(self.box_define[ii][0]), np.max(self.box_define[ii][0])
			ymin,ymax = np.min(self.box_define[ii][1]), np.max(self.box_define[ii][1])
			xmin,xmax,ymin,ymax = xmin-20.0, xmax+20.0, ymin-20.0, ymax+20.0
			xmin,xmax,ymin,ymax = np.max([xmin,-180.0]),np.min([xmax,180.0]),np.max([ymin,-90.0]),np.min([ymax,90.0])
			
			#___________________________________________________________________
			fig, ax = plt.figure(figsize=(13, 13)), plt.gca()
			map 	= Basemap(projection = 'cyl',resolution = 'c',
						llcrnrlon = xmin, urcrnrlon = xmax, llcrnrlat = ymin, urcrnrlat = ymax)
			mx,my 	= map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
			map.drawmapboundary(fill_color='0.9',linewidth=1.0)
			
			xlabels,ylabels=[0,0,0,1],[1,0,0,0]
			xticks,yticks = np.arange(0.,360.,10.), np.arange(-90.,90.,5.)
			map.drawparallels(yticks,labels=ylabels,fontsize=14)
			map.drawmeridians(xticks,labels=xlabels,fontsize=14)
			map.bluemarble()
			fesom_plot_lmask(map,mesh,ax,'none','r')
			ax.grid(color='k', linestyle='-', linewidth=0.5)
			
			#___________________________________________________________________
			patch=[]
			patch.append(Polygon(zip(self.box_define[ii][0],self.box_define[ii][1]),closed=True,clip_on=True) )
			ax.plot(self.box_define[ii][0]    ,self.box_define[ii][1] ,linestyle='None'   ,color='w',linewidth=2.0,marker='o',mfc='w',mec='k',axes=ax)
			ax.add_collection(PatchCollection(patch, alpha=1.0,facecolor='none',edgecolor='w',zorder=1,linewidth=2.0,hatch='/'))	
			
			#___________________________________________________________________
			# save figure
			if inputarray['save_fig']==True:
				print(' --> save figure: png')
				sfname = 'boxplot_'+self.box_define[ii][2]+'_position'+'.png'
				sdname = inputarray['save_figpath']
				if os.path.isdir(sdname)==False: os.mkdir(sdname)
				plt.savefig(sdname+sfname, \
							format='png', dpi=600, \
							bbox_inches='tight', pad_inches=0,\
							transparent=True,frameon=True)
			
			#___________________________________________________________________
			plt.show(block=False)
			fig.canvas.draw()
