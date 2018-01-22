import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sub_fesom_mesh 		import * 
from colormap_c2c			import *

def fesom_init_box():

	global box 
	box 			= fesom_box()
	box.box_refxy = []
	#line.line_refxy.append([])
	#line.line_ipxy  = s
	#line.line_ipdist= []
	return box
	
#+_____________________________________________________________________________+
class fesom_box:
	#+_________________________________________________________________________+
	def __init__(self):
		self.press 		= 'None'
		self.x 			= []
		self.y			= []
		self.xc 		= []
		self.yc			= []
		self.zoomfac    = 20.0
		self.b			= []
		self.k			= []
		self.fig        = []
		self.ax         = []
		self.cid_pressb = []
		self.cid_pressk = []
		self.xlim       = []
		self.ylim       = []
		self.box       = []
		self.boxx      = []
		self.boxy      = []
		self.box_refxy = []
		self.box_idx   = []
		self.box_data  = []
	#+_________________________________________________________________________+
	def connect(self,fig,ax,map):
		print(' stept in connect')
		self.figure = plt.gcf()
		self.ax     = plt.gca()
		self.map    = map
		#self.canvas = FigureCanvas(self.figure)
		self.xlim   = self.ax.get_xlim()
		self.ylim   = self.ax.get_ylim()
		self.text   = plt.figtext(0.5,0.05,'Interactive Selection',
							fontsize=24,
							fontweight='bold',
							horizontalalignment='center',
							bbox=dict(facecolor='w',edgecolor='k'),
							)
		# This line stops proceeding of the code...
		plt.show(block=True)
		
	#+_________________________________________________________________________+
	# change zoom factor to zoom in 
	def drawzoom_in(self):	
		print(' zoom_in')
		self.zoomfac = self.zoomfac/2.0
		
	#+_________________________________________________________________________+
	# change zoom factor to zoom out 
	def drawzoom_out(self):	
		print(' zoom_out')
		self.zoomfac = self.zoomfac*2.0
		
	#+_________________________________________________________________________+
	# zoom full out
	def drawzoom_outfull(self):	
		print(' zoom_fullout')
		self.ax.set_xlim(self.xlim[0],self.xlim[1])
		self.ax.set_ylim(self.ylim[0],self.ylim[1])
		#plt.xscale('linear')
		#plt.yscale('linear')
		self.figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# aplly actual zoom factor
	def zoom(self, k):
		self.k = k
		if self.k == '+':
			self.drawzoom_in()
		elif self.k == '-':
			self.drawzoom_out()
		self.ax.set_xlim(np.max([self.xc-self.zoomfac,-180.0]), np.min([self.xc+self.zoomfac,180.0]))
		self.ax.set_ylim(np.max([self.yc-self.zoomfac, -90.0]), np.min([self.yc+self.zoomfac, 90.0]))
		#plt.xscale('linear')
		#plt.yscale('linear')
		self.figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# update center position with left mouse click 
	def update_center(self, xc, yc):
		self.xc, self.yc, = xc, yc
		
	#+_________________________________________________________________________+
	# center new leftclick point
	def move_center(self):
		self.ax.set_xlim(np.max([self.xc-self.zoomfac,-180.0]), np.min([self.xc+self.zoomfac,180.0]))
		self.ax.set_ylim(np.max([self.yc-self.zoomfac, -90.0]), np.min([self.yc+self.zoomfac, 90.0]))
		#plt.xscale('linear')
		#plt.yscale('linear')
		self.figure.canvas.draw()
		
	#+_________________________________________________________________________+
	# BoxBuilder
	def Boxbuilder(self,x,y,mode='add'):
		if np.size(x)==0 or np.size(y)==0: 
			self.boxx = list(self.box.get_xdata())
			self.boxy = list(self.box.get_ydata())
			plt.xscale('linear')
			plt.yscale('linear')
			return
		if mode == 'add':
			self.boxx.append(np.float16(x))
			self.boxy.append(np.float16(y))
			if len(self.boxx)>=4:
				self.boxx = [np.min(self.boxx),np.max(self.boxx),np.max(self.boxx),np.min(self.boxx),np.min(self.boxx)]
				self.boxy = [np.min(self.boxy),np.min(self.boxy),np.max(self.boxy),np.max(self.boxy),np.min(self.boxy)]
				self.box.set_linestyle('-')
		elif mode == 'remove':
			self.boxx=[]
			self.boxy=[]
			self.box.set_linestyle('None')
			
		mboxx,mboxy = self.map(self.boxx,self.boxy)
		self.box.set_data(mboxx, mboxy)
		self.box.set_axes(self.ax)
		self.box.figure.canvas.draw()
		plt.xscale('linear')
		plt.yscale('linear')
		self.figure.canvas.draw()
	
	#+_________________________________________________________________________+
	# what should happen if a mouse button is clicked when cursor is over axis 
	def anybutton(self,event):
		#global self
		if event.inaxes:
			x, y, b, k = event.xdata, event.ydata, event.button, []
			#print('button=%d, xdata=%f, ydata=%f' % ( event.button, event.xdata, event.ydata))
			#+___left mouse button_____________________________________________+
			if event.button==1 : 
				xc, yc = x, y
				if self.press=='None':
					self.update_center(xc, yc)
					self.move_center()
				elif self.press=='Box':
					self.Boxbuilder(x,y)
			#+___middle mouse button___________________________________________+
			if event.button==2 : return
			if (event.xdata is None): return
			
			#+___right mouse button____________________________________________+
			if event.button==3 : 
				self.drawzoom_outfull()
		
	#+_________________________________________________________________________+
	# what should happen if a keboard button is pressed when cursor is over axis 
	def anykey(self,event):
		#global self
		if event.inaxes:
			x, y , b, k = event.xdata, event.ydata, [], event.key
			
			#___ZOOM IN [+]_____________________________________________________
			if k=='+':
				self.zoom(k)
				
			#___ZOOM OUT [-]____________________________________________________
			elif k=='-':
				self.zoom(k)
				
			#___FINISH INTERACTION [q]__________________________________________
			elif k=='q':
				self.drawzoom_outfull()
				
				# disconnects button and key events
				self.figure.canvas.mpl_disconnect(self.cid_pressb)
				self.figure.canvas.mpl_disconnect(self.cid_pressk)
				
				# This line allows proceeding of the code...
				plt.close(self.figure)
				
			#___BOX MODE [l]___________________________________________________
			elif k=='b':
				#---------------------------------------------------------------
				# activate line mode, start line 
				if self.press == 'None':
					print('draw box: [ON]')
					self.text.set_text('draw box: [ON]')
					plt.draw(),self.figure.canvas.draw()
					
					self.press = 'Box'
					self.box, = plt.plot([], [],color='w',linestyle='None',linewidth=2,marker='o',axes=self.ax)
					self.Boxbuilder([],[])
					#plt.xscale('linear'),
					#plt.yscale('linear')
					self.figure.canvas.draw()
				#---------------------------------------------------------------
				# switch off line mode, end line
				elif self.press =='Box': 
					print('draw line: [OFF]')
					self.text.set_text('draw line: [OFF]')
					plt.draw(),self.figure.canvas.draw()
					self.press = 'None'
					self.box.set_color([0.0,0.8,0.0])
					#self.line.set_color('w')
					self.box.set_linewidth(2.0)
					plt.xscale('linear')
					plt.yscale('linear')
					self.figure.canvas.draw()
					#ask for name 
					#self.text.set_text('enter name in command line')
					#plt.draw(),self.figure.canvas.draw()
					#name = raw_input("Enter Name of line: ")
					name = 'default'
					self.box_refxy.append([])
					self.box_refxy[len(self.box_refxy)-1]=[self.boxx[0:-1],self.boxy[0:-1],name]
					
			#___DELETE DRAW POINTS [d]__________________________________________
			elif k=='d':
				if self.press=='Box':
					self.Boxbuilder(x,y,mode='remove')
	
	#+_________________________________________________________________________+
	# interpolate line points on fesom data
	def select_pointsinbox(self,mesh):
		self.box_idx = []
		for ii in range(0,len(self.box_refxy)):
			self.box_idx.append([])
			self.box_idx[ii-1] = mesh.nodes_2d_xg >= np.min(self.box_refxy[ii][0])
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_xg <= np.max(box.box_refxy[ii][0]))
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_yg >= np.min(box.box_refxy[ii][1]))
			self.box_idx[ii-1] = np.logical_and(self.box_idx[ii-1],mesh.nodes_2d_yg <= np.max(box.box_refxy[ii][1]))
		self.box_idx[ii-1] = np.array(self.box_idx[ii-1])


#___PLOT FESOM2.0 DATA IN INDEX BOX OVER DEPTH AND TIME DATA____________________
# 
#
#_______________________________________________________________________________
def plot_index_t_x_z(mesh,data,box):
	from set_inputarray import inputarray
	
	fsize=14
	fig = plt.figure(figsize=(20, 13))
	ax1  = plt.gca()
	# duplicate x-axes

	#+_____________________________________________________________________________+
	cnumb= 25
	cmin = np.nanmin(data.value)
	cmax = np.nanmax(data.value)
	cref = cmin + (cmax-cmin)/2
	cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
	#cref =0.0

	#+_____________________________________________________________________________+
	# if anomaly data	
	if data.anom==True: 
		cmax,cmin,cref = np.nanmax(data.value), np.nanmin(data.value), 0.0
		data.cmap='blue2red'
	#+_____________________________________________________________________________+
	# if predefined color range	
	if len(data.crange)!=0:
		if len(data.crange)==2:
			cmin = np.float(data.crange[0])
			cmax = np.float(data.crange[1])
			cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
		elif len(data.crange)==3:
			cmin = np.float(data.crange[0])
			cmax = np.float(data.crange[1])
			cref = np.float(data.crange[2])
		else:
			print(' this colorrange definition is not supported !!!')
			print('data.crange=[cmin,cmax] or data.crange=[cmin,cmax,cref]')
	
	print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
	cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,data.cmap)
	print('clevel = ',clevel)
	
	do_drawedges=True
	if clevel.size>30: do_drawedges=False

	# overwrite discrete colormap
	#cmap0 = cmocean.cm.balance
	#do_drawedges=False

	#+_____________________________________________________________________________+
	# make pcolor or contour plot 
	depth = mesh.zlev[0:-1] + (mesh.zlev[1::]-mesh.zlev[0:-1])/2.0
	depth = -depth
	depthlim = np.sum(~np.isnan(data.value[0,:])).max()
	if depthlim==depth.shape: depthlim=depthlim-1
	yy,xx = np.meshgrid(depth,data.time)
	data_plot = data.value
	data_plot[data_plot<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
	data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps


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

	#+_____________________________________________________________________________+
	# set main axes
	ax1.set_xlim(data.time.min(),data.time.max())
	ax1.set_ylim(0,depth[depthlim-1])
	ax1.invert_yaxis()
	ax1.set_axis_bgcolor([0.25,0.25,0.25])
	ax1.tick_params(axis='both',which='major',direction='out',length=8)
	ax1.minorticks_on()
	ax1.tick_params(axis='both',which='minor',direction='out',length=4)
	ax1.set_xlabel('Time [years]',fontdict=dict(fontsize=12))
	ax1.set_ylabel('Depth [km]',fontdict=dict(fontsize=12))
	plt.title(data.descript+' - '+box.box_refxy[0][2]+'\n',fontdict= dict(fontsize=24),verticalalignment='bottom')

	#+_____________________________________________________________________________+
	# draw colorbar
	divider = make_axes_locatable(ax1)
	cax     = divider.append_axes("right", size="2.5%", pad=0.5)
	plt.clim(clevel[0],clevel[-1])

	cbar = plt.colorbar(hp,ax=ax1,cax=cax,ticks=clevel,drawedges=do_drawedges)
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

	#+_____________________________________________________________________________+
	plt.show(block=False)
	
	return(fig,ax1,cax)
	print('finish')

