# Patrick Scholz, 23.01.2018
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from sub_fesom_mesh 		import * 
from sub_fesom_plot 		import *
from sub_regriding_adapt 	import *
from colormap_c2c			import *

#+_____________________________________________________________________________+
class fesom_line:
	
	which_obj					= 'line'
	
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
	
	##____selection variable______________________
	#x, y, xc, yc				= [], [], [], []
	#b, k, cid_pressb, cid_pressk= [], [], [] ,[]
	#zoomfac, press				= 20.0, 'None'
	#fig, ax, xlim, ylim			= [], [], [], []
	#ptsx, ptsy, drawline 		= [], [], []
		
	#____line variable___________________________
	line_define, line_define_dist	 = [], []
	line_interp_p, line_interp_pm    = [], []
	line_interp_p_dist				 = []
	line_interp_pm_dist				 = []
	line_interp_pm_dr 				 = []
	line_interp_pm_nvec				 = []
	line_interp_pm_evec				 = []
	
	#___INIT LINE OBJECT_______________________________________________________#
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
	# LineBuilder
	def _Linebuilder_(self,x,y,mode='add'):
		if np.size(x)==0 or np.size(y)==0: 
			self._ptsx = list(self._drawline.get_xdata())
			self._ptsy = list(self._drawline.get_ydata())
			plt.xscale('linear')
			plt.yscale('linear')
			return
		if mode == 'add':
			self._ptsx.append(np.float16(x))
			self._ptsy.append(np.float16(y))
		elif mode == 'remove':
			self._ptsx.remove(self._ptsx[-1])
			self._ptsy.remove(self._ptsy[-1])
			
		mptsx,mptsy = self._map(self._ptsx,self._ptsy)
		self._drawline.set_data(mptsx, mptsy)
		self._drawline.set_axes(self._ax)
		self._drawline.figure.canvas.draw()
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
				elif self._press=='Line':
					self._Linebuilder_(x,y)
			#+___middle mouse button___________________________________________+
			if event.button==2 : return
			if (event.xdata is None): return
			
			#+___right mouse button____________________________________________+
			if event.button==3 : 
				self._drawzoom_out_full()
		
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
				
				# disconnects button and key events
				self._figure.canvas.mpl_disconnect(self._cid_pressb)
				self._figure.canvas.mpl_disconnect(self._cid_pressk)
				
				# This line allows proceeding of the code...
				plt.close(self._figure)
				
			#___LINE MODE [l]___________________________________________________
			elif k=='l':
				#---------------------------------------------------------------
				# activate line mode, start line 
				if self._press == 'None':
					print('draw line: [ON]')
					self._text.set_text('draw line: [ON]')
					plt.draw(),self._figure.canvas.draw()
					
					self._press = 'Line'
					self._drawline, = plt.plot([], [],color='w',linewidth=2,marker='o',axes=self._ax)
					self._Linebuilder_([],[])
					#plt.xscale('linear'),
					#plt.yscale('linear')
					self._figure.canvas.draw()
				#---------------------------------------------------------------
				# switch off line mode, end line
				elif self._press =='Line': 
					print('draw line: [OFF]')
					self._text.set_text('draw line: [OFF]')
					plt.draw(),self._figure.canvas.draw()
					self._press = 'None'
					self._drawline.set_color([0.0,0.8,0.0])
					#self._drawline.set_color('w')
					self._drawline.set_linewidth(2.0)
					plt.xscale('linear')
					plt.yscale('linear')
					self._figure.canvas.draw()
					#ask for name 
					#self._text.set_text('enter name in command line')
					#plt.draw(),self._figure.canvas.draw()
					#name = raw_input("Enter Name of line: ")
					if len(self._ptsx)!=0 :
						name = 'default'
						self.line_define.append([])
						self.line_define[len(self.line_define)-1]=[self._ptsx,self._ptsy,name]
					
			#___DELETE DRAW POINTS [d]__________________________________________
			elif k=='d':
				if self._press=='Line':
					self._Linebuilder_(x,y,mode='remove')
	
	#+_________________________________________________________________________+
	#|																		   |
	#+_________________________________________________________________________+
	# analyse slected lines and calculate interpolation points ip
	def analyse_lines(self,which='res',npoints=50,res=10):
		
		#+_____________________________________________________________________+
		self.line_interp_p   	 = []
		self.line_interp_p_dist  = []
		self.line_interp_pm_dr   = []
		self.line_interp_pm   	 = []
		self.line_interp_pm_dist = []
		self.line_interp_pm_nvec = []
		self.line_interp_pm_evec = [] 
		self.line_define_dist	 = []
		for ii in range(0,len(self.line_define)):
			
			#+_________________________________________________________________+
			self.line_interp_p.append([ [],[] ])
			self.line_interp_p_dist.append([])
			self.line_interp_p_dist[ii].append(0.0)
			self.line_interp_pm_dr.append([])
			self.line_interp_pm.append([ [],[] ])
			self.line_interp_pm_dist.append([])
			self.line_interp_pm_nvec.append([ [],[] ])
			self.line_interp_pm_evec.append([ [],[] ])
			self.line_define_dist.append([])
			for jj in range(0,len(self.line_define[ii][0])-1):
				
				#+_____________________________________________________________+
				P1    = [self.line_define[ii][0][jj],self.line_define[ii][1][jj]]
				P2    = [self.line_define[ii][0][jj+1],self.line_define[ii][1][jj+1]]
				#evec  = [self.line_define[ii][0][jj+1]-self.line_define[ii][0][jj],
						 #self.line_define[ii][1][jj+1]-self.line_define[ii][1][jj],]
				#+_____________________________________________________________+
				# unit vector of line
				evec  = [P2[0]-P1[0],
						 P2[1]-P1[1]]
				evec  = np.array(evec)
				evecn = (evec[0]**2+evec[1]**2)**0.5
				evec  = evec/evecn
				
				if which=='npoints':
					evecnl = np.linspace(0,evecn,npoints)
					loop_pts = npoints
				elif which=='res':
					#Rearth=6371.0
					evecnl = np.arange(0,evecn,res)
					loop_pts=evecnl.size
					
				#_______________________________________________________________
				# normal vector 
				nvec = np.array([-evec[1],evec[0]])
				
				#_______________________________________________________________
				for kk in range(0,loop_pts):
					# lon/lat coordinates of interpolation points
					self.line_interp_p[ii][0].append(P1[0]+evecnl[kk]*evec[0])
					self.line_interp_p[ii][1].append(P1[1]+evecnl[kk]*evec[1])
					
					# lon/lat coordinates of midinterpolation points
					if kk>0:
						aux_mid = evecnl[kk-1]+(evecnl[kk]-evecnl[kk-1])/2
						self.line_interp_pm[ii][0].append(P1[0]+aux_mid*evec[0])
						self.line_interp_pm[ii][1].append(P1[1]+aux_mid*evec[1])
						del aux_mid
					
					# calc distance of interpolation points from start point
					if len(self.line_interp_p[ii][1])>1:
						Rearth=6371.0
						xc,yc,zc = geo2cart([self.line_interp_p[ii][0][-2],self.line_interp_p[ii][0][-1]],
											[self.line_interp_p[ii][1][-2],self.line_interp_p[ii][1][-1]])
						dr = np.pi/180*Rearth*np.degrees(np.arccos( (xc[0]*xc[1]+yc[0]*yc[1]+zc[0]*zc[1])/(Rearth**2) ))
						self.line_interp_p_dist[ii].append(self.line_interp_p_dist[ii][-1]+dr)
						self.line_interp_pm_dr[ii].append(dr)
						
					# calc distance of mid-interpolation points from start point
					if len(self.line_interp_pm[ii][1])>0:
						Rearth=6371.0
						xc,yc,zc = geo2cart([self.line_interp_pm[ii][0][-1],self.line_interp_p[ii][0][-1]],
											[self.line_interp_pm[ii][1][-1],self.line_interp_p[ii][1][-1]])
						dr = np.pi/180*Rearth*np.degrees(np.arccos( (xc[0]*xc[1]+yc[0]*yc[1]+zc[0]*zc[1])/(Rearth**2) ))
						self.line_interp_pm_dist[ii].append(self.line_interp_p_dist[ii][-2]+dr)
						
						self.line_interp_pm_nvec[ii][0].append(nvec[0])
						self.line_interp_pm_nvec[ii][1].append(nvec[1])
						self.line_interp_pm_evec[ii][0].append(evec[0])
						self.line_interp_pm_evec[ii][1].append(evec[1])
						
					if kk==0:
						self.line_define_dist[ii].append(self.line_interp_p_dist[ii][-1])
			self.line_define_dist[ii].append(self.line_interp_p_dist[ii][-1])
	
	#+_________________________________________________________________________+
	#|																		   |
	#+_________________________________________________________________________+
	# interpolate line points on fesom data
	def interp_lines(self,mesh,usemidpts=True):
		
		#+_____________________________________________________________________+
		bckp_value,bckp_value2,bckp_value3 = [],[],[]
		bckp_value,self.value  = self.value,[]
		if len(self.value2)!=0 : bckp_value2,self.value2  = self.value2,[]
		if len(self.value3)!=0 : bckp_value3,self.value3  = self.value3,[]
		
		for ii in range(0,len(self.line_define)):
			
			#___________________________________________________________________
			if usemidpts==True:
				npts    = len(self.line_interp_pm[ii][0])
			else:
				npts    = len(self.line_interp_p[ii][0])
			nlev    = mesh.nlev-1
			
			#___________________________________________________________________
			self.value.append([])
			self.value[ii]=np.zeros((npts,nlev))
			if len(bckp_value2)!=0: 
				self.value2.append([]) 
				self.value2[ii]=np.zeros((npts,nlev))
			if len(bckp_value3)!=0: 
				self.value3.append([])
				self.value3[ii]=np.zeros((npts,nlev))
			
			#___________________________________________________________________
			for di in range(0,nlev):
				data_di   = np.copy(bckp_value[:,di])
				if bckp_value.shape[0]==mesh.n2dna:
					data_di[di>=mesh.nodes_2d_iz-1]=np.nan
				elif bckp_value.shape[0]==mesh.n2dea:
					data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
				#data_di[data_di==0]=np.nan
				
				which='node'
				if data_di.shape[0]==mesh.n2dea: which='elem'
				if usemidpts==True:
					distances, inds = create_indexes_and_distances(mesh,
													np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
													k=10, n_jobs=2,which=which)
					self.value[ii][:,di]=fesom2regular(data_di, mesh,
									np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
									distances=distances, inds=inds)
									#how='idist',k=10)
				else:
					distances, inds = create_indexes_and_distances(mesh,
													np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),\
													k=10, n_jobs=2,which=which)
					self.value[ii][:,di]=fesom2regular(data_di, mesh,
									np.array(self.line_interp_p[ii][0]),np.array(self.line_interp_p[ii][1]),
									distances=distances, inds=inds)
					
								
				#_______________________________________________________________
				# also interpolate second value
				if len(bckp_value2)!=0:
					data_di   = np.copy(bckp_value2[:,di])
					if bckp_value2.shape[0]==mesh.n2dna:
						data_di[di>=mesh.nodes_2d_iz-1]=np.nan
					elif bckp_value2.shape[0]==mesh.n2dea:
						data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
					#data_di[data_di==0]=np.nan
					
					if usemidpts==True:
						self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
										np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
										distances=distances, inds=inds)
										#how='idist',k=10)
					else:
						self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
									np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
									distances=distances, inds=inds)
									#how='idist',k=10)
									
				#_______________________________________________________________
				# also interpolate third value
				if len(bckp_value3)!=0:
					data_di   = np.copy(bckp_value3[:,di])
					if bckp_value3.shape[0]==mesh.n2dna:
						data_di[di>=mesh.nodes_2d_iz-1]=np.nan
					elif bckp_value3.shape[0]==mesh.n2dea:
						data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
					#data_di[data_di==0]=np.nan
					
					which='node'
					if data_di.shape[0]==mesh.n2dea: which='elem'
					if usemidpts==True:
						distances, inds = create_indexes_and_distances(mesh,
													np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
													k=10, n_jobs=2,which=which)
					
						self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
										np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
										distances=distances, inds=inds)
										#how='idist',k=10)
					else:
						distances, inds = create_indexes_and_distances(mesh,
													np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
													k=10, n_jobs=2,which=which)
					
						self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
									np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
									distances=distances, inds=inds)
									#how='idist',k=10)
			
			#___________________________________________________________________
			# copy metainfo from data object to line object
			#self.var 		= data.var
			#self.runid 		= data.runid
			#self.path 		= data.path
			
			#self.descript	= data.descript
			#self.sname 		= data.sname
			#self.lname 		= data.lname
			#self.unit 		= data.unit
			#self.str_time 	= data.str_time
			#self.str_dep 	= data.str_dep
			
			#self.cmap		= data.cmap
			#self.cnumb		= data.cnumb
			#self.crange		= data.crange
			
			#self.anom 		= data.anom
			
			#self.depth 		= data.depth
	
	#+_________________________________________________________________________+
	#|																		   |
	#+_________________________________________________________________________+
	# calculate fluxes through section
	def calc_flux(self, Tref=[], Sref=35.0):
		
		#_______________________________________________________________________
		# calculate volume flux trough section 
		if self.var == 'vec_uv':
			print(' --> calculate volume flux through section')
			self.valueflx = []
			for ii in range(0,len(self.value)):
				self.valueflx.append([])
				self.valueflx[ii]=np.zeros(self.value[ii].shape)
				#_______________________________________________________________
				# loop over zlevel
				for jj in range(0,self.value[ii].shape[1]):
					self.valueflx[ii][:,jj] = (self.value[ii][:,jj]*self.line_interp_pm_nvec[ii][0] + \
									  self.value2[ii][:,jj]*self.line_interp_pm_nvec[ii][1])*self.line_interp_pm_dr[ii]
				
				#_______________________________________________________________
				# loop over sample points
				dz = self.depth[1:]-self.depth[:-1]
				for jj in range(0,self.value[ii].shape[0]):
					self.valueflx[ii][jj,:] = self.valueflx[ii][jj,:] * -dz/1000
				
				#_______________________________________________________________
				self.unit='[Sv]'	
				self.sname='vtransp'	
				self.lname='Volume Transport'
				
		#_______________________________________________________________________
		# calculate heat flux trough section 
		if self.var == 'vec_tuv':
			print(' --> calculate heat flux through section: Tref='+str(Tref))
			rho0=1000.0 
			c_p =3996.0
			self.valueflx = []
			for ii in range(0,len(self.value)):
				self.valueflx.append([])
				self.valueflx[ii]=np.zeros(self.value[ii].shape)
				#_______________________________________________________________
				# loop over zlevel
				for jj in range(0,self.value[ii].shape[1]):
					if len(Tref)==0:
						# calc absolute heatflux 
						self.valueflx[ii][:,jj] =  (self.value[ ii][:,jj]*self.line_interp_pm_nvec[ii][0] + \
													self.value2[ii][:,jj]*self.line_interp_pm_nvec[ii][1]) * \
													self.value3[ii][:,jj] * \
													self.line_interp_pm_dr[ii]
												#(self.value3[ii][:,jj]+273.15) 
					else:
						# calc net heat flux with respect to reference temperature in deg C
						self.valueflx[ii][:,jj] =  (self.value[ ii][:,jj]*self.line_interp_pm_nvec[ii][0] + \
													self.value2[ii][:,jj]*self.line_interp_pm_nvec[ii][1]) * \
													(self.value3[ii][:,jj]-Tref) * \
													self.line_interp_pm_dr[ii]
												
				#(self.value3[ii][:,jj]-Tref) 
				#_______________________________________________________________
				# loop over sample points
				dz = self.depth[1:]-self.depth[:-1]
				for jj in range(0,self.value[ii].shape[0]):
					self.valueflx[ii][jj,:] = self.valueflx[ii][jj,:] * -dz * 1000 * c_p * rho0 * 1.0e-15
					
				#_______________________________________________________________
				self.unit='[PW]'	
				self.sname='hflux'	
				self.lname='Heat Flux'
			
		#_______________________________________________________________________
		# calculate liquit freshwater flux trough section 
		if self.var == 'vec_suv':
			print(' --> calculate liquit freshwater flux through section: Sref='+str(Sref))
			rho0=1000.0 
			self.valueflx = []
			for ii in range(0,len(self.value)):
				self.valueflx.append([])
				self.valueflx[ii]=np.zeros(self.value[ii].shape)
				#_______________________________________________________________
				# loop over zlevel
				for jj in range(0,self.value[ii].shape[1]):
						# calc absolute heatflux 
						self.valueflx[ii][:,jj] =  (self.value[ ii][:,jj]*self.line_interp_pm_nvec[ii][0] + \
													self.value2[ii][:,jj]*self.line_interp_pm_nvec[ii][1]) * \
													(Sref-self.value3[ii][:,jj])/Sref * \
													self.line_interp_pm_dr[ii]
												#(self.value3[ii][:,jj]+273.15) 
					
				#_______________________________________________________________
				# loop over sample points
				dz = self.depth[1:]-self.depth[:-1]
				for jj in range(0,self.value[ii].shape[0]):
					self.valueflx[ii][jj,:] = self.valueflx[ii][jj,:] * -dz * rho /1000
					
				#_______________________________________________________________
				self.unit='[10^-3 Sv]'	
				self.sname='lfwflux'	
				self.lname='Liquit Freshwater Flux'
				
	#+_________________________________________________________________________+
	#|																		   |
	#+_________________________________________________________________________+
	# calculate fluxes through section
	def data_anom(self, line, line2):
		self.descript     = line2.descript+'-'+line.descript
		if line.str_time!=line2.str_time:
			self.str_time = line2.str_time+'-'+line.str_time
		else:
			self.str_time = line.str_time
		self.str_dep = line.str_dep
		#_______________________________________________________________________
		self.value   = []
		self.value2  = []
		self.value3  = []
		self.valueflx= []
		if len(line2.value)!=0:
			for ii in range(0,len(line2.value)):
				self.value.append([])
				self.value[ii] = line2.value[ii]-line.value[ii]
				if len(line2.value2)!=0:
					self.value2.append([])
					self.value2[ii]  = line2.value2[ii]-line.value2[ii]
				if len(line2.value3)!=0:
					self.value3.append([])
					self.value3[ii]  = line2.value3[ii]-line.value3[ii]
				if len(line2.valueflx)!=0:
					self.valueflx.append([])
					self.valueflx[ii]= line2.valueflx[ii]-line.valueflx[ii]
			
		self.anom  = True
		self.crange= []
		self.cmap  = 'blue2red'
		self.depth = line.depth
		self.sname = line.sname 
		self.lname = line.lname
		self.unit  = line.unit
	#+_________________________________________________________________________+
	#|						 PLOT LINE OVER DEPTH 							   |
	#+_________________________________________________________________________+
	def plot_lines_dist_x_z(self):
		from set_inputarray import inputarray
	
		# loop over number of drawn section lines
		for ii in range(0,len(self.line_define)):
			
			#___________________________________________________________________
			fsize=14
			fig = plt.figure(figsize=(20, 13))
			ax1  = plt.gca()
			# duplicate x-axes
			ax2  = ax1.twiny()
			
			#___________________________________________________________________
			if self.var.find('vec')!=-1:
				data_plot = self.valueflx[ii]
			else:
				data_plot = self.value[ii]
			
			#___________________________________________________________________
			cnumb= 25
			cmin = np.nanmin(data_plot)
			cmax = np.nanmax(data_plot)
			cref = cmin + (cmax-cmin)/2
			cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
			#cref =0.0
			cmap = 'grads'
			
			if self.var.find('vec')!=-1: cref, cmap=0.0, 'grads'
			if self.anom==True: cref, cmap=0.0, 'blue2red'
			
			#cref, cmap=0.0, 'blue2red'
			
			cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,cmap)
			do_drawedges=True
			if clevel.size>30: do_drawedges=False
			
			# overwrite discrete colormap
			#cmap0 = cmocean.cm.balance
			#do_drawedges=False
			
			#___________________________________________________________________
			# make pcolor or contour plot 
			xxx = np.array(self.line_interp_pm_dist[ii])
			#xxx1= np.array(self.line_ipxy[ii][0])
			yyy = self.depth[0:-1] + (self.depth[1::]-self.depth[0:-1])/2.0
			yyy = -yyy
			yyylim = np.sum(~np.isnan(data_plot),axis=1).max()+1
			if yyylim==yyy.shape: yyylim=yyylim-1
			yy,xx = np.meshgrid(yyy,xxx)
			
			if inputarray['which_plot']=='pcolor':
				hp=plt.pcolormesh(xx,yy,data_plot,
							shading='flat',#flat
							antialiased=False,
							edgecolor='None',
							cmap=cmap0,
							vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
			else: 
				hp=plt.contourf(xx,yy,data_plot,levels=clevel,
							antialiased=False,
							cmap=cmap0,
							vmin=clevel[0], vmax=clevel[-1])
							#vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
				#plt.contour(xx,yy,data_plot,levels=clevel,
							#antialiased=True,
							#colors='k',
							#linewidths=0.5,
							#linestyles='solid',
							#vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
			hp.cmap.set_under([0.4,0.4,0.4])
			
			#___________________________________________________________________
			# plot position of section points
			plt.plot(self.line_define_dist[ii][0] ,0.0,        color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
			plt.plot(self.line_define_dist[ii][0] ,yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
			plt.plot(self.line_define_dist[ii][-1],0.0,        color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
			plt.plot(self.line_define_dist[ii][-1],yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
			if len(self.line_define_dist[ii])>2:
				for jj in range(1,len(self.line_define_dist[ii])-1):
					plt.plot(np.ones((2,))*self.line_define_dist[ii][jj],np.array([0.0,yyy[yyylim]]),color='k',linewidth=0.5,linestyle='--',marker=None)
					plt.plot(self.line_define_dist[ii][jj],0.0,        color='k',linewidth=0.5,marker="^",markersize=8,mfc='k', clip_on=False)
					plt.plot(self.line_define_dist[ii][jj],yyy[yyylim],color='k',linewidth=0.5,marker="v",markersize=8,mfc='k', clip_on=False)
			
			#___________________________________________________________________
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
			plt.title(self.descript+' - '+self.line_define[0][2]+'\n',fontdict= dict(fontsize=24),verticalalignment='bottom')
			
			#___________________________________________________________________
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
			for jj in range(0,len(self.line_define[ii][0])): 
				strlon='E' 
				if self.line_define[ii][0][jj]<0:strlon='W'
				strlat='N' ; 
				if self.line_define[ii][1][jj]<0:strlon='S'
				lonticklabels.append('{:2.2f}$^{{\\circ}}${:s}\n{:2.2f}$^{{\\circ}}${:s}'.format(np.abs(self.line_define[ii][0][jj]),strlon,np.abs(self.line_define[ii][1][jj]),strlat))
			ax2.set_xticks(self.line_define_dist[ii])
			ax2.set_xticklabels(lonticklabels,fontdict=dict(fontsize=10))
			ax2.set_xlabel('Longitude/Latitude [deg]',fontdict=dict(fontsize=12),verticalalignment='top')
			ax2.tick_params(axis='both',which='major',direction='out',length=8)
			
			#___________________________________________________________________
			# draw colorbar
			divider = make_axes_locatable(ax1)
			cax     = divider.append_axes("right", size="2.5%", pad=0.5)
			plt.clim(clevel[0],clevel[-1])
			
			cbar = plt.colorbar(hp,ax=[ax1,ax2],cax=cax,ticks=clevel,drawedges=do_drawedges)
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
			# bug fix workaround to proper draw secondary axes
			fig.canvas.draw()
			ax2.set_position(ax1.get_position())
			
			#___________________________________________________________________
			# save figure
			if inputarray['save_fig']==True:
				print(' --> save figure: png')
				str_times= self.str_time.replace(' ','').replace(':','')
				str_deps = self.str_dep.replace(' ','').replace(',','').replace(':','')
				sfname = 'lineplot_'+self.line_define[ii][2]+'_'+self.descript+'_'+self.sname+'_'+str_times+'_'+str_deps+'.png'
				sdname = inputarray['save_figpath']
				if os.path.isdir(sdname)==False: os.mkdir(sdname)
				plt.savefig(sdname+sfname, \
							format='png', dpi=600, \
							bbox_inches='tight', pad_inches=0,\
							transparent=True,frameon=True)
			
			#___________________________________________________________________
			plt.show(block=False)
		
		
	#+_________________________________________________________________________+
	#|						 PLOT LINE OVER DEPTH 							   |
	#+_________________________________________________________________________+
	def plot_lines_position(self,mesh):
		from set_inputarray import inputarray
		# draw position of line 
		for ii in range(0,len(self.line_define)):
			#___________________________________________________________________
			xmin,xmax = np.min(self.line_define[ii][0]), np.max(self.line_define[ii][0])
			ymin,ymax = np.min(self.line_define[ii][1]), np.max(self.line_define[ii][1])
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
			ax.plot(self.line_define[ii][0],self.line_define[ii][1] ,color='w',linewidth=2,marker='o',mfc='w',mec='k',axes=ax)
			ax.plot(self.line_define[ii][0][0],self.line_define[ii][1][0] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.0,0.8,0.0],mec='k',axes=ax)
			ax.plot(self.line_define[ii][0][-1],self.line_define[ii][1][-1] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.8,0.0,0.0],mec='k',axes=ax)
			
			ax.plot(self.line_interp_pm[ii][0],self.line_interp_pm[ii][1] ,
			color='r',linewidth=4,marker='+',linestyle='None',axes=ax,markersize=6,mfc='r')
			step=10
			ax.quiver(self.line_interp_pm[ii][0][0::step],self.line_interp_pm[ii][1][0::step],
					  self.line_interp_pm_nvec[ii][0][0::step],self.line_interp_pm_nvec[ii][1][0::step],color='r')
			
			#___________________________________________________________________
			# save figure
			if inputarray['save_fig']==True:
				print(' --> save figure: png')
				sfname = 'lineplot_'+self.line_define[ii][2]+'_position'+'.png'
				sdname = inputarray['save_figpath']
				if os.path.isdir(sdname)==False: os.mkdir(sdname)
				plt.savefig(sdname+sfname, \
							format='png', dpi=600, \
							bbox_inches='tight', pad_inches=0,\
							transparent=True,frameon=True)
			
			#___________________________________________________________________
			plt.show(block=False)
			fig.canvas.draw()