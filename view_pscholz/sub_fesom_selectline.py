import sys
import numpy as np
import matplotlib.pyplot as plt
from sub_fesom_mesh 		import * 
from sub_regriding_adapt 	import *

def fesom_init_line():

	global line 
	line 			= fesom_line()
	line.line_refxy = []
	#line.line_refxy.append([])
	#line.line_ipxy  = s
	#line.line_ipdist= []
	return line
	
#+_____________________________________________________________________________+
class fesom_line:
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
		self.line       = []
		self.linex      = []
		self.liney      = []
		self.line_refxy = []
		self.line_refdist= []
		self.line_ipxy  = []
		self.line_ipdist= []#
		self.line_data  = []
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
	# LineBuilder
	def Linebuilder(self,x,y,mode='add'):
		if np.size(x)==0 or np.size(y)==0: 
			self.linex = list(self.line.get_xdata())
			self.liney = list(self.line.get_ydata())
			plt.xscale('linear')
			plt.yscale('linear')
			return
		if mode == 'add':
			self.linex.append(np.float16(x))
			self.liney.append(np.float16(y))
		elif mode == 'remove':
			self.linex.remove(self.linex[-1])
			self.liney.remove(self.liney[-1])
			
		mlinex,mliney = self.map(self.linex,self.liney)
		self.line.set_data(mlinex, mliney)
		self.line.set_axes(self.ax)
		self.line.figure.canvas.draw()
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
				elif self.press=='Line':
					self.Linebuilder(x,y)
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
				
			#___LINE MODE [l]___________________________________________________
			elif k=='l':
				#---------------------------------------------------------------
				# activate line mode, start line 
				if self.press == 'None':
					print('draw line: [ON]')
					self.text.set_text('draw line: [ON]')
					plt.draw(),self.figure.canvas.draw()
					
					self.press = 'Line'
					self.line, = plt.plot([], [],color='w',linewidth=2,marker='o',axes=self.ax)
					self.Linebuilder([],[])
					#plt.xscale('linear'),
					#plt.yscale('linear')
					self.figure.canvas.draw()
				#---------------------------------------------------------------
				# switch off line mode, end line
				elif self.press =='Line': 
					print('draw line: [OFF]')
					self.text.set_text('draw line: [OFF]')
					plt.draw(),self.figure.canvas.draw()
					self.press = 'None'
					self.line.set_color([0.0,0.8,0.0])
					#self.line.set_color('w')
					self.line.set_linewidth(0.5)
					plt.xscale('linear')
					plt.yscale('linear')
					self.figure.canvas.draw()
					#ask for name 
					#self.text.set_text('enter name in command line')
					#plt.draw(),self.figure.canvas.draw()
					#name = raw_input("Enter Name of line: ")
					name = 'default'
					self.line_refxy.append([])
					self.line_refxy[len(self.line_refxy)-1]=[self.linex,self.liney,name]
					
			#___DELETE DRAW POINTS [d]__________________________________________
			elif k=='d':
				if self.press=='Line':
					self.Linebuilder(x,y,mode='remove')
		
	#+_________________________________________________________________________+
	# analyse slected lines and calculate interpolation points ip
	def analyse_lines(self,npoints=50):
		#self.line_ipxy  =[[],[]]
		#self.line_ipdist=[0.0]
		
		#+_____________________________________________________________________+
		self.line_ipxy   = []
		self.line_ipdist = []
		self.line_refdist= []
		for ii in range(0,len(self.line_refxy)):
			
			#+_________________________________________________________________+
			self.line_ipxy.append([ [],[] ])
			self.line_ipdist.append([])
			self.line_ipdist[ii].append(0.0)
			self.line_refdist.append([])
			for jj in range(0,len(self.line_refxy[ii][0])-1):
				
				#+_____________________________________________________________+
				P1    = [self.line_refxy[ii][0][jj],self.line_refxy[ii][1][jj]]
				evec  = [self.line_refxy[ii][0][jj+1]-self.line_refxy[ii][0][jj],
						 self.line_refxy[ii][1][jj+1]-self.line_refxy[ii][1][jj],]
				evec  = np.array(evec)
				evecn = (evec[0]**2+evec[1]**2)**0.5
				evec  = evec/evecn
				evecnl = np.linspace(0,evecn,npoints)
				
				#+_____________________________________________________________+
				for kk in range(0,npoints):
					# lon/lat coordinates of interpolation points
					self.line_ipxy[ii][0].append(P1[0]+evecnl[kk]*evec[0])
					self.line_ipxy[ii][1].append(P1[1]+evecnl[kk]*evec[1])
					# calc distance of interpolation points from start point
					if len(self.line_ipxy[ii][1])>1:
						#px1,py1,pz1=geo2cart(self.line_ipxy[0][-1],self.line_ipxy[1][-1])
						#px2,py2,pz2=geo2cart(self.line_ipxy[0][-2],self.line_ipxy[1][-2])
						#dphi = np.arccos((px1*px2+py1*py2+pz1*pz2)/Rearth/Rearth)*180.0/np.pi
						dphi = np.sqrt((self.line_ipxy[ii][0][-1]-self.line_ipxy[ii][0][-2])**2.0+
									   (self.line_ipxy[ii][1][-1]-self.line_ipxy[ii][1][-2])**2.0)
						Rearth=6371.0
						dr = np.pi*Rearth*dphi/180.0
						self.line_ipdist[ii].append(self.line_ipdist[ii][-1]+dr)
					if kk==0:
						self.line_refdist[ii].append(self.line_ipdist[ii][-1])
			self.line_refdist[ii].append(self.line_ipdist[ii][-1])
	#+_________________________________________________________________________+
	# interpolate line points on fesom data
	def interp_lines(self,mesh,data):
		
		#+_____________________________________________________________________+
		self.line_data  = [];
		for ii in range(0,len(self.line_refxy)):
			
			self.line_data.append([])
			npts    = len(line.line_ipxy[ii][0])
			linedata = np.zeros((npts,mesh.nlev-1))
			
			for di in range(0,mesh.nlev-1):
				data_di   = np.array(data.value[:,di])
				data_di[di>=mesh.nodes_2d_iz]=np.nan
				data_di[data_di==0]=np.nan
				
				linedata[:,di]=fesom2regular(data_di, mesh, 
								np.array(line.line_ipxy[ii][0]), 
								np.array(line.line_ipxy[ii][1]))
			
			self.line_data[ii] = linedata