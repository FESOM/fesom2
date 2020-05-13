# Patrick Scholz, 23.01.2018
import sys
import os
import numpy as np
import cmocean
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from sub_fesom_mesh         import * 
from sub_fesom_plot         import *
from sub_regriding_adapt     import *
from colormap_c2c            import *

#+_____________________________________________________________________________+
#|                                                                             |
#|                    *** FESOM2.0 LINE-SECTION CLASS ***                      |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_line:
    
    which_obj                    = 'line'
    
    #____data run variables______________________
    var                         = ''
    runid, path, descript       = 'fesom', '', ''
    
    #____data time variables_____________________
    year,  month, record, depth = [], [], [], []
    str_time, str_dep           = '', ''
    
    #____data projection variables_______________
    proj, proj_lon, proj_lat    = '', 0.0, 90.0
    
    #____data description info___________________
    sname, lname, unit          = '', '', ''
    
    #____data plot variables_____________________
    cmap, crange, cnumb         = 'grads', [], []
    which_plot                  = ''
    
    #____data variable___________________________
    value, value2, value3, time = [], [], [], [] 
    valueflx                    = []; # store calculated fluxes (volume, heat, salt)
    which_mean                  = 'monthly'
    anom                        = False
    
    #____line variable___________________________
    line_define, line_define_dist= [], []
    line_interp_p, line_interp_pm= [], []
    line_interp_p_dist           = []
    line_interp_pm_dist          = []
    line_interp_pm_dr            = []
    line_interp_pm_nvec          = []
    line_interp_pm_evec          = []
    bottom                       = []
    
    #____auxilary field for climatology___________
    lon,lat,depth,zlev,fname,zvec= [],[],[],[],[],''
    
    #+___INIT LINE OBJECT______________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def __init__(self):
        #____selection variable______________________
        self._x, self._y, self._xc, self._yc  = [], [], [], []
        self._b, self._k,                       = [] ,[]
        self._cid_pressb, self._cid_pressk      = [], [],
        self._zoomfac, self._press              = 20.0, 'None'
        self._fig, self._ax, self._map          = [], [], []
        self._xlim, self._ylim                  = [], [],
        self._ptsx, self._ptsy, self._drawline= [], [], []
        self._text                              = []
        self.line_define                      = []
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def _connect_(self,fig,ax,map):
        print(' --> start interactive line selection')
        self._figure = fig
        self._ax     = ax
        self._map    = map
        self._xlim   = self._ax.get_xlim()
        self._ylim   = self._ax.get_ylim()
        self._text   = plt.figtext(0.5,0.05,'Interactive Selection',
                            fontsize=24,
                            fontweight='bold',
                            horizontalalignment='center',
                            bbox=dict(facecolor='w',edgecolor='k'),
                            )
        
        self._cid_pressb = self._figure.canvas.mpl_connect('button_press_event', self._anybutton_)
        self._cid_pressk = self._figure.canvas.mpl_connect('key_press_event', self._anykey_)
        
        self._drawline, = plt.plot([], [],color='w',linewidth=2,marker='o')
        
        # This line stops proceeding of the code...
        plt.show(block=True)
    
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def _disconnect_(self):
        self._figure.canvas.mpl_disconnect(self._cid_pressb)
        self._figure.canvas.mpl_disconnect(self._cid_pressk)
        plt.show(block=False)
    
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # change zoom factor to zoom in 
    def _drawzoom_in_(self):    
        print(' zoom_in')
        self._zoomfac = self._zoomfac/2.0
    
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # change zoom factor to zoom out 
    def _drawzoom_out_(self):    
        print(' zoom_out')
        self._zoomfac = self._zoomfac*2.0
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # zoom full out
    def _drawzoom_outfull_(self):    
        print(' zoom_fullout')
        self._ax.set_xlim(self._xlim[0],self._xlim[1])
        self._ax.set_ylim(self._ylim[0],self._ylim[1])
        self._ax.figure.canvas.draw()
        
    #+_________________________________________________________________________+
    #|                                                                         |
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
        self._ax.figure.canvas.draw()
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # update center position with left mouse click 
    def _update_center_(self, xc, yc):
        self._xc, self._yc, = xc, yc
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # center new leftclick point
    def _move_center_(self):
        self._ax.set_xlim(np.max([self._xc-self._zoomfac,-180.0]), np.min([self._xc+self._zoomfac,180.0]))
        self._ax.set_ylim(np.max([self._yc-self._zoomfac, -90.0]), np.min([self._yc+self._zoomfac, 90.0]))
        self._ax.figure.canvas.draw()
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # LineBuilder
    def _Linebuilder_(self,x,y,mode='add'):
        if np.size(x)==0 or np.size(y)==0: 
            self._ptsx = list(self._drawline.get_xdata())
            self._ptsy = list(self._drawline.get_ydata())
            return
        if mode == 'add':
            self._ptsx.append(np.float16(x))
            self._ptsy.append(np.float16(y))
        elif mode == 'remove':
            self._ptsx.remove(self._ptsx[-1])
            self._ptsy.remove(self._ptsy[-1])
            
        mptsx,mptsy = self._map(self._ptsx,self._ptsy)
        self._drawline.set_data(mptsx, mptsy)
        self._drawline.figure.canvas.draw()
        self._drawline.axes.set_xscale('linear')
        self._drawline.axes.set_yscale('linear')
        self._drawline.figure.canvas.draw()
        
        
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # what should happen if a mouse button is clicked when cursor is over axis 
    def _anybutton_(self,event):
        
        if event.inaxes:
            x, y, b, k = event.xdata, event.ydata, event.button, []
            #print('button=%d, xdata=%f, ydata=%f' % ( event.button, event.xdata, event.ydata))
            #____left mouse button______________________________________________
            if event.button==1 : 
                xc, yc = x, y
                if self._press=='None':
                    self._update_center_(xc, yc)
                    self._move_center_()
                elif self._press=='Line':
                    self._Linebuilder_(x,y)
                    
            #____middle mouse button____________________________________________
            if event.button==2 : return
            if (event.xdata is None): return
            
            #____right mouse button_____________________________________________
            if event.button==3 : 
                self._drawzoom_outfull_()
        
    #+_________________________________________________________________________+
    #|                                                                         |
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
                self._disconnect_()
                
                # This line allows proceeding of the code...
                plt.close(self._figure)
                plt.close('all')
                
            #___LINE MODE [l]___________________________________________________
            elif k=='l':
                #---------------------------------------------------------------
                # activate line mode, start line 
                if self._press == 'None':
                    self._press = 'Line'
                    print('draw line: [ON]')
                    self._text.set_text('draw line: [ON]')
                    self._text.figure.canvas.draw()
                    
                    self._Linebuilder_([],[])
                    self._drawline.axes.set_xscale('linear')
                    self._drawline.axes.set_yscale('linear')
                    self._drawline.axes.figure.canvas.draw()
                    
                #---------------------------------------------------------------
                # switch off line mode, end line
                elif self._press =='Line': 
                    print('draw line: [OFF]')
                    self._text.set_text('draw line: [OFF]')
                    self._text.figure.canvas.draw()
                    self._press = 'None'
                    
                    self._drawline.set_color([0.0,0.8,0.0])
                    self._drawline.set_linewidth(2.0)
                    self._drawline.axes.set_xscale('linear')
                    self._drawline.axes.set_yscale('linear')
                    self._drawline.figure.canvas.draw()
                    
                    #ask for name 
                    #self._text.set_text('enter name in command line')
                    #self._text.figure.canvas.draw()
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
    #|                                                                         |
    #+_________________________________________________________________________+
    # analyse slected lines and calculate interpolation points ip
    def analyse_lines(self,which='res',npoints=50,res=10):
        
        #_______________________________________________________________________
        self.line_interp_p        = []
        self.line_interp_p_dist  = []
        self.line_interp_pm_dr   = []
        self.line_interp_pm        = []
        self.line_interp_pm_dist = []
        self.line_interp_pm_nvec = []
        self.line_interp_pm_evec = [] 
        self.line_define_dist     = []
        for ii in range(0,len(self.line_define)):
            
            #___________________________________________________________________
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
                
                #_______________________________________________________________
                P1    = [self.line_define[ii][0][jj],self.line_define[ii][1][jj]]
                P2    = [self.line_define[ii][0][jj+1],self.line_define[ii][1][jj+1]]
                #evec  = [self.line_define[ii][0][jj+1]-self.line_define[ii][0][jj],
                         #self.line_define[ii][1][jj+1]-self.line_define[ii][1][jj],]
                #_______________________________________________________________
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
                    
                    if kk>0:
                        # lon/lat coordinates of midinterpolation points
                        aux_mid = evecnl[kk-1]+(evecnl[kk]-evecnl[kk-1])/2
                        self.line_interp_pm[ii][0].append(P1[0]+aux_mid*evec[0])
                        self.line_interp_pm[ii][1].append(P1[1]+aux_mid*evec[1])
                        del aux_mid
                    
                        # calc distance of interpolation points from start point
                        Rearth=6371.0
                        xc,yc,zc = geo2cart([self.line_interp_p[ii][0][-2],self.line_interp_p[ii][0][-1]],
                                            [self.line_interp_p[ii][1][-2],self.line_interp_p[ii][1][-1]])
                        dr = np.pi/180*Rearth*np.degrees(np.arccos( (xc[0]*xc[1]+yc[0]*yc[1]+zc[0]*zc[1])/(Rearth**2) ))
                        self.line_interp_p_dist[ii].append(self.line_interp_p_dist[ii][-1]+dr)
                        self.line_interp_pm_dr[ii].append(dr)
                        
                        # calc distance of mid-interpolation points from start point
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
    #|                                                                         |
    #+_________________________________________________________________________+
    # interpolate line points on fesom data
    def interp_lines(self,mesh,usemidpts=True):
        
        #_______________________________________________________________________
        bckp_value,bckp_value2,bckp_value3 = [],[],[]
        bckp_value,self.value  = self.value,[]
        if len(self.value2)!=0 : bckp_value2,self.value2  = self.value2,[]
        if len(self.value3)!=0 : bckp_value3,self.value3  = self.value3,[]
        
        #_______________________________________________________________________
        # copy metainfo from data object to line object
        #self.zlev = mesh.zlev
        self.bottom = []
        
        #_______________________________________________________________________
        # loop over all defined cross section lines
        for ii in range(0,len(self.line_define)):
            
            #___________________________________________________________________
            if usemidpts==True:
                npts    = len(self.line_interp_pm[ii][0])
            else:
                npts    = len(self.line_interp_p[ii][0])
            nlev    = mesh.nlev-1
            
            #___________________________________________________________________
            self.value.append([])
            self.value[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value2)!=0: 
                self.value2.append([]) 
                self.value2[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value3)!=0: 
                self.value3.append([])
                self.value3[ii]=np.zeros((npts,nlev))*np.nan
            
            
            
            #___________________________________________________________________
            # precalculate distance ind indices for interpolation
            which='node'
            if bckp_value.shape[0]==mesh.n2dea: which='elem'
            if usemidpts==True:
                distances, inds = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
                                                k=3, n_jobs=2,which=which)
                                                #k=10, n_jobs=2,which=which)
            else:
                distances, inds = create_indexes_and_distances(mesh,
                                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),\
                                                    k=3, n_jobs=2,which=which)
                                                    #k=10, n_jobs=2,which=which)            
            if len(bckp_value3)!=0:
                which='node'
                if bckp_value3.shape[0]==mesh.n2dea: which='elem'
                if usemidpts==True:
                    distances3, inds3 = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
                                                k=3, n_jobs=2,which=which)
                                                #k=10, n_jobs=2,which=which)                    
                else:
                    distances3, inds3 = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                                k=3, n_jobs=2,which=which)
                                                #k=10, n_jobs=2,which=which)                    
                
            #___________________________________________________________________
            for di in range(0,nlev):
                data_di   = np.copy(bckp_value[:,di])
                if bckp_value.shape[0]==mesh.n2dna:
                    #data_di[di>=mesh.nodes_2d_izg-1]=np.nan
                    data_di[di>mesh.nodes_2d_izg-1]=np.nan
                elif bckp_value.shape[0]==mesh.n2dea:
                    data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    
                if usemidpts==True:
                    #aux = fesom2regular(data_di, mesh,
                                    #np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                    #distances=distances, inds=inds)
                    self.value[ii][:,di]=fesom2regular(data_di, mesh,
                                    np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                    distances=distances, inds=inds)
                                    #how='idist',k=10)
                else:
                    self.value[ii][:,di]=fesom2regular(data_di, mesh,
                                    np.array(self.line_interp_p[ii][0]),np.array(self.line_interp_p[ii][1]),
                                    distances=distances, inds=inds)
                #_______________________________________________________________
                # also interpolate second value
                if len(bckp_value2)!=0:
                    data_di   = np.copy(bckp_value2[:,di])
                    if bckp_value2.shape[0]==mesh.n2dna:
                        data_di[di>=mesh.nodes_2d_izg-1]=np.nan
                    elif bckp_value2.shape[0]==mesh.n2dea:
                        data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    if usemidpts==True:
                        self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
                                        np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                        distances=distances, inds=inds)
                    else:
                        self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                    distances=distances, inds=inds)
                #_______________________________________________________________
                # also interpolate third value
                if len(bckp_value3)!=0:
                    data_di   = np.copy(bckp_value3[:,di])
                    if bckp_value3.shape[0]==mesh.n2dna:
                        data_di[di>=mesh.nodes_2d_iz-1]=np.nan
                    elif bckp_value3.shape[0]==mesh.n2dea:
                        data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    if usemidpts==True:
                        self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
                                        np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                        distances=distances3, inds=inds3)
                    else:
                        self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                    distances=distances3, inds=inds3)
           
            
            #___________________________________________________________________
            self.zvec = np.abs(self.zlev[:-1])+(np.abs(self.zlev[1:])-np.abs(self.zlev[:-1]))/2
            
            #___________________________________________________________________
            # calculate bottom patch
            aux = np.ones(self.value[ii].shape,dtype='int16')
            aux[np.isnan(self.value[ii])]=0
            aux = aux.sum(axis=1)
            aux[aux!=0]=aux[aux!=0]-1
            #bottom = np.abs(self.zlev[aux])
            bottom = np.abs(self.zvec[aux])
            # smooth bottom patch
            filt=np.array([1,2,1]) #np.array([1,2,3,2,1])
            filt=filt/np.sum(filt)
            aux = np.concatenate( (np.ones((filt.size,))*bottom[0],bottom,np.ones((filt.size,))*bottom[-1] ) )
            aux = np.convolve(aux,filt,mode='same')
            bottom = aux[filt.size:-filt.size]
            self.bottom.append([])
            self.bottom[ii]=bottom
            
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # interpolate line points on fesom data
    def interp_lines_fesom14(self,mesh,usemidpts=True):
        
        #_______________________________________________________________________
        bckp_value,bckp_value2,bckp_value3 = [],[],[]
        bckp_value,self.value  = self.value,[]
        if len(self.value2)!=0 : bckp_value2,self.value2  = self.value2,[]
        if len(self.value3)!=0 : bckp_value3,self.value3  = self.value3,[]
        
        #_______________________________________________________________________
        # copy metainfo from data object to line object
        #self.zlev = mesh.zlev
        self.bottom = []
        
        #_______________________________________________________________________
        # loop over all defined cross section lines
        for ii in range(0,len(self.line_define)):
            
            #___________________________________________________________________
            if usemidpts==True:
                npts    = len(self.line_interp_pm[ii][0])
            else:
                npts    = len(self.line_interp_p[ii][0])
            nlev    = mesh.nlev-1
            
            #___________________________________________________________________
            self.value.append([])
            self.value[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value2)!=0: 
                self.value2.append([]) 
                self.value2[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value3)!=0: 
                self.value3.append([])
                self.value3[ii]=np.zeros((npts,nlev))*np.nan
            
            
            
            #___________________________________________________________________
            # precalculate distance ind indices for interpolation
            which='node'
            if bckp_value.shape[0]==mesh.e2d: which='elem'
            if usemidpts==True:
                distances, inds = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
                                                k=10, n_jobs=2,which=which)
            else:
                distances, inds = create_indexes_and_distances(mesh,
                                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),\
                                                    k=10, n_jobs=2,which=which)
            
            if len(bckp_value3)!=0:
                which='node'
                if bckp_value3.shape[0]==mesh.e2d: which='elem'
                if usemidpts==True:
                    distances3, inds3 = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_pm[ii][0]), np.array(self.line_interp_pm[ii][1]),
                                                k=10, n_jobs=2,which=which)
                else:
                    distances3, inds3 = create_indexes_and_distances(mesh,
                                                np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                                k=10, n_jobs=2,which=which)
                
            #___________________________________________________________________
            for di in range(0,nlev):
                data_di   = np.copy(bckp_value[:,di])
                #if bckp_value.shape[0]==mesh.n2dna:
                    #data_di[di>=mesh.nodes_2d_izg-1]=np.nan
                #elif bckp_value.shape[0]==mesh.n2dea:
                    #data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    
                if usemidpts==True:
                    #aux = fesom2regular(data_di, mesh,
                                    #np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                    #distances=distances, inds=inds)
                    
                    self.value[ii][:,di]=fesom2regular(data_di, mesh,
                                    np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                    distances=distances, inds=inds)
                                    #how='idist',k=10)
                else:
                    self.value[ii][:,di]=fesom2regular(data_di, mesh,
                                    np.array(self.line_interp_p[ii][0]),np.array(self.line_interp_p[ii][1]),
                                    distances=distances, inds=inds)
                #_______________________________________________________________
                # also interpolate second value
                if len(bckp_value2)!=0:
                    data_di   = np.copy(bckp_value2[:,di])
                    #if bckp_value2.shape[0]==mesh.n2dna:
                        #data_di[di>=mesh.nodes_2d_izg-1]=np.nan
                    #elif bckp_value2.shape[0]==mesh.n2dea:
                        #data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    if usemidpts==True:
                        self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
                                        np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                        distances=distances, inds=inds)
                    else:
                        self.value2[ii][:,di]=fesom2regular(data_di, mesh, 
                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                    distances=distances, inds=inds)
                #_______________________________________________________________
                # also interpolate third value
                if len(bckp_value3)!=0:
                    data_di   = np.copy(bckp_value3[:,di])
                    #if bckp_value3.shape[0]==mesh.n2dna:
                        #data_di[di>=mesh.nodes_2d_iz-1]=np.nan
                    #elif bckp_value3.shape[0]==mesh.n2dea:
                        #data_di[di>=np.concatenate((mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))-1]=np.nan
                    if usemidpts==True:
                        self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
                                        np.array(self.line_interp_pm[ii][0]),np.array(self.line_interp_pm[ii][1]),
                                        distances=distances3, inds=inds3)
                    else:
                        self.value3[ii][:,di]=fesom2regular(data_di, mesh, 
                                    np.array(self.line_interp_p[ii][0]), np.array(self.line_interp_p[ii][1]),
                                    distances=distances3, inds=inds3)
           
            
            #___________________________________________________________________
            self.zvec = np.abs(self.zlev)
            
            #___________________________________________________________________
            # calculate bottom patch
            aux = np.ones(self.value[ii].shape,dtype='int16')
            aux[np.isnan(self.value[ii])]=0
            aux = aux.sum(axis=1)
            aux[aux!=0]=aux[aux!=0]-1
            #bottom = np.abs(self.zlev[aux])
            bottom = np.abs(self.zvec[aux])
            # smooth bottom patch
            filt=np.array([1,2,1]) #np.array([1,2,3,2,1])
            filt=filt/np.sum(filt)
            aux = np.concatenate( (np.ones((filt.size,))*bottom[0],bottom,np.ones((filt.size,))*bottom[-1] ) )
            aux = np.convolve(aux,filt,mode='same')
            bottom = aux[filt.size:-filt.size]
            self.bottom.append([])
            self.bottom[ii]=bottom            
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # interpolate line points on fesom data
    def interp_lines_reg(self,radius_of_influence=100000, k=3, n_jobs=2 ):
        
        #_______________________________________________________________________
        bckp_value,bckp_value2,bckp_value3 = [],[],[]
        bckp_value,self.value  = self.value,[]
        if len(self.value2)!=0 : bckp_value2,self.value2  = self.value2,[]
        if len(self.value3)!=0 : bckp_value3,self.value3  = self.value3,[]
        self.bottom = []
        
        for ii in range(0,len(self.line_define)):
            
            #___________________________________________________________________
            npts    = len(self.line_interp_pm[ii][0])
            nlev    = len(self.zlev)
            
            #___________________________________________________________________
            self.value.append([])
            self.value[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value2)!=0: 
                self.value2.append([]) 
                self.value2[ii]=np.zeros((npts,nlev))*np.nan
            if len(bckp_value3)!=0: 
                self.value3.append([])
                self.value3[ii]=np.zeros((npts,nlev))*np.nan
            
            #___________________________________________________________________
            # precalculate distance ind indices for interpolation
            mlon,mlat         = np.meshgrid(self.lon, self.lat)
            # source coordinates
            xsrc, ysrc, zsrc = lon_lat_to_cartesian(mlon.flatten(), mlat.flatten())
            # destination coordinates
            xdst, ydst, zdst = lon_lat_to_cartesian(self.line_interp_pm[ii][0], self.line_interp_pm[ii][1])
            tree             = cKDTree(list(zip(xsrc, ysrc, zsrc)))
            distances, inds  = tree.query(list(zip(xdst, ydst, zdst)), k = k, n_jobs=n_jobs)
                
            #___________________________________________________________________
            # do horizontal interpoaltion over depth layers
            data_dst = np.zeros((npts,nlev))
            
            #   calc distance weighting
            distances_ma = np.ma.masked_greater(distances, radius_of_influence)
            w = 1.0 / distances_ma**2
            
            # loop over depth
            for di in range(0,nlev):
                
                #______________________________________________________________________
                data_src =  np.copy(bckp_value[di,:,:]).flatten()
                
                #______________________________________________________________________
                data_interp = np.ma.sum(w * data_src[inds], axis=1) / np.ma.sum(w, axis=1)
                #data_interpolated.shape = lons.shape
                data_interp = np.ma.masked_invalid(data_interp)
                #______________________________________________________________________
                self.value[ii][:,di] = np.ma.filled(data_interp,np.nan)
                
                #_______________________________________________________________
                # also interpolate second value
                if len(bckp_value2)!=0:
                    #______________________________________________________________________
                    data_src =  np.copy(bckp_value2[di,:,:]).flatten()   
                    
                    #______________________________________________________________________
                    data_interp = np.ma.sum(w * data_src[inds], axis=1) / np.ma.sum(w, axis=1)
                    #data_interpolated.shape = lons.shape
                    data_interp = np.ma.masked_invalid(data_interp)
                    #______________________________________________________________________
                    self.value2[ii][:,di] = np.ma.filled(data_interp,np.nan)
                    
                #_______________________________________________________________
                # also interpolate third value
                if len(bckp_value3)!=0:
                    #______________________________________________________________________
                    data_src =  np.copy(bckp_value3[di,:,:]).flatten()   
                    
                    #______________________________________________________________________
                    data_interp = np.ma.sum(w * data_src[inds], axis=1) / np.ma.sum(w, axis=1)
                    #data_interpolated.shape = lons.shape
                    data_interp = np.ma.masked_invalid(data_interp)
                    #______________________________________________________________________
                    self.value2[ii][:,di] = np.ma.filled(data_interp,np.nan)
         
            #___________________________________________________________________
            self.zvec = np.abs(self.zlev)
            
            #___________________________________________________________________
            # calculate bottom patch
            aux = np.ones(self.value[ii].shape,dtype='int16')
            aux[np.isnan(self.value[ii])]=0
            aux = aux.sum(axis=1)
            aux[aux!=0]=aux[aux!=0]-1
            bottom = np.abs(self.zvec[aux])
            # smooth bottom patch
            filt=np.array([1,2,1]) #np.array([1,2,3,2,1])
            filt=filt/np.sum(filt)
            aux = np.concatenate( (np.ones((filt.size,))*bottom[0],bottom,np.ones((filt.size,))*bottom[-1] ) )
            aux = np.convolve(aux,filt,mode='same')
            bottom = aux[filt.size:-filt.size]
            self.bottom.append([])
            self.bottom[ii]=bottom 

    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # interpolate line points on fesom data
    def interp_vert(self,levels):
        data_out = np.zeros((len(self.line_interp_pm[0][0]),levels.size))
        for di in range(0,len(levels)):
            # find upper and lower layer indices
            #idx_dwn = np.array(np.where( data.depth[di]<=abs(mesh.zlev))).squeeze()
            idx_dwn = np.array(np.where( abs(levels[di])<=abs(self.zlev))).squeeze()
            if len(idx_dwn)==0: idx_dwn=[len(self.zlev)-1]
            idx_dwn = idx_dwn[0]
            idx_up  = idx_dwn-1
            if idx_up<0: idx_up=0
            
            # linear vertical interpolant
            deltaz   = abs(self.zlev[idx_dwn])-abs(self.zlev[idx_up])
            #deltaz_i = abs(mesh.zlev[idx_dwn])-data.depth[di]
            deltaz_i = abs(self.zlev[idx_dwn])-abs(levels[di])
            
            # interpoalte verticaly and sum up
            if deltaz_i==0:
                auxval = self.value[0][:,idx_dwn]
            else:
                auxval = self.value[0][:,idx_dwn]-(self.value[0][:,idx_dwn]-self.value[0][:,idx_up])*deltaz_i/deltaz
                
            #aux_div[~np.isnan(auxval)]=aux_div[~np.isnan(auxval)]+1.0
            #auxval[np.isnan(auxval)]=0.0
            data_out[:,di] = auxval
        self.value[0] = data_out
        
    #+_________________________________________________________________________+
    #|                                                                         |
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
                dz = self.zlev[1:]-self.zlev[:-1]
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
                dz = self.zlev[1:]-self.zlev[:-1]
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
                dz = self.zlev[1:]-self.zlev[:-1]
                for jj in range(0,self.value[ii].shape[0]):
                    self.valueflx[ii][jj,:] = self.valueflx[ii][jj,:] * -dz * rho /1000
                    
                #_______________________________________________________________
                self.unit='[10^-3 Sv]'    
                self.sname='lfwflux'    
                self.lname='Liquit Freshwater Flux'
                
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # calculate fluxes through section
    def data_anom(self, line2, line):
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
        self.zlev = line.zlev
        self.sname = line.sname 
        self.lname = line.lname
        self.unit  = line.unit
    #+_________________________________________________________________________+
    #|                         PLOT LINE OVER DEPTH                                |
    #+_________________________________________________________________________+
    def plot_lines_dist_x_z(self,numb=[],figsize=[],do_subplot=[],do_output=True,
                            do_xlabel=True,do_ylabel=True,do_title=True,do_cbar=True,
                            allow_save=True,maxdep=[],which_lines=[1,1,1],do_2ndax=True,do_lname=False):
        from set_inputarray import inputarray
        fsize = 12
        
        if len(figsize)==0 : figsize=[12,6]
        #_______________________________________________________________________
        # numb...which cross section should be drawn in case several sections where
        # defined, numb = [] means all are plotted
        if isinstance(numb,int)==True : numb = [numb]
        if len(numb)==0 : numb=range(0,len(self.line_define))
        
        #_______________________________________________________________________
        # loop over number of drawn section lines
        for ii in numb:
            #_______________________________________________________________________
            # plot is not part of subplot
            if len(do_subplot)==0:
                fig,ax1 = plt.figure(figsize=figsize), plt.gca()
            else:
                fig,ax1 = do_subplot[0], do_subplot[1]
                fig.sca(ax1)
            
            if np.unique(self.line_interp_p[0][0]).size!=1 and \
               np.unique(self.line_interp_p[0][1]).size!=1: 
                if do_2ndax: ax2  = ax1.twiny()
            
            #___________________________________________________________________
            if self.var.find('vec')!=-1:
                data_plot = np.copy(self.valueflx[ii])
            else:
                data_plot = np.copy(self.value[ii])
            
            #___________________________________________________________________
            # set colorrange, colormap, number of colors
            cmin,cmax = np.nanmin(data_plot),np.nanmax(data_plot)
            cref = cmin + (cmax-cmin)/2
            cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
            if len(self.crange)!=0: cmin,cmax,cref = self.crange[0],self.crange[1],self.crange[2]
            
            cnumb=25
            if np.size(self.cnumb)!=0: cnumb = self.cnumb
            
            cmap = 'grads'
            if len(self.cmap)!=0: cmap = self.cmap
            
            if self.var.find('vec')!=-1: cref, self.cmap=0.0, 'grads'
            if self.anom==True: cref, self.cmap=0.0, 'blue2red'
            
            #___________________________________________________________________
            cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,cmap)
            do_drawedges=True
            if clevel.size>60: do_drawedges=False
            
            #___________________________________________________________________
            # limit data to colorrange
            data_plot[data_plot<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
            data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
            
            #___________________________________________________________________
            # set x-axis
            if np.unique(self.line_interp_p[0][0]).size==1:
                xvec,str_xlabel = np.array(self.line_interp_pm[ii][1]),'Latitude [deg]'
            elif np.unique(self.line_interp_p[0][1]).size==1:    
                xvec,str_xlabel = np.array(self.line_interp_pm[ii][0]),'Longitude [deg]'
            else:    
                xvec,str_xlabel = np.array(self.line_interp_pm_dist[ii]),'Distance from start point [km]'
            
            #___________________________________________________________________
            # set y-axis
            #if len(self.lon)==0:
                #yvec,str_ylabel = -(self.zlev[0:-1] + (self.zlev[1::]-self.zlev[0:-1])/2.0), 'Depth [m]'
                #ylim = np.sum(~np.isnan(data_plot),axis=1).max()-1
                #if ylim<yvec.shape[0]-1: ylim=ylim+1
            #else:
                #yvec,str_ylabel = self.zlev, 'Depth [m]'
                #ylim = np.sum(~np.isnan(data_plot),axis=1).max()-1
                #if ylim<yvec.shape[0]-1: ylim=ylim+1
            yvec,str_ylabel = self.zvec, 'Depth [m]'
            ylim = np.sum(~np.isnan(data_plot),axis=1).max()-1
            if ylim<yvec.shape[0]-1: ylim=ylim+1
            if np.isscalar(maxdep)==False: maxdep=yvec[ylim]
            
            #___________________________________________________________________
            # make pcolor or contour plot 
            yy,xx = np.meshgrid(yvec,xvec)
            if inputarray['which_plot']=='pcolor':
                #hp=plt.pcolormesh(xx,yy,data_plot,
                hp=ax1.pcolormesh(xx,yy,data_plot,
                            shading='gouraud',#flat
                            antialiased=False,
                            edgecolor='None',
                            cmap=cmap0,
                            clim=[clevel[0],clevel[-1]],
                            vmin=clevel[0],vmax=clevel[-1])
            else: 
                #hp=plt.contourf(xx,yy,data_plot,levels=clevel,
                hp=ax1.contourf(xx,yy,data_plot,levels=clevel,
                            antialiased=False,
                            cmap=cmap0,
                            vmin=clevel[0], vmax=clevel[-1])
                            #vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
                #plt.contour(xx,yy,data_plot,levels=clevel[clevel<=cref],
                if which_lines[0]==1:
                    ax1.contour(xx,yy,data_plot,levels=clevel[clevel<cref],
                                    antialiased=True,
                                    colors='k',
                                    linewidths=0.5,
                                    linestyles='-',
                                    vmin=clevel[0], vmax=clevel[-1])
                if which_lines[1]==1:
                    ax1.contour(xx,yy,data_plot,levels=clevel[clevel==cref],
                                    antialiased=True,
                                    colors='k',
                                    linewidths=1.0,
                                    linestyles='-',
                                    vmin=clevel[0], vmax=clevel[-1])
                if which_lines[2]==1:
                    ax1.contour(xx,yy,data_plot,levels=clevel[clevel>cref],
                                    antialiased=True,
                                    colors='k',
                                    linewidths=0.5,
                                    linestyles='--',
                                    vmin=clevel[0], vmax=clevel[-1])
            
            #___________________________________________________________________
            # plot bottom patch
            ax1.fill_between(xvec,self.bottom[ii], maxdep,color=[0.5,0.5,0.5])#,alpha=0.95)
            ax1.plot(xvec,self.bottom[ii],color='k')
            
            #___________________________________________________________________
            # plot position of section points
            if len(self.line_define_dist[ii])>2:
                plt.plot(xvec[0] ,0.0,       color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
                plt.plot(xvec[0] ,maxdep,color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.0,0.8,0.0], clip_on=False)
                plt.plot(xvec[-1],0.0,       color='k',linewidth=0.5,marker="^",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
                plt.plot(xvec[-1],maxdep,color='k',linewidth=0.5,marker="v",markersize=10,mfc=[0.8,0.0,0.0], clip_on=False)
                for jj in range(1,len(self.line_define_dist[ii])-1):
                    plt.plot(np.ones((2,))*self.line_define_dist[ii][jj],np.array([0.0,maxdep]),color='k',linewidth=0.5,linestyle='--',marker=None)
                    plt.plot(self.line_define_dist[ii][jj],0.0,        color='k',linewidth=0.5,marker="^",markersize=8,mfc='k', clip_on=False)
                    plt.plot(self.line_define_dist[ii][jj],maxdep,color='k',linewidth=0.5,marker="v",markersize=8,mfc='k', clip_on=False)
            
            #___________________________________________________________________
            # set main axes
            ax1.set_xlim(xvec.min(),xvec.max())
            ax1.set_ylim(0,maxdep)
            ax1.invert_yaxis()
            ax1.tick_params(axis='both',which='major',direction='out',length=8,labelsize=fsize)
            ax1.minorticks_on()
            ax1.tick_params(axis='both',which='minor',direction='out',length=4,labelsize=fsize)
            if do_xlabel: 
                ax1.set_xlabel(str_xlabel,fontdict=dict(fontsize=fsize))
            else:
                ax1.set_xticklabels([])
            if do_ylabel: 
                ax1.set_ylabel(str_ylabel,fontdict=dict(fontsize=fsize))
            else:
                ax1.set_yticklabels([])
            #if do_title:  ax1.set_title(self.descript+' '+self.line_define[0][2]+'\n',fontdict= dict(fontsize=fsize+2),verticalalignment='bottom')
            
            txtx, txty = xvec[0]+(xvec[-1]-xvec[0])*0.025,maxdep-(maxdep*0.025)
            if do_lname:
                ax1.text(txtx,txty,self.descript+' '+self.line_define[0][2] , fontsize=14, fontweight='bold',horizontalalignment='left')
            else: 
                ax1.text(txtx,txty,self.descript, fontsize=14, fontweight='bold',horizontalalignment='left')
            #___________________________________________________________________
            # set upper secondary x-axes
            if np.unique(self.line_interp_p[0][0]).size!=1 and \
               np.unique(self.line_interp_p[0][1]).size!=1 and do_2ndax:  
                ax2.set_xlim(xvec.min(),xvec.max())
                ax2.set_ylim(0,maxdep)
                ax2.invert_yaxis()
                ax1.invert_yaxis()
                #_______________________________________________________________
                lonticklabels=[]
                for jj in range(0,len(self.line_define[ii][0])): 
                    strlon='E' 
                    if self.line_define[ii][0][jj]<0:strlon='W'
                    strlat='N' ; 
                    if self.line_define[ii][1][jj]<0:strlat='S'
                    lonticklabels.append('{:2.2f}$^{{\\circ}}${:s}\n{:2.2f}$^{{\\circ}}${:s}'.format((self.line_define[ii][0][jj]),strlon,(self.line_define[ii][1][jj]),strlat))
                ax2.plot(self.line_define_dist[ii][0],0 ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.0,0.8,0.0],mec='k',clip_on=False)
                ax2.plot(self.line_define_dist[ii][1],0 ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.8,0.0,0.0],mec='k',clip_on=False)
                
                ax2.set_xticks(self.line_define_dist[ii])
                ax2.set_xticklabels(lonticklabels,fontdict=dict(fontsize=10))
                #ax2.set_xlabel('Longitude/Latitude [deg]',fontdict=dict(fontsize=fsize),verticalalignment='top')
                ax2.tick_params(axis='both',which='major',direction='out',length=8)
            
            #___________________________________________________________________
            # draw colorbar
            if do_cbar : 
                divider = make_axes_locatable(ax1)
                cax     = divider.append_axes("right", size="2.5%", pad=0.5)
                #plt.clim(clevel[0],clevel[-1])
                #ax1.clim(clevel[0],clevel[-1])
                for im in ax1.get_images():     im.set_clim(clevel1[0],clevel1[-1])
                
                #cbar = plt.colorbar(hp,ax=[ax1,ax2],cax=cax,ticks=clevel,drawedges=do_drawedges)
                cbar = plt.colorbar(hp,ax=[ax1],cax=cax,ticks=clevel,drawedges=do_drawedges,orientation='vertical')
                cbar.set_label(self.lname+' '+self.unit+'\n'+self.str_time, size=fsize)
                
                ncl = 10
                rotation = 0
                if cbar.orientation=='vertical': tickl = cbar.ax.get_yticklabels()
                else:                            tickl = cbar.ax.get_xticklabels()
                ncbar_l=len(tickl[:])
                idx_cref = np.where(clevel==cref)[0]
                idx_cref = np.asscalar(idx_cref)
                nstep = ncbar_l/ncl
                nstep = np.int(np.floor(nstep))
                if nstep==0: nstep=1
                idx = np.arange(0,len(tickl),1)
                idxb = np.ones((len(tickl),), dtype=bool)                
                idxb[idx_cref::nstep]  = False
                idxb[idx_cref::-nstep] = False
                idx = idx[idxb==True]
                for jj in list(idx):
                    tickl[jj]=''
                if cbar.orientation=='vertical':cbar.ax.set_yticklabels(tickl,fontsize=fsize,rotation=rotation)
                else:                           cbar.ax.set_xticklabels(tickl,fontsize=fsize,rotation=rotation)    
    
                
                #cl = plt.getp(cbar.ax, 'ymajorticklabels')
                #plt.setp(cl, fontsize=fsize)
                
                ## kickout some colormap labels if there are to many
                #ncbar_l=len(cbar.ax.get_yticklabels()[:])
                #idx_cref = np.where(clevel==cref)[0]
                #idx_cref = np.asscalar(idx_cref) 
                #nmax_cbar_l = 10
                #nstep = ncbar_l/nmax_cbar_l
                #nstep = np.int(np.floor(nstep))
                #if nstep==0:nstep=1
                #plt.setp(cbar.ax.get_yticklabels()[:], visible=False)
                #plt.setp(cbar.ax.get_yticklabels()[idx_cref::nstep], visible=True)
                #plt.setp(cbar.ax.get_yticklabels()[idx_cref::-nstep], visible=True)
            
            #___________________________________________________________________
            # bug fix workaround to proper draw secondary axes
            fig.canvas.draw()
            if np.unique(self.line_interp_p[0][0]).size!=1 and \
               np.unique(self.line_interp_p[0][1]).size!=1 and do_2ndax: ax2.set_position(ax1.get_position())
            
            #___________________________________________________________________
            # save figure
            if inputarray['save_fig']==True and allow_save==True:
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
            
            #___________________________________________________________________
            if do_cbar : 
                return(fig,ax1,cbar,clevel,hp)
            else:
                return(fig,ax1,clevel,hp)
        
        
    #+_________________________________________________________________________+
    #|                         PLOT LINE POSITION                                |
    #+_________________________________________________________________________+
    def plot_lines_position(self,mesh,cmap,clevel,cref,numb=[],figsize=[],
                            do_subplot=[],do_nvec=False, do_cbar=True, allow_save=True,do_grid=False):
        from set_inputarray import inputarray
        fsize=10
        
        if len(figsize)==0 : figsize=[12,6]
        #_______________________________________________________________________
        if isinstance(numb,int)==True : numb = [numb]
        if len(numb)==0 : numb=range(0,len(self.line_define))
        
        #_______________________________________________________________________
        # draw position of line
        for ii in numb:
            #_______________________________________________________________________
            # plot is not part of subplot
            if len(do_subplot)==0:
                fig,ax = plt.figure(figsize=figsize), plt.gca()
            else:
                fig,ax = do_subplot[0], do_subplot[1]
                fig.sca(ax)
                
            #___________________________________________________________________
            xmin,xmax = np.min(self.line_define[ii][0]), np.max(self.line_define[ii][0])
            ymin,ymax = np.min(self.line_define[ii][1]), np.max(self.line_define[ii][1])
            
            if xmin==xmax: xmin,xmax = xmin-(ymax-ymin)/2,xmin+(ymax-ymin)/2,  
            if ymin==ymax: ymin,ymax = ymin-(xmax-xmin)/2,ymin+(xmax-xmin)/2,  
            
            xmin,xmax,ymin,ymax = xmin-5.0, xmax+5.0, ymin-5.0, ymax+5.0
            xmin,xmax,ymin,ymax = np.max([xmin,-180.0]),np.min([xmax,180.0]),np.max([ymin,-90.0]),np.min([ymax,90.0])
            print(xmin,xmax,ymin,ymax)
            #___________________________________________________________________
            #fig, ax = plt.figure(figsize=(13, 13)), plt.gca()
            map     = Basemap(projection = 'cyl',resolution = 'c',
                        llcrnrlon = xmin, urcrnrlon = xmax, llcrnrlat = ymin, urcrnrlat = ymax)
            mx,my     = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
            
            #___________________________________________________________________
            tri     = Triangulation(mx, my,mesh.elem_2d_i)
            idxbox_e = mesh.nodes_2d_xg[mesh.elem_2d_i].max(axis=1)<xmin
            idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_xg[mesh.elem_2d_i].min(axis=1)>xmax)
            idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].max(axis=1)<ymin)
            idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].min(axis=1)>ymax)
            tri.set_mask(idxbox_e)
            #hp1        = plt.tripcolor(tri,-mesh.nodes_2d_zg,
                        #cmap=cmap, #cm.gist_ncar, #cm.nipy_spectral, #cmocean.cm.deep,
                        #antialiased=False,
                        #edgecolors='None',
                        #shading='gouraud')
            data_plot = -mesh.nodes_2d_zg
            data_plot[data_plot<clevel[ 0]] = clevel[ 0]+np.finfo(np.float32).eps
            data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
            hp1=ax.tricontourf(tri,data_plot,levels=clevel,antialiased=False,extend='both',cmap=cmap)
            if do_grid: ax.triplot(tri,color='k',linewidth=.8,alpha=0.8)
    
            #___________________________________________________________________
            step=np.array([1.0,2.0,5.0,10.0,15.0,30.0,45.0,60.0])
            nxyl = [5,5]
            idx = np.array(np.where(step>=(xmax-xmin)/nxyl[0])[0])
            if idx.size==0 : idx = np.array([step.size-1])
            stepx=step[idx[0]]
            idx = np.array(np.where(step>=(ymax-ymin)/nxyl[1])[0])
            if idx.size==0 : idx = np.array([step.size-1])
            stepy=step[idx[0]]
            xticks , yticks =np.arange(-180.,180.+1,stepx), np.arange(-90.,90.+1,stepy)
            xlabels,ylabels=[0,0,0,1],[1,0,0,0]
            map.drawparallels(yticks,labels=ylabels,fontsize=fsize)
            map.drawmeridians(xticks,labels=xlabels,fontsize=fsize)
            fesom_plot_lmask(map,mesh,ax,'0.6','k')
            #map.bluemarble()
            map.drawmapboundary(fill_color='0.9',linewidth=1.0)
            
            #___________________________________________________________________
            ax.plot(self.line_define[ii][0],self.line_define[ii][1] ,color='w',linewidth=4,marker='o',mfc='w',mec='k',axes=ax)
            ax.plot(self.line_define[ii][0][0],self.line_define[ii][1][0] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.0,0.8,0.0],mec='k',axes=ax)
            ax.plot(self.line_define[ii][0][-1],self.line_define[ii][1][-1] ,color='k',linewidth=0.5,marker='s',markersize=8,mfc=[0.8,0.0,0.0],mec='k',axes=ax)
            ax.plot(self.line_interp_pm[ii][0],self.line_interp_pm[ii][1] ,color='k',marker='.',linestyle='None',linewidth=0.5,markersize=4,axes=ax)
            
            #ax.plot(self.line_interp_pm[ii][0],self.line_interp_pm[ii][1] ,
            #color='r',linewidth=4,marker='+',linestyle='None',axes=ax,markersize=6,mfc='r')
            if do_nvec:
                # plot arrow in middle of each line section element
                for jj in range(0,len(self.line_define_dist[ii])-1):
                    middist = self.line_define_dist[ii][jj]+(self.line_define_dist[ii][jj+1]-self.line_define_dist[ii][jj])/2
                    mididx  = np.abs(np.array(self.line_interp_pm_dist[0]-middist)).argmin()
                    ax.quiver(self.line_interp_pm[ii][0][mididx],self.line_interp_pm[ii][1][mididx],
                        self.line_interp_pm_nvec[ii][0][mididx],self.line_interp_pm_nvec[ii][1][mididx],color='w',scale=10.0)
                plt.title(self.line_define[ii][2],fontdict=dict(fontsize=fsize*2),verticalalignment='bottom')
            
            #___________________________________________________________________
            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("right", size="2.5%", pad=0.1)
            if do_cbar:
                cbar = plt.colorbar(hp1,ax=ax,orientation='vertical',\
                            ticks=clevel,drawedges=True,extend='neither',\
                            extendrect=True,extendfrac='auto')
                ncl=10
                
                cbar.set_label('depth [m]', size=12)
                if cbar.orientation=='vertical': tickl = cbar.ax.get_yticklabels()
                else:                            tickl = cbar.ax.get_xticklabels()
                ncbar_l=len(tickl[:])
                idx_cref = np.where(clevel==cref)[0]
                idx_cref = np.asscalar(idx_cref)
                nstep = ncbar_l/ncl
                nstep = np.int(np.floor(nstep))
                if nstep==0: nstep=1
                idx = np.arange(0,len(tickl),1)
                idxb = np.ones((len(tickl),), dtype=bool)                
                idxb[idx_cref::nstep]  = False
                idxb[idx_cref::-nstep] = False
                idx = idx[idxb==True]
                for jj in list(idx):
                    tickl[jj]=''
                if cbar.orientation=='vertical':cbar.ax.set_yticklabels(tickl)
                else:                           cbar.ax.set_xticklabels(tickl)    
            plt.show(block=False)
            fig.canvas.draw()
            
            #___________________________________________________________________
            # save figure
            if inputarray['save_fig']==True and allow_save:
                print(' --> save figure: png')
                sfname = 'lineplot_'+self.line_define[ii][2]+'_position'+'.png'
                sdname = inputarray['save_figpath']
                if os.path.isdir(sdname)==False: os.mkdir(sdname)
                plt.savefig(sdname+sfname, \
                            format='png', dpi=600, \
                            bbox_inches='tight', pad_inches=0,\
                            transparent=True,frameon=True)
            
            #___________________________________________________________________
            if do_cbar : 
                return(fig,ax,cbar,clevel,hp1)
            else:
                return(fig,ax,clevel,hp1)
            
