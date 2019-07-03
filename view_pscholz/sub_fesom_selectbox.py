# Patrick Scholz, 23.01.2018
import sys
import os
import numpy as np
import copy as  cp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sub_fesom_mesh         import * 
from sub_fesom_plot         import *
from colormap_c2c            import *


#+_____________________________________________________________________________+
#|                                                                             |
#|                    *** FESOM2.0 INDEXBOX CLASS ***                          |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_box:
    
    which_obj                    = 'box'
    
    #____data run variables______________________
    var                         = ''
    runid, path, descript       = 'fesom', '', ''
    
    #____data time variables_____________________
    year,  month, record, depth,zlev = [], [], [], [],[]
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
    
    #____box variable____________________________
    box_define                  = []
    box_idx                     = []
    
    #____additional opttion______________________
    botlimit                    = []
    
    #+___INIT BOX OBJECT_______________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def __init__(self):
        #____selection variable______________________
        self._x, self._y, self._xc, self._yc  = [], [], [], []
        self._b, self._k,                     = [] ,[]
        self._cid_pressb, self._cid_pressk    = [], [],
        self._zoomfac, self._press            = 20.0, 'None'
        self._figure, self._ax, self._map     = [], [], []
        self._xlim, self._ylim                = [], [],
        self._ptsx, self._ptsy, self._drawline= [], [], []
        self._text                            = []
        self.box_define, self.box_idx         = [],[]
        
    #+_________________________________________________________________________+
    #|                                                                         |
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
        
        self._cid_pressb = self._figure.canvas.mpl_connect('button_press_event', self._anybutton_)
        self._cid_pressk = self._figure.canvas.mpl_connect('key_press_event', self._anykey_)
        
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
    # BoxBuilder
    def _Boxbuilder_(self,x,y,mode='add'):
        if np.size(x)==0 or np.size(y)==0: 
            self._ptsx = list(self._drawbox.get_xdata())
            self._ptsy = list(self._drawbox.get_ydata())
            return
        if mode == 'add':
            if len(self._ptsx)==2:
                self._drawbox.set_linestyle('-')
                #dist0 = np.sqrt((self._ptsx[0]-np.float16(x))**2+(self._ptsy[0]-np.float16(y))**2)
                #dist1 = np.sqrt((self._ptsx[1]-np.float16(x))**2+(self._ptsy[1]-np.float16(y))**2)
                dist0x = np.abs(self._ptsx[0]-np.float16(x))
                dist0y = np.abs(self._ptsy[0]-np.float16(y))
                dist1x = np.abs(self._ptsx[1]-np.float16(x))
                dist1y = np.abs(self._ptsy[1]-np.float16(y))
                if dist0x<=dist1x : self._ptsx[0]=np.float16(x)
                else:                self._ptsx[1]=np.float16(x)
                if dist0y<=dist1y : self._ptsy[0]=np.float16(y)
                else:                self._ptsy[1]=np.float16(y)
                    
                #self._ptsy.append(np.float16(y))
                #self._ptsx.remove(self._ptsx[0])
                #self._ptsy.remove(self._ptsy[0])
                #self._ptsx.append(np.float16(x))
                #self._ptsy.append(np.float16(y))
                
            else:
                self._ptsx.append(np.float16(x))
                self._ptsy.append(np.float16(y))
            
            #if len(self._ptsx)>=4:
                #self._ptsx = [np.min(self._ptsx),np.max(self._ptsx),np.max(self._ptsx),\
                              #np.min(self._ptsx),np.min(self._ptsx)]
                #self._ptsy = [np.min(self._ptsy),np.min(self._ptsy),np.max(self._ptsy),\
                              #np.max(self._ptsy),np.min(self._ptsy)]
                #self._drawbox.set_linestyle('-')
        elif mode == 'remove':
            self._ptsx=[]
            self._ptsy=[]
            self._drawbox.set_linestyle('None')
            
        auxboxx = [np.min(self._ptsx),np.max(self._ptsx),np.max(self._ptsx),\
                   np.min(self._ptsx),np.min(self._ptsx)]
        auxboxy = [np.min(self._ptsy),np.min(self._ptsy),np.max(self._ptsy),\
                  np.max(self._ptsy),np.min(self._ptsy)]    
        
        self._drawbox.set_data(auxboxx,auxboxy)
        self._drawbox.set_axes(self._ax)
        self._drawbox.figure.canvas.draw()
        self._drawbox.axes.set_xscale('linear')
        self._drawbox.axes.set_yscale('linear')
        self._drawbox.axes.figure.canvas.draw()
    
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # PolygonBuilder
    def _Polygonbuilder_(self,x,y,mode='add'):
        if np.size(x)==0 or np.size(y)==0: 
            self._ptsx = list(self._drawbox.get_xdata())
            self._ptsy = list(self._drawbox.get_ydata())
            return
        #_______________________________________________________________________
        if mode == 'add':
            self._ptsx.append(np.float16(x))
            self._ptsy.append(np.float16(y))
            if len(self._ptsx)>=2:
                self._drawbox.set_linestyle('-')
        #_______________________________________________________________________
        elif mode == 'remove':
            self._ptsx.remove(self._ptsx[-1])
            self._ptsy.remove(self._ptsy[-1])
            
        #_______________________________________________________________________
        # what happens if polygon is closed
        if np.sqrt((self._ptsx[0]-self._ptsx[-1])**2 + (self._ptsy[0]-self._ptsy[-1])**2)<=0.5 and len(self._ptsx)>=3:
            #___________________________________________________________________
            self._text.set_text('draw polygon: [OFF]')
            self._text.figure.canvas.draw()
            self._press = 'None'
            
            #___________________________________________________________________
            self._ptsx.remove(self._ptsx[-1])
            self._ptsy.remove(self._ptsy[-1])
            self._ptsx.append(self._ptsx[0])
            self._ptsy.append(self._ptsy[0])
            #mboxx,mboxy = self._map(self._ptsx,self._ptsy)
            #self._drawbox.set_data(mboxx, mboxy)
            self._drawbox.set_data(self._ptsx,self._ptsy)
            self._drawbox.set_color([0.0,0.8,0.0])
            self._drawbox.set_linewidth(2.0)
            self._drawbox.axes.set_xscale('linear')
            self._drawbox.axes.set_yscale('linear')
            self._drawbox.axes.figure.canvas.draw()
            
            #___________________________________________________________________
            #ask for name 
            #self._text.set_text('enter name in command line')
            #plt.draw(),self._figure.canvas.draw()
            #name = raw_input("Enter Name of line: ")
            name = 'default'
            self.box_define.append([])
            self.box_define[len(self.box_define)-1]=[self._ptsx,self._ptsy,name]
        else:    
            mboxx,mboxy = self._map(self._ptsx,self._ptsy)
            self._drawbox.set_data(mboxx, mboxy)
            self._drawbox.axes.set_xscale('linear')
            self._drawbox.axes.set_yscale('linear')
            self._drawbox.axes.figure.canvas.draw()
            
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # what should happen if a mouse button is clicked when cursor is over axis 
    def _anybutton_(self,event):
        #global self
        if event.inaxes:
            x, y, b, k = event.xdata, event.ydata, event.button, []
            #print('button=%d, xdata=%f, ydata=%f' % ( event.button, event.xdata, event.ydata))
            #____left mouse button______________________________________________
            if event.button==1 : 
                xc, yc = x, y
                if self._press=='None':
                    self._update_center_(xc, yc)
                    self._move_center_()
                elif self._press=='Box':
                    self._Boxbuilder_(x,y)
                elif self._press=='Polygon':
                    self._Polygonbuilder_(x,y)
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
                #plt.close('all')
                plt.show(block=False)
                
                plt.close(self._figure)
                plt.close('all')
                
                return
                
            #___BOX MODE [l]___________________________________________________
            elif k=='b' or k=='1':
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
                    print('draw box: [OFF]')
                    self._text.set_text('draw box: [OFF]')
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
                    self.box_define[len(self.box_define)-1]=[list(np.sort(self._ptsx)),list(np.sort(self._ptsy)),name]
                    
            #___BOX MODE [l]___________________________________________________
            elif k=='p' or k=='2':
                #---------------------------------------------------------------
                # activate line mode, start line 
                if self._press == 'None':
                    print('draw polygon: [ON]')
                    self._text.set_text('draw polygon: [ON]')
                    plt.draw(),self._figure.canvas.draw()
                    
                    self._press = 'Polygon'
                    self._drawbox, = plt.plot([], [],color='w',linestyle='None',linewidth=2,marker='o',axes=self._ax)
                    self._Polygonbuilder_([],[])
                    #plt.xscale('linear'),
                    #plt.yscale('linear')
                    self._figure.canvas.draw()
                #---------------------------------------------------------------
                ## switch off line mode, end line
                elif self._press =='Polygon': 
                    print('draw polygon: [OFF]')
                    self._text.set_text('draw polygon: [OFF]')
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
                if   self._press=='Box':
                    self._Boxbuilder_(x,y,mode='remove')
                elif self._press=='Polygon':
                    self._Polygonbuilder_(x,y,mode='remove')
    
    
    #+_________________________________________________________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    # select which node points are within selected box/polygon
    def select_pointsinbox(self,mesh):
        from matplotlib import path
        self.box_idx = []
        for ii in range(0,len(self.box_define)):
            if len(self.box_define[ii][0])>2:
                p = path.Path(list(zip(self.box_define[ii][0],self.box_define[ii][1])))
            else:
                auxboxx = [ self.box_define[ii][0][0],\
                            self.box_define[ii][0][1],\
                            self.box_define[ii][0][1],\
                            self.box_define[ii][0][0]]
                auxboxy = [ self.box_define[ii][1][0],\
                            self.box_define[ii][1][0],\
                            self.box_define[ii][1][1],\
                            self.box_define[ii][1][1]]
                p = path.Path(list(zip(auxboxx,auxboxy)))
            
            self.box_idx.append([])
            self.box_idx[ii] = p.contains_points(np.array(list(zip(mesh.nodes_2d_xg,mesh.nodes_2d_yg))))
            
            # exlude node points whos bottom depth is shallower than botlimit
            if len(self.botlimit)>0:
                self.box_idx[ii][np.abs(mesh.nodes_2d_zg)<np.abs(self.botlimit)] = False
        
    #+_________________________________________________________________________+
    #|                                                                         |
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
        for ii in range(0,len(self.box_define)):
            self.value[ii] = box2.value[ii]-box.value[ii]
        self.anom  = True
        self.crange= []
        self.cmap  = 'blue2red'
        self.depth = box.depth
        self.sname = box.sname 
        self.lname = box.lname
        self.unit  = box.unit
    
    
    #+___PLOT FESOM2.0 DATA IN INDEX BOX OVER DEPTH AND TIME DATA______________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def plot_index_t_x_z(self,numb=[],figsize=[],do_subplot=[],do_output=True):
        from set_inputarray import inputarray
        fsize=16
        #_______________________________________________________________________
        if isinstance(numb,int)==True : numb = [numb]
        if len(numb)==0 : numb=range(0,len(self.box_define))
        
        #_______________________________________________________________________
        for ii in numb:
            #fig = plt.figure(figsize=(12, 6))
            #ax1  = plt.gca()
            
            if len(figsize)==0 : figsize=[10,5]
            #___________________________________________________________________________
            # plot is not part of subplot
            if len(do_subplot)==0:
                fig = plt.figure(figsize=figsize)
                ax1 = plt.gca()
            else:
                fig=do_subplot[0]
                ax1=do_subplot[1]
                fig.sca(ax1)
            resolution = 'c'
            fsize = 14
            
            #___________________________________________________________________
            cnumb= 30
            cmin = np.nanmin(self.value[ii])
            cmax = np.nanmax(self.value[ii])
            cref = cmin + (cmax-cmin)/2
            cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
            #cref =0.0
            
            #___________________________________________________________________
            # if anomaly data    
            if self.anom==True: 
                cmax,cmin,cref = np.nanmax(self.value[ii]), np.nanmin(self.value[ii]), 0.0
                self.cmap='blue2red'
            #___________________________________________________________________
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
            
            if do_output==True: print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
            if do_output==True: print('[cnum]=',cnumb)
            cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,self.cmap)
            if do_output==True: print('clevel = ',clevel)
            
            do_drawedges=True
            if clevel.size>50: do_drawedges=False
            
            # overwrite discrete colormap
            #cmap0 = cmocean.cm.balance
            #do_drawedges=False
            
            #___________________________________________________________________
            # make pcolor or contour plot 
            #depth = self.depth[:-1] + (self.depth[1:]-self.depth[:-1])/2.0
            #depth = self.zlev[:-1] + (self.zlev[1:]-self.zlev[:-1])/2.0
            depth = -self.zlev
            depthlim = np.sum(~np.isnan(self.value[ii][0,:])).max()
            if depthlim==depth.shape: depthlim=depthlim-1
            yy,xx = np.meshgrid(depth,self.time)
            
            #___________________________________________________________________
            data_plot = np.copy(self.value[ii])
            data_plot[data_plot<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
            data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
            
            #___________________________________________________________________
            if inputarray['which_plot']=='pcolor':
                hp=ax1.pcolormesh(xx[:,0:depthlim],yy[:,0:depthlim],data_plot[:,0:depthlim],
                            shading='flat',#flat
                            antialiased=False,
                            edgecolor='None',
                            cmap=cmap0),
                            #vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
            else: 
                hp=ax1.contourf(xx[:,0:depthlim],yy[:,0:depthlim],data_plot[:,0:depthlim],levels=clevel,
                            antialiased=True,
                            cmap=cmap0,
                            vmin=clevel[0], vmax=clevel[-1])
                ax1.contour(xx[:,0:depthlim],yy[:,0:depthlim],data_plot[:,0:depthlim],levels=clevel,
                            antialiased=True,
                            colors='k',
                            linewidths=[0.25,0.1],
                            linestyles=['solid'],
                            vmin=np.nanmin(data_plot), vmax=np.nanmax(data_plot))
            #hp.cmap.set_under([0.4,0.4,0.4])
            
            #___________________________________________________________________
            # set main axes
            ax1.set_xlim(self.time.min(),self.time.max())
            ax1.set_ylim(0,depth[depthlim-1])
            ax1.invert_yaxis()
            #ax1.set_axis_bgcolor([0.25,0.25,0.25])
            ax1.tick_params(axis='both',which='major',direction='out',length=8,labelsize=fsize)
            ax1.minorticks_on()
            ax1.tick_params(axis='both',which='minor',direction='out',length=4,labelsize=fsize)
            ax1.set_xlabel('Time [years]',fontdict=dict(fontsize=fsize))
            ax1.set_ylabel('Depth [km]',fontdict=dict(fontsize=fsize))
            plt.title(self.descript+' - '+self.box_define[0][2],fontdict= dict(fontsize=fsize*2),verticalalignment='bottom')
            
            #___________________________________________________________________
            # draw colorbar
            #divider = make_axes_locatable(ax1)
            #cax     = divider.append_axes("right", size="2.5%", pad=0.5)
            #plt.clim(clevel[0],clevel[-1])
            him = ax1.get_images()
            for im in ax1.get_images():
                im.set_clim(clevel[0],clevel[-1])
            
            #cbar = plt.colorbar(hp,ax=ax1,cax=cax,ticks=clevel,drawedges=do_drawedges)
            cbar = plt.colorbar(hp,ax=ax1,ticks=clevel,drawedges=do_drawedges)
            cbar.set_label(self.lname+' '+self.unit+'\n'+self.str_time, size=fsize)
            
            cl = plt.getp(cbar.ax, 'ymajorticklabels')
            plt.setp(cl, fontsize=fsize)
            
            # kickout some colormap labels if there are to many
            ncbar_l=len(cbar.ax.get_yticklabels()[:])
            idx_cref = np.where(clevel==cref)[0]
            idx_cref = np.asscalar(idx_cref)
            nmax_cbar_l = 10
            nstep = ncbar_l/nmax_cbar_l
            nstep = np.int(np.floor(nstep))
            #plt.setp(cbar.ax.get_yticklabels()[:], visible=False)
            #plt.setp(cbar.ax.get_yticklabels()[idx_cref::nstep], visible=True)
            #plt.setp(cbar.ax.get_yticklabels()[idx_cref::-nstep], visible=True)
            if nstep==0:nstep=1
            tickl = cbar.ax.get_yticklabels()
            idx = np.arange(0,len(tickl),1)
            idxb = np.ones((len(tickl),), dtype=bool)                
            idxb[idx_cref::nstep]  = False
            idxb[idx_cref::-nstep] = False
            idx = idx[idxb==True]
            for ii in list(idx):
                tickl[ii]=''
            cbar.ax.set_yticklabels(tickl)
            
            # reposition colorbar and axes
            plt.tight_layout()
            pos_ax   = ax1.get_position()
            ax1.set_position([pos_ax.x0,pos_ax.y0,0.90-pos_ax.y0,pos_ax.height])
            pos_ax   = ax1.get_position()
            pos_cbar = cbar.ax.get_position()
            cbar.ax.set_position([pos_ax.x0+pos_ax.width+0.01,pos_cbar.y0, pos_cbar.width, pos_cbar.height])
            fig.canvas.draw()
            #___________________________________________________________________
            # save figure
            if inputarray['save_fig']==True:
                print(' --> save figure: png')
                str_times= self.str_time.replace(' ','').replace(':','')
                str_deps = self.str_dep.replace(' ','').replace(',','').replace(':','')
                sfname = 'boxplot_'+self.box_define[ii][2]+'_'+self.descript+'_'+self.sname+'_'+str_times+'_'+str_deps+'.png'
                sdname = inputarray['save_figpath']
                if os.path.isdir(sdname)==False: os.mkdir(sdname)
                plt.savefig(sdname+sfname, \
                            format='png', dpi=600, \
                            bbox_inches='tight', pad_inches=0,\
                            transparent=True,frameon=True)
            #___________________________________________________________________
            plt.show(block=False)
            print('finish')
        
        return(fig,ax1,cbar)
        
        
    #+___PLOT FESOM2.0 DATA IN INDEX BOX POSITION______________________________+
    #|                                                                         |
    #+_________________________________________________________________________+
    def plot_index_position(self,mesh,numb=[]):
        from set_inputarray import inputarray
        fsize=16
        #_______________________________________________________________________
        if isinstance(numb,int)==True : numb = [numb]
        if len(numb)==0 : numb=range(0,len(self.box_define))
        
        #_______________________________________________________________________
        # draw position of box 
        for ii in numb:
        
            #___________________________________________________________________
            xmin,xmax = np.min(self.box_define[ii][0]), np.max(self.box_define[ii][0])
            ymin,ymax = np.min(self.box_define[ii][1]), np.max(self.box_define[ii][1])
            xmin,xmax,ymin,ymax = xmin-20.0, xmax+20.0, ymin-20.0, ymax+20.0
            xmin,xmax,ymin,ymax = np.max([xmin,-180.0]),np.min([xmax,180.0]),np.max([ymin,-90.0]),np.min([ymax,90.0])
            
            #___________________________________________________________________
            figp, ax = plt.figure(figsize=(13, 13)), plt.gca()
            map     = Basemap(projection = 'cyl',resolution = 'c',
                        llcrnrlon = xmin, urcrnrlon = xmax, llcrnrlat = ymin, urcrnrlat = ymax)
            
            mx,my     = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
            
            
            map.bluemarble()
            fesom_plot_lmask(map,mesh,ax,'0.6')
            
            xlabels,ylabels=[0,0,0,1],[1,0,0,0]
            xticks,yticks = np.arange(0.,360.,10.), np.arange(-90.,90.,5.)
            map.drawparallels(yticks,labels=ylabels,fontsize=fsize)
            map.drawmeridians(xticks,labels=xlabels,fontsize=fsize)
            map.drawmapboundary(linewidth=1.0)
            
            #___________________________________________________________________
            patch=[]
            if len(self.box_define[ii][0])>2:
                ax.plot(self.box_define[ii][0]    ,self.box_define[ii][1] ,linestyle='None'   ,color='w',linewidth=2.0,marker='o',mfc='w',mec='k',axes=ax)
                patch.append(Polygon(list(zip(self.box_define[ii][0],self.box_define[ii][1])),closed=True,clip_on=True) )
            else:
                auxboxx = [ self.box_define[ii][0][0],\
                            self.box_define[ii][0][1],\
                            self.box_define[ii][0][1],\
                            self.box_define[ii][0][0],\
                            self.box_define[ii][0][0]]
                auxboxy = [ self.box_define[ii][1][0],\
                            self.box_define[ii][1][0],\
                            self.box_define[ii][1][1],\
                            self.box_define[ii][1][1],\
                            self.box_define[ii][1][0],]
                ax.plot(auxboxx    ,auxboxy ,linestyle='None'   ,color='w',linewidth=2.0,marker='o',mfc='w',mec='k',axes=ax)
                patch.append(Polygon(list(zip(auxboxx,auxboxy)),closed=True,clip_on=True) )
            ax.add_collection(PatchCollection(patch, alpha=1.0,facecolor='none',edgecolor='w',zorder=1,linewidth=2.0,hatch='/'))    
            
            # plot selected mesh points
            ax.plot(mesh.nodes_2d_xg[self.box_idx[ii]],mesh.nodes_2d_yg[self.box_idx[ii]],color='r',linestyle='None',marker='.',alpha=0.5)
            plt.title(self.box_define[ii][2],fontdict= dict(fontsize=fsize*2),verticalalignment='bottom')
            
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
            plt.show()
            figp.canvas.draw()
