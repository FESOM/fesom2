# Patrick, Scholz 02.09.2018
import numpy as np
import time
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches as Polygon
from colormap_c2c    import *
import matplotlib.path as mpltPath
from matplotlib.tri import Triangulation
from numba import jit, njit, prange
import xarray

#+___CALCULATE MERIDIONAL OVERTURNING FROM VERTICAL VELOCITIES_________________+
#| Global MOC, Atlantik MOC, Indo-Pacific MOC, Indo MOC                        |
#|                                                                             |
#+_____________________________________________________________________________+
def calc_xmoc(mesh,data,dlat=1.0,do_onelem=True,do_output=True,which_moc='gmoc',\
              in_elemidx=[], out_elemidx=False, usemeshdiag=[]):
    #do_onelem=False
    #_________________________________________________________________________________________________
    t1=time.time()
    if do_output==True: print('_____calc. '+which_moc.upper()+' from vertical velocities via meridional bins_____')
        
    #___________________________________________________________________________
    # calculate/use index for basin domain limitation
    
    if which_moc=='gmoc':
        in_elemidx=[]
    else:    
        tt1=time.time()
        if len(in_elemidx)==0:
            #___________________________________________________________________
            # amoc2 ... calculate amoc without arctic
            if which_moc=='amoc2':
                # for calculation of amoc mesh focus must be on 0 degree longitude
                if mesh.focus!=0:
                    mesh.focus=0
                    mesh.fesom_grid_rot_r2g(str_mode='focus')
                box_moc = [-100.0,36.0,-30.0,79.0]
                in_elemidx=calc_basindomain(mesh,box_moc,do_output=do_output)
            #___________________________________________________________________
            # amoc ... calculate amoc including arctic
            if which_moc=='amoc':
                # for calculation of amoc mesh focus must be on 0 degree longitude
                if mesh.focus!=0:
                    mesh.focus=0
                    mesh.fesom_grid_rot_r2g(str_mode='focus')
                box_moc = [[-100.0,36.0,-30.0,79.0],[-180.0,180.0,65.0,90.0]]
                in_elemidx=calc_basindomain(mesh,box_moc,do_output=do_output)   
            #___________________________________________________________________
            # pmoc ... calculate indo-pacific moc
            elif which_moc=='pmoc':
                # for calculation of pmoc mesh focus must be on 180 degree longitude
                if mesh.focus!=180:
                    mesh.focus=180
                    mesh.fesom_grid_rot_r2g(str_mode='focus')
                box_moc = [[115.0,240.0,40.0,65.0],[30.0,295.0,-30.0,40.0]]
                in_elemidx=calc_basindomain(mesh,box_moc,do_output=do_output)    
            #___________________________________________________________________
            # imoc ... calculate indian ocean moc
            elif which_moc=='imoc':
                box_moc = [48.0,77.0,9.0,32.0]
                in_elemidx=calc_basindomain(mesh,box_moc,do_output=do_output)    
            #fig = plt.figure(figsize=[6,3])
            #plt.triplot(mesh.nodes_2d_xg,mesh.nodes_2d_yg,mesh.elem_2d_i[in_elemidx,:],linewidth=0.2)
            #plt.axis('scaled')
            #plt.title('Basin limited domain')
            #plt.show()
            #fig.canvas.draw()
            #STOP
            
    #___________________________________________________________________________
    # do moc calculation either on nodes or on elements. In moment local moc 
    # calculation is just supported on elements
    if do_onelem==True:
        
        if len(usemeshdiag)==0:
            mat_2d_iz = np.concatenate(( mesh.elem0_2d_iz,mesh.elem0_2d_iz[mesh.pbndtri_2d_i]))
            
        else:    
            ncfile  = Dataset(os.path.join(usemeshdiag))
            elem0_2d_iz=ncfile.variables['nlevels'][:]-1
            mat_2d_iz = np.concatenate((elem0_2d_iz,elem0_2d_iz[mesh.pbndtri_2d_i]))
            
        # calc triangle area if not alredy exist
        if len(mesh.elem0_2d_area)==0: mesh.fesom_calc_triarea()
        
        # augment elemental triangle area so it fits periodic augmented mesh
        mat_2d_area = np.concatenate((mesh.elem0_2d_area,mesh.elem0_2d_area[mesh.pbndtri_2d_i]))
        mat_2d_area = mat_2d_area*1e6
        
        # mean over elements
        mat_mean = data.value[mesh.elem_2d_i,:].sum(axis=1)/3.0*1e-6        
        
        if len(in_elemidx)>0:
            mat_mean = mat_mean[in_elemidx,:]
            # create meridional bins
            triy = mesh.nodes_2d_yg[mesh.elem_2d_i[in_elemidx,:]]
            lat   = np.arange(np.floor(triy.min())+dlat/2, np.ceil(triy.max()), dlat)
            lat_i = (( mesh.nodes_2d_yg[mesh.elem_2d_i[in_elemidx,:]].sum(axis=1)/3.0-lat[0])/dlat).astype('int')
            mat_2d_area = mat_2d_area[in_elemidx]
            mat_2d_iz   = mat_2d_iz[in_elemidx]
            
        else:
            # binning latitudes
            triy = mesh.nodes_2d_yg[mesh.elem_2d_i]
            lat   = np.arange(np.floor(triy.min())+dlat/2, np.ceil(triy.max()), dlat)
            lat_i = (( mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0-lat[0])/dlat).astype('int')
        del triy
        
        # calculate area weighted mean
        mat_mean = np.multiply(mat_mean, mat_2d_area[:,np.newaxis])
                
        del mat_2d_area
    else:
        mat_2d_iz = mesh.nodes_2d_iz
        
        # keep in mind that node area info is changing over depth--> therefor load from file 
        ncfile  = Dataset(os.path.join(data.path, data.runid+'.mesh.diag.nc'))
        mat_2d_area=ncfile.variables['nod_area'][:,:].transpose()
        mat_2d_area = np.concatenate((mat_2d_area,mat_2d_area[mesh.pbndn_2d_i,:]),axis=0)
        
        # create meridional bins
        lat   = np.arange(-90+dlat/2, 90-dlat/2, dlat)
        lat_i = (( mesh.nodes_2d_yg-lat[0])/dlat).astype('int')
        
        # mean over elements
        mat_mean = data.value*1e-6
        
        # calculate area weighted mean
        mat_mean = np.multiply(mat_mean, mat_2d_area)
        del mat_2d_area
        
    #___________________________________________________________________________
    # This approach is five time faster than the original from dima at least for
    # COREv2 mesh but needs probaply a bit more RAM
    moc     = np.zeros([mesh.nlev,lat.size])
    bottom  = np.zeros([lat.size,])
    numbtri = np.zeros([lat.size,])
    topo    = np.float16(mesh.zlev[mat_2d_iz])
    
    # this is more or less requird so bottom patch looks aceptable
    if which_moc=='pmoc':
        topo[np.where(topo>-30.0)[0]]=np.nan  
    else:
        #topo[np.where(topo>-100.0)[0]]=np.nan  
        topo[np.where(topo>-30.0)[0]]=np.nan  
    
    # be sure ocean floor is setted to zero
    for di in range(0,mesh.nlev):
        mat_idx = np.where(di>=mat_2d_iz)[0]
        mat_mean[mat_idx,di]=0.0
        
    # loop over meridional bins
    for bini in range(lat_i.min(), lat_i.max()+1):
        numbtri[bini]= np.sum(lat_i==bini)
        moc[:, bini]=mat_mean[lat_i==bini,:].sum(axis=0)
        #bottom[bini] = np.nanmedian(topo[lat_i==bini])
        bottom[bini] = np.nanpercentile(topo[lat_i==bini],15)
        
    # kickout outer bins where eventually no triangles are found
    idx    = numbtri>0
    moc    = moc[:,idx]
    bottom = bottom[idx]
    lat    = lat[idx]
    
    # do cumulative summation to finally calculate moc
    if len(in_elemidx)>0:
        moc = np.fliplr(moc)
        moc = -moc.cumsum(axis=1)
        moc = np.fliplr(moc)
    else:
        moc = moc.cumsum(axis=1)
    
    #___________________________________________________________________________
    # smooth bottom line a bit 
    filt=np.array([1,2,3,2,1])
    #filt=np.array([1,2,1])
    filt=filt/np.sum(filt)
    aux = np.concatenate( (np.ones((filt.size,))*bottom[0],bottom,np.ones((filt.size,))*bottom[-1] ) )
    aux = np.convolve(aux,filt,mode='same')
    bottom = aux[filt.size:-filt.size]
    
    #___________________________________________________________________________
    t2=time.time()
    if do_output==True: print(' --> total time:{:.3f} s'.format(t2-t1))
    
    #___________________________________________________________________________
    if which_moc=='pmoc':
        # rotate mesh back to starting point
        if mesh.focus!=0:
           mesh.focus=0
           mesh.fesom_grid_rot_r2g(str_mode='focus')
    
    #___________________________________________________________________________
    # variable number of output fields if you also want to write out the basin limited domain index
    #if out_elemidx==True and which_moc!='gmoc':    
    if out_elemidx==True:       
        return(moc,lat,bottom,in_elemidx)
    else:
        return(moc,lat,bottom)
    
 
#+___PLOT MERIDIONAL OVERTRUNING CIRCULATION  _________________________________+
#|                                                                             |
#+_____________________________________________________________________________+
def plot_xmoc(lat, depth, moc, bottom=[], which_moc='gmoc',
              str_descript='', str_time='', figsize=[], crange=[], cnumb=20,
              do_clabel=True, do_subplot=[]):    
    
    #___________________________________________________________________________
    if len(figsize)==0: figsize=[13,6]
    
    #___________________________________________________________________________
    # plot is not part of subplot
    if len(do_subplot)==0:
        fig = plt.figure(figsize=figsize)
        ax1 = plt.gca()
    else:
        fig = do_subplot[0]
        ax1 = do_subplot[1]
    
    resolution = 'c'
    fsize = 10
    txtx, txty = lat[0]+5,depth[-1]+100
    
    #+_________________________________________________________________________+
    #| set minimum, maximum and reference values for the creation of the       |
    #| adjustable colormap                                                     |
    #+_________________________________________________________________________+
    #cnumb = 20; # minimum number of colors
    #cmin,cmax,cref  = -6,16,0 # [cmin, cmax, cref]  --> MLD2
    cmin,cmax,cref  = moc[np.where(depth<=-500)[0][0]::,:].min(),moc[np.where(depth<=-500)[0][0]::,:].max(),0.0 # [cmin, cmax, cref]  --> MLD2
    if len(crange)!=0: cmin,cmax,cref = crange[0],crange[1],0.0
    print(cmin,cmax,cref)
    cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,'blue2red')
    cbot = [0.5,0.5,0.5]
    do_drawedges=True
    if clevel.size>30: do_drawedges=False

    #+_________________________________________________________________________+
    #| plot AXES1                                                              |
    #+_________________________________________________________________________+
    plt.sca(ax1)
    data_plot = moc
    data_plot[data_plot<clevel[ 0]]  = clevel[ 0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    hp1=plt.contourf(lat,depth,data_plot,levels=clevel,extend='both',cmap=cmap0)
    CS=plt.contour(lat,depth,data_plot,levels=clevel,colors='k',linewidths=[0.5,0.25],antialised=True)
    if do_clabel: ax1.clabel(CS,CS.levels[np.where(CS.levels!=cref)],inline=1,inline_spacing=1, fontsize=6,fmt='%1.1f Sv')
    if len(bottom)>0:
        ax1.plot(lat,bottom,color='k')
        ax1.fill_between(lat, bottom, depth[-1],color=cbot,zorder=2)#,alpha=0.95)

    #___________________________________________________________________________
    for im in ax1.get_images():
        im.set_clim(clevel[0],clevel[-1])
    #___________________________________________________________________________   
    ax1.text(txtx,txty,str_descript , fontsize=14, fontweight='bold',horizontalalignment='left')
    ax1.grid(color='k', linestyle='-', linewidth=0.25,alpha=1.0)
    ax1.set_xlabel('Latitudes [deg]',fontsize=12)
    ax1.set_ylabel('Depth [m]',fontsize=12)
    
    #+_________________________________________________________________________+
    #| plot first colorbar                                                     |
    #+_________________________________________________________________________+
    cbar1 = plt.colorbar(hp1,ax=ax1,ticks=clevel,drawedges=True,extend='neither',extendrect=False,extendfrac='auto')
    cbar1.set_label('Global Meridional Overturning Circulation [Sv]'+'\n'+str_time, size=fsize+2)
    if which_moc=='gmoc':
        cbar1.set_label('Global Meridional Overturning Circulation [Sv]'+'\n'+str_time, size=fsize+2)
    elif which_moc=='amoc' or which_moc=='amoc2':
        cbar1.set_label('Atlantic Meridional Overturning Circulation [Sv]'+'\n'+str_time, size=fsize+2)
    elif which_moc=='pmoc':
        cbar1.set_label('Indo-Pacific Meridional Overturning Circulation [Sv]'+'\n'+str_time, size=fsize+2)
    elif which_moc=='imoc':
        cbar1.set_label('Indo Meridional Overturning Circulation [Sv]'+'\n'+str_time, size=fsize+2)

    ncbar_l=len(cbar1.ax.get_yticklabels()[:])
    idx_cref = np.where(clevel==cref)[0]
    idx_cref = np.asscalar(idx_cref)
    nmax_cbar_l = 10
    nstep = ncbar_l/nmax_cbar_l
    nstep = np.int(np.floor(nstep))
    if nstep==0: nstep=1
    
    #plt.setp(cbar1.ax.get_yticklabels()[:], visible=False)
    #plt.setp(cbar1.ax.get_yticklabels()[idx_cref::nstep], visible=True)
    #plt.setp(cbar1.ax.get_yticklabels()[idx_cref::-nstep], visible=True)
    #plt.show(block=False)    
    fig.canvas.draw()
    tickl = cbar1.ax.get_yticklabels()
    idx = np.arange(0,len(tickl),1)
    idxb = np.ones((len(tickl),), dtype=bool)                
    idxb[idx_cref::nstep]  = False
    idxb[idx_cref::-nstep] = False
    idx = idx[idxb==True]
    for ii in list(idx):
        tickl[ii]=''
    cbar1.ax.set_yticklabels(tickl)
    
    fig.canvas.draw()
    return(fig,ax1)
    

    
#+___PLOT MERIDIONAL OVERTRUNING CIRCULATION TIME-SERIES_______________________+
#|                                                                             |
#+_____________________________________________________________________________+
def plot_xmoc_tseries(time,moc_t,which_lat=['max'],which_moc='amoc',str_descript='',str_time='',figsize=[]):    
    import matplotlib.patheffects as path_effects
    from matplotlib.ticker import AutoMinorLocator

    if len(figsize)==0: figsize=[13,4]
    fig,ax= plt.figure(figsize=figsize),plt.gca()
    
    for ii in range(len(which_lat)):
        if which_lat[ii]=='max':
            str_label='max {:s}: 30°N<=lat<=45°N'.format(which_moc.upper(),which_lat[ii])
        else:
            str_label='max {:s} at: {:2.1f}°N'.format(which_moc.upper(),which_lat[ii])
        hp=ax.plot(time,moc_t[:,ii],\
                   linewidth=2,label=str_label,marker='o',markerfacecolor='w',\
                   path_effects=[path_effects.SimpleLineShadow(offset=(1.5,-1.5),alpha=0.3),path_effects.Normal()])
        # plot mean value with trinagle 
        plt.plot(time[0]-(time[-1]-time[0])*0.015,moc_t[:,ii].mean(),\
                 marker='<',markersize=8,markeredgecolor='k',markeredgewidth=0.5,\
                 color=hp[0].get_color(),zorder=3,clip_box=False,clip_on=False)
        
        # plot std. range
        plt.plot(time[0]-(time[-1]-time[0])*0.015,moc_t[:,ii].mean()+moc_t[:,ii].std(),\
                 marker='^',markersize=6,markeredgecolor='k',markeredgewidth=0.5,\
                 color=hp[0].get_color(),zorder=3,clip_box=False,clip_on=False)
        
        plt.plot(time[0]-(time[-1]-time[0])*0.015,moc_t[:,ii].mean()-moc_t[:,ii].std(),\
                 marker='v',markersize=6,markeredgecolor='k',markeredgewidth=0.5,\
                 color=hp[0].get_color(),zorder=3,clip_box=False,clip_on=False)
        
    ax.legend(loc='lower right', shadow=True,fancybox=True,frameon=True,mode='None')
    ax.set_xlabel('Time [years]',fontsize=12)
    ax.set_ylabel('{:s} in [Sv]'.format(which_moc.upper()),fontsize=12)
    minor_locator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.xaxis.set_minor_locator(minor_locator)
    plt.grid(which='major')
    plt.xticks(np.arange(1940,2015,5))
    plt.xlim(time[0]-(time[-1]-time[0])*0.015,time[-1]+(time[-1]-time[0])*0.015)    
    plt.show(block=False)
    
    fig.canvas.draw()
    return(fig,ax)


#+___CALCULATE BASIN LIMITED DOMAIN____________________________________________+
#| to calculate the regional moc (amoc,pmoc,imoc) the domain needs be limited  |
#| to corresponding basin. here the elemental index of the triangels in the    |
#| closed basin is calcualted                                                  |
#+_____________________________________________________________________________+
def calc_basindomain(mesh,box_moc,do_output=False):
    
    if do_output==True: print(' --> calculate basin limited domain',end='')
    t1=time.time()
    
    #___________________________________________________________________________
    # 1st. pre-limit ocean domain by pre defined boxes for atlantic, 
    # indo-pacific ... basin
    box_moc = np.matrix(box_moc)
    allbox_idx  =np.zeros((mesh.n2dea,),dtype=bool)
    for bi in range(0,box_moc.shape[0]):
        #_______________________________________________________________________
        box_idx = mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3.0<box_moc[bi,0]
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3.0>box_moc[bi,1])
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0<box_moc[bi,2])
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0>box_moc[bi,3])
        allbox_idx[np.where(box_idx==False)[0]]=True
    
    # kick out all the elements that are not located within predefined box
    allbox_elem = mesh.elem_2d_i[allbox_idx==True,:]
    
    #fig = plt.figure(figsize=[10,5])
    #plt.triplot(mesh.nodes_2d_xg,mesh.nodes_2d_yg,mesh.elem_2d_i[allbox_idx==True,:],linewidth=0.2)
    #plt.axis('scaled')
    #plt.title('Box limited domain')
    #plt.show()
    #fig.canvas.draw()
    
    #___________________________________________________________________________
    # 2nd. select main ocean basin cluster for either atlantic or pacific by finding 
    # the outer coastline of the main basin
    edge    = np.concatenate((allbox_elem[:,[0,1]], allbox_elem[:,[0,2]], allbox_elem[:,[1,2]]),axis=0)
    edge    = np.sort(edge,axis=1) 
    
    # python  sortrows algorythm --> matlab equivalent
    # twice as fast as list sorting
    sortidx = np.lexsort((edge[:,0],edge[:,1]))
    edge    = edge[sortidx,:].squeeze()
    edge    = np.array(edge)
    
    # list sorting
    #edge    = edge.tolist()
    #edge.sort()
    #edge    = np.array(edge)
        
    idx     = np.diff(edge,axis=0)==0
    idx     = np.all(idx,axis=1)
    idx     = np.logical_or(np.concatenate((idx,np.array([False]))),np.concatenate((np.array([False]),idx)))
        
    # all outer ocean edges that belong to preselected boxes defined domain
    bnde    = edge[idx==False,:]
    nbnde    = bnde.shape[0];
    
    # find most northern ocean edge and start from there
    maxi = np.argmax(mesh.nodes_2d_yg[bnde].sum(axis=1)/2)
    
    #fig = plt.figure(figsize=[10,5])
    #plt.plot(mesh.nodes_2d_xg[np.unique(bnde.flatten())],mesh.nodes_2d_yg[np.unique(bnde.flatten())],'*')
    #plt.plot(mesh.nodes_2d_xg[np.unique(bnde[maxi].flatten())],mesh.nodes_2d_yg[np.unique(bnde[maxi].flatten())],'*',color='red')
    #plt.axis('scaled')
    #plt.title('Outer boundary edges')
    #plt.show()
    #fig.canvas.draw()
    
    #___________________________________________________________________________
    # start with on outer coastline edge and find the next edge that is 
    # connected and so forth, like that build entire outer coastline of the main 
    # basin
    
    run_cont        = np.zeros((1,nbnde+1))*np.nan
    run_cont[0,:2]  = bnde[maxi,:] # initialise the first landmask edge
    #run_bnde        = bnde[1:,:] # remaining edges that still need to be distributed
    run_bnde        = np.delete(bnde, (maxi), axis=0)
    count_init      = 1;
    init_ind        = run_cont[0,0];
    ind_lc_s        = 0;
    
    ocebasin_polyg = []
    for ii in range(0,nbnde):
        #_______________________________________________________________________
        # search for next edge that contains the last node index from 
        # run_cont
        kk_rc = np.column_stack(np.where( run_bnde==np.int(run_cont[0,count_init]) ))
        kk_r  = kk_rc[:,0]
        kk_c  = kk_rc[:,1]
        count_init  = count_init+1
        if len(kk_c)==0 : break        
        
        #_______________________________________________________________________
        if kk_c[0] == 0 :
            run_cont[0,count_init] = run_bnde[kk_r[0],1]
        else:
            run_cont[0,count_init] = run_bnde[kk_r[0],0]
            
        #_______________________________________________________________________
        # if a land sea mask polygon is closed
        if  np.any(run_bnde[kk_r[0],:] == init_ind):
            count_init  = count_init+1
            
            aux_lx = mesh.nodes_2d_xg[np.int64(run_cont[0,0:count_init])];
            aux_ly = mesh.nodes_2d_yg[np.int64(run_cont[0,0:count_init])];
            aux_xy = np.zeros((count_init,2))
            aux_xy[:,0] = aux_lx
            aux_xy[:,1] = aux_ly
            ocebasin_polyg=aux_xy
            del aux_lx; del aux_ly; del aux_xy
            
            ind_lc_s = ind_lc_s+count_init+1;
            
            count_init = count_init+1
            aux_ind  = np.arange(0,run_bnde.shape[0],1)
            run_bnde = run_bnde[aux_ind!=kk_r[0],:]
            if np.size(run_bnde)==0:
                break
                
            #___________________________________________________________________
            run_cont        = np.zeros((1,nbnde))*np.nan
            run_cont[0,:2]  = run_bnde[0,:]
            run_bnde        = run_bnde[1:,:]
            count_init=1;
        else:
            aux_ind =np.arange(0,run_bnde.shape[0],1)
            run_bnde=run_bnde[aux_ind!=kk_r[0],:]
            
    #___________________________________________________________________________
    # check which preselected triangle centroids are within main ocean basin 
    # polygon 
    ptsc = list(zip(mesh.nodes_2d_xg[allbox_elem].sum(axis=1)/3,mesh.nodes_2d_yg[allbox_elem].sum(axis=1)/3))
    
    #___________________________________________________________________________
    # Option (1)
    # python Matplotlib mplPath seems to be faster at least by a factor of 2 when
    # compared to ray_tracing method
    #print(' >> use mpltPath ',end='')
    #path = mpltPath.Path(ocebasin_polyg)
    #inside_ocebasin = path.contains_points(ptsc)
    
    # Option (2)
    # determine points in polygon by ray tracing method
    #print(' >> use rtracing ',end='')
    #inside_ocebasin = [calc_ray_tracing(point[0], point[1], np.array(ocebasin_polyg)) for point in ptsc]
    #inside_ocebasin = np.array(inside_ocebasin)
    
    # Option (3)
    # determine points in polygon by parallel (numba optimized) ray tracing 
    # method --> considerable faster for large meshes
    print(' >> use rtracing parallel ',end='')
    inside_ocebasin = calc_ray_tracing_parallel(np.array(ptsc),np.array(ocebasin_polyg),np.zeros((len(ptsc),),dtype=bool))
    
    #___________________________________________________________________________
    # write out regional indices with respect to the global elemental array
    allbox_tidx = np.where(allbox_idx==True)[0]
    allbox_fin = allbox_tidx[inside_ocebasin==True]
    
    #fig = plt.figure(figsize=[10,5])
    #plt.triplot(mesh.nodes_2d_xg,mesh.nodes_2d_yg,mesh.elem_2d_i[allbox_fin,:],linewidth=0.2)
    #plt.axis('scaled')
    #plt.title('Basin limited domain')
    #plt.show()
    #fig.canvas.draw()
    
    #___________________________________________________________________________
    t2=time.time()
    print(" >> time: {:.3f} s".format(t2-t1))   
    
    return(allbox_fin)


#+___CALCULATE BASIN LIMITED DOMAIN___________________________________________________________________+
#| to calculate the regional moc (amoc,pmoc,imoc) the domain needs be limited to corresponding basin.
#| here the elemental index of the triangels in the closed basin is calcualted
#+____________________________________________________________________________________________________+
def calc_basindomain_slow(mesh,box_moc,do_output=False):
    
    if do_output==True: print('     --> calculate regional basin limited domain')
    box_moc = np.matrix(box_moc)
    for bi in range(0,box_moc.shape[0]):
        #_____________________________________________________________________________________________
        box_idx = mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3.0<box_moc[bi,0]
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3.0>box_moc[bi,1])
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0<box_moc[bi,2])
        box_idx = np.logical_or(box_idx,mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3.0>box_moc[bi,3])
        box_idx = np.where(box_idx==False)[0]
        box_elem2di = mesh.elem_2d_i[box_idx,:]

        #_____________________________________________________________________________________________
        # calculate edge indices of box limited domain
        edge_12     = np.sort(np.array(box_elem2di[:,[0,1]]),axis=1)
        edge_23     = np.sort(np.array(box_elem2di[:,[1,2]]),axis=1)
        edge_31     = np.sort(np.array(box_elem2di[:,[2,0]]),axis=1)
        edge_triidx = np.arange(0,box_elem2di.shape[0],1)

        #_____________________________________________________________________________________________
        # start with seed triangle
        seed_pts     = [box_moc[bi,0]+(box_moc[bi,1]-box_moc[bi,0])/2.0,box_moc[bi,2]+(box_moc[bi,3]-box_moc[bi,2])/2.0]
        seed_triidx  = np.argsort((mesh.nodes_2d_xg[box_elem2di].sum(axis=1)/3.0-seed_pts[0])**2 + (mesh.nodes_2d_yg[box_elem2di].sum(axis=1)/3.0-seed_pts[1])**2,axis=-0)[0]
        seed_elem2di = box_elem2di[seed_triidx,:]
        seed_edge    = np.concatenate((seed_elem2di[:,[0,1]], seed_elem2di[:,[1,2]], seed_elem2di[:,[2,0]]),axis=0)     
        seed_edge    = np.sort(seed_edge,axis=1) 
        
        # already delete seed triangle and coresbonding edges from box limited domain list
        edge_triidx = np.delete(edge_triidx,seed_triidx)
        edge_12     = np.delete(edge_12,seed_triidx,0)
        edge_23     = np.delete(edge_23,seed_triidx,0)
        edge_31     = np.delete(edge_31,seed_triidx,0)

        #_____________________________________________________________________________________________
        # do iterative search of which triangles are connected to each other and form cluster
        t1 = time.time()
        tri_merge_idx = np.zeros((box_elem2di.shape[0],),dtype='int')
        tri_merge_count = 0
        for ii in range(0,10000): 
            #print(ii,tri_merge_count,seed_edge.shape[0])
        
            # determine which triangles contribute to edge
            triidx12 = ismember_rows(seed_edge,edge_12)
            triidx23 = ismember_rows(seed_edge,edge_23)
            triidx31 = ismember_rows(seed_edge,edge_31)
        
            # calculate new seed edges
            seed_edge = np.concatenate((edge_23[triidx12,:],edge_31[triidx12,:],\
                                        edge_12[triidx23,:],edge_31[triidx23,:],\
                                        edge_12[triidx31,:],edge_23[triidx31,:]))
            
            # collect all found connected triagles    
            triidx = np.concatenate((triidx12,triidx23,triidx31))
            triidx = np.unique(triidx)
            
            # break out of iteration loop 
            if triidx.size==0: break 
                
            # add found trinagles to final domain list    
            tri_merge_idx[tri_merge_count:tri_merge_count+triidx.size]=edge_triidx[triidx]
            tri_merge_count = tri_merge_count+triidx.size
            
            # delete already found trinagles and edges from list
            edge_triidx = np.delete(edge_triidx,triidx)
            edge_12     = np.delete(edge_12,triidx,0)
            edge_23     = np.delete(edge_23,triidx,0)
            edge_31     = np.delete(edge_31,triidx,0)
    
            del triidx,triidx12,triidx23,triidx31
        
        tri_merge_idx = tri_merge_idx[:tri_merge_count-1]
        t2=time.time()
        if do_output==True: print('         elpased time:'+str(t2-t1)+'s')
        
        #_____________________________________________________________________________________________
        # calculate final domain limited trinagle cluster element index
        if bi==0:
            box_idx_fin = box_idx[tri_merge_idx]
        else:
            box_idx_fin = np.concatenate((box_idx_fin,box_idx[tri_merge_idx]))
        
    return(box_idx_fin)


#+___EQUIVALENT OF MATLAB ISMEMBER FUNCTION___________________________________________________________+
#|                                                                                                    |
#+____________________________________________________________________________________________________+
def ismember_rows(a, b):
    return np.flatnonzero(np.in1d(b[:,0], a[:,0]) & np.in1d(b[:,1], a[:,1]))


#+___RAY TRACING METHOD TO CHECK IF POINT IS IN POLYGON________________________+
#| see...https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-
#| checking-if-a-point-is-inside-a-polygon-in-python
#| 
#+_____________________________________________________________________________+
@jit(nopython=True)
def calc_ray_tracing(x,y,poly):
    n = len(poly)
    inside = False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside
#+___RAY TRACING METHOD PARALLEL TO CHECK IF POINT IS IN POLYGON_______________+
#| see...https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-|
#| checking-if-a-point-is-inside-a-polygon-in-python                           |
#|                                                                             |
#+_____________________________________________________________________________+
#@njit(parallel=True,nopython=True)
@njit(parallel=True)
def calc_ray_tracing_parallel(pts,poly,inside):
    n = len(poly)
    npts=len(pts)
    for j in prange(npts):
        x,y = pts[j]
        p1x,p1y = poly[0]
        for i in range(n+1):
            p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xints:
                            inside[j] = not inside[j]
            p1x,p1y = p2x,p2y

    return inside
