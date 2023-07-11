import matplotlib.pyplot as plt
import copy as cp
from sub_climatology  import *
from colormap_c2c    import *
from sub_regriding_adapt import *
from sub_fesom_plot  import *
from sub_fesom_data_netcdf4  import * 

#_______________________________________________________________________________
# compute bottom density
def compute_bottom_density(mesh,data,which_sigma4):
    aux1 = np.zeros([mesh.n2dna,])*np.nan
    for node in range(0,mesh.n2dna):
        #if data.value[node,mesh.nodes_2d_izg[node]-1]<which_sigma4:continue
        #if mesh.nodes_2d_yg[node]>60: continue    
        ## exclude mediteranian sea    
        #if mesh.nodes_2d_xg[node]>=-5.5 and mesh.nodes_2d_xg[node]<=35.5 and \
           #mesh.nodes_2d_yg[node]>=30.0 and mesh.nodes_2d_xg[node]<=46.0: continue
        aux1[node] = data.value[node,mesh.nodes_2d_izg[node]-1]
    data.value = aux1    
    data.str_dep = ', bottom'
    return(data)

def compute_bottom(mesh,data,nlayer=1):
    aux1 = np.zeros([mesh.n2dna,])*np.nan
    for node in range(0,mesh.n2dna): 
        if nlayer>1:
            aux1[node] = data.value[node,mesh.nodes_2d_izg[node]-1-nlayer+1:mesh.nodes_2d_izg[node]-1+1].mean(axis=0)
        else:    
            aux1[node] = data.value[node,mesh.nodes_2d_izg[node]-1]
    data.value = aux1    
    data.str_dep = ', bottom'
    return(data)

#_______________________________________________________________________________
# compute climatological anomaly 
def clim_data_anom(mesh,clim,data,clim_descript):    
    # calculate anomaly climatology to fesom2.0
    data11            = cp.deepcopy(clim)
    
    data_clim1        = clim_vinterp(data11,data.depth)
    mlon,mlat         = np.meshgrid(data11.lon, data11.lat)
    
    data11.str_dep    = data.str_dep
    distances1, inds1 = create_indexes_and_distances(mesh, mlon, mlat,k=5, n_jobs=2)
    data_fesom1       = fesom2regular(data.value, mesh, mlon, mlat, distances=distances1, inds=inds1, radius_of_influence=100000)
    data11.anom       = data_fesom1-data_clim1
    data11.descript   = data.descript+'-'+clim_descript
    data11.value      =[]
    return(data11)

#_______________________________________________________________________________
# setup colorbar 
def setup_colorbar(hp, ax, clevel, cref, str_cbar, ncl=6, orientation='horizontal',fsize=12,rotation=0):
    #___________________________________________________________________________
    cbar = plt.colorbar(hp,ax=ax,orientation=orientation,\
                         ticks=clevel,drawedges=True,extend='neither',\
                         extendrect=True,extendfrac='auto')
    cbar.set_label(str_cbar, fontsize=fsize)
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
    for ii in list(idx):
        tickl[ii]=''
    if cbar.orientation=='vertical':cbar.ax.set_yticklabels(tickl,fontsize=fsize,rotation=rotation)
    else:                           cbar.ax.set_xticklabels(tickl,fontsize=fsize,rotation=rotation)    
    
    return(cbar)

#_______________________________________________________________________________
# plot single moc panel
def plot_map_xmoc(fig, ax, mesh, lat, data, bottom, clevel, cmap, str_tl=[], \
                  str_yl=[], str_xl=[], str_txt=[], cbot= [0.5,0.5,0.5],fsize=12,str_enum=[]):
    fig.sca(ax)
    #___________________________________________________________________________
    data_plot = np.copy(data)
    data_plot[data_plot<clevel[ 0]]  = clevel[ 0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    hp=plt.contourf(lat,mesh.zlev,data_plot,levels=clevel,antialiased=False,extend='both',cmap=cmap)
    plt.contour(lat,mesh.zlev,data_plot,levels=clevel,colors='k',linewidths=[0.25],antialised=True)
    ax.plot(lat,bottom,color='k',linewidth=1)
    ax.fill_between(lat, bottom, mesh.zlev[-1],color=cbot,zorder=2)#,alpha=0.95)
    ax.set_aspect(np.abs(lat[-1]-lat[0])/np.abs(mesh.zlev[-1]-mesh.zlev[0])/2)
    for im in ax.get_images(): im.set_clim(clevel[0],clevel[-1])
    #___________________________________________________________________________    
    if len(str_tl)!=0 : ax.set_title(str_tl , fontsize=fsize)
    if len(str_yl)!=0 : 
        ax.set_ylabel(str_yl,fontsize=fsize)
    else:
        ax.set_yticklabels([])
    if len(str_xl)!=0 : ax.set_xlabel(str_xl,fontsize=fsize)   
    else:
        ax.set_xticklabels([])    
    if len(str_txt)!=0:
        txtx, txty = lat[0]+5, mesh.zlev[-1]+100
        ax.text(txtx,txty,str_txt , fontsize=fsize, fontweight='bold',horizontalalignment='left')
    if len(str_enum)!=0:    
        txtx, txty = lat[0]+2, bottom.min()+50
        ax.text(txtx,txty,str_enum , fontsize=fsize, fontweight='normal',horizontalalignment='left',verticalalignment='top')
        
    ax.axes.figure.canvas.draw()    
    return(ax,hp)

#_______________________________________________________________________________
# plot single panel of northern/southern stereographic projection 
def plot_map_stere(fig, ax, mesh, data, inputarray, setbndlat, clevel, cmap, \
                   nxyl=[6,6], pos_txt=[45,75], fsize=12, nargout=2, 
                   str_tl=[], str_yl=[], str_xl=[], str_txt=[], str_enum=[]):
    fig.sca(ax)
    resolution = 'c'
    #___________________________________________________________________________
    if setbndlat>0:
        map = Basemap(projection = 'npstere', resolution = 'c',lon_0 = 0, boundinglat = setbndlat,round= True)
    else:
        map = Basemap(projection = 'spstere', resolution = 'c',lon_0 = 180, boundinglat = setbndlat,round= True)
        
    #___________________________________________________________________________
    mx,my = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    tri      = Triangulation(mx, my,mesh.elem_2d_i)
    mask_tri = TriAnalyzer(tri).get_flat_tri_mask(0.00001)
    data_plot = np.copy(data.value)
    if np.any(np.isnan(data.value)):
        if data.value.size ==mesh.n2dna:
            mask_tri = np.logical_or(mask_tri,np.isnan(np.sum(data.value[mesh.elem_2d_i],axis=1)))
        elif data.value.size ==mesh.n2dea:    
            mask_tri = np.logical_or(mask_tri,np.isnan(data.value))
    
    if setbndlat>0:
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3<setbndlat)
    else:    
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3>setbndlat)
    tri.set_mask(mask_tri)    
    data_plot[data_plot<clevel[ 0]]  = clevel[ 0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    hp=ax.tricontourf(tri,data_plot,levels=clevel,antialiased=False,extend='both',cmap=cmap)

    #___________________________________________________________________________
    step=np.array([1.0,2.0,5.0,10.0,15.0,30.0,45.0,60.0])
    idx = np.array(np.where(step>=(inputarray['which_box'][1]-inputarray['which_box'][0])/nxyl[0])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepx=step[idx[0]]
    idx = np.array(np.where(step>=(inputarray['which_box'][3]-inputarray['which_box'][2])/nxyl[1])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepy=step[idx[0]]
    xticks , yticks =np.arange(-180.,180.+1,stepx), np.arange(-90.,90.+1,stepy)
    map.drawmeridians(xticks,labels=[],linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2,latmax=85)
    map.drawparallels(yticks,labels=[],linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2)
    
    #___________________________________________________________________________    
    map.drawmapboundary(fill_color='0.9',linewidth=1.0)
    fesom_plot_lmask(map,mesh,ax,'0.6')
    
    #___________________________________________________________________________
    for im in ax.get_images():     im.set_clim(clevel[0],clevel[-1])
    if len(str_tl)!=0 : ax.set_title(str_tl , fontsize=fsize)
    if len(str_yl)!=0 : ax.set_ylabel(str_yl, fontsize=fsize)
    if len(str_xl)!=0 : ax.set_xlabel(str_xl, fontsize=fsize)    
    if len(str_txt)!=0:
        txtx, txty = pos_txt[0],pos_txt[1]
        ax.text(txtx,txty,str_txt , fontsize=fsize, fontweight='bold',horizontalalignment='left')
    if len(str_enum)!=0:
        ax.text(0.025,0.1,str_enum,fontweight='normal',horizontalalignment='left',verticalalignment='top',fontsize=14,transform=ax.transAxes)   
    #___________________________________________________________________________
    ax.axes.figure.canvas.draw()
    if   nargout==2: return(ax,hp)
    elif nargout==3: return(ax,hp,tri)
    elif nargout==4: return(ax,hp,tri,map)
    else           : print(' --> nargout number is not supported')
    
#_______________________________________________________________________________
# plot single panel of cylindric proection on triangular fesom mesh
def plot_map_fesom(fig, ax, mesh, data, inputarray, clevel, cmap, xl=[0,0,0,1], \
                   yl=[1,0,0,0], nxyl=[6,6], fsize=12, which_lines=[0,0,0],
                   str_tl=[], str_yl=[], str_xl=[], str_txt=[], str_enum=[], nargout=2):
    fig.sca(ax)
    resolution = 'c'
    #___________________________________________________________________________
    map = Basemap(projection=data.proj,resolution='c', \
                  llcrnrlon=inputarray['which_box'][0], urcrnrlon=inputarray['which_box'][1], \
                  llcrnrlat=inputarray['which_box'][2], urcrnrlat=inputarray['which_box'][3],\
                  area_thresh=10000,ax=ax,suppress_ticks=True)
    
    #___________________________________________________________________________
    mx,my = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    tri      = Triangulation(mx, my,mesh.elem_2d_i)
    mask_tri = TriAnalyzer(tri).get_flat_tri_mask(0.00001)
    data_plot = np.copy(data.value)
    if np.any(np.isnan(data_plot)):
        mask_tri = np.logical_or(mask_tri,np.isnan(np.sum(data_plot[mesh.elem_2d_i],axis=1)))
        tri.set_mask(mask_tri)
    data_plot[data_plot<clevel[ 0]] = clevel[ 0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    hp=ax.tricontourf(tri,data_plot,levels=clevel,antialiased=False,extend='both',cmap=cmap)
    if which_lines[0]==1: ax.tricontour(tri,data_plot,levels=clevel[clevel<data.crange[2]],antialiased=True,colors='k', linewidths=0.2, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    if which_lines[1]==1: ax.tricontour(tri,data_plot,levels=clevel[clevel==data.crange[2]],antialiased=True,colors='k', linewidths=0.5, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    if which_lines[2]==1: ax.tricontour(tri,data_plot,levels=clevel[clevel>data.crange[2]],antialiased=True,colors='k', linewidths=0.2, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    #___________________________________________________________________________
    step=np.array([1.0,2.0,5.0,10.0,15.0,30.0,45.0,60.0])
    idx = np.array(np.where(step>=(inputarray['which_box'][1]-inputarray['which_box'][0])/nxyl[0])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepx=step[idx[0]]
    idx = np.array(np.where(step>=(inputarray['which_box'][3]-inputarray['which_box'][2])/nxyl[1])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepy=step[idx[0]]
    xticks , yticks =np.arange(-180.+stepx,180.,stepx), np.arange(-90.+stepy,90.-stepy+1,stepy)
    map.drawmeridians(xticks,labels=xl,linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2)
    map.drawparallels(yticks,labels=yl,linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2)
    
    #___________________________________________________________________________
    map.drawmapboundary(fill_color='0.8',linewidth=1.0)
    fesom_plot_lmask(map,mesh,ax,'0.6')
    
    #___________________________________________________________________________
    for im in ax.get_images():     im.set_clim(clevel1[0],clevel1[-1])
    if len(str_tl )!=0: ax.set_title(str_tl , fontsize=fsize)
    if len(str_yl )!=0: ax.set_ylabel(str_yl, fontsize=fsize) 
    if len(str_xl )!=0: ax.set_xlabel(str_xl, fontsize=fsize) 
    if len(str_txt)!=0:
        fac = 2.5
        txtx, txty = inputarray['which_box'][0]+fac,inputarray['which_box'][2]+fac
        ax.text(txtx,txty,str_txt , fontsize=fsize, fontweight='normal',horizontalalignment='left')
    if len(str_enum)!=0:    
        fac = 2.5
        txtx, txty = inputarray['which_box'][0]+fac,inputarray['which_box'][3]-fac
        ax.text(txtx,txty,str_enum, fontweight='normal',horizontalalignment='left',\
                verticalalignment='top',fontsize=fsize)
    #___________________________________________________________________________
    ax.axes.figure.canvas.draw()
    if   nargout==2: return(ax,hp)
    elif nargout==3: return(ax,hp,tri)
    elif nargout==4: return(ax,hp,tri,map)
    else           : print(' --> nargout number is not supported')
    
#_______________________________________________________________________________
# plot single panel of cylindric proection of climatoligical anomaly on regular grid
def plot_map_clima(fig, ax, mesh, data, inputarray, clevel, cmap, xl=[0,0,0,1], \
                   yl=[1,0,0,0], nxyl=[6,6], str_tl=[], str_yl=[], str_xl=[],str_txt=[], str_enum=[], fsize=12, which_lines=[0,0,0]):
    fig.sca(ax)
    resolution = 'c'
    #___________________________________________________________________________
    map = Basemap(projection=data.proj,resolution=resolution, \
                  llcrnrlon=inputarray['which_box'][0], urcrnrlon=inputarray['which_box'][1], \
                  llcrnrlat=inputarray['which_box'][2], urcrnrlat=inputarray['which_box'][3],\
                  area_thresh=10000,ax=ax,suppress_ticks=True)
    
    #___________________________________________________________________________
    if len(data.value)!=0 : data_plot = np.copy(data.value)
    elif len(data.anom)!=0: data_plot = np.copy(data.anom)
    data_plot = np.ma.masked_where(np.isnan(data_plot), data_plot)
    data_plot[data_plot<clevel[ 0]] = clevel[ 0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    hp=ax.contourf(data.lon,data.lat,data_plot,levels=clevel,antialiased=False,extend='both',cmap=cmap)
    if which_lines[0]==1: ax.contour(data.lon,data.lat,data_plot,levels=clevel[clevel<data.crange[2]],antialiased=True,colors='k', linewidths=0.2, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    if which_lines[1]==1: ax.contour(data.lon,data.lat,data_plot,levels=clevel[clevel==data.crange[2]],antialiased=True,colors='k', linewidths=0.5, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    if which_lines[2]==1: ax.contour(data.lon,data.lat,data_plot,levels=clevel[clevel>data.crange[2]],antialiased=True,colors='k', linewidths=0.2, linestyle='-',vmin=clevel[0], vmax=clevel[-1], alpha=0.25)
    #
    #___________________________________________________________________________
    step=np.array([5.0,10.0,15.0,30.0,45.0,60.0])
    idx = np.array(np.where(step>=(inputarray['which_box'][1]-inputarray['which_box'][0])/nxyl[0])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepx=step[idx[0]]
    idx = np.array(np.where(step>=(inputarray['which_box'][3]-inputarray['which_box'][2])/nxyl[1])[0])
    if idx.size==0 : idx = np.array([step.size-1])
    stepy=step[idx[0]]
    xticks , yticks =np.arange(-180.+stepx,180.,stepx), np.arange(-90.+stepy,90.-stepy+1,stepy)
    map.drawmeridians(xticks,labels=xl,linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2)
    map.drawparallels(yticks,labels=yl,linewidth=0.25,dashes=[1,1e-10],fontsize=fsize-2)
    
    #___________________________________________________________________________
    map.drawmapboundary(fill_color='0.9',linewidth=1.0)
    fesom_plot_lmask(map,mesh,ax,'0.6')
    
    #___________________________________________________________________________
    for im in ax.get_images(): im.set_clim(clevel[0],clevel[-1])
    if len(str_tl )!=0: ax.set_title(str_tl , fontsize=fsize)
    if len(str_yl )!=0: ax.set_ylabel(str_yl, fontsize=fsize)
    if len(str_xl )!=0: ax.set_xlabel(str_xl, fontsize=fsize) 
    if len(str_txt)!=0:
        fac = 2.5
        txtx, txty = inputarray['which_box'][0]+fac,inputarray['which_box'][2]+fac
        ax.text(txtx,txty,str_txt , fontsize=fsize, fontweight='normal',horizontalalignment='left')
    if len(str_enum)!=0:    
        fac = 2.5
        txtx, txty = inputarray['which_box'][0]+fac,inputarray['which_box'][3]-fac
        ax.text(txtx,txty,str_enum, fontweight='normal',horizontalalignment='left',\
                verticalalignment='top',fontsize=fsize)
    
    return(ax,hp)

#_______________________________________________________________________________
# plot alimatological anomaly panel matrix
#
# --> required infrastrucutre
#
# nr, nc         = 5, 6 # number of row and column of panel matrix
# var_columns    = [[None]]*nc # init list array list([ [],[],[],[],[],[] ....])
# var_columns[0] = 'tkmin=1.00e-6', '../results/new_linfs/cvmix_tke/'+str(which_cycle)+'/'
# var_columns[1] = 'tkmin=2.00e-6', '../results/new_linfs/cvmix_tke_tkmin2.00e-6/'+str(which_cycle)+'/'
# var_columns[2] = 'tkmin=4.00e-6', '../results/new_linfs/cvmix_tke_tkmin4.00e-6/'+str(which_cycle)+'/'
# var_columns[3] = 'tkmin=0.50e-6', '../results/new_linfs/cvmix_tke_tkmin0.50e-6/'+str(which_cycle)+'/'
# var_columns[4] = 'tkmin=0.25e-6', '../results/new_linfs/cvmix_tke_tkmin0.25e-6/'+str(which_cycle)+'/'
# var_columns[5] = 'tkmin=0.15e-6', '../results/new_linfs/cvmix_tke_tkmin0.15e-6/'+str(which_cycle)+'/'
#
# var_rows       = [[None]]*nr # init list array list([ [],[],[],[],[],[] ....])
# var_rows[0]    = np.arange(   0,  250+1,10)
# var_rows[1]    = np.arange( 250,  500+1,25)
# var_rows[2]    = np.arange( 500, 1000+1,25)
# var_rows[3]    = np.arange(1000, 2000+1,50)
# var_rows[4]    = np.arange(2000, 4000+1,50)
#
# which_year     = [1989,2009]
# which_mon      = list(range(1,12+1))
# which_dep      = []
#
def plot_map_clima_panel_matrix(inputarray, mesh, var_rows, var_columns, which_var, \
    whatisrow='depth', whatiscolumn='descript/path',\
    which_year=[1989,2009], which_mon=list(range(1,12+1)), which_dep=np.arange(0,250+1,10), \
    which_woa='woa2018', do_rmse=True, fsize=10, do_enumerate=True):

    nr, nc = len(var_rows), len(var_columns)
    
    #___________________________________________________________________________
    data1 = fesom_data(inputarray) # init fesom2.0 data object
    data1.var, data1.year, data1.month, data1.depth =  which_var, which_year, which_mon, which_dep

    #___________________________________________________________________________
    if data1.var=='temp' : cmin1,cmax1,cref1,cnumb1  = -4.0, 4.0, 0.0, 20 # [cmin, cmax, cref] --> salt+ 
    if data1.var=='salt' : cmin1,cmax1,cref1,cnumb1  = -0.5, 0.5, 0.0, 20 # [cmin, cmax, cref] --> salt
    cmap1,clevel1= colormap_c2c(cmin1,cmax1,cref1,cnumb1,'blue2red')
    do_drawedges=True

    #___________________________________________________________________________
    # LOAD WOA CLIMATOLOGY
    if which_woa=='woa2005': clim_descr, clim_path, clim_fname  = 'WOA05', '../view/woa2005/', 'woa2005TS.nc'
    if which_woa=='woa2018': clim_descr, clim_path, clim_fname  = 'WOA18', '/work/ollie/pscholz/INIT_HYDRO/woa2018/', 'woa2018TS_0.25deg.nc'
    clim  = clim_data(clim_path,clim_fname,data1.var)
    clim.descript = data1.descript+'-'+clim_descr
    clim.crange,clim.cmap, clim.cnumb      = data1.crange,data1.cmap, data1.cnumb
    clim.str_time, clim.str_dep            = data1.str_time, data1.str_dep
    clim.sname, clim.lname, clim.unit      = data1.sname, data1.lname, data1.unit
    clim.proj, clim.proj_lon, clim.proj_lat= data1.proj, data1.proj_lon, data1.proj_lat

    #___________________________________________________________________________
    # setup figures and axes
    if   nc!=1 and nr!=1: fig, ax = plt.subplots(nrows=nr, ncols=nc, figsize=[4*nc,3*nr])
    elif nc==1          : fig, ax = plt.subplots(nrows=nr, ncols=nc, figsize=[8*nc,4*nr])
    elif nr==1          : fig, ax = plt.subplots(nrows=nr, ncols=nc, figsize=[4*nc,8*nr])
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95,wspace=0.01, hspace=0.01)

    #___________________________________________________________________________
    # start row/column loop
    print(' --> total nr. of matrix panels: {:d} |'.format(nr*nc), end='')
    ni = 0
    for nri in range(0,nr):
        for nci in range(0,nc):
            ni = ni+1
            print('{:d}|'.format(ni), end='')
            #___________________________________________________________________
            # if subplot should be left empty
            if not var_columns[nci]: 
                ax[nri,nci].remove()
                continue
            
            #___________________________________________________________________
            # set row/column specific info
            if   whatisrow    =='depth':
                data1.depth = var_rows[nri]
            elif whatiscolumn =='depth':
                data1.depth = var_columns[nri]
            else: 
                print(' --> ERROR: this panel matrix configuration is not supported')
                exit
            
            if   whatisrow    =='descript/path':
                 data1.descript, data1.path = var_rows[nci][0],var_rows[nci][1]
            elif whatiscolumn =='descript/path':
                data1.descript, data1.path = var_columns[nci][0],var_columns[nci][1]
            else: 
                print(' --> ERROR: this panel matrix configuration is not supported')
                exit
            
            #___________________________________________________________________
            # load data
            fesom_load_data_horiz_netcdf4(mesh,data1,do_output=False)
            data11 = clim_data_anom(mesh,clim,data1,clim_descr)
            
            #___________________________________________________________________
            # compute RMSE
            meanbias1 = np.nansum(np.abs(data11.anom))/np.isnan(data11.anom).sum()
            
            #___________________________________________________________________
            # do axis writtings
            xl, yl, str_tl, str_yl, str_xl = [], [], [], [], []
            if nri == nr-1      : xl = [0,0,0,1]
            if nci == 0         : yl = [1,0,0,0]
            
            if  whatisrow =='depth':
                if nci == 0: str_yl=data11.str_dep[7::]+'\n\n'  
            elif whatisrow =='descript/path':
                if nci == 0: str_yl=data11.descript+'\n\n'  
                    
            if whatiscolumn =='depth':
                if nri == 0: str_tl=data11.str_dep[7::]
            elif whatiscolumn =='descript/path':
                if nri == 0: str_tl=data11.descript
                
            str_err=[]
            if do_rmse:
                if data1.var=='temp': str_err='\n Err={:2.2f}°C'.format(meanbias1)
                if data1.var=='salt': str_err='\n Err={:2.2f}psu'.format(meanbias1)   
            
            #___________________________________________________________________
            # do subplot
            if nc==1 and nr==1:
                dumax,hp = plot_map_clima(fig, ax, mesh, data11, inputarray, clevel1, cmap1, \
                                        xl=xl, yl=yl, str_tl=str_tl, str_yl=str_yl, str_xl=str_xl, \
                                        str_txt=str_err, fsize=fsize )
            
            elif nc!=1 and nr!=1:
                dumax,hp = plot_map_clima(fig, ax[nri,nci], mesh, data11, inputarray, clevel1, cmap1, \
                                        xl=xl, yl=yl, str_tl=str_tl, str_yl=str_yl, str_xl=str_xl, \
                                        str_txt=str_err, fsize=fsize )
            elif nc==1:
                dumax,hp = plot_map_clima(fig, ax[nri], mesh, data11, inputarray, clevel1, cmap1, \
                                        xl=xl, yl=yl, str_tl=str_tl, str_yl=str_yl, str_xl=str_xl, \
                                        str_txt=str_err, fsize=fsize )
            elif nr==1:
                dumax,hp = plot_map_clima(fig, ax[nci], mesh, data11, inputarray, clevel1, cmap1, \
                                        xl=xl, yl=yl, str_tl=str_tl, str_yl=str_yl, str_xl=str_xl, \
                                        str_txt=str_err, fsize=fsize )
            
            #___________________________________________________________________
            # enumerate a, b, c, ....
            fac = 2.5
            if do_enumerate:
                txtx, txty = inputarray['which_box'][0]+fac,inputarray['which_box'][3]-fac
                if nc*nr>26 : str_enum = chr(ord('`')+nci+1)+str(nri+1)+')'
                else        : str_enum = chr(ord('`')+ni)+')'
                if nc==1 and nr==1:
                    ax.text(txtx,txty,str_enum, \
                                    fontweight='normal',horizontalalignment='left',\
                                    verticalalignment='top',fontsize=fsize)
                elif nc!=1 and nr!=1:
                    ax[nri,nci].text(txtx,txty,str_enum, \
                                    fontweight='normal',horizontalalignment='left',\
                                    verticalalignment='top',fontsize=fsize)
                elif nc==1:
                    ax[nri].text(txtx,txty,str_enum, \
                                    fontweight='normal',horizontalalignment='left',\
                                    verticalalignment='top',fontsize=fsize)
                elif nr==1:
                    ax[nci].text(txtx,txty,str_enum, \
                                    fontweight='normal',horizontalalignment='left',\
                                    verticalalignment='top',fontsize=fsize)
    #___________________________________________________________________________
    # do colorbar
    str_unit = data1.unit
    if data1.var=='temp': str_unit = '[°C]'
    cbar1 = setup_colorbar(hp, ax, clevel1, cref1, data1.lname+' anomaly '+' '+str_unit+'\n '+data1.str_time, ncl=10, orientation='vertical')

    #___________________________________________________________________________
    # adapt panel and colorbar positions
    fig.canvas.draw()
    pos_cbar1 = cbar1.ax.get_position()
    if   nc==1 and nr==1: pos_ax    = ax.get_position()
    elif nc!=1 and nr!=1: pos_ax    = ax[0,0].get_position()
    else              : pos_ax    = ax[0].get_position()

    x0,y0  = 0.1, 0.1
    #xg, yg = 0.005, 0.005
    xg, yg = 0.01, 0.01
    w, h   = pos_ax.width, pos_ax.height

    fac    = 0.70/(w*nc+xg*(nc-1))
    w, h   = w*fac, h*fac

    print(w,h)
    for nri in range(0,nr):
        for nci in range(0,nc):
            if   nc==1 and nr==1: ax.set_position( [x0+(w+xg)*nci, y0+(h+yg)*(nr-nri-1), w, h])
            elif nc!=1 and nr!=1: ax[nri,nci].set_position( [x0+(w+xg)*nci, y0+(h+yg)*(nr-nri-1), w, h])
            elif nc==1          : ax[nri].set_position(     [x0+(w+xg)*nci, y0+(h+yg)*(nr-nri-1), w, h])
            elif nr==1          : ax[nci].set_position(     [x0+(w+xg)*nci, y0+(h+yg)*(nr-nri-1), w, h])
    cbar1.ax.set_position([x0+(w+xg)*nc, y0, pos_cbar1.width, nr*h+(nr-1)*yg])
    #cbar1.ax.set_aspect(1/(pos_cbar1.width/(nr*h+(nr-1)*yg)/20*4))
    cbar1.ax.set_aspect(1/(pos_cbar1.width/(nr*h+(nr-1)*yg)/20*4*(5-(nr-1))))
    fig.canvas.draw()        
    
    #___________________________________________________________________________
    # return value
    return(data11)


#_______________________________________________________________________________
# plot single panel of cylindric proection of climatoligical anomaly on regular grid
def save_plot(inputarray,data1,str_general,which_cycle=[],which_res=600, which_format='png', box=[], dname=[]):
    if inputarray['save_fig']==True:
        print(' --> save figure as '+which_format+': ',end='')
        # str_general = 'tunning_overflowtopo_okDS'
        str_descrpt = data1.descript.replace(' ','_') 
        str_times= data1.str_time.replace(' ','').replace(':','')
        str_dep= data1.str_dep[5::].replace(' ','').replace(':','').replace(',','')
        str_reg='regio'
        if (inputarray['which_box'][1]-inputarray['which_box'][0])==360 and  (inputarray['which_box'][3]-inputarray['which_box'][2])>160:
            str_reg='globe'
        sfname = str_general+'_'+data1.which_obj+'_'+str_descrpt
        if not which_cycle==False:
            sfname = sfname+'_'+'scycl:'+str(which_cycle)
        str_name = data1.sname  
        if len(str_name)==0: str_name = data1.var  
        sfname = sfname+'_'+str_reg+'_'+str_name+'_'+str_times+'_'+str_dep+'.png'
        print(sfname,end='')
        sdname = inputarray['save_figpath'] 
        if os.path.isdir(sdname)==False: os.makedirs(sdname, exist_ok=True)
        plt.savefig(sdname+sfname, format=which_format, dpi=which_res, bbox_inches='tight', pad_inches=0,transparent=True)#,frameon=True)
        print(' --> finished!')
