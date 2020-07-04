# Patrick Scholz, 23.01.2018
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from matplotlib.tri import Triangulation,TriAnalyzer
from mpl_toolkits.axes_grid1 import make_axes_locatable
from colormap_c2c import colormap_c2c
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

#___PLOT 2D FESOM DATA__________________________________________________________
# input : data dictionary: data.value, data.sname, data.lname, data.unit
#               data['levels']
#_______________________________________________________________________________
def fesom_plot2d_data(mesh,data,figsize=[],do_subplot=[],do_output=True,do_grid=False, 
                      which_orient='vertical', nmax_cbar_l=8 ):
    if do_output==True:
        print('')
        print('___PLOT 2D DATA____________________________________________')
    from set_inputarray import inputarray
    if (data.value.min() == data.value.max() ):
        print(' --> can\'t plot data are everywhere ', data.value.min())
        sys.exit(' --> STOP HERE!')
    elif (len(np.unique(data.value[~np.isnan(data.value)]))==1 ) :   
        print(' --> can\'t plot data are everywhere ', np.unique(data.value[~np.isnan(data.value)]) )
        sys.exit(' --> STOP HERE!')
        
    
    if len(figsize)==0 : figsize=[13,13]
    #___________________________________________________________________________
    # plot is not part of subplot
    if len(do_subplot)==0:
        fig = plt.figure(figsize=figsize)
        ax  = plt.gca()
    else:
        fig=do_subplot[0]
        ax =do_subplot[1]
        fig.sca(ax)
    resolution = 'c'
    fsize = 14
    
    #+_________________________________________________________________________+
    #| SET PROJECTION PARAMETERS                                               |
    #+_________________________________________________________________________+
    if (data.proj=='ortho' or data.proj=='stere'):
        map = Basemap(projection = data.proj, resolution = resolution,
                      lat_0      = data.proj_lat, lon_0 = data.proj_lon)
        ylabels= [0,0,0,0]
        xlabels= [1,1,1,1]
        xticks = np.arange(-180,180,20)
        yticks = np.arange(-90,90,15)
    elif data.proj=='npstere':
        setbndlat = np.max([0.0,inputarray['which_box'][2]])
        map = Basemap(projection = data.proj, resolution = resolution,
                  lon_0      = 0, boundinglat = setbndlat,
                  round      = True)
        ylabels=[0,0,0,0]
        xlabels=[1,1,1,1]
        xticks = np.arange(0.,360.,20.)
        yticks = np.arange(-90.,90.+1,10.)
    elif data.proj=='spstere':
        setbndlat = np.min([0.0,inputarray['which_box'][3]])
        map = Basemap(projection = data.proj, resolution = resolution,
                  boundinglat= setbndlat, lon_0      = 180,
                  round      = True)
        ylabels=[0,0,0,0]
        xlabels=[1,1,1,1]
        xticks = np.arange(0.,360.,20.)
        yticks = np.arange(-90.,90.+1,10.)
    else:   
        map = Basemap(projection = data.proj,resolution = resolution,
                llcrnrlon  = inputarray['which_box'][0], urcrnrlon  = inputarray['which_box'][1],
                llcrnrlat  = inputarray['which_box'][2], urcrnrlat  = inputarray['which_box'][3],
                )
        ylabels=[1,0,0,0]
        xlabels=[0,0,0,1]
        ticknr   = 8
        tickstep = np.array([0.5,1.0,2.0,2.5,5.0,10.0,15.0,20.0,30.0])
        idx      = (inputarray['which_box'][1]-inputarray['which_box'][0])/ticknr
        idx1      = np.array(np.where(tickstep>=idx))
        if idx1.size==0 : 
            xticks   = np.arange(-180.+tickstep[-1],180.-tickstep[-1]+1,tickstep[-1])
        else:
            xticks   = np.arange(-180.+tickstep[idx1[0,0]],180.-tickstep[idx1[0,0]]+1,tickstep[idx1[0,0]])
        del idx, idx1
        idx      = (inputarray['which_box'][3]-inputarray['which_box'][2])/ticknr
        idx1      = np.array(np.where(tickstep>=idx))
        if idx1.size==0 : 
            yticks   = np.arange(-90.+tickstep[-1],90.-tickstep[-1]+1,tickstep[-1])
        else:    
            yticks   = np.arange(-90.+tickstep[idx1[0,0]],90.-tickstep[idx1[0,0]]+1,tickstep[idx1[0,0]])
        del idx, idx1    
    #___________________________________________________________________________
    # go from geo coord. to projection coord.
    mx,my = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    
    #___________________________________________________________________________
    # make python triangulation object 
    tri      = Triangulation(mx, my,mesh.elem_2d_i)
    
    #___________________________________________________________________________
    # calculate index of which data points are within boi
    idxbox_e = mesh.nodes_2d_xg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][0]
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_xg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][1])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][2])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][3])
    idxbox_e = idxbox_e==False # true index for triangles that are within box 
    
    # case of node data
    if data.value.size ==mesh.n2dna:
        idx_box = np.ones((mesh.n2dna,),dtype='bool') 
    elif data.value.size ==mesh.n2dea:
        idx_box = np.ones((mesh.n2dea,),dtype='bool')
    
    if data.value.size ==mesh.n2dna:
        idxbox_n  = mesh.elem_2d_i[idxbox_e,:].flatten().transpose()
        idxbox_n  = np.array(idxbox_n).squeeze()
        idxbox_n  = np.unique(idxbox_n)
        idx_box   = np.zeros((mesh.n2dna,), dtype=bool)
        idx_box[idxbox_n]=True
        del idxbox_n 
        del idxbox_e
    # case of element data
    elif data.value.size ==mesh.n2dea:
        idx_box   = idxbox_e
        del idxbox_e
    
    #___________________________________________________________________________
    # make triangle mask arrays from flat triangles and nan trinagles
    mask_tri = TriAnalyzer(tri).get_flat_tri_mask(0.00001)
    
    if np.any(np.isnan(data.value)):
        if data.value.size ==mesh.n2dna:
            mask_tri = np.logical_or(mask_tri,np.isnan(np.sum(data.value[mesh.elem_2d_i],axis=1)))
        elif data.value.size ==mesh.n2dea:    
            mask_tri = np.logical_or(mask_tri,np.isnan(data.value))
        idx_box = np.logical_and(idx_box,np.isnan(data.value)==False)
    
    if data.proj=='npstere':
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3<setbndlat)
    elif data.proj=='spstere':    
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3>setbndlat)
    tri.set_mask(mask_tri)
    
    #+_________________________________________________________________________+
    #| set minimum, maximum and reference values for the creation of the       |
    #| adjustable colormap                                                     |
    #+_________________________________________________________________________+
    cnumb    = 20; # minimum number of colors
    if np.size(data.cnumb)!=0:
        cnumb = data.cnumb
    # predined color ranges
    if data.sname=='a_ice':
        cmax,cmin,cref = 100.0, 0.0, 50.0
        data.cmap='wbgyr'
    elif data.sname=='m_ice':
            cmax = np.nanmax(data.value[idx_box])
            cmin = np.nanmin(data.value[idx_box])
            cref = cmin + (cmax-cmin)/2
            cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
            data.cmap='wbgyr'    
    elif data.sname=='ssh':
        cmax = np.max(data.value[idx_box])
        cmin = np.min(data.value[idx_box])
        cref = 0.0
    elif data.sname=='salt' or data.sname=='sss':
        cmax = np.max(data.value[idx_box])
        cmin = 20.0
        cref = cmin + (cmax-cmin)/2
        cref = np.round_(cref, np.int32(np.floor(np.log10(np.abs(cref)))) ) 
    elif data.sname=='triresol' or data.sname=='triarea' or data.sname=='depth':
        cmax = np.max(data.value[idx_box])
        cmin = 0
        cref = cmin + (cmax-cmin)/2
        cref = np.around(cref, -np.int32(np.floor(np.log10(cref))) ) 
    elif data.sname=='u' or data.sname=='v' or data.sname=='w' or \
         data.sname=='bolus_u' or data.sname=='bolus_v'  or data.sname=='bolus_w':    
        cmax = np.nanmax(data.value[idx_box])
        cmin = np.nanmin(data.value[idx_box])
        cref=0
        data.cmap = 'blue2red'
    elif data.sname=='norm_uv':
        cmax = np.max(data.value[idx_box])
        cmin = 0.0
        cref = cmin + (cmax-cmin)/2
        cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))) ) 
        data.cmap='wbgyr'
    elif data.var.find('MLD')!=-1:
        #data.cmap='rygbw'
        data.cmap='wbgyr'
        cmin = 0.0
        cmax = np.max(data.value[idx_box])
        cref = cmin + (cmax-cmin)/2
        cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
    else:
        cmax = np.max(data.value[idx_box])
        cmin = np.min(data.value[idx_box])
        cref = cmin + (cmax-cmin)/2
        cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1) ) 
    #___________________________________________________________________________
    # if anomaly data
    if data.anom==True: 
        cmax,cmin,cref = np.max(data.value[idx_box]), np.min(data.value[idx_box]), 0.0
        data.cmap='blue2red'
    
    #___________________________________________________________________________
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
    
    #___________________________________________________________________________
    # make custom colormap
    if do_output==True: print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
    cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,data.cmap)
    if do_output==True:print('clevel = ',clevel)
    do_drawedges=True
    if clevel.size>50: do_drawedges=False
    
    #if data.sname=='a_ice' or data.sname=='m_ice' or data.sname.find('MLD')==0:
    if 'MLD' in data.var:
        if data.anom==False:
            # make no sea ice transparent
            from matplotlib.colors import ListedColormap
            auxcmap = cmap0(np.arange(cmap0.N))
            #if data.sname=='a_ice' or data.sname=='m_ice' :
                ##auxcmap[0,-1]=0.0
            #elif data.sname.find('MLD')==0 :
            if data.sname.find('MLD')==0 :    
                #auxcmap[-1,-1]=0.0
                auxcmap[0,-1]=0.0
            cmap0 = ListedColormap(auxcmap)
    
    #___________________________________________________________________________
    # limit minimum and maximum of data vector to setted colorrange avoid holes 
    # in the plot
    data_plot = np.copy(data.value)
    data_plot[data_plot<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
    data_plot[data_plot>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    
    data.levels = clevel
    #+_________________________________________________________________________+
    #| plot data on triangluar grid                                            |
    #+_________________________________________________________________________+
    # plot data defined on nodes 
    if data.value.size==mesh.n2dna:
        if data.which_plot=='pcolor':
            hp1=ax.tripcolor(tri,data_plot,                              
                antialiased=False,
                edgecolors='None',
                cmap=cmap0,
                shading='gouraud',
                clim=[clevel[0],clevel[-1]],
                vmin=clevel[0],vmax=clevel[-1])
                #shading='flat')
                #shading='gouraud')
            if do_grid==True: ax.triplot(tri,color='k',linewidth=.5,alpha=0.25)        
        elif data.which_plot=='contourf':
            hp1=ax.tricontourf(tri,data_plot,
                levels=clevel, 
                antialiased=False,
                extend='both',
                cmap=cmap0)
            if do_grid==True: ax.triplot(tri,color='k',linewidth=.5,alpha=0.25)
    # plot data defined on elements
    elif data.value.size==mesh.n2dea:
        hp1=ax.tripcolor(tri,data_plot,                          
                    antialiased=False,
                    cmap=cmap0,
                    clim=[clevel[0],clevel[-1]],
                    vmin=clevel[0],vmax=clevel[-1])
        
        if do_grid==True: ax.triplot(tri,color='k',linewidth=.5,alpha=0.15)
    
    #___________________________________________________________________________
    # arange zonal & meriodional gridlines and labels
    map.drawmapboundary(fill_color='0.9',linewidth=1.0)
    
    # label lon lat grid for ortho projection 
    if data.proj in ['ortho','stere']:
        map.drawparallels([0.0],
                linewidth=1.0,
                dashes=[1,1e-10],
                fontsize=fsize)
        map.drawmeridians([0.0],
                linewidth=1.0,
                dashes=[1,1e-10],
                fontsize=fsize)
        
        for i in np.arange(len(xticks)):
            xt,yt=xy=map(xticks[i],0)
            if xticks[i]>0 :
                plt.text(xt,yt,' {:.1f}$^{{\\circ}}$W'.format(xticks[i]),fontsize=12,fontweight='bold',verticalalignment='center',horizontalalignment='left')
            elif xticks[i]<0: 
                plt.text(xt,yt,' {:.1f}$^{{\\circ}}$W'.format(xticks[i]),fontsize=12,fontweight='bold',verticalalignment='center',horizontalalignment='left')
            else: 
                plt.text(xt,yt,' {:.1f}$^{{\\circ}}$'.format(xticks[i]),fontsize=12,fontweight='bold',verticalalignment='center',horizontalalignment='left')
                plt.plot(xt,yt,'o',linestyle='None',color='k')    
        for i in np.arange(len(yticks)):
            xt,yt=xy=map(0,yticks[i])
            if yticks[i]>0 :
                plt.text(xt,yt,' {:.1f}$^{{\\circ}}$N'.format(yticks[i]),fontsize=12,fontweight='bold',verticalalignment='center',horizontalalignment='left')
            elif yticks[i]<0: 
                plt.text(xt,yt,' {:.1f}$^{{\\circ}}$S'.format(yticks[i]),fontsize=12,fontweight='bold',verticalalignment='center',horizontalalignment='left')
                plt.plot(xt,yt,'o',linestyle='None',color='k')    
    # label lon lat grid for cycl projection             
    else:
        map.drawparallels(yticks,#np.arange(-90.,90.,30.),
            labels=ylabels,
            linewidth=0.25,
            dashes=[1,1e-10],
            fontsize=fsize)
        meridians = map.drawmeridians(xticks,#np.arange(0.,360.,30.),
            linewidth=0.25,
            labels=xlabels,
            dashes=[1,1e-10],
            fontsize=fsize)
            
        # try to rotate meridian labels
        if data.proj=='cyl':
            for m in meridians:
                try:
                    meridians[m][1][0].set_rotation(25)
                except:
                    pass
        
    #___________________________________________________________________________
    # draw land mask patch
    if inputarray['which_mask'] == 'fesom':
        fesom_plot_lmask(map,mesh,ax,'0.6')
    elif inputarray['which_mask'] == 'bluemarble':
        map.bluemarble()
        fesom_plot_lmask(map,mesh,ax,'none')
    elif inputarray['which_mask'] == 'etopo':
        map.etopo()
        fesom_plot_lmask(map,mesh,ax,'none')
    
    #___________________________________________________________________________
    # arange colorbar position and labels
    if (data.proj=='ortho'   or data.proj=='stere' or \
        data.proj=='npstere' or data.proj=='spstere'):
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="2.5%", pad=0.70)
    else:
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="2.5%", pad=0.1)
    
    # give each subplot its own colorrange
    #plt.clim(clevel[0],clevel[-1])
    him = ax.get_images()
    for im in ax.get_images():
        im.set_clim(clevel[0],clevel[-1])
    
    cbar = plt.colorbar(hp1,cax=cax,ticks=clevel,drawedges=do_drawedges, \
                        extend='neither',extendrect=False,extendfrac=None,\
                        orientation=which_orient)
    cbar.set_label(data.lname+' '+data.unit+'\n'+data.str_time+data.str_dep, size=fsize+2)
    cl = plt.getp(cbar.ax, 'ymajorticklabels')
    plt.setp(cl, fontsize=fsize)
    
    # kickout some colormap labels if there are to many
    ncbar_l=len(cbar.ax.get_yticklabels()[:])
    idx_cref = np.where(clevel==cref)[0]
    idx_cref = np.asscalar(idx_cref)
    
    #nmax_cbar_l = 12
    nstep = ncbar_l/nmax_cbar_l
    nstep = np.int(np.floor(nstep))
    if nstep==0:nstep=1
    
    fig.canvas.draw() # this is need so cbar.ax.get_yticklabels() always finds the labels
    if cbar.orientation=='vertical':
        tickl = cbar.ax.get_yticklabels()
    else:
        tickl = cbar.ax.get_xticklabels()
    
    idx = np.arange(0,len(tickl),1)
    idxb = np.ones((len(tickl),), dtype=bool)                
    idxb[idx_cref::nstep]  = False
    idxb[idx_cref::-nstep] = False
    idx = idx[idxb==True]
    for ii in list(idx):
        tickl[ii]=''
    if cbar.orientation=='vertical':    
        cbar.ax.set_yticklabels(tickl)
    else:    
        cbar.ax.set_xticklabels(tickl)
    #___________________________________________________________________________
    ax.tick_params(axis='both', direction='out')
    #ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
    #ax.get_yaxis().tick_left()
    plt.sca(ax)
    
    if data.proj=='npstere' or data.proj=='spstere' or data.proj=='ortho':
        plt.title(data.descript+'\n',fontdict= dict(fontsize=20),verticalalignment='baseline')
    else:    
        plt.title(data.descript,fontdict= dict(fontsize=20),verticalalignment='center')
#         plt.title(data.descript+'\n',fontdict= dict(fontsize=24),verticalalignment='center')
    #ax.set_xlabel(' ',fontsize=20,verticalalignment='top')
    #ax.set_ylabel(' ',fontsize=30,verticalalignment='bottom')
    #+_________________________________________________________________________+
    #| SAVE FIGURE                                                             |
    #+_________________________________________________________________________+
    if inputarray['save_fig']==True:
        print(' --> save figure: png')
        #print(fig.dpi)
        str_times= data.str_time.replace(' ','').replace(':','')
        str_deps = data.str_dep.replace(' ','').replace(',','').replace(':','')
        sfname = 'plot_'+data.descript+'_'+data.proj+'_'+data.sname+'_'+str_times+'_'+str_deps+'.png'
        sdname = inputarray['save_figpath']
        if os.path.isdir(sdname)==False: os.makedirs(sdname)
        plt.savefig(sdname+sfname, \
                    format='png', dpi=600, \
                bbox_inches='tight', pad_inches=0,\
                transparent=False,frameon=True)
    
    #___________________________________________________________________________
    plt.show(block=False)
    fig.canvas.draw()
    #___________________________________________________________________________
    return(fig,ax,map,cbar,hp1,tri)
    #return map
    
    
#___PLOT 2D VEC FESOM DATA__________________________________________________________
# input : data dictionary: data.value, data.sname, data.lname, data.unit
#               data['levels']
#_______________________________________________________________________________
def fesom_plot2dvec_data(mesh,data,figsize=[],do_subplot=[],do_output=True,scal_fac=1.0, nmax_cbar_l=8):
    if do_output==True: 
        print('')
        print('___PLOT 2D VECTOR DATA____________________________________________')
    from set_inputarray import inputarray
    
    if len(figsize)==0 : figsize=[13,13]
    #___________________________________________________________________________
    # plot is not part of subplot
    if len(do_subplot)==0:
        fig = plt.figure(figsize=figsize)
        ax  = plt.gca()
    else:
        fig=do_subplot[0]
        ax =do_subplot[1]
        plt.sca(ax)
    resolution = 'l'
    fsize = 14
    
    #+_________________________________________________________________________+
    #| SET PROJECTION PARAMETERS                                               |
    #+_________________________________________________________________________+
    if (data.proj=='ortho' or data.proj=='stere'):
        map = Basemap(projection = data.proj,
                  resolution = resolution,
                  lat_0      = data.proj_lat,
                  lon_0      = data.proj_lon)
        ylabels= [0,0,0,0]
        xlabels= [1,1,1,1]
        xticks = np.arange(-180,180,20)
        yticks = np.arange(-90,90,15)
    elif data.proj=='npstere':
        setbndlat = np.max([0.0,inputarray['which_box'][2]])
        map = Basemap(projection = data.proj,
                  resolution = resolution,
                  lon_0      = 0, 
                  boundinglat= setbndlat,
                  round      = True)
        ylabels=[0,0,0,0]
        xlabels=[1,1,1,1]
        xticks = np.arange(0.,360.,20.)
        yticks = np.arange(-90.,90.,15.)
    elif data.proj=='spstere':
        setbndlat = np.min([0.0,inputarray['which_box'][3]])
        map = Basemap(projection = data.proj,
                  resolution = resolution,
                  boundinglat= setbndlat,
                  lon_0      = 180,
                  round      = True)
        ylabels=[0,0,0,0]
        xlabels=[1,1,1,1]
        xticks = np.arange(0.,360.,20.)
        yticks = np.arange(-90.,90.,15.)
    else:    
        map = Basemap(projection = data.proj,
                  resolution = resolution,
                  llcrnrlon  = inputarray['which_box'][0],
                  urcrnrlon  = inputarray['which_box'][1],
                  llcrnrlat  = inputarray['which_box'][2],
                  urcrnrlat  = inputarray['which_box'][3])
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
    mx,my = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    tri = Triangulation(mx,my,mesh.elem_2d_i)
    
    #___________________________________________________________________________
    # calculate index of which data points are within boi
    idxbox_e = mesh.nodes_2d_xg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][0]
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_xg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][1])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][2])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][3])
    idxbox_e = idxbox_e==False # true  index for triangles that are within box 
    
    # case of node data
    if data.value.size ==mesh.n2dna:
        idx_box = np.ones((mesh.n2dna,),dtype='bool')   
    elif data.value.size ==mesh.n2dea:
        idx_box = np.ones((mesh.n2dea,),dtype='bool')
        
    if data.value.size ==mesh.n2dna:
        idxbox_n  = mesh.elem_2d_i[idxbox_e,:].flatten().transpose()
        idxbox_n  = np.array(idxbox_n).squeeze()
        idxbox_n  = np.unique(idxbox_n)
        idx_box   = np.zeros((mesh.n2dna,), dtype=bool)
        idx_box[idxbox_n]=True
        del idxbox_n 
        del idxbox_e
    # case of element data
    elif data.value.size ==mesh.n2dea:
        idx_box   = idxbox_e
        del idxbox_e
    
    #___________________________________________________________________________
    # make triangle mask arrays from flat triangles and nan trinagles
    mask_tri = TriAnalyzer(tri).get_flat_tri_mask(0.00001)
    if np.any(np.isnan(data.value)):
        if data.value.size ==mesh.n2dna:
            mask_tri = np.logical_or(mask_tri,np.isnan(np.sum(data.value[mesh.elem_2d_i],axis=1)))
        elif data.value.size ==mesh.n2dea:    
            mask_tri = np.logical_or(mask_tri,np.isnan(data.value))
        idx_box = np.logical_and(idx_box,np.isnan(data.value)==False)
    
    if data.proj=='npstere':
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3<setbndlat)
    elif data.proj=='spstere':    
        mask_tri = np.logical_or(mask_tri,np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3>setbndlat)
    tri.set_mask(mask_tri)
    
    #___________________________________________________________________________
    plt.sca(ax)
    #ax.tricontourf(tri,mesh.nodes_2d_zg,levels=np.array([-6000,-5500,-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-750,-500,-250,-150,-80,-40,-20]),\
                    #cmap=cm.gray,alpha=0.75,antialiased=False,extend='both',edgealpha=.75)
    data_plotu = np.copy(data.value)
    data_plotv = np.copy(data.value2)
    
    #___________________________________________________________________________
    if data_plotu.size==mesh.n2dea:
        mx_in = np.sum(mesh.nodes_2d_xg[mesh.elem_2d_i],axis=1)/3
        my_in = np.sum(mesh.nodes_2d_yg[mesh.elem_2d_i],axis=1)/3
    else:
        mx_in = mesh.nodes_2d_xg
        my_in = mesh.nodes_2d_yg
    mx_in = mx_in[idx_box]    
    my_in = my_in[idx_box]
    data_plotu = data_plotu[idx_box]    
    data_plotv = data_plotv[idx_box]
    
    #___________________________________________________________________________
    data_plotu,data_plotv,mx,my = map.rotate_vector(data_plotu,data_plotv, mx_in, my_in, returnxy=True)
    if data.proj=='spstere': data_plotu,data_plotv = -data_plotu, -data_plotv
    
    #___________________________________________________________________________
    # limit all arrays to plotable box 
    idx      = ~np.isnan(data_plotu)
    mx         = mx[idx]
    my         = my[idx]
    data_plotu  = data_plotu[idx]
    data_plotv  = data_plotv[idx]
    
    idx      = np.where(data_plotu!=0.0)
    mx         = mx[idx]
    my         = my[idx]
    data_plotu  = data_plotu[idx]
    data_plotv  = data_plotv[idx]
    
    norm        = np.sqrt(data_plotu**2.0+data_plotv**2.0)
    norm_orig   = np.copy(norm)
    
    #___________________________________________________________________________
    cnumb    = 10; #minimum number of colors
    data.cmap = 'wbgyr'
    cmin     = 0 ; 
    cmax     = np.max(norm)
    #cmax     = 0.25
    cref     = cmin + (cmax-cmin)/2
    cref     = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))) ) 
    
    #___________________________________________________________________________
    # if predefined color range --> data.crange
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
    cnumb = data.cnumb    
    if do_output==True: print('[cmin,cmax,cref] = ['+str(cmin)+', '+str(cmax)+', '+str(cref)+']')
    cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,data.cmap)
    if do_output==True: print('clevel = ',clevel)
    
    #___________________________________________________________________________
    # limit minimum and maximum of vector norm
    norm[norm<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
    norm[norm>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    
    #___________________________________________________________________________
    # shrink all vectors that are originaly longer than norm to the norm length
    data_plotu = data_plotu*norm/norm_orig
    data_plotv = data_plotv*norm/norm_orig
    
    #___________________________________________________________________________
    # eliminate all vectors that are close to 0
    #idx     = norm>=clevel[-1]*0.01
    #idx     = norm>=clevel[1]*0.5
    #mx         = mx[idx]
    #my         = my[idx]
    #data_plotu  = data_plotu[idx]
    #data_plotv  = data_plotv[idx]
    #norm        = norm[idx]
    #del idx 
    
    # normalize vector components
    #data_plotu  = data_plotu/clevel[-1]
    #data_plotv  = data_plotv/clevel[-1]
    
    # limit smalles size of vectors
    
    #limit = cmax/5
    #limit = cmax/10
    limit = clevel[-1]*0.15
    aux = (data_plotu**2.0+data_plotv**2.0)**0.5
    data_plotu[aux<limit]=data_plotu[aux<limit]/aux[aux<limit]*limit
    data_plotv[aux<limit]=data_plotv[aux<limit]/aux[aux<limit]*limit
    del aux
    
    #___________________________________________________________________________
    #hp1=ax.quiver(mx,my,data_plotu,data_plotv,norm,#  
    #hp1=map.quiver(mx,my,data_plotu,data_plotv,norm,#    
                    #scale_units='xy',scale=None, #np.max(norm)*0.01,#scale=np.max(norm)*1.15, #25,#scale_units='xy'/'width',scale=clevel[-1],
                    #edgecolor='k',
                    #linewidth=0.1,
                    #cmap = cmap0,#alpha=0.9,
                    #pivot='middle', #'middle',
                    #antialiased=True,
                    #minlength=0.0,
                    #width=0.005,
                    #headlength=3,
                    #headaxislength=3,
                    #headwidth=2.5)
    
    aux = np.sqrt((mx.max()-mx.min())**2 + (my.max()-my.min())**2)
    print(norm.max())
    
    
    
    #hp1=ax.quiver(mx,my,data_plotu,data_plotv,norm,#  
                    #scale_units='inches',scale=0.6, #np.max(norm)*0.01,#scale=np.max(norm)*1.15, #25,#scale_units='xy'/'width',scale=clevel[-1],
                    #edgecolor='k',
                    #linewidth=0.1,
                    #cmap = cmap0,#alpha=0.9,
                    #pivot='middle', #'middle',
                    #antialiased=True,
                    #minlength=0.0,
                    #width=0.005,
                    #headlength=2.5,
                    #headaxislength=2.5,
                    #headwidth=1.5)
    
    mx,my,data_plotu,data_plotv,norm = mx[::2], my[::2], data_plotu[::2], data_plotv[::2], norm[::2]
    
    hp1=map.quiver(mx,my,data_plotu,data_plotv,norm,#  
                    scale_units='xy',scale=1/aux*50*scal_fac*clevel[-1], #scale=1/aux*80*scal_fac*clevel[-1], #np.max(norm)*0.01,#scale=np.max(norm)*1.15, #25,#scale_units='xy'/'width',scale=clevel[-1],
                    edgecolor='k',
                    linewidth=0.1,
                    cmap = cmap0,#alpha=0.9,
                    pivot='middle', #'middle',
                    antialiased=True,
                    minlength=0.0,
                    width=0.005,
                    headlength=3,
                    headaxislength=3,
                    headwidth=2.0)
    
    hp1.set_clim([clevel[0],clevel[-1]])
    # draw reference vector in corner
    #hpref = plt.quiverkey(hp1, 0.85, 0.05,1.0 ,'${:.2f}\\ \\frac{{m}}{{s}}$'.format(clevel[-1]), labelpos='E',
            #coordinates='axes',fontproperties={"size":20, "weight":"extra bold"})
    #hpref = plt.quiverkey(hp1, 0.1, 0.05,1.0 ,'${:.2f}\\ \\frac{{m}}{{s}}$'.format(clevel[-1]), labelpos='E',
            #coordinates='axes',fontproperties={"size":20, "weight":"extra bold"})
    
    #hpref = ax.quiverkey(hp1, 0.05, 0.05,clevel[-1],'${:.2f}\\ \\frac{{m}}{{s}}$'.format(clevel[-1]), labelpos='E',
            #coordinates='axes',fontproperties={"size":20, "weight":"extra bold"})
    
    #___________________________________________________________________________
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
    
    # label lon lat grid for ortho projection 
    if (data.proj=='ortho'   or data.proj=='stere'):
        for i in np.arange(len(xticks)):
            plt.annotate(np.str(xticks[i]),xy=map(xticks[i],0),xycoords='data')
        for i in np.arange(len(yticks)):
            plt.annotate(np.str(yticks[i]),xy=map(0,yticks[i]),xycoords='data')
    
    #___________________________________________________________________________
    # draw land mask patch
    if inputarray['which_mask'] == 'fesom':
        fesom_plot_lmask(map,mesh,ax,'0.6')
    elif inputarray['which_mask'] == 'bluemarble':
        map.bluemarble()
        fesom_plot_lmask(map,mesh,ax,'none')
    elif inputarray['which_mask'] == 'etopo':
        map.etopo()
        fesom_plot_lmask(map,mesh,ax,'none')
    
    #___________________________________________________________________________
    # arange colorbar position and labels
    if (data.proj=='ortho'   or data.proj=='stere' or \
        data.proj=='npstere' or data.proj=='spstere'):
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="2.5%", pad=0.70)
    else:
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="2.5%", pad=0.1)
    # give each subplot its own colorrange
    #plt.clim(clevel[0],clevel[-1])
  
    cbar = plt.colorbar(hp1,cax=cax,ticks=clevel,drawedges=True,
            extend='neither',extendrect=True,extendfrac='auto')
    
    him = ax.get_images()
    for im in ax.get_images():
        im.set_clim(clevel[0],clevel[-1])
    
    
    cbar.set_label(data.lname+' '+data.unit+'\n'+data.str_time+data.str_dep, size=fsize+2)
    cl = plt.getp(cbar.ax, 'ymajorticklabels')
    plt.setp(cl, fontsize=fsize)
    
    # kickout some colormap labels if there are to many
    ncbar_l=len(cbar.ax.get_yticklabels()[:])
    idx_cref = np.where(clevel==cref)[0]
    idx_cref = np.asscalar(idx_cref)
    
    #nmax_cbar_l = 12
    nstep = ncbar_l/nmax_cbar_l
    nstep = np.int(np.floor(nstep))
    if nstep==0:nstep=1
    
    fig.canvas.draw() # this is need so cbar.ax.get_yticklabels() always finds the labels
    if cbar.orientation=='vertical':
        tickl = cbar.ax.get_yticklabels()
    else:
        tickl = cbar.ax.get_xticklabels()
    idx = np.arange(0,len(tickl),1)
    idxb = np.ones((len(tickl),), dtype=bool)                
    idxb[idx_cref::nstep]  = False
    idxb[idx_cref::-nstep] = False
    idx = idx[idxb==True]
    for ii in list(idx):
        tickl[ii]=''
    if cbar.orientation=='vertical':    
        cbar.ax.set_yticklabels(tickl)
    else:    
        cbar.ax.set_xticklabels(tickl)
    #___________________________________________________________________________
    ax.tick_params(axis='both', direction='out')
    plt.sca(ax)
    
    #plt.title(data.descript+'\n',fontdict= dict(fontsize=24),verticalalignment='center')
    if data.proj=='npstere' or data.proj=='spstere' or data.proj=='ortho':
        plt.title(data.descript+'\n',fontdict= dict(fontsize=24),verticalalignment='baseline')
    else:    
        plt.title(data.descript+'\n',fontdict= dict(fontsize=24),verticalalignment='center')
    #___________________________________________________________________________
    if inputarray['save_fig']==True:
        print(' --> save figure: png')
        #print(fig.dpi)
        sdname = inputarray['save_figpath']
        if os.path.isdir(sdname)==False: os.makedirs(sdname)
        plt.savefig(sdname+'plot_'+data.proj+'_'+data.sname+'.png', \
                format='png', dpi=600, \
                bbox_inches='tight', pad_inches=0,\
                transparent=True,frameon=True)
    
    #___________________________________________________________________________
    plt.show(block=False)
    
    #___________________________________________________________________________
    return(fig,ax,cbar,hp1)
    
    
#___PLOT 2D FESOM MESH IN GEO COORDINATES_______________________________________
#
#
#_______________________________________________________________________________
def fesom_plot2d_geomesh(mesh):
    print('')
    print('___PLOT 2D GEO MESH________________________________________')
    from set_inputarray import inputarray
    
    #___________________________________________________________________________
    fig = plt.figure()
    ax  = plt.gca()
    #___________________________________________________________________________
    map = Basemap(projection ='cyl',
          resolution = 'h',
          llcrnrlat  = -90,
          urcrnrlat  = 90,
          llcrnrlon  = -180,
          urcrnrlon  = 180)
    
    mx,my = map(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    
    #make python triangulation object 
    tri  = Triangulation(mx, my,mesh.elem_2d_i)
    #tri  = Triangulation(mx, my,mesh['pos_i'][mesh['elem_2d_i']])
    mask = TriAnalyzer(tri).get_flat_tri_mask(0.00001)
    #tri.set_mask(mask)
    
    #___________________________________________________________________________
    plt.triplot(tri,linewidth=0.1,color='k')
    #plt.tripcolor(tri,mesh['nodes_2d_xg']*0,
              #edgecolors='k',
          #facecolors=np.ones((mesh['n2dea'],1)),
          #linewidth=0.25)
    #___________________________________________________________________________
    if inputarray['which_mask'] == 'fesom':
        #if data.proj!='ortho':
        fesom_plot_lmask(map,mesh,ax,'0.6')
    elif inputarray['which_mask'] == 'bluemarble':
        map.bluemarble()
        #if data.proj!='ortho':
        fesom_plot_lmask(map,mesh,ax,'none')
    elif inputarray['which_mask'] == 'etopo':
        map.etopo()
        #if data.proj!='ortho':
        fesom_plot_lmask(map,mesh,ax,'none')
    
    #___________________________________________________________________________
    map.drawparallels(np.arange(-90.,90.,15.),
            linewidth=0.5,
            labels=[1,0,0,0],
            fontsize=8,
            dashes=[1,1e-10])
    map.drawmeridians(np.arange(0.,360.,30.),
            linewidth=0.5,
            labels=[0,0,0,1],
            fontsize=8,
            dashes=[1,1e-10])
    map.drawmapboundary(fill_color='1.0',linewidth=1.0)
    
    #___________________________________________________________________________
    plt.show(block=False)
    
    #___________________________________________________________________________
    if inputarray['save_fig']==True:
        print(' --> save figure')
        #fig.savefig('mesh_resolution.eps',dpi=300)
        sdname = inputarray['save_figpath']
        if os.path.isdir(sdname)==False: os.makedirs(sdname)
        plt.savefig(sdname+'mesh_geo_coord.png', \
                format='png', dpi=600, \
                bbox_inches='tight', pad_inches=0,\
                transparent=True,frameon=True)
    
    #___________________________________________________________________________
    return map

    
#___PLOT 2D FESOM MESH IN GEO COORDINATES_______________________________________
#
#
#_______________________________________________________________________________
def fesom_plot2d_rotmesh(mesh):
    print('')
    print('___PLOT 2D ROT MESH________________________________________')
    from set_inputarray import inputarray
    
    #___________________________________________________________________________
    fig = plt.figure()
    ax  = plt.gca()
    #___________________________________________________________________________
    map = Basemap(projection ='cyl',
          resolution = 'i',
          llcrnrlat  = -90,
          urcrnrlat  = 90,
          llcrnrlon  = np.min(mesh.nodes_2d_xg),
          urcrnrlon  = np.max(mesh.nodes_2d_xg))
    
    mx,my = map(mesh.nodes_2d_xr, mesh.nodes_2d_yr)
    
    #make python triangulation object 
    tri  = Triangulation(mx, my,mesh.elem_2d_i)
    mask = TriAnalyzer(tri).get_flat_tri_mask(0.1)
    tri.set_mask(mask)
    
    #___________________________________________________________________________
    plt.triplot(tri,linewidth=0.2,color='k')
    #plt.tripcolor(tri,mesh['nodes_2d_xg']*0,
              #edgecolors='k',
          #facecolors=np.ones((mesh['n2dea'],1)),
          #linewidth=0.25)
    
    #___________________________________________________________________________
    map.drawparallels(np.arange(-90.,90.,15.),
            linewidth=0.5,
            labels=[1,0,0,0],
            dashes=[1,1e-10])
    map.drawmeridians(np.arange(0.,360.,30.),
            linewidth=0.5,
            labels=[0,0,0,1],
            fontsize=12,
            dashes=[1,1e-10])
    map.drawmapboundary(fill_color='1.0',linewidth=1.0)
    
    #___________________________________________________________________________
    plt.show(block=False)
    
    #___________________________________________________________________________
    if inputarray['save_fig']==True:
        print(' --> save figure')
        sdname = inputarray['save_figpath']
        #fig.savefig('mesh_resolution.eps',dpi=300)
        if os.path.isdir(sdname)==False: os.makedirs(sdname)
        plt.savefig(sdname+'mesh_rot_coord.png', format='png', dpi=600)
    
    #___________________________________________________________________________
    return map

#___PLOT LAND MASK PATCH________________________________________________________
#
#
#_______________________________________________________________________________
def fesom_plot_lmask(map,mesh,ax,fcolor,ecolor='k'):
    
    #___________________________________________________________________________
    m_patches=[]
    npoly = np.shape(mesh.polygon_xy)
    #for ii in range(0,npoly[0]):
    for ii in range(0,npoly[0]):
        #_______________________________________________________________________
        aux_x = np.array(mesh.polygon_xy[ii][:,0])
        aux_y = np.array(mesh.polygon_xy[ii][:,1])
        #_______________________________________________________________________
        if (map.projection=='ortho' or map.projection=='stere'):
            # in case of ortho projection basemap has some problems. landmask
            # points that are outside the seeable ortho domain are set to 1e30. 
            # I try to set e points so that they lie on the outer boundary 
            # of the ortho domain
            if np.any(aux_y<map.boundinglat) and map.projparams['lat_0']==90:
                if np.all(aux_y<map.boundinglat):
                    continue
                ind = np.where(aux_y<map.boundinglat);
                aux_y[ind]=map.boundinglat;
            elif np.any(aux_y>map.boundinglat) and map.projparams['lat_0']==-90:
                if np.all(aux_y>map.boundinglat):
                    continue
                ind = np.where(aux_y>map.boundinglat);
                aux_y[ind]=map.boundinglat;
            else:
                # lon/lat coordinates of outer boundary
                outerbnd_lon=map.boundarylons
                outerbnd_lat=map.boundarylats

                # outer boundary in projection coordinates
                outerbnd_x,outerbnd_y = map(outerbnd_lon, outerbnd_lat)

                # landmask points in projection coord.
                mlmaskx,mlmasky = map(aux_x, aux_y)

                # seek for bad land mask points x or y = 1e30 in projection coord
                bad_xy    = np.where((mlmaskx>1e15)|(mlmasky>1e15))
                if np.size(np.array(bad_xy))==mlmaskx.size:
                    continue
                # look what are thes bad points in lon/lat coordinates
                bad_lon   = aux_x[bad_xy]
                if map.projparams['lat_0']>=0:
                    bad_lon[bad_lon<=np.min(outerbnd_lon)]=bad_lon[bad_lon<=np.min(outerbnd_lon)]+360
                elif map.projparams['lat_0']<0:
                    bad_lon[bad_lon>=np.max(outerbnd_lon)]=bad_lon[bad_lon>=np.max(outerbnd_lon)]-360

                # interpolate
                bad_lat   = np.interp(bad_lon,outerbnd_lon,outerbnd_lat,period=360)

                # set lat coordinates of bad landmask points to lat values of the
                # outer boundary domain
                aux_y[bad_xy]=bad_lat

            #fig = plt.figure()
            #plt.plot(aux_x,aux_y)
            #plt.show(block=False)
            #STOP
        #_______________________________________________________________________
        elif map.projection=='npstere':
            if np.all(aux_y<0):
                continue
            elif np.any(aux_y<map.boundinglat):
                aux_y[aux_y<map.boundinglat]=map.boundinglat
        #_______________________________________________________________________
        elif map.projection=='spstere':
            if np.any(aux_y==-90):
                ind = np.where(aux_y!=-90);
                aux_x = aux_x[ind]
                aux_y = aux_y[ind]
            elif np.any(aux_y>map.boundinglat):
                aux_y[aux_y>map.boundinglat]=map.boundinglat
        #_______________________________________________________________________
        mlmaskx,mlmasky = map(aux_x, aux_y)

        del aux_x; del aux_y

        mlmaskx  = np.reshape(mlmaskx,(mlmaskx.size,1))
        mlmasky  = np.reshape(mlmasky,(mlmasky.size,1))

        bad_x    = np.array(np.where(mlmaskx>1e15,0,1))
        bad_y    = np.array(np.where(mlmasky>1e15,0,1))
        #bad_x    = np.array(np.where(mlmaskx>1e29,0,1))
        #bad_y    = np.array(np.where(mlmasky>1e29,0,1))

        mlmaskx  = mlmaskx[np.logical_and(bad_x,bad_y)]
        mlmaskx  = np.reshape(mlmaskx,[mlmaskx.size,1])

        mlmasky  = mlmasky[np.logical_and(bad_x,bad_y)]
        mlmasky  = np.reshape(mlmasky,[mlmasky.size,1])

        mlmaskxy = np.concatenate((mlmaskx, mlmasky),axis=1)

        if mlmaskxy.size>0:
            m_patches.append( Polygon(np.array(mlmaskxy), 
                    closed=True,\
                    clip_on=True,\
                    capstyle='projecting') )

            del mlmaskx; del mlmasky; del mlmaskxy

            ax.add_collection(PatchCollection(    m_patches,        \
                            facecolor=fcolor,    \
                            edgecolor=ecolor,    \
                            linewidths=0.5,     \
                            zorder=1))
        
    
    
#___PLOT 3D FESOM GLOBE IN SPHERICAL COORDINATES________________________________
#
#
#_______________________________________________________________________________
def fesom_plot3d_earth(mesh,data):
    #___________________________________________________________________________
    from mpl_toolkits.mplot3d import Axes3D
    from mayavi import mlab
    from matplotlib.path import Path
    from tvtk.api import tvtk
    
    #___________________________________________________________________________
    zz     = np.abs(mesh.nodes_2d_z)
    zmax   = np.max(zz)
    zz[mesh.nodes_2d_i==1]=0
    pfac   = 0.075
    R      = 1-zz/zmax*pfac
    cart_x,cart_y,cart_z = geo2spherical(R, mesh.nodes_2d_xg, mesh.nodes_2d_yg)
    
    #___________________________________________________________________________
    tri      = Triangulation(mesh.nodes_2d_xg, mesh.nodes_2d_yg,mesh.elem_2d_i)
    
    #___________________________________________________________________________
    # mayavi 
    fig = mlab.figure(size=(1100, 1100), bgcolor=(0.25, 0.25, 0.25))
    fig.scene.interactor.interactor_style = tvtk.InteractorStyleImage()
    #fig.scene.interactor.interactor_style = tvtk.InteractorStyleTerrain()
    
    #___________________________________________________________________________
    cnumb=20
    cmax = np.max(data.value[np.isnan(data.value)==False])
    cmin = np.min(data.value[np.isnan(data.value)==False])
    cref = 0.0
    cmap0,clevel = colormap_c2c(cmin,cmax,cref,cnumb,data.cmap)
    print(clevel)
    
    data.value[data.value<clevel[0]]  = clevel[0]+np.finfo(np.float32).eps
    data.value[data.value>clevel[-1]] = clevel[-1]-np.finfo(np.float32).eps
    #___________________________________________________________________________
    # plot 3d spherical globe
    globe=mlab.triangular_mesh(cart_x, \
                 cart_y, \
             cart_z, \
             tri.triangles,\
             scalars=data.value,\
             colormap='blue-red')
    
    #___________________________________________________________________________
    # set colorbar
    cb = mlab.colorbar(title='Depth [m]', orientation='vertical', nb_colors=clevel.size-1,\
              label_fmt='%4.2f')
    cb.scalar_bar_representation.position = [0.9,0.2]
    cb.scalar_bar_representation.position2 = [0.05,0.7]
    cb.label_text_property.font_size=10
    
    globe.module_manager.scalar_lut_manager.reverse_lut=True
    #globe.module_manager.scalar_lut_manager.number_of_colors=512
    
    #___________________________________________________________________________
    # plot 3d spherical land sea mask
    xy_all,tri_all=calc_3d_lsmask(mesh)
    cart_x,cart_y,cart_z = geo2spherical(1.0, xy_all[:,0], xy_all[:,1])
    mlab.triangular_mesh(cart_x, \
             cart_y, \
             cart_z, \
             tri_all,\
             color=(0.8,0.8,0.8) )
               
    mlab = plot_3d_lmaskline(mlab,mesh)
    #___________________________________________________________________________
    # plot zonal & meridional 3d gridlines
    aux_lon  = np.arange(-180,180.5,0.5)
    aux_lat  = np.arange(-80,80.5,0.5)
    
    grid_step= 30
    grid_lon = np.arange(-180,180+grid_step,grid_step)
    grid_lat = np.arange(-90,90+grid_step,grid_step)
    calc_3d_lon_lat_gridline(mlab,grid_lon,grid_lat,aux_lon,aux_lat)

    # plot thicker 0 +-90 +-180 meridians and equator & zonal&meridional 
    # thicklines and labels
    fig.scene.disable_render = True # Super duper trick
    grid_lon_step = grid_lon
    grid_lat_step = grid_lat
    grid_lon      = np.array([0,180])
    grid_lat      = np.array([-80,0,80])
    calc_3d_lon_lat_maingrid(mlab,grid_lon,grid_lat,aux_lon,aux_lat,grid_lon_step,grid_lat_step)
    fig.scene.disable_render = False # Super duper trick
    
    #___________________________________________________________________________
    globe.scene.light_manager.light_mode = "vtk"

    return(globe)
    #return(globe,cb)
    #fig.scene.disable_render = False # Super duper trick
    mlab.show()
    
    
#___TRANSFORM 2D GEO COORDINATES TO 3D SPHERICAL CARTESIAN COORD._______________
#
#_______________________________________________________________________________
def geo2spherical(R,lon,lat):
    cart_x = R * np.cos(np.radians(lat)) * np.cos(np.radians(lon)) 
    cart_y = R * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
    cart_z = R * np.sin(np.radians(lat))
    return(cart_x,cart_y,cart_z) 
    
    
#___CALC. SPHERICAL 3d LAND SEA MASK____________________________________________
#
#_______________________________________________________________________________
def calc_3d_lsmask(mesh):
    
    from matplotlib.path import Path
    
    #___________________________________________________________________________
    npoly = np.shape(mesh.polygon_xy)
    for ii in range(0,npoly[0]):
        #_______________________________________________________________________
        x_min = np.floor(np.min(mesh.polygon_xy[ii][:,0]))  ; x_max = np.ceil(np.max(mesh.polygon_xy[ii][:,0])) ; 
        y_min = np.floor(np.min(mesh.polygon_xy[ii][:,1])) ; y_max = np.ceil(np.max( mesh.polygon_xy[ii][:,1])) ; 
        x_vec = np.arange(x_min,x_max,5)
        y_vec = np.arange(y_min,y_max,5)
        x_m,y_m = np.meshgrid(x_vec,y_vec);
        del x_vec; del y_vec
        x_m   = x_m.reshape((x_m.size,1));
        y_m   = y_m.reshape((y_m.size,1));

        #_______________________________________________________________________
        # check if regular points are within polygon
        IN    = Path(mesh.polygon_xy[ii]).contains_points(np.concatenate((x_m,y_m),axis=1))  # False
        x_m   = x_m[IN==True]
        y_m   = y_m[IN==True]
        del IN
        xy_m = np.concatenate((mesh.polygon_xy[ii],np.concatenate((x_m,y_m),axis=1)))
        del x_m; del y_m

        #_______________________________________________________________________
        # make delaunay triangulation with polygon points and additional points 
        # that are inside the polygon
        tri      = Triangulation(xy_m[:,0], xy_m[:,1])

        tri_x = np.sum(xy_m[tri.triangles,0],axis=1)/3
        tri_y = np.sum(xy_m[tri.triangles,1],axis=1)/3
        tri_x = np.reshape(tri_x,(tri_x.size,1))
        tri_y = np.reshape(tri_y,(tri_y.size,1))
        IN    = Path(mesh.polygon_xy[ii]).contains_points(np.concatenate((tri_x,tri_y),axis=1))  # False
        tri.triangles=tri.triangles[IN==True,:]

        #_______________________________________________________________________
        if ii==0:
            xy_all  = xy_m
            tri_all = tri.triangles
        else:
            tri_all = np.concatenate((tri_all,tri.triangles+xy_all.shape[0]),axis=0)
            xy_all  = np.concatenate((xy_all,xy_m),axis=0)
    
    #_______________________________________________________________________
    return(xy_all,tri_all)
    
    
#___CALC. LON/LAT 3D GRID LINES_________________________________________________
#
#_______________________________________________________________________________
def calc_3d_lon_lat_gridline(mlab,grid_lon,grid_lat,aux_lon,aux_lat):
    
    x = list()
    y = list()
    z = list()
    s = list()
    connections = list()
    
    # The index of the current point in the total amount of points
    index = 0
    N = aux_lat.shape[0]
    # Create each line one after the other in a loop
    for ii in range(0,grid_lon.shape[0]):
        cart_x,cart_y,cart_z = geo2spherical(1.0, aux_lat*0+grid_lon[ii], aux_lat)

        x.append(cart_x)
        y.append(cart_y)
        z.append(cart_z)
        s.append(aux_lat)
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
                    [np.arange(index,   index + N - 1.5),
                    np.arange(index + 1, index + N - .5)]
                    ).T)
        index += N

        N = aux_lon.shape[0]
    # Create each line one after the other in a loop
    for ii in range(0,grid_lat.shape[0]):
        cart_x,cart_y,cart_z = geo2spherical(1.0, aux_lon, aux_lon*0+grid_lat[ii])

        x.append(cart_x)
        y.append(cart_y)
        z.append(cart_z)
        s.append(aux_lon)
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
                    [np.arange(index,   index + N - 1.5),
                    np.arange(index + 1, index + N - .5)]
                    ).T)
        index += N


        # Now collapse all positions, scalars and connections in big arrays
        x = np.hstack(x)
        y = np.hstack(y)
        z = np.hstack(z)
        s = np.hstack(s)

        connections = np.vstack(connections)

        # Create the points
        src = mlab.pipeline.scalar_scatter(x, y, z, s)

        # Connect them
        src.mlab_source.dataset.lines = connections
        src.update()

        # The stripper filter cleans up connected lines
        lines = mlab.pipeline.stripper(src)

        # Finally, display the set of lines
        mlab.pipeline.surface(lines, color=(0,0,0), line_width=2, opacity=.4)
    
    return lines
    
    
#___CALC. LON/LAT MAIN THICK GRID LINES + THEIR TICKS___________________________
#
#_______________________________________________________________________________
def calc_3d_lon_lat_maingrid(mlab,grid_lon,grid_lat,aux_lon,aux_lat,grid_lon_step,grid_lat_step):
    
    x = list()
    y = list()
    z = list()
    s = list()
    connections = list()
    
    # The index of the current point in the total amount of points
    index = 0
    N = aux_lat.shape[0]
    # Create each line one after the other in a loop
    for ii in range(0,grid_lon.shape[0]):
        cart_x,cart_y,cart_z = geo2spherical(1.0, aux_lat*0+grid_lon[ii], aux_lat)

        x.append(cart_x)
        y.append(cart_y)
        z.append(cart_z)
        s.append(aux_lat)
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
                    [np.arange(index,   index + N - 1.5),
                    np.arange(index + 1, index + N - .5)]
                    ).T)
        index += N

        N = aux_lon.shape[0]
    # Create each line one after the other in a loop
    for ii in range(0,grid_lat.shape[0]):
        cart_x,cart_y,cart_z = geo2spherical(1.0, aux_lon, aux_lon*0+grid_lat[ii])

        x.append(cart_x)
        y.append(cart_y)
        z.append(cart_z)
        s.append(aux_lon)
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
                    [np.arange(index,   index + N - 1.5),
                    np.arange(index + 1, index + N - .5)]
                    ).T)
        index += N
    
    #___________________________________________________________________________
    # plot tick line
    N = 2
    ticklength = 0.05;
    for ii in range(0,grid_lon.shape[0]):
        for jj in range(0,grid_lat_step.shape[0]):
            cart_x,cart_y,cart_z = geo2spherical(np.array([1.0, 1.0+ticklength]), \
                                  np.repeat(grid_lon[ii],2), \
                                  np.repeat(grid_lat_step[jj],2))

            x.append(cart_x) ; y.append(cart_y) ; z.append(cart_z) ; s.append(cart_z)
            connections.append(np.vstack([np.arange(index,   index + N - 1.5),np.arange(index + 1, index + N - .5)]).T)
            index += N
        
        for ii in range(0,grid_lon_step.shape[0]):
            cart_x,cart_y,cart_z = geo2spherical(np.array([1.0, 1.0+ticklength]), \
                                  np.repeat(grid_lon_step[ii],2), \
                                  np.repeat(0.0,2))

            x.append(cart_x) ; y.append(cart_y) ; z.append(cart_z) ; s.append(cart_z)
            connections.append(np.vstack([np.arange(index,   index + N - 1.5),np.arange(index + 1, index + N - .5)]).T)
            index += N
    
    #___________________________________________________________________________
    # Now collapse all positions, scalars and connections in big arrays
    x = np.hstack(x)
    y = np.hstack(y)
    z = np.hstack(z)
    s = np.hstack(s)
    connections = np.vstack(connections)
    
    # Create the points
    src = mlab.pipeline.scalar_scatter(x, y, z, s)
    
    # Connect them
    src.mlab_source.dataset.lines = connections
    src.update()
    
    # The stripper filter cleans up connected lines
    lines = mlab.pipeline.stripper(src)
    
    # Finally, display the set of lines
    mlab.pipeline.surface(lines, color=(0,0,0), line_width=2, opacity=0.75)
    
    #___________________________________________________________________________
    # plot ticklabels
    for ii in range(0,grid_lon.shape[0]):
        for jj in range(0,grid_lat_step.shape[0]):
            cart_x,cart_y,cart_z = geo2spherical(1.0+ticklength+0.02, \
                                  grid_lon[ii], \
                                  grid_lat_step[jj])
            t=mlab.text3d(cart_x, cart_y, cart_z, \
                      np.str(grid_lat_step[jj]),\
                color=(0,0,0),\
                scale=0.04,opacity=0.75)
            t.vector_text.update()
            #mlab.text(cart_x, cart_y, 
                      #np.str(grid_lat_step[jj]),\
                #z=cart_z, \
                #width=0.02)
        for jj in range(0,grid_lon_step.shape[0]):
            cart_x,cart_y,cart_z = geo2spherical(1.0+ticklength+0.02, \
                                  grid_lon_step[jj], \
                                  0)
            t=mlab.text3d(cart_x, cart_y, cart_z, \
                      np.str(grid_lon_step[jj]),\
                    color=(0,0,0),\
                    scale=0.04,opacity=0.75)        
            t.vector_text.update()    
            #mlab.text(cart_x, cart_y, \
                      #np.str(grid_lon_step[jj]),\
                    #z=cart_z, \
                    #width=0.02)        
    
    return mlab
    
    
#___CALC. LON/LAT 3D GRID LINES_________________________________________________
#
#_______________________________________________________________________________
def plot_3d_lmaskline(mlab,mesh):
    
    x = list()
    y = list()
    z = list()
    s = list()
    connections = list()
    
    # The index of the current point in the total amount of points
    index = 0
    
    npoly = np.shape(mesh.polygon_xy)
    for ii in range(0,npoly[0]):
        N = mesh.polygon_xy[ii].shape[0]
        cart_x,cart_y,cart_z = geo2spherical(1.0, mesh.polygon_xy[ii][:,0], mesh.polygon_xy[ii][:,1])

        x.append(cart_x)
        y.append(cart_y)
        z.append(cart_z)
        s.append(cart_z)
        # This is the tricky part: in a line, each point is connected
        # to the one following it. We have to express this with the indices
        # of the final set of points once all lines have been combined
        # together, this is why we need to keep track of the total number of
        # points already created (index)
        connections.append(np.vstack(
                    [np.arange(index,   index + N - 1.5),
                    np.arange(index + 1, index + N - .5)]
                    ).T)
        index += N
    
    # Now collapse all positions, scalars and connections in big arrays
    x = np.hstack(x)
    y = np.hstack(y)
    z = np.hstack(z)
    s = np.hstack(s)
    connections = np.vstack(connections)
    
    # Create the points
    src = mlab.pipeline.scalar_scatter(x, y, z, s)
    
    # Connect them
    src.mlab_source.dataset.lines = connections
    
    # The stripper filter cleans up connected lines
    lines = mlab.pipeline.stripper(src)
    
    # Finally, display the set of lines
    mlab.pipeline.surface(lines, color=(0,0,0), line_width=1, opacity=.5)
    
    return mlab


#_______________________________________________________________________________
#
#_______________________________________________________________________________
def fesom_idxinbox(mesh,data1,inputarray):
    #___________________________________________________________________________
    # calculate index of which data points are within boi
    idxbox_e = mesh.nodes_2d_xg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][0]
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_xg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][1])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].max(axis=1)<inputarray['which_box'][2])
    idxbox_e = np.logical_or(idxbox_e,mesh.nodes_2d_yg[mesh.elem_2d_i].min(axis=1)>inputarray['which_box'][3])
    idxbox_e = idxbox_e==False # true index for triangles that are within box 
    
    # case of node data
    if mesh.n2dna in data1.value.shape:
        idxbox_n  = mesh.elem_2d_i[idxbox_e,:].flatten().transpose()
        idxbox_n  = np.array(idxbox_n).squeeze()
        idxbox_n  = np.unique(idxbox_n)
        idx_box   = np.zeros((mesh.n2dna,), dtype='bool')
        idx_box[idxbox_n]=True
        del idxbox_n 
        del idxbox_e
    # case of element data
    elif mesh.n2dea in data1.value.shape:
        idx_box = idxbox_e
        del idxbox_e
    
    #___________________________________________________________________________
    # make triangle mask arrays from flat triangles and nan trinagles
    if np.any(np.isnan(data1.value)):
        idx_box = np.logical_and(idx_box,np.isnan(data1.value)==False)
    
    return(idx_box)


#_______________________________________________________________________________
# select optimal color range by histogramm
#_______________________________________________________________________________
def fesom_choose_best_crange(in_data,in_weights,limit=0.99,fac=1.0,do_output=False,increace_dezimal=0):
    cmin, cmax = np.nanmin(in_data), np.nanmax(in_data)
    if do_output:print(' --> orig cmin,cmax:',cmin,cmax)
    
    in_data = in_data*fac
    
    binrange=[cmin, cmax];
    if cmin<0.0 and cmax>0.0 : binrange=[-max(np.abs(cmin),cmax),max(np.abs(cmin),cmax)]
    if len(in_weights)!=0:
        hist, binedge = np.histogram(in_data,range=(binrange[0], binrange[1]), bins=10000, weights=in_weights,density=True, normed=True) 
    else:
        hist, binedge = np.histogram(in_data,range=(binrange[0], binrange[1]), bins=10000,density=True, normed=True) 
        
    binedge_mid, cumsum, limit = binedge[:-1]+(binedge[1:]-binedge[:-1])/2, np.cumsum(hist*(binedge[1:]-binedge[:-1])), limit
    
    idx_min = np.where(cumsum<=1-limit)[0]
    if len(idx_min)==0: 
        idx_min=0
    else: 
        idx_min=idx_min[-1]
    
    idx_max = np.where(cumsum>=limit)[0]
    if len(idx_max)==0: 
        idx_max=len(binedge_mid)
    else:    
        idx_max=idx_max[0]
    cmin, cmax = binedge_mid[idx_min], binedge_mid[idx_max]
    
    if cmax!=0.0: cmax = np.around(cmax, -np.int32(np.floor(np.log10(np.abs(cmax)))-1+increace_dezimal) ) 
    if cmin!=0.0: cmin = np.around(cmin, -np.int32(np.floor(np.log10(np.abs(cmin)))-1+increace_dezimal) ) 
    
    if do_output:print(' --> best cmin,cmax:',cmin,cmax)
    
    return(cmin,cmax)

#_______________________________________________________________________________
# select optimal color range by histogramm
#_______________________________________________________________________________
def fesom_choose_best_cref(cmin,cmax,varname,do_rescale='auto',fac=0):
    cref = cmin + (cmax-cmin)/2
    if cref!=0.0 : cref = np.around(cref, -np.int32(np.floor(np.log10(np.abs(cref)))-1+fac)) 
    if varname in ['u','v','w','ssh','fw','fh'] or \
       any(x in varname for x in ['vec','anom','dvd']):
       if not do_rescale=='log10': cref=0.0
    if varname in ['Kv']: cref = np.around(cref,0)
    return(cref)    
