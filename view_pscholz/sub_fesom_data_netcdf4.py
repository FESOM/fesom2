# Patrick Scholz, 23.01.2018
import numpy as np
import numpy.matlib
import time
import os
from netCDF4 import Dataset, MFDataset
from set_inputarray import *
from sub_fesom_mesh import *
import matplotlib.pyplot as plt
import seawater as sw
global inputarray
    
#+_____________________________________________________________________________+
#|                                                                             |
#|                        *** FESOM DATA CLASS ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_data(object):
    
    which_obj                   = 'data'
    
    #____data run variables______________________
    var                         = ''
    runid, path, descript       = 'fesom', '', ''
    
    #____data time variables_____________________
    year,  month, record, depth = [], [], [], []
    str_time, str_dep           = '', ''
    
    #____data projection variables_______________
    proj, proj_lon, proj_lat    = '', 0.0, 90.0
    
    #____data description info___________________
    sname, lname, unit            = '', '', ''
    
    #____data plot varaibles_____________________
    cmap, crange, cnumb         = 'grads', [], []
    which_plot                  = ''
    
    #____data varaible___________________________
    value, value2, value3, time = [], [], [], []
    rescale                     = []
    which_mean                  = 'monthly'
    anom                        = False
    
    #___INIT DATA OBJECT________________________________________________________
    #
    #___________________________________________________________________________
    def __init__(self,inputarray):
        if len(inputarray)!=0:
            self.runid       = inputarray['data_id']
            self.path        = inputarray['data_dir1']
            self.proj        = inputarray['proj']
            self.proj_lon    = inputarray['proj_lon']
            self.proj_lat    = inputarray['proj_lat']
            self.which_plot  = inputarray['which_plot']
    
#___LOAD FESOM2.0 DATA AS TIME AVERAGED HORIZONTAL SLICE________________________
#
#_______________________________________________________________________________
def fesom_load_data_horiz_netcdf4(mesh,data,             \
                                  which_files='default', \
                                  do_tmean=True,         \
                                  do_loadloop=False,     \
                                  do_rescale='auto',     \
                                  do_interp_e2n=True,    \
                                  do_vecrot=True,        \
                                  do_output=True):
    #_______________________________________________________________________________
    # plot depth
    if data.var=='depth':
        data.value 	= -mesh.nodes_2d_zg
        data.sname, data.lname, data.unit, data.cmap = 'depth', 'Depth', 'm', 'wbgyr'
        return data
    #_______________________________________________________________________________
    # plot triangle resolution interpolated to node
    elif data.var=='triresol':
        if len(mesh.nodes_2d_resol)==0: mesh.fesom_calc_triresol()
        data.value 	= mesh.nodes_2d_resol
        data.sname, data.lname, data.unit, data.cmap = 'triresol', 'Resolution', 'km', 'rygbw'
        return data
    #_______________________________________________________________________________
    # plot triangle area interpolated to node
    elif data.var=='triarea':
        if len(mesh.nodes_2d_area)==0: mesh.fesom_calc_triarea()
        data.value 	= mesh.nodes_2d_area
        data.sname,data.lname, data.unit, data.cmap= 'triarea', 'Area', 'km^2', 'cmocean.cm.balance'
        return data
    
    #___________________________________________________________________________
    # number of years to average 
    nyi         = data.year[1]-data.year[0]+1
    # array with years
    ayi         = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    # number of month to average
    nmi         = np.array(data.month).size
    
    #___________________________________________________________________________
    # create multi-year file list to load 
    fname_list,\
    var_list    = do_multiyear_fname_list(data, which_files, do_output)
    
    #___________________________________________________________________________
    # load first one to get the dimensions and attributes
    fname_data  = Dataset(fname_list[0][0],'r')
     
    #___________________________________________________________________________
    # compute dimension contained in file
    nti,nsi,ndi = do_filedims(fname_data, nyi, do_output)
    
    #___________________________________________________________________________
    # get attributes of variable 
    data        = do_fileattr(fname_data, data, var_list)
    
    #___________________________________________________________________________
    # do pre-selection of depth levels in each file thats opened with open_mfdataset
    sel_levidx  = do_select_levidx(mesh, data, ndi, do_output)
        
    #___________________________________________________________________________
    # do pre-selection of time steps in each file thats opened with open_mfdataset
    sel_timeidx = do_select_timeidx(data, nti, nyi, nmi, do_loadloop, do_output)
    
    #___________________________________________________________________________
    # select time and depth range and apply time averagng
    if not do_loadloop:
        data    = do_load_mfdata(  mesh, data, fname_list, var_list, sel_timeidx, \
                                   sel_levidx, nsi, ndi, do_tmean, do_output,)
    else:    
        data    = do_load_dataloop(mesh, data, fname_list, var_list, sel_timeidx, \
                                   sel_levidx, nti, nsi, ndi, do_tmean, do_output,)
        
    #___________________________________________________________________________
    # do nessessary postprocesssing to rotate vectors, compute norm, set NaNs ...
    data        = do_postprocess(mesh, data, do_rescale, do_interp_e2n, do_vecrot, do_output)
    
    #___________________________________________________________________________
    # add augmented periodic boundary
    data        = do_pbnd_augmentation(mesh, data, nti, ndi, do_tmean)
    
    #___________________________________________________________________________
    # in case data object coms from FESOm1.4 than is loaded with pyfesom and it
    # is mesh.zlevs instead of mesh.zlev
    if len(data.depth)==0:
        try:
            mesh.zlevs
        except AttributeError:
            data.zlev = mesh.zlev
        else:     
            data.zlev = mesh.zlevs
    else:
        data.zlev = data.depth
    
    #___________________________________________________________________________
    # set up descriptiv variable names
    data        = do_descriptstr(data, nyi, nmi, nti, ndi)
    if do_output: print('')
    return(data)

#___COMBINE MULTI-YEAR FILE LOADING LISTS_______________________________________
#
#_______________________________________________________________________________
def do_multiyear_fname_list(data, which_files, do_output):
    
    # number of years to average 
    nyi = data.year[1]-data.year[0]+1
    # array with years
    ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    
    # initialise fname_list
    fname_list = [[],[],[]]
    if do_output==True:
        print(' --> path          :',data.path)
        print(' --> year/month/rec:',data.year,',',data.month,',',data.record)
        print(' --> depth         :',data.depth)
        print(' --> var           :',data.var, end='')
    # fill fname_list
    for yi in range(0, nyi):
        
        #_______________________________________________________________________
        # build filename list for data variables
        if any(x in data.var for x in ['norm','vec','ptemp','pdens','sigma']) or \
                    data.var in ['u','v','uice','vice'] :
            if   any(x in data.var for x in ['tuv']):
                var_list = ['u', 'v', 'temp']
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
                fname_list[2].append(data.path+'/'+do_fname_mask(which_files,var_list[2],data.runid,str(ayi[yi])))
                
            elif any(x in data.var for x in ['suv']):
                var_list = ['u', 'v', 'salt']
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
                fname_list[2].append(data.path+'/'+do_fname_mask(which_files,var_list[2],data.runid,str(ayi[yi])))
            
            elif any(x in data.var for x in ['uice','vice']):
                var_list = ['uice', 'vice',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['f_uwind','f_vwind','f_uvwind']):
                var_list = ['f_uwind', 'f_vwind',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['uwind','vwind','uvwind']):
                var_list = ['uwind', 'vwind',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['f_tx_sur','f_ty_sur','f_txy_sur']):
                var_list = ['f_tx_sur', 'f_ty_sur',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['tx_sur','ty_sur','txy_sur']):
                var_list = ['tx_sur', 'ty_sur',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['pdens','ptemp','sigma']):   
                var_list = ['temp', 'salt',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            elif any(x in data.var for x in ['pgf_x','pgf_y','pgf_xy','pgf']):
                var_list = ['pgf_x', 'pgf_y',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))    
            elif any(x in data.var for x in ['uv','u','v']):
                var_list = ['u', 'v',[]]
                fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
                fname_list[1].append(data.path+'/'+do_fname_mask(which_files,var_list[1],data.runid,str(ayi[yi])))
            else:
                print(' --> vector variable not defined')
            
        else:
            var_list = [data.var ,[],[]]
            fname_list[0].append(data.path+'/'+do_fname_mask(which_files,var_list[0],data.runid,str(ayi[yi])))
    
    if do_output==True: print(' --> var_list:',var_list)
    
    return(fname_list, var_list)

#___FILE NAME MASK______________________________________________________________
#
#_______________________________________________________________________________
def do_fname_mask(which_files,varname,runid,year):
    if   which_files=='default':
        fname = varname+'.'+runid+'.'+year+'.nc'
    elif which_files=='blowup_oce':
        fname = 'fesom.'+year+'.'+'oce'+'.'+'blowup'+'.nc'
    elif which_files=='blowup_ice':    
        fname = 'fesom.'+year+'.'+'ice'+'.'+'blowup'+'.nc'
    elif which_files=='restart_oce':
        fname = 'fesom.'+year+'.'+'oce'+'.'+'restart'+'.nc'
    elif which_files=='restart_ice':
        fname = 'fesom.'+year+'.'+'ice'+'.'+'restart'+'.nc'
        
    return(fname)    

#___COMPUTE DIMENSIONS CONTAINED IN FILES_______________________________________
#
#_______________________________________________________________________________
def do_filedims(fname_data,nyi,do_output):
    dimname     = list(fname_data.dimensions.keys())
    nti,nsi,ndi = 0,0,0
    for aux in dimname:
        #if aux=='time'  : nti = int(len(fname_data.dimensions[aux])/nyi)
        if aux=='time'  : nti = len(fname_data.dimensions[aux])
        if aux=='nod2'  : nsi = len(fname_data.dimensions[aux])
        if aux=='elem'  : nsi = len(fname_data.dimensions[aux])
        if aux=='nz'    : ndi = len(fname_data.dimensions[aux])
        if aux=='nz1'   : ndi = len(fname_data.dimensions[aux])
    del aux    
    if do_output==True: print(' --> nti/nsi/ndi   : [{:d}, {:d}, {:d}]'.format(nti,nsi,ndi))  
    
    return(nti,nsi,ndi)
    
#___LOAD MULTI-FILE DATA VIA MFDATASET__________________________________________
#
#_______________________________________________________________________________
def do_load_mfdata(mesh, data, fname_list, var_list, sel_timeidx, sel_levidx, \
                   nsi, ndi, do_tmean, do_output,):    
    if ndi!=0:
        
        #_______________________________________________________________________
        # select time+depth range + compute time mean
        if do_tmean:
            data.value = MFDataset(fname_list[0],'r').variables[var_list[0]][sel_timeidx,:,sel_levidx].mean(axis=0)
            if len(fname_list[1]): data.value2 = MFDataset(fname_list[1],'r').variables[var_list[1]][sel_timeidx,:,sel_levidx].mean(axis=0)
            if len(fname_list[2]): data.value3 = MFDataset(fname_list[2],'r').variables[var_list[2]][sel_timeidx,:,sel_levidx].mean(axis=0)
        else:
            data.value = MFDataset(fname_list[0],'r').variables[var_list[0]][sel_timeidx,:,sel_levidx]
            if len(fname_list[1]): data.value2 = MFDataset(fname_list[1],'r').variables[var_list[1]][sel_timeidx,:,sel_levidx]
            if len(fname_list[2]): data.value3 = MFDataset(fname_list[2],'r').variables[var_list[2]][sel_timeidx,:,sel_levidx]
        
        #_______________________________________________________________________
        # compute potential density & temperatur if selected
        if any(x in data.var for x in ['pdens','ptemp','sigma']):
            dep   = np.matlib.repmat(mesh.zmid[sel_levidx],nsi,1)
            lat   = np.matlib.repmat(mesh.nodes_2d_yg[0:mesh.n2dn],len(sel_levidx),1).transpose()
            press = sw.pres(dep,lat)
            press_ref = 0
            if   '0' in data.var : press_ref=0
            elif '1' in data.var : press_ref=1000
            elif '2' in data.var : press_ref=2000
            elif '3' in data.var : press_ref=3000
            elif '4' in data.var : press_ref=4000
            elif '5' in data.var : press_ref=5000
            if 'sigma' in data.var: data.lname = '$\sigma_{'+str(int(press_ref/1000))+'}$ '+data.lname
            del dep,lat
            if do_tmean:
                if 'ptemp' in data.var: data.value = sw.ptmp(data.value2,data.value,press,press_ref)
                if any(x in data.var for x in ['pdens','sigma']): data.value = sw.pden(data.value2,data.value,press,press_ref)-1000.025 
            else:
                for it in range(0,data.value.shape[0]):
                    if 'ptemp' in data.var: data.value[it,:,:] = sw.ptmp(data.value2[it,:,:],data.value[it,:,:],press,press_ref)
                    if any(x in data.var for x in ['pdens','sigma']): data.value[it,:,:] = sw.pden(data.value2[it,:,:],data.value[it,:,:],press,press_ref)-1000.025 
            fname_list[1]=[]
            
        #_______________________________________________________________________    
        # compute depth mean + linear interpolation to selected depth levels
        if do_tmean:
            data.value = do_zinterp(mesh, data.value, data.depth, ndi, data.value.shape[0], sel_levidx,do_output)
            if len(fname_list[1]): data.value2 = do_zinterp(mesh, data.value2, data.depth, ndi, data.value2.shape[0], sel_levidx,do_output)
            if len(fname_list[2]): data.value3 = do_zinterp(mesh, data.value3, data.depth, ndi, data.value3.shape[0], sel_levidx,do_output)
        else:
            data.value = do_zinterp(mesh, data.value, data.depth, ndi, data.value.shape[1], sel_levidx,do_output)
            if len(fname_list[1]): data.value2 = do_zinterp(mesh, data.value2, data.depth, ndi, data.value2.shape[1], sel_levidx,do_output)
            if len(fname_list[2]): data.value3 = do_zinterp(mesh, data.value3, data.depth, ndi, data.value3.shape[1], sel_levidx,do_output)
        
    # 2D data:
    else: 
        
        if do_tmean:
            data.value = MFDataset(fname_list[0],'r').variables[var_list[0]][sel_timeidx,:].mean(axis=0)
            if len(fname_list[1]): data.value2 = MFDataset(fname_list[1],'r').variables[var_list[1]][sel_timeidx,:].mean(axis=0)
            if len(fname_list[2]): data.value3 = MFDataset(fname_list[2],'r').variables[var_list[2]][sel_timeidx,:].mean(axis=0)
        else:
            data.value = MFDataset(fname_list[0],'r').variables[var_list[0]][sel_timeidx,:]
            if len(fname_list[1]): data.value2 = MFDataset(fname_list[1],'r').variables[var_list[1]][sel_timeidx,:]
            if len(fname_list[2]): data.value3 = MFDataset(fname_list[2],'r').variables[var_list[2]][sel_timeidx,:]
    
    # kickout single array dimension
    data.value = data.value.squeeze()
    if len(fname_list[1]): data.value2 = data.value2.squeeze()
    if len(fname_list[2]): data.value3 = data.value3.squeeze()
    
    return(data) 

#___LOAD MULTI-FILE DATA VIA MFDATASET__________________________________________
#
#_______________________________________________________________________________
def do_load_dataloop(mesh, data, fname_list, var_list, sel_timeidx, sel_levidx, \
                   nti, nsi, ndi, do_tmean, do_output,):    
    if ndi!=0:
        
        #_______________________________________________________________________
        # initialize data.value array
        if do_tmean:
            data.value = np.zeros((nsi,len(sel_levidx)),dtype='float32')
            if len(fname_list[1]): data.value2 = np.zeros((nsi,len(sel_levidx)),dtype='float32')
            if len(fname_list[1]): data.value3 = np.zeros((nsi,len(sel_levidx)),dtype='float32')
        else:
            data.value = np.zeros((nti,nsi,len(sel_levidx)))
            if len(fname_list[1]): data.value2 = np.zeros((nti*len(fname_list),nsi,len(sel_levidx)),dtype='float32')
            if len(fname_list[1]): data.value3 = np.zeros((nti*len(fname_list),nsi,len(sel_levidx)),dtype='float32')
        
        #_______________________________________________________________________
        # select time+depth range + compute time mean
        for it in range(0,len(fname_list[0])):
            if do_tmean:
                data.value = data.value + Dataset(fname_list[0][it],'r').variables[var_list[0]][sel_timeidx,:,sel_levidx].mean(axis=0)
                if len(fname_list[1]): data.value2 = data.value2 + Dataset(fname_list[1][it],'r').variables[var_list[1]][sel_timeidx,:,sel_levidx].mean(axis=0)
                if len(fname_list[2]): data.value3 = data.value3 + Dataset(fname_list[2][it],'r').variables[var_list[2]][sel_timeidx,:,sel_levidx].mean(axis=0)
            else:
                t_idx = sel_timeidx+nti*it
                data.value[t_idx,:,:] = data.value[t_idx,:,:] + Dataset(fname_list[0][it],'r').variables[var_list[0]][sel_timeidx,:,sel_levidx]
                if len(fname_list[1]): data.value2[t_idx,:,:] = data.value2[t_idx,:,:] + Dataset(fname_list[1][it],'r').variables[var_list[1]][sel_timeidx,:,sel_levidx]
                if len(fname_list[2]): data.value3[t_idx,:,:] = data.value3[t_idx,:,:] + Dataset(fname_list[2][it],'r').variables[var_list[2]][sel_timeidx,:,sel_levidx]
        
        # divide by loaded file number --> to final time mean
        if do_tmean:
            data.value = data.value/len(fname_list[0]) 
            if len(fname_list[1]): data.value2 = data.value2/len(fname_list[1]) 
            if len(fname_list[2]): data.value3 = data.value2/len(fname_list[2])
                
        #_______________________________________________________________________
        # compute potential density & temperatur if selected
        if any(x in data.var for x in ['pdens','ptemp','sigma']):
            dep      = np.matlib.repmat(mesh.zmid[sel_levidx],nsi,1)
            lat      = np.matlib.repmat(mesh.nodes_2d_yg[0:mesh.n2dn],nsi,1).transpose()
            press    = sw.pres(dep,lat)
            press_ref= 0
            if   '0' in data.var : press_ref=0
            elif '1' in data.var : press_ref=1000
            elif '2' in data.var : press_ref=2000
            elif '3' in data.var : press_ref=3000
            elif '4' in data.var : press_ref=4000
            elif '5' in data.var : press_ref=5000
            if press_ref!=0: data.lname = '$\sigma_{'+str(int(press_ref/1000))+'}$ '+data.lname
            del dep,lat
            if do_tmean:
                if 'ptemp' in data.var: data.value = sw.ptmp(data.value2,data.value,press,press_ref)
                if any(x in data.var for x in ['pdens','sigma']): data.value = sw.pden(data.value2,data.value,press,press_ref)-1000.025 
            else:
                for it in range(0,data.value.shape[0]):
                    if 'ptemp' in data.var: data.value[it,:,:] = sw.ptmp(data.value2[it,:,:],data.value[it,:,:],press,press_ref)
                    if any(x in data.var for x in ['pdens','sigma']): data.value[it,:,:] = sw.pden(data.value2[it,:,:],data.value[it,:,:],press,press_ref)-1000.025 
            fname_list[1]=[]
            
        #_______________________________________________________________________    
        # compute depth mean + linear interpolation to selected depth levels
        data.value = do_zinterp(mesh, data.value, data.depth, ndi, nsi, sel_levidx,do_output)
        if len(fname_list[1]): data.value2 = do_zinterp(mesh, data.value2, data.depth, ndi, nsi, sel_levidx,do_output)
        if len(fname_list[2]): data.value3 = do_zinterp(mesh, data.value3, data.depth, ndi, nsi, sel_levidx,do_output)
        
    # 2D data:
    else: 
        #_______________________________________________________________________
        # initialize data.value array
        if do_tmean:
            data.value = np.zeros((nsi,))
            if len(fname_list[1]): data.value2 = np.zeros((nsi,))
            if len(fname_list[1]): data.value3 = np.zeros((nsi,))
        else:
            data.value = np.zeros((nti,nsi))
            if len(fname_list[1]): data.value2 = np.zeros((nti*len(fname_list),nsi))
            if len(fname_list[1]): data.value3 = np.zeros((nti*len(fname_list),nsi))
        
        #_______________________________________________________________________
        # select time+depth range + compute time mean
        for it in range(0,len(fname_list[0])):
            if do_tmean:
                data.value = data.value + Dataset(fname_list[0][it],'r').variables[var_list[0]][sel_timeidx,:].mean(axis=0)
                if len(fname_list[1]): data.value2 = data.value2 + Dataset(fname_list[1][it],'r').variables[var_list[1]][sel_timeidx,:].mean(axis=0)
                if len(fname_list[2]): data.value3 = data.value3 + Dataset(fname_list[2][it],'r').variables[var_list[2]][sel_timeidx,:].mean(axis=0)
            else:
                t_idx = sel_timeidx+nti*it
                data.value[t_idx,:] = data.value[t_idx,:] + Dataset(fname_list[0][it],'r').variables[var_list[0]][sel_timeidx,:]
                if len(fname_list[1]): data.value2[t_idx,:] = data.value2[t_idx,:] + Dataset(fname_list[1][it],'r').variables[var_list[1]][sel_timeidx,:]
                if len(fname_list[2]): data.value3[t_idx,:] = data.value3[t_idx,:] + Dataset(fname_list[2][it],'r').variables[var_list[2]][sel_timeidx,:]
                
        # divide by loaded file number --> to final time mean
        if do_tmean:
            data.value = data.value/len(fname_list[0]) 
            if len(fname_list[1]): data.value2 = data.value2/len(fname_list[1]) 
            if len(fname_list[2]): data.value3 = data.value2/len(fname_list[2])
    
    return(data) 

#___WRITE FILE ATTRIBUTES TO DATA ARRAY_________________________________________
#
#_______________________________________________________________________________
def do_fileattr(fname_data,data,var_list):
    data.sname = data.var 
    for auxname, auxvariable in fname_data.variables.items(): 
        if auxname==var_list[0]:
            for attrname in auxvariable.ncattrs():
                if attrname=='description' : data.lname = getattr(auxvariable, attrname)
                if attrname=='units'       : data.unit  = '['+getattr(auxvariable, attrname)+']'
    if any(x in data.var for x in ['ptemp']): data.lname='potential temperature'
    if any(x in data.var for x in ['pdens','sigma']): data.lname,data.unit='potential density','[kg/m^3]'
    if any(x in data.var for x in ['tuv'  ]): data.lname,data.unit='temperature advection', '[°C*m/s]'
    if any(x in data.var for x in ['suv'  ]): data.lname,data.unit='salt advection', '[psu*m/s]'  
    if any(x in data.var for x in ['norm_uv'  ]): data.lname='norm of horizontal velocity'
    if any(x in data.var for x in ['norm_tuv','norm_suv']): data.lname='norm of horizontal '+data.lname
    if data.var == 'v' : data.lname='meridional velocity'
    if data.var == 'u' : data.lname='zonal velocity'
    return(data)

#___DO VERTICAL INTERPOLTION____________________________________________________
#
#_______________________________________________________________________________
def do_zinterp(mesh, idata, idepth, ndi, nsi, sel_levidx,do_output):
    # only do vertical interpolation when something is written in the data.depth
    if len(idepth)!=0:
        weight   = np.zeros((len(idepth),3))
        # div ... divisor for correct vertical averaging 
        div      = np.zeros((nsi,)) 
        # valid_lay ... is interpolated layer valid [1.0 yes, 0.0 no ] 
        valid_lay= np.zeros((nsi,len(idepth))) 
        
        ndimax = mesh.nodes_2d_izg.max()-1
        if ndi==mesh.nlev: ndimax = mesh.nodes_2d_izg.max()
        
        for zi in range(0,len(idepth)):
            if ndi!=mesh.nlev: 
                auxidx = np.searchsorted(-mesh.zmid,idepth[zi])   
            else: 
                auxidx = np.searchsorted(-mesh.zlev,idepth[zi])
            
            if auxidx>ndimax: auxidx = ndimax
            if auxidx>0: auxidx = auxidx-1
            
            weight[zi,0], weight[zi,1] = sel_levidx.index(auxidx), sel_levidx.index(auxidx+1)
            if ndi!=mesh.nlev: 
                depth = mesh.zmid
                if   nsi==mesh.n2dn: isvalid = mesh.nodes_2d_izg[:nsi]-1>=auxidx+1
                elif nsi==mesh.n2de: isvalid = mesh.elem0_2d_iz[:nsi]-1>=auxidx+1
                valid_lay[isvalid,zi] = valid_lay[isvalid,zi] + 1.0
            else             :
                depth = mesh.zlev
                if   nsi==mesh.n2dn: isvalid = mesh.nodes_2d_izg[:nsi]>=auxidx+1
                elif nsi==mesh.n2de: isvalid = mesh.elem0_2d_iz[:nsi]>=auxidx+1
                valid_lay[isvalid,zi] = valid_lay[isvalid,zi] + 1.0
            weight[zi,2] = (idepth[zi]+depth[auxidx])/(depth[auxidx]-depth[auxidx+1])    
            div[isvalid] = div[isvalid]+1.0
            del isvalid    
            
            if weight[zi,2]<0.0:
                if do_output: print(' --> WARNING !!!   : no vert. extrap. supported! reset '+str(idepth[zi])+'m to nearest valid layer '+str(abs(depth[auxidx]))+'m')
                weight[zi,2] = max([weight[zi,2],0.0])
            if weight[zi,2]>1.0:
                if do_output: print(' --> WARNING !!!   : no vert. extrap. supported! reset '+str(idepth[zi])+'m to nearest valid layer '+str(abs(depth[auxidx+1]))+'m')
                weight[zi,2] = min([weight[zi,2],1.0])    
        
        #_______________________________________________________________________
        # compute depth interpolation
        if idata.ndim==2:
            idata = (idata[:,weight[:,0].astype('int')]         \
                    + (                                           \
                        idata[:,weight[:,1].astype('int')]\
                        -idata[:,weight[:,0].astype('int')]\
                        )*weight[:,2])*valid_lay
            idata = idata.sum(axis=1)
            
            #___________________________________________________________________
            # divide by number of vlid layers --> to compute mean depth
            isvalid                        = div!=0
            idata[isvalid]                 = idata[isvalid]/div[isvalid]
            # set nodes where are no valid layers for interpolation to nan 
            idata[np.logical_not(isvalid)] = np.nan
        
        # compute depth interpolation in case no time mean was applied    
        elif idata.ndim==3:
            nt,ns,nd=idata.shape
            idata2 = np.zeros((nt,ns,len(idepth)))
            for it in range(0,nt):
                idata2[it,:,:] =  (np.transpose(idata[it,:,weight[:,0].astype('int')])         \
                                + np.transpose(                                           \
                                    idata[it,:,weight[:,1].astype('int')]\
                                    -idata[it,:,weight[:,0].astype('int')]\
                                )*weight[:,2])*valid_lay
                        
            idata = idata2.sum(axis=2)
            #_______________________________________________________________________
            # divide by number of vlid layers --> to compute mean depth
            isvalid = div!=0
            for it in range(0,nt):
                idata[it,isvalid]                 = idata[it,isvalid]/div[isvalid]
                # set nodes where are no valid layers for interpolation to nan 
                idata[it,np.logical_not(isvalid)] = np.nan
    
    return(idata)

#___DO NECCESSARY POSTPROCESSING________________________________________________
#
#_______________________________________________________________________________
def do_postprocess(mesh, data, do_rescale, do_interp_e2n, do_vecrot, do_output):

    # in case more than one variable has been loaded compute norm, vec, ...
    if any(x in data.var for x in ['tuv','suv']):
        # interpolate data3 value to triangle center like u and v
        data.value3 = data.value3[mesh.elem0_2d_i].mean(axis=1)
        data.value, data.value2, data.value3 = data.value*data.value3, data.value2*data.value3, []

    ## do vector rotation from ROT --> GEO in case of vector data   
    if do_vecrot==True:
        if any(x in data.var for x in ['norm','vec']) or \
        data.var in ['u','v','uice','vice']:
            data.value, data.value2 = fesom_vector_rot(mesh, data.value, data.value2,do_output=do_output)
        
    # compute norm of vector data
    if all(x in data.var for x in ['norm']):
        data.value, data.value2 = np.sqrt(data.value**2+data.value2**2),[]
    
    # only single velocity component is choosen
    if data.var in ['u','uice'] or any(x in data.var for x in ['pdens','sigma']): data.value2 = []
    if data.var in ['v','vice']: data.value, data.value2  = data.value2, []
    
    # interpolate elemental values to nodes
    if 'vec' not in data.var and do_interp_e2n==True:
        if data.value.ndim==1 and data.value.shape[0]==mesh.n2de:
            data.value = mesh.fesom_interp_e2n(np.array(data.value))
            
        elif data.value.ndim>1 and (data.value.shape[0]==mesh.n2de or data.value.shape[1]==mesh.n2dea):
            data.value = mesh.fesom_interp_e2n(np.array(data.value))
            
    # customise sea ice variables
    if data.var in ['a_ice']        : data.value = data.value * 100.0
    #if data.var in ['a_ice','m_ice']: data.value[data.value<=0.0]=np.nan
    
    # change sign of mixed layer depth from negativ to positve
    if 'mld' in data.var.lower(): data.value = -data.value
    
    # cutoff exponentials --> add therefore string to unit parameter
    if do_rescale=='auto':
        if np.nanmax(np.abs(data.value))<1e-2 and np.nanmax(np.abs(data.value))>0.0:
            scal = 10**(np.floor(np.log10(max(abs(np.nanmin(data.value)),abs(np.nanmax(data.value))))-1))
            data.value = data.value/scal
            data.rescale=scal
            if any(x in data.var for x in ['vec']): data.value2 = data.value2/scal
            data.unit  = ' $ \cdot 10^{'+str(int(np.log10(scal)))+'} $'+data.unit
    elif do_rescale=='log10':
        data.value[data.value!=0.0] = np.log10(data.value[data.value!=0.0])
        data.rescale='log10'
        if any(x in data.var for x in ['vec']): data.value2 = np.log10(data.value2)
        data.unit  = ' log10() '+data.unit
        
    return(data)

#___AUGMENT PERIODIC BOUNDARY VALUES____________________________________________
#
#_______________________________________________________________________________
def do_pbnd_augmentation(mesh, data, nti, ndi, do_tmean):
    if len(data.depth)!=0 or ndi==0:
        if do_tmean or nti==1: 
            if data.value.shape[0]==mesh.n2dn: 
                data.value = np.concatenate((data.value,data.value[mesh.pbndn_2d_i]))
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[mesh.pbndn_2d_i]))
            if data.value.shape[0]==mesh.n2de: 
                data.value = np.concatenate((data.value,data.value[mesh.pbndtri_2d_i]))
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[mesh.pbndtri_2d_i]))
        else:
            if data.value.shape[1]==mesh.n2dn: 
                data.value = np.concatenate((data.value,data.value[:,mesh.pbndn_2d_i]),axis=1)
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[:,mesh.pbndn_2d_i]),axis=1)
            if data.value.shape[1]==mesh.n2de: 
                data.value = np.concatenate((data.value,data.value[:,mesh.pbndtri_2d_i]),axis=1)
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[:,mesh.pbndtri_2d_i]),axis=1)
    else:
        if do_tmean or nti==1:
            if data.value.shape[0]==mesh.n2dn: 
                data.value = np.concatenate((data.value,data.value[mesh.pbndn_2d_i,:]))
                if any(x in data.var for x in ['vec']):
                    data.value2 = np.concatenate((data.value2,data.value2[mesh.pbndn_2d_i,:]))
            if data.value.shape[0]==mesh.n2de: 
                data.value = np.concatenate((data.value,data.value[mesh.pbndtri_2d_i,:]))
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[mesh.pbndtri_2d_i,:]))    
        else:
            if data.value.shape[1]==mesh.n2dn: 
                data.value = np.concatenate((data.value,data.value[:,mesh.pbndn_2d_i,:]),axis=1)
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[:,mesh.pbndn_2d_i,:]),axis=1)
            if data.value.shape[1]==mesh.n2de: 
                data.value = np.concatenate((data.value,data.value[:,mesh.pbndtri_2d_i,:]), axis=1)
                if any(x in data.var for x in ['vec']): 
                    data.value2 = np.concatenate((data.value2,data.value2[:,mesh.pbndtri_2d_i,:]), axis=1)    
    return(data)

#___DO DESCRIPTIVE TIME STRINGS_____________________________________________________
#
#_______________________________________________________________________________
def do_descriptstr(data,nyi,nmi,nti,ndi):
    # array with month
    ami = np.array(data.month)
    
    if len(data.record)==0:
        mon_list   = np.array(['J','F','M','A','M','J','J','A','S','O','N','D'])
        mon_list2  = np.array(['January','February','March','April','May','June','July','August','September','October','Novemver','December'])
        # year info 
        if nyi == 1:
            str_time_1 = 'y: '+str(data.year[0])
        else:    
            str_time_1 = 'y: '+str(data.year[0])+'-'+str(data.year[1])
            
        # month info 
        if nmi == 1 and nti>=12:
            str_time_1 = str_time_1+' ['+mon_list2[np.array(data.month)-1][0]+']'
        elif nmi>1 and nmi<12 and nti>=12:
            str_mon=''
            for ii in range(0,nmi):
                if ami[ii]==12 : str_mon=mon_list[ami[ii]-1]+str_mon
                else            : str_mon=str_mon+mon_list[ami[ii]-1]
            str_time_1 = str_time_1 + ' [' + str_mon + ']'
        data.str_time = str_time_1
    else: 
        data.str_time = 'rec_i: '+str(data.record[0])
    
    if ndi!=0 and len(data.zlev)!=0:
        if len(data.zlev)<=1:
            data.str_dep=', dep: '+str(data.zlev[0])+'m'
        else:
            data.str_dep=', dep: '+str(data.zlev[0])+'m-'+str(data.zlev[-1])+'m'
    return(data)
    
#___PRE SELECT INDICES FOR VERTICAL LEVEL SUBSET________________________________
#
#_______________________________________________________________________________
def do_select_levidx(mesh,data,ndi,do_output):
    sel_levidx = []
    if ndi!=0 : # deal with 3d data
        # select all depth layers
        if len(data.depth) == 0 :    
            sel_levidx = list(range(0,ndi))
        # select single depth level --> still need interpolation so two depth 
        # indices have to be selected   
        elif len(data.depth) == 1 :  
            # variable is on full depth levels
            if ndi==mesh.nlev:
                #ndimax     = mesh.zlev.size-1
                ndimax     = mesh.nodes_2d_izg.max()
                auxidx     = np.searchsorted(-mesh.zlev,data.depth[0])
                
            # variable is on mid-depth levels    
            else:
                #ndimax     = mesh.zmid.size-1
                ndimax     = mesh.nodes_2d_izg.max()-1
                auxidx     = np.searchsorted(-mesh.zmid,data.depth[0])
               
            if   auxidx>ndimax : sel_levidx = [ndimax-1,ndimax]       
            elif auxidx>=1     : sel_levidx = [auxidx-1,auxidx]
            else               : sel_levidx = [auxidx,auxidx+1]    
            
            del auxidx
        
        # select a range for depth levels that have to be interpolated
        elif len(data.depth) >= 2:
            sel_levidx=[]
            
            ndimax = mesh.nodes_2d_izg.max()-1
            if ndi==mesh.nlev: ndimax = mesh.nodes_2d_izg.max()
            
            for aux in data.depth:
                if ndi==mesh.nlev:
                    auxidx     = np.searchsorted(-mesh.zlev,aux)
                else:
                    auxidx     = np.searchsorted(-mesh.zmid,aux)
                if auxidx>ndimax and ndimax not in sel_levidx: sel_levidx.append(ndimax)    
                if auxidx>=1 and auxidx-1 not in sel_levidx: sel_levidx.append(auxidx-1)
                if (auxidx not in sel_levidx): sel_levidx.append(auxidx)
                if (auxidx==0 and 1 not in sel_levidx): sel_levidx.append(auxidx+1)
                
    if do_output==True: print(' --> sel_levidx    :',sel_levidx)            
    
    return(sel_levidx)
    
#___PRE SELECT INDICES FOR TEMPORAL SUBSET______________________________________
#
#_______________________________________________________________________________
def do_select_timeidx(data ,nti, nyi, nmi, do_loadloop, do_output):
    sel_timeidx=[]
    if len(data.record)!=0:
        if data.record[0]>nti: 
            print(' --> WARNING --> data.record indices to large, correct it to maximum value in the file')
            data.record=[1]
        sel_timeidx = [x-1 for x in data.record]
    else:
        
        # file contains annual data
        if   nti==1:
            sel_timeidx = [0]
            if nmi!=12: print(' --> WARNING --> no monthly information in file --> load annual instead')
            
        # file contains monthly data
        elif nti==12:                   
            sel_timeidx = [x-1 for x in data.month]
            
        # file contains 5 daily data    
        elif nti==73:                   
            fivedayinyear=np.arange(5, 365+1,5)
            daypermonth = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
            idx_5day2mon = np.zeros((73,))
            for ii in range(0,12): 
                idx_5day2mon[(fivedayinyear>daypermonth[ii]) & (fivedayinyear<=daypermonth[ii+1])]=ii+1
            del fivedayinyear, daypermonth
            # find out which 5day index belongs to selected monthly subset
            sel_timeidx=[i for i, e in enumerate(idx_5day2mon) if any(e==data.month) ]
            del idx_5day2mon  
            
        # file contains daily data    
        elif nti==365 or nti==366:
            if   nti==365: daypermonth, idx_day = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]), list(range(1,365+1))
            elif nti==366: daypermonth, idx_day = np.cumsum([0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]), list(range(1,366+1))
            # find out which day belongs to which month
            idx_day2mon = np.zeros((365,))
            for ii in range(0,12): 
                idx_day2mon[(idx_day>daypermonth[ii]) & (idx_day<=daypermonth[ii+1])]=ii+1
            del idx_day, daypermonth
            # find out which day index belongs to selected monthly subset
            sel_timeidx=[i for i, e in enumerate(idx_day2mon) if any(e==data.month) ]
            del idx_day2mon  

        else:
            print(' --> time frame not supported')
            #return
            sel_timeidx=list(range(0,nti))
            
    # expand time selction indices in case several years are loaded to the entire loaded file 
    if nyi>1 and not do_loadloop:
        orig_timeidx = np.copy(sel_timeidx)
        if nyi>1:
            for yi in range(1, nyi): sel_timeidx.extend([x+(nti*yi) for x in orig_timeidx])
        del orig_timeidx    
    
    if do_output==True: print(' --> sel_timeidx   :',sel_timeidx)
    
    return(sel_timeidx)

#___CALCULATE ANOMALY BETWEEN TWO DATA OBJECTS__________________________________
#
#_______________________________________________________________________________
def fesom_data_anom(data,data2):
    anom = fesom_data([])
    
    anom.which_obj = 'data'
    
    #____data run variables______________________
    if data.var!=data2.var:
        anom.var                        = data.var+'-'+data2.var
    else:
        anom.var                        = data.var
    anom.runid, anom.path               = data.runid, data.path
    anom.descript                       = data.descript+'-'+data2.descript
    
    #____data time variables_____________________
    anom.year,  anom.month, anom.record = data.year, data.month, data.record
    anom.depth                          = data.depth
    if data.str_time!=data2.str_time:
        anom.str_time = data.str_time+'-'+data2.str_time
    else:
        anom.str_time = data.str_time
    anom.str_dep                        = data.str_dep
    
    #____data projection variables_______________
    anom.proj                           = data.proj
    anom.proj_lon, anom.proj_lat        = data.proj_lon, data.proj_lat
    anom.cmap,anom.cnumb                = data.cmap,data.cnumb  
    
    #____data description info___________________
    #anom.sname, anom.lname, anom.unit   = data.sname, data.lname, data.unit
    if data.sname!=data2.sname:
        anom.sname = data.sname+'-'+data2.sname
    else:
        anom.sname = data.sname
    if data.lname!=data2.lname:
        anom.lname = data.lname+'-'+data2.lname
    else:
        anom.lname = data.lname
    anom.unit   = data.unit
    #____data plot varaibles_____________________
    #anom.cmap, anom.crange, anom.cnumb  = data.cmap, data.crange, data.cnumb
    anom.which_plot                     = data.which_plot 
    #if len(data.crange)!=0 :anom.crange=[]
    #____data varaible___________________________
    anom.value                          = data.value-data2.value
    if len(data.value2)!=0:
        anom.value2                     = data.value2-data2.value2
    if len(data.value3)!=0:
        anom.value3                     = data.value3-data2.value3
    anom.time                           = data.time
    anom.anom                           = True
    anom.value2                         = data.value2
    anom.which_mean                     = data.which_mean
    
    #___________________________________________________________________________
    return(anom)

#___COPY DATA OBJECTS___________________________________________________________
#
#_______________________________________________________________________________
def fesom_data_copy(data):
    copy = fesom_data([])
    
    #____data run variables______________________
    copy.var                            = data.var
    copy.runid, copy.path               = data.runid, data.path
    copy.descript                       = data.descript
    
    #____data time variables_____________________
    copy.year,  copy.month, copy.record = data.year, data.month, data.record
    copy.depth                          = data.depth
    copy.str_time                         = data.str_time
    copy.str_dep                        = data.str_dep
    
    #____data projection variables_______________
    copy.proj                           = data.proj
    copy.proj_lon, copy.proj_lat        = data.proj_lon, data.proj_lat
    
    #____data description info___________________
    copy.sname, copy.lname, copy.unit   = data.sname, data.lname, data.unit
    
    #____data plot varaibles_____________________
    #anom.cmap, anom.crange, anom.cnumb  = data.cmap, data.crange, data.cnumb
    copy.which_plot                     = data.which_plot 
    #if len(data.crange)!=0 :anom.crange=[]
    #____data varaible___________________________
    copy.value                          = data.value
    copy.time                           = data.time
    copy.anom                           = data.anom  
    copy.value2                         = data.value2
    copy.which_mean                     = data.which_mean
    
    #___________________________________________________________________________
    return(copy)
