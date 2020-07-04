# Patrick Scholz, 23.01.2018
import numpy as np
import time
import os
from netCDF4 import Dataset
from set_inputarray import *
from sub_fesom_mesh import *
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy.matlib
import seawater as sw
#from sub_fesom_mesh import fesom_vector_rot
global inputarray
    
#+_____________________________________________________________________________+
#|                                                                             |
#|                        *** FESOM DATA CLASS ***                             |
#|                                                                             |
#+_____________________________________________________________________________+
class fesom_data(object):
    
    which_obj                    = 'data'
    
    #____data run variables______________________
    var                         = ''
    runid, path, descript       = 'fesom', '', ''
    
    #____data time variables_____________________
    year,  month, record, depth,zlev = [], [], [], [],[]
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
    which_mean                  = 'monthly'
    anom                        = False
    
    
    #___INIT DATA OBJECT_______________________________________________________#
    #
    #__________________________________________________________________________#
    def __init__(self,inputarray):
        if len(inputarray)!=0:
            self.runid        = inputarray['data_id']
            self.path        = inputarray['data_dir1']
            self.proj        = inputarray['proj']
            self.proj_lon    = inputarray['proj_lon']
            self.proj_lat    = inputarray['proj_lat']
            self.which_plot    = inputarray['which_plot']
    
    
#___LOAD FESOM2.0 DATA AS TIME AVERAGED HORIZONTAL SLICE________________________
#
#_______________________________________________________________________________
def fesom_load_data_horiz(mesh,data,do_output=True):
    
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
    nyi = data.year[1]-data.year[0]+1
    
    # array with years
    ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    # number of month to average
    nmi = np.array(data.month).size
    # array with month
    ami = np.array(data.month)
    
    if do_output==True:
        print('')
        print('     -----+-----------------------------------+------------')
        print('     Year |               MON                 |')
        print('     -----+-----------------------------------+------------')
        print('     --> '+data.path)
        print('     --> '+data.var)
    #____START YEAR LOOP________________________________________________________
    aux_datavar = data.var
    for yi in range(0, nyi):
        if do_output==True: print('     {:4.0f} |'.format(ayi[yi]),end='')
        tstart = time.time();
        if data.var.find('norm')!=-1 or data.var.find('vec')!=-1 or data.var=='ptemp' or data.var=='pdens':
            if data.var.find('tuv')!=-1 or data.var.find('suv')!=-1:
                
                if   data.var.find('tuv')!=-1 : aux_datavar,aux_datavar2,aux_datavar3 = 'u','v','temp'
                elif data.var.find('suv')!=-1 : aux_datavar,aux_datavar2,aux_datavar3 = 'u','v','salt'
                #_______________________________________________________________________
                # build filename where data are located
                fname  =   data.path+'/'\
                            + aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                fname2 =   data.path+'/'\
                            + aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                fname3 =   data.path+'/'\
                            + aux_datavar3 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'        
                #_______________________________________________________________________
                # open netcdf file --> read NETCDF4
                ncid  = Dataset(fname , 'r') 
                ncid2 = Dataset(fname2, 'r') 
                ncid3 = Dataset(fname3, 'r')
                
            elif data.var.find('uice')!=-1 or data.var.find('vice')!=-1:    
                aux_datavar,aux_datavar2 = 'uice','vice'
                #_______________________________________________________________________
                # build filename where data are located
                fname  =   data.path+'/'\
                            + aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                fname2 =   data.path+'/'\
                            + aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                #_______________________________________________________________________
                # open netcdf file --> read NETCDF4
                ncid  = Dataset(fname, 'r') 
                ncid2 = Dataset(fname2, 'r') 
                
            else:
                aux_datavar,aux_datavar2 = 'u','v'
                if data.var =='ptemp' or data.var=='pdens':
                    aux_datavar,aux_datavar2 = 'temp','salt'
                #_______________________________________________________________________
                # build filename where data are located
                fname  =   data.path+'/'\
                            + aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                fname2 =   data.path+'/'\
                            + aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
                #_______________________________________________________________________
                # open netcdf file --> read NETCDF4
                ncid  = Dataset(fname, 'r') 
                ncid2 = Dataset(fname2, 'r') 
        else:
            #_______________________________________________________________________
            # build filename where data are located
            fname =   data.path+'/'\
                    + data.var + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            
            #_______________________________________________________________________
            # open netcdf file --> read NETCDF4
            ncid = Dataset(fname, 'r') 
        
        
        # read dimension of variable 
        ncdims=ncid.variables[aux_datavar].shape
        #print(ncdims)
        
        # number of dimension is 2d or 3d
        dim_num = len(ncdims)
        
        # check time dimension of loaded data aka number of records
        nrec  = ncdims[0]
        nsmple= ncdims[1]
        
        #_______________________________________________________________________
        # at first loaded year allocate matrices to calc means
        if yi==0:
            
            # if data are 3d and no loadable depth layers are given load all fesom 
            # layers
            if dim_num==3 and  (len(data.depth)==0 or np.all(data.depth==mesh.zlev)==True):
                # allocate array for mean data
                mean_data = np.zeros((nsmple,ncdims[2]),dtype='float16')
                #mean_data = np.zeros((nsmple,ncdims[2]),dtype='float32')
                if data.var.find('vec')!=-1:
                    ncdims2    = ncid2.variables[aux_datavar2].shape
                    nsmple2    = ncdims2[1]
                    mean_data2 = np.zeros((nsmple2,ncdims2[2]),dtype='float16')
                    #mean_data2 = np.zeros((nsmple2,ncdims2[2]),dtype='float32')
                if data.var.find('tuv')!=-1 or data.var.find('suv')!=-1:
                    ncdims3    = ncid3.variables[aux_datavar3].shape
                    nsmple3    = ncdims3[1]
                    mean_data3 = np.zeros((nsmple3,ncdims2[2]),dtype='float16')    
                    #mean_data3 = np.zeros((nsmple3,ncdims2[2]),dtype='float32')    
            
            else:
                # allocate array for mean data
                mean_data = np.zeros((nsmple,),dtype='float16')
                #mean_data = np.zeros((nsmple,),dtype='float32')
                if data.var.find('vec')!=-1:
                    nsmple2    = ncid2.variables[aux_datavar2].shape[1]
                    mean_data2 = np.zeros((nsmple2,),dtype='float16')
                    #mean_data2 = np.zeros((nsmple2,),dtype='float32')
                if data.var.find('tuv')!=-1 or data.var.find('suv')!=-1:
                    nsmple3    = ncid3.variables[aux_datavar3].shape[1]
                    mean_data3 = np.zeros((nsmple3,),dtype='float16')
                    #mean_data3 = np.zeros((nsmple3,),dtype='float32')
            #___________________________________________________________________
            # setup variable name, description and unit
            ncvar = ncid.variables[aux_datavar]
            #data.sname=ncvar.name
            data.sname=data.var
            for attrname in ncvar.ncattrs():
                if attrname=='description':
                    data.lname=ncvar.getncattr(attrname)
                elif attrname=='units':
                    data.unit='['+ncvar.getncattr(attrname)+']'
            if data.var =='ptemp': data.lname='potential temperature'
            if data.var =='pdens': data.lname,data.unit='potential density','[kg/m^3]'
            if data.var.find('tuv')!=-1: data.lname='temperature advection'
            if data.var.find('suv')!=-1: data.lname='salt advection'
            
            #___________________________________________________________________
            # calc 3d depth and latitude for calculation of potential density 
            # and temperature
            if data.var =='ptemp' or data.var=='pdens':
                depth3d = np.zeros(ncid.variables[aux_datavar][:].shape)
                lat3d   = np.zeros(ncid.variables[aux_datavar][:].shape)
                for nreci in range(0,nrec):
                    depth3d[nreci,:,:] = np.matlib.repmat((mesh.zlev[0:-1]+(mesh.zlev[1::]-mesh.zlev[0:-1])/2),nsmple,1)
                    lat3d[nreci,:,:] = np.matlib.repmat(mesh.nodes_2d_yg[0:mesh.n2dn],ncid.variables[aux_datavar][:].shape[2],1).transpose()
                press3d = sw.pres(depth3d,lat3d)
                del depth3d, lat3d
            
        
        #_______________________________________________________________________
        #ncval = np.array(ncid.variables[aux_datavar][:])
        ncval = np.copy(ncid.variables[aux_datavar][:])
        if data.var.find('norm')!=-1 or data.var.find('vec')!=-1:
            ncval2 = np.copy(ncid2.variables[aux_datavar2][:])
            
        if data.var.find('tuv')!=-1 or data.var.find('tuv')!=-1:
            ncval3 = np.copy(ncid3.variables[aux_datavar3][:])
            if data.var.find('norm')!=-1:
                ncval3 = ncval3[:,mesh.elem0_2d_i,:]
                ncval3 = ncval3.mean(axis=2)
                ncval  = ncval*ncval3
                ncval2 = ncval2*ncval3
                del ncval3
                ncid3.close()
            
        if data.var.find('norm')!=-1:
            ncval = np.sqrt(ncval[:]**2+ncval2[:]**2)
            del ncval2
            ncid2.close()
            
        elif data.var=='ptemp':    
            #ncval2 = np.array(ncid2.variables[aux_datavar2][:])
            ncval2 = np.copy(ncid2.variables[aux_datavar2][:])
            ncval = sw.ptmp(ncval2,ncval,press3d)
            del ncval2
            
        elif data.var=='pdens':    
            #ncval2 = np.array(ncid2.variables[aux_datavar2][:])
            ncval2 = np.copy(ncid2.variables[aux_datavar2][:])
            ncval = sw.pden(ncval2,ncval,press3d,0)-1000.025 
            del ncval2
            
        #_______________________________________________________________________
        # check if a single time slice record is selcted
        # no ?  --> than use data.month to select time subset
        if len(data.record)==0:
            #___________________________________________________________________
            if nrec==1:
                if do_output==True: print('annual |',end='')
            elif nrec==12:
                if do_output==True: print('monthly|',end='')
            #elif nrec==73:
                #print('5daily |'),
            elif nrec==365:
                if do_output==True: print('daily  |',end='')
            else:    
                print(' --> error: temporal length of data unclear or not supported!!!|'),
                break
            #___________________________________________________________________
            # select monthly or seasonal subset of data
            if nmi<12 and nrec>=12:
                if do_output==True: 
                    for ii in data.month: print('{:d}|'.format(ii),end='')
                #_______________________________________________________________
                # select from monthly data
                if nrec==12 : 
                    ncval=fesom_time_depth_mean(mesh,ncval,np.array(data.month)-1,data.depth)
                    if data.var.find('vec')!=-1: 
                        ncval2=fesom_time_depth_mean(mesh,ncval2,np.array(data.month)-1,data.depth)
                    if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:
                        ncval3=fesom_time_depth_mean(mesh,ncval3,np.array(data.month)-1,data.depth)
                #___________________________________________________________________
                # select from daily data
                elif nrec==365 : 
                    # find out which day belongs to selected subset
                    idx_sel = sel_timesubset_daily(data)
                    ncval=fesom_time_depth_mean(mesh,ncval,idx_sel==True,data.depth)
                    if data.var.find('vec')!=-1: 
                        ncval2=fesom_time_depth_mean(mesh,ncval2,idx_sel==True,data.depth)
                    if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:    
                        ncval3=fesom_time_depth_mean(mesh,ncval3,idx_sel==True,data.depth)
                    del idx_sel
                    
            #_______________________________________________________________________
            # select full data
            else:    
                if nrec>=12: 
                    if do_output==True: 
                        for ii in data.month: print('{:d}|'.format(ii)),
                ncval=fesom_time_depth_mean(mesh,ncval,[],data.depth)
                if data.var.find('vec')!=-1: 
                    ncval2=fesom_time_depth_mean(mesh,ncval2,[],data.depth)
                if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:    
                    ncval3=fesom_time_depth_mean(mesh,ncval3,[],data.depth)
        #_______________________________________________________________________
        # yes ?  --> use data.record to select time slice
        else:
            if np.max(data.record)-1>nrec:
                print(' --> error: this record number doesnt exist in that file')
                break
            print('select single record time slice|',end='')
            for ii in data.record: print('{:d}|'.format(ii)),
            ncval=fesom_time_depth_mean(mesh,ncval,data.record[0]-1,data.depth)
            if data.var.find('vec')!=-1: 
                ncval2=fesom_time_depth_mean(mesh,ncval2,data.record[0]-1,data.depth)
            if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:    
                ncval3=fesom_time_depth_mean(mesh,ncval3,data.record[0]-1,data.depth)
        #_______________________________________________________________________
        # close netcdf file
        ncid.close()
        if data.var.find('vec')!=-1: ncid2.close()
        if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:ncid3.close()
        #_______________________________________________________________________
        # now average over yearly subset
        mean_data = mean_data + ncval/nyi
        del ncval 
        if data.var.find('vec')!=-1: 
            mean_data2 = mean_data2 + ncval2/nyi
            del ncval2 
        if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:    
            mean_data3 = mean_data3 + ncval3/nyi
            del ncval3 
        
        #_______________________________________________________________________
        tend = time.time();        
        if do_output==True: print(' --> t={:2.2f}s'.format(tend-tstart))
    #___________________________________________________________________________
    # in case of ice data set no ice to nan
    if (data.var=='a_ice'):
        mean_data=mean_data*100.0
    if (data.var=='a_ice' or data.var=='m_ice'):
        #mean_data[mean_data==0.0]=np.nan
        mean_data[mean_data<=0.0]=np.nan
    
    #____ROTATE VECTORES FROM ROT TO GEO________________________________________
    # in future check if rotation is s till neccessary -> we want to/should write out 
    # the data already rotated ?!
    if data.var.find('vec')==0:
        mean_data,mean_data2=fesom_vector_rot(mesh,mean_data,mean_data2)
        
    #___________________________________________________________________________
    # augment periodic boundary of data.value
    if nsmple==mesh.n2dn:
        # augment node data
        data.value = np.concatenate((mean_data,mean_data[mesh.pbndn_2d_i]))
        del mean_data
    else:
        # augment elem data
        data.value = np.concatenate((mean_data,mean_data[mesh.pbndtri_2d_i]))
        del mean_data
        # delete trinangle that are seen as abnormal 
        if len(mesh.abnormtri_2d_i)!=0:
            data.value = np.delete(data.value,mesh.abnormtri_2d_i)
    
    #___________________________________________________________________________
    # augment periodic boundary of data.value2
    if data.var.find('vec')==0: 
        if nsmple2==mesh.n2dn:
        # augment node data
            data.value2 = np.concatenate((mean_data2,mean_data2[mesh.pbndn_2d_i]))
            del mean_data2
        else:
            # augment elem data
            data.value2 = np.concatenate((mean_data2,mean_data2[mesh.pbndtri_2d_i]))
            del mean_data2
            # delete trinangle that are seen as abnormal 
            if len(mesh.abnormtri_2d_i)!=0:
                data.value2 = np.delete(data.value2,mesh.abnormtri_2d_i)
    
    #___________________________________________________________________________
    # augment periodic boundary of data.value2
    if data.var.find('vec_tuv')!=-1 or data.var.find('vec_suv')!=-1:    
        if nsmple3==mesh.n2dn:
        # augment node data
            data.value3 = np.concatenate((mean_data3,mean_data3[mesh.pbndn_2d_i]))
            del mean_data3
        else:
            # augment elem data
            data.value3 = np.concatenate((mean_data3,mean_data3[mesh.pbndtri_2d_i]))
            del mean_data3
            # delete trinangle that are seen as abnormal 
            if len(mesh.abnormtri_2d_i)!=0:
                data.value3 = np.delete(data.value3,mesh.abnormtri_2d_i)
    
    #___________________________________________________________________________
    if len(data.depth)==0:
        try:
            mesh.zlevs
        except AttributeError:
            data.zlev = mesh.zlev
        else:     
            data.zlev = mesh.zlevs
        #data.depth = mesh.zlev
    else:
        data.zlev = data.depth
    
    #___________________________________________________________________________
    # set up descriptiv variable names
    if len(data.record)==0:
        mon_list   = np.array(['J','F','M','A','M','J','J','A','S','O','N','D'])
        mon_list2  = np.array(['January','February','March','April','May','June','July','August','September','October','Novemver','December'])
        # year info 
        if nyi == 1:
            str_time_1 = 'y: '+str(data.year[0])
        else:    
            str_time_1 = 'y: '+str(data.year[0])+'-'+str(data.year[1])
            
        # month info 
        if nmi == 1 and nrec>=12:
            str_time_1 = str_time_1+' ['+mon_list2[np.array(data.month)-1][0]+']'
        elif nmi>1 and nmi<12 and nrec>=12:
            str_mon=''
            for ii in range(0,nmi):
                if ami[ii]==12 : str_mon=mon_list[ami[ii]-1]+str_mon
                else            : str_mon=str_mon+mon_list[ami[ii]-1]
            str_time_1 = str_time_1 + ' [' + str_mon + ']'
        data.str_time = str_time_1
    else: 
        data.str_time = 'rec_i: '+str(data.record[0])
    
    if dim_num==3 and len(data.zlev)!=0:
        if len(data.zlev)<=1:
            data.str_dep=', dep: '+str(data.zlev[0])+'m'
        else:
            data.str_dep=', dep: '+str(data.zlev[0])+'m-'+str(data.zlev[-1])+'m'
    #___________________________________________________________________________
    return data
    
    
#___LOAD FESOM2.0 DATA OVER TIME IN INDEX BOX___________________________________
#
#_______________________________________________________________________________
def fesom_load_data_overtime(mesh,data,do_output=True):
    #___________________________________________________________________________
    # number of years to average 
    nyi = data.year[1]-data.year[0]+1
    
    # array with years
    ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    # number of month to average
    nmi = np.array(data.month).size
    # array with month
    ami = np.array(data.month)
    
    ##___________________________________________________________________________
    ## 3d lsmask
    #lsmask3d = np.zeros((mesh.n2dn,mesh.nlev-1))
    #for ii in range(0,mesh.n2dn):
        #lsmask3d[ii,0:mesh.nodes_2d_iz[ii]-1]=1
    
    if do_output==True:
        print('     -----+-----------------------------------+------------')
        print('     Year |               MON                 |')
        print('     -----+-----------------------------------+------------')
    print('     --> '+data.path)
    print('     --> '+data.var)
    #____START YEAR LOOP________________________________________________________
    aux_datavar = data.var
    count_ti=0
    for yi in range(0, nyi):
        if do_output==True: print('     {:4.0f} |'.format(ayi[yi]),end=''),
        
        tstart = time.time();
        if data.var.find('norm')!=-1 or data.var=='ptemp':
            aux_datavar,aux_datavar2 = 'u','v'
            if data.var =='ptemp':
                aux_datavar,aux_datavar2 = 'temp','salt'
            #_______________________________________________________________________
            # build filename where data are located
            fname  =   data.path+'/'\
                        + aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            fname2 =   data.path+'/'\
                        + aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            #_______________________________________________________________________
            # open netcdf file --> read NETCDF4
            ncid  = Dataset(fname, 'r') 
            ncid2 = Dataset(fname2, 'r')
            
        elif data.var.find('tuv')!=-1 or data.var.find('suv')!=-1:
            if   data.var.find('tuv')!=-1 : aux_datavar,aux_datavar2,aux_datavar3 = 'u','v','temp'
            elif data.var.find('suv')!=-1 : aux_datavar,aux_datavar2,aux_datavar3 = 'u','v','salt'
            #_______________________________________________________________________
            # build filename where data are located
            fname  =   data.path+'/'\
                        + aux_datavar + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            fname2 =   data.path+'/'\
                        + aux_datavar2 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            fname3 =   data.path+'/'\
                        + aux_datavar3 + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'        
            #_______________________________________________________________________
            # open netcdf file --> read NETCDF4
            ncid  = Dataset(fname, 'r') 
            ncid2 = Dataset(fname2, 'r') 
            ncid3 = Dataset(fname3, 'r') 
        else:
            #_______________________________________________________________________
            # build filename where data are located
            fname =   data.path+'/'\
                    + data.var + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            
            #_______________________________________________________________________
            # open netcdf file --> read NETCDF4
            ncid = Dataset(fname, 'r') 
        
        # read dimension of variable 
        ncdims=ncid.variables[aux_datavar].shape
        
        # number of dimension is 2d or 3d
        dim_num = len(ncdims)
        
        # check time dimension of loaded data aka number of records
        nrec  = ncdims[0]
        nsmple= ncdims[1]
        nlev  = 0
        if dim_num == 3: nlev=ncdims[2]
        
        
        #_______________________________________________________________________
        # at first loaded year allocate matrices to calc means
        if yi==0:
            nti = np.min([nmi,nrec])*nyi
            
            #___________________________________________________________________
            # 3d lsmask
            if dim_num == 3:
                lsmask3d = np.zeros((mesh.n2dn,nlev))
                if mesh.nlev==nlev:
                    for ii in range(0,mesh.n2dn):
                        lsmask3d[ii,0:mesh.nodes_2d_iz[ii]]=1
                else:    
                    for ii in range(0,mesh.n2dn):
                        lsmask3d[ii,0:mesh.nodes_2d_iz[ii]-1]=1
            
            #___________________________________________________________________
            # if data are 3d and no loadable depth layers are given load all fesom 
            # layers
            if dim_num==3 and (len(data.depth)==0 or np.all(data.depth==mesh.zlev)==True):
                # allocate array for mean box index data
                if data.which_obj=='box':
                    data.value = []
                    for ii in range(0,len(data.box_define)):
                        data.value.append([])
                        data.value[ii] = np.zeros((nti,ncdims[2]),dtype='float32')
                else:
                    data.value = np.zeros((nti,ncdims[2]),dtype='float32')
            else:
                # allocate array for mean box index  data
                if data.which_obj=='box':
                    data.value = []
                    for ii in range(0,len(data.box_define)):
                        data.value.append([])
                        data.value[ii] = np.zeros((nti,),dtype='float32')
                else:
                    data.value = np.zeros((nti,),dtype='float32')
                
            data.time  = np.zeros((nti,),dtype='float32')
            
            #___________________________________________________________________
            # setup variable name, description and unit
            ncvar = ncid.variables[aux_datavar]
            #data.sname=ncvar.name
            data.sname=data.var
            for attrname in ncvar.ncattrs():
                if attrname=='description':
                    data.lname=ncvar.getncattr(attrname)
                elif attrname=='units':
                    data.unit='['+ncvar.getncattr(attrname)+']'
            if data.var =='ptemp': data.lname='potential temperature'
            if data.var.find('tuv')!=-1: data.lname='temperature advection'
            if data.var.find('suv')!=-1: data.lname='salt advection'
            
        #_______________________________________________________________________
        ncval = np.copy(ncid.variables[aux_datavar][:])
        if data.var.find('norm')!=-1:
            ncval2 = np.copy(ncid2.variables[aux_datavar2][:])
            
        if data.var.find('tuv')!=-1 or data.var.find('tuv')!=-1:
            ncval3 = np.copy(ncid3.variables[aux_datavar][:])
            ncval  = ncval*ncval3
            ncval2 = ncval2*ncval3
            del ncval3
            ncid3.close()
            
        if data.var.find('norm')!=-1:
            ncval = np.sqrt(ncval[:]**2+ncval2[:]**2)
            del ncval2
            ncid2.close()
            
        elif data.var=='ptemp':    
            import numpy.matlib
            import seawater as sw
            ncval2  = np.copy(ncid2.variables[aux_datavar2][:])
            depth3d = np.zeros(ncval.shape)
            press3d = np.zeros(ncval.shape)
            for nreci in range(0,nrec):
                depth3d[nreci,:,:] = np.matlib.repmat((mesh.zlev[0:-1]+(mesh.zlev[1::]-mesh.zlev[0:-1])/2),nsmple,1)
                press3d[nreci,:,:] = np.matlib.repmat(mesh.nodes_2d_yg[0:mesh.n2dn],ncval.shape[2],1).transpose()
            press3d = sw.pres(depth3d,press3d)
            ncval   = sw.ptmp(ncval2,ncval,press3d)
        
        #_______________________________________________________________________
        # put in nan 3d landsea mask
        ncval[:,lsmask3d==0]=np.nan
        
        #_______________________________________________________________________
        # check if a single time slice record is selcted
        # no ?  --> than use data.month to select time subset
        if len(data.record)==0:
            #___________________________________________________________________
            if nrec==1:
                if do_output==True and yi==0 : print('annual |',end=''),
            elif nrec==12:
                if do_output==True and yi==0 : print('monthly|',end=''),
            #elif nrec==73:
                #print('5daily |'),
            elif nrec==365:
                if do_output==True and yi==0 : print('daily  |',end=''),
            else:    
                print(' --> error: temporal length of data unclear or not supported!!!|'),
                break
            
            #___________________________________________________________________
            # select data like they are stored in the file 
            if         data.which_mean=='None':
                if data.which_obj=='box':
                    for ii in range(0,len(data.box_define)):
                        ncval2 = np.nanmean(ncval[:,np.where(data.box_idx[ii][0:mesh.n2dn])[0],:],axis=1)
                        data.value[ii][count_ti:count_ti+nrec,:] = ncval2
                data.time[count_ti:count_ti+nrec   ] = ayi[yi]+np.arange(0.0,nrec,1.0)/nrec
                count_ti = count_ti+nrec
                
            #___________________________________________________________________
            # select data by data.month calculate monthly mean from
            elif     data.which_mean=='monthly' and nrec!=1:
                if data.which_obj=='box':
                    for ii in range(0,len(data.box_define)):
                        #_______________________________________________________________
                        # calc mean over points in box
                        ncval2 = np.nanmean(ncval[:,np.where(data.box_idx[ii][0:mesh.n2dn])[0],:],axis=1)
                        #_______________________________________________________________
                        # select from monthly data
                        if   nrec==12 : 
                            # select month subset
                            ncval2 = fesom_time_depth_mean(mesh,ncval2,np.array(data.month)-1,[],do_tmean=False) 
                        # select from daily data
                        elif nrec==365 : 
                            # find out which day belongs to selected subset
                            idx_sel = sel_timesubset_daily(data)
                            ncval2=fesom_time_depth_mean(mesh,ncval2,idx_sel==True,[],do_tmean=False)
                        #_______________________________________________________________
                        data.value[ii][count_ti:count_ti+ncval.shape[0],:] = ncval2[:,:]
                        
                data.time[count_ti:count_ti+ncval.shape[0]   ] = ayi[yi]+np.arange(0.0,ncval.shape[0],1.0)/ncval.shape[0]
                count_ti = count_ti+ncval.shape[0]
                
            #___________________________________________________________________
            # select data by data.month calculate seasonal mean from
            elif     data.which_mean=='seasonal' and nrec!=1:
                if data.which_obj=='box':
                    for ii in range(0,len(data.box_define)):
                        #_______________________________________________________________
                        # calc mean over points in box
                        ncval2 = np.nanmean(ncval[:,np.where(data.box_idx[ii][0:mesh.n2dn])[0],:],axis=1)
                        #_______________________________________________________________
                        # select from monthly data
                        if nrec==12 : 
                            # select month subset
                            ncval2 = fesom_time_depth_mean(mesh,ncval2,np.array(data.month)-1,[],do_tmean=True)
                        #_______________________________________________________________
                        # select from daily data
                        elif nrec==365 : 
                            # find out which day belongs to selected subset
                            idx_sel = sel_timesubset_daily(data)
                            ncval2=fesom_time_depth_mean(mesh,ncval2,idx_sel==True,[],do_tmean=True)
                            
                        #_______________________________________________________________
                        data.value[ii][count_ti:count_ti+1,:] = ncval2[:,:]
                data.time[count_ti:count_ti+1   ] = ayi[yi]
                count_ti = count_ti+1
            
            #___________________________________________________________________
            # calculate annual mean from daily, monthly or annual data
            elif     data.which_mean=='annual' or nrec==1:
                if data.which_obj=='box':
                    for ii in range(0,len(data.box_define)):
                        #_______________________________________________________________
                        # calc mean over points in box
                        ncval2 = np.nanmean(ncval[:,np.where(data.box_idx[ii][0:mesh.n2dn])[0],:],axis=1)
                        #_______________________________________________________________
                        data.value[ii][count_ti:count_ti+1,:] = ncval2.mean(axis=0)
                data.time[count_ti:count_ti+1   ] = ayi[yi]
                count_ti = count_ti+1
        
        #_______________________________________________________________________
        # close netcdf file
        ncid.close()
        
        
        #_______________________________________________________________________
        tend = time.time();        
        if do_output==True: print(' --> t={:2.2f}s'.format(tend-tstart))
    #___________________________________________________________________________
    # cut of time_data
    if data.which_obj=='box':
        for ii in range(0,len(data.box_define)):
            data.value[ii] = data.value[ii][0:count_ti,:]
            
    data.time  = data.time[0:count_ti]
    if len(data.zlev)==0:
        if nlev==mesh.nlev:
            data.zlev = mesh.zlev
        elif nlev==mesh.nlev-1:    
            data.zlev = mesh.zmid
    #if len(data.depth)==0: data.depth = mesh.zlev
    
    #___________________________________________________________________________
    # set up descriptiv variable names
    if len(data.record)==0:
        mon_list   = np.array(['J','F','M','A','M','J','J','A','S','O','N','D'])
        mon_list2  = np.array(['January','February','March','April','May','June','July','August','September','October','Novemver','December'])
        # year info 
        if nyi == 1:
            str_time_1 = 'y: '+str(data.year[0])
        else:    
            str_time_1 = 'y: '+str(data.year[0])+'-'+str(data.year[1])
            
        # month info 
        if nmi == 1 and nrec>=12:
            str_time_1 = str_time_1+' ['+mon_list2[np.array(data.month)-1][0]+']'
        elif nmi>1 and nmi<12 and nrec>=12:
            str_mon=''
            for ii in range(0,nmi):
                if ami[ii]==12 : str_mon=mon_list[ami[ii]-1]+str_mon
                else            : str_mon=str_mon+mon_list[ami[ii]-1]
            str_time_1 = str_time_1 + ' [' + str_mon + ']'
        data.str_time = str_time_1
    else: 
        data.str_time = 'rec_i: '+str(data.record[0])
    
    #___________________________________________________________________________
    return(data)
    
    
#___LOAD FESOM2.0 DATA AS TIME AVERAGED HORIZONTAL SLICE________________________
#
#_______________________________________________________________________________
def fesom_load_blowup(mesh,data,do_output=True,which_file=['blowup','oce']):
    
    #___________________________________________________________________________
    # number of years to average 
    nyi = data.year[1]-data.year[0]+1
    
    # array with years
    ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    # number of month to average
    nmi = np.array(data.month).size
    # array with month
    ami = np.array(data.month)
    
    #____START YEAR LOOP________________________________________________________
    fname = data.path+'/'+'fesom.'+str(ayi[0])+'.'+which_file[1]+'.'+which_file[0]+'.nc'
    print('     --> '+fname)
    print('     --> '+data.var)
    
    #_______________________________________________________________________
    # open netcdf file --> read NETCDF4
    ncid = Dataset(fname, 'r') 
    
    
    # read dimension of variable 
    ncdims=ncid.variables[data.var].shape
    #print(ncdims)
    
    # number of dimension is 2d or 3d
    dim_num = len(ncdims)
    
    # check time dimension of loaded data aka number of records
    nrec  = ncdims[0]
    nsmple= ncdims[1]
    
    #_______________________________________________________________________
    # allocate array for mean data
    if dim_num==3 and  (len(data.depth)==0 or np.all(data.depth==mesh.zlev)==True):
        # allocate array for mean data
        mean_data = np.zeros((nsmple,ncdims[2]),dtype='float32')
    else:
        # allocate array for mean data
        mean_data = np.zeros((nsmple,),dtype='float32')
    
    #___________________________________________________________________
    # setup variable name, description and unit
    ncvar = ncid.variables[data.var]
    #data.sname=ncvar.name
    data.sname=data.var
    for attrname in ncvar.ncattrs():
        if attrname=='description':
            data.lname=ncvar.getncattr(attrname)
        elif attrname=='units':
            data.unit='['+ncvar.getncattr(attrname)+']'
    if data.var =='ptemp': data.lname='potential temperature'
    if data.var.find('tuv')!=-1: data.lname='temperature advection'
    if data.var.find('suv')!=-1: data.lname='salt advection'
    
    #_______________________________________________________________________
    ncval = np.copy(ncid.variables[data.var][:])
    ncval = fesom_time_depth_mean(mesh,ncval,data.record[0]-1,data.depth,do_tmean=False)
    
    #_______________________________________________________________________
    # close netcdf file
    ncid.close()
    
    #_______________________________________________________________________
    # now average over yearly subset
    mean_data = mean_data + ncval/nyi
    
    #___________________________________________________________________________
    # in case of ice data set no ice to nan
    if (data.var=='a_ice'):
        mean_data=mean_data*100.0
    #if (data.var=='a_ice' or data.var=='m_ice'):
        #mean_data[mean_data==0.0]=np.nan
    
    #___________________________________________________________________________
    # augment periodic boundary of data.value
    if nsmple==mesh.n2dn:
        # augment node data
        data.value = np.concatenate((mean_data,mean_data[mesh.pbndn_2d_i]))
        del mean_data
    else:
        # augment elem data
        data.value = np.concatenate((mean_data,mean_data[mesh.pbndtri_2d_i]))
        del mean_data
        # delete trinangle that are seen as abnormal 
        if len(mesh.abnormtri_2d_i)!=0:
            data.value = np.delete(data.value,mesh.abnormtri_2d_i)
    
    #___________________________________________________________________________
    return data
    
    
#___LOAD FESOM2.0 DATA AS 3D FOR VERY BIG MESHES LIKE FRON OR STORM_____________
# the time slice in each file are loaded and averaged serialy
#
#_______________________________________________________________________________
def fesom_load_data3d_4bm(mesh,data,do_output=True):
    
    #___________________________________________________________________________
    # number of years to average 
    nyi = data.year[1]-data.year[0]+1
    
    # array with years
    ayi = np.linspace(data.year[0],data.year[1],nyi,dtype='uint16')
    # number of month to average
    nmi = np.array(data.month).size
    # array with month
    ami = np.array(data.month)
    
    if do_output==True:
        print('')
        print('     -----+-----------------------------------+------------')
        print('     Year |               MON                 |')
        print('     -----+-----------------------------------+------------')
        print('     --> '+data.path)
        print('     --> '+data.var)
    #____START YEAR LOOP________________________________________________________
    for yi in range(0, nyi):
        if do_output==True: print('     {:4.0f} |'.format(ayi[yi]),end=''),
        tstart = time.time();
        
        #_______________________________________________________________________
        # build filename where data are located
        fname =   data.path+'/'\
               + data.var + '.' + data.runid + '.' + str(ayi[yi]) + '.nc'
            
        #_______________________________________________________________________
        # open netcdf file --> read NETCDF4
        ncid = Dataset(fname, 'r')
        
        # read dimension of variable 
        ncdims=ncid.variables[data.var].shape
        
        # number of dimension is 2d or 3d
        dim_num = len(ncdims)
        
        # check time dimension of loaded data aka number of records
        nrec  = ncdims[0]
        nsmple= ncdims[1]
        
        #_______________________________________________________________________
        # at first loaded year allocate matrices to calc means
        if yi==0:
            
            # if data are 3d and no loadable depth layers are given load all fesom 
            # layers
            if dim_num==3 and  (len(data.depth)==0 or np.all(data.depth==mesh.zlev)==True):
                # allocate array for mean data
                data.value = np.zeros((nsmple,ncdims[2]),dtype='float32')
            else:
                # allocate array for mean data
                data.value = np.zeros((nsmple,),dtype='float32')
            #___________________________________________________________________
            # setup variable name, description and unit
            data.sname=data.var
            for attrname in ncid.variables[data.var].ncattrs():
                if attrname=='description':
                    data.lname=ncid.variables[data.var].getncattr(attrname)
                elif attrname=='units':
                    data.unit='['+ncid.variables[data.var].getncattr(attrname)+']'
            if data.var =='ptemp': data.lname='potential temperature'
            if data.var.find('tuv')!=-1: data.lname='temperature advection'
            if data.var.find('suv')!=-1: data.lname='salt advection'
            
        #_______________________________________________________________________
        if nrec>1:
            for ii in range(nrec):
                if do_output==True: print('{:02d}|'.format(ii+1),end=''),
                data.value = data.value + ncid.variables[data.var][ii,:]/nyi/nrec
        else:
            data.value = data.value + ncid.variables[data.var][:].mean(axis=0)/nyi
        #_______________________________________________________________________
        # close netcdf file
        ncid.close()
        
        #_______________________________________________________________________
        tend = time.time();        
        if do_output==True: print(' --> t={:2.2f}s'.format(tend-tstart))
    #___________________________________________________________________________
    # in case of ice data set no ice to nan
    if (data.var=='a_ice'):
        mean_data=mean_data*100.0
    
    #___________________________________________________________________________
    # augment periodic boundary of data.value
    if nsmple==mesh.n2dn:
        # augment node data
        data.value = np.concatenate((data.value,data.value[mesh.pbndn_2d_i]))
    else:
        # augment elem data
        data.value = np.concatenate((data.value,data.value[mesh.pbndtri_2d_i]))
        # delete trinangle that are seen as abnormal 
        if len(mesh.abnormtri_2d_i)!=0:
            data.value = np.delete(data.value,mesh.abnormtri_2d_i)
    
    #___________________________________________________________________________
    if len(data.depth)==0:
        data.depth = mesh.zlev
    
    #___________________________________________________________________________
    # set up descriptiv variable names
    if len(data.record)==0:
        mon_list   = np.array(['J','F','M','A','M','J','J','A','S','O','N','D'])
        mon_list2  = np.array(['January','February','March','April','May','June','July','August','September','October','Novemver','December'])
        # year info 
        if nyi == 1:
            str_time_1 = 'y: '+str(data.year[0])
        else:    
            str_time_1 = 'y: '+str(data.year[0])+'-'+str(data.year[1])
            
        # month info 
        if nmi == 1 and nrec>=12:
            str_time_1 = str_time_1+' ['+mon_list2[np.array(data.month)-1][0]+']'
        elif nmi>1 and nmi<12 and nrec>=12:
            str_mon=''
            for ii in range(0,nmi):
                if ami[ii]==12 : str_mon=mon_list[ami[ii]-1]+str_mon
                else            : str_mon=str_mon+mon_list[ami[ii]-1]
            str_time_1 = str_time_1 + ' [' + str_mon + ']'
        data.str_time = str_time_1
    else: 
        data.str_time = 'rec_i: '+str(data.record[0])
    
    if dim_num==3 and len(data.depth)!=0:
        if len(data.depth)<=1:
            data.str_dep=', dep: '+str(data.depth[0])+'m'
        else:
            data.str_dep=', dep: '+str(data.depth[0])+'m-'+str(data.depth[-1])+'m'
    #___________________________________________________________________________
    #return data
    
    
    
#___DO VERTICAL INTERPOLATE AVERAGE OVER CERTAIN LAYERS_________________________
#
#_______________________________________________________________________________
def fesom_vinterp(data_in,mesh,levels):
     
    #___________________________________________________________________________
    # do vertical linear interpolation of certain layer
    dims = data_in.shape
    if   len(dims)==2: nsmple, nlev     = dims[0], dims[1]
    elif len(dims)==3: nt, nsmple, nlev = dims[0], dims[1], dims[2]
    
    #___________________________________________________________________________
    #data_zlev=mesh.zlev[0:-1]+(mesh.zlev[1::]-mesh.zlev[0:-1])/2
    if nlev==mesh.nlev:
        # data are on levels
        data_zlev = mesh.zlev
    elif nlev==mesh.nlev-1:
        ##data are on mid depth levels
        data_zlev=mesh.zlev[0:-1]+(mesh.zlev[1::]-mesh.zlev[0:-1])/2
    else:
        print(' --> error: wrong number of levels for vertical interpolation')
        
    #___________________________________________________________________________
    if   len(dims)==2 : data_out, aux_div = np.zeros((nsmple,),dtype='float32'), np.zeros((nsmple,),dtype='float32')
    elif len(dims)==3 : data_out, aux_div = np.zeros((nti,nsmple),dtype='float32'), np.zeros((nti,nsmple,),dtype='float32')
    
    #___________________________________________________________________________
    for di in range(0,len(levels)):
        #_______________________________________________________________________
        # find upper and lower layer indices
        #idx_dwn = np.array(np.where( data.depth[di]<=abs(mesh.zlev))).squeeze()
        #idx_dwn = np.array(np.where( levels[di]<=abs(mesh.zlev))).squeeze()
        idx_dwn = np.array(np.where( levels[di]<=abs(data_zlev))).squeeze()
        idx_dwn = idx_dwn[0]
        idx_up  = idx_dwn-1
        if idx_up<0: idx_up=0
        
        #_______________________________________________________________________
        # linear vertical interpolant
        #deltaz   = abs(mesh.zlev[idx_dwn])-abs(mesh.zlev[idx_up])
        #deltaz_i = abs(mesh.zlev[idx_dwn])-levels[di]
        deltaz   = abs(data_zlev[idx_dwn])-abs(data_zlev[idx_up])
        deltaz_i = abs(data_zlev[idx_dwn])-levels[di]
        
        #_______________________________________________________________________
        # data_in is defined on nodes temp, salt, ssh ...
        if nsmple==mesh.n2dn:
            if deltaz_i==0 or idx_dwn==idx_up:
                if   len(dims)==2:    
                    auxval = data_in[:,idx_dwn]
                    aux_div[mesh.nodes_2d_iz>=idx_dwn]=aux_div[mesh.nodes_2d_iz>=idx_dwn]+1.0
                else : 
                    auxval = data_in[:,:,idx_dwn]
                    aux_div[:,mesh.nodes_2d_iz>=idx_dwn]=aux_div[:,mesh.nodes_2d_iz>=idx_dwn]+1.0
                
            else:
                if   len(dims)==2:
                    auxval = data_in[:,idx_dwn]-(data_in[:,idx_dwn]-data_in[:,idx_up])*deltaz_i/deltaz
                    auxval[mesh.nodes_2d_iz<idx_dwn]=0.0
                    aux_div[mesh.nodes_2d_iz>=idx_dwn]=aux_div[mesh.nodes_2d_iz>=idx_dwn]+1.0
                else:
                    auxval = data_in[:,:,idx_dwn]-(data_in[:,:,idx_dwn]-data_in[:,:,idx_up])*deltaz_i/deltaz
                    auxval[:,mesh.nodes_2d_iz<idx_dwn]=0.0
                    aux_div[:,mesh.nodes_2d_iz>=idx_dwn]=aux_div[:,mesh.nodes_2d_iz>=idx_dwn]+1.0
                    
            data_out = data_out + auxval
        #_______________________________________________________________________
        #data_in is defined on elements u,v....
        elif nsmple==mesh.n2de:
            if deltaz_i==0 or idx_dwn==idx_up:
                if   len(dims)==2:
                    auxval = data_in[:,idx_dwn]
                    aux_div[mesh.elem0_2d_iz>=idx_dwn]=aux_div[mesh.elem0_2d_iz>=idx_dwn]+1.0
                else:
                    auxval = data_in[:,:,idx_dwn]
                    aux_div[:,mesh.elem0_2d_iz>=idx_dwn]=aux_div[:,mesh.elem0_2d_iz>=idx_dwn]+1.0
                    
            else:
                if   len(dims)==2:
                    auxval     = data_in[:,idx_dwn]-(data_in[:,idx_dwn]-data_in[:,idx_up])*deltaz_i/deltaz
                    auxval[mesh.elem0_2d_iz<idx_dwn]=0
                    aux_div[mesh.elem0_2d_iz>=idx_dwn]=aux_div[mesh.elem0_2d_iz>=idx_dwn]+1.0
                else:
                    auxval     = data_in[:,:,idx_dwn]-(data_in[:,:,idx_dwn]-data_in[:,:,idx_up])*deltaz_i/deltaz
                    auxval[:,mesh.elem0_2d_iz<idx_dwn]=0
                    aux_div[:,mesh.elem0_2d_iz>=idx_dwn]=aux_div[:,mesh.elem0_2d_iz>=idx_dwn]+1.0
            data_out = data_out + auxval
            
    #___________________________________________________________________________
    # do mean over averaged layers
    data_out[aux_div!=0.0] = data_out[aux_div!=0.0]/aux_div[aux_div!=0.0]
    data_out[aux_div==0.0      ]=np.nan
    data_out[np.isinf(data_out)]=np.nan
    
    #fig=plt.figure()
    #tri = Triangulation(mesh.nodes_2d_xg,mesh.nodes_2d_yg,mesh.elem_2d_i)
    #data_plot = np.concatenate((data_out,data_out[mesh.pbndn_2d_i]))
    #plt.tripcolor(tri,data_plot,edgecolors='None')
    #plt.colorbar()
    #plt.show
    #fig.canvas.draw()
    #STOP
    #___________________________________________________________________________
    return(data_out)
    
    
#___SELECT TIME SUBSET OF DAILY DATA TO CALCULATE SEASONAL MEAN_________________
#
#_______________________________________________________________________________
def sel_timesubset_daily(data):
    
    daypermonth=np.array([0,31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    daypermonth=daypermonth.cumsum()
    idx_day2mon = np.zeros((365,))
    idx_day     = np.arange(1,365+1,1)
    for ii in range(0,12): 
        idx_day2mon[(idx_day>daypermonth[ii]) & (idx_day<=daypermonth[ii+1])]=ii+1
    del idx_day
    del daypermonth
                
    # find out which day belongs to selected subset
    idx_sel = np.in1d(idx_day2mon,data.month)
    del idx_day2mon    
    
    return(idx_sel)
    
    
#___SELECT SUBSET AND DO TEMPORAL AND VERTICAL AVERAGE__________________________
#
#_______________________________________________________________________________
def fesom_time_depth_mean(mesh,ncval,sel_time,levels=[],do_tmean=True,do_zmean=True):
    
    dim_num = len(ncval.shape)
    if np.isscalar(sel_time)==True:
        nselt = 0
    else:
        nselt = len(sel_time)
    #___________________________________________________________________________
    # do no time selection
    #if len(sel_time)==0:
    if nselt==0:
        if dim_num==2:
            # case 2d variable
            #ncval= ncval[:,:]
            if np.isscalar(sel_time)==True:
                ncval= ncval[sel_time,:]
            elif do_tmean==True and np.isscalar(sel_time)==False: 
                ncval=ncval.mean(axis=0)
            
        elif dim_num==3:
            # case 3d variable
            #ncval= ncval[:,:,:]
            if np.isscalar(sel_time)==True:
                ncval= ncval[sel_time,:,:]
            elif do_tmean==True and np.isscalar(sel_time)==False: 
                ncval=ncval.mean(axis=0)
            # do vertical interpolation of certain layer if data.depth is not empty
            #if len(data.depth)!=0 : ncval=fesom_vinterp(ncval,mesh)
            if do_zmean==True:
                #if len(levels)!=0 : 
                if len(levels)!=0 and np.all(levels==mesh.zlev)==False:
                    ncval=fesom_vinterp(ncval,mesh,levels)
        else:
            print(' --> error: more than 3 variable dimensions are not supported' )
    
    #___________________________________________________________________________
    # select certain time slice and average
    else:
        if dim_num==2:
            # case 3d variable
            ncval= ncval[sel_time,:]
            if do_tmean==True: ncval=ncval.mean(axis=0)
        elif dim_num==3:
            # case 3d variable
            ncval= ncval[sel_time,:,:]
            if do_tmean==True: ncval=ncval.mean(axis=0)
            # do vertical interpolation of certain layer if data.depth is not empty
            #if len(data.depth)!=0 : ncval=fesom_vinterp(ncval,mesh)
            if do_zmean==True:
                #if len(levels)!=0 : 
                if len(levels)!=0 and np.all(levels==mesh.zlev)==False:
                    ncval=fesom_vinterp(ncval,mesh,levels)
        else:
            print(' --> error: more than 3 variable dimensions are not supported' )
    return(ncval)
    
    
#___CALCULATE ANOMALY BETWEEN TWO DATA OBJECTS__________________________________
#
#_______________________________________________________________________________
def fesom_data_anom(data,data2):
    anom = fesom_data([])
    
    #____data run variables______________________
    if data.var!=data2.var:
        anom.var                            = data2.var+'-'+data.var
    else:
        anom.var                            = data.var
    anom.runid, anom.path               = data.runid, data.path
    anom.descript                       = data2.descript+'-'+data.descript
    
    #____data time variables_____________________
    anom.year,  anom.month, anom.record = data.year, data.month, data.record
    anom.depth                          = data.depth
    if data.str_time!=data2.str_time:
        anom.str_time = data2.str_time+'-'+data.str_time
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
        anom.sname = data2.sname+'-'+data.sname
    else:
        anom.sname = data.sname
    if data.lname!=data2.lname:
        anom.lname = data2.lname+'-'+data.lname
    else:
        anom.lname = data.lname
    anom.unit   = data.unit
    #____data plot varaibles_____________________
    #anom.cmap, anom.crange, anom.cnumb  = data.cmap, data.crange, data.cnumb
    anom.which_plot                     = data.which_plot 
    #if len(data.crange)!=0 :anom.crange=[]
    #____data varaible___________________________
    anom.value                          = data2.value-data.value
    if len(data.value2)!=0:
        anom.value2                          = data2.value2-data.value2
    if len(data.value3)!=0:
        anom.value3                          = data2.value3-data.value3
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
