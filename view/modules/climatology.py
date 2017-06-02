# This file is part of pyfesom
#
################################################################################
#
# Original code by Dmitry Sidorenko, 2013
#
# Modifications:
#   Nikolay Koldunov, 2016
#          - change to netCDF4 
#          - change scipy griddata interpolation to KDTree for speed
# TODO
# Add seasonal climatology
################################################################################
from regriding import fesom2regular, create_indexes_and_distances
import numpy as np
import scipy as sc
from numpy import nanmean
from netCDF4 import Dataset
import os


class climatology(object):
    '''
    Class that contains information from ocean 1 degree annual climatology.
    Presently options are WOA2005 and PHC3.0

    Parameters
    ----------
    path : str
        Path to the directory with climatology files
    climname : str
        Name of the climatology ('woa05' or 'phc')

    Returns
    -------
    object with climatology fields

    '''
    def __init__(self, path, climname='woa05'):
        if climname=='woa05':
            ncfile = Dataset(os.path.join(path, 'woa2005TS.nc'))
            self.T = np.copy(ncfile.variables['t00an1'][0,:,:,:])
            x=np.copy(ncfile.variables['lon'][:])
            x[x>180]=x[x>180]-360
            ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
            x=np.sort(x)
            self.x=x
            self.y=ncfile.variables['lat'][:]
            self.z=ncfile.variables['depth'][:]
            self.T[:,:,:]=self.T[:,:,ind]
            self.S=np.copy(ncfile.variables['s00an1'][0,:,:,:])
            self.S[:,:,:]=self.S[:,:,ind]
            ncfile.close()
            self.Tyz=nanmean(self.T, 2)
            self.Syz=nanmean(self.S, 2)
            self.T = np.ma.masked_greater(self.T,1000)
            self.S = np.ma.masked_greater(self.S,1000)

        if climname=='phc':
            ncfile = Dataset(os.path.join(path, 'phc3.0_annual.nc'))
            self.T = np.copy(ncfile.variables['temp'][:,:,:])
            x=np.copy(ncfile.variables['lon'][:])
            x[x>180]=x[x>180]-360
            ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
            x=np.sort(x)
            self.x=x
            self.y=ncfile.variables['lat'][:]
            self.z=ncfile.variables['depth'][:]
            self.T[:,:,:]=self.T[:,:,ind]
            self.S=np.copy(ncfile.variables['salt'][:,:,:])
            self.S[:,:,:]=self.S[:,:,ind]
            ncfile.close()
            self.Tyz=nanmean(self.T, 2)
            self.Syz=nanmean(self.S, 2)


def fesom_2_clim(data, depth, mesh, climatology, verbose=True, radius_of_influence=100000):
    '''
    Interpolation of fesom data to grid of the climatology.

    Parameters
    ----------
    data  : array
        1d array of FESOM 2d data slice
    depth : depth of the slice
    mesh  : mesh object
        FESOM mesh object        
    climatology: climatology object
        FESOM climatology object

    Returns
    -------
    iz : the index of the closest depth in climatology to the input depth
    xx : 2d array longitudes
    yy : 2d array latitudes
    out_data : 2d array
       array with data interpolated to climatology level

    '''
    xx,yy = np.meshgrid(climatology.x, climatology.y)
    out_data=np.copy(climatology.T)
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)
    
    #import pdb
    #pdb.set_trace()
    iz=abs(abs(climatology.z)-abs(depth)).argmin()
    print('the model depth is: ', depth, '; the closest depth in climatology is: ', climatology.z[iz])
    out_data=fesom2regular(data, mesh, xx, yy, distances=distances, inds=inds, radius_of_influence=radius_of_influence)
    out_data[np.isnan(climatology.T[iz,:,:])]=np.nan
    return iz, xx, yy, out_data