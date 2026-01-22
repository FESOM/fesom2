# This file is part of pyfesom
#
################################################################################
#
# Original code by Nikolay Koldunov, 2016
# adapted by Patrick Scholz, 2017 for my plot stuff
#
################################################################################

from scipy.spatial import cKDTree
import numpy as np

def lon_lat_to_cartesian(lon, lat, R = 6371000):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R. Taken from http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z

def create_indexes_and_distances(mesh, lons, lats, k=1, n_jobs=2,which='node'):
    '''
    Creates KDTree object and query it for indexes of points in FESOM mesh that are close to the
    points of the target grid. Also return distances of the original points to target points.

 
    Parameters
    ----------
    mesh : fesom_mesh object
        pyfesom mesh representation
    lons/lats : array
        2d arrays with target grid values.
    k : int
        k-th nearest neighbors to return.
    n_jobs : int, optional
        Number of jobs to schedule for parallel processing. If -1 is given
        all processors are used. Default: 1.

    Returns
    -------
    distances : array of floats
        The distances to the nearest neighbors. 
    inds : ndarray of ints
        The locations of the neighbors in data.

    '''
    try:
        mesh.x2
    except AttributeError:
        if which=='node':
            xs, ys, zs = lon_lat_to_cartesian(mesh.nodes_2d_xg, mesh.nodes_2d_yg)
        elif which=='elem':
            xs, ys, zs = lon_lat_to_cartesian(mesh.nodes_2d_xg[mesh.elem_2d_i].sum(axis=1)/3, 
                                            mesh.nodes_2d_yg[mesh.elem_2d_i].sum(axis=1)/3)
    else:
        if which=='node':
            xs, ys, zs = lon_lat_to_cartesian(mesh.x2, mesh.y2)
        elif which=='elem':
            xs, ys, zs = lon_lat_to_cartesian(mesh.x2[mesh.elem].sum(axis=1)/3, 
                                              mesh.y2[mesh.elem].sum(axis=1)/3)
    
    xt, yt, zt = lon_lat_to_cartesian(lons.flatten(), lats.flatten())
    
    tree = cKDTree(list(zip(xs, ys, zs)))
    distances, inds = tree.query(list(zip(xt, yt, zt)), k = k, n_jobs=n_jobs)
    
    return distances, inds

def fesom2regular(data, mesh, lons, lats, distances=None, \
                  inds=None, how='nn', k=10, radius_of_influence=100000, n_jobs = 2 ):
    '''
    Interpolates data from FESOM mesh to target (usually regular) mesh.

    Parameters
    ----------
    data : array
        1d array that represents FESOM data at one 
    mesh : fesom_mesh object
        pyfesom mesh representation
    lons/lats : array
        2d arrays with target grid values.
    distances : array of floats, optional
        The distances to the nearest neighbors.
    inds : ndarray of ints, optional
        The locations of the neighbors in data.
    how : str
       Interpolation method. Options are 'nn' (nearest neighbor) and 'idist' (inverce distance)
    k : int
        k-th nearest neighbors to use. Only used when how==idist
    radius_of_influence : int
        Cut off distance in meters.
    n_jobs : int, optional
        Number of jobs to schedule for parallel processing. If -1 is given
        all processors are used. Default: 1.
    
    Returns
    -------
    data_interpolated : 2d array
        array with data interpolated to the target grid.
    '''
    which='node'
    try:
        mesh.x2
    except AttributeError:
        if data.shape[0]==mesh.n2dea: which='elem'
    else:
        if data.shape[0]==mesh.e2d: which='elem'
    
    #print distances
    if (distances is None) or (inds is None):
        
        if how=='nn':
            distances, inds = create_indexes_and_distances(mesh, lons, lats,\
                                                           k=1, n_jobs=n_jobs,which=which)
        elif how=='idist':
            distances, inds = create_indexes_and_distances(mesh, lons, lats,\
                                                           k=k, n_jobs=n_jobs,which=which)

    if distances.ndim == 1:
        #distances_ma = np.ma.masked_greater(distances, radius_of_influence)
        data_interpolated = data[inds]

        data_interpolated[distances>=radius_of_influence] = np.nan
        
        data_interpolated = data_interpolated.reshape(lons.shape)
        data_interpolated = np.ma.masked_invalid(data_interpolated)
        data_interpolated = np.ma.filled(data_interpolated,np.nan)
    else:
        distances_ma = np.ma.masked_greater(distances, radius_of_influence)
        #distances_ma = np.ma.masked_greater(distances, np.nan)
        
        w = 1.0 / distances_ma**2
        data_interpolated = np.ma.sum(w * data[inds], axis=1) / np.ma.sum(w, axis=1)
        data_interpolated.shape = lons.shape
        data_interpolated = np.ma.masked_invalid(data_interpolated)
        data_interpolated = np.ma.filled(data_interpolated,np.nan)
        
        #data_interpolated[distances[:,1]>=radius_of_influence] = np.nan
        #data_interpolated[np.where(distances.min(axis=1)>=radius_of_influence)[0]] = np.nan
    
    return data_interpolated

def fesom2clim(data, depth, mesh, climatology, verbose=True, radius_of_influence=100000):
    '''
    Interpolation of fesom data to grid of the climatology.
    !NOTE! there is no interpolation to the climatology level, we interpolate
    the model level that is closest to level in climatology. 

    Parameters
    ----------
    data  : array
        1d array of FESOM 2d data slice
    depth : int 
        depth of the slice
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

def regular2regular(data, ilons, ilats, olons, olats, distances=None, \
                  inds=None, how='nn', k=10, radius_of_influence=100000, n_jobs = 2):
    '''
    Interpolates from regular to regular grid. 
    It's a wraper around `fesom2regular` that creates an object that 
    mimic fesom mesh class and contain only coordinates and flatten the data.

    Parameters
    ----------
    data : array
        1d or 2d array that represents gridded data at one level.
    ilons : arrea
    mesh : fesom_mesh object
        pyfesom mesh representation
    lons/lats : array
        2d arrays with target grid values.
    distances : array of floats, optional
        The distances to the nearest neighbors.
    inds : ndarray of ints, optional
        The locations of the neighbors in data.
    how : str
       Interpolation method. Options are 'nn' (nearest neighbor) and 'idist' (inverce distance)
    k : int
        k-th nearest neighbors to use. Only used when how==idist
    radius_of_influence : int
        Cut off distance in meters.
    n_jobs : int, optional
        Number of jobs to schedule for parallel processing. If -1 is given
        all processors are used. Default: 1.
    
    Returns
    -------
    data_interpolated : 2d array
        array with data interpolated to the target grid.

    '''
    fmesh = namedtuple('mesh', 'x2 y2')
    mesh = fmesh(x2=ilons.ravel(), y2=ilats.ravel())

    data = data.ravel()

    data_interpolated = fesom2regular(data, mesh, lons=olons, lats=olats,\
                                      distances=distances, \
                                      inds=inds, how=how,\
                                      k=k, radius_of_influence=radius_of_influence,\
                                      n_jobs = n_jobs)
    
    
    return data_interpolated

def regular2clim(data, ilons, ilats, izlevs, climatology, levels=None, verbose=True):
    '''
    Interpolation of data on the regular grid to climatology for the set of levels.

    Parameters
    ----------
    data : array
        3d array of data on regular grid for one time step
    ilons : array
        1d or 2d array of longitudes
    ilats : array
        1d or 2d array of latitudes
    izlevs : array
        depths of the source data
    climatology: climatology object
        FESOM climatology object
    levels : list like
        list of levels to interpolate. At present you can use only
        standard levels of WOA. If levels are not specified, all standard WOA
        levels will be used. 


    Returns
    -------
    xx : 2d array
        longitudes
    yy : 2d array
        latitudes
    out_data : 2d array
       array with data interpolated to climatology level

    '''
    if ((ilons.ndim == 1) & (ilats.ndim == 1)):
        ilons, ilats = np.meshgrid(ilons, ilats)
    elif (ilons.ndim == 2) & (ilats.ndim == 2):
        pass
    else:
        raise ValueError('Wrong dimentions for ilons and ilats')

    fmesh = namedtuple('mesh', 'x2 y2 zlevs')
    mesh = fmesh(x2=ilons.ravel(), y2=ilats.ravel(), zlevs=izlevs)

    #data = data.ravel()

    if levels is None:
        levels = climatology.z
    else:
        levels = np.array(levels)
        check = np.in1d(levels, climatology.z)
        if False in check:
            raise ValueError('Not all of the layers that you specify are WOA2005 layers')

    xx,yy = np.meshgrid(climatology.x, climatology.y)
    out_data=np.zeros((levels.shape[0],climatology.T.shape[1], climatology.T.shape[2]))
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)

    for dep_ind in range(len(levels)):
        wdep=levels[dep_ind]
        if verbose:
            print('interpolating to level: {}'.format(str(wdep)))
        dep_up=[z for z in abs(mesh.zlevs) if z<=wdep][-1]
        if verbose:
            print('Upper level: {}'.format(str(dep_up)))
        dep_lo=[z for z in abs(mesh.zlevs) if z>wdep][0]
        if verbose:
            print('Lower level: {}'.format(str(dep_lo)))
        i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
        i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
        dind_up=(abs(mesh.zlevs-dep_up)).argmin()
        dind_lo=(abs(mesh.zlevs-dep_lo)).argmin()
        data2=i_up*data[dind_up,:,:]
        data2=data2+i_lo*data[dind_lo,:,:]
        #zz[dep_ind,:,:] = pf.fesom2regular(data2, mesh, xx,yy)
        out_data[dep_ind,:,:] = regular2regular(data2, ilons, ilats,\
                                                xx, yy,\
                                                distances=distances, inds=inds)
    depth_indexes = [np.where(climatology.z==i)[0][0] for i in levels]
    out_data = np.ma.masked_where(climatology.T[depth_indexes,:,:].mask, out_data)
    return xx, yy, out_data

def clim2regular(climatology, param, olons, olats, \
                 levels=None, verbose=True, radius_of_influence=100000):
    '''
    Interpolation of data on the regular grid to climatology for the set of levels.

    Parameters
    ----------
    climatology: climatology object
        FESOM climatology object
    param : str
        name of the parameter to interpolate. Only 'T' for temperature and
        'S' for salinity are currently supported.
    olons : array
        1d or 2d array of longitudes to interpolate climatology on.
    olats : array
        1d or 2d array of longitudes to interpolate climatology on.
    levels : list like
        list of levels to interpolate to.
    

    Returns
    -------
    xx : 2d array
        longitudes
    yy : 2d array
        latitudes
    out_data : 2d array
       array with climatology data interpolated to desired levels
    '''
    clons,clats = np.meshgrid(climatology.x, climatology.y)

    fmesh = namedtuple('mesh', 'x2 y2 zlevs')
    mesh = fmesh(x2=clons.ravel(), y2=clats.ravel(), zlevs=climatology.z)

    #data = data.ravel()

    if levels is None:
        raise ValueError('provide levels for interpolation')
    else:
        levels = np.array(levels)   
        
    if olons.ndim == 1:
        xx,yy = np.meshgrid(olons, olats)
    elif olons.ndim == 2:
        xx = olons
        yy = olats
    else:
        raise ValueError('Wrong dimentions for olons')


    out_data=np.zeros((levels.shape[0],olons.shape[0], olons.shape[1]))
    distances, inds = create_indexes_and_distances(mesh, xx, yy,\
                                                k=10, n_jobs=2)

    if param == 'T':
        data = climatology.T[:]
    elif param == 'S':
        data = climatology.S[:]
    else:
        raise ValueError('Wrong variable {}'.format(param))


    for dep_ind in range(len(levels)):
        wdep=levels[dep_ind]
        if verbose:
            print('interpolating to level: {}'.format(str(wdep)))
        dep_up=[z for z in abs(mesh.zlev) if z<=wdep][-1]
        if verbose:
            print('Upper level: {}'.format(str(dep_up)))
        dep_lo=[z for z in abs(mesh.zlev) if z>wdep][0]
        if verbose:
            print('Lower level: {}'.format(str(dep_lo)))
        i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
        i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
        dind_up=(abs(mesh.zlevs-dep_up)).argmin()
        dind_lo=(abs(mesh.zlevs-dep_lo)).argmin()
        data2=i_up*data[dind_up,:,:]
        data2=data2+i_lo*data[dind_lo,:,:]
        #zz[dep_ind,:,:] = pf.fesom2regular(data2, mesh, xx,yy)
        out_data[dep_ind,:,:] = regular2regular(data2, clons, clats,\
                                                xx, yy,\
                                                distances=distances, inds=inds,\
                                                radius_of_influence=radius_of_influence)
    #depth_indexes = [np.where(climatology.z==i)[0][0] for i in levels]
    #out_data = np.ma.masked_where(climatology.T[depth_indexes,:,:].mask, out_data)
    return xx, yy, out_data
