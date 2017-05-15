# This file is part of pyfesom
#
################################################################################
#
# Original code by Nikolay Koldunov, 2016
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

def create_indexes_and_distances(mesh, lons, lats, k=1, n_jobs=2, ):
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
    xs, ys, zs = lon_lat_to_cartesian(mesh.x2, mesh.y2)
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

    '''
    #print distances
    if (distances is None) or (inds is None):
        
        if how=='nn':
            distances, inds = create_indexes_and_distances(mesh, lons, lats,\
                                                           k=1, n_jobs=n_jobs)
        elif how=='idist':
            distances, inds = create_indexes_and_distances(mesh, lons, lats,\
                                                           k=k, n_jobs=n_jobs)

    if distances.ndim == 1:
        #distances_ma = np.ma.masked_greater(distances, radius_of_influence)
        data_interpolated = data[inds]

        data_interpolated[distances>=radius_of_influence] = np.nan
        
        data_interpolated = data_interpolated.reshape(lons.shape)
        data_interpolated = np.ma.masked_invalid(data_interpolated)
    else:
        distances_ma = np.ma.masked_greater(distances, radius_of_influence)
        
        w = 1.0 / distances_ma**2
        data_interpolated = np.ma.sum(w * data[inds], axis=1) / np.ma.sum(w, axis=1)
        data_interpolated.shape = lons.shape
        data_interpolated = np.ma.masked_invalid(data_interpolated)
    
    return data_interpolated

