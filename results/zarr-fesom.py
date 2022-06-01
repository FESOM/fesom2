#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dask.array as da
from zarr.storage import DirectoryStore
import zarr 
import numpy as np
import xarray as xr


# In[2]:


def dataset_from_fesom_zarr(dataset_path):
    store = DirectoryStore(dataset_path)
    root = zarr.open_group(store)
    variable_dict = {} 
    for variable in root.group_keys():
        variable_group = root[variable]
        # below variable_group[k] corresponsds to group of arrays output by a pe.
        # find array dims from first pe's arrays ! TODO: later put it in variable group attrs.
        dims_list=variable_group[0].attrs['_ARRAY_DIMENSIONS']
        # find index of nod2d it has to be on left end to be able to concattinated using dask
        nod2d_index = dims_list.index('nod2d')
        nod2d_on_left = True if nod2d_index==0 else False
        darr_conc_list = [da.from_zarr(variable_group[k]).swapaxes(0,1) if not nod2d_on_left else
                          da.from_zarr(variable_group[k]) for k in variable_group.array_keys()] 
            
        darr_conc = da.concatenate(darr_conc_list)
        # reswap dims after concatination if nod2d was not on left
        darr_conc = darr_conc.swapaxes(0,1) if not nod2d_on_left else darr_conc
        variable_dict.update({variable:(dims_list,darr_conc)})
    return xr.Dataset(variable_dict,coords={})

# above function contains dask arrays, no data is completely loaded yet, 
# below check_bounds 1. indirectly checks loading all the data 2.checks bounds 
def check_bounds(xr_dataset):
    from pprint import pprint
    var_bounds = {}
    for variable in xr_dataset.data_vars.keys():
        var_bounds.update({variable:{'min': xr_dataset[variable].min().values,
                                     'max': xr_dataset[variable].max().values}})
    pprint(var_bounds)


# In[3]:


fesom_da=dataset_from_fesom_zarr('test2.zarr')
fesom_da


# In[4]:


# if prefered to have time at left end of dims of variables, e.g, sst(time,lev,nod2d)
fesom_da.transpose("time",...)


# In[5]:


# check loading of all data and bounds of variables


# In[6]:


check_bounds(fesom_da)

