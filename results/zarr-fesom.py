#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dask.array as da
from zarr.storage import DirectoryStore
import zarr 
import numpy as np
import xarray as xr


# In[12]:


def dataset_from_fesom_zarr(dataset_path):
    store = DirectoryStore(dataset_path)
    root = zarr.open_group(store)
    variables=list(root.group_keys())
    variable_group = root[variables[0]]
    darr_conc=da.concatenate([da.from_zarr(variable_group[k], 
                                           chunks=variable_group[k].chunks) for k in variable_group.array_keys()])
    return xr.Dataset({variables[0]:(('nod2d','time'),darr_conc)},coords={})


# In[13]:


fesom_da=dataset_from_fesom_zarr('test2.zarr')
fesom_da


# In[11]:


fesom_da.to_netcdf('test2.nc')


# In[9]:


fesom_da.sst.values


# In[10]:


fesom_da.sst.values.min(),fesom_da.sst.values.max()


# In[ ]:


# still need to add coordinate information

