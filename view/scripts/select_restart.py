# -*- coding: utf-8 -*-

# This script will select one time step from the restart file 
# and put it in to a new restart file with only one time step.
# This is needed to save space when copy very large restarts.
# In this example the time step 6 (in python counting starts with zero)
# will be selected. 
# Note that there might be incompatability between different versions
# of netCDF libraries, so in case of problems make sure you compile 
# FESOM2 and xarray with the same version of netCDF
 
import xarray as xr

# input files
ocean = '/work/ollie/nkolduno/output_CORE_test1/fesom.1948.oce.restart.nc'
ice   = '/work/ollie/nkolduno/output_CORE_test1/fesom.1948.ice.restart.nc'

# output files
ocean2 = '/work/ollie/nkolduno/output_CORE_test2/fesom.1948.oce.restart.nc'
ice2   = '/work/ollie/nkolduno/output_CORE_test2/fesom.1948.ice.restart.nc'

fo = xr.open_dataset(ocean)
fi = xr.open_dataset(ice)

# Select time slice. We use slice to keep the time dimention
fos = fo.isel(time=slice(5,6), drop=False)
fis = fi.isel(time=slice(5,6), drop=False)

fos.to_netcdf(ocean2, unlimited_dims='time', format='NETCDF4_CLASSIC')
fis.to_netcdf(ice2  , unlimited_dims='time', format='NETCDF4_CLASSIC')
