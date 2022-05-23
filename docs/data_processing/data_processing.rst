.. _chap_data_processing:

Data pre/post processing
************************

netCDF files for initial conditions
===================================

The netCDF files have to satisfy the following criteria:

- should have DIMENSIONS named ``lon/longitude/LON`` , ``lat/latitude/LAT`` and ``depth``
- should have VARIABLES  named ``lon/longitude/LON`` , lat/latitude/LAT and ``depth``
- ``lon/lat`` dimentions should be one dimentional (e.g ``lon(lon)``)
- each variable with initial conditions should have only three dimentions (e.g. ``temp(depth, lat, lon)``)
- The fields should start from ``0th`` meridian and longitudes should have values from ``0 to 360``
- the missing values should have values larger than ``1e11``

The file that would be read potentially without problems can be created with the following python code (variables lat,lon_reshaped, depth, salt_reshaped, temp_reshaped should be prepeared from the original file):

.. code-block:: python

    from netCDF4 import Dataset

    fw = Dataset('woa18_netcdf.nc', 'w', )

    fw.createDimension('time', 1)
    fw.createDimension('lat', lat.shape[0])
    fw.createDimension('lon', lon_reshaped.shape[0])
    fw.createDimension('depth', depth.shape[0])

    latitude  = fw.createVariable('lat', 'd', ('lat',))
    latitude[:] = lat[:]

    longitude = fw.createVariable('lon', 'd', ('lon',))
    longitude[:] = lon_reshaped[:]

    ddepth = fw.createVariable('depth', 'd', ('depth',))
    ddepth[:] = depth[:]

    salinity = fw.createVariable('salt','d', ('depth', 'lat', 'lon'), fill_value= 1e+20)
    salinity[:] = salt_reshaped[:]
    salinity.missing_value = 1e+20

    temperature = fw.createVariable('temp','d', ('depth', 'lat', 'lon'), fill_value= 1e+20)
    temperature[:] = temp_reshaped[:]
    temperature.missing_value = 1e+20

We will try to provide convertion instructions in the form of jupyter notebooks to all files with initial conditions.


Convert grid to netCDF that CDO understands
===========================================

We are going to use ``spheRlab`` for conversion. You have to have **R** already installed.

Clone ``spheRlab``:

::

    git clone https://github.com/FESOM/spheRlab.git spheRlab

Build package:

::

    cd spheRlab/
    R CMD build spheRlab

Make sure you have cdo installed (``cdo -V``) and launch R (type ``R``).

Install the package:

::

    R>install.packages("spheRlab_1.1.0.tar.gz",repos=NULL)

If you don't have netCDF library installed, you also have to do:

::

    R>install.packages("ncdf4")

Load libraries:

::

    R>library(spheRlab)
    R>library(ncdf4)

You can get help (for any function) by typing, e.g.:

::

    R>?sl.grid.writeCDO

Define path to the mesh:

::

    R>meshpath="/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/"

Read the grid in to R structure (the arguments rot etc. might be different for different meshes, but this is the standard):

For rotated meshes:

::

    R>grid = sl.grid.readFESOM(griddir=meshpath,rot=TRUE,rot.invert=TRUE,rot.abg=c(50,15,-90))

For unrotated meshes:

::

    R>grid = sl.grid.readFESOM(griddir=meshpath,rot=FALSE,rot.invert=FALSE,rot.abg=c(0,0,0), threeD=FALSE)

Define path to the output file:

::

    R>ofile = paste0(meshpath, "sl.grid.CDO", sep = "")

Directrly write netCDF file with mesh description:

::

    R>sl.grid.writeCDO(grid, ofile=ofile, netcdf=TRUE, depth=FALSE)

Conservative remapping with cdo (interpolate topography to mesh)
----------------------------------------------------------------
::

    $bash> export MESHPATH=/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/
    $bash> export DATAPATH=/work/ollie/dsidoren/ETOPO5/etopo5_lonlat.nc
    $bash> cdo remapycon,$MESHPATH/sl.grid.CDO.nc -selname,topo $DATAPATH $MESHPATH/topo.nc


