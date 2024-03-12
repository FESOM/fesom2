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

Interpolate initial conditions from existing run
================================================

The easiest way to bring your results from one mesh to another (for example to use low resolution spinup to initialise high resolution run), and continue computations is to use the results of source simulation as initial condition for target simulation. The most common case is to interpolate sea ice conditions, that usually need few years to spinup, and default behavour of FESOM2, that puts 2 meter thick ice everywere temperature is below freesing is not very optimal, espetially for 21st century runs.

Interpolation of sea ice
------------------------

Preparation of interpolated files is `described in this notebook <https://github.com/FESOM/pyfesom2/blob/master/notebooks/extra_notebooks/initial_conditions_sea_ice_interpolation.ipynb>`_

Resulting files should be located in `ClimateDataPath`, that is set in `namelist.config`.

In the `nemelist.tra` one should edit the `&tracer_init2d` section, that by default looks like this::

    &tracer_init2d                                      
    n_ic2d   = 3                                        
    idlist   = 1, 2, 3                                  
    filelist = 'a_ice2.nc', 'm_ice2.nc', 'm_snow2.nc' 
    varlist  = 'a_ice', 'm_ice', 'm_snow'
    ini_ice_from_file=.false.
    /

You should set `ini_ice_from_file` to `.true.`, and provide list of files with interpolated ice area (`a_ice`), sea ice thickness (`m_ice`) and snow thickness over sea ice (`m_snow`).

Interpolation of TS
-------------------

Preparation of interpolated files is `described in this notebook <https://github.com/FESOM/pyfesom2/blob/master/notebooks/extra_notebooks/initial_conditions_TS_interpolation.ipynb>`_

Resulting files should be located in `ClimateDataPath`, that is set in `namelist.config`.

In the `nemelist.tra` one should edit the `&tracer_init3d` section, that by default looks like this::

    &tracer_init3d                           
    n_ic3d   = 2                             
    idlist   = 2, 1                           
    filelist = 'phc3.0_winter.nc', 'phc3.0_winter.nc' 
    varlist  = 'salt', 'temp'                 
    t_insitu = .true.
    /

If you interpolate from model data, most probably you have data in potential temperature, so switch `t_insitu` to `.false.`, and provide files with your initial conditions. If you follow the notebook example, and create two separate files for temperature and salinity, your configuration might look like this::

    &tracer_init3d                            
    n_ic3d   = 2                              
    idlist   = 2, 1                           
    filelist = 'phc3.0_winter.nc', 'phc3.0_winter.nc' 
    varlist  = 'salt', 'temp'                 
    t_insitu = .true.
    /

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


