## Initial conditions in FESOM2

Initial conditions are loaded from netCDF file (or files) that contain values on regular grid.

### In the namelist.oce

One have to specify several parameters in the `namelist.oce`. The section of the `namelist.oce` might look like this:

```bash
&oce_init3d                               ! initial conditions for tracers
n_ic3d   = 2                              ! number of tracers to initialize
idlist   = 1, 0                           ! their IDs (0 is temperature, 1 is salinity, etc.). The reading order is defined here!
filelist = 'phc3.0_winter.nc', 'phc3.0_winter.nc' ! list of files in ClimateDataPath to read (one file per tracer), same order as idlist
varlist  = 'salt', 'temp'                 ! variables to read from specified files
t_insitu = .true.                         ! if T is insitu it will be converted to potential after reading it
```

* **n_ic3d** - how many tracers you want to initialise. As a minimum you have to initialise temperature and salinity. They have reserved id's 0 and 1 respectivelly (see `idlist` below). In this example only two tracers (temperature and salinity) are initialised.
* **idlist** - IDs of tracers. In this variable you define the reading order. Here we on purpose change the order of temperature and salinity to demonstrate that the order can be arbitrary. Once again remember that ID 0 and 1 are reserved for temperature and salinity respectively.
* **filelist** - coma separated list of files (each in qoutation marks) that contain initial conditions (see the section below about requirements to the file format). The path to the folder with this files is defined in `namelist.config` (`ClimateDataPath` variable). In this case the file `phc3.0_winter.nc` is the same for temperature and salinity since it contains both variables.
* **varlist** - names of the variables in the netCDF files specified above. Note again the order of the variables in the example, it can be arbitrary and in this case temperature comes after salinity. 
* **t_insitu** - most of climatologies are distributed with *in situ* temperature, while model needs potential temperature. This flag allows to do the conversion (UNESCO equation) on the fly.

### netCDF files for initial conditions

The netCDF files have to satisfy the following criteria:

* shoudl have dimentions named lon/longitude/LON , lat/latitude/LAT and depth
* should have variables  named lon/longitude/LON , lat/latitude/LAT and depth
* lon/lat dimentions should be one dimentional (e.g `lon(lon)`)
* each variable with initial conditions should have only three dimentions (e.g. `temp(depth, lat, lon)`)
* The fields should start from 0th meridian and longitudes should have values from 0 to 360
* the missing values should have values larger than 1e11

The file that would be read potentially without problems can be created with the following python code (variables `lat`,`lon_reshaped`, `depth`, `salt_reshaped`, `temp_reshaped` should be prepeared from the original file):

```python
from netCDF4 import Dataset

fw = Dataset('woa18_netcdf.nc', 'w', )
    
fw.createDimension('time', 366)
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
```

We will try to provide convertion instructions in the form of jupyter notebooks to all files with initial conditions.

