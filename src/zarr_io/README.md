# Zarr-IO notes

## To checkout branch and pull zarr-fortran submodule
git clone https://github.com/FESOM/fesom2.git fesom_zarr
cd fesom_zarr;
git checkout suvi_zarr_output
git submodule update --init

... usual compilation steps using `./configure *options*`
## new io options

To facilitate some benchmarking with zarr-io with existing netcdf-io. There is a option to disable/enable netcdf-io and zarr-io for output using
boolean flags in namelist.io,&io_ctl  `lnetcdf_io` , `lzarr_io`

`namelist.io`
```
....
....
!using just one variable to be able to compare netcdf and zarr and of same kind(8)
&nml_list
io_list =  'sst       ',1, 's', 8,
/
&io_ctl
lnetcdf_io = .false.
lzarr_io = .true.
/

example namelist.io file is in `results/namelist.io`
```

Note: 
1. Currently all IO is of real, kind=8, inline with kind used by a variable in memory, so for benchmarking it may be better to set netcdf io to also use the same kind like in namelist.io example above. While this saves some model memory will increase storage size.(This can  be fixed easily in future) 
2. output is in `test2.zarr` (a hard coded name) in workdir/results dir. Post-processing code to load the data in Python will follow soon.
3. current zarr-io implementation is syncronous IO (but parallel), unlike async-IO of netcdf using threads in fesom2. zarr-io subroutine output_zarr is called in fvom_main.F90. To add a new zarr output variable, see line 56 in  src/zarr_io/io_zarrdata.F90, adding another call to write_array should be enough. ( This can also be fixed easily in future in lines of what was done for hecuba/cassandra io) 

