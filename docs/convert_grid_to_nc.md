Convert grid to netCDF that CDO understands
===========================================

We are going to use spheRlab for conversion. You have to have R already installed.

Clone spheRlab

```bash
git clone https://github.com/FESOM/spheRlab.git spheRlab
```
Build package:
```bash
cd spheRlab/
R CMD build spheRlab
```
Make sure you have cdo installed (`cdo -V`) and launch R (type `R`).
Install the package:
```R
R>install.packages("spheRlab_1.1.0.tar.gz",repos=NULL)
```
If you don;t have netCDF library installed, you also have to do:

```R
R>install.packages("ncdf4")
```

Load libraries:
```R
R>library(spheRlab)
R>library(ncdf4)
```
You can get help (for any function) by typing, e.g.:
```R
R>?sl.grid.writeCDO
```

Define path to the mesh
```R
R>meshpath="/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/"
```
Read the grid in to R structure (the arguments `rot` etc. might be different for different meshes, but this is the standard):
For rotated meshes:
```R
R>grid = sl.grid.readFESOM(griddir=meshpath,rot=TRUE,rot.invert=TRUE,rot.abg=c(50,15,-90))
```
For unrotated meshes:
```R
R>grid = sl.grid.readFESOM(griddir=meshpath,rot=FALSE,rot.invert=FALSE,rot.abg=c(0,0,0), threeD=FALSE)
```
Define path to the output file:
```R
R>ofile = paste0(meshpath, "sl.grid.CDO", sep = "")
```
Directrly write netCDF file with mesh description:
```R
R>sl.grid.writeCDO(grid, ofile=ofile, netcdf=TRUE, depth=FALSE)
```

Conservative remapping with cdo (interpolate topography to mesh)
---------------------------------------------------------------
```bash
$bash> export MESHPATH=/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/
$bash> export DATAPATH=/work/ollie/dsidoren/ETOPO5/etopo5_lonlat.nc
$bash> cdo remapycon,$MESHPATH/sl.grid.CDO.nc -selname,topo $DATAPATH $MESHPATH/topo.nc
```