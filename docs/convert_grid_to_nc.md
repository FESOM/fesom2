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
Load library:
```R
R>library(spheRlab)
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
```R
R>grid = sl.grid.readFESOM(griddir=meshpath,rot=TRUE,rot.invert=TRUE,rot.abg=c(50,15,-90))
```
Define path to the output file:
```R
R>ofile = paste0(meshpath, "sl.grid.CDO")
```
Write ascii mesh:
```R
R>sl.grid.writeCDO(grid,ofile=ofile)
```
Convert ascii to NetCDF via system call:
```R
R>system(paste0("cdo -f nc const,0,",ofile," ",ofile,".nc"))
```
You could end here. However, `sl.grid.readFESOM` has generated other useful information on the grid that can be added to the NetCDF mesh file now as follows:
```R
R>install.packages("ncdf4") # only needed if not yet installed
R>sl.grid.addinfo(grid,ncdf.file.in=paste0(ofile,".nc"))
```

Conservative remapping with cdo (interpolate topography to mesh)
---------------------------------------------------------------
```bash
$bash> export MESHPATH=/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/
$bash> export DATAPATH=/work/ollie/dsidoren/ETOPO5/etopo5_lonlat.nc
$bash> cdo remapycon,$MESHPATH/sl.grid.CDO.nc -selname,topo $DATAPATH $MESHPATH/topo.nc
```