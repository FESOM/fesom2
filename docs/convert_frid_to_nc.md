Convert grid to netCDF that CDO understands
===========================================

We are going to use spheRlab for conversion. You have to have R already installed.

Clone spheRlab

```bash
$bash> git clone https://github.com/FESOM/spheRlab.git spheRlab
```
Buiild package:
```bash
$bash>cd spheRlab/
$bash>R CMD build spheRlab
```
Make sure you have cdo installed (`cdo -V`) and launch R (type `R`)
#make sure that you have cdo installed:
#cdo -V
R>install.packages("spheRlab_1.1.0.tar.gz",repos=NULL)
R>library(spheRlab)
R>?sl.grid. #will give you help
R>meshpath="/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/"
R>grid = sl.grid.readFESOM(griddir=meshpath,rot=TRUE,rot.invert=TRUE,rot.abg=c(50,15,-90))
R>ofile = paste0(meshpath, "sl.grid.CDO")
R>sl.grid.writeCDO(grid,ofile=ofile)
R>system(paste0("cdo -f nc const,0,",ofile," ",ofile,".nc"))
R>sl.grid.addinfo(grid,ncdf.file.in=paste0(ofile,".nc")) #This step might be ignored, it might also not work if netcdf libraries are not loaded.

$bash> export MESHPATH=/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/
$bash> export DATAPATH=/work/ollie/dsidoren/ETOPO5/etopo5_lonlat.nc
$bash> cdo remapycon,$MESHPATH/sl.grid.CDO.nc -selname,topo $DATAPATH $MESHPATH/topo.nc