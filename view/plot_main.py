# import standard python packages
import sys
import numpy as np
# import basemap
from mpl_toolkits.basemap import Basemap
# import FESOM packages
sys.path.append("./modules/")
from load_mesh_data import *
sys.path.append("/home/h/hbkdsido/utils/seawater-1.1/")
import seawater as sw
from fesom_plot_tools import *
import cmocean.cm as cmo

# set the paths to the mesh & results
result_path ='../results/zlevel/withpice/'
runid='fesom'
str_id, depth	='sst', 0
# specify records and year to read
records, year=np.linspace(55,55,1).astype(int), 1948

# read the mesh
# set the path to the mesh
meshpath  ='/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/'
alpha, beta, gamma=[50, 15, -90]
try:
	mesh
except NameError:
	print("mesh will be loaded")
	mesh=load_mesh(meshpath, abg=[alpha, beta, gamma], usepickle = False)
else:
	print("mesh with this name already exists and will not be loaded")

#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, records, year, mesh, result_path, runid)

# Plot the field
# set the label for the colorbar & contour intervals for ftriplot
cbartext, cont = '$^\circ$C', [-3., 25., .1]
cmap=cmo.balance
fig = plt.figure(figsize=(12,8))
# ftriplot is defined in fesom_plot_tools.py
[im, map, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]), oce='global', cmap=cmap)
cbar.set_label(cbartext, fontsize=22)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
cbar.ax.tick_params(labelsize=22)


