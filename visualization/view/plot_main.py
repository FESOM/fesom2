
# coding: utf-8

# In[ ]:

get_ipython().magic('matplotlib notebook')
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')


# In[ ]:

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
result_path ='../results/'
runid='fesom'
str_id, depth	='u', 100
# specify records and year to read
records, year=np.linspace(1,1,1).astype(int), 1948


# In[ ]:

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


# In[ ]:

# read the model result from fesom.XXXX.oce.nc
data=read_fesom_slice(str_id, records, year, mesh, result_path, runid, ilev=ilev)


# In[ ]:

# Example 1: plot the 2D slice from data at depth
# set the label for the colorbar & contour intervals for ftriplot
cbartext, cont = '$^\circ$C', [-3., 25., .1]
cmap=cmo.balance
fig = plt.figure(figsize=(12,8))
# ftriplot is defined in fesom_plot_tools.py
data[data==0]=np.nan
[im, map, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]), oce='global', cmap=cmap)
cbar.set_label(cbartext, fontsize=22)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
cbar.ax.tick_params(labelsize=22)


# In[ ]:

# Example 2: comparing to climatology
# read the climatology first
from climatology import *
phc=climatology('phc3/', climname='phc')


# In[ ]:

# interpolate the data onto the climatology grig
[iz, xx, yy, zz]=fesom_2_clim(data, depth, mesh, phc)


# In[ ]:

# plot the difference to climatology
cbartext=('$^\circ$ C' if (str_id=='temp') else 'psu')
cont = [-2., 2., .1]
# compute the difference to climatology (only T & S are supported)
dd=zz[:,:]-(phc.T[iz,:,:] if (str_id=='temp') else phc.S[iz,:,:])

fig = plt.figure(figsize=(12,8))
[im, map, cbar]=wplot_xy(xx,yy,dd,np.arange(cont[0], cont[1]+cont[2], cont[2]), cmap=cmo.balance, do_cbar=True)
cbar.set_label(cbartext, fontsize=22)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
cbar.ax.tick_params(labelsize=22)


# In[ ]:

# Example 3: plot the norm of velocity (given on elements)
# read the ocean velocities
u=read_fesom_slice('u', records, year, mesh, result_path, runid, ilev=ilev)
v=read_fesom_slice('v', records, year, mesh, result_path, runid, ilev=ilev)


# In[ ]:

# interpolate them onto nodes
unodes=np.zeros(shape=mesh.n2d)
vnodes=np.zeros(shape=mesh.n2d)

var_elem=u*mesh.voltri
for i in range(mesh.e2d):
    unodes[mesh.elem[i,:]]=unodes[mesh.elem[i,:]]+[var_elem[i], var_elem[i], var_elem[i]]
unodes=unodes/mesh.lump2/3.

var_elem=v*mesh.voltri
for i in range(mesh.e2d):
    vnodes[mesh.elem[i,:]]=vnodes[mesh.elem[i,:]]+[var_elem[i], var_elem[i], var_elem[i]]
vnodes=vnodes/mesh.lump2/3.


data=np.hypot(unodes, vnodes)


# In[ ]:

# Example 3: plot the velocity norm
# set the label for the colorbar & contour intervals for ftriplot
cbartext, cont = 'm/s', [0.,.5, 0.01]
cmap=cmo.amp
fig = plt.figure(figsize=(12,8))
# ftriplot is defined in fesom_plot_tools.py
data[data==0]=np.nan
[im, map, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]), oce='global', cmap=cmap)
cbar.set_label(cbartext, fontsize=22)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
cbar.ax.tick_params(labelsize=22)

