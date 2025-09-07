#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(0,0,1).astype(int), [2010, 2010]

lev=10
#######################################
str_id = 'ssh'
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-1, 1, .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
str_id	= 'temp'
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'degree C', [-2, 25., .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
str_id ='salt'
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'psu', [29., 37., .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.show(block=False)
#######################################
lev=11
str_id ='temp'
months, years=np.linspace(0,11,12).astype(int), [1948, 1948]
data1=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)

months, years=np.linspace(0,11,12).astype(int), [2010, 2010]
data2=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data=data2-data1
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'degree C', [-5., 5., .01]
#	read the model result from fesom.XXXX.oce.nc
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
lev=25
str_id ='salt'
months, years=np.linspace(0,11,12).astype(int), [1948, 1948]
data1=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)

months, years=np.linspace(0,11,12).astype(int), [2010, 2010]
data2=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data=data2-data1
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'psu', [-.5, .5, .01]
#	read the model result from fesom.XXXX.oce.nc
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
