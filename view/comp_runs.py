#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(0,11,12).astype(int), [1965, 1965]

lev=20
#######################################
str_id, depth	='ssh', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-.1, .1, .001]
#	read the model result from fesom.XXXX.oce.nc
result_path='../results/CORE2_0001/'
data1	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
result_path='../results/CORE2_005/'
data2	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
data=data2-data1
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
str_id, depth	='temp', 100
#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'degree C', [-2, 25., .1]
cbartext, cont	= 'degree C', [-2., 2., .01]
#	read the model result from fesom.XXXX.oce.nc
result_path='../results/CORE2_0001/'
data1	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
result_path='../results/CORE2_005/'
data2	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data=data2-data1
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
str_id, depth	='salt', 100
#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'psu', [29., 37., .1]
cbartext, cont	= 'psu', [-1., 1., .01]
#	read the model result from fesom.XXXX.oce.nc
result_path='../results/CORE2_0001/'
data1	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
result_path='../results/CORE2_005/'
data2	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
data=data2-data1
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.show(block=False)



