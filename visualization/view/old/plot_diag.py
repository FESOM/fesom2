#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(0,11,12).astype(int), [1952, 1952]

#######################################
str_id, depth	='ssh', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-2., 2., .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+', '+'year ' +str(years[0]))
#######################################
str_id, depth	='temp', 10
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'deg C', [-2., 30., .1]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev)
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
#######################################
str_id, depth	='salt', 10
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'PSU', [29., 37., .1]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev)
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
#######################################
str_id, depth	=['u', 'v'], 0
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm/s', [0., .5, .01]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data1	=read_fesom_3del(str_id[0], months, years, mesh, result_path, runid, '.oce.nc', lev)
data2	=read_fesom_3del(str_id[1], months, years, mesh, result_path, runid, '.oce.nc', lev)
data=np.sqrt(np.square(data1)+np.square(data2))
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]), data_on_elem=1)
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))

############################# PLOT TRENDS #############################
str_id, depth	='ssh', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-.2, .2, .01]
#	read the model result from fesom.XXXX.oce.nc
data1	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
data2	=read_fesom_2d(str_id, months, [years[0]-1, years[1]-1], mesh, result_path, runid, '.oce.nc')
data=data1-data2
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
data[17848-1]=5
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' trend, '+'year ' +str(years[0]))
#######################################
str_id, depth	='temp', 10
#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'm', [-2, 30., .1]
cbartext, cont	= 'deg C', [-1., 1., .01]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data1	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev)
data2	=read_fesom_3d(str_id, months, [years[0]-1, years[1]-1], mesh, result_path, runid, '.oce.nc', lev)
data=data1-data2
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
plt.title('FVOM '+str(str_id)+' trend at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
#######################################
str_id, depth	='salt', 10
#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'm', [29., 37., .1]
cbartext, cont	= 'PSU', [-1., 1., .01]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data1	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev)
data2	=read_fesom_3d(str_id, months, [years[0]-1, years[1]-1], mesh, result_path, runid, '.oce.nc', lev)
data=data1-data2
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' trend at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
plt.show()

