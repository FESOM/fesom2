#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(0,0,1).astype(int), [1949, 1949]

#######################################
str_id, depth	='ssh', 0
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-2., 2., .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc', ncfile='fvsom.restart.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
#data[122363-1]=5
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+', '+'year ' +str(years[0]))
#######################################
str_id, depth	='temp', 9
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'deg C', [-2., 30., .1]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
lev=0
data	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev, ncfile='fvsom.restart.oce.nc')
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
str_id, depth	='salt', 100. #60
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'PSU', [30., 37., .1]
#cbartext, cont	= 'PSU', [34.37, 34.8, .01]
#cbartext, cont	= 'PSU', [0., 5., .01]
#cbartext, cont	= 'PSU', [-10, 0.1, .1]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
lev=0
data	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev, ncfile='fvsom.restart.oce.nc')
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
#######################################
str_id, lev	='w', 0

#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'm', [-.1, .1, .01]
cbartext, cont	= 'm', [-.00001, .00001, .000001]
cbartext, cont	= 'm', [-.001, .001, .0001]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev, ncfile='fvsom.restart.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
