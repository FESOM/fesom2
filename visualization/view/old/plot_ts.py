#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(11,11,1).astype(int), [1948, 1948]
months, years=np.linspace(11,11,1).astype(int), [1948, 1948]

#######################################
str_id, depth	='ssh', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-1, 1, .1]
contours=np.arange(cont[0], cont[1]+cont[2], cont[2])
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
im=plt.tricontourf(mesh.x2, mesh.y2, mesh.elem, data, levels=contours, extend='both')
cbar=plt.colorbar(im,orientation="horizontal")
plt.title(str(str_id)+' snapshot in year ' +str(years[0]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
str_id, depth	='temp', 200
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'degC', [3., 10., .1]
contours=np.arange(cont[0], cont[1]+cont[2], cont[2])
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', 25)
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
im=plt.tricontourf(mesh.x2, mesh.y2, mesh.elem, data, levels=contours, extend='both')
cbar=plt.colorbar(im,orientation="horizontal")
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title(str(str_id)+' snapshot at depth '+ str(depth)+ 'm in year ' +str(years[0]))
plt.show(block=False)
