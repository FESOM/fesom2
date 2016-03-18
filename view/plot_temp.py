#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
years=[1975, 1975]
months=np.linspace(5,5,1).astype(int)
str_id = 'temp'
lev=10 # the depth
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'degree C', [-2, 27., .1]
#######################################
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



w=woa2005('woa2005/')
str_id	= 'temp'
cbartext, cont	=  '$^\circ C$', [-2, 2, .01]
depth=100.
wlev=np.argmin(abs(w.z-depth))
xx,yy,zz=fesom_2_woa2005(str_id, months, years, mesh, result_path, runid, '.oce.nc', w)
#	plot the model difference to climatology
fig = plt.figure(figsize=(8,5))
[im, cbar]=wplot_xy(xx,yy,zz[wlev,:,:]-w.S[wlev,:,:], np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_ticks([round(i,2) for i in np.linspace(cont[0], cont[1], 5)])
cbar.set_label(cbartext)


str_id	= 'salt'
cbartext, cont	= 'psu', [-.5, .5, .01]
depth=100.
wlev=np.argmin(abs(w.z-depth))
xx,yy,zz=fesom_2_woa2005(str_id, months, years, mesh, result_path, runid, '.oce.nc', w)
#	plot the model difference to climatology
fig = plt.figure(figsize=(8,5))
[im, cbar]=wplot_xy(xx,yy,zz[wlev,:,:]-w.S[wlev,:,:], np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_ticks([round(i,2) for i in np.linspace(cont[0], cont[1], 5)])
cbar.set_label(cbartext)

plt.show(block=False)




