#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(1,1,1).astype(int), [1948, 1948]

#######################################
str_id, depth	='w', 10
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= ' ', [0., 1., .01]
#cbartext, cont	= 'PSU', [0., 5., .01]
#	read the model result from fesom.XXXX.oce.nc
arr=abs(abs(mesh.zlev)-abs(depth))
lev=[i for i, j in enumerate(arr) if j == min(arr)]
data	=read_fesom_3d(str_id, months, years,   mesh, result_path, runid, '.oce.nc', lev, ncfile='fvsom.oce.nc')
data[data==0]=np.nan
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
plt.title('FVOM '+str(str_id)+' at '+str(depth)+'m'+' , '+'year ' +str(years[0]))
#######################################

