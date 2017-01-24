#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months, years=np.linspace(0,0,1).astype(int), [0, 0]

#######################################
str_id, depth	='ssh', 0
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [0., .1, .005]
#	read the model result from fesom.XXXX.oce.nc
#data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc', ncfile='fvsom.1969_1974.ssh.timvar.nc')
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc', ncfile='agu02.1981_1986.ssh.timvar.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
#data[122363-1]=5
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#plt.title('FVOM '+str(str_id)+', '+'year ' +str(years[0]))
plt.show(block=False)
