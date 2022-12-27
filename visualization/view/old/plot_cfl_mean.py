#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#######################################
#months, years=np.linspace(171,171,1).astype(int), [1948, 1948]
#str_id, lev	='cfl_w', 10

months, years=np.linspace(0,11,12).astype(int), [1954, 1954]
str_id, lev	='w', 10

#	set the label for the colorbar & contour intervals for ftriplot
#cbartext, cont	= 'm', [-.1, .1, .01]
cbartext, cont	= 'm', [-.00001, .00001, .000001]
cbartext, cont	= 'm', [-.0001, .0001, .000001]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc', lev)
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
#cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
