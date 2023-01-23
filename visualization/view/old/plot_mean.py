#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
exec(open("include.py").read())
# 	define the paths where the mesh & data are stored
exec(open("init.py").read())

######################################
str_id, depth	='hice', 100
#	read the model result from fesom.XXXX.oce.nc
months, years=np.linspace(0,2,3).astype(int), [1948, 1948]
h_march=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.ice.nc') #, ncfile='fvsom.restart.oce.nc')

months, years=np.linspace(8,10,3).astype(int), [1948, 1948]
h_sept=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.ice.nc')#, ncfile='fvsom.restart.oce.nc')

#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [0, 3, .05]
# plot the northern ice
contours_abs=[0, 5]
fig = plt.figure(figsize=(11.5,5))
fig.subplots_adjust(wspace=0,hspace=0)
# left subplot (ice in March)
ax1 = fig.add_subplot(1,2,1)

[im, cbar]=ftriplot(mesh, h_march, np.arange(contours_abs[0], contours_abs[1]+.1, .1), oce='np', do_cbar=False, mlabels=[1,0,0,0], extend='max')

# right subplot (ice in September)
ax2 = fig.add_subplot(1,2,2)
[im, cbar]=ftriplot(mesh, h_sept, np.arange(contours_abs[0], contours_abs[1]+.1, .1), oce='np', do_cbar=False, extend='max')

cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
cbar=fig.colorbar(im, cax, orientation='vertical')
cbar.set_ticks([round(i,1) for i in np.linspace(contours_abs[0], contours_abs[1], 7)])
cbar.set_label('m')

#######################################
str_id, depth	='hice', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [0, 3, .05]
# plot the southern ice
contours_abs=[0, 3]
fig = plt.figure(figsize=(11.5,5))
fig.subplots_adjust(wspace=0,hspace=0)
# left subplot (ice in March)
ax1 = fig.add_subplot(1,2,1)
h=np.copy(h_march)
h[mesh.y2>0]=np.nan
[im, cbar]=ftriplot(mesh, h, np.arange(contours_abs[0], contours_abs[1]+.02, .05), oce='sp', do_cbar=False, mlabels=[1,0,0,0], extend='max')


ax2 = fig.add_subplot(1,2,2)
h=np.copy(h_sept)
h[mesh.y2>0]=np.nan
[im, cbar]=ftriplot(mesh, h, np.arange(contours_abs[0], contours_abs[1]+.1, .1), oce='sp', do_cbar=False, mlabels=[1,0,0,0], extend='max')

cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
cbar=fig.colorbar(im, cax, orientation='vertical')
cbar.set_ticks([round(i,1) for i in np.linspace(contours_abs[0], contours_abs[1], 7)])
cbar.set_label('m')

#######################################
str_id, depth	='ssh', 100
months, years=np.linspace(0,11,12).astype(int), [1948, 1948]
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [-2, 2, .1]
#	read the model result from fesom.XXXX.oce.nc
data	=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.oce.nc')
#       plot the stuff
fig = plt.figure(figsize=(8,5))
#       ftriplot is defined in fesom_plot_tools.py
[im, cbar]=ftriplot(mesh, data, np.arange(cont[0], cont[1]+cont[2], cont[2]))
cbar.set_label(cbartext)
cbar.set_ticks([round(i,4) for i in np.linspace(cont[0], cont[1], 5)])
#######################################
