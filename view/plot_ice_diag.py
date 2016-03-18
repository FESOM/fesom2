#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")

str_id, years='area', [1971, 1975]
hemi='np' #or 'sp'
# set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [0., 1., .01] #concentration
# read March result from fesom.XXXX.ice.nc
months=np.linspace(2,2,1).astype(int)
h_march=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.ice.nc')

# read September result from fesom.XXXX.ice.nc
months=np.linspace(8,8,1).astype(int)
h_sept=read_fesom_2d(str_id, months, years, mesh, result_path, runid, '.ice.nc')
if (hemi=='sp'):
        h_sept[mesh.y2>0.]=0.

############ Plot the ice ############
fig = plt.figure(figsize=(11.5,5))
fig.subplots_adjust(wspace=0,hspace=0)
# left subplot (ice in March)
ax1 = fig.add_subplot(1,2,1)
[im, cbar]=ftriplot(mesh, h_march, np.arange(cont[0], cont[1]+cont[2], cont[2]), oce=hemi, do_cbar=False, mlabels=[1,0,0,0], extend='max')

# right subplot (ice in September)
ax2 = fig.add_subplot(1,2,2)
[im, cbar]=ftriplot(mesh, h_sept, np.arange(cont[0], cont[1]+cont[2], cont[2]), oce=hemi, do_cbar=False, extend='max')

cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
cbar=fig.colorbar(im, cax, orientation='vertical')
cbar.set_ticks([round(i,1) for i in np.linspace(contours_abs[0], contours_abs[1], 7)])
cbar.set_label(cbartext)
