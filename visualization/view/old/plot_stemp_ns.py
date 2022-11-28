#sys.exit("!stop!")
str_id, depth	='temp', 100
#	set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	= 'm', [0, 3, .05]
#	read the model result from fesom.XXXX.oce.nc
months, years=np.linspace(2,2,1).astype(int), [1950, 1950]
h_march=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc')

months, years=np.linspace(8,8,1).astype(int), [1950, 1950]
h_sept=read_fesom_3d(str_id, months, years, mesh, result_path, runid, '.oce.nc')

# plot the northern ice
contours_abs=[-3, 3]
fig = plt.figure(figsize=(11.5,5))
fig.subplots_adjust(wspace=0,hspace=0)
# left subplot (ice in March)
ax1 = fig.add_subplot(1,2,1)

[im, cbar]=ftriplot(mesh, h_march, np.arange(contours_abs[0], contours_abs[1]+.1, .1), oce='np', do_cbar=False, mlabels=[1,0,0,0], extend='both')

# right subplot (ice in September)
ax2 = fig.add_subplot(1,2,2)
[im, cbar]=ftriplot(mesh, h_sept, np.arange(contours_abs[0], contours_abs[1]+.1, .1), oce='np', do_cbar=False, extend='both')

cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
cbar=fig.colorbar(im, cax, orientation='vertical')
cbar.set_ticks([round(i,1) for i in np.linspace(contours_abs[0], contours_abs[1], 7)])
cbar.set_label('m')
