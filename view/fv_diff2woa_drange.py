#  plot the difference to WOA2005 climatology integrated over different depth layers
#  modify include.py and init.py before running this script
#  include required system & FESOM modules
execfile("include.py")
# define the paths where the mesh & data are stored
execfile("init.py")
# set the years & months for the average
years=[1992, 2007]
months=np.linspace(0,11,12).astype(int)
str_id = 'temp'
str_id = 'salt'
# set the label for the colorbar & contour intervals for ftriplot
cbartext, cont	=  '$^\circ C$', [-2., 2., .01]
cbartext, cont	=  'psu', [-.5, .5, .002]
# Read WOA2005
w=woa2005('woa2005/')
w.T_insitu=np.copy(w.T) #backup the insitu temperature
for i in range(len(w.z)):
	w.T[i,:,:]=sw.temppot0(w.S[i,:,:], w.T[i,:,:], w.z[i])

# read FESOM2.0 output on WOA2005 mesh
xx,yy,zz=fesom_2_woa2005(str_id, months, years, mesh, result_path, runid, '.oce.nc', w)

# Do the plotting
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,20))
depth=[0., 200.]
dmap=0.
dz=0.
# compute the depth averaged difference over depth[0] and deth[1]
#####################################
for wlev in range(1, len(w.z)):
	if (abs(w.z[wlev-1])>=depth[0] and abs(w.z[wlev]) <= depth[1]):
		aux_up=(zz[wlev-1,:,:]-(w.T[wlev-1,:,:] if (str_id=='temp') else w.S[wlev-1,:,:]))
		aux_lo=(zz[wlev,:,:]  -(w.T[wlev,:,:]   if (str_id=='temp') else w.S[wlev,:,:]))
		dz_loc=w.z[wlev]-w.z[wlev-1]
		dmap=dmap+0.5*dz_loc*(aux_up+aux_lo)
		dz=dz+dz_loc
dmap=dmap/dz
# plot the difference
ax1 = plt.subplot(4,1,1)
[im, map]=wplot_xy(xx,yy,dmap, np.arange(cont[0], cont[1]+cont[2], cont[2]), do_cbar=False)
# repeat the same for other layers
#####################################
depth=[200., 500.]
dmap=0.
dz=0.
for wlev in range(1, len(w.z)):
	if (abs(w.z[wlev-1])>=depth[0] and abs(w.z[wlev]) <= depth[1]):
		aux_up=(zz[wlev-1,:,:]-(w.T[wlev-1,:,:] if (str_id=='temp') else w.S[wlev-1,:,:]))
		aux_lo=(zz[wlev,:,:]  -(w.T[wlev,:,:]   if (str_id=='temp') else w.S[wlev,:,:]))
		dz_loc=w.z[wlev]-w.z[wlev-1]
		dmap=dmap+0.5*dz_loc*(aux_up+aux_lo)
		dz=dz+dz_loc
dmap=dmap/dz

ax2 = plt.subplot(4,1,2)
[im, map]=wplot_xy(xx,yy,dmap, np.arange(cont[0], cont[1]+cont[2], cont[2]), do_cbar=False)
#####################################
depth=[500., 1000.]
dmap=0.
dz=0.
for wlev in range(1, len(w.z)):
	if (abs(w.z[wlev-1])>=depth[0] and abs(w.z[wlev]) <= depth[1]):
		aux_up=(zz[wlev-1,:,:]-(w.T[wlev-1,:,:] if (str_id=='temp') else w.S[wlev-1,:,:]))
		aux_lo=(zz[wlev,:,:]  -(w.T[wlev,:,:]   if (str_id=='temp') else w.S[wlev,:,:]))
		dz_loc=w.z[wlev]-w.z[wlev-1]
		dmap=dmap+0.5*dz_loc*(aux_up+aux_lo)
		dz=dz+dz_loc
dmap=dmap/dz
ax3 = plt.subplot(4,1,3)
[im, map]=wplot_xy(xx,yy,dmap, np.arange(cont[0], cont[1]+cont[2], cont[2]), do_cbar=False)
#####################################
depth=[1000., 1500.]
dmap=0.
dz=0.
for wlev in range(1, len(w.z)):
	if (abs(w.z[wlev-1])>=depth[0] and abs(w.z[wlev]) <= depth[1]):
		aux_up=(zz[wlev-1,:,:]-(w.T[wlev-1,:,:] if (str_id=='temp') else w.S[wlev-1,:,:]))
		aux_lo=(zz[wlev,:,:]  -(w.T[wlev,:,:]   if (str_id=='temp') else w.S[wlev,:,:]))
		dz_loc=w.z[wlev]-w.z[wlev-1]
		dmap=dmap+0.5*dz_loc*(aux_up+aux_lo)
		dz=dz+dz_loc
dmap=dmap/dz
ax4 = plt.subplot(4,1,4)
[im, map]=wplot_xy(xx,yy,dmap, np.arange(cont[0], cont[1]+cont[2], cont[2]), do_cbar=False)
#####################################
# plot the colorbar
fig.subplots_adjust(bottom=0.2) 
cbar_ax = fig.add_axes([0.255, 0.15, 0.5, 0.02])
cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', ticks=np.arange(cont[0], cont[1]+cont[2], np.diff(cont[0:2])/4.))
cbar.set_label(cbartext)
plt.show(block=False)

