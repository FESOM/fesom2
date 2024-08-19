#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#######################################
from scipy.io import netcdf
import matplotlib.ticker

ncfile = netcdf.netcdf_file(result_path+'fvom_1965_1967_transports_5deg.nc', 'r')	
lons    = ncfile.variables['lons'].data
lats    = ncfile.variables['lats'].data
deps    = ncfile.variables['deps'].data
T	= ncfile.variables['time'].data
moc	= ncfile.variables['moc'].data
boc	= ncfile.variables['boc'].data


# Meridional overturning (in Z coordinate)
data=moc.mean(axis=0)
contours_abs=[-20, 20]
d_abs=1
data[data<-900]=np.nan
data[data<=contours_abs[0]]=contours_abs[0]
data[data>=contours_abs[1]]=contours_abs[1]

yy, zz = np.meshgrid(lats, deps)
fig = plt.figure(figsize=(8,5))
im=plt.contourf(yy, zz, data, levels=np.arange(contours_abs[0],contours_abs[1]+d_abs,d_abs))
cbar=plt.colorbar(orientation='horizontal', aspect=35, pad=0.11)#shrink=0.5)
tick_locs   = np.linspace(contours_abs[0],contours_abs[1],5)
cbar.locator   = matplotlib.ticker.FixedLocator(tick_locs)
cbar.update_ticks()
cbar.set_label('Sv')
fig.gca().set_xlabel('lat, degree',fontsize=12)
fig.gca().set_ylabel('depth, [m]',fontsize=12)
plt.grid()
#plt.title('Meridional overturning streamfunction')
