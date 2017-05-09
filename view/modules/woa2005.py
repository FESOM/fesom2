from scipy.io import netcdf
import numpy as np
import scipy as sc
from load_mesh_data import *
from scipy.interpolate import griddata
from scipy.stats.stats import nanmean
class woa2005:
	"""existing instances are: x,y,z,T,S
	example: w=woa2005(woa05_path)"""
	def __init__(self, woa05_path):
		ncfile= netcdf.netcdf_file(woa05_path+'t00an1.nc', 'r')
		self.T= np.copy(ncfile.variables['t00an1'].data[0,:,:,:])
		x=np.copy(ncfile.variables['lon'].data)
		x[x>180]=x[x>180]-360
		ind=[i[0] for i in sorted(enumerate(x), key=lambda x:x[1])]
		x=np.sort(x)
		self.x=x
		self.y=ncfile.variables['lat'].data
		self.z=ncfile.variables['depth'].data
		ncfile.close()
		self.T[:,:,:]=self.T[:,:,ind]
		ncfile=netcdf.netcdf_file(woa05_path+'s00an1.nc', 'r')
		self.S=np.copy(ncfile.variables['s00an1'].data[0,:,:,:])
		self.S[:,:,:]=self.S[:,:,ind]
		ncfile.close()
		self.T[self.T>90]=np.nan
		self.S[self.S>90]=np.nan
		self.Tyz=nanmean(self.T, 2)
		self.Syz=nanmean(self.S, 2)

def fesom_2_woa2005(str_id, months, years, mesh, result_path, runid, woa05):
	xx,yy = np.meshgrid(woa05.x, woa05.y)
	zz=np.copy(woa05.T)
# Preload the 3D whole field if fits into the memory
	data = read_fesom_3d_full(str_id, months, years, mesh, result_path, runid)
	for dep_ind in range(len(woa05.z)):
		print('interpolating level: ', dep_ind)
		wdep=woa05.z[dep_ind]
# find indicies and depth of the nodes below and above the target			
		iz_up=[(i,z) for (i,z) in enumerate(abs(mesh.zbars)) if z <= wdep]
		iz_lo=[(i,z) for (i,z) in enumerate(abs(mesh.zbars)) if z >  wdep]
		if (not iz_up):
			iz_up=[(0, 0.)]
		if (not iz_lo):
			iz_lo=[(len(mesh.zbars), wdep)]	# the wdep and not mesh.zbars[-1] is chosen in order to 
                                                        # extrapolate the last column with the weight of 1
		dep_up=iz_up[-1][1]
		dep_lo=iz_lo[0][1]
		i_up=iz_up[-1][0]
		i_lo=iz_lo[0][0]
		c_up=1.-abs(wdep-dep_up)/(dep_lo-dep_up)
		c_lo=1.-abs(wdep-dep_lo)/(dep_lo-dep_up)
		print "i_up, i_lo=", i_up, i_lo
		print "c_up, c_lo=", c_up, c_lo
# Slow but less memory consuming technique
#		data_up	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, ext, i_up)
#		data_lo	=read_fesom_3d(str_id, months, years, mesh, result_path, runid, ext, i_lo)
#		data2=c_up*data_up+c_lo*data_lo
# We preloaded the whole 3D field 'data' instead
		data2=c_up*data[:, i_up]+c_lo*data[:, i_lo]
		zz[dep_ind,:,:] = griddata((mesh.x2,mesh.y2), data2, (xx,yy), method='linear')
	
	zz[np.isnan(woa05.T)]=np.nan
	return xx, yy, zz
