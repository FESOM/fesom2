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

def fesom_2_woa2005(data, mesh, woa05):
	xx,yy = np.meshgrid(woa05.x, woa05.y)
	zz=np.copy(woa05.T)
	for dep_ind in range(len(woa05.z)):
		print 'interpolating level: ', dep_ind
		wdep=woa05.z[dep_ind]	
		dep_up=[z for z in abs(mesh.zlevs) if z<=wdep][-1]
		dep_lo=[z for z in abs(mesh.zlevs) if z>wdep][0]
		i_up=1-abs(wdep-dep_up)/(dep_lo-dep_up)
		i_lo=1-abs(wdep-dep_lo)/(dep_lo-dep_up)
		data2=i_up*fesom2depth(dep_up, data, mesh)
		data2=data2+i_lo*fesom2depth(dep_lo, data, mesh)		
		zz[dep_ind,:,:] = griddata((mesh.x2,mesh.y2), data2, (xx,yy), method='linear')
	zz[np.isnan(woa05.T)]=np.nan
	return xx, yy, zz
