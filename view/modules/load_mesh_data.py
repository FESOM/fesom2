class fesom_mesh:
	"""existing instances are: path, n2d, e2d, nlev, zlevs, x2, y2, elem, n32, no_cyclic_elem, alpha, beta, gamma"""
	def __init__(self):
		self.path=''
		self.n2d=0
		self.n3d=0
		self.e2d=0
		self.nlev=0
		self.zlevs= []
		self.x2= []
		self.y2= []
		self.elem= []
		self.n32=[]
		self.no_cyclic_elem=[]
		self.voltri=[]
		self.alpha=0
		self.beta=0
		self.gamma=0		
	def __str__(self):
		return "mesh path=%s" % self.path

def read_fesom_mesh(path, alpha, beta, gamma, read_diag=True): 
	import matplotlib.pyplot as plt
	from scalar_r2g import scalar_r2g
	import numpy as np
	from scipy.io import netcdf

	mesh=fesom_mesh()
	mesh.path=path
	mesh.alpha=alpha
	mesh.beta=beta
	mesh.gamma=gamma	

	nod2dfile=mesh.path+'nod2d.out'
	elm2dfile=mesh.path+'elem2d.out'
	aux3dfile=mesh.path+'aux3d.out'

	file_content = np.loadtxt(nod2dfile, comments='!', skiprows=1)
	mesh.x2=file_content[:,1]
	mesh.y2=file_content[:,2]
	mesh.n2d=len(mesh.x2)
	
	file_content = np.loadtxt(elm2dfile, comments='!', skiprows=1)
	mesh.elem=file_content.astype(int)-1
	mesh.e2d=np.shape(mesh.elem)[0]

	with open(aux3dfile) as f:
		mesh.nlev=int(f.next())
		mesh.zlev=np.array([f.next().rstrip() for x in xrange(mesh.nlev)]).astype(float)

	###########################################	
	#here we compute the volumes of the triangles
	#this should be moved into fesom generan mesh output netcdf file
	#
	r_earth=6367500.0
	rad=np.pi/180
	edx=mesh.x2[mesh.elem]
	edy=mesh.y2[mesh.elem]
	ed=np.array([edx, edy])

	jacobian2D=ed[:, :, 1]-ed[:, :, 0]
	jacobian2D=np.array([jacobian2D, ed[:, :, 2]-ed[:, :, 0]])
	for j in range(2):
		mind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val > 355]
		pind = [i for (i, val) in enumerate(jacobian2D[j,0,:]) if val < -355]
		jacobian2D[j,0,mind]=jacobian2D[j,0,mind]-360
		jacobian2D[j,0,pind]=jacobian2D[j,0,pind]+360

	jacobian2D=jacobian2D*r_earth*rad

	for k in range(2):
		jacobian2D[k,0,:]=jacobian2D[k,0,:]*np.cos(edy.mean(1)*rad)

	mesh.voltri=np.zeros(mesh.e2d)
	for k in range(mesh.e2d):
		mesh.voltri[k]=abs(np.linalg.det(jacobian2D[:,:,k]))/2
	# compute the 2D lump operator
	cnt=np.array((0,)*mesh.n2d)
	mesh.lump2=np.array((0.,)*mesh.n2d)
	for i in range(3):
		for j in range(mesh.e2d):
			n=mesh.elem[j,i]
			#cnt[n]=cnt[n]+1
			mesh.lump2[n]=mesh.lump2[n]+mesh.voltri[j]
	mesh.lump2=mesh.lump2/3.

	#here we read the 3D cluster volumes
	if (read_diag):
		f = netcdf.netcdf_file(mesh.path+'fesom.initial.mesh.diag.nc', 'r')
		mesh.cluster_vol3=f.variables['cluster_vol'].data
		mesh.cluster_vol2=f.variables['cluster_area'].data
		f.close()
	else:
		mesh.cluster_vol3=0
		mesh.cluster_vol2=0
	#we should rotate the mesh to the geographical coordinates
	(mesh.x2,mesh.y2)=scalar_r2g(mesh.alpha,mesh.beta,mesh.gamma,mesh.x2,mesh.y2)
	d=mesh.x2[mesh.elem].max(axis=1)-mesh.x2[mesh.elem].min(axis=1)
	mesh.no_cyclic_elem = [i for (i, val) in enumerate(d) if val < 100]
	return mesh

def read_fesom_3d(str_id, months, years, mesh, result_path, runid, ext, ind, how='mean', ncfile=''): 
	import numpy as np
	from scipy.io import netcdf
	str_id		=str_id
	ext		=ext
	runid		=runid
	years		=years
	months		=months
	result_path	=result_path

	y			=years[0]
	data3		=np.zeros(shape=(mesh.n2d))
	while y<=years[1]:
		print ['reading year '+str(y)+':']
		if ncfile=='':
			ncfile =result_path+runid+'.'+str(y)+ext
		else:
			ncfile =result_path+ncfile

		f = netcdf.netcdf_file(ncfile, 'r')
		if how=='mean':
			data3 = data3+f.variables[str_id].data[months,:,ind].mean(axis=0)
		elif how=='max':
			data3 = data3+f.variables[str_id].data[months,:,ind].max(axis=0)
		f.close()
		y=y+1
	data3=data3/(years[1]-years[0]+1)
	return data3

def read_fesom_3del(str_id, months, years, mesh, result_path, runid, ext, ind, how='mean', ncfile=''): 
	import numpy as np
	from scipy.io import netcdf
	str_id		=str_id
	ext		=ext
	runid		=runid
	years		=years
	months		=months
	result_path	=result_path

	y		=years[0]
	data3		=np.zeros(shape=(mesh.e2d))
	while y<=years[1]:
		print ['reading year '+str(y)+':']
		if ncfile=='':
			ncfile =result_path+runid+'.'+str(y)+ext
		else:
			ncfile =result_path+ncfile

		f = netcdf.netcdf_file(ncfile, 'r')
		if how=='mean':
			data3 = data3+f.variables[str_id].data[months,:,ind].mean(axis=0)
		elif how=='max':
			data3 = data3+f.variables[str_id].data[months,:,ind].max(axis=0)
		f.close()
		y=y+1
	data3=data3/(years[1]-years[0]+1)
	return data3

def read_fesom_2d(str_id, months, years, mesh, result_path, runid, ext, how='mean', ncfile=''): 
	import numpy as np
	from scipy.io import netcdf
	y=years[0]
	data2=np.zeros(shape=(mesh.n2d))
	while y<=years[1]:
		print ['reading year '+str(y)+':']
		if ncfile=='':
			ncfile =result_path+runid+'.'+str(y)+ext
		else:
			ncfile =result_path+ncfile

		f = netcdf.netcdf_file(ncfile, 'r')
		if how=='mean':
			data2 = data2+f.variables[str_id].data[months,:].mean(axis=0)
		elif how=='max':
			data2 = data2+f.variables[str_id].data[months,:].max(axis=0)
		f.close()
		y=y+1
	data2=data2/(years[1]-years[0]+1)
	return data2

def fesom2depth(depth, data3, mesh): 
	import numpy as np
	data2=np.zeros(shape=(mesh.n2d))
	dind=(abs(mesh.zlevs-depth)).argmin()
	ind_depth=mesh.n32[:,dind]-1
	ind_noempty=[i for i in range(mesh.n2d) if ind_depth[i]>=0]
	ind_empty=[i for i in range(mesh.n2d) if ind_depth[i]<0]
	data2[ind_noempty]=data3[ind_depth[ind_noempty]]
	data2[ind_empty]=np.nan
	return data2
