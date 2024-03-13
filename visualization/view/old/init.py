# set the paths to the mesh & results
meshpath  ='/work/ollie/dsidoren/input/fesom2.0/meshes/mesh_CORE2_final/'
result_path ='../results/core2_nsw/'
runid='fesom'

# read the mesh
alpha, beta, gamma=[50, 15, -90] # rotated mesh
#alpha, beta, gamma=[0., 0., 0.]  # not rotated mesh
try:
	mesh
except NameError:
	print("mesh will be loaded")
	mesh=read_fesom_mesh(meshpath, [alpha, beta, gamma])
else:
	print("mesh with this name already exists and will not be loaded")
