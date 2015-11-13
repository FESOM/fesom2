# set the paths to the mesh & results
result_path='../results/'
meshpath  ='../input/mesh/mesh_CORE2_final/'

#result_path='../results/expl/'
#meshpath  ='../input/mesh/mesh_ref/'

runid='fvsom'

# read the mesh
alpha, beta, gamma=[50, 15, -90]
try:
	mesh
except NameError:
	print "mesh will be loaded"
	mesh=read_fesom_mesh(meshpath, alpha, beta, gamma,read_diag=False)
else:
	print "mesh with this name already exists and will not be loaded"
