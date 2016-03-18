# set the paths to the mesh & results
result_path='../results/'
meshpath  ='/gfs2/work/hbkdsido/input/mesh/glob_15km/'
meshpath  ='/gfs2/work/hbkdsido/input/mesh/mesh_CORE2_final/'
runid='fvsom'

# read the mesh
alpha, beta, gamma=[50, 15, -90] # rotated mesh
#alpha, beta, gamma=[0., 0., 0.]  # not rotated mesh
try:
	mesh
except NameError:
	print "mesh will be loaded"
	mesh=read_fesom_mesh(meshpath, alpha, beta, gamma,read_diag=False)
else:
	print "mesh with this name already exists and will not be loaded"
