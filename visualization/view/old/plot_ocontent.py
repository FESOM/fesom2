#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#	set the years & months for the average
months=np.linspace(0,11,12).astype(int)
years=[1948, 1977]


str_id ='salt'
dimZ=len(mesh.zlev)-1
layer_content=np.zeros([len(years), dimZ])
layer_volume =np.zeros(dimZ)

for I, Y in enumerate(years):
	for K in np.arange(1, dimZ, 1):	
		dz=np.diff(mesh.zlev[[K-1, K]])
		data	=read_fesom_3d(str_id, months, [Y, Y], mesh, result_path, runid, '.oce.nc', K)
		data[data==0]=np.nan
		el_data=data[mesh.elem].sum(axis=1)/3.
		layer_content[I, K]=np.nansum(el_data*mesh.voltri)*dz
		if (I==0):
			layer_volume[K]=sum(v for i, v in enumerate(mesh.voltri) if not np.isnan(el_data[i]))*dz


fig = plt.figure(figsize=(7, 7))

for I, Y in enumerate(years):
	plt.plot(layer_content[I,:]/layer_volume, mesh.zlev[0:dimZ],label='FVOM '+str(Y))

import sys
sys.exit()

agu03=np.loadtxt('agu03.dat')
if str_id =='temp':
	plt.plot(agu03[:, 0], -agu03[:, 2], label='FESOM 2010')

if str_id =='salt':
	plt.plot(agu03[:, 1], -agu03[:, 2], label='FESOM 2010')

plt.grid()

if str_id =='temp': 
	fig.gca().set_xlabel('degree C', fontsize=14)
if str_id =='salt':
	fig.gca().set_xlabel('PSU', fontsize=14)

fig.gca().set_ylabel('m',fontsize=14)
plt.legend(loc=4)
plt.show(block=False)
