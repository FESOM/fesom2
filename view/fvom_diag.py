#	modify include.py and init.py before running this script
# 	include required system & FESOM modules
execfile("include.py")
# 	define the paths where the mesh & data are stored
execfile("init.py")
#######################################

d=mesh.x2[mesh.elem].max(axis=1)-mesh.x2[mesh.elem].min(axis=1)
ind=(d < 2)

n=1190238-1
elvls=np.loadtxt(mesh.path+'elvls.out')

ind=[];
#indexies for the box
#box=[-20, 40, 40, 80]
box=[mesh.x2[n]-3, mesh.x2[n]+3, mesh.y2[n]-3, mesh.y2[n]+3]

i1=mesh.x2[mesh.elem].max(axis=1)>=box[0]
i2=mesh.x2[mesh.elem].max(axis=1)<=box[1]
i3=mesh.y2[mesh.elem].max(axis=1)>=box[2]
i4=mesh.y2[mesh.elem].max(axis=1)<=box[3]

ind=i1&i2&i3&i4


fig = plt.figure()
#ax=plt.triplot(mesh.x2, mesh.y2, mesh.elem[ind,:], color='b', facecolors=elvls[ind])
#elvls[elvls>25]=25
#elvls[elvls<15]=15
ax=plt.tripcolor(mesh.x2, mesh.y2, mesh.elem[ind,:], color='w', facecolors=elvls[ind])
plt.grid()
fig.gca().set_xlim(box[0:2])
fig.gca().set_ylim(box[2:4])
plt.grid()
#cbar = plt.colorbar(ax, ticks=np.arange(elvls[ind].min(),elvls[ind].max()+1,1).astype(int), orientation='horizontal')
cbar = plt.colorbar(ax, orientation='horizontal') #ticks=np.arange(5,12,1).astype(int)
plt.plot(mesh.x2[n], mesh.y2[n], markersize=10, marker='o')
plt.show(block=False)

