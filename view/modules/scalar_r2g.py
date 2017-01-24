def scalar_r2g(al, be, ga, rlon, rlat):
	import numpy as np
	import math as mt
#	from numpy import *
#	from math import *
	result = []
	rad=mt.pi/180
	al=al*rad
	be=be*rad
	ga=ga*rad
	rotate_matrix=np.zeros(shape=(3,3))
	rotate_matrix[0,0]=np.cos(ga)*np.cos(al)-np.sin(ga)*np.cos(be)*np.sin(al)
	rotate_matrix[0,1]=np.cos(ga)*np.sin(al)+np.sin(ga)*np.cos(be)*np.cos(al)
	rotate_matrix[0,2]=np.sin(ga)*np.sin(be)
	rotate_matrix[1,0]=-np.sin(ga)*np.cos(al)-np.cos(ga)*np.cos(be)*np.sin(al)
	rotate_matrix[1,1]=-np.sin(ga)*np.sin(al)+np.cos(ga)*np.cos(be)*np.cos(al)
	rotate_matrix[1,2]=np.cos(ga)*np.sin(be)
	rotate_matrix[2,0]=np.sin(be)*np.sin(al)
	rotate_matrix[2,1]=-np.sin(be)*np.cos(al)
	rotate_matrix[2,2]=np.cos(be)
	rotate_matrix=np.linalg.pinv(rotate_matrix)
	rlat=rlat*rad
	rlon=rlon*rad	#Rotated Cartesian coordinates:
	xr=np.cos(rlat)*np.cos(rlon)
	yr=np.cos(rlat)*np.sin(rlon)
	zr=np.sin(rlat)	#Geographical Cartesian coordinates:
	xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr
	yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr
	zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr		#Geographical coordinates:
	lat=[mt.asin(val) for (i, val) in enumerate(zg)]
	lon=[mt.atan2(yg[i],xg[i]) for i in range(len(xg))]
	a = [i for (i, val) in enumerate((abs(xg)+abs(yg))) if val ==0]
	#.astype(float32)
	if a: lon[a]=0
	lat = [lat[i]/rad for i in range(len(lat))]
	lon = [lon[i]/rad for i in range(len(lon))]
	lat=np.array(lat)
	lon=np.array(lon)
	return (lon,lat)
