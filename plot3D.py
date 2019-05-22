from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def compute_poly_AI(X,Y):
	##Area Calc
	A = 0
	for x in range(1,len(X)-2):
		A = A + (X[x]*Y[x+1] - X[x+1]*Y[x])
	A = 0.5*A

	##Centroid Calc
	Cx = 0
	Cy = 0
	for x in range(1,len(X)-2):
		Cx=Cx + (X[x] + X[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
		Cy=Cy + (Y[x] + Y[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
	Cx = Cx/(6*A)
	Cy = Cy/(6*A)
	
	Ix =0
	Iy=0
	Ixy=0
	for x in range(1,len(X)-2):
		Ix = Ix +(Y[x]**2 + Y[x]*Y[x+1] +Y[x+1]**2)*(X[x]*Y[x+1] - X[x+1]*Y[x])
		Iy = Iy + (X[x]**2 + X[x]*X[x+1] +X[x+1]**2)*(X[x]*Y[x+1] - X[x+1]*Y[x])
		Ixy = Ixy + (X[x]*Y[x+1] + 2*X[x]*Y[x] + 2*X[x+1]*Y[x+1] + X[x+1]*Y[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
	Ix =Ix/12
	Iy = Iy/12
	Ixy = Ixy/24
	return A,Ix,Iy,Ixy,Cx,Cy

def rotate_around_point(x, y, radians, x_pt,y_pt):

	adjusted_x = np.subtract(x,x_pt)
	adjusted_y = np.subtract(y,y_pt)

	qx = x_pt + np.cos(radians) * adjusted_x + np.sin(radians) * adjusted_y
	qy = y_pt + -np.sin(radians) * adjusted_x + np.cos(radians) * adjusted_y
	return qx, qy
	
def plotter(X,Y,twist,layers):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.set_aspect("equal")

	Z = np.linspace(0,1,layers)
	twist_per_loop = 2*np.pi*0.1/layers

	for Z in Z:
		ax.plot3D(X,Y,Z, 'b')
		A,Ix,Iy,Ixy,Cx,Cy = compute_poly_AI(X,Y)
		X,Y  = rotate_around_point(X,Y,twist_per_loop*Z,Cx,Cy)

	xx,yy,zz = np.meshgrid(X,Y,Z)

	# ax.set_xlim3d(-1, 1)
	# ax.set_ylim3d(-1, 1)
	# ax.set_zlim3d(0, 1)
	ax.axis('equal')
	plt.show()