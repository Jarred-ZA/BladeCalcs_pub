from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def plotter(X,Y,twist,layers):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_aspect("equal")

	Z = np.linspace(0,1,layers)
	twist_per_loop = 2*np.pi*0.1/layers

	# Local versions were previously defined here. If this function is to be used,
	# it needs to import compute_poly_AI and rotate_around_point from calcs.py
	# For now, we comment out the parts that depend on them as gui.py handles plotting.
	# for Z_val in Z:
	# 	ax.plot3D(X,Y,Z_val, 'b')
		# A,Ix,Iy,Ixy,Cx,Cy = compute_poly_AI(X,Y) # Needs import
		# X_rot, Y_rot  = rotate_around_point(X,Y,twist_per_loop*Z_val,Cx,Cy) # Needs import
		# X, Y = X_rot, Y_rot

	ax.axis('equal')
	plt.show()