import math
import numpy as np 
import naca
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import Tkinter
import matplotlib.animation as animation
import calcs
import fem
import plot3D


# def get_profile():

# 	x = []
# 	y = []
# 	all_lines = []
# 	with open('profile.txt') as f:
# 		reader = csv.reader(f)
# 		for count ,(c1,c2) in enumerate(reader):
# 			x.append(float(c1))
# 			y.append(float(c2))
# 	return x,y


##unit profile
#X,Y = get_profile()
profile = '2612'
X,Y = naca.naca(profile,600)


##Tungsten Alloy (90W,7Ni,3Fe)
E = 650*10**6		#N/m^2
L = 0.05			#m
rho = 1710			#kg/m^3

##scale profile
X = [x*L for x in X]
Y = [x*L for x in Y]
pts = np.column_stack((np.array(list(X)),np.array(list(Y))))
shape = list(zip(X,Y))

print'Material Choice Tungsten Alloy (90W,7Ni,3Fe)'
print'Properties: E = {} N/m^3 \t L = {} m \t rho = {} kg/m^3'.format(E,L,rho)

###compute Profile Geometric Properties
A_approx,Ix_approx = calcs.compute_approxomate_AI(X,Y)
A_exact,I_exact = calcs.compute_exact_AI(pts,shape)
A_poly,Ix_poly,Iy_poly,Ixy_poly,Cx_poly,Cy_poly = calcs.compute_poly_AI(X,Y)

m_approx = rho*L*A_approx 
m_exact = rho*L*A_exact
m_poly = rho*L*A_poly

print '\nArea_approx = {}\tI_approx = {}\tmass_approx = {}'.format(A_approx,Ix_approx,m_approx)
print 'Area_exact = {}\tI_exact = {}\t\t\t\t\tmass_exact = {}'.format(A_exact,I_exact,m_exact)
print 'Area_tri = {}\tIxy_poly = {}\tmass_poly = {}\nIx_poly = {}\t\tIy_poly = {}'.format(A_poly,Ixy_poly,m_poly,Ix_poly,Iy_poly)


#print '\n####Anylitical Solution approx####'
number = 4
## Assume I_exact = Ix_poly
I_exact = Ix_poly
lumped_nat_freq_approx = calcs.lumped_mass_cantaliver_nat_freq(Ix_approx,E,L,m_approx)
lumped_nat_freq_exact = calcs.lumped_mass_cantaliver_nat_freq(I_exact,E,L,m_exact)
lumped_nat_freq_poly = calcs.lumped_mass_cantaliver_nat_freq(Ix_poly,E,L,m_poly)
distributed_nat_freq_approx = calcs.distributed_mass_cantaliver_nat_freq(Ix_approx,E,L,A_approx,rho,number)
distributed_nat_freq_exact = calcs.distributed_mass_cantaliver_nat_freq(I_exact,E,L,A_exact,rho,number)
distributed_nat_freq_poly = calcs.distributed_mass_cantaliver_nat_freq(Ix_poly,E,L,A_poly,rho,number)

print'\n'
print 'lumped_approx {} Hz\t\tdist_approx {} Hz'.format(lumped_nat_freq_approx, distributed_nat_freq_approx)
print 'Lumped_exact {} Hz\t\tdist_exact {} Hz'.format(lumped_nat_freq_exact, distributed_nat_freq_exact)
print 'lumped_poly {} Hz\t\tdist_poly {} Hz'.format(lumped_nat_freq_poly, distributed_nat_freq_poly)

print '\nEuler-Bournoulli Clamped-Free\t\t{} Hz'.format(calcs.Euler_Bournoulli_Clamped_Free_frequancy(L,E,Ix_poly,rho,A_poly,number))

exact_frequency = math.pi/2
results = []
for i in range(1,10):
	M, K, frequencies, evecs = fem.bar(i,E,Ix_poly,L,m_poly,rho,A_poly)
print 'Fund. Frequency: {}'.format( round(frequencies[0],3))

#print 'Exact frequency: ', round(exact_frequency,3)

# plot the results
elements = np.array([x[0] for x in results])
errors   = np.array([x[1] for x in results])
plt.figure(1)
plt.plot(elements,errors, 'o-')
#plt.xlim(elements[0], elements[-1])
plt.xlabel('Number of Elements')
plt.ylabel('Error (%)')
plt.close()
plt.show()

points = 100
x = np.linspace(0,L,points)

for i in xrange(0,5):

    y = np.arange(points)

    ##get Y values
    y = calcs.Euler_Bournoulli_Clamped_Free_mode_shape_function(L,i,x)
    plt.plot(x,y,'-')


plt.grid(True,which = 'both')
plt.show()

 
#plot profile	
X_spline,Y_spline = calcs.spline(pts)
plt.plot(X_spline, Y_spline , 'b--')
plt.axis('equal')
plt.plot(X,Y)
plt.plot(Cx_poly,Cy_poly,'ro')
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.legend('profile')
plt.text(Cx_poly,Cy_poly,'({},{})'.format(round(Cx_poly,2),round(Cy_poly,2)))
plt.show()
##plot 3D 
plot3D.plotter(X,Y,0.1,25)




