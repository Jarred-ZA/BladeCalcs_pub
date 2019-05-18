import math
import numpy as np 
import naca
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.linalg import eigh
import fsolve_example
import Tkinter
import matplotlib.animation as animation



def lumped_mass_cantaliver_nat_freq(I,E,L,m):
	k = 3*E*I/(L**3)
	nat_freq = (math.sqrt(k/m))/(2*math.pi)
	return	nat_freq

def distributed_mass_cantaliver_nat_freq(I,E,L,A,rho,number):
	k = E*I/(rho*A) 
	nat_freq = [] 
	#for i in xrange(1,number):
		#nat_freq.append(((		(2*number -1)		*math.pi/2)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((0.59686*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((1.49418*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((2.50025*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((3.49999*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	return nat_freq

def get_profile():

	x = []
	y = []
	all_lines = []
	with open('profile.txt') as f:
		reader = csv.reader(f)
		for count ,(c1,c2) in enumerate(reader):
			x.append(float(c1))
			y.append(float(c2))
	return x,y

class Display(object):
    def __init__(self):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.h = []
        self.label = []
        self.fig, self.ax = self.plt.subplots()
        self.plt.axis('equal')
        self.plt.xlabel('x')
        self.plt.ylabel('y')
        self.ax.grid(True)
    def plot(self, X, Y, Cx,Cy,label1,label2):
        h, = self.plt.plot(X, Y, '-', linewidth = 1)
        self.h.append(h)
     	h, = self.plt.plot(Cx, Cy, 'ro', linewidth = 1)
        self.h.append(h)
        self.label.append(label1)
        self.label.append(label2)

    def show(self):
        self.plt.axis((0,max(X))+self.plt.axis()[2:])
        self.ax.legend(self.h, self.label)
        self.plt.show()

def compute_approxomate_AI(X,Y):
	## compute aproxomate area and bending inertia
	c = (max(X) - min(X))
	t = (max(Y)- min(Y))
	h = (max(Y)+ min(Y))/2

	tou = t/c
	eps = h/c

	##Assume Proportial factors
	K_A = 0.6
	K_I = 0.036

	A_approx = K_A*(c**2)*tou
	I_approx = K_I*(c**4)*tou*((tou**2) + (eps**2))
	return A_approx,I_approx

def PolyArea(x,y):
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def compute_exact_AI(X,Y):
	A = PolyArea(X,Y)
	I = 100
	return A,I

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

def bar(num_elems,E,I,L,m,rho,A):
	restrained_dofs = [0]

	k_element = E*I/(rho*A)
	m_element = m/num_elems

	# element mass and stiffness matrices for a bar
	m = np.array([[2,1],[1,2]]) / (6. * num_elems)
	k = np.array([[1,-1],[-1,1]]) * float(num_elems)

	m = m*m_element
	k = k*k_element

	# construct global mass and stiffness matrices
	M = np.zeros((num_elems+1,num_elems+1))
	K = np.zeros((num_elems+1,num_elems+1))

	# assembly of elements
	for i in range(num_elems):
		M_temp = np.zeros((num_elems+1,num_elems+1))
		K_temp = np.zeros((num_elems+1,num_elems+1))
		M_temp[i:i+2,i:i+2] = m
		K_temp[i:i+2,i:i+2] = k
		M += M_temp
		K += K_temp

	# remove the fixed degrees of freedom
	for dof in restrained_dofs:
		for i in [0,1]:
			M = np.delete(M, dof, axis=i)
			K = np.delete(K, dof, axis=i)

	# eigenvalue problem
	evals, evecs = eigh(K,M)
	frequencies = np.sqrt(evals)
	return M, K, frequencies, evecs

def Euler_Bournoulli_Clamped_Free_frequancy(E,I,rho,A,number):
	number = number +1
	k = math.sqrt((E*I)/(rho*A) )
	wn = []
	freq = []
	for x in xrange(1,number):
		wn.append(float(k*(((2*x -1)*math.pi/2)/L)**2 ))
	length_wn = len(wn)
	for x in range(length_wn):
	     freq.append(float(wn[int(x)] * (1/(2*math.pi))  ))
 	return (freq)

def Euler_Bournoulli_Clamped_Free_mode_shape_function(mode_number,x):
    ##work out weighted freq
    f = lambda x : np.cos(x)*np.cosh(x) + 1

    Bnl = fsolve_example.Bnl(f, 20 )

    Bnl_mode = 0
    if mode_number <= 4:
    	Bnl_mode = Bnl[mode_number]
    if mode_number > 4:	
    	Bnl_mode = float((((2*(mode_number+1) -1)*math.pi/2)))

    #print 'Bnl {}'.format(Bnl_mode)
    #work out characteristic freq
    Bn = float(Bnl_mode/L)
    #print 'Bn {}'.format(Bn)
    #work out sigma n
    sigma = 1
    if mode_number == 1:
    	sigma = 0.7341


	if mode_number == 2:
		sigma = 1.0185


    if mode_number == 3:
		sigma = 0.9995

    x = np.cosh(Bn*x) - np.cos(Bn*x) - sigma* ( np.sinh(Bn*x)  - np.sin(Bn*x) )
    return x

##unit profile
#X,Y = get_profile()
profile = '2612'
X,Y = naca.naca(profile,600)

##Tungsten Alloy (90W,7Ni,3Fe)
E = 650*10**6		#N/m^2
L = 0.5			#m
rho = 1710			#kg/m^3

##scale profile
X = [x*L for x in X]
Y = [x*L for x in Y]

print'Material Choice Tungsten Alloy (90W,7Ni,3Fe)'
print'Properties: E = {} N/m^3 \t L = {} m \t rho = {} kg/m^3'.format(E,L,rho)

###compute Profile Geometric Properties
A_approx,Ix_approx = compute_approxomate_AI(X,Y)
A_exact,I_exact = compute_exact_AI(X,Y)
A_poly,Ix_poly,Iy_poly,Ixy_poly,Cx_poly,Cy_poly = compute_poly_AI(X,Y)

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
lumped_nat_freq_approx = lumped_mass_cantaliver_nat_freq(Ix_approx,E,L,m_approx)
lumped_nat_freq_exact = lumped_mass_cantaliver_nat_freq(I_exact,E,L,m_exact)
lumped_nat_freq_poly = lumped_mass_cantaliver_nat_freq(Ix_poly,E,L,m_poly)
distributed_nat_freq_approx = distributed_mass_cantaliver_nat_freq(Ix_approx,E,L,A_approx,rho,number)
distributed_nat_freq_exact = distributed_mass_cantaliver_nat_freq(I_exact,E,L,A_exact,rho,number)
distributed_nat_freq_poly = distributed_mass_cantaliver_nat_freq(Ix_poly,E,L,A_poly,rho,number)

print'\n'
print 'lumped_approx {} Hz\t\tdist_approx {} Hz'.format(lumped_nat_freq_approx, distributed_nat_freq_approx)
print 'Lumped_exact {} Hz\t\tdist_exact {} Hz'.format(lumped_nat_freq_exact, distributed_nat_freq_exact)
print 'lumped_poly {} Hz\t\tdist_poly {} Hz'.format(lumped_nat_freq_poly, distributed_nat_freq_poly)

print '\nEuler-Bournoulli Clamped-Free\t\t{} Hz'.format(Euler_Bournoulli_Clamped_Free_frequancy(E,Ix_poly,rho,A_poly,number))

exact_frequency = math.pi/2
results = []
for i in range(1,10):
	M, K, frequencies, evecs = bar(i,E,Ix_poly,L,m_poly,rho,A_poly)
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
    y = Euler_Bournoulli_Clamped_Free_mode_shape_function(i,x)
    plt.plot(x,y,'-')

# top = Tkinter.Tk()
# # Code to add widgets will go here...
# top.mainloop()


plt.grid(True,which = 'both')
plt.show()

 
#plot profile	
#d = Display()
plt.figure(2)
plt.plot(X,Y)
plt.plot(Cx_poly,Cy_poly,'ro')
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.legend('profile')
plt.text(Cx_poly,Cy_poly,'({},{})'.format(round(Cx_poly,2),round(Cy_poly,2)))
plt.show()




