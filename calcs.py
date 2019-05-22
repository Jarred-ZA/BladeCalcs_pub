import math
import numpy as np 
from scipy.interpolate import splprep, splev
from scipy.optimize import fsolve

def lumped_mass_cantaliver_nat_freq(I,E,L,m):
	k = 3*E*I/(L**3)
	nat_freq = (math.sqrt(k/m))/(2*math.pi)
	return	nat_freq

def distributed_mass_cantaliver_nat_freq(I,E,L,A,rho,number):
	k = E*I/(rho*A) 
	nat_freq = [] 
	nat_freq.append(((0.59686*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((1.49418*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((2.50025*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	nat_freq.append(((3.49999*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
	return nat_freq

def compute_exact_AI(pts,shape):
	A = PolyArea(pts)
	Ixx, Iyy, Ixy = inertia(shape)
	return A,Ixx

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

def PolyArea(pts):
	x = pts[:,0]
	y = pts[:,1]
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def Euler_Bournoulli_Clamped_Free_frequancy(L,E,I,rho,A,number):
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

def Euler_Bournoulli_Clamped_Free_mode_shape_function(L,mode_number,x):
    ##work out weighted freq
    f = lambda x : np.cos(x)*np.cosh(x) + 1

    Bnl = Bnler(f, 20 )

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

def spline(pts):
	tck, u = splprep(pts.T, u=None, s=0.0, per=1) 
	u_new = np.linspace(u.min(), u.max(), 1000)
	X_spline, Y_spline = splev(u_new, tck, der=0)
	return X_spline,Y_spline

def remove_duplicates(lst):
	lst.sort()
	newlst = []
	newlst.append(lst[1])
	i = 0
	while i <  len(lst) - 1:
		if (lst[i] - newlst[-1]) > 0.1:
			newlst.append(lst[i])
		i += 1

	return newlst

def Bnler(f,Guess_range):
	

	sol =  []
	Guess = []
	for x in xrange(1,Guess_range):
		sol.append(fsolve(f,x)[0])
		Guess.append(x)

	sol = remove_duplicates(sol)
	return sol

def area(pts):
  # 'Area of cross-section.'
  
  if pts[0] != pts[-1]:
    pts = pts + pts[:1]
  x = [ c[0] for c in pts ]
  y = [ c[1] for c in pts ]
  s = 0
  for i in range(len(pts) - 1):
    s += x[i]*y[i+1] - x[i+1]*y[i]
  return s/2

def centroid(pts):
  # 'Location of centroid.'
  
  if pts[0] != pts[-1]:
    pts = pts + pts[:1]
  x = [ c[0] for c in pts ]
  y = [ c[1] for c in pts ]
  sx = sy = 0
  a = area(pts)
  for i in range(len(pts) - 1):
    sx += (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
    sy += (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
  return sx/(6*a), sy/(6*a)

def inertia(pts):
  # 'Moments and product of inertia about centroid.'
  
  if pts[0] != pts[-1]:
    pts = pts + pts[:1]
  x = [ c[0] for c in pts ]
  y = [ c[1] for c in pts ]
  sxx = syy = sxy = 0
  a = area(pts)
  cx, cy = centroid(pts)
  for i in range(len(pts) - 1):
    sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
    syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
    sxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i])*(x[i]*y[i+1] - x[i+1]*y[i])
  return sxx/12 - a*cy**2, syy/12 - a*cx**2, sxy/24 - a*cx*cy

def principal(Ixx, Iyy, Ixy):
  # 'Principal moments of inertia and orientation.'
  
  avg = (Ixx + Iyy)/2
  diff = (Ixx - Iyy)/2      # signed
  I1 = avg + sqrt(diff**2 + Ixy**2)
  I2 = avg - sqrt(diff**2 + Ixy**2)
  theta = atan2(-Ixy, diff)/2
  return I1, I2, theta
