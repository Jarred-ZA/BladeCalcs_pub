import math
import numpy as np 
from scipy.interpolate import splprep, splev
from scipy.optimize import fsolve
from math import sqrt, atan2
from scipy.linalg import eigh

def lumped_mass_cantaliver_nat_freq(I, E, L, m):
    """Calculate natural frequency of a cantilever beam with lumped mass.
    
    Args:
        I: Moment of inertia
        E: Young's modulus
        L: Length of beam
        m: Mass
        
    Returns:
        Natural frequency in Hz
    """
    k = 3*E*I/(L**3)
    nat_freq = (math.sqrt(k/m))/(2*math.pi)
    return nat_freq

def distributed_mass_cantaliver_nat_freq(I, E, L, A, rho, number):
    """Calculate natural frequencies of a cantilever beam with distributed mass.
    
    Args:
        I: Moment of inertia
        E: Young's modulus
        L: Length of beam
        A: Cross-sectional area
        rho: Density
        number: Number of modes to calculate
        
    Returns:
        List of natural frequencies in Hz
    """
    k = E*I/(rho*A) 
    nat_freq = [] 
    nat_freq.append(((0.59686*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
    nat_freq.append(((1.49418*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
    nat_freq.append(((2.50025*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
    nat_freq.append(((3.49999*math.pi)**2)/((L**2)*math.sqrt(k))*((1/(2*math.pi))))
    return nat_freq

def compute_exact_AI(pts, shape):
    """Calculate exact area and moment of inertia.
    
    Args:
        pts: Array of points
        shape: List of (x,y) coordinates
        
    Returns:
        Tuple of (area, moment of inertia about x-axis)
    """
    A = PolyArea(pts)
    Ixx, Iyy, Ixy = inertia(shape)
    return A, Ixx

def compute_poly_AI(X, Y):
    """Calculate area and moments of inertia using polygon method.
    
    Args:
        X: List of x-coordinates
        Y: List of y-coordinates
        
    Returns:
        Tuple of (area, Ix, Iy, Ixy, Cx, Cy)
    """
    # Area calculation
    A = 0
    for x in range(1, len(X)-2):
        A = A + (X[x]*Y[x+1] - X[x+1]*Y[x])
    A = 0.5*A

    # Centroid calculation
    Cx = 0
    Cy = 0
    for x in range(1, len(X)-2):
        Cx = Cx + (X[x] + X[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
        Cy = Cy + (Y[x] + Y[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
    Cx = Cx/(6*A) if A != 0 else 0 # Avoid division by zero
    Cy = Cy/(6*A) if A != 0 else 0 # Avoid division by zero
    
    # Moments of inertia calculation
    Ix = 0
    Iy = 0
    Ixy = 0
    for x in range(1, len(X)-2):
        Ix = Ix + (Y[x]**2 + Y[x]*Y[x+1] + Y[x+1]**2)*(X[x]*Y[x+1] - X[x+1]*Y[x])
        Iy = Iy + (X[x]**2 + X[x]*X[x+1] + X[x+1]**2)*(X[x]*Y[x+1] - X[x+1]*Y[x])
        Ixy = Ixy + (X[x]*Y[x+1] + 2*X[x]*Y[x] + 2*X[x+1]*Y[x+1] + X[x+1]*Y[x+1])*(X[x]*Y[x+1] - X[x+1]*Y[x])
    Ix = Ix/12
    Iy = Iy/12
    Ixy = Ixy/24
    return A, Ix, Iy, Ixy, Cx, Cy

def compute_approxomate_AI(X, Y):
    """Calculate approximate area and moment of inertia.
    
    Args:
        X: List of x-coordinates
        Y: List of y-coordinates
        
    Returns:
        Tuple of (approximate area, approximate moment of inertia)
    """
    c = (max(X) - min(X))
    t = (max(Y) - min(Y))
    h = (max(Y) + min(Y))/2

    tou = t/c
    eps = h/c

    # Assume proportional factors
    K_A = 0.6
    K_I = 0.036

    A_approx = K_A*(c**2)*tou
    I_approx = K_I*(c**4)*tou*((tou**2) + (eps**2))
    return A_approx, I_approx

def PolyArea(pts):
    """Calculate area of polygon using shoelace formula.
    
    Args:
        pts: Array of points
        
    Returns:
        Area of polygon
    """
    x = pts[:,0]
    y = pts[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def Euler_Bournoulli_Clamped_Free_frequancy(L, E, I, rho, A, number):
    """Calculate natural frequencies of a clamped-free Euler-Bernoulli beam.
    
    Args:
        L: Length of beam
        E: Young's modulus
        I: Moment of inertia
        rho: Density
        A: Cross-sectional area
        number: Number of modes to calculate
        
    Returns:
        List of natural frequencies in Hz
    """
    number = number + 1
    k = math.sqrt((E*I)/(rho*A))
    wn = []
    freq = []
    for x in range(1, number):
        wn.append(float(k*(((2*x -1)*math.pi/2)/L)**2))
    length_wn = len(wn)
    for x in range(length_wn):
        freq.append(float(wn[int(x)] * (1/(2*math.pi))))
    return freq

def Euler_Bournoulli_Clamped_Free_mode_shape_function(L, mode_number, x):
    """Calculate mode shape function for a clamped-free Euler-Bernoulli beam.
    
    Args:
        L: Length of beam
        mode_number: Mode number
        x: Position along beam
        
    Returns:
        Mode shape value at position x
    """
    # Work out weighted frequency
    f = lambda x: np.cos(x)*np.cosh(x) + 1
    Bnl = Bnler(f, 20)
    
    Bnl_mode = 0
    if mode_number <= 4:
        Bnl_mode = Bnl[mode_number]
    if mode_number > 4:
        Bnl_mode = float((((2*(mode_number+1) -1)*math.pi/2)))
    
    # Work out characteristic frequency
    Bn = float(Bnl_mode/L)
    
    # Work out sigma n
    sigma = 1
    if mode_number == 1:
        sigma = 0.7341
    if mode_number == 2:
        sigma = 1.0185
    if mode_number == 3:
        sigma = 0.9995
    
    x = np.cosh(Bn*x) - np.cos(Bn*x) - sigma*(np.sinh(Bn*x) - np.sin(Bn*x))
    return x

def spline(pts):
    """Create a spline interpolation of points.
    
    Args:
        pts: Array of points
        
    Returns:
        Tuple of (x-coordinates, y-coordinates) of spline
    """
    tck, u = splprep(pts.T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), 1000)
    X_spline, Y_spline = splev(u_new, tck, der=0)
    return X_spline, Y_spline

def remove_duplicates(lst):
    """Remove duplicate values from list.
    
    Args:
        lst: List of values
        
    Returns:
        List with duplicates removed
    """
    lst.sort()
    newlst = []
    newlst.append(lst[1])
    i = 0
    while i < len(lst) - 1:
        if (lst[i] - newlst[-1]) > 0.1:
            newlst.append(lst[i])
        i += 1
    return newlst

def Bnler(f, Guess_range):
    """Solve for Bn values.
    
    Args:
        f: Function to solve
        Guess_range: Range of guesses
        
    Returns:
        List of solutions
    """
    sol = []
    Guess = []
    for x in range(1, Guess_range):
        sol.append(fsolve(f, x)[0])
        Guess.append(x)
    
    sol = remove_duplicates(sol)
    return sol

def area(pts):
    """Calculate area of cross-section.
    
    Args:
        pts: List of (x,y) coordinates
        
    Returns:
        Area of cross-section
    """
    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    s = 0
    for i in range(len(pts) - 1):
        s += x[i]*y[i+1] - x[i+1]*y[i]
    return s/2

def centroid(pts):
    """Calculate location of centroid.
    
    Args:
        pts: List of (x,y) coordinates
        
    Returns:
        Tuple of (x-coordinate, y-coordinate) of centroid
    """
    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    sx = sy = 0
    a = area(pts)
    for i in range(len(pts) - 1):
        sx += (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sy += (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
    return sx/(6*a), sy/(6*a)

def inertia(pts):
    """Calculate moments and product of inertia about centroid.
    
    Args:
        pts: List of (x,y) coordinates
        
    Returns:
        Tuple of (Ixx, Iyy, Ixy)
    """
    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    sxx = syy = sxy = 0
    a = area(pts)
    cx, cy = centroid(pts)
    for i in range(len(pts) - 1):
        sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        sxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i])*(x[i]*y[i+1] - x[i+1]*y[i])
    return sxx/12 - a*cy**2, syy/12 - a*cx**2, sxy/24 - a*cx*cy

def principal(Ixx, Iyy, Ixy):
    """Calculate principal moments of inertia and orientation.
    
    Args:
        Ixx: Moment of inertia about x-axis
        Iyy: Moment of inertia about y-axis
        Ixy: Product of inertia
        
    Returns:
        Tuple of (I1, I2, theta)
    """
    avg = (Ixx + Iyy)/2
    diff = (Ixx - Iyy)/2      # signed
    I1 = avg + sqrt(diff**2 + Ixy**2)
    I2 = avg - sqrt(diff**2 + Ixy**2)
    theta = atan2(-Ixy, diff)/2
    return I1, I2, theta

def rotate_around_point(x, y, radians, x_pt, y_pt):
    """Rotate points (x, y) around a point (x_pt, y_pt) by radians."""
    adjusted_x = np.subtract(x, x_pt)
    adjusted_y = np.subtract(y, y_pt)

    cos_rad = np.cos(radians)
    sin_rad = np.sin(radians)

    qx = x_pt + cos_rad * adjusted_x - sin_rad * adjusted_y # Corrected rotation formula
    qy = y_pt + sin_rad * adjusted_x + cos_rad * adjusted_y # Corrected rotation formula
    return qx, qy

def run_calculations(X, Y, E, L, rho):
    # Get blade profile points
    # X, Y = get_blade_profile()
    
    # --- Calculate Properties ---
    A_approx, I_approx, mass_approx = calculate_approximate_properties(X, Y, L, rho)
    A_exact, I_exact, mass_exact = calculate_exact_properties(X, Y, L, rho)
    A_poly, Ix_poly, Iy_poly, Ixy_poly, mass_poly = calculate_polynomial_properties(X, Y, L, rho)
    
    # --- Calculate Frequencies ---
    lumped_approx = calculate_lumped_frequency(E, I_approx, L, mass_approx)
    dist_approx = calculate_distributed_frequency(E, I_approx, L, rho, A_approx)
    
    lumped_exact = calculate_lumped_frequency(E, I_exact, L, mass_exact)
    dist_exact = calculate_distributed_frequency(E, I_exact, L, rho, A_exact)
    
    lumped_poly = calculate_lumped_frequency(E, Ix_poly, L, mass_poly)
    dist_poly = calculate_distributed_frequency(E, Ix_poly, L, rho, A_poly)
    
    # --- Format Results String --- 
    results = """
========================================
        BLADE CALCULATION RESULTS
========================================

INPUT PARAMETERS:
-----------------
  Young's Modulus (E) : {E:.3e} N/m^2 
  Length (L)          : {L:.4f} m
  Density (rho)       : {rho:.1f} kg/m^3

AREA & INERTIA (Cross-Section):
---------------------------------
  Method      | Area (m^2)   | Ix (m^4)     | Iy (m^4)     | Ixy (m^4)    | Mass (kg)
  ------------|--------------|--------------|--------------|--------------|----------
  Approximate | {A_approx:12.6e} | {I_approx:12.6e} |    N/A       |    N/A       | {mass_approx:8.6f}
  Exact       | {A_exact:12.6e} | {I_exact:12.6e} |    N/A       |    N/A       | {mass_exact:8.6f}
  Polynomial  | {A_poly:12.6e} | {Ix_poly:12.6e} | {Iy_poly:12.6e} | {Ixy_poly:12.6e} | {mass_poly:8.6f}

NATURAL FREQUENCIES (Hz):
-------------------------
  Method      | Lumped Mass  | Distributed Mass (Modes 1-4)
  ------------|--------------|--------------------------------------------------
  Approximate | {lumped_approx:12.2f} | {dist_approx[0]:.2f}, {dist_approx[1]:.2f}, {dist_approx[2]:.2f}, {dist_approx[3]:.2f}
  Exact       | {lumped_exact:12.2f} | {dist_exact[0]:.2f}, {dist_exact[1]:.2f}, {dist_exact[2]:.2f}, {dist_exact[3]:.2f}
  Polynomial  | {lumped_poly:12.2f} | {dist_poly[0]:.2f}, {dist_poly[1]:.2f}, {dist_poly[2]:.2f}, {dist_poly[3]:.2f}

========================================
"""
    return results

def calculate_approximate_properties(X, Y, L, rho):
    # Simple rectangular approximation
    max_y = np.max(np.abs(Y))
    A = 2 * max_y * L
    I = (2 * max_y)**3 * L / 12
    mass = A * L * rho
    return A, I, mass

def calculate_exact_properties(X, Y, L, rho):
    # Use the actual profile shape
    A = np.trapz(Y, X) * L
    I = np.trapz(Y**3, X) * L / 12
    mass = A * L * rho
    return A, I, mass

def calculate_polynomial_properties(X, Y, L, rho):
    # Use polynomial approximation
    A, Ix, Iy, Ixy, Cx, Cy = compute_poly_AI(X, Y)
    mass_per_unit_length = A * rho
    total_mass = mass_per_unit_length * L
    
    # The moments of inertia calculated by compute_poly_AI are area moments (m^4).
    # They don't need to be multiplied by L again here.
    # Mass moments of inertia would be different.
    # Let's return area moments and total mass.
    return A, Ix, Iy, Ixy, total_mass

def calculate_lumped_frequency(E, I, L, m):
    # Natural frequency for lumped mass system
    k = 3 * E * I / L**3
    return np.sqrt(k / m) / (2 * np.pi)

def calculate_distributed_frequency(E, I, L, rho, A):
    # Natural frequencies for distributed mass system
    beta = np.array([1.875, 4.694, 7.855, 10.996])
    return beta**2 * np.sqrt(E * I / (rho * A * L**4)) / (2 * np.pi)
