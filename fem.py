import numpy as np
from scipy.linalg import eigh


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

def beam(num_elems,E,I,L,m,rho,A):
	restrained_dofs = [1, 0, -2, -1]

	l = L
	Cm = rho*A
	Ck = E*I

	# element mass and stiffness matrices
	m = np.array([[156, 22*l, 54, -13*l],
				  [22*l, 4*l*l, 13*l, -3*l*l],
				  [54, 13*l, 156, -22*l],
				  [-13*l, -3*l*l, -22*l, 4*l*l]]) * Cm * l / 420

	k = np.array([[12, 6*l, -12, 6*l],
				  [6*l, 4*l*l, -6*l, 2*l*l],
				  [-12, -6*l, 12, -6*l],
				  [6*l, 2*l*l, -6*l, 4*l*l]]) * Ck / l**3

	# construct global mass and stiffness matrices
	M = np.zeros((2*num_elems+2,2*num_elems+2))
	K = np.zeros((2*num_elems+2,2*num_elems+2))

	# for each element, change to global coordinates
	for i in range(num_elems):
		M_temp = np.zeros((2*num_elems+2,2*num_elems+2))
		K_temp = np.zeros((2*num_elems+2,2*num_elems+2))
		M_temp[2*i:2*i+4, 2*i:2*i+4] = m
		K_temp[2*i:2*i+4, 2*i:2*i+4] = k
		M += M_temp
		K += K_temp

	# remove the fixed degrees of freedom
	for dof in restrained_dofs:
		for i in [0,1]:
			M = np.delete(M, dof, axis=i)
			K = np.delete(K, dof, axis=i)

	evals, evecs = eigh(K,M)
	frequencies = np.sqrt(evals)
	return M, K, frequencies, evecs