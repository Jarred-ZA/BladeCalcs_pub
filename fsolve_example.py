import numpy as np
from scipy.optimize import fsolve


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

def Bnl(f,Guess_range):
	

	sol =  []
	Guess = []
	for x in xrange(1,Guess_range):
		sol.append(fsolve(f,x)[0])
		Guess.append(x)

	sol = remove_duplicates(sol)
	return sol


