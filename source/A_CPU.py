import numpy as np
import time

def compute_A_CPU(X, X_grid, Phi, A, d_l, nt, nspat, finest):

#--------- Unaltered 

	tic     = time.time()
	A_unalt = np.zeros((nt, nt))
	for i in range(nt):
		for j in range(nt):
			a_sum = 0
			for k in range(nspat):
				a_sum += X[k,i] * Phi[k,j]
			A_unalt[i,j] = a_sum
	time_A_unalt = time.time() - tic

#--------- Implemented

	tic   = time.time()
	A_imp = np.zeros((nt,nt))
	G     = np.zeros(nspat, dtype=int)
	i     = 0
	for n in range(nspat):
		if i < nspat:
			X_grid_max = X_grid[i,0]
			for j in range(1,nt):
				if X_grid[i,j] > X_grid_max:
					X_grid_max = X_grid[i,j]
			G[i]  = d_l[X_grid_max]
			i    += G[i]
		else:
			break

	for i in range(nt):
		for j in range(nt):
			a_sum = 0
			k     = 0
			for n in range(nspat):
				if k < nspat:
					a_sum += G[k] * X[k,i] * Phi[k,j]
					k += G[k]
				else:
					break
			A_imp[i,j] = a_sum
	time_A_imp = time.time() - tic

#--------- Theoretical Computation

	theory_mult = 0
	theory_add  = 0

	for i in range(nt):
		for j in range(nt):
			k = 0
			for n in range(nt):
				if k < nspat:
					theory_mult += 2
					theory_add  += 1
					k += G[k]
				else:
					break

	unalt_mult = nspat * nt**2
	unalt_add  = nspat * (nt-1) * nt

#--------- Time for individual operations

	B = 0
	tic = time.time()
	for i in range(unalt_add):
		B = B + 1
	time_add_unalt = time.time() - tic

	B = 1.0000001
	tic = time.time()
	for i in range(unalt_mult):
		B = B*1.0000001
	time_mult_unalt = time.time() - tic

	time_A_theory = (theory_mult/unalt_mult)*time_mult_unalt + (theory_add/unalt_add)*time_add_unalt

#--------- Check matrices for correctness

	if np.max(abs(np.subtract(A_imp, A))) < 1e-8:
		print('The implemented A is correct')
	else:
		print('The implemented A is incorrect')

	if np.max(abs(np.subtract(A_unalt, A))) < 1e-8:
		print('The unaltered A is correct')
	else:
		print('The unaltered A is incorrect')

	return time_A_unalt, time_A_imp



