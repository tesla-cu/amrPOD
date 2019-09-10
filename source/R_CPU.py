import numpy as np
import time
#from GridGen import GridGen

def compute_R_CPU(X, X_grid, R, d_l, nt, nspat):

#--------- Unaltered 

	tic     = time.time()
	R_unalt = np.zeros((nt, nt))
	for i in range(nt):
		for j in range(i+1):
			r_sum = 0
			for k in range(nspat):
				r_sum += X[k,i] * X[k,j]
			R_unalt[i,j] = r_sum
			R_unalt[j,i] = r_sum
	time_unalt = time.time() - tic

#--------- Implemented Computation

	tic   = time.time()
	R_imp = np.zeros((nt, nt))
	for i in range(nt):
		for j in range(i+1):
			r_sum = 0
			k = 0
			if i == j:
				for n in range(nspat):
					if k < nspat: 		
						d_val = d_l[X_grid[k,i]]
						r_sum += d_val * X[k,j] * X[k,i]
						k += d_val
					else:
						break
			else:
				for n in range(nspat):
					if k < nspat:
						if X_grid[k,i] > X_grid[k,j]:
							d_val = d_l[X_grid[k,i]]
						else:
							d_val = d_l[X_grid[k,j]]
						r_sum += d_val * X_grid[k, j] * X[k,i]
						k += d_val
					else:
						break
			R_imp[i,j] = r_sum
			R_imp[j,i] = r_sum
	time_imp = time.time() - tic

	#--------- Theoretical Computation

	# theory_mult = 0
	# theory_add  = 0
	# for i in range(nt):
	# 	for j in range(i+1):
	# 		if i == j:
	# 			r_sum = 0
	# 			k = 0
	# 			for n in range(nspat):
	# 				if k < nspat: 		
	# 					d_val = d_l[X_grid[k,i]]
	# 					theory_mult += 2
	# 					theory_add  += 1
	# 					k += d_val
	# 				else:
	# 					break
	# 		else:
	# 			r_sum = 0
	# 			k = 0
	# 			for n in range(nspat):
	# 				if k < nspat:
	# 					if X_grid[k,i] > X_grid[k,j]:
	# 						d_val = d_l[X_grid[k,i]]
	# 					else:
	# 						d_val = d_l[X_grid[k,j]]
	# 					theory_mult += 2
	# 					theory_add  += 1
	# 					k += d_val
	# 				else:
	# 					break

	# unalt_mult = nspat*nt*(nt+1)
	# unalt_add  = (nspat-1)*nt*(nt+1)

#--------- Time for individual operations

	# B = 0
	# tic = time.time()
	# for i in range(unalt_add):
	# 	B = B + 1
	# time_add_unalt = time.time() - tic

	# B = 1.0000001
	# tic = time.time()
	# for i in range(unalt_mult):
	# 	B = B*1.0000001
	# time_mult_unalt = time.time() - tic

	# time_theory = (theory_mult/unalt_mult)*time_mult_unalt + (theory_add/unalt_add)*time_add_unalt

#--------- Check matrices for correctness

	if np.max(abs(np.subtract(R_imp, R))) < 1e-8:
		print('The implemented R is correct')
	else:
		print('The implemented R is incorrect')

	if np.max(abs(np.subtract(R_unalt, R))) < 1e-8:
		print('The unaltered R is correct')
	else:
		print('The unaltered R is incorrect')

	# return time_imp, time_unalt, time_theory
	return time_imp, time_unalt

	

