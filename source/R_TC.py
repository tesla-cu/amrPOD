# ================================================= #
# Code:        Computing R in POD                   #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np

# =========================================================================== #
# Function to compute the covariance R in POD using a standard matrix 
# operation technique and the new algorithm leveraging AMR repetitions
#
# Inputs:
# - X      : snapshot matrix
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - R      : covariance matrix computed using matrix operations 
#              - if R == False, then we do not check the correctness
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
# - wt_art : weighting of arithmetic operations
# - wt_acc : weighting of accessing operations
# - wt_asn : weighting of assignment operations
# - wt_log : weighting of logical operations
# - wt_fun : weighting of function calls
#
# Outputs:
# - im : num of operations to compute R using implemented algorithm
# - un : num of operations to compute R using unaltered algorithm
# =========================================================================== #
def compute_R_TC(X, X_grid, R, d_l, nt, nspat, finest, \
	             wt_art, wt_acc, wt_asn, wt_log, wt_fun):

	# Initialize Operation Counts ---------------------------------------------

	# Unaltered Algorithm
	un     = 0 # Total operation counts
	un_art = 0 # Arithmetic operations
	un_acc = 0 # Access operations
	un_asn = 0 # Assignment operations
	un_fun = 0 # Function call
	
	# Implemented Algorithm
	im 	   = 0 # Total operation counts
	im_art = 0 # Arithmetic operations
	im_acc = 0 # Access operations
	im_asn = 0 # Assignment operations
	im_log = 0 # Logical operations
	im_fun = 0 # Function call

	# =========================================================================
	# Unaltered Computation 
	# =========================================================================

	R_un = np.empty((nt, nt))
	un_asn += wt_asn
	un_fun += wt_fun

	for m in range(nt):
		un_art += wt_art
		un_asn += wt_asn

		un_art += wt_art
		for n in range(m+1):
			un_art += wt_art
			un_asn += wt_asn

			r_sum = 0
			un_asn += wt_asn

			for i in range(nspat):
				un_art += wt_art
				un_asn += wt_asn

				r_sum += X[i,m] * X[i,n]
				un_acc += 2*wt_acc
				un_art += 2*wt_art
				un_asn += wt_asn

			R_un[m,n] = r_sum
			un_acc += wt_acc
			un_asn += wt_asn

			R_un[n,m] = r_sum
			un_acc += wt_acc
			un_asn += wt_asn

	# =========================================================================
	# Implemented Computation 
	# =========================================================================

	R_im = np.empty((nt, nt))
	im_asn += wt_asn
	im_fun += wt_fun

	d_f1 = d_l[finest-1]
	im_asn += wt_asn
	im_art += wt_art
	im_acc += wt_acc

	for m in range(nt):
		im_art += wt_art
		im_asn += wt_asn

		im_art += wt_art
		for n in range(m+1):
			im_art += wt_art
			im_asn += wt_asn

			r_sum  = 0
			im_asn += wt_asn

			i = 0
			im_asn += wt_asn

			im_log += wt_log
			if m == n:

				for ii in range(nspat): 
					im_art += wt_art
					im_asn += wt_asn

					im_log += wt_log
					if i < nspat:

						im_log += wt_log
						im_acc += wt_acc
						if X_grid[i,m] == finest:

							im_art += wt_art
							for j in range(i,i+d_f1):
								im_art += wt_art
								im_asn += wt_asn

								r_sum += X[j,n]*X[j,m]
								im_art += 2*wt_art
								im_acc += 2*wt_acc
								im_asn += wt_asn

							i += d_f1
							im_art += wt_art
							im_asn += wt_asn
						else:

							d_val  = d_l[X_grid[i,m]]
							im_acc += 2*wt_acc
							im_asn += wt_asn

							r_sum  += d_val*X[i,n]*X[i,m]
							im_acc += 2*wt_acc
							im_asn += wt_asn
							im_art += 3*wt_art

							i      += d_val
							im_asn += wt_asn
							im_art += wt_art

					else:
						break

			else:
				for ii in range(nspat):
					im_art += wt_art
					im_asn += wt_asn
					
					im_log += wt_log
					if i < nspat:

						im_log += 2*wt_log
						im_acc += 2*wt_acc
						if X_grid[i,m] == finest or X_grid[i,n] == finest:

							im_art += wt_art
							for j in range(i,i+d_f1):
								im_art += wt_art
								im_asn += wt_asn

								r_sum += X[j,n]*X[j,m]
								im_art += 2*wt_art
								im_acc += 2*wt_acc
								im_asn += wt_asn

							i += d_f1
							im_art += wt_art
							im_asn += wt_asn
						else:
						
							im_log += wt_log
							im_acc += 2*wt_acc
							if X_grid[i,m] > X_grid[i,n]:
								d_val  = d_l[X_grid[i,m]]
								im_acc += 2*wt_acc
								im_asn += wt_asn

							else:
								d_val  = d_l[X_grid[i,n]]
								im_acc += 2*wt_acc
								im_asn += wt_asn

							r_sum  += d_val*X[i,n]*X[i,m]
							im_acc += 2*wt_acc
							im_art += 3*wt_art
							im_asn += wt_asn

							i      += d_val
							im_art += wt_art
							im_asn += wt_asn

					else:
						break

			R_im[m,n] = r_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

			R_im[n,m] = r_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

	# =========================================================================
	# Check Correctness of Matrices
	# =========================================================================

	# Check if we should check for correctness
	if type(R) != bool:
	
		# Compute relative error for each cell
		err_im = np.max(abs(np.subtract(R_im, R)) / abs(R))
		err_un = np.max(abs(np.subtract(R_un, R)) / abs(R))

		if err_im < 1e-6:
			print('The implemented R is correct')
		else:
			print('The implemented R is incorrect')

		if err_un < 1e-6:
			print('The unaltered R is correct')
		else:
			print('The unaltered R is incorrect')

	# Sum operations from im and un 
	un = un_asn + un_acc + un_art + un_fun
	im = im_asn + im_acc + im_art + im_fun + im_log

	# Return op counts of implemented and unaltered algorithm
	return im, un
