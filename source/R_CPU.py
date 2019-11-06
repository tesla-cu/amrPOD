# ================================================= #
# Code:        Computing R in POD                   #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
import time

# ================================================================= #
# Function to compute the covariance R in POD using a standard 
# matrix operation technique and the new algorithm leveraging AMR
# repetitions
#
# Inputs:
# - X      : snapshot matrix
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - R      : covariance matrix computed using matrix operations 
#            (this is used as a check we did the computation right)
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
#
# Outputs:
# - time_im : CPU time to compute R using implemented algorithm
# - time_un : CPU time to compute R using unaltered algorithm
# ================================================================= #
def compute_R_CPU(X, X_grid, R, d_l, nt, nspat):

	# ========== Unaltered Computation ============================ #

	# Initialize timer
	tic  = time.time()

	# Initialize R matrix for unaltered computation  
	R_un = np.zeros((nt, nt))

	# Compute R matrix with unaltered algorithm
	for m in range(nt):      # iterate over nt rows
		for n in range(m+1): # iterate up to and including diagnoal

			# Initialize temporary variable to store sum of an 
			# element of R
			r_sum = 0

			# Compute inner product
			for i in range(nspat):
				r_sum += X[i,m] * X[i,n]

			# Using symmetry, assign values of element of R
			R_un[m,n] = r_sum
			R_un[n,m] = r_sum

	# Compute total cpu time
	time_un = time.time() - tic

	# ========== Implemented Computation ========================== #

	# Initialize timer
	tic   = time.time()

	# Initialize R matrix for computation of implemented algorithm
	R_im = np.zeros((nt, nt))

	# Compute R matrix with implemented algorithm
	for m in range(nt):      # iterate over nt rows
		for n in range(m+1): # iterate up to and including diagnoal

			# Initialize temporary variable to store sum of an 
			# element of R
			r_sum = 0

			# Initialize index of spatial location
			i = 0

			# If this is a diagonal element, we know the grid level
			# is the same between snapshots and we do not need to 
			# check which grid level is the higher because they are
			# the same
			if m == n:
				# Compute inner product by weighting and skipping 
				# cells that are repeated
				for ii in range(nspat): # dummy iteration
					if i < nspat: 		
						d_val = d_l[X_grid[i,m]]     # # of repeats
						r_sum += d_val*X[i,n]*X[i,m] # weight computation
						i += d_val                   # skip repeats
					else:
						break

			# Off-diagonal element
			else:
				# Compute inner product by weighting and skipping 
				# cells that are repeated
				for ii in range(nspat): # dummy iteration
					if i < nspat:
						# Determine which grid level is higher and 
						# use this as weight and # of repeats
						if X_grid[i,m] > X_grid[i,n]:
							d_val = d_l[X_grid[i,m]]
						else:
							d_val = d_l[X_grid[i,n]]
						r_sum += d_val*X[i,n]*X[i,m] # weight computation
						i += d_val                   # skip repeats
					else:
						break

			# Using symmetry, assign values of element of R
			R_im[m,n] = r_sum
			R_im[n,m] = r_sum

	# Compute total cpu time
	time_im = time.time() - tic

	# ========== Check Correctness of Matrices ==================== #
	
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

	# Return CPU time of implemented and unaltered algorithm
	return time_im, time_un
