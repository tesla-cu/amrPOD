# ================================================= #
# Code:        Computing A in POD                   #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
import time

# =========================================================================== #
# Function to compute the temporal coefficients A in POD using a  standard 
# matrix operation technique and the new algorithm leveraging AMR repetitions
#
# Inputs:
# - X      : snapshot matrix
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - Phi    : spatial mode matrix computed using matrix operations 
# - A      : temporal coefficent matrix computed using matrix 
#            operations 
#              - if A == False, then we do not check the correctness
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
# - finest : finest AMR grid level
#
# Outputs:
# - time_im : CPU time to compute A using implemented algorithm
# - time_un : CPU time to compute A using unaltered algorithm
# =========================================================================== #
def compute_A_CPU(X, X_grid, Phi, A, d_l, nt, nspat, finest):

	# =========================================================================
	# Unaltered Computation 
	# =========================================================================

	# Initialize timer
	tic = time.time()

	# Initialize A matrix for unaltered computation
	A_un = np.empty((nt, nt))

	# Compute A matrix with unaltered algorithm
	for m in range(nt):     # iterate over nt rows
		for n in range(nt): # iterate over nt columns

			# Initialize temporary variable to store sum of an 
			# element of A
			a_sum = 0

			# Compute inner product
			for i in range(nspat):
				a_sum += X[i,m] * Phi[i,n]

			# Assign value of element of A
			A_un[m,n] = a_sum

	# Compute total cpu time
	time_un = time.time() - tic

	# =========================================================================
	# Implemented Computation
	# =========================================================================

	# Initialize timer
	tic = time.time()

	# Initialize A matrix for computation of implemented algorithm
	A_im = np.empty((nt, nt))

	# Initialize matrix to store maximum grid level
	G = np.empty((nspat), dtype=int)

	d_f1 = d_l[finest-1]

	# Initialize index of spatial location
	i = 0

	# Find the finest cell for all spatial locations
	for ii in range(nspat): # dummy loop
		if i < nspat:       # exit loop if we are at the end

			# Initialize max grid level
			X_grid_max = X_grid[i,0]

			# Find the max grid level for a spatial location
			for m in range(1,nt):

				# Check if current cell is bigger than current max
				if X_grid[i,m] > X_grid_max:
					X_grid_max = X_grid[i,m]

					# If this is the finest, no point in continuing
					# looking for finer cells
					if X_grid_max == finest:
						break

			if X_grid_max == finest:
				G[i] = 1
				i    += d_f1
			else:
				G[i] =  d_l[X_grid_max] # get # of repeats
				i    += G[i]            # skip cells that are repeated
		else:
			break

	# Compute elements of A
	for m in range(nt):     # iterate over rows
		for n in range(nt): # iterate over columns

			# Initialize temporary variable to store sum of an 
			# element of A
			a_sum = 0

			# Initialize index of spatial location
			i = 0

			# Compute value of one element in A
			for ii in range(nspat): # dummy loop
				if i < nspat:       # exit loop if at end
					if G[i] == 1:
						for j in range(i,i+d_f1):
							a_sum += X[j,m]*Phi[j,n]
						i += d_f1
					else:
						a_sum += G[i]*X[i,m]*Phi[i,n] # weight computation
						i += G[i]                     # skip repeats
				else:
					break

			# Assign values of element of R
			A_im[m,n] = a_sum

	# Compute total cpu time
	time_im = time.time() - tic

	# =========================================================================
	# Check Correctness of Matrices
	# =========================================================================
	
	# Check if we should check for correctness
	if type(A) != bool:

		# Compute relative error for each cell
		err_im = np.max(abs(np.subtract(A_im, A)) / abs(A))
		err_un = np.max(abs(np.subtract(A_un, A)) / abs(A))

		if err_im < 1e-6:
			print('The implemented A is correct')
		else:
			print('The implemented A is incorrect')

		if err_un < 1e-6:
			print('The unaltered A is correct')
		else:
			print('The unaltered A is incorrect')

	# Return CPU time of implemented and unaltered algorithm
	return time_im, time_un
