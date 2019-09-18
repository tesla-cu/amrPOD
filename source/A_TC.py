# ================================================= #
# Code:        Computing A in POD                   #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np

# ================================================================= #
# Function to compute the temporal coefficients A in POD using a 
# standard matrix operation technique and the new algorithm 
# leveraging AMR repetitions
#
# Inputs:
# - X      : snapshot matrix
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - Phi    : spatial mode matrix computed using matrix operations 
# - A      : temporal coefficent matrix computed using matrix 
#            operations (this is used as a check we did the 
#            computation right)
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
# - finest : finest AMR grid level
#
# Outputs:
# - time_im : CPU time to compute A using implemented algorithm
# - time_un : CPU time to compute A using unaltered algorithm
# ================================================================= #
def compute_A_TC(X, X_grid, Phi, A, d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log):


	# ========== Define Unaltered Op Counts ======================= #
	un        = 0          # Total operation counts
	un_art    = 0			# Arithmetic operations
	un_acc    = 0			# Access operations
	un_asn    = 0			# Assignment operations
	
	# ========== Define Implemented Op Counts ===================== #
	im 	      = 0			# Total operation counts
	im_art    = 0			# Arithmetic operations
	im_acc    = 0			# Access operations
	im_asn    = 0			# Assignment operations
	im_log    = 0			# Logical operations


	# ========== Unaltered Computation ============================ #

	# Initialize A matrix for unaltered computation
	A_un = np.zeros((nt, nt))
	un_asn += wt_asn

	# Compute A matrix with unaltered algorithm
	for m in range(nt):     # iterate over nt rows
		un_art += wt_art
		un_acc += wt_acc

		for n in range(nt): # iterate over nt columns
			un_art += wt_art
			un_acc += wt_acc

			# Initialize temporary variable to store sum of an 
			# element of A
			a_sum  = 0
			un_asn += wt_asn

			# Compute inner product
			for i in range(nspat):
				un_art += wt_art
				un_acc += wt_acc

				a_sum  += X[i,m] * Phi[i,n]
				un_acc += 2*wt_acc
				un_art += 2*wt_art
				un_asn += wt_asn

			# Assign value of element of A
			A_un[m,n] = a_sum
			un_acc    += wt_acc
			un_asn    += wt_asn

	# ========== Implemented Computation ========================== #

	# Initialize R matrix for computation of implemented algorithm
	A_im   = np.zeros((nt, nt))
	im_asn += wt_asn

	# Initialize matrix to store maximum grid level
	G      = np.zeros((nspat), dtype=int)
	im_asn += wt_asn	

	# Initialize index of spatial location
	i      = 0
	im_asn += wt_asn

	# Find the finest cell for all spatial locations
	for ii in range(nspat): # dummy loop
		im_art += wt_art
		im_acc += wt_acc

		im_log += wt_log
		if i < nspat:       # exit loop if we are at the end

			# Initialize max grid level
			X_grid_max = X_grid[i,0]
			im_acc     += wt_acc
			im_asn     += wt_asn

			# Find the max grid level for a spatial location
			for m in range(1,nt):
				im_art += wt_art
				im_acc += wt_acc

				# Check if current cell is bigger than current max
				im_log += wt_log
				im_acc += wt_acc
				if X_grid[i,m] > X_grid_max:

					X_grid_max = X_grid[i,m]
					im_acc     += wt_acc
					im_asn     += wt_asn

					# If this is the finest, no point in continuing
					# looking for finer cells
					im_log += wt_log
					if X_grid_max == finest:
							break

			G[i]   =  d_l[X_grid_max] # get # of repeats
			im_acc += 2*wt_acc
			im_asn += wt_asn

			i      += G[i]            # skip cells that are repeated
			im_art += wt_art
			im_acc += wt_acc
			im_asn += wt_asn

		else:
			break

	# Compute elements of A
	for m in range(nt):     # iterate over rows
		im_art += wt_art
		im_acc += wt_acc

		for n in range(nt): # iterate over columns
			im_art += wt_art
			im_acc += wt_acc

			# Initialize temporary variable to store sum of an 
			# element of A
			a_sum  = 0
			im_asn += wt_asn

			# Initialize index of spatial location
			i      = 0
			im_asn += wt_asn

			# Compute value of one element in A
			for ii in range(nspat): # dummy loop
				im_art += wt_art
				im_acc += wt_acc

				im_log += wt_log
				if i < nspat:       # exit loop if at end

					a_sum  += G[i]*X[i,m]*Phi[i,n] # weight computation
					im_acc += 3*wt_acc
					im_art += 3*wt_art
					im_asn += wt_asn

					i      += G[i]                     # skip repeats
					im_acc += wt_acc
					im_asn += wt_asn
					im_art += wt_art

				else:
					break

			# Assign values of element of R
			A_im[m,n] = a_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

	# ========== Check Correctness of Matrices ==================== #

	if np.max(abs(np.subtract(A_im, A))) < 1e-8:
		print('The implemented A is correct')
	else:
		print('The implemented A is incorrect')

	if np.max(abs(np.subtract(A_un, A))) < 1e-8:
		print('The unaltered A is correct')
	else:
		print('The unaltered A is incorrect')

	# ========== Sum operations from im and un =================== #

	un = un_asn + un_acc + un_art
	im = im_asn + im_acc + im_art + im_log

	# Return op counts of implemented and unaltered algorithm
	return im, un

	