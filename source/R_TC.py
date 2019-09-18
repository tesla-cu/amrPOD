# ================================================= #
# Code:        Computing R in POD                   #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np

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
def compute_R_TC(X, X_grid, R, d_l, nt, nspat, wt_art, wt_acc, wt_asn, wt_log):

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
	# Initialize R matrix for unaltered computation  
	R_un = np.zeros((nt, nt))
	un_asn += wt_asn

	# Compute R matrix with unaltered algorithm

	for m in range(nt):      # iterate over nt rows
		un_art += wt_art
		un_acc  += wt_acc

		for n in range(m+1): # iterate up to and including diagnoal
			un_art += wt_art
			un_acc  += wt_acc

			# Initialize temporary variable to store sum of an 
			# element of R
			r_sum    = 0
			un_asn += wt_asn

			# Compute inner product
			for i in range(nspat):
				un_art += wt_art
				un_acc  += wt_acc

				r_sum += X[i,m] * X[i,n]
				un_acc  += 2*wt_acc
				un_art  += 2*wt_art
				un_asn  += wt_asn

			# Using symmetry, assign values of element of R
			R_un[m,n] = r_sum
			un_acc += wt_acc
			un_asn += wt_asn

			R_un[n,m] = r_sum
			un_acc += wt_acc
			un_asn += wt_asn

	# ========== Implemented Computation ========================== #

	# Initialize R matrix for computation of implemented algorithm
	R_im = np.zeros((nt, nt))
	im_asn += wt_asn

	# Compute R matrix with implemented algorithm
	for m in range(nt):      # iterate over nt rows
		im_art += wt_art
		im_acc += wt_acc

		for n in range(m+1): # iterate up to and including diagnoal
			im_art += wt_art
			im_acc += wt_acc

			# Initialize temporary variable to store sum of an 
			# element of R
			r_sum  = 0
			im_asn += wt_asn

			# Initialize index of spatial location
			i = 0
			im_asn += wt_asn

			# If this is a diagonal element, we know the grid level
			# is the same between snapshots and we do not need to 
			# check which grid level is the higher because they are
			# the same
			im_log += wt_log
			if m == n:
				# Compute inner product by weighting and skipping 
				# cells that are repeated
				for ii in range(nspat): # dummy iteration
					im_art += wt_art
					im_acc += wt_acc

					im_log += wt_log
					if i < nspat:

						d_val  = d_l[X_grid[i,m]]     # # of repeats
						im_acc += 2*wt_acc
						im_asn += wt_asn

						r_sum  += d_val*X[i,n]*X[i,m] # weight computation
						im_acc += 2*wt_acc
						im_asn += wt_asn
						im_art += 3*wt_art

						i      += d_val                   # skip repeats
						im_asn += wt_asn
						im_art += wt_art

					else:
						break

			# Off-diagonal element
			else:
				# Compute inner product by weighting and skipping 
				# cells that are repeated

				for ii in range(nspat): # dummy iteration
					im_art += wt_art
					im_acc += wt_acc
					
					im_log += wt_log
					if i < nspat:
						
						# Determine which grid level is higher and 
						# use this as weight and # of repeats
						im_log += wt_log
						im_acc += wt_acc
						if X_grid[i,m] > X_grid[i,n]:

							d_val  = d_l[X_grid[i,m]]
							im_acc += 2*wt_acc
							im_asn += wt_asn

						else:
							d_val  = d_l[X_grid[i,n]]
							im_acc += 2*wt_acc
							im_asn += wt_asn

						r_sum  += d_val*X_grid[i,n]*X[i,m] # weight computation
						im_acc += 2*wt_acc
						im_art += 3*wt_art
						im_asn += wt_asn

						i      += d_val                        # skip repeats
						im_art += wt_art
						im_asn += wt_asn

					else:
						break

			# Using symmetry, assign values of element of R
			R_im[m,n] = r_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

			R_im[n,m] = r_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

	# ========== Check Correctness of Matrices ==================== #

	if np.max(abs(np.subtract(R_im, R))) < 1e-8:
		print('The implemented R is correct')
	else:
		print('The implemented R is incorrect')

	if np.max(abs(np.subtract(R_un, R))) < 1e-8:
		print('The unaltered R is correct')
	else:
		print('The unaltered R is incorrect')

	# ========== Sum operations from im and un =================== #

	un = un_asn + un_acc + un_art
	im = im_asn + im_acc + im_art + im_log

	# Return op counts of implemented and unaltered algorithm
	return im, un
