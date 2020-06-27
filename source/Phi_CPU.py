# ================================================= #
# Code:        Computing Phi in POD                 #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
import time

# =========================================================================== #
# Function to compute the POD spatial modes in using a standard  matrix 
# operation technique and the new algorithm leveraging AMR repetitions
#
# Inputs:
# - X      : snapshot matrix
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - Psi    : matrix containing eigenvectors of R
# - Lambda : matrix containing eigenvalues of R
# - method : integer (1 or 2) specifying the method to use for the 
#            new algorithm
# - Phi    : spatial mode matrix computed using matrix operations 
#              - if Phi == False, then we do not check the correctness
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
# - finest : finest AMR grid level
#
# Outputs:
# - time_im : CPU time to compute Phi using implemented algorithm
# - time_un : CPU time to compute Phi using unaltered algorithm
# =========================================================================== #
def compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest):

	# =========================================================================
	# Unaltered Computation 
	# =========================================================================

	# Initialize timer
	tic = time.time()

	# Initialize Phi matrix for unaltered computation
	Phi_un = np.empty((nspat,nt))

	# Compute Phi matrix with unaltered algorithm
	for i in range(nspat):  # iterate through rows of Phi
		for m in range(nt):  # iterate through columns of Phi

			# Initialize temporary variable to store sum of an 
			# element of Phi
			phi_sum = 0

			# Compute inner product
			for n in range(nt):
				phi_sum += X[i,n] * Psi[n,m]

			# Assign value of element of Phi
			Phi_un[i,m] = phi_sum/np.sqrt(Lambda[m,m])

	# Compute total cpu time
	time_un = time.time() - tic
	
	# =========================================================================
	# Implemented Computation - Method 1
	# =========================================================================

	if method == 1:

		# Initialize timer
		tic = time.time()

		# Initialize Phi matrix for implemented computation
		Phi_im = np.empty((nspat, nt))

		# Define number of repetitions for a coarse cell
		d_0 = d_l[0]

		# Iterate on the cells that could be new coarse cells
		for i in range(0, nspat, d_0):

			# Initialize matrix to store the maximum grid level at 
			# a particular cell location
			G = np.zeros((d_0), dtype=int)

			idx = 0 # index used for iterating within a coarse cell
			j = i   # index of global matrix

			# Find maximum value for each spatial location
			for jj in range(d_0): # dummy loop
				if idx < d_0:     # iterate within coarse cell

					# Initialize max grid level
					X_grid_max = X_grid[j,0]

					# Find the max grid level for a spatial location
					for m in range(1,nt):

						# Check if current cell is bigger than current max
						if X_grid[j,m] > X_grid_max:
							X_grid_max = X_grid[j,m] #SHOULD THESE BE j'S OR i'S??

							# If this is the finest, no point in 
							# continuing looking for finer cells
							if X_grid_max == finest:
								break

					g_val  =  d_l[X_grid_max] # # of repetitions of cell
					G[idx] =  g_val # store value
					idx    += g_val # skip repeats in local coarse cell
					j      += g_val # skip repeats in global cell

				else:
					break
	
			# Compute elements of Phi
			for m in range(nt): 
				idx = 0 # index used for iterating within a coarse cell
				j   = i # index of global matrix

				# Iterate within the coarse cell
				for ii in range(d_0): # dummy loop
					if idx < d_0:

						# Initialize temporary variable to store sum of an 
						# element of Phi
						phi_sum = 0

						# Compute one element in Phi
						for n in range(nt):
							phi_sum += X[j,n] * Psi[n,m]

						g_val = G[idx] # # of repetitions of cell
						

						# Assign value of Phi after dividing by sqrt(lamb_ii)
						# Phi_im[j:j+g_val, m] = phi_sum / np.sqrt(Lambda[m,m])
						# Assign value of Phi after dividing by sqrt(lamb_ii)
						phi_sum = phi_sum / np.sqrt(Lambda[m,m])
						for k in range(j,j+g_val):
							Phi_im[k,m] = phi_sum


						idx += g_val # skip repeats in local coarse cell
						j += g_val # skip repeats in global cell
						
					else:
						break
		# Compute total cpu time
		time_im = time.time() - tic

	# =========================================================================
	# Implemented Computation - Method 2
	# =========================================================================

	elif method == 2:

		# Initialize timer
		tic = time.time()

		# Initialize Phi matrix for implemented computation
		Phi_im = np.empty((nspat, nt))

		# Determine number of levels and finest
		nlev  = len(d_l)
		finest = nlev - 1

		# Specify number of repeats for coarsest and one level
		d_0  = d_l[0]
		d_1  = d_l[1]
		d_f1 = d_l[finest-1]

		# Iterate on the cells that could be new coarse cells
		for i in range(0, nspat, d_0):

			# Initialize matrices to store locations of each level and
			# the number of levels
			G_mat = np.zeros((nlev, d_1, nt), dtype=int)
			nl    = np.zeros((nlev, d_1),     dtype=int)

			# Iterate over columns of X (different times)
			for n in range(nt):

				# Get level of that grid cell
				lvl = X_grid[i,n]

				# Determine if the lvl is the coarsest, if so, 
				# add this cell to G and add another value to nl
				if lvl == 0:
					G_mat[0,0,nl[0,0]] = n
					nl[0,0] += 1
					
				else:
					# Check if the finest is greater than 1. If it is,
					# we will need to look at higher levels
					if finest > 1:
						# Set local and global variables inside coarse cell
						idx = 0
						j = i

						# Iterate over a coarse cell
						for jj in range(d_1):
							if idx < d_1: # dummy loop

								# Get level of this cell
								lvl = X_grid[j,n]

								# Check finest since finest must occur in blocks
								# of d_l we don't look where all finest would be
								if lvl == finest:
									# Tabulate cell and skip corresponding number
									G_mat[finest, idx, nl[finest, idx]] = n
									nl[finest, idx] += 1
									idx += 1
									j += d_f1

								# If not the finest, we need to check each level 
								# from 1 to f-1 iteratively
								else:
									# Tabulate the current level
									G_mat[lvl, idx, nl[lvl, idx]] = n
									nl[lvl, idx] += 1
									idx += d_l[lvl+1]
									j += d_l[lvl]
							else:
								break

					# If not, this cell must be a level one so we 
					# add it to the counts
					else:
						G_mat[1,0,nl[1,0]] = n
						nl[1,0] += 1

			# Compute values of Phi for the elements contained in 
			# coarse cell
			for n in range(nt):

				# Initialize matrix to store contribution of each level
				H = np.zeros((d_0, nlev))

				# Check if we have any l=0 cells
				if nl[0,0] > 0:

					# Initialize temporary variable to store sum of  
					# l=0 contributions to an element of Phi
					l_sum = 0

					# Compute contributions of l=0 cells
					for m in range(nl[0,0]):
						k = G_mat[0,0,m]
						l_sum += X[i,k] * Psi[k,n]

					# Assign this l=0 contribution to H to be summed
					# later
					for m in range(d_0):
						H[m,0] = l_sum

				# Check if we have any cells that aren't l=0
				if nl[0,0] < nt:

					# Check if the finest is greater than 2.
					if finest > 2:

						# Get contributions from l=1 up to l=finest-2
						for l in range(1,finest-1): # 
							idx = 0

							# Iterate on cells that could have new contributions
							for j in range(i, i+d_0, d_l[l]):

								# Check if we have any cells
								if nl[l, idx] > 0:

									# Initialize temporary variable to store sum of  
									# clvl contributions to an element of Phi
									l_sum = 0

									# Compute contribution of current level
									for m in range(nl[l, idx]):
										k = G_mat[l, idx, m]
										l_sum += X[j,k] * Psi[k,n]

									# Assign contribution to H
									for m in range(idx*d_f1, idx*d_f1 + d_l[l]):
										H[m, l] = l_sum

								idx += d_l[l + 1]


					# Check if the finest is greater than 1. This is 
					# where we get contributions from finest-1 and 
					# finest
					if finest > 1:

						idx = 0

						# Compute finest and finest-1 level contribution
						for j in range(i, i + d_0, d_f1):

							# Check if we have any cells at the finest-1 level
							if nl[finest-1, idx] > 0:

								l = finest-1

								# Initialize temporary variable to store sum of  
								# l contributions to an element of Phi
								l_sum = 0

								# Compute contribution of current level
								for m in range(nl[l, idx]):
									k = G_mat[l, idx, m]
									l_sum += X[j,k] * Psi[k,n]

								# Assign contribution to H
								for m in range(idx*d_l[l], (idx+1)*d_l[l]):
									H[m, l] = l_sum

							# Compute finest level contribution
							if nl[finest, idx] > 0:

								# Compute contribution of finest level
								for k in range(j, j+d_f1):

									# Initialize temporary variable to store sum of  
									# l contributions to an element of Phi
									l_sum = 0

									# Compute contribution of current level
									for m in range(nl[finest, idx]):
										p = G_mat[finest, idx, m]
										l_sum += X[k,p] * Psi[p,n]

									# Assign contribution to H
									H[k-j+d_f1*idx, finest] = l_sum

							idx += 1

					# If not, we know the last cells must be l=1
					else:
						# Compute contributions of l=1 cells
						for k in range(i, i+d_0):

							# Initialize temporary variable to store sum of  
							# l=0 contributions to an element of Phi
							l_sum = 0

							# Compute a particular contribution of
							# l=1 forr an element of Phi
							for m in range(nl[1,0]):
								p = G_mat[1, 0, m]
								l_sum += X[k,p] * Psi[p,n]

							# Assign this l=1 contribution to H to be
							# summed later
							H[k-i,1] = l_sum

				# Iterate over spatial locations within coarse cell
				# and compute elements of Phi by summing rows of H
				for m in range(d_0):

					# Initialize temporary variable to store sum of  
					# H rows, which is used for computing elements
					# of Phi
					H_sum = 0
					for l in range(nlev):
						H_sum += H[m,l]

					# Assign value of element of Phi
					Phi_im[m+i, n] = H_sum/np.sqrt(Lambda[n,n])

		# Compute total cpu time
		time_im = time.time() - tic

	# =========================================================================
	# Check Correctness of Matrices
	# =========================================================================
	
	# Check if we should check for correctness
	if type(Phi) != bool:

		# Compute relative error for each cell
		err_im = np.max(abs(np.subtract(Phi_im, Phi)) / abs(Phi))
		err_un = np.max(abs(np.subtract(Phi_un, Phi)) / abs(Phi))

		if err_im < 1e-6:
			print('The implemented Phi is correct')
		else:
			print('The implemented Phi is incorrect')

		if err_un < 1e-6:
			print('The unaltered Phi is correct')
		else:
			print('The unaltered Phi is incorrect')

	# Return CPU time of implemented and unaltered algorithm
	return time_im, time_un





