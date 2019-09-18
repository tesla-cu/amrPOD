# ================================================= #
# Code:        Computing Phi in POD                 #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
import time

# ================================================================= #
# Function to compute the POD spatial modes in using a standard 
# matrix operation technique and the new algorithm leveraging AMR
# repetitions
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
#            (this is used as a check we did the computation right)
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - nt     : number of time steps
# - nspat  : number of spatial locations
# - finest : finest AMR grid level
#
# Outputs:
# - time_im : CPU time to compute Phi using implemented algorithm
# - time_un : CPU time to compute Phi using unaltered algorithm
# ================================================================= #
def compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest):

	# ========== Unaltered Computation ============================ #

	# Initialize timer
	tic = time.time()

	# Initialize Phi matrix for unaltered computation
	Phi_un = np.zeros((nspat,nt))

	# Compute Phi matrix with unaltered algorithm
	for i in range(nspat):  # iterate through rows of Phi
		for m in range(nt):  # iterate through columns of Phi

			# Initialize temporary variable to store sum of an 
			# element of Phi
			phi_sum = 0

			# Compute inner product
			for k in range(nt):
				phi_sum += X[i,k] * Psi[k,m]

			# Assign value of element of Phi
			Phi_un[i,m] = phi_sum/np.sqrt(Lambda[m,m])

	# Compute total cpu time
	time_un = time.time() - tic
	
	# ========== Implemented Computation - Method 1 =============== #

	if method == 1:

		# Initialize timer
		tic = time.time()

		# Initialize Phi matrix for implemented computation
		Phi_im = np.zeros((nspat, nt))

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

	# ========== Implemented Computation - Method 2 =============== #

	elif method == 2:

		# Initialize timer
		tic = time.time()

		# Initialize Phi matrix for implemented computation
		Phi_im = np.zeros((nspat, nt))

		# Determine number of levels and finest
		nlev  = len(d_l)
		finest = nlev - 1

		# Specify number of repeats for coarsest and one level
		d_0 = d_l[0]
		d_1 = d_l[1]

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
					# we will need to perform recursion
					if finest > 1:
						G_mat, nl = find_lvl_indices(X_grid, i, 0, n, 1, finest, d_l, G_mat, nl)
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
					# Check if the finest is greater than 1. If it is,
					# we will need to perform recursion
					if finest > 1:
						H = compute_H(X_grid, Psi, i, 0, n, nt, 1, finest, d_l, G_mat, nl, H)
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
					for l in range(finest+1):
						H_sum += H[m,l]

					# Assign value of element of Phi
					Phi_im[m+i, n] = H_sum/np.sqrt(Lambda[n,n])

		# Compute total cpu time
		time_im = time.time() - tic

	if np.max(abs(np.subtract(Phi_im, Phi))) < 1e-8:
		print('The implemented Phi is correct')
	else:
		print('The implemented Phi is incorrect')

	if np.max(abs(np.subtract(Phi_un, Phi))) < 1e-8:
		print('The unaltered Phi is correct')
	else:
		print('The unaltered Phi is incorrect')

	return time_im, time_un




# ================================================================= #
# Function used by compute_Phi_CPU to tabulate all cells according 
# to grid level. This is a recursive algorithm that determines the
# exact location to find these cells without looking in unncessary
# locations.
#
# Inputs:
# - X_grid : companion matrix to snapshot matrix that stores grid 
#            levels instead of solution values
# - i      : the global index in X that we are looking at
# - idx    : index of the cell within a given coarse cell
# - n      : snapshot or column of X we are considering
# - clvl   : current level that we consider the coarsest cell
# - finest : finest AMR grid level
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - G_mat  : matrix to store locations of each level
# - nl     : vector to store number of cells for each level
#
# Outputs:
# - G_mat  : (as above) but updated with clvl locations
# - nl     : (as above) but updated with number of cells at clvl
# ================================================================= #
def find_lvl_indices(X_grid, i, idx, n, clvl, finest, d_l, G_mat, nl):

	# If this current level is 2 or less more than the finest, we 
	# can skip cells  
	if clvl < finest-1:

		# Iterate on a cell with repeats corresponding to clvl
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			
			# Get level of that grid cell
			lvl = X_grid[j,n]

			# If that grid level is the current level, tabulate its
			# location and add to nl
			if lvl == clvl:
				G_mat[clvl, idx, nl[clvl, idx]] = n
				nl[clvl, idx] += 1

			# If not, we do recursion to look in the appropriate cells
			else:
				G_mat, nl = find_lvl_indices(X_grid, j, idx, n, clvl+1, finest, d_l, G_mat, nl)

			# Move index forward for idx within one finer cell
			idx += d_l[clvl+1]

	# If the clvl is the finest or one coarser, we know exactly 
	# where these cells are located
	else:

		# Iterate on a cell with repeats corresponding to clvl
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

			# Get level of that grid cell
			lvl = X_grid[j,n]

			# If that grid level is the current level, tabulate its
			# location and add to nl
			if lvl == clvl:
				G_mat[clvl, idx, nl[clvl, idx]] = n
				nl[clvl, idx] += 1
			# If not current level, it must be the finest. Tabulate.
			else:
				G_mat[finest, idx, nl[finest, idx]] = n
				nl[finest, idx] += 1

			# Move index forward for number of repeats of finest, 
			# which is always 1
			idx += 1

	# Return updated G_mat and nl
	return G_mat, nl

# ================================================================= #
# Function used by compute_Phi_CPU to tabulate all cells according 
# to grid level. This is a recursive algorithm that determines the
# exact location to find these cells without looking in unncessary
# locations.
#
# Inputs:
# - X      : snapshot matrix
# - Psi    : matrix containing eigenvectors of R
# - i      : the global index in X that we are looking at
# - idx    : index of the cell within a given coarse cell
# - n      : snapshot or column of X we are considering
# - clvl   : current level that we consider the coarsest cell
# - finest : finest AMR grid level
# - d_l    : number of repeated cells for a given level l (called 
#            c_\ell^d in the paper)
# - G_mat  : matrix to store locations of each level
# - nl     : vector to store number of cells for each level
# - H      : matrix storing contributions to an element of Phi for
#            each level
# 
# Outputs:
# - H      : (as above) but updated with clvl contributions to H
# ================================================================= #
def compute_H(X, Psi, i, idx, n, nt, clvl, finest, d_l, G_mat, nl, H):
	
	# If this current level is 2 or less more than the finest, we 
	# can skip cells  
	if clvl < finest-1:

		# Iterate on a cell with repeats corresponding to clvl
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

			# Check if we have any cells at the current level
			if nl[clvl, idx] > 0:

				# Initialize temporary variable to store sum of  
				# clvl contributions to an element of Phi
				l_sum = 0

				# Compute contribution of current level
				for m in range(nl[clvl, idx]):
					k = G_mat[clvl, idx, m]
					l_sum += X[j,k] * Psi[k,n]

				# Temporary index needed for assigning 
				# contribution to correct cells in H
				l_idx = int((j-i)/d_l[clvl])

				# Assign contribution to H
				for m in range(l_idx*d_l[clvl], (l_idx+1)*d_l[clvl]):
					H[m, clvl] = l_sum

			# Determine if we need to look at a finer level
			nccells = 0
			for l in range(clvl+1):
				nccells += nl[l, idx]

			# If this is true, we need to check contributions of 
			# finer cells
			if nccells < nt:
				H = compute_H(X, Psi, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H)

			# Move index forward for idx within one finer cell
			idx += d_l[clvl + 1]
	
	# If this current level is the finest or one coarse, we know 
	# where the rest of the cells reside
	else:

		# Iterate on a cell with repeats corresponding to clvl
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

			# Check if we have any cells at the current level
			if nl[clvl, idx] > 0:

				# Initialize temporary variable to store sum of  
				# clvl contributions to an element of Phi
				l_sum = 0

				# Compute contribution of current level
				for m in range(nl[clvl, idx]):
					k = G_mat[clvl, idx, m]
					l_sum += X[j,k] * Psi[k,n]

				# Assign contribution to H
				for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
					H[m, clvl] = l_sum

			# Check if we have any cells at the finest level
			if nl[finest, idx] > 0:

				# Compute contribution of current level
				for k in range(j, j+d_l[clvl]):

					# Initialize temporary variable to store sum of  
					# clvl contributions to an element of Phi
					l_sum = 0

					# Compute contribution of current level
					for m in range(nl[finest, idx]):
						p = G_mat[finest, idx, m]
						l_sum += X[k,p] * Psi[p,n]

					# Assign contribution to H
					H[k-j+d_l[finest-1] * (idx), finest] = l_sum

			# Move index forward for number of repeats of finest, 
			# which is always 1
			idx += 1

	return H



# Old method 1
# if method == 1:

# 	Phi_imp = np.zeros((nspat, nt))
# 	G       = np.zeros((nspat), dtype=int)
# 	i       = 0

# 	for n in range(nspat):
# 		if i < nspat:
# 			X_grid_max = X_grid[i,0]
# 			for j in range(1,nt):
# 				if X_grid[i,j] > X_grid_max:
# 					X_grid_max = X_grid[i,j]
# 			G[i]  = d_l[X_grid_max]
# 			i     += G[i]
# 		else:
# 			break

# 	for i in range(nt):
# 		j = 0
# 		for n in range(nspat):
# 			if j < nspat:
# 				phi_sum = 0
# 				for k in range(nt):
# 					phi_sum += X[j,k] * Psi[k,i]
# 				Phi_imp[j:j+G[j], i] = phi_sum / np.sqrt(Lambda[i,i])
# 				j = j + G[j]
	

