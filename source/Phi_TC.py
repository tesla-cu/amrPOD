# ================================================= #
# Code:        Computing Phi in POD                 #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np

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
def compute_Phi_TC(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun):

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
	im_fun    = 0			# Function call

	# ========== Unaltered Computation ============================ #

	# Initialize Phi matrix for unaltered computation
	Phi_un = np.zeros((nspat,nt))
	un_asn += wt_asn

	# Compute Phi matrix with unaltered algorithm
	for i in range(nspat):  # iterate through rows of Phi
		un_art += wt_art
		un_acc += wt_acc

		for m in range(nt):  # iterate through columns of Phi
			un_art += wt_art
			un_acc += wt_acc

			# Initialize temporary variable to store sum of an 
			# element of Phi
			phi_sum = 0
			un_asn  += wt_asn

			# Compute inner product
			for k in range(nt):
				un_art += wt_art
				un_acc += wt_acc

				phi_sum += X[i,k] * Psi[k,m]
				un_acc  += 2*wt_acc
				un_art  += 2*wt_art
				un_asn  += wt_asn

			# Assign value of element of Phi
			Phi_un[i,m] = phi_sum/np.sqrt(Lambda[m,m])
			un_acc      += 2*wt_acc
			un_art      += 2*wt_art
			un_asn      += wt_asn
	
	# ========== Implemented Computation - Method 1 =============== #

	if method == 1:

		# Initialize Phi matrix for implemented computation
		Phi_im = np.zeros((nspat, nt))
		im_asn += wt_asn

		# Define number of repetitions for a coarse cell
		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		# Iterate on the cells that could be new coarse cells
		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_acc += wt_acc

			# Initialize matrix to store the maximum grid level at 
			# a particular cell location
			G 	   = np.zeros((d_0), dtype=int)
			im_asn += wt_asn

			idx    = 0 # index used for iterating within a coarse cell
			im_asn += wt_asn

			j      = i   # index of global matrix
			im_asn += wt_asn

			# Find maximum value for each spatial location
			for jj in range(d_0): # dummy loop
				im_art += wt_art
				im_acc += wt_acc

				im_log += wt_log
				if idx < d_0:     # iterate within coarse cell

					# Initialize max grid level
					X_grid_max = X_grid[j,0]
					im_acc     += wt_acc
					im_asn     += wt_asn

					# Find the max grid level for a spatial location
					for m in range(1,nt):
						im_art += wt_art
						im_acc += wt_acc

						# Check if current cell is bigger than current max
						im_acc += wt_acc
						im_log += wt_log
						if X_grid[j,m] > X_grid_max:

							X_grid_max = X_grid[j,m] #SHOULD THESE BE j'S OR i'S??
							im_acc     += wt_acc
							im_asn     += wt_asn

							# If this is the finest, no point in 
							# continuing looking for finer cells
							im_log += wt_log
							if X_grid_max == finest:
								break

					g_val  =  d_l[X_grid_max] # # of repetitions of cell
					im_acc += wt_acc
					im_asn += wt_asn

					G[idx] =  g_val # store value
					im_acc += wt_acc
					im_asn += wt_asn

					idx    += g_val # skip repeats in local coarse cell
					im_art += wt_art
					im_asn += wt_asn

					j      += g_val # skip repeats in global cell
					im_art += wt_art
					im_asn += wt_asn

				else:
					break
	
			# Compute elements of Phi
			for m in range(nt):
				im_art += wt_art
				im_acc += wt_acc

				idx    = 0 # index used for iterating within a coarse cell
				im_asn += wt_asn

				j      = i # index of global matrix
				im_asn += wt_asn

				# Iterate within the coarse cell
				for ii in range(d_0): # dummy loop
					im_art += wt_art
					im_acc += wt_acc

					im_log += wt_log
					if idx < d_0:

						# Initialize temporary variable to store sum of an 
						# element of Phi
						phi_sum = 0
						im_asn  += wt_asn

						# Compute one element in Phi
						for n in range(nt):
							im_art += wt_art
							im_acc += wt_acc

							phi_sum += X[j,n] * Psi[n,m]
							im_acc  += 2*wt_acc
							im_art  += 2*wt_art
							im_asn  += wt_asn

						g_val  = G[idx] # # of repetitions of cell
						im_acc += wt_acc
						im_asn += wt_asn

						# # Assign value of Phi after dividing by sqrt(lamb_ii)
						# Phi_im[j:j+g_val, m] = phi_sum / np.sqrt(Lambda[m,m])
						# im_acc += 2*wt_acc
						# im_art += 3*wt_art
						# im_asn += wt_asn

						# Assign value of Phi after dividing by sqrt(lamb_ii)

						phi_sum = phi_sum / np.sqrt(Lambda[m,m])
						im_acc  += wt_acc
						im_art  += 2*wt_art
						im_asn  += wt_asn

						im_art  += wt_art
						for k in range(j,j+g_val):
							im_art += wt_art
							im_acc += wt_acc

							Phi_im[k,m] = phi_sum
							im_acc += wt_acc
							im_asn += wt_asn

						idx    += g_val # skip repeats in local coarse cell
						im_asn += wt_asn
						im_art += wt_art

						j 	   += g_val # skip repeats in global cell
						im_asn += wt_asn
						im_art += wt_art
						
					else:
						break

	# ========== Implemented Computation - Method 2 =============== #

	elif method == 2:

		# Initialize Phi matrix for implemented computation
		Phi_im = np.zeros((nspat, nt))
		im_asn += wt_asn

		# Determine number of levels and finest
		nlev   = len(d_l)
		im_acc += wt_acc
		im_asn += wt_asn

		finest = nlev - 1
		im_art += wt_art
		im_asn += wt_asn

		# Specify number of repeats for coarsest and one level
		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		d_1 = d_l[1]
		im_acc += wt_acc
		im_asn += wt_asn

		# Iterate on the cells that could be new coarse cells
		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_acc += wt_acc

			# Initialize matrices to store locations of each level and
			# the number of levels
			G_mat  = np.zeros((nlev, d_1, nt), dtype=int)
			im_asn += wt_asn

			nl     = np.zeros((nlev, d_1),     dtype=int)
			im_asn += wt_asn

			# Iterate over columns of X (different times)
			for n in range(nt):
				im_art += wt_art
				im_acc += wt_acc

				# Get level of that grid cell
				lvl    = X_grid[i,n]
				im_acc += wt_acc
				im_asn += wt_asn

				# Determine if the lvl is the coarsest, if so, 
				# add this cell to G and add another value to nl
				im_log += wt_log
				if lvl == 0:

					nl[0,0] += 1
					im_acc  += wt_acc
					im_art  += wt_art
					im_asn  += wt_asn

					G_mat[0,0,nl[0,0]-1] = n
					im_acc += wt_acc
					im_art += wt_art
					im_asn += wt_asn

				else:
					# Check if the finest is greater than 1. If it is,
					# we will need to perform recursion

					im_log += wt_log
					if finest > 1:
						G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun = find_lvl_indices(X_grid, i, 0, n, 1, finest, d_l, G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
						im_fun    += wt_fun
						im_asn    += 2*wt_asn

					# If not, this cell must be a level one so we 
					# add it to the counts
					else:
						nl[1,0] += 1
						im_art  += wt_art
						im_acc  += wt_acc
						im_asn  += wt_asn

						G_mat[1,0,nl[1,0]-1] = n
						im_art += wt_art
						im_acc += 2*wt_acc
						im_asn += wt_asn

			# Compute values of Phi for the elements contained in 
			# coarse cell
			for n in range(nt):
				im_art += wt_art
				im_acc += wt_acc

				# Initialize matrix to store contribution of each level
				H 	   = np.zeros((d_0, nlev))
				im_asn += wt_asn

				# Check if we have any l=0 cells
				im_log += wt_log
				im_acc += wt_acc
				if nl[0,0] > 0:

					# Initialize temporary variable to store sum of  
					# l=0 contributions to an element of Phi
					l_sum  = 0
					im_asn += wt_asn

					# Compute contributions of l=0 cells
					im_acc += wt_acc
					for m in range(nl[0,0]):
						im_art += wt_art
						im_acc += wt_acc

						k      = G_mat[0,0,m]
						im_acc += wt_acc
						im_asn += wt_asn

						l_sum  += X[i,k] * Psi[k,n]
						im_acc += 2*wt_acc
						im_art += 2*wt_art
						im_asn += wt_asn

					# Assign this l=0 contribution to H to be summed
					# later
					for m in range(d_0):
						im_art += wt_art
						im_acc += wt_acc

						H[m,0] = l_sum
						im_acc += wt_acc
						im_asn += wt_asn

				# Check if we have any cells that aren't l=0
				im_log += wt_log
				im_acc += wt_acc
				if nl[0,0] < nt:
					# Check if the finest is greater than 1. If it is,
					# we will need to perform recursion
					im_log += wt_log
					if finest > 1:

						H, im_art, im_acc, im_asn, im_log, im_fun = compute_H(X_grid, Psi, i, 0, n, nt, 1, finest, d_l, G_mat, nl, H, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
						im_fun += wt_fun
						im_asn += wt_asn
					# If not, we know the last cells must be l=1
					else:
						# Compute contributions of l=1 cells
						im_acc += wt_acc
						for k in range(i, i+d_0):
							im_art += wt_art
							im_acc += wt_acc

							# Initialize temporary variable to store sum of  
							# l=0 contributions to an element of Phi
							l_sum  = 0
							im_asn += wt_asn

							# Compute a particular contribution of
							# l=1 forr an element of Phi
							im_acc += wt_acc
							for m in range(nl[1,0]):
								im_art += wt_art
								im_acc += wt_acc

								p      = G_mat[1, 0, m]
								im_acc += wt_acc
								im_asn += wt_asn

								l_sum  += X[k,p] * Psi[p,n]
								im_acc += 2*wt_acc
								im_art += 2*wt_art
								im_asn += wt_asn 

							# Assign this l=1 contribution to H to be
							# summed later
							H[k-i,1] = l_sum
							im_art   += wt_art
							im_acc   += wt_acc
							im_asn   += wt_asn

				# Iterate over spatial locations within coarse cell
				# and compute elements of Phi by summing rows of H
				for m in range(d_0):
					im_art += wt_art
					im_acc += wt_acc

					# Initialize temporary variable to store sum of  
					# H rows, which is used for computing elements
					# of Phi
					H_sum  = 0
					im_asn += wt_asn

					im_art += wt_art
					for l in range(finest+1):
						im_art += wt_art
						im_acc += wt_acc

						H_sum  += H[m,l]
						im_acc += wt_acc
						im_art += wt_art
						im_asn += wt_asn 

					# Assign value of element of Phi
					Phi_im[m+i, n] = H_sum/np.sqrt(Lambda[n,n])
					im_acc += 2*wt_acc
					im_art += 3*wt_art
					im_asn += wt_asn

	if np.max(abs(np.subtract(Phi_im, Phi))) < 1e-8:
		print('The implemented Phi is correct')
	else:
		print('The implemented Phi is incorrect')

	if np.max(abs(np.subtract(Phi_un, Phi))) < 1e-8:
		print('The unaltered Phi is correct')
	else:
		print('The unaltered Phi is incorrect')

	# ========== Sum operations from im and un =================== #

	un = un_asn + un_acc + un_art 
	im = un_art + im_acc + im_asn + im_log + im_fun

	# Return op counts of implemented and unaltered algorithm
	return im, un 


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
def find_lvl_indices(X_grid, i, idx, n, clvl, finest, d_l, G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun):

	# If this current level is 2 or less more than the finest, we 
	# can skip cells  
	im_art += wt_art
	im_log += wt_log
	if clvl < finest-1:

		# Iterate on a cell with repeats corresponding to clvl
		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_acc += wt_acc
			
			# Get level of that grid cell
			lvl    = X_grid[j,n]
			im_acc += wt_acc
			im_asn += wt_asn

			# If that grid level is the current level, tabulate its
			# location and add to nl
			im_log += wt_log
			if lvl == clvl:

				G_mat[clvl, idx, nl[clvl, idx]] = n
				im_acc += 2*wt_acc
				im_asn += wt_asn

				nl[clvl, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

			# If not, we do recursion to look in the appropriate cells
			else:
				G_mat, nl = find_lvl_indices(X_grid, j, idx, n, clvl+1, finest, d_l, G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
				im_fun    += wt_fun
				im_asn    += 2*wt_asn

			# Move index forward for idx within one finer cell
			idx    += d_l[clvl+1]
			im_art += 2*wt_art
			im_acc += wt_acc
			im_asn += wt_asn

	# If the clvl is the finest or one coarser, we know exactly 
	# where these cells are located
	else:

		# Iterate on a cell with repeats corresponding to clvl
		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_acc += wt_acc

			# Get level of that grid cell
			lvl    = X_grid[j,n]
			im_acc += wt_acc
			im_asn += wt_asn 

			# If that grid level is the current level, tabulate its
			# location and add to nl
			im_log += wt_log
			if lvl == clvl:

				G_mat[clvl, idx, nl[clvl, idx]] = n
				im_acc += 2*wt_acc
				im_asn += wt_asn

				nl[clvl, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

			# If not current level, it must be the finest. Tabulate.
			else:
				nl[finest, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

				G_mat[finest, idx, nl[finest, idx]-1] = n
				im_art += wt_art
				im_acc += 2*wt_acc
				im_asn += wt_asn

			# Move index forward for number of repeats of finest, 
			# which is always 1
			idx    += 1
			im_art += wt_art
			im_asn += wt_asn

	# Return updated G_mat and nl
	return G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun

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
def compute_H(X, Psi, i, idx, n, nt, clvl, finest, d_l, G_mat, nl, H, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun):
	
	# If this current level is 2 or less more than the finest, we 
	# can skip cells  
	im_log += wt_log
	im_art += wt_art
	if clvl < finest-1:

		# Iterate on a cell with repeats corresponding to clvl
		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_acc += wt_acc

			# Check if we have any cells at the current level
			im_log += wt_log
			im_acc += wt_acc
			if nl[clvl, idx] > 0:

				# Initialize temporary variable to store sum of  
				# clvl contributions to an element of Phi
				l_sum  = 0
				im_asn += wt_asn

				# Compute contribution of current level
				im_acc += wt_acc
				for m in range(nl[clvl, idx]):
					im_art += wt_art
					im_acc += wt_acc

					k = G_mat[clvl, idx, m]
					im_acc += wt_acc
					im_asn += wt_asn

					l_sum  += X[j,k] * Psi[k,n]
					im_art += 2*wt_art
					im_acc += 2*wt_acc
					im_asn += wt_asn

				# Temporary index needed for assigning 
				# contribution to correct cells in H
				# l_idx  = int((j-i)/d_l[clvl])
				# im_art += 2*wt_art
				# im_acc += wt_acc
				# im_asn += wt_asn

				# Assign contribution to H
				im_art += 3*wt_art
				im_acc += 2*wt_acc

				# for m in range(l_idx*d_l[clvl], (l_idx+1)*d_l[clvl]):
				for m in range(idx*d_l[finest-1], idx*d_l[finest-1] + d_l[clvl]):
					im_art += wt_art
					im_acc += wt_acc

					H[m, clvl] = l_sum
					im_acc     += wt_acc
					im_asn     += wt_asn

			# Determine if we need to look at a finer level
			nccells = 0
			im_asn  += wt_asn

			im_art  += wt_art
			for l in range(clvl+1):
				im_art += wt_art
				im_acc += wt_acc

				nccells += nl[l, idx]
				im_acc  += wt_acc
				im_art  += wt_art
				im_asn  += wt_asn

			# If this is true, we need to check contributions of 
			# finer cells
			im_log += wt_log
			if nccells < nt:

				H, im_art, im_acc, im_asn, im_log, im_fun = compute_H(X, Psi, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
				im_fun += wt_fun
				im_asn += wt_asn

			# Move index forward for idx within one finer cell
			idx    += d_l[clvl + 1]
			im_acc += wt_acc
			im_art += 2*wt_art
			im_asn += wt_asn

	# If this current level is the finest or one coarse, we know 
	# where the rest of the cells reside
	else:

		# Iterate on a cell with repeats corresponding to clvl
		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_acc += wt_acc

			# Check if we have any cells at the current level
			im_log += wt_log
			im_acc += wt_acc
			if nl[clvl, idx] > 0:

				# Initialize temporary variable to store sum of  
				# clvl contributions to an element of Phi
				l_sum  = 0
				im_asn += wt_asn

				# Compute contribution of current level
				im_acc += wt_acc
				for m in range(nl[clvl, idx]):
					im_art += wt_art
					im_acc += wt_acc

					k = G_mat[clvl, idx, m]
					im_acc += wt_acc
					im_asn += wt_asn

					l_sum  += X[j,k] * Psi[k,n]
					im_acc += 2*wt_acc
					im_art += 2*wt_art
					im_asn += wt_asn

				# Assign contribution to H
				im_acc += 2*wt_acc
				im_art += 3*wt_art
				for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
					im_art += wt_art
					im_acc += wt_acc

					H[m, clvl] = l_sum
					im_acc += wt_acc
					im_asn += wt_asn

			# Check if we have any cells at the finest level
			im_log += wt_log
			im_acc += wt_acc
			if nl[finest, idx] > 0:

				# Compute contribution of current level
				im_art += wt_art
				im_acc += wt_acc
				for k in range(j, j+d_l[clvl]):
					im_art += wt_art
					im_acc += wt_acc

					# Initialize temporary variable to store sum of  
					# clvl contributions to an element of Phi
					l_sum  = 0
					im_asn += wt_asn

					# Compute contribution of current level
					im_acc += wt_acc
					for m in range(nl[finest, idx]):
						im_art += wt_art
						im_acc += wt_acc

						p = G_mat[finest, idx, m]
						im_acc += wt_acc
						im_asn += wt_asn

						l_sum  += X[k,p] * Psi[p,n]
						im_art += 2*wt_art
						im_acc += 2*wt_acc
						im_asn += wt_asn

					# Assign contribution to H
					H[k-j+d_l[finest-1] * (idx), finest] = l_sum
					im_acc += 2*wt_acc
					im_art += 4*wt_art
					im_asn += wt_asn

			# Move index forward for number of repeats of finest, 
			# which is always 1
			idx    += 1
			im_art += wt_art
			im_asn += wt_asn

	return H, im_art, im_acc, im_asn, im_log, im_fun



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
	

