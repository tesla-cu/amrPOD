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
# repetitions that recursively finds unique cells.
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
# - wt_art : weighting of arithmetic operations
# - wt_acc : weighting of accessing operations
# - wt_asn : weighting of assignment operations
# - wt_log : weighting of logical operations
# - wt_fun : weighting of function calls
#
# Outputs:
# - time_im : CPU time to compute Phi using implemented algorithm
# - time_un : CPU time to compute Phi using unaltered algorithm
# ================================================================= #
def compute_Phi_TC(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun):

	# ========== Initialize Operation Counts ====================== #

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

	# ========== Unaltered Computation ============================ #

	Phi_un = np.zeros((nspat,nt))
	un_asn += wt_asn
	un_fun += wt_fun

	for i in range(nspat):
		un_art += wt_art
		un_asn += wt_asn

		for m in range(nt):
			un_art += wt_art
			un_asn += wt_asn

			phi_sum = 0
			un_asn  += wt_asn

			for k in range(nt):
				un_art += wt_art
				un_asn += wt_asn

				phi_sum += X[i,k] * Psi[k,m]
				un_acc  += 2*wt_acc
				un_art  += 2*wt_art
				un_asn  += wt_asn

			Phi_un[i,m] = phi_sum/np.sqrt(Lambda[m,m])
			un_acc      += 2*wt_acc
			un_art      += 2*wt_art
			un_asn      += wt_asn
			un_fun      += wt_fun
	
	# ========== Implemented Computation - Method 1 =============== #

	if method == 1:

		Phi_im = np.zeros((nspat, nt))
		im_asn += wt_asn
		im_fun += im_fun

		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_asn += wt_asn

			G 	   = np.zeros((d_0), dtype=int)
			im_asn += wt_asn
			im_fun += wt_fun

			idx    = 0 
			im_asn += wt_asn

			j      = i  
			im_asn += wt_asn

			for jj in range(d_0): 
				im_art += wt_art
				im_asn += wt_asn

				im_log += wt_log
				if idx < d_0:

					X_grid_max = X_grid[j,0]
					im_acc     += wt_acc
					im_asn     += wt_asn

					for m in range(1,nt):
						im_art += wt_art
						im_asn += wt_asn

						im_acc += wt_acc
						im_log += wt_log
						if X_grid[j,m] > X_grid_max:

							X_grid_max = X_grid[j,m]
							im_acc     += wt_acc
							im_asn     += wt_asn

							im_log += wt_log
							if X_grid_max == finest:
								break

					g_val  =  d_l[X_grid_max] 
					im_acc += wt_acc
					im_asn += wt_asn

					G[idx] =  g_val 
					im_acc += wt_acc
					im_asn += wt_asn

					idx    += g_val 
					im_art += wt_art
					im_asn += wt_asn

					j      += g_val 
					im_art += wt_art
					im_asn += wt_asn

				else:
					break
	
			for m in range(nt):
				im_art += wt_art
				im_asn += wt_asn

				idx    = 0
				im_asn += wt_asn

				j      = i 
				im_asn += wt_asn

				for ii in range(d_0):
					im_art += wt_art
					im_asn += wt_asn

					im_log += wt_log
					if idx < d_0:

						phi_sum = 0
						im_asn  += wt_asn

						for n in range(nt):
							im_art += wt_art
							im_asn += wt_asn

							phi_sum += X[j,n] * Psi[n,m]
							im_acc  += 2*wt_acc
							im_art  += 2*wt_art
							im_asn  += wt_asn

						g_val  = G[idx]
						im_acc += wt_acc
						im_asn += wt_asn

						phi_sum = phi_sum / np.sqrt(Lambda[m,m])
						im_acc  += wt_acc
						im_art  += 2*wt_art
						im_asn  += wt_asn
						im_fun  += wt_fun

						im_art  += wt_art
						for k in range(j,j+g_val):
							im_art += wt_art
							im_asn += wt_asn

							Phi_im[k,m] = phi_sum
							im_acc += wt_acc
							im_asn += wt_asn

						idx    += g_val 
						im_asn += wt_asn
						im_art += wt_art

						j 	   += g_val
						im_asn += wt_asn
						im_art += wt_art
						
					else:
						break

	# ========== Implemented Computation - Method 2 =============== #

	elif method == 2:

		Phi_im = np.zeros((nspat, nt))
		im_asn += wt_asn
		im_fun += im_fun

		nlev   = finest + 1
		im_art += wt_art
		im_asn += wt_asn

		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		d_1 = d_l[1]
		im_acc += wt_acc
		im_asn += wt_asn

		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_asn += wt_asn

			G_mat  = np.zeros((nlev, d_1, nt), dtype=int)
			im_asn += wt_asn
			im_fun += wt_fun

			nl     = np.zeros((nlev, d_1), dtype=int)
			im_asn += wt_asn
			im_fun += wt_fun

			for n in range(nt):
				im_art += wt_art
				im_asn += wt_asn

				lvl    = X_grid[i,n]
				im_acc += wt_acc
				im_asn += wt_asn

				im_log += wt_log
				if lvl == 0:

					G_mat[0,0,nl[0,0]] = n
					im_acc += 2*wt_acc
					im_asn += wt_asn

					nl[0,0] += 1
					im_acc  += wt_acc
					im_art  += wt_art
					im_asn  += wt_asn


				else:

					im_log += wt_log
					if finest > 1:
						G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun = find_lvl_indices(X_grid, i, 0, n, 1, finest, d_l, G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
						im_fun    += wt_fun
						im_asn    += 2*wt_asn

					else:
						G_mat[1,0,nl[1,0]] = n
						im_acc += 2*wt_acc
						im_asn += wt_asn

						nl[1,0] += 1
						im_art  += wt_art
						im_acc  += wt_acc
						im_asn  += wt_asn


			for n in range(nt):
				im_art += wt_art
				im_asn += wt_asn

				H 	   = np.zeros((d_0, nlev))
				im_asn += wt_asn
				im_fun += wt_fun

				im_log += wt_log
				im_acc += wt_acc
				if nl[0,0] > 0:

					l_sum  = 0
					im_asn += wt_asn

					im_acc += wt_acc
					for m in range(nl[0,0]):
						im_art += wt_art
						im_asn += wt_asn

						k      = G_mat[0,0,m]
						im_acc += wt_acc
						im_asn += wt_asn

						l_sum  += X[i,k] * Psi[k,n]
						im_acc += 2*wt_acc
						im_art += 2*wt_art
						im_asn += wt_asn

					for m in range(d_0):
						im_art += wt_art
						im_asn += wt_asn

						H[m,0] = l_sum
						im_acc += wt_acc
						im_asn += wt_asn

				im_log += wt_log
				im_acc += wt_acc
				if nl[0,0] < nt:

					im_log += wt_log
					if finest > 1:

						H, im_art, im_acc, im_asn, im_log, im_fun = compute_H(X_grid, Psi, i, 0, n, nt, 1, finest, d_l, G_mat, nl, H, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
						im_fun += wt_fun
						im_asn += wt_asn

					else:

						im_art += wt_art
						for k in range(i, i+d_0):
							im_art += wt_art
							im_asn += wt_asn

							l_sum  = 0
							im_asn += wt_asn

							im_acc += wt_acc
							for m in range(nl[1,0]):
								im_art += wt_art
								im_asn += wt_asn

								p      = G_mat[1, 0, m]
								im_acc += wt_acc
								im_asn += wt_asn

								l_sum  += X[k,p] * Psi[p,n]
								im_acc += 2*wt_acc
								im_art += 2*wt_art
								im_asn += wt_asn 

							H[k-i,1] = l_sum
							im_art   += wt_art
							im_acc   += wt_acc
							im_asn   += wt_asn

				for m in range(d_0):
					im_art += wt_art
					im_asn += wt_asn

					H_sum  = 0
					im_asn += wt_asn

					im_art += wt_art
					for l in range(finest+1):
						im_art += wt_art
						im_asn += wt_asn

						H_sum  += H[m,l]
						im_acc += wt_acc
						im_art += wt_art
						im_asn += wt_asn 

					Phi_im[m+i, n] = H_sum/np.sqrt(Lambda[n,n])
					im_acc += 2*wt_acc
					im_art += 3*wt_art
					im_asn += wt_asn
					im_fun += wt_fun

	# ========== Check Correctness of Matrices ==================== #

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

	# ========== Sum operations from im and un =================== #

	un = un_asn + un_acc + un_art + un_fun
	im = un_art + im_acc + im_asn + im_fun + im_log

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

	im_art += wt_art
	im_log += wt_log
	if clvl < finest-1:

		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_asn += wt_asn
			
			lvl    = X_grid[j,n]
			im_acc += wt_acc
			im_asn += wt_asn

			im_log += wt_log
			if lvl == clvl:

				G_mat[clvl, idx, nl[clvl, idx]] = n
				im_acc += 2*wt_acc
				im_asn += wt_asn

				nl[clvl, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

			else:
				G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun = find_lvl_indices(X_grid, j, idx, n, clvl+1, finest, d_l, G_mat, nl, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
				im_fun    += wt_fun
				im_asn    += 2*wt_asn

			idx    += d_l[clvl+1]
			im_art += 2*wt_art
			im_acc += wt_acc
			im_asn += wt_asn

	else:

		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_asn += wt_asn

			lvl    = X_grid[j,n]
			im_acc += wt_acc
			im_asn += wt_asn 

			im_log += wt_log
			if lvl == clvl:

				G_mat[clvl, idx, nl[clvl, idx]] = n
				im_acc += 2*wt_acc
				im_asn += wt_asn

				nl[clvl, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

			else:
				G_mat[finest, idx, nl[finest, idx]] = n
				im_acc += 2*wt_acc
				im_asn += wt_asn

				nl[finest, idx] += 1
				im_acc += wt_acc
				im_art += wt_art
				im_asn += wt_asn

			idx    += 1
			im_art += wt_art
			im_asn += wt_asn

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
	
	im_log += wt_log
	im_art += wt_art
	if clvl < finest-1:

		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_asn += wt_asn

			im_log += wt_log
			im_acc += wt_acc
			if nl[clvl, idx] > 0:

				l_sum  = 0
				im_asn += wt_asn

				im_acc += wt_acc
				for m in range(nl[clvl, idx]):
					im_art += wt_art
					im_asn += wt_asn

					k = G_mat[clvl, idx, m]
					im_acc += wt_acc
					im_asn += wt_asn

					l_sum  += X[j,k] * Psi[k,n]
					im_art += 2*wt_art
					im_acc += 2*wt_acc
					im_asn += wt_asn

				im_art += 5*wt_art
				im_acc += 3*wt_acc
				for m in range(idx*d_l[finest-1], idx*d_l[finest-1] + d_l[clvl]):
					im_art += wt_art
					im_asn += wt_asn

					H[m, clvl] = l_sum
					im_acc     += wt_acc
					im_asn     += wt_asn

			nccells = 0
			im_asn  += wt_asn

			im_art  += wt_art
			for l in range(clvl+1):
				im_art += wt_art
				im_asn += wt_asn

				nccells += nl[l, idx]
				im_acc  += wt_acc
				im_art  += wt_art
				im_asn  += wt_asn

			im_log += wt_log
			if nccells < nt:

				H, im_art, im_acc, im_asn, im_log, im_fun = compute_H(X, Psi, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H, im_art, im_acc, im_asn, im_log, im_fun, wt_art, wt_acc, wt_asn, wt_log, wt_fun)
				im_fun += wt_fun
				im_asn += wt_asn

			idx    += d_l[clvl + 1]
			im_acc += wt_acc
			im_art += 2*wt_art
			im_asn += wt_asn

	else:

		im_art += 3*wt_art
		im_acc += 2*wt_acc
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			im_art += wt_art
			im_asn += wt_asn

			im_log += wt_log
			im_acc += wt_acc
			if nl[clvl, idx] > 0:

				l_sum  = 0
				im_asn += wt_asn

				im_acc += wt_acc
				for m in range(nl[clvl, idx]):
					im_art += wt_art
					im_asn += wt_asn

					k = G_mat[clvl, idx, m]
					im_acc += wt_acc
					im_asn += wt_asn

					l_sum  += X[j,k] * Psi[k,n]
					im_acc += 2*wt_acc
					im_art += 2*wt_art
					im_asn += wt_asn

				im_acc += 2*wt_acc
				im_art += 3*wt_art
				for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
					im_art += wt_art
					im_asn += wt_asn

					H[m, clvl] = l_sum
					im_acc += wt_acc
					im_asn += wt_asn

			im_log += wt_log
			im_acc += wt_acc
			if nl[finest, idx] > 0:

				im_art += wt_art
				im_acc += wt_acc
				for k in range(j, j+d_l[clvl]):
					im_art += wt_art
					im_asn += wt_asn

					l_sum  = 0
					im_asn += wt_asn

					im_acc += wt_acc
					for m in range(nl[finest, idx]):
						im_art += wt_art
						im_asn += wt_asn

						p = G_mat[finest, idx, m]
						im_acc += wt_acc
						im_asn += wt_asn

						l_sum  += X[k,p] * Psi[p,n]
						im_art += 2*wt_art
						im_acc += 2*wt_acc
						im_asn += wt_asn

					H[k-j+(d_l[finest-1]*idx), finest] = l_sum
					im_acc += 2*wt_acc
					im_art += 4*wt_art
					im_asn += wt_asn

			idx    += 1
			im_art += wt_art
			im_asn += wt_asn

	return H, im_art, im_acc, im_asn, im_log, im_fun

