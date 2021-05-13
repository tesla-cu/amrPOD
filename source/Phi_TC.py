# ================================================= #
# Code:        Computing Phi in POD                 #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
from numba import njit

# =========================================================================== #
# Function to compute the POD spatial modes in using a standard  matrix 
# operation technique and the new algorithm leveraging AMR repetitions that
# iteratively finds unique cells.
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
# - wt_art : weighting of arithmetic operations
# - wt_acc : weighting of accessing operations
# - wt_asn : weighting of assignment operations
# - wt_log : weighting of logical operations
# - wt_fun : weighting of function calls
#
# Outputs:
# - im : num operations to compute Phi using implemented algorithm
# - un : num operations to compute Phi using unaltered algorithm
# =========================================================================== #
@njit
def compute_Phi_TC(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest,\
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

	Phi_un = np.empty((nspat,nt))
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

			for n in range(nt):
				un_art += wt_art
				un_asn += wt_asn

				phi_sum += X[i,n] * Psi[n,m]
				un_acc  += 2*wt_acc
				un_art  += 2*wt_art
				un_asn  += wt_asn

			Phi_un[i,m] = phi_sum/np.sqrt(Lambda[m,m])
			un_acc      += 2*wt_acc
			un_art      += 2*wt_art
			un_asn      += wt_asn
			un_fun      += wt_fun
	
	# =========================================================================
	# Implemented Computation - Method 1
	# =========================================================================

	if method == 1:

		Phi_im = np.empty((nspat, nt))
		im_asn += wt_asn
		im_fun += wt_fun

		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_asn += wt_asn

			G 	   = np.zeros((d_0))
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

						g_val  = int(G[idx])
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

	# =========================================================================
	# Implemented Computation - Method 2 
	# =========================================================================

	elif method == 2:

		Phi_im = np.empty((nspat, nt))
		im_asn += wt_asn
		im_fun += wt_fun

		nlev   = finest + 1
		im_art += wt_art
		im_asn += wt_asn

		d_0    = d_l[0]
		im_acc += wt_acc
		im_asn += wt_asn

		d_1 = d_l[1]
		im_acc += wt_acc
		im_asn += wt_asn

		d_f1 = d_l[finest-1]
		im_asn += wt_asn
		im_art += wt_art
		im_acc += wt_acc

		for i in range(0, nspat, d_0):
			im_art += wt_art
			im_asn += wt_asn

			G_mat  = np.zeros((nlev, d_1, nt))
			im_asn += wt_asn
			im_fun += wt_fun

			nl     = np.zeros((nlev, d_1))
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

					G_mat[0,0,int(nl[0,0])] = n
					im_acc += 2*wt_acc
					im_asn += wt_asn

					nl[0,0] += 1
					im_acc  += wt_acc
					im_art  += wt_art
					im_asn  += wt_asn

				else:

					im_log += wt_log
					if finest > 1:

						idx = 0
						im_asn += wt_asn

						j = i
						im_asn += wt_asn

						for jj in range(d_1):
							im_art += wt_art
							im_asn += wt_asn

							im_log += wt_log
							if idx < d_1:

								lvl = X_grid[j,n]
								im_acc += wt_acc
								im_asn += wt_asn

								im_log += wt_log
								if lvl == finest:

									G_mat[finest, idx, int(nl[finest, idx])] = n
									im_acc += 2*wt_acc
									im_asn += wt_asn

									nl[finest, idx] += 1
									im_acc += wt_acc
									im_art += wt_art
									im_asn += wt_asn

									idx += 1
									im_art += wt_art
									im_asn += wt_asn

									j += d_f1
									im_art += wt_art
									im_asn += wt_asn

								else:
									G_mat[lvl, idx, int(nl[lvl, idx])] = n
									im_acc += 2*wt_acc
									im_asn += wt_asn

									nl[lvl, idx] += 1
									im_acc += wt_acc
									im_art += wt_art
									im_asn += wt_asn

									idx += d_l[lvl+1]
									im_art += 2*wt_art
									im_acc += wt_acc
									im_asn += wt_asn

									j += d_l[lvl]
									im_acc += wt_acc
									im_art += wt_art
									im_asn += wt_asn

							else:
								break

					else:
						G_mat[1,0,int(nl[1,0])] = n
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

						k      = int(G_mat[0,0,m])
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
					if finest > 2:

						im_art += wt_art
						for l in range(1,finest-1): 
							im_art += wt_art
							im_asn += wt_asn

							idx = 0
							im_asn += wt_asn

							im_art += wt_art
							im_acc += wt_acc
							for j in range(i, i+d_0, d_l[l]):
								im_art += wt_art
								im_asn += wt_asn

								im_log += wt_log
								im_acc += wt_acc
								if nl[l, idx] > 0:

									l_sum = 0
									im_asn += wt_asn

									im_acc += wt_acc
									for m in range(nl[l, idx]):
										im_art += wt_art
										im_asn += wt_asn

										k = int(G_mat[l, idx, m])
										im_acc += wt_acc
										im_asn += wt_asn

										l_sum += X[j,k] * Psi[k,n]
										im_acc += 2*wt_acc
										im_art += 2*wt_art
										im_asn += wt_asn

									im_art += 3*wt_art
									im_acc += wt_acc
									for m in range(idx*d_f1, idx*d_f1 + d_l[l]):
										im_art += wt_art
										im_asn += wt_asn

										H[m, l] = l_sum
										im_acc += wt_acc
										im_asn += wt_asn

								idx += d_l[l + 1]
								im_acc += wt_acc
								im_art += 2*wt_art
								im_asn += wt_asn

					im_log += wt_log
					if finest > 1:

						idx = 0
						im_asn += wt_asn

						im_art += wt_art
						for j in range(i, i + d_0, d_f1):
							im_art += wt_art
							im_asn += wt_asn

							im_acc += wt_acc
							im_art += wt_art
							im_log += wt_log
							if nl[finest-1, idx] > 0:

								l = finest-1
								im_art += wt_art
								im_asn += wt_asn

								l_sum = 0
								im_asn += wt_asn

								im_acc += wt_acc
								for m in range(int(nl[l, idx])):
									im_art += wt_art
									im_asn += wt_asn

									k = int(G_mat[l, idx, m])
									im_acc += wt_acc
									im_asn += wt_asn

									l_sum += X[j,k] * Psi[k,n]
									im_acc += 2*wt_acc
									im_asn += wt_asn
									im_art += 2*wt_art


								im_art += 3*wt_art
								im_acc += 2*wt_acc
								for m in range(idx*d_l[l], (idx+1)*d_l[l]):
									im_art += wt_art
									im_asn += wt_asn

									H[m, l] = l_sum
									im_acc += wt_acc
									im_asn += wt_asn

							im_log += wt_log
							im_acc += wt_acc
							if nl[finest, idx] > 0:

								im_art += wt_art
								for k in range(j, j+d_f1):
									im_art += wt_art
									im_asn += wt_asn

									l_sum = 0
									im_asn += wt_asn

									im_acc += wt_acc
									for m in range(int(nl[finest, idx])):
										im_art += wt_art
										im_asn += wt_asn

										p = int(G_mat[finest, idx, m])
										im_acc += wt_acc
										im_asn += wt_asn

										l_sum += X[k,p] * Psi[p,n]
										im_acc += 2*wt_acc
										im_asn += wt_asn
										im_art += 2*wt_art

									H[k-j+d_f1*idx, finest] = l_sum
									im_art += 3*wt_art
									im_acc += wt_acc
									im_asn += wt_asn

							idx += 1
							im_art += wt_art
							im_asn += wt_asn

					else:

						im_art += wt_art
						for k in range(i, i+d_0):
							im_art += wt_art
							im_asn += wt_asn

							l_sum  = 0
							im_asn += wt_asn

							im_acc += wt_acc
							for m in range(int(nl[1,0])):
								im_art += wt_art
								im_asn += wt_asn

								p      = int(G_mat[1, 0, m])
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

					for l in range(nlev):
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

	# =========================================================================
	# Check Correctness of Matrices
	# =========================================================================

	# Check if we should check for correctness
	if Phi is not None:

		# Compute relative error for each cell
		err_im = np.max(np.abs(np.subtract(Phi_im, Phi)) / np.abs(Phi))
		err_un = np.max(np.abs(np.subtract(Phi_un, Phi)) / np.abs(Phi))

		if err_im < 1e-6:
			print('The implemented Phi is correct')
		else:
			print('The implemented Phi is incorrect')

		if err_un < 1e-6:
			print('The unaltered Phi is correct')
		else:
			print('The unaltered Phi is incorrect')

	# Sum operations from im and un 
	un = un_asn + un_acc + un_art + un_fun
	im = un_art + im_acc + im_asn + im_fun + im_log

	# Return op counts of implemented and unaltered algorithm
	return im, un 
