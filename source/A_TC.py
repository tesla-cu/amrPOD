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
# - wt_art : weighting of arithmetic operations
# - wt_acc : weighting of accessing operations
# - wt_asn : weighting of assignment operations
# - wt_log : weighting of logical operations
# - wt_fun : weighting of function calls
#
# Outputs:
# - im : num of operations to compute A using implemented algorithm
# - un : num of operations to compute A using unaltered algorithm
# ================================================================= #
def compute_A_TC(X, X_grid, Phi, A, d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun):

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

	A_un = np.zeros((nt, nt))
	un_asn += wt_asn
	un_fun += wt_fun

	for m in range(nt):
		un_art += wt_art
		un_asn += wt_asn

		for n in range(nt):
			un_art += wt_art
			un_asn += wt_asn

			a_sum  = 0
			un_asn += wt_asn

			for i in range(nspat):
				un_art += wt_art
				un_asn += wt_asn

				a_sum  += X[i,m] * Phi[i,n]
				un_acc += 2*wt_acc
				un_art += 2*wt_art
				un_asn += wt_asn

			A_un[m,n] = a_sum
			un_acc    += wt_acc
			un_asn    += wt_asn

	# ========== Implemented Computation ========================== #

	A_im   =  np.zeros((nt, nt))
	im_asn += wt_asn
	im_fun += wt_fun

	G      =  np.zeros((nspat), dtype=int)
	im_asn += wt_asn
	im_fun += wt_fun

	i       = 0
	im_asn += wt_asn

	d_f1 = d_l[finest-1]
	im_asn += wt_asn
	im_art += wt_art
	im_acc += im_acc
	
	for ii in range(nspat):
		im_art += wt_art
		im_asn += wt_asn

		im_log += wt_log
		if i < nspat:

			X_grid_max = X_grid[i,0]
			im_acc     += wt_acc
			im_asn     += wt_asn

			for m in range(1,nt):
				im_art += wt_art
				im_asn += wt_asn

				im_log += wt_log
				im_acc += wt_acc
				if X_grid[i,m] > X_grid_max:

					X_grid_max = X_grid[i,m]
					im_acc     += wt_acc
					im_asn     += wt_asn

					im_log += wt_log
					if X_grid_max == finest:
						break

			im_log += wt_log
			if X_grid_max == finest:
				G[i] = 1
				im_acc += wt_acc
				im_asn += wt_asn

				i    += d_f1
				im_art += wt_art
				im_asn += wt_asn
			else:
				G[i]   =  d_l[X_grid_max]
				im_acc += 2*wt_acc
				im_asn += wt_asn

				i      += G[i]
				im_art += wt_art
				im_acc += wt_acc
				im_asn += wt_asn

		else:
			break

	for m in range(nt):
		im_art += wt_art
		im_asn += wt_asn

		for n in range(nt): 
			im_art += wt_art
			im_asn += wt_asn

			a_sum  = 0
			im_asn += wt_asn

			i      = 0
			im_asn += wt_asn

			for ii in range(nspat):
				im_art += wt_art
				im_asn += wt_asn

				im_log += wt_log
				if i < nspat:

					im_log += wt_log
					if G[i] == 1:

						im_art += wt_art
						for j in range(i,i+d_f1):
							im_art += wt_art
							im_asn += wt_asn

							a_sum += X[j,m]*Phi[j,n]
							im_acc += 2*wt_acc
							im_art += 2*wt_art
							im_asn += wt_asn

						i += d_f1
					else:
						a_sum  += G[i]*X[i,m]*Phi[i,n]
						im_acc += 3*wt_acc
						im_art += 3*wt_art
						im_asn += wt_asn

						i      += G[i]
						im_acc += wt_acc
						im_asn += wt_asn
						im_art += wt_art

				else:
					break

			A_im[m,n] = a_sum
			im_acc    += wt_acc
			im_asn    += wt_asn

	# ========== Check Correctness of Matrices ==================== #

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

	# ========== Sum operations from im and un =================== #

	un = un_asn + un_acc + un_art + un_fun
	im = im_asn + im_acc + im_art + im_fun + im_log

	# Return op counts of implemented and unaltered algorithm
	return im, un

	