import numpy as np

 # X_grid - X matrix with all grid levels
 # i      - spatial location
 # idx    - location within lvl=0 block
 # n      - time step
 # clvl   - current lvl within recursion
 # finest - finest lvl 
 # d_l    - array of c 
 # G_mat  - matrix holding elements that contribute to particular lvl within coarse cell
 # nl     - matrix holding number of cells for each lvl within coarse cell

def find_lvl_indices_TC(X_grid, i, idx, n, clvl, finest, d_l, G_mat, nl, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func):
	
	Phi_count_imp_arith   += 1
	Phi_count_imp_logtest += 1
	if clvl < finest-1:
		

		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			Phi_count_imp_arith  += 1
			Phi_count_imp_access += 1

			lvl = X_grid[j,n]
			Phi_count_imp_access += 1
			Phi_count_imp_assign += 1

			Phi_count_imp_logtest += 1
			if lvl == clvl:
				
				nl[clvl, idx] += 1
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1
				Phi_count_imp_arith  += 1

				G_mat[clvl, idx, nl[clvl, idx]-1] = n
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 2
				Phi_count_imp_assign += 1

			else:
				G_mat, nl, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func = \
				find_lvl_indices(X_grid, j, idx, n, clvl+1, finest, d_l, G_mat, nl)
				Phi_count_imp_func   += 1
				Phi_count_imp_assign += 2


			idx += d_l[clvl+1]
			Phi_count_imp_arith  += 2
			Phi_count_imp_access += 1
			Phi_count_imp_assign += 1

	else:
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			Phi_count_imp_arith  += 1
			Phi_count_imp_access += 1

			lvl = X_grid[j,n]
			Phi_count_imp_access += 1
			Phi_count_imp_assign += 1

			Phi_count_imp_logtest += 1
			if lvl == clvl:

				nl[clvl, idx] += 1
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1

				G_mat[clvl, idx, nl[clvl, idx]-1] = n
				Phi_count_imp_access += 2
				Phi_count_imp_arith  += 1
				Phi_count_imp_assign += 1

			else:
				nl[finest, idx] += 1
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1

				G_mat[finest, idx, nl[finest, idx]-1] = n
				Phi_count_imp_access += 2
				Phi_count_imp_arith  += 1
				Phi_count_imp_assign += 1

			idx += 1
			Phi_count_imp_arith  += 1
			Phi_count_imp_assign += 1

	return G_mat, nl, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func


