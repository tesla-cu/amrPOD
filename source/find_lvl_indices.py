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

def find_lvl_indices(X_grid, i, idx, n, clvl, finest, d_l, G_mat, nl):

	print(i)
	print(X_grid[i,:])

	if clvl < finest-1:
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			#print(j)
			lvl = X_grid[j,n]

			if lvl == clvl:
				nl[clvl, idx] += 1
				G_mat[clvl, idx, nl[clvl, idx]-1] = n

			else:
				G_mat, nl = find_lvl_indices(X_grid, j, idx, n, clvl+1, finest, d_l, G_mat, nl)
			idx += d_l[clvl+1]

	else:
		for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):
			lvl = X_grid[j,n]
			if lvl == clvl:
				# print(nl)
				nl[clvl, idx] += 1
				G_mat[clvl, idx, nl[clvl, idx]-1] = n
			else:
				nl[finest, idx] += 1
				G_mat[finest, idx, nl[finest, idx]-1] = n

			idx += 1

	return G_mat, nl