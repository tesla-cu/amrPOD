import numpy as np
import time
# from Compute_H import compute_H
# from find_lvl_indices import find_lvl_indices


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

def compute_H(X, Psi, i, idx, n, nt, clvl, finest, d_l, G_mat, nl, H):
    
    if clvl < finest-1:
        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

            if nl[clvl, idx] > 0:

                l_sum = 0
                for m in range(nl[clvl, idx]):
                    k = G_mat[clvl, idx, m]
                    l_sum += X[j,k] * Psi[k,n]

                l_idx = int([j-i]/d_l[clvl])
                for m in range(l_idx*d_l[clvl], (l_idx+1)*d_l[clvl]):
                    H[m, clvl] = l_sum

            nccells = 0
            for l in range(clvl+1):
                nccells += nl[l, idx]

            if nccells < nt:
                H = compute_H(X, Psi, j, idx, n, nt, clvl+1, finest, d_l, G_mat, nl, H)
            idx += d_l[clvl + 1]
            
    else:

        for j in range(i, i+d_l[clvl-1]-1, d_l[clvl]):

            if nl[clvl, idx] > 0:
                l_sum = 0
                for m in range(nl[clvl, idx]):

                    k = G_mat[clvl, idx, m]
                    l_sum += X[j,k] * Psi[k,n]

                for m in range((idx)*d_l[clvl], (idx+1)*d_l[clvl]):
                    H[m, clvl] = l_sum

            if nl[finest, idx] > 0:

                for k in range(j, j+d_l[clvl]):
                    l_sum = 0
                    for m in range(nl[finest, idx]):
                        p = G_mat[finest, idx, m]
                        l_sum += X[k,p] * Psi[p,n]

                    H[k-j+d_l[finest-1] * (idx), finest] = l_sum

            idx += 1

    return H

def compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat, finest):

#--------- Unaltered 

	tic = time.time()
	Phi_unalt = np.zeros((nspat,nt))
	for i in range(nt):
		for j in range(nspat):
			phi_sum = 0
			for k in range(nt):
				phi_sum += X[j,k] * Psi[k,i]
			Phi_unalt[j,i] = phi_sum/np.sqrt(Lambda[i,i])
	time_Phi_unalt = time.time() - tic

	tic = time.time()

#------------------
#--------- Method 1
#------------------ 
	
	if method == 1:

		Phi_imp = np.zeros((nspat, nt))
		G       = np.zeros((nspat), dtype=int)
		i       = 0

		for n in range(nspat):
			if i < nspat:
				X_grid_max = X_grid[i,0]
				for j in range(1,nt):
					if X_grid[i,j] > X_grid_max:
						X_grid_max = X_grid[i,j]
				G[i]  = d_l[X_grid_max]
				i     += G[i]
			else:
				break

#--------- Compute Phi

		for i in range(nt):
			j = 0
			for n in range(nspat):
				if j < nspat:
					phi_sum = 0
					for k in range(nt):
						phi_sum += X[j,k] * Psi[k,i]
					Phi_imp[j:j+G[j], i] = phi_sum / np.sqrt(Lambda[i,i])
					j = j + G[j]

#------------------
#--------- Method 5
#------------------ 

	elif method == 5:

		Phi_imp = np.zeros((nspat, nt))
		finest = len(d_l) - 1

		for i in range(0, nspat, d_l[0]):
			G_mat = np.zeros((finest+1, d_l[1], nt), dtype=int)
			nl    = np.zeros((finest+1, d_l[1]), dtype=int)

			for n in range(nt):

				lvl = X_grid[i,n]
				if lvl == 0:
					nl[0,0] += 1
					G_mat[0,0,nl[0,0]-1] = n

				else:
					if finest > 1:
						G_mat, nl = find_lvl_indices(X_grid, i, 0, n, 1, finest, d_l, G_mat, nl)
					else:
						nl[1,0] += 1
						G_mat[1,0,nl[1,0]-1] = n

	
			for n in range(nt):

				H = np.zeros((d_l[0], finest+1))

				if nl[0,0] > 0:
					l_sum = 0
					for m in range(nl[0,0]):
						k = G_mat[0,0,m]
						l_sum += X[i,k] * Psi[k,n]

					for m in range(d_l[0]):
						H[m,0] = l_sum

				if nl[0,0] < nt:
					if finest > 1:
						H = compute_H(X_grid, Psi, i, 0, n, nt, 1, finest, d_l, G_mat, nl, H)
					else:

						for k in range(i, i+d_l[0]):
							l_sum = 0
							for m in range(nl[1,0]):
								p = G_mat[1, 0, m]
								l_sum += X[k,p] * Psi[p,n]

							H[k-i,1] = l_sum

				for m in range(d_l[0]):
					H_sum = 0
					for l in range(finest+1):
						H_sum += H[m,l]
					Phi_imp[m+i, n] = H_sum/np.sqrt(Lambda[n,n])


	time_Phi_imp = time.time() - tic

	if np.max(abs(np.subtract(Phi_imp, Phi))) < 1e-8:
		print('The implemented Phi is correct')
	else:
		print('The implemented Phi is incorrect')

	if np.max(abs(np.subtract(Phi_unalt, Phi))) < 1e-8:
		print('The unaltered Phi is correct')
	else:
		print('The unaltered Phi is incorrect')

	return time_Phi_unalt, time_Phi_imp











	

