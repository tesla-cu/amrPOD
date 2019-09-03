import numpy as np
import time
from Compute_H import compute_H
from find_lvl_indices import find_lvl_indices

def compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat):

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

# #------------------
# #--------- Method 4
# #------------------ 

# 	elif method == 4:

# 		idx_arr = np.arange(0, nt)
# 		Phi_imp = np.zeros((nspat, nt))
# 		rat     = int(d_l[0]/d_l[1])

# 		for i in range(0, nspat, d_l[0]):
# 			G0_mat = np.zeros((nt), dtype=int)
# 			G1_mat = np.zeros((rat, nt), dtype=int)
# 			G2_mat = np.zeros((rat*rat, nt), dtype=int)
# 			G3_mat = np.zeros((rat*rat, nt), dtype=int)

# 			nl0    = -1
# 			nl1	   = -1*np.ones((rat), dtype=int)
# 			nl2    = -1*np.ones((rat*rat), dtype=int)
# 			nl3    = -1*np.ones((rat*rat), dtype=int)

# 			for n in range(nt):
# 				il1 = 0
# 				il2 = 0
# 				lvl = X_grid[i,n]

# 				if lvl == 0:
# 					nl0 += 1
# 					G0_mat[nl0] = n
# 				else:
# 					for j in range(i, i+d_l[0]-1, d_l[1]): # Check -1
# 						lvl = X_grid[j,n]

# 						if lvl == 1:
# 							nl1[il1] += 1
# 							G1_mat[il1, nl1[il1]] = n

# 							il1 += 1
# 							il2 += d_l[2]

# 						else:
# 							for k in range(j, j+d_l[1]-1, d_l[2]):
# 								lvl = X_grid[k,n]

# 								if lvl == 2:
# 									nl2[il2] += 1
# 									G2_mat[il2, nl2[il2]] = n

# 								else:
# 									nl3[il2] += 1
# 									#print(il2)
# 									#print('\n', nl3.shape)
# 									print(nl3[il2])
# 									print('\n', G3_mat.shape)
# 									G3_mat[il2, nl3[il2]] = n

# 							il2 += 1

# 			for m in range(nt):
# 				il1 = 0
# 				il2 = 0
# 				H   = np.zeros((d_l[0], 4))

# 				if nl0 > 0:
# 					l0_sum = 0
# 					for n in range(nl0):
# 						idx = G0_mat[n]
# 						l0_sum += X[i,idx] * Psi[idx,m]
# 					for n in range(d_l[0]):
# 						H[n,0] = l0_sum

# 				if nl0 < nt:
# 					for j in range(i, i+d_l[0]-1, d_l[1]):
# 						if nl1[il1] > 0:
# 							l1_sum = 0
# 							for n in range(nl1[il1]):
# 								idx = G1_mat[il1,n]
# 								l1_sum += X[j,idx] * Psi[idx,m]
# 							for n in range(j-i+1, j-i+d_l[1]):
# 								H[n,1] = l1_sum

# 						if nl0 + nl1[il1] < nt:
# 							for k in range(j, j+d_l[1]-1, d_l[2]):
# 								if nl2[il2] > 0:
# 									l2_sum = 0
# 									for n in range(nl2[il2]):
# 										idx = G2_mat[il2,n]
# 										l2_sum += X[k,idx] * Psi[idx,m]
# 									for n in range(k-i+1, k-i+d_l[2]):
# 										H[n,2] = l2_sum

# 								if nl3[il2] > 0:
# 									for l in range(k, k+d_l[2]-1):
# 										l3_sum = 0
# 										for n in range(nl3[il2]):
# 											idx = G3_mat[il2, n]
# 											l3_sum += X[l,idx] * Psi[idx,m]
# 										H[l-i+1,4] = l3_sum

# 								il2 += 1
# 						else:
# 							il2 += d_l[2]

# 						il1 += 1

# 				for n in range(d_l[0]):
# 					H_sum = 0
# 					for l in range(4):
# 						H_sum += H[n,l]
# 					Phi_imp[n+i-1,m] = H_sum/np.sqrt(Lambda[m,m])

#------------------
#--------- Method 5
#------------------ 

	elif method == 5:

		Phi_imp = np.zeros((nspat, nt))
		rat = d_l[0]/d_l[1]
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

	if np.max(abs(np.subtract(Phi_imp, Phi))) < 1e-2:
		print('The implemented Phi is correct')
	else:
		print('The implemented Phi is incorrect')

	if np.max(abs(np.subtract(Phi_unalt, Phi))) < 1e-2:
		print('The unaltered Phi is correct')
	else:
		print('The unaltered Phi is incorrect')

	return time_Phi_unalt, time_Phi_imp











	

