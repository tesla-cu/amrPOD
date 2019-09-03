import numpy as np
import time

from find_lvl_indices_TC import find_lvl_indices_TC
from Compute_H_TC import compute_H_TC

def compute_Phi_TC(X_grid, method, d_l, nt, nspat, finest):

	Phi_count_unalt        = 0
	Phi_count_unalt_arith  = 0
	Phi_count_unalt_access = 0
	Phi_count_unalt_assign = 0

	Phi_count_imp          = 0
	Phi_count_imp_arith    = 0
	Phi_count_imp_access   = 0
	Phi_count_imp_assign   = 0
	Phi_count_imp_logtest  = 0
	Phi_count_imp_func     = 0		

	#--------- Unaltered 

	#Phi_unalt = np.zeros((nspat,nt))
	Phi_count_unalt_assign += 1

	for i in range(nt):
		Phi_count_unalt_arith += 1
		Phi_count_unalt_access  += 1

		for j in range(nspat):
			Phi_count_unalt_arith += 1
			Phi_count_unalt_access  += 1			

			#phi_sum = 0
			Phi_count_unalt_assign += 1

			for k in range(nt):
				Phi_count_unalt_arith += 1
				Phi_count_unalt_access  += 1

				#phi_sum += X[j,k] * Psi[k,i]
				Phi_count_unalt_access += 2
				Phi_count_unalt_arith  += 2
				Phi_count_unalt_assign += 1

			#Phi_unalt[j,i] = phi_sum/np.sqrt(Lambda[i,i])
			Phi_count_unalt_access += 2
			Phi_count_unalt_arith  += 2
			Phi_count_unalt_assign += 1

#------------------
#--------- Method 1
#------------------ 
	
	if method == 1:

		#Phi_imp = np.zeros((nspat, nt))
		Phi_count_imp_assign += 1

		G       = np.zeros((nspat), dtype=int)
		Phi_count_imp_assign += 1

		i       = 0
		Phi_count_imp_assign += 1

		for n in range(nspat):
			Phi_count_imp_arith  += 1
			Phi_count_imp_access += 1

			Phi_count_imp_logtest += 1
			if i < nspat:

				X_grid_max = X_grid[i,0]
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1

				Phi_count_imp_logtest += 1
				if X_grid_max < finest:

					for j in range(1,nt):
						Phi_count_imp_arith  += 1
						Phi_count_imp_access += 1

						Phi_count_imp_logtest += 1
						Phi_count_imp_access  += 1
						if X_grid[i,j] > X_grid_max:

							X_grid_max = X_grid[i,j]
							Phi_count_imp_access += 1
							Phi_count_imp_assign += 1

							Phi_count_imp_logtest += 1
							if X_grid_max == finest:

								break

				G[i]  = d_l[X_grid_max]
				Phi_count_imp_access += 2
				Phi_count_imp_assign += 1

				i     += G[i]
				Phi_count_imp_access += 1
				Phi_count_imp_arith  += 1
				Phi_count_imp_assign += 1

			else:
				break

#--------- Compute Phi

		for i in range(nt):
			Phi_count_imp_arith  += 1
			Phi_count_imp_access += 1

			j = 0
			Phi_count_imp_assign += 1

			for n in range(nspat):
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 1

				Phi_count_imp_logtest += 1 
				if j < nspat:

					phi_sum = 0
					Phi_count_imp_assign += 1

					for k in range(nt):
						Phi_count_imp_arith  += 1
						Phi_count_imp_access += 1

						#phi_sum += X[j,k] * Psi[k,i]
						Phi_count_imp_access += 2
						Phi_count_imp_arith  += 2
						Phi_count_imp_assign += 1

					#Phi_imp[j:j+G[j], i] = phi_sum / np.sqrt(Lambda[i,i])
					Phi_count_imp_access += 3
					Phi_count_imp_arith  += 3
					Phi_count_imp_assign += 1

					j += G[j]
					Phi_count_imp_arith  += 1
					Phi_count_imp_access += 1
					Phi_count_imp_assign += 1

#------------------
#--------- Method 5
#------------------ 
	
	elif method == 5:

		#Phi_imp = np.zeros((nspat, nt))
		Phi_count_imp_assign += 1

		#rat = d_l[0]/d_l[1]

		finest = len(d_l) - 1
		Phi_count_imp_access += 1
		Phi_count_imp_arith  += 1
		Phi_count_imp_assign += 1

		for i in range(0, nspat, d_l[0]):
			Phi_count_imp_arith  += 1
			Phi_count_imp_access += 1

			G_mat = np.zeros((finest+1, d_l[1], nt), dtype=int)
			Phi_count_imp_access += 1
			Phi_count_imp_assign += 1

			nl    = np.zeros((finest+1, d_l[1]), dtype=int)
			Phi_count_imp_access += 1
			Phi_count_imp_assign += 1

			for n in range(nt):
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 1

				lvl = X_grid[i,n]
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1

				Phi_count_imp_logtest += 1
				if lvl == 0:

					nl[0,0] += 1
					Phi_count_imp_access += 1
					Phi_count_imp_arith  += 1
					Phi_count_imp_assign += 1

					G_mat[0,0,nl[0,0]-1] = n
					Phi_count_imp_access += 2
					Phi_count_imp_arith  += 1
					Phi_count_imp_assign += 1

				else:

					Phi_count_imp_logtest += 1
					if finest > 1:

						G_mat, nl, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func = \
						find_lvl_indices_TC(X_grid, i, 0, n, 1, finest, d_l, G_mat, nl, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func)
						Phi_count_imp_func   += 1
						Phi_count_imp_assign += 2

					else:
						nl[1,0] += 1
						Phi_count_imp_arith    += 1
						Phi_count_imp_access += 1
						Phi_count_imp_assign   += 1

						G_mat[1,0,nl[1,0]-1] = n
						Phi_count_imp_access   += 2
						Phi_count_imp_arith    += 1
						Phi_count_imp_assign   += 1
	
			for n in range(nt):
				Phi_count_imp_arith  += 1
				Phi_count_imp_access += 1

				H = np.zeros((d_l[0], finest+1))
				Phi_count_imp_access += 1
				Phi_count_imp_assign += 1

				Phi_count_imp_access  += 1
				Phi_count_imp_logtest += 1
				if nl[0,0] > 0:

					l_sum = 0
					Phi_count_imp_assign += 1

					for m in range(nl[0,0]):
						Phi_count_imp_arith  += 1
						Phi_count_imp_access += 1

						k = G_mat[0,0,m]
						Phi_count_imp_access += 1
						Phi_count_imp_assign   += 1

						#l_sum += X[i,k] * Psi[k,n]
						Phi_count_imp_access += 2
						Phi_count_imp_arith  += 2
						Phi_count_imp_assign += 1

					for m in range(d_l[0]):
						Phi_count_imp_arith  += 1
						Phi_count_imp_access += 1

						#H[m,0] = l_sum
						Phi_count_imp_access += 1
						Phi_count_imp_assign += 1

				Phi_count_imp_logtest += 1
				Phi_count_imp_access  += 1
				if nl[0,0] < nt:

					Phi_count_imp_logtest += 1
					if finest > 1:

						H, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func = \
						compute_H_TC(X_grid, i, 0, n, nt, 1, finest, d_l, G_mat, nl, H, Phi_count_imp_arith, Phi_count_imp_access, Phi_count_imp_assign, Phi_count_imp_logtest, Phi_count_imp_func)

						Phi_count_imp_assign += 1
						Phi_count_imp_func   += 1

					else:

						for k in range(i, i+d_l[0]):
							Phi_count_imp_arith  += 1
							Phi_count_imp_access += 1

							# l_sum = 0
							Phi_count_imp_assign += 1

							for m in range(nl[1,0]):
								Phi_count_imp_arith  += 1
								Phi_count_imp_access += 1

								p = G_mat[1, 0, m]
								Phi_count_imp_access += 1
								Phi_count_imp_assign += 1

								#l_sum += X[k,p] * Psi[p,n]
								Phi_count_imp_arith  += 2
								Phi_count_imp_access += 2
								Phi_count_imp_assign += 1

							#H[k-i,1] = l_sum
							Phi_count_imp_access += 1
							Phi_count_imp_arith  += 1
							Phi_count_imp_assign += 1

				for m in range(d_l[0]):
					Phi_count_imp_arith  += 1
					Phi_count_imp_access += 1

					#H_sum = 0
					Phi_count_imp_assign += 1

					for l in range(finest+1):
						Phi_count_imp_arith  += 1
						Phi_count_imp_access += 1

						#H_sum += H[m,l]
						Phi_count_imp_access += 1
						Phi_count_imp_arith  += 1
						Phi_count_imp_assign += 1

					#Phi_imp[m+i, n] = H_sum/np.sqrt(Lambda[n,n])
					Phi_count_imp_access += 2
					Phi_count_imp_arith  += 3
					Phi_count_imp_assign += 1



	Phi_count_unalt = Phi_count_unalt_access + Phi_count_unalt_assign + Phi_count_unalt_arith
	Phi_count_imp   = Phi_count_imp_access + Phi_count_imp_assign + Phi_count_imp_arith + Phi_count_imp_logtest + Phi_count_imp_func

	return Phi_count_imp, Phi_count_unalt 