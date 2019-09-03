import numpy as np
import time

def compute_R_TC(X_grid, R, d_l, nt, nspat):


	R_count_unalt        = 0
	R_count_unalt_arith  = 0
	R_count_unalt_access = 0
	R_count_unalt_assign = 0

	R_count_imp          = 0
	R_count_imp_arith    = 0
	R_count_imp_access   = 0
	R_count_imp_assign   = 0
	R_count_imp_logtest  = 0	

#--------- Unaltered 

	#R_unalt = np.zeros(nt, nt)
	R_count_unalt_assign += 1

	for i in range(nt):
		R_count_unalt_arith  += 1
		R_count_unalt_access += 1

		for j in range(i+1):
			R_count_unalt_arith  += 1
			R_count_unalt_access += 1

			#r_sum = 0
			R_count_unalt_assign += 1

			for k in range(nspat):
				R_count_unalt_arith  += 1
				R_count_unalt_access += 1

				#r_sum += X[k,i] * X[k,j]
				R_count_unalt_arith  += 2
				R_count_unalt_access += 2
				R_count_unalt_assign += 1

			#R_unalt[i,j] = r_sum
			R_count_unalt_access += 1
			R_count_unalt_assign += 1			

			#R_unalt[j,i] = r_sum
			R_count_unalt_access += 1
			R_count_unalt_assign += 1	

#--------- Implemented Computation

	#R_imp = np.zeros(nt, nt)
	R_count_imp_assign += 1

	for i in range(nt):
		R_count_imp_arith  += 1
		R_count_imp_access += 1

		for j in range(i+1):
			R_count_imp_arith  += 1
			R_count_imp_access += 1

			R_count_imp_logtest += 1
			if i == j:

				#r_sum = 0
				R_count_imp_assign += 1

				k = 0
				R_count_imp_assign += 1

				for n in range(nspat):
					R_count_imp_arith  += 1
					R_count_imp_access += 1

					R_count_imp_logtest += 1
					if k < nspat:

						c_val = d_l[X_grid[k,i]]
						R_count_imp_access += 2
						R_count_imp_assign += 1

						#r_sum += c_val * X[k,j] * X[k,i]
						R_count_imp_access += 2
						R_count_imp_arith += 3
						R_count_imp_assign += 1

						k += c_val
						R_count_imp_assign += 1
						R_count_imp_arith += 1

					else:
						break
				#R_imp[i,j] = r_sum
				R_count_imp_assign += 1

				#R_imp[j,i] = r_sum
				R_count_imp_assign += 1

			else:
				#r_sum = 0
				R_count_imp_assign += 1

				k = 0
				R_count_imp_assign += 1

				for n in range(nspat):
					R_count_imp_arith += 1
					R_count_imp_access += 1

					R_count_imp_logtest += 1
					if k < nspat:

						R_count_imp_access += 2
						R_count_imp_logtest += 1
						if X_grid[k,i] > X_grid[k,j]:

							c_min = d_l[X_grid[k,i]]
							R_count_imp_access += 2 
							R_count_imp_assign += 1

						else:
							c_min = d_l[X_grid[k,j]]
							R_count_imp_access += 2
							R_count_imp_assign += 1

						#r_sum += c_min * X_grid[k, j] * X[k,i]
						R_count_imp_arith += 3
						R_count_imp_access +=2
						R_count_imp_assign += 1

						k += c_min
						R_count_imp_arith += 1
						R_count_imp_assign += 1

					else:
						break
				#R_imp[i,j] = r_sum
				R_count_imp_assign += 1
				R_count_imp_access += 1

				#R_imp[j,i] = r_sum
				R_count_imp_assign += 1
				R_count_imp_access += 1


	R_count_unalt = R_count_unalt_access + R_count_unalt_assign + R_count_unalt_arith
	R_count_imp   = R_count_imp_access + R_count_imp_assign + R_count_imp_arith + R_count_imp_logtest

	return R_count_imp, R_count_unalt
