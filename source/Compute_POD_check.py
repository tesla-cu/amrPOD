import numpy as np
import matplotlib.pyplot as plt 
from numpy import linalg as LA

from GenGrid import GenGrid

from R_CPU   import compute_R_CPU
from R_TC    import compute_R_TC

# from Phi_CPU import compute_Phi_CPU
# from Phi_TC  import compute_Phi_TC

from Phi_CPU_iter import compute_Phi_CPU
from Phi_TC_iter  import compute_Phi_TC

from A_CPU   import compute_A_CPU
from A_TC    import compute_A_TC

def Compute_POD_check(nx, ny, nz, finest, nt, amr_datadir):
 
 	# ---------- Helpful quantities derived from user inputs --------
	nspat = nx*ny*nz    
	nlev  = finest + 1
	ndim  = 0            # num dimensions
	if nx > 1: ndim += 1 
	if ny > 1: ndim += 1 
	if nz > 1: ndim += 1
	levels = np.arange(0, nlev)
	
	# ---------- Arrays for repetition ------------------------------
	c_l   = np.zeros((nlev), dtype=int)
	d_l   = np.zeros((nlev), dtype=int)
	for i in range(nlev):
		c_l[i]    = 2**(finest-i)
		d_l[i]    = (2**ndim)**(finest-i)

	# ---------- Load or generate data ------------------------------
	X      = np.zeros((nspat,nt))
	X_grid = np.zeros((nspat,nt), dtype=int)
	for n in range(nt):

		# Load data from file
		grid = np.fromfile(amr_datadir + 'grid_level%05d.bin' % n).astype(int)
		grid = np.reshape(grid, [nx, ny, nz])
		data = np.fromfile(amr_datadir + 'density%05d.bin' % n)
		data = np.reshape(data, [nx, ny, nz])

		# Perform reshaping procedure
		# 1D, no reshaping required
		if ndim == 1:
			grid_1D = np.squeeze(grid)
			data_1D = np.squeeze(data)

		# 2D reshaping procedure, see text for details
		elif ndim == 2:
			grid_1D = np.squeeze(grid)
			data_1D = np.squeeze(data)
			for c in c_l:
				nxr = grid_1D.shape[0]

				grid_1D = np.transpose(grid_1D, ( 1,  0))
				grid_1D = np.reshape(  grid_1D, (-1,  c,  nxr))
				grid_1D = np.transpose(grid_1D, ( 1,  0,  2))
				grid_1D = np.reshape(  grid_1D, ( c, -1))

				data_1D = np.transpose(data_1D, ( 1,  0))
				data_1D = np.reshape(  data_1D, (-1,  c,  nxr))
				data_1D = np.transpose(data_1D, ( 1,  0,  2))
				data_1D = np.reshape(  data_1D, ( c, -1))

		# 3D reshaping procedure, see text for details
		elif ndim == 3:
			grid_1D = grid
			data_1D = data
			for c in c_l:
				nxr = grid_1D.shape[0]
				nyr = grid_1D.shape[1]

				grid_1D = np.transpose(grid_1D, ( 2,  1,  0))
				grid_1D = np.reshape(  grid_1D, (-1,  c,  nyr, nxr))
				grid_1D = np.transpose(grid_1D, ( 1,  0,  2,   3))
				grid_1D = np.reshape(  grid_1D, ( c, -1,  c,   nxr))
				grid_1D = np.transpose(grid_1D, ( 0,  2,  1,   3))
				grid_1D = np.reshape(  grid_1D, ( c,  c, -1))

				data_1D = np.transpose(data_1D, ( 2,  1,  0))
				data_1D = np.reshape(  data_1D, (-1,  c,  nyr, nxr))
				data_1D = np.transpose(data_1D, ( 1,  0,  2,   3))
				data_1D = np.reshape(  data_1D, ( c, -1,  c,   nxr))
				data_1D = np.transpose(data_1D, ( 0,  2,  1,   3))
				data_1D = np.reshape(  data_1D, ( c,  c, -1))

		# Assign new reshaped data to corresponding X matrix
		X[:,n]      = data_1D
		X_grid[:,n] = grid_1D

	# ---------- Compute grid information from X_grid ---------------
	l_comp  = np.zeros((nlev)) # computed level fractions
	lc_comp = np.zeros((nlev)) # computed level constant fractions

	# Compute l_comp
	for n in range(nt):
		for l in levels:
			l_comp[l] += np.sum(X_grid[:,n] == l)/nspat
	l_comp = l_comp/nt
	# print("l_comp = ", l_comp)

	# Compute lc_comp
	for l in levels:
		for i in range(nspat):
			if np.all(X_grid[i,:] == l):
				lc_comp[l] += 1
	lc_comp = lc_comp/nspat
	# print("lc_comp = ", lc_comp)

	# ---------- Calculate POD with matrix operations ---------------
	X_tp        = np.transpose(X)
	R           = np.matmul(X_tp, X)
	Lambda, Psi = LA.eig(R)
	idx_eig     = np.argsort(Lambda) # sort eigenvalues
	Lambda      = Lambda[idx_eig]
	Psi         = Psi[:,idx_eig]
	Phi         = np.matmul(X,Psi)
	Phi         = np.matmul(Phi, np.diag(1/np.sqrt(Lambda)))
	A           = np.matmul(X_tp, Phi)
	Lambda      = np.diag(Lambda) # make this a matrix

	# ---------- Calculate POD with iterative operations ------------
	R_imp,  R_unalt  = compute_R_CPU(X, X_grid, R, d_l, nt, nspat)
	P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, Phi, d_l, nt, nspat, finest)
	P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, Phi, d_l, nt, nspat, finest)
	A_imp,  A_unalt  = compute_A_CPU(X, X_grid, Phi, A, d_l, nt, nspat, finest)

	# ---------- Reshape back to original shape ---------------------

	# Create array that specifies a dimension of reshape
	c_lr       = np.zeros((nlev), dtype=int) # c_l in reverse
	c_lr[0:-1] = c_l[-2::-1] # reversed elements of c_l excluding last
	c_lr[-1]   = nx          # last iteration must be nx

	# Iterate through all snapshots
	for n in range(nt):
		
		phi_1D = Phi[:,n]
		phi_1D = np.expand_dims(phi_1D, axis=0)

		# Perform reshaping procedure
		# 1D, no reshaping required
		if ndim == 1:
			print('1D')
		# 2D reshaping procedure, see text for details
		elif ndim == 2:
			for c in c_lr:
				nxr = phi_1D.shape[0]

				phi_1D = np.reshape(  phi_1D, (nxr, -1, c))
				phi_1D = np.transpose(phi_1D, ( 1,  0,  2))
				phi_1D = np.reshape(  phi_1D, (-1, c))
				phi_1D = np.transpose(phi_1D, ( 1,  0))

		# 3D reshaping procedure, see text for details
		elif ndim == 3:
			for c in c_lr:
				nxr = phi_1D.shape[0]
				nyr = phi_1D.shape[1]

				phi_1D = np.reshape(  phi_1D, ( nxr, nyr,  -1,   c))
				phi_1D = np.transpose(phi_1D, ( 0,  2,  1,   3))
				phi_1D = np.reshape(  phi_1D, ( nxr, -1,  c,   c))
				phi_1D = np.transpose(phi_1D, ( 1,  0,  2,   3))
				phi_1D = np.reshape(  phi_1D, (-1,  c,  c))
				phi_1D = np.transpose(phi_1D, ( 2,  1,  0))

		phi_1D = np.reshape(phi_1D, (nspat))
		Phi[:,n] = phi_1D

	return R, Phi, A