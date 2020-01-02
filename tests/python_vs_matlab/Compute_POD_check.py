import numpy as np
import matplotlib.pyplot as plt 
from numpy import linalg as LA

from GenGrid import GenGrid
from Reshape_AMR import Reshape_AMR

from R_CPU   import compute_R_CPU
from R_TC    import compute_R_TC

from Phi_CPU import compute_Phi_CPU
from Phi_TC  import compute_Phi_TC

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

	# ---------- Load data ------------------------------------------
	X      = np.zeros((nspat,nt))
	X_grid = np.zeros((nspat,nt), dtype=int)
	for n in range(nt):

		# Load data from file
		grid = np.fromfile(amr_datadir + 'grid_level%05d.bin' % n).astype(int)
		data = np.fromfile(amr_datadir + 'density%05d.bin' % n)

		# Perform reshaping
		# X[:,n]      = Reshape_AMR(nx, ny, nz, finest, data, 'forward')
		# X_grid[:,n] = Reshape_AMR(nx, ny, nz, finest, grid, 'forward')

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
	R_imp,  R_unalt  = compute_R_CPU  (X, X_grid,                 R,   d_l, nt, nspat, finest)
	P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, Phi, d_l, nt, nspat, finest)
	P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, Phi, d_l, nt, nspat, finest)
	A_imp,  A_unalt  = compute_A_CPU  (X, X_grid, Phi,            A,   d_l, nt, nspat, finest)

	# ---------- Reshape back to original shape ---------------------

	# Iterate through all snapshots
	for n in range(nt):
		Phi[:,n] = Reshape_AMR(nx, ny, nz, finest, Phi[:,n], 'reverse')

	# ---------- Return data ----------------------------------------
	return R, Phi, A
