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

def Compute_POD(gen_grid, nx, ny, nz, finest, l_fracs, lc_fracs, nt, TC_CPU='CPU', amr_datadir=None):
 
	nspat = nx*ny*nz    
	nlev  = finest + 1
	c_l   = np.zeros((nlev), dtype=int)
	d_l   = np.zeros((nlev), dtype=int)

	# ---------- Helpful quantities derived from user inputs --------
	ndim  = 0          # num dimensions
	if nx > 1: ndim += 1 
	if ny > 1: ndim += 1 
	if nz > 1: ndim += 1
	levels = np.arange(0, nlev)
	
	for i in range(nlev):
		c_l[i]    = 2**(finest-i)
		d_l[i]    = (2**ndim)**(finest-i)

	# ---------- Define operation weighting for counting ------------
	wt_art    = 1       # Arithmetic
	wt_acc    = 1       # Memory accessing
	wt_asn    = 1       # Variable assignment
	wt_log    = 1       # Logical test
	wt_fun    = 1       # Function call

	# ---------- Load or generate data ------------------------------
	X      = np.zeros((nspat,nt))
	X_grid = np.zeros((nspat,nt), dtype=int)
	for n in range(nt):

		if gen_grid:
			grid = GenGrid(nx, ny, nz, c_l, d_l, l_fracs, lc_fracs)
			data = grid + 1.5
		else:
			grid = np.fromfile(amr_datadir + 'grid_level%05d.bin' % n).astype(int)
			data = np.fromfile(amr_datadir + 'density%05d.bin' % n)

		# Perform reshaping
		X[:,n]      = Reshape_AMR(nx, ny, nz, finest, data, 'forward')
		X_grid[:,n] = Reshape_AMR(nx, ny, nz, finest, grid, 'forward')

	# ---------- Compute grid information from X_grid ---------------
	l_comp  = np.zeros((nlev)) # computed level fractions
	lc_comp = np.zeros((nlev)) # computed level constant fractions

	# Iterate through all time steps
	for n in range(nt):

		# Find fraction of grid at a particular level. Add fractions
		# for each time step then average
		for l in levels:
			l_comp[l] += np.sum(X_grid[:,n] == l)/nspat

	# Take average
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
	idx_eig     = np.argsort(Lambda) # sort evals and evecs
	Lambda      = Lambda[idx_eig]
	Psi         = Psi[:,idx_eig]
	Phi         = np.matmul(X,Psi)
	Phi         = np.matmul(Phi, np.diag(1/np.sqrt(Lambda)))
	A           = np.matmul(X_tp, Phi)
	Lambda      = np.diag(Lambda) # make this a matrix

	# ---------- Compute time complexity of each operation ----------
	if TC_CPU == 'TC':

		R_imp,  R_unalt  = compute_R_TC(X, X_grid, R, d_l, nt, nspat, \
			wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		P1_imp, P1_unalt = compute_Phi_TC(X, X_grid, Psi, Lambda, 1, Phi, d_l, nt, nspat, finest, \
			wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		P2_imp, P2_unalt = compute_Phi_TC(X, X_grid, Psi, Lambda, 2, Phi, d_l, nt, nspat, finest, \
			wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		A_imp,  A_unalt  = compute_A_TC(X, X_grid, Phi, A, d_l, nt, nspat, finest, \
			wt_art, wt_acc, wt_asn, wt_log, wt_fun)

	# ---------- Compute CPU time of each operation -----------------
	elif TC_CPU == 'CPU':

		R_imp,  R_unalt  = compute_R_CPU(X, X_grid, R, d_l, nt, nspat)
		P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, Phi, d_l, nt, nspat, finest)
		P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, Phi, d_l, nt, nspat, finest)
		A_imp,  A_unalt  = compute_A_CPU(X, X_grid, Phi, A, d_l, nt, nspat, finest)

	else:
		print("Input must be either 'CPU' or 'TC'")
		sys.exit()

	# ---------- Reshape back to original shape ---------------------

	# Iterate through all snapshots
	for n in range(nt):
		Phi[:,n] = Reshape_AMR(nx, ny, nz, finest, Phi[:,n], 'reverse')

	return R_imp, R_unalt, P1_imp, P1_unalt, P2_imp, P2_unalt, A_imp, A_unalt

