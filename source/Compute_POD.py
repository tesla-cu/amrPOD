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

# =========================================================================== #
# Function to generate data on the computational cost to compute POD
#
# Inputs:
# - gen_grid    : are we generating the grid 
# - nx          : number of cells in the x direction
# - ny          : number of cells in the y direction
# - nz          : number of cells in the z direction
# - finest      : finest level of AMR
# - l_fracs     : array of the fraction of the domain at a particular level
# - lc_fracs    : array of the fraction of the domain at a particular level that 
#                 is held constant
# - nt          : number of time steps
# - TC_CPU      : are we doing operation counts (TC) or CPU time (CPU)
# - amr_datadir : directory storing AMR data if we aren't generating it
#
# Outputs:
# - R_imp    : either TC or CPU of implemented algorithm to compute R
# - R_unalt  : either TC or CPU of unaltered algorithm to compute R
# - P1_imp   : either TC or CPU of implemented algorithm to compute Phi, M1
# - P1_unalt : either TC or CPU of unaltered algorithm to compute Phi, M1
# - P2_imp   : either TC or CPU of implemented algorithm to compute Phi, M2
# - P2_unalt : either TC or CPU of unaltered algorithm to compute Phi, M2
# - A_imp    : either TC or CPU of implemented algorithm to compute A
# - A_unalt  : either TC or CPU of unaltered algorithm to compute A
# =========================================================================== #
def Compute_POD(gen_grid, nx, ny, nz, finest, l_fracs, lc_fracs, nt, \
	TC_CPU='CPU', amr_datadir=None):
 
	# Helpful quantities derived from user inputs -----------------------------
	nspat = nx*ny*nz   # number of spatial points
	nlev  = finest + 1 # number of AMR levels
	ndim  = 0          # num dimensions
	if nx > 1: ndim += 1 
	if ny > 1: ndim += 1 
	if nz > 1: ndim += 1
	levels = np.arange(0, nlev)
	
	c_l   = np.empty((nlev), dtype=int) # num reps in one dimension
	d_l   = np.empty((nlev), dtype=int) # num reps in all dimensions
	for i in range(nlev):
		c_l[i]    = 2**(finest-i)
		d_l[i]    = (2**ndim)**(finest-i)

	# Define operation weighting for counting ---------------------------------
	wt_art    = 1       # Arithmetic
	wt_acc    = 1       # Memory accessing
	wt_asn    = 1       # Variable assignment
	wt_log    = 1       # Logical test
	wt_fun    = 1       # Function call

	# Load or generate data ---------------------------------------------------
	X      = np.empty((nspat,nt))
	X_grid = np.empty((nspat,nt), dtype=int)

	# Iterate over each snapshot
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

	# Compute grid information from X_grid ------------------------------------
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

	# Initialize dummy matrices
	Psi    = np.ones((nt,nt))
	Lambda = np.ones((nt,nt))
	Phi    = np.ones((nspat,nt))

	# =========================================================================
	# Begin POD operations
	# =========================================================================

	# Compute time complexity of each operation
	if TC_CPU == 'TC':

		R_imp,  R_unalt  = compute_R_TC  (X, X_grid,                 False, \
			d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		P1_imp, P1_unalt = compute_Phi_TC(X, X_grid, Psi, Lambda, 1, False, \
			d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		P2_imp, P2_unalt = compute_Phi_TC(X, X_grid, Psi, Lambda, 2, False, \
			d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun)

		A_imp,  A_unalt  = compute_A_TC  (X, X_grid, Phi,            False, \
			d_l, nt, nspat, finest, wt_art, wt_acc, wt_asn, wt_log, wt_fun)

	# Compute CPU time of each operation
	elif TC_CPU == 'CPU':

		R_imp,  R_unalt  = compute_R_CPU  (X, X_grid,                 False, \
			d_l, nt, nspat, finest)

		P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, False, \
			d_l, nt, nspat, finest)

		P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, False, \
			d_l, nt, nspat, finest)

		A_imp,  A_unalt  = compute_A_CPU  (X, X_grid, Phi,            False, \
			d_l, nt, nspat, finest)

	else:
		print("Input must be either 'CPU' or 'TC'")
		sys.exit()

	return R_imp, R_unalt, P1_imp, P1_unalt, P2_imp, P2_unalt, A_imp, A_unalt

