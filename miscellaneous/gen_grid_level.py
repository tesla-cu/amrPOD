import numpy as np

# Directory storing extracted data
datadir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/data/AMR_sim/'

# Set data parameters
dim  = 3   # dimension of extracted data
f    = 4   # finest AMR level
nx   = 512 # nx of extracted data
ny   = 512 # ny of extracted data
nz   = 512 # nz of extracted data
nt   = 2   # number of time steps
nlev = f+1 # number of AMR levels

# Compute c_l
c_l = np.empty((nlev), dtype=int)
for i in range(nlev):
	c_l[i] = 2**(f-i)

# Write out grid level data for each snapshot
for n in range(1,nt):
	print('writing out grid level data for snapshot %i' % n, flush=True)

	# Load z velocity data
	data = np.fromfile(datadir + 'z_velocity%05d.bin' % n)
	data = np.reshape(data, (nx,ny,nz))

	grid_level    = np.full((nx,ny,nz), f, dtype=np.float64)

	# Iterate reverse through all c_l except finest where c_f=1
	for l_inv, c in enumerate(c_l[-2::-1]):

		# Current level we are checking
		l = f - l_inv - 1

		# Make temp slightly smaller grid level mat to store diffs in cells
		grid_level_tmp = grid_level[:,:,0:nz-c+1]

		# Check if cells c-1 are identical
		diff = data[:,:,c-1:nz] - data[:,:,0:nz-c+1]

		# Find cells that are identical
		grid_level_tmp[diff == 0.0] = l

		# Assign correct cells the grid level
		for i in range(c):
			grid_level[:,:,i::c] = grid_level_tmp[:,:,0::c]

	grid_level.tofile(datadir + 'grid_level%05d.bin' % n)






