# ================================================= #
# Code:        Reshape for AMR data                 #
# Authors:     Michael Meehan and Sam Simons-Wellin #
# Institution: University of Colorado Boulder       #
# Year:        2019                                 #
# ================================================= #

import numpy as np
import sys

# ================================================================= #
# Function to reshape AMR data that forces repeated cells to stay
# continguous when reshaped into a 1D array.
#
# Inputs:
# - nx        : number of cells in x direction
# - ny        : number of cells in y direction
# - nz        : number of cells in z direction
# - finest    : finest AMR grid level
# - data      : 1D array of data simply reshaped
# - direction : either 'forward' or 'reverse'
#               - forward - put simply reshaped into strategically
#                           reshaped
#               - reverse - put strategically reshaped into simply
#                           reshaped
#
# Outputs:
# - data : 1D array with repeated cells now contiguous
# ================================================================= #
def Reshape_AMR(nx, ny, nz, finest, data, direction):

	# ---------- Helpful quantities derived from user inputs --------
	nspat = nx*ny*nz    
	nlev  = finest + 1
	ndim  = 0            # num dimensions
	if nx > 1: ndim += 1 
	if ny > 1: ndim += 1 
	if nz > 1: ndim += 1

	# Array for repetition 
	c_l   = np.zeros((nlev), dtype=int)
	for i in range(nlev):
		c_l[i]    = 2**(finest-i)

	# ========== Forward direction reshaping ========================
	if direction[0].lower() == 'f':

		# Reshape data back into original form 
		data = np.reshape(data, [nx, ny, nz])

		# ---------- Perform reshaping procedure --------------------

		# 1D, no reshaping required
		if ndim == 1:
			data_1D = np.squeeze(data)

		# 2D reshaping procedure, see text for details
		elif ndim == 2:
			data_1D = np.squeeze(data)
			for c in c_l:
				nxr = data_1D.shape[0]

				data_1D = np.transpose(data_1D, ( 1,  0))
				data_1D = np.reshape(  data_1D, (-1,  c,  nxr))
				data_1D = np.transpose(data_1D, ( 1,  0,  2))
				data_1D = np.reshape(  data_1D, ( c, -1))

		# 3D reshaping procedure, see text for details
		elif ndim == 3:
			data_1D = data
			for c in c_l:
				nxr = data_1D.shape[0]
				nyr = data_1D.shape[1]

				data_1D = np.transpose(data_1D, ( 2,  1,  0))
				data_1D = np.reshape(  data_1D, (-1,  c,  nyr, nxr))
				data_1D = np.transpose(data_1D, ( 1,  0,  2,   3))
				data_1D = np.reshape(  data_1D, ( c, -1,  c,   nxr))
				data_1D = np.transpose(data_1D, ( 0,  2,  1,   3))
				data_1D = np.reshape(  data_1D, ( c,  c, -1))


	# ========== Reverse direction reshaping ========================
	elif direction[0].lower() == 'r':

		# Create array that specifies a dimension of reshape
		c_lr       = np.zeros((nlev), dtype=int) # c_l in reverse
		c_lr[0:-1] = c_l[-2::-1] # reversed elements of c_l excluding last
		c_lr[-1]   = nx          # last iteration must be nx

		# 1D, no reshaping required
		if ndim == 1:
			data_1D = data

		# 2D reshaping procedure, see text for details
		elif ndim == 2:
			data_2D = np.expand_dims(data, axis=0)
			for c in c_lr:
				nxr = data_2D.shape[0]

				data_2D = np.reshape(  data_2D, (nxr, -1, c))
				data_2D = np.transpose(data_2D, ( 1,  0,  2))
				data_2D = np.reshape(  data_2D, (-1, c))
				data_2D = np.transpose(data_2D, ( 1,  0))

			data_1D = np.reshape(data_2D, (nspat))

		# 3D reshaping procedure, see text for details
		elif ndim == 3:
			data_3D = np.expand_dims(data, axis=0)
			for c in c_lr:
				nxr = data_3D.shape[0]
				nyr = data_3D.shape[1]

				data_3D = np.reshape(  data_3D, ( nxr, nyr,  -1,   c))
				data_3D = np.transpose(data_3D, ( 0,  2,  1,   3))
				data_3D = np.reshape(  data_3D, ( nxr, -1,  c,   c))
				data_3D = np.transpose(data_3D, ( 1,  0,  2,   3))
				data_3D = np.reshape(  data_3D, (-1,  c,  c))
				data_3D = np.transpose(data_3D, ( 2,  1,  0))

			data_1D = np.reshape(data_3D, (nspat))

	else:
		print("Reshaping direction must be either 'forward' or 'reverse'!")
		sys.exit()

	# ---------- Return newly reshaped data -------------------------
	return np.squeeze(data_1D)














