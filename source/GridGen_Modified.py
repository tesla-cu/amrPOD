import numpy as np
import matplotlib.pyplot as plt 
import random as rand

import numpy as np
import matplotlib.pyplot as plt 
import random as rand

def GridGen_Modified(nx, ny, nz, nt, c_l, d_l, l_arr, lc_arr):

        # ---------- Arguments ------------------------------------------
        # nx     - num grid cells in x-dir at finest resolution
        # ny     - num grid cells in y-dir at finest resolution
        # nz     - num grid cells in z-dir at finest resolution
        # nt     - num time steps
        # c_l    - 
        # d_l    - 
        # l_arr  - 
        # lc_arr - 

        # Note, set lc_arr[0] in argument

        # ---------- Error checking 

        if sum(l_arr) != 1:
                print('Error: l_arr does not sum to 1.0!')
                exit()

        # error for l<lc and other things

        # # ---------- Helpful quantities derived from user inputs
        finest = len(l_arr)-1    # finest lvl of AMR
        nspat = int(nx*ny*nz)   # num spatial points
        nlev  = finest + 1 # num levels
        ndim  = 0          # num dimensions
        if nx > 1: ndim += 1 
        if ny > 1: ndim += 1 
        if nz > 1: ndim += 1 

        # ---------- Initialization of arrays and matrices
        X_grid = np.zeros((nspat, nt), dtype=int) # matrix for POD operations
        nccells_lc_arr  = np.zeros((nlev), dtype=int)
        levels = np.arange(0, nlev)

        # ---------- Compute quantities used to generate grids  

        # Compute the number of coarse cells we can tag for lc's
        nc_cells  = nspat / d_l[0]

        # Compute the number of coarse cells we need to tag for the lc
        # Note: This is a cumulative sum because once we are done tagging the
        #       lc of a coarser level, we then move to taking the current lc 
        for l in levels:
                nccells_lc_arr[l] = sum(nc_cells*lc_arr[0:l+1])

        # ========== Begin operation to compute grids =======================

        # Iterate in time computing new grids with each snapshot
        for n in range(nt):
                # Tag cells for lc
                count = 0
                grid = np.zeros((nx, ny, nz))
                for i in range(0, nx, c_l[0]):
                        for j in range(0, ny, c_l[0]):
                                for k in range(0, nz, c_l[0]):
                                        count += 1
                                        for l in levels:
                                                if count <= nccells_lc_arr[l]:
                                                        grid[i:i+c_l[0], j:j+c_l[0], k:k+c_l[0]] = l
                                                        break

                # Find the indices that we can start at without disrupting lc
                count = 0
                flag = False
                for i in range(0, nx, c_l[0]):
                        for j in range(0, ny, c_l[0]):
                                for k in range(0, nz, c_l[0]):
                                        if count < nccells_lc_arr[-1]:
                                                count += 1
                                        else:
                                                flag = True
                                                break
                                if flag:
                                        break
                        if flag:
                                break
                start_i = i
                start_j = j
                start_k = k

                # Find the number of valid cells we can tag and number of cells
                # we need to tag.
                # Note: We say we can only tag cells of one level coarser, so 
                #       only a l=2 can be tagged if it is a l=1 cell, etc. So,the
                #       number of cells valid for l=2 is the number of levels 
                #       tagged for l=1. A level that has finer grid cells must
                #       account for the additional cells at those finer levels
                nc_valid    = np.zeros(finest)
                nc_valid[0] = nspat * (1 - sum(lc_arr)) # must not dirupt lc's
                nc_tag      = np.zeros(finest)
                nc_tag[0]   = nspat * (sum(l_arr[1:nlev] - lc_arr[1:nlev]))
                for i in range(1,finest):
                        nc_valid[i]   = nc_tag[i-1]
                        nc_tag[i]  = nspat * sum(l_arr[i+1:nlev] - lc_arr[i+1:nlev])

                # Strategy: iterate through levels up to finest-1 starting at 
                # indices after lc tagged. Define a probability based on number
                # of valid cells and number of cells we need to tag. Draw a 
                # random number, and based on that relative to ratio of valid
                # to tagged, decide if we tag and split that cell
                for l in range(finest):

                        # Iterate through x-dimension
                        for i in range(start_i, nx, c_l[l]):
                                if i >= start_i and i < start_i + c_l[0]:
                                        j = start_j
                                else:
                                        j = 0

                                # Iterate through y-dimension
                                while j < ny:

                                        if i >= start_i and i < start_i + c_l[0] and \
                                           j >= start_j and j < start_j + c_l[0]:
                                                k = start_k
                                        else:
                                                k = 0

                                        # Iterate through z-dimension
                                        while k < nz:
                                                rval = rand.random()

                                                # Ensure probability is a valid value (not inf/nan)
                                                if nc_valid[l] == 0 and nc_tag[l] != 0:
                                                        prob = 1
                                                elif nc_valid[l] == 0 and nc_tag[l] == 0:
                                                        prob = 0 
                                                else:
                                                        prob = (nc_tag[l]/nc_valid[l]) 

                                                if l == 0:
                                                        if rval < prob:   
                                                                grid[i:i+c_l[l], j:j+c_l[l], k:k+c_l[l]] = l+1
                                                                nc_tag[l] -= d_l[l]
                                                        nc_valid[l] -= d_l[l]
                                                        
                                                else:
                                                        if grid[i,j,k] > l-1:
                                                                if rval < prob:
                                                                        grid[i:i+c_l[l], j:j+c_l[l], k:k+c_l[l]] = l+1
                                                                        nc_tag[l] -= d_l[l]
                                                                nc_valid[l] -= d_l[l]

                                                k += c_l[l]
                                        j += c_l[l]


                # 1D, no reshaping required
                if ndim == 1:
                        grid_1D = np.squeeze(grid)

                # 2D reshaping procedure, see text for details
                elif ndim == 2:
                        grid_1D = np.squeeze(grid)
                        for c in c_l:
                                nxr = grid_1D.shape[0]
                                grid_1D = np.transpose(grid_1D, ( 1,  0))
                                grid_1D = np.reshape(  grid_1D, (-1,  c,  nxr))
                                grid_1D = np.transpose(grid_1D, ( 1,  0,  2))
                                grid_1D = np.reshape(  grid_1D, ( c, -1))

                # 3D reshaping procedure, see text for details
                elif ndim == 3:
                        grid_1D = grid
                        for c in c_l:
                                nxr = grid_1D.shape[0]
                                nyr = grid_1D.shape[1]
                                grid_1D = np.transpose(grid_1D, ( 2,  1,  0))
                                grid_1D = np.reshape(  grid_1D, (-1,  c,  nyr, nxr))
                                grid_1D = np.transpose(grid_1D, ( 1,  0,  2,   3))
                                grid_1D = np.reshape(  grid_1D, ( c, -1,  c,   nxr))
                                grid_1D = np.transpose(grid_1D, ( 0,  2,  1,   3))
                                grid_1D = np.reshape(  grid_1D, ( c,  c, -1))

                else:
                        print("No reshaping procedure for ndim = %i" % ndim)
                        exit()

                X_grid[:,n] = grid_1D
        

        return X_grid



# def GridGen_Modified(nx, ny, nz, nt, nspat, ndim, nlev, finest, c_l, d_l):

#       # ---------- User inputs --------------------------------------------
#       lc_0   = 1/16 # fraction of grid cells fixed at lvl 0
#       lc_1   = 1/16
#       lc_2   = 1/16
#       lc_3   = 1/16
#       l_1    = 4/16 # fraction of grid cells at lvl 1
#       l_2    = 4/16
#       l_3    = 4/16

#       # ---------- Error checking 

#       # error for l<lc and other things

#       # ---------- Initialization of arrays and matrices
#       X_grid = np.zeros((nspat, nt), dtype=int) # matrix for POD operations
#       lc_arr = np.zeros((nlev))
#       l_arr  = np.zeros((nlev))
#       l_exp  = np.zeros(nlev) # computed (experiment) level fractions
#       lc_exp = np.zeros(nlev) # computed (experiment) level constant fractions
#       nccells_lc_arr  = np.zeros((nlev), dtype=int)
#       levels = np.arange(0, nlev)

#       # ---------- Compute quantities used to generate grids

#       # Find fractional values of lc_arr and l_arr
#       for l in levels:
#               exec('lc_arr[%i] = lc_%i' % (l,l))
#               if l>0: exec('l_arr[%i]  = l_%i'  % (l,l))
#       l_arr[0] = 1 - sum(l_arr)

#       # Compute the number of coarse cells we can tag for lc's
#       nc_cells  = nspat / c_l[0]**ndim

#       # Compute the number of coarse cells we need to tag for the lc
#       # Note: This is a cumulative sum because once we are done tagging the
#       #       lc of a coarser level, we then move to taking the current lc 
#       for l in levels:
#               nccells_lc_arr[l] = sum(nc_cells*lc_arr[0:l+1])

#       # ========== Begin operation to compute grids =======================

#       # Iterate in time computing new grids with each snapshot
#       for n in range(nt):
#               # Tag cells for lc
#               count = 0
#               grid = np.zeros((nx, ny, nz))
#               for i in range(0, nx, c_l[0]):
#                       for j in range(0, ny, c_l[0]):
#                               for k in range(0, nz, c_l[0]):
#                                       count += 1
#                                       for l in levels:
#                                               if count <= nccells_lc_arr[l]:
#                                                       grid[i:i+c_l[0], j:j+c_l[0], k:k+c_l[0]] = l
#                                                       break

#               # Find the indices that we can start at without disrupting lc
#               count = 0
#               flag = False
#               for i in range(0, nx, c_l[0]):
#                       for j in range(0, ny, c_l[0]):
#                               for k in range(0, nz, c_l[0]):
#                                       if count < nccells_lc_arr[-1]:
#                                               count += 1
#                                       else:
#                                               flag = True
#                                               break
#                               if flag:
#                                       break
#                       if flag:
#                               break
#               start_i = i
#               start_j = j
#               start_k = k

#               # Find the number of valid cells we can tag and number of cells
#               # we need to tag.
#               # Note: We say we can only tag cells of one level coarser, so 
#               #       only a l=2 can be tagged if it is a l=1 cell, etc. So,the
#               #       number of cells valid for l=2 is the number of levels 
#               #       tagged for l=1. A level that has finer grid cells must
#               #       account for the additional cells at those finer levels
#               nc_valid    = np.zeros(finest)
#               nc_valid[0] = nspat * (1 - sum(lc_arr)) # must not dirupt lc's
#               nc_tag      = np.zeros(finest)
#               nc_tag[0]   = nspat * (sum(l_arr[1:nlev] - lc_arr[1:nlev]))
#               for i in range(1,finest):
#                       nc_valid[i]   = nc_tag[i-1]
#                       nc_tag[i]  = nspat * sum(l_arr[i+1:nlev] - lc_arr[i+1:nlev])

#               # Strategy: iterate through levels up to finest-1 starting at 
#               # indices after lc tagged. Define a probability based on number
#               # of valid cells and number of cells we need to tag. Draw a 
#               # random number, and based on that relative to ratio of valid
#               # to tagged, decide if we tag and split that cell
#               for l in range(finest):

#                       # Iterate through x-dimension
#                       for i in range(start_i, nx, c_l[l]):
#                               if i >= start_i and i < start_i + c_l[0]:
#                                       j = start_j
#                               else:
#                                       j = 0

#                               # Iterate through y-dimension
#                               while j < ny:

#                                       if i >= start_i and i < start_i + c_l[0] and \
#                                          j >= start_j and j < start_j + c_l[0]:
#                                               k = start_k
#                                       else:
#                                               k = 0

#                                       # Iterate through z-dimension
#                                       while k < nz:
#                                               rval = rand.random()

#                                               # Ensure probability is a valid value (not inf/nan)
#                                               if nc_valid[l] == 0 and nc_tag[l] != 0:
#                                                       prob = 1
#                                               elif nc_valid[l] == 0 and nc_tag[l] == 0:
#                                                       prob = 0 
#                                               else:
#                                                       prob = (nc_tag[l]/nc_valid[l]) 

#                                               if l == 0:
#                                                       if rval < prob:   
#                                                               grid[i:i+c_l[l], j:j+c_l[l], k:k+c_l[l]] = l+1
#                                                               nc_tag[l] -= d_l[l]
#                                                       nc_valid[l] -= d_l[l]
                                                        
#                                               else:
#                                                       if grid[i,j,k] > l-1:
#                                                               if rval < prob:
#                                                                       grid[i:i+c_l[l], j:j+c_l[l], k:k+c_l[l]] = l+1
#                                                                       nc_tag[l] -= d_l[l]
#                                                               nc_valid[l] -= d_l[l]

#                                               k += c_l[l]
#                                       j += c_l[l]


#               # 1D, no reshaping required
#               if ndim == 1:
#                       grid_1D = np.squeeze(grid)

#               # 2D reshaping procedure, see text for details
#               elif ndim == 2:
#                       grid_1D = np.squeeze(grid)
#                       for c in c_l:
#                               nxr = grid_1D.shape[0]
#                               grid_1D = np.transpose(grid_1D, ( 1,  0))
#                               grid_1D = np.reshape(  grid_1D, (-1,  c,  nxr))
#                               grid_1D = np.transpose(grid_1D, ( 1,  0,  2))
#                               grid_1D = np.reshape(  grid_1D, ( c, -1))

#               # 3D reshaping procedure, see text for details
#               elif ndim == 3:
#                       grid_1D = grid
#                       for c in c_l:
#                               nxr = grid_1D.shape[0]
#                               nyr = grid_1D.shape[1]
#                               grid_1D = np.transpose(grid_1D, ( 2,  1,  0))
#                               grid_1D = np.reshape(  grid_1D, (-1,  c,  nyr, nxr))
#                               grid_1D = np.transpose(grid_1D, ( 1,  0,  2,   3))
#                               grid_1D = np.reshape(  grid_1D, ( c, -1,  c,   nxr))
#                               grid_1D = np.transpose(grid_1D, ( 0,  2,  1,   3))
#                               grid_1D = np.reshape(  grid_1D, ( c,  c, -1))

#               else:
#                       print("No reshaping procedure for ndim = %i" % ndim)
#                       exit()

#               X_grid[:,n] = grid_1D

#       # ========== Compute l and lc from data =============================

#       # Place holder for 
#       grid_lc = X_grid[:,0]

#       # Iterate through all time steps
#       for n in range(nt):

#               # Find fraction of grid that is lc based on if there is a
#               # difference in level from grid_lc. 
#               #
#               # Strategy: if there is a difference, mark that grid cell in 
#               # grid_lc as -1, so not there will also be a difference for the
#               # rest of time. At the end, we will just look at the final grid_lc
#               # to see how many cells of each level there are
#               grid_diff = X_grid[:,n] - grid_lc
#               grid_lc[grid_diff != 0] = -1

#               # Find fraction of grid at a particular level. 
#               #
#               # Strategy: simply add fractions for each time step then average
#               for l in levels:
#                       l_exp[l] += np.sum(X_grid[:,n] == l)/nspat

#       # Take average
#       l_exp = l_exp/nt

#       # Compute lc_exp from grid_lc
#       for l in levels:
#               lc_exp[l] = np.sum(grid_lc == l)/nspat

#       # ---------- Output results
#       print('\n The defined level percentages are:')
#       print(l_arr)
#       print('\n The measured level percentages are:')
#       print(l_exp)
#       print('\n The defined constant level percentages are:')
#       print(lc_arr)
#       print('\n The measured constant level percentages are:')
#       print(lc_exp)

#       return X_grid
