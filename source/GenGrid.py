import numpy as np
import random as rand

# =========================================================================== #
# Function to generate synthetic AMR data based on user inputs
#
# Inputs:
# - nx     : number of cells in the x direction
# - ny     : number of cells in the x direction
# - nz     : number of cells in the x direction
# - nt     : number of time steps
# - c_l    : number of repeated computations in a particular dimension
# - d_l    : total number of repeated computations for a cell
# - l_arr  : array of the fraction of the domain at a particular level
# - lc_arr : array of the fraction of the domain at a particular level that 
#            is held constand
#
# Outputs:
# - grid : generated grid with specified parameters
# =========================================================================== #
def GenGrid(nx, ny, nz, c_l, d_l, l_arr, lc_arr):

    # ------------ Helpful quantities derived from user inputs -------------- #
    finest = len(l_arr)-1 # finest lvl of AMR
    nspat = int(nx*ny*nz) # num spatial points
    nlev  = finest + 1    # number of  levels
    ndim  = 0             # number of dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1 

    # -------------------------- Error checking ----------------------------- #

    # Check sum(l_arr) is close to 1
    if not np.isclose(np.sum(l_arr), 1.0, atol=1e-12):
        print('Error: l_arr must sum to 1.0!')
        exit()

    # Check lc is less than or equal to l for a given level
    if any(l_arr - lc_arr < 0.0):
        print('Error: lc cannot be greater than l!')
        exit()

    # Check we have enough resolution and enough cells to generate the grid
    if not (nx/d_l[-2]).is_integer() and not (nx/d_l[0]).is_integer() and \
       not (ny/d_l[-2]).is_integer() and not (ny/d_l[0]).is_integer() and  \
       not (nz/d_l[-2]).is_integer() and not (nz/d_l[0]).is_integer():
        print('Error: AMR blocks are not evenly divisible by nspat!')
        exit()
    

    # --------------- Initialization of arrays and matrices ----------------- #
    nccells_lc_arr  = np.zeros((nlev), dtype=int)
    levels = np.arange(0, nlev)

    # ------------- Compute quantities used to generate grids --------------- #

    # Compute the number of coarse cells we can tag for lc's
    nc_cells  = nspat/d_l[0]

    # Compute the number of coarse cells we need to tag for the lc
    # Note: This is a cumulative sum because once we are done tagging the
    #       lc of a coarser level, we then move to taking the current lc 
    for l in levels:
        nccells_lc_arr[l] = sum(nc_cells*lc_arr[0:l+1])

    # ================= Begin operation to compute grids ==================== #

    # ------------------------- Tag cells for lc ---------------------------- #
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

    # ----- Find the indices that we can start at without disrupting lc ----- # 
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

    # Grid is done being generated
    return grid

