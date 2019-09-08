import numpy as np
import matplotlib.pyplot as plt 
from numpy import linalg as LA
import matplotlib.pyplot as plt

from GenGrid import GenGrid

from R_CPU   import   compute_R_CPU
from R_TC    import   compute_R_TC

from Phi_CPU import   compute_Phi_CPU
from Phi_TC  import   compute_Phi_TC

from A_CPU   import   compute_A_CPU
from A_TC    import   compute_A_TC

def Compute_POD(gen_grid, nx, ny, nz, finest, l_fracs, lc_fracs, nt, TC_CPU='TC', amr_datadir=None):
 
    nspat = nx*ny*nz    
    nlev  = finest + 1
    c_l   = np.zeros((nlev), dtype=int)
    d_l   = np.zeros((nlev), dtype=int)

    # ---------- Helpful quantities derived from user inputs
    ndim  = 0          # num dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1
    levels = np.arange(0, nlev)
    
    for i in range(nlev):
        c_l[i]    = 2**(finest-i)        
        d_l[i]    = (2**ndim)**(finest-i)

    # print('Computing for nt = %i, l_0=%0.8f, l_1 = %0.8f, l_2 = %0.8f' % (nt, l_fracs[0],l_fracs[1], l_fracs[2]))

    #--------- Load or generate data
    X      = np.zeros((nspat,nt))
    X_grid = np.zeros((nspat,nt), dtype=int)
    for n in range(nt):

        if gen_grid:
            grid = GenGrid(nx, ny, nz, c_l, d_l, l_fracs, lc_fracs)
            data = grid
        else:
            print('need to write this code!!')
            sys.exit()

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

    #--------- Compute grid information from X_grid
    l_comp  = np.zeros((nlev)) # computed level fractions
    lc_comp = np.zeros((nlev)) # computed level constant fractions

    # Place holder for 
    grid_lc = np.zeros((nlev))

    # Iterate through all time steps
    for n in range(nt):

        # Find fraction of grid that is lc based on if there is a
        # difference in level from grid_lc. 
        #
        # Strategy: if there is a difference, mark that grid cell in 
        # grid_lc as -1, so not there will also be a difference for the
        # rest of time. At the end, we will just look at the final grid_lc
        # to see how many cells of each level there are
        #grid_diff = X_grid[:,n] - grid_lc
        #grid_lc[grid_diff != 0] = -1

        # Find fraction of grid at a particular level. 
        #
        # Strategy: simply add fractions for each time step then average
        for l in levels:
            l_comp[l] += np.sum(X_grid[:,n] == l)/nspat
        #print(l_comp)

    # Take average
    l_comp = l_comp/nt

    # Compute lc_comp from grid_lc
    for l in levels:
        for i in range(nspat):
            if np.all(X_grid[i,:] == l):
                lc_comp[l] += 1
    lc_comp = lc_comp/nspat
    # print('lc_comp= ', lc_comp)

    #--------- Calculate POD with matrix operations
    X_tp             = np.transpose(X)
    R                = np.matmul(X_tp, X)
    Lambda, Psi      = LA.eig(R)
    Lambda_inv_sqrt  = 1/np.sqrt(Lambda)
    Lambda_diag      = np.diag(Lambda_inv_sqrt)
    Phi              = np.matmul(X,Psi)
    Phi              = np.matmul(Phi, Lambda_diag)
    A                = np.matmul(X_tp, Phi)
    Lambda           = np.diag(Lambda)

    #--------- Get CPU time and time complexity for R
    # time_R_imp[int_arr], time_R_unalt[int_arr], time_R_theory[int_arr] = compute_R_CPU(X, X_grid, R, d_l, nt, nspat)
    TC_R_imp, TC_R_unalt = compute_R_TC(X_grid, R, d_l, nt, nspat)

    #--------- Get CPU time and time complexity for Phi Method 1
    method = 1
    # time_Phi_1_imp[int_arr], time_Phi_1_unalt[int_arr] = compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat)
    TC_Phi_1_imp, TC_Phi_1_unalt     = compute_Phi_TC(X_grid, method, d_l, nt, nspat, finest)

    #--------- Get CPU time and time complexity for Phi Method 4
    # method = 4
    # time_Phi_4_imp[int_arr], time_Phi_4_unalt[int_arr] = compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat)
    # TC_Phi_4_imp, TC_Phi_4_unalt    = compute_Phi_TC(X_grid, method, d_l, nt, nspat)

    #--------- Get CPU time and time complexity for Phi Method 4
    method = 5
    # time_Phi_5_imp, time_Phi_5_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, method, Phi, d_l, nt, nspat)
    TC_Phi_5_imp, TC_Phi_5_unalt   = compute_Phi_TC(X_grid, method, d_l, nt, nspat, finest)


    #--------- Get CPU time and time complexity for A
    # time_A_unalt[int_arr], time_A_imp[int_arr], time_A_theory[int_arr] = compute_A_CPU(X, X_grid, Phi, A, d_l, nt, nspat)
    TC_A_imp, TC_A_unalt = compute_A_TC(X_grid, R, d_l, nt, nspat, finest)

    return TC_R_imp, TC_R_unalt, TC_Phi_1_imp, TC_Phi_1_unalt, TC_Phi_5_imp, TC_Phi_5_unalt, TC_A_imp, TC_A_unalt
