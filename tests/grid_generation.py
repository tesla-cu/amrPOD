import numpy as np
import os
import sys
sys.path.insert(1, '../source/')
from GenGrid import GenGrid

# ========== Main function ==========================================
if __name__ == '__main__':

    print('starting script to write synthetic AMR grids ...')

    # ---------- User defined inputs --------------------------------
    nx     = 24                   # x spatial points                  
    ny     = 24                   # y spatial points
    nz     = 24                   # z spatial points
    finest = 3                    # finest level of AMR in the domain
    nt     = 50                   # spanning nt
    ls     = np.array([15/27, 1/27, 1/27, 10/27]) # amount of l0 and l1
    # lcs    = np.array([1/4, 1/4]) # amount of lc0 and lc1
    lcs    = np.zeros((finest+1)) # amount of lc0 and lc1

    # Direction where /code/ lives
    basedir = '../../'

    # Directory where we want to store data
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory where we will write grid level data
    checkdir = datadir + 'check_grid_generation/'
    if not os.path.exists(checkdir):
        os.mkdir(checkdir)

    # ---------- Helpful quantities derived from user inputs --------
    nlev = finest + 1
    c_l  = np.zeros((nlev), dtype=int)
    d_l  = np.zeros((nlev), dtype=int)

    ndim  = 0          # num dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1
    nspat = nx*ny*nz
    
    for i in range(nlev):
        c_l[i] = 2**(finest-i)
        d_l[i] = (2**ndim)**(finest-i)

    # ---------- Generate data ----------------------------------------
    X_grid = np.zeros((nspat, nt), dtype=int)
    for n in range(nt):
        grid = GenGrid(nx, ny, nz, c_l, d_l, ls, lcs)
        X_grid[:,n] = np.reshape(grid, (-1))

    # ---------- Compute grid information from X_grid ---------------
    l_comp  = np.zeros((nlev)) # computed level fractions
    lc_comp = np.zeros((nlev)) # computed level constant fractions

    # Iterate through all time steps
    for n in range(nt):

        # Find fraction of grid at a particular level. Add fractions
        # for each time step then average
        for l in range(nlev):
            l_comp[l] += np.sum(X_grid[:,n] == l)/nspat

    # Take average
    l_comp = l_comp/nt
    print("l_comp = ", l_comp)

    # Compute lc_comp
    for l in range(nlev):
        for i in range(nspat):
            if np.all(X_grid[i,:] == l):
                lc_comp[l] += 1
    lc_comp = lc_comp/nspat
    print("lc_comp = ", lc_comp)

    # Get error amounts
    l_err  = (ls  - l_comp) /ls
    lc_err = (lcs - lc_comp)/lcs

    sim_info = open(checkdir + '/sim_info.txt', 'w')
    sim_info.write("finest:   %i\n" % finest)
    sim_info.write("nt:       %i\n" % nt)
    sim_info.write("nx:       %i\n" % nx)
    sim_info.write("ny:       %i\n" % ny)
    sim_info.write("nz:       %i\n" % nz)
    [sim_info.write("l%i:       %0.8f\n" % (i,l))  for i,l  in enumerate(ls)]
    [sim_info.write("l_comp%i:  %0.8f\n" % (i,l))  for i,l  in enumerate(l_comp)]
    [sim_info.write("l_err%i:   %0.8e\n" % (i,l))  for i,l  in enumerate(l_err)]
    [sim_info.write("lc%i:      %0.8f\n" % (i,lc)) for i,lc in enumerate(lcs)]
    [sim_info.write("lc_comp%i: %0.8f\n" % (i,lc)) for i,lc in enumerate(lc_comp)]
    [sim_info.write("lc_err%i:  %0.8e\n" % (i,lc)) for i,lc in enumerate(lc_err)]
    sim_info.close()

    print('\tdone checking synthetic AMR grids')

