import numpy as np
import os
import sys
sys.path.insert(1, '../source/')
from GenGrid import GenGrid

# =============================================================================
# Main function 
# =============================================================================

if __name__ == '__main__':

    print('starting script to write synthetic AMR grids ...')

    # User defined inputs -----------------------------------------------------

    # Grid information
    nx     = 128                   # x spatial points                  
    ny     = 128                   # y spatial points
    nz     = 64                    # z spatial points
    nt     = 5                     # spanning nt
    # finest = 1                    # finest level of AMR in the domain
    # ls     = np.array([1/2, 1/2]) # amount of l0, l1, ...
    # lcs    = np.array([1/4, 1/4]) # amount of lc0, lc1, ...
    finest = 3                    # finest level of AMR in the domain
    ls     = np.array([1/2, 1/4, 1/8, 1/8]) # amount of l0, l1, ...
    lcs    = np.array([1/4, 1/8, 0/8, 0/8]) # amount of lc0, lc1, ...

    # Direction where /code/ lives
    basedir = '../../'

    # Directory where we want to store data
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory where we will write grid level data
    amr_datadir = datadir + 'AMR_fortran_test/'
    if not os.path.exists(amr_datadir):
        os.mkdir(amr_datadir)

    # Helpful quantities derived from user inputs -----------------------------
    nlev = finest + 1
    c_l  = np.zeros((nlev), dtype=int)
    d_l  = np.zeros((nlev), dtype=int)

    ndim  = 0          # num dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1
    
    for i in range(nlev):
        c_l[i] = 2**(finest-i)
        d_l[i] = (2**ndim)**(finest-i)

    # Generate data -----------------------------------------------------------
    for n in range(nt):
        grid = GenGrid(nx, ny, nz, c_l, d_l, ls, lcs)
        grid = np.transpose(grid, (2, 1, 0)) # need to load into fortran
        grid.tofile(amr_datadir + 'grid_level%05d.bin' % n)

    sim_info = open(amr_datadir + '/sim_info.txt', 'w')
    sim_info.write("finest:   %i\n" % finest)
    sim_info.write("nt:       %i\n" % nt)
    sim_info.write("nx:       %i\n" % nx)
    sim_info.write("ny:       %i\n" % ny)
    sim_info.write("nz:       %i\n" % nz)
    [sim_info.write("l%i:       %0.8f\n" % (i,l))  for i,l  in enumerate(ls)]
    [sim_info.write("lc%i:      %0.8f\n" % (i,lc)) for i,lc in enumerate(lcs)]
    sim_info.close()

    print('\tdone writing synthetic AMR grids')

