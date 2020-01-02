import numpy as np
import os
import sys
from numpy import linalg as LA

sys.path.insert(1, '../../source/')

from GenGrid     import GenGrid
from Reshape_AMR import Reshape_AMR

from R_CPU   import compute_R_CPU
from R_TC    import compute_R_TC

from Phi_CPU import compute_Phi_CPU
from Phi_TC  import compute_Phi_TC

from A_CPU   import compute_A_CPU
from A_TC    import compute_A_TC


if __name__ == '__main__':

    print('starting script to perform POD on AMR grids ...')

    # =========================================================================
    # User defined inputs 
    # =========================================================================

    nx     = 512   # x spatial points                  
    ny     = 512   # y spatial points
    nz     = 1     # z spatial points
    finest = 3     # finest level of AMR in the domain
    nt     = 10    # number of time steps

    # Direction where /code/ lives
    basedir = '../../../'

    # Directory where AMR data is stored
    amr_datadir ='/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/'+\
        'data/slice/x0.000-0.000_y-1.000-1.000_z0.000-2.000_t40.0000-42.0000/'

    # Directory where we want to store computed data
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory that describes where we want to store computed R, Phi, and A
    studydir = datadir + 'python_vs_matlab/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # =========================================================================
    # Begin POD computation
    # =========================================================================

    # Helpful quantities derived from user inputs 
    nspat = nx*ny*nz    
    nlev  = finest + 1
    ndim  = 0            # num dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1
    levels = np.arange(0, nlev)
    
    # Arrays for repetition
    c_l   = np.zeros((nlev), dtype=int)
    d_l   = np.zeros((nlev), dtype=int)
    for i in range(nlev):
        c_l[i]    = 2**(finest-i)
        d_l[i]    = (2**ndim)**(finest-i)

    # Load data
    X      = np.zeros((nspat,nt))
    X_grid = np.zeros((nspat,nt), dtype=int)
    for n in range(nt):

        # Load data from file
        grid = np.fromfile(amr_datadir + 'grid_level%05d.bin' % n).astype(int)
        data = np.fromfile(amr_datadir + 'z_velocity%05d.bin' % n)

        # Perform reshaping
        X[:,n]      = Reshape_AMR(nx, ny, nz, finest, data, 'forward')
        X_grid[:,n] = Reshape_AMR(nx, ny, nz, finest, grid, 'forward')

    # -------------------------------------------------------------------------
    # Calculate POD with matrix operations
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Calculate POD with iterative operations 
    # -------------------------------------------------------------------------
    R_imp,  R_unalt  = compute_R_CPU  (X, X_grid,                 R,   d_l, \
        nt, nspat, finest)
    P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, Phi, d_l, \
        nt, nspat, finest)
    P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, Phi, d_l, \
        nt, nspat, finest)
    A_imp,  A_unalt  = compute_A_CPU  (X, X_grid, Phi,            A,   d_l, \
        nt, nspat, finest)

    # Reshape back to original shape
    for n in range(nt):
        Phi[:,n] = Reshape_AMR(nx, ny, nz, finest, Phi[:,n], 'reverse')

    # =========================================================================
    # Save data 
    # =========================================================================
    
    print('writing out text files ...')    
    np.savetxt(studydir + "/R_py.txt",   R)
    np.savetxt(studydir + "/Phi_py.txt", Phi)
    np.savetxt(studydir + "/A_py.txt",    A)
    print('done performing POD on AMR grids.')




