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

    # User defined inputs -----------------------------------------------------
    nx          = 128   # x spatial points                  
    ny          = 128   # y spatial points
    nz          = 64    # z spatial points
    nt          = 5     # spanning nt
    finest      = 5     # finest level of AMR in the domain
    ls          = np.array([4/16, 4/16, 2/16, 2/16, 2/16, 2/16])
    lcs         = np.array([1/16, 1/16, 0/16, 0/16, 0/16, 2/16]) 
    # finest      = 3     # finest level of AMR in the domain
    # ls          = np.array([6/16, 4/16, 4/16, 2/16])
    # lcs         = np.array([2/16, 1/16, 0/16, 0/16]) 

    # Direction where /code/ livesc
    basedir = '../../../'

    # Directory where we want to store data on speed up
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory that describes the study we are looking at
    studydir = datadir + 'reshaping/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # Helpful quantities derived from user inputs -----------------------------
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

    # Generate data ----------------------------------------------------------
    X      = np.zeros((nspat,nt))
    X_grid = np.zeros((nspat,nt), dtype=int)
    for n in range(nt):

        # Generate data
        grid = GenGrid(nx, ny, nz, c_l, d_l, ls, lcs)
        data = grid + 1.5 # force data to not align with grid 

        # Reshape into 1D array
        grid_1D = np.reshape(grid, (nspat))
        data_1D = np.reshape(data, (nspat))

        # Assign new data to corresponding X matrix
        X_grid[:,n] = grid_1D
        X[:,n]      = data_1D

    # Compute grid information from X_grid ------------------------------------
    l_comp  = np.zeros((nlev)) # computed level fractions
    lc_comp = np.zeros((nlev)) # computed level constant fractions

    # Compute l_comp
    for n in range(nt):
        for l in levels:
            l_comp[l] += np.sum(X_grid[:,n] == l)/nspat
    l_comp = l_comp/nt
    print("l_comp = ", l_comp)

    # Compute lc_comp
    for l in levels:
        for i in range(nspat):
            if np.all(X_grid[i,:] == l):
                lc_comp[l] += 1
    lc_comp = lc_comp/nspat
    print("lc_comp = ", lc_comp)

    # =========================================================================
    # Calculate POD without reshaping procedure
    # =========================================================================
    X_tp        = np.transpose(X)
    R_nr        = np.matmul(X_tp, X)
    Lambda, Psi = LA.eig(R_nr)
    idx_eig     = np.argsort(Lambda) # sort eigenvalues
    Lambda      = Lambda[idx_eig]
    Psi         = Psi[:,idx_eig]
    Phi         = np.matmul(X,Psi)
    Phi_nr      = np.matmul(Phi, np.diag(1/np.sqrt(Lambda)))
    A_nr        = np.matmul(X_tp, Phi_nr)
    Lambda      = np.diag(Lambda) # make this a matrix

    # =========================================================================
    # Calculate POD with reshaping procedure
    # =========================================================================

    # Perform reshaping
    for n in range(nt):
        X[:,n]      = Reshape_AMR(nx, ny, nz, finest, X[:,n],      'forward')
        X_grid[:,n] = Reshape_AMR(nx, ny, nz, finest, X_grid[:,n], 'forward')


    # Calculate POD with reshaped data
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

    # Calculate POD with iterative operations
    R_imp,  R_unalt  = compute_R_CPU  (X, X_grid,                 R,   d_l, \
        nt, nspat, finest)
    P1_imp, P1_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 1, Phi, d_l, \
        nt, nspat, finest)
    P2_imp, P2_unalt = compute_Phi_CPU(X, X_grid, Psi, Lambda, 2, Phi, d_l, \
        nt, nspat, finest)
    A_imp,  A_unalt  = compute_A_CPU  (X, X_grid, Phi,            A,   d_l, \
        nt, nspat, finest)

    R_imp,  R_unalt  = compute_R_TC   (X, X_grid,                 R,   d_l, \
        nt, nspat, finest, 1, 1, 1, 1, 1)
    P1_imp, P1_unalt = compute_Phi_TC (X, X_grid, Psi, Lambda, 1, Phi, d_l, \
        nt, nspat, finest, 1, 1, 1, 1, 1)
    P2_imp, P2_unalt = compute_Phi_TC (X, X_grid, Psi, Lambda, 2, Phi, d_l, \
        nt, nspat, finest, 1, 1, 1, 1, 1)
    A_imp,  A_unalt  = compute_A_TC   (X, X_grid, Phi,            A,   d_l, \
        nt, nspat, finest, 1, 1, 1, 1, 1)

    # Reshape back to original shape
    for n in range(nt):
        Phi[:,n] = Reshape_AMR(nx, ny, nz, finest, Phi[:,n], 'reverse')

    # =========================================================================
    # Compare reshaped and non-reshaped
    # =========================================================================

    # Compute relative error
    R_err   = np.max( abs(R    -     R_nr)    /     R_nr)
    Phi_err = np.max((abs(Phi) - abs(Phi_nr)) / abs(Phi_nr))
    A_err   = np.max((abs(A)   - abs(A_nr))   / abs(A_nr)) 

    print('Error in R:   %0.8e' % R_err)
    print('Error in Phi: %0.8e' % Phi_err)
    print('Error in A:   %0.8e' % A_err)

    error_info = open(studydir + '/error_info.txt', 'w')
    error_info.write('Error in R:   %0.8e\n' % R_err)
    error_info.write('Error in Phi: %0.8e\n' % Phi_err)
    error_info.write('Error in A:   %0.8e\n' % A_err)
    error_info.close()

    