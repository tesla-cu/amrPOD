import numpy as np
import os
import sys
from numpy import linalg as LA
sys.path.insert(1, '../source/')
from Compute_POD_correct import Compute_POD_correct

# ========== Main function ==========================================
if __name__ == '__main__':

    print('starting script to perform POD on AMR grids ...')

    # ---------- User defined inputs --------------------------------
    nx          = 128   # x spatial points                  
    ny          = 128    # y spatial points
    nz          = 64   # z spatial points
    finest      = 5     # finest level of AMR in the domain
    nt          = 5     # spanning nt
    ls          = np.array([4/16, 4/16, 2/16, 2/16, 2/16, 2/16])
    lcs         = np.array([1/16, 1/16, 0/16, 0/16, 0/16, 2/16]) 

    # Direction where /code/ livesc
    # basedir = '/Users/samsimonswellin/desktop/'
    basedir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/'

    # Directory where AMR data is stored
    # amr_datadir = '/Users/samsimonswellin/desktop/' + \
    #     'x0.000-0.000_y-1.000-1.000_z0.000-2.000_t40.0000-42.0000/'
    amr_datadir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/' + \
        'slice/x0.000-0.000_y-1.000-1.000_z0.000-2.000_t40.0000-42.0000/'

    # Directory where we was to store data on speed up
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory that describes the study we are looking at
    studydir = datadir + 'check_internal/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # ---------- Compute POD ----------------------------------------
    R_err, Phi_err, A_err = Compute_POD_correct(nx, ny, nz, finest, nt, ls, lcs)

    print('Error in R:   %0.8e' % R_err)
    print('Error in Phi: %0.8e' % Phi_err)
    print('Error in A:   %0.8e' % A_err)

    error_info = open(studydir + '/error_info.txt', 'w')
    error_info.write('Error in R:   %0.8e\n' % R_err)
    error_info.write('Error in Phi: %0.8e\n' % Phi_err)
    error_info.write('Error in A:   %0.8e\n' % A_err)
    error_info.close()


