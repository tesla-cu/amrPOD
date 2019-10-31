import numpy as np
import os
import sys
from numpy import linalg as LA
sys.path.insert(1, '../source/')
from Compute_POD_check import Compute_POD_check



# This function simply enables us to do a parallel run where we have
# multiple arguments for the function
def parallel_run(args):
   return Compute_POD(*args)


if __name__ == '__main__':

    print('starting script to perform POD on AMR grids ...')

    # ---------- User defined inputs -------------------------------- 
    # gen_grid    = False  # are we generating synthetic data?
    # compute_tc  = True  # are we computing the time complexity?
    # compute_cpu = True # are we computing the cpu time?
    nx          = 512   # x spatial points                  
    ny          = 512    # y spatial points
    nz          = 1     # z spatial points
    finest      = 3     # finest level of AMR in the domain
    nt          = 10     # number of time steps
    # lc_fracs     = np.array([1/16, 0/16, 0/16])

    # Variables we will iterate through
    # xvals = nt_arr # what will be plotted on the x-axis
    # yvals = l1_arr # what will be plotted on the y-axis

    # Direction where /code/ lives
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
    studydir = datadir + 'check_MAT/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # Directory that holds all txt files of POD variables
    # txtdir = studydir + 'txt_files/'
    # if not os.path.exists(txtdir):
    #     os.mkdir(txtdir)

    # Directory that holds all images generated from txt files
    # imgdir = studydir + 'images/'
    # if not os.path.exists(imgdir):
    #     os.mkdir(imgdir)


    R, Phi, A = Compute_POD_check(nx, ny, nz, finest, nt,  amr_datadir)

    # ========== Save data
    print('writing out text files ...')
    np.savetxt(studydir + "/R_py.txt",   R)
    np.savetxt(studydir + "/Phi_py.txt", Phi)
    np.savetxt(studydir + "/A_py.txt",    A)
    print('done performing POD on AMR grids.')




