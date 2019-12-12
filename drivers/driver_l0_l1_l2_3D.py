import numpy as np
import matplotlib.pyplot as plt 
import multiprocessing as mp
import os
import sys
from numpy import linalg as LA
sys.path.insert(1, '../source/')
sys.path.insert(1, '../plotting/')
from Compute_POD import Compute_POD
from plot_single import plot_single

# This function simply enables us to do a parallel run where we have
# multiple arguments for the function
def parallel_run(args):
   return Compute_POD(*args)

# ========== Main function ==========================================
if __name__ == '__main__':

    print('starting script to perform POD on AMR grids ...')

    # ---------- User defined inputs -------------------------------- 
    gen_grid    = True  # are we generating synthetic data?
    compute_tc  = True  # are we computing the time complexity?
    compute_cpu = False # are we computing the cpu time?
    nx          = 24   # x spatial points                  
    ny          = 24   # y spatial points
    nz          = 24    # z spatial points
    finest      = 3     # finest level of AMR in the domain
    nsample     = 8    # number of samples for each parameter set
    nt          = 50     # spanning nt
    l0_arr      = np.arange(0.0, 10/27, 1/27) # spanning l0
    l1_arr      = np.arange(0.0, 10/27, 1/27) # spanning l1
    l2_arr      = np.arange(0.0, 10/27, 1/27) # spanning l2
    lcs         = np.zeros((finest+1)) # fraction of grid that stays constant in time

    # Direction where /code/ livesc
    # basedir = '/Users/samsimonswellin/desktop/'
    # basedir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/'
    basedir = '../../'

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
    studydir = datadir + 'l0_l1_l2_3D/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # Directory that holds all txt files of TC and CPU
    npydir = studydir + 'npy_files/'
    if not os.path.exists(npydir):
        os.mkdir(npydir)

    # ---------- Useful things based on inputs ----------------------

    # Variables we will iterate through
    xvals = l0_arr
    yvals = l1_arr
    zvals = l2_arr

    # Useful things we will need
    ndim  = 0          # num dimensions
    if nx > 1: ndim += 1 
    if ny > 1: ndim += 1 
    if nz > 1: ndim += 1 
    nlev = finest+1
    d_l  = np.zeros((nlev), dtype=int)
    for i in range(nlev):
        d_l[i]    = (2**ndim)**(finest-i)

    
    # Initialize matrices to store data
    #    - dim 1: xvar
    #    - dim 2: yvar
    #    - dim 3: zvar
    #    - dim 4: the different algorithms
    #        - index 0: R
    #        - index 1: Phi, method 1
    #        - index 2: Phi, method 2
    #        - index 3: A
    #    - avg is the average of all the samples
    #    - rms is the rms of all the samples
    if compute_tc:
        TC_avg_imp   = np.zeros((len(xvals), len(yvals), len(zvals), 4), dtype=int)
        TC_avg_unalt = np.zeros((len(xvals), len(yvals), len(zvals), 4), dtype=int)
        TC_rms_imp   = np.zeros((len(xvals), len(yvals), len(zvals), 4), dtype=int)
        TC_rms_unalt = np.zeros((len(xvals), len(yvals), len(zvals), 4), dtype=int)

    if compute_cpu:
        CPU_avg_imp   = np.zeros((len(xvals), len(yvals), len(zvals), 4))
        CPU_avg_unalt = np.zeros((len(xvals), len(yvals), len(zvals), 4))
        CPU_rms_imp   = np.zeros((len(xvals), len(yvals), len(zvals), 4))
        CPU_rms_unalt = np.zeros((len(xvals), len(yvals), len(zvals), 4))


    # Start parallel processing
    nthread = mp.cpu_count()
    print('starting pool with %i threads ...' % nthread)
    pool = mp.Pool(processes=nthread)

    # ========== Begin computation of various parameters ============

    for ixval, xval in enumerate(xvals):
        for iyval, yval in enumerate(yvals):
            for izval, zval in enumerate(zvals):
                l0 = xval # less confusing
                l1 = yval # less confusing
                l2 = zval # less confusing

                print('l0 = %0.6f, l1 = %0.6f, l2 = %0.6f' % (l0, l1, l2))

                # Since we have l in our vars, we need to get ls
                ls = np.array([l0, l1, l2, 1.0-l0-l1-l2])

                # ---------- Compute time complexity --------------------
                if compute_tc:
                    # Generate tuple that stores all relevant parameters
                    #    - note, we need list of these equal to number of samples
                    if gen_grid:
                        arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt, \
                            'TC')
                    else:
                        arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt, \
                            'TC', amr_datadir)
                    arg_list = []
                    [arg_list.append(arg_tuple) for i in range(nsample)]

                    # Farm out nsample to each processor
                    res_tuple = pool.map(parallel_run, arg_list)

                    for i in range(8):
                        res = np.zeros((nsample), dtype=int)

                        for j in range(nsample):
                            res[j] = res_tuple[j][i]

                        res_avg = np.mean(res)
                        res_rms = np.std(res-res_avg)

                        if   i==0: 
                            TC_avg_imp[ixval, iyval, izval, 0]   = res_avg
                            TC_rms_imp[ixval, iyval, izval, 0]   = res_rms
                        elif i==1: 
                            TC_avg_unalt[ixval, iyval, izval, 0] = res_avg
                            TC_rms_unalt[ixval, iyval, izval, 0] = res_rms
                        elif i==2: 
                            TC_avg_imp[ixval, iyval, izval, 1]   = res_avg
                            TC_rms_imp[ixval, iyval, izval, 1]   = res_rms
                        elif i==3: 
                            TC_avg_unalt[ixval, iyval, izval, 1] = res_avg
                            TC_rms_unalt[ixval, iyval, izval, 1] = res_rms
                        elif i==4: 
                            TC_avg_imp[ixval, iyval, izval, 2]   = res_avg
                            TC_rms_imp[ixval, iyval, izval, 2]   = res_rms
                        elif i==5: 
                            TC_avg_unalt[ixval, iyval, izval, 2] = res_avg
                            TC_rms_unalt[ixval, iyval, izval, 2] = res_rms
                        elif i==6: 
                            TC_avg_imp[ixval, iyval, izval, 3]   = res_avg
                            TC_rms_imp[ixval, iyval, izval, 3]   = res_rms
                        elif i==7: 
                            TC_avg_unalt[ixval, iyval, izval, 3] = res_avg
                            TC_rms_unalt[ixval, iyval, izval, 3] = res_rms

                # ---------- Compute CPU time ---------------------------
                if compute_cpu:
                    # Generate tuple that stores all relevant parameters
                    #    - note, we need list of these equal to number of samples
                    if gen_grid:
                        arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt, \
                            'CPU')
                    else:
                        arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt, \
                            'CPU', amr_datadir)
                    arg_list = []
                    [arg_list.append(arg_tuple) for i in range(nsample)]

                    # Farm out nsample to each processor
                    res_tuple = pool.map(parallel_run, arg_list)

                    for i in range(8):
                        res = np.zeros((nsample))

                        for j in range(nsample):
                            res[j] = res_tuple[j][i]

                        res_avg = np.mean(res)
                        res_rms = np.std(res-res_avg)

                        if   i==0: 
                            CPU_avg_imp[ixval, iyval, izval, 0]   = res_avg
                            CPU_rms_imp[ixval, iyval, izval, 0]   = res_rms
                        elif i==1: 
                            CPU_avg_unalt[ixval, iyval, izval, 0] = res_avg
                            CPU_rms_unalt[ixval, iyval, izval, 0] = res_rms
                        elif i==2: 
                            CPU_avg_imp[ixval, iyval, izval, 1]   = res_avg
                            CPU_rms_imp[ixval, iyval, izval, 1]   = res_rms
                        elif i==3: 
                            CPU_avg_unalt[ixval, iyval, izval, 1] = res_avg
                            CPU_rms_unalt[ixval, iyval, izval, 1] = res_rms
                        elif i==4: 
                            CPU_avg_imp[ixval, iyval, izval, 2]   = res_avg
                            CPU_rms_imp[ixval, iyval, izval, 2]   = res_rms
                        elif i==5: 
                            CPU_avg_unalt[ixval, iyval, izval, 2] = res_avg
                            CPU_rms_unalt[ixval, iyval, izval, 2] = res_rms
                        elif i==6: 
                            CPU_avg_imp[ixval, iyval, izval, 3]   = res_avg
                            CPU_rms_imp[ixval, iyval, izval, 3]   = res_rms
                        elif i==7: 
                            CPU_avg_unalt[ixval, iyval, izval, 3] = res_avg
                            CPU_rms_unalt[ixval, iyval, izval, 3] = res_rms

    pool.close()
    pool.join()

    # ========== Save data ==========================================

    # Write out data from this run
    sim_info = open(npydir + '/sim_info.txt', 'w')
    sim_info.write("gen_grid: %s\n" % gen_grid)
    sim_info.write("finest:   %i\n" % finest)
    sim_info.write("nsample:  %i\n" % nsample)
    sim_info.write("nx:       %i\n" % nx)
    sim_info.write("ny:       %i\n" % ny)
    sim_info.write("nz:       %i\n" % nz)
    sim_info.write("nt:       %i\n" % nt)
    [sim_info.write("lc%i:      %0.8f\n" % (i,lc)) for i,lc in enumerate(lcs)]
    sim_info.write("xvar:     %s\n" % "l0")
    sim_info.write("x_0:      %0.8f\n" % xvals[0])
    sim_info.write("x_inc:    %0.8f\n" % np.diff(xvals[0:2]))
    sim_info.write("x_end:    %0.8f\n" % xvals[-1])
    sim_info.write("yvar:     %s\n" % "l1")
    sim_info.write("y_0:      %0.8f\n" % yvals[0])
    sim_info.write("y_inc:    %0.8f\n" % np.diff(yvals[0:2]))
    sim_info.write("y_end:    %0.8f\n" % yvals[-1])
    sim_info.write("zvar:     %s\n" % "l2")
    sim_info.write("z_0:      %0.8f\n" % zvals[0])
    sim_info.write("z_inc:    %0.8f\n" % np.diff(zvals[0:2]))
    sim_info.write("z_end:    %0.8f\n" % zvals[-1])
    sim_info.close()

    # Write out time complexity data
    if compute_tc:
        np.save(npydir + "/TC_R_avg_imp.npy",    TC_avg_imp[:,:,:,0])
        np.save(npydir + "/TC_R_avg_unalt.npy",  TC_avg_unalt[:,:,:,0])
        np.save(npydir + "/TC_R_rms_imp.npy",    TC_rms_imp[:,:,:,0])
        np.save(npydir + "/TC_R_rms_unalt.npy",  TC_rms_unalt[:,:,:,0])

        np.save(npydir + "/TC_P1_avg_imp.npy",   TC_avg_imp[:,:,:,1])
        np.save(npydir + "/TC_P1_avg_unalt.npy", TC_avg_unalt[:,:,:,1])
        np.save(npydir + "/TC_P1_rms_imp.npy",   TC_rms_imp[:,:,:,1])
        np.save(npydir + "/TC_P1_rms_unalt.npy", TC_rms_unalt[:,:,:,1])

        np.save(npydir + "/TC_P2_avg_imp.npy",   TC_avg_imp[:,:,:,2])
        np.save(npydir + "/TC_P2_avg_unalt.npy", TC_avg_unalt[:,:,:,2])
        np.save(npydir + "/TC_P2_rms_imp.npy",   TC_rms_imp[:,:,:,2])
        np.save(npydir + "/TC_P2_rms_unalt.npy", TC_rms_unalt[:,:,:,2])

        np.save(npydir + "/TC_A_avg_imp.npy",    TC_avg_imp[:,:,:,3])
        np.save(npydir + "/TC_A_avg_unalt.npy",  TC_avg_unalt[:,:,:,3])
        np.save(npydir + "/TC_A_rms_imp.npy",    TC_rms_imp[:,:,:,3])
        np.save(npydir + "/TC_A_rms_unalt.npy",  TC_rms_unalt[:,:,:,3])

    # Write out CPU time data
    if compute_cpu:
        np.save(npydir + "/CPU_R_avg_imp.npy",    CPU_avg_imp[:,:,:,0])
        np.save(npydir + "/CPU_R_avg_unalt.npy",  CPU_avg_unalt[:,:,:,0])
        np.save(npydir + "/CPU_R_rms_imp.npy",    CPU_rms_imp[:,:,:,0])
        np.save(npydir + "/CPU_R_rms_unalt.npy",  CPU_rms_unalt[:,:,:,0])

        np.save(npydir + "/CPU_P1_avg_imp.npy",   CPU_avg_imp[:,:,:,1])
        np.save(npydir + "/CPU_P1_avg_unalt.npy", CPU_avg_unalt[:,:,:,1])
        np.save(npydir + "/CPU_P1_rms_imp.npy",   CPU_rms_imp[:,:,:,1])
        np.save(npydir + "/CPU_P1_rms_unalt.npy", CPU_rms_unalt[:,:,:,1])

        np.save(npydir + "/CPU_P2_avg_imp.npy",   CPU_avg_imp[:,:,:,2])
        np.save(npydir + "/CPU_P2_avg_unalt.npy", CPU_avg_unalt[:,:,:,2])
        np.save(npydir + "/CPU_P2_rms_imp.npy",   CPU_rms_imp[:,:,:,2])
        np.save(npydir + "/CPU_P2_rms_unalt.npy", CPU_rms_unalt[:,:,:,2])

        np.save(npydir + "/CPU_A_avg_imp.npy",    CPU_avg_imp[:,:,:,3])
        np.save(npydir + "/CPU_A_avg_unalt.npy",  CPU_avg_unalt[:,:,:,3])
        np.save(npydir + "/CPU_A_rms_imp.npy",    CPU_rms_imp[:,:,:,3])
        np.save(npydir + "/CPU_A_rms_unalt.npy",  CPU_rms_unalt[:,:,:,3])

