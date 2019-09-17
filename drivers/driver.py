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


if __name__ == '__main__':

    print('starting script to perform POD on AMR grids ...')

    # ---------- User defined inputs -------------------------------- 
    gen_grid    = True  # are we generating synthetic data?
    compute_tc  = False  # are we computing the time complexity?
    compute_cpu = True # are we computing the cpu time?
    nx          = 32    # x spatial points                  
    ny          = 1    # y spatial points
    nz          = 1     # z spatial points
    finest      = 4     # finest level of AMR in the domain
    nsample     = 1     # number of samples for each parameter set
    # nt_arr      = np.arange(101, 102, 5)     # spanning nt
    nt_arr      = np.arange(5, 6, 5)     # spanning nt
    l1_arr      = np.arange(4*0.0625, 7*0.0625, .25) # spanning l1
    lcs         = np.zeros((finest+1)) # fraction of grid that stays constant in time
    # lc_fracs     = np.array([1/16, 0/16, 0/16])

    # Variables we will iterate through
    xvals = nt_arr # what will be plotted on the x-axis
    yvals = l1_arr # what will be plotted on the y-axis

    
    # lc_fracs     = np.array([1/16, 0/16, 0/16])
    # l0_frac_arr  = np.zeros(np.size(rc_arr))
    # l1_frac_arr  = np.zeros(np.size(rc_arr))
    # l2_frac_arr  = np.zeros(np.size(rc_arr))
    # l_frac_data  = np.zeros((np.size(l1_frac_arr), nlev))

    # Direction where /code/ lives
    basedir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/'

    # Directory where AMR data is stored
    amr_datadir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/' + \
        'slice/x0.000-0.000_y-1.000-1.000_z0.000-2.000_t40.0000-42.0000/'

    # Directory where we was to store data on speed up
    datadir = basedir + 'data/'
    if not os.path.exists(datadir):
        os.mkdir(datadir)

    # Directory that describes the study we are looking at
    studydir = datadir + 'nt_l1_small/'
    if not os.path.exists(studydir):
        os.mkdir(studydir)

    # Directory that holds all txt files of TC and CPU
    txtdir = studydir + 'txt_files/'
    if not os.path.exists(txtdir):
        os.mkdir(txtdir)

    # Directory that holds all images generated from txt files
    imgdir = studydir + 'images/'
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)

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
    #    - dim 1: what will be plotted on the xaxis
    #    - dim 2: what will be plotted on the yaxis
    #    - dim 3: the different algorithms
    #        - index 0: R
    #        - index 1: Phi, method 1
    #        - index 2: Phi, method 2
    #        - index 3: A
    #    - avg is the average of all the samples
    #    - rms is the rms of all the samples
    if compute_tc:
        TC_avg_imp   = np.zeros((len(xvals), len(yvals), 4), dtype=int)
        TC_avg_unalt = np.zeros((len(xvals), len(yvals), 4), dtype=int)
        TC_rms_imp   = np.zeros((len(xvals), len(yvals), 4), dtype=int)
        TC_rms_unalt = np.zeros((len(xvals), len(yvals), 4), dtype=int)

    if compute_cpu:
        CPU_avg_imp   = np.zeros((len(xvals), len(yvals), 4))
        CPU_avg_unalt = np.zeros((len(xvals), len(yvals), 4))
        CPU_rms_imp   = np.zeros((len(xvals), len(yvals), 4))
        CPU_rms_unalt = np.zeros((len(xvals), len(yvals), 4))

    # TC_R_imp_avg       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_R_unalt_avg     = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_R_imp_rms       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_R_unalt_rms     = np.zeros((np.size(rc_arr), np.size(nt_arr)))

    # TC_Phi_1_imp_avg   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_1_unalt_avg = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_1_imp_rms   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_1_unalt_rms = np.zeros((np.size(rc_arr), np.size(nt_arr)))

    # TC_Phi_5_imp_avg   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_5_unalt_avg = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_5_imp_rms   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_Phi_5_unalt_rms = np.zeros((np.size(rc_arr), np.size(nt_arr)))

    # TC_A_imp_avg       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_A_unalt_avg     = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_A_imp_rms       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
    # TC_A_unalt_rms     = np.zeros((np.size(rc_arr), np.size(nt_arr)))

    # Start parallel processing
    # nthread = mp.cpu_count()
    nthread = 1
    print('starting pool with %i threads ...' % nthread)
    pool = mp.Pool(processes=nthread)

    for ixval, xval in enumerate(xvals):
        for iyval, yval in enumerate(yvals):
            nt = xval # less confusing
            l1 = yval # less confusing

            print('nt = %i, l1 = %0.4f' % (nt, l1))

            # Since we have l in our vars, we need to ls
            # ls = np.array([1.0-l1, l1])
            # ls = np.array([1.0-3*l1, l1, l1, l1/2, l1/2])
            # ls = np.array([1.0-3*l1, l1, l1, l1])

            ls = np.array([0.5, 0.25, 0.125, 0.0625, 0.0625])
            # ls = np.array([0.5, 0.25, 0.125, 0.125])

            # ls = np.array([1.0-l1-0.125-0.0625, l1, 0.0625, 0.0625, 0.0625])
            # l_fracs       = l_frac_data[irc_arr,:] 

            # ========== Compute time complexity
            if compute_tc:
                # Generate tuple that stores all relevant parameters
                #    - note, we need list of these equal to number of samples
                if gen_grid:
                    arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt)
                else:
                    arg_tuple = (gen_grid, nx, ny, nz, finest, ls, lcs, nt, \
                        'TC', amr_datadir)
                arg_list = []
                [arg_list.append(arg_tuple) for i in range(nsample)]
                # print(arg_list)

                # Farm out nsample to each processor
                res_tuple = pool.map(parallel_run, arg_list)
                # print('res_tuple = ', res_tuple)
                # print('comp ratio = ', rc_arr[irc_arr])

                for i in range(8):
                    res = np.zeros((nsample), dtype=int)

                    for j in range(nsample):
                            res[j] = res_tuple[j][i]

                    res_avg = np.mean(res)
                    res_rms = np.std(res-res_avg)
                    # print('res_avg = ', res_avg)
                    # print('res_rms =', res_rms)

                    if   i==0: 
                        TC_avg_imp[ixval, iyval, 0]   = res_avg
                        TC_rms_imp[ixval, iyval, 0]   = res_rms
                    elif i==1: 
                        TC_avg_unalt[ixval, iyval, 0] = res_avg
                        TC_rms_unalt[ixval, iyval, 0] = res_rms
                    elif i==2: 
                        TC_avg_imp[ixval, iyval, 1]   = res_avg
                        TC_rms_imp[ixval, iyval, 1]   = res_rms
                    elif i==3: 
                        TC_avg_unalt[ixval, iyval, 1] = res_avg
                        TC_rms_unalt[ixval, iyval, 1] = res_rms
                    elif i==4: 
                        TC_avg_imp[ixval, iyval, 2]   = res_avg
                        TC_rms_imp[ixval, iyval, 2]   = res_rms
                    elif i==5: 
                        TC_avg_unalt[ixval, iyval, 2] = res_avg
                        TC_rms_unalt[ixval, iyval, 2] = res_rms
                    elif i==6: 
                        TC_avg_imp[ixval, iyval, 3]   = res_avg
                        TC_rms_imp[ixval, iyval, 3]   = res_rms
                    elif i==7: 
                        TC_avg_unalt[ixval, iyval, 3] = res_avg
                        TC_rms_unalt[ixval, iyval, 3] = res_rms

            # ========== Compute CPU time
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
                # print(arg_list)

                # Farm out nsample to each processor
                res_tuple = pool.map(parallel_run, arg_list)
                # print('res_tuple = ', res_tuple)
                # print('comp ratio = ', rc_arr[irc_arr])

                for i in range(8):
                    res = np.zeros((nsample))

                    for j in range(nsample):
                            res[j] = res_tuple[j][i]

                    res_avg = np.mean(res)
                    res_rms = np.std(res-res_avg)
                    # print('res_avg = ', res_avg)
                    # print('res_rms =', res_rms)

                    if   i==0: 
                        CPU_avg_imp[ixval, iyval, 0]   = res_avg
                        CPU_rms_imp[ixval, iyval, 0]   = res_rms
                    elif i==1: 
                        CPU_avg_unalt[ixval, iyval, 0] = res_avg
                        CPU_rms_unalt[ixval, iyval, 0] = res_rms
                    elif i==2: 
                        CPU_avg_imp[ixval, iyval, 1]   = res_avg
                        CPU_rms_imp[ixval, iyval, 1]   = res_rms
                    elif i==3: 
                        CPU_avg_unalt[ixval, iyval, 1] = res_avg
                        CPU_rms_unalt[ixval, iyval, 1] = res_rms
                    elif i==4: 
                        CPU_avg_imp[ixval, iyval, 2]   = res_avg
                        CPU_rms_imp[ixval, iyval, 2]   = res_rms
                    elif i==5: 
                        CPU_avg_unalt[ixval, iyval, 2] = res_avg
                        CPU_rms_unalt[ixval, iyval, 2] = res_rms
                    elif i==6: 
                        CPU_avg_imp[ixval, iyval, 3]   = res_avg
                        CPU_rms_imp[ixval, iyval, 3]   = res_rms
                    elif i==7: 
                        CPU_avg_unalt[ixval, iyval, 3] = res_avg
                        CPU_rms_unalt[ixval, iyval, 3] = res_rms

    pool.close()
    pool.join()

    # ========== Save data

    # Write out data from this run
    sim_info = open(txtdir + '/sim_info.txt', 'w')
    sim_info.write("gen_grid: %s\n" % gen_grid)
    sim_info.write("finest:   %i\n" % finest)
    sim_info.write("nsample:  %i\n" % nsample)
    sim_info.write("nx:       %i\n" % nx)
    sim_info.write("ny:       %i\n" % ny)
    sim_info.write("nz:       %i\n" % nz)
    sim_info.write("xval:     %s\n" % "nt")
    [sim_info.write("lc%i:      %0.8f\n" % (i,lc)) for i,lc in enumerate(lcs)]
    sim_info.write("x_0:      %i\n" % xvals[0])
    sim_info.write("x_inc:    %i\n" % np.diff(xvals[0:2]))
    sim_info.write("x_end:    %i\n" % xvals[-1])
    sim_info.write("yval:     %s\n" % "l1")
    sim_info.write("y_0:      %i\n" % yvals[0])
    sim_info.write("y_inc:    %i\n" % np.diff(yvals[0:2]))
    sim_info.write("y_end:    %i\n" % yvals[-1])
    sim_info.close()

    # np.savetxt(txtdir + "/rc_array.txt", rc_arr)
    # np.savetxt(txtdir + "/nt_array.txt", nt_arr)

    # Write out time complexity data
    if compute_tc:
        np.savetxt(txtdir + "/TC_R_avg_imp.txt",    TC_avg_imp[:,:,0])
        np.savetxt(txtdir + "/TC_R_avg_unalt.txt",  TC_avg_unalt[:,:,0])
        np.savetxt(txtdir + "/TC_R_rms_imp.txt",    TC_rms_imp[:,:,0])
        np.savetxt(txtdir + "/TC_R_rms_unalt.txt",  TC_rms_unalt[:,:,0])

        np.savetxt(txtdir + "/TC_P1_avg_imp.txt",   TC_avg_imp[:,:,1])
        np.savetxt(txtdir + "/TC_P1_avg_unalt.txt", TC_avg_unalt[:,:,1])
        np.savetxt(txtdir + "/TC_P1_rms_imp.txt",   TC_rms_imp[:,:,1])
        np.savetxt(txtdir + "/TC_P1_rms_unalt.txt", TC_rms_unalt[:,:,1])

        np.savetxt(txtdir + "/TC_P2_avg_imp.txt",   TC_avg_imp[:,:,2])
        np.savetxt(txtdir + "/TC_P2_avg_unalt.txt", TC_avg_unalt[:,:,2])
        np.savetxt(txtdir + "/TC_P2_rms_imp.txt",   TC_rms_imp[:,:,2])
        np.savetxt(txtdir + "/TC_P2_rms_unalt.txt", TC_rms_unalt[:,:,2])

        np.savetxt(txtdir + "/TC_A_avg_imp.txt",    TC_avg_imp[:,:,3])
        np.savetxt(txtdir + "/TC_A_avg_unalt.txt",  TC_avg_unalt[:,:,3])
        np.savetxt(txtdir + "/TC_A_rms_imp.txt",    TC_rms_imp[:,:,3])
        np.savetxt(txtdir + "/TC_A_rms_unalt.txt",  TC_rms_unalt[:,:,3])

    # Write out CPU time data
    if compute_cpu:
        np.savetxt(txtdir + "/CPU_R_avg_imp.txt",    CPU_avg_imp[:,:,0])
        np.savetxt(txtdir + "/CPU_R_avg_unalt.txt",  CPU_avg_unalt[:,:,0])
        np.savetxt(txtdir + "/CPU_R_rms_imp.txt",    CPU_rms_imp[:,:,0])
        np.savetxt(txtdir + "/CPU_R_rms_unalt.txt",  CPU_rms_unalt[:,:,0])

        np.savetxt(txtdir + "/CPU_P1_avg_imp.txt",   CPU_avg_imp[:,:,1])
        np.savetxt(txtdir + "/CPU_P1_avg_unalt.txt", CPU_avg_unalt[:,:,1])
        np.savetxt(txtdir + "/CPU_P1_rms_imp.txt",   CPU_rms_imp[:,:,1])
        np.savetxt(txtdir + "/CPU_P1_rms_unalt.txt", CPU_rms_unalt[:,:,1])

        np.savetxt(txtdir + "/CPU_P2_avg_imp.txt",   CPU_avg_imp[:,:,2])
        np.savetxt(txtdir + "/CPU_P2_avg_unalt.txt", CPU_avg_unalt[:,:,2])
        np.savetxt(txtdir + "/CPU_P2_rms_imp.txt",   CPU_rms_imp[:,:,2])
        np.savetxt(txtdir + "/CPU_P2_rms_unalt.txt", CPU_rms_unalt[:,:,2])

        np.savetxt(txtdir + "/CPU_A_avg_imp.txt",    CPU_avg_imp[:,:,3])
        np.savetxt(txtdir + "/CPU_A_avg_unalt.txt",  CPU_avg_unalt[:,:,3])
        np.savetxt(txtdir + "/CPU_A_rms_imp.txt",    CPU_rms_imp[:,:,3])
        np.savetxt(txtdir + "/CPU_A_rms_unalt.txt",  CPU_rms_unalt[:,:,3])

    
    # Plot these various quantities
    plot_single(txtdir, imgdir, compute_tc, compute_cpu, xvals, yvals, 'nt', 'l1')


    # Code if we want to span rc
    """
    rc_arr       = np.arange(12,1,-2)
    n            = np.size(rc_arr)
    l_frac_0     = .9
    tol          = 10000

    #-------- Solve for other level fracs based on specified r_c
    for i in range(n):
        flag = 1
        while flag < tol:
            a       = np.array([ [1, 1], [ 1/d_l[1], 1/d_l[2] ]])
            b       = np.array([1 - l_frac_0, 1/rc_arr[i] - l_frac_0/d_l[0]])
            l1l2    = LA.solve(a,b)
            if l1l2[0] >=0 and l1l2[0] <= 1 and l1l2[1] >= 0 and l1l2[1] <= 1:
                break
            else:
                flag     += 1
                l_frac_0 -= d_l[0]/(nx*ny*nz)

        l1_frac_arr[i] = l1l2[0]
        l2_frac_arr[i] = l1l2[1]
        l0_frac_arr[i] = 1 - l1l2[0] - l1l2[1]
        l_frac_data[i,:] = np.array([l0_frac_arr[i],l1_frac_arr[i],l2_frac_arr[i]])
    """


