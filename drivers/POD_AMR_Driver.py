import numpy as np
import matplotlib.pyplot as plt 
import multiprocessing as mp
import os
import sys
from numpy import linalg as LA
from Compute_POD import Compute_POD
from Plot_Varying_rc_nt import Plot_Varying_rc_nt

def parallel_run(args):
   return Compute_POD(*args)

if __name__ == '__main__':

        generate_data              = True
        nx                         = 32                                                        
        ny                         = 32
        nz                         = 1
        finest                     = 2
        nsample                    = 2
        nlev                       = finest+1
        nt_arr                     = np.arange(10, 31, 10)        
        ndim_driver                = 2
        d_l                        = np.zeros((nlev), dtype=int)
        for i in range(nlev):
            d_l[i]    = (2**ndim_driver)**(finest-i) 

        rc_arr       = np.arange(12,1,-2)
        n            = np.size(rc_arr)
        l_frac_0     = .9
        tol          = 10000
        lc_fracs     = np.array([1/16, 0/16, 0/16])
        l0_frac_arr  = np.zeros(np.size(rc_arr))
        l1_frac_arr  = np.zeros(np.size(rc_arr))
        l2_frac_arr  = np.zeros(np.size(rc_arr))
        l_frac_data  = np.zeros((np.size(l1_frac_arr), nlev))

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

        
        basedir = os.getcwd()        
        datadir = basedir + ('/Desktop')
        if not os.path.exists(datadir):
            os.mkdir(datadir)

        studydir = datadir + ('/Varying_rc_Exp/')   # Folder containing all data from this case
        if not os.path.exists(studydir):
            os.mkdir(studydir)

        txtdir = studydir + ('/Data/')              # Folder containing all .txt files produced
        if not os.path.exists(txtdir):
            os.mkdir(txtdir)        

        TC_R_imp_avg       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_R_unalt_avg     = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_R_imp_rms       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_R_unalt_rms     = np.zeros((np.size(rc_arr), np.size(nt_arr)))

        TC_Phi_1_imp_avg   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_1_unalt_avg = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_1_imp_rms   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_1_unalt_rms = np.zeros((np.size(rc_arr), np.size(nt_arr)))

        TC_Phi_5_imp_avg   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_5_unalt_avg = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_5_imp_rms   = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_Phi_5_unalt_rms = np.zeros((np.size(rc_arr), np.size(nt_arr)))

        TC_A_imp_avg       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_A_unalt_avg     = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_A_imp_rms       = np.zeros((np.size(rc_arr), np.size(nt_arr)))
        TC_A_unalt_rms     = np.zeros((np.size(rc_arr), np.size(nt_arr)))

        # Start parallel processing
        nthread = mp.cpu_count()
        nthread = 2
        print('starting pool with %i threads ...' % nthread)
        pool = mp.Pool(processes=nthread)

        for int_arr, nt in enumerate(nt_arr):
                for irc_arr, rc in enumerate(rc_arr):
                        l_fracs       = l_frac_data[irc_arr,:]                        
                        print('lc_fracs_Driver = ', lc_fracs)
                        
                        # Get list of tuples length of nsamp 
                        arg_tuple = (generate_data, nx, ny, nz, finest, l_fracs, lc_fracs, nt)
                        arg_list = []
                        [arg_list.append(arg_tuple) for i in range(nsample)]

                        # Farm out nsample to each processor
                        res_tuple = pool.map(parallel_run, arg_list)
                        print('res_tuple = ', res_tuple)
                        print('comp ratio = ', rc_arr[irc_arr])

                        for i in range(8):
                                res = np.zeros((nsample), dtype=int)

                                for j in range(nsample):
                                        res[j] = res_tuple[j][i]

                                res_avg = np.mean(res)
                                res_rms = np.std(res-res_avg)
                                print('res_avg = ', res_avg)
                                print('res_rms =', res_rms)

                                if i==0: 
                                        TC_R_imp_avg[irc_arr, int_arr]       = res_avg
                                        TC_R_imp_rms[irc_arr, int_arr]       = res_rms
                                elif i==1: 
                                        TC_R_unalt_avg[irc_arr, int_arr]     = res_avg
                                        TC_R_unalt_rms[irc_arr, int_arr]     = res_rms
                                elif i==2: 
                                        TC_Phi_1_imp_avg[irc_arr, int_arr]   = res_avg
                                        TC_Phi_1_imp_rms[irc_arr, int_arr]   = res_rms
                                elif i==3: 
                                        TC_Phi_1_unalt_avg[irc_arr, int_arr] = res_avg
                                        TC_Phi_1_unalt_rms[irc_arr, int_arr] = res_rms
                                elif i==4: 
                                        TC_Phi_5_imp_avg[irc_arr, int_arr]   = res_avg
                                        TC_Phi_5_imp_rms[irc_arr, int_arr]   = res_rms
                                elif i==5: 
                                        TC_Phi_5_unalt_avg[irc_arr, int_arr] = res_avg
                                        TC_Phi_5_unalt_rms[irc_arr, int_arr] = res_rms
                                elif i==6: 
                                        TC_A_imp_avg[irc_arr, int_arr]       = res_avg
                                        TC_A_imp_rms[irc_arr, int_arr]       = res_rms
                                elif i==7: 
                                        TC_A_unalt_avg[irc_arr, int_arr]     = res_avg
                                        TC_A_unalt_rms[irc_arr, int_arr]     = res_rms

        pool.close()
        pool.join()

        #-------- Save data
        np.savetxt(txtdir + "/TC_R_imp_avg.txt", TC_R_imp_avg)
        np.savetxt(txtdir + "/TC_R_unalt_avg.txt", TC_R_unalt_avg)
        np.savetxt(txtdir + "/TC_R_imp_rms.txt", TC_R_imp_rms)
        np.savetxt(txtdir + "/TC_R_unalt_rms.txt", TC_R_unalt_rms)

        np.savetxt(txtdir + "/TC_Phi_1_imp_avg.txt", TC_Phi_1_imp_avg)
        np.savetxt(txtdir + "/TC_Phi_1_unalt_avg.txt", TC_Phi_1_unalt_avg)
        np.savetxt(txtdir + "/TC_Phi_1_imp_rms.txt", TC_Phi_1_imp_rms)
        np.savetxt(txtdir + "/TC_Phi_1_unalt_rms.txt", TC_Phi_1_unalt_rms)

        np.savetxt(txtdir + "/TC_Phi_5_imp_avg.txt", TC_Phi_5_imp_avg)
        np.savetxt(txtdir + "/TC_Phi_5_unalt_avg.txt", TC_Phi_5_unalt_avg)        
        np.savetxt(txtdir + "/TC_Phi_5_imp_rms.txt", TC_Phi_5_imp_rms)
        np.savetxt(txtdir + "/TC_Phi_5_unalt_rms.txt", TC_Phi_5_unalt_rms)

        np.savetxt(txtdir + "/TC_A_imp_avg.txt", TC_A_imp_avg)
        np.savetxt(txtdir + "/TC_A_unalt_avg.txt", TC_A_unalt_avg)
        np.savetxt(txtdir + "/TC_A_imp_rms.txt", TC_A_imp_rms)
        np.savetxt(txtdir + "/TC_A_unalt_rms.txt", TC_A_unalt_rms)

        np.savetxt(txtdir + "/rc_array.txt", rc_arr)
        np.savetxt(txtdir + "/nt_array.txt", nt_arr)

        sim_info = open(txtdir + '/sim_info.txt', 'w')
        sim_info.write("Data Generated: %s\n" % generate_data)
        sim_info.write("finest: %i\n" % finest)
        sim_info.write("samples: %i\n" % nsample)
        sim_info.write("nx: %i\n" % nx)
        sim_info.write("ny: %i\n" % ny)
        sim_info.write("nz: %i\n" % nz)
        sim_info.write("nt_0: %i\n" % nt_arr[0])
        sim_info.write("nt_end: %i\n" % nt_arr[-1])

        
        # ------ Plot Results 
        plot_res = 600                      #image resolution in dpi

        plotdir  = studydir + ('Plots/')    # Folder containing all .png images 
        if not os.path.exists(plotdir):
            os.mkdir(plotdir)

        fig_R, fig_R_imp_rms, fig_R_unalt_rms, fig_Phi_1, fig_Phi_1_imp_rms, fig_Phi_1_unalt_rms, fig_Phi_5, fig_Phi_5_imp_rms, \
        fig_Phi_5_unalt_rms, fig_A, fig_A_imp_rms, fig_A_unalt_rms = Plot_Varying_rc_nt(txtdir)

        fig_R.savefig(plotdir + 'R_efficiency.png', dpi = plot_res)
        fig_R_imp_rms.savefig(plotdir + 'R_imp_rms.png', dpi = plot_res)
        fig_R_unalt_rms.savefig(plotdir + 'R_unalt_rms.png', dpi = plot_res)

        fig_Phi_1.savefig(plotdir + 'Phi_1_efficiency.png', dpi = plot_res)
        fig_Phi_1_imp_rms.savefig(plotdir + 'Phi_1_imp_rms.png', dpi = plot_res)
        fig_Phi_1_unalt_rms.savefig(plotdir + 'Phi_1_unalt_rms.png', dpi = plot_res)

        fig_Phi_5.savefig(plotdir + 'Phi_5_efficiency.png', dpi = plot_res)
        fig_Phi_5_imp_rms.savefig(plotdir + 'Phi_5_imp_rms.png', dpi = plot_res)
        fig_Phi_5_unalt_rms.savefig(plotdir + 'Phi_5_unalt_rms.png', dpi = plot_res)

        fig_A.savefig(plotdir + 'A_efficiency.png', dpi = plot_res)
        fig_A_imp_rms.savefig(plotdir + 'A_imp_rms.png', dpi = plot_res)
        fig_A_unalt_rms.savefig(plotdir + 'A_unalt_rms.png', dpi = plot_res)



