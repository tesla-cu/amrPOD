import numpy as np
import matplotlib.pyplot as plt 
import multiprocessing as mp
import os
import sys
from numpy import linalg as LA
from Compute_POD import Compute_POD
from Plot_Varying_level_comp import Plot_Varying_level_comp

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
        nt                         = 25  
        ndim_driver                = 2
        d_l                        = np.zeros((nlev), dtype=int)
        
        for i in range(nlev):
            d_l[i]    = (2**ndim_driver)**(finest-i) 
        
        l_frac_0     = np.arange(.25, .5, .025)
        l_frac_1     = np.arange(.25, .5, .025)
        rc_arr       = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))

        lc_fracs_arr             = np.zeros((4, 3))
        lc_fracs_arr[1,0]        = 4/16
        lc_fracs_arr[2,1]        = 4/16
        lc_fracs_arr[3,0:2]      = 4/16

        basedir = os.getcwd()        
        datadir = basedir + ('/Desktop')
        if not os.path.exists(datadir):
            os.mkdir(datadir)

        studydir = datadir + ('/Varying_level_comp_Exp/')   # Folder containing all data from this case
        if not os.path.exists(studydir):
            os.mkdir(studydir)

        txtdir = studydir + ('/Data/')              # Folder containing all .txt files produced
        if not os.path.exists(txtdir):
            os.mkdir(txtdir)    
        
        for k in range(len(lc_fracs_arr)):

                lc_fracs = lc_fracs_arr[k,:]

                TC_R_imp_avg       = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_R_unalt_avg     = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_R_imp_rms       = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_R_unalt_rms     = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))

                TC_Phi_1_imp_avg   = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_1_unalt_avg = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_1_imp_rms   = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_1_unalt_rms = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))

                TC_Phi_5_imp_avg   = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_5_unalt_avg = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_5_imp_rms   = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_Phi_5_unalt_rms = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))

                TC_A_imp_avg       = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_A_unalt_avg     = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_A_imp_rms       = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))
                TC_A_unalt_rms     = np.zeros((np.size(l_frac_0), np.size(l_frac_1)))

                # Start parallel processing
                nthread = mp.cpu_count()
                nthread = 2
                print('starting pool with %i threads ...' % nthread)
                pool = mp.Pool(processes=nthread)

                for il_frac_0, var_frac_0 in enumerate(l_frac_0):
                        for il_frac_1, var_frac_1 in enumerate(l_frac_1):

                                var_frac_2    = 1 - var_frac_0 - var_frac_1    
                                l_fracs       = np.array([var_frac_0, var_frac_1, var_frac_2])

                                # compute current compression ratio
                                rc_arr[il_frac_0, il_frac_1] = 1/np.sum(np.divide(l_fracs, d_l))

                                # Get list of tuples length of nsamp 
                                arg_tuple = (generate_data, nx, ny, nz, finest, l_fracs, lc_fracs, nt)
                                arg_list = []
                                [arg_list.append(arg_tuple) for i in range(nsample)]

                                # Farm out nsample to each processor
                                res_tuple = pool.map(parallel_run, arg_list)
                                
                                print('res_tuple = ', res_tuple)
                                print('comp ratio = ', rc_arr[il_frac_0, il_frac_1])
                                print('lc_fracs_Driver = ', lc_fracs)

                                for i in range(8):
                                        res = np.zeros((nsample), dtype=int)
                                        for j in range(nsample):
                                                res[j] = res_tuple[j][i]
                                        res_avg = np.mean(res)
                                        res_rms = np.std(res-res_avg)
                                        print('res_avg = ', res_avg)
                                        print('res_rms =', res_rms)

                                        if i==0: 
                                                TC_R_imp_avg[il_frac_0, il_frac_1]       = res_avg
                                                TC_R_imp_rms[il_frac_0, il_frac_1]       = res_rms
                                        elif i==1: 
                                                TC_R_unalt_avg[il_frac_0, il_frac_1]     = res_avg
                                                TC_R_unalt_rms[il_frac_0, il_frac_1]     = res_rms
                                        elif i==2: 
                                                TC_Phi_1_imp_avg[il_frac_0, il_frac_1]   = res_avg
                                                TC_Phi_1_imp_rms[il_frac_0, il_frac_1]   = res_rms
                                        elif i==3: 
                                                TC_Phi_1_unalt_avg[il_frac_0, il_frac_1] = res_avg
                                                TC_Phi_1_unalt_rms[il_frac_0, il_frac_1] = res_rms
                                        elif i==4: 
                                                TC_Phi_5_imp_avg[il_frac_0, il_frac_1]   = res_avg
                                                TC_Phi_5_imp_rms[il_frac_0, il_frac_1]   = res_rms
                                        elif i==5: 
                                                TC_Phi_5_unalt_avg[il_frac_0, il_frac_1] = res_avg
                                                TC_Phi_5_unalt_rms[il_frac_0, il_frac_1] = res_rms
                                        elif i==6: 
                                                TC_A_imp_avg[il_frac_0, il_frac_1]       = res_avg
                                                TC_A_imp_rms[il_frac_0, il_frac_1]       = res_rms
                                        elif i==7: 
                                                TC_A_unalt_avg[il_frac_0, il_frac_1]     = res_avg
                                                TC_A_unalt_rms[il_frac_0, il_frac_1]     = res_rms

                pool.close()
                pool.join()

                #-------- Save data
                np.savetxt(txtdir + "/TC_R_imp_avg" + str(k) + ".txt", TC_R_imp_avg)
                np.savetxt(txtdir + "/TC_R_unalt_avg" + str(k) + ".txt", TC_R_unalt_avg)
                np.savetxt(txtdir + "/TC_R_imp_rms" + str(k) + ".txt", TC_R_imp_rms)
                np.savetxt(txtdir + "/TC_R_unalt_rms" + str(k) + ".txt", TC_R_unalt_rms)

                np.savetxt(txtdir + "/TC_Phi_1_imp_avg" + str(k) + ".txt", TC_Phi_1_imp_avg)
                np.savetxt(txtdir + "/TC_Phi_1_unalt_avg" + str(k) + ".txt", TC_Phi_1_unalt_avg)
                np.savetxt(txtdir + "/TC_Phi_1_imp_rms" + str(k) + ".txt", TC_Phi_1_imp_rms)
                np.savetxt(txtdir + "/TC_Phi_1_unalt_rms" + str(k) + ".txt", TC_Phi_1_unalt_rms)

                np.savetxt(txtdir + "/TC_Phi_5_imp_avg" + str(k) + ".txt", TC_Phi_5_imp_avg)
                np.savetxt(txtdir + "/TC_Phi_5_unalt_avg" + str(k) + ".txt", TC_Phi_5_unalt_avg)        
                np.savetxt(txtdir + "/TC_Phi_5_imp_rms" + str(k) + ".txt", TC_Phi_5_imp_rms)
                np.savetxt(txtdir + "/TC_Phi_5_unalt_rms" + str(k) + ".txt", TC_Phi_5_unalt_rms)

                np.savetxt(txtdir + "/TC_A_imp_avg" + str(k) + ".txt", TC_A_imp_avg)
                np.savetxt(txtdir + "/TC_A_unalt_avg" + str(k) + ".txt", TC_A_unalt_avg)
                np.savetxt(txtdir + "/TC_A_imp_rms" + str(k) + ".txt", TC_A_imp_rms)
                np.savetxt(txtdir + "/TC_A_unalt_rms" + str(k) + ".txt", TC_A_unalt_rms)
                np.savetxt(txtdir + "/rc_array" + str(k) + ".txt", rc_arr)

                np.savetxt(txtdir + "/l_frac_0_arr" + str(k) + ".txt", l_frac_0)
                np.savetxt(txtdir + "/l_frac_1_arr" + str(k) + ".txt", l_frac_1)
                np.savetxt(txtdir + "/l_c_fracs" + str(k) + ".txt", lc_fracs_arr[k,:])


                sim_info = open(txtdir + '/sim_info.txt', 'w')
                sim_info.write("Data Generated: %s\n" % generate_data)
                sim_info.write("finest: %i\n" % finest)
                sim_info.write("samples: %i\n" % nsample)
                sim_info.write("nx: %i\n" % nx)
                sim_info.write("ny: %i\n" % ny)
                sim_info.write("nz: %i\n" % nz)
                sim_info.write("nt: %i\n" % nt)

                # ------ Plot Results 
                plot_res = 600                      #image resolution in dpi

                plotdir  = studydir + ('Plots/')    # Folder containing all .png images 
                if not os.path.exists(plotdir):
                    os.mkdir(plotdir)

                fig_R, fig_R_imp_rms, fig_R_unalt_rms, fig_Phi_1, fig_Phi_1_imp_rms, fig_Phi_1_unalt_rms, fig_Phi_5, fig_Phi_5_imp_rms, \
                fig_Phi_5_unalt_rms, fig_A, fig_A_imp_rms, fig_A_unalt_rms = Plot_Varying_level_comp(txtdir, k)

                fig_R.savefig(plotdir + "R_efficiency{k}.png".format(k=k), dpi = plot_res)
                fig_R_imp_rms.savefig(plotdir + "R_imp_rms{k}.png".format(k=k), dpi = plot_res)
                fig_R_unalt_rms.savefig(plotdir + "R_unalt_rms{k}.png".format(k=k), dpi = plot_res)

                fig_Phi_1.savefig(plotdir + "Phi_1_efficiency{k}.png".format(k=k), dpi = plot_res)
                fig_Phi_1_imp_rms.savefig(plotdir + "Phi_1_imp_rms{k}.png".format(k=k), dpi = plot_res)
                fig_Phi_1_unalt_rms.savefig(plotdir + "Phi_1_unalt_rms{k}.png".format(k=k), dpi = plot_res)

                fig_Phi_5.savefig(plotdir + "Phi_5_efficiency{k}.png".format(k=k), dpi = plot_res)
                fig_Phi_5_imp_rms.savefig(plotdir + "Phi_5_imp_rms{k}.png".format(k=k), dpi = plot_res)
                fig_Phi_5_unalt_rms.savefig(plotdir + "Phi_5_unalt_rms{k}.png".format(k=k), dpi = plot_res)

                fig_A.savefig(plotdir + "A_efficiency{k}.png".format(k=k), dpi = plot_res)
                fig_A_imp_rms.savefig(plotdir + "A_imp_rms{k}.png".format(k=k), dpi = plot_res)
                fig_A_unalt_rms.savefig(plotdir + "A_unalt_rms{k}.png".format(k=k), dpi = plot_res)

