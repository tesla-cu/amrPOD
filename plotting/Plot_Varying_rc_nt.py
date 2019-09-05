import numpy as np
import matplotlib.pyplot as plt 
import os
import sys

def Plot_Varying_rc_nt(txtdir):

	TC_R_imp_avg       = np.loadtxt(txtdir + "/TC_R_imp_avg.txt")
	TC_R_unalt_avg     = np.loadtxt(txtdir + "/TC_R_unalt_avg.txt")
	TC_R_imp_rms       = np.loadtxt(txtdir + "/TC_R_imp_rms.txt")
	TC_R_unalt_rms     = np.loadtxt(txtdir + "/TC_R_unalt_rms.txt")

	TC_Phi_1_imp_avg   = np.loadtxt(txtdir + "/TC_Phi_1_imp_avg.txt")
	TC_Phi_1_unalt_avg = np.loadtxt(txtdir + "/TC_Phi_1_unalt_avg.txt")
	TC_Phi_1_imp_rms   = np.loadtxt(txtdir + "/TC_Phi_1_imp_rms.txt")
	TC_Phi_1_unalt_rms = np.loadtxt(txtdir + "/TC_Phi_1_unalt_rms.txt")

	TC_Phi_5_imp_avg   = np.loadtxt(txtdir + "/TC_Phi_5_imp_avg.txt")
	TC_Phi_5_unalt_avg = np.loadtxt(txtdir + "/TC_Phi_5_unalt_avg.txt")        
	TC_Phi_5_imp_rms   = np.loadtxt(txtdir + "/TC_Phi_5_imp_rms.txt")
	TC_Phi_5_unalt_rms = np.loadtxt(txtdir + "/TC_Phi_5_unalt_rms.txt")

	TC_A_imp_avg       = np.loadtxt(txtdir + "/TC_A_imp_avg.txt")
	TC_A_unalt_avg	   = np.loadtxt(txtdir + "/TC_A_unalt_avg.txt")
	TC_A_imp_rms       = np.loadtxt(txtdir + "/TC_A_imp_rms.txt")
	TC_A_unalt_rms     = np.loadtxt(txtdir + "/TC_A_unalt_rms.txt")

	rc_arr             = np.loadtxt(txtdir + "/rc_array.txt")
	nt_arr             = np.loadtxt(txtdir + "/nt_array.txt")

	cont_levs          = 100

	# ----------------- Plot R
	fig_R    	       = plt.figure()
	TC_R_Norm          = np.divide(TC_R_imp_avg, TC_R_unalt_avg)
	plt.contourf(nt_arr, rc_arr, TC_R_Norm, cont_levs)
	cbar 			   = plt.colorbar()
	cbar.set_label('${T_{imp}}/{T_{unalt}}$')
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('R Efficiency')
	plt.set_cmap('bwr')

	fig_R_imp_rms      = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_R_imp_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('R Implemented RMS')
	plt.set_cmap('bwr')

	fig_R_unalt_rms     = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_R_unalt_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('R Unaltered RMS')
	plt.set_cmap('bwr')

	# ----------------- Plot Phi 1
	fig_Phi_1	           = plt.figure()
	TC_Phi_1_Norm          = np.divide(TC_Phi_1_imp_avg, TC_Phi_1_unalt_avg)
	plt.contourf(nt_arr, rc_arr, TC_Phi_1_Norm, cont_levs)
	cbar 			       = plt.colorbar()
	cbar.set_label('${T_{imp}}/{T_{unalt}}$')
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{1}$ Efficiency')
	plt.set_cmap('bwr')

	fig_Phi_1_imp_rms      = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_Phi_1_imp_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{1}$ Implemented RMS')
	plt.set_cmap('bwr')

	fig_Phi_1_unalt_rms     = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_Phi_1_unalt_rms, cont_levs)
	plt.colorbar()
	plt.xticks(range(len(nt_arr)), nt_arr)
	plt.yticks(range(len(rc_arr)), rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{1}$ Unaltered RMS')
	plt.set_cmap('bwr')

	# ----------------- Plot Phi 5
	fig_Phi_5	           = plt.figure()
	TC_Phi_5_Norm          = np.divide(TC_Phi_5_imp_avg, TC_Phi_5_unalt_avg)
	plt.contourf(nt_arr, rc_arr, TC_Phi_5_Norm, cont_levs)
	cbar 			       = plt.colorbar()
	cbar.set_label('${T_{imp}}/{T_{unalt}}$')
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{5}$ Efficiency')
	plt.set_cmap('bwr')

	fig_Phi_5_imp_rms      = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_Phi_5_imp_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{5}$ Implemented RMS')
	plt.set_cmap('bwr')

	fig_Phi_5_unalt_rms     = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_Phi_5_unalt_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('$Phi_{1}$ Unaltered RMS')
	plt.set_cmap('bwr')

	# ----------------- Plot A
	fig_A	           = plt.figure()
	TC_A_Norm          = np.divide(TC_A_imp_avg, TC_A_unalt_avg)
	plt.contourf(nt_arr, rc_arr, TC_A_Norm, cont_levs)
	cbar 			   = plt.colorbar()
	cbar.set_label('${T_{imp}}/{T_{unalt}}$')
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('A Efficiency')
	plt.set_cmap('bwr')

	fig_A_imp_rms      = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_A_imp_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('A Implemented RMS')
	plt.set_cmap('bwr')

	fig_A_unalt_rms     = plt.figure()
	plt.contourf(nt_arr, rc_arr, TC_A_unalt_rms, cont_levs)
	plt.colorbar()
	plt.xticks(nt_arr)
	plt.yticks(rc_arr)
	plt.xlabel('$n_{t}$')
	plt.ylabel('$r_{c}$')
	plt.title('A Unaltered RMS')
	plt.set_cmap('bwr')

	plt.close('all')
	
	return fig_R, fig_R_imp_rms, fig_R_unalt_rms, fig_Phi_1, fig_Phi_1_imp_rms, fig_Phi_1_unalt_rms, fig_Phi_5, fig_Phi_5_imp_rms, fig_Phi_5_unalt_rms, fig_A, fig_A_imp_rms, fig_A_unalt_rms










