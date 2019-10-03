import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import os
import sys

def gen_fig(num, den, imgdir, figfile, xaxis, yaxis, xlabel, ylabel):
	# num - numerator of images
	# den - denomenator for normalization

	img = np.transpose(np.divide(num, den))
	# img = np.divide(num, den)
	Xaxis, Yaxis = np.meshgrid(xaxis, yaxis, copy=False, indexing='ij')
	cont_levs = 100
	xmin = np.min(xaxis)
	xmax = np.max(xaxis)
	ymin = np.min(yaxis)
	ymax = np.max(yaxis)


	fig = plt.figure(clear=True)
	ax = fig.add_subplot(1,1,1)

	# im = ax.contourf(Xaxis, Yaxis, img, cont_levs, origin='lower', \
	# 	extent=[xmin,xmax,ymin,ymax], cmap='bwr')
	im = ax.imshow(img, aspect='auto', origin='lower', \
		extent=[xmin,xmax,ymin,ymax], cmap='bwr')
	cbar = ax.figure.colorbar(im, ax=ax)

	ax.set_xticks(xaxis)
	ax.set_yticks(yaxis)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	# cbar = plt.colorbar()
	# cbar.set_label('${T_{imp}}/{T_{unalt}}$')
	# plt.xticks(xaxis)
	# plt.yticks(yaxis)
	# plt.xlabel(xlabel)
	# plt.ylabel(ylabel)
	# plt.title(figfile)

	# ax.set_xticks(xaxis)
	# ax.set_yticks(yaxis)

	# ax.set_xticks(xaxis-2.5, minor=True)
	# ax.set_yticks(yaxis-.03125, minor=True)

	# ax.xaxis.set_major_locator(plt.NullLocator())
	# ax.xaxis.set_major_formatter(plt.NullFormatter())
	# ax.yaxis.set_major_locator(plt.NullLocator())
	# ax.yaxis.set_major_formatter(plt.NullFormatter())

	

	# ax.set_yticks([0.5,1.5], minor=True)

	# ax.set_xticklabels(['10','15','20'], minor=True)
	# ax.set_yticklabels(['0','1'])

    # Turn spines off and create white grid.
	# for edge, spine in ax.spines.items():
	# 	spine.set_visible(False)

	

	# plt.setp(ax.get_xticklabels(), ha="right")

	# Hide major tick labels
	# ax.set_xticklabels('')
	# ax.set_yticklabels('')
	# ax.tick_params(axis='both', which	='minor')

	# # Customize minor tick labels
	# ax.set_xticks(xaxis + (np.diff(xaxis[0:2])/2.0),      minor=True)
	# ax.set_xticklabels(['10','15','20'], minor=True)


	# # Hide major tick labels
	# ax.xaxis.set_major_formatter(ticker.NullFormatter())

	# # Customize minor tick labels
	# ax.xaxis.set_minor_locator(ticker.FixedLocator(xaxis))
	# ax.xaxis.set_minor_formatter(ticker.FixedFormatter(['10','15','20']))



	# ax.tick_params(ha='right')
	# for label in ax.get_xticklabels():
	# 	label.set_horizontalalignment('right')
	# plt.set_cmap('bwr')
	plt.savefig(imgdir + figfile, dpi=300)
	plt.close()





def plot_single(txtdir, imgdir, compute_tc, compute_cpu, xaxis, yaxis, xlabel, ylabel):

	# ========== Load data
	if compute_tc:
		TC_R_avg_imp   = np.loadtxt(txtdir + "/TC_R_avg_imp.txt")
		TC_R_avg_unalt = np.loadtxt(txtdir + "/TC_R_avg_unalt.txt")
		TC_R_rms_imp   = np.loadtxt(txtdir + "/TC_R_rms_imp.txt")
		TC_R_rms_unalt = np.loadtxt(txtdir + "/TC_R_rms_unalt.txt")

		TC_P1_avg_imp   = np.loadtxt(txtdir + "/TC_P1_avg_imp.txt")
		TC_P1_avg_unalt = np.loadtxt(txtdir + "/TC_P1_avg_unalt.txt")
		TC_P1_rms_imp   = np.loadtxt(txtdir + "/TC_P1_rms_imp.txt")
		TC_P1_rms_unalt = np.loadtxt(txtdir + "/TC_P1_rms_unalt.txt")

		TC_P2_avg_imp   = np.loadtxt(txtdir + "/TC_P2_avg_imp.txt")
		TC_P2_avg_unalt = np.loadtxt(txtdir + "/TC_P2_avg_unalt.txt")
		TC_P2_rms_imp   = np.loadtxt(txtdir + "/TC_P2_rms_imp.txt")
		TC_P2_rms_unalt = np.loadtxt(txtdir + "/TC_P2_rms_unalt.txt")

		TC_A_avg_imp   = np.loadtxt(txtdir + "/TC_A_avg_imp.txt")
		TC_A_avg_unalt = np.loadtxt(txtdir + "/TC_A_avg_unalt.txt")
		TC_A_rms_imp   = np.loadtxt(txtdir + "/TC_A_rms_imp.txt")
		TC_A_rms_unalt = np.loadtxt(txtdir + "/TC_A_rms_unalt.txt")

	if compute_cpu:
		CPU_R_avg_imp   = np.loadtxt(txtdir + "/CPU_R_avg_imp.txt")
		CPU_R_avg_unalt = np.loadtxt(txtdir + "/CPU_R_avg_unalt.txt")
		CPU_R_rms_imp   = np.loadtxt(txtdir + "/CPU_R_rms_imp.txt")
		CPU_R_rms_unalt = np.loadtxt(txtdir + "/CPU_R_rms_unalt.txt")

		CPU_P1_avg_imp   = np.loadtxt(txtdir + "/CPU_P1_avg_imp.txt")
		CPU_P1_avg_unalt = np.loadtxt(txtdir + "/CPU_P1_avg_unalt.txt")
		CPU_P1_rms_imp   = np.loadtxt(txtdir + "/CPU_P1_rms_imp.txt")
		CPU_P1_rms_unalt = np.loadtxt(txtdir + "/CPU_P1_rms_unalt.txt")

		CPU_P2_avg_imp   = np.loadtxt(txtdir + "/CPU_P2_avg_imp.txt")
		CPU_P2_avg_unalt = np.loadtxt(txtdir + "/CPU_P2_avg_unalt.txt")
		CPU_P2_rms_imp   = np.loadtxt(txtdir + "/CPU_P2_rms_imp.txt")
		CPU_P2_rms_unalt = np.loadtxt(txtdir + "/CPU_P2_rms_unalt.txt")

		CPU_A_avg_imp   = np.loadtxt(txtdir + "/CPU_A_avg_imp.txt")
		CPU_A_avg_unalt = np.loadtxt(txtdir + "/CPU_A_avg_unalt.txt")
		CPU_A_rms_imp   = np.loadtxt(txtdir + "/CPU_A_rms_imp.txt")
		CPU_A_rms_unalt = np.loadtxt(txtdir + "/CPU_A_rms_unalt.txt")

    # ========== Plot data

    # ---------- Plot R
	if compute_tc:
		gen_fig(TC_R_avg_imp,  TC_R_avg_unalt,  imgdir, 'TC_R_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(TC_R_rms_imp,  TC_R_avg_imp,    imgdir, 'TC_R_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)
	if compute_cpu:
		gen_fig(CPU_R_avg_imp, CPU_R_avg_unalt, imgdir, 'CPU_R_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(CPU_R_rms_imp, CPU_R_avg_imp,   imgdir, 'CPU_R_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)

	# ---------- Plot Phi 1
	if compute_tc:
		gen_fig(TC_P1_avg_imp,  TC_P1_avg_unalt,  imgdir, 'TC_P1_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(TC_P1_rms_imp,  TC_P1_avg_imp,    imgdir, 'TC_P1_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)
	if compute_cpu:
		gen_fig(CPU_P1_avg_imp, CPU_P1_avg_unalt, imgdir, 'CPU_P1_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(CPU_P1_rms_imp, CPU_P1_avg_imp,   imgdir, 'CPU_P1_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)

	# ---------- Plot Phi 2
	if compute_tc:
		gen_fig(TC_P2_avg_imp,  TC_P2_avg_unalt,  imgdir, 'TC_P2_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(TC_P2_rms_imp,  TC_P2_avg_imp,    imgdir, 'TC_P2_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)
	if compute_cpu:
		gen_fig(CPU_P2_avg_imp, CPU_P2_avg_unalt, imgdir, 'CPU_P2_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(CPU_P2_rms_imp, CPU_P2_avg_imp,   imgdir, 'CPU_P2_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)

	# ---------- Plot A
	if compute_tc:
		gen_fig(TC_A_avg_imp,  TC_A_avg_unalt,  imgdir, 'TC_A_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(TC_A_rms_imp,  TC_A_avg_imp,    imgdir, 'TC_A_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)
	if compute_cpu:
		gen_fig(CPU_A_avg_imp, CPU_A_avg_unalt, imgdir, 'CPU_A_imp_unalt', \
			xaxis, yaxis, xlabel, ylabel)
		gen_fig(CPU_A_rms_imp, CPU_A_avg_imp,   imgdir, 'CPU_A_rms_avg', \
			xaxis, yaxis, xlabel, ylabel)


	"""

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

	"""








