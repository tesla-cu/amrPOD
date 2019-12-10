import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import sys

def gen_fig(num, den, imgdir, figfile, xaxis, yaxis, xlabel, ylabel):
	# num - numerator of images
	# den - denomenator for normalization

	# img = np.transpose(np.divide(num, den))
	img = np.divide(num, den)
	Xaxis, Yaxis = np.meshgrid(xaxis, yaxis, copy=False, indexing='ij')
	cont_levs = 100
	xmin = np.min(xaxis)
	xmax = np.max(xaxis)
	ymin = np.min(yaxis)
	ymax = np.max(yaxis)


	fig = plt.figure(clear=True)
	ax = fig.add_subplot(1,1,1)

	im = ax.contourf(Xaxis, Yaxis, img, cont_levs, origin='lower', \
		extent=[xmin,xmax,ymin,ymax], cmap='bwr')
	# im = ax.imshow(img, aspect='auto', origin='lower', \
	# 	extent=[xmin,xmax,ymin,ymax], cmap='bwr')
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








