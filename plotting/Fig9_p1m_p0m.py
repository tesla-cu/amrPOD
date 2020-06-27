
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import interpolate

from scipy.interpolate import RegularGridInterpolator as rgi

def Fig9_p1m_p0m(datadir, imgdir):

    print('making figure 9 ...')

    # Set up figure
    fig = plt.figure()

    # Load data
    txtdir = datadir + 'lc0_lc1/txt_files/'

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

    # Data information
    lc0min = 0/64  # x min
    lc0max = 32/64 # x max
    lc0inc = 1/64  # x increment
    lc1min = 0/64  # y min
    lc1max = 32/64 # y max
    lc1inc = 1/64  # y increment

    lc0 = np.arange(lc0min, lc0max+lc0inc, lc0inc)
    lc1 = np.arange(lc1min, lc1max+lc1inc, lc1inc)

    Lc0, Lc1 = np.meshgrid(lc0, lc1, indexing='ij')


    # Plot figure -------------------------------------------------------------
    grid = AxesGrid(fig, (0.08,0.13,0.81,0.80),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = True,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")
    

    for i in range(4):

        ax = grid[i].axes

        if i == 0:
            num = TC_R_avg_imp
            den = TC_R_avg_unalt
            ax.set_title(r'$\mathbf{R}$')
        elif i == 1:
            num = TC_P1_avg_imp
            den = TC_P1_avg_unalt
            ax.set_title(r'$\mathbf{\Phi}$ -- Method 1')
        elif i == 2:
            num = TC_P2_avg_imp
            den = TC_P2_avg_unalt
            ax.set_title(r'$\mathbf{\Phi}$ -- Method 2')
        elif i == 3:
            num = TC_A_avg_imp
            den = TC_A_avg_unalt
            ax.set_title(r'$\mathbf{A}$')

        # Get ratio of AMR over standard operations
        ratio = np.log10(num/den)
        print(np.max(ratio))
        print(np.min(ratio))

        # Display contoured image
        im = ax.contourf(Lc0, Lc1, ratio, 100, origin='lower', \
            extent=[lc0min,lc0max,lc1min,lc1max], cmap='bwr', \
            vmin=-0.09, vmax=0.09)

        # Colorbar information
        cbar = ax.cax.colorbar(im, format='%0.0e')
        ax.cax.set_ylabel(r'$\log(\overline{T}_a/T_s)$ [ops]')
        cbar.ax.set_ylim(-0.09,0.09)
        cbar.ax.set_yticks(np.linspace(-0.08,0.08,5))
        cbar.ax.set_yticklabels(['-0.08','-0.04','0.00','0.04','0.08'])

        ax.set_xticks(np.linspace(1/8,3/8,3))
        ax.set_yticks(np.linspace(1/8,3/8,3))
        ax.set_xticklabels(['1/8','2/8','3/8'])
        ax.set_yticklabels(['1/8','2/8','3/8'])
        ax.set_xlabel(r'$\tilde{p}_0^m$')
        ax.set_ylabel(r'$\tilde{p}_1^m$')

    # Save image
    fig.set_size_inches(6.5,1.8,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'Fig9_p1m_p0m.png', dpi=600)



    print('\tdone with figure 9')