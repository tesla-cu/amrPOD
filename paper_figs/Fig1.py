import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def Fig1(datadir, imgdir):

    print('making figure 1 ...')

    # Set up figure
    fig = plt.figure()

    # Load data
    txtdir = datadir + 'l1_nt/txt_files/'

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
    ntmin = 2    # x min
    ntmax = 60   # x max
    ntinc = 2    # x increment
    l1min = 0.0  # y min
    l1max = 3/4  # y max
    l1inc = 1/64 # y increment

    nt = np.arange(ntmin, ntmax+ntinc, ntinc)
    l1 = np.arange(l1min, l1max+l1inc, l1inc)

    Nt, L1 = np.meshgrid(nt, l1, indexing='ij')

    # Top half of figure ------------------------------------------------------
    grid = AxesGrid(fig, (0.07,0.53,0.84,0.39),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
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
        # ratio = num/den
        print(np.max(ratio))
        print(np.min(ratio))

        # Display contoured image
        im = ax.contourf(Nt, L1, ratio, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
            vmin=-0.25, vmax=0.25)
        # im = ax.contourf(Nt, L1, ratio, 100, origin='lower', \
        #     extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
        #     vmin=0.5, vmax=1.5)

        # Colorbar information
        cbar = ax.cax.colorbar(im)
        ax.cax.set_ylabel(r'$\log(\overline{T}_a/\overline{T}_s)$ (ops)')
        cbar.ax.set_ylim(-0.25,0.25)
        cbar.ax.set_yticks(np.linspace(-0.2,0.2,5))
        cbar.ax.set_yticklabels(['-0.2','-0.1','0.0','0.1','0.2'])
        # cbar.ax.set_yticklabels([])
        # ax.cax.set_ylabel(r'$\overline{T}_a/\overline{T}_s$ (ops)')
        # cbar.ax.set_ylim(0.5,1.5)
        # cbar.ax.set_yticks(np.linspace(0.6,1.4,5))
        # cbar.ax.set_yticklabels(['0.6','0.8','1.0','1.2','1.4'])

        # Label information
        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(1/8,5/8,3))
        ax.set_xticklabels([])
        ax.set_yticklabels(['1/8','3/8','5/8'])
        ax.set_xlabel([])
        ax.set_ylabel('$p_1$')

    # Bottom half of figure ---------------------------------------------------
    grid = AxesGrid(fig, (0.07,0.11,0.84,0.39),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    for i in range(4):

        ax = grid[i].axes

        if i == 0:
            num = TC_R_rms_imp
            den = TC_R_avg_imp
        elif i == 1:
            num = TC_P1_rms_imp
            den = TC_P1_avg_imp
        elif i == 2:
            num = TC_P2_rms_imp
            den = TC_P2_avg_imp
        elif i == 3:
            num = TC_A_rms_imp
            den = TC_A_avg_imp

        # Get ratio of fluctuations over average
        ratio = (num/den)*1000
        print(np.max(ratio))
        print(np.min(ratio))

        # Display contoured image
        # im = ax.contourf(Nt, L1, ratio, 100, origin='lower', \
        #     extent=[ntmin,ntmax,l1min,l1max], cmap='Reds', \
        #     vmin=0.0, vmax=1.6)
        im = ax.contourf(Nt, L1, ratio, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='Reds', \
            vmin=0.0, vmax=8.0)

        # Colorbar information
        cbar = ax.cax.colorbar(im)
        ax.cax.set_ylabel(r'RMS$({T}_a)/\overline{T}_a \times 10^3$ (ops)')
        cbar.ax.set_ylim(0,8)
        cbar.ax.set_yticks(np.linspace(0,8,5))
        cbar.ax.set_yticklabels(['0.0','2.0','4.0','6.0','8.0'])


        # Label information
        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(1/8,5/8,3))
        ax.set_xticklabels(['10','30','50'])
        ax.set_yticklabels(['1/8','3/8','5/8'])
        ax.set_xlabel('$n_t$')
        ax.set_ylabel('$p_1$')


    # Save figure
    fig.set_size_inches(6.5,3.4,forward=True)
    plt.savefig(imgdir + 'fig1.png', dpi=300)
    

    print('\tdone with figure 1')