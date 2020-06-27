import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.interpolate import interp1d


def Fig8_p1_nt_3D(datadir, imgdir):

    print('making figure 8 ...')

    # Set up figure
    fig = plt.figure()

    # Load 2D data
    txtdir_2D = datadir + 'l1_nt/txt_files/'
    TC_R_avg_imp_2D    = np.loadtxt(txtdir_2D + "/TC_R_avg_imp.txt")
    TC_R_avg_unalt_2D  = np.loadtxt(txtdir_2D + "/TC_R_avg_unalt.txt")
    TC_P1_avg_imp_2D   = np.loadtxt(txtdir_2D + "/TC_P1_avg_imp.txt")
    TC_P1_avg_unalt_2D = np.loadtxt(txtdir_2D + "/TC_P1_avg_unalt.txt")
    TC_P2_avg_imp_2D   = np.loadtxt(txtdir_2D + "/TC_P2_avg_imp.txt")
    TC_P2_avg_unalt_2D = np.loadtxt(txtdir_2D + "/TC_P2_avg_unalt.txt")
    TC_A_avg_imp_2D    = np.loadtxt(txtdir_2D + "/TC_A_avg_imp.txt")
    TC_A_avg_unalt_2D  = np.loadtxt(txtdir_2D + "/TC_A_avg_unalt.txt")

    # Load 3D data
    txtdir_3D = datadir + 'l1_nt_3D/txt_files/'
    TC_R_avg_imp_3D    = np.loadtxt(txtdir_3D + "/TC_R_avg_imp.txt")
    TC_R_avg_unalt_3D  = np.loadtxt(txtdir_3D + "/TC_R_avg_unalt.txt")
    TC_P1_avg_imp_3D   = np.loadtxt(txtdir_3D + "/TC_P1_avg_imp.txt")
    TC_P1_avg_unalt_3D = np.loadtxt(txtdir_3D + "/TC_P1_avg_unalt.txt")
    TC_P2_avg_imp_3D   = np.loadtxt(txtdir_3D + "/TC_P2_avg_imp.txt")
    TC_P2_avg_unalt_3D = np.loadtxt(txtdir_3D + "/TC_P2_avg_unalt.txt")
    TC_A_avg_imp_3D    = np.loadtxt(txtdir_3D + "/TC_A_avg_imp.txt")
    TC_A_avg_unalt_3D  = np.loadtxt(txtdir_3D + "/TC_A_avg_unalt.txt")

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

    # Top image set up
    grid1 = AxesGrid(fig, (0.08,0.53,0.84,0.39),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    # Bottom image set up
    grid2 = AxesGrid(fig, (0.08,0.11,0.82,0.39),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="none")


    for i in range(4):

        # Top half of figure --------------------------------------------------
        ax = grid1[i].axes

        if i == 0:
            num_2D = TC_R_avg_imp_2D
            den_2D = TC_R_avg_unalt_2D
            num_3D = TC_R_avg_imp_3D
            den_3D = TC_R_avg_unalt_3D
            ax.set_title(r'$\mathbf{R}$')
        elif i == 1:
            num_2D = TC_P1_avg_imp_2D
            den_2D = TC_P1_avg_unalt_2D
            num_3D = TC_P1_avg_imp_3D
            den_3D = TC_P1_avg_unalt_3D
            ax.set_title(r'$\mathbf{\Phi}$ -- Method 1')
        elif i == 2:
            num_2D = TC_P2_avg_imp_2D
            den_2D = TC_P2_avg_unalt_2D
            num_3D = TC_P2_avg_imp_3D
            den_3D = TC_P2_avg_unalt_3D
            ax.set_title(r'$\mathbf{\Phi}$ -- Method 2')
        elif i == 3:
            num_2D = TC_A_avg_imp_2D
            den_2D = TC_A_avg_unalt_2D
            num_3D = TC_A_avg_imp_3D
            den_3D = TC_A_avg_unalt_3D
            ax.set_title(r'$\mathbf{A}$')

        # Get ratios of AMR over standard operations
        ratio_2D = np.log10(np.sqrt(num_2D/den_2D))
        ratio_3D = np.log10(np.sqrt(num_3D/den_3D))
        print(np.max(ratio_3D))
        print(np.min(ratio_3D))

        # Display contoured image
        im = ax.contourf(Nt, L1, ratio_3D, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
            vmin=-0.25, vmax=0.25)
        # im = ax.contourf(Nt, L1, ratio, 100, origin='lower', \
        #     extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
        #     vmin=0.5, vmax=1.5)

        # Colorbar information
        cbar = ax.cax.colorbar(im)
        ax.cax.set_ylabel(r'$\log(\overline{T}_a/T_s)$ (ops)')
        cbar.ax.set_ylim(-0.25,0.1)
        cbar.ax.set_yticks(np.linspace(-0.2,0.1,4))
        cbar.ax.set_yticklabels(['-0.2','-0.1','0.0','0.1'])
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

        # Bottom half of figure -----------------------------------------------
        ax = grid2[i].axes

        # Interpolate data for each nt onto l1 of interest
        yr1_2D = np.zeros((len(nt))) # for 2D, 1/5
        yr2_2D = np.zeros((len(nt))) # for 2D, 2/5
        yr3_2D = np.zeros((len(nt))) # for 2D, 3/5
        yr1_3D = np.zeros((len(nt))) # for 3D, 1/5
        yr2_3D = np.zeros((len(nt))) # for 3D, 2/5
        yr3_3D = np.zeros((len(nt))) # for 3D, 3/5

        # Interpolate data for each nt
        for j in range(len(nt)):

            f_2D = interp1d(l1, ratio_2D[j,:], kind='linear')
            f_3D = interp1d(l1, ratio_3D[j,:], kind='linear')

            yr1_2D[j] = f_2D(1/5)
            yr2_2D[j] = f_2D(2/5)
            yr3_2D[j] = f_2D(3/5)
            yr1_3D[j] = f_3D(1/5)
            yr2_3D[j] = f_3D(2/5)
            yr3_3D[j] = f_3D(3/5)

        # Plot interpolated data
        ax.plot(nt, yr1_2D, 'r-',  label='$p_1 = 1/5$')
        ax.plot(nt, yr2_2D, 'b-',  label='$p_1 = 2/5$')
        ax.plot(nt, yr3_2D, 'g-',  label='$p_1 = 3/5$')
        ax.plot(nt, yr1_3D, 'r--')
        ax.plot(nt, yr2_3D, 'b--')
        ax.plot(nt, yr3_3D, 'g--')

        # Unused in figure, only used for legend
        ax.plot(30, 0, 'k-',  label='2D')
        ax.plot(30, 0, 'k--', label='3D')

        # Turn grid on
        ax.grid('on', linestyle=':')

        # Label information
        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(-0.1,0.1,3))
        ax.set_xticklabels(['10', '30', '50'])
        ax.set_yticklabels(['-0.1','0','0.1'])
        ax.set_xlabel('$n_t$')
        ax.set_ylabel(r'$\log(\overline{T}_a/T_s)$ (ops)')
        ax.set_xlim(np.min(nt),np.max(nt))
        if i == 3:
            ax.legend(loc='right', bbox_to_anchor=(1.5, 0.3), ncol=2, \
                frameon=True, facecolor='white', framealpha=1.0)


    # Save figure
    fig.set_size_inches(6.5,3.2,forward=True)
    plt.savefig(imgdir + 'Fig8_p1_nt_3D.png', dpi=600)    

    print('\tdone with figure 8')