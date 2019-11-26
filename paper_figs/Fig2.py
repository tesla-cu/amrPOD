
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
# from scipy import interpolate
from scipy.interpolate import interp1d

# from scipy.interpolate import RegularGridInterpolator as rgi

def Fig2(datadir, imgdir):

    print('making figure 2 ...')

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

    grid1 = AxesGrid(fig, (0.08,0.54,0.80,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    grid2 = AxesGrid(fig, (0.08,0.15,0.783,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="none")

    # Get locations on fixed compression ratio
    # ax = grid1[1].axes
    # Rim = ax.contour(L2, L1, R, [6,8,10,12], origin='lower', \
    #     colors=('m','b','g','y'))
    # Rim = ax.contour(L2, L1, L0, [0.2,0.4,0.6,0.8], origin='lower', \
    #     colors=('m','b','g','y'))
    # ax.clabel(Rim, Rim.levels, inline=True, inline_spacing=30, \
    #     fmt=r'$r = %i$', fontsize=8, rightside_up=True, \
    #     manual=[(0.1,0.1),(0.2,0.2),(0.3,0.3),(0.4,0.4)])
    # xr6  = Rim.allsegs[0][0][:,0]
    # yr6  = Rim.allsegs[0][0][:,1]
    # xr8  = Rim.allsegs[1][0][:,0]
    # yr8  = Rim.allsegs[1][0][:,1]
    # xr10 = Rim.allsegs[2][0][:,0]
    # yr10 = Rim.allsegs[2][0][:,1]
    # xr12 = Rim.allsegs[3][0][:,0]
    # yr12 = Rim.allsegs[3][0][:,1]

    # yr1 = np.linspace(1/5, 1/5, len(nt))
    # yr2 = np.linspace(2/5, 2/5, len(nt))
    # yr3 = np.linspace(3/5, 3/5, len(nt))

    # xyr6  = np.stack((xr6, yr6),  axis=1)
    # xyr8  = np.stack((xr8, yr8),  axis=1)
    # xyr10 = np.stack((xr10,yr10), axis=1)
    # xyr12 = np.stack((xr12,yr12), axis=1)
    
    # Xr6,  Yr6  = np.meshgrid(xr6,  yr6,  indexing='ij')
    # Xr8,  Yr8  = np.meshgrid(xr8,  yr8,  indexing='ij')
    # Xr10, Yr10 = np.meshgrid(xr10, yr10, indexing='ij')
    # Xr12, Yr12 = np.meshgrid(xr12, yr12, indexing='ij')

    

    for i in range(4):

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

        # img = num_3D/den_3D
        # img = np.exp(num_3D/den_3D)
        img_2D = np.log10(np.sqrt(num_2D/den_2D))
        img_3D = np.log10(np.sqrt(num_3D/den_3D))
        print(np.max(img_3D))
        print(np.min(img_3D))

        # im = ax.contourf(L2, L1, img_3D, 100, origin='lower', \
        #     extent=[l2min,l2max,l1min,l1max], cmap='bwr', \
        #     vmin=-0.8, vmax=0.2)
        im = ax.contourf(Nt, L1, img_3D, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
            vmin=-0.3, vmax=0.3)


        im.set_clim(-0.3, 0.3)
        # im.set_clim(0, 2)
        cbar = ax.cax.colorbar(im)
        ax.cax.set_ylabel('$\log(TC_{im}/TC_{un})$')
        # ax.cax.set_ylabel('$TC_{im}/TC_{un}$')

        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(1/8,5/8,3))
        ax.set_xticklabels([])
        ax.set_yticklabels(['1/8','3/8','5/8'])
        ax.set_ylabel('$\ell_1$')

        ax = grid2[i].axes

        # Interpolate data for each nt onto l1 of interest
        yr1_2D = np.zeros((len(nt)))
        yr2_2D = np.zeros((len(nt)))
        yr3_2D = np.zeros((len(nt)))
        yr1_3D = np.zeros((len(nt)))
        yr2_3D = np.zeros((len(nt)))
        yr3_3D = np.zeros((len(nt)))

        for j in range(len(nt)):

            f_2D = interp1d(l1, img_2D[j,:], kind='linear')
            f_3D = interp1d(l1, img_3D[j,:], kind='linear')

            yr1_2D[j] = f_2D(1/5)
            yr2_2D[j] = f_2D(2/5)
            yr3_2D[j] = f_2D(3/5)
            yr1_3D[j] = f_3D(1/5)
            yr2_3D[j] = f_3D(2/5)
            yr3_3D[j] = f_3D(3/5)


        # Get data along each compression ratio line
        # f = rgi((l2, l1), img, method='linear')
        # fr6  = f(xyr6)
        # fr8  = f(xyr8)
        # fr10 = f(xyr10)
        # fr12 = f(xyr12)

        # ax.plot(nt, yr1_2D, 'r-',  label='2D, $\ell_1 = 1/5$')
        # ax.plot(nt, yr2_2D, 'b-',  label='2D, $\ell_1 = 2/5$')
        # ax.plot(nt, yr3_2D, 'g-',  label='2D, $\ell_1 = 3/5$')
        # ax.plot(nt, yr1_3D, 'r--', label='3D, $\ell_1 = 1/5$')
        # ax.plot(nt, yr2_3D, 'b--', label='3D, $\ell_1 = 2/5$')
        # ax.plot(nt, yr3_3D, 'g--', label='3D, $\ell_1 = 3/5$')
        ax.plot(nt, yr1_2D, 'r-',  label='$\ell_1 = 1/5$')
        ax.plot(nt, yr2_2D, 'b-',  label='$\ell_1 = 2/5$')
        ax.plot(nt, yr3_2D, 'g-',  label='$\ell_1 = 3/5$')
        ax.plot(nt, yr1_3D, 'r--')
        ax.plot(nt, yr2_3D, 'b--')
        ax.plot(nt, yr3_3D, 'g--')

        ax.plot(30, 0, 'k-',  label='2D')
        ax.plot(30, 0, 'k--', label='3D')

        ax.grid('on', linestyle=':')



        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(-0.1,0.1,3))
        ax.set_xticklabels(['10', '30', '50'])
        ax.set_yticklabels(['-0.1','0','0.1'])
        ax.set_xlabel('$n_t$')
        ax.set_ylabel('$\log(TC_{im}/TC_{un})$')
        ax.set_xlim(np.min(nt),np.max(nt))
        if i == 3:
            ax.legend(loc='right', bbox_to_anchor=(1.8, 0.5))


    # print('saving image ...')
    fig.set_size_inches(6.5,3.0,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'fig2.png', dpi=300)


    # img = np.transpose(np.divide(num, den))
    # img = np.divide(num, den)
    # Xaxis, Yaxis = np.meshgrid(xaxis, yaxis, copy=False, indexing='ij')
    # cont_levs = 100
    # xmin = np.min(xaxis)
    # xmax = np.max(xaxis)
    # ymin = np.min(yaxis)
    # ymax = np.max(yaxis)


    # fig = plt.figure(clear=True)
    # ax = fig.add_subplot(1,1,1)

    

    print('\tdone with figure 2')