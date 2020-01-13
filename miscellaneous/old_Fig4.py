
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import interpolate

from scipy.interpolate import RegularGridInterpolator as rgi

def Fig4(datadir, imgdir):

    print('making figure 4 ...')

    # Set up figure
    fig = plt.figure()
    


    # Load data
    txtdir = datadir + 'l1_l2/txt_files/'

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
    l2min = 0/64  # x min
    l2max = 32/64 # x max
    l2inc = 1/64  # x increment
    l1min = 0/64  # y min
    l1max = 32/64 # y max
    l1inc = 1/64  # y increment

    l2 = np.arange(l2min, l2max+l2inc, l2inc)
    l1 = np.arange(l1min, l1max+l1inc, l1inc)

    L2, L1 = np.meshgrid(l2, l1, indexing='ij')

    # Make 2D array of compression ratio
    L0 = 1.0 - L1 - L2
    R  = 16*L0 + 4*L1 + L2
    print(R)

    # R = 16*L0**2 + 4*L1**2 + L2**2 + \
    #      4*L1*L2  + L1*L2   + L0*L2
    # print(R2)

    grid1 = AxesGrid(fig, (0.08,0.54,0.80,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = True,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    grid2 = AxesGrid(fig, (0.08,0.12,0.80,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.1,
        aspect = False,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="none")

    # Get locations on fixed compression ratio
    ax = grid1[1].axes
    Rim = ax.contour(L2, L1, R, [6,8,10,12], origin='lower', \
        colors=('m','b','g','y'))
    # Rim = ax.contour(L2, L1, L0, [0.2,0.4,0.6,0.8], origin='lower', \
    #     colors=('m','b','g','y'))
    ax.clabel(Rim, Rim.levels, inline=True, inline_spacing=30, \
        fmt=r'$r = %i$', fontsize=8, rightside_up=True, \
        manual=[(0.1,0.1),(0.2,0.2),(0.3,0.3),(0.4,0.4)])
    xr6  = Rim.allsegs[0][0][:,0]
    yr6  = Rim.allsegs[0][0][:,1]
    xr8  = Rim.allsegs[1][0][:,0]
    yr8  = Rim.allsegs[1][0][:,1]
    xr10 = Rim.allsegs[2][0][:,0]
    yr10 = Rim.allsegs[2][0][:,1]
    xr12 = Rim.allsegs[3][0][:,0]
    yr12 = Rim.allsegs[3][0][:,1]

    # xr6  = np.linspace(0, 0.5, 9)
    # xr8  = np.linspace(0, 0.5, 9)
    # xr10 = np.linspace(0, 0.5, 9)
    # xr12 = np.linspace(0, 0.5, 9)
    # yr6  = np.linspace(1/8, 1/8, 9)
    # yr8  = np.linspace(2/8, 2/8, 9)
    # yr10 = np.linspace(3/8, 3/8, 9)
    # yr12 = np.linspace(4/8, 4/8, 9)

    xyr6  = np.stack((xr6, yr6),  axis=1)
    xyr8  = np.stack((xr8, yr8),  axis=1)
    xyr10 = np.stack((xr10,yr10), axis=1)
    xyr12 = np.stack((xr12,yr12), axis=1)
    
    Xr6,  Yr6  = np.meshgrid(xr6,  yr6,  indexing='ij')
    Xr8,  Yr8  = np.meshgrid(xr8,  yr8,  indexing='ij')
    Xr10, Yr10 = np.meshgrid(xr10, yr10, indexing='ij')
    Xr12, Yr12 = np.meshgrid(xr12, yr12, indexing='ij')

    

    for i in range(4):

        ax = grid1[i].axes

        if i == 0:
            num = TC_R_avg_imp
            den = TC_R_avg_unalt
            ax.set_title('Computing R')
        elif i == 1:
            num = TC_P1_avg_imp
            den = TC_P1_avg_unalt
            ax.set_title('Computing P1')
        elif i == 2:
            num = TC_P2_avg_imp
            den = TC_P2_avg_unalt
            ax.set_title('Computing P2')
        elif i == 3:
            num = TC_A_avg_imp
            den = TC_A_avg_unalt
            ax.set_title('Computing A')

        img = num/den
        # img = np.exp(num/den)
        # img = np.log10(np.sqrt(num/den))
        print(np.max(img))
        print(np.min(img))

        # im = ax.contourf(L2, L1, img, 100, origin='lower', \
        #     extent=[l2min,l2max,l1min,l1max], cmap='bwr', \
        #     vmin=-0.8, vmax=0.2)
        im = ax.contourf(L2, L1, img, 100, origin='lower', \
            extent=[l2min,l2max,l1min,l1max], cmap='bwr', \
            vmin=0, vmax=2)


        # im.set_clim(-0.8, 0.8)
        im.set_clim(0, 2)
        cbar = ax.cax.colorbar(im)
        # ax.cax.set_ylabel('$\log(TC_{im}/TC_{un})$')
        ax.cax.set_ylabel('$TC_{im}/TC_{un}$')



        # ax.set_xticks([ntmin,(ntmax+ntmin)/2,ntmax])
        # ax.set_yticks([l1min,(l1max+l1min)/2,l1max])
        ax.set_xticks(np.linspace(1/8,3/8,3))
        ax.set_yticks(np.linspace(1/8,3/8,3))
        ax.set_xticklabels([])
        ax.set_yticklabels(['1/8','2/8','3/8'])
        ax.set_ylabel('$\ell_1$')


        ax = grid2[i].axes

        # Get data along each compression ratio line
        f = rgi((l2, l1), img, method='linear')
        fr6  = f(xyr6)
        fr8  = f(xyr8)
        fr10 = f(xyr10)
        fr12 = f(xyr12)

        ax.plot(xr6,  fr6,  'm')
        ax.plot(xr8,  fr8,  'b')
        ax.plot(xr10, fr10, 'g')
        ax.plot(xr12, fr12, 'y')
        ax.grid('on')



        # ax.set_xticks([ntmin,(ntmax+ntmin)/2,ntmax])
        # ax.set_yticks([l1min,(l1max+l1min)/2,l1max])
        ax.set_xticks(np.linspace(1/8,3/8,3))
        # ax.set_yticks(np.linspace(1/8,3/8,3))
        ax.set_xticklabels(['1/8','2/8','3/8'])
        # ax.set_yticklabels(['1/8','2/8','3/8'])
        ax.set_xlabel('$\ell_2$')
        # ax.set_ylabel('$\ell_1$')
        ax.set_xlim(0,0.5)


    # print('saving image ...')
    fig.set_size_inches(6.5,3.25,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'fig4.png', dpi=300)


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

    

    print('\tdone with figure 4')