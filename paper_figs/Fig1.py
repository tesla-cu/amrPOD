
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

    grid = AxesGrid(fig, (0.08,0.54,0.80,0.35),
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

        img = np.log10(num/den)
        print(np.max(img))
        print(np.min(img))


        
        # ax.title = i

        # im = ax.contourf(Nt, L1, img, 100, origin='lower', cmap='bwr')
        im = ax.contourf(Nt, L1, img, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='bwr', \
            vmin=-0.3, vmax=0.3)
        im.set_clim(-0.3, 0.3)
        # cbar = ax.figure.colorbar(im, ax=ax)
        cbar = ax.cax.colorbar(im)
        ax.cax.set_ylabel('$\log(TC_{im}/TC_{un})$')

        # ax.set_xticks([ntmin,(ntmax+ntmin)/2,ntmax])
        # ax.set_yticks([l1min,(l1max+l1min)/2,l1max])
        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(1/8,5/8,3))
        ax.set_xticklabels([])
        ax.set_yticklabels(['1/8','3/8','5/8'])
        ax.set_ylabel('$\ell_1$')


    grid = AxesGrid(fig, (0.08,0.12,0.80,0.35),
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

        img = (num/den)*100
        print(np.max(img))
        print(np.min(img))


        
        # ax.title = i

        # im = ax.contourf(Nt, L1, img, 100, origin='lower', cmap='bwr')
        im = ax.contourf(Nt, L1, img, 100, origin='lower', \
            extent=[ntmin,ntmax,l1min,l1max], cmap='Reds', \
            vmin=0.0, vmax=1.6)
        im.set_clim(0.0, 1.6)
        # cbar = ax.figure.colorbar(im, ax=ax)
        cbar = ax.cax.colorbar(im)
        # ax.cax.set_ylabel('rms$(TC_{im})/$rms$(TC_{un})$ X 100')
        ax.cax.set_ylabel('Fluctuations')

        # ax.set_xticks([ntmin,(ntmax+ntmin)/2,ntmax])
        # ax.set_yticks([l1min,(l1max+l1min)/2,l1max])
        ax.set_xticks(np.linspace(10,50,3))
        ax.set_yticks(np.linspace(1/8,5/8,3))
        ax.set_xticklabels(['10','30','50'])
        ax.set_yticklabels(['1/8','3/8','5/8'])
        ax.set_xlabel('$n_t$')
        ax.set_ylabel('$\ell_1$')


    # print('saving image ...')
    fig.set_size_inches(6.5,3.25,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'fig1.png', dpi=300)


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

    

    print('\tdone with figure 1')