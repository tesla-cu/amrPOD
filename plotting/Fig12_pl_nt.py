import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def Fig12_pl_nt(datadir, imgdir):

    print('making figure 12 ...')

    # Set up figure
    fig = plt.figure()
    cmap_acc = mpl.cm.get_cmap('Accent')
    # lvlC0    = np.tile(cmap_acc(.125), (3,1))
    # lvlC1    = np.tile(cmap_acc(.375), (3,1))
    # lvlC2    = np.tile(cmap_acc(.625), (3,1))
    clrs = [cmap_acc(.125), cmap_acc(.375), cmap_acc(.625)]
    lins = ['-', ':']
    mrks = ['o', 's', '^']

    # Set up data directories
    f1dir = datadir + 'AMR_gen/f1_3D/'
    f2dir = datadir + 'AMR_gen/f2_3D/'

    # Data information
    ntmin = 10
    ntmax = 60
    ntinc = 10

    nt = np.arange(ntmin, ntmax+ntinc, ntinc)

    # Image set up
    grid = AxesGrid(fig, (0.10,0.17,0.88,0.72),
        nrows_ncols = (1, 2),
        axes_pad    = 0.1,
        aspect      = False,
        label_mode  = "L",
        share_all   = True,
        cbar_mode   = "none")

    for i in range(2):

        if i==0:
            ps = np.zeros((2, len(nt)))
            pms = np.zeros((2, len(nt)))
        elif i==1:
            ps = np.zeros((3, len(nt)))
            pms = np.zeros((3, len(nt)))

        # Plot grid composition -----------------------------------------------
        ax = grid[i].axes

        for j, n in enumerate(nt):

            if i==0:
                fold = (f1dir + 'CPU_timing_nt%i/' % n)
            elif i==1:
                fold = (f2dir + 'CPU_timing_nt%i/' % n)

            ps[:,j]  = np.loadtxt(fold + 'p.txt')
            pms[:,j] = np.loadtxt(fold + 'pm.txt')

        # print(pms)
        for j in range(ps.shape[0]):

            ax.plot(nt, ps[j,:],  linestyle=lins[0], color=clrs[j], \
                marker=mrks[j], markersize=4, linewidth=1, \
                label=('$\overline{p}_%i$'%j))
            ax.plot(nt, pms[j,:], linestyle=lins[1], color=clrs[j], \
                marker=mrks[j], markersize=4, linewidth=1, \
                label=('$p_%i^m$'%j))

            ax.grid('on', linestyle=':')

            # Label information
            ax.set_xticks(nt)
            ax.set_xticklabels(['10','20','30','40','50','60'])
            ax.set_xlabel(r'$n_t$')
            ax.set_ylabel('Grid Composition')

        if i==0:
            ax.set_title(r'(a) $f=1$', loc='left')
        elif i==1:
            ax.set_title(r'(b) $f=2$', loc='left')
    
    # Set legend with ordered elements
    handles, labels = ax.get_legend_handles_labels()
    order = [0,2,4,1,3,5]
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
        loc=3, frameon=True, facecolor='white', framealpha=1.0, \
        ncol=2, bbox_to_anchor=(-0.8,0.05))

    # Save figure -------------------------------------------------------------
    fig.set_size_inches(4.5, 2.0,forward=True)
    plt.savefig(imgdir + 'Fig12_pl_nt.png', dpi=300)    

    print('\tdone with figure 12')