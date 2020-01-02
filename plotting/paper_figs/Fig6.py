import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def comp_avg_rms(arr):

    avg = np.mean(arr)
    rms = np.std(arr - avg)
    return avg, rms


def Fig6(datadir, imgdir):

    print('making figure 6 ...')

    # Set up figure
    fig = plt.figure()

    # Set up data directories
    f1dir = datadir + 'AMR_sim/f1_3D/'
    f2dir = datadir + 'AMR_sim/f2_3D/'

    # Data information
    ntmin = 10
    ntmax = 60
    ntinc = 10

    nt = np.arange(ntmin, ntmax+ntinc, ntinc)

    # Image set up
    grid = AxesGrid(fig, (0.10,0.10,0.84,0.82),
        nrows_ncols = (2, 3),
        axes_pad    = 0.2,
        aspect      = False,
        label_mode  = "L",
        share_all   = True,
        direction   = "column",
        cbar_mode   = "none")

    # Initialize data to store CPU time
    #  - index 1: nt
    #  - index 2: num of finest levels
    #  - index 3: first is avg, second is rms
    R_AMR    = np.zeros((len(nt), 2, 2))
    R_std    = np.zeros((len(nt), 2, 2))
    Phi_AMR1 = np.zeros((len(nt), 2, 2))
    Phi_AMR2 = np.zeros((len(nt), 2, 2))
    Phi_std  = np.zeros((len(nt), 2, 2))
    A_AMR    = np.zeros((len(nt), 2, 2))
    A_std    = np.zeros((len(nt), 2, 2))

    for i, n in enumerate(nt):
        for j in range(2):
            if j==0:
                fold = (f1dir + 'CPU_timing_nt%i/' % n)
            elif j==1:
                fold = (f2dir + 'CPU_timing_nt%i/' % n)

            # CPU times for R
            CPUs = np.loadtxt(fold + 'R_CPU_AMR.txt')
            R_AMR[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'R_CPU_standard.txt')
            R_std[i,j,:] = comp_avg_rms(CPUs)

            # CPU times for Phi
            CPUs = np.loadtxt(fold + 'Phi_CPU_AMR1.txt')
            Phi_AMR1[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'Phi_CPU_AMR2.txt')
            Phi_AMR2[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'Phi_CPU_standard.txt')
            Phi_std[i,j,:] = comp_avg_rms(CPUs)

            # CPU times for A
            CPUs = np.loadtxt(fold + 'A_CPU_AMR.txt')
            A_AMR[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'A_CPU_standard.txt')
            A_std[i,j,:] = comp_avg_rms(CPUs)


    # Plot R CPU times for f1 -------------------------------------------------
    ax = grid[0].axes

    # ax.errorbar(nt, R_std[:,0,0], yerr=2*R_std[:,0,1], linestyle='-', \
    #     elinewidth=1, capsize=5, color='k', label='$\overline{T}_s$')
    # ax.errorbar(nt, R_AMR[:,0,0], yerr=2*R_AMR[:,0,1], linestyle='--', \
    #     elinewidth=1, capsize=5, color='b', label='$\overline{T}_a$')
    ax.plot(nt, R_std[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, R_AMR[:,0,0], linestyle='--', \
        color='b', label='$\overline{T}_a$')
    ax.set_title(r'(a) $\bf{R}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')
    ax.legend(loc=4, frameon=True, facecolor='white', framealpha=1.0, \
        bbox_to_anchor=(1.1,-0.1))

    # Plot R CPU times for f2 -------------------------------------------------
    ax = grid[1].axes

    # ax.errorbar(nt, R_std[:,1,0], yerr=2*R_std[:,1,1], linestyle='-', \
    #     elinewidth=1, capsize=5, color='k', label='$\overline{T}_s$')
    # ax.errorbar(nt, R_AMR[:,1,0], yerr=2*R_AMR[:,1,1], linestyle='--', \
    #     elinewidth=1, capsize=5, color='b', label='$\overline{T}_a$')
    ax.plot(nt, R_std[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, R_AMR[:,1,0], linestyle='--', \
        color='b', label='$\overline{T}_a$')
    ax.set_title(r'(d) $\bf{R}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')

    # Plot Phi CPU times for f1 -----------------------------------------------
    ax = grid[2].axes

    # ax.errorbar(nt, Phi_std[:,0,0], yerr=2*Phi_std[:,0,1], linestyle='-', \
    #     elinewidth=1, capsize=5, color='k', label='$\overline{T}_s$')
    # ax.errorbar(nt, Phi_AMR1[:,0,0], yerr=2*Phi_AMR1[:,0,1], linestyle='--', \
    #     elinewidth=1, capsize=5, color='b', label='$\overline{T}_a$ -- M1')
    # ax.errorbar(nt, Phi_AMR2[:,0,0], yerr=2*Phi_AMR2[:,0,1], linestyle=':', \
    #     elinewidth=1, capsize=5, color='orange',label='$\overline{T}_a$ -- M2')
    ax.plot(nt, Phi_std[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, Phi_AMR1[:,0,0], linestyle='--', \
        color='b', label='$\overline{T}_a$ -- M1')
    ax.plot(nt, Phi_AMR2[:,0,0], linestyle=':', \
        color='orange', label='$\overline{T}_a$ -- M2')
    ax.set_title(r'(b) $\bf{\Phi}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')
    ax.legend(loc=4, frameon=True, facecolor='white', framealpha=1.0, \
        bbox_to_anchor=(1.1,-0.1))

    # Plot Phi CPU times for f2 -----------------------------------------------
    ax = grid[3].axes

    # ax.errorbar(nt, Phi_std[:,1,0], yerr=2*Phi_std[:,1,1], linestyle='-', \
    #     elinewidth=1, capsize=5, color='k', label='$\overline{T}_s$')
    # ax.errorbar(nt, Phi_AMR1[:,1,0], yerr=2*Phi_AMR1[:,1,1], linestyle='--', \
    #     elinewidth=1, capsize=5, color='b', label='$\overline{T}_a$ -- M1')
    # ax.errorbar(nt, Phi_AMR2[:,1,0], yerr=2*Phi_AMR2[:,1,1], linestyle=':', \
    #     elinewidth=1, capsize=5, color='orange',label='$\overline{T}_a$ -- M2')
    ax.plot(nt, Phi_std[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, Phi_AMR1[:,1,0], linestyle='--', \
        color='b', label='$\overline{T}_a$ -- M1')
    ax.plot(nt, Phi_AMR2[:,1,0], linestyle=':', \
        color='orange',label='$\overline{T}_a$ -- M2')
    ax.set_title(r'(e) $\bf{\Phi}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')

    # Plot A CPU times for f1 -------------------------------------------------
    ax = grid[4].axes

    ax.plot(nt, A_std[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, A_AMR[:,0,0], linestyle='--', \
        color='b', label='$\overline{T}_a$')
    ax.set_title(r'(c) $\bf{A}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')
    ax.legend(loc=4, frameon=True, facecolor='white', framealpha=1.0, \
        bbox_to_anchor=(1.1,-0.1))

    # Plot A CPU times for f2 -------------------------------------------------
    ax = grid[5].axes

    ax.plot(nt, A_std[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_s$')
    ax.plot(nt, A_AMR[:,1,0], linestyle='--', \
        color='b', label='$\overline{T}_a$')
    ax.set_title(r'(f) $\bf{A}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','30','40','50','60'])
    ax.set_xlabel(r'$n_t$')
    ax.set_ylabel(r'$\overline{T}$ [CPU (s)]')

    # Save figure -------------------------------------------------------------
    fig.set_size_inches(6.5,3.75,forward=True)
    plt.savefig(imgdir + 'fig6.png', dpi=300)    

    print('\tdone with figure 6')