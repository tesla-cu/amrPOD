import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def comp_avg_rms(arr):

    avg = np.mean(arr)
    rms = np.std(arr - avg)
    return avg, rms

def Fig13_CPU(datadir, imgdir):

    print('making figure 13 ...')

    # Set up figure
    fig = plt.figure()

    # Set up data directories
    f1dir_CPU = datadir + 'AMR_gen/f1_3D/intel_timing/'
    f2dir_CPU = datadir + 'AMR_gen/f2_3D/intel_timing/'
    f3dir_CPU = datadir + 'AMR_gen/f3_3D/intel_timing/'

    f1dir_TC = datadir + 'AMR_gen/f1_3D/TC_nt/txt_files/'
    f2dir_TC = datadir + 'AMR_gen/f2_3D/TC_nt/txt_files/'
    f3dir_TC = datadir + 'AMR_gen/f3_3D/TC_nt/txt_files/'

    # Data information
    ntmin = 10
    ntmax = 80
    ntinc = 10
    nt = np.arange(ntmin, ntmax+ntinc, ntinc)

    # Image set up
    grid = AxesGrid(fig, (0.085,0.075,0.90,0.89),
        nrows_ncols = (3, 3),
        axes_pad    = 0.2,
        aspect      = False,
        label_mode  = "L",
        share_all   = True,
        direction   = "column",
        cbar_mode   = "none")

    # Initialize data to store CPU time and TC counts
    #  - index 1: nt
    #  - index 2: num of finest levels
    #  - index 3: first is avg, second is rms
    R_AMR_CPU    = np.zeros((len(nt), 3, 2))
    R_std_CPU    = np.zeros((len(nt), 3, 2))
    Phi_AMR1_CPU = np.zeros((len(nt), 3, 2))
    Phi_AMR2_CPU = np.zeros((len(nt), 3, 2))
    Phi_std_CPU  = np.zeros((len(nt), 3, 2))
    A_AMR_CPU    = np.zeros((len(nt), 3, 2))
    A_std_CPU    = np.zeros((len(nt), 3, 2))

    R_AMR_TC     = np.zeros((len(nt), 3))
    R_std_TC     = np.zeros((len(nt), 3))
    Phi_AMR1_TC  = np.zeros((len(nt), 3))
    Phi_AMR2_TC  = np.zeros((len(nt), 3))
    Phi_std_TC   = np.zeros((len(nt), 3))
    A_AMR_TC     = np.zeros((len(nt), 3))
    A_std_TC     = np.zeros((len(nt), 3))

    # Load CPU times
    for i, n in enumerate(nt):
        for j in range(3):
            if j==0:
                fold = (f1dir_CPU + 'CPU_timing_nt%i/' % n)
            elif j==1:
                fold = (f2dir_CPU + 'CPU_timing_nt%i/' % n)
            elif j==2:
                fold = (f3dir_CPU + 'CPU_timing_nt%i/' % n)

            # CPU times for R
            CPUs = np.loadtxt(fold + 'R_CPU_AMR.txt')
            R_AMR_CPU[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'R_CPU_standard.txt')
            R_std_CPU[i,j,:] = comp_avg_rms(CPUs)

            # CPU times for Phi
            CPUs = np.loadtxt(fold + 'Phi_CPU_AMR1.txt')
            Phi_AMR1_CPU[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'Phi_CPU_AMR2.txt')
            Phi_AMR2_CPU[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'Phi_CPU_standard.txt')
            Phi_std_CPU[i,j,:] = comp_avg_rms(CPUs)

            # CPU times for A
            CPUs = np.loadtxt(fold + 'A_CPU_AMR.txt')
            A_AMR_CPU[i,j,:] = comp_avg_rms(CPUs)
            CPUs = np.loadtxt(fold + 'A_CPU_standard.txt')
            A_std_CPU[i,j,:] = comp_avg_rms(CPUs)

    # Load TC counts
    for j in range(3):
        if j==0:
            fold = f1dir_TC
        elif j==1:
            fold = f2dir_TC
        elif j==2:
            fold = f3dir_TC

        # TC counts for R
        TCs = np.loadtxt(fold + 'TC_R_avg_imp.txt')
        R_AMR_TC[:,j] = TCs[:-2]
        TCs = np.loadtxt(fold + 'TC_R_avg_unalt.txt')
        R_std_TC[:,j] = TCs[:-2]

        # TC counts for Phi
        TCs = np.loadtxt(fold + 'TC_P1_avg_imp.txt')
        Phi_AMR1_TC[:,j] = TCs[:-2]
        TCs = np.loadtxt(fold + 'TC_P2_avg_imp.txt')
        Phi_AMR2_TC[:,j] = TCs[:-2]
        TCs = np.loadtxt(fold + 'TC_P1_avg_unalt.txt')
        Phi_std_TC[:,j] = TCs[:-2]

        # TC counts for A
        TCs = np.loadtxt(fold + 'TC_A_avg_imp.txt')
        A_AMR_TC[:,j] = TCs[:-2]
        TCs = np.loadtxt(fold + 'TC_A_avg_unalt.txt')
        A_std_TC[:,j] = TCs[:-2]

    # Normalize the CPU times and TC counts
    for j in range(3):
        R_AMR_CPU   [:,j,:] = R_AMR_CPU   [:,j,:] / R_std_CPU  [0,j,0]
        Phi_AMR1_CPU[:,j,:] = Phi_AMR1_CPU[:,j,:] / Phi_std_CPU[0,j,0]
        Phi_AMR2_CPU[:,j,:] = Phi_AMR2_CPU[:,j,:] / Phi_std_CPU[0,j,0]
        A_AMR_CPU   [:,j,:] = A_AMR_CPU   [:,j,:] / A_std_CPU  [0,j,0]

        R_std_CPU  [:,j,:] = R_std_CPU  [:,j,:] / R_std_CPU  [0,j,0]
        Phi_std_CPU[:,j,:] = Phi_std_CPU[:,j,:] / Phi_std_CPU[0,j,0]
        A_std_CPU  [:,j,:] = A_std_CPU  [:,j,:] / A_std_CPU  [0,j,0]

        R_AMR_TC   [:,j] = R_AMR_TC   [:,j] / R_std_TC  [0,j]
        Phi_AMR1_TC[:,j] = Phi_AMR1_TC[:,j] / Phi_std_TC[0,j]
        Phi_AMR2_TC[:,j] = Phi_AMR2_TC[:,j] / Phi_std_TC[0,j]
        A_AMR_TC   [:,j] = A_AMR_TC   [:,j] / A_std_TC  [0,j]

        R_std_TC  [:,j] = R_std_TC  [:,j] / R_std_TC  [0,j]
        Phi_std_TC[:,j] = Phi_std_TC[:,j] / Phi_std_TC[0,j]
        A_std_TC  [:,j] = A_std_TC  [:,j] / A_std_TC  [0,j]

    # Plot R CPU times for f1 -------------------------------------------------
    ax = grid[0].axes
    ax.plot([0,0], [0,0], linestyle='-', color='k', linewidth=5.0, label='CPU')
    ax.plot([0,0], [0,0], linestyle='-', color='b', linewidth=5.0, label='ops')
    ax.plot(nt, R_std_TC[:,0], linestyle='-', color='b')
    ax.plot(nt, R_AMR_TC[:,0], linestyle='--', color='b')
    ax.plot(nt, R_std_CPU[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, R_AMR_CPU[:,0,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(a) $\bf{R}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)
    ax.legend(loc=4, bbox_to_anchor=(1.05,-0.15), frameon=True, facecolor='white', framealpha=1.0)

    # Plot R CPU times for f2 -------------------------------------------------
    ax = grid[1].axes
    ax.plot(nt, R_std_TC[:,1], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, R_AMR_TC[:,1], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$')
    ax.plot(nt, R_std_CPU[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, R_AMR_CPU[:,1,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(d) $\bf{R}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Plot R CPU times for f3 -------------------------------------------------
    ax = grid[2].axes
    ax.plot(nt, R_std_TC[:,2], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, R_AMR_TC[:,2], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$')
    ax.plot(nt, R_std_CPU[:,2,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, R_AMR_CPU[:,2,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(g) $\bf{R}$, $f=3$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Plot Phi CPU times for f1 -----------------------------------------------
    ax = grid[3].axes
    ax.plot([0,0], [0,0], linestyle='-', color='k', linewidth=5.0, label='CPU')
    ax.plot([0,0], [0,0], linestyle='-', color='b', linewidth=5.0, label='ops')
    ax.plot(nt, Phi_std_TC[:,0], linestyle='-', color='b')
    ax.plot(nt, Phi_AMR1_TC[:,0], linestyle='--', color='b')
    ax.plot(nt, Phi_AMR2_TC[:,0], linestyle=':', color='b')
    ax.plot(nt, Phi_std_CPU[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, Phi_AMR1_CPU[:,0,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$ (M1)')
    ax.plot(nt, Phi_AMR2_CPU[:,0,0], linestyle=':', \
        color='k', label='$\overline{T}_\mathrm{AMR}$ (M2)')
    ax.set_title(r'(b) $\bf{\Phi}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)
    ax.legend(loc=4, bbox_to_anchor=(1.05,-0.15), frameon=True, facecolor='white', framealpha=1.0)

    # Plot Phi CPU times for f2 -----------------------------------------------
    ax = grid[4].axes
    ax.plot(nt, Phi_std_TC[:,1], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, Phi_AMR1_TC[:,1], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$ (M1)')
    ax.plot(nt, Phi_AMR2_TC[:,1], linestyle=':', \
        color='b', label='$\overline{T}_\mathrm{AMR}$ (M2)')
    ax.plot(nt, Phi_std_CPU[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, Phi_AMR1_CPU[:,1,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$ (M1)')
    ax.plot(nt, Phi_AMR2_CPU[:,1,0], linestyle=':', \
        color='k',label='$\overline{T}_\mathrm{AMR}$ (M2)')
    ax.set_title(r'(e) $\bf{\Phi}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Plot Phi CPU times for f3 -----------------------------------------------
    ax = grid[5].axes
    ax.plot(nt, Phi_std_TC[:,2], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, Phi_AMR1_TC[:,2], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$ (M1)')
    ax.plot(nt, Phi_AMR2_TC[:,2], linestyle=':', \
        color='b', label='$\overline{T}_\mathrm{AMR}$ (M2)')
    ax.plot(nt, Phi_std_CPU[:,2,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, Phi_AMR1_CPU[:,2,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$ (M1)')
    ax.plot(nt, Phi_AMR2_CPU[:,2,0], linestyle=':', \
        color='k',label='$\overline{T}_\mathrm{AMR}$ (M2)')
    ax.set_title(r'(h) $\bf{\Phi}$, $f=3$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Plot A CPU times for f1 -------------------------------------------------
    ax = grid[6].axes
    ax.plot([0,0], [0,0], linestyle='-', color='k', linewidth=5.0, label='CPU')
    ax.plot([0,0], [0,0], linestyle='-', color='b', linewidth=5.0, label='ops')
    ax.plot(nt, A_std_TC[:,0], linestyle='-', color='b')
    ax.plot(nt, A_AMR_TC[:,0], linestyle='--', color='b')
    ax.plot(nt, A_std_CPU[:,0,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, A_AMR_CPU[:,0,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(c) $\bf{A}$, $f=1$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)
    ax.legend(loc=4, bbox_to_anchor=(1.05,-0.15), frameon=True, facecolor='white', framealpha=1.0)

    # Plot A CPU times for f2 -------------------------------------------------
    ax = grid[7].axes
    ax.plot(nt, A_std_TC[:,1], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, A_AMR_TC[:,1], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$')
    ax.plot(nt, A_std_CPU[:,1,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, A_AMR_CPU[:,1,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(f) $\bf{A}$, $f=2$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Plot A CPU times for f3 -------------------------------------------------
    ax = grid[8].axes
    ax.plot(nt, A_std_TC[:,2], linestyle='-', \
        color='b', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, A_AMR_TC[:,2], linestyle='--', \
        color='b', label='$\overline{T}_\mathrm{AMR}$')
    ax.plot(nt, A_std_CPU[:,2,0], linestyle='-', \
        color='k', label='$\overline{T}_\mathrm{std}$')
    ax.plot(nt, A_AMR_CPU[:,2,0], linestyle='--', \
        color='k', label='$\overline{T}_\mathrm{AMR}$')
    ax.set_title(r'(i) $\bf{A}$, $f=3$', loc='left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid('on', linestyle=':')
    
    # Label information
    ax.set_xticks(nt)
    ax.set_xticklabels(['10','20','','40','','60','','80'])
    ax.set_xlabel(r'$N_\mathrm{t}$')
    ax.set_ylabel(r'$\overline{T}_\mathrm{AMR}/T_\mathrm{std}(N_\mathrm{t}=10)$')
    ax.set_xlim(10,80)

    # Save figure -------------------------------------------------------------
    fig.set_size_inches(6.5,5.5,forward=True)
    plt.savefig(imgdir + 'Fig13_CPU.png', dpi=600)    

    print('\tdone with figure 13')