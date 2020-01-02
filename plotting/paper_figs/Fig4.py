
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import interpolate

from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.optimize import curve_fit

# =========================================================================== #
# Equation to fit:
#
#   l0_coeff*l0 + l1_coeff*l1 + l2_coeff*l2
# =========================================================================== #
def fit_equation_linear(l, l0_coeff, l1_coeff, l2_coeff):
    l0,l1,l2 = l
    return l0_coeff*l0 + \
           l1_coeff*l1 + \
           l2_coeff*l2

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0sq_coeff*l0^2  + l1sq_coeff*l1^2  + l2sq_coeff*l2^2 + ...               #
#   l0l1_coeff*l0*l1 + l0l2_coeff*l0*l2 + l1l2_coeff*l2^2                     #
# =========================================================================== #
def fit_equation_quadratic(l, l0sq_coeff, l1sq_coeff, l2sq_coeff, \
                              l0l1_coeff, l0l2_coeff, l1l2_coeff):
    l0,l1,l2 = l
    return   l0sq_coeff*l0**2 +   l1sq_coeff*l1**2 +   l2sq_coeff*l2**2 + \
           2*l0l1_coeff*l0*l1 + 2*l0l2_coeff*l0*l2 + 2*l1l2_coeff*l1*l2

# =========================================================================== #
# Figure 4
# =========================================================================== #
def Fig4(datadir, imgdir):

    print('making figure 4 ...')

    # Set up figure
    fig = plt.figure()
    
    # Load data
    l_txtdir  = datadir + 'l1_l0/txt_files/'
    lc_txtdir = datadir + 'lc1l1_lc0l0/txt_files/'

    TC_R_avg_imp    = np.loadtxt(l_txtdir  + "/TC_R_avg_imp.txt")
    TC_R_avg_unalt  = np.loadtxt(l_txtdir  + "/TC_R_avg_unalt.txt")
    TC_R_rms_imp    = np.loadtxt(l_txtdir  + "/TC_R_rms_imp.txt")
    TC_R_rms_unalt  = np.loadtxt(l_txtdir  + "/TC_R_rms_unalt.txt")

    TC_P1_avg_imp   = np.loadtxt(lc_txtdir + "/TC_P1_avg_imp.txt")
    TC_P1_avg_unalt = np.loadtxt(lc_txtdir + "/TC_P1_avg_unalt.txt")
    TC_P1_rms_imp   = np.loadtxt(lc_txtdir + "/TC_P1_rms_imp.txt")
    TC_P1_rms_unalt = np.loadtxt(lc_txtdir + "/TC_P1_rms_unalt.txt")

    TC_P2_avg_imp   = np.loadtxt(l_txtdir  + "/TC_P2_avg_imp.txt")
    TC_P2_avg_unalt = np.loadtxt(l_txtdir  + "/TC_P2_avg_unalt.txt")
    TC_P2_rms_imp   = np.loadtxt(l_txtdir  + "/TC_P2_rms_imp.txt")
    TC_P2_rms_unalt = np.loadtxt(l_txtdir  + "/TC_P2_rms_unalt.txt")

    TC_A_avg_imp    = np.loadtxt(lc_txtdir + "/TC_A_avg_imp.txt")
    TC_A_avg_unalt  = np.loadtxt(lc_txtdir + "/TC_A_avg_unalt.txt")
    TC_A_rms_imp    = np.loadtxt(lc_txtdir + "/TC_A_rms_imp.txt")
    TC_A_rms_unalt  = np.loadtxt(lc_txtdir + "/TC_A_rms_unalt.txt")

    # Data information
    l0min = 0/64  # x min
    l0max = 32/64 # x max
    l0inc = 1/64  # x increment
    l1min = 0/64  # y min
    l1max = 32/64 # y max
    l1inc = 1/64  # y increment

    l0 = np.arange(l0min, l0max+l0inc, l0inc)
    l1 = np.arange(l1min, l1max+l1inc, l1inc)

    L0, L1 = np.meshgrid(l0, l1, indexing='ij')

    # Get corresponding l2 values
    L2 = 1.0 - L0 - L1

    # Set up grids for plotting
    grid1 = AxesGrid(fig, (0.06,0.55,0.88,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.29,
        aspect = True,
        label_mode = "all",
        share_all = False,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    grid2 = AxesGrid(fig, (0.06,0.12,0.88,0.35),
        nrows_ncols = (1, 4),
        axes_pad = 0.18,
        aspect = True,
        label_mode = "all",
        share_all = False,
        cbar_location="right",
        cbar_mode="each",
        cbar_size="7%",
        cbar_pad="3%")


    for i in range(4):

        # Top half of figure --------------------------------------------------
        ax = grid1[i].axes

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

        ratio = num/den
        # ratio = np.exp(num/den)
        # ratio = np.log10(np.sqrt(num/den))
        print(np.max(ratio))
        print(np.min(ratio))

        im = ax.contourf(L0, L1, ratio, 100, origin='lower', \
            extent=[l0min,l0max,l1min,l1max], cmap='bwr', \
            vmin=0.3, vmax=1.7)

        # Colorbar information
        cbar = ax.cax.colorbar(im, format='%0.1f')
        ax.cax.set_ylabel(r'$\overline{T}_a/T_s$ [ops]')
        cbar.ax.set_ylim(0.32,1.27)
        cbar.ax.set_yticks(np.linspace(0.4,1.2,5))


        # Label inforrmation
        ax.set_xticks(np.linspace(1/8,3/8,3))
        ax.set_yticks(np.linspace(1/8,3/8,3))
        ax.set_xticklabels([])
        ax.set_yticklabels([]) # add text later
        ax.set_xlabel('')
        if i==0:
            ax.set_ylabel(r'$p_1$', labelpad=20)
            ax.text(-0.02, 1/8, '1/8', ha='right', va='center')
            ax.text(-0.02, 2/8, '2/8', ha='right', va='center')
            ax.text(-0.02, 3/8, '3/8', ha='right', va='center')
        if i==2:
            ax.set_ylabel(r'$p_1$', labelpad=-3.5)
        elif i==1 or i==3:
            ax.set_ylabel(r'$p_1 = p_0^m$', labelpad=-3.5)


        # Bottom half of figure -----------------------------------------------
        
        ax = grid2[i].axes

        l0_1D    = np.reshape(L0,    (-1))
        l1_1D    = np.reshape(L1,    (-1))
        l2_1D    = np.reshape(L2,    (-1))
        ratio_1D = np.reshape(ratio, (-1))

        # best_vals, covar = \
        #     curve_fit(fit_equation_quadratic,(l0_1D,l1_1D,l2_1D), ratio_1D)

        # print('best_vals: {}'.format(best_vals))
        # print('error: {}'.format(np.sqrt(np.diag(covar))))

        # best_vals, covar = \
        #     curve_fit(fit_equation_linear,   (l0_1D,l1_1D,l2_1D), ratio_1D)

        # print('best_vals: {}'.format(best_vals))
        # print('error: {}'.format(np.sqrt(np.diag(covar))))

        ratio_fit = np.zeros((len(ratio_1D)))
        if i==0:
            best_vals, covar = \
                curve_fit(fit_equation_quadratic,(l0_1D,l1_1D,l2_1D), ratio_1D)

            for j in range(len(ratio_1D)):
                ratio_fit[j] = fit_equation_quadratic(\
                    (l0_1D[j], l1_1D[j], l2_1D[j]), \
                    best_vals[0], best_vals[1], best_vals[2], \
                    best_vals[3], best_vals[4], best_vals[5])

        elif i==1 or i==2 or i==3:
            best_vals, covar = \
                curve_fit(fit_equation_linear,   (l0_1D,l1_1D,l2_1D), ratio_1D)

            for j in range(len(ratio_1D)):
                ratio_fit[j] = fit_equation_linear(\
                    (l0_1D[j], l1_1D[j], l2_1D[j]), \
                    best_vals[0], best_vals[1], best_vals[2])

        # print('best_vals: {}'.format(best_vals))
        # print('error: {}'.format(np.sqrt(np.diag(covar))))

        ratio_fit = np.reshape(ratio_fit, ratio.shape)

        error = np.log(np.abs(ratio_fit-ratio)/ratio)
        # max_error = np.max(error)
        # min_error = np.min(error)
        # if np.isinf(min_error):
        #     min
        # print(min_error)


        im = ax.contourf(L0, L1, error, 100, origin='lower',\
            extent=[l0min,l0max,l1min,l1max], cmap='cool_r')

        # Colorbar information
        cbar = ax.cax.colorbar(im, format='%0.1f')
        # cbar.ax.set_ylim(min_error, max_error)
        clims = np.array(im.get_clim())
        cbar.ax.set_yticks(clims)
        cbar.ax.set_yticklabels(['',''])
        if i==3:
            ax.cax.set_ylabel('log(error) [ops]')
        ax.text(0.52, 0.51, '%0.1f'%clims[0], ha='left', va='bottom')
        ax.text(0.52,-0.015, '%0.1f'%clims[1], ha='left', va='top')


        # Label inforrmation
        ax.set_xticks(np.linspace(1/8,3/8,3))
        ax.set_yticks(np.linspace(1/8,3/8,3))
        ax.set_xticklabels(['1/8','2/8','3/8'])
        ax.set_yticklabels([]) # add text later
        if i==0:
            ax.set_xlabel(r'$p_1$')
            ax.set_ylabel(r'$p_1$', labelpad=20)
            ax.text(-0.02, 1/8, '1/8', ha='right', va='center')
            ax.text(-0.02, 2/8, '2/8', ha='right', va='center')
            ax.text(-0.02, 3/8, '3/8', ha='right', va='center')
        if i==2:
            ax.set_xlabel(r'$p_1$')
            ax.set_ylabel(r'$p_1$', labelpad=-3.5)
        elif i==1 or i==3:
            ax.set_xlabel(r'$p_1 = p_0^m$')
            ax.set_ylabel(r'$p_1 = p_0^m$', labelpad=-3.5)

        

        # Get data along each compression ratio line
        # f = rgi((l2, l1), ratio, method='linear')
        # fr6  = f(xyr6)
        # fr8  = f(xyr8)
        # fr10 = f(xyr10)
        # fr12 = f(xyr12)

        # ax.plot(xr6,  fr6,  'm')
        # ax.plot(xr8,  fr8,  'b')
        # ax.plot(xr10, fr10, 'g')
        # ax.plot(xr12, fr12, 'y')
        # ax.grid('on')



        # ax.set_xticks([ntmin,(ntmax+ntmin)/2,ntmax])
        # ax.set_yticks([l1min,(l1max+l1min)/2,l1max])
        # ax.set_xticks(np.linspace(1/8,3/8,3))
        # ax.set_yticks(np.linspace(1/8,3/8,3))
        # ax.set_xticklabels(['1/8','2/8','3/8'])
        # ax.set_yticklabels(['1/8','2/8','3/8'])
        # ax.set_xlabel('$\ell_2$')
        # ax.set_ylabel('$\ell_1$')
        # ax.set_xlim(0,0.5)


    # print('saving image ...')
    fig.set_size_inches(6.5,3.2,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'fig4.png', dpi=300)

    

    print('\tdone with figure 4')