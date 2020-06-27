
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import interpolate

from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.optimize import curve_fit



# =========================================================================== #
# Equation to fit:
#
#   l0_coeff*l0^l0_exp + l1_coeff*l1^l1_exp + l2_coeff*l2^l2_exp
# =========================================================================== #
def fit_equation(l, l0_coeff, l1_coeff, l2_coeff, l3_coeff, \
                    l0_exp,   l1_exp,   l2_exp,   l3_exp):
    l0,l1,l2,l3 = l
    return l0_coeff*l0**l0_exp + \
           l1_coeff*l1**l1_exp + \
           l2_coeff*l2**l2_exp + \
           l3_coeff*l3**l3_exp

# =========================================================================== #
# Equation to fit:
#
#   l0_coeff*l0 + l1_coeff*l1 + l2_coeff*l2
# =========================================================================== #
def fit_equation2(l, l0_coeff, l1_coeff, l2_coeff, l3_coeff):
    l0,l1,l2,l3 = l
    return l0_coeff*l0 + \
           l1_coeff*l1 + \
           l2_coeff*l2 + \
           l3_coeff*l3

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0sq_coeff*l0^2  + l1sq_coeff*l1^2  + l2sq_coeff*l2^2 + l3sq_coeff*l3^2 + #
#   l0l1_coeff*l0*l1 + l0l2_coeff*l0*l2 + l1l2_coeff*l2^2                     #
# =========================================================================== #
def fit_equation3(l, l0sq_coeff, l1sq_coeff, l2sq_coeff, l3sq_coeff, \
                     l0l1_coeff, l0l2_coeff, l0l3_coeff, \
                     l1l2_coeff, l1l3_coeff, l2l3_coeff):
    l0,l1,l2,l3 = l
    return l0sq_coeff*l0**2 + l1sq_coeff*l1**2 + l2sq_coeff*l2**2 + l3sq_coeff*l3**2 + \
           l0l1_coeff*l0*l1 + l0l2_coeff*l0*l2 + l0l3_coeff*l0*l3 + \
           l1l2_coeff*l1*l2 + l1l3_coeff*l1*l3 + l2l3_coeff*l2*l3


def find_reduced_dim_3D(datadir, imgdir):

    print('doing data reduced for 3D ...')

    # Set up figure
    fig = plt.figure()

    # Load data
    npydir = datadir + 'l1_l2_l3_3D/npy_files/'

    TC_R_avg_imp   = np.load(npydir + "/TC_R_avg_imp.npy")
    TC_R_avg_unalt = np.load(npydir + "/TC_R_avg_unalt.npy")
    TC_R_rms_imp   = np.load(npydir + "/TC_R_rms_imp.npy")
    TC_R_rms_unalt = np.load(npydir + "/TC_R_rms_unalt.npy")

    TC_P1_avg_imp   = np.load(npydir + "/TC_P1_avg_imp.npy")
    TC_P1_avg_unalt = np.load(npydir + "/TC_P1_avg_unalt.npy")
    TC_P1_rms_imp   = np.load(npydir + "/TC_P1_rms_imp.npy")
    TC_P1_rms_unalt = np.load(npydir + "/TC_P1_rms_unalt.npy")

    TC_P2_avg_imp   = np.load(npydir + "/TC_P2_avg_imp.npy")
    TC_P2_avg_unalt = np.load(npydir + "/TC_P2_avg_unalt.npy")
    TC_P2_rms_imp   = np.load(npydir + "/TC_P2_rms_imp.npy")
    TC_P2_rms_unalt = np.load(npydir + "/TC_P2_rms_unalt.npy")

    TC_A_avg_imp   = np.load(npydir + "/TC_A_avg_imp.npy")
    TC_A_avg_unalt = np.load(npydir + "/TC_A_avg_unalt.npy")
    TC_A_rms_imp   = np.load(npydir + "/TC_A_rms_imp.npy")
    TC_A_rms_unalt = np.load(npydir + "/TC_A_rms_unalt.npy")

    # Data information
    l1min = 0/64  # x min
    l1max = 16/64 # x max
    l1inc = 8/64  # x increment
    l2min = 0/64  # y min
    l2max = 16/64 # y max
    l2inc = 8/64  # y increment
    l3min = 0/64  # y min
    l3max = 16/64 # y max
    l3inc = 8/64  # y increment

    l1 = np.arange(l1min, l1max+l1inc, l1inc)
    l2 = np.arange(l2min, l2max+l2inc, l2inc)
    l3 = np.arange(l3min, l3max+l3inc, l3inc)

    L1, L2, L3 = np.meshgrid(l1, l2, l3, indexing='ij')
    L0 = 1.0 - L1 - L2 - L3

    for i in range(4):

        if i == 0:
            num = TC_R_avg_imp
            den = TC_R_avg_unalt
        elif i == 1:
            num = TC_P1_avg_imp
            den = TC_P1_avg_unalt
        elif i == 2:
            num = TC_P2_avg_imp
            den = TC_P2_avg_unalt
        elif i == 3:
            num = TC_A_avg_imp
            den = TC_A_avg_unalt

        ratio = num/den
        # ratio = np.exp(num/den)
        # ratio = np.log10(np.sqrt(num/den))
        print(np.max(ratio))
        print(np.min(ratio))


        l0_1D    = np.reshape(L0,    (-1))
        l1_1D    = np.reshape(L1,    (-1))
        l2_1D    = np.reshape(L2,    (-1))
        l3_1D    = np.reshape(L3,    (-1))
        ratio_1D = np.reshape(ratio, (-1))


        best_vals, covar = curve_fit(fit_equation, (l0_1D,l1_1D,l2_1D,l3_1D), \
            ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------')

        best_vals, covar = curve_fit(fit_equation2, (l0_1D,l1_1D,l2_1D,l3_1D), \
            ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------')

        best_vals, covar = curve_fit(fit_equation3, (l0_1D,l1_1D,l2_1D,l3_1D),\
            ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('==============================================================')

