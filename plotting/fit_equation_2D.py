
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0_coeff*l0^l0_exp + l1_coeff*l1^l1_exp + l2_coeff*l2^l2_exp              #
# =========================================================================== #
def fit_exponential_equation(l, l0_coeff, l1_coeff, l2_coeff, \
                                l0_exp,   l1_exp,   l2_exp):
    l0,l1,l2 = l
    return l0_coeff*l0**l0_exp + \
           l1_coeff*l1**l1_exp + \
           l2_coeff*l2**l2_exp

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0_coeff*l0 + l1_coeff*l1 + l2_coeff*l2                                   #
# =========================================================================== #
def fit_linear_equation(l, l0_coeff, l1_coeff, l2_coeff):
    l0,l1,l2 = l
    return l0_coeff*l0 + \
           l1_coeff*l1 + \
           l2_coeff*l2

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0sq_coeff*l0^2  + l1sq_coeff*l1^2  + l2sq_coeff*l2^2                     #
#   l0l1_coeff*l0*l1 + l0l2_coeff*l0*l2 + l1l2_coeff*l2^2                     #
# =========================================================================== #
def fit_quadratic_equation(l, l0sq_coeff, l1sq_coeff, l2sq_coeff, \
                              l0l1_coeff, l0l2_coeff, l1l2_coeff):
    l0,l1,l2 = l
    return   l0sq_coeff*l0**2 +   l1sq_coeff*l1**2 +   l2sq_coeff*l2**2 + \
           2*l0l1_coeff*l0*l1 + 2*l0l2_coeff*l0*l2 + 2*l1l2_coeff*l1*l2


def fit_equation_2D(datadir, imgdir):

    print('fitting equations for 2D data ...')

    # Set up figure
    fig = plt.figure()

    # Load data
    txtdir_l  = datadir + 'l1_l0/txt_files/'
    txtdir_lc = datadir + 'lc1l1_lc0l0/txt_files/'

    TC_R_avg_imp   = np.loadtxt(txtdir_l + "/TC_R_avg_imp.txt")
    TC_R_avg_unalt = np.loadtxt(txtdir_l + "/TC_R_avg_unalt.txt")
    TC_R_rms_imp   = np.loadtxt(txtdir_l + "/TC_R_rms_imp.txt")
    TC_R_rms_unalt = np.loadtxt(txtdir_l + "/TC_R_rms_unalt.txt")

    TC_P1_avg_imp   = np.loadtxt(txtdir_lc + "/TC_P1_avg_imp.txt")
    TC_P1_avg_unalt = np.loadtxt(txtdir_lc + "/TC_P1_avg_unalt.txt")
    TC_P1_rms_imp   = np.loadtxt(txtdir_lc + "/TC_P1_rms_imp.txt")
    TC_P1_rms_unalt = np.loadtxt(txtdir_lc + "/TC_P1_rms_unalt.txt")

    TC_P2_avg_imp   = np.loadtxt(txtdir_l + "/TC_P2_avg_imp.txt")
    TC_P2_avg_unalt = np.loadtxt(txtdir_l + "/TC_P2_avg_unalt.txt")
    TC_P2_rms_imp   = np.loadtxt(txtdir_l + "/TC_P2_rms_imp.txt")
    TC_P2_rms_unalt = np.loadtxt(txtdir_l + "/TC_P2_rms_unalt.txt")

    TC_A_avg_imp   = np.loadtxt(txtdir_lc + "/TC_A_avg_imp.txt")
    TC_A_avg_unalt = np.loadtxt(txtdir_lc + "/TC_A_avg_unalt.txt")
    TC_A_rms_imp   = np.loadtxt(txtdir_lc + "/TC_A_rms_imp.txt")
    TC_A_rms_unalt = np.loadtxt(txtdir_lc + "/TC_A_rms_unalt.txt")

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
    L2 = 1 - L1 - L0

    # File to store best fit values
    best_fit_vals_2D = open(imgdir + '/best_fit_vals_2D.txt', 'w')

    print('===================================================' \
          '===================================================')
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
        ratio_1D = np.reshape(ratio, (-1))

        """
        best_vals, covar = curve_fit(fit_exponential_equation, \
            (l0_1D,l1_1D,l2_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------------------')

        best_vals, covar = curve_fit(fit_linear_equation, \
            (l0_1D,l1_1D,l2_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------------------')

        best_vals, covar = curve_fit(fit_quadratic_equation, \
            (l0_1D,l1_1D,l2_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('===================================================' \
              '===================================================')
        """

        if i==0:
            best_vals, covar = curve_fit(fit_quadratic_equation, \
                (l0_1D,l1_1D,l2_1D), ratio_1D)
            print('best_vals: {}'.format(best_vals))
            print('error: {}'.format(np.sqrt(np.diag(covar))))

        else:
            best_vals, covar = curve_fit(fit_linear_equation, \
                (l0_1D,l1_1D,l2_1D), ratio_1D)
            print('best_vals: {}'.format(best_vals))
            print('error: {}'.format(np.sqrt(np.diag(covar))))

        print('===================================================' \
              '===================================================')

        if   i==0: best_fit_vals_2D.write("R info:\n")
        elif i==1: best_fit_vals_2D.write("Phi - M1 info:\n")
        elif i==2: best_fit_vals_2D.write("Phi - M2 info:\n")
        elif i==3: best_fit_vals_2D.write("A info:\n")

        best_fit_vals_2D.write("    best values: ")
        [best_fit_vals_2D.write("%0.8f " % b) for b in best_vals]
        best_fit_vals_2D.write("\n    one standard deviation error: ")
        [best_fit_vals_2D.write("%0.8e " % e) for e in np.sqrt(np.diag(covar))]
        best_fit_vals_2D.write("\n\n")

    best_fit_vals_2D.close()

