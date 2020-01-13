import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0_coeff*l0^l0_exp + l1_coeff*l1^l1_exp + l2_coeff*l2^l2_exp +            #
#   l3_coeff*l3^l3_exp                                                        #
# =========================================================================== #
def fit_exponential_equation(l, l0_coeff, l1_coeff, l2_coeff, l3_coeff, \
                                l0_exp,   l1_exp,   l2_exp,   l3_exp):
    l0,l1,l2,l3 = l
    return l0_coeff*l0**l0_exp + \
           l1_coeff*l1**l1_exp + \
           l2_coeff*l2**l2_exp + \
           l3_coeff*l3**l3_exp

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0_coeff*l0 + l1_coeff*l1 + l2_coeff*l2 + l3_coeff*l3                     #
# =========================================================================== #
def fit_linear_equation(l, l0_coeff, l1_coeff, l2_coeff, l3_coeff):
    l0,l1,l2,l3 = l
    return l0_coeff*l0 + \
           l1_coeff*l1 + \
           l2_coeff*l2 + \
           l3_coeff*l3

# =========================================================================== #
# Equation to fit:                                                            #
#                                                                             #
#   l0sq_coeff*l0^2  + l1sq_coeff*l1^2  + l2sq_coeff*l2^2 + l3sq_coeff*l3^2 + #
#   2*l0l1_coeff*l0*l1 + 2*l0l2_coeff*l0*l2 + 2*l1l2_coeff*l2^2 +             #
#   2*l0l1_coeff*l0*l1 + 2*l0l2_coeff*l0*l2 + 2*l1l2_coeff*l2^2               #
# =========================================================================== #
def fit_quadratic_equation(l, l0sq_coeff, l1sq_coeff, l2sq_coeff, l3sq_coeff, \
                              l0l1_coeff, l0l2_coeff, l0l3_coeff, \
                              l1l2_coeff, l1l3_coeff, l2l3_coeff):
    l0,l1,l2,l3 = l
    return   l0sq_coeff*l0**2 +   l1sq_coeff*l1**2 +   l2sq_coeff*l2**2 + l3sq_coeff*l3**2 + \
           2*l0l1_coeff*l0*l1 + 2*l0l2_coeff*l0*l2 + 2*l0l3_coeff*l0*l3 + \
           2*l1l2_coeff*l1*l2 + 2*l1l3_coeff*l1*l3 + 2*l2l3_coeff*l2*l3


def fit_equation_3D(datadir, imgdir):

    print('fitting equations for 3D data ...')

    # Set up figure
    fig = plt.figure()

    # Load data
    # npydir_l  = datadir + 'l0_l1_l2_3D/npy_files/'
    # npydir_lc = datadir + 'lc0l0_lc1l1_lc2l2_3D/npy_files/'
    # npydir_l  = datadir + 'lc0l0_lc1l1_lc2l2_3D_small24/npy_files/'
    # npydir_lc = datadir + 'lc0l0_lc1l1_lc2l2_3D_small24/npy_files/'
    npydir_l  = datadir + 'lc0l0_lc1l1_lc2l2_3D_small16/npy_files/'
    npydir_lc = datadir + 'lc0l0_lc1l1_lc2l2_3D_small16/npy_files/'

    TC_R_avg_imp   = np.load(npydir_l + "/TC_R_avg_imp.npy")
    TC_R_avg_unalt = np.load(npydir_l + "/TC_R_avg_unalt.npy")
    TC_R_rms_imp   = np.load(npydir_l + "/TC_R_rms_imp.npy")
    TC_R_rms_unalt = np.load(npydir_l + "/TC_R_rms_unalt.npy")

    TC_P1_avg_imp   = np.load(npydir_lc + "/TC_P1_avg_imp.npy")
    TC_P1_avg_unalt = np.load(npydir_lc + "/TC_P1_avg_unalt.npy")
    TC_P1_rms_imp   = np.load(npydir_lc + "/TC_P1_rms_imp.npy")
    TC_P1_rms_unalt = np.load(npydir_lc + "/TC_P1_rms_unalt.npy")

    TC_P2_avg_imp   = np.load(npydir_l + "/TC_P2_avg_imp.npy")
    TC_P2_avg_unalt = np.load(npydir_l + "/TC_P2_avg_unalt.npy")
    TC_P2_rms_imp   = np.load(npydir_l + "/TC_P2_rms_imp.npy")
    TC_P2_rms_unalt = np.load(npydir_l + "/TC_P2_rms_unalt.npy")

    TC_A_avg_imp   = np.load(npydir_lc + "/TC_A_avg_imp.npy")
    TC_A_avg_unalt = np.load(npydir_lc + "/TC_A_avg_unalt.npy")
    TC_A_rms_imp   = np.load(npydir_lc + "/TC_A_rms_imp.npy")
    TC_A_rms_unalt = np.load(npydir_lc + "/TC_A_rms_unalt.npy")

    # Data information
    # l0min = 0/27 # x min
    # l0max = 9/27 # x max
    # l0inc = 1/27 # x increment
    # l1min = 0/27 # y min
    # l1max = 9/27 # y max
    # l1inc = 1/27 # y increment
    # l2min = 0/27 # y min
    # l2max = 9/27 # y max
    # l2inc = 1/27 # y increment

    # l0min = 0/27 # x min
    # l0max = 9/27 # x max
    # l0inc = 3/27 # x increment
    # l1min = 0/27 # y min
    # l1max = 9/27 # y max
    # l1inc = 3/27 # y increment
    # l2min = 0/27 # y min
    # l2max = 9/27 # y max
    # l2inc = 3/27 # y increment

    l0min = 0/64 # x min
    l0max = 20/64 # x max
    l0inc = 4/64 # x increment
    l1min = 0/64 # y min
    l1max = 20/64 # y max
    l1inc = 4/64 # y increment
    l2min = 0/64 # y min
    l2max = 20/64 # y max
    l2inc = 4/64 # y increment

    l0 = np.arange(l0min, l0max+l0inc, l0inc)
    l1 = np.arange(l1min, l1max+l1inc, l1inc)
    l2 = np.arange(l2min, l2max+l2inc, l2inc)

    L0, L1, L2 = np.meshgrid(l0, l1, l2, indexing='ij')
    L3 = 1 - L2 - L1 - L0

    # File to store best fit values
    best_fit_vals_3D = open(imgdir + '/best_fit_vals_3D.txt', 'w')

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
        l3_1D    = np.reshape(L3,    (-1))
        ratio_1D = np.reshape(ratio, (-1))

        """
        best_vals, covar = curve_fit(fit_exponential_equation, \
            (l0_1D,l1_1D,l2_1D,l3_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------------------')

        best_vals, covar = curve_fit(fit_linear_equation, \
            (l0_1D,l1_1D,l2_1D,l3_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('------------------------------------------------------')

        best_vals, covar = curve_fit(fit_quadratic_equation, \
            (l0_1D,l1_1D,l2_1D,l3_1D), ratio_1D)
        print('best_vals: {}'.format(best_vals))
        print('error: {}'.format(np.sqrt(np.diag(covar))))
        print('covar: {}'.format(covar))

        print('===================================================' \
              '===================================================')
        """
        
        
        # """
        if i==0:
            best_vals, covar = curve_fit(fit_quadratic_equation, \
                (l0_1D,l1_1D,l2_1D,l3_1D), ratio_1D)
            print('best_vals: {}'.format(best_vals))
            print('error: {}'.format(np.sqrt(np.diag(covar))))

        else:
            best_vals, covar = curve_fit(fit_linear_equation, \
                (l0_1D,l1_1D,l2_1D,l3_1D), ratio_1D)
            print('best_vals: {}'.format(best_vals))
            print('error: {}'.format(np.sqrt(np.diag(covar))))


        print('===================================================' \
              '===================================================')
        # """


        if   i==0: best_fit_vals_3D.write("R info:\n")
        elif i==1: best_fit_vals_3D.write("Phi - M1 info:\n")
        elif i==2: best_fit_vals_3D.write("Phi - M2 info:\n")
        elif i==3: best_fit_vals_3D.write("A info:\n")

        best_fit_vals_3D.write("    best values: ")
        [best_fit_vals_3D.write("%0.8f " % b) for b in best_vals]
        best_fit_vals_3D.write("\n    one standard deviation error: ")
        [best_fit_vals_3D.write("%0.8e " % e) for e in np.sqrt(np.diag(covar))]
        best_fit_vals_3D.write("\n\n")

        

    best_fit_vals_3D.close()

