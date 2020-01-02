import sys
import os
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from Fig0 import Fig0
from Fig1 import Fig1
from Fig2 import Fig2
from Fig3 import Fig3
from Fig4 import Fig4
from Fig5 import Fig5

print("Python version:     {}".format(sys.version))
print("matplotlib version: {}".format(mpl.__version__))
print("numpy version:      {}".format(np.__version__))
print("scipy version:      {}".format(sp.__version__))
print("yt version:         {}".format(yt.__version__))

basedir = '../../'
datadir = basedir + 'data/'
imgdir  = basedir + 'images/'

# mpl_fname = mpl.matplotlib_fname()

# print('Setting up rcParams')

# Fig0(datadir, imgdir)
# Fig1(datadir, imgdir)
# Fig2(datadir, imgdir)
# Fig3(datadir, imgdir)
Fig4(datadir, imgdir)
# Fig5(datadir, imgdir)



