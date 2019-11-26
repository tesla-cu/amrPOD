import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from Fig1 import Fig1
from Fig2 import Fig2
from Fig3 import Fig3

print("Python version:     {}".format(sys.version))
print("matplotlib version: {}".format(mpl.__version__))
print("numpy version:      {}".format(np.__version__))

basedir = '../../'
datadir = basedir + 'data/'
imgdir  = basedir + 'images/'

# mpl_fname = mpl.matplotlib_fname()

# print('Setting up rcParams')

Fig1(datadir, imgdir)
Fig2(datadir, imgdir)
Fig3(datadir, imgdir)



