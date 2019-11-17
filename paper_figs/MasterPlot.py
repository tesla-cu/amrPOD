import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np



from Fig1 import Fig1

print("Python version:     {}".format(sys.version))
print("matplotlib version: {}".format(mpl.__version__))
print("numpy version:      {}".format(np.__version__))



basedir = '/home/mike/Research/POD_AMR/'
datadir = basedir + 'data/'
imgdir  = basedir + 'images/'

# mpl_fname = mpl.matplotlib_fname()

# print('Setting up rcParams')

Fig1(datadir, imgdir)


