import sys
import os
import yt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from Fig1_reshape   import Fig1_reshape
from Fig6_synthetic import Fig6_synthetic
from Fig7_p1_nt     import Fig7_p1_nt
from Fig8_p1_nt_3D  import Fig8_p1_nt_3D
from Fig9_p1m_p0m   import Fig9_p1m_p0m
from Fig10_l2_error import Fig10_l2_error

from Tab2_fit2D import Tab2_fit2D
from Tab3_fit3D import Tab3_fit3D

from Fig11_genuine import Fig11_genuine
from Fig12_pl_nt   import Fig12_pl_nt
from Fig13_CPU     import Fig13_CPU

print("python version:     {}".format(sys.version))
print("matplotlib version: {}".format(mpl.__version__))
print("numpy version:      {}".format(np.__version__))
print("scipy version:      {}".format(sp.__version__))
print("yt version:         {}".format(yt.__version__))

basedir = '../../'
datadir = basedir + 'data/'
imgdir  = basedir + 'images/'

Fig1_reshape  (imgdir)
Fig6_synthetic(datadir, imgdir)
Fig7_p1_nt    (datadir, imgdir)
Fig8_p1_nt_3D (datadir, imgdir)
Fig9_p1m_p0m  (datadir, imgdir)
Fig10_l2_error(datadir, imgdir)

Tab2_fit2D(datadir, imgdir)
Tab3_fit3D(datadir, imgdir)

Fig11_genuine(datadir, imgdir)
Fig12_pl_nt  (datadir, imgdir)
Fig13_CPU    (datadir, imgdir)



