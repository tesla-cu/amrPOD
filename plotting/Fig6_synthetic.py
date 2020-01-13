
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import interpolate

from scipy.interpolate import RegularGridInterpolator as rgi

def Fig6_synthetic(datadir, imgdir):

    print('making figure 6 ...')

    # Set up figure
    fig = plt.figure()

    # Directory where synthetic data is stored
    syndir = datadir + 'AMR_syn/'

    # Generate discrete colormap based on 'Accent' for grid level images
    cmap_acc = matplotlib.cm.get_cmap('Accent')
    cmap = matplotlib.colors.ListedColormap(\
        [cmap_acc(0.125), cmap_acc(0.375)])

    # Plot figure -------------------------------------------------------------
    grid = AxesGrid(fig, (0.06,0.07,0.89,0.91),
        nrows_ncols = (1, 5),
        axes_pad = 0.1,
        aspect = True,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")
    

    for i in range(5):

        ax = grid[i].axes

        # Load data
        grid_data = np.reshape(np.fromfile(syndir + 'grid_level%05d.bin' % i),\
            (64,64))

        # Display grid level data
        image = ax.imshow(grid_data.T, origin='lower', extent=[1,64,1,64], \
            cmap=cmap, vmin=0, vmax=2)

        # Colorbar information
        cbar = ax.cax.colorbar(image, format='%i')
        ax.cax.set_ylabel(r'$\ell$')
        cbar.ax.set_ylim(0,2)
        cbar.ax.set_yticks([0.5,1.5])
        cbar.ax.set_yticklabels(['0','1'])

        # Label information
        ax.set_title(r'$t=%i$' % i)
        ax.set_xticks(np.linspace(8,56,4))
        ax.set_yticks(np.linspace(8,56,4))
        ax.set_xticklabels(['8','24','40','56'])
        ax.set_yticklabels(['8','24','40','56'])
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')

    # Save image
    fig.set_size_inches(6.5,1.7,forward=True) # figure size must be set here
    plt.savefig(imgdir + 'Fig6_synthetic.png', dpi=300)

    print('\tdone with figure 6')