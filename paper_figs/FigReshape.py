import numpy as np 
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as patches

def FigReshape(imgdir):

    # Set up figure
    fig = plt.figure()
    # Define grid size in image. 
    # If this changes, colspan and rowspan will need to change.
    gridspec.GridSpec(256, 256)     
    d = .1 # dot diameter for 3 continuation dots
    # define colormap
    cmap_acc = mpl.cm.get_cmap('Accent')
    lvlC0    = np.tile(cmap_acc(.125), (3,1))
    lvlC1    = np.tile(cmap_acc(.375), (3,1))
    lvlC2    = np.tile(cmap_acc(.625), (3,1))
    lvlC3    = np.tile(cmap_acc(.875), (4,1))

    # Adjust color map to include proper number of colors
    for i in range(4):
        if i < 3:
            lvlC0[i,1] = lvlC0[i,1] - (i*.1)
            lvlC1[i,1] = lvlC1[i,1] - (i*.1)
            lvlC2[i,0] = lvlC2[i,0] - (i*.1)
        if i < 4:
            lvlC3[i,0] = lvlC3[i,0] - (i*.1)

    cmap1    = mpl.colors.ListedColormap([lvlC0[0,:], lvlC0[1,:], lvlC0[2,:],\
                                          lvlC1[0,:], lvlC1[1,:], lvlC1[2,:],\
                                          lvlC2[0,:], lvlC2[1,:], lvlC2[2,:],\
                                          lvlC3[0,:], lvlC3[1,:], lvlC3[2,:], lvlC3[3,:]])

    # Initialize X
    X = np.zeros((16, 16))

    # Level zero
    X[0:8,  8:16] = 0/13
    X[0:8,  0:8]  = 1/13
    X[8:16, 0:8]  = 2/13

    # Level one
    X[8:12,  8:12]  = 3/13
    X[12:16, 12:16] = 4/13
    X[8:12,  12:16] = 5/13

    # Level two
    X[14:16, 8:10]  = 6/13
    X[14:16, 10:12] = 7/13
    X[12:14, 8:10]  = 8/13

    # Level three
    X[12:13, 10:11] = 9/13
    X[13:14, 10:11] = 10/13
    X[13:14, 11:12] = 11/13
    X[12:13, 11:12] = 12/13

    # Set up of reshaping procedure
    nx     = 16
    ny     = 16
    nspat  = int(nx*ny)
    finest = 3
    nlev   = finest + 1
    ndim   = 2
    c_l    = np.zeros((nlev), dtype=int)
    for i in range(nlev):
        c_l[i] = 2**(finest-i)
    data_1D = np.squeeze(X)

    plt.subplot2grid((256,256), (0, 0), colspan=64, rowspan=64) 
    plt.imshow(data_1D.T, cmap=cmap1, origin='lower')
    plt.box(False)

    # Plot original grid
    ax = plt.gca()
    ax.set_xticks(np.arange(0, nx, 1))
    ax.set_yticks(np.arange(0, ny, 1))
    ax.set_xticks(np.arange(-.5, nx, 1), minor=True);
    ax.set_yticks(np.arange(-.5, ny, 1), minor=True);
    ax.grid(which='minor', color='k', linestyle='--', linewidth=.25)
    ax.tick_params(labelbottom=False)
    ax.tick_params(labelleft=False)
    ax.tick_params(which='minor', axis='y', colors='white')
    ax.tick_params(which='minor', axis='x', colors='white')
    ax.arrow(16, 7.5, 6, 0, head_width=.5, head_length=1, fc='k', ec='k', clip_on=False)

    nxr = nx

    for k, c in enumerate(c_l):
        
        # Reshape data
        data_1D = np.transpose(data_1D, ( 1,  0))
        data_1D = np.reshape(  data_1D, (-1,  c,  nxr))
        data_1D = np.transpose(data_1D, ( 1,  0,  2))
        data_1D = np.reshape(  data_1D, ( c, -1))

        # Get dimensions of the current reshaped matrix
        nxr = data_1D.shape[0]
        nyr = nspat/nxr
        
        if k == 0:
            plt.subplot2grid((256,256), (0, 80), colspan=32, rowspan=128)
            plt.imshow(data_1D.T, cmap=cmap1, origin='lower')
            plt.box(False)

            ax = plt.gca()
            ax.set_xticks(np.arange(0, nxr, 1))
            ax.set_yticks(np.arange(0, nyr, 1))
            ax.set_xticks(np.arange(-.5, nxr, 1), minor=True);
            ax.set_yticks(np.arange(-.5, nyr, 1), minor=True);
            ax.grid(which='minor', color='k', linestyle='--', linewidth=.25)
            ax.tick_params(labelbottom=False)
            ax.tick_params(labelleft=False)
            ax.tick_params(which='minor', axis='y', colors='white')
            ax.tick_params(which='minor', axis='x', colors='white')
            ax.arrow(8, 23.5, 6, 0, head_width=.5, head_length=1, fc='k', ec='k', clip_on=False)

        if k == 1:
            plt.subplot2grid((256,256), (0, 134), colspan=16, rowspan=200)
            plt.imshow(data_1D.T, cmap=cmap1, origin='lower')
            plt.box(False)

            ax = plt.gca()
            ax.set_xticks(np.arange(0, nxr, 1))
            ax.set_yticks(np.arange(0, nyr, 1))
            ax.set_xticks(np.arange(-.5, nxr, 1), minor=True);
            ax.set_yticks(np.arange(-.5, nyr, 1), minor=True);
            ax.grid(which='minor', color='k', linestyle='--', linewidth=.25)
            ax.tick_params(labelbottom=False)
            ax.tick_params(labelleft=False)
            ax.tick_params(which='minor', axis='y', colors='white')
            ax.tick_params(which='minor', axis='x', colors='white')
            ax.tick_params(which='major', axis='x', colors='white')
            ax.set_ylim((13.5,63.5))
            ax.arrow(4, 55.5, 6, 0, head_width=.5, head_length=1, fc='k', ec='k', clip_on=False)
            for i in range(3):
                circ = plt.Circle((1.5, 12.5-i), d, color='k', clip_on=False)
                ax.add_patch(circ)

        if k == 2:
            plt.subplot2grid((256,256), (0, 175), colspan=8, rowspan=220)
            plt.imshow(data_1D.T, cmap=cmap1, origin='lower')
            plt.box(False)

            ax = plt.gca()
            ax.set_xticks(np.arange(0, nxr, 1))
            ax.set_yticks(np.arange(0, nyr, 1))
            ax.set_xticks(np.arange(-.5, nxr, 1), minor=True);
            ax.set_yticks(np.arange(-.5, nyr, 1), minor=True);
            ax.grid(which='minor', color='k', linestyle='--', linewidth=.15)
            ax.tick_params(labelbottom=False)
            ax.tick_params(labelleft=False)
            ax.tick_params(which='minor', axis='y', colors='white')
            ax.tick_params(which='minor', axis='x', colors='white') 
            ax.tick_params(which='major', axis='x', colors='white')
            ax.set_ylim((72.5,127.5))
            ax.arrow(2, 119.5, 6, 0, head_width=.5, head_length=1, fc='k', ec='k', clip_on=False)
            for i in range(3):
                circ = plt.Circle((.5, 71.5-i), d, color='k', clip_on=False)
                ax.add_patch(circ)      
        
        
        if k==3:

            plt.subplot2grid((256,256), (0, 210), colspan=4, rowspan=260)
            plt.imshow(data_1D.T, cmap=cmap1, origin='lower')
            plt.box(False)

            ax = plt.gca()
            ax.set_xticks(np.arange(-.5, 1, 1))
            ax.set_yticks(np.arange(0, 256.5, 1))
            ax.set_xticks(np.arange(-.5, 1, 1), minor=True);
            ax.set_yticks(np.arange(-.5, 256.5, 1), minor=True);
            ax.grid(which='minor', color='k', linestyle='--', linewidth=.15)
            ax.tick_params(labelbottom=False)
            ax.tick_params(labelleft=False)
            ax.tick_params(which='minor', axis='y', colors='white')
            ax.tick_params(which='minor', axis='x', colors='white') 
            ax.tick_params(which='major', axis='x', colors='white')
            ax.set_ylim((190.5,255.5))
            for i in range(3):
                circ = plt.Circle((0, 189.5-i), d, color='k', clip_on=False)
                ax.add_patch(circ)

    fig.set_size_inches(6.5,3.25,forward=True)
    plt.savefig(imgdir + 'fig_reshape.png', dpi=300)
