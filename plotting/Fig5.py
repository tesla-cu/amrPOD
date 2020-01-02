import numpy as np
import yt
import matplotlib
import matplotlib.pyplot as plt
 
from mpl_toolkits.axes_grid1 import AxesGrid

def Fig5(datadir, imgdir):

    yt.funcs.mylog.setLevel(50)

    print('making figure 5 ...')

    # Set up figure
    fig = plt.figure()

    # Directories storing simulation and synthetic data
    simdir = datadir + 'AMR_sim/'

    # Units from PeleLM simulation
    units_override = {"length_unit": (1.0, "m"),
                      "time_unit": (1.0, "s"),
                      "mass_unit": (1.0, "kg"),
                      "temperature_unit": (1.0, "K")}

    # Generate discrete colormap based on 'Accent' for grid level images
    cmap_acc = matplotlib.cm.get_cmap('Accent')
    cmap = matplotlib.colors.ListedColormap(\
        [cmap_acc(0.125), cmap_acc(0.375), cmap_acc(0.625), cmap_acc(0.875)])

    # Top half of figure ------------------------------------------------------
    grid1 = AxesGrid(fig, (0.055,0.53,0.88,0.40),
        nrows_ncols = (1, 5),
        axes_pad = 0.1,
        aspect = True,
        label_mode = "L",
        share_all = True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="7%",
        cbar_pad="3%")

    grid2 = AxesGrid(fig, (0.055,0.11,0.88,0.40),
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

        if i == 0:
            ds = yt.load(simdir + 'plt27130', units_override=units_override)
        elif i == 1:
            ds = yt.load(simdir + 'plt27130', units_override=units_override)
        elif i == 2:
            ds = yt.load(simdir + 'plt27130', units_override=units_override)
        elif i == 3:
            ds = yt.load(simdir + 'plt27130', units_override=units_override)
        elif i == 4:
            ds = yt.load(simdir + 'plt27130', units_override=units_override)

        # Top half of figure --------------------------------------------------

        sd = yt.SlicePlot(ds,'x','density', origin="native", \
            center=((0.0,0.0,0.375),'m'), width=(0.75,'m'))
        sd.set_cmap(field='density', cmap='dusk')
        sd.annotate_grids(linewidth=0.25, min_level=1)
        sd.set_axes_unit('cm')
        sd.set_xlabel('')
        sd.set_ylabel('$z$ [cm]')
        sd.set_font_size(8)
        sd.set_minorticks('all', 'off')
        sd.set_unit('density','kg/m**3')
        sd.set_log('density',False)
        sd.set_colorbar_label(field='density',label='Density [kg/m$^3$]')

        # Redraw on combined figure
        plot = sd.plots['density']
        plot.figure = fig
        plot.axes = grid1[i].axes
        plot.cax = grid1.cbar_axes[i]
        sd._setup_plots()

        # Set details after redrawing
        sd.plots['density'].axes.set_xticks([-30,-15,0,15,30])
        sd.plots['density'].axes.set_xticklabels([])
        sd.plots['density'].axes.set_yticks([0,15,30,45,60])
        if i == 0:
            sd.plots['density'].axes.set_title(r'$t=0$')
        elif i == 1:
            sd.plots['density'].axes.set_title(r'$t=\tau/4$')
        elif i == 2:
            sd.plots['density'].axes.set_title(r'$t=\tau/2$')
        elif i == 3:
            sd.plots['density'].axes.set_title(r'$t=3\tau/4$')
        elif i == 4:
            sd.plots['density'].axes.set_title(r'$t=\tau$')
        sd.plots['density'].cb.set_ticks([0.2,0.5,0.8,1.1])
        sd.plots['density'].cb.ax.tick_params(length=2, direction='in')

        # Bottom half of figure -----------------------------------------------
        sg = yt.SlicePlot(ds,'x','grid_level', origin="native", \
            center=((0.0,0.0,0.375),'m'), width=(0.75,'m'))
        sg.set_cmap(field='grid_level', cmap=cmap)
        sg.set_axes_unit('cm')
        sg.set_xlabel('$x$ [cm]')
        sg.set_ylabel('$z$ [cm]')
        sg.set_font_size(8)
        sg.set_minorticks('all', 'off')
        sg.set_log('grid_level',False)
        sg.set_zlim('grid_level', 0.5, 4.5)
        sg.set_colorbar_label(field='grid_level',label='$\ell$')

        # Redraw on combined figure
        plot = sg.plots['grid_level']
        plot.figure = fig
        plot.axes = grid2[i].axes
        plot.cax = grid2.cbar_axes[i]
        sg._setup_plots()

        # Set details after redrawing
        sg.plots['grid_level'].axes.set_xticks([-30,-15,0,15,30])
        sg.plots['grid_level'].axes.set_xticklabels(\
            ['-30','-15','0','15','30'])
        sg.plots['grid_level'].axes.set_yticks([0,15,30,45,60])
        # sg.plots['grid_level'].cb.set_ticks(np.linspace(0.5,3.5,4))
        sg.plots['grid_level'].cb.set_ticks(np.linspace(1,4,4))
        sg.plots['grid_level'].cb.set_ticklabels(['0','1','2','3'])
        sg.plots['grid_level'].cb.ax.tick_params(length=2, direction='in')
        sg.plots['grid_level'].cb.ax.set_ylim(0.5,4.5)

    # Save figure
    fig.set_size_inches(6.5,2.75,forward=True)
    plt.savefig(imgdir + 'fig5.png', dpi=300)
    

    print('\tdone with figure 5')