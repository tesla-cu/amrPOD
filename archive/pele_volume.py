import sys
import yt
import numpy as np
import os
import glob
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

yt.enable_parallelism()

if yt.is_root():
    print('starting python script to extract slice data')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set up
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get the data directory
basedir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/'
os.chdir(basedir)
# basedir = os.getcwd() 

# Get image directory to save images
# imagedir = '/Users/mike/Research/InternalResearchPapers/AMR_POD/images/'

# Parameters set by user for extraction
starttime = 39.999   # simulation start time
endtime = 40.0201     # simulation end time
nskip = 1         # reduce number of images by factor
finest = None      # finest AMR level, either specify level or None for auto
# region of slice: xlo, xhi, ylo, yhi, zlo, zhi
region = np.array([-0.5, 0.5, -0.25, 0.25, 0.0, 1.0]) # if 2D plane, do NOT choose exact edge
# region = np.array([-1.0, 1.0, -1.0, 1.0, 0.0, 2.0]) # if 2D plane, do NOT choose exact edge
variables = ['z_velocity'] # variables for extraction
yt.funcs.mylog.setLevel(30) # eliminate output from yt load

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Find index of the start and end times
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create the array of plt directories and time of plt file
all_plt_files = glob.glob("plt?????")
all_ts = yt.DatasetSeries(all_plt_files, parallel=True)

# Initialize storage array for all times and plt files
time_storage = {}

# Iterate through plt's to find time indices
# this method leaves None's and index numbers for plt's not in time range
for sto, ds in all_ts.piter(storage=time_storage):
    if ds.current_time >= starttime and ds.current_time <= endtime:
        sto.result = float(ds.current_time)
        sto.result_id = str(ds)

# Convert the storage dictionary values to time and plt file names
time1   = np.array(list(time_storage.values()))
time    = [x for x in time1 if x is not None]
numplt1 = np.array(list(time_storage.keys()))
numplt  = [x for x in numplt1 if x.startswith("plt")]

# Sort these
numplt.sort()
time.sort()

# Reduce these to number we are skipping
numplt = numplt[0::nskip]
time   = time[0::nskip]

# Print out plt number and exact starting time
startplt  = numplt[0]
endplt    = numplt[-1]
starttime = time[0]
endtime   = time[-1]
if yt.is_root():
    print('start plt:  %s'    % startplt)
    print('end plt:    %s'    % endplt)
    print('start time: %0.6f' % starttime)
    print('end time:   %0.6f' % endtime)

# Define number of time steps and a string to the time
nt     = len(numplt)
ttext  = 't%0.4f-%0.4f' % (time[0], time[-1])
dt     = (endtime-starttime)/(nt-1)

# Load arbitrary time step to aquire info about the data
ds = yt.load(numplt[-1])

# Determine finest if set to None
if finest == None:
    reg = ds.r[:, :, :]
    finest = int(reg.max('grid_level'))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Redefine bounds of region  
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define domain xlo, xhi, ..., nx, ...
dxlo = float(ds.domain_left_edge[0])
dxhi = float(ds.domain_right_edge[0])
dylo = float(ds.domain_left_edge[1])
dyhi = float(ds.domain_right_edge[1])
dzlo = float(ds.domain_left_edge[2])
dzhi = float(ds.domain_right_edge[2])
dnx  = int(ds.domain_dimensions[0])*(2**finest)
dny  = int(ds.domain_dimensions[1])*(2**finest)
dnz  = int(ds.domain_dimensions[2])*(2**finest)

# Define user specified bounds
xlo = region[0]
xhi = region[1]
ylo = region[2]
yhi = region[3]
zlo = region[4]
zhi = region[5]

# Redefine coordinates if user is above the domain
if xlo < dxlo:
    xlo = dxlo
if xhi > dxhi:
    xhi = dxhi
if ylo < dylo:
    ylo = dylo
if yhi > dyhi:
    yhi = dyhi
if zlo < dzlo:
    zlo = dzlo
if zhi > dzhi:
    zhi = dzhi

# Define text for saving
xtext = 'x%0.3f-%0.3f' % (xlo, xhi)
ytext = 'y%0.3f-%0.3f' % (ylo, yhi)
ztext = 'z%0.3f-%0.3f' % (zlo, zhi)

# Define folder to save data
volumedir = basedir + 'volume/'
coor_volumedir = '%s_%s_%s_%s_f%i/' % (xtext, ytext, ztext, ttext, finest)

# Directory to store all slice data
if not os.path.exists(volumedir):
    os.mkdir(volumedir)

# Directory for specific region and time frame
if not os.path.exists(volumedir + coor_volumedir):
    os.mkdir(volumedir + coor_volumedir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Specify number of grid points 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Specifying resolution
reg = ds.r[:, :, :]
if xlo != xhi:
    dx = (dxhi - dxlo)/dnx
    nx = int(dnx*(xhi - xlo)/(dxhi - dxlo))
    if yt.is_root():
        print('    dx = %0.11f' % dx, flush=True)
        print('    nx = %i'     % nx, flush=True)
elif xlo==xhi:
    nx = 1

if ylo != yhi:
    dy = (dyhi - dylo)/dny
    ny = int(dny*(yhi - ylo)/(dyhi - dylo))
    if yt.is_root():
        print('    dy = %0.11f' % dy, flush=True)
        print('    ny = %i'     % ny, flush=True)
elif ylo==yhi:
    ny = 1

if zlo != zhi:
    dz = (dzhi - dzlo)/dnz
    nz = int(dnz*(zhi - zlo)/(dzhi - dzlo))
    if yt.is_root():
        print('    dz = %0.11f' % dz, flush=True)
        print('    nz = %i'     % nz, flush=True)
elif zlo==zhi:
    nz = 1

del reg

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Specify x, y, z coordinates
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Identify if the slice is along a particular dimension, and specify dimensions
slice_dim = -1
if xlo == xhi:
    slice_dim = 0
elif ylo == yhi:
    slice_dim = 1
elif zlo == zhi:
    slice_dim = 2
else:
    xf = np.linspace(dxlo+(dx/2),dxhi-(dx/2),dnx)
    yf = np.linspace(dylo+(dy/2),dyhi-(dy/2),dny)
    zf = np.linspace(dzlo+(dz/2),dzhi-(dz/2),dnz)

if yt.is_root():
    if xlo == xhi:
        x    = xlo
        y    = np.linspace(ylo+(dy/2), yhi-(dy/2), ny)
        z    = np.linspace(zlo+(dz/2), zhi-(dz/2), nz)
        x    = np.full((ny,nz), x)
        y, z = np.meshgrid(y, z, indexing='ij')
    elif ylo == yhi:
        x    = np.linspace(xlo+(dx/2), xhi-(dx/2), nx)
        y    = ylo
        z    = np.linspace(zlo+(dz/2), zhi-(dz/2), nz)
        y    = np.full((nx,nz), y)
        x, z = np.meshgrid(x, z, indexing='ij')
    elif zlo == zhi:
        x    = np.linspace(xlo+(dx/2), xhi-(dx/2), nx)
        y    = np.linspace(ylo+(dy/2), yhi-(dy/2), ny)
        z    = zlo
        z    = np.full((nx,ny), z)
        x, y = np.meshgrid(x, y, indexing='ij')
    else:
        x       = np.linspace(xlo+(dx/2),xhi-(dx/2),nx)
        y       = np.linspace(ylo+(dy/2),yhi-(dy/2),ny)
        z       = np.linspace(zlo+(dz/2),zhi-(dz/2),nz)
        x, y, z = np.meshgrid(x, y, z, indexing='ij')


    fh = open(volumedir + coor_volumedir + 'x.bin', "bw")
    x.tofile(fh)
    fh = open(volumedir + coor_volumedir + 'y.bin', "bw")
    y.tofile(fh)
    fh = open(volumedir + coor_volumedir + 'z.bin', "bw")
    z.tofile(fh)

    del x, y, z

# Reassign nx, ny, nz to one array (makes it easier later)
nx = np.array([nx, ny, nz])
# nx = np.array([1,64,64])

# Output slice information
if yt.is_root():
    print('\nregion for extraction is ...', flush=True)
    print('    xlo = %0.3f' % xlo, flush=True)
    print('    xhi = %0.3f' % xhi, flush=True)
    print('    ylo = %0.3f' % ylo, flush=True)
    print('    yhi = %0.3f' % yhi, flush=True)
    print('    zlo = %0.3f' % zlo, flush=True)
    print('    zhi = %0.3f' % zhi, flush=True)
    print('\nnumber of cells in each direction is ...', flush=True)
    print('    nx = %i' % nx[0], flush=True)
    print('    ny = %i' % nx[1], flush=True)
    print('    nz = %i' % nx[2], flush=True)
    print('\nvariables to extract ...', flush=True)
    # for var in variables:
        # print('    %s' % var, flush=True)
    print('\nthe slice will be along the dimension: %i' % slice_dim, flush=True)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract data from region  
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create time series
ts = yt.DatasetSeries(numplt, parallel=True)

# Iterate through all plt files 
for ds in ts.piter():
    
    count = np.argmin(abs(np.array(time)-float(ds.current_time)))

    print('time: %0.6f sec' % ds.current_time, flush=True)

    # Create yt data using ds.covering grid
    # dd = ds.smoothed_covering_grid(level=finest, left_edge=[dxlo,dylo,dzlo], dims=ds.domain_dimensions*(2**finest))
    dd = ds.covering_grid(level=finest, left_edge=[dxlo,dylo,dzlo], dims=ds.domain_dimensions*(2**finest))
    # dd = ds.r[xlo:xhi:complex(0,nx[0]), ylo:yhi:complex(0,nx[1]), zlo:zhi:complex(0,nx[2])]
    if slice_dim == 0:
        ixloc = int(dnx*(2**finest)/2)
    elif slice_dim == 1:
        iyloc = int(dny*(2**finest)/2)
    elif slice_dim == 2:
        izloc = int(dnz*(2**finest)/2)

    # Find indices of coordinates if not on boundary edges
    if xlo > dxlo: ixlo = np.argmin(abs(xf-(dx/2)-xlo))
    else:          ixlo = 0
    if xhi < dxhi: ixhi = np.argmin(abs(xf-(dx/2)-xhi))
    else:          ixhi = dnx
    if ylo > dylo: iylo = np.argmin(abs(yf-(dy/2)-ylo))
    else:          iylo = 0
    if yhi < dyhi: iyhi = np.argmin(abs(yf-(dy/2)-yhi))
    else:          iyhi = dny
    if zlo > dzlo: izlo = np.argmin(abs(zf-(dz/2)-zlo))
    else:          izlo = 0
    if zhi < dzhi: izhi = np.argmin(abs(zf-(dz/2)-zhi))
    else:          izhi = dnz    

    # Iterate through all variables to extract data
    for var in variables:

        # Pull out variable data one at a time, ds.r method
        if slice_dim == -1 or slice_dim == 1:
            varmat = dd[var]
        else:
            varmat = np.transpose(dd[var])

        varmat = np.array(varmat[ixlo:ixhi, iylo:iyhi, izlo:izhi])

        # Save some images on first iteration to see what region we are looking at
        if count == 0 and slice_dim >= 0:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.imshow(np.array(varmat), origin='lower')
            fig.savefig(volumedir + coor_volumedir + var + str(count).zfill(5) + '.png')
        
        # Save data
        fh = open(volumedir + coor_volumedir + var + str(count).zfill(5) + '.bin', "bw")
        varmat.tofile(fh)

        print(np.max(varmat), flush=True)
        print(np.min(varmat), flush=True)

        del varmat

    ds.index.clear_all_data()


# Now we want to move all the files we created into its own folder
if yt.is_root():
    print('finished python script to extract slice data', flush=True)






