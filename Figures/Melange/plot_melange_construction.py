
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def read_grid_components(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                  model_name + '_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    ds.close()

    return (XC, YC, Depth)

def read_iceberg_fields(iceberg_dir, rows, cols, Nr):

    bergMask = np.fromfile(os.path.join(iceberg_dir, 'bergMask.bin'), '>f8').reshape((rows, cols))
    #bergConc = np.fromfile(os.path.join(iceberg_dir, 'bergConc.bin'), '>f8').reshape((rows, cols))
    bergMaskNums = np.fromfile(os.path.join(iceberg_dir, 'bergMaskNums.bin'), '>f8').reshape((rows, cols))
    driftMask = np.fromfile(os.path.join(iceberg_dir, 'driftMask.bin'), '>f8').reshape((rows, cols))
    barrierMask = np.fromfile(os.path.join(iceberg_dir, 'barrierMask.bin'), '>f8').reshape((rows, cols))
    meltMask = np.fromfile(os.path.join(iceberg_dir, 'meltMask.bin'), '>f8').reshape((rows, cols))
    numBergsPerCell = np.fromfile(os.path.join(iceberg_dir, 'numBergsPerCell.bin'), '>f8').reshape((rows, cols))
    openFrac = np.fromfile(os.path.join(iceberg_dir, 'openFrac.bin'), '>f8').reshape((Nr, rows, cols))
    totalBergArea = np.fromfile(os.path.join(iceberg_dir, 'totalBergArea.bin'), '>f8').reshape((Nr, rows, cols))

    iceberg_lengths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_length.bin'), '>f8').reshape((rows*cols,500))
    iceberg_widths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_widths.bin'), '>f8').reshape((rows*cols,500))
    iceberg_depths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_depths.bin'), '>f8').reshape((rows*cols,500))

    berg_components = {'bergMask': bergMask,
                        'bergMaskNums': bergMaskNums,
                        'driftMask': driftMask,
                        'barrierMask': barrierMask,
                        'meltMask': meltMask,
                        'numBergsPerCell': numBergsPerCell,
                        'openFrac': openFrac,
                        'totalBergArea': totalBergArea,
                        'iceberg_lengths': iceberg_lengths,
                        'iceberg_widths': iceberg_widths,
                        'iceberg_depths': iceberg_depths}

    return(berg_components)


def plot_iceberg_components(project_dir, model, XC, YC, Depth, berg_components):

    Cols, Rows = np.meshgrid(np.arange(XC.shape[1]), np.arange(XC.shape[0]))


    plot_width_small = 3
    plot_width_big = int(4*plot_width_small/3)
    plot_height = 4
    v_spacing = 1

    fig = plt.figure(figsize=(14, 10))

    gs = GridSpec(3*plot_height+v_spacing, 4*plot_width_small, figure=fig,
                    wspace=0.7, hspace=0.7, left=0.07, right=0.97, top=0.95, bottom=0.05)

    ax = fig.add_subplot(gs[:plot_height, :plot_width_small])
    bergMask = berg_components['bergMask']
    C = plt.pcolormesh(Cols, Rows, bergMask, cmap='Blues_r')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Iceberg Mask')
    plt.ylabel('Model Rows')
    plt.gca().set_xticklabels([])

    ax = fig.add_subplot(gs[:plot_height, plot_width_small:2*plot_width_small])
    openFrac = berg_components['openFrac']
    C = plt.pcolormesh(Cols, Rows, openFrac[0,:,:], cmap='Blues')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Open Fraction (surface)')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    ax = fig.add_subplot(gs[:plot_height, 2*plot_width_small:3*plot_width_small])
    bergMaskNums = berg_components['bergMaskNums']
    bergMaskNums = np.ma.masked_where(bergMaskNums == 0, bergMaskNums)
    C = plt.pcolormesh(Cols, Rows, bergMaskNums, cmap='viridis')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Iceberg Mask Numbers')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    ax = fig.add_subplot(gs[:plot_height, 3*plot_width_small:4*plot_width_small])
    driftMask = berg_components['driftMask']
    C = plt.pcolormesh(Cols, Rows, driftMask, cmap='Greys', vmin=0, vmax=1)
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Drift Mask')
    ax.set_yticklabels([])
    plt.gca().set_xticklabels([])

    ################################################################################

    ax = fig.add_subplot(gs[plot_height:2*plot_height,:plot_width_small])
    meltMask = berg_components['meltMask']
    C = plt.pcolormesh(Cols, Rows, meltMask, cmap='Blues')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Melt Mask')
    plt.ylabel('Model Rows')
    plt.xlabel('Model Columns')

    ax = fig.add_subplot(gs[plot_height:2*plot_height,plot_width_small:2*plot_width_small])
    numBergsPerCell = berg_components['numBergsPerCell']
    numBergsPerCell = np.ma.masked_where(numBergsPerCell == 0, numBergsPerCell)
    C = plt.pcolormesh(Cols, Rows, numBergsPerCell, cmap='turbo')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Number of Bergs Per Cell')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')

    ax = fig.add_subplot(gs[plot_height:2*plot_height,2*plot_width_small:3*plot_width_small])
    totalBergArea = berg_components['totalBergArea']
    C = plt.pcolormesh(Cols, Rows, totalBergArea[0, :, :], cmap='turbo')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Total Berg Area (surface)')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')

    ax = fig.add_subplot(gs[plot_height:2*plot_height,3*plot_width_small:4*plot_width_small])
    barrierMask = berg_components['barrierMask']
    C = plt.pcolormesh(Cols, Rows, barrierMask, cmap='Greys', vmin=0, vmax=1)
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Barrier Mask')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')

    ################################################################################

    ax = fig.add_subplot(gs[2*plot_height+v_spacing:3 * plot_height+v_spacing, :plot_width_big, ])
    iceberg_lengths = berg_components['iceberg_lengths'].ravel()
    iceberg_lengths=iceberg_lengths[iceberg_lengths > 0]
    C = plt.hist(iceberg_lengths, bins=20, color='silver', edgecolor='black')
    plt.title('Iceberg Lengths')
    plt.ylabel('Count')
    plt.xlabel('Length (m)')

    ax = fig.add_subplot(gs[2 * plot_height+v_spacing:3 * plot_height+v_spacing, plot_width_big:2*plot_width_big ])
    iceberg_widths = berg_components['iceberg_widths'].ravel()
    iceberg_widths = iceberg_widths[iceberg_widths > 0]
    C = plt.hist(iceberg_widths, bins=20, color='silver', edgecolor='black')
    plt.title('Iceberg Widths')
    plt.xlabel('Width (m)')

    ax = fig.add_subplot(gs[2 * plot_height+v_spacing:3 * plot_height+v_spacing, 2*plot_width_big:3*plot_width_big, ])
    iceberg_depths = berg_components['iceberg_depths'].ravel()
    iceberg_depths = iceberg_depths[iceberg_depths > 0]
    C = plt.hist(iceberg_depths, bins=20, color='silver', edgecolor='black')
    plt.title('Iceberg Depths')
    plt.xlabel('Depth (m)')

    output_file = os.path.join(project_dir, 'Figures', 'Melange',
                               model+'_iceberg_components.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/MITgcm/' \
              'configurations/downscale_darwin/'

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/'

model = 'L2_Upernavik'

if model == 'L2_Upernavik':
    iceberg_dir = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'input', 'iceberg')
    XC, YC, Depth = read_grid_components(config_dir, 'L2_Upernavik')
    Nr = 61
if model == 'verify':
    iceberg_dir = os.path.join(config_dir,'L2','L2_Upernavik','input_verify')
    XC = np.zeros((10,30))
    YC = np.zeros((10,30))
    Depth = np.fromfile(os.path.join(config_dir,'L2','L2_Upernavik','input_verify','bathymetry.bin'), '>f8').reshape((10,30))
    Nr = 12

rows, cols = XC.shape

berg_components = read_iceberg_fields(iceberg_dir, rows, cols, Nr)

plot_iceberg_components(project_dir, model, XC, YC, Depth, berg_components)