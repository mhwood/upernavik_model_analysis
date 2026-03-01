
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

    precision = '>f4'

    bergMask = np.fromfile(os.path.join(iceberg_dir, 'bergMask.bin'), precision).reshape((rows, cols))
    #bergConc = np.fromfile(os.path.join(iceberg_dir, 'bergConc.bin'), precision).reshape((rows, cols))
    bergMaskNums = np.fromfile(os.path.join(iceberg_dir, 'bergMaskNums.bin'), precision).reshape((rows, cols))
    driftMask = np.fromfile(os.path.join(iceberg_dir, 'driftMask.bin'), precision).reshape((rows, cols))
    barrierMask = np.fromfile(os.path.join(iceberg_dir, 'barrierMask.bin'), precision).reshape((rows, cols))
    meltMask = np.fromfile(os.path.join(iceberg_dir, 'meltMask.bin'), precision).reshape((rows, cols))
    numBergsPerCell = np.fromfile(os.path.join(iceberg_dir, 'numBergsPerCell.bin'), precision).reshape((rows, cols))
    openFrac = np.fromfile(os.path.join(iceberg_dir, 'openFrac.bin'), precision).reshape((Nr, rows, cols))
    totalBergArea = np.fromfile(os.path.join(iceberg_dir, 'totalBergArea.bin'), precision).reshape((Nr, rows, cols))

    iceberg_lengths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_length.bin'), precision).reshape((rows*cols,500))
    iceberg_widths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_widths.bin'), precision).reshape((rows*cols,500))
    iceberg_depths = np.fromfile(os.path.join(iceberg_dir, 'icebergs_depths.bin'), precision).reshape((rows*cols,500))

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

    fig = plt.figure(figsize=(12, 10))

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
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[:plot_height, plot_width_small:2*plot_width_small])
    openFrac = berg_components['openFrac']
    C = plt.pcolormesh(Cols, Rows, openFrac[0,:,:], cmap='Blues')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Open Fraction (surface)')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[:plot_height, 2*plot_width_small:3*plot_width_small])
    bergMaskNums = berg_components['bergMaskNums']
    bergMaskNums = np.ma.masked_where(bergMaskNums == 0, bergMaskNums)
    C = plt.pcolormesh(Cols, Rows, bergMaskNums, cmap='viridis')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Iceberg Mask Numbers')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[:plot_height, 3*plot_width_small:4*plot_width_small])
    driftMask = berg_components['driftMask']
    C = plt.pcolormesh(Cols, Rows, driftMask, cmap='Greys', vmin=0, vmax=1)
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Drift Mask')
    ax.set_yticklabels([])
    plt.gca().set_xticklabels([])
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ################################################################################

    ax = fig.add_subplot(gs[plot_height:2*plot_height,:plot_width_small])
    meltMask = berg_components['meltMask']
    C = plt.pcolormesh(Cols, Rows, meltMask, cmap='Blues')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Melt Mask')
    plt.ylabel('Model Rows')
    plt.xlabel('Model Columns')
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[plot_height:2*plot_height,plot_width_small:2*plot_width_small])
    numBergsPerCell = berg_components['numBergsPerCell']
    numBergsPerCell = np.ma.masked_where(numBergsPerCell == 0, numBergsPerCell)
    C = plt.pcolormesh(Cols, Rows, numBergsPerCell, cmap='turbo')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Number of Bergs Per Cell')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[plot_height:2*plot_height,2*plot_width_small:3*plot_width_small])
    totalBergArea = berg_components['totalBergArea']
    C = plt.pcolormesh(Cols, Rows, totalBergArea[0, :, :], cmap='turbo')
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Total Berg Area (surface)')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')
    plt.gca().set_xlim([250, np.shape(XC)[1]])

    ax = fig.add_subplot(gs[plot_height:2*plot_height,3*plot_width_small:4*plot_width_small])
    barrierMask = berg_components['barrierMask']
    C = plt.pcolormesh(Cols, Rows, barrierMask, cmap='Greys', vmin=0, vmax=1)
    plt.contour(Cols, Rows, Depth, levels=[0.1], colors='k', linewidths=0.5)
    plt.colorbar(C)
    plt.title('Barrier Mask')
    ax.set_yticklabels([])
    plt.xlabel('Model Columns')
    plt.gca().set_xlim([250, np.shape(XC)[1]])

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

    output_file = os.path.join(project_dir, 'Figures','Ocean', 'Melange',
                               model+'_iceberg_components.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

def print_berg_stats(berg_components):

    numBergsPerCell = berg_components['numBergsPerCell']
    bergMaskNums = berg_components['bergMaskNums']

    upernavik_rows = np.arange(100)
    kakivfaat_rows = np.arange(150, 250)
    ussing_braeer_rows = np.arange(250, 350)

    upernavik_berg_nums = np.arange(np.min(bergMaskNums[upernavik_rows,:][bergMaskNums[upernavik_rows,:]>0]),
                                    np.max(bergMaskNums[upernavik_rows,:])+1).astype(int)
    kakivfaat_berg_nums = np.arange(np.min(bergMaskNums[kakivfaat_rows,:][bergMaskNums[kakivfaat_rows,:]>0]),
                                    np.max(bergMaskNums[kakivfaat_rows,:])+1).astype(int)
    ussing_braeer_berg_nums = np.arange(np.min(bergMaskNums[ussing_braeer_rows,:][bergMaskNums[ussing_braeer_rows,:]>0]),
                                    np.max(bergMaskNums[ussing_braeer_rows,:])+1).astype(int)

    iceberg_widths = berg_components['iceberg_widths']#.ravel()
    iceberg_depths = berg_components['iceberg_depths']#.ravel()
    iceberg_lengths = berg_components['iceberg_lengths']#.ravel()

    upernavik_widths = iceberg_widths[upernavik_berg_nums.ravel(), :]
    upernavik_depths = iceberg_depths[upernavik_berg_nums.ravel(), :]
    upernavik_lengths = iceberg_lengths[upernavik_berg_nums.ravel(), :]
    median_upernavik_width = np.median(upernavik_widths[upernavik_widths > 0])
    median_upernavik_depth = np.median(upernavik_depths[upernavik_depths > 0])
    median_upernavik_length = np.median(upernavik_lengths[upernavik_lengths > 0])

    kakivfaat_widths = iceberg_widths[kakivfaat_berg_nums.ravel(), :]
    kakivfaat_depths = iceberg_depths[kakivfaat_berg_nums.ravel(), :]
    kakivfaat_lengths = iceberg_lengths[kakivfaat_berg_nums.ravel(), :]
    median_kakivfaat_width = np.median(kakivfaat_widths[kakivfaat_widths > 0])
    median_kakivfaat_depth = np.median(kakivfaat_depths[kakivfaat_depths > 0])
    median_kakivfaat_length = np.median(kakivfaat_lengths[kakivfaat_lengths > 0])

    ussing_braeer_widths = iceberg_widths[ussing_braeer_berg_nums.ravel(), :]
    ussing_braeer_depths = iceberg_depths[ussing_braeer_berg_nums.ravel(), :]
    ussing_braeer_lengths = iceberg_lengths[ussing_braeer_berg_nums.ravel(), :]
    median_ussing_braeer_width = np.median(ussing_braeer_widths[ussing_braeer_widths > 0])
    median_ussing_braeer_depth = np.median(ussing_braeer_depths[ussing_braeer_depths > 0])
    median_ussing_braeer_length = np.median(ussing_braeer_lengths[ussing_braeer_lengths > 0])

    print('Upernavik Stats:')
    print('    Total bergs: ', np.sum(numBergsPerCell[upernavik_rows,:]))
    print('     median width: ', median_upernavik_width)
    print('     median depth: ', median_upernavik_depth)
    print('     median length: ', median_upernavik_length)

    print('Kakivfaat Stats:')
    print('    Total bergs: ', np.sum(numBergsPerCell[kakivfaat_rows,:]))
    print('     median width: ', median_kakivfaat_width)
    print('     median depth: ', median_kakivfaat_depth)
    print('     median length: ', median_kakivfaat_length)

    print('Ussing Braeer Stats:')
    print('    Total bergs: ', np.sum(numBergsPerCell[ussing_braeer_rows,:]))
    print('     median width: ', median_ussing_braeer_width)
    print('     median depth: ', median_ussing_braeer_depth)
    print('     median length: ', median_ussing_braeer_length)

    print('Total Bergs', np.sum(numBergsPerCell))


# config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/MITgcm/' \
#               'configurations/downscale_darwin/'
#
# project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/'

config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin/'
# config_d'/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/Downscaled_Darwin/darwin3/configurations/downscaled_ecco_v5_darwin/'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

model = 'L2_Upernavik'


iceberg_dir = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'input', 'iceberg')
XC, YC, Depth = read_grid_components(config_dir, 'L2_Upernavik')
Nr = 61

rows, cols = XC.shape

berg_components = read_iceberg_fields(iceberg_dir, rows, cols, Nr)

print_berg_stats(berg_components)

plot_iceberg_components(project_dir, model, XC, YC, Depth, berg_components)