
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from pyproj import Transformer
from matplotlib.patches import Rectangle


def reproject_polygon(polygon_array, inputCRS, outputCRS, x_column=0, y_column=1, run_test=True):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(polygon_array)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def read_model_domain(config_dir, model_name):
    domain_file = os.path.join(config_dir, 'nc_grids', model_name+'_grid.nc')
    ds = nc4.Dataset(domain_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    poitns = np.column_stack([XC.ravel(), YC.ravel()])
    reprojected_points = reproject_polygon(poitns, outputCRS=3413, inputCRS=4326)
    X = reprojected_points[:, 0].reshape(np.shape(XC))
    Y = reprojected_points[:, 1].reshape(np.shape(YC))

    return(X, Y, Depth)

def get_glacier_domain_3413(glacier_name):

    if glacier_name == 'Upernavik SS':
        min_x = -308352
        max_x = -298154
        min_y = -1856426
        max_y = -1848510
    elif glacier_name == 'Upernavik S':
        min_x = -308546
        max_x = -299999
        min_y = -1851618
        max_y = -1842707
    elif glacier_name == 'Upernavik C':
        min_x = -308570
        max_x = -299052
        min_y = -1840012
        max_y = -1832096
    elif glacier_name == 'Upernavik N':
        min_x = -304758
        max_x = -298469
        min_y = -1832387
        max_y = -1824836
    elif glacier_name == 'Upernavik NW':
        min_x = -308910
        max_x = -303593
        min_y = -1830663
        max_y = -1823670
    elif glacier_name == 'Akullikassaap':
        min_x = -318477
        max_x = -312392
        min_y = -1827410
        max_y = -1821233
    elif glacier_name == 'Akullikassaap E':
        min_x = -313815
        max_x = -313961
        min_y = -1832509
        max_y = -1821534
    elif glacier_name=='Nunatakassaap':
        min_x = -325932
        max_x = -318599
        min_y = -1808082
        max_y = -1797252
    elif glacier_name == 'Kakivfaat':
        min_x = -330424
        max_x = -319643
        min_y = -1780279
        max_y = -1766730
    elif glacier_name == 'Kakivfaat S':
        min_x = -319521
        max_x = -314252
        min_y = -1788996
        max_y = -1782950
    elif glacier_name == 'Qeqertarsuup':
        min_x = -331953
        max_x = -325859
        min_y = -1764715
        max_y = -1755974
    elif glacier_name == 'Ussing Braeer':
        min_x = -329477
        max_x = -320298
        min_y = -1734290
        max_y = -1725840
    elif glacier_name == 'Ussing Braeer N':
        min_x = -330399
        max_x = -322386
        min_y = -1724772
        max_y = -1715885
    elif glacier_name == 'Cornell':
        min_x = -334187
        max_x = -326126
        min_y = -1691069
        max_y = -1682182
    elif glacier_name == 'Cornell N':
        min_x = -333750
        max_x = -325834
        min_y = -1682765
        max_y = -1676840
    else:
        raise ValueError('Glacier name not recognized for domain extents')

    # if x extent smaller than y extent, extend left
    x_extent = max_x - min_x
    y_extent = max_y - min_y
    if x_extent < y_extent:
        difference = y_extent - x_extent
        min_x = min_x - difference / 2
        max_x = max_x + difference / 2
    if y_extent < x_extent:
        difference = x_extent - y_extent
        min_y = min_y - difference / 2
        max_y = max_y + difference / 2

    return (min_x, max_x, min_y, max_y)

def read_bedmachine_subset_around_glacier(min_x, max_x, min_y, max_y):
    bm5_file = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Bathymetry/BedMachineGreenland-v5.nc'

    ds = nc4.Dataset(bm5_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    bed = ds.variables['bed'][:, :]
    surface = ds.variables['surface'][:, :]
    mask = ds.variables['mask'][:, :]
    ds.close()

    x_indices = np.where((x >= min_x) & (x <= max_x))[0]
    y_indices = np.where((y >= min_y) & (y <= max_y))[0]
    bed = bed[y_indices, :][:, x_indices]
    surface = surface[y_indices, :][:, x_indices]
    mask = mask[y_indices, :][:, x_indices]

    x_subset = x[x_indices]
    y_subset = y[y_indices]

    return (x_subset, y_subset, bed, surface, mask)

def get_glacier_iceplume_row_col_extent(glacier_name):
    if glacier_name == 'Upernavik SS':
        min_row = 22
        max_row = 28
        min_col = 415
        max_col = 423
    elif glacier_name == 'Upernavik S':
        min_row = 33
        max_row = 42
        min_col = 415
        max_col = 420
    elif glacier_name == 'Upernavik C':
        min_row = 55
        max_row = 61
        min_col = 418
        max_col = 422
    elif glacier_name == 'Upernavik N':
        min_row = 69
        max_row = 77
        min_col = 420
        max_col = 426
    elif glacier_name == 'Upernavik NW':
        min_row = 74
        max_row = 79
        min_col = 411
        max_col = 417
    elif glacier_name == 'Akullikassaap':
        min_row = 79
        max_row = 82
        min_col = 392
        max_col = 396
    elif glacier_name == 'Akullikassaap E':
        min_row = 67
        max_row = 72
        min_col = 400
        max_col = 404
    elif glacier_name=='Nunatakassaap':
        min_row = 120
        max_row = 132
        min_col = 378
        max_col = 385
    elif glacier_name == 'Kakivfaat':
        min_row = 175
        max_row = 195
        min_col = 370
        max_col = 380
    elif glacier_name == 'Kakivfaat S':
        min_row = 157
        max_row = 162
        min_col = 393
        max_col = 395
    elif glacier_name == 'Qeqertarsuup':
        min_row = 206
        max_row = 215
        min_col = 366
        max_col = 370
    elif glacier_name == 'Ussing Braeer':
        min_row = 268
        max_row = 273
        min_col = 373
        max_col = 381
    elif glacier_name == 'Ussing Braeer N':
        min_row = 285
        max_row = 295
        min_col = 368
        max_col = 378
    elif glacier_name == 'Cornell':
        min_row = 352
        max_row = 367
        min_col = 362
        max_col = 370
    elif glacier_name == 'Cornell N':
        min_row = 367
        max_row = 372
        min_col = 365
        max_col = 369
    else:
        raise ValueError('Glacier name not recognized for iceplume row/col extents')

    return (min_row, max_row, min_col, max_col)

def read_iceplume_mask(config_dir, min_row, max_row, min_col, max_col):

    mask_file = os.path.join(config_dir,'L2','L2_Upernavik','input', 'iceplume',
                             'L2_iceplume_mask.bin')
    n_rows = 375
    n_cols = 450

    iceplume_mask = np.fromfile(mask_file, dtype='>f4').reshape((n_rows, n_cols))
    #iceplume_mask_subset = iceplume_mask[min_row:max_row+1, min_col:max_col+1]

    return iceplume_mask

def read_iceplume_Qsg_timeseries(config_dir, row, col):
    mask_file = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'input', 'iceplume','L2_Qsg_2016')
    n_rows = 375
    n_cols = 450
    Qsg = np.fromfile(mask_file, dtype='>f4').reshape((366, n_rows, n_cols))
    Qsg_timeseries = Qsg[:, row, col]
    return(Qsg_timeseries)

def plot_glacier_implementation(config_dir, glacier_name,
                                min_x, max_x, min_y, max_y,
                                bm_x, bm_y, bed, surface, mask,
                                model_X, model_Y, model_Depth,
                                min_row, max_row, min_col, max_col, iceplume_mask):

    figure_file = os.path.join(config_dir,'L2','L2_Upernavik','plots','init_files', 'iceplume',
                               glacier_name.replace(' ', '_') +'_iceplume_implementation.png')

    plot_width = 12
    plot_height = 12
    timeseries_height = 4
    legend_width = 3

    fig = plt.figure(figsize=(12,8))

    gs = GridSpec(plot_height+timeseries_height+2,2*plot_width+legend_width, figure=fig,
                  left =0.09, right=0.95, wspace=0.3, bottom =0.07, top=0.95)

    ax1 = fig.add_subplot(gs[:plot_height,:plot_width])
    ax1.pcolormesh(bm_x, bm_y, bed, cmap='Blues_r', vmin=-300, vmax=0)

    land_mask = surface>0
    land_mask = np.ma.masked_where(~land_mask, land_mask)
    ax1.pcolormesh(bm_x, bm_y, land_mask, cmap='Greys')
    ax1.pcolormesh(bm_x, bm_y, land_mask, cmap='BrBG', alpha=0.5)

    ice_mask = mask==2
    ice_mask = np.ma.masked_where(~ice_mask, ice_mask)
    ax1.pcolormesh(bm_x, bm_y, ice_mask, cmap='Greys')

    ax1.set_title(glacier_name + ' BedMachine Domain')
    ax1.set_xlim(min_x, max_x)
    ax1.set_ylim(min_y, max_y)
    ax1.set_xlabel('Distance East (m, 3413)')
    ax1.set_ylabel('Distance North (m, 3413)')

    ################################################################################

    ax2 = fig.add_subplot(gs[:plot_height,plot_width:2*plot_width])
    ax2.pcolormesh(model_X, model_Y, -1*model_Depth, cmap='Blues_r', vmin=-300, vmax=0)

    model_land_mask = model_Depth<=0
    model_land_mask = np.ma.masked_where(~model_land_mask, model_land_mask)
    ax2.pcolormesh(model_X, model_Y, model_land_mask, cmap='BrBG',alpha=0.5)

    ice_rows, ice_cols = np.where(iceplume_mask!=0)
    for i in range(len(ice_rows)):
        row = ice_rows[i]
        col = ice_cols[i]
        x0 = model_X[row, col] - 250
        y0 = model_Y[row, col] - 250
        width = 500
        height = 500
        rect_outline = Rectangle((x0, y0), width, height,
                         linewidth=1, edgecolor='k', facecolor='none', zorder=99)
        ax2.add_patch(rect_outline)

        if row>=min_row and row<=max_row and col>=min_col and col<=max_col:
            mask_val = iceplume_mask[row, col]

            # print(x0, y0)
            # print(min_row + row, min_col + col)

            if mask_val == 1:
                color = '#ff7f0e'
            elif mask_val == -1:
                color = '#d62728'
            elif mask_val == 6:
                color = '#2ca02c'
            elif mask_val == -6:
                color = '#9467bd'
            else:
                raise ValueError('Unrecognized iceplume mask value')
            rect = Rectangle((x0, y0), width, height,
                             linewidth=1, edgecolor='k', facecolor=color, zorder=99)
            ax2.add_patch(rect)

    # c = ax2.pcolormesh(np.arange(min_col, max_col+1), np.arange(min_row, max_row+1),
    #                       iceplume_mask)

    ax2.set_title(glacier_name + ' Model Domain with Iceplume Mask')
    ax2.set_xlim(min_x, max_x)
    ax2.set_ylim(min_y, max_y)
    ax2.set_yticklabels([])
    ax2.set_xlabel('Distance East (m, 3413)')

    ################################################################################

    ax3 = fig.add_subplot(gs[3:-3-timeseries_height-2,2*plot_width:2*plot_width+legend_width])

    colors = ['#ff7f0e', '#d62728', '#2ca02c', '#9467bd']
    labels = ['Vertical\nMelt Cell (1)', 'Horizontal \nMelt Cell (-1)',
                'Vertical\nPlume Cell (6)', 'Horizontal\nPlume Cell (-6)']

    for i in range(len(colors)):
        rect = Rectangle((0, i+0.1), 1, 0.8, facecolor=colors[i], edgecolor='k', linewidth=1)
        ax3.add_patch(rect)
        ax3.text(1.2, i + 0.5, labels[i], va='center', fontsize=10)
    ax3.set_ylim(-1, 5)
    ax3.set_xlim(0, 5)
    ax3.axis('off')

    ################################################################################

    ax4 = fig.add_subplot(gs[-timeseries_height:, :2*plot_width])

    for i in range(len(ice_rows)):
        row = ice_rows[i]
        col = ice_cols[i]
        if row>=min_row and row<=max_row and col>=min_col and col<=max_col:
            mask_val = iceplume_mask[row, col]
            if mask_val == 6 or mask_val == -6:
                Qsg_timeseries = read_iceplume_Qsg_timeseries(config_dir, row, col)

                dec_yrs = np.arange(2016, 2017, 1/366)
                ax4.plot(dec_yrs, Qsg_timeseries, 'k-', label=f'Cell at row {row}, col {col}')

    ax4.set_ylabel('Subglacial Discharge\n(m$^3$/s)')
    ax4.set_title(glacier_name + ' Iceplume Subglacial Discharge Timeseries (2016)')
    ax4.legend()

    plt.savefig(figure_file, dpi=300)
    plt.close()

config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/' \
              'Greenland Model Analysis/Fjord/Upernavik'

glacier_names = ['Cornell N', 'Cornell','Ussing Braeer N','Ussing Braeer',
                 'Qeqertarsuup','Kakivfaat','Kakivfaat S','Nunatakassaap',
                 'Upernavik NW', 'Akullikassaap', 'Akullikassaap E',
                 'Upernavik N','Upernavik C','Upernavik S','Upernavik SS']

for glacier_name in glacier_names:
    print('Processing glacier:', glacier_name)

    model_X, model_Y, model_Depth = read_model_domain(config_dir, 'L2_Upernavik')

    min_x, max_x, min_y, max_y = get_glacier_domain_3413(glacier_name)

    bm_x, bm_y, bed, surface, mask = read_bedmachine_subset_around_glacier(min_x, max_x, min_y, max_y)

    min_row, max_row, min_col, max_col = get_glacier_iceplume_row_col_extent(glacier_name)

    iceplume_mask = read_iceplume_mask(config_dir, min_row, max_row, min_col, max_col)

    plot_glacier_implementation(config_dir, glacier_name,
                                min_x, max_x, min_y, max_y,
                                bm_x, bm_y, bed, surface, mask,
                                model_X, model_Y, model_Depth,
                                min_row, max_row, min_col, max_col, iceplume_mask)





