


import os
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from matplotlib.gridspec import GridSpec
from pyproj import Transformer

# ignore RuntimeWarnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1, run_test = True):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_model_domain(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                    model_name+'_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    ds.close()

    return(XC, YC, Depth)

def read_west_greenland_chl_data(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Observations','Chlorophyll',
                                  'West_Greenland_Chl_Climatology.nc'))
    chl = ds.variables['chlor_a'][:]
    depth = ds.variables['depth'][:]
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()

    x = x / 1000
    y = y / 1000

    return(chl, depth, x, y)

def read_upernavik_chl_data(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Ocean','Chlorophyll',
                                  'L2_Upernavik_Chl_Climatology.nc'))
    chl = ds.variables['chlor_a'][:]
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()

    x = x/1000
    y = y/1000

    return(chl, x, y)

def read_west_greenland_MODIS_data(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Imagery', 'MODIS',
                                  'West_Greenland_MOIDS_Imagery.nc'))

    R = ds.variables['band_4'][:]
    G = ds.variables['band_3'][:]
    B = ds.variables['band_1'][:]
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]

    R = R.reshape((rows, cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B, G, R], axis=2)
    image = np.flipud(image)

    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()

    minval = np.percentile(image[~np.isnan(image)], 2)
    maxval = np.percentile(image[~np.isnan(image)], 95)
    image = np.clip(image, minval, maxval)
    # print(minval, maxval)

    image = ((image - minval) / (maxval - minval)) * 255
    image = image.astype(int)

    x = x / 1000
    y = y / 1000

    return(image, x, y)

def read_upernavik_MODIS_data(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Imagery', 'MODIS',
                                  'L2_Upernavik_MODIS_20220720_3413.nc'))

    R = ds.variables['band_4'][:]
    G = ds.variables['band_3'][:]
    B = ds.variables['band_1'][:]
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]

    R = R.reshape((rows, cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B, G, R], axis=2)
    image = np.flipud(image)

    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()

    minval = np.percentile(image[~np.isnan(image)], 2)
    maxval = np.percentile(image[~np.isnan(image)], 95)
    image = np.clip(image, minval, maxval)

    # plt.imshow(image)
    # plt.show()

    # image = ((image - minval) / (maxval - minval)) * 255
    # image = image.astype(int)

    # plt.imshow(image)
    # plt.show()

    x = x / 1000
    y = y / 1000

    return(image, x, y)

def read_upernavik_landsat_data(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Imagery', 'Landsat',
                                  'Upernavik Fjord Landsat.nc'))

    R = ds.variables['band_4'][:]
    G = ds.variables['band_3'][:]
    B = ds.variables['band_2'][:]
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]

    R = R.reshape((rows, cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B, G, R], axis=2)
    image = np.flipud(image)

    x = ds.variables['X'][:,:]
    y = ds.variables['Y'][:,:]
    ds.close()

    minval = np.percentile(image, 2)
    maxval = np.percentile(image, 95)
    image = np.clip(image, minval, maxval)

    image = ((image - minval) / (maxval - minval)) * 255
    image = image.astype(int)

    x = x / 1000
    y = y / 1000

    return(image)

def get_domain_extents(config_dir, domain_name):

    if domain_name=='West Greenland':
        min_x = -943593
        max_x = -111839
        min_y = -3363131
        max_y = -1095405

    elif domain_name=='Upernavik':
        XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
        points = np.column_stack([XC.ravel(), YC.ravel()])
        reprojected_points = reproject_polygon(points, outputCRS=3413, inputCRS=4326)
        X = reprojected_points[:,0]
        Y = reprojected_points[:,1]
        min_x = np.min(X)
        max_x = np.max(X)
        min_y = np.min(Y)
        max_y = np.max(Y)

    elif domain_name=='Upernavik Fjord':

        min_x = -389568
        max_x = -284882
        min_y = -1864778
        max_y = -1803830

    # extents = (min_x, max_x, min_y, max_y)

    extents = (min_x/1000, max_x/1000, min_y/1000, max_y/1000)

    return(extents)

def add_lonlat_annotations_to_ax(ax, X, Y, Lon, Lat, side, locations, region):

    if side=='left':
        X_subset = X[:,0]
        Y_subset = Y[:,0]
        Lon_subset = Lon[:,0]
        Lat_subset = Lat[:,0]
    elif side=='bottom':
        X_subset = X[0, :]
        Y_subset = Y[0, :]
        Lon_subset = Lon[0, :]
        Lat_subset = Lat[0, :]
    elif side == 'top':
        X_subset = X[-1, :]
        Y_subset = Y[-1, :]
        Lon_subset = Lon[-1, :]
        Lat_subset = Lat[-1, :]
    elif side == 'right':
        X_subset = X[:,-1]
        Y_subset = Y[:,-1]
        Lon_subset = Lon[:,-1]
        Lat_subset = Lat[:,-1]

    for location in locations:

        if side in ['left', 'right']:
            index = (np.abs(Lat_subset - location)).argmin()
            if region=='WG':
                shift = -5 if side=='left' else 5
                txt =  f'{int(round(Lat_subset[index]))}°N'
            else:
                shift = -1 if side == 'left' else 1
                txt = f'{Lat_subset[index]:.1f}°N'

            ax.text(X_subset[index]+shift, Y_subset[index], txt,
                    verticalalignment='center',
                    horizontalalignment='right' if side=='left' else 'left',
                    fontsize=10, color='k')
        elif side in ['top', 'bottom']:
            index = (np.abs(-1 * Lon_subset - location)).argmin()
            if region == 'WG':
                shift = -5 if side == 'bottom' else 5
                txt = f'{int(-1*round(Lon_subset[index]))}°W'
            else:
                shift = -1 if side == 'bottom' else 1
                txt = f'{int(-1*round(Lon_subset[index]))}°W'

            ax.text(X_subset[index], Y_subset[index]+shift, txt,
                    verticalalignment='top' if side=='bottom' else 'bottom',
                    horizontalalignment='center',
                    fontsize=10, color='k')



def create_map_with_chl_figure(config_dir, project_dir):

    XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
    points = np.column_stack([XC.ravel(), YC.ravel()])
    reprojected_points = reproject_polygon(points, outputCRS=3413, inputCRS=4326)
    X = reprojected_points[:,0].reshape(XC.shape)
    Y = reprojected_points[:,1].reshape(YC.shape)

    print(np.min(X), np.max(X))
    print(np.min(Y), np.max(Y))

    X = X / 1000
    Y = Y / 1000
    bbox = np.array([[np.min(X), np.min(Y)],
                     [np.max(X), np.min(Y)],
                     [np.max(X), np.max(Y)],
                     [np.min(X), np.max(Y)],
                     [np.min(X), np.min(Y)]])

    extents = get_domain_extents(config_dir, 'Upernavik Fjord')
    fjord_bbox = np.array([[extents[0], extents[2]],
                           [extents[1], extents[2]],
                           [extents[1], extents[3]],
                           [extents[0], extents[3]],
                           [extents[0], extents[2]]])


    fig = plt.figure(figsize=(10, 8))

    left_width = 10
    right_width = 8
    top_left_height = 6
    bottom_left_height = 4
    spacing = 1
    right_height = top_left_height + bottom_left_height + spacing


    gs = GridSpec(top_left_height + bottom_left_height + spacing,
                  left_width + right_width + spacing, figure=fig,
                  left = 0.1, right = 0.95, top = 0.95, bottom = 0.1)

    #############################################################################
    ax1 = fig.add_subplot(gs[0:top_left_height, 0:left_width])

    extents = get_domain_extents(config_dir, 'Upernavik')

    chl_upernavik, x_chl_upernavik, y_chl_upernavik = read_upernavik_chl_data(project_dir)
    modis_image, modis_x, modis_y = read_upernavik_MODIS_data(project_dir)
    modis_extent = (np.min(modis_x), np.max(modis_x),
                    np.min(modis_y), np.max(modis_y))
    ax1.imshow(modis_image, extent=modis_extent)
    chl_upernavik = np.ma.masked_where(Depth<0.1, chl_upernavik)
    C = plt.pcolormesh(x_chl_upernavik, y_chl_upernavik, chl_upernavik,
                       cmap='turbo', vmin=0, vmax=5)
    ax1.contour(X, Y, Depth, levels=[0.1], colors='k', linewidths=0.5)
    ax1.plot(fjord_bbox[:, 0], fjord_bbox[:, 1], 'k-', linewidth=2)

    ax1.contour(X, Y, Depth, levels=[0.1], colors='k', linewidths=0.5)

    ax1.contour(X, Y, -1*XC, levels=[55, 60], colors='k', linewidths=0.5, alpha=0.5)
    ax1.contour(X, Y, YC, levels=[72.5, 73, 73.5, 74], colors='k', linewidths=0.5, alpha=0.5)
    ax1.set_xlim(extents[0], extents[1])
    ax1.set_ylim(extents[2], extents[3])

    add_lonlat_annotations_to_ax(ax1, X, Y, XC, YC, side='left', locations=[72.5, 73, 73.5], region='L2')
    add_lonlat_annotations_to_ax(ax1, X, Y, XC, YC, side='bottom', locations=[55, 60], region='L2')
    add_lonlat_annotations_to_ax(ax1, X, Y, XC, YC, side='top', locations=[55, 60], region='L2')
    add_lonlat_annotations_to_ax(ax1, X, Y, XC, YC, side='right', locations=[73, 73.5, 74], region='L2')

    ax1.text(0.98, 0.98, 'a)', transform=ax1.transAxes,
             verticalalignment='top', horizontalalignment='right',
             fontsize=12, color='k',
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))

    ax1.set_xticks([])
    ax1.set_yticks([])


    #############################################################################

    ax2 = fig.add_subplot(gs[:, left_width + spacing:])
    extents = get_domain_extents(config_dir, 'West Greenland')

    modis_image, modis_x, modis_y = read_west_greenland_MODIS_data(project_dir)
    modis_extent = (np.min(modis_x), np.max(modis_x),
                    np.min(modis_y), np.max(modis_y))
    ax2.imshow(modis_image, extent=modis_extent)

    chl, depth, x_chl, y_chl = read_west_greenland_chl_data(project_dir)
    chl = np.ma.masked_where(depth<0.5, chl)
    C = plt.pcolormesh(x_chl, y_chl, chl, cmap='turbo', vmin=0, vmax=5)

    # y_chl_flipped = np.flipud(y_chl)
    X_WG, Y_WG = np.meshgrid(x_chl, y_chl)
    points_WG = np.column_stack([X_WG.ravel()*1000, Y_WG.ravel()*1000])
    lonlat_WG = reproject_polygon(points_WG, outputCRS=4326, inputCRS=3413)
    Lon_WG = lonlat_WG[:,0].reshape(X_WG.shape)
    Lat_WG = lonlat_WG[:,1].reshape(X_WG.shape)

    ax2.plot(bbox[:,0], bbox[:,1], 'k-', linewidth=2)

    ax2.contour(X_WG, Y_WG, depth, levels=[0.1], colors='k', linewidths=0.5)

    ax2.contour(X_WG, Y_WG, -1*Lon_WG, levels=[50, 60, 70, 80], colors='k', linewidths=0.5, alpha=0.5)
    ax2.contour(X_WG, Y_WG, Lat_WG, levels=[60,65,70,75], colors='k',linewidths=0.5, alpha=0.5)

    ax2.set_xlim(extents[0], extents[1])
    ax2.set_ylim(extents[2], extents[3])

    add_lonlat_annotations_to_ax(ax2, X_WG, Y_WG, Lon_WG, Lat_WG, side='left', locations=[60,65,70,75], region='WG')
    add_lonlat_annotations_to_ax(ax2, X_WG, Y_WG, Lon_WG, Lat_WG, side='bottom', locations=[50, 60], region='WG')
    add_lonlat_annotations_to_ax(ax2, X_WG, Y_WG, Lon_WG, Lat_WG, side='top', locations=[60, 70, 80], region='WG')
    add_lonlat_annotations_to_ax(ax2, X_WG, Y_WG, Lon_WG, Lat_WG, side='right', locations=[60,65,70,75], region='WG')

    # add a b) annotation letter in the upper right corner with transparent white background
    ax2.text(0.98, 0.99, 'b)', transform=ax2.transAxes,
             verticalalignment='top', horizontalalignment='right',
             fontsize=12, color='k',
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))


    ax2.set_xticks([])
    ax2.set_yticks([])

    #############################################################################
    ax3 = fig.add_subplot(gs[top_left_height + spacing:,
                          0:left_width])
    landsat_image = read_upernavik_landsat_data(project_dir)
    ax3.imshow(landsat_image)#,
               #extent=(extents[0], extents[1], extents[2], extents[3]))

    ax3.text(0.98, 0.98, 'c)', transform=ax3.transAxes,
             verticalalignment='top', horizontalalignment='right',
             fontsize=12, color='k',
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))


    plt.savefig(os.path.join(project_dir, 'Figures', 'L2_Upernavik_domains_with_chl.png'),
                dpi=300)
    plt.close(fig)


config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/' \
             'Downscale_Darwin/MITgcm/configurations/downscale_darwin/'

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

create_map_with_chl_figure(config_dir, project_dir)