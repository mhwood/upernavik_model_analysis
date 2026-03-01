
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
from osgeo import gdal
import shapefile
from pyproj import Transformer

def read_fjord_geometry(fjord_file):

    ds = nc4.Dataset(fjord_file)
    x = ds.variables['local_x'][:]/1e3
    y = ds.variables['local_y'][:]/1e3
    Lon = ds.variables['Lon'][:, :]
    Lat = ds.variables['Lat'][:, :]
    Depth = ds.variables['Bed'][:, :]
    ds.close()

    return(Lon, Lat)

def read_landsat_image(file_path):

    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]

    transform = ds.GetGeoTransform()
    extents = [transform[0], transform[0] + transform[1] * cols,
               transform[3] + transform[5] * rows, transform[3]]

    ds = None

    return(R, extents)



    # x_resolution = transform[1]
    # y_resolution = transform[5]

def get_melange_region_bounds(location):

    if location=='Upernavik_Fjord':
        min_x = -342356
        min_y = -1859037
        # max_x = -287777
        # max_y = -1808695
        distance = 48000
        max_x = min_x + distance
        max_y = min_y + distance
    if location=='Ussing_Braeer':
        min_x = -342541
        min_y = -1737970
        distance = 30000
        # max_x = -311555
        # max_y = -1703856
        max_x = min_x + distance
        max_y = min_y + distance
    if location=='Kakivfaat':
        min_x = -339699
        min_y = -1781649
        # max_x = -312216
        # max_y = -1748370
        distance = 25000
        max_x = min_x + distance
        max_y = min_y + distance

    # make sure the subset is square
    if max_x-min_x > max_y-min_y:
        print('Expanded by ' + str(max_x-min_x) + ' m in y direction to make square')
        max_y = min_y + (max_x - min_x)
    elif max_y-min_y > max_x-min_x:
        print('Expanded by ' + str(max_y-min_y) + ' m in x direction to make square')
        max_x = min_x + (max_y - min_y)

    x_3413 = np.arange(min_x, max_x + 50.0, 50.0)
    y_3413 = np.arange(min_y, max_y + 50.0, 50.0)
    X_3413, Y_3413 = np.meshgrid(x_3413, y_3413)

    points_3413 = np.column_stack([X_3413.ravel(), Y_3413.ravel()])

    points_4326 = reproject_polygon(points_3413, inputCRS=3413, outputCRS=4326)
    Lon_4326 = points_4326[:, 0].reshape(X_3413.shape)
    Lat_4326 = points_4326[:, 1].reshape(X_3413.shape)

    points_32621 = reproject_polygon(points_4326, inputCRS=4326, outputCRS=32621)
    X_32621 = points_32621[:, 0].reshape(X_3413.shape)
    Y_32621 = points_32621[:, 1].reshape(X_3413.shape)

    return Lon_4326, Lat_4326, X_3413, Y_3413, X_32621, Y_32621

def read_bedmachine_mask(X, Y):

    bedmachine_file = '/Users/mike/Documents/Research/Data Repository/Greenland/Bathymetry/' \
                      'BedMachine/BedMachineGreenland-v5.nc'

    ds = nc4.Dataset(bedmachine_file)
    bm_x = ds.variables['x'][:]
    bm_y = ds.variables['y'][:]
    bm_mask = ds.variables['mask'][:]
    ds.close()

    buffer = 1000
    x_indices = np.logical_and(bm_x >= np.min(X)-buffer, bm_x <= np.max(X)+buffer)
    y_indices = np.logical_and(bm_y >= np.min(Y)-buffer, bm_y <= np.max(Y)+buffer)

    bm_x = bm_x[x_indices]
    bm_y = bm_y[y_indices]
    bm_X, bm_Y = np.meshgrid(bm_x, bm_y)
    bm_mask = bm_mask[np.ix_(y_indices, x_indices)]
    # print(np.shape(bm_x), np.shape(bm_y), np.shape(bm_mask))

    mask = griddata(np.column_stack([bm_X.ravel(), bm_Y.ravel()]), bm_mask.ravel(), (X, Y), method='nearest')
    return(mask)

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1,run_test = True):

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
    elif inputCRS == 4326 and str(outputCRS)[:3] == '326':
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
    elif inputCRS == 3413 and str(outputCRS)[:3] == '326':
        x2, y2 = transformer.transform(polygon_array[:,x_column], polygon_array[:,y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def subset_landsat_to_new_grid(grid, extent, dest_X_32621, dest_Y_32621):

    src_x = np.arange(extent[0], extent[1], 30)
    src_y = np.arange(extent[2], extent[3], 30)
    src_y = np.flip(src_y)
    src_X, src_Y = np.meshgrid(src_x, src_y)

    # print(np.shape(src_X), np.shape(src_Y), np.shape(grid))

    ll_dist = (src_X-np.min(dest_X_32621))**2 + (src_Y-np.min(dest_Y_32621))**2
    ll_row, ll_col = np.where(ll_dist==np.min(ll_dist))
    ll_row = ll_row[0]
    ll_col = ll_col[0]

    ul_dist = (src_X - np.max(dest_X_32621)) ** 2 + (src_Y -  np.min(dest_Y_32621)) ** 2
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))
    ul_row = ul_row[0]
    ul_col = ul_col[0]

    ur_dist = (src_X - np.min(dest_X_32621)) ** 2 + (src_Y - np.max(dest_Y_32621)) ** 2
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))
    ur_row = ur_row[0]
    ur_col = ur_col[0]

    lr_dist = (src_X - np.max(dest_X_32621)) ** 2 + (src_Y - np.max(dest_Y_32621)) ** 2
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    lr_row = lr_row[0]
    lr_col = lr_col[0]

    min_row = np.min([ll_row, lr_row, ur_row, ul_row])-100
    max_row = np.max([ll_row, lr_row, ur_row, ul_row])+100

    min_col = np.min([ll_col, lr_col, ur_col, ul_col])-100
    max_col = np.max([ll_col, lr_col, ur_col, ul_col])+100

    # print(ll_row, lr_row, ur_row, ul_row)
    # print(min_row, max_row)
    # print(ll_col, lr_col, ur_col, ul_col)
    # print(min_col, max_col)

    src_X = src_X[min_row:max_row, :]
    src_Y = src_Y[min_row:max_row, :]
    grid = grid[min_row:max_row, :]

    src_X = src_X[:, min_col:max_col]
    src_Y = src_Y[:, min_col:max_col]
    grid = grid[:, min_col:max_col]

    # print(np.shape(src_X), np.shape(src_Y), np.shape(grid))
    # print(np.shape(dest_X_32621), np.shape(dest_Y_32621))

    interp_grid = griddata(np.column_stack([src_X.ravel(), src_Y.ravel()]),
                           grid.ravel(), (dest_X_32621, dest_Y_32621),
                           method='nearest')

    return(interp_grid)

def stack_landsat_subsets_to_nc(output_file, Lon, Lat, X, Y, grids, bedmachine_mask):

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('rows',np.shape(Lon)[0])
    ds.createDimension('cols',np.shape(Lat)[1])

    x = ds.createVariable('X','f4',('rows','cols'))
    x[:, :] = X

    x = ds.createVariable('Y', 'f4', ('rows', 'cols'))
    x[:, :] = Y

    x = ds.createVariable('Lon', 'f4', ('rows', 'cols'))
    x[:, :] = Lon

    x = ds.createVariable('Lat', 'f4', ('rows', 'cols'))
    x[:, :] = Lat

    x = ds.createVariable('band_2', 'f4', ('rows', 'cols'))
    x[:, :] = grids[0]

    x = ds.createVariable('band_3', 'f4', ('rows', 'cols'))
    x[:, :] = grids[1]

    x = ds.createVariable('band_4', 'f4', ('rows', 'cols'))
    x[:, :] = grids[2]

    x = ds.createVariable('bedmachine_mask', 'f4', ('rows', 'cols'))
    x[:, :] = bedmachine_mask

    ds.close()



project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/' \
                 'Fjord/Upernavik'

fjord_file = project_folder + '/Data/Upernavik Fjord Geometry.nc'

# plt.plot(dest_X_32621[:,0], dest_Y_32621[:,0])
# plt.plot(dest_X_32621[:,-1], dest_Y_32621[:,-1])
# plt.plot(dest_X_32621[0, :], dest_Y_32621[0, :])
# plt.plot(dest_X_32621[-1, :], dest_Y_32621[-1, :])
# plt.show()

file_names = ['LC09_L2SP_017008_20250717_20250723_02_T1_SR_B2.TIF',
              'LC09_L2SP_017008_20250717_20250723_02_T1_SR_B3.TIF',
              'LC09_L2SP_017008_20250717_20250723_02_T1_SR_B4.TIF']

locations = ['Kakivfaat']#'Kakivfaat']#, , 'Upernavik_Fjord']'Ussing_Braeer'

for location in locations:

    Lon_4326, Lat_4326, X_3413, Y_3413, X_32621, Y_32621 = get_melange_region_bounds(location)

    grids = []

    bedmachine_mask = read_bedmachine_mask(X_3413, Y_3413)

    # plt.imshow(bedmachine_mask, origin='lower')
    # plt.show()

    for file_name in file_names:

        file_path = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/Data/Observations/Imagery/Landsat/'+file_name
        print(file_name)

        grid, extent = read_landsat_image(file_path)

        # plt.plot(dest_X_32621[:,0], dest_Y_32621[:,0])
        # plt.plot(dest_X_32621[:,-1], dest_Y_32621[:,-1])
        # plt.plot(dest_X_32621[0, :], dest_Y_32621[0, :])
        # plt.plot(dest_X_32621[-1, :], dest_Y_32621[-1, :])
        #
        # plt.plot([extent[0], extent[1]], [extent[2], extent[2]], 'k-')
        # plt.plot([extent[0], extent[1]], [extent[3], extent[3]], 'k-')
        # plt.plot([extent[0], extent[0]], [extent[2], extent[3]], 'k-')
        # plt.plot([extent[1], extent[1]], [extent[2], extent[3]], 'k-')
        # plt.show()

        interp_grid = subset_landsat_to_new_grid(grid, extent, X_32621, Y_32621)

        grids.append(interp_grid)

        if 'B2.TIF' in file_name:
            plt.imshow(interp_grid, origin='lower')
            plt.show()

    output_file = project_folder+'/Data/Observations/Imagery/Landsat/'+location+' Landsat Imagery.nc'
    stack_landsat_subsets_to_nc(output_file, Lon_4326, Lat_4326, X_3413, Y_3413, grids, bedmachine_mask)








