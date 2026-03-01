
import os
import numpy as np
import netCDF4 as nc4
from scipy.interpolate import griddata
from pyproj import Transformer


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
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                  model_name + '_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    ds.close()

    return (XC, YC, Depth)

def get_domain_extents(config_dir, domain_name):
    if domain_name == 'West Greenland':
        min_x = -943593
        max_x = -111839
        min_y = -3363131
        max_y = -1095405

    elif domain_name == 'Upernavik':
        XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
        points = np.column_stack([XC.ravel(), YC.ravel()])
        reprojected_points = reproject_polygon(points, outputCRS=3413, inputCRS=4326)
        X = reprojected_points[:, 0]
        Y = reprojected_points[:, 1]
        min_x = np.min(X)
        max_x = np.max(X)
        min_y = np.min(Y)
        max_y = np.max(Y)

    extents = (min_x, max_x, min_y, max_y)
    return (extents)

def interpolate_chl_to_model_grid(Lon, Lat, X, Y):

    chl_file = '/Users/mhwood/Documents/Research/Data Repository/Global/Chlorophyll/' \
               'OC-CCI/AQUA_MODIS.20020801_20250831.L3m.MC.CHL.chlor_a.4km.nc'

    ds = nc4.Dataset(chl_file, 'r')

    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    chl = ds.variables['chlor_a'][:]
    ds.close()

    min_lon_index = np.argmin(np.abs(lon - np.min(Lon))) - 1
    max_lon_index = np.argmin(np.abs(lon - np.max(Lon))) + 1
    max_lat_index = np.argmin(np.abs(lat - np.min(Lat))) - 1
    min_lat_index = np.argmin(np.abs(lat - np.max(Lat))) + 1

    lon_subset = lon[min_lon_index:max_lon_index+1]
    lat_subset = lat[min_lat_index:max_lat_index+1]
    chl_subset = chl[min_lat_index:max_lat_index+1, min_lon_index:max_lon_index+1]

    Lon_subset, Lat_subset = np.meshgrid(lon_subset, lat_subset)
    points = np.column_stack((Lon_subset.ravel(), Lat_subset.ravel()))
    reprojected_points = reproject_polygon(points, inputCRS=4326, outputCRS=3413)

    chl_interpolated = griddata(reprojected_points, chl_subset.ravel(), (X, Y),
                                method='linear')
    return chl_interpolated

def interpolate_GEBCO_to_model_grid(Lon, Lat, X, Y):

    gebco_file = '/Users/mhwood/Documents/Research/Data Repository/Global/Bathymetry/GEBCO_2023.nc'

    ds = nc4.Dataset(gebco_file, 'r')

    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    gebco = ds.variables['elevation'][:]
    ds.close()

    min_lon_index = np.argmin(np.abs(lon - np.min(Lon))) - 1
    max_lon_index = np.argmin(np.abs(lon - np.max(Lon))) + 1
    min_lat_index = np.argmin(np.abs(lat - np.min(Lat))) - 1
    max_lat_index = np.argmin(np.abs(lat - np.max(Lat))) + 1

    # print(min_lon_index, max_lon_index, min_lat_index, max_lat_index)

    lon_subset = lon[min_lon_index:max_lon_index+1]
    lat_subset = lat[min_lat_index:max_lat_index+1]
    gebco = gebco[min_lat_index:max_lat_index+1, min_lon_index:max_lon_index+1]

    skip = 5
    lon_subset = lon_subset[::skip]
    lat_subset = lat_subset[::skip]
    gebco = gebco[::skip, ::skip]

    Lon_subset, Lat_subset = np.meshgrid(lon_subset, lat_subset)
    points = np.column_stack((Lon_subset.ravel(), Lat_subset.ravel()))
    reprojected_points = reproject_polygon(points, inputCRS=4326, outputCRS=3413)

    gebco_interpolated = griddata(reprojected_points, gebco.ravel(), (X, Y),
                                method='linear')

    gebco_interpolated = -gebco_interpolated  # Convert to depth (positive down)
    gebco_interpolated[gebco_interpolated < 0] = 0  # Set any negative depths to zero

    return gebco_interpolated

def create_global_subset(config_dir, project_dir):

    print('- Creating West Greenland Chlorophyll and Bathymetry Subset...')

    extent = get_domain_extents(config_dir, 'West Greenland')

    resolution = 1000
    x = np.arange(extent[0], extent[1]+resolution, resolution)
    y = np.arange(extent[2], extent[3]+resolution, resolution)
    X, Y = np.meshgrid(x, y)

    points_4326 = reproject_polygon(np.column_stack([X.ravel(), Y.ravel()]),
                                      inputCRS=3413, outputCRS=4326)
    Lon = points_4326[:, 0].reshape(X.shape)
    Lat = points_4326[:, 1].reshape(Y.shape)

    print('    - Interpolating Bathymetry Data...')
    bathy_interpolated = interpolate_GEBCO_to_model_grid(Lon, Lat, X, Y)

    print('    - Interpolating Chlorophyll Data...')
    chl_interpolated = interpolate_chl_to_model_grid(Lon, Lat, X, Y)


    subset_file = os.path.join(project_dir, 'Data', 'West_Greenland_Chl_Climatology.nc')

    ds = nc4.Dataset(subset_file, 'w', format='NETCDF4')
    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))

    chl_var = ds.createVariable('chlor_a', 'f4', ('y', 'x'))
    bathy_var = ds.createVariable('depth', 'f4', ('y', 'x'))
    x_var = ds.createVariable('x', 'f4', ('x',))
    y_var = ds.createVariable('y', 'f4', ('y',))

    chl_var[:, :] = chl_interpolated
    bathy_var[:, :] = bathy_interpolated
    x_var[:] = x
    y_var[:] = y
    ds.close()

def create_model_domain_subset(config_dir, project_dir):

    XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
    points = np.column_stack([XC.ravel(), YC.ravel()])
    reprojected_points = reproject_polygon(points, outputCRS=3413, inputCRS=4326)
    X = reprojected_points[:, 0].reshape(XC.shape)
    Y = reprojected_points[:, 1].reshape(YC.shape)
    x = X[0,:]
    y = Y[:,0]

    Lon = XC
    Lat = YC

    chl_interpolated = interpolate_chl_to_model_grid(Lon, Lat, X, Y)

    subset_file = os.path.join(project_dir, 'Data',
                               'L2_Upernavik_Chl_Climatology.nc')

    ds = nc4.Dataset(subset_file, 'w', format='NETCDF4')
    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))

    chl_var = ds.createVariable('chlor_a', 'f4', ('y', 'x'))
    x_var = ds.createVariable('x', 'f4', ('x',))
    y_var = ds.createVariable('y', 'f4', ('y',))

    chl_var[:, :] = chl_interpolated
    x_var[:] = x
    y_var[:] = y
    ds.close()

config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'MITgcm/configurations/downscale_darwin/'


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

create_global_subset(config_dir, project_dir)

# create_model_domain_subset(config_dir, project_dir)