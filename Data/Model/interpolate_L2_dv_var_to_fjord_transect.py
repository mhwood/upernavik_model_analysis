
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
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

def read_polygons_from_nc(project_dir, fjord_name):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Models','Fjord Transects',
                                  fjord_name, fjord_name+' Fjord Geometry Transect.nc'))

    distance = ds.variables['distance'][:]
    surface = ds.variables['surface'][:]
    surface*=-1
    # distance *=1e-3

    transect_x = ds.variables['transect_x'][:]
    transect_y = ds.variables['transect_y'][:]

    ice_polygon_x = ds.variables['ice_polygon_x'][:]
    ice_polygon_y = ds.variables['ice_polygon_y'][:]
    ice_polygon = np.column_stack([ice_polygon_x, ice_polygon_y])
    # ice_polygon[:,0]*=1e-3
    ice_polygon[:,1]*=-1

    bed_polygon_x = ds.variables['bathy_polygon_x'][:]
    bed_polygon_y = ds.variables['bathy_polygon_y'][:]
    bed_polygon = np.column_stack([bed_polygon_x, bed_polygon_y])
    # bed_polygon[:,0]*=1e-3
    bed_polygon[:,1]*=-1

    return(distance, transect_x, transect_y, surface, ice_polygon, bed_polygon)

def read_fjord_data_from_nc(config_dir, fjord_name, experiment, year, month, var_name):

    file_name = os.path.join(config_dir,'L2','L2_Upernavik', 'results_'+experiment, 'dv', '_'.join(fjord_name.split()),
                                var_name, var_name + '_' + str(year) + f'{month:02d}.nc')
    ds = nc4.Dataset(file_name)
    depth = ds.variables['depths'][:]
    longitude = ds.variables['longitude'][:]
    latitude = ds.variables['latitude'][:]
    grid = ds.variables[var_name][:, :, :]
    grid = np.nanmean(grid, axis=0)
    ds.close()

    return(depth, longitude, latitude, grid)

def interpolate_var_to_transect(depth, longitude, latitude, grid,
                                transect_x, transect_y, distance, surface):

    points = np.column_stack((longitude, latitude))
    reprojected_points = reproject_polygon(points, 4326, 3413, x_column=0, y_column=1)
    x = reprojected_points[:,0]
    y = reprojected_points[:,1]

    num_depths = len(depth)
    num_transect_points = len(transect_x)
    transect_var = np.full((num_depths, num_transect_points), np.nan)

    dist_thresh = 500 * np.sqrt(2) # 500 m grid cells
    for n in range(num_transect_points):
        dist = np.sqrt((x - transect_x[n])**2 + (y - transect_y[n])**2)
        local_indices = np.where(dist <= dist_thresh)[0]
        # if n==20:
        #     for ind in local_indices:
        #         plt.plot(grid[:,ind], depth, 'b-', markersize=1)
        #     plt.show()

        # if n==0:
        #     print(longitude[local_indices], latitude[local_indices])

        # print(n, len(transect_x), len(local_indices), local_indices, [transect_x[n], transect_y[n]])
        for d in range(num_depths):
            points = np.column_stack((x, y))
            values = grid[d, :].ravel()
            points = points[local_indices, :]
            values = values[local_indices]
            transect_point = np.array([[transect_x[n], transect_y[n]]])
            # if n==0:
                # print(d, depth[d], points, values, transect_point)
            if  len(values)<3: # np.any(values==0) or
                continue
            else:
                interp_value = griddata(points, values, transect_point, method='nearest')
                transect_var[d, n] = interp_value

    return(transect_var)

def save_interpolated_transect_to_nc(project_dir, fjord_name, experiment, distance, depth, transect_var, var_name, year, month):
    output_file = os.path.join(project_dir, 'Data','Models','Fjord Transects',fjord_name,
                               '_'.join(fjord_name.split())+f'_Fjord_{experiment}_{var_name}_Transect_{year}{month:02d}.nc')
    ds = nc4.Dataset(output_file, 'w', format='NETCDF4')

    ds.createDimension('distance', len(distance))
    ds.createDimension('depth', len(depth))

    distance_var = ds.createVariable('distance', 'f4', ('distance',))
    depth_var = ds.createVariable('depth', 'f4', ('depth',))
    transect_var_nc = ds.createVariable(var_name, 'f4', ('depth', 'distance'))

    distance_var[:] = distance
    depth_var[:] = depth
    transect_var_nc[:, :] = transect_var
    ds.close()

# project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'
project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'


# config_dir = '/Volumes/eqipsermia/downscale_darwin/L2_Disko_Bay'
# config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
#              'MITgcm/configurations/downscale_darwin/'
config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'



year = 2020
month = 11


XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
points = np.column_stack([XC.ravel(), YC.ravel()])
reprojected_points = reproject_polygon(points, 4326, 3413, x_column=0, y_column=1)
X = reprojected_points[:,0].reshape(XC.shape)
Y = reprojected_points[:,1].reshape(YC.shape)

for fjord_name in ['Upernavik N','Kakivfaat','Ussing Braeer']:
    distance, transect_x, transect_y, surface, ice_polygon, bed_polygon = read_polygons_from_nc(project_dir, fjord_name)

    for experiment in ['baseline','baseline_melange','baseline_iceplume','baseline_melange_iceplume']:
        for var_name in ['THETA','SALT','PTRACE28']:
            depth, longitude, latitude, grid = read_fjord_data_from_nc(config_dir, fjord_name, experiment, year, month, var_name)

            transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
            x, y = transformer.transform(latitude.ravel(), longitude.ravel())

            transect_var = interpolate_var_to_transect(depth, longitude, latitude, grid,
                                                        transect_x, transect_y, distance, surface)

            save_interpolated_transect_to_nc(project_dir, fjord_name, experiment, distance, depth, transect_var, var_name, year, month)





# # Plotting for verification
# plt.figure(figsize=(8,6))
# plt.contourf(distance, -depth, transect_var, levels=20, cmap='viridis')
# plt.colorbar(label=var_name)
# # plt.plot(ice_polygon[:,0], ice_polygon[:,1], 'k-', linewidth=2, label='Ice Shelf')
# plt.plot(bed_polygon[:,0], -1*bed_polygon[:,1], 'k--', linewidth=2, label='Bathy')
# # plt.plot(distance, surface, 'w-', linewidth=2, label='Surface')
# plt.ylim([-1000, 0])
# plt.xlabel('Distance along transect (km)')
# plt.ylabel('Depth (m)')
# plt.title(f'{var_name} along Upernavik N Fjord Transect - {year}-{month:02d}')
# # plt.legend()
# plt.show()