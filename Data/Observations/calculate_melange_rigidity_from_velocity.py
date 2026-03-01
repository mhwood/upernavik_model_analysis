
import os
import numpy as np
import matplotlib.pyplot as plt
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


def read_melange_fraction_from_file(project_dir, region):

    file_path = os.path.join(project_dir, 'Data', 'Observations','Melange', region+'_Landsat_8_Melange_Fraction.nc')

    ds = nc4.Dataset(file_path, 'r')
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    melange_fraction = ds.variables['melange_fraction'][:,:]
    ds.close()

    melange_fraction[melange_fraction < 0.3] = 0

    X, Y = np.meshgrid(x, y)
    upper_left_indices = np.logical_and((X < -340000), (Y > -1.8e6))
    melange_fraction[upper_left_indices] = 0

    return x, y, melange_fraction

def read_velocity_from_file(project_dir, X, Y, Depth):

    file_path = os.path.join(project_dir, 'Data', 'Observations','Melange', 'melange_CW_Feb2026.nc')

    ds = nc4.Dataset(file_path, 'r')
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    velocity = ds.variables['velocity'][:,:,:]
    time = ds.variables['time'][:]
    ds.close()

    min_date = 365-3 # Jan 1 of 2016 in this file
    max_date = min_date + 365*4 +366 # Jan 1 of 2021 in this file
    # time = time[(time >= min_date) & (time <= max_date)]
    velocity = velocity[(time >= min_date) & (time <= max_date),:,:]
    print(np.shape(velocity))

    # choose 10 random timeserpes and subset the velocit grid on those
    # random_time_indices = np.random.choice(velocity.shape[0], size=10, replace=False)
    # velocity = velocity[random_time_indices,:,:]

    # vel_mean = np.nanmean(velocity, axis=0)
    # plt.pcolormesh(x, y, vel_mean)
    # plt.colorbar(label='Mean Velocity (m/yr)')
    # plt.contour(X, Y, Depth, levels=[1], colors='k', linewidths=2)
    # plt.contour(X, Y, Depth, levels=[1], colors='w', linewidths=1)
    # plt.xlim([np.min(x), np.max(x)])
    # plt.ylim([np.min(y), np.max(y)])
    # plt.show()

    return x, y, velocity

def compute_melange_rigidity(region,melange_x, melange_y, melange_fraction,
                           velocity_x, velocity_y, velocity):

    melange_X, melange_Y = np.meshgrid(melange_x, melange_y)
    velocity_X, velocity_Y = np.meshgrid(velocity_x, velocity_y)

    rigidity_mask = np.zeros(np.shape(melange_fraction))
    N = 0

    for time in range(velocity.shape[0]):
        if time%10==0:
            print(f'Processing time step {time+1}/{velocity.shape[0]} for '+region)
        velocity_slice = velocity[time,:,:]
        velocity_slice = ~np.isnan(velocity_slice)
        velocity_slice = velocity_slice.astype(float)
        velocity_slice = griddata(np.column_stack([velocity_X.ravel(), velocity_Y.ravel()]),
                                  velocity_slice.ravel(), (melange_X, melange_Y), method='nearest')
        # if time==5:
        #     plt.pcolormesh(melange_x, melange_y, velocity_slice)
        #     plt.colorbar(label='Velocity Slice')
        #     plt.show()
        rigidity_mask += velocity_slice
        N += 1

    rigidity_mask /= N

    # plt.pcolormesh(melange_x, melange_y, rigidity_mask)
    # plt.colorbar(label='Rigidity Mask')
    # plt.show()

    return rigidity_mask

def save_rigidity_mask_as_nc(project_dir, region, rigidity_mask, rigidity_mask_fraction, x, y):
    output_file = os.path.join(project_dir, 'Data', 'Observations', 'Melange',
                               region+'_Melange_Rigidity_Mask.nc')

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))
    x_var = ds.createVariable('x', 'f4', ('x',))
    x_var[:] = x
    y_var = ds.createVariable('y', 'f4', ('y',))
    y_var[:] = y
    melange_var = ds.createVariable('rigidity_mask', 'f4', ('y', 'x',))
    melange_var[:, :] = rigidity_mask
    melange_fraction_var = ds.createVariable('rigidity_mask_fraction', 'f4', ('y', 'x',))
    melange_fraction_var[:, :] = rigidity_mask_fraction
    ds.close()


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'


config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/Downscaled_Darwin/darwin3/configurations/downscaled_ecco_v5_darwin'

XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')
points = np.column_stack([XC.ravel(), YC.ravel()])
reprojected_points = reproject_polygon(points, outputCRS=3413, inputCRS=4326)
X = reprojected_points[:, 0].reshape(XC.shape)
Y = reprojected_points[:, 1].reshape(YC.shape)

print('Reading velocity data from file...')
velocity_x, velocity_y, velocity = read_velocity_from_file(project_dir, X, Y, Depth)

for region in ['Upernavik_Fjord']:#, 'Kakivfaat', 'Ussing_Braeer']:
    print('Processing region: '+region)

    melange_x, melange_y, melange_fraction = read_melange_fraction_from_file(project_dir, region)

    rigidity_mask_fraction = compute_melange_rigidity(region, melange_x, melange_y, melange_fraction,
                                             velocity_x, velocity_y, velocity)

    rigidity_mask = np.zeros_like(rigidity_mask_fraction)
    rigidity_mask[rigidity_mask_fraction > 0.2] = 1
    rigidity_mask[melange_fraction == 0] = 0

    save_rigidity_mask_as_nc(project_dir, region, rigidity_mask, rigidity_mask_fraction, melange_x, melange_y)

    # plt.subplot(1, 3, 1)
    # plt.pcolormesh(melange_x, melange_y, melange_fraction)
    # plt.colorbar(label='Melange Fraction')
    #
    # plt.subplot(1, 3, 2)
    # plt.pcolormesh(melange_x, melange_y, rigidity_mask_fraction)
    # plt.colorbar(label='Melange Fraction')
    #
    # plt.subplot(1,3,3)
    # plt.pcolormesh(melange_x, melange_y, rigidity_mask)
    # plt.colorbar(label='Rigidity Mask')
    # plt.show()