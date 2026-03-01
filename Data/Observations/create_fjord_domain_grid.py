
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import netCDF4 as nc4
from pyproj import Transformer

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
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def create_XY_3413():

    botttom_left_x = -372026
    botttom_left_y = -1823386

    angle = -29 #degrees

    resolution = 50

    n_rows = 800
    n_cols = 1700

    X = np.zeros((n_rows, n_cols))
    Y = np.zeros((n_rows, n_cols))

    for row in range(n_rows):
        x_left = botttom_left_x + row * resolution * np.cos(np.deg2rad(1-angle))
        y_left = botttom_left_y + row * resolution * np.sin(np.deg2rad(1-angle))
        for col in range(n_cols):
            X[row, col] = x_left + col * resolution * np.cos(np.deg2rad(angle))
            Y[row, col] = y_left + col * resolution * np.sin(np.deg2rad(angle))

    # plt.plot(X[0,:], Y[0,:], 'k-')
    # plt.plot(X[-1, :], Y[-1, :], 'g-')
    # plt.show()

    local_x = resolution*np.arange(n_cols)
    local_y = resolution*np.arange(n_rows)

    return(local_x, local_y, X,Y)

def read_bedmachine_subset(dest_X, dest_Y):

    bed_machine_file = '/Users/mike/Documents/Research/Data Repository/' \
                       'Greenland/Bathymetry/BedMachine/BedMachineGreenland-v5.nc'
    ds = nc4.Dataset(bed_machine_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    bed = ds.variables['bed'][:, :]
    thickness = ds.variables['thickness'][:, :]
    surface = ds.variables['surface'][:, :]
    ds.close()

    x_indices = np.logical_and(x<np.max(dest_X)+500, x>np.min(dest_X)-500)
    y_indices = np.logical_and(y < np.max(dest_Y) + 500, y > np.min(dest_Y) - 500)

    x = x[x_indices]
    y = y[y_indices]
    X, Y = np.meshgrid(x,y)

    bed = bed[y_indices, :]
    bed = bed[:, x_indices]

    thickness = thickness[y_indices, :]
    thickness = thickness[:, x_indices]

    surface = surface[y_indices, :]
    surface = surface[:, x_indices]

    bed[thickness>0]=0
    bed[surface > 0] = 0

    interp_bed = griddata(np.column_stack([X.ravel(), Y.ravel()]), bed.ravel(), (dest_X, dest_Y))

    return(interp_bed)

def store_fjord_domain_to_grid(output_file, local_x, local_y, X, Y, lon, lat, Bed):
    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('rows', np.shape(X)[0])
    ds.createDimension('cols', np.shape(Y)[1])

    x = ds.createVariable('X', 'f4', ('rows', 'cols'))
    x[:, :] = X

    x = ds.createVariable('Y', 'f4', ('rows', 'cols'))
    x[:, :] = Y

    x = ds.createVariable('Lon', 'f4', ('rows', 'cols'))
    x[:, :] = lon

    x = ds.createVariable('Lat', 'f4', ('rows', 'cols'))
    x[:, :] = lat

    x = ds.createVariable('Bed', 'f4', ('rows', 'cols'))
    x[:, :] = Bed

    x = ds.createVariable('local_x', 'f4', ('cols', ))
    x[:] = local_x

    x = ds.createVariable('local_y', 'f4', ('rows', ))
    x[:] = local_y

    ds.close()

print('Generating the X and Y coordinates of the grid')
local_x, local_y, X, Y = create_XY_3413()

print('Reprojecting to lon/lat coords')
lon_lat = reproject_polygon(np.column_stack([X.ravel(), Y.ravel()]),
                            3413, 4326)
lon = lon_lat[:,0].reshape(np.shape(X))
lat = lon_lat[:,1].reshape(np.shape(X))

print('Generating the bed subset')
Bed = read_bedmachine_subset(X, Y)

# C = plt.pcolormesh(local_x, local_y, Bed)
# plt.colorbar(C)
# plt.show()

project_folder = '/Users/mike/Documents/Research/Projects/Iceberg Modeling'
output_file = project_folder+'/Data/Glacier/Upernavik Fjord Geometry.nc'

store_fjord_domain_to_grid(output_file, local_x, local_y, X, Y, lon, lat, Bed)















