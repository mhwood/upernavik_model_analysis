
import os
import argparse
import sys
import numpy as np
import netCDF4 as nc4
from pyproj import Transformer
from scipy.interpolate import griddata
# from osgeo import gdal
# from osgeo import osr


def read_extent_from_chl_chlimatology_nc(project_dir):

    ds = nc4.Dataset(os.path.join(project_dir,'Data','West_Greenland_Chl_Climatology.nc'))
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    ds.close()
    X, Y = np.meshgrid(x, y)

    return(X, Y)


def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon


def read_RGB_bands_from_nc(nc_file_path):

    ds = nc4.Dataset(nc_file_path)
    Lon = ds.variables['Longitude'][:, :]
    Lat = ds.variables['Latitude'][:, :]
    band_1 = ds.variables['sur_refl_b01'][:, :, :]
    band_1 = band_1[0, :, :]
    band_3 = ds.variables['sur_refl_b03'][:, :, :]
    band_3 = band_3[0, :, :]
    band_4 = ds.variables['sur_refl_b04'][:, :, :]
    band_4 = band_4[0, :, :]

    ds.close()

    return(Lon, Lat, band_1, band_3, band_4)


def read_MODIS_points_to_domain(modis_path,min_lon,max_lon,min_lat,max_lat,reproject_to_polar):

    file_names = ['h16v01.ncml.nc4','h16v02.ncml.nc4','h16v00.ncml.nc4',
                  'h15v01.ncml.nc4','h15v02.ncml.nc4',
                  'h17v01.ncml.nc4',]

    resolution_buffer = 0.1

    points_started = False

    for file_name in file_names:
        print('   - Reading '+file_name)
        nc_file_path = modis_path+'/'+file_name

        Lon, Lat, band_1, band_3, band_4 = read_RGB_bands_from_nc(nc_file_path)

        points = np.column_stack([np.ravel(Lon), np.ravel(Lat)])
        values_1 = np.reshape(band_1,(np.size(band_1),1)) # red
        values_3 = np.reshape(band_3,(np.size(band_3),1)) # green
        values_4 = np.reshape(band_4,(np.size(band_4),1)) # blue

        indices_lon = np.logical_and(points[:,0]>=min_lon-resolution_buffer,
                                     points[:,0]<=max_lon+resolution_buffer)
        indices_lat = np.logical_and(points[:,1]>=min_lat-resolution_buffer,
                                     points[:,1]<=max_lat+resolution_buffer)
        indices = np.logical_and(indices_lon, indices_lat)

        if np.any(indices):

            points_subset = points[indices, :]
            band_1_subset = values_1[indices]
            band_3_subset = values_3[indices]
            band_4_subset = values_4[indices]

            if not points_started:
                points_started = True
                all_points = points_subset
                all_band_1_points = band_1_subset
                all_band_3_points = band_3_subset
                all_band_4_points = band_4_subset
            else:
                all_points = np.vstack([all_points,points_subset])
                all_band_1_points = np.vstack([all_band_1_points,band_1_subset])
                all_band_3_points = np.vstack([all_band_3_points,band_3_subset])
                all_band_4_points = np.vstack([all_band_4_points,band_4_subset])

    if reproject_to_polar:
        all_points = reproject_points(all_points, inputCRS=4326, outputCRS=3413)

    return(all_points,all_band_1_points,all_band_3_points,all_band_4_points)


def interpolate_points_to_grid(points, X, Y, band_1_points, band_3_points, band_4_points):

    print('    - Interpolating band 1')
    band_1_grid = griddata(points,band_1_points.ravel(),(X,Y))
    print('    - Interpolating band 3')
    band_3_grid = griddata(points,band_3_points.ravel(),(X,Y))
    print('    - Interpolating band 4')
    band_4_grid = griddata(points,band_4_points.ravel(),(X,Y))

    return(band_1_grid, band_3_grid, band_4_grid)


def write_data_to_nc(output_file, x, y, band_1_grid, band_3_grid, band_4_grid):
    ds = nc4.Dataset(output_file, 'w', format='NETCDF4')

    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))

    x_var = ds.createVariable('x', 'f4', ('x',))
    y_var = ds.createVariable('y', 'f4', ('y',))

    band_1_var = ds.createVariable('band_1', 'f4', ('y', 'x'))
    band_3_var = ds.createVariable('band_3', 'f4', ('y', 'x'))
    band_4_var = ds.createVariable('band_4', 'f4', ('y', 'x'))

    x_var[:] = x
    y_var[:] = y

    band_1_var[:, :] = band_1_grid
    band_3_var[:, :] = band_3_grid
    band_4_var[:, :] = band_4_grid

    ds.close()


def create_tif_file(project_dir, modis_dir):

    buffer = 0.1

    reproject_to_polar = True
    resolution = 250

    # reproject_to_polar = False
    # resolution = 0.001

    # step 1: read in the geometry
    print(' - Reading in the geometry')
    X, Y = read_extent_from_chl_chlimatology_nc(project_dir)
    x = X[0, :]
    y = Y[:, 0]

    points = np.column_stack([X.ravel(), Y.ravel()])
    reprojected_points = reproject_points(points, inputCRS=3413, outputCRS=4326)
    Lon = reprojected_points[:, 0].reshape(X.shape)
    Lat = reprojected_points[:, 1].reshape(Y.shape)

    min_lon = np.min(Lon)
    max_lon = np.max(Lon)
    min_lat = np.min(Lat)
    max_lat = np.max(Lat)

    print('   - Domain Lon: '+str(min_lon)+' to '+str(max_lon))
    print('   - Domain Lat: '+str(min_lat)+' to '+str(max_lat))

    # step 3: read in the modis points
    print(' - Reading in the MODIS data')
    points, band_1_points, band_3_points, band_4_points = \
        read_MODIS_points_to_domain(modis_dir, min_lon, max_lon, min_lat, max_lat, reproject_to_polar)

    # step 4: interpolate the modis points onto the grid
    print(' - Interpolating the points onto the domain')
    band_1_grid, band_3_grid, band_4_grid = \
        interpolate_points_to_grid(points, X, Y, band_1_points, band_3_points, band_4_points)

    # step 6: output the files to nc
    print(' - Outputting the bands to nc')
    output_file = os.path.join(project_dir, 'Data','Observations', 'West_Greenland_MOIDS_Imagery.nc')
    write_data_to_nc(output_file, x, y, band_1_grid, band_3_grid, band_4_grid)


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

modis_dir = '/Users/mhwood/Documents/Research/Data Repository/Greenland/MODIS'

create_tif_file(project_dir, modis_dir)

