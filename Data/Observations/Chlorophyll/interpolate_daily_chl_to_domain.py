
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from pyproj import Transformer
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

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

def read_model_domain(config_dir, model_name):
    domain_file = os.path.join(config_dir, 'nc_grids', model_name+'_grid.nc')
    ds = nc4.Dataset(domain_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()
    return(XC, YC, Depth)

def interpolate_obs_to_domain(chl_dir, year, XC, YC, X, Y, Depth):

    first_file = True

    dec_yrs = []
    chl_grids = []

    for month in range(1,13):
        print('Adding in data for month '+str(month)+' in year '+str(year))
        if month in [1, 3, 5, 7, 8, 10, 12]:
            n_days = 31
        elif month in [4, 6, 9, 11]:
            n_days = 30
        else:
            if year % 4 == 0:
                n_days = 29
            else:
                n_days = 28
        for day in range(1,n_days+1):
            file_path = os.path.join(chl_dir, str(year),
                                     f'ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-{year}{month:02d}{day:02d}-fv6.0.nc')
            if os.path.exists(file_path):
                ds = nc4.Dataset(file_path)
                chl = ds.variables['chlor_a'][:, :, :]
                if first_file:
                    lon = ds.variables['lon'][:]
                    lat = ds.variables['lat'][:]
                ds.close()

                chl = chl[0, :, :]

                if first_file:

                    ll_buffer = 1
                    lon_indices = np.where((lon >= np.min(XC) - ll_buffer) &
                                              (lon <= np.max(XC) + ll_buffer))[0]
                    lat_indices = np.where((lat >= np.min(YC) - ll_buffer) &
                                                (lat <= np.max(YC) + ll_buffer))[0]
                    lon = lon[lon_indices]
                    lat = lat[lat_indices]

                    lon_grid, lat_grid = np.meshgrid(lon, lat)
                    points = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
                    reprojected_points = reproject_points(points, inputCRS=4326, outputCRS=3413)
                    X_obs = reprojected_points[:, 0].reshape(lon_grid.shape)
                    Y_obs = reprojected_points[:, 1].reshape(lat_grid.shape)

                    first_file = False

                chl = chl[lat_indices, :][:, lon_indices]
                chl_points = np.column_stack([X_obs.ravel(), Y_obs.ravel()])
                chl_values = chl.ravel()
                chl_interpolated = griddata(chl_points, chl_values,
                                            (X, Y), method='linear')

                dec_yr = YMD_to_DecYr(year, month, day)
                dec_yrs.append(dec_yr)
                chl_grids.append(chl_interpolated)

                # plt.pcolormesh(X,Y,chl_interpolated, cmap='turbo', vmin=0, vmax=5)
                # plt.show()

    dec_yrs = np.array(dec_yrs)
    chl_grids = np.array(chl_grids)
    return(dec_yrs, chl_grids)

def output_to_nc(output_file, dec_yrs, chl_grids, XC, YC, X, Y):

    ds = nc4.Dataset(output_file, 'w', format='NETCDF4')

    ds.createDimension('time', len(dec_yrs))
    ds.createDimension('y', YC.shape[0])
    ds.createDimension('x', XC.shape[1])

    time_var = ds.createVariable('time', 'f4', ('time',))
    y_var = ds.createVariable('Y', 'f4', ('y','x'))
    x_var = ds.createVariable('X', 'f4', ('y','x'))
    lon_var = ds.createVariable('lon', 'f4', ('y', 'x'))
    lat_var = ds.createVariable('lat', 'f4', ('y', 'x'))

    chl_var = ds.createVariable('chlor_a', 'f4', ('time', 'y', 'x'))

    time_var[:] = dec_yrs
    y_var[:] = Y
    x_var[:] = X
    lon_var[:] = XC
    lat_var[:] = YC

    chl_var[:, :, :] = chl_grids

    ds.close()

project_dir = '/Users/mhwood/Documents/Research/Projects/' \
              'Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/' \
             'Downscaled_Darwin/darwin3/configurations/downscaled_ecco_v5_darwin'


chl_dir = '/Users/mhwood/Desktop/Chl'


XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')

points = np.column_stack([XC.ravel(), YC.ravel()])
reprojected_points = reproject_points(points, inputCRS=4326, outputCRS=3413)
X = reprojected_points[:, 0].reshape(XC.shape)
Y = reprojected_points[:, 1].reshape(YC.shape)

for year in range(2021,2022):
    dec_yrs, chl_grid = interpolate_obs_to_domain(chl_dir, year, XC, YC, X, Y, Depth)

    output_file = os.path.join(project_dir, 'Data','Observations','Chlorophyll',
                               'Upernavik_Chl_Observations_'+str(year)+'.nc')
    output_to_nc(output_file, dec_yrs, chl_grid, XC, YC, X, Y)