

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def iter_number_to_date(iter_number,seconds_per_iter,start_year=1992):

    total_seconds = iter_number*seconds_per_iter
    date = dt.datetime(start_year,1,1) + dt.timedelta(seconds=total_seconds)
    return(date)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    drF = ds.variables['drF'][:]
    hFaC = ds.variables['HFacC'][:, :, :]
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, Z, Depth, hFaC)


def compute_Chl_map(project_dir, year, month):

    start_day_index = 0
    end_day_index = 0
    for m in range(1,13):
        if m in [1,3,5,7,8,10,12]:
            days_in_month = 31
        elif month in [4,6,9,11]:
            days_in_month = 30
        elif m == 2:
            if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
                days_in_month = 29
            else:
                days_in_month = 28

        if m < month:
            start_day_index += days_in_month
            end_day_index += days_in_month
        if m==month:
            end_day_index += days_in_month
            break

    print('Start day index:', start_day_index, 'End day index:', end_day_index, 'for month:', month, 'and year:', year)

    file_path = os.path.join(project_dir, 'Data','Observations','Chlorophyll',
                             'Upernavik_Chl_Observations_'+str(year)+'.nc')
    ds = nc4.Dataset(file_path)
    grid = ds.variables['chlor_a'][:, :, :]
    grid = grid[start_day_index:end_day_index, :, :]
    grid[grid>100] = np.nan
    grid = np.nanmean(grid, axis=0)
    Lon = ds.variables['lon'][:]
    Lat = ds.variables['lat'][:]
    X = ds.variables['X'][:]
    Y = ds.variables['Y'][:]
    ds.close()

    chl_map = grid

    return(Lon, Lat, X, Y, chl_map)


def write_map_to_nc(project_dir, year, month, Lon, Lat, X, Y, chl_map):

    output_file = os.path.join(project_dir, 'Data', 'Observations', 'Chlorophyll',
                               'Chl_observations_mean_map_'+str(year)+'{:02d}'.format(month)+'.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('cols', X.shape[1])
    ds.createDimension('rows', X.shape[0])

    ivar = ds.createVariable('longitude','f4',('rows','cols'))
    ivar[:] = Lon

    vvar = ds.createVariable('latitude', 'f4', ('rows','cols'))
    vvar[:] = Lat

    xvar = ds.createVariable('X', 'f4', ('rows','cols'))
    xvar[:] = X

    yvar = ds.createVariable('Y', 'f4', ('rows','cols'))
    yvar[:] = Y

    mvar = ds.createVariable(var_name, 'f4', ('rows','cols'))
    mvar[:] = chl_map


    ds.close()


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

model_name = 'L2_Upernavik'
var_name = 'Chl'

year = 2021


experiments = ['control']
for experiment in experiments:
    for month in range(1,12):

        Lon, Lat, X, Y, chl_map = compute_Chl_map(project_dir, year, month)

        C = plt.pcolormesh(chl_map)
        plt.colorbar(C)
        plt.title('Chlorophyll Mean Map - '+experiment+ ' - '+str(year)+'/{:02d}'.format(month))
        plt.show()

        write_map_to_nc(project_dir, year, month, Lon, Lat, X, Y, chl_map)





