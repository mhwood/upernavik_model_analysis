

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


def compute_Chl_map(results_dir, var_name, year, month, depth_index, hFaC):


    chl_grid = np.zeros_like(hFaC[0,:,:])

    for chl_number in range(1,6):
        var_name = 'Chl0'+str(chl_number)

        file_path = os.path.join(results_dir, 'daily_mean', var_name, var_name + '_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[var_name][:, :, :, :]
        grid = grid[:,depth_index,:,:]
        grid = np.mean(grid, axis=0)
        chl_grid += grid
        ds.close()

    chl_map = chl_grid/5

    return(chl_map)


def write_map_to_nc(project_dir, year, month, depth, experiment, var_name, XC, YC, chl_map):

    output_file = os.path.join(project_dir, 'Data', 'Models', 'Chlorophyll',
                               var_name +'_'+experiment+ '_mean_map_'+str(year)+'{:02d}'.format(month)+'_'+str(int(depth))+'m.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('cols', XC.shape[1])
    ds.createDimension('rows', XC.shape[0])

    ivar = ds.createVariable('longitude','f4',('rows','cols'))
    ivar[:] = XC

    vvar = ds.createVariable('latitude', 'f4', ('rows','cols'))
    vvar[:] = YC

    mvar = ds.createVariable(var_name, 'f4', ('rows','cols'))
    mvar[:] = chl_map

    ds.depth = depth

    ds.close()


config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/' \
             'Downscale_Darwin/darwin3/configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'


model_name = 'L2_Upernavik'
var_name = 'Chl'

XC, YC, Z, Depth, hFaC = read_grid_geometry_from_nc(config_dir, model_name)

year = 2016
month = 8
depth_index = 6

depth = Z[depth_index]

experiments = ['control']
for experiment in experiments:
    for month in range(6,10):

        print('  - Workin on year '+experiment+ ' month '+str(month))
        results_dir = os.path.join(config_dir, 'L2', model_name, 'results_' + experiment)

        chl_map = compute_Chl_map(results_dir, var_name, year, month, depth_index, hFaC)

        C = plt.pcolormesh(chl_map)
        plt.colorbar(C)
        plt.title('Chlorophyll Mean Map - '+experiment+ ' - '+str(year)+'/{:02d}'.format(month)+' - Depth: '+str(int(depth))+'m')
        plt.show()

        write_map_to_nc(project_dir, year, month, depth, experiment, var_name, XC, YC, chl_map)





