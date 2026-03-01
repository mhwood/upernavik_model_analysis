

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


def compute_timeseries(results_dir, metric, year, month, hFaC):

    # min row, max row, min col, max col
    bounds = {'Upernavik':[0, 100, 350, 450],
              'Kakivfaat':[150, 230, 270, 450],
              'Ussing Braeer':[230, 375, 270, 450]}

    all_month_timeseries = {}
    for glacier in list(bounds.keys()):
        all_month_timeseries[glacier] = {}

    for var_name in ['fwflx', 'htflx', 'mltrt']:

        subset_dict = {'fwflx':'BRGfwFlx',
                       'htflx':'BRGhtFlx',
                       'mltrt':'BRGmltRt'}

        file_path = os.path.join(results_dir, 'daily_mean', subset_dict[var_name],
                                 subset_dict[var_name] + '_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_path)
        grid = ds.variables[subset_dict[var_name]][:, :, :, :]
        iterations = ds.variables['iterations'][:]
        ds.close()

        for glacier in list(bounds.keys()):
            print('         - Working on '+var_name+' for glacier '+glacier)
            min_row, max_row, min_col, max_col = bounds[glacier]
            timeseries = np.zeros((np.shape(grid)[0],))

            for t in range(np.shape(grid)[0]):
                level_set = grid[t, :, :, :]
                level_set = level_set[:,min_row:max_row, min_col:max_col]
                if var_name in ['fwflx', 'htflx']:
                    timeseries[t] = np.sum(level_set)
                else:
                    if metric=='mean':
                        timeseries[t] = np.mean(level_set[level_set!=0])
                    elif metric=='median':
                        timeseries[t] = np.median(level_set[level_set!=0])

            all_month_timeseries[glacier][var_name] = timeseries

    return(iterations,all_month_timeseries)


def write_timeseries_to_nc(project_dir, metric, experiment, dec_yrs, all_timeseries):

    output_file = os.path.join(project_dir, 'Data', 'Models', 'Melange',
                               'Melange_'+experiment+ '_'+metric+'_timeseries.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time',len(dec_yrs))

    ivar = ds.createVariable('time','f4',('time',))
    ivar[:] = dec_yrs

    for glacier in list(all_timeseries.keys()):

        grp = ds.createGroup(glacier)

        for var_name in list(all_timeseries[glacier].keys()):
            timeseries = all_timeseries[glacier][var_name]

            vvar = grp.createVariable(var_name, 'f4', ('time',))
            vvar[:] = timeseries

    ds.close()


config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/' \
             'Downscale_Darwin/darwin3/configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'


model_name = 'L2_Upernavik'
metric = 'mean'

XC, YC, Z, Depth, hFaC = read_grid_geometry_from_nc(config_dir, model_name)

years = [2016]

experiments = ['melange']
for experiment in experiments:

    print('  - Workin on year '+experiment)
    results_dir = os.path.join(config_dir, 'L2', model_name, 'results_' + experiment)

    n_days = 0
    for year in years:
        if year%4==0:
            n_days+=366
        else:
            n_days+=365
        if year==2016:
            n_days-=31 # no Jan this year

    dec_yrs = np.zeros((n_days,))

    all_timeseries = {'Upernavik':{'fwflx':np.zeros((n_days,)), 'htflx':np.zeros((n_days,)), 'mltrt':np.zeros((n_days,))},
                      'Kakivfaat':{'fwflx':np.zeros((n_days,)), 'htflx':np.zeros((n_days,)), 'mltrt':np.zeros((n_days,))},
                      'Ussing Braeer':{'fwflx':np.zeros((n_days,)), 'htflx':np.zeros((n_days,)), 'mltrt':np.zeros((n_days,))}}

    counter = 0
    for year in years:
        if year==2016:
            start_month = 2
        else:
            start_month = 1
        for month in range(start_month,13):
            print('     - Reading in month '+str(month)+' in year '+str(year))

            month_iterations, all_month_timeseries = compute_timeseries(results_dir, metric, year, month, hFaC)
            # plt.plot(month_iterations, month_timeseries)
            # plt.show()

            for glacier in list(all_month_timeseries.keys()):
                for var_name in list(all_month_timeseries[glacier].keys()):
                    month_timeseries = all_month_timeseries[glacier][var_name]
                    # print(glacier, var_name, np.mean(month_timeseries))
                    all_timeseries[glacier][var_name][counter:counter+np.shape(month_timeseries)[0]] = month_timeseries

            dec_yrs[counter:counter+np.shape(month_timeseries)[0]] = [YMD_to_DecYr(year,month,day) for day in range(1,np.shape(month_timeseries)[0]+1)]
            # iterations[counter:counter+np.shape(month_timeseries)[0]] = month_iterations

            counter += np.shape(month_timeseries)[0]

    # plt.plot(dec_yrs[:counter], timeseries[:counter])
    # plt.show()

    write_timeseries_to_nc(project_dir, metric, experiment, dec_yrs, all_timeseries)





