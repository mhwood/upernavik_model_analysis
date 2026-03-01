

import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_fluxCO2_timeseries(config_dir, experiment, year):

    dec_yrs = []
    fluxCO2 = []

    for month in range(2,9):
        file_name = os.path.join(config_dir, 'L2', 'L2_Upernavik','results_'+experiment,'daily_mean','fluxCO2',
                                    'fluxCO2_' + str(year) + '{:02d}'.format(month) + '.nc')
        ds = nc4.Dataset(file_name)
        grid = ds.variables['fluxCO2'][:, :, :]
        ds.close()

        # print(month, np.shape(grid))

        for day in range(1, np.shape(grid)[0]+1):
            dec_yr = YMD_to_DecYr(year, month, day)
            dec_yrs.append(dec_yr)
            # print(day)
            # if month==8 and day==7:
            #     C = plt.pcolormesh(grid[day-1, :, :])
            #     plt.colorbar(C)
            #     plt.show()
            edge_buffer = 10
            daily_flux = np.nansum(np.array(grid[day-1, edge_buffer:-edge_buffer, edge_buffer:]))
            if daily_flux<-550:
                C = plt.pcolormesh(grid[day - 1, :, :])
                plt.colorbar(C)
                plt.show()
            fluxCO2.append(daily_flux)

    dec_yrs = np.array(dec_yrs)
    fluxCO2 = np.array(fluxCO2)

    timeseries = np.column_stack((dec_yrs, fluxCO2))

    return timeseries


config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/' \
              'Greenland Model Analysis/Fjord/Upernavik'

timeseries_control_2016 = read_fluxCO2_timeseries(config_dir, experiment='control', year=2016)

plt.plot(timeseries_control_2016[:,0], timeseries_control_2016[:,1])
plt.show()










