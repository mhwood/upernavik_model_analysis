
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_PFG_data_from_dv_output(config_dir, experiment):

    loc_index = 13
    depth_index = 6

    years = [2016]
    year_months = []
    for year in years:
        if year==2016:
            start_month = 2
            end_month = 9
        else:
            start_month = 1
            end_month = 13
        for month in range(start_month, end_month):
            year_months.append((year, month))

    dv_dir = os.path.join(config_dir, 'results_'+experiment, 'dv', 'CTD')

    all_timeseries = {}
    dec_yrs = []
    ptracers = list(range(29, 34))

    for ptrace in ptracers:
        all_data = []
        for (year, month) in year_months:
            file_name = os.path.join(dv_dir, f'PTRACE{ptrace:02d}',
                                     f'PTRACE{ptrace:02d}_{year:04d}{month:02d}.nc')
            ds = nc4.Dataset(file_name)
            pfg = ds.variables[f'PTRACE{ptrace:02d}'][:,:,:]
            ds.close()

            t, k, n = np.where(pfg == np.max(pfg))
            print(np.min(pfg), np.max(pfg), t, k, n)

            pfg = pfg[:, :, loc_index]
            all_data.append(pfg[:])

            if ptrace == ptracers[0]:
                file_dec_yrs = []
                for day in range(1, pfg.shape[0]+1):
                    dec_yr = YMD_to_DecYr(year, month, day)
                    file_dec_yrs.append(dec_yr)
                dec_yrs.extend(file_dec_yrs)
        all_data = np.concatenate(all_data, axis=0)

        all_timeseries[f'PTRACE{ptrace:02d}'] = all_data[:,depth_index]

    dec_yrs = np.array(dec_yrs)

    return(all_timeseries, dec_yrs)


project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

# config_dir = '/Volumes/eqipsermia/downscale_darwin/L2_Disko_Bay'
config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'MITgcm/configurations/downscale_darwin/L2/L2_Upernavik'

experiment = 'control'

all_timeseries, dec_yrs = read_PFG_data_from_dv_output(config_dir, experiment)

for ptrace in range(29,34):
    plt.plot(dec_yrs, all_timeseries['PTRACE'+str(ptrace)], label='PTRACE'+str(ptrace))
    # plt.colorbar()
plt.legend()
plt.show()

