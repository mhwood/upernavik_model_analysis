

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4



def read_annual_chl_obs_grid(project_dir, year):
    chl_file = os.path.join(project_dir, 'Data','Observations', 'Chlorophyll',
                            'Upernavik_Chl_Observations_'+str(year)+'.nc')
    ds = nc4.Dataset(chl_file)
    chl = ds.variables['chlor_a'][:, :,:]
    X = ds.variables['X'][:,:]
    Y = ds.variables['Y'][:,:]
    dec_yrs = ds.variables['time'][:]
    ds.close()
    return chl, X, Y, dec_yrs

def read_annual_chl_model_timeseries(project_dir, experiment, year, chl_number):
    chl_file = os.path.join(project_dir, 'Data','Models', 'Chlorophyll',
                            'Chl'+'{:02d}'.format(chl_number)+'_'+experiment+'_median_timeseries.nc')
    ds = nc4.Dataset(chl_file)
    chl = ds.variables['Chl'+'{:02d}'.format(chl_number)][:]
    dec_yrs = ds.variables['time'][:]
    ds.close()
    timeseries = np.column_stack((dec_yrs, chl))
    timeseries = timeseries[timeseries[:,0]!=0,:]
    return timeseries

def compute_chl_timeseries(chl, dec_yrs, window_size=5):

    chl_timeseries = np.zeros((len(dec_yrs),))

    for t in range(len(dec_yrs)):
        chl_snapshot = np.array(chl[t, :, :])
        chl_snapshot[chl_snapshot>20]=np.nan
        chl_mean = np.nanmean(chl_snapshot)
        chl_timeseries[t] = chl_mean
    return chl_timeseries





project_dir = '/Users/mhwood/Documents/Research/Projects/' \
              'Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

chl_obs, X, Y, dec_yrs_2021 = read_annual_chl_obs_grid(project_dir, year=2021)
chl_obs_timeseries_2021 = compute_chl_timeseries(chl_obs, dec_yrs_2021, window_size=5)

chl_model_timeseries = read_annual_chl_model_timeseries(project_dir, experiment='baseline', year=2021, chl_number=2)

plt.figure(figsize=(10,5))

# plt.plot(dec_yrs_2021-2021, chl_obs_timeseries_2021, 'o')
plt.plot(chl_model_timeseries[:,0]-2021, chl_model_timeseries[:,1], 'o')

plt.xlabel('Decimal Year')
plt.ylabel('Mean Chlorophyll (mg m$^{-3}$)')
plt.show()
