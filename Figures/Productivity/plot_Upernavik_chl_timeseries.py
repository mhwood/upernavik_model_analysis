

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


def compute_chl_obs_timeseries(chl, dec_yrs, window_size=5):

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

chl, X, Y, dec_yrs_2016 = read_annual_chl_obs_grid(project_dir, year=2016)
chl_timeseries_2016 = compute_chl_obs_timeseries(chl, dec_yrs_2016, window_size=5)

chl, X, Y, dec_yrs_2017 = read_annual_chl_obs_grid(project_dir, year=2017)
chl_timeseries_2017 = compute_chl_obs_timeseries(chl, dec_yrs_2017, window_size=5)

chl, X, Y, dec_yrs_2018 = read_annual_chl_obs_grid(project_dir, year=2018)
chl_timeseries_2018 = compute_chl_obs_timeseries(chl, dec_yrs_2018, window_size=5)

chl, X, Y, dec_yrs_2019 = read_annual_chl_obs_grid(project_dir, year=2019)
chl_timeseries_2019 = compute_chl_obs_timeseries(chl, dec_yrs_2019, window_size=5)

plt.figure(figsize=(10,5))
plt.plot(dec_yrs_2016-2016, chl_timeseries_2016, '-o')
plt.plot(dec_yrs_2017-2017, chl_timeseries_2017, '-o')
plt.plot(dec_yrs_2018-2018, chl_timeseries_2018, '-o')
plt.plot(dec_yrs_2019-2019, chl_timeseries_2019, '-o')

plt.xlabel('Decimal Year')
plt.ylabel('Mean Chlorophyll (mg m$^{-3}$)')
plt.show()
