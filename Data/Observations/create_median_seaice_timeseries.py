

import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
import datetime


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_annual_mean_seaice_conc_timeseries(project_dir, year):

    data_file = os.path.join(project_dir,'Data','Observations','Sea Ice','Upernavik Sea Ice '+str(year)+'.nc')
    ds = nc4.Dataset(data_file)
    seaice = ds.variables['seaice_conc'][:,:,:]
    ds.close()

    timeseries = np.zeros((np.shape(seaice)[0],))
    for t in range(len(timeseries)):
        seaice_subset = seaice[t,:,:]
        timeseries[t] = np.median(seaice_subset[seaice_subset<=1])

    dec_yrs = np.linspace(year, year+1, np.shape(seaice)[0]+1)
    dec_yrs = dec_yrs[:-1]

    timeseries = np.column_stack([dec_yrs, timeseries])

    # plt.plot(timeseries[:,0], timeseries[:,1])
    # plt.show()

    return(timeseries)

def write_annual_timeseries_to_nc(project_dir, year, timeseries):

    ds = nc4.Dataset(os.path.join(project_dir,'Data','Observations','Sea Ice',
                                  'Upernavik Sea Ice Concentration Timeseries '+str(year)+'.nc'),'w')

    ds.createDimension('dec_yrs',np.shape(timeseries)[0])

    d = ds.createVariable('dec_yrs','f4',('dec_yrs',))
    d[:] = timeseries[:,0]

    s = ds.createVariable('median_conc', 'f4', ('dec_yrs',))
    s[:] = timeseries[:, 1]

    ds.close()

def write_timeseries_to_nc(project_dir, full_timeseries):

    ds = nc4.Dataset(os.path.join(project_dir,'Data','Observations','Sea Ice','Upernavik Sea Ice Concentration Timeseries.nc'),'w')

    ds.createDimension('dec_yrs',np.shape(full_timeseries)[0])

    d = ds.createVariable('dec_yrs','f4',('dec_yrs',))
    d[:] = full_timeseries[:,0]

    s = ds.createVariable('median_conc', 'f4', ('dec_yrs',))
    s[:] = full_timeseries[:, 1]

    ds.close()

    a=1



project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

first_file = True

for year in range(2016,2021):
    print('Reading year '+str(year))

    annual_timeseries = read_annual_mean_seaice_conc_timeseries(project_dir, year)

    # write_annual_timeseries_to_nc(project_dir, year, annual_timeseries)

    if first_file:
        full_timeseries = annual_timeseries
        first_file = False
    else:
        full_timeseries = np.vstack([full_timeseries, annual_timeseries])

write_timeseries_to_nc(project_dir, full_timeseries)

plt.plot(full_timeseries[:, 0], full_timeseries[:, 1])
plt.show()



