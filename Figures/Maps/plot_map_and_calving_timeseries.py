

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import interp1d
from matplotlib.gridspec import GridSpec
import datetime

def read_flux_timeseries_from_nc(project_dir,glacier_names):

    file_path = os.path.join(project_dir, 'Data', 'Glacier', 'Upernavik Ice Flux Timeseries.nc')
    ds = nc4.Dataset(file_path)
    dec_yrs = ds.variables['dec_yr'][:]

    all_timeseries = []
    for g in range(len(glacier_names)):
        all_timeseries.append(np.column_stack([dec_yrs,
                                               ds.variables[glacier_names[g]][:]]))
    ds.close()

    for t in range(len(all_timeseries)):
        if t==0:
            total_flux_timeseries = np.copy(all_timeseries[t])
        else:
            total_flux_timeseries[:,1] += all_timeseries[t][:,1]
    total_flux_timeseries = total_flux_timeseries[total_flux_timeseries[:,0]>=2016,:]
    total_flux_timeseries = total_flux_timeseries[total_flux_timeseries[:, 0] < 2021, :]
    total_flux_timeseries[:,1] *= 1e12/917

    cumulative_flux_timeseries = np.copy(total_flux_timeseries)
    cumulative_flux_timeseries[:,1]=0
    for year in range(2016,2021):
        indices = np.where(np.floor(cumulative_flux_timeseries[:,0])==year)[0]
        for index in indices:
            if index!=indices[0]:
                cumulative_flux_timeseries[index,1] = cumulative_flux_timeseries[index-1,1] +\
                    total_flux_timeseries[index,1]*(total_flux_timeseries[index,0]-total_flux_timeseries[index-1,0])
            else:
                cumulative_flux_timeseries[index, 1] = 0

    return(all_timeseries, cumulative_flux_timeseries)

def read_calving_timeseries_from_nc(project_dir,glacier_names):
    years = [2016, 2017, 2018, 2019, 2020]

    file_path = os.path.join(project_dir, 'Data', 'Glacier', 'Upernavik Calving Timeseries.nc')
    ds = nc4.Dataset(file_path)

    all_timeseries = []
    for g in range(len(glacier_names)):
        grp = ds.groups[glacier_names[g]]
        all_timeseries.append(np.column_stack([grp.variables['dec_yr'][:],
                                               grp.variables['volume'][:]]))
        print(glacier_names[g],'{:.2e}'.format(np.sum(grp.variables['volume'][:])/5))
    ds.close()

    month_timeseries_template = np.zeros(((years[-1] - years[0] + 1) * 12, 3))
    counter = 0
    for year in range(years[0], years[-1] + 1):
        for month in range(1, 13):
            if month < 12:
                month_timeseries_template[counter, 0] = YMD_to_DecYr(year, month, 1)
                month_timeseries_template[counter, 1] = YMD_to_DecYr(year, month + 1, 1)
            else:
                month_timeseries_template[counter, 0] = YMD_to_DecYr(year, month, 1)
                month_timeseries_template[counter, 1] = YMD_to_DecYr(year + 1, 1, 1)
            counter += 1

    monthly_calving_timeseries = []
    for t in range(len(all_timeseries)):
        monthly_timeseries = np.copy(month_timeseries_template)
        for m in range(np.shape(monthly_timeseries)[0]):
            indices = np.logical_and(all_timeseries[t][:, 0] >= monthly_timeseries[m, 0],
                                     all_timeseries[t][:, 0] < monthly_timeseries[m, 1])
            monthly_timeseries[m, 2] = np.sum(all_timeseries[t][indices,1])
        monthly_calving_timeseries.append(monthly_timeseries)


    total_calving_timeseries = np.zeros((len(years)*365,2))
    total_calving_timeseries[:,0] = np.linspace(years[0],years[-1]+1,len(years)*365)
    for t in range(len(all_timeseries)):
        for d in range(np.shape(all_timeseries[t])[0]):
            total_calving_timeseries[np.argmin(np.abs(total_calving_timeseries[:,0]-all_timeseries[t][d,0])),1] += all_timeseries[t][d,1]
    print('Total', '{:.2e}'.format(np.sum(total_calving_timeseries[:,1]) / 5))

    cumulative_calving_timeseries = np.copy(total_calving_timeseries)
    cumulative_calving_timeseries[:,1]=0
    for year in years:
        indices = np.where(np.floor(cumulative_calving_timeseries[:,0])==year)[0]
        for index in indices:
            if index!=indices[0]:
                cumulative_calving_timeseries[index,1] = cumulative_calving_timeseries[index-1,1] +\
                    total_calving_timeseries[index,1]
            else:
                cumulative_calving_timeseries[index, 1] = 0

    return(all_timeseries, cumulative_calving_timeseries, monthly_calving_timeseries)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def model_timeseries_to_dec_yrs(timesteps):
    dec_yrs = np.zeros_like(timesteps)
    for t in range(len(timesteps)):
        date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=timesteps[t])
        dec_yrs[t] = YMD_to_DecYr(date.year, date.month, date.day)
    return(dec_yrs)

def collect_monthly_calving_timrseries(config_dir, calving_location_numbers):

    n_timesteps_in_schedule = 10000
    n_years_in_schedule = 8

    start_year = 2016
    end_year = 2020

    month_timeseries_template = np.zeros(((end_year-start_year+1)*12, 3))
    counter = 0
    for year in range(start_year, end_year+1):
        for month in range(1, 13):
            if month<12:
                month_timeseries_template[counter, 0] = YMD_to_DecYr(year, month,1)
                month_timeseries_template[counter, 1] = YMD_to_DecYr(year, month+1, 1)
            else:
                month_timeseries_template[counter, 0] = YMD_to_DecYr(year, month, 1)
                month_timeseries_template[counter, 1] = YMD_to_DecYr(year+1, 1, 1)
            counter +=1

    monthly_calving_timeseries = []
    total_calving_timeseries = np.zeros(((end_year-start_year+1)*365,2))
    for year in range(start_year,end_year+1):
        total_calving_timeseries[365*(year-start_year):365*(year-start_year+1),0] = np.linspace(year,year+1-1/365,365)

    for g in range(len(calving_location_numbers)):
        monthly_timeseries = np.copy(month_timeseries_template)
        for n in calving_location_numbers[g]:
            calving_file = os.path.join(config_dir,'input','calving_schedules',
                                        'calving_schedule_'+'{:03d}'.format(n))
            calving_schedule = np.fromfile(calving_file,'>f8').reshape((3,n_timesteps_in_schedule*n_years_in_schedule)).T
            calving_schedule = calving_schedule[calving_schedule[:,0]!=0,:]
            calving_schedule[:,0] = model_timeseries_to_dec_yrs(calving_schedule[:,0])
            volume = 1.62*calving_schedule[:,1]*calving_schedule[:,1]*calving_schedule[:,2]

            # make a cumulative timeseries of calving
            for c in range(len(calving_schedule)):
                if calving_schedule[c,0]>=start_year and calving_schedule[c,0]<=end_year+1:
                    index = np.argmin(np.abs(total_calving_timeseries[:,0]-calving_schedule[c,0]))
                    total_calving_timeseries[index,1]+=volume[c]

            for m in range(np.shape(monthly_timeseries)[0]):
                indices = np.logical_and(calving_schedule[:,0]>=monthly_timeseries[m,0],
                                         calving_schedule[:,0]<monthly_timeseries[m,1])
                monthly_timeseries[m,2] = np.sum(volume[indices])
        monthly_calving_timeseries.append(monthly_timeseries)

    cumulative_calving_timeseries = np.copy(total_calving_timeseries)
    cumulative_calving_timeseries[:,1]=0
    for year in range(start_year,end_year+1):
        indices = np.where(np.floor(cumulative_calving_timeseries[:,0])==year)[0]
        for index in indices:
            if index!=indices[0]:
                cumulative_calving_timeseries[index,1] = cumulative_calving_timeseries[index-1,1] +\
                    total_calving_timeseries[index,1]*(total_calving_timeseries[index,0]-total_calving_timeseries[index-1,0])
            else:
                cumulative_calving_timeseries[index, 1] = 0

    # plt.plot(monthly_timeseries[:,0], monthly_timeseries[:,2])
    # plt.show()

    return(monthly_calving_timeseries, cumulative_calving_timeseries)

def plot_map_and_timeseires(project_dir, glacier_names,
                            flux_timeseries, cumulative_flux_timeseries,
                            monthly_calving_volume_timeseries, cumulative_calving_timeseries):

    colors = ['red','orange','green','blue','purple']

    fig = plt.figure(figsize=(8, 8))

    gs1 = GridSpec(4,2, left=0.1, right=0.92, bottom=0.05, top=0.95, hspace=0.05)

    ax2 = fig.add_subplot(gs1[1, :])
    min_flux = 0
    max_flux = 8
    for g in range(len(glacier_names)):
        ax2.plot(flux_timeseries[g][:,0], flux_timeseries[g][:,1],'-',
                 color=colors[g],label=glacier_names[g].split('_')[-1])
    ax2.legend(ncols=5)
    ax2.set_xlim([2016,2021])
    ax2.set_ylim([min_flux,max_flux])
    ax2.set_ylabel('Mass Flux (Gt/yr)')
    ax2t=ax2.twinx()
    ax2t.set_ylabel('Volume Flux (km$^3$/yr)')
    ax2t.set_ylim([0,max_flux*1e12/917/1e9])
    ax2.set_xticklabels([])

    ax3 = fig.add_subplot(gs1[2,:])
    for g in range(len(glacier_names)):
        timeseries = monthly_calving_volume_timeseries[g]
        if g==0:
            bottom=np.zeros((np.shape(timeseries)[0],))
        else:
            bottom+=monthly_calving_volume_timeseries[g-1][:,2]/1e9
        # print(g, np.shape(timeseries))
        ax3.bar(timeseries[:,0], timeseries[:,2]/1e9,
                width = timeseries[:,1]-timeseries[:,0], align='edge',
                color = colors[g], bottom=bottom, edgecolor='k',linewidth=0.5)
    ax3.set_xlim([2016, 2021])
    ax3.set_ylabel('Monthly Calving\nVolume (km$^3$)')
    ax3.set_xticklabels([])

    ax4 = fig.add_subplot(gs1[3,:])
    cumulative_flux_timeseries[cumulative_flux_timeseries[:,1]==0, 1] = np.nan
    cumulative_calving_timeseries[cumulative_calving_timeseries[:, 1] == 0, 1] = np.nan
    ax4.plot(cumulative_flux_timeseries[:,0], cumulative_flux_timeseries[:,1]/1e9,'-',
             label='Total Ice Flux',color='grey')
    ax4.plot(cumulative_calving_timeseries[:,0],cumulative_calving_timeseries[:,1]/1e9,'k-',
             label='Total Calving Flux')
    ax4.legend()
    ax4.set_ylabel('Annually Integrated\nFlux (km$^3$)')
    ax4.grid(linestyle='--',linewidth=0.5, alpha=0.5)
    ax4.set_xlim([2016, 2021])

    output_file = os.path.join(project_dir,'Figures','Glacier','Glacier Map and Calving Timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)
    a=1

project_dir = '/Users/mike/Documents/Research/Projects/Iceberg Modeling'

glacier_names = ['Upernavik_Isstrom_NW',
                 'Upernavik_Isstrom_N',
                 'Upernavik_Isstrom_C',
                 'Upernavik_Isstrom_S',
                 'Upernavik_Isstrom_SS']

calving_location_numbers = [np.arange(1,2),
                            np.arange(11,12),
                            np.arange(21,22),
                            np.arange(31,32),
                            np.arange(41,42)]

flux_timeseries, cumulative_flux_timeseries = read_flux_timeseries_from_nc(project_dir,glacier_names)

_, cumulative_calving_timeseries, monthly_calving_timeseries = read_calving_timeseries_from_nc(project_dir,glacier_names)


plot_map_and_timeseires(project_dir, glacier_names,
                        flux_timeseries, cumulative_flux_timeseries,
                        monthly_calving_timeseries, cumulative_calving_timeseries)