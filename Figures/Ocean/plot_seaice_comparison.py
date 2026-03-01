
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import netCDF4 as nc4
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def retrieve_model_seaice_timeseries(config_dir, results_dir, glacier_locations):

    model_seaice_timeseries = {}
    glacier_location_index = {'NG':5,'MF':3,'CS':11}
    glacier_location_lon_lat = {}

    min_year = 2016
    max_year = 2020
    n_timesteps = 0
    for year in range(min_year,max_year+1):
        if year%4==0:
            n_timesteps += 366*4
        else:
            n_timesteps += 365*4

    for glacier_location in glacier_locations:
        model_seaice_timeseries[glacier_location] = np.zeros((n_timesteps,2))

    var_dir = os.path.join(config_dir,'L2','L2_Upernavik',results_dir,'dv','CTD','AREA')
    steps_counted = 0
    first_file = True
    for year in range(min_year, max_year+1):
        for month in range(1,13):

            if month in [1,3,5,7,8,10,12]:
                n_days = 31
            elif month in [4,6,9,11]:
                n_days = 30
            else:
                if year%4==0:
                    n_days=29
                else:
                    n_days=28

            for day in range(1,n_days+1):
                for h in range(4):
                    for glacier_location in glacier_locations:
                        model_seaice_timeseries[glacier_location][steps_counted+4*(day-1)+h,0] = YMD_to_DecYr(year,month,day,h*6)

            file_name = 'AREA_'+str(year)+'{:02d}'.format(month)+'.nc'
            if file_name in os.listdir(var_dir):
                ds = nc4.Dataset(os.path.join(var_dir, file_name))
                var_grid = ds.variables['AREA'][:, :]
                longitude = ds.variables['longitude'][:]
                latitude = ds.variables['latitude'][:]
                ds.close()

                if first_file:
                    for glacier_location in glacier_locations:
                        glacier_location_lon_lat[glacier_location] = [longitude[glacier_location_index[glacier_location]],
                                                                      latitude[glacier_location_index[glacier_location]]]
                    first_file=False

                # for i in range(len(longitude)):
                #     print(i,longitude[i], latitude[i])
                #     plt.plot(longitude[i], latitude[i], 'k.')
                #     plt.text(longitude[i], latitude[i], str(i), color='k')
                # plt.show()

                for glacier_location in glacier_locations:
                    model_seaice_timeseries[glacier_location][steps_counted:steps_counted+n_days*4,1] = var_grid[:, glacier_location_index[glacier_location]]

            steps_counted += n_days*4

    return(model_seaice_timeseries, glacier_location_lon_lat)

def read_model_grid_from_nc(config_dir):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Upernavik_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()
    return(XC,YC)

def retrieve_obs_seaice_timeseries(project_dir, glacier_locations, glacier_location_lon_lat,
                                   XC, YC):

    obs_seaice_timeseries = {}
    sample_locations = {}
    for glacier_location in glacier_locations:
        location = glacier_location_lon_lat[glacier_location]
        dist = ((XC-location[0])**2 + (YC-location[1])**2)**0.5
        row, col = np.where(dist==np.min(dist))
        sample_locations[glacier_location] = [row[0], col[0]]

    min_year = 2016
    max_year = 2020
    n_timesteps = 0
    for year in range(min_year,max_year+1):
        if year%4==0:
            n_timesteps += 366
        else:
            n_timesteps += 365

    for glacier_location in glacier_locations:
        obs_seaice_timeseries[glacier_location] = np.zeros((n_timesteps,2))

    var_dir = os.path.join(project_dir,'Data','Observations','Sea Ice')
    steps_counted = 0
    for year in range(min_year, max_year+1):

        if year%4==0:
            year_days = 366
        else:
            year_days = 365

        file_name = 'Upernavik Sea Ice ' + str(year) +'.nc'
        if file_name in os.listdir(var_dir):
            ds = nc4.Dataset(os.path.join(var_dir, file_name))
            var_grid = ds.variables['seaice_conc'][:, :, :]
            ds.close()

            for glacier_location in glacier_locations:
                obs_seaice_timeseries[glacier_location][steps_counted:steps_counted + year_days, 1] = \
                    var_grid[:,sample_locations[glacier_location][0], sample_locations[glacier_location][1]]

        for month in range(1,13):

            if month in [1,3,5,7,8,10,12]:
                n_days = 31
            elif month in [4,6,9,11]:
                n_days = 30
            else:
                if year%4==0:
                    n_days=29
                else:
                    n_days=28

            for day in range(1,n_days+1):
                for glacier_location in glacier_locations:
                    obs_seaice_timeseries[glacier_location][steps_counted+day-1,0] = YMD_to_DecYr(year,month,day)

            steps_counted += n_days

    return(obs_seaice_timeseries)

def plot_seaice_comparison(project_dir, glacier_locations, glacier_location_names,
                           model_seaice_timeseries, obs_seaice_timeseries):

    colors = ['red','orange','green','blue','purple']

    min_year = 2016
    max_year = 2017

    fig = plt.figure(figsize=(10, 7))

    gs1 = GridSpec(3,1, left=0.1, right=0.92, bottom=0.05, top=0.95, hspace=0.1)

    for g in range(len(glacier_locations)):
        ax2 = fig.add_subplot(gs1[g, :])
        ax2.plot(model_seaice_timeseries[glacier_locations[g]][:,0],
                 model_seaice_timeseries[glacier_locations[g]][:,1])
        ax2.plot(obs_seaice_timeseries[glacier_locations[g]][:, 0],
                 obs_seaice_timeseries[glacier_locations[g]][:, 1])
        ax2.set_xlim([min_year, max_year+1])
        ax2.set_ylim([-0.05, 1.05])
        if g==0:
            ax2.set_title('Upernavik Fjord Sea Ice Concentration')
        if g==1:
            ax2.set_ylabel('Sea Ice Concentration (m$^2$/m$^2$')
        ax2.text(2016.02,0,glacier_location_names[g],ha='left',va='bottom')
        ax2.grid(linestyle='--',linewidth=0.5, alpha=0.5)
        if g<2:
            ax2.set_xticklabels([])

        ax2.set_xticks(np.arange(min_year + 15 / 365.25, max_year + 1, 30 / 365.25))
        labels = []
        for year in range(min_year, max_year+1):
            labels+=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
            for i in range(12):
                time = year + i * 30 / 365.25
                plt.plot([time, time], [-0.05,1.05], 'k--', linewidth=0.5, alpha=0.5)
        if g==2:
            ax2.set_xticklabels(labels)

    output_file = os.path.join(project_dir,'Figures','Model','Upernavik Sea Ice Timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)
    a=1

def write_output_to_nc():
    a = 

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/' \
             'Downscale_Greenland/MITgcm/configurations/downscaled_greenland'

project_dir = '/Users/mike/Documents/Research/Projects/Iceberg Modeling'

results_dir = 'results_control'

glacier_locations = ['NG','MF','CS']
glacier_location_names = ['Near-Glacier','Mid-Fjord','Fjord-Mouth','Continental Shelf']

model_seaice_timeseries, glacier_location_lon_lat = retrieve_model_seaice_timeseries(config_dir, results_dir, glacier_locations)

XC, YC = read_model_grid_from_nc(config_dir)

obs_seaice_timeseries = retrieve_obs_seaice_timeseries(project_dir, glacier_locations,
                                                       glacier_location_lon_lat, XC, YC)

plot_seaice_comparison(project_dir, glacier_locations, glacier_location_names,
                           model_seaice_timeseries, obs_seaice_timeseries)