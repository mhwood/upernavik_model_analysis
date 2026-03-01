


import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta, date
import gsw
from scipy.interpolate import interp1d

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    this_date = datetime(year,month,day,hour,minute,second)
    start = date(this_date.year, 1, 1).toordinal()
    year_length = date(this_date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(this_date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_float_profiles(project_folder, float_ID):

    dec_yrs = np.zeros((1000,))
    profiles = []

    counter = 0
    for file_name in sorted(os.listdir(os.path.join(project_folder,'Data','Ocean','Float Profiles', float_ID))):
        if file_name[0]!='.' and file_name[-3:]=='.nc':
            ds = nc4.Dataset(os.path.join(project_folder,'Data','Ocean','Float Profiles', float_ID,file_name))
            depth = ds.variables['depth'][:]
            temp = ds.variables['potential_temperature'][:]
            salt = ds.variables['practical_salinity'][:]
            ds.close()
            profiles.append(np.column_stack([depth, temp, salt]))

            year = int(file_name.split('_')[0][:4])
            month = int(file_name.split('_')[0][4:6])
            day = int(file_name.split('_')[0][6:8])
            dec_yrs[counter] = YMD_to_DecYr(year,month,day)

            counter+=1

    dec_yrs = dec_yrs[dec_yrs!=0]

    return(profiles, dec_yrs)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def read_L2_CTD_dv_float_output(config_dir):

    months = np.arange(2,12).tolist()

    for month in months:
        theta_file = os.path.join(config_dir,'L2','L2_Upernavik','results_baseline', 'dv', 'CTD',
                                   'THETA', 'THETA_2016'+'{:02d}'.format(month)+'.nc')
        ds = nc4.Dataset(theta_file)
        depth = ds.variables['depths'][:]
        Theta = ds.variables['THETA'][:, :, :]
        longitude = ds.variables['longitude'][:]
        latitude = ds.variables['latitude'][:]
        ds.close()

        salt_file = os.path.join(config_dir,'L2','L2_Upernavik','results_baseline', 'dv', 'CTD',
                                 'SALT', 'SALT_2016'+'{:02d}'.format(month)+'.nc')
        ds = nc4.Dataset(salt_file)
        Salt = ds.variables['SALT'][:, :, :]
        ds.close()

        dec_yrs = []
        for day in range(1,np.shape(Theta)[0]+1):
            dec_yr = YMD_to_DecYr(2016, month, day)
            dec_yrs.append(dec_yr)
        dec_yrs = np.array(dec_yrs)

        if month == months[0]:
            theta_timeseries = Theta
            salt_timeseries = Salt
            all_dec_yrs = dec_yrs
        else:
            theta_timeseries = np.concatenate([theta_timeseries, Theta], axis=0)
            salt_timeseries = np.concatenate([salt_timeseries, Salt], axis=0)
            all_dec_yrs = np.concatenate([all_dec_yrs, dec_yrs], axis=0)

    lon = -57.89070
    lat = 73.257333
    dist = great_circle_distance(lon, lat, longitude.ravel(), latitude.ravel())
    min_dist_index = np.where(dist==np.min(dist))[0][0]
    print('Closest model grid point is %.1f m away from float location' % np.min(dist))

    theta_timeseries = theta_timeseries[:, :, min_dist_index]
    salt_timeseries = salt_timeseries[:, :, min_dist_index]

    return(depth, all_dec_yrs, theta_timeseries, salt_timeseries)

def format_axes(ax, min_dec_yr, max_dec_yr, letter='a'):
    ax.set_ylim([999, 0])
    ax.set_xlim([min_dec_yr, max_dec_yr])
    # ax.plot([np.min(dec_yrs_F10052), np.min(dec_yrs_F10052)], [500, 0], 'k--')
    ax.set_xticklabels([])

    month_dec_yrs = []
    month_labels = []
    year = int(np.floor(min_dec_yr))
    for yr in range(year, int(np.ceil(max_dec_yr))):
        for month in range(1,13):
            dec_yr = YMD_to_DecYr(yr,month,1)
            if dec_yr>=min_dec_yr and dec_yr<=max_dec_yr:
                month_dec_yrs.append(dec_yr)
                if month==1:
                    month_labels.append("'"+str(yr)[-2:])
                else:
                    month_labels.append(datetime(yr,month,1).strftime('%b')[0])

    ax.set_xticks(month_dec_yrs)
    ax.set_xticklabels(month_labels, rotation=0, ha='center', fontsize=8)

    for year in range(int(np.ceil(min_dec_yr)), int(np.ceil(max_dec_yr))):
        ax.plot([year, year], [500, 0], 'k--', linewidth=0.5, alpha=0.5)

    # for i in range(12):
    #     time = year + i * 30 / 365.25
    #     plt.plot([time, time], [500, 0], 'k--', linewidth=0.5, alpha=0.5)

    ax.text(year+0.98, 10, letter+')', ha='right', va='top',
             fontsize=16, color='black')

def plot_float_data(project_folder,
                    profiles_F9186, dec_yrs_F9186,
                    model_depth, model_dec_yrs,
                    theta_timeseries_control, salt_timeseries_control
                    ):

    fig = plt.figure(figsize=(10,10))

    plot_rows = 5
    plot_cols = 5
    gs = GridSpec(plot_rows*5+2, plot_cols*2, left=0.1, right=0.95, bottom=0.05, top=0.95, hspace=0.3)

    ##########################################################################################

    ax11 = fig.add_subplot(gs[:plot_rows, :plot_cols])
    tmin = -1.5
    tmax = 4
    for p in range(len(profiles_F9186)):
        t0 = dec_yrs_F9186[p]
        t1 = dec_yrs_F9186[p]+10/365.25
        time = np.column_stack([t0, t1])
        profile = np.column_stack([profiles_F9186[p][:,1], profiles_F9186[p][:,1]])
        depth = profiles_F9186[p][:,0]
        ax11.pcolormesh(time,depth,profile,cmap='turbo', vmin=tmin, vmax=tmax)

    ax11.text(2021 + 10 / 365.25, 475, '2023 (Argo Float Profiles)', color='black', ha='left', va='bottom',
              bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    # ax11.text(np.min(dec_yrs_F10052)-10/365.25, 475, 'Float F9186', color='black', ha='right', va='bottom',
    #         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    # ax11.text(np.min(dec_yrs_F10052) + 10 / 365.25, 475, 'Float F10052', color='black', ha='left', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    ax11.set_title('Temperature')
    format_axes(ax11, 2021, 2022, letter='a')


    ##########################################################################################

    ax12 = fig.add_subplot(gs[:plot_rows, plot_cols:])
    smin = 32
    smax = 35
    for p in range(len(profiles_F9186)):
        t0 = dec_yrs_F9186[p]
        t1 = dec_yrs_F9186[p] + 10 / 365.25

        time = np.column_stack([t0, t1])
        profile = np.column_stack([profiles_F9186[p][:, 2], profiles_F9186[p][:, 2]])
        depth = profiles_F9186[p][:, 0]
        ax12.pcolormesh(time, depth, profile, cmap='viridis', vmin=smin, vmax=smax)

    # ax12.text(np.min(dec_yrs_F10052) - 10 / 365.25, 750, 'Float F9186', color='black', ha='right', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    # ax12.text(np.min(dec_yrs_F10052) + 10 / 365.25, 750, 'Float F10052', color='black', ha='left', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax12.set_title('Salinity')
    ax12.set_yticklabels([])
    format_axes(ax12, 2021, 2022, letter='b')

    ##########################################################################################
    # 2008

    ax21 = fig.add_subplot(gs[plot_rows:2*plot_rows, :plot_cols])
    # print(np.shape(dec_yrs_2008), np.shape(d))
    ax21.pcolormesh(model_dec_yrs, model_depth, theta_timeseries_control.T, cmap='turbo', vmin=tmin, vmax=tmax)
    format_axes(ax21, min_dec_yr=2016, max_dec_yr=2017, letter='c')
    ax21.text(2016 + 10 / 365.25, 475, '2008 (Model)', color='black', ha='left', va='bottom',
              bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))

    ax22 = fig.add_subplot(gs[plot_rows:2 * plot_rows, plot_cols:])
    ax22.pcolormesh(model_dec_yrs, model_depth, salt_timeseries_control.T, cmap='viridis', vmin=smin, vmax=smax)
    format_axes(ax22, min_dec_yr=2016, max_dec_yr=2017, letter='d')
    ax22.set_yticklabels([])

    # ##########################################################################################
    # # 2012
    #
    # ax31 = fig.add_subplot(gs[2*plot_rows:3 * plot_rows, :plot_cols])
    # ax31.pcolormesh(dec_yrs_2012, model_depth, Theta_2012, cmap='turbo', vmin=tmin, vmax=tmax)
    # format_axes(ax31, year=2012, letter='e')
    # ax31.text(2012 + 10 / 365.25, 475, '2012 (Model)', color='black', ha='left', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    #
    # ax31.set_ylabel('Depth (m)')
    #
    # ax32 = fig.add_subplot(gs[2*plot_rows:3 * plot_rows, plot_cols:])
    # ax32.pcolormesh(dec_yrs_2012, model_depth, Salt_2012, cmap='viridis', vmin=smin, vmax=smax)
    # format_axes(ax32, year=2012, letter='f')
    # ax32.set_yticklabels([])
    #
    # ##########################################################################################
    # # 2017
    #
    # ax41 = fig.add_subplot(gs[3*plot_rows:4 * plot_rows, :plot_cols])
    # ax41.pcolormesh(dec_yrs_2017, model_depth, Theta_2017, cmap='turbo', vmin=tmin, vmax=tmax)
    # format_axes(ax41, year=2017,letter='g')
    # ax41.text(2017 + 10 / 365.25, 475, '2017 (Model)', color='black', ha='left', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    #
    # ax42 = fig.add_subplot(gs[3*plot_rows:4 * plot_rows, plot_cols:])
    # ax42.pcolormesh(dec_yrs_2017, model_depth, Salt_2017, cmap='viridis', vmin=smin, vmax=smax)
    # format_axes(ax42, year=2017,letter='h')
    # ax42.set_yticklabels([])
    #
    # ##########################################################################################
    # # 2019
    #
    # ax51 = fig.add_subplot(gs[4*plot_rows:5 * plot_rows, :plot_cols])
    # ax51.pcolormesh(dec_yrs_2019, model_depth, Theta_2019, cmap='turbo', vmin=tmin, vmax=tmax)
    # format_axes(ax51, year=2019,letter='i')
    # ax51.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    # ax51.text(2019 + 10 / 365.25, 475, '2019 (Model)', color='black', ha='left', va='bottom',
    #           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
    #
    # ax52 = fig.add_subplot(gs[4*plot_rows:5 * plot_rows, plot_cols:])
    # ax52.pcolormesh(dec_yrs_2019, model_depth, Salt_2019, cmap='viridis', vmin=smin, vmax=smax)
    # format_axes(ax52, year=2019,letter='j')
    # ax52.set_yticklabels([])
    # ax52.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])

    ##########################################################################################

    ax11c = fig.add_subplot(gs[-1, 1:4])
    y = np.array([0, 1])
    x = np.linspace(tmin,tmax,100)
    X,Y = np.meshgrid(x,y)
    ax11c.pcolormesh(x,y,X, cmap='turbo')
    ax11c.set_yticks([])
    ax11c.set_xlabel('Temperature ($^{\circ}$C)')

    ax12c = fig.add_subplot(gs[-1, plot_cols+1:plot_cols+4])
    y = np.array([0, 1])
    x = np.linspace(smin, smax, 100)
    X, Y = np.meshgrid(x, y)
    ax12c.pcolormesh(x, y, X, cmap='viridis', vmin=smin, vmax=smax)
    ax12c.set_yticks([])
    ax12c.set_xlabel('Salinity (psu)')


    output_file = os.path.join(project_folder,'Figures','Ocean','Upernavik Float Comparison.png')
    plt.savefig(output_file)
    plt.close(fig)




project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/MITgcm/' \
             'configurations/downscale_darwin'


profiles_F9186, dec_yrs_F9186 = read_float_profiles(project_folder, float_ID = 'F9186')

# iceplume = False
# dec_yrs_2008, depth, Theta_2008 = read_model_data(project_folder, 2008, 'Theta', iceplume)
# dec_yrs_2012, _, Theta_2012 = read_model_data(project_folder, 2012, 'Theta', iceplume)
# dec_yrs_2017, _, Theta_2017 = read_model_data(project_folder, 2017, 'Theta', iceplume)
# dec_yrs_2019, _, Theta_2019 = read_model_data(project_folder, 2019, 'Theta', iceplume)
# 
# _, model_depth, Salt_2008 = read_model_data(project_folder, 2008, 'Salt', iceplume)
# _, _, Salt_2012 = read_model_data(project_folder, 2012, 'Salt', iceplume)
# _, _, Salt_2017 = read_model_data(project_folder, 2017, 'Salt', iceplume)
# _, _, Salt_2019 = read_model_data(project_folder, 2019, 'Salt', iceplume)

model_depth, model_dec_yrs, theta_timeseries_control, salt_timeseries_control =\
    read_L2_CTD_dv_float_output(config_dir)

plot_float_data(project_folder,
                profiles_F9186, dec_yrs_F9186,
                model_depth, model_dec_yrs,
                theta_timeseries_control, salt_timeseries_control)



