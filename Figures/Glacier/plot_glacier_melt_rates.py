
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime as dt
from matplotlib.gridspec import GridSpec

def read_model_domain(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                  model_name + '_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    drF = ds.variables['drF'][:]
    ds.close()

    Z_bottom = np.cumsum(drF)
    # print(Z_bottom)
    Z_top = Z_bottom - drF
    Z = 0.5 * (Z_top + Z_bottom)

    return (XC, YC, Depth, Z)

def get_glacier_locations():
    glacier_locations = {'Upernavik SS': {'latitude': 72.79292, 'longitude': -54.27084},
                         'Upernavik S': {'latitude': 72.84225, 'longitude': -54.39244},
                         'Upernavik C': {'latitude': 72.93910, 'longitude': -54.38256},
                         'Upernavik N': {'latitude': 73.01690, 'longitude': -54.38086},
                         'Upernavik NW': {'latitude': 73.02248, 'longitude': -54.52449},
                         'Nunatakassaap': {'latitude': 73.22435, 'longitude': -55.18467},
                         'Kakivfaat': {'latitude': 73.49262, 'longitude': -55.49989},
                         'Qeqertarsuup': {'latitude': 73.57625, 'longitude': -55.58233},
                         'Ussing Braeer': {'latitude': 73.84979, 'longitude': -55.63595},
                         'Ussing Braeer N': {'latitude': 73.94423, 'longitude': -55.76761},
                         'Cornell': {'latitude': 74.23043, 'longitude': -56.11842},
                         'Cornell N': {'latitude': 74.27986, 'longitude': -56.12629}}
    return glacier_locations

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

def Date_to_DecYr(date):

    start_of_year = dt.datetime(date.year, 1, 1, 0, 0, 0)
    start_of_next_year = dt.datetime(date.year + 1, 1, 1, 0, 0, 0)

    year_length = (start_of_next_year - start_of_year).total_seconds()
    seconds_into_year = (date - start_of_year).total_seconds()

    decimal_fraction = seconds_into_year / year_length
    dec_yr = date.year + decimal_fraction

    return dec_yr

def read_glacier_Qsg(config_dir, glaciers):

    glacier_locations = get_glacier_locations()

    XC, YC, Depth, Z = read_model_domain(config_dir, 'L2_Upernavik')

    Qsg_file = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'input', 'iceplume', 'L2_Qsg_modified_2020')
    Qsg = np.fromfile(Qsg_file, '>f4').reshape(366, 375, 450)

    dec_yrs = []
    for i in range(Qsg.shape[0]):
        date = dt.datetime(2020, 1, 1) + dt.timedelta(days=i)
        dec_yrs.append(Date_to_DecYr(date))
    dec_yrs = np.array(dec_yrs)

    all_glacier_Qsg_timeseries = {}

    for glacier in glaciers:

        distance = great_circle_distance(glacier_locations[glacier]['longitude'], glacier_locations[glacier]['latitude'],
                                        XC, YC)
        glacier_row, glacier_col = np.where(distance == np.min(distance))
        glacier_Qsg = Qsg[:, glacier_row, glacier_col]
        glacier_Qsg_timeseries = np.column_stack((dec_yrs, glacier_Qsg))
        all_glacier_Qsg_timeseries[glacier] = glacier_Qsg_timeseries

    return all_glacier_Qsg_timeseries

def read_Theta_Salt_near_plume_location(config_dir, glaciers):
    glacier_locations = get_glacier_locations()

    XC, YC, Depth, Z = read_model_domain(config_dir, 'L2_Upernavik')

    theta_file = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'results_baseline', 'monthly_snapshot','Theta', 'Theta_202011.nc')
    ds = nc4.Dataset(theta_file)
    Theta = ds.variables['Theta'][:, :, :]
    ds.close()

    salt_file = os.path.join(config_dir, 'L2', 'L2_Upernavik', 'results_baseline', 'monthly_snapshot','Salt', 'Salt_202011.nc')
    ds = nc4.Dataset(salt_file)
    Salt = ds.variables['Salt'][:, :, :]
    ds.close()

    for glacier in glaciers:

        distance = great_circle_distance(glacier_locations[glacier]['longitude'], glacier_locations[glacier]['latitude'],
                                        XC, YC)
        glacier_row, glacier_col = np.where(distance == np.min(distance))

        salt_profile = Salt[0,:, glacier_row[0], glacier_col[0]]
        theta_profile = Theta[0,:, glacier_row[0], glacier_col[0]]
        depth = Depth[glacier_row[0], glacier_col[0]]

        fig = plt.figure(figsize=(9,4))
        plt.subplot(1,2,1)
        plt.plot(salt_profile, Z, 'b-')
        plt.plot([np.min(salt_profile), np.max(salt_profile)], [depth, depth], 'k--')
        plt.title(glacier)
        plt.xlim([np.min(salt_profile[salt_profile>0]), np.max(salt_profile[salt_profile>0])])
        plt.ylim([depth+50, 0])

        plt.subplot(1,2,2)
        plt.plot(theta_profile, Z, 'r-')
        plt.plot([np.min(theta_profile), np.max(theta_profile)], [depth, depth], 'k--')
        plt.xlim([np.min(theta_profile[theta_profile>0]), np.max(theta_profile[theta_profile>0])])
        plt.ylim([depth+50, 0])

        plt.show()



    a=1

def read_melt_timeseries_from_nc(config_dir, experiment):

    glacier_melt_rates = {}

    ds = nc4.Dataset(os.path.join(config_dir,'L2','L2_Upernavik',
                                  'results_'+experiment, 'dv', 'iceplume',
                                  'ICEFRNTM', 'ICEFRNTM_202011.nc'))
    depth = ds.variables['depths'][:]
    longitude = ds.variables['longitude'][:]
    latitude = ds.variables['latitude'][:]
    melt_rate = ds.variables['ICEFRNTM'][:, :, :]
    ds.close()

    non_zero_locations = np.sum(melt_rate[0,:,:], axis=0) != 0

    melt_rate = melt_rate[:, :, non_zero_locations]
    longitude = longitude[non_zero_locations]
    latitude = latitude[non_zero_locations]

    glacier_locations = get_glacier_locations()

    dec_yrs = []
    for i in range(melt_rate.shape[0]):
        date = dt.datetime(2020, 11, 1) + dt.timedelta(hours=i*6)
        dec_yrs.append(Date_to_DecYr(date))
    dec_yrs = np.array(dec_yrs)
    # print(dec_yrs)

    for i in range(np.shape(melt_rate)[2]):
        for glacier in glacier_locations.keys():
            distance = great_circle_distance(glacier_locations[glacier]['longitude'], glacier_locations[glacier]['latitude'],
                                            longitude[i], latitude[i])
            if distance < 1500:
                print(f'Glacier: {glacier}, Distance: {distance} m')
                glacier_melt_rates[glacier] = melt_rate[:, :, i]

    return(depth, dec_yrs, glacier_melt_rates)

def plot_glacier_melt_rates(project_dir, glacier, all_glacier_Qsg_timeseries,
                            dec_yrs, depth,
                            melt_rates_iceplume, melt_rates_melange_iceplume):

    #max_depth_index = np.sum(np.max(melt_rates_iceplume[glacier].T, axis=1) > 0)
    max_depth_index = 42

    m_min = 0
    m_max = np.max(melt_rates_iceplume[glacier])
    cmap = 'turbo'

    d_min = -0.5
    d_max = 0.5
    d_cmap = 'seismic'

    fig = plt.figure(figsize=(8,10))

    gs = GridSpec(5,21, left=0.1, right=0.9, top=0.95, bottom=0.05)

    ax = fig.add_subplot(gs[0, :-2])
    ax.plot(all_glacier_Qsg_timeseries[glacier][:,0], all_glacier_Qsg_timeseries[glacier][:,1], 'k-')
    ax.set_ylabel('Subglacial Discharge (m$^3$/s)')
    ax.set_title(f'{glacier} Subglacial Discharge and Ice Front Melt Rates')


    ax = fig.add_subplot(gs[1,:-2])
    ax.pcolormesh(dec_yrs, depth, melt_rates_iceplume[glacier].T, vmin=m_min, vmax=m_max, cmap=cmap)
    ax.set_ylim([depth[max_depth_index], 0])

    ax = fig.add_subplot(gs[2, :-2])
    ax.pcolormesh(dec_yrs, depth, melt_rates_melange_iceplume[glacier].T, vmin=m_min, vmax=m_max, cmap=cmap)
    ax.set_ylim([depth[max_depth_index], 0])

    axc = fig.add_subplot(gs[1:3, -1])
    cx = np.array([0, 1])
    cy = np.linspace(m_min, m_max, 100)
    CX, CY = np.meshgrid(cx, cy)
    c = axc.pcolormesh(CX, CY, CY, vmin=m_min, vmax=m_max, cmap=cmap)
    axc.set_xticks([])
    axc.set_ylabel('Melt Rate (m/day)')
    # put axis on right
    axc.yaxis.set_label_position("right")
    axc.yaxis.tick_right()

    ax = fig.add_subplot(gs[3, :-2])
    ax.pcolormesh(dec_yrs, depth, melt_rates_melange_iceplume[glacier].T -melt_rates_iceplume[glacier].T,
                  vmin=d_min, vmax=d_max, cmap=d_cmap)
    ax.set_ylim([depth[max_depth_index], 0])

    axc = fig.add_subplot(gs[3, -1])
    cx = np.array([0, 1])
    cy = np.linspace(d_min, d_max, 100)
    CX, CY = np.meshgrid(cx, cy)
    c = axc.pcolormesh(CX, CY, CY, vmin=d_min, vmax=d_max, cmap=d_cmap)
    axc.set_xticks([])
    axc.set_ylabel('Melt Rate Difference (m/day)')
    # put axis on right
    axc.yaxis.set_label_position("right")
    axc.yaxis.tick_right()

    max_melt_rate_timeseries = np.max(melt_rates_iceplume[glacier], axis=1)
    max_melt_rate_timeseries_melange = np.max(melt_rates_melange_iceplume[glacier], axis=1)

    ax = fig.add_subplot(gs[4, :-2])
    ax.plot(dec_yrs, max_melt_rate_timeseries, 'b-', label='Ice Plume')
    ax.plot(dec_yrs, max_melt_rate_timeseries_melange, 'r-', label='Melange + Ice Plume')


    plt.savefig(os.path.join(project_dir, 'Figures', 'Ocean', 'Ice Front', "_".join(glacier.split())+'_melt_rates.png'), dpi=300)
    plt.close(fig)





config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
             'configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

glaciers = ['Upernavik SS',  'Upernavik C', 'Upernavik N',# 'Upernavik NW', 'Upernavik S',
            'Nunatakassaap', 'Kakivfaat', 'Qeqertarsuup', 'Ussing Braeer',
            'Ussing Braeer N', 'Cornell', 'Cornell N']
glaciers = ['Upernavik N']

# read_Theta_Salt_near_plume_location(config_dir, glaciers)

all_glacier_Qsg_timeseries = read_glacier_Qsg(config_dir, glaciers)

depth, dec_yrs, melt_rates_iceplume = read_melt_timeseries_from_nc(config_dir, 'baseline_iceplume')

_, _, melt_rates_melange_iceplume = read_melt_timeseries_from_nc(config_dir, 'baseline_melange_iceplume')

for glacier in glaciers:

    plot_glacier_melt_rates(project_dir, glacier, all_glacier_Qsg_timeseries,
                            dec_yrs, depth,
                            melt_rates_iceplume, melt_rates_melange_iceplume)




