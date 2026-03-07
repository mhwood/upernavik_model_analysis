

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
from matplotlib.gridspec import GridSpec
import shutil

def read_grid_geometry_from_nc(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc'))
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    drF = ds.variables['drF'][:]
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, drF, Z, Depth)


def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def iter_number_to_dec_yr(iter_number,seconds_per_iter=60):

    total_seconds = iter_number*seconds_per_iter
    date = datetime.datetime(1992,1,1) + datetime.timedelta(seconds=total_seconds)

    dec_yr = YMD_to_DecYr(date.year,date.month,date.day)
    # print(date)
    return(dec_yr)


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


def read_melt_timeseries_from_nc(config_dir, results_dir, year, glaciers, glacier_to_center_locations):

    if year %4==0:
        n_timesteps = 366
    else:
        n_timesteps = 365

    glacier_melt_grids = {}
    for glacier in glaciers:
        glacier_melt_grids[glacier] = np.zeros((n_timesteps*4, 61))
    dec_yrs = np.zeros((n_timesteps*4,))

    glacier_to_point = {'NW': 0, 'N': 0, 'C': 0, 'S': 0, 'SS': 0}

    first_file = True
    points_counted = 0

    for month in range(1,13):
        file_name = os.path.join(config_dir,'L2','L2_Upernavik',results_dir,'dv','iceplume','ICEFRNTM',
                                 'ICEFRNTM_'+str(year)+'{:02d}'.format(month)+'.nc')

        if month in [1,3,5,7,8,10,12]:
            n_days = 31
        elif month in [4,6,9,11]:
            n_days = 30
        else:
            if year%4==0:
                n_days = 29
            else:
                n_days = 28

        counter = 0
        for day in range(1,n_days+1):
            for hour in [0, 6, 12, 18]:
                dec_yrs[points_counted+counter] = YMD_to_DecYr(year, month, day, hour)
                counter += 1

        if os.path.exists(file_name):
            ds = nc4.Dataset(file_name)
            melt = ds.variables['ICEFRNTM'][:, :, :]
            iterations = ds.variables['iterations'][:]
            if first_file:
                longitude = ds.variables['longitude'][:]
                latitude = ds.variables['latitude'][:]
                melt_sum = np.sum(melt[0, :, :], axis=0)
                point_indices = np.arange(len(longitude))
                longitude = longitude[melt_sum != 0]
                latitude = latitude[melt_sum != 0]
                point_indices = point_indices[melt_sum != 0]
                for glacier in glaciers:
                    location = glacier_to_center_locations[glacier]
                    distances = great_circle_distance(location[0], location[1], longitude, latitude)
                    glacier_to_point[glacier] = point_indices[np.argmin(distances)]
                    print(glacier, glacier_to_point[glacier], np.min(distances))
                first_file = False
            ds.close()

            for glacier in glaciers:
                point = glacier_to_point[glacier]
                glacier_melt_grids[glacier][points_counted:points_counted+4*n_days, :] = melt[:, :, point]

        points_counted += n_days*4

    return(dec_yrs, glacier_melt_grids)


def plot_melt_timeseries(project_dir, dec_yrs, Z, glacier_melt_grids):

    min_year = 2016
    max_year = 2017

    vmin = 0
    vmax = 10

    dmin = -1
    dmax = 1

    max_depth = 850
    # melange_draft = 240

    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    x_ticks = []
    x_grid_locations = []
    x_tick_labels = []
    for year in range(min_year,max_year):
        for month in range(1,13):
            x_ticks.append(YMD_to_DecYr(year,month,15))
            x_tick_labels.append(month_labels[month-1])
            if month in [1,3,5,7,8,10,12]:
                x_grid_locations.append(YMD_to_DecYr(year,month,31))
            elif month in [4,6,9,11]:
                x_grid_locations.append(YMD_to_DecYr(year, month, 30))
            else:
                if year%4==0:
                    x_grid_locations.append(YMD_to_DecYr(year, month, 29))
                else:
                    x_grid_locations.append(YMD_to_DecYr(year, month, 28))



    fig = plt.figure(figsize=(8,10), dpi=300)

    gs = GridSpec(6, 23, top=0.96, bottom=0.05, left = 0.08, right=0.91)

    ax1 = fig.add_subplot(gs[0,:-2])
    melt_plume_plot = glacier_melt_grids['NW'].T#np.ma.masked_where(melt_plume == 0, melt_plume)
    ax1.pcolormesh(dec_yrs, Z, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax1.set_xlim([min_year, max_year])
    ax1.set_ylim([max_depth,0])
    ax1.set_ylabel('Depth (m)')
    ax1.grid(linestyle='--',linewidth=0.5,alpha=0.5)
    ax1.set_title('Glacier Submarine Melt Rate')
    # ax1.set_xticklabels([])
    ax1.text(2016.1, 30, 'a)', ha='left', va='top', fontsize=12, color='white')
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_tick_labels)
    for loc in range(len(x_grid_locations)):
        ax1.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth,0], '-', linewidth=0.5, color='silver')

    ax1c = fig.add_subplot(gs[1:4,-1])
    # ax1c.tick_params(axis='both', which='major', labelsize=12)
    cx = np.array([0, 1])
    cy = np.arange(vmin, vmax+0.01,0.01)
    CX, CY = np.meshgrid(cx, cy)
    ax1c.pcolormesh(CX, CY, CY, cmap='turbo')
    ax1c.set_xticks([])
    ax1c.set_ylabel('Melt Rate (m/day)')
    ax1c.yaxis.tick_right()
    ax1c.yaxis.set_label_position("right")

    ax2 = fig.add_subplot(gs[1, :-2])
    melt_plume_plot = glacier_melt_grids['N'].T  # np.ma.masked_where(melt_plume == 0, melt_plume)
    ax2.pcolormesh(dec_yrs, Z, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax2.set_xlim([min_year, max_year])
    ax2.set_ylim([max_depth, 0])
    ax2.set_ylabel('Depth (m)')
    ax2.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    # ax1.set_xticklabels([])
    ax2.text(2016.1, 30, 'a)', ha='left', va='top', fontsize=12, color='white')
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_tick_labels)
    for loc in range(len(x_grid_locations)):
        ax2.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax3 = fig.add_subplot(gs[2, :-2])
    melt_plume_plot = glacier_melt_grids['C'].T  # np.ma.masked_where(melt_plume == 0, melt_plume)
    ax3.pcolormesh(dec_yrs, Z, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax3.set_xlim([min_year, max_year])
    ax3.set_ylim([max_depth, 0])
    ax3.set_ylabel('Depth (m)')
    ax3.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    # ax1.set_xticklabels([])
    ax3.text(2016.1, 30, 'a)', ha='left', va='top', fontsize=12, color='white')
    ax3.set_xticks(x_ticks)
    ax3.set_xticklabels(x_tick_labels)
    for loc in range(len(x_grid_locations)):
        ax3.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax4 = fig.add_subplot(gs[3, :-2])
    melt_plume_plot = glacier_melt_grids['S'].T  # np.ma.masked_where(melt_plume == 0, melt_plume)
    ax4.pcolormesh(dec_yrs, Z, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax4.set_xlim([min_year, max_year])
    ax4.set_ylim([max_depth, 0])
    ax4.set_ylabel('Depth (m)')
    ax4.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    # ax1.set_xticklabels([])
    ax4.text(2016.1, 30, 'a)', ha='left', va='top', fontsize=12, color='white')
    ax4.set_xticks(x_ticks)
    ax4.set_xticklabels(x_tick_labels)
    for loc in range(len(x_grid_locations)):
        ax4.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax5 = fig.add_subplot(gs[4, :-2])
    melt_plume_plot = glacier_melt_grids['SS'].T  # np.ma.masked_where(melt_plume == 0, melt_plume)
    ax5.pcolormesh(dec_yrs, Z, melt_plume_plot, vmin=vmin, vmax=vmax, cmap='turbo')
    ax5.set_xlim([min_year, max_year])
    ax5.set_ylim([max_depth, 0])
    ax5.set_ylabel('Depth (m)')
    ax5.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    # ax1.set_xticklabels([])
    ax5.text(2016.1, 30, 'a)', ha='left', va='top', fontsize=12, color='white')
    ax5.set_xticks(x_ticks)
    ax5.set_xticklabels(x_tick_labels)
    for loc in range(len(x_grid_locations)):
        ax5.plot([x_grid_locations[loc], x_grid_locations[loc]], [max_depth, 0], '-', linewidth=0.5, color='silver')

    ax5 = fig.add_subplot(gs[5, :-2])
    plt.plot(dec_yrs, np.max(glacier_melt_grids['NW'][:, :],axis=1), label='NW', color='red')
    plt.plot(dec_yrs, np.max(glacier_melt_grids['N'][:, :], axis=1), label='N', color='orange')
    plt.plot(dec_yrs, np.max(glacier_melt_grids['C'][:, :], axis=1), label='C', color='green')
    plt.plot(dec_yrs, np.max(glacier_melt_grids['S'][:, :], axis=1), label='S', color='blue')
    plt.plot(dec_yrs, np.max(glacier_melt_grids['SS'][:, :], axis=1), label='SS', color='purple')
    plt.legend()
    ax5.set_ylim([5, 8.5])
    plt.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax5.set_xticks(x_ticks)
    ax5.set_xticklabels(x_tick_labels)

    output_file = os.path.join(project_dir, 'Upernavik Modeled Ice Front Melt.jpg')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_melt_rate_comparison():
    config_dir = '/Volumes/kullorsuaq/Research/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
                 'configurations/downscale_greenland'

    project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

    XC, YC, drF, Z, Depth = read_grid_geometry_from_nc(config_dir, model_name='L2_Upernavik')

    glaciers = ['NW','N','C','S','SS']
    glacier_to_center_locations = {'NW': [-54.52218, 73.01875],
                            'N': [-54.37945, 73.01643,],
                            'C': [-54.38338,72.93889],
                            'S': [-54.39157,72.84235],
                            'SS': [-54.27122,72.79309]}

    dec_yrs, glacier_melt_grids =\
        read_melt_timeseries_from_nc(config_dir, 'results_iceplume_iceberg', 2016, glaciers, glacier_to_center_locations)

    print(dec_yrs)

    plot_melt_timeseries(project_dir, dec_yrs, Z, glacier_melt_grids)


if __name__ == '__main__':
    plot_melt_rate_comparison()
