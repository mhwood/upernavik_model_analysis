


import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)


def read_iceplume_grid(config_dir):
    grid = np.fromfile(os.path.join(config_dir, 'L2', 'L2_Upernavik', 'input','iceplume', 'L2_Qsg_2021'), '>f4')
    grid = grid.reshape((365,375,450))
    return(grid)

def sample_glacier_timeseries(Qsg):

    glacier_timeseries = {}

    glaciers = ['Upernavik Isstrom SS', 'Upernavik Isstrom S',
                'Upernavik Isstrom C', 'Upernavik Isstrom N', 'Upernavik Isstrom NW',
                'Nunatakassaap', 'Kakivfaat Sermia', 'Qeqertarsuup Sermia',
                'Ussing Braeer', 'Ussing Braeer N', 'Cornell']

    nonzero_grid = np.sum(Qsg, axis=0)
    rows, cols = np.where(nonzero_grid > 0)

    for g in range(len(glaciers)):
        glacier_timeseries[glaciers[g]] = Qsg[:,rows[g],cols[g]]

    return(glacier_timeseries)
def compute_ticks_and_labels(years):

    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

    x_ticks = []
    x_grid_locations = []

    for year in years:
        for month in range(1, 13):
            x_ticks.append(YMD_to_DecYr(year, month, 15))
            if month in [1, 3, 5, 7, 8, 10, 12]:
                x_grid_locations.append(YMD_to_DecYr(year, month, 31))
            elif month in [4, 6, 9, 11]:
                x_grid_locations.append(YMD_to_DecYr(year, month, 30))
            else:
                if year % 4 == 0:
                    x_grid_locations.append(YMD_to_DecYr(year, month, 29))
                else:
                    x_grid_locations.append(YMD_to_DecYr(year, month, 28))
        x_ticks = np.array(x_ticks)
    return(x_ticks, month_labels, x_grid_locations)


def plot_subglacial_discharge_timeseries(project_dir, glacier_timeseries):

    x_ticks, x_tick_labels, x_grid_locations = compute_ticks_and_labels([2021])
    dec_yrs = np.linspace(2021, 2022, 365)

    Qsg_min = -50
    Qsg_max = 1000

    fig = plt.figure(figsize=(8,8))

    plt.subplot(3,1,1)
    ax0 = plt.gca()
    ax0.plot(dec_yrs, glacier_timeseries['Ussing Braeer N'], color='green', label='Ussing Braeer N')
    ax0.plot(dec_yrs, glacier_timeseries['Ussing Braeer'], color='blue', label='Ussing Braeer')
    ax0.set_xlim([2021, 2022])
    ax0.set_ylim([Qsg_min, Qsg_max])
    ax0.set_xticks(x_ticks)
    ax0.set_xticklabels(x_tick_labels, fontsize=10)
    for loc in range(len(x_grid_locations)):
        ax0.plot([x_grid_locations[loc], x_grid_locations[loc]], [Qsg_min, Qsg_max], '-', linewidth=0.5, color='silver')
    ax0.set_title('Subglacial Discharge')
    ax0.legend(loc=2)

    plt.subplot(3, 1, 2)
    ax0 = plt.gca()
    ax0.plot(dec_yrs, glacier_timeseries['Qeqertarsuup Sermia'], color='green', label='Qeqertarsuup Sermia')
    ax0.plot(dec_yrs, glacier_timeseries['Kakivfaat Sermia'], color='blue', label='Kakivfaat Sermia')
    ax0.set_xlim([2021, 2022])
    ax0.set_ylim([Qsg_min, Qsg_max])
    ax0.set_xticks(x_ticks)
    ax0.set_xticklabels(x_tick_labels, fontsize=10)
    for loc in range(len(x_grid_locations)):
        ax0.plot([x_grid_locations[loc], x_grid_locations[loc]], [Qsg_min, Qsg_max], '-', linewidth=0.5, color='silver')
    plt.ylabel('(m$^3$/s)')
    ax0.legend(loc=2)

    plt.subplot(3, 1, 3)
    ax0 = plt.gca()
    ax0.plot(dec_yrs, glacier_timeseries['Upernavik Isstrom NW'], color='red',  label='Upernavik Isstrom NW')
    ax0.plot(dec_yrs, glacier_timeseries['Upernavik Isstrom N'], color='orange', label='Upernavik Isstrom N')
    ax0.plot(dec_yrs, glacier_timeseries['Upernavik Isstrom C'], color='gold', label='Upernavik Isstrom C')
    ax0.plot(dec_yrs, glacier_timeseries['Upernavik Isstrom S'], color='green', label='Upernavik Isstrom S')
    ax0.plot(dec_yrs, glacier_timeseries['Upernavik Isstrom SS'], color='blue', label='Upernavik Isstrom SS')
    ax0.set_xlim([2021, 2022])
    ax0.set_ylim([Qsg_min, Qsg_max])
    ax0.set_xticks(x_ticks)
    ax0.set_xticklabels(x_tick_labels, fontsize=10)
    for loc in range(len(x_grid_locations)):
        ax0.plot([x_grid_locations[loc], x_grid_locations[loc]], [Qsg_min, Qsg_max], '-', linewidth=0.5, color='silver')
    ax0.legend(loc=2)


    plt.savefig(os.path.join(project_dir, 'Figures', 'Glacier', 'Upernavik_subglacial_discharge_timeseries.png'), dpi=300)
    plt.close(fig)


config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'


Qsg = read_iceplume_grid(config_dir)

glacier_timeseries = sample_glacier_timeseries(Qsg)

plot_subglacial_discharge_timeseries(project_dir, glacier_timeseries)

