

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shapefile
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
import datetime as dt

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def read_seaice_timeseries(data_folder):
    ds = nc4.Dataset(os.path.join(data_folder,'Sea Ice','Upernavik Sea Ice Concentration Timeseries.nc'))
    dec_yrs = ds.variables['dec_yrs'][:]
    seaice = ds.variables['median_conc'][:]
    ds.close()
    timeseries = np.column_stack([dec_yrs, seaice])
    return(timeseries)

def read_modeled_seaice_timeseries(project_dir, experiment):


    file_name = os.path.join(project_dir,'Data','Models','Sea Ice','SIarea_'+experiment+'_median_timeseries.nc')
    ds = nc4.Dataset(file_name)
    siarea = ds.variables['SIarea'][:]
    dec_yrs = ds.variables['time'][:]
    ds.close()

    seaice_timeseries = np.column_stack([dec_yrs, siarea])

    return(seaice_timeseries)

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

def plot_timeseries(project_dir, obs_seaice_timeseries,
                    experiments, all_model_timeseries):

    years = [2016]
    x_ticks, x_tick_labels, x_grid_locations = compute_ticks_and_labels(years)

    fig = plt.figure(figsize=(8,6))

    gs1 = GridSpec(4, 1, left=0.13, right=0.97, bottom=0.05, top=0.95, hspace=0.05)

    letters = ['a)', 'b)', 'c)', 'd)']

    for e, experiment in enumerate(experiments):

        ax0 = fig.add_subplot(gs1[e, :])

        if experiment!='':
            seaice_timeseries = all_model_timeseries[experiment]
            ax0.plot(seaice_timeseries[:, 0], seaice_timeseries[:, 1], '-', color='b', label=experiment)
        ax0.plot(obs_seaice_timeseries[:, 0], obs_seaice_timeseries[:, 1], '-', color='k', label='Observations')


        seaice_min = -0.1
        seaice_max = 1.25
        ax0.set_ylim([seaice_min,seaice_max])
        ax0.set_xlim([2016, 2017])
        ax0.text(2016.01,0.98*seaice_max,letters[e]+' '+experiment, ha='left', va='top', fontsize=12, color='black')

        if e==0:
            ax0.set_title('Median Sea Ice Concentration in Upernavik')
            ax0.legend()
        if e==1:
            ax0.set_ylabel('Sea Ice Concentration\nm$^2$/m$^2$')
        if e==3:
            ax0.set_xticks(x_ticks)
            ax0.set_xticklabels(x_tick_labels, fontsize=10)
        else:
            ax0.set_xticks([])
            ax0.set_xticklabels([])
        for loc in range(len(x_grid_locations)):
            ax0.plot([x_grid_locations[loc], x_grid_locations[loc]], [seaice_min,seaice_max], '-', linewidth=0.5,
                     color='silver')

        # add a label for the experiment with a white outline
        ax0.text(0.98, 0.02*seaice_max, experiment, ha='right', va='bottom',
                 fontsize=12, color='black')



    output_file = os.path.join(project_dir, 'Figures','Ocean','Sea Ice', 'Upernavik_Sea_Ice_Concentration_Timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

location = 'Upernavik'

data_folder = os.path.join(project_dir, 'Data', 'Observations')

seaice_obs_timeseries = read_seaice_timeseries(data_folder)

experiments = ['control','control_old','','']
all_model_timeseries = {}
for experiment in experiments:
    if experiment != '':
        seaice_timeseries_control = read_modeled_seaice_timeseries(project_dir, experiment)
        all_model_timeseries[experiment] = seaice_timeseries_control


plot_timeseries(project_dir, seaice_obs_timeseries,
                experiments, all_model_timeseries)


