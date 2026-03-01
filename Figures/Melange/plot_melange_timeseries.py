
import os
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from matplotlib.gridspec import GridSpec

def read_melange_timeseries_from_nc(project_dir, experiment):
    file_name = os.path.join(project_dir, 'Data', 'Models', 'Melange',
                               'Melange_' + experiment + '_mean_timeseries.nc')

    ds = nc4.Dataset(file_name)
    dec_yrs = ds.variables['time'][:]
    all_timeseries = {}
    for glacier in list(ds.groups.keys()):
        all_timeseries[glacier] = {}
        for var_name in list(ds.groups[glacier].variables.keys()):
            all_timeseries[glacier][var_name] = ds.groups[glacier].variables[var_name][:]
    ds.close()

    return(dec_yrs, all_timeseries)

def plot_melange_timeseries(project_dir, glacier, all_timeseries_melange, dec_yrs):

    fig = plt.figure(figsize=(8, 8))

    gs = GridSpec(3, 1, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    timeseries_melange = all_timeseries_melange[glacier]['fwflx']
    ax1.plot(dec_yrs, timeseries_melange, label='melange', color='dodgerblue')
    ax1.set_title(glacier+' Melange Melt Timeseries')
    ax1.set_ylabel('Freshwater Flux\n(kg/s, into ocean)')
    ax1.set_xlim([2016,2017])
    plt.legend()

    ax2 = fig.add_subplot(gs[1, 0])
    timeseries_melange = all_timeseries_melange[glacier]['htflx']/1e6
    ax2.plot(dec_yrs, timeseries_melange, label='melange', color='dodgerblue')
    ax2.set_ylabel('Heat Flux\n(MW/m$^2$, out of ocean)')
    ax2.set_xlim([2016, 2017])

    ax3 = fig.add_subplot(gs[2, 0])
    timeseries_melange = all_timeseries_melange[glacier]['mltrt']
    ax3.plot(dec_yrs, timeseries_melange, label='melange', color='dodgerblue')
    ax3.set_ylabel('Melt Rate\n(m/s)')
    ax3.set_xlim([2016, 2017])

    output_file = os.path.join(project_dir, 'Figures','Ocean', 'Melange', glacier+'_melange_timeseries.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

dec_yrs, all_timeseries_melange = read_melange_timeseries_from_nc(project_dir, 'melange')

plot_melange_timeseries(project_dir, 'Upernavik', all_timeseries_melange, dec_yrs)

