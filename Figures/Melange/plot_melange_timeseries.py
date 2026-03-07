
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

def plot_melange_timeseries(project_dir, glacier, dec_yrs_melange, all_timeseries_melange,
                        dec_yrs_melange_iceplume, all_timeseries_melange_iceplume):

    okabe_ito_colors = ["#E69F00","#56B4E9","#009E73","#F5C710","#0072B2","#D55E00","#CC79A7","#999999","#000000"]
    colors = [okabe_ito_colors[0], okabe_ito_colors[5], okabe_ito_colors[4], okabe_ito_colors[2]]

    fig = plt.figure(figsize=(8, 8))

    gs = GridSpec(3, 1, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    timeseries_melange = all_timeseries_melange[glacier]['fwflx']
    ax1.plot(dec_yrs_melange, timeseries_melange, label='melange', color=colors[1])
    timeseries_melange_iceplume = all_timeseries_melange_iceplume[glacier]['fwflx']
    ax1.plot(dec_yrs_melange_iceplume, timeseries_melange_iceplume, label='melange + ice plume', color=colors[3])
    ax1.set_title(glacier+' Melange Melt Timeseries')
    ax1.set_ylabel('Freshwater Flux\n(kg/s, into ocean)')
    ax1.set_xlim([2020.8,2022])
    plt.legend()

    ax2 = fig.add_subplot(gs[1, 0])
    timeseries_melange = all_timeseries_melange[glacier]['htflx']/1e6
    ax2.plot(dec_yrs_melange, timeseries_melange, label='melange', color=colors[1])
    timeseries_melange_iceplume = all_timeseries_melange_iceplume[glacier]['htflx']/1e6
    ax2.plot(dec_yrs_melange_iceplume, timeseries_melange_iceplume, label='melange + ice plume', color=colors[3])
    ax2.set_ylabel('Heat Flux\n(MW/m$^2$, out of ocean)')
    ax2.set_xlim([2020.8,2022])

    ax3 = fig.add_subplot(gs[2, 0])
    timeseries_melange = all_timeseries_melange[glacier]['mltrt']
    ax3.plot(dec_yrs_melange, timeseries_melange, label='melange', color=colors[1])
    timeseries_melange_iceplume = all_timeseries_melange_iceplume[glacier]['mltrt']
    ax3.plot(dec_yrs_melange_iceplume, timeseries_melange_iceplume, label='melange + ice plume', color=colors[3])
    ax3.set_ylabel('Melt Rate\n(m/s)')
    ax3.set_xlim([2020.8,2022])

    output_file = os.path.join(project_dir, 'Figures','Ocean', 'Melange', glacier+'_melange_timeseries.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

dec_yrs_melange , all_timeseries_melange = read_melange_timeseries_from_nc(project_dir, 'baseline_melange')

dec_yrs_melange_iceplume, all_timeseries_melange_iceplume = read_melange_timeseries_from_nc(project_dir, 'baseline_melange_iceplume')

plot_melange_timeseries(project_dir, 'Upernavik',
                        dec_yrs_melange, all_timeseries_melange,
                        dec_yrs_melange_iceplume, all_timeseries_melange_iceplume)

