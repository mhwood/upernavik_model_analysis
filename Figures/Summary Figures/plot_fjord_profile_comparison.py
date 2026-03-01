
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
from pyproj import Transformer


def read_model_transect_mean_profile(project_dir, fjord_name, experiment, year, month, var_name):

    file_name = os.path.join(project_dir, 'Data','Models','Fjord Transects',fjord_name,
                               '_'.join(fjord_name.split())+f'_Fjord_{experiment}_{var_name.upper()}_Transect_{year}{month:02d}.nc')
    ds = nc4.Dataset(file_name)
    depth = ds.variables['depth'][:]
    grid = ds.variables[var_name.upper()][:, :]
    ds.close()

    grid[grid==0] = np.nan
    profile = np.nanmean(grid, axis=1)
    profile = np.column_stack([depth, profile])

    return(profile)

def plot_profile_comparisons(project_dir, fjord_names, tracer_names, profile_dict):

    print(profile_dict.keys())

    fig = plt.figure(figsize=(12, 8))

    gs = GridSpec(len(fjord_names), 5, figure=fig, wspace=0.3, hspace=0.3)

    for i, fjord_name in enumerate(fjord_names):

        fjord_dict = profile_dict[fjord_name]

        for j, tracer_name in enumerate(tracer_names):
            ax = fig.add_subplot(gs[i, j])
            tracer_profiles = fjord_dict[tracer_name]

            for experiment in tracer_profiles.keys():
                if experiment!='control':
                    print(experiment)
                    ax.plot(tracer_profiles[experiment][:, 1]-tracer_profiles['control'][:, 1],
                            tracer_profiles[experiment][:, 0], label=experiment)

            if i==0:
                ax.set_title(f'{tracer_name}')
            if i==2:
                ax.set_xlabel('Tracer Concentration')
            if j==0:
                ax.set_ylabel(fjord_name+'\nDepth (m)')
            else:
                ax.set_yticklabels([])
            if i==0 and j==0:
                ax.legend()
            ax.set_ylim([800,0])
            # ax.invert_yaxis()
            # ax.legend()

    output_file = os.path.join(project_dir, 'Figures', 'Ocean', 'Fjord_Profile_Comparison.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    a=1

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/'

tracer_names = ['THETA','SALT','WVEL','PTRACE02','PTRACE28']

fjord_names = ['Ussing Braeer','Kakivfaat', 'Upernavik N']

profile_dict = {}

for fjord_name in fjord_names:

    fjord_tracer_dict = {}

    for tracer_name in tracer_names:

        tracer_dict = {}

        model_profile_control = read_model_transect_mean_profile(project_dir, fjord_name,
                                                                 'baseline', 2020, 11, tracer_name)

        model_profile_melange = read_model_transect_mean_profile(project_dir, fjord_name,
                                                                 'baseline_melange', 2020, 11, tracer_name)

        model_profile_iceplume = read_model_transect_mean_profile(project_dir, fjord_name,
                                                                  'baseline_iceplume', 2020, 11, tracer_name)

        model_profile_melange_iceplume = read_model_transect_mean_profile(project_dir,fjord_name,
                                                                          'baseline_melange_iceplume', 2020, 11, tracer_name)

        tracer_dict['control'] = model_profile_control
        tracer_dict['melange'] = model_profile_melange
        tracer_dict['iceplume'] = model_profile_iceplume
        tracer_dict['melange_iceplume'] = model_profile_melange_iceplume

        fjord_tracer_dict[tracer_name] = tracer_dict

    profile_dict[fjord_name] = fjord_tracer_dict


plot_profile_comparisons(project_dir, fjord_names, tracer_names, profile_dict)

