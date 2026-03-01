
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
from pyproj import Transformer

def reproject_polygon(polygon_array, inputCRS, outputCRS, x_column=0, y_column=1, run_test=True):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(polygon_array)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def read_model_domain(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                  model_name + '_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    ds.close()

    return (XC, YC, Depth)

def read_polygons_from_nc(project_dir):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Models','Fjord Transects',
                                  'Upernavik N Fjord Geometry Transect.nc'))

    distance = ds.variables['distance'][:]
    surface = ds.variables['surface'][:]
    surface*=-1
    distance *=1e-3

    transect_x = ds.variables['transect_x'][:]
    transect_y = ds.variables['transect_y'][:]

    ice_polygon_x = ds.variables['ice_polygon_x'][:]
    ice_polygon_y = ds.variables['ice_polygon_y'][:]
    ice_polygon = np.column_stack([ice_polygon_x, ice_polygon_y])
    ice_polygon[:,0]*=1e-3
    ice_polygon[:,1]*=-1

    bed_polygon_x = ds.variables['bathy_polygon_x'][:]
    bed_polygon_y = ds.variables['bathy_polygon_y'][:]
    bed_polygon = np.column_stack([bed_polygon_x, bed_polygon_y])
    bed_polygon[:,0]*=1e-3
    bed_polygon[:,1]*=-1

    return(distance, transect_x, transect_y, surface, ice_polygon, bed_polygon)

def read_model_transect_data(project_dir, experiment, years, month, var_name):

    year = 2020

    annual_model_grids = {}

    file_name = os.path.join(project_dir, 'Data','Models','Fjord Transects',
                               f'Upernavik_N_Fjord_{experiment}_{var_name.upper()}_Transect_{year}{month:02d}.nc')
    ds = nc4.Dataset(file_name)
    depth = ds.variables['depth'][:]
    distance = ds.variables['distance'][:]
    grid = ds.variables[var_name.upper()][:, :]
    ds.close()

    annual_model_grids[year] = grid

    return(depth, distance, annual_model_grids)

def plot_transects(project_dir, years, tracer_name,
                   distance, surface, ice_polygon, bed_polygon,
                   model_depth_control, model_distance_control, annual_model_grids_control,
                   model_depth_melange, model_distance_melange, annual_model_grids_melange,
                   model_depth_iceplume, model_distance_iceplume, annual_model_grids_iceplume,
                   model_depth_melange_iceplume, model_distance_melange_iceplume, annual_model_grids_melange_iceplume):

    fig = plt.figure(figsize=(9, 10))

    plot_rows = 6
    plot_cols = 12
    cb_cols = 1
    spacing = 1

    gs = GridSpec(4*plot_rows+spacing, plot_cols + cb_cols + spacing,
                  top=0.97, bottom=0.07, left=0.12, right=0.90, hspace=0.6, wspace=0.1)

    ###########################################################
    # Metadata

    if tracer_name == 'PTRACE02':
        wvel_min = 0
        wvel_max = 12
        cmap = 'RdPu'
        dmin = -3
        dmax = 3

    letters = ['a','b','c','d']

    ###########################################################
    # Plot wvel directly from control model

    ax2 = fig.add_subplot(gs[:plot_rows, :-cb_cols-spacing])
    # po4 = po4_dict[years[y]]
    bed = Polygon(bed_polygon, facecolor='grey', edgecolor='k')
    ice = Polygon(ice_polygon, facecolor='white', edgecolor='k')

    ax2.pcolormesh(distance, model_depth_control, annual_model_grids_control[years[0]],
                   cmap=cmap,vmin=wvel_min,vmax=wvel_max)
    ax2.plot(distance, surface, 'k-', linewidth=0.5)
    ax2.add_patch(ice)
    ax2.add_patch(bed)
    ax2.set_ylim([np.max(bed_polygon[:,1]),-50])
    ax2.set_xlim([np.min(distance), np.max(distance)])

    ax2.set_title('Control Model Vertical Velocity')

    ax2.set_xticklabels([])
    # ax2.set_yticklabels([])
    ax2.text(0.05 * np.max(distance), np.max(ice_polygon[:, 1]), r"$\bf{"+letters[0]+"}$", color='white',
             fontsize=14, zorder=99, ha='left', va='bottom')
    ax2.text(0.98, 0.02, letters[0] + ')', transform=ax2.transAxes,
             verticalalignment='bottom', horizontalalignment='right',
             fontsize=12, color='k')

    ###########################################################
    # Plot wvel colorbar

    ax1c = fig.add_subplot(gs[:plot_rows, -cb_cols:])
    y = np.linspace(wvel_min, wvel_max, 100)
    x = np.array([0, 1])
    X, Y = np.meshgrid(x, y)
    ax1c.pcolormesh(X, Y, Y, vmin=wvel_min, vmax=wvel_max, cmap='RdPu')
    ax1c.set_ylabel('Vertical Velocity (mm/s)')
    ax1c.set_xticks([])
    ax1c.yaxis.tick_right()
    ax1c.yaxis.set_label_position("right")

    ###########################################################
    # Plot wvel directly from control model

    for y in range(1,4):
        ax3 = fig.add_subplot(gs[y * plot_rows+spacing:(y + 1) * plot_rows+spacing, :-cb_cols-spacing])
        bed = Polygon(bed_polygon, facecolor='grey', edgecolor='k')
        ice = Polygon(ice_polygon, facecolor='white', edgecolor='k')

        if y == 2:
            annual_model_grids = annual_model_grids_iceplume
            model_depth = model_depth_iceplume
            grid = annual_model_grids_iceplume[years[0]]
            diff_grid = grid - annual_model_grids_control[years[0]]
            diff_grid[diff_grid<dmin] = dmin
            diff_grid[diff_grid>dmax] = dmax

            ax3.contourf(distance, model_depth, diff_grid,
                         levels=np.linspace(dmin, dmax, 100),
                           cmap='seismic')#, vmin=wvel_min, vmax=wvel_max)

            print("Plotting iceplume model")
            print(np.nanmin(annual_model_grids[years[0]]), np.nanmax(annual_model_grids[years[0]]))
        elif y == 1:
            annual_model_grids = annual_model_grids_melange
            model_depth = model_depth_melange
            grid = annual_model_grids[years[0]]
            diff_grid = grid - annual_model_grids_control[years[0]]
            diff_grid[diff_grid < dmin] = dmin
            diff_grid[diff_grid > dmax] = dmax

            ax3.contourf(distance, model_depth, diff_grid,
                             levels=np.linspace(dmin, dmax, 100),
                               cmap='seismic')
        elif y == 3:
            annual_model_grids = annual_model_grids_melange_iceplume
            model_depth = model_depth_melange_iceplume
            grid = annual_model_grids[years[0]]
            diff_grid = grid - annual_model_grids_control[years[0]]
            diff_grid[diff_grid < dmin] = dmin
            diff_grid[diff_grid > dmax] = dmax

            ax3.contourf(distance, model_depth, diff_grid,
                             levels=np.linspace(dmin, dmax, 100),
                               cmap='seismic')
        ax3.plot(distance, surface, 'k-', linewidth=0.5)
        ax3.add_patch(ice)
        ax3.add_patch(bed)
        ax3.set_ylim([np.max(bed_polygon[:, 1]), -50])
        ax3.set_xlim([np.min(distance), np.max(distance)])
        if y == 1:
            ax3.set_title('Anomalies Relative to Control')
            ax3.set_ylabel('Depth (m)')
        if y == 3:
            ax3.set_xlabel('Distance Along Fjord (km)')
        else:
            ax3.set_xticklabels([])
        # ax3.set_yticklabels([])
        ax3.text(0.98, 0.02, letters[y] + ')', transform=ax3.transAxes,
                 verticalalignment='bottom', horizontalalignment='right',
                 fontsize=12, color='k')

    ax2c = fig.add_subplot(gs[int(1.5*plot_rows)+spacing:int(3.5*plot_rows)+spacing, -cb_cols:])
    y = np.linspace(dmin, dmax, 100)
    x = np.array([0, 1])
    X, Y = np.meshgrid(x, y)
    ax2c.pcolormesh(X, Y, Y, vmin=dmin, vmax=dmax, cmap='seismic')
    ax2c.set_ylabel('Anomaly ()')
    ax2c.set_xticks([])
    ax2c.yaxis.tick_right()
    ax2c.yaxis.set_label_position("right")


    output_file = os.path.join(project_dir, 'Figures', 'Ocean', 'Upernavik N '+tracer_name+' Transect Comparison.png')
    plt.savefig(output_file)
    plt.close(fig)


project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'MITgcm/configurations/downscale_darwin'

XC, YC, Depth = read_model_domain(config_dir, 'L2_Upernavik')

tracer_name = 'PTRACE02'

years = [2020]

distance, transect_x, transect_y, surface, ice_polygon, bed_polygon = read_polygons_from_nc(project_dir)


model_depth_control, model_distance_control, annual_model_grids_control = \
    read_model_transect_data(project_dir, 'baseline', 2020, 11, tracer_name)

model_depth_melange, model_distance_melange, annual_model_grids_melange = \
    read_model_transect_data(project_dir, 'baseline_melange', 2020, 11, tracer_name)

model_depth_iceplume, model_distance_iceplume, annual_model_grids_iceplume = \
    read_model_transect_data(project_dir, 'baseline_iceplume', 2020, 11, tracer_name)

model_depth_melange_iceplume, model_distance_melange_iceplume, annual_model_grids_melange_iceplume = \
    read_model_transect_data(project_dir, 'baseline_melange_iceplume', 2020, 11, tracer_name)


# plt.imshow(annual_model_grids_iceplume[2016])
# plt.show()

plot_transects(project_dir, years, tracer_name,
               distance, surface, ice_polygon, bed_polygon,
               model_depth_control, model_distance_control, annual_model_grids_control,
               model_depth_melange, model_distance_melange, annual_model_grids_melange,
               model_depth_iceplume, model_distance_iceplume, annual_model_grids_iceplume,
               model_depth_melange_iceplume, model_distance_melange_iceplume, annual_model_grids_melange_iceplume)


