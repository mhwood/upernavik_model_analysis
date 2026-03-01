
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
from pyproj import Transformer

def read_polygons_from_nc(project_dir):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data','Ocean','Fjord Transects',
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

def read_AXCTD_transect_data(project_dir):

    years = [2016, 2018, 2019, 2020]

    file_name = os.path.join(project_dir, 'Data', 'Ocean', 'AXCTDs', 'Upernavik_N_Fjord_AXCTD_Transects.nc')

    ds = nc4.Dataset(file_name)

    depth = ds.variables['depth'][:]
    distance = ds.variables['distance'][:]

    annual_grids = {}
    for year in years:
        drp = ds[str(year)]
        theta_var = drp.variables['potential_temperature'][:, :]
        salt_var = drp.variables['practical_salinity'][:, :]
        ctd_depths = drp.variables['ctd_max_depths'][:]
        ctd_distances = drp.variables['ctd_distances'][:]
        annual_grids[year] = {'Theta': theta_var[:, :], 'Salt': salt_var[:, :],
                              'ctd_depths': ctd_depths[:], 'ctd_distances': ctd_distances[:]}

    ds.close()
    return(depth, distance, annual_grids)

def read_model_transect_data(project_dir, years, month, var_name):

    year = 2016

    annual_model_grids = {}

    file_name = os.path.join(project_dir, 'Data','Ocean','Fjord Transects',
                               f'Upernavik_N_Fjord_{var_name.upper()}_Transect_{year}{month:02d}.nc')
    ds = nc4.Dataset(file_name)
    depth = ds.variables['depth'][:]
    distance = ds.variables['distance'][:]
    grid = ds.variables[var_name.upper()][:, :]
    ds.close()

    annual_model_grids[year] = grid

    return(depth, distance, annual_model_grids)

def plot_transects(project_dir, years, var_name,
                   axctd_depth, model_dpeth, distance, surface, ice_polygon, bed_polygon,
                   annual_AXCTD_grids, annual_model_grids):

    fig = plt.figure(figsize=(8, 10))

    plot_rows = 6
    cb_rows = 1

    gs = GridSpec(4*plot_rows+cb_rows+2, 2,
                  top=0.97, bottom=0.07, left=0.12, right=0.95, hspace=0.6, wspace=0.1)

    ###########################################################
    # AXCTDs

    theta_min = -1
    theta_max = 2

    letters = ['a','c','e','g']

    for y in range(len(years)):
        ax1 = fig.add_subplot(gs[y*plot_rows:(y+1)*plot_rows, 0])
        grid = annual_AXCTD_grids[years[y]][var_name]
        bed = Polygon(bed_polygon, facecolor='grey', edgecolor='k')
        ice = Polygon(ice_polygon, facecolor='white', edgecolor='k')

        ax1.pcolormesh(distance,axctd_depth,grid,cmap='turbo',vmin=theta_min,vmax=theta_max)
        ax1.plot(distance, surface, 'k-', linewidth=0.5)
        ax1.add_patch(ice)
        ax1.add_patch(bed)
        ax1.set_ylim([np.max(bed_polygon[:,1]),-50])
        ax1.set_xlim([np.min(distance), np.max(distance)])
        if y==0:
            ax1.set_title('Observations')
        if y==3:
            ax1.set_xlabel('Distance Along Transect (km)')
        else:
            ax1.set_xticklabels([])
        if y==1:
            ax1.set_ylabel('Depth (m)')
        ax1.text(0.98, 0.02, letters[y]+') '+str(years[y]), transform=ax1.transAxes,
                 verticalalignment='bottom', horizontalalignment='right',
                 fontsize=12, color='k')

    ###########################################################
    # Model

    po4_min = -3
    po4_max = 3

    letters = ['b','d','f','h']

    for y in range(1):
        ax2 = fig.add_subplot(gs[y*plot_rows:(y+1)*plot_rows, 1])
        # po4 = po4_dict[years[y]]
        bed = Polygon(bed_polygon, facecolor='grey', edgecolor='k')
        ice = Polygon(ice_polygon, facecolor='white', edgecolor='k')

        ax2.pcolormesh(distance, model_depth, annual_model_grids[years[y]],
                       cmap='turbo',vmin=theta_min,vmax=theta_max)
        ax2.plot(distance, surface, 'k-', linewidth=0.5)
        ax2.add_patch(ice)
        ax2.add_patch(bed)
        ax2.set_ylim([np.max(bed_polygon[:,1]),-50])
        ax2.set_xlim([np.min(distance), np.max(distance)])
        if y == 0:
            ax2.set_title('Model')
        if y==3:
            ax2.set_xlabel('Distance Along Transect (km)')
        else:
            ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.text(0.05 * np.max(distance), np.max(ice_polygon[:, 1]), r"$\bf{"+letters[y]+"}$", color='white',
                 fontsize=14, zorder=99, ha='left', va='bottom')
        ax2.text(0.98, 0.02, letters[y] + ')', transform=ax2.transAxes,
                 verticalalignment='bottom', horizontalalignment='right',
                 fontsize=12, color='k')

    ax1c = fig.add_subplot(gs[-1,:])
    x = np.linspace(-1, 2, 100)
    y=np.array([0,1])
    X,Y=np.meshgrid(x,y)
    ax1c.pcolormesh(X,Y,X,vmin=theta_min,vmax=theta_max,cmap='turbo')
    ax1c.set_xlabel('Temperature (°C)')
    ax1c.set_yticks([])

    # ax2c = fig.add_subplot(gs[-1, 1])
    # x = np.linspace(po4_min, po4_max, 100)
    # y = np.array([0, 1])
    # X, Y = np.meshgrid(x, y)
    # ax2c.pcolormesh(X, Y, X, vmin=po4_min, vmax=po4_max, cmap='seismic')
    # ax2c.set_xlabel('$\Delta$[PO$_4$] ($\mu M$)')
    # ax2c.set_yticks([])


    output_file = os.path.join(project_dir, 'Figures', 'Ocean', 'Upernavik N '+var_name+' Transect.png')
    plt.savefig(output_file)
    plt.close(fig)




project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

# config_dir = '/Volumes/eqipsermia/downscale_darwin/L2_Disko_Bay'
config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'MITgcm/configurations/downscale_darwin/L2/L2_Upernavik'

distance, transect_x, transect_y, surface, ice_polygon, bed_polygon = read_polygons_from_nc(project_dir)

years = [2016,2018,2019,2020]

axctd_depth, _, annual_AXCTD_grids = read_AXCTD_transect_data(project_dir)

model_depth, _, annual_model_grids = read_model_transect_data(project_dir, 2016, 8, 'Theta')

plot_transects(project_dir, years, 'Theta',
               axctd_depth, model_depth, distance, surface, ice_polygon, bed_polygon,
               annual_AXCTD_grids, annual_model_grids)


