
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pyproj import Transformer

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1,run_test = True):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_grid_components(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir, 'nc_grids',
                                  model_name + '_grid.nc'))
    XC = ds.variables['XC'][:]
    YC = ds.variables['YC'][:]
    Depth = ds.variables['Depth'][:]
    ds.close()

    return (XC, YC, Depth)

def read_chl_fields(project_dir, experiment, year, month, depth):

    if experiment == 'observations':
        chl_file = os.path.join(project_dir, 'Data', 'Observations',
        'Chlorophyll', 'Chl_observations_mean_map_{}{:02d}.nc'.format(year, month))
    else:
        chl_file = os.path.join(project_dir, 'Data', 'Models',
        'Chlorophyll', 'Chl_{}_mean_map_{:04d}{:02d}_{}m.nc'.format(experiment, year, month, depth))

    ds = nc4.Dataset(chl_file)
    chl = ds.variables['Chl'][:]
    ds.close()

    return(chl)

def read_chl_timeseries(project_dir, experiment):

    if experiment == 'observations':
        chl_file = os.path.join(project_dir, 'Data', 'Observations',
        'Chlorophyll', 'Chl_observations_mean_timeseries.nc')
    else:
        chl_file = os.path.join(project_dir, 'Data', 'Models',
        'Chlorophyll', 'Chl_{}_mean_timeseries.nc'.format(experiment))

    ds = nc4.Dataset(chl_file)
    time = ds.variables['time'][:]
    chl_ts = ds.variables['Chl'][:]
    ds.close()

    timeseries = np.column_stack([time, chl_ts])

    timeseries = timeseries[timeseries[:,0]!=0,:]

    return(timeseries)

def plot_chl_comparison(project_dir, year, month, X, Y, Depth,
                        Chl_observations, Chl_control, Chl_control_timeseries):

    if month==8:
        month_name = 'August'

    fig = plt.figure(figsize=(10, 9))

    plot_height = 3
    plot_width = 4
    v_spacing = 1
    colorbar_height = 1
    timeseries_height = 3

    vmin = 0
    vmax = 2

    gs = GridSpec(2*plot_height + timeseries_height + v_spacing, 3*plot_width,
                  left = 0.08, right = 0.95, top = 0.92, bottom = 0.05, hspace = 0.6)

    ax1 = fig.add_subplot(gs[0:plot_height, 0:plot_width])
    plt.pcolormesh(X, Y, Chl_observations, cmap='turbo', vmin=vmin, vmax=vmax)
    ax1.contour(X, Y, Depth, levels=[0.5], colors='k', linewidths=0.5)
    ax1.set_title('Observations')
    ax1.set_xticks([])
    ax1.set_yticks([])

    ###########################################################################

    titles = ['Observations', 'Control','Melange Only','Discharge Only', 'Melange + Discharge']
    for m in range(1,5):
        if m==1:
            ax = fig.add_subplot(gs[0:plot_height, plot_width:2*plot_width])
            chl_grid = Chl_control
        if m==2:
            ax = fig.add_subplot(gs[0:plot_height, 2*plot_width:3*plot_width])
            ax1.set_title('Observations')
        if m==3:
            ax = fig.add_subplot(gs[plot_height:2*plot_height, plot_width:2*plot_width])
        if m==4:
            ax = fig.add_subplot(gs[plot_height:2*plot_height, 2*plot_width:3*plot_width])

        if m<2:
            plt.pcolormesh(X, Y, chl_grid, cmap='turbo', vmin=vmin, vmax=vmax)
        ax.contour(X, Y, Depth, levels=[0.5], colors='k', linewidths=0.5)
        if m!=1:
            ax.set_title(titles[m])
        else:
            ax.set_title('Mean '+month_name+' Chlorophyll-a Concentration\n\nControl')

        ax.set_xticks([])
        ax.set_yticks([])

    ###########################################################################

    axC = fig.add_subplot(gs[plot_height+v_spacing:plot_height+v_spacing+colorbar_height,
                          0:plot_width])

    cx = np.linspace(vmin, vmax, 100)
    cy = np.array([0,1])
    CX, CY = np.meshgrid(cx, cy)
    C = axC.pcolormesh(CX, CY, CX, cmap='turbo', vmin=vmin, vmax=vmax)
    axC.set_yticks([])
    axC.set_title('Chlorophyll-a\nConcentration')
    axC.set_xlabel('mg/m$^3$')

    ###########################################################################

    axT = fig.add_subplot(gs[2*plot_height:2*plot_height+timeseries_height,
                          0:3*plot_width])

    axT.plot(Chl_control_timeseries[:,0], Chl_control_timeseries[:,1], label='Control')

    axT.legend()

    axT.set_xlim([2016,2017])
    axT.set_title('Mean Chlorophyll-a Timeseries in Model Domain')
    axT.set_ylabel('Chlorophyll-a\nConcentration (mg/m$^3$)')
    v_range = vmax - vmin
    axT.set_ylim([vmin-0.1*v_range, vmax+0.1*v_range])

    output_file = os.path.join(project_dir, 'Figures', 'Chlorophyll',
                               'Chl_spatial_comparison_{}{:02d}.png'.format(year, month))
    plt.savefig(output_file, dpi=300)
    plt.close(fig)


config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/MITgcm/' \
              'configurations/downscale_darwin/'

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/'

year = 2016
month = 8
depth = 12

XC, YC, Depth = read_grid_components(config_dir, 'L2_Upernavik')

points = reproject_polygon(np.column_stack([XC.ravel(), YC.ravel()]), inputCRS=4326, outputCRS=3413)
X = points[:,0].reshape(XC.shape)
Y = points[:,1].reshape(YC.shape)

Chl_observations = read_chl_fields(project_dir, 'observations', year, month, depth)
Chl_control = read_chl_fields(project_dir, 'control', year, month, depth)

Chl_control_timeseries = read_chl_timeseries(project_dir, 'control')

print('Observations', np.nanmin(Chl_observations), np.nanmax(Chl_observations))
print('Control', np.min(Chl_control), np.max(Chl_control))

plot_chl_comparison(project_dir, year, month, X, Y, Depth,
                    Chl_observations, Chl_control, Chl_control_timeseries)

