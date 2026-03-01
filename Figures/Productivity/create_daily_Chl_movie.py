
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata
from pyproj import Proj, Transformer

def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:2] == '34' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()

    points = np.column_stack([XC.ravel(), YC.ravel()])
    points = reproject_polygon(points, 4326, 3413)
    X_3413 = np.reshape(points[:, 0], np.shape(XC))
    Y_3413 = np.reshape(points[:, 1], np.shape(YC))

    return(XC, YC, X_3413, Y_3413, Depth)

def get_dates():

    year = 2021
    start_month = 6
    end_month = 7

    date_strs = []
    for month in range(start_month, end_month+1):
        days_in_month = 31
        if month in [4,6,9,11]:
            days_in_month = 30
        if month ==2:
            if (year%4==0 and year%100!=0) or year%400==0:
                days_in_month = 29
            else:
                days_in_month = 28
        for day in range(1,days_in_month+1):
            date_str = f'{year}{month:02d}{day:02d}'
            date_strs.append(date_str)

    return date_strs

def plot_panel(output_dir, output_file, chl_dir, date_str,
               XC, YC, X_3413, Y_3413, Depth):

    chl_file = os.path.join(chl_dir, date_str[:4],
                            f'ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-{date_str}-fv6.0.nc')
    ds = nc4.Dataset(chl_file)
    chl = ds.variables['chlor_a'][:,:,:]
    lon = ds.variables['lon'][:]
    lat = ds.variables['lat'][:]
    ds.close()
    chl = chl[0, :, :]

    # plt.pcolormesh(chl)
    # plt.show()

    # negative_lon_indices = lon < 0
    # lon = np.concatenate((lon[~negative_lon_indices], lon[negative_lon_indices] + 360))
    # chl = np.hstack((chl[:, ~negative_lon_indices], chl[:, negative_lon_indices]))
    # chl = np.ma.masked_where(chl >100, chl)

    # plt.plot(lon)
    # plt.show()

    # plt.pcolormesh(lon, lat, chl)
    # plt.show()

    min_lon = np.min(XC)-1
    max_lon = np.max(XC)+1
    min_lat = np.min(YC)-1
    max_lat = np.max(YC)+1

    lon_indices = np.logical_and(lon >= min_lon, lon <= max_lon)
    lat_indices = np.logical_and(lat >= min_lat, lat <= max_lat)

    lon = lon[lon_indices]
    lat = lat[lat_indices]
    chl = chl[np.ix_(lat_indices, lon_indices)]
    chl[chl>100] = np.nan

    Lon,Lat = np.meshgrid(lon, lat)
    reprojected_points = reproject_polygon(np.column_stack([Lon.ravel(), Lat.ravel()]), 4326, 3413)
    X_chl = np.reshape(reprojected_points[:, 0], np.shape(Lon))
    Y_chl = np.reshape(reprojected_points[:, 1], np.shape(Lat))
    chl = griddata((X_chl.ravel(), Y_chl.ravel()), chl.ravel(), (X_3413, Y_3413), method='linear')

    fig = plt.figure(figsize=(8,6))
    plt.style.use('dark_background')

    chlmin = np.log10(0.05)
    chlmax = np.log10(10)

    gs = GridSpec(1, 21, left = 0.1, right = 0.9, bottom=0.1, top=0.9, wspace=0.05)

    ax = fig.add_subplot(gs[:, :-2])
    plt.pcolormesh(X_3413, Y_3413, np.log10(chl), vmin=chlmin, vmax=chlmax, cmap='turbo')

    plt.contour(X_3413, Y_3413, Depth, levels=[0], colors='white', linewidths=0.5)

    # plt.arrow(obs_point[0]-0.5, obs_point[1]-0.5, 0.5, 0.5, head_width=0.2, head_length=0.2, fc='white', ec='white')
    # plt.text(obs_point[0]-0.51, obs_point[1]-0.51, 'Ship Position on 19 June 2022\n (27°13.3\'N, 178°11.2\'E)',
    #          color='white', fontsize=8, verticalalignment='top', horizontalalignment='right')

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Upernavik Chl-a Estimate (OC-CCI v6) - '+date_str[:4]+'-'+date_str[4:6]+'-'+date_str[6:8])

    cbar_ax = fig.add_subplot(gs[:, -1:])

    x = np.array([0, 1])
    y = np.linspace(chlmin, chlmax, 100)
    X, Y = np.meshgrid(x, y)
    C = cbar_ax.pcolormesh(X, Y, Y, cmap='turbo', vmin=chlmin, vmax=chlmax)
    cbar_ax.set_xticks([])
    # put ticks and label on right
    cbar_ax.yaxis.set_ticks_position('right')
    cbar_ax.set_ylabel('Chl-a Concentration (mg m$^{-3}$)')
    # put the y label on the right
    cbar_ax.yaxis.set_label_position("right")
    cbar_ax.set_yticks(np.log10([0.05, 1, 0.5, 1, 5, 10]))
    cbar_ax.set_yticklabels(['0.05', '1', '0.5', '1', '5', '10'])


    plt.savefig(os.path.join(output_dir,'Panels', output_file))
    plt.close(fig)

def compile_panels_to_movie(output_dir):
    pwd = os.getcwd()

    panels_dir = os.path.join(output_dir,'Panels')

    os.chdir(panels_dir)

    output_name = 'Upernavik_Chl_Summer_2021.mp4'

    os.system("ffmpeg -r 5 -i chl_panel_%03d.png -vcodec mpeg4 -b 3M -y " + output_name)
    os.rename(output_name, os.path.join('..', output_name))

    os.chdir(pwd)

chl_dir = '/Users/mhwood/Desktop/Chl'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/' \
             'Downscale_Greenland/MITgcm/configurations/downscale_greenland'


XC, YC, X_3413, Y_3413, Depth = read_grid_geometry_from_nc(config_dir, 'L2_Upernavik')

output_dir = os.path.join(project_dir, 'Figures', 'Ocean', 'Chlorophyll')

date_strs = get_dates()

counter = 1
for date_str in date_strs:
    output_file = f'chl_panel_{counter:03d}.png'
    if output_file not in []:#os.listdir(os.path.join(output_dir,'Panels')):
        print('Creating panel for ' + date_str)
        plot_panel(output_dir, output_file, chl_dir, date_str,
                   XC, YC, X_3413, Y_3413, Depth)
    counter+=1

compile_panels_to_movie(output_dir)



