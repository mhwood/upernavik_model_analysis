
import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
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

def read_annual_data_stack_to_model_grid(data_dir, year, XC, YC, X_3413, Y_3413, Depth):

    first_file = True
    if year%4==0:
        seaice_stack = np.zeros((366,np.shape(X_3413)[0], np.shape(X_3413)[1]))
    else:
        seaice_stack = np.zeros((365, np.shape(X_3413)[0], np.shape(X_3413)[1]))

    counter = 0

    for month in range(1, 13):
        print('  - Reading in data in month '+str(month)+' in year '+str(year))
        if month in [1, 3, 5, 7, 8, 10, 12]:
            n_days = 31
        elif month in [4, 6, 9, 11]:
            n_days = 30
        else:
            if year % 4 == 0:
                n_days = 29
            else:
                n_days = 28

        yr = str(year)
        mo = '{:02d}'.format(month)

        for day in range(1,n_days+1):
            # try:
            dy = '{:02d}'.format(day)
            # if year>=2008:
            #     file_path = os.path.join(data_dir, 'seaice_conc_daily_nh_'+yr+mo+dy+'_f17_v04r00.nc')
            # else:
            #     file_path = os.path.join(data_dir, 'seaice_conc_daily_nh_' + yr + mo + dy + '_f13_v04r00.nc')
            file_path = os.path.join(data_dir, 'sic_psn25_'+yr+mo+dy+'_F17_v05r00.nc')
            print('       - Reading '+'sic_psn25_'+yr+mo+dy+'_F17_v05r00.nc')

            ds = nc4.Dataset(file_path)
            # if first_file:
            #     x = ds.variables['xgrid'][:]
            #     y = ds.variables['ygrid'][:]
            # seaice = ds.variables['nsidc_nt_seaice_conc'][:,:,:]
            if first_file:
                x = ds.variables['x'][:]
                y = ds.variables['y'][:]
            seaice = ds.variables['cdr_seaice_conc'][:,:,:]
            ds.close()
            seaice = np.array(seaice[0, :, :])

            if first_file:
                X_3411,Y_3411 = np.meshgrid(x,y)
                # plt.pcolormesh(X_3411, Y_3411, seaice)
                # plt.show()

                points = np.column_stack([X_3411.ravel(), Y_3411.ravel()])
                points = reproject_polygon(points, 3411, 3413)
                X_3413_seaice = np.reshape(points[:, 0],np.shape(X_3411))
                Y_3413_seaice = np.reshape(points[:, 1], np.shape(Y_3411))

                first_file = False

            # fill the inland points to avoid interpolation issues?
            points_nonnan = np.column_stack([X_3411.ravel(),Y_3411.ravel()])
            seaice_nonnan = np.ravel(seaice)
            non_zero_locations = seaice_nonnan<2
            points_nonnan = points_nonnan[non_zero_locations,:]
            seaice_nonnan = seaice_nonnan[non_zero_locations]
            seaice_nearest = griddata(points_nonnan, seaice_nonnan, (X_3411, Y_3411),
                                      method='nearest')

            # plt.subplot(1, 2, 1)
            # C = plt.pcolormesh(seaice, vmin=0, vmax=1)
            # plt.colorbar(C, orientation='horizontal')
            # plt.subplot(1, 2, 2)
            # C=plt.pcolormesh(seaice_nearest)
            # plt.colorbar(C, orientation='horizontal')
            # plt.show()

            seaice_3413= griddata(np.column_stack([X_3413_seaice.ravel(),Y_3413_seaice.ravel()]),
                                  seaice_nearest.ravel(), (X_3413, Y_3413))
            seaice_3413[Depth<=0]=254

            seaice_stack[counter,:,:] = seaice_3413

            # plt.subplot(1, 2, 1)
            # C = plt.pcolormesh(X_3413_seaice, Y_3413_seaice, seaice,vmin=0,vmax=1)
            # plt.gca().set_xlim([np.min(X_3413), np.max(X_3413)])
            # plt.gca().set_ylim([np.min(Y_3413), np.max(Y_3413)])
            # plt.colorbar(C, orientation='horizontal')
            # plt.subplot(1, 2, 2)
            # C=plt.pcolormesh(X_3413, Y_3413, seaice_3413,vmin=0,vmax=1)
            # plt.colorbar(C, orientation='horizontal')
            # plt.show()
            # except:
            #     print('   There was an issue with this file')

            counter +=1

    return(seaice_stack)

def write_data_to_annual_nc(output_file,X,Y,seaice):

    ds= nc4.Dataset(output_file,'w')
    tdim = ds.createDimension('time',np.shape(seaice)[0])
    ydim = ds.createDimension('y',np.shape(seaice)[1])
    xdim = ds.createDimension('x', np.shape(seaice)[2])

    x = ds.createVariable('x','f4',('x',))
    x[:] = X[0,:]

    y = ds.createVariable('y', 'f4', ('y',))
    y[:] = Y[:,0]

    d = ds.createVariable('days', 'f4', ('time',))
    d[:] = np.arange(1,np.shape(seaice)[0]+1)

    s = ds.createVariable('seaice_conc', 'f4', ('time','y', 'x'))
    s[:,:,:] = seaice

    ds.close()

def read_data_from_annual_nc(input_file):

    ds = nc4.Dataset(input_file, 'r')
    seaice = ds.variables['seaice_conc'][:,:,:]
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    days = ds.variables['days'][:]
    ds.close()

    return seaice, x, y, days


# config_dir = '/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/darwin3/' \
#              'configurations/downscale_darwin'
config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/' \
             'Downscale_Greenland/MITgcm/configurations/downscale_greenland'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'
# project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Disko Bay'

data_dir = '/Volumes/CoOL/Data_Repository/Arctic/Sea_Ice/daily/v5/'

XC, YC, X_3413, Y_3413, Depth = read_grid_geometry_from_nc(config_dir, model_name='L2_Upernavik')

for year in range(2021,2023):

    file_name = 'Upernavik Sea Ice '+str(year)+'.nc'

    if file_name not in os.listdir(os.path.join(project_dir,'Data','Observations')):
        print('Stacking data in year ' + str(year))

        seaice_stack = read_annual_data_stack_to_model_grid(data_dir, year, XC, YC, X_3413, Y_3413, Depth)

        output_file = os.path.join(project_dir,'Data','Observations',file_name)
        # print(output_file)

        write_data_to_annual_nc(output_file,X_3413, Y_3413,seaice_stack)

    seaice, x, y, days = read_data_from_annual_nc(os.path.join(project_dir,'Data','Observations',file_name))

    seaice[seaice>1]=1
    seaice = np.array(seaice)

    short_file_name = 'Upernavik_Sea_Ice_'+str(year)
    bin_file = os.path.join(project_dir,'Data','Observations',short_file_name)
    seaice.ravel('C').astype('>f4').tofile(bin_file)

