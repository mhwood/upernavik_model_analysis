
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


def read_model_grid(config_dir, model_name):
    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',model_name+'_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()
    return(XC, YC)


def read_L2_CTD_locations(config_dir, L2_model_name, L2_XC, L2_YC):

    dv_file = os.path.join(config_dir,'L2',L2_model_name,'input','dv','CTD_mask.bin')

    dv_grid = np.fromfile(dv_file,'>f4').reshape(np.shape(L2_XC))

    rows, cols = np.where(dv_grid!=0)

    locations = np.zeros((len(rows),2))
    for i in range(len(rows)):
        locations[i,0] = L2_XC[rows[i],cols[i]]
        locations[i,1] = L2_YC[rows[i],cols[i]]

    return(locations)


def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)


def read_L2_timeseries_at_CTD_locations(config_dir, L2_model_name, locations, L2_XC, L2_YC):

    # find the indices of the locations in L2
    rows = np.zeros((np.shape(locations)[0],)).astype(int)
    cols = np.zeros((np.shape(locations)[0],)).astype(int)
    for ll in range(np.shape(locations)[0]):
        dist = great_circle_distance(locations[ll,0], locations[ll,1], L2_XC.ravel(), L2_YC.ravel())
        dist = np.reshape(dist,np.shape(L2_XC))
        row,col = np.where(dist==np.min(dist))
        rows[ll]=row[0]
        cols[ll]=col[0]

    print(rows)

    Nr = 51

    output_dir = os.path.join(config_dir,'L2',L2_model_name,'results','state_3D_mon_mean')
    output_files = os.listdir(output_dir)

    # loop through the files and stack em into a big grid
    L2_timeseries_theta = np.zeros((np.shape(locations)[0], 7*12, Nr))
    L2_timeseries_salt = np.zeros((np.shape(locations)[0], 7 * 12, Nr))
    iterations = np.zeros((7*12,))
    counter = -1
    for year in range(2015,2017):
        for month in range(1,13):
            counter+=1
            datestr = str(year)+'{:02d}'.format(month)
            for file_name in output_files:
                if file_name[0]!='.' and file_name[-3:]=='.nc':
                    if file_name.split('.')[1]==datestr:
                        print('    - Reading from '+file_name)
                        ds = nc4.Dataset(os.path.join(output_dir,file_name))
                        iters = ds.variables['iterations'][:]
                        Theta = ds.variables['Theta'][:, : , :, :]
                        Salt = ds.variables['Salt'][:, :, :, :]
                        if counter==10:
                            depths = ds.variables['depths'][:]
                        ds.close()

                        iterations[counter] = iters[0]

                        for ll in range(len(rows)):
                            L2_timeseries_theta[ll,counter,:] = Theta[0,:, rows[ll], cols[ll]]
                            L2_timeseries_salt[ll, counter, :] = Salt[0, :, rows[ll], cols[ll]]

    return(iterations, depths, L2_timeseries_theta, L2_timeseries_salt)


def store_timeseries_as_nc(project_dir, iterations, depths, L2_timeseries_theta, L2_timeseries_salt, L2_CTD_locations):

    output_file = os.path.join(project_dir, 'Data', 'Ocean', 'Modeling', 'L2_state_at_L2_CTDs.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('iteration', len(iterations))
    ds.createDimension('depth', len(depths))
    ds.createDimension('loc_num', np.shape(L2_CTD_locations)[0])

    t = ds.createVariable('iterations', 'f4', ('iteration',))
    t[:] = iterations

    t = ds.createVariable('depths', 'f4', ('depth',))
    t[:] = depths

    t = ds.createVariable('loc_num', 'f4', ('loc_num',))
    t[:] = np.arange(np.shape(L2_CTD_locations)[0])+1

    t = ds.createVariable('Theta', 'f4', ('loc_num','iteration','depth'))
    t[:,:,:] = L2_timeseries_theta

    t = ds.createVariable('Salt', 'f4', ('loc_num', 'iteration', 'depth'))
    t[:,:,:] = L2_timeseries_salt

    t = ds.createVariable('longitude', 'f4', ('loc_num',))
    t[:] = L2_CTD_locations[:,0]

    t = ds.createVariable('latitude', 'f4', ('loc_num',))
    t[:] = L2_CTD_locations[:, 1]

    ds.close()



config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
             'configurations/downscaled_greenland'


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

L2_model_name = 'L2_Upernavik'

L2_XC, L2_YC = read_model_grid(config_dir, L2_model_name)

L2_CTD_locations = read_L2_CTD_locations(config_dir, L2_model_name, L2_XC, L2_YC)

iterations, depths, L2_timeseries_theta, L2_timeseries_salt = \
    read_L2_timeseries_at_CTD_locations(config_dir, L2_model_name, L2_CTD_locations, L2_XC, L2_YC)


# store_timeseries_as_nc(project_dir, iterations, depths, L2_timeseries_theta, L2_timeseries_salt, L2_CTD_locations)