
import os
import netCDF4 as nc4
import numpy as np
from pyproj import Transformer
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

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

def read_polygons_from_nc(project_dir):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Upernavik N Fjord Geometry Transect.nc'))

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

def retrieve_model_obs_pairs(ctd_dir, project_folder):

    mapping_file = os.path.join(project_folder, 'Data', 'Ocean',  'AXCTDs',
                                'AXCTD List at Model Points.csv')
    f = open(mapping_file, 'r')
    lines = f.readlines()
    f.close()
    lines.pop(0)

    obs_profiles = {}

    for line in lines:
        line = line.split(',')
        file_name = line[0]

        obs_file_path = os.path.join(ctd_dir, file_name[:4], 'CTD_'+file_name+'.nc')
        ds = nc4.Dataset(obs_file_path)
        depth = ds.variables['depth'][:]
        theta_obs = ds.variables['potential_temperature'][:]
        salt_obs = ds.variables['practical_salinity'][:]
        lon = ds.longitude
        lat = ds.latitude
        ds.close()

        x, y = reproject_polygon(np.array([[lon, lat]]), 4326, 3413)[0]

        profile = np.column_stack((depth, theta_obs, salt_obs))

        obs_profiles[file_name] = {'profile': profile, 'year': int(file_name[:4]), 'x': x, 'y': y}

    return(obs_profiles)

def create_interpolated_transects(obs_profiles, year, distance, transect_x, transect_y, depth):

    Distance, Depth = np.meshgrid(distance, depth)
    first_profile = True

    ctd_depths = []
    ctd_distances = []

    for file_name in list(obs_profiles.keys()):
        profile = obs_profiles[file_name]['profile']
        profile_year = obs_profiles[file_name]['year']
        x = obs_profiles[file_name]['x']
        y = obs_profiles[file_name]['y']

        if profile_year== year:
            # Find closest point on transect
            dists = np.sqrt((transect_x - x)**2 + (transect_y - y)**2)
            closest_index = np.argmin(dists)
            transect_distance = dists[closest_index]

            if transect_distance < 10000:
                print(f'    - Year: {year}, File: {file_name}, Transect Distance: {transect_distance:.2f} m')
                print('    - Closest Transect Point Distance: {:.2f} m'.format(distance[closest_index]))

                if first_profile:
                    all_points = np.column_stack([distance[closest_index]*np.ones((np.shape(profile)[0],)),
                                                  profile[:,0]])
                    all_theta = profile[:,1:2]
                    all_salt = profile[:,2:3]
                    first_profile = False
                else:
                    new_points = np.column_stack([distance[closest_index]*np.ones((np.shape(profile)[0],)),
                                                  profile[:,0]])
                    all_points = np.vstack([all_points, new_points])
                    all_theta = np.vstack([all_theta, profile[:,1:2]])
                    all_salt = np.vstack([all_salt, profile[:,2:3]])

                ctd_depths.append(np.max(profile[:,0]))
                ctd_distances.append(distance[closest_index])

    Theta_grid = griddata(all_points, all_theta.ravel(), (Distance, Depth), method='linear')
    Salt_grid = griddata(all_points, all_salt.ravel(), (Distance, Depth), method='linear')

    return(Theta_grid, Salt_grid, ctd_depths, ctd_distances)

def store_interpolated_transects_to_nc(project_folder, distance, depth,
                                       Theta_grids, Salt_grids, all_ctd_depths, all_ctd_distances):

    output_dir = os.path.join(project_folder, 'Data', 'Ocean', 'AXCTDs')
    file_name = os.path.join(output_dir, 'Upernavik_N_Fjord_AXCTD_Transects.nc')

    ds = nc4.Dataset(file_name, 'w', format='NETCDF4')

    ds.createDimension('distance', len(distance))
    ds.createDimension('depth', len(depth))

    distance_var = ds.createVariable('distance', 'f4', ('distance',))
    depth_var = ds.createVariable('depth', 'f4', ('depth',))

    distance_var[:] = distance
    depth_var[:] = depth

    for year in Theta_grids.keys():

        drp = ds.createGroup(str(year))

        drp.createDimension('ctd_profiles', len(all_ctd_depths[year]))


        theta_var = drp.createVariable('potential_temperature', 'f4', ('depth', 'distance'), fill_value=np.nan)
        salt_var = drp.createVariable('practical_salinity', 'f4', ('depth', 'distance'), fill_value=np.nan)
        ctd_depths_var = drp.createVariable('ctd_max_depths', 'f4', ('ctd_profiles',))
        ctd_distances_var = drp.createVariable('ctd_distances', 'f4', ('ctd_profiles',))


        theta_var[:, :] = Theta_grids[year]
        salt_var[:, :] = Salt_grids[year]
        ctd_depths_var[:] = all_ctd_depths[year]
        ctd_distances_var[:] = all_ctd_distances[year]

    ds.close()


def create_AXCTD_transects(project_folder, ctd_dir):

    model_name = 'L2_Upernavik'

    distance, transect_x, transect_y, surface, ice_polygon, bed_polygon = read_polygons_from_nc(project_folder)

    depth = np.arange(0, np.max(bed_polygon[:,1]))

    obs_profiles = retrieve_model_obs_pairs(ctd_dir, project_folder)

    years = [2016,2018,2019,2020]

    Theta_grids = {}
    Salt_grids = {}
    all_ctd_depths = {}
    all_ctd_distances = {}

    for year in years:
        print('Processing year '+str(year))
        Theta_grid, Salt_grid, ctd_depths, ctd_distances =\
                create_interpolated_transects(obs_profiles, year, distance, transect_x, transect_y, depth)

        Theta_grids[year] = Theta_grid
        Salt_grids[year] = Salt_grid
        all_ctd_depths[year] = ctd_depths
        all_ctd_distances[year] = ctd_distances

    store_interpolated_transects_to_nc(project_folder, distance, depth,
                                        Theta_grids, Salt_grids, all_ctd_depths, all_ctd_distances)















project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

ctd_dir = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_CTDs/Processed'

create_AXCTD_transects(project_folder, ctd_dir)







