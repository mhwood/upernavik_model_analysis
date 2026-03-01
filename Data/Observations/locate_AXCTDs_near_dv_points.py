

import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc4


def read_model_geometry(config_dir):

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids','L2_Upernavik_grid.nc'))
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    return(XC,YC,Depth)


def read_model_dv_locations(config_dir, XC, YC):

    file_name = os.path.join(config_dir,'L2','L2_Upernavik','input','dv','CTD_mask.bin')

    grid = np.fromfile(file_name, '>f4').reshape(np.shape(XC))
    N = int(np.max(grid))

    ctd_locations = np.zeros((N,2))
    for i in range(1,N+1):
        rows, cols = np.where(grid==i)
        ctd_locations[i-1, 0] = XC[rows[0], cols[0]]
        ctd_locations[i-1, 1] = YC[rows[0], cols[0]]

    return(ctd_locations)


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


def find_axctds_near_dv_points(axctd_folder, model_ctd_locations, XC, YC):

    threshold = 5000

    ref_file = os.path.join(axctd_folder,'OMG_CTD_Locations.csv')
    f = open(ref_file)
    lines = f.read()
    f.close()

    lines = lines.split('\n')
    lines.pop(0)

    file_ids = []
    lons = []
    lats = []

    for line in lines:
        line = line.split(',')
        if len(line)>5:
            file_ids.append(line[0])
            lons.append(float(line[-2]))
            lats.append(float(line[-1]))

    lons = np.array(lons)
    lats = np.array(lats)

    file_id_lists = {}

    # plt.plot(model_ctd_locations[:,0], model_ctd_locations[:,1], 'go')
    # plt.plot(lons, lats, 'k.')
    # plt.gca().set_xlim([np.min(XC), np.max(XC)])
    # plt.gca().set_ylim([np.min(YC), np.max(YC)])
    # plt.show()

    for c in range(np.shape(model_ctd_locations)[0]):
        lon = model_ctd_locations[c,0]
        lat = model_ctd_locations[c,1]
        dist = great_circle_distance(lon, lat, lons, lats)
        locations = np.where(dist<threshold)[0]
        # print(locations)
        if len(locations)>0:
            location_files = []
            for ll in locations:
                location_files.append(file_ids[ll])
            file_id_lists[c+1] = location_files
            # print('Found '+str(len(location_files))+' for locations '+str(c+1))

    return(file_id_lists)


def write_file_id_lists_to_reference(project_folder, model_ctd_locations, file_id_lists):

    output_file = os.path.join(project_folder,'Data','Ocean','AXCTDs', 'AXCTD List at Model Points.csv')

    output = 'AXCTD File,Model Point,Longitude,Latitude'
    for c in range(np.shape(model_ctd_locations)[0]):
        if c+1 in list(file_id_lists.keys()):
            file_list = file_id_lists[c+1]
            for file_name in file_list:
                output+='\n'+file_name+','+str(c+1)+','+str(model_ctd_locations[c,0])+','+str(model_ctd_locations[c, 1])

    f = open(output_file, 'w')
    f.write(output)
    f.close()




project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

axctd_folder = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_CTDs'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/' \
             'configurations/downscaled_greenland'



XC, YC, Depth = read_model_geometry(config_dir)

model_ctd_locations = read_model_dv_locations(config_dir, XC, YC)

file_id_lists = find_axctds_near_dv_points(axctd_folder, model_ctd_locations,XC, YC)
print(file_id_lists)

write_file_id_lists_to_reference(project_folder, model_ctd_locations, file_id_lists)



