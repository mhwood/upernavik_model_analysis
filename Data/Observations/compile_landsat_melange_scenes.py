
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata


def get_melange_region_bounds(location):

    if location=='Upernavik_Fjord':
        min_x = -342356
        min_y = -1859037
        distance = 48000
        max_x = min_x + distance
        max_y = min_y + distance
    if location=='Ussing_Braeer':
        min_x = -342541
        min_y = -1737970
        distance = 30000
        max_x = min_x + distance
        max_y = min_y + distance
    if location=='Kakivfaat':
        min_x = -339699
        min_y = -1781649
        distance = 25000
        max_x = min_x + distance
        max_y = min_y + distance

    x_3413 = np.arange(min_x, max_x + 50.0, 50.0)
    y_3413 = np.arange(min_y, max_y + 50.0, 50.0)
    X_3413, Y_3413 = np.meshgrid(x_3413, y_3413)

    return X_3413, Y_3413

def grid_cell_ID_to_cell_bounds_function(cell_ID):

    row = int(cell_ID.split('_')[1])
    col = int(cell_ID.split('_')[3])

    x_0 = -652925
    y_0 = -3384425

    min_x = x_0+col*50000
    max_x = min_x + 50000
    min_y = y_0 + row * 50000
    max_y = min_y + 50000

    extents =[min_x,min_y,max_x,max_y]
    return(extents)

def check_rectangle_overlap(extents_1, extents_2):
    min_x1, min_y1, max_x1, max_y1 = extents_1
    min_x2, min_y2, max_x2, max_y2 = extents_2

    if (min_x1 < max_x2 and max_x1 > min_x2 and
        min_y1 < max_y2 and max_y1 > min_y2):
        return True
    else:
        return False

def get_cell_ID_list_from_coordinates(X, Y):

    min_scene_x = np.min(X)
    max_scene_x = np.max(X)
    min_scene_y = np.min(Y)
    max_scene_y = np.max(Y)
    scene_extents = [min_scene_x, min_scene_y, max_scene_x, max_scene_y]

    min_row = 1e22
    max_row = 0
    min_col = 1e22
    max_col = 0

    cell_ID_list = []
    for r in range(0,50):
        for c in range(0,50):
            cell_ID = 'r_{:02d}_c_{:02d}'.format(r,c)
            cell_extents = grid_cell_ID_to_cell_bounds_function(cell_ID)
            if check_rectangle_overlap(scene_extents, cell_extents):
                cell_ID_list.append(cell_ID)
                if r < min_row:
                    min_row = r
                if r > max_row:
                    max_row = r
                if c < min_col:
                    min_col = c
                if c > max_col:
                    max_col = c

    cell_ID_list_gridded = []
    for r in range(min_row, max_row+1):
        row_IDs = []
        for c in range(min_col, max_col+1):
            cell_ID = 'r_{:02d}_c_{:02d}'.format(r,c)
            if cell_ID in cell_ID_list:
                row_IDs.append(cell_ID)
        cell_ID_list_gridded.append(row_IDs)

    return(cell_ID_list, cell_ID_list_gridded)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    drF = ds.variables['drF'][:]
    hFaC = ds.variables['HFacC'][:, :, :]
    ds.close()
    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2
    return(XC, YC, Z, Depth, hFaC)

def get_dict_of_cell_scene_dates(cell_dir, cell_IDs):

    date_file_list_dict = {}
    for year in range(2016,2021):
        for month in range(1,13):
            n_days = 31
            if month in [4,6,9,11]:
                n_days = 30
            if month == 2:
                if (year%4==0 and year%100!=0) or (year%400==0):
                    n_days = 29
                else:
                    n_days = 28
            for day in range(1,n_days+1):
                date_str = '{:04d}{:02d}{:02d}'.format(year,month,day)
                date_file_list_dict[date_str] = []

    for cell_id in cell_IDs:
        scenes_dir = os.path.join(cell_dir, cell_id,'Imagery','Data','Scenes')
        if os.path.exists(scenes_dir):

            for file_name in os.listdir(scenes_dir):
                if 'Landsat_8' in file_name and file_name[0]!='.':
                    ds = nc4.Dataset(os.path.join(scenes_dir,file_name))
                    dates = list(ds.groups.keys())
                    for date in dates:
                        date_str = date.split('_')[3]
                        if int(date_str[:4]) >= 2016 and int(date_str[:4]) <= 2020:
                            date_file_list_dict[date_str].append(file_name)

    return(date_file_list_dict)

def compile_scene_composite_mask(cell_dir, cell_IDs_gridded):

    full_mask = 2*np.ones( (len(cell_IDs_gridded)*500, len(cell_IDs_gridded[0])*500) )
    full_x = np.zeros( (len(cell_IDs_gridded[0])*500, ) )
    full_y = np.zeros( (len(cell_IDs_gridded)*500, ) )

    for r, cell_IDs_row in enumerate(cell_IDs_gridded):
        for c, cell_id in enumerate(cell_IDs_row):

            mask_file = os.path.join(cell_dir, cell_id, 'Masks', cell_id+' BedMachineV5 Mask.nc')

            if os.path.exists(mask_file):
                ds = nc4.Dataset(mask_file)
                mask = ds.variables['mask'][:,:]
                x = ds.variables['x'][:]
                y = ds.variables['y'][:]
                ds.close()

                full_mask[r*500:(r+1)*500, c*500:(c+1)*500] = mask
                full_x[c*500:(c+1)*500] = x
                full_y[r*500:(r+1)*500] = y

    return(full_mask, full_x, full_y)

def compile_date_scene_composite(cell_dir, cell_IDs_gridded, date_str, file_list):

    full_scene_band_3 = np.zeros( (len(cell_IDs_gridded)*1000, len(cell_IDs_gridded[0])*1000) )
    full_scene_band_QA = np.zeros( (len(cell_IDs_gridded)*1000, len(cell_IDs_gridded[0])*1000) )
    full_scene_x = np.zeros( (len(cell_IDs_gridded[0])*1000, ) )
    full_scene_y = np.zeros( (len(cell_IDs_gridded)*1000, ) )

    scene_IDs = []

    for r, cell_IDs_row in enumerate(cell_IDs_gridded):
        for c, cell_id in enumerate(cell_IDs_row):

            for file_name in file_list:
                if file_name.startswith(cell_id):

                    scene_file = os.path.join(cell_dir, cell_id,'Imagery','Data','Scenes',
                                              file_name)

                    ds = nc4.Dataset(scene_file)
                    scene_x = ds.variables['x'][:]
                    scene_y = ds.variables['y'][:]
                    grps = list(ds.groups.keys())
                    for grp in grps:
                        if date_str in grp:
                            scene = ds.groups[grp]
                            if grp not in scene_IDs:
                                scene_IDs.append(grp)
                            scene_3 = scene.variables['band_3'][:,:]
                            scene_QA = scene.variables['band_QA'][:,:]
                            apply_mask = np.logical_and(np.logical_and(scene_3>0, scene_3<65535), ~np.isnan(scene_3))
                            # print(f'         - Applying scene in {100*np.sum(apply_mask)/np.size(apply_mask)}% of valid pixels')

                            full_scene_band_3[r*1000:(r+1)*1000, c*1000:(c+1)*1000][apply_mask] = scene_3[apply_mask]
                            full_scene_band_QA[r*1000:(r+1)*1000, c*1000:(c+1)*1000][apply_mask] = scene_QA[apply_mask]
                            full_scene_x[c*1000:(c+1)*1000] = scene_x
                            full_scene_y[r*1000:(r+1)*1000] = scene_y

                            # if date_str=='20200904':
                            #     plt.subplot(1,2,1)
                            #     plt.pcolormesh(scene_3)
                            #     plt.subplot(1, 2, 2)
                            #     plt.pcolormesh(full_scene_band_3)
                            #     plt.show()

                    ds.close()

    return(full_scene_band_3, full_scene_band_QA, full_scene_x, full_scene_y, scene_IDs)

def write_scene_composites_to_nc(project_dir, region, scene_mask, scene_dict):

    output_file = os.path.join(project_dir, 'Data', 'Observations','Imagery',
                               region+'_Landsat_8_Scenes.nc')


    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('x',len(scene_x))
    ds.createDimension('y',len(scene_y))

    xvar = ds.createVariable('x','f4',('x',))
    xvar[:] = scene_x

    yvar = ds.createVariable('y','f4',('y',))
    yvar[:] = scene_y

    mask_var = ds.createVariable('mask', 'f4', ('y', 'x'))
    mask_var[:, :] = scene_mask

    for date_str in scene_dict.keys():
        grp = ds.createGroup(date_str)
        scene_band_3 = scene_dict[date_str]['band_3']
        scene_band_QA = scene_dict[date_str]['band_QA']

        band_3_var = grp.createVariable('band_3', 'f4', ('y','x'))
        band_3_var[:,:] = scene_band_3

        band_QA_var = grp.createVariable('band_QA', 'f4', ('y','x'))
        band_QA_var[:,:] = scene_band_QA

        grp.scenes = ','.join(scene_dict[date_str]['scene_IDs'])

    ds.close()



config_dir = '/Users/mhwood/Documents/Research/Projects/Ocean_Modelling/Projects/Downscaled_Darwin/' \
             'darwin3/configurations/downscaled_ecco_v5_darwin/'


cell_dir = '/Volumes/backups/Research Backup/Projects/Greenland Dynamics/Data/Grids/'

XC, YC, Z, Depth, hFaC = read_grid_geometry_from_nc(config_dir, 'L2_Upernavik')

region = 'Upernavik_Fjord'

for region in ['Upernavik_Fjord']:#, 'Kakivfaat']: #'Upernavik_Fjord',
    print('Working on region: '+region)

    X_3413, Y_3413 = get_melange_region_bounds(region)

    print('    - Determining available scenes...')
    # cell_IDs = []
    # cell_IDs_gridded = []
    # for r in range(30,32):
    #     rows = []
    #     for c in range(6,8):
    #         cell_IDs.append('r_{:02d}_c_{:02d}'.format(r,c))
    #         rows.append('r_{:02d}_c_{:02d}'.format(r,c))
    #     cell_IDs_gridded.append(rows)

    cell_IDs, cell_IDs_gridded = get_cell_ID_list_from_coordinates(X_3413, Y_3413)
    print('      - cell IDs:',cell_IDs_gridded)


    date_file_list_dict = get_dict_of_cell_scene_dates(cell_dir, cell_IDs)

    # print(date_file_list_dict)

    # max_date = ''
    # max_dates = 0
    # for date in date_file_list_dict.keys():
    #     n_files = len(date_file_list_dict[date])
    #     if n_files > max_dates:
    #         max_dates = n_files
    #         max_date = date
    #
    # max_dates = 1
    print('    - Checking which files have the correct scene information')
    date_file_list_dict_small = {}
    for date in date_file_list_dict.keys():
        n_files = len(date_file_list_dict[date])
        if n_files>0:
            # print('Date with max scenes: '+date+' with '+str(n_files)+' scenes')
            date_file_list_dict_small[date] = date_file_list_dict[date]
    # print('    - Using all dates with at least '+str(max_dates)+' scenes for composites.')

    print(date_file_list_dict_small)

    # if max_dates==0:
    #     print('No scenes found for region '+region+'. Skipping to next region.')
    #     continue
    # else:
    print('    - Compiling scene composites...')
    scene_dict = {}
    for d, date_str in enumerate(list(date_file_list_dict_small.keys())):
        print('      - Working on date '+date_str+' ('+str(d)+' of '+str(len(date_file_list_dict_small))+' scenes)')
        scene_band_3, scene_band_QA, scene_x, scene_y, scene_IDs = compile_date_scene_composite(cell_dir, cell_IDs_gridded,
                                                                                    date_str,
                                                                                    date_file_list_dict_small[date_str])
        scene_dict[date_str] = {}
        scene_dict[date_str]['band_3'] = scene_band_3
        scene_dict[date_str]['band_QA'] = scene_band_QA
        scene_dict[date_str]['x'] = scene_x
        scene_dict[date_str]['y'] = scene_y
        scene_dict[date_str]['scene_IDs'] = scene_IDs
        print('            - Min and max of band 3 in composite: '+str(np.min(scene_band_3))+' and '+str(np.max(scene_band_3)))

        # if date_str=='20200904':
        #     C = plt.pcolormesh(scene_band_3, vmin=7000, vmax=20000)
        #     plt.colorbar(C)
        #     plt.show()

    print('    - Compiling scene mask and interpolating...')
    scene_mask, mask_x, mask_y = compile_scene_composite_mask(cell_dir, cell_IDs_gridded)

    # C = plt.pcolormesh(scene_mask)
    # plt.colorbar(C)
    # plt.show()

    mask_X, mask_Y = np.meshgrid(mask_x, mask_y)
    X, Y = np.meshgrid(scene_x, scene_y)

    scene_mask = griddata( np.column_stack((mask_X.ravel(), mask_Y.ravel())),
                           scene_mask.ravel(),
                           (X, Y),
                           method='nearest')

    # C = plt.pcolormesh(scene_mask)
    # plt.colorbar(C)
    # plt.show()

    print('    - Writing scene composites to netCDF...')
    project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'
    write_scene_composites_to_nc(project_dir, region, scene_mask, scene_dict)





