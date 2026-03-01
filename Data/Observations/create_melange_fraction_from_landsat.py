
import os
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
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


    return x_3413, y_3413

def create_melange_mask_from_landsat(project_dir, region):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Observations', 'Imagery',
                                               region+'_Landsat_8_Scenes.nc') )
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    mask = ds.variables['mask'][:, :]

    melange_mask_sum = np.zeros(np.shape(mask))
    melange_mask_count = np.zeros(np.shape(mask))

    # if region=='Upernavik_Fjord':
    #     date_strs = ['20160531','20170813','20190929','20200915']
    # else:
    date_strs = list(ds.groups.keys())

    for grp in date_strs:
        group = ds.groups[grp]
        band_3 = group.variables['band_3'][:, :]
        band_3 = np.array(band_3)
        band_3[band_3<0] = 0
        band_3[band_3>65535] = 65535

        # band_QA = group.variables['band_QA'][:, :]
        # # plt.pcolormesh(band_QA)
        # # plt.title(grp)
        # # plt.show()
        # cloud_mask = band_QA<7000

        # rescale from landsat 8 bits to 0-1
        band_3 = band_3 / 65535
        melange_indices = band_3 > 0.2
        melange_mask_sum[melange_indices] += 1
        nonzero_indices = band_3 > 0
        melange_mask_count[nonzero_indices] += 1

    ds.close()

    melange_mask_fraction = np.zeros(np.shape(mask))
    valid_indices = melange_mask_count > 0
    melange_mask_fraction[valid_indices] = melange_mask_sum[valid_indices] / melange_mask_count[valid_indices]

    melange_mask_fraction[mask!=0] = 0

    return(melange_mask_fraction, mask, x, y)

def save_melange_mask_as_nc(project_dir, region, melange_mask_fraction, mask, x, y):
    output_file = os.path.join(project_dir, 'Data', 'Observations', 'Melange',
                               region+'_Landsat_8_Melange_Fraction.nc')

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('x', len(x))
    ds.createDimension('y', len(y))
    x_var = ds.createVariable('x', 'f4', ('x',))
    x_var[:] = x
    y_var = ds.createVariable('y', 'f4', ('y',))
    y_var[:] = y
    mask_var = ds.createVariable('mask', 'f4', ('y', 'x',))
    mask_var[:, :] = mask
    melange_var = ds.createVariable('melange_fraction', 'f4', ('y', 'x',))
    melange_var[:, :] = melange_mask_fraction
    ds.close()


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

for region in ['Upernavik_Fjord','Kakivfaat','Ussing_Braeer']:

    x_3413, y_3413 = get_melange_region_bounds(region)
    X_3413, Y_3413 = np.meshgrid(x_3413, y_3413)

    melange_mask_fraction, mask, x, y = create_melange_mask_from_landsat(project_dir, region)
    # melange_mask_fraction[melange_mask_fraction<0.3] = 0

    X, Y = np.meshgrid(x, y)
    points = np.column_stack([X.ravel(), Y.ravel()])
    melange_mask_fraction = griddata(points, melange_mask_fraction.ravel(), (X_3413, Y_3413), method='nearest')
    mask = griddata(points, mask.ravel(), (X_3413, Y_3413), method='nearest')

    save_melange_mask_as_nc(project_dir, region, melange_mask_fraction, mask, x_3413, y_3413)

    area_grid = np.ones(np.shape(mask))* (x[1]-x[0]) * (y[1]-y[0])
    total_melange_area = np.sum(area_grid[melange_mask_fraction>0.5])
    print(f'{region} Melange Area: {total_melange_area/1e6:.2f} km^2')

    # plt.subplot(1,2,1)
    # C = plt.pcolormesh(x_3413, y_3413, melange_mask_fraction, vmin=0, vmax=1, cmap='turbo')
    # plt.contour(x_3413, y_3413, mask, levels=[0.01], colors='w', linewidths=0.5)
    # plt.colorbar(C, label='Melange Fraction')
    #
    # plt.subplot(1, 2, 2)
    # C = plt.pcolormesh(x_3413, y_3413, melange_mask_fraction>0.5, vmin=0, vmax=1, cmap='turbo')
    # plt.contour(x_3413, y_3413, mask, levels=[0.01], colors='w', linewidths=0.5)
    # plt.colorbar(C, label='Melange Fraction')
    # plt.show()
