
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


def read_landsat_8_scene(project_dir, date_str):

    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Observations', 'Imagery',
                                               'Upernavik_Fjord_Landsat_8_Scenes.nc') )
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    mask = ds.variables['mask'][:, :]

    group = ds.groups[date_str]
    band_3 = group.variables['band_3'][:, :]
    band_3 = np.array(band_3)

    band_3 = band_3 / 65535
    ds.close()

    return(band_3, mask, x, y)

def plot_landsat_8_scene_example(project_dir, date_str, band_3, mask, x, y):

    output_file = os.path.join(project_dir, 'Figures','Ocean','Melange', f'Upernavik_Landsat8_Band3_{date_str}.png')


    plt.figure(figsize=(8,8))
    plt.pcolormesh(x, y, band_3, cmap='gray', vmin=0, vmax=0.5)
    # plt.colorbar(label='Landsat 8 Band 3 Reflectance')
    plt.contour(x, y, mask, levels=[0.5], colors='yellow', linewidths=0.5)
    plt.title(f'Upernavik Fjord - Landsat 8 Band 3 on {date_str}')
    plt.xlabel('Distance East (km, EPSG: 3413)')
    plt.ylabel('Distance North (km, EPSG: 3413)')
    #plt.axis('equal')

    plt.savefig(output_file, dpi=300)
    plt.close()

def read_melange_fraction(project_dir):
    ds = nc4.Dataset(os.path.join(project_dir, 'Data', 'Observations', 'Imagery',
                                               'Upernavik_Fjord_Landsat_8_Melange_Fraction.nc') )
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    mask = ds.variables['mask'][:, :]
    melange_fraction = ds.variables['melange_fraction'][:, :]
    ds.close()

    return(melange_fraction, mask, x, y)

def plot_melange_fraction_example(project_dir, melange_fraction, band_3, mask, x, y):

    output_file = os.path.join(project_dir, 'Figures','Ocean','Melange',
                               f'Upernavik_Landsat8_Melange_Fraction.png')

    plt.figure(figsize=(10, 8))
    plt.pcolormesh(x, y, band_3, cmap='gray', vmin=0, vmax=0.5)
    melange_plot = np.ma.masked_where(melange_fraction<0.3, melange_fraction)
    C = plt.pcolormesh(x, y, melange_plot, cmap='turbo', vmin=0.3, vmax=1)
    plt.colorbar(C, label='Melange Fraction')
    plt.contour(x, y, mask, levels=[0.5], colors='w', linewidths=0.2)
    plt.title(f'Upernavik Fjord - Landsat 8 Band 3 on {date_str}')
    plt.xlabel('Distance East (km, EPSG: 3413)')
    plt.ylabel('Distance North (km, EPSG: 3413)')
    # plt.axis('equal')

    plt.savefig(output_file, dpi=300)
    plt.close()

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

date_str = '20200915'
band_3, mask, x, y = read_landsat_8_scene(project_dir, date_str)

x /= 1000
y /= 1000

# plot_landsat_8_scene_example(project_dir, date_str, band_3, mask, x, y)

melange_fraction, _, _, _ = read_melange_fraction(project_dir)

plot_melange_fraction_example(project_dir, melange_fraction, band_3, mask, x, y)

