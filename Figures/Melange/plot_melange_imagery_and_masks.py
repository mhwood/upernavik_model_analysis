
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import netCDF4 as nc4
from matplotlib.colors import ListedColormap, BoundaryNorm

def read_melange_image(project_dir, location_dict, location):
    image_file = os.path.join(project_dir, 'Data','Observations','Imagery','Landsat',
                             '_'.join(location.split())+' Landsat Imagery.nc')
    with nc4.Dataset(image_file, 'r') as ds:
        band_2 = ds.variables['band_2'][:,:]
        band_3 = ds.variables['band_3'][:,:]
        band_4 = ds.variables['band_4'][:,:]
        image = np.stack([band_4, band_3, band_2], axis=2)
        y = ds.variables['Y'][:,:]
        x = ds.variables['X'][:,:]
        image_mask = ds.variables['bedmachine_mask'][:,:]

    # brighten the rgb image by normalizing and cli
    image = image.astype(float)
    image = (image - image.min()) / (image.max() - image.min())

    # image[image>0.9] = .0.9
    # image = (image - imagemin()) / (image.max() - image.min())

    brightness_factor = 0.7
    image = np.clip(image * (1 + brightness_factor), 0, 1)


    location_dict[location]['image'] = image
    location_dict[location]['image_mask'] = image_mask
    location_dict[location]['image_x'] = x
    location_dict[location]['image_y'] = y

    return location_dict

def read_melange_mask(project_dir, location_dict, location):
    mask_file = os.path.join(project_dir, 'Data','Observations','Melange',
                             '_'.join(location.split())+'_Landsat_8_Melange_Fraction.nc')
    with nc4.Dataset(mask_file, 'r') as ds:
        mask = ds.variables['mask'][:]
        fraction = ds.variables['melange_fraction'][:]
        y = ds.variables['y'][:]
        x = ds.variables['x'][:]
    fraction[fraction<0.5] = 0

    location_dict[location]['fraction_mask'] = mask
    location_dict[location]['fraction'] = fraction
    location_dict[location]['fraction_x'] = x
    location_dict[location]['fraction_y'] = y

    return location_dict

def read_rigidity_mask(project_dir, location_dict, location):
    mask_file = os.path.join(project_dir, 'Data','Observations','Melange',
                             '_'.join(location.split())+'_Melange_Rigidity_Mask.nc')
    with nc4.Dataset(mask_file, 'r') as ds:
        # mask = ds.variables['mask'][:]
        rigidity = ds.variables['rigidity_mask_fraction'][:]
        y = ds.variables['y'][:]
        x = ds.variables['x'][:]
    if location=='Upernavik Fjord':
        threshold = 0.5
    else:
        threshold = 0.5
    rigidity[rigidity<threshold] = 0
    rigidity[rigidity>=threshold] = 1
    rigidity[location_dict[location]['fraction']==0] = 0

    # location_dict[location]['rigidity_mask'] = mask
    location_dict[location]['rigidity'] = rigidity
    location_dict[location]['rigidity_x'] = x
    location_dict[location]['rigidity_y'] = y

    return location_dict

def plot_melange_mask(project_dir, locations, location_dict):

    fig = plt.figure(figsize=(7.5, 8))

    plot_height = 6
    v_spacing = 1
    colorbar_height = 1

    fraction_vmin = 0.5
    fraction_vmax = 1
    fraction_cmap = 'turbo'

    rigid_color = 'firebrick'
    drifting_color = 'steelblue'
    other_color = 'white'
    rigid_cmap = ListedColormap([other_color, drifting_color, rigid_color])

    gs = GridSpec(3*plot_height + colorbar_height + v_spacing,
                  3, left = 0.05, right=0.97, top=0.93, bottom=0.08)

    for row in range(len(locations)):

        image = location_dict[locations[row]]['image']
        image_mask = location_dict[locations[row]]['image_mask']
        image_x = location_dict[locations[row]]['image_x']
        image_y = location_dict[locations[row]]['image_y']
        fraction = location_dict[locations[row]]['fraction']
        fraction_x = location_dict[locations[row]]['fraction_x']
        fraction_y = location_dict[locations[row]]['fraction_y']
        fraction_mask = location_dict[locations[row]]['fraction_mask']
        rigidity = location_dict[locations[row]]['rigidity']
        rigidity_x = location_dict[locations[row]]['rigidity_x']
        rigidity_y = location_dict[locations[row]]['rigidity_y']
        # rigidity_mask = location_dict[locations[row]]['rigidity_mask']

        ######################################################################
        # Satellite imagery
        ax = fig.add_subplot(gs[plot_height*row:plot_height*(row+1), 0])
        ax.set_ylabel(locations[row])

        # ax.pcolormesh(image_x, image_y, image, cmap='gray', vmin=vmin, vmax=vmax)
        ax.imshow(image, extent=(image_x.min(), image_x.max(), image_y.min(), image_y.max()),
                  origin='lower')
        # ax.contour(image_x, image_y, image_mask, levels=[1.0],
        #            colors='yellow', linewidths=1)
        # ax.contour(fraction_x, fraction_y, fraction_mask, levels=[0.01], colors='yellow', linewidths=0.5)

        if row==0:
            ax.set_title('Landsat 8\nReference Imagery', fontsize=12)

        ax.set_xticks([])
        ax.set_yticks([])

        ######################################################################
        # Melange Fraction
        ax = fig.add_subplot(gs[plot_height*row:plot_height*(row+1), 1])

        ax.imshow(image, extent=(image_x.min(), image_x.max(), image_y.min(), image_y.max()),
                  origin='lower', alpha=0.5)

        # add a white mask on the image to mask out the ocean
        ocean_mask = np.ma.masked_where(fraction_mask!=0, np.ones_like(fraction_mask))
        ax.imshow(ocean_mask, extent=(fraction_x.min(), fraction_x.max(), fraction_y.min(), fraction_y.max()),
                  origin='lower', cmap='gray', vmin=0, vmax=1)

        fraction = np.ma.masked_where(fraction_mask!=0, fraction)
        fraction = np.ma.masked_where(fraction<0.5, fraction)
        # ax.pcolormesh(fraction_x, fraction_y, fraction, cmap='turbo', vmin=fraction_vmin, vmax=fraction_vmax)
        ax.imshow(fraction, extent=(fraction_x.min(), fraction_x.max(), fraction_y.min(), fraction_y.max()),
                  origin='lower', cmap='turbo', vmin=fraction_vmin, vmax=fraction_vmax)

        ax.contour(fraction_x, fraction_y, fraction_mask, levels=[0.01], colors='k', linewidths=0.5)

        # make a scale bar that is 10 km long in the upper right
        if locations[row]=='Upernavik Fjord':
            length = 10000
        if locations[row]=='Kakivfaat':
            length = 5000
        if locations[row]=='Ussing Braeer':
            length = 5000
        scale_y = (fraction_y.max() - fraction_y.min()) * 0.85 + fraction_y.min()
        h_spacing = (fraction_x.max() - fraction_x.min()) * 0.05
        scale_min_x = fraction_x.max() - h_spacing - length
        scale_max_x = fraction_x.max() - h_spacing
        ax.plot([scale_min_x, scale_max_x], [scale_y, scale_y], color='k', linewidth=2)
        ax.text((scale_min_x + scale_max_x)/2, scale_y + (fraction_y.max() - fraction_y.min()) * 0.02,
                str(int(length/1000))+' km',
                ha='center', va='bottom', fontsize=12)

        if row==0:
            ax.set_title('Melange Fraction', fontsize=12)

        ax.set_xticks([])
        ax.set_yticks([])

        ######################################################################
        # Melange Rigidity
        ax = fig.add_subplot(gs[plot_height*row:plot_height*(row+1), 2])

        ax.imshow(image, extent=(image_x.min(), image_x.max(), image_y.min(), image_y.max()),
                  origin='lower', alpha=0.5)

        # add a white mask on the image to mask out the ocean
        # ocean_mask = np.ma.masked_where(fraction_mask != 0, np.ones_like(fraction_mask))
        ax.imshow(ocean_mask, extent=(fraction_x.min(), fraction_x.max(), fraction_y.min(), fraction_y.max()),
                  origin='lower', cmap='gray', vmin=0, vmax=1)

        norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], rigid_cmap.N)

        rigidity_plot = np.zeros_like(rigidity)
        rigidity_plot[rigidity==0] = 1
        rigidity_plot[rigidity==1] = 2
        rigidity_plot[fraction_mask!=0] = 0
        rigidity_plot[fraction<0.5] = 0
        rigidity_plot = np.ma.masked_where(rigidity_plot==0, rigidity_plot)

        plt.imshow(rigidity_plot, norm=norm, cmap=rigid_cmap, extent=(rigidity_x.min(), rigidity_x.max(), rigidity_y.min(), rigidity_y.max()),
                   origin='lower')

        # ax.pcolormesh(rigidity_x, rigidity_y, rigidity, cmap='Blues', vmin=0, vmax=1)
        ax.contour(fraction_x, fraction_y, fraction_mask, levels=[0.01], colors='k', linewidths=0.5)

        if row==0:
            ax.set_title('Melange Rigidity', fontsize=12)

        ax.set_xticks([])
        ax.set_yticks([])

    #################################################################
    # Colorbars
    # ax = fig.add_subplot(gs[-colorbar_height:, 0])
    # cx = np.linspace(0, 1, 100)
    # cy = np.array([0,1])
    # cX, cY = np.meshgrid(cx, cy)
    # c = ax.pcolormesh(cX, cY, cX, cmap='gray', vmin=0, vmax=1)
    # ax.set_xlabel('Relative Pixel Brightness')
    # ax.set_yticks([])

    ax = fig.add_subplot(gs[-colorbar_height:, 1])
    cx = np.linspace(fraction_vmin, fraction_vmax, 100)
    cy = np.array([0, 1])
    cX, cY = np.meshgrid(cx, cy)
    c = ax.pcolormesh(cX, cY, cX, cmap=fraction_cmap, vmin=fraction_vmin, vmax=fraction_vmax)
    ax.set_xlabel('Melange Fraction')
    ax.set_yticks([])

    ax = fig.add_subplot(gs[-colorbar_height:, 2])
    cx = np.linspace(0, 2, 3)
    cy = np.array([0, 1])
    cX, cY = np.meshgrid(cx, cy)
    c = ax.pcolormesh(cX, cY, cX, cmap=rigid_cmap)
    # ax.set_xlabel('Rigidity Mask')
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(['Not\nMelange', 'Drifting\nMelange', 'Rigid\nMelange'])
    ax.set_yticks([])

    output_file = os.path.join(project_dir, 'Figures','Ocean','L2_Upernavik_Melange_Imagery_and_Masks.png')
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

locations = ['Ussing Braeer', 'Kakivfaat', 'Upernavik Fjord']

location_dict = {'Ussing Braeer':{}, 'Kakivfaat':{}, 'Upernavik Fjord':{}}

#######################
# Read in imagery

for location in locations:
    location_dict = read_melange_image(project_dir, location_dict, location)
    location_dict = read_melange_mask(project_dir, location_dict, location)
    location_dict = read_rigidity_mask(project_dir, location_dict, location)


plot_melange_mask(project_dir, locations, location_dict)







