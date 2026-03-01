

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shapefile
from pyproj import Transformer
from scipy.interpolate import griddata

def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)

def read_transect_from_shapefile(project_dir, fjord_name):
    file_name = project_dir+'/Map/Shapefiles/'+'_'.join(fjord_name.split())+'_short_with_ice'
    sf = shapefile.Reader(file_name)
    points_4326 = np.array(sf.shapes()[0].points)
    points_4326 = series_to_N_points(points_4326,1000)

    transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
    points_x, points_y = transformer.transform(points_4326[:, 1].ravel(), points_4326[:, 0].ravel())
    points = np.column_stack([points_x, points_y])

    distance = np.zeros((np.shape(points)[0],))
    for d in range(len(distance)-1):
        distance[d+1] = distance[d]+ ((points[d+1,0]-points[d,0])**2 + (points[d+1,1]-points[d,1])**2)**0.5

    return(points, points_4326, distance)

def sample_bedmachine_on_shapefile(distance, transect):
    bedmachine_file = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Bathymetry/BedMachineGreenland-v5.nc'
    ds = nc4.Dataset(bedmachine_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]

    min_x_index = np.argmin(np.abs(x - np.min(transect[:, 0])))
    max_x_index = np.argmin(np.abs(x - np.max(transect[:, 0])))
    max_y_index = np.argmin(np.abs(y - np.min(transect[:, 1])))
    min_y_index = np.argmin(np.abs(y - np.max(transect[:, 1])))

    buffer=5
    x = x[min_x_index-buffer:max_x_index+buffer]
    y = y[min_y_index - buffer:max_y_index + buffer]

    surface = ds.variables['surface'][:,:]
    surface = surface[min_y_index - buffer:max_y_index + buffer, min_x_index-buffer:max_x_index+buffer]

    bed = ds.variables['bed'][:, :]
    bed = bed[min_y_index - buffer:max_y_index + buffer, min_x_index - buffer:max_x_index + buffer]

    ds.close()


    X, Y = np.meshgrid(x,y)

    # C = plt.pcolormesh(X, Y, bed)
    # plt.colorbar(C)
    # plt.show()

    surface_transect = griddata(np.column_stack([X.ravel(), Y.ravel()]), surface.ravel(), (transect[:,0], transect[:,1]))

    X, Y = np.meshgrid(x, y)
    bed_transect = griddata(np.column_stack([X.ravel(), Y.ravel()]), bed.ravel(),
                                (transect[:, 0], transect[:, 1]))

    ice_base = -surface_transect*917/(1024-917)
    for i in range(len(ice_base)):
        if ice_base[i]<bed_transect[i]:
            ice_base[i] = bed_transect[i]

    # plt.plot(distance, surface_transect)
    # plt.plot(distance, ice_base)
    # plt.plot(distance, bed_transect)
    # plt.show()

    return(surface_transect, ice_base, bed_transect)

def create_plotting_polygons(distance, surface, ice_base, bed):

    bathy_top = np.column_stack([distance, bed])
    bathy_bottom = np.flipud(np.copy(bathy_top))
    bathy_bottom[:, 1] = np.min(bed)-100
    bathy_outline = np.vstack([bathy_top, bathy_bottom])

    ice_top = np.column_stack([distance, surface])
    ice_bottom = np.column_stack([distance, ice_base])
    ice_indices = surface>0
    ice_top = ice_top[ice_indices,:]
    ice_bottom = ice_bottom[ice_indices,:]
    ice_bottom = np.flipud(ice_bottom)
    ice_outline = np.vstack([ice_top, ice_bottom])

    # plt.plot(bathy_outline[:,0], bathy_outline[:,1], 'k-', linewidth=3)
    # plt.plot(ice_outline[:,0], ice_outline[:,1],'b-')
    # plt.show()

    return(bathy_outline, ice_outline)

def save_transect_as_nc(project_dir, fjord_name, transect, transect_4326, distance,
                        surface, ice_base, bed,
                        bathy_outline, ice_outline):

    ds = nc4.Dataset(os.path.join(project_dir,'Data','Models','Fjord Transects',fjord_name,
                                  fjord_name+' Fjord Geometry Transect.nc'),'w')

    ds.createDimension('distance',len(distance))
    ds.createDimension('bed_outline_len',np.shape(bathy_outline)[0])
    ds.createDimension('ice_outline_len', np.shape(ice_outline)[0])

    t = ds.createVariable('transect_x','f4',('distance',))
    t[:] = transect[:,0]
    t = ds.createVariable('transect_y', 'f4', ('distance',))
    t[:] = transect[:, 1]

    t = ds.createVariable('transect_lon', 'f4', ('distance',))
    t[:] = transect_4326[:, 0]
    t = ds.createVariable('transect_lat', 'f4', ('distance',))
    t[:] = transect_4326[:, 1]

    t = ds.createVariable('distance', 'f4', ('distance',))
    t[:] = distance

    t = ds.createVariable('bed', 'f4', ('distance',))
    t[:] = bed

    t = ds.createVariable('surface', 'f4', ('distance',))
    t[:] = surface

    t = ds.createVariable('ice_base', 'f4', ('distance',))
    t[:] = ice_base

    t = ds.createVariable('bathy_polygon_x', 'f4', ('bed_outline_len',))
    t[:] = bathy_outline[:, 0]
    t = ds.createVariable('bathy_polygon_y', 'f4', ('bed_outline_len',))
    t[:] = bathy_outline[:, 1]

    t = ds.createVariable('ice_polygon_x', 'f4', ('ice_outline_len',))
    t[:] = ice_outline[:, 0]
    t = ds.createVariable('ice_polygon_y', 'f4', ('ice_outline_len',))
    t[:] = ice_outline[:, 1]

    ds.close()


project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

for fjord_name in ['Upernavik N','Kakivfaat','Ussing Braeer']:

    transect, transect_4326, distance = read_transect_from_shapefile(project_dir, fjord_name)

    surface, ice_base, bed = sample_bedmachine_on_shapefile(distance, transect)

    bathy_outline, ice_outline = create_plotting_polygons(distance, surface, ice_base, bed)

    save_transect_as_nc(project_dir, fjord_name, transect, transect_4326, distance,
                            surface, ice_base, bed,
                            bathy_outline, ice_outline)





