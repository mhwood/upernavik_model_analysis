import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import shapefile
import sys
from bisect import bisect_left
from matplotlib.path import Path
from scipy.interpolate import griddata

def read_grid_geometry_from_nc(config_dir, model_name):
    nc_file = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    dXC = ds.variables['dxC'][:, :]
    dYC = ds.variables['dyC'][:, :]
    rA = ds.variables['rA'][:, :]
    Depth = ds.variables['Depth'][:, :]
    hFacC = np.array(ds.variables['HFacC'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return (XC, YC, dXC, dYC, Depth, rA, hFacC, delR)

def read_melange_polygons_from_shapefile(config_dir, model_name):
    shapefile_path = os.path.join(config_dir, 'L2', model_name, 'input', 'melange_polygons', model_name+'_melange_polygons.shp')
    sf = shapefile.Reader(shapefile_path)
    shapes = sf.shapes()

    melange_polygons = {}
    for shape in shapes:
        points = np.array(shape.points)
        fjord = sf.record(shapes.index(shape))[0]
        melange_polygons[fjord] = points

    return melange_polygons

def read_melange_fraction_and_rigidiy(project_dir, region):
    file_path = os.path.join(project_dir, 'Data', 'Observations','Melange', region+'_Landsat_8_Melange_Fraction.nc')
    ds = nc4.Dataset(file_path, 'r')
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    melange_fraction = ds.variables['melange_fraction'][:, :]
    ds.close()
    melange_fraction[melange_fraction < 0.5] = 0

    file_path = os.path.join(project_dir, 'Data', 'Observations','Melange', region+'_Melange_Rigidity_Mask.nc')
    ds = nc4.Dataset(file_path, 'r')
    melange_rigidity = ds.variables['rigidity_mask_fraction'][:, :]
    ds.close()

    threshold=0.5
    melange_rigidity[melange_rigidity < threshold] = 0
    melange_rigidity[melange_rigidity >= threshold] = 1
    melange_rigidity[melange_fraction==0] = 0

    return melange_fraction, melange_rigidity, x, y

def create_masks(XC, YC, Depth, hFacC, # melange_polygons,
                 ussing_braeer_melange_fraction, ussing_braeer_melange_rigidity,
                 ussing_braeer_x_melange, ussing_braeer_y_melange,
                 kakivfaat_melange_fraction, kakivfaat_melange_rigidity,
                 kakivfaat_x_melange, kakivfaat_y_melange,
                 upernavik_melange_fraction, upernavik_melange_rigidity,
                 upernavik_x_melange, upernavik_y_melange):

    points_4326 = np.column_stack((XC.flatten(), YC.flatten()))
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    points_3413 = np.array(transformer.transform(points_4326[:, 0], points_4326[:, 1]))
    X_3413 = points_3413[0, :].reshape(XC.shape)
    Y_3413 = points_3413[1, :].reshape(YC.shape)

    bergMask = np.zeros_like(np.array(Depth))
    bergConc = np.zeros_like(np.array(Depth))
    driftMask = np.zeros_like(np.array(Depth))
    meltMask = np.zeros_like(np.array(Depth))
    barrierMask = np.zeros_like(np.array(Depth))

    # Iceberg mask
    for location in ['Upernavik_Fjord','Kakivfaat','Ussing_Braeer']:

        if location == 'Upernavik_Fjord':
            melange_fraction = upernavik_melange_fraction
            rigidity = upernavik_melange_rigidity
            x_melange = upernavik_x_melange
            y_melange = upernavik_y_melange
        elif location == 'Kakivfaat':
            melange_fraction = kakivfaat_melange_fraction
            rigidity = kakivfaat_melange_rigidity
            x_melange = kakivfaat_x_melange
            y_melange = kakivfaat_y_melange
        elif location == 'Ussing_Braeer':
            melange_fraction = ussing_braeer_melange_fraction
            rigidity = ussing_braeer_melange_rigidity
            x_melange = ussing_braeer_x_melange
            y_melange = ussing_braeer_y_melange

        X_melange, Y_melange = np.meshgrid(x_melange, y_melange)
        fraction_on_domain = griddata((X_melange.flatten(), Y_melange.flatten()), melange_fraction.flatten(),
                                      (X_3413.flatten(), Y_3413.flatten()), method='linear').reshape(XC.shape)
        rigidity_on_domain = griddata((X_melange.flatten(), Y_melange.flatten()), rigidity.flatten(),
                                        (X_3413.flatten(), Y_3413.flatten()), method='linear').reshape(XC.shape)
        rigidity_on_domain[rigidity_on_domain<0.5] = 0
        rigidity_on_domain[rigidity_on_domain>=0.5] = 1

        points = np.vstack((X_3413.flatten(), Y_3413.flatten())).T

        # path = Path(polygon)
        # inside = path.contains_points(points)
        # inside_mask = inside.reshape(np.shape(XC))

        inside_mask = XC>-1e22

        fraction_mask = fraction_on_domain >= 0.5

        bergMask[inside_mask & fraction_mask] = 1
        bergConc[inside_mask & fraction_mask] = fraction_on_domain[inside_mask & fraction_mask]

        barrierMask[inside_mask & fraction_mask] = rigidity_on_domain[inside_mask & fraction_mask]
        driftMask[inside_mask & fraction_mask] = 1-rigidity_on_domain[inside_mask & fraction_mask]

        # # transition berg conc from min to max over the fjord length
        # min_col, max_col = np.min(np.where(inside_mask)[1]), np.max(np.where(inside_mask)[1])
        # fjord_length = max_col - min_col
        # for col in range(min_col, max_col + 1):
        #     col_mask = inside_mask[:,col]
        #     if np.any(col_mask):
        #         conc = min_conc + (max_conc - min_conc) * (col - min_col) / fjord_length
        #         bergConc[col_mask,col] = conc

    # apply land and shallow masks
    bergMask[hFacC[0, :, :] == 0] = 0  # Mask out land points
    bergMask[Depth < 50] = 0  # Mask out shallow points
    bergConc[hFacC[0, :, :] == 0] = 0  # Mask out land points
    bergConc[Depth < 50] = 0  # Mask out shallow points

    barrierMask[hFacC[0, :, :] == 0] = 0  # Mask out land points
    barrierMask[Depth < 50] = 0  # Mask out shallow points
    driftMask[hFacC[0, :, :] == 0] = 0  # Mask out land points
    driftMask[Depth < 50] = 0  # Mask out shallow points

    # # Drift mask
    # driftMask[np.logical_and(bergConc > 0, bergConc <= 80)] = 1  # Low conc bergs can drift
    #
    # # Barrier mask
    # barrierMask[bergConc >= 80] = 1  # High conc bergs act as barriers

    # # Drift mask
    print('Only drifting in these files for testing')
    driftMask[bergConc > 0] = 1  # All bergs drift
    barrierMask[bergConc > 0] = 0  # No bergs act as barriers, for testing

    # Melt mask, only let bergs melt in this region (make melt water, these don't change size)
    meltMask[bergMask==1] = 1  # Allow focus on blocking effect only

    return bergMask, bergConc, driftMask, meltMask, barrierMask

def find_closest_indices(sorted_A, sorted_B):
    closest_indices = []
    for a in sorted_A:
        pos = bisect_left(sorted_B, a)  # Find position in B where a would fit
        # Compare neighbors to find the closest
        if pos == 0:
            closest_indices.append(0)
        elif pos == len(sorted_B):
            closest_indices.append(len(sorted_B) - 1)
        else:
            before = pos - 1
            after = pos
            closest_indices.append(before if abs(sorted_B[before] - a) <= abs(sorted_B[after] - a) else after)
    return closest_indices

def populate_masks_with_icebergs(bergMask, dXC, dYC, hFacC,rA, delR):

    # putting these manually for now but should be done
    # with dXC and dYC
    deltaX = 500
    deltaY = 500

    # This function is from Paul Summers

    iceBergDepth = 400  # max iceberg depth [meters], used for ICEBERG package
    iceCoverage = 60  # % of ice cover in melange, stay under 90% ideally
    bergType = 1  # 1 = block 2 = cone (not implemented)
    alpha = 1.9 * 2  # slope of inverse power law size frequency distribution
    scaling = 1  # 1 = Sulak 2017 2 = Barker 2004
    maxBergDepth = iceBergDepth  # (m) - set to zero if 'prescribing' max iceberg width, set at top here
    minBergDepth = 40  # (m)
    maxBergWidth = 0  # (m) - set to zero if 'prescribing' max iceberg depth
    minBergWidth = 40  # (m)

    hfacThreshold = .95

    # Iceberg concentration (# of each surface cell that is filled in plan view)
    bergConc = np.zeros_like(bergMask)
    bergConc[bergMask == 1] = iceCoverage

    desiredBergArea = np.sum(bergConc / 100.0 * rA)
    bergMaskArea = np.sum(bergMask * rA)
    print('Area where bergs live: ' + str(bergMaskArea) + ' m^2')
    print('Desired berg area: ' + str(desiredBergArea) + ' m^2')
    print('Ratio: ' + str(desiredBergArea / bergMaskArea * 100) + '%')

    if (scaling == 1):  # then use Sulak17 volume-area scaling volume = 6.0*area^1.30
        # assumes volume = L*W*D and W = L/1.62 (Dowdeswell et al 1992)
        if (maxBergWidth == 0):
            maxBergWidth = 0.0642449 * maxBergDepth ** (5 / 3)
            # minBergWidth = 0.0642449*minBergDepth**(5/3)
        elif (maxBergDepth == 0):
            maxBergDepth = 5.19155 * maxBergWidth ** (5 / 3)
            minBergDepth = 5.19155 * minBergWidth ** (5 / 3)
    elif (scaling == 2):  # Then use Barker04 width-depth relationship
        # Depth = 2.91*Width^0.71
        if (maxBergWidth == 0):
            maxBergWidth = (100 * 10 ** (58 / 71) * maxBergDepth ** (100 / 71)) / (291 * 291 ** (29 / 71))
            # minBergWidth = (100*10**(58/71)*minBergDepth**(100/71)) / (291*291**(29/71))
        elif (maxBergDepth == 0):
            maxBergDepth = 2.91 * maxBergWidth ^ 0.71
            minBergDepth = 2.91 * minBergWidth ^ 0.71

    numberOfBergs = 50  # low start, immediately doubled by scheme below, so guess low, high guesses (300%+) can cause to fail
    bergTopArea = 0
    areaResidual = 1
    # Generate the Inverse Power Law cumulative distribution function
    # over the range minBergWidth-maxBergWidth with a slope of alpha.
    print('Making bergs, this can take a few loops...')
    loop_count = 1

    np.random.seed(3)
    print('random seed set, not really random anymore')

    while (np.abs(areaResidual) > .005):  # Create random power dist of bergs, ensure correct surface area
        numberOfBergs = round(numberOfBergs * (1 + areaResidual))
        print('\tnumberOfBergs: ' + str(numberOfBergs))
        x_width = np.arange(minBergWidth, maxBergWidth, (maxBergWidth - minBergWidth) / (numberOfBergs * 1e2))
        inversePowerLawPDF_width = ((alpha - 1) / minBergWidth) * (x_width / minBergWidth) ** (-alpha)
        # Get the CDF numerically
        inversePowerLawCDF_width = np.cumsum(inversePowerLawPDF_width)
        # Normalize
        inversePowerLawCDF_width = inversePowerLawCDF_width / inversePowerLawCDF_width[-1]

        # Generate number_of_bergs uniformly distributed random numbers.
        uniformlyDistributedRandomNumbers = np.random.uniform(0, 1, numberOfBergs)

        inversePowerLawDistNumbers_width = np.zeros(uniformlyDistributedRandomNumbers.size)

        nearestIndex_width = [0] * uniformlyDistributedRandomNumbers.size

        nearestIndex_width = find_closest_indices(uniformlyDistributedRandomNumbers, inversePowerLawCDF_width)

        inversePowerLawDistNumbers_width = x_width[nearestIndex_width]
        inversePowerLawDistNumbers_length = inversePowerLawDistNumbers_width / 1.12  # Widths are bigger

        randScale = np.random.normal(6, 1.22, numberOfBergs)
        randPower = np.random.normal(0.3, 0.016, numberOfBergs)
        inversePowerLawDistNumbers_depth = randScale * (
                inversePowerLawDistNumbers_width * inversePowerLawDistNumbers_length) ** randPower * (920 / 1025)
        inversePowerLawDistNumbers_depth[
            inversePowerLawDistNumbers_depth < 5.123] = 5.123  # doest like round numbers, cap low end of bergs

        tooWide = np.count_nonzero(
            inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold))  # disallow completely full cells
        tooLong = np.count_nonzero(inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold))
        inversePowerLawDistNumbers_width[
            inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(
            hfacThreshold) - .01  # Max width is grid cell (assumed square)
        inversePowerLawDistNumbers_length[
            inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(
            hfacThreshold) - .01  # Max length is grid cell (assumed square)
        if (tooLong + tooWide > 0):
            print('\t\tBergs clipped: %i for width, %i for length' % (tooWide, tooLong))

        bergTopArea = sum(inversePowerLawDistNumbers_width * inversePowerLawDistNumbers_length)
        areaResidual = (desiredBergArea - bergTopArea) / desiredBergArea
        print('\t\t%.2f %% Bergs' % (bergTopArea / bergMaskArea * 100))
        print('\t\tareaResidual %.2f %%' % (areaResidual * 100))
        loop_count += 1
    print('====== Success! Found our bergs =====')
    print('Width min/mean/max: %f/%f/%f [m]' % (
        np.min(inversePowerLawDistNumbers_width), np.mean(inversePowerLawDistNumbers_width),
        np.max(inversePowerLawDistNumbers_width)))
    print('Depth min/mean/max: %f/%f/%f [m]' % (
        np.min(inversePowerLawDistNumbers_depth), np.mean(inversePowerLawDistNumbers_depth),
        np.max(inversePowerLawDistNumbers_depth)))
    print('Total Berg Area %f' % bergTopArea)
    print('Total Berg fract: %.2f %%' % (bergTopArea / bergMaskArea * 100))

    print('Shape of berg arrays: ' + str(np.shape(inversePowerLawDistNumbers_depth)))

    # Now we sort these berg into cell, randomly
    bergMaski = 0  # Bad name, but this is the count of cells that will recieve bergs
    bergDict = {}

    bergMaskNums = np.zeros_like(bergMask)  # This assigns each cell a unique number for indexing
    for j in range(np.shape(bergMask)[0]):
        for i in range(np.shape(bergMask)[1]):
            if (bergMask[j, i] == 1):
                # print('i,j, bergmask',i,j,bergMask[j,i])
                bergMaski = 1 + bergMaski  # Needs to start at 1, as non-bergs will be 0
                bergMaskNums[j, i] = bergMaski  # Assign Mask Nums, not random as we'll randomly place bergs in cells
                bergDict[bergMaski] = [j, i]  # This lets us do 1-D loops for the whole grid
    print('%i cells with bergs' % bergMaski)

    # Sort my bergs
    sorted_indices = np.argsort(
        -inversePowerLawDistNumbers_depth)  # Sort backwards to get descending from big to small bergs
    sorted_depth = inversePowerLawDistNumbers_depth[sorted_indices]
    sorted_width = inversePowerLawDistNumbers_width[sorted_indices]
    sorted_length = inversePowerLawDistNumbers_length[sorted_indices]
    assignedCell = np.random.randint(0, bergMaski, [numberOfBergs])  # In this script, every berg has a home

    # Array for bergs
    bergsPerCellLimit = 500
    icebergs_depths = np.zeros([bergMaski, bergsPerCellLimit])
    icebergs_widths = np.zeros([bergMaski, bergsPerCellLimit])
    icebergs_length = np.zeros([bergMaski, bergsPerCellLimit])  # careful, not plural as to length match

    np.random.seed(2)
    assignedCell = np.random.randint(0, bergMaski, [numberOfBergs])  # every Berg has a spot

    icebergs_per_cell = np.zeros([bergMaski], dtype=np.int16)
    icebergs_area_per_cell = np.zeros([bergMaski])

    for i in range(numberOfBergs):
        j = assignedCell[i]
        # print('looking at mask number',j,'at berg',i)
        # print('Berg number', icebergs_per_cell[j],'in this cell')
        # print('Width, Length',sorted_width[i], sorted_length[i])
        bergArea = sorted_width[i] * sorted_length[i]
        loopLimiter = 0
        while (bergArea > (
                deltaX * deltaY * (bergConc[bergDict[j + 1][0], bergDict[j + 1][1]] / 100) - icebergs_area_per_cell[
            j])):  # if above 'full', pick random new cell, accept with decreasing probability
            j_old = j
            j = np.random.randint(0, bergMaski)
            loopLimiter += 1
            if ((bergArea + icebergs_area_per_cell[j]) / (
                    deltaX * deltaY) < hfacThreshold - .01):  # only consider accepting if under 95
                odds = np.abs(np.random.normal(0, .5,
                                               1))  # randomly accepts those that are big in overfull cells, but at decreasing frequency
                overFull = ((bergArea + icebergs_area_per_cell[j]) / (deltaX * deltaY) * 100 - bergConc[
                    bergDict[j + 1][0], bergDict[j + 1][1]])
                if (odds > overFull):
                    # print('accepting overfull')
                    assignedCell[i] = j  # if we it a shuffling critera
                    break
            if (loopLimiter > bergMaski * 20):  # eventually we have to force some in
                indexesAllowed = np.where((deltaX * deltaY * hfacThreshold - icebergs_area_per_cell) > bergArea)
                # print(indexesAllowed)
                if (len(indexesAllowed[0]) == 0):
                    print('ERROR: no cells with room for this berg, reduce size or increase hfacThreshold')
                    sys.exit()
                randi = np.random.randint(0, len(indexesAllowed[0]))
                j = indexesAllowed[0][randi]
                assignedCell[i] = j  # if we it a shuffling critera, must line up for calculation below
                # print('\t Randomly missed, will force into cell with room: %i' % j)
                if ((np.min(icebergs_area_per_cell) + bergArea) / (deltaX * deltaY) > hfacThreshold):
                    print('WARNING cell very full: %.2f%%' % (
                            (np.min(icebergs_area_per_cell) + bergArea) * 100 / (deltaX * deltaY)))
                break

        icebergs_depths[j, icebergs_per_cell[j]] = sorted_depth[i]
        icebergs_widths[j, icebergs_per_cell[j]] = sorted_width[i]
        icebergs_length[j, icebergs_per_cell[j]] = sorted_length[i]
        icebergs_per_cell[j] += 1
        # icebergs_area_per_cell[j] = np.sum(icebergs_widths[j,:]*icebergs_length[j,:])
        icebergs_area_per_cell[j] += bergArea
    print('Bergs per cell and filled faction at surface for spot check')
    print(icebergs_per_cell)
    print(np.round(icebergs_area_per_cell / (deltaX * deltaY), 2))
    print('Max fill is: %.2f%%' % (np.nanmax(icebergs_area_per_cell / (deltaX * deltaY)) * 100))

    # All bergs now sorted
    openFrac = np.zeros_like(np.array(hFacC))
    SA = np.zeros_like(np.array(hFacC))
    SA[:, :, :] = np.nan

    # This loop knows where all bergs are already, different from searching for all bergs across entire grid
    numBergsPerCell = np.zeros(np.shape(rA), dtype=np.int64)
    sum_z = np.cumsum(delR)
    for i in range(bergMaski):
        bergCount = icebergs_per_cell[i]
        numBergsPerCell[bergDict[i + 1][0], bergDict[i + 1][1]] = bergCount
        if (bergCount > 0):
            lengths = icebergs_length[i, icebergs_length[i, :] > 0]  # return only non-zeros
            widths = icebergs_widths[i, icebergs_widths[i, :] > 0]  # return only non-zeros
            depths = icebergs_depths[i, icebergs_depths[i, :] > 0]  # return only non-zeros
            for k in range(np.shape(hFacC)[0]):
                cellVolume = deltaX * deltaY * delR[k]
                d_bot = sum_z[k]  # bottom of depth bin
                d_top = sum_z[k] - delR[k]
                volume1 = delR[k] * lengths[depths > d_bot] * widths[depths > d_bot]
                SA1 = delR[k] * 2 * (lengths[depths > d_bot] + widths[depths > d_bot])
                partialFill = (depths < d_bot) & (depths > d_top)
                # partial fill
                volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
                # partial sides
                SA2 = (depths[partialFill] - d_top) * 2 * (lengths[partialFill] + widths[partialFill])
                # bottom
                SA3 = lengths[partialFill] * widths[partialFill]
                # print(np.sum(volume1), np.sum(volume2))
                openFrac[k, bergDict[i + 1][0], bergDict[i + 1][1]] = 1 - (
                            (np.sum(volume1) + np.sum(volume2)) / cellVolume)
                SA[k, bergDict[i + 1][0], bergDict[i + 1][1]] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
        elif (bergCount == 0):
            openFrac[:, bergDict[i + 1][0], bergDict[i + 1][1]] = 1
            SA[:, bergDict[i + 1][0], bergDict[i + 1][1]] = 0

    icebergs_depths2D = np.zeros([bergsPerCellLimit, np.shape(rA)[0], np.shape(rA)[1]])
    icebergs_widths2D = np.zeros([bergsPerCellLimit, np.shape(rA)[0], np.shape(rA)[1]])
    icebergs_length2D = np.zeros([bergsPerCellLimit, np.shape(rA)[0], np.shape(rA)[1]])

    for k in range(bergMaski):
        j = bergDict[k + 1][0]
        i = bergDict[k + 1][1]
        icebergs_depths2D[:, j, i] = icebergs_depths[k, :]
        icebergs_widths2D[:, j, i] = icebergs_widths[k, :]
        icebergs_length2D[:, j, i] = icebergs_length[k, :]

    # # if(forceDraft): #this lets you resize all the icebergs to fit a desired mélange draft
    # #     mX=np.load(run_config['run_dir']+'/input/melangeX.npy')
    # #     mH=np.load(run_config['run_dir']+'/input/melangeH.npy')
    # #     x_fd = np.arange(deltaX/2,deltaX*(nx+.5),deltaX)
    # #     depthHelper = icebergs_depths2D.copy()
    # #     depthHelper[depthHelper == 0] = np.nan
    # #     depthMedian = np.nanmedian(depthHelper,axis=[0,1])
    # #     del depthHelper
    # #     for i in range(nx):
    # #         for j in range(ny):
    # #             if(~np.isnan(depthMedian[i])):
    # #                 icebergs_depths2D[:,j,i] = icebergs_depths2D[:,j,i]/depthMedian[i]*np.interp(x,mX,mH,0,0)[i]
    # #                 # pass
    #
    # #     for i in range(nx):
    # #         for j in range(ny):
    # #             numberOfBergs = len(icebergs_length2D[icebergs_length2D[:,j,i] > 0])
    # #             if(numberOfBergs > 0):
    # #                 lengths = icebergs_length2D[icebergs_length2D[:,j,i] > 0,j,i] #return only non-zeros
    # #                 widths = icebergs_widths2D[icebergs_widths2D[:,j,i] > 0,j,i] #return only non-zeros
    # #                 depths = icebergs_depths2D[icebergs_depths2D[:,j,i] > 0,j,i] #return only non-zeros
    # #                 for k in range(nz):
    # #                     cellVolume = deltaX*deltaY*dz[k]
    # #                     d_bot = sum_z[k] #bottom of depth bin
    # #                     d_top = sum_z[k] - dz[k]
    # #                     volume1 = dz[k] * lengths[depths > d_bot] * widths[depths > d_bot]
    # #                     SA1 = dz[k]*2*(lengths[depths > d_bot] + widths[depths > d_bot])
    # #                     partialFill = (depths < d_bot) & (depths > d_top)
    # #                     #partial fill
    # #                     volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
    # #                     #partial sides
    # #                     SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill])
    # #                     #bottom
    # #                     SA3 = lengths[partialFill] * widths[partialFill]
    # #                     #print(np.sum(volume1), np.sum(volume2))
    # #                     openFrac[k,j,i] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
    # #                     SA[k,j,i] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
    # #             elif(bergCount == 0):
    # #                 openFrac[:,j,i] = 1
    # #                 SA[:,j,i] = 0

    return(bergMaskNums, numBergsPerCell, openFrac, SA, bergConc,
           icebergs_depths2D, icebergs_widths2D, icebergs_length2D,
           inversePowerLawDistNumbers_depth, inversePowerLawDistNumbers_width, inversePowerLawDistNumbers_length)

def create_L2_iceberg_files(config_dir, L2_model_name, print_level):
    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils', 'init_file_creation'))

    XC, YC, dXC, dYC, Depth, rA, hFacC, delR = read_grid_geometry_from_nc(config_dir, L2_model_name)

    # melange_polygons = read_melange_polygons_from_shapefile(config_dir, L2_model_name)

    ussing_braeer_melange_fraction, ussing_braeer_melange_rigidity, ussing_braeer_x_melange, ussing_braeer_y_melange =\
        read_melange_fraction_and_rigidiy(project_dir, 'Ussing_Braeer')

    kakivfaat_melange_fraction, kakivfaat_melange_rigidity, kakivfaat_x_melange, kakivfaat_y_melange =\
        read_melange_fraction_and_rigidiy(project_dir, 'Kakivfaat')

    upernavik_melange_fraction, upernavik_melange_rigidity, upernavik_x_melange, upernavik_y_melange =\
        read_melange_fraction_and_rigidiy(project_dir, 'Upernavik_Fjord')

    iceberg_mask, iceberg_conc, drift_mask, melt_mask, barrier_mask = create_masks(XC, YC, Depth, hFacC, # melange_polygons,
                                                                                   ussing_braeer_melange_fraction, ussing_braeer_melange_rigidity,
                                                                                   ussing_braeer_x_melange, ussing_braeer_y_melange,
                                                                                   kakivfaat_melange_fraction, kakivfaat_melange_rigidity,
                                                                                   kakivfaat_x_melange, kakivfaat_y_melange,
                                                                                   upernavik_melange_fraction, upernavik_melange_rigidity,
                                                                                   upernavik_x_melange, upernavik_y_melange)

    precision = '>f4'

    # rows = np.arange(375)
    # cols = np.arange(450)
    # Cols, Rows = np.meshgrid(cols, rows)
    # plt.subplot(1,2,1)
    # plt.pcolormesh(Cols, Rows, iceberg_mask)
    # plt.contour(Cols, Rows, Depth, levels=[1], colors='k')
    # plt.title('Iceberg Mask')
    # plt.subplot(1,2,2)
    # C = plt.pcolormesh(Cols, Rows, iceberg_conc)
    # plt.contour(Cols, Rows, Depth, levels=[1], colors='k')
    # plt.title('Iceberg Concentration (%)')
    # plt.colorbar(C)
    # plt.show()

    iceberg_mask_file = os.path.join(config_dir,'L2',L2_model_name,'input','iceberg','bergMask.bin')
    iceberg_mask.ravel('C').astype(precision).tofile(iceberg_mask_file)
    iceberg_conc_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'iceberg', 'bergConc.bin')
    iceberg_conc.ravel('C').astype(precision).tofile(iceberg_conc_file)
    drift_mask_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'iceberg', 'driftMask.bin')
    drift_mask.ravel('C').astype(precision).tofile(drift_mask_file)
    melt_mask_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'iceberg', 'meltMask.bin')
    melt_mask.ravel('C').astype(precision).tofile(melt_mask_file)
    barrier_mask_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'iceberg', 'barrierMask.bin')
    barrier_mask.ravel('C').astype(precision).tofile(barrier_mask_file)

    bergMaskNums, numBergsPerCell, openFrac, SA, bergConc,\
        icebergs_depths2D, icebergs_widths2D, icebergs_length2D,\
        inversePowerLawDistNumbers_depth, inversePowerLawDistNumbers_width, inversePowerLawDistNumbers_length = \
        populate_masks_with_icebergs(iceberg_mask, dXC, dYC, hFacC, rA, delR)



    # write_bin('brg_tracerMask.bin', brgTracerMask)

    bergMaskNums_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'bergMaskNums.bin')
    bergMaskNums.ravel('C').astype(precision).tofile(bergMaskNums_file)

    numBergsPerCell_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'numBergsPerCell.bin')
    numBergsPerCell.ravel('C').astype(precision).tofile(numBergsPerCell_file)

    openFrac_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'openFrac.bin')
    openFrac.ravel('C').astype(precision).tofile(openFrac_file)

    totalBergArea_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'totalBergArea.bin')
    SA.ravel('C').astype(precision).tofile(totalBergArea_file)

    bergConc_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'bergConc.bin')
    bergConc.ravel('C').astype(precision).tofile(bergConc_file)

    icebergs_depths2D_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                     'iceberg', 'icebergs_depths.bin')
    icebergs_depths2D.ravel('C').astype(precision).tofile(icebergs_depths2D_file)

    icebergs_widths2D_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                          'iceberg', 'icebergs_widths.bin')
    icebergs_widths2D.ravel('C').astype(precision).tofile(icebergs_widths2D_file)

    icebergs_lengths2D_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                          'iceberg', 'icebergs_length.bin')
    icebergs_length2D.ravel('C').astype(precision).tofile(icebergs_lengths2D_file)

    inversePowerLawDistNumbers_depth_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                            'iceberg', 'inversePowerLawDistNumbers_depth.bin')
    inversePowerLawDistNumbers_depth.astype(precision).tofile(inversePowerLawDistNumbers_depth_file)

    inversePowerLawDistNumbers_width_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                            'iceberg', 'inversePowerLawDistNumbers_width.bin')
    inversePowerLawDistNumbers_width.astype(precision).tofile(inversePowerLawDistNumbers_width_file)

    inversePowerLawDistNumbers_length_file = os.path.join(config_dir, 'L2', L2_model_name, 'input',
                                            'iceberg', 'inversePowerLawDistNumbers_length.bin')
    inversePowerLawDistNumbers_length.astype(precision).tofile(inversePowerLawDistNumbers_length_file)



config_dir = '/Volumes/upernavik/Research/Ocean_Modeling/Projects/Downscale_Darwin/' \
             'darwin3/configurations/downscale_darwin/'

project_dir = '/Users/mhwood/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

create_L2_iceberg_files(config_dir, L2_model_name='L2_Upernavik', print_level=10)


