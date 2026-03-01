
import os
# import requests
import numpy as np
import matplotlib.pyplot as plt
# import xarray as xr
import gsw
import json
import argparse
import netCDF4 as nc4
from datetime import datetime

def read_alamo_data_from_json(file_path):
    f = open(file_path, 'r')
    xx = json.load(f)[0]
    f.close()

    lat = []
    lon = []
    depth = []
    date = []
    potential_temperature = []
    practical_salinity = []
    error_found = False

    dives = xx['dives'][0]

    if 'ascending' in dives['science'].keys():
        if len(dives['trajectory']['gps'])>0:
            lon = dives['trajectory']['gps'][0]['lon']
            lat = dives['trajectory']['gps'][0]['lat']
            date = np.datetime64(dives['trajectory']['gps'][0]['datetime'][:-1])
            pressure = dives['science']['ascending']['binned']['pressure']
            temperature = dives['science']['ascending']['binned']['temperature']
            practical_salinity = dives['science']['ascending']['binned']['salinity']
        else:
            error_found = True
    else:
        error_found = True

    if lon==0 or lat==0 or lon==None or lat==None:
        error_found=True

    if not error_found:
        depth = -1*gsw.conversions.z_from_p(pressure, lat)
        absolute_salinity = gsw.conversions.SA_from_SP(practical_salinity, pressure, lon, lat)
        potential_temperature = gsw.conversions.pt_from_t(absolute_salinity, temperature,pressure, 0)

        # plt.subplot(1,2,1)
        # plt.plot(potential_temperature,depth)
        # plt.title(str(lon)+', '+str(lat))
        # plt.subplot(1, 2, 2)
        # plt.plot(practical_salinity, depth)
        # plt.title(str(date))
        # plt.show()
    else:

        print('Error found in file '+file_path)

    return(lon, lat, date, depth, potential_temperature, practical_salinity, error_found)

def output_file_as_nc(output_file, lon, lat,
                      year,month,day,hour,minute,
                      binned_depth,binned_potential_temperature,binned_practical_salinity,source,project_name):

    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('records', len(binned_depth))
    dep = ds.createVariable('depth', 'f4', 'records')
    dep[:] = binned_depth
    salt = ds.createVariable('practical_salinity', 'f4', 'records')
    salt[:] = binned_practical_salinity
    ptemp = ds.createVariable('potential_temperature', 'f4', 'records')
    ptemp[:] = binned_potential_temperature

    ds.longitude = lon
    ds.latitude = lat
    ds.year = year
    ds.month = month
    ds.day = day
    ds.hour = hour
    ds.minute = minute
    ds.source = source
    ds.project_name = project_name

    ds.close()


def clean_and_bin_alamo_data(data_dir, project_dir, float_ID, printing):

    json_files = []
    float_data_path = os.path.join(data_dir,  float_ID)
    for file_name in os.listdir(float_data_path):
        if file_name[-5:] == '.json' and file_name not in ['Dive_038.json','Dive_043.json']:
            json_files.append(file_name)

    json_files = sorted(json_files)

    metadata_output ='File_ID,Float_ID,Year,Month,Day,Hour,Minute,Second,Longitude,Latitude'

    for file_name in json_files:
        if printing:
            print('    Binning file ' + file_name +' ('+str(json_files.index(file_name)+1)+' of '+str(len(json_files))+')')
        file_path = os.path.join(float_data_path, file_name)
        lon, lat, date, depth, potential_temperature, practical_salinity, error_found = read_alamo_data_from_json(file_path)
        if not error_found:
            print('   Lon: '+str(lon)+', Lat: '+str(lat))
            year = date.astype(object).year
            month = date.astype(object).month
            day = date.astype(object).day
            hour = date.astype(object).hour
            minute = date.astype(object).minute
            second = date.astype(object).second
            file_ID = str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)
            file_ID +='_'+ '{:02d}'.format(hour)+'{:02d}'.format(minute)+'{:02d}'.format(second)
            source = file_name
            project_name = 'Oceans Melting Greenland'
            output_file = os.path.join(os.path.join(project_dir,str(float_ID),
                                                    file_ID+'_F'+str(float_ID)+'.nc'))
            output_file_as_nc(output_file, lon, lat,
                              year, month, day, hour, minute,
                              depth, potential_temperature, practical_salinity, source,
                              project_name)
            metadata_output += '\n'+file_ID+',' + str(float_ID)+','+str(year)+','+str(month)+','+str(day)+','
            metadata_output += str(hour)+','+str(minute)+','+str(second)+','+str(lon)+','+str(lat)

    metadata_output_file = os.path.join(os.path.join(project_dir,str(float_ID)+'_Locations.csv'))
    f = open(metadata_output_file,'w')
    f.write(metadata_output)
    f.close()


def read_apex_data_from_csv(file_path):
    f = open(file_path)
    lines = f.read()
    f.close()
    lines=lines.split('\n')

    start_line_found = False
    for ll in range(len(lines)):
        line = lines[ll]
        if not start_line_found:
            if 'CP stopped' in line or 'CP already finished' in line:
                start_line_found=True
                start_line = ll+1

    lat = 0
    lon = 0
    depth = []
    date = []
    practical_salinity = []
    pressure = []
    potential_temperature = []
    error_found = False
    date_found = False

    if start_line_found:
        print(file_path)
        temperature = []
        for ll in range(start_line,len(lines)):
            line = lines[ll].split(',')
            if line[0]=='LGR_CP_PTSCI':
                pressure.append(float(line[2]))
                temperature.append(float(line[3]))
                practical_salinity.append(float(line[4]))
            if line[0]=='GPS':
                date = line[1]
                lat = float(line[2])
                lon = float(line[3])
                date_found = True

        # if np.abs(lon)>0:

        pressure = np.array(pressure)
        temperature = np.array(temperature)
        practical_salinity = np.array(practical_salinity)

        depth = -1*gsw.conversions.z_from_p(pressure, lat)
        absolute_salinity = gsw.conversions.SA_from_SP(practical_salinity, pressure, lon, lat)
        potential_temperature = gsw.conversions.pt_from_t(absolute_salinity, temperature,pressure, 0)

            # plt.subplot(1,2,1)
            # plt.plot(potential_temperature,depth)
            # plt.title(str(lon)+', '+str(lat))
            # plt.subplot(1, 2, 2)
            # plt.plot(practical_salinity, depth)
            # plt.title(str(date))
            # plt.show()
        # else:
        #     error_found = True
        #     print('Longitude was not found after CP stopped in ' + file_path)
    else:
        error_found = True
        print('No CP stopped message was found in '+file_path)

    if not date_found:
        date = lines[0].split(',')[1]

    return(lon, lat, date, depth, potential_temperature, practical_salinity, error_found)


def clean_and_bin_apex_data(data_dir, project_dir, float_ID, printing):

    # if 'Floats' not in os.listdir(os.path.join(project_dir,'Data')):
    #     os.mkdir(os.path.join(project_dir,'Data','Floats'))
    # if 'F' + str(float_ID) not in os.listdir(os.path.join(project_dir,'Data','Floats')):
    #     os.mkdir(os.path.join(project_dir,'Data','Floats','F' + str(float_ID)))
    # if 'Floats' not in os.listdir(os.path.join(project_dir,'Metadata')):
    #     os.mkdir(os.path.join(project_dir,'Metadata','Floats'))

    science_files = []
    float_data_path = os.path.join(data_dir, 'APEX', float_ID, 'Processed')
    for file_name in os.listdir(float_data_path):
        if file_name[-len('.science_log.csv'):] == '.science_log.csv':
            science_files.append(file_name)

    science_files = sorted(science_files)

    print(science_files)

    metadata_output ='File_ID,Float_ID,Year,Month,Day,Hour,Minute,Second,Longitude,Latitude'

    for file_name in science_files:
        # if printing:
        #     print('    Binning file ' + file_name +' ('+str(science_files.index(file_name)+1)+' of '+str(len(science_files))+')')
        file_path = os.path.join(float_data_path, file_name)
        lon, lat, date, depth, potential_temperature, practical_salinity, error_found = read_apex_data_from_csv(file_path)
        if not error_found:
            # print('   Lon: '+str(lon)+', Lat: '+str(lat))
            year = int(date.split('T')[0][:4])
            month = int(date.split('T')[0][4:6])
            day = int(date.split('T')[0][6:8])
            hour = int(date.split('T')[1][:2])
            minute = int(date.split('T')[1][2:4])
            second = int(date.split('T')[1][4:6])

            file_ID = str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)
            file_ID +='_'+ '{:02d}'.format(hour)+'{:02d}'.format(minute)+'{:02d}'.format(second)
            source = file_name
            project_name = 'Oceans Melting Greenland'
            output_file = os.path.join(os.path.join(project_dir,str(float_ID),
                                                    file_ID+'_F'+str(float_ID)+'.nc'))
            output_file_as_nc(output_file, lon, lat,
                              year, month, day, hour, minute,
                              depth, potential_temperature, practical_salinity, source,
                              project_name)
            metadata_output += '\n'+file_ID+',F' + str(float_ID)+','+str(year)+','+str(month)+','+str(day)+','
            metadata_output += str(hour)+','+str(minute)+','+str(second)+','+str(lon)+','+str(lat)

    metadata_output_file = os.path.join(os.path.join(project_dir,str(float_ID)+'_Locations.csv'))
    f = open(metadata_output_file,'w')
    f.write(metadata_output)
    f.close()


def bin_ctd_data(depth,potential_temperature,practical_salinity):

    max_depth = np.max(depth)
    binned_depth = np.arange(1,int(max_depth) + 1) #ignore the 1m surface layer - data is not trustworthy here
    binned_potential_temperature = np.zeros_like(binned_depth).astype(float)
    binned_practical_salinity = np.zeros_like(binned_depth).astype(float)
    binning_count = np.zeros_like(binned_depth).astype(int)
    for di in range(len(binned_depth) - 1):
        depth_indices = np.logical_and(depth >= binned_depth[di], depth < binned_depth[di + 1])
        if np.sum(depth_indices)>0:
            binned_potential_temperature[di] = np.mean(potential_temperature[depth_indices])
            binned_practical_salinity[di] = np.mean(practical_salinity[depth_indices])
            binning_count[di] = np.sum(depth_indices)
        # else:
        #     print('    No data for depth '+str(binned_depth[di]))

    # remove bins that didnt have any data
    has_data = binning_count > 0
    binned_depth = binned_depth[has_data]
    binned_potential_temperature = binned_potential_temperature[has_data]
    binned_practical_salinity = binned_practical_salinity[has_data]

    #filter out some ridiculous values
    salt_indices = np.logical_and(binned_practical_salinity>10,binned_practical_salinity<45)
    binned_depth = binned_depth[salt_indices]
    binned_potential_temperature = binned_potential_temperature[salt_indices]
    binned_practical_salinity = binned_practical_salinity[salt_indices]

    return(binned_depth,binned_potential_temperature,binned_practical_salinity)



def clean_and_bin_float_data(data_dir, project_dir, float_ID, float_type, printing):

    ########################################################################################################################
    # Use the Alamo read functions, if its an Alamo float

    if float_type == 'Alamo':
        clean_and_bin_alamo_data(data_dir, project_dir, float_ID, printing)

    ########################################################################################################################
    # Use the APEX read functions, if its an APEX float

    elif float_type == 'APEX':
        clean_and_bin_apex_data(data_dir, project_dir, float_ID, printing)

    else:
        raise ValueError('Float type '+float_type+' is not recognized - must be Alamo or APEX')



# data_dir = '/Users/mike/Documents/Research/Projects/Disko Bay/Data/Floats/Raw'
# project_dir = '/Users/mike/Documents/Research/Projects/Disko Bay/Data/Floats/Processed'

data_dir = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_Floats'
project_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik/Data/Ocean/Float Profiles'


float_ID = 'F9186'
float_type='APEX'

# float_ID = 'F9313'
# float_type='Alamo'

printing = True

clean_and_bin_float_data(data_dir, project_dir, float_ID, float_type, printing)