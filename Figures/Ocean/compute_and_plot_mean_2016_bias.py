

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d
import datetime

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

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
    return(XC, YC, Z, Z_top, Z_bottom, Depth, hFaC)

def read_L2_bias_grids(project_folder):
    output_file = os.path.join(project_folder, 'Data','Ocean', 'Modeling', 'L1_state_bias_vs_OMG.nc')
    ds = nc4.Dataset(output_file)
    depth = ds.variables['depths'][:]
    Theta = ds.variables['Theta'][:, :]
    Salt = ds.variables['Salt'][:, :]
    ds.close()
    return(depth, Theta, Salt)

def read_L2_CTD_dv_output(config_dir):

    months = np.arange(2,12).tolist()

    for month in months:
        theta_file = os.path.join(config_dir,'L2','L2_Upernavik','results_baseline', 'dv', 'CTD',
                                   'THETA', 'THETA_2016'+'{:02d}'.format(month)+'.nc')
        ds = nc4.Dataset(theta_file)
        depth = ds.variables['depths'][:]
        Theta = ds.variables['THETA'][:, :, :]
        longitude = ds.variables['longitude'][:]
        latitude = ds.variables['latitude'][:]
        ds.close()

        salt_file = os.path.join(config_dir,'L2','L2_Upernavik','results_baseline', 'dv', 'CTD',
                                 'SALT', 'SALT_2016'+'{:02d}'.format(month)+'.nc')
        ds = nc4.Dataset(salt_file)
        Salt = ds.variables['SALT'][:, :, :]
        ds.close()

        dec_yrs = []
        for day in range(1,np.shape(Theta)[0]+1):
            dec_yr = YMD_to_DecYr(2016, month, day)
            dec_yrs.append(dec_yr)
        dec_yrs = np.array(dec_yrs)

        if month == months[0]:
            theta_timeseries = Theta
            salt_timeseries = Salt
            all_dec_yrs = dec_yrs
        else:
            theta_timeseries = np.concatenate([theta_timeseries, Theta], axis=0)
            salt_timeseries = np.concatenate([salt_timeseries, Salt], axis=0)
            all_dec_yrs = np.concatenate([all_dec_yrs, dec_yrs], axis=0)

    return(depth, all_dec_yrs, theta_timeseries, salt_timeseries)

def retrieve_model_obs_pairs(ctd_dir, project_folder, model_depth, model_dec_yrs, theta_timeseries, salt_timeseries):

    mapping_file = os.path.join(project_folder, 'Data', 'Ocean',  'AXCTDs',
                                'AXCTD List at Model Points.csv')
    f = open(mapping_file, 'r')
    lines = f.readlines()
    f.close()
    lines.pop(0)

    model_profiles = []
    obs_profiles = []

    for line in lines:
        line = line.split(',')
        file_name = line[0]
        model_idx = int(line[1])

        if file_name[:4]=='2016':
            # print(file_name)
            obs_file_path = os.path.join(ctd_dir, '2016', 'CTD_'+file_name+'.nc')
            ds = nc4.Dataset(obs_file_path)
            depth = ds.variables['depth'][:]
            theta_obs = ds.variables['potential_temperature'][:]
            salt_obs = ds.variables['practical_salinity'][:]
            ds.close()

            obs_profile = np.column_stack([depth, theta_obs, salt_obs])
            obs_dec_yr = YMD_to_DecYr(int(file_name[0:4]), int(file_name[4:6]), int(file_name[6:8]))
            obs_profiles.append(obs_profile)

            model_time_idx = np.argmin(np.abs(model_dec_yrs - obs_dec_yr))
            theta_model_profile = theta_timeseries[model_time_idx, :, model_idx]
            salt_model_profile = salt_timeseries[model_time_idx, :, model_idx]
            model_profile = np.column_stack([model_depth, theta_model_profile, salt_model_profile])
            model_profiles.append(model_profile)


    return(model_profiles, obs_profiles)

def compute_bias_profiles(model_profiles, obs_profiles):

    bias_profiles = []
    for i in range(len(model_profiles)):
        model_profile = np.array(model_profiles[i])
        obs_profile = np.array(obs_profiles[i])

        # remove 0's
        model_profile = model_profile[(model_profile[:,1]!=0) & (model_profile[:,2]!=0), :]
        obs_profile = obs_profile[(obs_profile[:,1]!=0) & (obs_profile[:,2]!=0), :]

        obs_theta_interp = interp1d(obs_profile[:,0], obs_profile[:,1], bounds_error=False, fill_value=np.nan)
        obs_salt_interp = interp1d(obs_profile[:,0], obs_profile[:,2], bounds_error=False, fill_value=np.nan)

        obs_theta_aligned = obs_theta_interp(model_profile[:,0])
        obs_salt_aligned = obs_salt_interp(model_profile[:,0])

        theta_bias = model_profile[:,1] - obs_theta_aligned
        salt_bias = model_profile[:,2] - obs_salt_aligned

        bias_profile = np.column_stack([model_profile[:,0], theta_bias, salt_bias])
        bias_profiles.append(bias_profile)

    return(bias_profiles)

def compute_mean_profiles(profile_set, Z, Z_top, Z_bottom):

    mean_theta_profile = np.zeros((np.shape(Z)[0],))
    count_theta_profile = np.zeros((np.shape(Z)[0],))
    mean_salt_profile = np.zeros((np.shape(Z)[0],))
    count_salt_profile = np.zeros((np.shape(Z)[0],))

    for z in range(len(Z)):

        for d in range(len(profile_set)):

            # get the profile
            profile_grid = profile_set[d]
            non_zero_indices = profile_grid[:,1]!=0
            profile_grid=profile_grid[non_zero_indices,:]

            # subset profile to depth range
            depth_indices = np.logical_and(profile_grid[:,0]>=Z_top[z], profile_grid[:,0]<=Z_bottom[z])

            # add the points to the profile
            if np.sum(depth_indices)>0:
                non_zero_indices = depth_indices
                mean_theta_profile[z] += np.sum(profile_grid[depth_indices,1])
                count_theta_profile[z] += np.sum(depth_indices)
                mean_salt_profile[z] += np.sum(profile_grid[depth_indices,2])
                count_salt_profile[z] += np.sum(depth_indices)

    for d in range(len(mean_theta_profile)):
        if count_theta_profile[d]>0:
            mean_theta_profile[d] = mean_theta_profile[d]/count_theta_profile[d]
        else:
            mean_theta_profile[d] = mean_theta_profile[d-1]

        if count_salt_profile[d]>0:
            mean_salt_profile[d] = mean_salt_profile[d]/count_salt_profile[d]
        else:
            mean_salt_profile[d] = mean_salt_profile[d-1]

    return(mean_theta_profile, mean_salt_profile)

def compute_mean_bias_profile(Z, mean_theta_profile_model, mean_salt_profile_model, mean_theta_profile_obs, mean_salt_profile_obs):
    theta_bias_profile = mean_theta_profile_model - mean_theta_profile_obs
    salt_bias_profile = mean_salt_profile_model - mean_salt_profile_obs

    # taper from 100% at 100m to 50% at surface
    for z in range(len(Z)):
        if Z[z]<100:
            taper_factor = 0.5 + 0.5*(Z[z]/100)
            theta_bias_profile[z] = theta_bias_profile[z]*taper_factor
            salt_bias_profile[z] = salt_bias_profile[z]*taper_factor

    return(theta_bias_profile, salt_bias_profile)


def plot_profiles_and_bias_grids(project_folder, model_profiles, obs_profiles,bias_profiles,
                                 mean_theta_profile_model, mean_salt_profile_model,
                                 mean_theta_profile_obs, mean_salt_profile_obs,
                                 theta_bias_profile, salt_bias_profile):

    fig = plt.figure(figsize=(10, 7))

    gs2 = GridSpec(2, 3, left=0.11, right=0.98, bottom=0.1, top=0.92, wspace=0.3, hspace=0.3)


    ##########################################################
    # model theta profiles

    ax = fig.add_subplot(gs2[0,0])
    for i in range(len(model_profiles)):
        model_profile = model_profiles[i]
        depth_profile = model_profile[:,0]
        theta_profile = model_profile[:,1]
        non_nan_indices = np.logical_and(~np.isnan(theta_profile), theta_profile!=0)
        C = ax.plot(theta_profile[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)
    ax.plot(mean_theta_profile_model, Z, 'r-', linewidth=2)
    ax.set_ylim([1000, 0])
    ax.set_xlim([-2,6])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_title('a) Baseline Model')
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Temperature ($^{\circ}$C)')

    ##########################################################
    # obs theta profiles
    ax = fig.add_subplot(gs2[0,1])
    for i in range(len(obs_profiles)):
        obs_profile = obs_profiles[i]
        depth_profile = obs_profile[:, 0]
        theta_profile = obs_profile[:, 1]
        non_nan_indices = np.logical_and(~np.isnan(theta_profile), theta_profile!=0)
        C = ax.plot(theta_profile[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)
    ax.plot(mean_theta_profile_obs, Z, 'r-', linewidth=2)
    ax.set_ylim([1000, 0])
    ax.set_xlim([-2, 6])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_title('b) OMG AXCTD Measurements')
    ax.set_xlabel('Temperature ($^{\circ}$C)')

    ##########################################################
    # theta bias

    ax = fig.add_subplot(gs2[0,2])
    for i in range(len(bias_profiles)):
        bias_profile = bias_profiles[i]
        depth_profile = bias_profile[:,0]
        theta_bias = bias_profile[:,1]
        non_nan_indices = ~np.isnan(theta_bias)
        C = ax.plot(theta_bias[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)
    ax.plot(theta_bias_profile,Z, 'r-')
    # ax.set_xlim([np.min(X), np.max(X)])
    ax.set_ylim([1000, 0])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_title('c) Model Bias \n (Model - Observations) ')
    ax.set_xlabel('Temperature ($^{\circ}$C)')

    ##########################################################
    # model salt profiles

    ax = fig.add_subplot(gs2[1,0])
    for i in range(len(model_profiles)):
        model_profile = model_profiles[i]
        depth_profile = model_profile[:,0]
        salt_profile = model_profile[:,2]
        non_nan_indices = np.logical_and(~np.isnan(salt_profile), salt_profile!=0)
        C = ax.plot(salt_profile[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)
    ax.plot(mean_salt_profile_model, Z, 'b-', linewidth=2)
    ax.set_ylim([1000, 0])
    ax.set_xlim([32,36])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Salinity (psu)')

    ##########################################################
    # obs salt profiles
    ax = fig.add_subplot(gs2[1,1])
    for i in range(len(obs_profiles)):
        obs_profile = obs_profiles[i]
        depth_profile = obs_profile[:, 0]
        salt_profile = obs_profile[:, 2]
        non_nan_indices = np.logical_and(~np.isnan(salt_profile), salt_profile!=0)
        C = ax.plot(salt_profile[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)
    ax.plot(mean_salt_profile_obs, Z, 'b-', linewidth=2)
    ax.set_ylim([1000, 0])
    ax.set_xlim([32, 36])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_xlabel('Salinity (psu)')



    ##########################################################
    # salt bias

    ax = fig.add_subplot(gs2[1, 2])
    for i in range(len(bias_profiles)):
        bias_profile = bias_profiles[i]
        depth_profile = bias_profile[:,0]
        salt_bias = bias_profile[:,2]
        non_nan_indices = ~np.isnan(salt_bias)
        C = ax.plot(salt_bias[non_nan_indices], depth_profile[non_nan_indices], '-', color='silver', linewidth=1)

    ax.plot(salt_bias_profile, Z, 'b-')
    # ax.set_xlim([np.min(X), np.max(X)])
    ax.set_ylim([1000, 0])
    ax.grid(linestyle='--', linewidth=0.5, alpha=0.5)
    # ax.set_title('b) Salinity Bias Correction')
    ax.set_xlabel('Salinity (psu)')

    output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'L2_Bias_Correction_Profile.png')
    # output_file = os.path.join(project_folder, 'Figures', 'Ocean', 'L2_Bias_Correction_Profile.pdf')
    plt.savefig(output_file)
    plt.close()
    a=1

def save_bias_profile_to_nc(project_folder, Z, theta_bias_profile, salt_bias_profile):

    output_file = os.path.join(project_folder, 'Data', 'Ocean','AXCTDs', 'L2_baseline_bias_2016_mean.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('depth', len(Z))

    t = ds.createVariable('depths', 'f4', ('depth',))
    t[:] = Z

    t = ds.createVariable('Theta', 'f4', ('depth', ))
    t[:] = theta_bias_profile

    t = ds.createVariable('Salt', 'f4', ('depth', ))
    t[:] = salt_bias_profile

    ds.close()
    a=1


project_folder = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

config_dir = '/Users/mike/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_Darwin/MITgcm/' \
             'configurations/downscale_darwin'

ctd_dir = '/Users/mike/Documents/Research/Data Repository/Greenland/Ocean Properties/OMG_CTDs/Processed'

# step 0: read the model grid
XC, YC, Z, Z_top, Z_bottom, Depth, hFaC = read_grid_geometry_from_nc(config_dir, model_name='L2_Upernavik')
# Z = Z[:51]

model_depth, dec_yrs, theta_timeseries, salt_timeseries = read_L2_CTD_dv_output(config_dir)

model_profiles, obs_profiles = retrieve_model_obs_pairs(ctd_dir, project_folder,
                                                        model_depth, dec_yrs, theta_timeseries, salt_timeseries)

bias_profiles = compute_bias_profiles(model_profiles, obs_profiles)

mean_theta_profile_model, mean_salt_profile_model = compute_mean_profiles(model_profiles, Z, Z_top, Z_bottom)

mean_theta_profile_obs, mean_salt_profile_obs = compute_mean_profiles(obs_profiles, Z, Z_top, Z_bottom)

theta_bias_profile, salt_bias_profile = compute_mean_bias_profile(Z, mean_theta_profile_model, mean_salt_profile_model,
                                                                    mean_theta_profile_obs, mean_salt_profile_obs)

plot_profiles_and_bias_grids(project_folder, model_profiles, obs_profiles, bias_profiles,
                             mean_theta_profile_model, mean_salt_profile_model,
                             mean_theta_profile_obs, mean_salt_profile_obs,
                             theta_bias_profile, salt_bias_profile)


save_bias_profile_to_nc(project_folder, Z, theta_bias_profile, salt_bias_profile)


