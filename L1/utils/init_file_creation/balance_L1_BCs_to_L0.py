
import os
import shutil
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import datetime
import argparse
import ast

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = datetime.datetime(year,month,day,hour,minute,second)
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def read_flux_timeseries(config_dir, L1_model_name, boundaries, field_name, level='L1', balanced = True):

    dec_yr = []

    # make the dec yr array
    for year in range(1992,2030):
        for month in range(1,13):
            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            for d in range(1,nDays+1):
                dec_yr.append(YMD_to_DecYr(year,month,d))

    time = np.array(dec_yr)

    if level=='L1' and balanced:
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L1_timeseries','balanced',L1_model_name+'_'+field_name+'_boundary_flux.nc')
    if level=='L1' and not balanced:
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L1_timeseries','unbalanced',L1_model_name+'_'+field_name+'_boundary_flux.nc')
    if level=='L0':
        file_path = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                     'timeseries','L0_timeseries','L0_'+L1_model_name[3:]+'_'+field_name+'_boundary_flux.nc')

    all_timeseries = []
    for boundary in boundaries:
        ds = nc4.Dataset(file_path)
        grp = ds.groups[boundary]
        timeseries = grp.variables['integrated_timeseries'][:]
        ds.close()
        all_timeseries.append(timeseries)

    time = time[:len(timeseries)]

    return(time,all_timeseries)

def calculate_integrated_flux_anomaly(time,timeseries):

    anomaly_timeseries = timeseries - np.mean(timeseries)

    timestep = 24*60*60 # 1 day

    integrated_anomaly = np.zeros_like(timeseries)
    for i in range(len(timeseries)-1):
        integrated_anomaly[i+1] = integrated_anomaly[i] + anomaly_timeseries[i]*timestep

    return(integrated_anomaly)

def read_grid_information(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    rA = np.array(ds.variables['rA'][:, :])
    dxG = np.array(ds.variables['dxG'][:, :])
    dyG = np.array(ds.variables['dyG'][:, :])
    hFacS = np.array(ds.variables['HFacS'][:, :, :])
    hFacW = np.array(ds.variables['HFacW'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return (rA, dxG, dyG, hFacS, hFacW, delR)

def create_dest_file_list(mask_name, var_name, start_year, final_year):

    dest_files = []
    start_date = datetime.datetime(start_year, 1, 1)
    final_date = datetime.datetime(final_year, 12, 31)
    years = []
    for year in range(1992,2022):
        for month in range(1, 13):
            test_date = datetime.datetime(year, month, 15)
            if test_date >= start_date and test_date <= final_date:
                dest_files.append(str(year) + '{:02d}'.format(month))
                if year not in years:
                    years.append(year)

    if var_name in ['UVEL','VVEL','UICE','VICE']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    dest_files_per_year = {}
    for year in years:
        year_files = []
        for dest_file in dest_files:
            if str(year) in dest_file:
                year_files.append('L1_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix)
        dest_files_per_year[year] = year_files

    return(years, dest_files_per_year)

def read_monthly_tracer_files(tracer_dir, file_list,year,mask_name,n_rows,n_cols,Nr):

    if year%4==0:
        total_timesteps = 366
    else:
        total_timesteps = 365

    if mask_name in ['north','south']:
        total_grid = np.zeros((total_timesteps, Nr, n_cols))
    if mask_name in ['east','west']:
        total_grid = np.zeros((total_timesteps, Nr, n_rows))

    timesteps_counted = 0
    for file_name in sorted(file_list):
        grid = np.fromfile(os.path.join(tracer_dir,file_name), '>f4')
        if mask_name in ['north', 'south']:
            timesteps = int(np.size(grid) / (n_cols * Nr))
            grid = np.reshape(grid, (timesteps, Nr, n_cols))
        if mask_name in ['east', 'west']:
            timesteps = int(np.size(grid) / (n_rows * Nr))
            grid = np.reshape(grid, (timesteps, Nr, n_rows))
        total_grid[timesteps_counted:timesteps_counted+timesteps,:,:] = grid
        timesteps_counted += timesteps

    return(total_grid)


def balance_flux_on_boundary(boundary, L0_timeseries, L1_timeseries, L1_vel_timeseries):

    # consider the original concentrations on the boundary (e.g. mol/m3 of NO3, or C)
    # we spatially integrated this across the boundary to get concentration per unit time (e.g. C * m2 * m/s)
    # now, we find the difference between the anomalies (C * m2 * m/s) and divide by the total flux (m3/s)
    #       to get a correction to the concentrations so there is no anomalies resulting from interpolation
    L0_anomaly_timeseries = L0_timeseries# - np.mean(L1_timeseries)
    L1_anomaly_timeseries = L1_timeseries# - np.mean(L1_timeseries)

    correction = (L0_anomaly_timeseries-L1_anomaly_timeseries)/L1_vel_timeseries
    print(correction)

    # plt.title(boundary)
    # plt.subplot(4,1,1)
    # plt.plot(L0_anomaly_timeseries)
    # plt.plot(L1_anomaly_timeseries)
    # plt.subplot(4, 1, 2)
    # plt.plot(L0_anomaly_timeseries - L1_anomaly_timeseries)
    # plt.subplot(4, 1, 3)
    # plt.plot(L1_vel_timeseries)
    # plt.subplot(4,1,4)
    # plt.plot(correction)
    # plt.show()

    balanced_flux = 1
    return(balanced_flux)



def balance_bc_field(config_dir, L1_model_name, boundaries, var_name, start_year, final_year, print_level):

    if print_level >= 1:
        print('    - Balancing the tracer fluxes into the model domain to correct for interpolation biases')

    # this is the dir where the obcs output will be stored
    if 'obcs' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs'))
    if 'balanced' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'balanced'))
    output_dir = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'balanced')

    for boundary in boundaries:
        if boundary not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, boundary))
        if var_name!='VEL':
            if var_name not in os.listdir(os.path.join(output_dir, boundary)):
                os.mkdir(os.path.join(output_dir, boundary, var_name))

    # read in the grid information
    rA, dxG, dyG, hFacS, hFacW, delR = read_grid_information(config_dir, L1_model_name)
    Nr = len(delR)

    # step 1: read the flux timeseries for a given year
    L1_time, all_L1_balanced_velocity_timeseries = read_flux_timeseries(config_dir, L1_model_name, boundaries, 'VEL', balanced=True)
    L1_time, all_L1_unbalanced_timeseries = read_flux_timeseries(config_dir, L1_model_name, boundaries, var_name, balanced=False)
    L0_time, all_L0_timeseries = read_flux_timeseries(config_dir, L1_model_name, boundaries, var_name, level='L0')



    for b in range(len(boundaries)):
        boundary = boundaries[b]

        _, dest_files_per_year = create_dest_file_list(boundary, var_name, start_year, final_year)

        for year in range(start_year, final_year+1):

            start_index = 0
            end_index = 366
            # make the dec yr array
            for yr in range(1992, 2030):
                if yr < year:
                    if yr % 4 == 0:
                        start_index += 366
                    else:
                        start_index += 365
                    if (yr + 1) % 4 == 0:
                        end_index += 365
                    else:
                        end_index += 366

            L0_timeseries = all_L0_timeseries[b][start_index:end_index]
            L1_timeseries = all_L1_unbalanced_timeseries[b][start_index:end_index]
            L1_vel_timeseries = all_L1_balanced_velocity_timeseries[b][start_index:end_index]

            correction = (L0_timeseries - L1_timeseries) / L1_vel_timeseries

            obcs_folder = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs')
            # # read in the var file
            # var_file = os.path.join(obcs_folder, 'unbalanced', boundary, var_name,
            #                               'L1_BC_'+boundary+'_'+var_name+'_' + str(year))
            # var_grid = np.fromfile(var_file, '>f4')
            # n_timesteps = end_index-start_index
            # if boundary in ['north','south']:
            #     N = np.shape(rA)[0]
            # if boundary in ['east','west']:
            #     N = np.shape(rA)[1]
            # var_grid = np.reshape(var_grid, (n_timesteps, Nr, N)).astype(float)

            n_rows = np.shape(rA)[0]
            n_cols = np.shape(rA)[1]
            obcs_folder = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs')
            tracer_dir = os.path.join(obcs_folder,'unbalanced',boundary,var_name)
            var_grid = read_monthly_tracer_files(tracer_dir, dest_files_per_year[year], year, boundary, n_rows, n_cols, Nr)

            for timestep in range(np.shape(var_grid)[0]):
                time_slice = var_grid[timestep, :, :]
                time_slice[time_slice!=0] += correction[timestep]
                var_grid[timestep,:,:] = time_slice

            output_file = os.path.join(obcs_folder, 'balanced', boundary, var_name,
                                          'L1_BC_'+boundary+'_'+var_name+'_' + str(year))
            var_grid.ravel('C').astype('>f4').tofile(output_file)






def balance_bc_fields(config_dir, L1_model_name, boundaries, proc_id, start_year, final_year, print_level):

    if L1_model_name in ['L1_W_Greenland','L1_mac_delta']:
        var_names = ['THETA','SALT','VEL','UICE','VICE','HSNOW','HEFF','AREA']
    else:
        var_names = ['THETA', 'SALT', 'VEL']
    for i in range(1,32):
        var_names.append('PTRACE'+'{:02d}'.format(i))

    var_name_list = []
    for var_name in var_names:
        var_name_list.append(var_name)

    var_name = var_name_list[proc_id % len(var_name_list)]

    print('Creating the L1 timeseries for ' + var_name + ' to cover years ' +
          str(start_year)+' to ' + str(final_year))

    print('  Running the balancing routine:')
    balance_bc_field(config_dir, L1_model_name, boundaries, var_name, start_year, final_year, print_level)









if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--L1_model_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="L1_model_name",
                        type=str, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name
    start_year = args.start_year
    final_year = args.final_year

    balance_bc_fields(config_dir, L1_model_name, start_year, final_year, print_level=4)
