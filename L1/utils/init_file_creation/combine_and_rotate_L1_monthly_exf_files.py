
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def get_dest_file_list(var_name, n_rows_L1, n_cols_L1, year):

    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0

    print(' WARNING: Months set to 2 for testing in combine_and_rotate_L1_monthly_exf')

    for month in range(1, 2):
        if month in [1, 3, 5, 7, 8, 10, 12]:
            nDays = 31
        elif month in [4, 6, 9, 11]:
            nDays = 30
        else:
            if year % 4 == 0:
                nDays = 29
            else:
                nDays = 28
        if var_name in ['UWIND','VWIND']:
            dest_file = 'L1_exf_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '_rotated.bin'
        else:
            dest_file = 'L1_exf_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '.bin'
        dest_files.append(dest_file)
        nTimesteps = nDays*4
        total_timesteps += nTimesteps
        dest_file_shapes[dest_file] = (nTimesteps, n_rows_L1, n_cols_L1)

    return(dest_files,dest_file_shapes,total_timesteps)

def stack_monthly_exf_files_to_one(config_dir, config_name, var_name, AngleCS, AngleSN, dest_files, dest_file_shapes,total_timesteps,print_level):

    rows = dest_file_shapes[dest_files[0]][1]
    cols = dest_file_shapes[dest_files[0]][2]

    output_grid = np.zeros((total_timesteps,rows,cols))
    timesteps_added = 0
    for dest_file in dest_files:

        # ymd = dest_file.split('_')[3]#[:8]
        # if int(ymd[4:6])<3:
        #     var_grid = np.zeros(dest_file_shapes[dest_file])
        #     print('WARNING - FILLING IN VALUES WITH ZEROS FOR TESTING')
        #     if print_level >= 2:
        #         print('        - Adding timesteps from file ' + dest_file + ' to levels ' + str(
        #             timesteps_added) + ' to ' + str(timesteps_added + np.shape(var_grid)[0]))
        # elif int(ymd[4:6])==3 and int(ymd[6:8])<30:
        #     var_grid = np.zeros(dest_file_shapes[dest_file])
        #     print('WARNING - FILLING IN VALUES WITH ZEROS FOR TESTING')
        #     if print_level >= 2:
        #         print('        - Adding timesteps from file ' + dest_file + ' to levels ' + str(
        #             timesteps_added) + ' to ' + str(timesteps_added + np.shape(var_grid)[0]))
        # else:

        if var_name in ['UWIND','VWIND']:
            if var_name=='UWIND':
                u_dest_file = dest_file
            else:
                u_dest_file = dest_file.replace('VWIND','UWIND')
            if var_name=='VWIND':
                v_dest_file = dest_file
            else:
                v_dest_file = dest_file.replace('UWIND','VWIND')
            u_var_grid = np.fromfile(os.path.join(config_dir, 'L1', config_name, 'input', 'exf', 'UWIND', u_dest_file), '>f4')
            u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
            v_var_grid = np.fromfile(os.path.join(config_dir, 'L1', config_name, 'input', 'exf', 'VWIND', v_dest_file),'>f4')
            v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])

            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(np.shape(var_grid)[0]):
                if var_name=='UWIND':
                    var_grid[timestep, :, :] = AngleCS*u_var_grid[timestep, :, :] + AngleSN*v_var_grid[timestep, :, :]
                if var_name=='VWIND':
                    var_grid[timestep, :, :] = -1*AngleSN*u_var_grid[timestep, :, :] + AngleCS*v_var_grid[timestep, :, :]
        else:
            var_grid = np.fromfile(os.path.join(config_dir,'L1',config_name, 'input','exf',var_name,dest_file),'>f4')
            var_grid = np.reshape(var_grid,dest_file_shapes[dest_file])

        if print_level >= 2:
            print('        - Adding timesteps from file ' + dest_file + ' to levels ' + str(
                timesteps_added) + ' to ' + str(timesteps_added + np.shape(var_grid)[0]))

        if print_level >= 4:
            print('                - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                - There are '+str(np.sum(np.isnan(var_grid[0,:,:])))+' nan values')


        output_grid[timesteps_added:timesteps_added+np.shape(var_grid)[0],:,:] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    return(output_grid)

def combine_and_rotate_L1_monthly_exf_files(config_dir, config_name, proc_id,start_year,final_year,print_level):

    var_name_list = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP','RUNOFF','IRONDUST','ATMOSCO2']
    var_name = var_name_list[proc_id % len(var_name_list)]

    grid_file = os.path.join(config_dir,'nc_grids',config_name+'_grid.nc')
    ds = nc4.Dataset(grid_file)
    AngleCS = ds.variables['AngleCS'][:,:]
    AngleSN = ds.variables['AngleSN'][:, :]
    n_rows_L1 = np.shape(AngleCS)[0]
    n_cols_L1 = np.shape(AngleCS)[1]
    ds.close()

    years = np.arange(start_year,final_year+1).tolist()

    if print_level >=1:
        print('    - Combining monthly exf files for ' + var_name)

    for year in years:

        if print_level >= 2:
            print('        - Combining files in year ' + str(year))

        if print_level >=3:
            print('            - Getting a list of monthly files')
        dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(var_name, n_rows_L1, n_cols_L1, year)

        if print_level >= 1:
            print('    - Stacking all of the monthly files into a big global file')
        output_grid = stack_monthly_exf_files_to_one(config_dir, config_name, var_name, AngleCS, AngleSN, dest_files,
                                                   dest_file_shapes, total_timesteps, print_level)

        output_file = os.path.join(config_dir, 'L1',config_name, 'input', 'exf', 'L1_exf_' + var_name + '_'+str(year))
        if print_level >= 1:
            print('    - Outputting to ' + output_file)
        output_grid.ravel('C').astype('>f4').tofile(output_file)
