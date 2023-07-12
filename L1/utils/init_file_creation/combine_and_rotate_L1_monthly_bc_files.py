
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def get_dest_file_list(mask_name, var_name, Nr, n_rows_L1,n_cols_L1, year):

    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0

    for month in range(1, 13):

        if month in [1, 3, 5, 7, 8, 10, 12]:
            nTimesteps = 31
        elif month in [4, 6, 9, 11]:
            nTimesteps = 30
        else:
            if year % 4 == 0:
                nTimesteps = 29
            else:
                nTimesteps = 28

        if var_name in ['UVEL','VVEL','UICE','VICE']:
            dest_file = 'L1_BC_'+mask_name+'_'+ var_name + '.' + str(year) + '{:02d}'.format(month) + '_rotated.bin'
        else:
            dest_file = 'L1_BC_'+mask_name+'_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '.bin'
        dest_files.append(dest_file)

        total_timesteps += nTimesteps
        if mask_name in ['west','east']:
            dest_file_shapes[dest_file] = (nTimesteps, Nr, n_rows_L1, 1)
        if mask_name in ['north','south']:
            dest_file_shapes[dest_file] = (nTimesteps, Nr, 1, n_cols_L1)

    return(dest_files,dest_file_shapes,total_timesteps)

def stack_monthly_bc_files_to_one(config_dir, L1_model_name, mask_name, var_name, AngleCS, AngleSN, dest_files, dest_file_shapes,total_timesteps,print_level,read_from_unbalanced):

    depth_levels = dest_file_shapes[dest_files[0]][1]
    rows = dest_file_shapes[dest_files[0]][2]
    cols = dest_file_shapes[dest_files[0]][3]

    # the 2 is added because we will duplicate the first and last field
    output_grid = np.zeros((total_timesteps,depth_levels,rows,cols))
    timesteps_added = 0
    for dest_file in dest_files:

        # ymd = dest_file.split('.')[1]  # [:8]
        # print(dest_file,ymd)
        # if int(ymd[4:6]) < 3:
        #     var_grid = np.zeros(dest_file_shapes[dest_file])
        #     print('WARNING - FILLING IN VALUES WITH ZEROS FOR TESTING')
        #     if print_level >= 2:
        #         print('        - Adding timesteps from file ' + dest_file + ' to levels ' + str(
        #             timesteps_added) + ' to ' + str(timesteps_added + np.shape(var_grid)[0]))
        # elif int(ymd[4:6]) == 3 and int(ymd[6:8]) < 30:
        #     var_grid = np.zeros(dest_file_shapes[dest_file])
        #     print('WARNING - FILLING IN VALUES WITH ZEROS FOR TESTING')
        #     if print_level >= 2:
        #         print('        - Adding timesteps from file ' + dest_file + ' to levels ' + str(
        #             timesteps_added) + ' to ' + str(timesteps_added + np.shape(var_grid)[0]))
        # else:

        if var_name in ['UVEL','VVEL','UICE','VICE']:
            if 'VEL' in var_name:
                if var_name=='UVEL':
                    u_dest_file = dest_file
                else:
                    u_dest_file = dest_file.replace('VVEL','UVEL')
                if var_name=='VVEL':
                    v_dest_file = dest_file
                else:
                    v_dest_file = dest_file.replace('UVEL','VVEL')

                if read_from_unbalanced:
                    u_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','unbalanced',mask_name, 'UVEL', u_dest_file), '>f4')
                    v_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','unbalanced',mask_name, 'VVEL', v_dest_file),'>f4')
                else:
                    u_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', mask_name, 'UVEL',u_dest_file), '>f4')
                    v_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', mask_name, 'VVEL',v_dest_file), '>f4')
                u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
                v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])

            if 'ICE' in var_name:
                if var_name=='UICE':
                    u_dest_file = dest_file
                else:
                    u_dest_file = dest_file.replace('VICE','UICE')
                if var_name=='VICE':
                    v_dest_file = dest_file
                else:
                    v_dest_file = dest_file.replace('UICE','VICE')
                if read_from_unbalanced:
                    u_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','unbalanced',mask_name, 'UICE', u_dest_file), '>f4')
                    v_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','unbalanced',mask_name, 'VICE', v_dest_file),'>f4')
                else:
                    u_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', mask_name, 'UICE',u_dest_file), '>f4')
                    v_var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', mask_name, 'VICE',v_dest_file), '>f4')
                u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
                v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])

            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(np.shape(var_grid)[0]):
                for k in range(np.shape(var_grid)[1]):
                    if var_name=='UICE' or var_name=='UVEL':
                        var_grid[timestep, k, :, :] = AngleCS*u_var_grid[timestep, k, :, :] + AngleSN*v_var_grid[timestep, k, :, :]
                    if var_name=='VICE' or var_name=='VVEL':
                        var_grid[timestep, k, :, :] = -1*AngleSN*u_var_grid[timestep, k, :, :] + AngleCS*v_var_grid[timestep, k, :, :]

        else:
            if read_from_unbalanced:
                var_grid = np.fromfile(os.path.join(config_dir,'L1',L1_model_name, 'input','obcs','unbalanced',mask_name,var_name,dest_file),'>f4')
            else:
                var_grid = np.fromfile(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', mask_name, var_name,dest_file), '>f4')
            var_grid = np.reshape(var_grid,dest_file_shapes[dest_file])

        if print_level>=2:
            print('        - Adding timesteps from file '+dest_file+' to levels '+str(timesteps_added)+' to '+str(timesteps_added+np.shape(var_grid)[0]))
        if print_level >= 4:
            print('                - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                - There are '+str(np.sum(np.isnan(var_grid[0,:,:])))+' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added+np.shape(var_grid)[0],:,:,:] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    return(output_grid)

def combine_and_rotate_L1_monthly_bcs(config_dir, L1_model_name, boundaries, var_names,
                                      start_year, final_year, print_level, velocity_only = False,
                                      read_from_unbalanced = True, write_to_balanced = False):
    
    # if L1_model_name in ['L1_W_Greenland', 'L1_mac_delta']:
    #     var_names = ['THETA', 'SALT', 'UVEL', 'VVEL', 'UICE', 'VICE', 'HSNOW', 'HEFF', 'AREA']
    # else:
    #     var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    # for i in range(1, 32):
    #     var_names.append('PTRACE' + '{:02d}'.format(i))

    # var_name_list = []
    # mask_name_list = []
    # for var_name in var_names:
    #     for boundary in boundaries:
    #         if velocity_only:
    #             if var_name in ['UVEL','VVEL']:
    #                 var_name_list.append(var_name)
    #                 mask_name_list.append(boundary)
    #         else:
    #             var_name_list.append(var_name)
    #             mask_name_list.append(boundary)

    # var_name = var_name_list[proc_id % len(var_name_list)]
    # mask_name = mask_name_list[proc_id % len(var_name_list)]

    # grid_file = os.path.join(config_dir,'nc_grids',L1_model_name+'_grid.nc')
    # ds = nc4.Dataset(grid_file)
    # AngleCS = ds.variables['AngleCS'][:,:]
    # AngleSN = ds.variables['AngleSN'][:, :]
    # n_rows_L1 = np.shape(AngleCS)[0]
    # n_cols_L1 = np.shape(AngleCS)[1]
    # drF = ds.variables['drF'][:]
    # Nr = len(drF)
    # ds.close()

    for mask_name in boundaries:
        for var_name in var_names:

            grid_file = os.path.join(config_dir, 'nc_grids', L1_model_name + '_grid.nc')
            ds = nc4.Dataset(grid_file)
            AngleCS = ds.variables['AngleCS'][:, :]
            AngleSN = ds.variables['AngleSN'][:, :]
            n_rows_L1 = np.shape(AngleCS)[0]
            n_cols_L1 = np.shape(AngleCS)[1]
            drF = ds.variables['drF'][:]
            Nr = len(drF)
            ds.close()

    # print('There is an error at the boundary in the AngleCS/SN fields so taking')
    # print('these fields from one row/col on the interior')
    # if mask_name=='north':
    #     AngleCS = AngleCS[-1:,:]
    #     AngleSN = AngleSN[-1:,:]
    # if mask_name=='south':
    #     AngleCS = AngleCS[:1,:]
    #     AngleSN = AngleSN[:1,:]
    # if mask_name=='west':
    #     AngleCS = AngleCS[:,:1]
    #     AngleSN = AngleSN[:,:1]
    # if mask_name=='east':
    #     AngleCS = AngleCS[:,-1:]
    #     AngleSN = AngleSN[:,-1:]
            if mask_name=='north':
                AngleCS = AngleCS[-2:-1,:]
                AngleSN = AngleSN[-2:-1,:]
                AngleCS[0, 0] = AngleCS[0, 1]
                AngleCS[0, -1] = AngleCS[0, -2]
                AngleSN[0, 0] = AngleSN[0, 1]
                AngleSN[0, -1] = AngleSN[0, -2]
            if mask_name=='south':
                AngleCS = AngleCS[1:2,:]
                AngleSN = AngleSN[1:2,:]
                AngleCS[0, 0] = AngleCS[0, 1]
                AngleCS[0, -1] = AngleCS[0, -2]
                AngleSN[0, 0] = AngleSN[0, 1]
                AngleSN[0, -1] = AngleSN[0, -2]
            if mask_name=='west':
                AngleCS = AngleCS[:,1:2]
                AngleSN = AngleSN[:,1:2]
                AngleCS[0,0] = AngleCS[1,0]
                AngleCS[-1,0] = AngleCS[-2,0]
                AngleSN[0, 0] = AngleSN[1, 0]
                AngleSN[-1, 0] = AngleSN[-2, 0]
            if mask_name=='east':
                AngleCS = AngleCS[:,-2:-1]
                AngleSN = AngleSN[:,-2:-1]
                AngleCS[0, 0] = AngleCS[1, 0]
                AngleCS[-1, 0] = AngleCS[-2, 0]
                AngleSN[0, 0] = AngleSN[1, 0]
                AngleSN[-1, 0] = AngleSN[-2, 0]

            if var_name in ['ETAN','UICE','VICE','AREA','HSNOW','HEFF']:
                Nr = 1

            years = np.arange(start_year,final_year+1).tolist()

            if print_level >=1:
                print('    - Combining monthly BC files for ' + var_name)

            for year in years:

                if write_to_balanced:
                    if var_name=='UVEL' and mask_name=='east':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'unbalanced',mask_name,var_name,
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='UVEL' and mask_name=='west':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'unbalanced',mask_name,var_name,
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='VVEL' and mask_name=='north':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'unbalanced',mask_name,var_name,
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='VVEL' and mask_name=='south':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'unbalanced',mask_name,var_name,
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    else:
                        if var_name not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'balanced',mask_name)):
                            os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'balanced',mask_name,var_name))
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'balanced',mask_name,var_name,
                                                   'L1_BC_' + mask_name + '_' + var_name + '_' + str(year))
                else:
                    if var_name=='UVEL' and mask_name=='east':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'annual',
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='UVEL' and mask_name=='west':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'annual',
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='VVEL' and mask_name=='north':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'annual',
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    elif var_name=='VVEL' and mask_name=='south':
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'annual',
                                                   'L1_BC_'+mask_name+'_' + var_name + '_' + str(year)+'_unbalanced')
                    else:
                        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'annual',
                                                   'L1_BC_' + mask_name + '_' + var_name + '_' + str(year))

                if not os.path.exists(output_file):

                    if print_level >= 2:
                        print('        - Combining files in year ' + str(year))

                    if print_level >= 3:
                        print('            - Getting a list of monthly files')
                    dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(mask_name, var_name, Nr, n_rows_L1, n_cols_L1, year)

                    if print_level >= 3:
                        print('            - Stacking all of the monthly files into a big global file')
                    output_grid = stack_monthly_bc_files_to_one(config_dir, L1_model_name, mask_name, var_name,
                                                              AngleCS, AngleSN, dest_files, dest_file_shapes, total_timesteps, print_level, read_from_unbalanced)

                    if print_level >= 2:
                        print('        - Outputting to ' + output_file)
                    output_grid.ravel('C').astype('>f4').tofile(output_file)

                else:

                    if print_level >= 2:
                        print('        - Skipping '+str(var_name)+' file in year ' + str(year)+' because it already exists')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--L1_model_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="L1_model_name",
                        type=str, required=True)

    parser.add_argument("-i", "--proc_id", action="store",
                        help="The id of the process to run.", dest="proc_id",
                        type=int, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year

    combine_and_rotate_L1_monthly_bcs(config_dir, L1_model_name, proc_id,
                                    start_year, final_year, print_level=4)