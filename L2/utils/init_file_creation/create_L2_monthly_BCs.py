
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys

def create_src_dest_dicts_from_ref(config_dir, L2_model_name, mask_name, var_name, start_year, final_year, start_month, final_month):
    dest_files = []
    start_date = datetime(start_year, start_month, 1)
    final_date = datetime(final_year, final_month, 28)
    for year in range(1992,2022):
        for month in range(1, 13):
            test_date = datetime(year, month, 15)
            if test_date >= start_date and test_date <= final_date:
                dest_files.append(str(year) + '{:02d}'.format(month))

    f = open(os.path.join(config_dir,'L2',L2_model_name, 'input', L2_model_name+'_BC_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UVEL','VVEL','UICE','VICE']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    for dest_file in dest_files:
        dest_files_out.append('L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(mask_name+'_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets['L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    if var_name in ['UVEL','UICE']:
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name in ['VVEL','VICE']:
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(XC,YC,AngleCS,AngleSN,mask,delR)

def read_mask_reference_from_nc_dict(nc_dict_file,model_name,mask_name):

    if model_name=='L2_Santa_Barbara':
        if mask_name=='north':
            read_boundary='south'
        elif mask_name=='south':
            read_boundary='north'
        else:
            read_boundary = mask_name
    else:
        read_boundary = mask_name

    ds = nc4.Dataset(nc_dict_file)
    if model_name == 'L2_Disko_Bay':
        grp = ds.groups['L3_'+read_boundary]
    if model_name == 'L2_Santa_Barbara':
        grp = ds.groups['L2_SB_' + read_boundary]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]
    ds.close()
    return(source_rows,source_cols)

def read_L1_boundary_variable_points(L1_run_dir, L2_model_name, boundary, var_name,
                                     source_files, source_file_read_indices,
                                     boundary_Nr, rows, cols,
                                     L1_XC, L1_YC, L1_AngleCS, L1_AngleSN,
                                     L1_wet_grid_3D, print_level):
    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0] + 1
        # n_Timesteps += index_set[1] - index_set[0]

    print('           + the L1 grid for this file will have ' + str(n_Timesteps) + ' timesteps')

    # make a blank grid of zeros
    points = np.zeros((np.size(rows), 2))
    hfac_points = np.zeros((boundary_Nr,np.size(rows)))
    if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
        values = np.zeros((n_Timesteps, 1, np.size(rows)))
    else:
        values = np.zeros((n_Timesteps, boundary_Nr, np.size(rows)))

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level >= 4:
            print('                - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]) + ' from file ' + source_file)

        start_file_index = index_set[0]
        end_file_index = index_set[1]

        if print_level >= 4:
            print('                - Storing at points ' + str(index_counter) + ' to ' + str(
            index_counter + (end_file_index - start_file_index)) + ' in the grid')

        N = len(rows)


        if L2_model_name in ['L2_Disko_Bay']:
            if var_name in ['UVEL','VVEL']:
                u_var_file = os.path.join(L1_run_dir, 'dv', L2_model_name,'L2_'+boundary, 'L2_' + boundary + '_BC_mask_' + 'UVEL.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L1_run_dir, 'dv', L2_model_name,'L2_'+boundary, 'L2_' + boundary + '_BC_mask_' +'VVEL.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            elif var_name in ['UICE','VICE']:
                u_var_file = os.path.join(L1_run_dir, 'dv', L2_model_name,'L2_'+boundary, 'L2_' + boundary + '_BC_mask_' +'UICE.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L1_run_dir, 'dv', L2_model_name,'L2_'+boundary, 'L2_' + boundary + '_BC_mask_' +'VICE.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            else:
                var_file = os.path.join(L1_run_dir, 'dv', L2_model_name,'L2_'+boundary, 'L2_'+ boundary + '_BC_mask_' +var_name +'.' + file_suffix)
                var_grid = np.fromfile(var_file, dtype='>f4')
        if L2_model_name == 'L2_Santa_Barbara':
            if var_name in ['UVEL','VVEL']:
                u_var_file = os.path.join(L1_run_dir, 'dv','L2_SB_'+boundary, 'L2_SB_' + boundary + '_UVEL.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L1_run_dir, 'dv', 'L2_SB_'+boundary, 'L2_SB_' + boundary + '_VVEL.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            elif var_name in ['UICE','VICE']:
                u_var_file = os.path.join(L1_run_dir, 'dv', 'L2_SB_'+boundary, 'L2_SB_' + boundary + '_UICE.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L1_run_dir, 'dv', 'L2_SB_'+boundary, 'L2_SB_' + boundary + '_VICE.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            else:
                var_file = os.path.join(L1_run_dir, 'dv', 'L2_SB_'+boundary, 'L2_SB_'+ boundary + '_' +var_name +'.' + file_suffix)
                var_grid = np.fromfile(var_file, dtype='>f4')

        if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
            if var_name in ['UICE','VICE']:
                timesteps_in_file = int(np.size(u_var_grid) / (N))
                u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
                u_var_grid = u_var_grid[start_file_index:end_file_index+1, :]
                # u_var_grid = u_var_grid[start_file_index:end_file_index, :]
                v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
                v_var_grid = v_var_grid[start_file_index:end_file_index+1, :]
                # v_var_grid = v_var_grid[start_file_index:end_file_index, :]
            else:
                timesteps_in_file = int(np.size(var_grid) / (N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, N))
                var_grid = var_grid[start_file_index:end_file_index+1, :]
                # var_grid = var_grid[start_file_index:end_file_index, :]

            for n in range(N):
                points[n, 0] = L1_XC[rows[n], cols[n]]
                points[n, 1] = L1_YC[rows[n], cols[n]]
                hfac_points[:, n] = L1_wet_grid_3D[:, rows[n], cols[n]]
                if var_name in ['UICE', 'VICE']:
                    angle_cos = L1_AngleCS[rows[n], cols[n]]
                    angle_sin = L1_AngleSN[rows[n], cols[n]]
                    if var_name=='UICE':
                        zonal_velocity = angle_cos * u_var_grid[:, n] - angle_sin * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1, 0, n] = zonal_velocity
                        # values[index_counter:index_counter + (end_file_index - start_file_index),n] = zonal_velocity
                    if var_name=='VICE':
                        meridional_velocity = angle_sin * u_var_grid[:, n] + angle_cos * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1,0, n] = meridional_velocity
                        # values[index_counter:index_counter + (end_file_index - start_file_index),n] = meridional_velocity
                else:
                    values[index_counter:index_counter + (end_file_index - start_file_index)+1, 0, n] = var_grid[:, n]
                    # values[index_counter:index_counter + (end_file_index - start_file_index), n] = var_grid[:,n]
        else:
            if var_name in ['UVEL','VVEL']:
                timesteps_in_file = int(np.size(u_var_grid) / (boundary_Nr * N))
                u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, boundary_Nr, N))
                u_var_grid = u_var_grid[start_file_index:end_file_index+1, :, :]
                # u_var_grid = u_var_grid[start_file_index:end_file_index, :, :]
                v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, boundary_Nr, N))
                v_var_grid = v_var_grid[start_file_index:end_file_index+1, :, :]
                # v_var_grid = v_var_grid[start_file_index:end_file_index, :, :]
            else:
                timesteps_in_file = int(np.size(var_grid) / (boundary_Nr * N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, boundary_Nr, N))
                var_grid = var_grid[start_file_index:end_file_index+1, :, :]
                # var_grid = var_grid[start_file_index:end_file_index, :, :]
            for n in range(N):
                points[n, 0] = L1_XC[rows[n], cols[n]]
                points[n, 1] = L1_YC[rows[n], cols[n]]
                hfac_points[:, n] = L1_wet_grid_3D[:, rows[n], cols[n]]
                if var_name in ['UVEL', 'VVEL']:
                    angle_cos = L1_AngleCS[rows[n], cols[n]]
                    angle_sin = L1_AngleSN[rows[n], cols[n]]
                    if var_name=='UVEL':
                        zonal_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                        for k in range(np.shape(zonal_velocity)[1]):
                            zonal_velocity[:,k] = angle_cos * u_var_grid[:, k, n] - angle_sin * v_var_grid[:, k, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1, :boundary_Nr,n] = zonal_velocity
                        # values[index_counter:index_counter + (end_file_index - start_file_index), :boundary_Nr, n] = zonal_velocity
                    if var_name=='VVEL':
                        meridional_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                        for k in range(np.shape(meridional_velocity)[1]):
                            meridional_velocity[:,k] = angle_sin * u_var_grid[:, k, n] + angle_cos * v_var_grid[:, k, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1,:boundary_Nr,n] = meridional_velocity
                        # values[index_counter:index_counter + (end_file_index - start_file_index), :boundary_Nr, n] = meridional_velocity
                else:
                    values[index_counter:index_counter + (end_file_index - start_file_index)+1, :boundary_Nr, n] = var_grid[:, :, n]
                    # values[index_counter:index_counter + (end_file_index - start_file_index), :boundary_Nr, n] = var_grid[:, :, n]

        index_counter += (end_file_index - start_file_index)+1
        # index_counter += (end_file_index - start_file_index)

    return(points,values,hfac_points)

def downscale_L1_bc_field_to_L2(df, mask_name, var_name,
                                L2_points, L2_values,
                                L2_wet_grid_3D_points, L2_wet_grid_on_L2_3D_points,
                                L2_XC_subset, L2_YC_subset, L2_wet_grid_subset,print_level):

    nTimestepsOut = np.shape(L2_values)[0]

    if mask_name=='south' or mask_name=='north':
        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            L2_boundary_var = np.zeros((nTimestepsOut, 1, np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))
        else:
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_wet_grid_subset)[0], np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))

    if mask_name=='west' or mask_name=='east':
        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            L2_boundary_var = np.zeros((nTimestepsOut, 1, np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))
        else:
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_wet_grid_subset)[0], np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))

    if print_level>=5:
        print('                -  Variable shapes entering downscale routine:')
        print('                    -  L2_point: '+str(np.shape(L2_points)))
        print('                    -  L2_values: ' + str(np.shape(L2_values)))
        print('                    -  L2_wet_grid: ' + str(np.shape(L2_wet_grid_3D_points)))
        # print('          -  L2_wet_grid_on_L2_subset: ' + str(np.shape(L2_wet_grid_on_L2_3D_points)))
        print('                    -  L2_XC_subset: ' + str(np.shape(L2_XC_subset)))
        print('                    -  L2_YC_subset: ' + str(np.shape(L2_YC_subset)))
        print('                    -  L2_wet_grid_subset: ' + str(np.shape(L2_wet_grid_subset)))

    for timestep in range(nTimestepsOut):
        if timestep%50==0:
            if print_level >= 5:
                print('                  - Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+mask_name)
                print('                  - The bounds at this timestep are '+str(np.min(L2_values[timestep, :, :]))+' to '+str(np.max(L2_values[timestep, :, :])))

        if var_name in ['UVEL','VVEL','THETA','SALT','UICE','VICE']:
            downscaled_field = df.downscale_3D_points(L2_points, L2_values[timestep, :, :],
                                                       L2_wet_grid_3D_points,
                                                       L2_XC_subset, L2_YC_subset, L2_wet_grid_subset,
                                                       remove_zeros=True, printing = False)
        else:
            downscaled_field = df.downscale_3D_points_with_zeros(L2_points, L2_values[timestep, :, :],
                                                                 L2_wet_grid_3D_points,
                                                                 L2_XC_subset, L2_YC_subset, L2_wet_grid_subset,
                                                                 printing=False)
        L2_boundary_var[timestep, :, :, :] = downscaled_field

    # if var_name == 'ETAN':
    #     plt.subplot(2,1,1)
    #     C = plt.imshow(L1_boundary_var_subset[0,:,:])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2,1,2)
    #     plt.imshow(L2_boundary_var[0,:,:])
    #     plt.show()
    # else:
    #     plt.subplot(2, 1, 1)
    #     C = plt.imshow(L1_boundary_var_subset[0, 0, :, :])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2, 1, 2)
    #     plt.imshow(L2_boundary_var[0, 0, :, :])
    #     plt.show()

    return(L2_boundary_var)


########################################################################################################################


def create_bc_field(config_dir, L1_model_name, L2_model_name, mask_name,
                    boundary_var_name, dest_files, source_file_read_dict, source_file_read_index_sets, print_level, write_to_unbalanced):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    # this is the dir where the obcs output will be stored
    if 'obcs' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','obcs'))
    if write_to_unbalanced:
        if 'unbalanced' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input','obcs')):
            os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','obcs','unbalanced'))
        output_dir = os.path.join(config_dir,'L2',L2_model_name,'input','obcs','unbalanced')
    else:
        output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'obcs')

    if mask_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,mask_name))
    if boundary_var_name not in os.listdir(os.path.join(output_dir,mask_name)):
        os.mkdir(os.path.join(output_dir,mask_name,boundary_var_name))

    if print_level >= 1:
        print('    - Reading in the geometry of the L1 domain')
    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_wet_grid_3D, delR_in = read_grid_geometry_from_nc(config_dir, L1_model_name, boundary_var_name)
    Nr_in = len(delR_in)

    if print_level >= 1:
        print('    - Reading in the geometry of the L2 domain')
    L2_XC, L2_YC, L2_AngleCS, L2_AngleSN, L2_wet_grid_3D, delR_out = read_grid_geometry_from_nc(config_dir, L2_model_name, boundary_var_name)
    Nr_out = len(delR_out)

    if mask_name=='north':
        if boundary_var_name=='VVEL':
            L2_XC_subset = L2_XC[-2:-1, :]
            L2_YC_subset = L2_YC[-2:-1, :]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:,-2:-1,:]
        else:
            L2_XC_subset = L2_XC[-1:, :]
            L2_YC_subset = L2_YC[-1:, :]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, -1:, :]
    if mask_name=='south':
        if boundary_var_name=='VVEL':
            L2_XC_subset = L2_XC[1:2, :]
            L2_YC_subset = L2_YC[1:2, :]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:,1:2,:]
        else:
            L2_XC_subset = L2_XC[:1, :]
            L2_YC_subset = L2_YC[:1, :]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, :1, :]
    if mask_name=='east':
        if boundary_var_name == 'UVEL':
            L2_XC_subset = L2_XC[:, -2:-1]
            L2_YC_subset = L2_YC[:, -2:-1]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, :, -2:-1]
        else:
            L2_XC_subset = L2_XC[:, -1:]
            L2_YC_subset = L2_YC[:, -1:]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, :, -1:]
    if mask_name=='west':
        if boundary_var_name=='UVEL':
            L2_XC_subset = L2_XC[:, 1:2]
            L2_YC_subset = L2_YC[:, 1:2]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, :,1:2]
        else:
            L2_XC_subset = L2_XC[:, :1]
            L2_YC_subset = L2_YC[:, :1]
            L2_wet_grid_3D_subset = L2_wet_grid_3D[:, :, :1]

    if boundary_var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
        L2_wet_grid_3D_subset = L2_wet_grid_3D_subset[:1,:,:]

    print('    - Reading the mask to reference the variable to the llc grid')
    if L2_model_name=='L2_Santa_Barbara':
        nc_dict_file = os.path.join(config_dir,'L1',L1_model_name, 'input','L1_dv_mask_reference_dict_SB.nc')
    else:
        nc_dict_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'L1_dv_mask_reference_dict.nc')
    source_rows, source_cols = read_mask_reference_from_nc_dict(nc_dict_file, L2_model_name, mask_name)

    for dest_file in dest_files:
        if dest_file not in os.listdir(os.path.join(output_dir,mask_name,boundary_var_name)):

            print('    - Downscaling the timesteps to be stored in file ' + str(dest_file))
            source_files = source_file_read_dict[dest_file]

            if len(source_files)>0:
                # try:

                source_file_read_indices = source_file_read_index_sets[dest_file]

                if print_level >= 4:
                    print('                - Reading in the L1 diagnostics_vec output')
                L1_run_dir = os.path.join(config_dir,'L1',L1_model_name, 'run')
                L1_boundary_points, L1_boundary_values, L1_boundary_point_hFacC = \
                    read_L1_boundary_variable_points(L1_run_dir, L2_model_name, mask_name, boundary_var_name,
                                                     source_files, source_file_read_indices,
                                                     Nr_in, source_rows, source_cols,
                                                     L1_XC, L1_YC, L1_AngleCS, L1_AngleSN,
                                                     L1_wet_grid_3D,
                                                     print_level)

                if print_level >= 4:
                    print('                - L1 points shape: '+str(np.shape(L1_boundary_points)))
                    print('                - L1 values shape: ' + str(np.shape(L1_boundary_values)))
                    print('                     - NaN values: '+str(np.sum(np.isnan(L1_boundary_values))))
                    print('                - L1 hfac shape: ' + str(np.shape(L1_boundary_point_hFacC)))


                    # plt.imshow(L1_boundary_values[10,:,:])
                    # plt.show()

                # plt.plot(L2_XC_subset,L2_YC_subset,'g.')
                # plt.plot(L1_boundary_points[:,0],L1_boundary_points[:,1],'k.')
                # plt.show()

                if boundary_var_name in ['THETA', 'SALT', 'UVEL', 'VVEL'] or 'PTRACE' in boundary_var_name:
                    if Nr_out!=Nr_in:
                        L1_boundary_values, L1_boundary_point_hFacC = df.interpolate_var_points_timeseries_to_new_depth_levels(
                            L1_boundary_values, L1_boundary_point_hFacC, delR_in, delR_out)

                L1_wet_grid_3D_points = np.copy(L1_boundary_point_hFacC)
                L1_wet_grid_on_L2_3D_points = np.copy(L2_wet_grid_3D_subset) # temporary

                print('        - Downscaling the output to the new boundary')
                L2_boundary_var = downscale_L1_bc_field_to_L2(df, mask_name, boundary_var_name,
                                                              L1_boundary_points, L1_boundary_values,
                                                              L1_wet_grid_3D_points, L1_wet_grid_on_L2_3D_points,
                                                              L2_XC_subset, L2_YC_subset, L2_wet_grid_3D_subset, print_level)

                # if 'north' in mask_name or 'south' in mask_name:
                #     plt.subplot(1,2,1)
                #     plt.imshow(L2_wet_grid_3D_subset[:,0,:]==0,vmin=-0.1,vmax=0.1,cmap='seismic_r')
                #     plt.title('Mask')
                #     plt.subplot(1, 2, 2)
                #     plt.imshow(L2_boundary_var[0, :, 0, :],cmap='viridis')#,vmin=-0.1,vmax=0.1
                #     plt.title('L2 (nan values: '+str(np.sum(np.isnan(L2_boundary_var[0,:,0,:]))))
                #     plt.show()
                #
                # if 'east' in mask_name or 'west' in mask_name:
                #     plt.subplot(1, 2, 1)
                #     plt.imshow(L2_wet_grid_3D_subset[:, :, 0] == 0, vmin=-0.1, vmax=0.1, cmap='seismic_r')
                #     plt.title('Mask')
                #     plt.subplot(1, 2, 2)
                #     plt.imshow(L2_boundary_var[0, :, :, 0], cmap='viridis')  # ,vmin=-0.1,vmax=0.1
                #     plt.title('L2 (nan values: ' + str(np.sum(np.isnan(L2_boundary_var[0, :, :, 0]))))
                #     plt.show()

                if mask_name not in os.listdir(output_dir):
                    os.mkdir(os.path.join(output_dir, mask_name))
                if boundary_var_name not in os.listdir(os.path.join(output_dir, mask_name)):
                    os.mkdir(os.path.join(output_dir, mask_name, boundary_var_name))

                if print_level>=4:
                    print('                 - Output shape: '+str(np.shape(L2_boundary_var)))

                # if boundary_var_name in ['UVEL','VVEL','UICE','VICE']:
                #     output_file = os.path.join(output_dir, mask_name, boundary_var_name, dest_file[:-4]+'_rotated.bin')
                # else:
                output_file = os.path.join(output_dir, mask_name, boundary_var_name, dest_file)

                if print_level>=4:
                    print('                 - Output file: '+str(np.shape(L2_boundary_var)))

                L2_boundary_var.ravel(order='C').astype('>f4').tofile(output_file)

            else:
                print('    - Skipping ' + str(dest_file) + ' because it was already created')

def create_bc_fields_via_interpolation(config_dir, L1_model_name, L2_model_name, boundaries, proc_id,
                                       start_year, final_year, start_month,
                                       final_month, print_level,write_to_unbalanced = True):

    if L2_model_name in ['L2_Disko_Bay']:
        var_names = ['THETA','SALT','UVEL','VVEL','UICE','VICE','HSNOW','HEFF','AREA']
    else:
        var_names = ['THETA', 'SALT', 'UVEL', 'VVEL']
    for i in range(1,32):
        var_names.append('PTRACE'+'{:02d}'.format(i))

    var_name_list = []
    mask_name_list = []
    for var_name in var_names:
        for boundary in boundaries:
            var_name_list.append(var_name)
            mask_name_list.append(boundary)

    var_name = var_name_list[proc_id % len(var_name_list)]
    mask_name = mask_name_list[proc_id % len(var_name_list)]

    print('Creating the bc field for ' + var_name + ' on mask ' +mask_name+ ' to cover year/month/days ' +
          str(start_year)+'/'+str(start_month) + ' to ' +
          str(final_year)+'/'+str(final_month))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L2_model_name, mask_name, var_name, start_year, final_year, start_month, final_month)

    print('  Running the Downscale routine:')
    create_bc_field(config_dir, L1_model_name,  L2_model_name, mask_name, var_name, dest_files,
                    source_file_read_dict, source_file_read_index_sets, print_level, write_to_unbalanced)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO files are stored.", dest="ecco_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--L2_model_name", action="store",
                        help="The directory where the L2, L2, and L2 configurations are stored.", dest="L2_model_name",
                        type=str, required=True)

    parser.add_argument("-b", "--boundaries", action="store",
                        help="List of boundaries for this model (e.g. north, south, east, and/or west)",
                        dest="boundaries",
                        type=str, required=True, nargs='+')

    parser.add_argument("-i", "--proc_id", action="store",
                        help="The id of the process to run.", dest="proc_id",
                        type=int, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-s", "--start_month", action="store",
                        help="The start month.", dest="start_month",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    parser.add_argument("-f", "--final_month", action="store",
                        help="The final ymonth.", dest="final_month",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir
    L2_model_name = args.L2_model_name
    boundaries = args.boundaries
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month

    create_bc_fields_via_interpolation(config_dir, ecco_dir, L2_model_name, boundaries, proc_id,
                                       start_year, final_year, start_month,
                                       final_month, print_level=4)

