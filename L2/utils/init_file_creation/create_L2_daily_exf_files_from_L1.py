
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys

def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    if var_name=='UWIND':
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name=='VWIND':
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    hfac = hfac[:1,:,:]
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(XC,YC,AngleCS,AngleSN,mask)

def create_annual_list_of_files(year,var_name):
    if year % 4 == 0:
        n_days = 366
    else:
        n_days = 365

    day = 0
    month = 1

    output_files = []
    start_indices = []
    end_indices = []

    for d in range(n_days):
        day+=1
        if var_name not in ['UWIND','VWIND']:
            output_file = 'L2_exf_'+var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.bin'
        else:
            output_file = 'L2_exf_' + var_name + '_' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(day) + '_rotated.bin'
        output_files.append(output_file)
        start_indices.append(d*4)
        end_indices.append((d+1)*4)

        if month in [1,3,5,7,8,10,12]:
            if day==31:
                day = 0
                month += 1
        elif month in [4,6,9,11]:
            if day == 30:
                day = 0
                month += 1
        else:
            if year%4==0:
                if day==29:
                    day = 0
                    month += 1
            else:
                if day==28:
                    day = 0
                    month += 1

    return(output_files,start_indices,end_indices)

def read_exf_variable_from_L1(L1_exf_dir, var_name, year,
                              L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_wet_grid, print_level):

    L1_file_name = 'L1_exf_'+var_name+'_'+str(year)

    n_rows = np.shape(L1_XC)[0]
    n_cols = np.shape(L1_YC)[1]

    if var_name not in ['UWIND','VWIND']:
        if print_level>=3:
            print('            - Reading '+L1_file_name)
        exf_field = np.fromfile(os.path.join(L1_exf_dir,L1_file_name),'>f4')
        n_timesteps = int(np.size(exf_field) / (n_rows * n_cols))
        exf_field = np.reshape(exf_field, (n_timesteps, n_rows, n_cols))
    else:
        if var_name=='UWIND':
            u_file = L1_file_name
            v_file = L1_file_name.replace('UWIND','VWIND')
        if var_name=='VWIND':
            v_file = L1_file_name
            u_file = L1_file_name.replace('VWIND', 'UWIND')

        if print_level>=3:
            print('            - Reading '+u_file)
        u_exf_field = np.fromfile(os.path.join(L1_exf_dir, u_file), '>f4')
        n_timesteps = int(np.size(u_exf_field) / (n_rows * n_cols))
        u_exf_field = np.reshape(u_exf_field, (n_timesteps, n_rows, n_cols))

        if print_level>=3:
            print('            - Reading '+v_file)
        v_exf_field = np.fromfile(os.path.join(L1_exf_dir, v_file), '>f4')
        v_exf_field = np.reshape(v_exf_field, (n_timesteps, n_rows, n_cols))

        exf_field = np.zeros_like(u_exf_field)
        for i in range(np.shape(u_exf_field)[0]):
            if var_name=='UWIND':
                exf_field[i,:,:] = L1_AngleCS * u_exf_field[i,:,:] - L1_AngleSN * v_exf_field[i,:,:]
            if var_name == 'VWIND':
                exf_field[i,:,:] = L1_AngleSN * u_exf_field[i,:,:] + L1_AngleCS * v_exf_field[i,:,:]

        del u_exf_field
        del v_exf_field

    if print_level>=4:
        print('                - the L1 grid for this file will have '+str(n_timesteps)+' timesteps')

    ##############################################################

    points = np.column_stack([L1_XC.ravel(),L1_YC.ravel()])
    values = np.zeros((n_timesteps,np.size(L1_XC)))
    for timestep in range(n_timesteps):
        values[timestep,:] = exf_field[timestep,:,:].ravel()
    mask = L1_wet_grid.ravel()

    points = points[mask!=0,:]
    values = values[:,mask!=0]
    mask = mask[mask!=0]

    return(points,values,mask)

def downscale_L1_exf_field_to_L2(df, var_name,
                                 L1_points, L1_values,
                                 L1_wet_grid_points, L1_wet_grid_on_L2_points,
                                 L2_XC, L2_YC, L2_wet_grid, print_level):

    nTimestepsOut = np.shape(L1_values)[0]

    L2_exf_var = np.zeros((nTimestepsOut, np.shape(L2_YC)[0], np.shape(L2_YC)[1]))

    if print_level >= 5:
        print('                -  Variable shapes entering downscale routine:')
        print('                    -  L1_point: '+str(np.shape(L1_points)))
        print('                    -  L1_values: ' + str(np.shape(L1_values)))
        print('                    -  L1_wet_grid: ' + str(np.shape(L1_wet_grid_points)))
        # print('          -  L1_wet_grid_on_L2_subset: ' + str(np.shape(L1_wet_grid_on_L2_3D_points)))
        print('                    -  L2_XC: ' + str(np.shape(L2_XC)))
        print('                    -  L2_YC: ' + str(np.shape(L2_YC)))
        print('                    -  L2_wet_grid: ' + str(np.shape(L2_wet_grid)))

    for timestep in range(nTimestepsOut):
        # if timestep%50==0:
        if print_level >= 5:
            print('                    - Working on timestep '+str(timestep)+' of '+str(nTimestepsOut))

        downscaled_field = df.downscale_2D_points_with_zeros(L1_points, L1_values[timestep, :],
                                                  L1_wet_grid_points,
                                                  L2_XC, L2_YC, L2_wet_grid[0,:,:])
        L2_exf_var[timestep, :, :] = downscaled_field

    # plt.imshow(L2_exf_var[0,:,:],origin='lower')
    # plt.show()

    return(L2_exf_var)


########################################################################################################################


def create_exf_field(config_dir, L1_model_name, L2_model_name,
                     var_name, years, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','exf'))
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input','exf')

    if var_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,var_name))

    if print_level>=2:
        print('        - Reading in the geometry of the L1 domain')

    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_wet_grid = read_grid_geometry_from_nc(config_dir, L1_model_name, var_name)
    # # for face in faces:
    # #     plt.subplot(1,2,1)
    # #     C = plt.imshow(L1_XC_faces[face],origin='lower')
    # #     plt.colorbar(C)
    # #     plt.subplot(1, 2, 2)
    # #     C = plt.imshow(L1_YC_faces[face], origin='lower')
    # #     plt.colorbar(C)
    # #     plt.show()

    L2_XC, L2_YC, _, _, L2_wet_grid = read_grid_geometry_from_nc(config_dir,L2_model_name,var_name)

    for year in years:
        if print_level >= 2:
            print('        - Downscaling the timesteps for year ' + str(year))

        output_files, start_indices, end_indices = create_annual_list_of_files(year, var_name)

        all_files_completed = True
        for output_file_name in output_files:
            if output_file_name not in os.listdir(os.path.join(output_dir,var_name)):
                all_files_completed = False

        if not all_files_completed:

            if print_level >= 3:
                print('            - Reading in the L1 file')

            L1_exf_dir = os.path.join(config_dir,'L1',L1_model_name,'input','exf')
            L1_points, L1_values, L1_wet_grid_points = read_exf_variable_from_L1(L1_exf_dir, var_name, year,
                                                                                 L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_wet_grid,
                                                                                 print_level)


            for ff in range(len(output_files)):
                output_file_name = output_files[ff]
                if output_file_name not in os.listdir(os.path.join(output_dir, var_name)):
                    start_index = start_indices[ff]
                    end_index = end_indices[ff]

                    if print_level >= 4:
                        print('                - Creating file '+output_file_name+' (indices '+str(start_index)+' to '+str(end_index)+')')

                    L1_values_subset = L1_values[start_index:end_index,:]

                    L1_wet_grid_on_L2_points = np.copy(L1_wet_grid_points)

                    if print_level >= 5:
                        print('                    - Downscaling the output to the new boundary')
                    L2_exf_var = downscale_L1_exf_field_to_L2(df, var_name,
                                                              L1_points, L1_values_subset,
                                                              L1_wet_grid_points, L1_wet_grid_on_L2_points,
                                                              L2_XC, L2_YC, L2_wet_grid, print_level)

                    output_file = os.path.join(output_dir, var_name, output_file_name)
                    L2_exf_var.ravel(order='C').astype('>f4').tofile(output_file)
                else:
                    if print_level >= 4:
                        print('                - Skipping file '+output_file_name+' (already created)')

        else:
            if print_level >= 3:
                print('            - Skipping ' + str(year) + ' - all daily files completed for this year')

def create_exf_fields_via_interpolation(config_dir, L2_model_name, L1_model_name, proc_id,
                                        start_year, final_year, print_level):

    var_name_list = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP','RUNOFF','ATMOSCO2','IRONDUST']
    var_name = var_name_list[proc_id % len(var_name_list)]

    if print_level>=1:
        print('    - Creating the exf fields for ' + var_name + ' to cover years ' +
          str(start_year) + ' to ' + str(final_year))
    years = np.arange(start_year,final_year+1).tolist()

    if print_level >= 1:
        print('    - Running the downscale routine formatted as grid:')
    create_exf_field(config_dir, L1_model_name, L2_model_name,
                     var_name, years, print_level)

