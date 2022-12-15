
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
#from MITgcmutils import mds
#from scipy.interpolate import griddata
import sys
from datetime import datetime
import ast

def read_grid_geometry(config_dir,model_name,var_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]

    if var_name in ['VWIND']:
        hFac = 'S'
    elif var_name in ['UWIND']:
        hFac = 'W'
    else:
        hFac = 'C'

    HFac = ds.variables['HFac' + hFac][:, :, :]
    if hFac == 'W':
        HFac = HFac[:, :, :-1]
    if hFac == 'S':
        HFac = HFac[:, :-1, :]

    delR = ds.variables['drF'][:]

    ds.close()

    return(XC, YC, AngleCS, AngleSN, HFac, delR)

def read_interpolation_grid(config_dir, model_name, var_name):

    interpolation_grid_file = model_name+'_interpolation_grid.nc'

    if interpolation_grid_file not in os.listdir(os.path.join(config_dir,'nc_grids')):

        raise ValueError('Need to generate the interpolation grid using one of the pickup scripts')

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',interpolation_grid_file))
    if var_name in ['Vvel','SIvice','siVICE','VWIND','VSTRESS']:
        grp = ds.groups['S']
    elif var_name in ['Uvel','SIuice','siUICE','UWIND','USTRESS']:
        grp = ds.groups['W']
    else:
        grp = ds.groups['C']

    interpolation_type_grid = grp.variables['interp_type'][:,:,:]
    source_rows_grid = grp.variables['source_rows'][:, :, :]
    source_cols_grid = grp.variables['source_cols'][:, :, :]
    source_levels_grid = grp.variables['source_levels'][:,:,:]

    ds.close()

    return(interpolation_type_grid, source_rows_grid, source_cols_grid, source_levels_grid)

def create_src_dest_dicts_from_ref(config_dir, L1_model_name, var_name, start_year, final_year, start_month, final_month):
    dest_files = []
    start_date = datetime(start_year, start_month, 1)
    final_date = datetime(final_year, final_month, 28)
    for year in range(1992,2022):
        for month in range(1, 13):
            test_date = datetime(year, month, 15)
            if test_date >= start_date and test_date <= final_date:
                dest_files.append(str(year) + '{:02d}'.format(month))

    f = open(os.path.join(config_dir,'L1',L1_model_name, 'input', L1_model_name+'_exf_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UWIND', 'VWIND']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    for dest_file in dest_files:

        dest_files_out.append('L1_exf_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L1_exf_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets['L1_exf_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def read_L0_surface_variable_points(L0_run_dir, model_name, var_name,
                                     source_files, source_file_read_indices,
                                     faces, rows, cols,
                                     ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                     print_level):
    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0]

    # make a blank grid of zeros
    points = np.zeros((np.size(faces), 2))
    hfac_points = np.zeros((1,np.size(faces)))
    values = np.zeros((n_Timesteps, np.size(faces)))

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

        N = len(faces)

        prefix = '_'.join(model_name.split('_')[:2])
        boundary = 'surf'
        if var_name in ['UWIND','VWIND']:
            u_var_file = os.path.join(L0_run_dir, 'dv',model_name+'_' + boundary + '_UWIND.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',model_name+'_' + boundary + '_VWIND.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        elif var_name in ['USTRESS','VSTRESS']:
            u_var_file = os.path.join(L0_run_dir, 'dv',model_name+'_' + boundary + '_USTRESS.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',model_name+'_' + boundary + '_VSTRESS.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        else:
            var_file = os.path.join(L0_run_dir, 'dv',model_name+'_'+ boundary + '_' + var_name +'.' + file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')

        if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
            timesteps_in_file = int(np.size(u_var_grid) / N)
            u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
            u_var_grid = u_var_grid[start_file_index:end_file_index, :]
            v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
            v_var_grid = v_var_grid[start_file_index:end_file_index, :]
        else:
            timesteps_in_file = int(np.size(var_grid) / N)
            var_grid = np.reshape(var_grid, (timesteps_in_file, N))
            var_grid = var_grid[start_file_index:end_file_index, :]

        for n in range(N):
            if faces[n] != 0:
                points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                hfac_points[0, n] = ecco_hFacC_faces[faces[n]][0, rows[n], cols[n]]
                if var_name in ['UWIND','VWIND','USTRESS', 'VSTRESS']:
                    angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                    angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                    if var_name=='UWIND' or var_name=='USTRESS':
                        zonal_velocity = angle_cos * u_var_grid[:, n] - angle_sin * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index), n] = zonal_velocity
                    if var_name=='VWIND' or var_name=='VSTRESS':
                        meridional_velocity = angle_sin * u_var_grid[:, n] + angle_cos * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index), n] = meridional_velocity
                else:
                    values[index_counter:index_counter + (end_file_index - start_file_index), n] = var_grid[:, n]


        index_counter += (end_file_index - start_file_index)

    return(points,values,hfac_points)

def create_exf_field(config_dir, ecco_dir, model_name, var_name,
                     dest_files, source_file_read_dict, source_file_read_index_sets, print_level, use_interpolation_grids = True):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if 'exf' not in os.listdir(os.path.join(config_dir, 'L1',model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf'))
    if var_name not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name))

    if print_level >= 1:
        print('    - Creating the ' + var_name + ' exf files for the ' + model_name + ' model from ECCOv5 output data')

    if print_level >= 1:
        print('    - Reading in the model geometry')
    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_hfac, delR = read_grid_geometry(config_dir, model_name, var_name)
    L1_mask = np.copy(L1_hfac)
    L1_mask[L1_mask>0]=1
    L1_mask = L1_mask[:1,:,:]

    # step 1: get the ecco faces geometry
    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    llc = 270
    ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
        ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    L0_run_dir = os.path.join(config_dir, 'L0', 'run')

    if print_level >= 3:
        print('            - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[model_name+'_surf']
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()

    ####################################################################################################################
    # Loop through the monthly files

    for dest_file in dest_files:
        if dest_file not in []:#os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name)):

            year = int(dest_file.split('.')[1][:4])
            month = int(dest_file.split('.')[1][4:6])

            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            n_timesteps = nDays * 4

            if print_level >= 3:
                print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))

            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            output_grid = np.zeros((n_timesteps, np.shape(L1_XC)[0], np.shape(L1_XC)[1]))

            L0_surface_points, L0_surface_values, L0_surface_hFacC = \
                read_L0_surface_variable_points(L0_run_dir, model_name, var_name,
                                            source_files, source_file_read_indices,
                                            faces, rows, cols,
                                            ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces,
                                            ecco_hFacC_faces,
                                            print_level)

            L0_surface_values = np.reshape(L0_surface_values, (
                np.shape(L0_surface_values)[0], 1, np.shape(L0_surface_values)[1]))
            L0_surface_mask = np.copy(L0_surface_hFacC)
            L0_surface_mask[L0_surface_mask > 0] = 1

            if model_name == 'L1_mac_delta':
                L0_surface_points[:,0] += 360

            # plt.plot(L0_surface_points[:,0],L0_surface_points[:,1],'b.')
            # plt.plot(L1_XC,L1_YC,'k.')
            # plt.show()

            if use_interpolation_grids:
                L1_interpolation_mask, L1_source_rows, L1_source_cols, L1_source_levels = \
                    read_interpolation_grid(config_dir, model_name, var_name)

            for timestep in range(np.shape(L0_surface_values)[0]):
                if timestep%50 == 0:
                    print('        - Downscaling timesteps '+str(timestep)+' to '+str(np.min([timestep+50,np.shape(L0_surface_values)[0]])))

                if print_level >= 5:
                    if timestep == 0:
                        print('                - L0_surface_points shape: ' + str(np.shape(L0_surface_points)))
                        print('                - L0_surface_values shape: ' + str(np.shape(L0_surface_values)))
                        print('                - L0_surface_values range: '+str(np.min(L0_surface_values))+' to '+str(np.max(L0_surface_values)))
                        print('                - L0_surface_mask: ' + str(np.shape(L0_surface_mask)))
                        print('                - L1_XC shape: ' + str(np.shape(L1_XC)))
                        print('                - L1_YC shape: ' + str(np.shape(L1_YC)))
                        print('                - L1_mask shape: ' + str(np.shape(L1_mask)))
                        print('                - output_grid shape: ' + str(np.shape(output_grid)))

                if use_interpolation_grids:
                    interp_field = df.downscale_3D_points_with_interpolation_mask(L0_surface_points,
                                                                                  L0_surface_values[timestep, :, :],
                                                                                  L0_surface_mask,
                                                                                  L1_XC, L1_YC, L1_mask,
                                                                                  L1_interpolation_mask, L1_source_rows,
                                                                                  L1_source_cols, L1_source_levels,
                                                                                  printing=False,
                                                                                  testing=False)
                else:
                    interp_field = df.downscale_3D_points_with_zeros(L0_surface_points,
                                                                     L0_surface_values[timestep, :, :],
                                                                     L0_surface_mask,
                                                                     L1_XC, L1_YC, L1_mask,
                                                                     mean_vertical_difference=0,
                                                                     fill_downward=True,
                                                                     printing=False)

                output_grid[timestep, :, :] = interp_field[0, :, :]

            # plot_grid = output_grid[0, :, :]
            # plt.imshow(output_grid[0, :, :],origin='lower')
            # plt.title(var_name+' (min = '+str(np.min(plot_grid[plot_grid!=0]))+' (max = '+str(np.max(plot_grid[plot_grid!=0]))+')')
            # plt.show()

            output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name, dest_file)
            output_grid.ravel(order='C').astype('>f4').tofile(output_file)


def create_L1_exf_fields(config_dir, ecco_dir, L1_model_name, proc_id,
                      start_year, final_year, start_month, final_month, print_level, use_interpolation_grids = True):

    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP','RUNOFF','ATMOSCO2','IRONDUST']
    var_name = var_names[proc_id % len(var_names)]

    print('Creating the exf field for ' + var_name + ' to cover year/month/days ' +
          str(start_year)+'/'+str(start_month) + ' to ' +
          str(final_year)+'/'+str(final_month))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L1_model_name, var_name, start_year, final_year, start_month, final_month)

    print('  Running the Downscale routine:')
    create_exf_field(config_dir, ecco_dir,  L1_model_name, var_name,
                     dest_files, source_file_read_dict, source_file_read_index_sets, print_level,use_interpolation_grids)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO files are stored.", dest="ecco_dir",
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
    L1_model_name = args.L1_model_name
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month

    create_L1_exf_fields(config_dir, ecco_dir, L1_model_name, proc_id,
                         start_year, final_year, start_month,
                         final_month, print_level=4)


