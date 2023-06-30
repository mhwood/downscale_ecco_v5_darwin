
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import shapefile
import netCDF4 as nc4
import argparse
import ast
import sys

def create_dest_file_list(mask_name, var_name, start_year, final_year):

    dest_files = []
    start_date = datetime(start_year, 1, 1)
    final_date = datetime(final_year, 12, 31)
    years = []
    for year in range(1992,2022):
        for month in range(1, 13):
            test_date = datetime(year, month, 15)
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

def read_grid_information(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    Depth = np.array(ds.variables['Depth'][:, :])
    dxG = np.array(ds.variables['dxG'][:, :])
    dyG = np.array(ds.variables['dyG'][:, :])
    hFacS = np.array(ds.variables['HFacS'][:, :, :])
    hFacW = np.array(ds.variables['HFacW'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return (Depth, dxG, dyG, hFacS, hFacW, delR)

def read_annual_velocity_file(velocity_file_path,mask_name,n_rows,n_cols,Nr):

    grid = np.fromfile(velocity_file_path,'>f4')
    if mask_name in ['north','south']:
        timesteps = int(np.size(grid)/(n_cols*Nr))
        grid = np.reshape(grid, (timesteps, Nr, n_cols))
    if mask_name in ['east','west']:
        timesteps = int(np.size(grid)/(n_rows*Nr))
        grid = np.reshape(grid, (timesteps, Nr, n_rows))

    return(grid)

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

def write_timeseries_to_nc(output_file,var_name,timeseries):
    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time',len(timeseries))
    var = ds.createVariable(var_name,'f4',('time',))
    var[:] = timeseries
    ds.close()

def output_L1_fluxes_to_nc(output_file, Depth, dxG, dyG, hFacS, hFacW, drF,
                           boundaries, all_mean_flux_grids, all_integrated_flux_timeseries):

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('depth_cells',len(drF))
    ds.createDimension('rows',np.shape(Depth)[0])
    ds.createDimension('cols', np.shape(Depth)[1])
    ds.createDimension('rows_p1', np.shape(Depth)[0]+1)
    ds.createDimension('cols_p1', np.shape(Depth)[1]+1)
    ds.createDimension('time',np.size(all_integrated_flux_timeseries[0]))

    depth_var = ds.createVariable('Depth','f4',('rows','cols'))
    depth_var[:, :] = Depth

    dxg_var = ds.createVariable('dxG', 'f4', ('rows_p1', 'cols'))
    dxg_var[:, :] = dxG

    dyg_var = ds.createVariable('dyG', 'f4', ('rows', 'cols_p1'))
    dyg_var[:, :] = dyG

    drf_var = ds.createVariable('drF', 'f4', ('depth_cells', ))
    drf_var[:] = drF

    # hFacW_var = ds.createVariable('hFacW', 'f4', ('depth_cells','rows','cols_p1'))
    # hFacW_var[:, :, :] = hFacW
    #
    # hFacS_var = ds.createVariable('hFacS', 'f4', ('depth_cells','rows_p1','cols'))
    # hFacS_var[:, :, :] = hFacS

    for b in range(len(boundaries)):
        grp = ds.createGroup(boundaries[b])
        if boundaries[b] in ['north','south']:
            mean_flux_var = grp.createVariable('mean_flux','f4',('depth_cells','cols'))
        if boundaries[b] in ['east','west']:
            mean_flux_var = grp.createVariable('mean_flux','f4',('depth_cells','rows'))
        mean_flux_var[:, :] = all_mean_flux_grids[b]

        tvar = grp.createVariable('integrated_timeseries','f4',('time',))
        tvar[:] = all_integrated_flux_timeseries[b]

    ds.close()




########################################################################################################################


def create_timeseries(config_dir, L1_model_name, balanced, boundaries,
                      boundary_var_name, start_year, end_year, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import flux_functions as ff

    if print_level >= 1:
        print('    - Reading in the '+L1_model_name+' geometry')

    if balanced:
        balanced_dir = 'balanced'
    else:
        balanced_dir = 'unbalanced'

    Depth, dxG, dyG, hFacS, hFacW, drF = read_grid_information(config_dir, L1_model_name)
    n_rows = np.shape(Depth)[0]
    n_cols = np.shape(Depth)[1]
    Nr = len(drF)

    all_mean_flux_grids = []
    all_integrated_flux_timeseries = []

    for mask_name in boundaries:
        print('    - Calculating fluxes on the '+mask_name+' boundary')
        years, dest_files_per_year = create_dest_file_list(mask_name, boundary_var_name, start_year, end_year)

        timestep_counter = 0
        if mask_name in ['north', 'south']:
            mean_flux_grid = np.zeros((Nr,n_cols))
        if mask_name in ['west', 'east']:
            mean_flux_grid = np.zeros((Nr,n_rows))

        total_timesteps = 0
        for year in years:
            if year%4==0:
                total_timesteps+=366
            else:
                total_timesteps+=365
        integrated_flux_timeseries = np.zeros((total_timesteps,))

        for year in years:
            print('        - Adding data from year '+str(year))
            if mask_name in ['north','south']:
                velocity_file = 'L1_BC_'+mask_name+'_VVEL_'+str(year)
                velocity_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                  'balanced', mask_name, 'VVEL', velocity_file)
            if mask_name in ['west','east']:
                velocity_file = 'L1_BC_'+mask_name+'_UVEL_'+str(year)
                velocity_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                  'balanced', mask_name, 'UVEL', velocity_file)

            # read in the velocity grid
            velocity_grid = read_annual_velocity_file(velocity_file_path,mask_name,n_rows,n_cols,Nr)
            if print_level>4:
                print('            - Velocity shape: '+str(np.shape(velocity_grid)))

            if not balanced:
                if boundary_var_name!='VEL':
                    # read in the monthly tracer files to the annual file
                    tracer_dir = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                      balanced_dir, mask_name, boundary_var_name)
                    tracer_grid = read_monthly_tracer_files(tracer_dir, dest_files_per_year[year], year, mask_name, n_rows, n_cols, Nr)
                else:
                    tracer_grid = np.copy(velocity_grid)
                    tracer_grid[velocity_grid!=0] = 1
            else:
                if boundary_var_name!='VEL':
                    tracer_file = 'L1_BC_' + mask_name + '_'+boundary_var_name+'_' + str(year)
                    tracer_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                      'balanced', mask_name, boundary_var_name, tracer_file)
                    tracer_grid = read_annual_velocity_file(tracer_file_path,mask_name,n_rows,n_cols,Nr)
                else:
                    tracer_grid = np.copy(velocity_grid)
                    tracer_grid[velocity_grid != 0] = 1

            # calculate the fluxes for this grid
            flux_grid = ff.calculate_flux_into_domain(mask_name, boundary_var_name,
                                       velocity_grid, tracer_grid,
                                       dxG, dyG, drF, hFacS, hFacW, apply_density = False)

            mean_flux_grid += np.sum(flux_grid,axis=0)
            integrated_flux_timeseries[timestep_counter:timestep_counter+np.shape(flux_grid)[0]] = np.sum(np.sum(flux_grid,axis=1),axis=1)
            timestep_counter += np.shape(flux_grid)[0]

        all_mean_flux_grids.append(mean_flux_grid/timestep_counter)
        all_integrated_flux_timeseries.append(integrated_flux_timeseries)

    # plt.plot(integrated_flux_timeseries)
    # plt.show()

    output_file = os.path.join(config_dir,'L1',L1_model_name,'input','obcs','timeseries','L1_timeseries',balanced_dir,
                               L1_model_name+'_'+boundary_var_name+'_boundary_flux.nc')
    output_L1_fluxes_to_nc(output_file, Depth, dxG, dyG, hFacS, hFacW, drF,
                           boundaries, all_mean_flux_grids, all_integrated_flux_timeseries)

def create_annual_timeseries(config_dir, L1_model_name, balanced, boundaries,
                      boundary_var_name, start_year, end_year, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import flux_functions as ff

    if print_level >= 1:
        print('    - Reading in the '+L1_model_name+' geometry')

    if balanced:
        balanced_dir = 'balanced'
    else:
        balanced_dir = 'unbalanced'

    output_dir = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'timeseries', 'L1_timeseries', balanced_dir)


    Depth, dxG, dyG, hFacS, hFacW, drF = read_grid_information(config_dir, L1_model_name)
    n_rows = np.shape(Depth)[0]
    n_cols = np.shape(Depth)[1]
    Nr = len(drF)

    for year in range(start_year, end_year + 1):

        all_mean_flux_grids = []
        all_integrated_flux_timeseries = []

        if str(year) not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, str(year)))
        output_file_name = L1_model_name + '_' + boundary_var_name + '_boundary_flux_' + str(year) + '.nc'

        if output_file_name not in os.listdir(os.path.join(output_dir, str(year))):

            for mask_name in boundaries:
                print('    - Calculating fluxes on the '+mask_name+' boundary')

                _, dest_files_per_year = create_dest_file_list(mask_name, boundary_var_name, year, year+1)

                timestep_counter = 0
                if mask_name in ['north', 'south']:
                    mean_flux_grid = np.zeros((Nr,n_cols))
                if mask_name in ['west', 'east']:
                    mean_flux_grid = np.zeros((Nr,n_rows))

                if year%4==0:
                    total_timesteps=366
                else:
                    total_timesteps=365
                integrated_flux_timeseries = np.zeros((total_timesteps,))

                print('        - Adding data from year '+str(year)+' on the '+mask_name+' boundary')
                if mask_name in ['north','south']:
                    velocity_file = 'L1_BC_'+mask_name+'_VVEL_'+str(year)
                    velocity_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                      'balanced', mask_name, 'VVEL', velocity_file)
                if mask_name in ['west','east']:
                    velocity_file = 'L1_BC_'+mask_name+'_UVEL_'+str(year)
                    velocity_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                      'balanced', mask_name, 'UVEL', velocity_file)

                # read in the velocity grid
                velocity_grid = read_annual_velocity_file(velocity_file_path,mask_name,n_rows,n_cols,Nr)
                if print_level>3:
                    print('            - Velocity shape: '+str(np.shape(velocity_grid)))

                if not balanced:
                    if boundary_var_name!='VEL':
                        # read in the monthly tracer files to the annual file
                        tracer_dir = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                          balanced_dir, mask_name, boundary_var_name)
                        tracer_grid = read_monthly_tracer_files(tracer_dir, dest_files_per_year[year], year, mask_name, n_rows, n_cols, Nr)
                    else:
                        # tracers are 1 for the vel grid because there is no scalar factor when calculating the fluxes
                        tracer_grid = np.copy(velocity_grid)
                        tracer_grid[velocity_grid!=0] = 1
                else:
                    if boundary_var_name!='VEL':
                        tracer_file = 'L1_BC_' + mask_name + '_'+boundary_var_name+'_' + str(year)
                        tracer_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                                          'balanced', mask_name, boundary_var_name, tracer_file)
                        tracer_grid = read_annual_velocity_file(tracer_file_path,mask_name,n_rows,n_cols,Nr)
                    else:
                        tracer_grid = np.copy(velocity_grid)
                        tracer_grid[velocity_grid != 0] = 1

                # calculate the fluxes for this grid
                flux_grid = ff.calculate_flux_into_domain(mask_name, boundary_var_name,
                                           velocity_grid, tracer_grid,
                                           dxG, dyG, drF, hFacS, hFacW, apply_density = False)

                mean_flux_grid += np.sum(flux_grid,axis=0)
                integrated_flux_timeseries[timestep_counter:timestep_counter+np.shape(flux_grid)[0]] = np.sum(np.sum(flux_grid,axis=1),axis=1)
                timestep_counter += np.shape(flux_grid)[0]

                # plt.plot(integrated_flux_timeseries)
                # plt.show()

                all_mean_flux_grids.append(mean_flux_grid/timestep_counter)
                all_integrated_flux_timeseries.append(integrated_flux_timeseries)

            output_file = os.path.join(output_dir, str(year), output_file_name)
            output_L1_fluxes_to_nc(output_file, Depth, dxG, dyG, hFacS, hFacW, drF,
                                   boundaries, all_mean_flux_grids, all_integrated_flux_timeseries)


def create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
                         start_year, final_year, balanced, print_level):

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

    print('  Running the timeseries routine:')
    create_annual_timeseries(config_dir, L1_model_name, balanced, boundaries, var_name, start_year, final_year, print_level)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--L1_model_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="L1_model_name",
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
    L1_model_name = args.L1_model_name
    boundaries = args.boundaries
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month

    create_L1_timeseries(config_dir, L1_model_name, boundaries, proc_id,
                                       start_year, final_year, start_month,
                                       final_month, print_level=4)

