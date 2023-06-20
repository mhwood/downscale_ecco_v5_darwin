
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

def create_src_dest_dicts_from_ref(config_dir, L1_model_name, mask_name, var_name, start_year, final_year):
    dest_files = []
    start_date = datetime(start_year, 1, 1)
    final_date = datetime(final_year, 12, 31)
    for year in range(1992,2022):
        for month in range(1, 13):
            test_date = datetime(year, month, 15)
            if test_date >= start_date and test_date <= final_date:
                dest_files.append(str(year) + '{:02d}'.format(month))

    f = open(os.path.join(config_dir,'L1',L1_model_name, 'input', L1_model_name+'_BC_dest_ref.txt'))
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
        dest_files_out.append('L0_flux_'+mask_name+'_'+var_name+'.'+dest_file+'.nc')
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(mask_name+'_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L0_flux_'+mask_name+'_'+var_name+'.'+dest_file+'.nc'] = source_files
        source_file_read_index_sets['L0_flux_'+mask_name+'_'+var_name+'.'+dest_file+'.nc'] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def get_boundary_Nr_from_manual_list(L1_model_name, boundary):
    if L1_model_name=='L1_W_Greenland':
        if boundary=='north':
            boundary_Nr = 26
        if boundary=='west':
            boundary_Nr = 27
        if boundary=='south':
            boundary_Nr = 50
        if boundary=='east':
            boundary_Nr = 50
    elif L1_model_name=='L1_GOM':
        if boundary=='north':
            boundary_Nr = 5
        if boundary=='south':
            boundary_Nr = 50
        if boundary=='east':
            boundary_Nr = 50
    elif L1_model_name=='L1_mac_delta':
        if boundary=='north':
            boundary_Nr = 50
        if boundary=='west':
            boundary_Nr = 50
        if boundary=='east':
            boundary_Nr = 24
    elif L1_model_name=='L1_East_Pacific':
        if boundary=='north':
            boundary_Nr = 50
        if boundary=='west':
            boundary_Nr = 50
        if boundary=='south':
            boundary_Nr = 50
    else:
        raise ValueError('Need to manually set boundary Nr for this model')

    return(boundary_Nr)

def read_mask_reference_from_nc_dict(nc_dict_file,model_name,mask_name):
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups[model_name+'_'+mask_name]
    source_faces = grp.variables['source_faces'][:]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]
    ds.close()
    return(source_faces, source_rows,source_cols)

def read_flux_ecco_geometry_to_faces(ef, ecco_dir,llc):
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    DXG_faces = {}
    DYG_faces = {}
    XC_faces = {}
    YC_faces = {}
    for i in [1, 2, 3, 4, 5]:
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   3 * llc)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3 * llc,
                                                   llc)
        DXG_face = grid_dict['DXG'].T
        DYG_face = grid_dict['DYG'].T
        DXG_faces[i] = DXG_face
        DYG_faces[i] = DYG_face

        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    hFacS_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init', 'hFacS.data')
    hFacS_faces = ef.read_ecco_field_to_faces(hFacS_path, llc, dim=3)

    hFacW_path = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'input_init', 'hFacW.data')
    hFacW_faces = ef.read_ecco_field_to_faces(hFacW_path, llc, dim=3)

    # plt.imshow(hFacC_faces[5][5,:,:])
    # plt.show()

    return(XC_faces, YC_faces, DXG_faces, DYG_faces, hFacS_faces, hFacW_faces)


def read_flux_points_from_shapefile(config_dir,L1_model_name,mask_name):

    shp_path = os.path.join(config_dir,'L1',L1_model_name,'grid',L1_model_name+'_L0_'+mask_name+'_boundary')
    sf = shapefile.Reader(shp_path)
    records = sf.records()
    points = np.zeros((len(records),2))
    for r in range(len(records)):
        points[r, 0] = records[r][3]
        points[r, 1] = records[r][4]
    return(points)

def read_L0_boundary_variable_points(L0_run_dir, model_name, boundary, var_name,
                                     source_files,source_file_read_indices,
                                     llc, boundary_Nr, Nr, faces, rows, cols,
                                     ecco_XC_faces, ecco_YC_faces, ecco_DXG_faces, ecco_DYG_faces, ecco_hFacS_faces, ecco_hFacW_faces,
                                     print_level):
    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0] + 1
        # n_Timesteps += index_set[1] - index_set[0]

    print('           + the L0 grid for this file will have ' + str(n_Timesteps) + ' timesteps')

    # make a blank grid of zeros
    points = np.zeros((np.size(faces), 4))
    hfac_points = np.zeros((Nr,np.size(faces)))
    if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
        values = np.zeros((n_Timesteps, 1, np.size(faces)))
        vel_values = np.zeros((n_Timesteps, 1, np.size(faces)))
    else:
        values = np.zeros((n_Timesteps, Nr, np.size(faces)))
        vel_values = np.zeros((n_Timesteps, Nr, np.size(faces)))

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

        if model_name == 'L1_W_Greenland':
            if boundary == 'west':
                vel_var_file = os.path.join(L0_run_dir, 'dv', model_name,
                                            model_name + '_' + boundary + '_UVEL.' + file_suffix)
                scalar_for_zonal_or_meridional_orientation = -1
                ecco_hFac_faces = ecco_hFacW_faces
            if boundary == 'south':
                vel_var_file = os.path.join(L0_run_dir, 'dv', model_name,
                                            model_name + '_' + boundary + '_UVEL.' + file_suffix)
                scalar_for_zonal_or_meridional_orientation = -1
                ecco_hFac_faces = ecco_hFacW_faces
            if boundary == 'east':
                vel_var_file = os.path.join(L0_run_dir, 'dv', model_name,
                                            model_name + '_' + boundary + '_VVEL.' + file_suffix)
                scalar_for_zonal_or_meridional_orientation = 1
                ecco_hFac_faces = ecco_hFacS_faces
            if boundary == 'north':
                vel_var_file = os.path.join(L0_run_dir, 'dv', model_name,
                                            model_name + '_' + boundary + '_VVEL.' + file_suffix)
                scalar_for_zonal_or_meridional_orientation = -1
                ecco_hFac_faces = ecco_hFacS_faces
        else:
            # this page has the tile map: https://ecco-v4-python-tutorial.readthedocs.io/ECCO_v4_Loading_the_ECCOv4_native_model_grid_parameters.html
            raise ValueError('Need to specify the velocity directions for '+model_name+' explicitly')
        vel_var_grid = np.fromfile(vel_var_file, dtype='>f4')
        vel_var_grid = vel_var_grid * scalar_for_zonal_or_meridional_orientation

        if var_name!='VEL':
            var_file = os.path.join(L0_run_dir, 'dv', model_name, model_name+'_'+ boundary + '_' + var_name +'.' + file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
        else:
            # the variable (concentration) is treated as a scalar - for velocity just use 1 as the scalar
            var_grid = np.copy(vel_var_grid)
            var_grid[vel_var_grid!=0]=1

        if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
            timesteps_in_file = int(np.size(var_grid) / (N))
            var_grid = np.reshape(var_grid, (timesteps_in_file, N))
            var_grid = var_grid[start_file_index:end_file_index+1, :]
            vel_var_grid = np.reshape(vel_var_grid, (timesteps_in_file, boundary_Nr, N))
            vel_var_grid = vel_var_grid[start_file_index:end_file_index + 1, :, :]
            for n in range(N):
                if faces[n]!=0:
                    points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 2] = ecco_DXG_faces[faces[n]][rows[n], cols[n]]
                    points[n, 3] = ecco_DYG_faces[faces[n]][rows[n], cols[n]]
                    hfac_points[:, n] = ecco_hFac_faces[faces[n]][:, rows[n], cols[n]]
                    values[index_counter:index_counter + (end_file_index - start_file_index) + 1, 0, n] = var_grid[:, n]
                    vel_values[index_counter:index_counter + (end_file_index - start_file_index) + 1, 0, n] = vel_var_grid[:, 0, n]
        else:
            timesteps_in_file = int(np.size(var_grid) / (boundary_Nr * N))
            var_grid = np.reshape(var_grid, (timesteps_in_file, boundary_Nr, N))
            var_grid = var_grid[start_file_index:end_file_index+1, :, :]
            vel_var_grid = np.reshape(vel_var_grid, (timesteps_in_file, boundary_Nr, N))
            vel_var_grid = vel_var_grid[start_file_index:end_file_index + 1, :, :]
                # var_grid = var_grid[start_file_index:end_file_index, :, :]
            for n in range(N):
                if faces[n] != 0:
                    points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 2] = ecco_DXG_faces[faces[n]][rows[n], cols[n]]
                    points[n, 3] = ecco_DYG_faces[faces[n]][rows[n], cols[n]]
                    hfac_points[:, n] = ecco_hFac_faces[faces[n]][:, rows[n], cols[n]]
                    values[index_counter:index_counter + (end_file_index - start_file_index)+1, :boundary_Nr, n] = var_grid[:, :, n]
                    vel_values[index_counter:index_counter + (end_file_index - start_file_index) + 1, :boundary_Nr,n] = vel_var_grid[:, :, n]

        index_counter += (end_file_index - start_file_index)+1
        # index_counter += (end_file_index - start_file_index)

    return(points,values,vel_values,hfac_points)

def subset_points_to_flux_boundary(boundary_XC_YC,points,values,vel_values,hfac):
    points_subset = np.zeros((np.shape(boundary_XC_YC)[0],4))
    values_subset = np.zeros((np.shape(values)[0],np.shape(values)[1],np.shape(boundary_XC_YC)[0]))
    vel_values_subset = np.zeros((np.shape(vel_values)[0], np.shape(vel_values)[1], np.shape(boundary_XC_YC)[0]))
    hfac_subset = np.zeros((np.shape(hfac)[0], np.shape(boundary_XC_YC)[0]))

    for r in range(np.shape(boundary_XC_YC)[0]):
        dist = (points[:,0]-boundary_XC_YC[r,0])**2 + (points[:,1]-boundary_XC_YC[r,1])**2
        index = np.argmin(dist)
        if dist[index]>0.000001:
            raise ValueError('Point was not found in dv mask - see subset_points_to_flux_boundary')
        points_subset[r,:] = points[index,:]
        values_subset[:,:,r] = values[:,:,index]
        vel_values_subset[:, :, r] = vel_values[:, :, index]
        hfac_subset[:,r] = hfac[:,index]

    return(points_subset,values_subset,vel_values_subset,hfac_subset)

def calculate_boundary_flux(boundary,points,values,vel_values,hfac, delR, boundary_Nr):

    timeseries = np.zeros((np.shape(values)[0],))
    if boundary in ['west','east']:
        width = points[:,2]
    if boundary in ['south','north']:
        width = points[:,3]

    flux_grid = np.zeros((np.shape(values)[0],len(delR),len(width)))

    for t in range(len(timeseries)):
        for i in range(np.shape(values)[2]):
            for d in range(boundary_Nr):
                #     [units/s]     []      * [m]    * [m]   * [units/m3]  * [m/s]
                flux_grid[t,d,i] = hfac[d,i]*width[i]*delR[d]*values[t,d,i]*vel_values[t,d,i]

    if boundary in ['north','east']:
        flux_grid *= -1

    return(flux_grid)

def write_timeseries_to_nc(output_file,var_name,timeseries):
    ds = nc4.Dataset(output_file,'w')
    ds.createDimension('time',len(timeseries))
    var = ds.createVariable(var_name,'f4',('time',))
    var[:] = timeseries
    ds.close()

def output_L0_fluxes_to_nc(output_file, delR,
                           boundaries, all_boundary_points, all_boundary_hfacs, all_mean_flux_grids, all_integrated_flux_timeseries):

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('depth_cells',len(delR))
    ds.createDimension('time',np.size(all_integrated_flux_timeseries[0]))

    drf_var = ds.createVariable('drF', 'f4', ('depth_cells', ))
    drf_var[:] = delR

    # hFacW_var = ds.createVariable('hFacW', 'f4', ('depth_cells','rows','cols_p1'))
    # hFacW_var[:, :, :] = hFacW
    #
    # hFacS_var = ds.createVariable('hFacS', 'f4', ('depth_cells','rows_p1','cols'))
    # hFacS_var[:, :, :] = hFacS

    for b in range(len(boundaries)):
        grp = ds.createGroup(boundaries[b])

        mean_flux_grid = all_mean_flux_grids[b]

        grp.createDimension('depth',np.shape(mean_flux_grid)[0])
        grp.createDimension('n_points', np.shape(mean_flux_grid)[1])

        XC_var = grp.createVariable('XC', 'f4', ('n_points',))
        XC_var[:] = all_boundary_points[b][:,0]

        YC_var = grp.createVariable('YC', 'f4', ('n_points',))
        YC_var[:] = all_boundary_points[b][:, 1]

        DXC_var = grp.createVariable('DXG', 'f4', ('n_points',))
        DXC_var[:] = all_boundary_points[b][:, 2]

        DXC_var = grp.createVariable('DYG', 'f4', ('n_points',))
        DXC_var[:] = all_boundary_points[b][:, 3]

        hfac_var = grp.createVariable('hfac', 'f4', ('depth', 'n_points'))
        hfac_var[:, :] = all_boundary_hfacs[b]

        mean_flux_var = grp.createVariable('mean_flux','f4',('depth','n_points'))
        mean_flux_var[:, :] = mean_flux_grid

        tvar = grp.createVariable('integrated_timeseries','f4',('time',))
        tvar[:] = all_integrated_flux_timeseries[b]

    ds.close()


########################################################################################################################


def create_timeseries(config_dir, ecco_dir, L1_model_name, boundaries,
                      boundary_var_name, start_year, final_year, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    # this is the dir where the obcs output will be stored
    if 'obcs' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'input','obcs'))
    if 'timeseries' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'input','obcs')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','timeseries'))
    if 'L0_timeseries' not in os.listdir(os.path.join(config_dir, 'L1', L1_model_name, 'input','obcs','timeseries')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs','timeseries','L0_timeseries'))
    output_dir = os.path.join(config_dir,'L1',L1_model_name,'input','obcs','L0_timeseries')

    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    llc = 270
    # ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
    #     ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    ecco_XC_faces, ecco_YC_faces, ecco_DXG_faces, ecco_DYG_faces, ecco_hFacS_faces, ecco_hFacW_faces = \
        read_flux_ecco_geometry_to_faces(ef, ecco_dir, llc)

    Nr_in = 50
    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    all_mean_flux_grids = []
    all_integrated_flux_timeseries = []
    all_boundary_points = []
    all_boundary_hfacs = []

    for mask_name in boundaries:

        dest_files, source_file_read_dict, source_file_read_index_sets = \
            create_src_dest_dicts_from_ref(config_dir, L1_model_name, mask_name, boundary_var_name, start_year, final_year)

        boundary_Nr = get_boundary_Nr_from_manual_list(L1_model_name,mask_name)

        print('    - Reading the mask to reference the variable to the llc grid')
        nc_dict_file = os.path.join(config_dir,'L0','input','L0_dv_mask_reference_dict.nc')
        source_faces, source_rows, source_cols = read_mask_reference_from_nc_dict(nc_dict_file, L1_model_name, mask_name)

        L0_boundary_points_manually_selected = read_flux_points_from_shapefile(config_dir,L1_model_name,mask_name)

        timestep_counter = 0
        mean_flux_grid = np.zeros((len(delR), len(L0_boundary_points_manually_selected)))

        total_timesteps = 0
        for year in range(start_year,final_year+1):
            if year % 4 == 0:
                total_timesteps += 366
            else:
                total_timesteps += 365
        integrated_flux_timeseries = np.zeros((total_timesteps,))

        for dest_file in dest_files:

            print('    - Integrating the fluxes for ' + str(dest_file))
            source_files = source_file_read_dict[dest_file]

            if len(source_files)>0:

                source_file_read_indices = source_file_read_index_sets[dest_file]

                if print_level >= 4:
                    print('                - Reading in the L0 diagnostics_vec output')
                L0_run_dir = os.path.join(config_dir,'L0','run')
                L0_boundary_points, L0_boundary_values, L0_boundary_vel_values, L0_boundary_points_hFac = \
                    read_L0_boundary_variable_points(L0_run_dir, L1_model_name, mask_name, boundary_var_name,
                                                     source_files, source_file_read_indices,
                                                     llc, boundary_Nr, Nr_in, source_faces, source_rows, source_cols,
                                                     ecco_XC_faces, ecco_YC_faces, ecco_DXG_faces, ecco_DYG_faces,
                                                     ecco_hFacS_faces, ecco_hFacW_faces,
                                                     print_level)

                L0_boundary_points_subset, L0_boundary_values_subset, L0_boundary_vel_values_subset, L0_boundary_points_hFac_subset = \
                    subset_points_to_flux_boundary(L0_boundary_points_manually_selected,
                                                   L0_boundary_points, L0_boundary_values, L0_boundary_vel_values, L0_boundary_points_hFac)

                flux_grid = calculate_boundary_flux(mask_name,L0_boundary_points_subset, L0_boundary_values_subset,
                                                          L0_boundary_vel_values_subset, L0_boundary_points_hFac_subset, delR, boundary_Nr)

                mean_flux_grid += np.sum(flux_grid, axis=0)
                integrated_flux_timeseries[timestep_counter:timestep_counter + np.shape(flux_grid)[0]] = np.sum(np.sum(flux_grid, axis=1), axis=1)
                timestep_counter += np.shape(flux_grid)[0]

                # output_file = os.path.join(output_dir,mask_name,boundary_var_name,dest_file)
                # write_timeseries_to_nc(output_file, boundary_var_name, flux_timeseries)

                # plt.plot(L0_boundary_points[:,0],L0_boundary_points[:,1],'k.')
                # plt.plot(L0_boundary_points_subset[:, 0], L0_boundary_points_subset[:, 1], 'g.')
                # plt.show()
        all_boundary_points.append(L0_boundary_points_subset)
        all_boundary_hfacs.append(L0_boundary_points_hFac_subset)

        all_mean_flux_grids.append(mean_flux_grid / timestep_counter)
        all_integrated_flux_timeseries.append(integrated_flux_timeseries)

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs', 'timeseries', 'L0_timeseries',
                               'L0_'+L1_model_name[3:] + '_' + boundary_var_name + '_boundary_flux.nc')
    output_L0_fluxes_to_nc(output_file, delR, # Depth, drF,
                           boundaries, all_boundary_points, all_boundary_hfacs,
                           all_mean_flux_grids, all_integrated_flux_timeseries)

def create_L0_timeseries(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
                         start_year, final_year, print_level):

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

    print('Creating the L0 timeseries for ' + var_name + ' to cover years ' + str(start_year)+' to ' + str(final_year))

    print('  Running the timeseries routine:')
    create_timeseries(config_dir, ecco_dir,  L1_model_name, boundaries, var_name, start_year, final_year, print_level)


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
    L1_model_name = args.L1_model_name
    boundaries = args.boundaries
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month

    create_L0_timeseries(config_dir, ecco_dir, L1_model_name, boundaries, proc_id,
                                       start_year, final_year, start_month,
                                       final_month, print_level=4)

