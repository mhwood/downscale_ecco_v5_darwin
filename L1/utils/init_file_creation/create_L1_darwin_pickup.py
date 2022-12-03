
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_geometry(config_dir,model_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    delR = ds.variables['drF'][:]
    ds.close()

    return(XC, YC, AngleCS, AngleSN, delR)

def read_grid_mask(config_dir,model_name,var_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)

    if var_name in ['Vvel','GvNm1', 'GvNm2']:
        hFac = 'S'
    elif var_name in ['Uvel','GuNm1', 'GuNm2']:
        hFac = 'W'
    else:
        hFac = 'C'

    HFac = ds.variables['HFac' + hFac][:, :, :]
    if hFac == 'W':
        HFac = HFac[:, :, :-1]
    if hFac == 'S':
        HFac = HFac[:, :-1, :]

    ds.close()

    mask = np.copy(HFac)
    mask[mask>0]=1

    return(mask)

def read_interpolation_grid(df,config_dir, model_name, var_name,
                            ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain, XC, YC, print_level):

    interpolation_grid_file = model_name+'_interpolation_grid.nc'

    if interpolation_grid_file not in os.listdir(os.path.join(config_dir,'nc_grids')):

        if print_level >= 3:
            print('            - Creating the interpolation grids for this model (hasnt been done yet!)')

        domain_wet_cells_3D_C = read_grid_mask(config_dir,model_name,'Theta')
        domain_wet_cells_3D_S = read_grid_mask(config_dir, model_name, 'Vvel')
        domain_wet_cells_3D_W = read_grid_mask(config_dir, model_name, 'Uvel')

        interpolation_type_grid_C, source_row_grid_C, source_col_grid_C, source_level_grid_C, \
        interpolation_type_grid_S, source_row_grid_S, source_col_grid_S, source_level_grid_S, \
        interpolation_type_grid_W, source_row_grid_W, source_col_grid_W, source_level_grid_W = \
            df.create_interpolation_grids(ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain,
                               XC, YC, domain_wet_cells_3D_C, domain_wet_cells_3D_S, domain_wet_cells_3D_W)

        df.write_interpolation_grid_to_nc(config_dir, model_name,
                                       interpolation_type_grid_C, source_row_grid_C, source_col_grid_C, source_level_grid_C,
                                       interpolation_type_grid_S, source_row_grid_S, source_col_grid_S, source_level_grid_S,
                                       interpolation_type_grid_W, source_row_grid_W, source_col_grid_W, source_level_grid_W)

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

def interpolate_ecco_wetgrid_to_domain(XC, YC, delR, ecco_XC, ecco_YC, ecco_delR, ecco_wet_cells):

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    ecco_Z_bottom = np.cumsum(ecco_delR)
    ecco_Z_top = np.concatenate([np.array([0]), ecco_Z_bottom[:-1]])
    ecco_Z = (ecco_Z_bottom + ecco_Z_top) / 2

    ##########
    # interpolate vertically first

    ecco_wet_cells_interpolated = np.zeros((len(delR),np.shape(ecco_wet_cells)[1], np.shape(ecco_wet_cells)[2]))
    for i in range(len(delR)):
        if np.any(ecco_Z>Z[i]):
            ecco_Z_subset = ecco_Z[ecco_Z>Z[i]]
            index = np.argmin(np.abs(ecco_Z-np.min(ecco_Z_subset)))
        else:
            index = len(ecco_delR)-1
        ecco_wet_cells_interpolated[i,:,:] = ecco_wet_cells[index,:,:]

    ##########
    # then interpolate onto the tile

    ecco_wet_cells_on_tile_domain = np.zeros((len(delR),np.shape(XC)[0],np.shape(XC)[1]))

    points = np.hstack([np.reshape(ecco_XC, (np.size(ecco_XC), 1)),
                        np.reshape(ecco_YC, (np.size(ecco_YC), 1))])

    for i in range(len(delR)):
        values = np.reshape(ecco_wet_cells_interpolated[i,:,:], (np.size(ecco_wet_cells_interpolated[i,:,:]), 1))
        grid = griddata(points, values, (XC, YC), method='nearest')[:, :, 0]
        ecco_wet_cells_on_tile_domain[i,:,:] = grid

    return(ecco_wet_cells_on_tile_domain)

def stack_grids_to_pickup(interp_grids, var_names):

    for i in range(len(interp_grids)):
        print('       - Stacking field for '+var_names[i])
        var_grid = interp_grids[i]
        if i==0:
            pickup_grid = var_grid
        else:
            pickup_grid = np.concatenate([pickup_grid, var_grid], axis=0)

    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata, Nr):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+",\n"
    output += " " + "{:5d}".format(Nr) + ",    1," + "{:5d}".format(Nr) + "\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [          "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [    "+str(subset_metadata['nflds'][0])+" ];\n"
    output += " fldList = {\n "
    for var_name in subset_metadata['fldlist']:
        output += "'"+var_name
        for i in range(8-len(var_name)):
            output+= " "
        output+="' "
    output += "\n };"

    f = open(output_file+'.meta','w')
    f.write(output)
    f.close()

def create_L1_darwin_pickup_file(config_dir, model_name,
                               ecco_dir, llc,ordered_ecco_tiles, ordered_ecco_tile_rotations,
                               parent_model_pickup_iteration, print_level, use_interpolation_grids = True):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if print_level>=1:
        print('    - Creating the pickup file for the '+model_name+' model from ECCO data')

    if print_level>=1:
        print('    - Reading in the L1 tile geometry')
    # step 0: get the model domain
    XC, YC, _, _, delR = read_grid_geometry(config_dir, model_name)
    Nr = len(delR)

    if print_level>=1:
        print('    - Reading in the ECCO tile geometry')
    # step 1: get the ecco faces geometry
    # ecco_AngleCS, ecco_AngleSN, ecco_hfacC,
    ecco_XC, ecco_YC, _, _, ecco_hfacC, ecco_delR = \
        ef.read_ecco_grid_geometry(ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    # C = plt.imshow(ecco_AngleSN, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    if print_level>=1:
        print('    - Reading in the ECCO pickup file')
    pickup_file = 'pickup_darwin.'+'{:010d}'.format(parent_model_pickup_iteration)
    pickup_file_path = os.path.join(config_dir,'L0','run_darwin',pickup_file)
    var_names, var_grids, global_metadata = ef.read_ecco_darwin_pickup_to_stiched_grid(pickup_file_path, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    for vn in range(len(var_names)):
        print(var_names[vn],np.shape(var_grids[vn]))

    # C = plt.imshow(var_grids[0][0,:,:], origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.colorbar(C)
    # plt.show()

    ecco_wet_cells = np.copy(ecco_hfacC)
    ecco_wet_cells[ecco_wet_cells>0]=1

    if use_interpolation_grids:
        L1_interpolation_mask, L1_source_rows, L1_source_cols, L1_source_levels = \
            read_interpolation_grid(df, config_dir, model_name, var_name,
                                    ecco_XC, ecco_YC, ecco_wet_cells, ecco_wet_cells_on_tile_domain,
                                    XC, YC, print_level)

    # make some bins where the tiles will be stored
    interp_grids = []

    ####################################################################################################################

    if print_level >= 3:
        print('            - Interpolating the ECCO wetgrid onto the domain')
    ecco_wet_cells_on_tile_domain = interpolate_ecco_wetgrid_to_domain(XC, YC, delR,
                                                                       ecco_XC, ecco_YC, ecco_delR, ecco_wet_cells)

    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            ecco_grid = var_grids[var_names.index(var_name)]
            if print_level >= 3:
                print('            - Downscaling ' + var_name)

            domain_wet_cells_3D = read_grid_mask(config_dir,model_name,var_name)

            # plt.subplot(1, 3, 1)
            # plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
            # plt.subplot(1, 3, 2)
            # plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
            # plt.subplot(1, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.show()

            mean_vertical_difference = 0
            subset_copy = np.copy(ecco_grid)

            # if var_name.lower() not in ['etan', 'detahdt', 'etah']:
            #     ecco_grid, ecco_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
            #         ecco_grid, ecco_wet_cells, ecco_delR, delR)
            # print('     - Skipping the vertical interpolation')

            # # plt.subplot(1,2,1)
            # # plt.imshow(subset_copy[:,10,:])
            # # plt.subplot(1, 2, 2)
            # # plt.imshow(ecco_grid[:, 10, :])
            # # plt.show()

            domain_wet_cells_3D_for_interpolation = domain_wet_cells_3D

            if print_level >= 4:
                printing = True
            else:
                printing = False

            if use_interpolation_grids:
                interp_field = df.downscale_3D_field_with_interpolation_mask(ecco_XC, ecco_YC,
                                                                             ecco_grid, ecco_wet_cells,
                                                                             ecco_wet_cells_on_tile_domain,
                                                                             XC, YC, domain_wet_cells_3D_for_interpolation,
                                                                             L1_interpolation_mask, L1_source_rows,
                                                                             L1_source_cols, L1_source_levels,
                                                                             mean_vertical_difference=0,
                                                                             fill_downward=True,
                                                                             remove_zeros=True, printing=printing,
                                                                             testing=True)
            else:
                interp_field = df.downscale_3D_field(ecco_XC, ecco_YC,
                                                     ecco_grid, ecco_wet_cells,
                                                     ecco_wet_cells_on_tile_domain,
                                                     XC, YC, domain_wet_cells_3D_for_interpolation,
                                                     mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                     printing=printing, testing=True)

            if np.sum(np.isnan(interp_field))>0:
                raise ValueError('Found nans in the pickup grid...')
                # print('Setting '+str(np.sum(np.isnan(interp_field)))+' values to 0 in this grid')
                # interp_field[np.isnan(interp_field)] = 0

            # # if var_name=='Theta':
            # plt.subplot(2, 3, 1)
            # plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 4)
            # plt.imshow(ecco_grid[0, :, :], origin='lower')
            # plt.subplot(2, 3, 5)
            # plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.title(var_name)
            # plt.show()

            interp_grids.append(interp_field)

    if print_level >= 1:
        print('    - Stacking the interpolated fields into a pickup grid')
    pickup_grid = stack_grids_to_pickup(interp_grids,var_names)

    print(np.shape(pickup_grid))

    if print_level >= 1:
        print('    - Outputting the compact pickup grid to the input directory')
    pickup_metadata = dict(global_metadata)
    output_dir = os.path.join(config_dir, 'L1', model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_darwin.' + '{:010d}'.format(4*parent_model_pickup_iteration))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [4*int(pickup_metadata['timestepnumber'][0])]

    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata, Nr)
