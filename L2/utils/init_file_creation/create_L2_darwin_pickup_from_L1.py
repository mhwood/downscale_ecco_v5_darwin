
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.interpolate import griddata
from MITgcmutils import mds
import netCDF4 as nc4
import argparse
import ast
import sys


def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return(XC,YC,AngleCS,AngleSN,delR)

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

def interpolate_L1_wetgrid_to_L2_domain(L2_XC, L2_YC, delR_out, L1_XC, L1_YC, delR_in, L1_wet_cells):

    Z_bottom = np.cumsum(delR_out)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    L1_Z_bottom = np.cumsum(delR_in)
    L1_Z_top = np.concatenate([np.array([0]), L1_Z_bottom[:-1]])
    L1_Z = (L1_Z_bottom + L1_Z_top) / 2

    ##########
    # interpolate vertically first

    L1_wet_cells_interpolated = np.zeros((len(delR_out),np.shape(L1_wet_cells)[1], np.shape(L1_wet_cells)[2]))
    for i in range(len(delR_out)):
        if np.any(L1_Z>Z[i]):
            L1_Z_subset = L1_Z[L1_Z>Z[i]]
            index = np.argmin(np.abs(L1_Z-np.min(L1_Z_subset)))
        else:
            index = len(delR_out)-1
        L1_wet_cells_interpolated[i,:,:] = L1_wet_cells[index,:,:]

    ##########
    # then interpolate onto the tile

    L1_wet_cells_on_tile_domain = np.zeros((len(delR_out),np.shape(L2_XC)[0],np.shape(L2_XC)[1]))

    points = np.hstack([np.reshape(L1_XC, (np.size(L1_XC), 1)),
                        np.reshape(L1_YC, (np.size(L1_YC), 1))])

    for i in range(len(delR_out)):
        values = np.reshape(L1_wet_cells_interpolated[i,:,:], (np.size(L1_wet_cells_interpolated[i,:,:]), 1))
        grid = griddata(points, values, (L2_XC, L2_YC), method='nearest')[:, :, 0]
        L1_wet_cells_on_tile_domain[i,:,:] = grid

    return(L1_wet_cells_on_tile_domain)

def read_interpolation_grid(df,config_dir, model_name, var_name,
                            L1_XC, L1_YC, L1_wet_cells, L1_wet_cells_on_L2_domain, L2_XC, L2_YC, print_level):

    interpolation_grid_file = model_name+'_interpolation_grid.nc'

    if interpolation_grid_file not in os.listdir(os.path.join(config_dir,'nc_grids')):

        if print_level >= 3:
            print('            - Creating the interpolation grids for this model (hasnt been done yet!)')

        L2_domain_wet_cells_3D_C = read_grid_mask(config_dir,model_name,'Theta')
        L2_domain_wet_cells_3D_S = read_grid_mask(config_dir, model_name, 'Vvel')
        L2_domain_wet_cells_3D_W = read_grid_mask(config_dir, model_name, 'Uvel')

        interpolation_type_grid_C, source_row_grid_C, source_col_grid_C, source_level_grid_C, \
        interpolation_type_grid_S, source_row_grid_S, source_col_grid_S, source_level_grid_S, \
        interpolation_type_grid_W, source_row_grid_W, source_col_grid_W, source_level_grid_W = \
            df.create_interpolation_grids(L1_XC, L1_YC, L1_wet_cells, L1_wet_cells_on_L2_domain,
                               L2_XC, L2_YC, L2_domain_wet_cells_3D_C, L2_domain_wet_cells_3D_S, L2_domain_wet_cells_3D_W)

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

def read_pickup_file_to_compact(pickup_file_path, Nr):

    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    has_Nr = True

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def read_pickup_grid(pickup_file_path, Nr):

    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path, Nr)

    # for vn in range(len(var_names)):
    #     print(var_names[vn],np.shape(compact_var_grids[vn]))

    return(var_names, compact_var_grids, global_metadata)

def read_grid_mask_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    if var_name in ['Uvel', 'GuNm1', 'GuNm2']:
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name in ['Vvel', 'GvNm1', 'GvNm2']:
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(mask)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
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
    output += " nrecords = [         "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [   "+str(subset_metadata['nflds'][0])+" ];\n"
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

########################################################################################################################


def create_L2_darwin_pickup_file(config_dir, L1_model_name, L1_iteration, L2_model_name, print_level):

    if print_level>=1:
        print('    - Creating the darwin pickup file for the ' + L2_model_name + ' model from the '+L1_model_name+' model')

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, delR_in = read_grid_geometry_from_nc(config_dir, L1_model_name)
    Nr_in = len(delR_in)

    L1_Hfac = read_grid_mask_from_nc(config_dir, L1_model_name, 'Theta')
    L1_wet_cells = np.copy(L1_Hfac)
    L1_wet_cells[L1_wet_cells > 0] = 1

    L2_XC, L2_YC, L2_AngleCS, L2_AngleSN, delR_out = read_grid_geometry_from_nc(config_dir,L2_model_name)
    Nr_out = len(delR_out)

    if print_level >= 3:
        print('            - Interpolating the L1 wetgrid onto the L2 domain')
    L1_wet_cells_on_L2_domain = interpolate_L1_wetgrid_to_L2_domain(L2_XC, L2_YC, delR_out,
                                                                    L1_XC, L1_YC, delR_in, L1_wet_cells)

    pickup_file = 'pickup_darwin.' + '{:010d}'.format(L1_iteration)
    pickup_file_path = os.path.join(config_dir, 'L1', L1_model_name, 'run', pickup_file)
    var_names, var_grids, pickup_metadata = read_pickup_grid(pickup_file_path, Nr_in)
    print(var_names)

    if print_level >= 1:
        print('    - Downscaling the darwin pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            var_grid = var_grids[0]
            print(np.shape(var_grid))
            if print_level >= 2:
                print('        - Downscaling ' + var_name)

            L2_domain_wet_cells_3D = read_grid_mask_from_nc(config_dir,L2_model_name,var_name)

            # mean_vertical_difference = 0
            # subset_copy = np.copy(var_grid)

            L1_wet_cells = np.copy(L1_Hfac)
            L1_wet_cells[L1_wet_cells > 0] = 1

            if print_level >= 3:
                print('            - var_grid shape: ' + str(np.shape(var_grid)))
                print('            - L1_wet_cells shape: ' + str(np.shape(L1_wet_cells)))

            if Nr_in!=Nr_out:
                if print_level >= 3:
                    print('            - Interpolating to new depth levels ')
                var_grid, L1_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
                    var_grid, L1_wet_cells, delR_in, delR_out)

            L2_interpolation_mask, L2_source_rows, L2_source_cols, L2_source_levels = \
                read_interpolation_grid(df, config_dir, L2_model_name, var_name,
                                        L1_XC, L1_YC, L1_wet_cells, L1_wet_cells_on_L2_domain,
                                        L2_XC, L2_YC, print_level)

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            L1_wet_cells_on_L2_domain_3D = np.copy(L2_domain_wet_cells_3D)

            if print_level >= 3:
                print('            - L1_XC shape: '+str(np.shape(L1_XC)))
                print('            - L1_YC shape: ' + str(np.shape(L1_YC)))
                print('            - var_grid shape: ' + str(np.shape(var_grid)))
                print('            - L1_wet_cells shape: ' + str(np.shape(L1_wet_cells)))
                print('            - L2_XC shape: ' + str(np.shape(L2_XC)))
                print('            - L2_YC shape: ' + str(np.shape(L2_YC)))
                print('            - L2 wet cells shape: ' + str(np.shape(L2_domain_wet_cells_3D)))

            if print_level >= 4:
                printing = True
            else:
                printing = False

            # interp_field = df.downscale_3D_field(L1_XC, L1_YC,
            #                                      var_grid, L1_wet_cells,
            #                                      L1_wet_cells_on_domain_3D,
            #                                      L2_XC, L2_YC, L2_domain_wet_cells_3D,
            #                                      mean_vertical_difference=0, fill_downward=True, remove_zeros=False,
            #                                      printing=True, testing=True)

            interp_field = df.downscale_3D_field_with_interpolation_mask(L1_XC, L1_YC,
                                                                         var_grid, L1_wet_cells,
                                                                         L1_wet_cells_on_L2_domain,
                                                                         L2_XC, L2_YC, L2_domain_wet_cells_3D,
                                                                         L2_interpolation_mask, L2_source_rows,
                                                                         L2_source_cols, L2_source_levels,
                                                                         mean_vertical_difference=0, fill_downward=True,
                                                                         remove_zeros=True, printing=printing,
                                                                         testing=False)

            if print_level >= 3:
                print('            - Field output shape: '+str(np.shape(interp_field)))

            if np.sum(np.isnan(interp_field)) > 0:
                print('Setting ' + str(np.sum(np.isnan(interp_field))) + ' values to 0 in this grid')
                interp_field[np.isnan(interp_field)] = 0

            # plt.subplot(2, 3, 1)
            # plt.imshow(L1_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(L1_wet_cells_on_domain_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(L2_domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 2, 3)
            # C = plt.imshow(var_grid[0, :, :], origin='lower')
            # plt.colorbar(C)
            # plt.subplot(2, 2, 4)
            # C = plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.colorbar(C)
            # plt.show()

            # interp_field[L2_domain_wet_cells_3D[:np.shape(interp_field)[0], :, :] == 0] = 0
        else:
            interp_field = np.zeros((len(delR_out), np.shape(L2_XC)[0], np.shape(L2_XC)[1]))

        # if var_name.lower() in ['etan','etah']:
        #     interp_field[interp_field!=0]+=2

        interp_grids.append(interp_field)

    # plt.subplot(2,2,1)
    # plt.imshow(old_grids[0][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.subplot(2, 2, 2)
    # plt.imshow(old_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 3)
    # plt.imshow(interp_grids[0][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 4)
    # plt.imshow(interp_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.show()

    if print_level >= 2:
        print('        - Stacking the grids for output')
    pickup_grid = stack_grids_to_pickup(interp_grids)
    if print_level >= 2:
        print('        - Output shape: '+str(np.shape(pickup_grid)))

    if print_level >= 1:
        print('    - Outputting the compact pickup grid to the input directory')
    # pickup_metadata = dict(global_metadata)
    output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_darwin.' + '{:010d}'.format(5 * L1_iteration))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [5 * int(pickup_metadata['timestepnumber'][0])]

    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata, Nr_out)



