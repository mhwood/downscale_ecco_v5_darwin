
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
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

def read_pickup_file_to_compact(pickup_file_path, Nr):

    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
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

def rotate_grids_to_natural_grids(uvel, vvel, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        for k in range(np.shape(uvel)[0]):
            zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
            meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
        return (zonal_velocity, meridional_velocity)

    zonal_uvel, meridional_vvel = rotate_velocity_vectors_to_natural(AngleCS, AngleSN,
                                                                     uvel, vvel)

    # plt.subplot(2,2,1)
    # C = plt.imshow(var_grids[uvel_grid_index][0,:,:],origin='lower',cmap='seismic',vmin=-0.4,vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(var_grids[vvel_grid_index][0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(zonal_uvel[0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(meridional_vvel[0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.show()

    return(zonal_uvel, meridional_vvel)

def rotate_directional_fields_to_natural(var_grids,var_names,L1_AngleCS, L1_AngleSN):
    uvel_grid = var_grids[var_names.index('Uvel')]
    vvel_grid = var_grids[var_names.index('Vvel')]
    natural_uvel_grid, natural_vvel_grid = rotate_grids_to_natural_grids(uvel_grid, vvel_grid, L1_AngleCS, L1_AngleSN)
    var_grids[var_names.index('Uvel')] = natural_uvel_grid
    var_grids[var_names.index('Vvel')] = natural_vvel_grid

    gunm1_grid = var_grids[var_names.index('GuNm1')]
    gvnm1_grid = var_grids[var_names.index('GvNm1')]
    natural_gunm1_grid, natural_gvnm1_grid = rotate_grids_to_natural_grids(gunm1_grid, gvnm1_grid, L1_AngleCS,
                                                                           L1_AngleSN)
    var_grids[var_names.index('GuNm1')] = natural_gunm1_grid
    var_grids[var_names.index('GvNm1')] = natural_gvnm1_grid

    gunm2_grid = var_grids[var_names.index('GuNm2')]
    gvnm2_grid = var_grids[var_names.index('GvNm2')]
    natural_gunm2_grid, natural_gvnm2_grid = rotate_grids_to_natural_grids(gunm2_grid, gvnm2_grid, L1_AngleCS,
                                                                           L1_AngleSN)
    var_grids[var_names.index('GuNm2')] = natural_gunm2_grid
    var_grids[var_names.index('GvNm2')] = natural_gvnm2_grid
    return(var_grids,var_names)


def rotate_interpolated_grids_to_domain(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    uvel_grid_index = var_names.index('Uvel')
    vvel_grid_index = var_names.index('Vvel')
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   var_grids[uvel_grid_index], var_grids[vvel_grid_index])

    gunm1_grid_index = var_names.index('GuNm1')
    gvnm1_grid_index = var_names.index('GvNm1')
    gunm1, gvnm1 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm1_grid_index], var_grids[gvnm1_grid_index])

    gunm2_grid_index = var_names.index('GuNm2')
    gvnm2_grid_index = var_names.index('GvNm2')
    gunm2, gvnm2 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm2_grid_index], var_grids[gvnm2_grid_index])

    var_grids[uvel_grid_index] = uvel
    var_grids[vvel_grid_index] = vvel

    var_grids[gunm1_grid_index] = gunm1
    var_grids[gvnm1_grid_index] = gvnm1

    var_grids[gunm2_grid_index] = gunm2
    var_grids[gvnm2_grid_index] = gvnm2

    return(var_grids)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
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


def create_L2_pickup_file(config_dir, L1_model_name, L1_iteration, L2_model_name, zero_etan, print_level):

    if print_level>=1:
        print('    - Creating the pickup file for the ' + L2_model_name + ' model from the '+L1_model_name+' model')

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, delR_in = read_grid_geometry_from_nc(config_dir, L1_model_name)
    Nr_in = len(delR_in)

    L1_Hfac = read_grid_mask_from_nc(config_dir, L1_model_name, 'Theta')
    L1_wet_cells = np.copy(L1_Hfac)
    L1_wet_cells[L1_wet_cells > 0] = 1

    L2_XC, L2_YC, L2_AngleCS, L2_AngleSN, delR_out = read_grid_geometry_from_nc(config_dir,L2_model_name)
    Nr_out = len(delR_out)

    pickup_file = 'pickup.'+'{:010d}'.format(L1_iteration)
    pickup_file_path = os.path.join(config_dir,'L1',L1_model_name,'run',pickup_file)
    var_names, var_grids, pickup_metadata = read_pickup_grid(pickup_file_path, Nr_in)

    var_grids, var_names = rotate_directional_fields_to_natural(var_grids, var_names, L1_AngleCS, L1_AngleSN)

    if print_level >= 1:
        print('    - Downscaling the pickup grids')
    interp_grids = []
    output_var_names = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            var_grid = var_grids[vn]
            if print_level >= 2:
                print('        - Downscaling ' + var_name)

            domain_wet_cells_3D = read_grid_mask_from_nc(config_dir,L2_model_name,var_name)

            # mean_vertical_difference = 0
            # subset_copy = np.copy(var_grid)

            L1_wet_cells = np.copy(L1_Hfac)
            L1_wet_cells[L1_wet_cells > 0] = 1

            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                if Nr_in!=Nr_out:
                    var_grid, L1_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
                        var_grid, L1_wet_cells, delR_in, delR_out)
            else:
                domain_wet_cells_3D = domain_wet_cells_3D[:1,:,:]

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            L1_wet_cells_on_domain_3D = np.copy(domain_wet_cells_3D)

            if print_level >= 3:
                print('            - L1_XC shape: '+str(np.shape(L1_XC)))
                print('            - L1_YC shape: ' + str(np.shape(L1_YC)))
                print('            - var_grid shape: ' + str(np.shape(var_grid)))
                print('            - L1_wet_cells shape: ' + str(np.shape(L1_wet_cells)))
                print('            - L2_XC shape: ' + str(np.shape(L2_XC)))
                print('            - L2_YC shape: ' + str(np.shape(L2_YC)))
                print('            - L1_XC shape: ' + str(np.shape(domain_wet_cells_3D)))

            interp_field = df.downscale_3D_field(L1_XC, L1_YC,
                                                 var_grid, L1_wet_cells,
                                                 L1_wet_cells_on_domain_3D,
                                                 L2_XC, L2_YC, domain_wet_cells_3D,
                                                 mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                 printing=True, testing=False)
            if print_level >= 3:
                print('            - Field output shape: '+str(np.shape(interp_field)))

            if np.sum(np.isnan(interp_field)) > 0:
                print('Setting ' + str(np.sum(np.isnan(interp_field))) + ' values to 0 in this grid')
                interp_field[np.isnan(interp_field)] = 0

            if zero_etan and var_name.lower() == 'etan':
                etan_field = np.copy(interp_field)
                indices = etan_field!=0
                etan_field[indices] -= np.mean(etan_field[indices])
                interp_field = etan_field
                # interp_field*=0

            if zero_etan and var_name.lower() == 'etah':
                etan_field = np.copy(interp_field)
                indices = etan_field != 0
                etan_field[indices] -= np.mean(etan_field[indices])
                interp_field = etan_field
                # interp_field*=0

            # plt.subplot(2, 3, 1)
            # plt.imshow(L1_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(L1_wet_cells_on_domain_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 2, 3)
            # C = plt.imshow(var_grid[0, :, :], origin='lower')
            # plt.colorbar(C)
            # plt.subplot(2, 2, 4)
            # C = plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.colorbar(C)
            # plt.show()

            # interp_field[domain_wet_cells_3D[:np.shape(interp_field)[0], :, :] == 0] = 0
        else:
            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                interp_field = np.zeros((len(delR_out), np.shape(L2_XC)[0], np.shape(L2_XC)[1]))
            else:
                interp_field = np.zeros((1, np.shape(L2_XC)[0], np.shape(L2_XC)[1]))

        # if var_name.lower() in ['etan','etah']:
        #     interp_field[interp_field!=0]+=2

        interp_grids.append(interp_field)

    # old_grids = [np.copy(interp_grids)[0],np.copy(interp_grids)[1]]

    interp_grids = rotate_interpolated_grids_to_domain(var_names, interp_grids, L2_AngleCS, L2_AngleSN)

    # plt.subplot(2,2,1)
    # plt.imshow(old_grids[0][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.subplot(2, 2, 2)
    # plt.imshow(old_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 3)
    # plt.imshow(interp_grids[0][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 4)
    # plt.imshow(interp_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.show()

    pickup_grid = stack_grids_to_pickup(interp_grids)

    output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup.' + '{:010d}'.format(L1_iteration*5))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [int(1)]
    pickup_metadata['nFlds'] = [11]
    pickup_metadata['nrecords'] = [int((len(var_names) - 3) * len(delR_out) + 3)]
    pickup_metadata['fldlist'] = var_names
    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)



