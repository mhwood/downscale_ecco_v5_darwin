
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

def read_L1_grid_tile_geometry(config_dir,model_name, sNx, sNy,
                               ordered_nonblank_tiles, ordered_nonblank_rotations,
                               faces, ordered_tiles_faces_dict):

    stitched_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict), sNx*len(ordered_nonblank_tiles[1])))
    stitched_YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))
    stitched_AngleCS_grid = np.zeros((sNy * len(ordered_tiles_faces_dict), sNx * len(ordered_nonblank_tiles[1])))
    stitched_AngleSN_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))
    stitched_HFac_grid = np.zeros((1, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    HFac = ds.variables['HFacC'][:,:,:]
                    HFac = HFac[:1, :, :]
                    ds.close()

                    for i in range(ordered_nonblank_rotations[r][c]):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)
                        HFac = np.rot90(HFac, axes=(1,2))

                    stitched_XC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = XC
                    stitched_YC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = YC
                    stitched_AngleCS_grid[r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = AngleCS
                    stitched_AngleSN_grid[r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = AngleSN
                    stitched_HFac_grid[:, r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = HFac

    return(stitched_XC_grid, stitched_YC_grid, stitched_AngleCS_grid, stitched_AngleSN_grid, stitched_HFac_grid)

def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    ds.close()
    return(XC,YC,AngleCS,AngleSN)

def read_seaice_pickup_file_to_compact(pickup_file_path):

    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def read_seaice_pickup(pickup_file_path):

    var_names,row_bounds,compact_var_grids,global_metadata = read_seaice_pickup_file_to_compact(pickup_file_path)

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
    hfac = hfac[:1,:,:]
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

    uvel_grid = var_grids[var_names.index('siUICE')]
    vvel_grid = var_grids[var_names.index('siVICE')]
    natural_uvel_grid, natural_vvel_grid = rotate_grids_to_natural_grids(uvel_grid, vvel_grid, L1_AngleCS, L1_AngleSN)
    var_grids[var_names.index('siUICE')] = natural_uvel_grid
    var_grids[var_names.index('siVICE')] = natural_vvel_grid

    return(var_grids,var_names)

def rotate_interpolated_grids_to_domain(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    uvel_grid_index = var_names.index('siUICE')
    vvel_grid_index = var_names.index('siVICE')
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   var_grids[uvel_grid_index], var_grids[vvel_grid_index])

    var_grids[uvel_grid_index] = uvel
    var_grids[vvel_grid_index] = vvel

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

def write_seaice_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [    "+str(subset_metadata['nrecords'][0])+" ];\n"
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

########################################################################################################################


def create_L2_seaice_pickup_file(config_dir, L1_model_name, L1_iteration, L2_model_name, print_level):

    print('    - Creating the pickup file for the ' + L2_model_name + ' model')

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the output will be stored
    # output_dir = os.path.join(config_dir,'L2',L2_model_name,'input')

    print('    - Reading in the geometry of the L1 domain')
    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN = read_grid_geometry_from_nc(config_dir, L1_model_name)

    L1_Hfac = read_grid_mask_from_nc(config_dir, L1_model_name, 'Theta')
    L1_wet_cells = np.copy(L1_Hfac)
    L1_wet_cells[L1_wet_cells > 0] = 1

    L2_XC, L2_YC, L2_AngleCS, L2_AngleSN = read_grid_geometry_from_nc(config_dir,L2_model_name)

    pickup_file = 'pickup_seaice.'+'{:010d}'.format(L1_iteration)
    pickup_file_path = os.path.join(config_dir,'L1',L1_model_name,'run',pickup_file)
    var_names, var_grids, pickup_metadata = read_seaice_pickup(pickup_file_path)

    # plt.imshow(var_grids[5][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.show()

    var_grids, var_names = rotate_directional_fields_to_natural(var_grids, var_names, L1_AngleCS, L1_AngleSN)

    # plt.imshow(var_grids[5][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.show()

    print('    - Downscaling the pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            var_grid = var_grids[vn]
            print('      - Downscaling ' + var_name)
            print(np.shape(var_grid))

            domain_wet_cells = read_grid_mask_from_nc(config_dir,L2_model_name,var_name)

            # mean_vertical_difference = 0
            # subset_copy = np.copy(var_grid)

            L1_wet_cells = np.copy(L1_Hfac)
            L1_wet_cells[L1_wet_cells > 0] = 1

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            print('    var_grid shape', np.shape(var_grid))
            print('    L1_wet_cells shape: ',np.shape(L1_wet_cells))
            print('    domain_wet_cells shape: ', np.shape(domain_wet_cells))

            L1_wet_cells_on_domain = np.copy(domain_wet_cells)


            interp_field = df.downscale_3D_seaice_field(L1_XC, L1_YC,
                                                 var_grid, L1_wet_cells,
                                                 L1_wet_cells_on_domain,
                                                 L2_XC, L2_YC, domain_wet_cells,
                                                 mean_vertical_difference=0, fill_downward=True,
                                                 printing=True)

            if np.sum(np.isnan(interp_field)) > 0:
                print('Setting ' + str(np.sum(np.isnan(interp_field))) + ' values to 0 in this grid')
                interp_field[np.isnan(interp_field)] = 0

            # plt.subplot(2, 3, 1)
            # plt.imshow(L1_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(L1_wet_cells_on_domain[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 2, 3)
            # plt.imshow(var_grid[0, :, :], origin='lower')
            # plt.subplot(2, 2, 4)
            # plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.show()

            # interp_field[domain_wet_cells_3D[:np.shape(interp_field)[0], :, :] == 0] = 0
        else:
            interp_field = np.zeros((1, np.shape(L2_XC)[0], np.shape(L2_XC)[1]))

        interp_grids.append(interp_field)

    # old_grids = [np.copy(interp_grids)[0],np.copy(interp_grids)[1]]

    interp_grids = rotate_interpolated_grids_to_domain(var_names, interp_grids, L2_AngleCS, L2_AngleSN)

    # # plt.subplot(2,2,1)
    # # plt.imshow(old_grids[0][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # # plt.subplot(2, 2, 2)
    # # plt.imshow(old_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # # plt.subplot(2, 2, 3)
    # # plt.imshow(interp_grids[0][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # # plt.subplot(2, 2, 4)
    # # plt.imshow(interp_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # # plt.show()

    pickup_grid = stack_grids_to_pickup(interp_grids)

    output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.' + '{:010d}'.format(L1_iteration*5))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [int(1)]
    pickup_metadata['nFlds'] = [6]
    pickup_metadata['nrecords'] = [6]
    pickup_metadata['fldlist'] = var_names
    write_seaice_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)



